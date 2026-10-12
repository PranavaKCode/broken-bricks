#!/usr/bin/env python
"""Stage A: corpus flags, PubTator3 gene-paper pairs with abstract-level verification (amendment A2), attention tables,
pilot reproduction (preregistration D0, D1', S1, S2, S4, S8, V1).

Reads : raw/pubmed/abstracts.parquet, raw/pubtator/gene2pubtator3_AD.tsv (awk-filtered), raw/ncbi/*gene_info.gz, raw/ncbi/gene2pubmed.gz
Writes: work/corpus.parquet, work/human_genes.parquet, work/pairs_human.parquet, work/gene_attention.parquet, work/A_gene_year*.parquet, out/results/*.csv
"""
import os, re, gzip, io, json, sys, datetime, collections, unicodedata
from pathlib import Path
import numpy as np
import polars as pl

LIT = Path(os.path.expanduser("~/orcd/pool/mantis-films/research/lit"))
RAW, WORK, RES = LIT / "raw", LIT / "work", LIT / "out" / "results"
WORK.mkdir(exist_ok=True, parents=True); RES.mkdir(exist_ok=True, parents=True)
LAST_FULL_YEAR = 2025

# ---- D0 corpus ---------------------------------------------------------------------------------
corp = pl.read_parquet(RAW / "pubmed" / "abstracts.parquet").sort("pmid")
AMY = re.compile(r"\bamyloid", re.I)
AMY2 = re.compile(r"\bamyloid|\bA[bB]eta\b|\bA-?beta\b|\bAβ|\bbeta-?amyloid", re.I)
TREM = re.compile(r"\bTREM-?2\b|triggering receptor expressed on myeloid cells 2", re.I)
texts = (corp["title"].fill_null("") + " " + corp["abstract"].fill_null("")).to_list()
corp = corp.with_columns(
    pl.Series("amyloid", [bool(AMY.search(t)) for t in texts]),
    pl.Series("amyloid_broad", [bool(AMY2.search(t)) for t in texts]),
    pl.Series("trem2_flag", [bool(TREM.search(t)) for t in texts]),
)
corp = corp.with_columns(
    pl.when(pl.col("year_epub").is_not_null() & (pl.col("year_epub") < pl.col("year"))).then(pl.col("year_epub")).otherwise(pl.col("year")).alias("year_earliest"))
corp.write_parquet(WORK / "corpus.parquet", compression="zstd")
print("corpus", corp.shape, "years", corp["year"].min(), corp["year"].max(), "null year", corp["year"].null_count(), flush=True)


# ---- genes -----------------------------------------------------------------------------------
def read_gene_info(path):
    df = pl.read_csv(path, separator="\t", comment_prefix=None, null_values=["-"], quote_char=None,
                     schema_overrides={"#tax_id": pl.Int64, "GeneID": pl.Int64, "chromosome": pl.Utf8, "Modification_date": pl.Utf8}, infer_schema_length=10000)
    return df.rename({"#tax_id": "tax_id"})

hum = read_gene_info(str(RAW / "ncbi" / "Homo_sapiens.gene_info.gz")).filter(pl.col("tax_id") == 9606)
hum = hum.with_columns(
    pl.col("dbXrefs").fill_null("").str.extract(r"Ensembl:(ENSG\d+)", 1).alias("ensembl_id"),
    pl.col("dbXrefs").fill_null("").str.extract(r"HGNC:HGNC:(\d+)", 1).alias("hgnc_id"),
).select(["GeneID", "Symbol", "Synonyms", "ensembl_id", "hgnc_id", "type_of_gene", "description"]).rename({"GeneID": "gene_id", "Symbol": "symbol"})
hum.write_parquet(WORK / "human_genes.parquet")
print("human genes", hum.shape, hum["type_of_gene"].value_counts().sort("count", descending=True).head(6).to_dicts(), flush=True)
human_ids = set(hum["gene_id"].to_list())

# ---- PubTator3 pairs ---------------------------------------------------------------------------------
pt = pl.read_csv(RAW / "pubtator" / "gene2pubtator3_AD.tsv", separator="\t", has_header=False, quote_char=None,
                 new_columns=["pmid", "type", "concept", "mentions", "resource"],
                 schema_overrides={"pmid": pl.Int64, "type": pl.Utf8, "concept": pl.Utf8, "mentions": pl.Utf8, "resource": pl.Utf8})
print("pubtator AD rows", pt.shape, flush=True)
pt = pt.filter(pl.col("pmid").is_in(corp["pmid"].to_list()))
pt = pt.with_columns(pl.col("concept").str.split(";").alias("gid")).explode("gid")
pt = pt.with_columns(pl.col("gid").cast(pl.Int64, strict=False).alias("gene_id")).drop_nulls("gene_id")
pairs_all = pt.group_by(["pmid", "gene_id"]).agg(pl.col("mentions").fill_null("").str.join("|").alias("mentions"), pl.col("resource").fill_null("").str.join("||").alias("resource"))
print("pairs: all species", pairs_all.height, flush=True)

# ---- amendment A2: abstract-level verification of each (PMID, gene) pair ---------------------------------
tx = dict(zip(corp["pmid"].to_list(), (t.lower() for t in texts)))
GREEK = str.maketrans({"α": "alpha", "β": "beta", "γ": "gamma", "δ": "delta", "ε": "epsilon", "κ": "kappa", "ω": "omega", "μ": "mu", "λ": "lambda",
                       "σ": "sigma", "θ": "theta", "ζ": "zeta", "η": "eta", "τ": "tau", "ρ": "rho", "π": "pi"})
NONALNUM = re.compile(r"[^0-9a-z]+")
def norm(s):
    return " " + NONALNUM.sub(" ", unicodedata.normalize("NFKC", s).lower().translate(GREEK)).strip() + " "
ntx = {}
def get_norm(pmid):
    v = ntx.get(pmid)
    if v is None:
        n = norm(tx.get(pmid, ""))
        v = (n, n.replace(" ", ""))
        ntx[pmid] = v
    return v
def verified_strict(pmid, mentions):      # A1: case-insensitive whole-token substring
    t = tx.get(pmid, "")
    for m in set(mentions.split("|")):
        m = m.strip().lower()
        if not m: continue
        i = t.find(m)
        while i >= 0:
            a_ = t[i - 1] if i > 0 else " "; b_ = t[i + len(m)] if i + len(m) < len(t) else " "
            if not a_.isalnum() and not b_.isalnum(): return True
            i = t.find(m, i + 1)
    return False
def verified_norm(pmid, mentions):        # A2: normalised, token boundaries, plus a space-free pass for strings >= 5 characters
    nt, nt_ns = get_norm(pmid)
    for m in set(mentions.split("|")):
        m = m.strip()
        if not m: continue
        nm = norm(m)
        if nm.strip() and nm in nt: return True
        ns = nm.replace(" ", "")
        if len(ns) >= 5 and ns in nt_ns: return True
    return False
pa_ = pairs_all.to_dict(as_series=False)
v_norm = [verified_norm(p, m) for p, m in zip(pa_["pmid"], pa_["mentions"])]
v_strict = [verified_strict(p, m) for p, m in zip(pa_["pmid"], pa_["mentions"])]
pairs_all = pairs_all.with_columns(pl.Series("abs_verified", v_norm), pl.Series("abs_verified_strict", v_strict)).join(corp.select(["pmid", "year"]), on="pmid", how="left")
pairs_h = pairs_all.filter(pl.col("gene_id").is_in(list(human_ids)))
print("human pairs", pairs_h.height, "A2-verified share", float(pairs_h["abs_verified"].mean()), "strict share", float(pairs_h["abs_verified_strict"].mean()), flush=True)
pairs_h.write_parquet(WORK / "pairs_human.parquet", compression="zstd")
ver = pairs_h["abs_verified"].to_list()
tm = pairs_h.filter(pl.col("mentions").str.strip_chars("|").str.len_chars() > 0).with_columns((pl.col("year") // 5 * 5).alias("year5"))
rec = tm.group_by("year5").agg(pl.len().alias("n_pairs_with_mention_text"), pl.col("abs_verified").mean().alias("share_verified_A2"), pl.col("abs_verified_strict").mean().alias("share_verified_A1_strict")).sort("year5")
rec.write_csv(RES / "verifier_share_by_period.csv"); print(rec.filter(pl.col("year5") >= 1975), flush=True)

# ---- S4 ortholog-by-symbol merge (mouse, rat Gene IDs whose upper-cased symbol equals a human symbol), abstract-verified only ----
sym2gid = {s.upper(): g for s, g in zip(hum["symbol"].to_list(), hum["gene_id"].to_list())}
orth = []
for sp, tx_id in (("Mus_musculus", 10090), ("Rattus_norvegicus", 10116)):
    gi = read_gene_info(str(RAW / "ncbi" / f"{sp}.gene_info.gz")).filter(pl.col("tax_id") == tx_id)
    gi = gi.select(["GeneID", "Symbol"]).with_columns(pl.col("Symbol").str.to_uppercase().alias("U"))
    gi = gi.filter(pl.col("U").is_in(list(sym2gid.keys())))
    orth.append(pl.DataFrame(dict(gene_id=gi["GeneID"].to_list(), human_gene_id=[sym2gid[u] for u in gi["U"].to_list()])))
orth = pl.concat(orth)
nonh = pairs_all.filter(pl.col("abs_verified") & ~pl.col("gene_id").is_in(list(human_ids))).join(orth, on="gene_id", how="inner")
s4 = pl.concat([pairs_h.filter(pl.col("abs_verified")).select(["pmid", "gene_id"]), nonh.select(["pmid", pl.col("human_gene_id").alias("gene_id")])]).unique()
s4 = s4.join(corp.select(["pmid", "year"]), on="pmid")
print("S4 ortholog pairs", s4.height, "(human abstract-verified only", int(pairs_h["abs_verified"].sum()), ")", flush=True)

# ---- S1 alternative attention definitions --------------------------------------------------------
ids_ad = set(corp["pmid"].to_list())
g2p = collections.Counter(); g2p_year = collections.Counter()
yr = dict(zip(corp["pmid"].to_list(), corp["year"].to_list()))
with gzip.open(RAW / "ncbi" / "gene2pubmed.gz", "rt") as fh:
    next(fh)
    for line in fh:
        if not line.startswith("9606\t"):
            continue
        t, g, p = line.rstrip("\n").split("\t")[:3]
        p = int(p)
        if p in ids_ad:
            g2p[int(g)] += 1
            g2p_year[(int(g), yr[p])] += 1
print("gene2pubmed AD pairs", sum(g2p.values()), flush=True)
tok = re.compile(r"[A-Za-z0-9]+(?:-[A-Za-z0-9]+)*")
symset = {s.lower(): g for s, g in zip(hum["symbol"].to_list(), hum["gene_id"].to_list())}
tiab = collections.Counter(); tiab_year = collections.Counter()
for p, t in zip(corp["pmid"].to_list(), texts):
    seen = {symset[w] for w in (x.lower() for x in tok.findall(t)) if w in symset}
    y = yr[p]
    for g in seen:
        tiab[g] += 1; tiab_year[(g, y)] += 1
print("tiab symbol hits", sum(tiab.values()), flush=True)

# ---- gene attention tables ---------------------------------------------------------------------------
pv = pairs_h.filter(pl.col("abs_verified"))
pv.group_by(["gene_id", "year"]).agg(pl.len().alias("n")).sort(["gene_id", "year"]).write_parquet(WORK / "A_gene_year.parquet")          # primary (A2)
pairs_h.group_by(["gene_id", "year"]).agg(pl.len().alias("n")).write_parquet(WORK / "A_gene_year_literal.parquet")                           # literal D1
pairs_h.filter(pl.col("abs_verified_strict")).group_by(["gene_id", "year"]).agg(pl.len().alias("n")).write_parquet(WORK / "A_gene_year_strict.parquet")
s4.group_by(["gene_id", "year"]).agg(pl.len().alias("n")).write_parquet(WORK / "A_gene_year_s4.parquet")
tot = pv.group_by("gene_id").agg(pl.len().alias("A_g"))
tot_lit = pairs_h.group_by("gene_id").agg(pl.len().alias("A_g_literal"))
tot_str = pairs_h.filter(pl.col("abs_verified_strict")).group_by("gene_id").agg(pl.len().alias("A_g_strict"))
tot_post = pv.filter((pl.col("year") >= 2010) & (pl.col("year") <= LAST_FULL_YEAR)).group_by("gene_id").agg(pl.len().alias("A_2010_2025"))
tot_s4 = s4.group_by("gene_id").agg(pl.len().alias("A_g_s4"))
ga = hum.join(tot, on="gene_id", how="left").join(tot_lit, on="gene_id", how="left").join(tot_str, on="gene_id", how="left").join(tot_post, on="gene_id", how="left").join(tot_s4, on="gene_id", how="left")
ga = ga.with_columns(
    pl.col("gene_id").map_elements(lambda g: g2p.get(g, 0), return_dtype=pl.Int64).alias("A_g_gene2pubmed"),
    pl.col("gene_id").map_elements(lambda g: tiab.get(g, 0), return_dtype=pl.Int64).alias("A_g_tiab"),
).with_columns([pl.col(c).fill_null(0) for c in ["A_g", "A_g_literal", "A_g_strict", "A_2010_2025", "A_g_s4"]])
ga.write_parquet(WORK / "gene_attention.parquet")
pl.DataFrame([dict(gene_id=g, year=y, n=n) for (g, y), n in g2p_year.items()]).write_parquet(WORK / "A_gene_year_gene2pubmed.parquet")
pl.DataFrame([dict(gene_id=g, year=y, n=n) for (g, y), n in tiab_year.items()]).write_parquet(WORK / "A_gene_year_tiab.parquet")
print(ga.sort("A_g", descending=True).head(15).select(["symbol", "A_g", "A_g_literal", "A_g_strict", "A_g_s4", "A_g_gene2pubmed", "A_g_tiab"]), flush=True)

# ---- V1b: TREM2 per year under different counting rules (diagnostic; does not replace L-V2) ---------------------------
tre = ga.filter(pl.col("symbol") == "TREM2")["gene_id"][0]
def peryear(path, name):
    return pl.read_parquet(path).filter(pl.col("gene_id") == tre).select(["year", "n"]).rename({"n": name})
vt = pl.DataFrame(dict(year=list(range(2008, 2018))))
vt = (vt.join(peryear(WORK / "A_gene_year.parquet", "trem2_pubtator_abstract_verified"), on="year", how="left")
        .join(peryear(WORK / "A_gene_year_literal.parquet", "trem2_pubtator_literal"), on="year", how="left")
        .join(corp.group_by("year").agg(pl.col("trem2_flag").sum().alias("trem2_regex_tiab_pubdate_year")), on="year", how="left")
        .join(corp.group_by("year_earliest").agg(pl.col("trem2_flag").sum().alias("trem2_regex_tiab_earliest_year")).rename({"year_earliest": "year"}), on="year", how="left")
        .fill_null(0).sort("year"))
vt.write_csv(RES / "v1b_trem2_by_year_variants.csv"); print(vt, flush=True)

# ---- validation V1 (amyloid shares, TREM2 counts) -------------------------------------------------------
by = corp.group_by("year").agg(pl.len().alias("n_papers"), pl.col("amyloid").sum().alias("n_amyloid"),
                               pl.col("amyloid_broad").sum().alias("n_amyloid_broad"), pl.col("trem2_flag").sum().alias("n_trem2_flag")).sort("year")
by = by.with_columns((pl.col("n_amyloid") / pl.col("n_papers")).alias("amyloid_share"), (pl.col("n_amyloid_broad") / pl.col("n_papers")).alias("amyloid_broad_share"))
by.write_csv(RES / "v1_corpus_by_year.csv")
pilot = {1990: 0.138, 2000: 0.245, 2010: 0.334, 2020: 0.329, 2025: 0.310}
rows = []
for y, p in pilot.items():
    r = by.filter(pl.col("year") == y).to_dicts()[0]
    rows.append(dict(check=f"amyloid_share_{y}", pilot=str(p), observed=float(r["amyloid_share"]), observed_broad=float(r["amyloid_broad_share"]), n_papers=int(r["n_papers"]),
                     abs_diff_pp=100 * (r["amyloid_share"] - p), passed=abs(r["amyloid_share"] - p) <= 0.03))
pool = by.filter((pl.col("year") >= 2010) & (pl.col("year") <= 2025))
pooled = pool["n_amyloid"].sum() / pool["n_papers"].sum()
rows.append(dict(check="amyloid_share_pooled_2010_2025", pilot="0.30-0.35", observed=float(pooled), observed_broad=float(pool["n_amyloid_broad"].sum() / pool["n_papers"].sum()),
                 n_papers=int(pool["n_papers"].sum()), abs_diff_pp=None, passed=0.30 <= pooled <= 0.35))
t12 = int(by.filter(pl.col("year") == 2012)["n_trem2_flag"].item()); t13 = int(by.filter(pl.col("year") == 2013)["n_trem2_flag"].item())
rows.append(dict(check="trem2_flag_2012 (pass: 2-8)", pilot="4", observed=float(t12), observed_broad=None, n_papers=None, abs_diff_pp=None, passed=2 <= t12 <= 8))
rows.append(dict(check="trem2_flag_2013 (pass: 30-50)", pilot="39", observed=float(t13), observed_broad=None, n_papers=None, abs_diff_pp=None, passed=30 <= t13 <= 50))
rows.append(dict(check="trem2_ratio_2013_over_2012 (pass: >=5)", pilot="9.75", observed=float(t13 / max(t12, 1)), observed_broad=None, n_papers=None, abs_diff_pp=None, passed=(t13 / max(t12, 1)) >= 5))
v1 = pl.DataFrame(rows, infer_schema_length=None)
v1.write_csv(RES / "v1_pilot_reproduction.csv"); print(v1, flush=True)
pl.DataFrame([dict(item=k, value=str(v)) for k, v in dict(
    n_corpus=corp.height, n_pairs_human=pairs_h.height, a2_verified_share=float(np.mean(ver)), a1_strict_share=float(pairs_h["abs_verified_strict"].mean()), n_pairs_s4=s4.height,
    n_genes_with_A=int((ga["A_g"] > 0).sum()), retrieved=datetime.datetime.now(datetime.timezone.utc).isoformat()).items()]).write_csv(RES / "stageA_summary.csv")
print("STAGE_A_DONE", flush=True)
