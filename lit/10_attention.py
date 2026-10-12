#!/usr/bin/env python
"""Stage A: corpus flags, PubTator3 gene-paper pairs, attention tables (preregistration D0, D1, S1, S2, S4, S8).

Reads : raw/pubmed/abstracts.parquet, raw/pubtator/gene2pubtator3_AD.tsv (awk-filtered), raw/ncbi/*gene_info.gz,
        raw/ncbi/gene2pubmed.gz
Writes: work/corpus.parquet, work/pairs_human.parquet, work/gene_attention.parquet, results/*.csv
"""
import os, re, gzip, io, json, sys, datetime, collections
from pathlib import Path
import numpy as np
import polars as pl

LIT = Path(os.path.expanduser("~/orcd/pool/mantis-films/research/lit"))
RAW, WORK, RES = LIT / "raw", LIT / "work", LIT / "out" / "results"
WORK.mkdir(exist_ok=True, parents=True); RES.mkdir(exist_ok=True, parents=True)
LAST_FULL_YEAR = 2025

# ---- D0 corpus ---------------------------------------------------------------------------------
corp = pl.read_parquet(RAW / "pubmed" / "abstracts.parquet")
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
                     schema_overrides={"#tax_id": pl.Int64, "GeneID": pl.Int64, "chromosome": pl.Utf8, "Modification_date": pl.Utf8},
                     infer_schema_length=0 if False else 10000)
    return df.rename({"#tax_id": "tax_id"})

hum = read_gene_info(str(RAW / "ncbi" / "Homo_sapiens.gene_info.gz")).filter(pl.col("tax_id") == 9606)
hum = hum.with_columns(
    pl.col("dbXrefs").fill_null("").str.extract(r"Ensembl:(ENSG\d+)", 1).alias("ensembl_id"),
    pl.col("dbXrefs").fill_null("").str.extract(r"HGNC:HGNC:(\d+)", 1).alias("hgnc_id"),
).select(["GeneID", "Symbol", "Synonyms", "ensembl_id", "hgnc_id", "type_of_gene", "description"]).rename({"GeneID": "gene_id", "Symbol": "symbol"})
hum.write_parquet(WORK / "human_genes.parquet")
print("human genes", hum.shape, hum["type_of_gene"].value_counts().sort("count", descending=True).head(6).to_dicts(), flush=True)
human_ids = set(hum["gene_id"].to_list())

# ---- D1 pairs (PubTator3, literal) -------------------------------------------------------------
pt = pl.read_csv(RAW / "pubtator" / "gene2pubtator3_AD.tsv", separator="\t", has_header=False, quote_char=None,
                 new_columns=["pmid", "type", "concept", "mentions", "resource"],
                 schema_overrides={"pmid": pl.Int64, "type": pl.Utf8, "concept": pl.Utf8, "mentions": pl.Utf8, "resource": pl.Utf8})
print("pubtator AD rows", pt.shape, pt["type"].value_counts().to_dicts(), pt["resource"].value_counts().to_dicts(), flush=True)
print("sample concepts with ';':", pt.filter(pl.col("concept").str.contains(";")).head(3).to_dicts(), flush=True)
print("sample concepts non-numeric:", pt.filter(~pl.col("concept").str.contains(r"^[0-9;]+$")).head(5).to_dicts(), flush=True)
pt = pt.filter(pl.col("pmid").is_in(corp["pmid"]))
pt = pt.with_columns(pl.col("concept").str.split(";").alias("gid")).explode("gid")
pt = pt.with_columns(pl.col("gid").cast(pl.Int64, strict=False).alias("gene_id")).drop_nulls("gene_id")
pairs_all = pt.group_by(["pmid", "gene_id"]).agg(pl.col("mentions").str.join("|").alias("mentions"))
pairs_h = pairs_all.filter(pl.col("gene_id").is_in(list(human_ids)))
print("pairs: all", pairs_all.height, "human", pairs_h.height, flush=True)

# S8: mention string occurs in title/abstract (whole-token, case-insensitive)
tx = dict(zip(corp["pmid"].to_list(), (t.lower() for t in texts)))
def verified(pmid, mentions):
    t = tx.get(pmid, "")
    for m in set(mentions.split("|")):
        m = m.strip().lower()
        if not m:
            continue
        i = t.find(m)
        while i >= 0:
            a = t[i - 1] if i > 0 else " "
            b = t[i + len(m)] if i + len(m) < len(t) else " "
            if not a.isalnum() and not b.isalnum():
                return True
            i = t.find(m, i + 1)
    return False
ph = pairs_h.to_dict(as_series=False)
ver = [verified(p, m) for p, m in zip(ph["pmid"], ph["mentions"])]
pairs_h = pairs_h.with_columns(pl.Series("abs_verified", ver))
pairs_h = pairs_h.join(corp.select(["pmid", "year"]), on="pmid", how="left")
pairs_h.write_parquet(WORK / "pairs_human.parquet", compression="zstd")
print("S8 verified share", float(np.mean(ver)), flush=True)

# ---- S4 ortholog-by-symbol merge (mouse, rat Gene IDs whose upper-cased symbol equals a human symbol) ------
sym2gid = {s.upper(): g for s, g in zip(hum["symbol"].to_list(), hum["gene_id"].to_list())}
orth = []
for sp, tx_id in (("Mus_musculus", 10090), ("Rattus_norvegicus", 10116)):
    gi = read_gene_info(str(RAW / "ncbi" / f"{sp}.gene_info.gz")).filter(pl.col("tax_id") == tx_id)
    gi = gi.select(["GeneID", "Symbol"]).with_columns(pl.col("Symbol").str.to_uppercase().alias("U"))
    gi = gi.filter(pl.col("U").is_in(list(sym2gid.keys())))
    orth.append(pl.DataFrame(dict(gene_id=gi["GeneID"].to_list(), human_gene_id=[sym2gid[u] for u in gi["U"].to_list()])))
orth = pl.concat(orth)
nonh = pairs_all.filter(~pl.col("gene_id").is_in(list(human_ids))).join(orth, on="gene_id", how="inner")
s4 = pl.concat([pairs_h.select(["pmid", pl.col("gene_id")]), nonh.select(["pmid", pl.col("human_gene_id").alias("gene_id")])]).unique()
s4 = s4.join(corp.select(["pmid", "year"]), on="pmid")
print("S4 ortholog pairs", s4.height, "(human-only", pairs_h.height, ")", flush=True)

# ---- S1 alternative attention definitions --------------------------------------------------------
ids_ad = set(corp["pmid"].to_list())
g2p = collections.Counter()
g2p_year = collections.Counter()
yr = dict(zip(corp["pmid"].to_list(), corp["year"].to_list()))
with gzip.open(RAW / "ncbi" / "gene2pubmed.gz", "rt") as fh:
    next(fh)
    for line in fh:
        t, g, p = line.rstrip("\n").split("\t")[:3]
        if t == "9606":
            p = int(p)
            if p in ids_ad:
                g2p[int(g)] += 1
                g2p_year[(int(g), yr[p])] += 1
print("gene2pubmed AD pairs", sum(g2p.values()), flush=True)

# title/abstract symbol search (case-insensitive token match of official symbol)
tok = re.compile(r"[A-Za-z0-9]+(?:-[A-Za-z0-9]+)*")
symset = {s.lower(): g for s, g in zip(hum["symbol"].to_list(), hum["gene_id"].to_list())}
tiab = collections.Counter(); tiab_year = collections.Counter()
for p, t in zip(corp["pmid"].to_list(), texts):
    seen = {symset[w] for w in (x.lower() for x in tok.findall(t)) if w in symset}
    y = yr[p]
    for g in seen:
        tiab[g] += 1; tiab_year[(g, y)] += 1
print("tiab symbol hits", sum(tiab.values()), flush=True)

# ---- gene attention table ---------------------------------------------------------------------------
yrs = list(range(1975, 2027))
a_y = pairs_h.group_by(["gene_id", "year"]).agg(pl.len().alias("n")).sort(["gene_id", "year"])
a_y.write_parquet(WORK / "A_gene_year.parquet")
s4y = s4.group_by(["gene_id", "year"]).agg(pl.len().alias("n"))
s4y.write_parquet(WORK / "A_gene_year_s4.parquet")
tot = pairs_h.group_by("gene_id").agg(pl.len().alias("A_g"), pl.col("abs_verified").sum().alias("A_g_abs"))
tot_post = pairs_h.filter(pl.col("year") >= 2010).filter(pl.col("year") <= LAST_FULL_YEAR).group_by("gene_id").agg(pl.len().alias("A_2010_2025"))
tot_s4 = s4.group_by("gene_id").agg(pl.len().alias("A_g_s4"))
ga = hum.join(tot, on="gene_id", how="left").join(tot_post, on="gene_id", how="left").join(tot_s4, on="gene_id", how="left")
ga = ga.with_columns(
    pl.col("gene_id").map_elements(lambda g: g2p.get(g, 0), return_dtype=pl.Int64).alias("A_g_gene2pubmed"),
    pl.col("gene_id").map_elements(lambda g: tiab.get(g, 0), return_dtype=pl.Int64).alias("A_g_tiab"),
).with_columns([pl.col(c).fill_null(0) for c in ["A_g", "A_g_abs", "A_2010_2025", "A_g_s4"]])
ga.write_parquet(WORK / "gene_attention.parquet")
pl.DataFrame([dict(gene_id=g, year=y, n=n) for (g, y), n in g2p_year.items()]).write_parquet(WORK / "A_gene_year_gene2pubmed.parquet")
pl.DataFrame([dict(gene_id=g, year=y, n=n) for (g, y), n in tiab_year.items()]).write_parquet(WORK / "A_gene_year_tiab.parquet")
print(ga.sort("A_g", descending=True).head(15).select(["symbol", "A_g", "A_g_abs", "A_g_s4", "A_g_gene2pubmed", "A_g_tiab"]), flush=True)

# ---- validation V1 (amyloid shares, TREM2 counts) -------------------------------------------------------
by = corp.group_by("year").agg(pl.len().alias("n_papers"), pl.col("amyloid").sum().alias("n_amyloid"),
                               pl.col("amyloid_broad").sum().alias("n_amyloid_broad"), pl.col("trem2_flag").sum().alias("n_trem2_flag")).sort("year")
by = by.with_columns((pl.col("n_amyloid") / pl.col("n_papers")).alias("amyloid_share"), (pl.col("n_amyloid_broad") / pl.col("n_papers")).alias("amyloid_broad_share"))
by.write_csv(RES / "v1_corpus_by_year.csv")
pilot = {1990: 0.138, 2000: 0.245, 2010: 0.334, 2020: 0.329, 2025: 0.310}
rows = []
for y, p in pilot.items():
    r = by.filter(pl.col("year") == y).to_dicts()[0]
    rows.append(dict(check=f"amyloid_share_{y}", pilot=p, observed=r["amyloid_share"], observed_broad=r["amyloid_broad_share"], n_papers=r["n_papers"],
                     abs_diff_pp=100 * (r["amyloid_share"] - p), within_3pp=abs(r["amyloid_share"] - p) <= 0.03))
pool = by.filter((pl.col("year") >= 2010) & (pl.col("year") <= 2025))
pooled = pool["n_amyloid"].sum() / pool["n_papers"].sum()
rows.append(dict(check="amyloid_share_pooled_2010_2025", pilot="0.30-0.35", observed=pooled, observed_broad=pool["n_amyloid_broad"].sum() / pool["n_papers"].sum(),
                 n_papers=pool["n_papers"].sum(), abs_diff_pp=None, within_3pp=0.30 <= pooled <= 0.35))
t12 = by.filter(pl.col("year") == 2012)["n_trem2_flag"].item(); t13 = by.filter(pl.col("year") == 2013)["n_trem2_flag"].item()
rows.append(dict(check="trem2_flag_2012", pilot=4, observed=t12, observed_broad=None, n_papers=None, abs_diff_pp=None, within_3pp=2 <= t12 <= 8))
rows.append(dict(check="trem2_flag_2013", pilot=39, observed=t13, observed_broad=None, n_papers=None, abs_diff_pp=None, within_3pp=30 <= t13 <= 50))
rows.append(dict(check="trem2_ratio_2013_over_2012", pilot=39 / 4, observed=t13 / max(t12, 1), observed_broad=None, n_papers=None, abs_diff_pp=None, within_3pp=(t13 / max(t12, 1)) >= 5))
v1 = pl.DataFrame(rows, infer_schema_length=None)
v1.write_csv(RES / "v1_pilot_reproduction.csv")
print(v1, flush=True)
pl.DataFrame([dict(item=k, value=v) for k, v in dict(
    n_corpus=corp.height, n_pairs_human=pairs_h.height, s8_verified_share=float(np.mean(ver)), n_pairs_s4=s4.height,
    n_genes_with_A=int((ga["A_g"] > 0).sum()), retrieved=datetime.datetime.utcnow().isoformat() + "Z").items()]).write_csv(RES / "stageA_summary.csv")
print("STAGE_A_DONE", flush=True)
