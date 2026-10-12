#!/usr/bin/env python
"""Stage D: assemble per-paper and per-gene metadata with map coordinates, clusters, macro-topics (exploratory).
Writes work/papers_meta.parquet, work/genes_meta.parquet, out/results/macro_topics.csv, out/results/macro_topic_share_by_year.csv,
and the Mantis floor CSVs (floors/)."""
import os, re, json, datetime
from pathlib import Path
import numpy as np
import polars as pl

SEED = 20261011
LIT = Path(os.path.expanduser("~/orcd/pool/mantis-films/research/lit"))
RAW, WORK, RES = LIT / "raw", LIT / "work", LIT / "out" / "results"
EMB = WORK / "embed"
FL = LIT / "out" / "floors"; FL.mkdir(parents=True, exist_ok=True)
CAP = 150_000
def log(*a): print(*a, flush=True)

corp = pl.read_parquet(WORK / "corpus.parquet").sort("pmid")
ids = pl.read_parquet(EMB / "ids_union.parquet")
Z = np.load(EMB / "pca64.npy"); U = np.load(EMB / "umap2.npy"); lab = np.load(EMB / "leiden_labels.npy")
cl = pl.read_csv(EMB / "cluster_labels.csv")
npap = corp.height
assert Z.shape[0] == ids.height and (ids["id"].to_numpy()[:npap] == corp["pmid"].to_numpy()).all()
gtex = pl.read_parquet(EMB / "gene_texts.parquet")
gm = pl.read_csv(RES / "gene_master.csv", infer_schema_length=20000)
gm_all = gm

# ---- macro-topics: KMeans on cluster centroids (paper members only), weighted by cluster size --------------------
from sklearn.cluster import KMeans
nc = lab.max() + 1
Zp, lp = Z[:npap], lab[:npap]
cent = np.zeros((nc, 64)); cnt = np.bincount(lp, minlength=nc)
np.add.at(cent, lp, Zp); ok = cnt > 0; cent[ok] /= cnt[ok][:, None]
K = 8
km = KMeans(n_clusters=K, n_init=20, random_state=SEED).fit(cent[ok], sample_weight=cnt[ok])
macro_of_cluster = np.full(nc, -1); macro_of_cluster[np.where(ok)[0]] = km.labels_
macro = macro_of_cluster[lab]            # for every point (paper or gene)
from sklearn.feature_extraction.text import TfidfVectorizer
paper_text = (corp["title"].fill_null("") + " " + corp["abstract"].fill_null("")).to_list()
rng = np.random.default_rng(SEED)
docs = []
for m in range(K):
    mem = np.where(macro[:npap] == m)[0]
    take = rng.choice(mem, min(800, len(mem)), replace=False)
    docs.append(" ".join(paper_text[i] for i in take))
vec = TfidfVectorizer(stop_words="english", ngram_range=(1, 2), min_df=1, max_df=0.5, sublinear_tf=True, token_pattern=r"(?u)\b[A-Za-z][A-Za-z0-9\-]{2,}\b")
Mt = vec.fit_transform(docs); terms = np.array(vec.get_feature_names_out())
mac_rows = []
for m in range(K):
    top = terms[np.argsort(-Mt.getrow(m).toarray().ravel())[:6]]
    mac_rows.append(dict(macro=m, n_papers=int((macro[:npap] == m).sum()), n_clusters=int((macro_of_cluster == m).sum()), label=", ".join(top[:3]), keywords="; ".join(top)))
mac = pl.DataFrame(mac_rows)
mac.write_csv(RES / "macro_topics.csv"); log(mac)
yr = corp["year"].to_numpy()
rows = []
for y in range(1975, 2027):
    sel = yr == y
    if sel.sum() == 0: continue
    for m in range(K):
        rows.append(dict(year=y, macro=m, label=mac_rows[m]["label"], n_papers=int(sel.sum()), n=int(((macro[:npap] == m) & sel).sum()), share=float(((macro[:npap] == m) & sel).sum() / sel.sum())))
pl.DataFrame(rows).write_csv(RES / "macro_topic_share_by_year.csv")

# ---- papers meta ------------------------------------------------------------------------------------------------
pairs = pl.read_parquet(WORK / "pairs_human.parquet").select(["pmid", "gene_id"])
sym_of = dict(zip(gm_all["gene_id"].to_list(), gm_all["symbol"].to_list()))
G_of = dict(zip(gm_all["gene_id"].to_list(), gm_all["G_genetic"].to_list()))
dom_of = {g: (d.split("|") if d else []) for g, d in zip(gm_all["gene_id"].to_list(), gm_all["biodomains"].fill_null("").to_list())}
pg = pairs.group_by("pmid").agg(pl.col("gene_id").alias("gids"))
pgd = {p: g for p, g in zip(pg["pmid"].to_list(), pg["gids"].to_list())}
genes_s, ngenes, maxG, dom, topG = [], [], [], [], []
for p in corp["pmid"].to_list():
    gs = pgd.get(p, [])
    syms = sorted({sym_of.get(g, str(g)) for g in gs})
    genes_s.append(";".join(syms[:30])); ngenes.append(len(gs))
    gv = [(G_of.get(g, 0.0) or 0.0, sym_of.get(g, "")) for g in gs]
    maxG.append(max([v for v, _ in gv], default=0.0))
    topG.append(max(gv, default=(0.0, ""))[1] if gv and max(gv)[0] > 0 else "")
    w = {}
    for g in gs:
        d = dom_of.get(g, [])
        for x in d: w[x] = w.get(x, 0.0) + 1.0 / len(d)
    dom.append(max(sorted(w), key=lambda k: w[k]) if w else "Unassigned")
cl_label = dict(zip(cl["cluster"].to_list(), cl["label"].to_list()))
papers = corp.select(["pmid", "year", "journal", "title", "abstract", "amyloid", "record_type"]).with_columns(
    pl.Series("genes", genes_s), pl.Series("n_genes", ngenes), pl.Series("ot_genetic_max", maxG), pl.Series("top_genetic_gene", topG), pl.Series("biodomain", dom),
    pl.Series("cluster", lab[:npap]), pl.Series("cluster_label", [cl_label[c] for c in lab[:npap]]), pl.Series("macro", macro[:npap]), pl.Series("macro_label", [mac_rows[m]["label"] for m in macro[:npap]]),
    pl.Series("coordinate1", U[:npap, 0]), pl.Series("coordinate2", U[:npap, 1]))
for k in range(64): papers = papers.with_columns(pl.Series(f"v{k}", Z[:npap, k]))
papers.write_parquet(WORK / "papers_meta.parquet", compression="zstd")

# ---- genes meta -----------------------------------------------------------------------------------------------------
gi = {g: i for i, g in enumerate(ids["id"].to_numpy()[npap:])}
gt = gtex.join(gm_all.select(["gene_id", "symbol", "ensembl_id", "A_g", "A_2010_2025", "A_g_abs_verified", "G_genetic", "treatad_genetics_score", "y0_gwas", "E_g", "n_biodomains",
                              "biodomains", "strong_G07", "sparse_G07_A100"]), on="gene_id", how="left")
und = pl.read_csv(RES / "under_explored_G07_ranked.csv").select(["gene_id", pl.col("rank").alias("under_explored_rank")])
gt = gt.join(und, on="gene_id", how="left")
order = [gi[g] for g in gt["gene_id"].to_list()]
gt = gt.with_columns(pl.Series("cluster", lab[npap:][order]), pl.Series("cluster_label", [cl_label[c] for c in lab[npap:][order]]), pl.Series("macro", macro[npap:][order]),
                     pl.Series("coordinate1", U[npap:, 0][order]), pl.Series("coordinate2", U[npap:, 1][order]))
for k in range(64): gt = gt.with_columns(pl.Series(f"v{k}", Z[npap:, k][order]))
gt.write_parquet(WORK / "genes_meta.parquet", compression="zstd")
log("papers", papers.shape, "genes", gt.shape)

# ---- Mantis floors (CSV with coordinate1/coordinate2 + 64 vector columns) ---------------------------------------------------
vcols = [f"v{k}" for k in range(64)]
def rnd(df):
    return df.with_columns([pl.col(c).round(5) for c in vcols + ["coordinate1", "coordinate2"]])
# floor 3: genes
g3 = gt.with_columns(
    pl.col("symbol").alias("title"),
    (pl.col("protein_name").fill_null("") + ". " + pl.col("text").str.slice(0, 600)).alias("text_snippet"),
    pl.lit(3).alias("floor_number"),
    pl.col("biodomains").fill_null("").str.split("|").list.first().fill_null("").alias("biodomain_first"),
).select(["floor_number", "gene_id", "title", "symbol", "ensembl_id", "uniprot", "protein_name", "text_snippet", "text_source", "G_genetic", "treatad_genetics_score", "A_g", "A_2010_2025",
          "y0_gwas", "E_g", "n_biodomains", "biodomains", "biodomain_first", "strong_G07", "sparse_G07_A100", "under_explored_rank", "cluster", "cluster_label", "coordinate1", "coordinate2"] + vcols)
g3 = g3.rename({"G_genetic": "ot_genetic_score", "A_g": "ad_paper_count", "A_2010_2025": "ad_paper_count_2010_2025", "y0_gwas": "first_gwas_year", "biodomains": "biodomain", "biodomain_first": "biodomain_primary"})
rnd(g3).write_csv(FL / "floor3_genes.csv"); log("floor3", g3.shape)
# floor 4: literature (full and capped)
p4 = papers.with_columns(pl.lit(4).alias("floor_number"), pl.col("abstract").str.slice(0, 500).alias("abstract_snippet"),
                         ("https://pubmed.ncbi.nlm.nih.gov/" + pl.col("pmid").cast(pl.Utf8) + "/").alias("url")).select(
    ["floor_number", "pmid", "title", "abstract_snippet", "url", "year", "journal", "genes", "n_genes", "biodomain", "amyloid", "ot_genetic_max", "top_genetic_gene", "cluster", "cluster_label",
     "macro_label", "coordinate1", "coordinate2"] + vcols).rename({"amyloid": "amyloid_flag", "ot_genetic_max": "ot_genetic_score_max_gene"})
p4 = p4.with_columns(pl.col("amyloid_flag").cast(pl.Int8))
rnd(p4).write_csv(FL / "floor4_literature_full.csv"); log("floor4 full", p4.shape)
if p4.height > CAP:
    # stratified by year x amyloid flag, proportional allocation, fixed seed (preregistration 2.2)
    strata = p4.with_columns((pl.col("year").cast(pl.Utf8) + "_" + pl.col("amyloid_flag").cast(pl.Utf8)).alias("stratum"))
    rng = np.random.default_rng(SEED)
    keep = []
    frac = CAP / p4.height
    rows = []
    for s, df in strata.group_by("stratum"):
        n = len(df); k = int(round(n * frac))
        if n > 0 and k == 0: k = 1 if rng.random() < n * frac else 0
        idx = rng.choice(n, min(k, n), replace=False)
        keep.append(df["pmid"].to_numpy()[idx])
        rows.append(dict(stratum=s[0], n_full=n, n_sample=len(idx)))
    keep = np.concatenate(keep)
    # trim or top up to exactly CAP
    if len(keep) > CAP: keep = rng.choice(keep, CAP, replace=False)
    elif len(keep) < CAP:
        rest = np.setdiff1d(p4["pmid"].to_numpy(), keep); keep = np.concatenate([keep, rng.choice(rest, CAP - len(keep), replace=False)])
    p4c = p4.filter(pl.col("pmid").is_in(keep))
    pl.DataFrame(rows).sort("stratum").write_csv(RES / "floor_sampling.csv")
    rnd(p4c).write_csv(FL / "floor4_literature_150k.csv"); log("floor4 capped", p4c.shape)
    chk = pl.DataFrame([dict(item="full_n", value=p4.height), dict(item="sample_n", value=p4c.height), dict(item="seed", value=SEED),
                        dict(item="amyloid_share_full", value=float(p4["amyloid_flag"].mean())), dict(item="amyloid_share_sample", value=float(p4c["amyloid_flag"].mean())),
                        dict(item="year_min_sample", value=int(p4c["year"].min())), dict(item="year_max_sample", value=int(p4c["year"].max())),
                        dict(item="n_with_G_gt0_gene_full", value=int((p4["ot_genetic_score_max_gene"] > 0).sum())), dict(item="n_with_G_gt0_gene_sample", value=int((p4c["ot_genetic_score_max_gene"] > 0).sum()))])
    chk.write_csv(RES / "floor_sampling_check.csv")
log("FLOORS_DONE")
