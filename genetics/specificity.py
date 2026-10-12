#!/usr/bin/env python
"""Day 3, A7 item 4: cell-type specificity covariates for MAGMA gene-property analysis, from the SEA-AD pseudobulk objects.

For each of the 13 pseudobulk objects (Immune = Microglia-PVM, Astrocyte, Oligodendrocyte, OPC, Sst, Pvalb, Vip, Lamp5, L23-IT, L4-IT, L5-IT,
Endothelial, VLMC-Perivascular): MTG rows of aged donors (Neurotypical reference False), raw counts (`.X` holds integer counts) summed over rows
and donors, CPM over all 36,601 genes. Genes with total CPM >= 1 over the 13 types; specificity = CPM_t / sum_t CPM; covariate = rank-based
inverse normal transform of the specificity across genes, per cell type. Genes mapped to NCBI Gene IDs (genemap.py, A3 rule) and kept if they are in
NCBI37.3.gene.loc; if several rows map to one Gene ID the row with the larger total CPM is kept.
"""
import glob, json, os, sys, time
import numpy as np, pandas as pd, scipy.stats as st, scipy.sparse as sp, anndata as ad
sys.path.insert(0, "/home/pranavak/orcd/pool/mantis-films/research/broken-bricks/genetics")
from genemap import GeneMap

R = "/home/pranavak/orcd/pool/mantis-films/research/seaad"
OUT = f"{R}/results/day3"; os.makedirs(OUT, exist_ok=True)
TYPES = {"Immune": "Microglia", "Astrocyte": "Astrocyte", "Oligodendrocyte": "Oligodendrocyte", "OPC": "OPC", "Sst": "Sst", "Pvalb": "Pvalb",
         "Vip": "Vip", "Lamp5": "Lamp5", "L23-IT": "L23IT", "L4-IT": "L4IT", "L5-IT": "L5IT", "Endothelial": "Endothelial", "VLMC-Perivascular": "VLMC"}
T0 = time.perf_counter()
sums = {}; nrows = {}; nuc = {}; var0 = None
for key, name in TYPES.items():
    f = glob.glob(f"{R}/data/s3/Multiregion_2026/pseudobulk_objects/SEAAD_{key}_HIP_MEC_*pseudobulked.2026-06-22.h5ad")[0]
    a = ad.read_h5ad(f, backed="r")
    obs = a.obs[["Brain Region", "Neurotypical reference", "Number of nuclei"]]
    m = ((obs["Brain Region"].astype(str) == "MTG") & (obs["Neurotypical reference"].astype(str) == "False")).values
    idx = np.where(m)[0]
    v = a.var[["gene_ids"]].copy(); v["symbol"] = a.var.index.values
    if var0 is None: var0 = v
    else: assert (var0["gene_ids"].values == v["gene_ids"].values).all(), key
    tot = np.zeros(a.n_vars)
    for i in range(0, len(idx), 400):
        blk = idx[i:i + 400]
        xb = a.X[blk[0]:blk[-1] + 1]
        xb = xb.toarray() if sp.issparse(xb) else np.asarray(xb)
        tot += xb[blk - blk[0]].sum(axis=0)
    a.file.close()
    sums[name] = tot; nrows[name] = int(len(idx)); nuc[name] = int(obs["Number of nuclei"].values[m].sum())
    print(f"{name}: MTG aged rows {len(idx)} nuclei {nuc[name]} total counts {tot.sum():.3e}  {time.perf_counter() - T0:.0f}s", flush=True)

names = list(TYPES.values())
C = pd.DataFrame(sums, index=var0["gene_ids"].values)[names]
cpm = C.div(C.sum(axis=0), axis=1) * 1e6
tot = cpm.sum(axis=1)
keep = tot >= 1
spec = cpm[keep].div(tot[keep], axis=0)
meta = var0.loc[keep.values].copy()
meta["total_cpm"] = tot[keep].values
gm = GeneMap()
meta = gm.map_table(meta.rename(columns={"gene_ids": "ensembl_id"}), "ensembl_id", "symbol")
meta.index = spec.index
ok = meta["in_magma_loc"].values
spec, meta = spec[ok], meta[ok]
# one row per Gene ID: larger total CPM
meta["spec_row"] = np.arange(len(meta))
best = meta.sort_values("total_cpm", ascending=False).drop_duplicates("entrez")
spec, meta = spec.loc[best.index], best
n = len(spec)
cov = pd.DataFrame({c: st.norm.ppf((st.rankdata(spec[c].values) - 0.5) / n) for c in names})
cov.insert(0, "GENE", meta["entrez"].astype(int).values)
cov.to_csv(f"{OUT}/spec_covar.txt", sep="\t", index=False)
out = pd.concat([meta[["ensembl_id", "symbol", "entrez", "total_cpm"]].reset_index(drop=True), spec.reset_index(drop=True).add_prefix("spec_")], axis=1)
out.to_csv(f"{OUT}/specificity_table.csv", index=False)
cpm_top = {c: spec[c].sort_values(ascending=False).head(8).index.map(dict(zip(spec.index, meta["symbol"]))).tolist() for c in ["Microglia", "Oligodendrocyte", "Sst"]}
json.dump({"genes_in_covariate_file": int(n), "rows_mtg_aged": nrows, "nuclei_mtg_aged": nuc, "top_specific_genes_check": cpm_top,
           "seconds": round(time.perf_counter() - T0, 1)}, open(f"{OUT}/specificity_meta.json", "w"), indent=1)
print("genes in covariate file", n, "sanity top specific:", cpm_top, flush=True)
