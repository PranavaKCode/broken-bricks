#!/usr/bin/env python
"""Final check 4 (part 1): pooled-neuron cell-type specificity covariates for MAGMA gene-property analysis.

Same inputs and gene universe as genetics/specificity.py (13 SEA-AD pseudobulk objects, MTG rows of aged donors, raw counts summed over rows and
donors, genes with total CPM >= 1 over the 13 types, mapped to NCBI Gene IDs present in NCBI37.3.gene.loc, one row per Gene ID).
New: neuronal types are POOLED by summing their raw counts before CPM, giving groups instead of 13 separate types:
   Neuron_all = Sst + Pvalb + Vip + Lamp5 + L23IT + L4IT + L5IT  (all seven neuronal pseudobulk objects)
   Exc3 = L23IT + L4IT + L5IT          Inh4 = Sst + Pvalb + Vip + Lamp5
   Microglia, Astrocyte, Oligodendrocyte, OPC, Endothelial, VLMC (unchanged)
Two specificity families, each = CPM_group / sum of CPM over its own set of groups, then rank-based inverse normal transform across genes:
   N7 : Neuron_all + the 6 non-neuronal types          (covariates N7_Neuron_all, N7_Microglia, ...)
   N8 : Exc3 + Inh4 + the 6 non-neuronal types         (covariates N8_Exc3, N8_Inh4, N8_Microglia, ...)
Writes results/final_checks/spec_covar_pooled.txt (GENE + covariates) and spec_pooled_meta.json
"""
import glob, json, os, sys, time
import numpy as np, pandas as pd, scipy.stats as st, scipy.sparse as sp, anndata as ad
sys.path.insert(0, "/home/pranavak/orcd/pool/mantis-films/research/broken-bricks/genetics")
from genemap import GeneMap

R = "/home/pranavak/orcd/pool/mantis-films/research/seaad"
OUT = f"{R}/results/final_checks"; os.makedirs(OUT, exist_ok=True)
TYPES = {"Immune": "Microglia", "Astrocyte": "Astrocyte", "Oligodendrocyte": "Oligodendrocyte", "OPC": "OPC", "Sst": "Sst", "Pvalb": "Pvalb",
         "Vip": "Vip", "Lamp5": "Lamp5", "L23-IT": "L23IT", "L4-IT": "L4IT", "L5-IT": "L5IT", "Endothelial": "Endothelial", "VLMC-Perivascular": "VLMC"}
NEURONAL = ["Sst", "Pvalb", "Vip", "Lamp5", "L23IT", "L4IT", "L5IT"]
NONNEU = ["Microglia", "Astrocyte", "Oligodendrocyte", "OPC", "Endothelial", "VLMC"]
T0 = time.perf_counter()
sums = {}; nrows = {}; nuc = {}; var0 = None
for key, name in TYPES.items():
    f = glob.glob(f"{R}/data/s3/Multiregion_2026/pseudobulk_objects/SEAAD_{key}_HIP_MEC_*pseudobulked.2026-06-22.h5ad")[0]
    a = ad.read_h5ad(f, backed="r")
    obs = a.obs[["Brain Region", "Neurotypical reference", "Number of nuclei"]]
    m = ((obs["Brain Region"].astype(str) == "MTG") & (obs["Neurotypical reference"].astype(str) == "False")).values
    idx = np.where(m)[0]
    v = a.var[["gene_ids"]].copy(); v["symbol"] = a.var.index.values
    if var0 is None:
        var0 = v
    else:
        assert (var0["gene_ids"].values == v["gene_ids"].values).all(), key
    tot = np.zeros(a.n_vars)
    for i in range(0, len(idx), 400):
        blk = idx[i:i + 400]
        xb = a.X[blk[0]:blk[-1] + 1]
        xb = xb.toarray() if sp.issparse(xb) else np.asarray(xb)
        tot += xb[blk - blk[0]].sum(axis=0)
    a.file.close()
    sums[name] = tot; nrows[name] = int(len(idx)); nuc[name] = int(obs["Number of nuclei"].values[m].sum())
    print(f"{name}: MTG aged rows {len(idx)} nuclei {nuc[name]} total counts {tot.sum():.6e}  {time.perf_counter() - T0:.0f}s", flush=True)

C13 = pd.DataFrame(sums, index=var0["gene_ids"].values)[list(TYPES.values())]
cpm13 = C13.div(C13.sum(axis=0), axis=1) * 1e6
tot13 = cpm13.sum(axis=1)
keep = tot13 >= 1                      # identical gene filter to genetics/specificity.py


def pool(cols):
    return C13[cols].sum(axis=1)


G = pd.DataFrame({"Neuron_all": pool(NEURONAL), "Exc3": pool(["L23IT", "L4IT", "L5IT"]), "Inh4": pool(["Sst", "Pvalb", "Vip", "Lamp5"])})
for c in NONNEU:
    G[c] = C13[c]
Gcpm = G.div(G.sum(axis=0), axis=1) * 1e6

fam = {"N7": ["Neuron_all"] + NONNEU, "N8": ["Exc3", "Inh4"] + NONNEU}
spec = {}
for fn, cols in fam.items():
    s = Gcpm[cols][keep.values]
    s = s.div(s.sum(axis=1), axis=0)
    for c in cols:
        spec[f"{fn}_{c}"] = s[c]
spec = pd.DataFrame(spec)
meta = var0.loc[keep.values].copy()
meta["total_cpm"] = tot13[keep.values].values
gm = GeneMap()
meta = gm.map_table(meta.rename(columns={"gene_ids": "ensembl_id"}), "ensembl_id", "symbol")
meta.index = spec.index
ok = meta["in_magma_loc"].values
spec, meta = spec[ok], meta[ok]
meta["spec_row"] = np.arange(len(meta))
best = meta.sort_values("total_cpm", ascending=False).drop_duplicates("entrez")
spec, meta = spec.loc[best.index], best
n = len(spec)
cov = pd.DataFrame({c: st.norm.ppf((st.rankdata(spec[c].values) - 0.5) / n) for c in spec.columns})
cov.insert(0, "GENE", meta["entrez"].astype(int).values)
cov.to_csv(f"{OUT}/spec_covar_pooled.txt", sep="\t", index=False)

# gene-universe check against the original covariate file
orig = pd.read_csv(f"{R}/results/day3/spec_covar.txt", sep="\t")
same = set(orig["GENE"]) == set(cov["GENE"])
# sanity: classic neuronal and glial markers
sym = dict(zip(meta["entrez"].astype(int), meta["symbol"]))
cov_s = cov.set_index("GENE")
chk = {}
for g in ["SNAP25", "SYT1", "RBFOX3", "GAD1", "SLC17A7", "CSF1R", "P2RY12", "GFAP", "MBP", "PLP1"]:
    ids = [k for k, v in sym.items() if v == g]
    if ids:
        chk[g] = {c: round(float(cov_s.loc[ids[0], c]), 2) for c in ["N7_Neuron_all", "N7_Microglia", "N7_Astrocyte", "N7_Oligodendrocyte"]}
top = {c: [sym[g] for g in cov_s[c].sort_values(ascending=False).head(10).index] for c in ["N7_Neuron_all", "N7_Microglia", "N8_Exc3", "N8_Inh4"]}
json.dump({"genes_in_covariate_file": int(n), "same_gene_set_as_original": bool(same), "n_orig": int(len(orig)), "rows_mtg_aged": nrows, "nuclei_mtg_aged": nuc,
           "pooled_total_counts": {c: float(G[c].sum()) for c in G.columns}, "marker_check_INT": chk, "top10_specific": top,
           "seconds": round(time.perf_counter() - T0, 1)}, open(f"{OUT}/spec_pooled_meta.json", "w"), indent=1)
print("genes", n, "same gene set as original:", same, flush=True)
print(json.dumps({"marker_check_INT": chk, "top10": top}, indent=1), flush=True)
