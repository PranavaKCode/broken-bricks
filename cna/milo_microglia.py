#!/usr/bin/env python
"""Broken Bricks Day 1: Milo (pertpy) cross-check of the microglia CNA result (amendment A4.5).

Same cells, same graph recipe (k = 30 on X_scVI, random_state 0), design ~ age + sex + pmi + method_frac + CPS_Local
(continuous covariates z-scored across donors), tested coefficient = CPS_Local (last term), pydeseq2 solver.
Outputs milo_microglia_nhoods.csv, milo_microglia_cells.parquet, milo_microglia_summary.json
"""
import argparse, json, os, sys, time, warnings
os.environ.setdefault("LOKY_MAX_CPU_COUNT", "4")   # pydeseq2/joblib would otherwise start one worker per node core and OOM
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import scanpy as sc
import pertpy as pt
from scipy.stats import spearmanr, pearsonr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *

ap = argparse.ArgumentParser()
ap.add_argument("--indir", default=OUT); ap.add_argument("--outdir", default=OUT)
ap.add_argument("--prop", type=float, default=0.1)
args = ap.parse_args()
t0 = time.perf_counter()

a = ad.read_h5ad(f"{H5AD_DIR}/Microglia-PVM_MTG.h5ad", backed="r")
obs = a.obs.copy(); Z = np.asarray(a.obsm["X_scVI"], dtype=np.float32); names = np.asarray(a.obs_names); a.file.close()
donors = pd.read_csv(DONOR_TABLE, index_col="donor_id")
obs["donor_id"] = obs["donor_id"].astype(str)
keep = (obs["donor_id"].isin(donors.index) & donors["CPS_Local"].reindex(obs["donor_id"]).notna().values).values
obs = obs.loc[keep].copy(); Z = Z[keep]; names = names[keep]
ids = pd.Index(sorted(obs["donor_id"].unique()))
mf = obs["assay"].astype(str).str.contains("multiome").groupby(obs["donor_id"]).mean().reindex(ids)
d = pd.DataFrame({"age": donors["age"].reindex(ids), "sex": donors["sex_male"].reindex(ids).astype(float),
                  "pmi": donors["pmi"].reindex(ids), "method_frac": mf, "CPS_Local": donors["CPS_Local"].reindex(ids)}, index=ids)
for c in ["age", "pmi", "CPS_Local"]:
    d[c] = (d[c] - d[c].mean()) / d[c].std(ddof=1)
cell_cov = d.reindex(obs["donor_id"].values); cell_cov.index = pd.Index(names).astype(str)
adata = ad.AnnData(X=sp.csr_matrix((len(obs), 1), dtype=np.float32), obs=cell_cov.copy())
adata.obs["donor_id"] = obs["donor_id"].values
adata.obs["Supertype"] = obs["Supertype"].astype(str).values
adata.obsm["X_scVI"] = Z
sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors=NEIGHBORS, random_state=0)
print("cells", adata.n_obs, "donors", len(ids), flush=True)

milo = pt.tl.Milo()
mdata = milo.load(adata)
milo.make_nhoods(mdata["rna"], prop=args.prop, seed=0)
mdata = milo.count_nhoods(mdata, sample_col="donor_id")
print("neighborhoods", mdata["milo"].n_vars, flush=True)
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    milo.da_nhoods(mdata, design="~ age + sex + pmi + method_frac + CPS_Local", solver="pydeseq2")
nh = mdata["milo"].var.copy()
nhoods = mdata["rna"].obsm["nhoods"].tocsc()           # cells x neighborhoods (binary)
n_cells_nh = np.asarray(nhoods.sum(axis=0)).ravel()
nh["n_cells"] = n_cells_nh

# CNA primary result for the same cells
cna_cells = pd.read_parquet(f"{args.indir}/microglia__all_cells.parquet").set_index("obs_name")
nc = cna_cells["ncorr_CPS_Local"].reindex(adata.obs_names).values.astype(float)
fd = cna_cells["fdr_CPS_Local"].reindex(adata.obs_names).values.astype(float)
ncv = np.nan_to_num(nc, nan=0.0)
sum_nc = nhoods.T.dot(ncv); cnt = np.asarray(nhoods.sum(axis=0)).ravel()
nh["cna_mean_ncorr"] = sum_nc / np.maximum(cnt, 1)
nh["cna_frac_fdr10"] = np.asarray(nhoods.T.dot(((fd <= 0.10) & ~np.isnan(nc)).astype(float))).ravel() / np.maximum(cnt, 1)

# cell-level Milo score = mean logFC over the neighborhoods containing the cell
lf = nh["logFC"].fillna(0).values
member = np.asarray(nhoods.sum(axis=1)).ravel()
cell_lfc = np.where(member > 0, nhoods.dot(lf) / np.maximum(member, 1), np.nan)
sig = (nh["SpatialFDR"] < 0.10).values
msig = np.asarray(nhoods[:, np.where(sig)[0]].sum(axis=1)).ravel() if sig.any() else np.zeros(adata.n_obs)
cell_lfc_sig = np.where(msig > 0, nhoods[:, np.where(sig)[0]].dot(lf[sig]) / np.maximum(msig, 1), np.nan) if sig.any() else np.full(adata.n_obs, np.nan)
cells = pd.DataFrame({"obs_name": adata.obs_names.values, "Supertype": adata.obs["Supertype"].values,
                      "milo_logFC_mean": cell_lfc.astype(np.float32), "milo_logFC_sig_only": cell_lfc_sig.astype(np.float32),
                      "n_nhoods": member.astype(np.int16), "ncorr_CPS_Local": nc.astype(np.float32)})

ok = nh["logFC"].notna()
rho_nh = spearmanr(nh.loc[ok, "logFC"], nh.loc[ok, "cna_mean_ncorr"])
r_nh = pearsonr(nh.loc[ok, "logFC"], nh.loc[ok, "cna_mean_ncorr"])
m = ~np.isnan(cell_lfc) & ~np.isnan(nc)
rho_cell = spearmanr(cell_lfc[m], nc[m])
S = nh[(nh["SpatialFDR"] < 0.10) & ok]
agree_sig = float((np.sign(S["logFC"]) == np.sign(S["cna_mean_ncorr"])).mean()) if len(S) else None
S5 = nh[(nh["SpatialFDR"] < 0.05) & ok]
cna_hit = nh[(nh["cna_frac_fdr10"] >= 0.5) & ok]
agree_cna = float((np.sign(cna_hit["logFC"]) == np.sign(cna_hit["cna_mean_ncorr"])).mean()) if len(cna_hit) else None
# supertype-level: mean Milo logFC vs mean ncorr
stm = cells.groupby("Supertype").agg(milo=("milo_logFC_mean", "mean"), cna=("ncorr_CPS_Local", "mean"), n=("obs_name", "size"))
summary = {
    "design": "~ age + sex + pmi + method_frac + CPS_Local (z-scored continuous), solver pydeseq2, prop=%.2f" % args.prop,
    "n_cells": int(adata.n_obs), "n_donors": int(len(ids)), "n_neighborhoods": int(len(nh)), "n_tested": int(ok.sum()),
    "median_cells_per_neighborhood": float(np.median(n_cells_nh)),
    "n_nhoods_spatialFDR_lt_0.10": int(len(S)), "n_pos": int((S["logFC"] > 0).sum()), "n_neg": int((S["logFC"] < 0).sum()),
    "n_nhoods_spatialFDR_lt_0.05": int(len(S5)),
    "spearman_nhood_logFC_vs_mean_ncorr": [float(rho_nh[0]), float(rho_nh[1])], "pearson_nhood": [float(r_nh[0]), float(r_nh[1])],
    "spearman_cell_milo_vs_ncorr": [float(rho_cell[0]), float(rho_cell[1])],
    "sign_agreement_among_milo_significant_nhoods": agree_sig,
    "n_nhoods_with_half_or_more_cells_cna_fdr10": int(len(cna_hit)), "sign_agreement_among_cna_hit_nhoods": agree_cna,
    "supertype_level_spearman": float(spearmanr(stm["milo"], stm["cna"])[0]) if len(stm) > 3 else None,
    "supertype_table": stm.round(4).reset_index().to_dict(orient="records"),
    "seconds": round(time.perf_counter() - t0, 1),
}
os.makedirs(args.outdir, exist_ok=True)
nh.reset_index().rename(columns={"index": "nhood"}).to_csv(f"{args.outdir}/milo_microglia_nhoods.csv", index=False)
cells.to_parquet(f"{args.outdir}/milo_microglia_cells.parquet", index=False)
json.dump(summary, open(f"{args.outdir}/milo_microglia_summary.json", "w"), indent=1)
print(json.dumps({k: v for k, v in summary.items() if k != "supertype_table"}, indent=1))
