#!/usr/bin/env python
"""Final check 2: the Day-1 CNA run (cna/run_cna.py) repeated with donor technical covariates added to the covariate set.

Identical to cna/run_cna.py in every setting (kNN k = 30 on X_scVI, random_state 0, Nnull = 10,000, seed 0, CPS_Local only, all 84 aged donors or the
no-SAD variant) except for the covariate matrix. cna standardises every covariate internally, so scales do not matter.

models (covariates beyond the Day-1 base set age, sex, PMI, method fraction):
  base   none (reproduces the Day-1 primary run; a pipeline control)
  sat    + sequencing saturation (nuclei-weighted donor mean of the library-level GEX_sequencing_saturation; donors without a value get the median)
  tech3  + sequencing saturation + RIN (donor mean) + log10 total final nuclei of the donor (all cell classes)
  neun   + NeuN-positive fraction of the donor (FANS), a sorting parameter that is not a severity proxy
  rin    + RIN only          nuc  + log10 total final nuclei only
  satperm  + sequencing saturation PERMUTED across donors (rng 20261011): a negative control that costs one degree of freedom and carries no information

usage: run_cna_tech.py LABEL MODEL [--variant all|nosad] [--nnull 10000]
Technical table: results/day2/donor_tech.csv (bricks/tech_table.py).
"""
import argparse, json, os, sys, time, warnings, platform
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import scanpy as sc
import cna

sys.path.insert(0, "/home/pranavak/orcd/pool/mantis-films/research/broken-bricks/cna")
from common import *

ap = argparse.ArgumentParser()
ap.add_argument("label"); ap.add_argument("model", choices=["base", "sat", "tech3", "neun", "rin", "nuc", "satperm"])
ap.add_argument("--variant", default="all", choices=["all", "nosad"])
ap.add_argument("--nnull", type=int, default=NNULL)
ap.add_argument("--seed", type=int, default=SEED)
ap.add_argument("--outdir", default=f"{ROOT}/results/final_checks/cna_tech")
args = ap.parse_args()
h5name, display = SUBCLASSES[args.label]
os.makedirs(args.outdir, exist_ok=True)
stem = f"{args.outdir}/{args.label}__{args.variant}__{args.model}"
T = {}; T0 = time.perf_counter()

a = ad.read_h5ad(f"{H5AD_DIR}/{h5name}", backed="r")
obs = a.obs.copy()
Z = np.asarray(a.obsm["X_scVI"], dtype=np.float32)
U = np.asarray(a.obsm["X_umap"], dtype=np.float32)
names = np.asarray(a.obs_names)
n0 = a.n_obs
a.file.close()

donors = pd.read_csv(DONOR_TABLE, index_col="donor_id")
tech = pd.read_csv(f"{ROOT}/results/day2/donor_tech.csv", index_col="donor_id")
obs["donor_id"] = obs["donor_id"].astype(str)
has_cps = obs["donor_id"].isin(donors.index) & donors["CPS_Local"].reindex(obs["donor_id"]).notna().values
ref_flag = obs["Neurotypical reference"].astype(str) == "True"
assert (has_cps.values == ~ref_flag.values).all(), "aged-by-CSV != not-reference-flag"
keep = has_cps.values.copy()
dropped_donors = []
if args.variant == "nosad":
    sad = set(donors.index[donors["severe"] == 1])
    assert len(sad) == 11, len(sad)
    keep &= ~obs["donor_id"].isin(sad).values
    dropped_donors = sorted(sad)
obs = obs.loc[keep].copy(); Z = Z[keep]; U = U[keep]; names = names[keep]
ids = pd.Index(sorted(obs["donor_id"].unique()))
cnt = obs["donor_id"].value_counts()
mf = obs["assay"].astype(str).str.contains("multiome").groupby(obs["donor_id"]).mean().reindex(ids)
cov = pd.DataFrame({"age": donors["age"].reindex(ids), "sex": donors["sex_male"].reindex(ids).astype(float),
                    "pmi": donors["pmi"].reindex(ids), "method_frac": mf}, index=ids)
assert not cov.isna().any().any()
t = tech.reindex(ids)
missing = {}
if args.model in ("sat", "tech3"):
    sat = t["seq_saturation"]; missing["seq_saturation"] = int(sat.isna().sum())
    cov["seq_saturation"] = sat.fillna(sat.median())
if args.model == "tech3":
    cov["rin"] = t["rin_mean"].fillna(t["rin_mean"].median()); missing["rin"] = int(t["rin_mean"].isna().sum())
    cov["log10_nuclei"] = np.log10(t["n_nuclei"]); missing["n_nuclei"] = int(t["n_nuclei"].isna().sum())
if args.model == "rin":
    cov["rin"] = t["rin_mean"].fillna(t["rin_mean"].median()); missing["rin"] = int(t["rin_mean"].isna().sum())
if args.model == "nuc":
    cov["log10_nuclei"] = np.log10(t["n_nuclei"]); missing["n_nuclei"] = int(t["n_nuclei"].isna().sum())
if args.model == "satperm":
    sat = t["seq_saturation"]; vals = sat.fillna(sat.median()).values
    cov["seq_saturation_permuted"] = vals[np.random.default_rng(20261011).permutation(len(vals))]
if args.model == "neun":
    cov["neun_frac"] = t["neun_frac"].fillna(t["neun_frac"].median()); missing["neun_frac"] = int(t["neun_frac"].isna().sum())
assert not cov.isna().any().any()
cov.index = cov.index.astype(str)
print(args.label, args.variant, args.model, "donors", len(ids), "cells", len(obs), "covariates", list(cov.columns), "imputed", missing, flush=True)

adata = ad.AnnData(X=sp.csr_matrix((len(obs), 1), dtype=np.float32), obs=obs[["donor_id", "Supertype"]].copy())
adata.obs_names = pd.Index(names).astype(str)
adata.obs["donor_id"] = adata.obs["donor_id"].astype(str)
adata.obs["Supertype"] = adata.obs["Supertype"].astype(str)
adata.obsm["X_scVI"] = Z
sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors=NEIGHBORS, random_state=0)

ph = "CPS_Local"
y = donors[ph].reindex(ids).astype(float); y.index = y.index.astype(str)
t1 = time.perf_counter()
with warnings.catch_warnings(record=True) as wl:
    warnings.simplefilter("always")
    r = cna.tl.association(adata, y, "donor_id", covs=cov.copy(), Nnull=args.nnull, seed=args.seed, key_added=f"cna_{ph}", return_full=True)
dt = round(time.perf_counter() - t1, 1)
nc = adata.obs[f"cna_{ph}"].astype(float); fdr = adata.obs[f"cna_{ph}_fdr"].astype(float)
s5 = (fdr <= 0.05) & nc.notna(); s10 = (fdr <= 0.10) & nc.notna()
cells = pd.DataFrame({"obs_name": adata.obs_names.values, "donor_id": adata.obs["donor_id"].values, "Supertype": adata.obs["Supertype"].values,
                      "ncorr": nc.values.astype(np.float32), "fdr": fdr.values.astype(np.float32)})
st = []
for stype, g in adata.obs.assign(nc=nc, s5=s5, s10=s10).groupby("Supertype"):
    st.append({"supertype": stype, "n_cells": int(len(g)), "mean_ncorr": float(g["nc"].mean()),
               "n_pos_fdr10": int(((g.nc > 0) & g.s10).sum()), "n_neg_fdr10": int(((g.nc < 0) & g.s10).sum())})
res = {"global_p": float(r.p), "p_floor": 1.0 / (args.nnull + 1), "k_pcs": int(r.k), "ks": [int(x) for x in r.ks], "r2": float(r.r2),
       "null_r2_mean": float(r.nullr2_mean), "null_r2_std": float(r.nullr2_std), "n_donors": int(len(r.nam)),
       "n_neighborhoods_kept": int(r.kept.sum()), "n_cells": int(len(adata)), "n_cells_fdr5": int(s5.sum()), "n_cells_fdr10": int(s10.sum()),
       "pos_fdr10": int(((nc > 0) & s10).sum()), "neg_fdr10": int(((nc < 0) & s10).sum()),
       "pos_fdr5": int(((nc > 0) & s5).sum()), "neg_fdr5": int(((nc < 0) & s5).sum()),
       "assoc_seconds": dt, "warnings": sorted({str(w.message)[:160] for w in wl}), "supertypes": st}
print(f"[{ph}] P={r.p:.5f} k={r.k} r2={r.r2:.3f} fdr10 cells={int(s10.sum())} pos/neg@10%={res['pos_fdr10']}/{res['neg_fdr10']} ({dt} s)", flush=True)
out = {"label": args.label, "display": display, "variant": args.variant, "model": args.model, "covariates": list(cov.columns), "imputed_donors": missing,
       "n_cells_used": int(len(adata)), "n_donors": int(len(ids)), "dropped_donors": dropped_donors, "neighbors": NEIGHBORS, "nnull": args.nnull, "seed": args.seed,
       "total_s": round(time.perf_counter() - T0, 1), "versions": {"python": platform.python_version(), "cna": "0.2.6", "scanpy": sc.__version__},
       "result": res}
json.dump(out, open(stem + ".json", "w"), indent=1)
cells.to_parquet(stem + "_cells.parquet", index=False)
print("wrote", stem + ".json", flush=True)
