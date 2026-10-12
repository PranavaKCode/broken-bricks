#!/usr/bin/env python
"""Broken Bricks Day 1: CNA on one SEA-AD MTG subclass, following preregistration.md section 1.2.

  sample = donor; phenotype = CPS_Local (then CPS_Local_ABeta, CPS_Local_pTau) from Global_and_Local_CPS.20260105.csv
  covariates = age, sex (male=1), PMI, method fraction (share of the donor's nuclei of this subclass from 10x multiome; C2)
  kNN graph on X_scVI, n_neighbors = 30, rebuilt on the retained cells
  reference donors dropped (they have no CPS in the CSV)

variants (primary = all):
  all    all 84 aged donors (PRIMARY: PC1, H1)
  nosad  without the 11 Severely Affected Donors (preregistered "with and without" run), kNN rebuilt
  ge50   only donors with >= 50 nuclei of this subclass (sensitivity), kNN rebuilt; CPS_Local only

usage: run_cna.py LABEL VARIANT [--nnull 10000] [--phenos a,b,c]
"""
import argparse, json, os, sys, time, warnings, platform
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import scanpy as sc
import cna

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *

ap = argparse.ArgumentParser()
ap.add_argument("label"); ap.add_argument("variant", choices=["all", "nosad", "ge50"])
ap.add_argument("--nnull", type=int, default=NNULL)
ap.add_argument("--phenos", default="")
ap.add_argument("--seed", type=int, default=SEED)
ap.add_argument("--outdir", default=OUT)
args = ap.parse_args()
h5name, display = SUBCLASSES[args.label]
phenos = [p for p in args.phenos.split(",") if p] or (["CPS_Local"] if args.variant == "ge50" else PHENOS)
os.makedirs(args.outdir, exist_ok=True)
stem = f"{args.outdir}/{args.label}__{args.variant}"
T = {}; T0 = time.perf_counter()
def tic(name, t0): T[name] = round(time.perf_counter() - t0, 2); print(f"[time] {name}: {T[name]} s", flush=True)

# ---- 1. obs + X_scVI + X_umap only (backed) ----
t = time.perf_counter()
a = ad.read_h5ad(f"{H5AD_DIR}/{h5name}", backed="r")
obs = a.obs.copy()
Z = np.asarray(a.obsm["X_scVI"], dtype=np.float32)
U = np.asarray(a.obsm["X_umap"], dtype=np.float32)
names = np.asarray(a.obs_names)
n0 = a.n_obs
a.file.close()
tic("read_obs_scvi_umap", t)
print(args.label, "cells", n0, "X_scVI", Z.shape, "X_umap", U.shape, flush=True)

# ---- 2. donors: drop reference donors; join CPS from the CSV-derived donor table ----
t = time.perf_counter()
donors = pd.read_csv(DONOR_TABLE, index_col="donor_id")
obs["donor_id"] = obs["donor_id"].astype(str)
has_cps = obs["donor_id"].isin(donors.index) & donors["CPS_Local"].reindex(obs["donor_id"]).notna().values
ref_flag = obs["Neurotypical reference"].astype(str) == "True"
assert (has_cps.values == ~ref_flag.values).all(), "aged-by-CSV != not-reference-flag"
keep = has_cps.values.copy()
cnt_all = obs.loc[keep, "donor_id"].value_counts()
dropped_donors = []
if args.variant == "nosad":
    sad = set(donors.index[donors["severe"] == 1])
    assert len(sad) == 11, len(sad)
    keep &= ~obs["donor_id"].isin(sad).values
    dropped_donors = sorted(sad)
elif args.variant == "ge50":
    small = sorted(cnt_all.index[cnt_all < 50])
    dropped_donors = small
    keep &= ~obs["donor_id"].isin(small).values
    if not small:
        json.dump({"label": args.label, "variant": "ge50", "identical_to": "all", "note": "all donors have >= 50 nuclei"},
                  open(stem + ".json", "w"), indent=1)
        print("ge50 identical to all; wrote stub"); sys.exit(0)
obs = obs.loc[keep].copy(); Z = Z[keep]; U = U[keep]; names = names[keep]
ids = pd.Index(sorted(obs["donor_id"].unique()))
cnt = obs["donor_id"].value_counts()
print("variant", args.variant, "donors", len(ids), "dropped", dropped_donors, "cells", len(obs),
      "cells/donor min/median/max", cnt.min(), cnt.median(), cnt.max(), flush=True)
mf = obs["assay"].astype(str).str.contains("multiome").groupby(obs["donor_id"]).mean().reindex(ids)
cov = pd.DataFrame({"age": donors["age"].reindex(ids), "sex": donors["sex_male"].reindex(ids).astype(float),
                    "pmi": donors["pmi"].reindex(ids), "method_frac": mf}, index=ids)
assert not cov.isna().any().any()
cov.index = cov.index.astype(str)
tic("filter_join", t)

adata = ad.AnnData(X=sp.csr_matrix((len(obs), 1), dtype=np.float32), obs=obs[["donor_id", "Supertype"]].copy())
adata.obs_names = pd.Index(names).astype(str)
adata.obs["donor_id"] = adata.obs["donor_id"].astype(str)
adata.obs["Supertype"] = adata.obs["Supertype"].astype(str)
adata.obsm["X_scVI"] = Z

# ---- 3. kNN ----
t = time.perf_counter()
sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors=NEIGHBORS, random_state=0)
tic("neighbors_k%d" % NEIGHBORS, t)

# ---- 4. CNA ----
res = {}
cells = pd.DataFrame({"obs_name": adata.obs_names.values, "donor_id": adata.obs["donor_id"].values,
                      "Supertype": adata.obs["Supertype"].values,
                      "umap_1": U[:, 0], "umap_2": U[:, 1]})
for ph in phenos:
    y = donors[ph].reindex(ids).astype(float); y.index = y.index.astype(str)
    t = time.perf_counter()
    with warnings.catch_warnings(record=True) as wl:
        warnings.simplefilter("always")
        r = cna.tl.association(adata, y, "donor_id", covs=cov.copy(), Nnull=args.nnull, seed=args.seed,
                               key_added=f"cna_{ph}", return_full=True)
    dt = round(time.perf_counter() - t, 1); T[f"cna_{ph}"] = dt
    nc = adata.obs[f"cna_{ph}"].astype(float); fdr = adata.obs[f"cna_{ph}_fdr"].astype(float)
    cells[f"ncorr_{ph}"] = nc.values.astype(np.float32); cells[f"fdr_{ph}"] = fdr.values.astype(np.float32)
    s5 = (fdr <= 0.05) & nc.notna(); s10 = (fdr <= 0.10) & nc.notna()
    st = []
    for stype, g in adata.obs.assign(nc=nc, s5=s5, s10=s10).groupby("Supertype"):
        st.append({"supertype": stype, "n_cells": int(len(g)), "mean_ncorr": float(g["nc"].mean()),
                   "n_pos_fdr5": int(((g.nc > 0) & g.s5).sum()), "n_neg_fdr5": int(((g.nc < 0) & g.s5).sum()),
                   "n_pos_fdr10": int(((g.nc > 0) & g.s10).sum()), "n_neg_fdr10": int(((g.nc < 0) & g.s10).sum())})
    res[ph] = {
        "global_p": float(r.p), "p_floor": 1.0 / (args.nnull + 1), "k_pcs": int(r.k), "ks": [int(x) for x in r.ks],
        "r2": float(r.r2), "null_r2_mean": float(r.nullr2_mean), "null_r2_std": float(r.nullr2_std),
        "n_donors": int(len(r.nam)), "n_neighborhoods_kept": int(r.kept.sum()), "n_cells": int(len(adata)),
        "n_cells_ncorr_nan": int(nc.isna().sum()),
        "fdr5_threshold": None if r.fdr_5p_t is None else float(r.fdr_5p_t),
        "fdr10_threshold": None if r.fdr_10p_t is None else float(r.fdr_10p_t),
        "n_cells_fdr5": int(s5.sum()), "n_cells_fdr10": int(s10.sum()),
        "pos_fdr5": int(((nc > 0) & s5).sum()), "neg_fdr5": int(((nc < 0) & s5).sum()),
        "pos_fdr10": int(((nc > 0) & s10).sum()), "neg_fdr10": int(((nc < 0) & s10).sum()),
        "ncorr_min": float(nc.min()), "ncorr_max": float(nc.max()), "ncorr_mean": float(nc.mean()),
        "assoc_seconds": dt, "warnings": sorted({str(w.message)[:160] for w in wl}),
        "supertypes": st,
    }
    print(f"[{ph}] P={r.p:.5f} k={r.k} r2={r.r2:.3f} fdr5/10 cells={int(s5.sum())}/{int(s10.sum())} "
          f"pos/neg@10%={res[ph]['pos_fdr10']}/{res[ph]['neg_fdr10']}  ({dt} s)", flush=True)
    del r

T["total_s"] = round(time.perf_counter() - T0, 1)
out = {"label": args.label, "display": display, "variant": args.variant, "h5ad": f"{H5AD_DIR}/{h5name}",
       "n_cells_file": int(n0), "n_cells_used": int(len(adata)), "n_donors": int(len(ids)),
       "dropped_donors": dropped_donors, "neighbors": NEIGHBORS, "nnull": args.nnull, "seed": args.seed,
       "covariates": list(cov.columns), "cells_per_donor": {"min": int(cnt.min()), "median": float(cnt.median()), "max": int(cnt.max())},
       "method_frac_mean": float(mf.mean()), "timings_s": T,
       "versions": {"python": platform.python_version(), "cna": "0.2.6",
                    "scanpy": sc.__version__, "anndata": ad.__version__, "numpy": np.__version__, "pandas": pd.__version__},
       "results": res}
json.dump(out, open(stem + ".json", "w"), indent=1)
cells.to_parquet(stem + "_cells.parquet", index=False)
print("wrote", stem + ".json", " total", T["total_s"], "s", flush=True)
