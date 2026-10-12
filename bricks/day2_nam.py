#!/usr/bin/env python
"""Broken Bricks Day 2 (cluster part): phenotype-free neighborhood abundance matrix (NAM) PCs for one SEA-AD MTG subclass.

Implements amendments A6 (items 1, 2, 4) and A8 (techNAM sensitivity): graph, full NAM, NAM PCs (primary: no covariates; covNAM: Day-1 covariates
residualised; foldPCA: PCA refitted inside each training fold with held-out donors projected) and 20 random split-half NAMs.
Nothing here touches CPS_Local: the phenotype enters only in the local script `day2_bricks.py`.

usage: day2_nam.py LABEL VARIANT [--nsplits 20] [--kmax 20] [--validate]
  VARIANT: all   = every aged donor with >= 50 nuclei of this subclass
           nosad = additionally without the 11 Severely Affected Donors
output: <OUT>/<label>__<variant>.npz  (+ .json with timings), graph cache <WORK>/graph_<label>__<variant>.npz
"""
import argparse, json, os, sys, time, platform
import numpy as np, pandas as pd, scipy.sparse as sp, scipy.stats as st
import anndata as ad, scanpy as sc

sys.path.insert(0, "/home/pranavak/orcd/pool/mantis-films/research/broken-bricks/cna")
from common import ROOT, H5AD_DIR, DONOR_TABLE, SUBCLASSES, NEIGHBORS

OUT = f"{ROOT}/results/day2"
WORK = f"{ROOT}/work/day2"
SEED = 20261011
MINCELLS = 50

ap = argparse.ArgumentParser()
ap.add_argument("label"); ap.add_argument("variant", choices=["all", "nosad"])
ap.add_argument("--nsplits", type=int, default=20)
ap.add_argument("--kmax", type=int, default=20)
ap.add_argument("--batch", type=int, default=4, help="half-NAMs diffused together")
ap.add_argument("--validate", action="store_true")
args = ap.parse_args()
os.makedirs(OUT, exist_ok=True); os.makedirs(WORK, exist_ok=True)
h5name, display = SUBCLASSES[args.label]
stem = f"{OUT}/{args.label}__{args.variant}"
T = {}; T0 = time.perf_counter()
def tic(name, t0): T[name] = round(time.perf_counter() - t0, 2); print(f"[time] {name}: {T[name]} s", flush=True)

# ---------------- 1. cells and donors ----------------
t = time.perf_counter()
a = ad.read_h5ad(f"{H5AD_DIR}/{h5name}", backed="r")
obs = a.obs[["donor_id", "assay", "Neurotypical reference", "Supertype"]].copy()
Z = np.asarray(a.obsm["X_scVI"], dtype=np.float32)
names = np.asarray(a.obs_names)
a.file.close()
donors = pd.read_csv(DONOR_TABLE, index_col="donor_id")
obs["donor_id"] = obs["donor_id"].astype(str)
aged = obs["donor_id"].isin(donors.index) & ~(obs["Neurotypical reference"].astype(str) == "True")
assert aged.sum() > 0
cnt_all = obs.loc[aged, "donor_id"].value_counts()
keep_d = set(cnt_all.index[cnt_all >= MINCELLS])
sad = set(donors.index[donors["severe"] == 1]); assert len(sad) == 11
if args.variant == "nosad":
    keep_d -= sad
keep = aged.values & obs["donor_id"].isin(keep_d).values
obs = obs.loc[keep].copy(); Z = Z[keep]
ids = sorted(keep_d)
D = len(ids); did = {d: i for i, d in enumerate(ids)}
sidx = obs["donor_id"].map(did).values.astype(np.int64)
n = len(obs)
counts = np.bincount(sidx, minlength=D)
print(args.label, args.variant, "donors", D, "cells", n, "min/median/max per donor", counts.min(), np.median(counts), counts.max(), flush=True)
all84 = sorted(donors.index)
assert len(all84) == 84
perm = np.random.default_rng(SEED).permutation(84)
fold84 = np.empty(84, int); fold84[perm] = np.arange(84) % 5          # fold of the i-th donor in sorted-ID order
fold = np.array([fold84[all84.index(d)] for d in ids])
mf = obs["assay"].astype(str).str.contains("multiome").groupby(obs["donor_id"]).mean().reindex(ids).values
cov = pd.DataFrame({"age": donors["age"].reindex(ids).values, "sex": donors["sex_male"].reindex(ids).astype(float).values,
                    "pmi": donors["pmi"].reindex(ids).values, "method_frac": mf}, index=ids)
assert not cov.isna().any().any()
tech = pd.read_csv(f"{ROOT}/results/day2/donor_tech.csv", index_col="donor_id").reindex(ids)
cov_tech = cov.copy()                                   # techNAM sensitivity: Day-1 covariates + the six technical covariates of A6 item 9 (A8)
cov_tech["log10_n_nuclei"] = np.log10(tech["n_nuclei"].values); cov_tech["rin"] = tech["rin_mean"].values
cov_tech["neun_frac"] = tech["neun_frac"].values; cov_tech["seq_saturation"] = tech["seq_saturation"].values
assert not cov_tech.isna().any().any()
tic("read_filter", t)

# ---------------- 2. graph (cached) ----------------
t = time.perf_counter()
gfile = f"{WORK}/graph_{args.label}__{args.variant}.npz"
if os.path.exists(gfile):
    conn = sp.load_npz(gfile).tocsr()
    print("graph loaded from cache", gfile, flush=True)
else:
    adata = ad.AnnData(X=sp.csr_matrix((n, 1), dtype=np.float32))
    adata.obs_names = pd.Index(names[keep]).astype(str)
    adata.obsm["X_scVI"] = Z
    sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors=NEIGHBORS, random_state=0)
    conn = adata.obsp["connectivities"].tocsr()
    sp.save_npz(gfile, conn)
    del adata
tic("graph", t)
conn = conn.astype(np.float64)
colsums = np.asarray(conn.sum(axis=0)).ravel() + 1.0

# ---------------- 3. diffusion helpers (replicates cna 0.2.6 _nam) ----------------
def diffuse_step(s, dtype):
    sn = s / colsums[:, None].astype(dtype)
    return conn_d[dtype].dot(sn) + sn
conn_d = {np.float64: conn, np.float32: conn.astype(np.float32)}

def full_nam():
    S = np.zeros((n, D)); S[np.arange(n), sidx] = 1.0
    C = S.sum(axis=0)
    s = S; prev = np.inf; nsteps = None
    for i in range(15):
        s = diffuse_step(s, np.float64)
        medkurt = np.median(st.kurtosis(s / C, axis=1))
        if prev - medkurt < 3 and i + 1 >= 3:
            nsteps = i + 1; break
        prev = medkurt
    if nsteps is None: nsteps = 15
    return (s / C).T, nsteps                      # D x n

def half_nams(halfidx_list, nsteps):
    """halfidx_list: list of boolean masks (n,) of cells in the half-sample; returns list of (D x n) NAMs"""
    nb = len(halfidx_list)
    S = np.zeros((n, D * nb), dtype=np.float32)
    Ch = np.zeros((nb, D))
    for b, m in enumerate(halfidx_list):
        rows = np.where(m)[0]
        S[rows, b * D + sidx[rows]] = 1.0
        Ch[b] = np.bincount(sidx[rows], minlength=D)
    s = S
    for _ in range(nsteps):
        s = diffuse_step(s, np.float32)
    out = []
    for b in range(nb):
        out.append((s[:, b * D:(b + 1) * D] / Ch[b].astype(np.float32)).T)
    return out, Ch

def nam_pcs(NAM, kmax, covs=None):
    """cna _resid_nam / svd_nam with batches=None: center, (residualise covs), standardise (ddof=1), SVD of the Gram matrix."""
    X = np.asarray(NAM, dtype=np.float64)
    X = X - X.mean(axis=0)
    if covs is not None:
        Cv = np.asarray(covs, dtype=np.float64)
        Cv = (Cv - Cv.mean(axis=0)) / Cv.std(axis=0, ddof=1)
        M = np.eye(len(X)) - Cv @ np.linalg.solve(Cv.T @ Cv, Cv.T)
        X = M @ X
    sd = X.std(axis=0, ddof=1)
    ok = sd > 0
    X = X[:, ok] / sd[ok]
    X = X - X.mean(axis=0)
    G = X @ X.T
    sv, U = np.linalg.eigh(G)
    o = np.argsort(-sv); sv, U = sv[o], U[:, o]
    sv = np.clip(sv, 1e-12, None)
    return U[:, :kmax] * np.sqrt(sv[:kmax]), sv, U[:, :kmax]

def fold_pcs(NAM, kmax):
    """PCA refitted on the training donors of each fold; held-out donors projected. returns (5, D, kmax)"""
    X0 = np.asarray(NAM, dtype=np.float64)
    out = np.zeros((5, D, kmax))
    for f in range(5):
        tr = np.where(fold != f)[0]; te = np.where(fold == f)[0]
        mu = X0[tr].mean(axis=0); sd = X0[tr].std(axis=0, ddof=1)
        ok = sd > 0
        Xtr = (X0[tr][:, ok] - mu[ok]) / sd[ok]          # columns are already centred on the training donors
        G = Xtr @ Xtr.T
        sv, U = np.linalg.eigh(G); o = np.argsort(-sv); sv, U = sv[o], U[:, o]
        sv = np.clip(sv, 1e-12, None)
        k = min(kmax, len(tr) - 1)
        V = Xtr.T @ U[:, :k] / np.sqrt(sv[:k])
        out[f][tr, :k] = U[:, :k] * np.sqrt(sv[:k])
        Xte = (X0[te][:, ok] - mu[ok]) / sd[ok]
        out[f][te, :k] = Xte @ V
    return out

# ---------------- 4. full NAM and PCs ----------------
t = time.perf_counter()
NAM, nsteps = full_nam()
tic("full_nam", t); print("nsteps", nsteps, "NAM shape", NAM.shape, flush=True)
t = time.perf_counter()
pcs_full, sv_full, U_full = nam_pcs(NAM, args.kmax)
pcs_cov, sv_cov, _ = nam_pcs(NAM, args.kmax, covs=cov.values)
pcs_tech, sv_tech, _ = nam_pcs(NAM, args.kmax, covs=cov_tech.values)
varexp = sv_full / sv_full.sum()
if args.validate:
    from cna.tools._nam import _resid_nam
    df = pd.DataFrame(NAM, index=pd.Index(ids, name="donor_id"))
    r = _resid_nam(df, None, None)
    Uc = r.namresid_sampleXpc.values[:, :args.kmax]
    cors = [abs(np.corrcoef(Uc[:, j], U_full[:, j])[0, 1]) for j in range(args.kmax)]
    print("VALIDATE vs cna _resid_nam: |corr| of first PCs", np.round(cors[:6], 6), "min over", args.kmax, "=", round(min(cors), 6), flush=True)
    T["validate_min_abs_corr"] = float(min(cors))
fpc = fold_pcs(NAM, args.kmax)
tic("pcs", t)
del NAM

# ---------------- 5. split halves ----------------
halfpcs = np.zeros((args.nsplits, 2, D, args.kmax), dtype=np.float64)
halfpcs_tech = np.zeros((args.nsplits, 2, D, args.kmax), dtype=np.float64)
halfn = np.zeros((args.nsplits, 2, D), dtype=np.int32)
masks = []
for s in range(args.nsplits):
    rng = np.random.default_rng(SEED + s)
    mA = np.zeros(n, bool)
    for d in range(D):
        idx = np.where(sidx == d)[0]
        p = rng.permutation(idx)
        mA[p[: (len(p) + 1) // 2]] = True
    masks.append((mA, ~mA))
flat = [(s, h, masks[s][h]) for s in range(args.nsplits) for h in range(2)]
t = time.perf_counter()
for i in range(0, len(flat), args.batch):
    chunk = flat[i:i + args.batch]
    nams, Ch = half_nams([m for (_, _, m) in chunk], nsteps)
    for (s, h, _), N_, c in zip(chunk, nams, Ch):
        halfpcs[s, h] = nam_pcs(N_, args.kmax)[0]
        halfpcs_tech[s, h] = nam_pcs(N_, args.kmax, covs=cov_tech.values)[0]
        halfn[s, h] = c
    print(f"[halves] {i + len(chunk)}/{len(flat)}  {round(time.perf_counter() - t, 1)} s", flush=True)
    del nams
tic("halves", t)

T["total_s"] = round(time.perf_counter() - T0, 1)
np.savez_compressed(stem + ".npz", donors=np.array(ids), n_nuclei=counts, fold=fold, nsteps=nsteps, pcs_full=pcs_full, sv_full=sv_full,
                    pcs_cov=pcs_cov, pcs_tech=pcs_tech, sv_tech=sv_tech, foldpcs=fpc, halfpcs=halfpcs, halfpcs_tech=halfpcs_tech, halfn=halfn, method_frac=mf)
json.dump({"label": args.label, "display": display, "variant": args.variant, "n_cells": int(n), "n_donors": int(D),
           "dropped_lt50": sorted(set(cnt_all.index[cnt_all < MINCELLS])), "nsteps": int(nsteps), "kmax": args.kmax,
           "nsplits": args.nsplits, "seed": SEED, "neighbors": NEIGHBORS, "varexp_first10": [float(x) for x in varexp[:10]],
           "timings_s": T, "versions": {"python": platform.python_version(), "numpy": np.__version__, "scanpy": sc.__version__,
                                        "anndata": ad.__version__, "pandas": pd.__version__}},
          open(stem + ".json", "w"), indent=1)
print("wrote", stem + ".npz", "total", T["total_s"], "s", flush=True)
