#!/usr/bin/env python
"""Final check 6 (optional): independent recomputation of one subclass's split-half reliability of the brick score B from nuclei.

Independent of bricks/day2_nam.py + bricks/day2_bricks.py in these respects: the full NAM comes from cna's own `_nam`, the half-sample diffusion from
cna's own `diffuse` (same graph, same number of steps), the PCA is a plain SVD of the standardised NAM, B is scikit-learn RidgeCV (inner leave-one-out)
in 5 folds by donor, and the 20 random half-splits use different seeds (777 + s, not 20261011 + s). Everything follows amendment A6 items 1-4:
  donors with >= 50 nuclei of the subclass, kNN k = 30 on X_scVI (random_state 0), k = 15 PCs divided by the SD of PC1, alphas logspace(-4, 4, 33),
  fold = position mod 5 of the permutation default_rng(20261011).permutation(84) over the 84 donors in sorted-ID order, B = out-of-fold prediction,
  split-half r = Pearson r of the two half vectors of B, median over splits.
usage: reliability_check.py LABEL [--nsplits 20]      writes results/final_checks/reliability_check/<label>.json
"""
import argparse, json, os, sys, time
import numpy as np, pandas as pd, scipy.sparse as sp
import anndata as ad, scanpy as sc
from sklearn.linear_model import RidgeCV
from cna.tools._nam import _nam, diffuse

sys.path.insert(0, "/home/pranavak/orcd/pool/mantis-films/research/broken-bricks/cna")
from common import ROOT, H5AD_DIR, DONOR_TABLE, SUBCLASSES, NEIGHBORS

ap = argparse.ArgumentParser(); ap.add_argument("label"); ap.add_argument("--nsplits", type=int, default=20); args = ap.parse_args()
OUT = f"{ROOT}/results/final_checks/reliability_check"; os.makedirs(OUT, exist_ok=True)
T0 = time.perf_counter()
h5name, display = SUBCLASSES[args.label]
a = ad.read_h5ad(f"{H5AD_DIR}/{h5name}", backed="r")
obs = a.obs[["donor_id", "Neurotypical reference"]].copy(); Z = np.asarray(a.obsm["X_scVI"], dtype=np.float32); names = np.asarray(a.obs_names); a.file.close()
donors = pd.read_csv(DONOR_TABLE, index_col="donor_id")
obs["donor_id"] = obs["donor_id"].astype(str)
aged = obs["donor_id"].isin(donors.index) & ~(obs["Neurotypical reference"].astype(str) == "True")
cnt = obs.loc[aged, "donor_id"].value_counts()
ids = sorted(cnt.index[cnt >= 50]); D = len(ids)
keep = aged.values & obs["donor_id"].isin(ids).values
obs = obs.loc[keep].copy(); Z = Z[keep]; names = names[keep]
print(args.label, "donors", D, "nuclei", len(obs), flush=True)

adata = ad.AnnData(X=sp.csr_matrix((len(obs), 1), dtype=np.float32), obs=obs[["donor_id"]].copy())
adata.obs_names = pd.Index(names).astype(str); adata.obs["donor_id"] = pd.Categorical(adata.obs["donor_id"], categories=ids)
adata.obsm["X_scVI"] = Z
sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors=NEIGHBORS, random_state=0)

# full NAM with cna's own stopping rule
NAM = _nam(adata, "donor_id", sids=ids)                       # DataFrame D x n
S_full = pd.get_dummies(adata.obs["donor_id"])[ids]
# how many steps did it take? re-derive from the kurtosis rule by calling with nsteps unset and counting via cna's generator
from cna.tools._nam import diffuse_stepwise
import scipy.stats as st
C = S_full.sum(axis=0); prev = np.inf; nsteps = None
for i, s in enumerate(diffuse_stepwise(adata, S_full)):
    mk = np.median(st.kurtosis(s / C, axis=1))
    if prev - mk < 3 and i + 1 >= 3:
        nsteps = i + 1; break
    prev = mk
print("diffusion steps", nsteps, flush=True)

def svd_scores(N, k=15):
    X = np.asarray(N, float); X = X - X.mean(axis=0); sd = X.std(axis=0, ddof=1); ok = sd > 0
    X = X[:, ok] / sd[ok]; X = X - X.mean(axis=0)
    U, s, _ = np.linalg.svd(X @ X.T)
    return U[:, :k] * np.sqrt(s[:k])

all84 = sorted(donors.index); perm = np.random.default_rng(20261011).permutation(84)
fold84 = np.empty(84, int); fold84[perm] = np.arange(84) % 5
fold = np.array([fold84[all84.index(d)] for d in ids])
y = donors["CPS_Local"].reindex(ids).values.astype(float)
ALPHAS = np.logspace(-4, 4, 33)

def oof_B(P):
    P = P / P[:, 0].std(ddof=1)
    B = np.zeros(D)
    for f in range(5):
        tr, te = fold != f, fold == f
        m = RidgeCV(alphas=ALPHAS).fit(P[tr], y[tr]); B[te] = m.predict(P[te])
    return B

from scipy.stats import spearmanr
B_full = oof_B(svd_scores(NAM.values))
rho = spearmanr(B_full, y)
print(f"full-sample out-of-fold rho {rho[0]:.3f}", flush=True)

rs = []
sidx = adata.obs["donor_id"].cat.codes.values
for s_ in range(args.nsplits):
    rng = np.random.default_rng(777 + s_)
    mA = np.zeros(len(adata), bool)
    for d in range(D):
        idx = np.where(sidx == d)[0]; p = rng.permutation(idx); mA[p[: (len(p) + 1) // 2]] = True
    Bs = []
    for m in (mA, ~mA):
        Sh = S_full.values.astype(np.float32) * m[:, None]
        Ch = Sh.sum(axis=0)
        sh = diffuse(adata, Sh, nsteps)
        Nh = (sh / Ch).T
        Bs.append(oof_B(svd_scores(Nh)))
    rs.append(float(np.corrcoef(Bs[0], Bs[1])[0, 1]))
    print(f"split {s_}: r = {rs[-1]:.3f}", flush=True)
res = {"label": args.label, "display": display, "n_donors": D, "n_nuclei": int(len(adata)), "nsteps": int(nsteps), "full_oof_rho": float(rho[0]),
       "split_half_r_median": float(np.median(rs)), "split_half_r_q25": float(np.percentile(rs, 25)), "split_half_r_q75": float(np.percentile(rs, 75)),
       "split_half_r_min": float(np.min(rs)), "r_all": rs, "seed_base": 777, "seconds": round(time.perf_counter() - T0, 1)}
json.dump(res, open(f"{OUT}/{args.label}.json", "w"), indent=1)
print(json.dumps({k: v for k, v in res.items() if k != "r_all"}), flush=True)
