#!/usr/bin/env python
"""Broken Bricks Day 3, item 1: red/blue gene programs per SEA-AD MTG subclass (amendment A7 item 1).

For one subclass: <= 50,000 aged nuclei (all if fewer; otherwise donor-stratified, min(n_donor, 100) per donor guaranteed, remainder
proportional, seed 20261011), raw counts from the CELLxGENE h5ad `raw.X`, log-normalised expression = log1p(count / total x 1e4),
genes detected in >= 1% of the sampled nuclei. Per-nucleus ncorr = Day-1 CPS_Local run (variant `all` primary; `nosad` sensitivity).

Primary statistic ("donor-blocked Spearman"): per donor the mean log-normalised expression and the mean ncorr over the donor's
sampled nuclei, then Spearman across donors that have >= 50 nuclei of the subclass. Sensitivity: within-donor Spearman, Fisher-z
averaged (weights n-3) over donors with >= 50 sampled nuclei (only for variant all).

usage: programs.py LABEL [all|nosad]
output: <ROOT>/results/day3/programs/<label>__<variant>_genes.csv (all tested genes with statistics and the program flags)
"""
import json, os, sys, time
import numpy as np, pandas as pd, scipy.sparse as sp, scipy.stats as st, h5py
sys.path.insert(0, "/home/pranavak/orcd/pool/mantis-films/research/broken-bricks/cna")
from common import ROOT, H5AD_DIR, SUBCLASSES

label = sys.argv[1]; variant = sys.argv[2] if len(sys.argv) > 2 else "all"
OUT = f"{ROOT}/results/day3/programs"; os.makedirs(OUT, exist_ok=True)
SEED = 20261011; CAP = 50000; MINDON = 50; DETECT = 0.01; TOP = 200; CHUNK = 20000
h5name, display = SUBCLASSES[label]
path = f"{H5AD_DIR}/{h5name}"
T0 = time.perf_counter()

cells = pd.read_parquet(f"{ROOT}/results/day1/{label}__{variant}_cells.parquet", columns=["obs_name", "donor_id", "ncorr_CPS_Local"])
cells = cells[np.isfinite(cells["ncorr_CPS_Local"])].reset_index(drop=True)
n_by_d = cells["donor_id"].value_counts()
elig = set(n_by_d.index[n_by_d >= MINDON])
cells = cells[cells["donor_id"].isin(elig)].reset_index(drop=True)
print(label, variant, "cells with ncorr (eligible donors)", len(cells), "donors", len(elig), flush=True)

# ---- donor-stratified sample
rng = np.random.default_rng(SEED)
if len(cells) <= CAP:
    sel = np.arange(len(cells))
else:
    cnt = cells["donor_id"].value_counts()
    guar = np.minimum(cnt.values, 100).astype(int)
    rest = CAP - guar.sum()
    room = cnt.values - guar
    q = room * rest / room.sum()
    base = np.minimum(np.floor(q).astype(int), room)
    short = rest - base.sum()
    order = np.argsort(-(q - np.floor(q)), kind="stable")
    for i in order:
        if short <= 0: break
        if base[i] < room[i]: base[i] += 1; short -= 1
    take = dict(zip(cnt.index, guar + base))
    parts = []
    for d, g in cells.groupby("donor_id"):
        idx = g.index.values
        parts.append(idx if take[d] >= len(idx) else rng.choice(idx, take[d], replace=False))
    sel = np.sort(np.concatenate(parts))
cs = cells.loc[sel].reset_index(drop=True)
print("sampled nuclei", len(cs), flush=True)

# ---- read raw counts for the sampled nuclei (chunked CSR reads)
import anndata as ad
a = ad.read_h5ad(path, backed="r")
obsn = pd.Index(a.obs_names.astype(str)); var = a.raw.var.copy(); a.file.close()
rows = obsn.get_indexer(cs["obs_name"].astype(str))
assert (rows >= 0).all()
order_rows = np.argsort(rows); rows_sorted = rows[order_rows]
nvar = len(var)
blocks = []
with h5py.File(path, "r") as h:
    g = h["raw/X"]; indptr = g["indptr"][:]
    ntot = len(indptr) - 1
    for r0 in range(0, ntot, CHUNK):
        r1 = min(r0 + CHUNK, ntot)
        lo, hi = np.searchsorted(rows_sorted, [r0, r1])
        if hi == lo: continue
        a0, a1 = indptr[r0], indptr[r1]
        M = sp.csr_matrix((g["data"][a0:a1], g["indices"][a0:a1], indptr[r0:r1 + 1] - a0), shape=(r1 - r0, nvar))
        blocks.append(M[rows_sorted[lo:hi] - r0])
        print(f"  read rows {r0}-{r1}  kept {hi - lo}  {time.perf_counter() - T0:.0f}s", flush=True)
C = sp.vstack(blocks).tocsr()
inv = np.empty_like(order_rows); inv[order_rows] = np.arange(len(order_rows))
C = C[inv]                                                  # back to the order of `cs`
tot = np.asarray(C.sum(axis=1)).ravel()
det = np.asarray((C > 0).sum(axis=0)).ravel() / C.shape[0]
Xn = sp.diags((1e4 / np.maximum(tot, 1)).astype(np.float32)) @ C.astype(np.float32)
Xn = Xn.tocsr(); Xn.data = np.log1p(Xn.data)
mean_le = np.asarray(Xn.mean(axis=0)).ravel()
print("normalised; total counts median", np.median(tot), flush=True)

genes_ok = np.where(det >= DETECT)[0]
print("genes detected in >=1% of nuclei:", len(genes_ok), "of", nvar, flush=True)

# ---- donor-blocked Spearman
did = pd.Index(sorted(cs["donor_id"].unique())); D = len(did)
dcode = did.get_indexer(cs["donor_id"])
Ind = sp.csr_matrix((np.ones(len(cs)), (dcode, np.arange(len(cs)))), shape=(D, len(cs)))
nper = np.asarray(Ind.sum(axis=1)).ravel()
Em = (sp.diags(1.0 / nper) @ Ind @ Xn[:, genes_ok]).toarray()          # D x genes : donor mean log-normalised expression
ncm = (sp.diags(1.0 / nper) @ Ind) @ cs["ncorr_CPS_Local"].values       # donor mean ncorr
Rg = st.rankdata(Em, axis=0); Rn = st.rankdata(ncm)
Rg = (Rg - Rg.mean(0)) / Rg.std(0, ddof=1).clip(1e-12); Rn = (Rn - Rn.mean()) / Rn.std(ddof=1)
rho = (Rg * Rn[:, None]).sum(0) / (D - 1)
# approximate P (t distribution) for ranking transparency only
tt = rho * np.sqrt((D - 2) / np.clip(1 - rho ** 2, 1e-12, None)); pval = 2 * st.t.sf(np.abs(tt), D - 2)

res = pd.DataFrame({"ensembl_id": var.index.values[genes_ok], "symbol": var["feature_name"].astype(str).values[genes_ok],
                    "biotype": var["feature_type"].astype(str).values[genes_ok], "detect_frac": det[genes_ok], "mean_logexpr": mean_le[genes_ok],
                    "rho_donor_blocked": rho, "p_donor_blocked_approx": pval})

# ---- within-donor Spearman (sensitivity, variant all only)
if variant == "all":
    zs, ws = np.zeros(len(genes_ok)), 0.0
    Xg = Xn[:, genes_ok].tocsr()
    for d in range(D):
        ix = np.where(dcode == d)[0]
        if len(ix) < MINDON: continue
        Xd = Xg[ix].toarray()
        rx = st.rankdata(Xd, axis=0); ry = st.rankdata(cs["ncorr_CPS_Local"].values[ix])
        rx = (rx - rx.mean(0)) / rx.std(0, ddof=1).clip(1e-12); ry = (ry - ry.mean()) / ry.std(ddof=1)
        r = np.clip((rx * ry[:, None]).sum(0) / (len(ix) - 1), -0.999, 0.999)
        zs += np.arctanh(r) * (len(ix) - 3); ws += len(ix) - 3
    res["rho_within_donor"] = np.tanh(zs / ws)

# ---- program flags: top 200 positive (red), 200 most negative with rho < 0 (blue); ties by symbol
def top_sets(df, col):
    d = df.sort_values([col, "symbol"], ascending=[False, True])
    red = d[d[col] > 0].head(TOP)["ensembl_id"].tolist()
    d2 = df.sort_values([col, "symbol"], ascending=[True, True])
    blue = d2[d2[col] < 0].head(TOP)["ensembl_id"].tolist()
    return red, blue
red, blue = top_sets(res, "rho_donor_blocked")
res["program"] = np.where(res.ensembl_id.isin(red), "red", np.where(res.ensembl_id.isin(blue), "blue", ""))
if variant == "all":
    rw, bw = top_sets(res, "rho_within_donor")
    res["program_within_donor"] = np.where(res.ensembl_id.isin(rw), "red", np.where(res.ensembl_id.isin(bw), "blue", ""))
    pc = res[res.biotype == "protein_coding"]
    rp, bp = top_sets(pc, "rho_donor_blocked")
    res["program_protein_coding"] = np.where(res.ensembl_id.isin(rp), "red", np.where(res.ensembl_id.isin(bp), "blue", ""))
res.to_csv(f"{OUT}/{label}__{variant}_genes.csv", index=False)
meta = {"label": label, "variant": variant, "n_nuclei_sampled": int(len(cs)), "n_donors": int(D), "n_genes_tested": int(len(genes_ok)),
        "n_nuclei_available": int(len(cells)), "nuclei_per_donor_sampled_min_med_max": [int(nper.min()), float(np.median(nper)), int(nper.max())],
        "top_red_rho": float(res.loc[res.program == "red", "rho_donor_blocked"].min()), "top_blue_rho": float(res.loc[res.program == "blue", "rho_donor_blocked"].max()),
        "seconds": round(time.perf_counter() - T0, 1)}
json.dump(meta, open(f"{OUT}/{label}__{variant}_meta.json", "w"), indent=1)
print(meta, flush=True)
print("red head:", res[res.program == "red"].sort_values("rho_donor_blocked", ascending=False).symbol.head(15).tolist())
print("blue head:", res[res.program == "blue"].sort_values("rho_donor_blocked").symbol.head(15).tolist())
