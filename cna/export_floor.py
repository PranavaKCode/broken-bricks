#!/usr/bin/env python
"""Export the Mantis "cells" floor (protocol 06 section 4.3 Day 4, floor 1) as a slim h5ad.

Cells: all aged Microglia-PVM nuclei + N nuclei of each of astro, oligo, OPC, Sst, L2/3 IT, Pvalb with
N = floor((150000 - n_microglia) / 6), allocated across donors proportionally to the donor's nuclei (largest remainder),
seed 20261011.  Reference donors are not in the CNA results and so are not in the floor.

obs: subclass, donor_id, supertype, leiden (global int code of subclass+supertype), CPS / CPS_ABeta / CPS_pTau (donor-level, from
the CPS CSV), ncorr_* and fdr_* from the primary CNA run, redblue_* flags (+1 red = FDR<=level and ncorr>0, -1 blue, 0 otherwise),
apoe, apoe4_carrier, cognitive_status, dementia, adnc, braak, age_at_death, sex, severe_donor.
obsm: X_umap (SEA-AD native UMAPs scaled into a 20x20 box and laid out on a grid), X_umap_subclass (native), X_pca (= X_scVI, 10-d;
the latent space is PER SUBCLASS, so distances across subclasses are not meaningful).
X: placeholder (one all-zero column) until the Day-3 program genes exist.
"""
import argparse, json, os, sys
import numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *

ap = argparse.ArgumentParser()
ap.add_argument("--indir", default=OUT); ap.add_argument("--outdir", default=OUT)
ap.add_argument("--cap", type=int, default=150000)
args = ap.parse_args()
SEED = 20261011
TYPES = ["microglia", "astrocyte", "oligodendrocyte", "opc", "sst", "l23it", "pvalb"]
GRID = {"microglia": (0, 0), "astrocyte": (1, 0), "oligodendrocyte": (2, 0), "opc": (3, 0),
        "sst": (0, 1), "pvalb": (1, 1), "l23it": (2, 1)}
BOX, GAP = 20.0, 6.0
donors = pd.read_csv(DONOR_TABLE, index_col="donor_id")

cells = {l: pd.read_parquet(f"{args.indir}/{l}__all_cells.parquet") for l in TYPES}
n_micro = len(cells["microglia"])
n_each = (args.cap - n_micro) // 6
print("microglia", n_micro, "per other type", n_each, "total", n_micro + 6 * n_each, flush=True)

def allocate(counts, target):
    """largest-remainder proportional allocation, capped at counts"""
    if target >= counts.sum():
        return counts.copy()
    q = counts * target / counts.sum()
    base = np.floor(q).astype(int); rem = target - base.sum()
    order = np.argsort(-(q - base), kind="stable")
    base[order[:rem]] += 1
    return np.minimum(base, counts)

parts, samp_rows, umaps, pcas = [], [], [], []
for i, l in enumerate(TYPES):
    df = cells[l].reset_index(drop=True)
    target = n_micro if l == "microglia" else n_each
    rng = np.random.default_rng([SEED, i])
    cnt = df.groupby("donor_id").size()
    alloc = pd.Series(allocate(cnt.values, min(target, len(df))), index=cnt.index)
    sel = []
    for d, g in df.groupby("donor_id"):
        k = int(alloc[d])
        sel.append(g.index.values if k >= len(g) else rng.choice(g.index.values, k, replace=False))
        samp_rows.append({"subclass": SUBCLASSES[l][1], "donor_id": d, "n_available": int(len(g)), "n_selected": k})
    sel = np.sort(np.concatenate(sel))
    assert len(sel) == min(target, len(df)), (l, len(sel), target)
    sub = df.loc[sel].copy(); sub["subclass_label"] = l
    # X_scVI for the selected cells (read from the h5ad, matched by obs_name)
    a = ad.read_h5ad(f"{H5AD_DIR}/{SUBCLASSES[l][0]}", backed="r")
    pos = a.obs_names.get_indexer(sub["obs_name"].values)
    assert (pos >= 0).all()
    Z = np.asarray(a.obsm["X_scVI"], dtype=np.float32)[pos]
    a.file.close()
    # native UMAP scaled into a BOX x BOX square (aspect preserved) and moved to its grid slot
    U = sub[["umap_1", "umap_2"]].values.astype(np.float32)
    lo, hi = cells[l][["umap_1", "umap_2"]].min().values, cells[l][["umap_1", "umap_2"]].max().values
    sc_ = BOX / float((hi - lo).max())
    gx, gy = GRID[l]
    G = (U - lo) * sc_ + np.array([gx * (BOX + GAP), -gy * (BOX + GAP)], dtype=np.float32)
    parts.append(sub); umaps.append((G.astype(np.float32), U)); pcas.append(Z)
    print(l, "selected", len(sub), "of", len(df), flush=True)

obs = pd.concat(parts, ignore_index=True)
Ug = np.vstack([u[0] for u in umaps]); Un = np.vstack([u[1] for u in umaps]); Zp = np.vstack(pcas)
D = donors.reindex(obs["donor_id"].values)
out = pd.DataFrame(index=pd.Index(obs["obs_name"].astype(str)))
out["subclass"] = obs["subclass_label"].map({k: v[1] for k, v in SUBCLASSES.items()}).values
out["donor_id"] = obs["donor_id"].values
out["supertype"] = obs["Supertype"].values
codes = {k: i for i, k in enumerate(sorted(set(zip(out["subclass"], out["supertype"]))))}
out["leiden"] = [codes[k] for k in zip(out["subclass"], out["supertype"])]
out["CPS"] = D["CPS_Local"].values.astype(np.float32)
out["CPS_ABeta"] = D["CPS_Local_ABeta"].values.astype(np.float32)
out["CPS_pTau"] = D["CPS_Local_pTau"].values.astype(np.float32)
for ph, short in [("CPS_Local", "CPS"), ("CPS_Local_ABeta", "ABeta"), ("CPS_Local_pTau", "pTau")]:
    nc = obs[f"ncorr_{ph}"].values.astype(np.float32); fd = obs[f"fdr_{ph}"].values.astype(np.float32)
    out[f"ncorr_{short}"] = nc
    out[f"fdr_{short}"] = fd
    for lev, tag in [(0.10, "fdr10"), (0.05, "fdr5")]:
        flag = np.where(np.isnan(nc) | (fd > lev), 0, np.sign(nc)).astype(np.int8)
        out[f"redblue_{tag}_{short}"] = flag
out["apoe"] = D["APOE Genotype"].astype(str).values
out["apoe4_carrier"] = D["apoe4_carrier"].astype(int).values
out["cognitive_status"] = D["Cognitive Status"].astype(str).values
out["dementia"] = D["dementia"].astype(int).values
out["adnc"] = D["Overall AD neuropathological Change"].astype(str).values
out["braak"] = D["Braak"].astype(str).values
out["age_at_death"] = D["age"].astype(np.float32).values
out["sex"] = D["Sex"].astype(str).values
out["severe_donor"] = D["severe"].astype(int).values

a = ad.AnnData(X=sp.csr_matrix((len(out), 1), dtype=np.float32), obs=out)
a.var_names = pd.Index(["placeholder"])
a.obsm["X_umap"] = Ug; a.obsm["X_umap_subclass"] = Un; a.obsm["X_pca"] = Zp
a.uns["description"] = ("Broken Bricks floor 1 (cells): SEA-AD MTG, aged donors only. X is a placeholder. X_pca is the per-subclass scVI latent "
                        "(10-d); X_umap is the SEA-AD native per-subclass UMAP placed on a grid.")
a.uns["leiden_codes"] = json.dumps({f"{k[0]}|{k[1]}": v for k, v in codes.items()})
os.makedirs(args.outdir, exist_ok=True)
path = f"{args.outdir}/floor1_cells.h5ad"
a.write_h5ad(path, compression="gzip")
pd.DataFrame(samp_rows).to_csv(f"{args.outdir}/floor1_sampling.csv", index=False)
summ = {"n_cells": int(len(out)), "per_subclass": out["subclass"].value_counts().to_dict(), "n_donors": int(out["donor_id"].nunique()),
        "n_supertypes": int(len(codes)), "cap": args.cap, "n_each_large": int(n_each), "seed": SEED,
        "file_bytes": os.path.getsize(path), "obs_columns": list(out.columns), "obsm": list(a.obsm.keys())}
json.dump(summ, open(f"{args.outdir}/floor1_cells_summary.json", "w"), indent=1)
print(json.dumps(summ, indent=1))
