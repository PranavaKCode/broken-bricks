#!/usr/bin/env python
"""Mantis "donors" floor (protocol 06 section 4.3 Day 4, floor 2): 84 rows, one per aged SEA-AD MTG donor.

Input : results/bricks/all/donor_bricks.csv (tier columns from day2_bricks.py)
Output: results/bricks/donors_floor.csv  (also copied to floors/floor2_donors.csv, not committed) and results/bricks/donor_bricks.csv

Coordinates (A6 item 12, with the tier choice stated in A10): first two principal components of the matrix of residual z-scores Z of the
tier `allcols` (all 11 subclass columns; the preregistered filter keeps only three, all neuronal; missing = 0), with a UMAP of the same matrix
(n_neighbors 10, min_dist 0.3, seed 20261011) in `umap1`/`umap2`. vector_columns v0..v10 = the 11 brick scores B, z-scored over donors,
missing = column median (a raw embedding for Mantis; do not read distances in it as biology beyond what B means).
Labels: `discordance_label` (preregistered rule, tier allcols: all donors "not discordant"), `profile_label` (exploratory description of the
single shared residual axis and the glia-versus-neuron contrast; thresholds |z| > 1 on factor1 and |delta| > 1.5).
"""
import os, sys, json
import numpy as np, pandas as pd
import umap

HERE = os.path.dirname(os.path.abspath(__file__))
R = os.path.join(HERE, "..", "results", "bricks")
LAB = ["microglia", "astrocyte", "oligodendrocyte", "opc", "sst", "pvalb", "vip", "lamp5", "l23it", "l4it", "l5it"]
SEED = 20261011
d = pd.read_csv(os.path.join(R, "all", "donor_bricks.csv"), index_col="donor_id")
Z = d[[f"Z_{l}" for l in LAB]]
X = Z.fillna(0).values
Xc = X - X.mean(0)
U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
pc = U[:, :2] * S[:2]
sign = np.sign(Vt[0].sum()); pc[:, 0] *= sign if sign != 0 else 1       # PC1 positive = sicker than CPS predicts across cell types
um = umap.UMAP(n_neighbors=10, min_dist=0.3, random_state=SEED, n_jobs=1).fit_transform(Xc)
B = d[[f"B_{l}" for l in LAB]]
Bz = (B - B.mean()) / B.std(ddof=1); Bz = Bz.fillna(Bz.median())
f = d.copy()
f.insert(0, "floor_number", 2); f.insert(1, "title", f.index)
f.insert(2, "coordinate1", pc[:, 0]); f.insert(3, "coordinate2", pc[:, 1]); f.insert(4, "umap1", um[:, 0]); f.insert(5, "umap2", um[:, 1])
for j, l in enumerate(LAB): f[f"v{j}"] = Bz[f"B_{l}"].values
shared = f["factor1__allcols"]
delta = f["delta__allcols"]
prof = np.where(shared > 1, "sicker than CPS predicts in most cell types", np.where(shared < -1, "healthier than CPS predicts in most cell types", "near CPS expectation"))
con = np.where(delta > 1.5, "; glia relatively sicker than neurons", np.where(delta < -1.5, "; neurons relatively sicker than glia", ""))
f["profile_label"] = [a + b for a, b in zip(prof, con)]
f["discordance_label"] = f["discordance_label__allcols"]
f["discordant"] = f["discordant__allcols"]
cols_first = ["floor_number", "title", "coordinate1", "coordinate2", "umap1", "umap2"] + [f"v{j}" for j in range(len(LAB))]
f = f[cols_first + [c for c in f.columns if c not in cols_first]]
f.to_csv(os.path.join(R, "donors_floor.csv"))
d.to_csv(os.path.join(R, "donor_bricks.csv"))
os.makedirs(os.path.join(HERE, "..", "floors"), exist_ok=True)
f.to_csv(os.path.join(HERE, "..", "floors", "floor2_donors.csv"))
json.dump({"n_rows": int(len(f)), "pc_var_explained": [float(x) for x in S[:2] ** 2 / (S ** 2).sum()], "n_missing_scores_total": int(B.isna().sum().sum()),
           "profile_counts": f["profile_label"].value_counts().to_dict()}, open(os.path.join(R, "donors_floor_summary.json"), "w"), indent=1)
print(f.shape, f["profile_label"].value_counts().to_dict())
