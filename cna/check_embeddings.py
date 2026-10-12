#!/usr/bin/env python
"""Are X_scVI / X_umap of the subclass files in one shared coordinate system, or per-subclass?  (decides how the cells floor is built)"""
import sys, os, json
import numpy as np, anndata as ad
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *
rows = {}
for label, (h5, disp) in SUBCLASSES.items():
    a = ad.read_h5ad(f"{H5AD_DIR}/{h5}", backed="r")
    Z = np.asarray(a.obsm["X_scVI"]); U = np.asarray(a.obsm["X_umap"])
    rows[label] = {"n": int(Z.shape[0]), "scvi_mean": Z.mean(0).round(2).tolist(), "scvi_std": Z.std(0).round(2).tolist(),
                   "umap_min": U.min(0).round(1).tolist(), "umap_max": U.max(0).round(1).tolist(), "umap_mean": U.mean(0).round(1).tolist()}
    a.file.close()
    print(label, json.dumps(rows[label]), flush=True)
# separability of subclasses in X_scVI: nearest-centroid accuracy on a 2,000-cell subsample per subclass
rng = np.random.default_rng(0)
samp = {}
for label, (h5, disp) in SUBCLASSES.items():
    a = ad.read_h5ad(f"{H5AD_DIR}/{h5}", backed="r")
    idx = np.sort(rng.choice(a.n_obs, 2000, replace=False))
    samp[label] = np.asarray(a.obsm["X_scVI"])[idx]; a.file.close()
cent = {k: v.mean(0) for k, v in samp.items()}
keys = list(cent); C = np.stack([cent[k] for k in keys])
acc = {k: float(np.mean([keys[int(np.argmin(((x - C) ** 2).sum(1)))] == k for x in samp[k]])) for k in keys}
print("nearest-centroid accuracy per subclass (chance ~0.09):", json.dumps(acc))
json.dump({"stats": rows, "nearest_centroid_acc": acc}, open(f"{OUT}/check_embeddings.json", "w"), indent=1)
