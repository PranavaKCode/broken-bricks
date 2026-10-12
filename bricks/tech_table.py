#!/usr/bin/env python
"""Donor-level technical covariates for the Day-2 technical-covariate regression (A6 item 9).

Reads the SEA-AD MTG final-nuclei metadata (library-level fields are repeated on every nucleus row) and writes one row per donor:
  n_nuclei, rin_mean, neun_frac, frac_multiome_all (method), seq_saturation (nuclei-weighted mean of library-level
  GEX_sequencing_saturation over libraries that have a value), seq_sat_missing_frac (share of the donor's nuclei in libraries with no value),
  n_libraries, umis_median, genes_median
"""
import sys, pandas as pd, numpy as np
sys.path.insert(0, "/home/pranavak/orcd/pool/mantis-films/research/broken-bricks/cna")
from common import ROOT
f = f"{ROOT}/data/s3/MTG/RNAseq/SEAAD_MTG_RNAseq_final-nuclei_metadata.2026-06-22.csv"
cols = ["Donor ID", "library_prep", "method", "RIN", "NeuN positive fraction on FANS", "GEX_sequencing_saturation",
        "Number of UMIs", "Genes detected"]
d = pd.read_csv(f, usecols=cols)
print(d.shape, flush=True)
g = d.groupby("Donor ID")
out = pd.DataFrame({"n_nuclei": g.size(), "rin_mean": g["RIN"].mean(), "neun_frac": g["NeuN positive fraction on FANS"].mean(),
                    "frac_multiome_all": g["method"].apply(lambda s: (s.astype(str) == "10xMulti").mean()),
                    "umis_median": g["Number of UMIs"].median(), "genes_median": g["Genes detected"].median(),
                    "n_libraries": g["library_prep"].nunique()})
lib = d.groupby(["Donor ID", "library_prep"]).agg(n=("method", "size"), sat=("GEX_sequencing_saturation", "first")).reset_index()
def wmean(x):
    m = x["sat"].notna()
    return pd.Series({"seq_saturation": np.average(x.loc[m, "sat"], weights=x.loc[m, "n"]) if m.any() else np.nan,
                      "seq_sat_missing_frac": x.loc[~m, "n"].sum() / x["n"].sum()})
out = out.join(lib.groupby("Donor ID").apply(wmean, include_groups=False))
out.index.name = "donor_id"
out.to_csv(f"{ROOT}/results/day2/donor_tech.csv")
print(out.describe().T.to_string())
print("donors with no seq saturation:", int(out["seq_saturation"].isna().sum()))
