#!/usr/bin/env python
"""Final check 3 (part 1): per-donor nucleus counts by Class and Subclass from the SEA-AD MTG final-nuclei metadata (1.15 GB CSV).

Writes results/final_checks/donor_subclass_counts.csv (long format: donor_id, neurotypical_reference, Class, Subclass, used_in_analysis, n)
and class_value_counts.json. The regressions against CPS_Local are in class_share_regress.py (local, pandas + statsmodels).
"""
import json, os, sys
import pandas as pd

R = "/home/pranavak/orcd/pool/mantis-films/research/seaad"
F = f"{R}/data/s3/MTG/RNAseq/SEAAD_MTG_RNAseq_final-nuclei_metadata.2026-06-22.csv"
OUT = f"{R}/results/final_checks"; os.makedirs(OUT, exist_ok=True)
d = pd.read_csv(F, usecols=["Donor ID", "Neurotypical reference", "Class", "Subclass", "Used in analysis"],
                dtype={"Donor ID": "category", "Neurotypical reference": "category", "Class": "category", "Subclass": "category", "Used in analysis": "category"})
print(d.shape, flush=True)
print(d["Used in analysis"].value_counts(dropna=False).to_dict(), flush=True)
print(d["Class"].value_counts(dropna=False).to_dict(), flush=True)
print(d["Subclass"].value_counts(dropna=False).to_dict(), flush=True)
g = d.groupby(["Donor ID", "Neurotypical reference", "Class", "Subclass", "Used in analysis"], observed=True, dropna=False).size().rename("n").reset_index()
g.columns = ["donor_id", "neurotypical_reference", "Class", "Subclass", "used_in_analysis", "n"]
g.to_csv(f"{OUT}/donor_subclass_counts.csv", index=False)
json.dump({"n_rows": int(len(d)), "used_in_analysis": {str(k): int(v) for k, v in d["Used in analysis"].value_counts(dropna=False).items()},
           "class": {str(k): int(v) for k, v in d["Class"].value_counts(dropna=False).items()},
           "subclass": {str(k): int(v) for k, v in d["Subclass"].value_counts(dropna=False).items()},
           "n_donors": int(d["Donor ID"].nunique())}, open(f"{OUT}/class_value_counts.json", "w"), indent=1)
print("wrote", len(g), "rows", flush=True)
