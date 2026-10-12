#!/usr/bin/env python
"""Collect MAGMA gene-property and set-analysis outputs into tidy CSVs (runs on ORCD or locally; pandas only).

Input : <MAGMA>/<tag>/celltype.gsa.out, <MAGMA>/<tag>/<label>.{marginal,conditional}.gsa.out
Output: results/day3/h3_celltype_property.csv, h3_program_sets.csv, h3_summary.json  (written next to the inputs' parent: <ROOT>/results/day3)
Decision rules (amendment A7 items 4-5):
  positive control: Microglia gene-property one-sided P < 1e-4
  open question   : Microglia_red conditional on Microglia specificity, one-sided P < 0.05
  prediction      : no neuronal primary program (14 sets) has Holm-adjusted conditional P < 0.05
"""
import json, os, sys
import numpy as np, pandas as pd
from statsmodels.stats.multitest import multipletests
R = sys.argv[1] if len(sys.argv) > 1 else "/home/pranavak/orcd/pool/mantis-films/research/seaad/results/day3"
M = f"{R}/magma"
COVAR = {"microglia": "Microglia", "astrocyte": "Astrocyte", "oligodendrocyte": "Oligodendrocyte", "opc": "OPC", "sst": "Sst", "pvalb": "Pvalb",
         "vip": "Vip", "lamp5": "Lamp5", "l23it": "L23IT", "l4it": "L4IT", "l5it": "L5IT"}
NEURONAL = ["Sst", "Pvalb", "Vip", "Lamp5", "L23IT", "L4IT", "L5IT"]

def read_gsa(path):
    meta = []; hdr = None; rows = []
    for line in open(path):
        if line.startswith("#"):
            meta.append(line.strip("# \n")); continue
        parts = line.split()
        if not parts: continue
        if hdr is None: hdr = parts; continue
        rows.append(parts[:len(hdr)])
    df = pd.DataFrame(rows, columns=hdr)
    for c in ["NGENES", "BETA", "BETA_STD", "SE", "P"]:
        if c in df: df[c] = pd.to_numeric(df[c])
    return df, meta

prop, sets = [], []
tags = [t for t in ["primary", "neff", "noapoe"] if os.path.isdir(f"{M}/{t}")]
for tag in tags:
    f = f"{M}/{tag}/celltype.gsa.out"
    if os.path.exists(f):
        df, meta = read_gsa(f); df.insert(0, "tag", tag); df["neglog10_p"] = -np.log10(df["P"]); prop.append(df)
        direction = [m for m in meta if "DIRECTION" in m.upper() or "ONE-SIDED" in m.upper()]
    for lab, cv in COVAR.items():
        for kind in ["marginal", "conditional"]:
            f = f"{M}/{tag}/{lab}.{kind}.gsa.out"
            if not os.path.exists(f): continue
            df, meta = read_gsa(f)
            df = df[df["TYPE"].astype(str).str.upper() == "SET"].copy() if "TYPE" in df else df
            df.insert(0, "tag", tag); df.insert(1, "subclass", cv); df.insert(2, "analysis", kind)
            df["neglog10_p"] = -np.log10(df["P"]); sets.append(df)
prop = pd.concat(prop) if prop else pd.DataFrame(); sets = pd.concat(sets) if sets else pd.DataFrame()
if len(sets):
    sets["color"] = sets["VARIABLE"].str.extract(r"_(red|blue)")[0]
    sets["kind"] = sets["VARIABLE"].str.replace(r"^[A-Za-z0-9]+_(red|blue)_?", "", regex=True).replace("", "primary")
    sets["neuronal"] = sets["subclass"].isin(NEURONAL)
    sets["holm_p_within_neuronal_primary"] = np.nan
    for (tag, an), g in sets[(sets.kind == "primary") & sets.neuronal].groupby(["tag", "analysis"]):
        sets.loc[g.index, "holm_p_within_neuronal_primary"] = multipletests(g["P"].values, method="holm")[1]
    sets = sets.reset_index(drop=True)
prop.to_csv(f"{R}/h3_celltype_property.csv", index=False); sets.to_csv(f"{R}/h3_program_sets.csv", index=False)

S = {"tags": tags}
if len(prop):
    for tag in tags:
        p = prop[(prop.tag == tag)]
        mg = p[p.VARIABLE == "Microglia"]
        if len(mg): S[f"{tag}_microglia_property"] = {"beta": float(mg.BETA.iloc[0]), "se": float(mg.SE.iloc[0]), "p": float(mg.P.iloc[0]), "ngenes": int(mg.NGENES.iloc[0])}
    S["positive_control_pass"] = bool(S.get("primary_microglia_property", {}).get("p", 1) < 1e-4)
if len(sets):
    for tag in tags:
        q = sets[(sets.tag == tag) & (sets.kind == "primary")]
        r = q[(q.VARIABLE == "Microglia_red") & (q.analysis == "conditional")]
        m = q[(q.VARIABLE == "Microglia_red") & (q.analysis == "marginal")]
        if len(r): S[f"{tag}_microglia_red_conditional"] = {"beta": float(r.BETA.iloc[0]), "se": float(r.SE.iloc[0]), "p": float(r.P.iloc[0]), "ngenes_in_set": int(r.NGENES.iloc[0])}
        if len(m): S[f"{tag}_microglia_red_marginal"] = {"beta": float(m.BETA.iloc[0]), "se": float(m.SE.iloc[0]), "p": float(m.P.iloc[0])}
        neu = q[q.neuronal & (q.analysis == "conditional")]
        if len(neu):
            S[f"{tag}_neuronal_min_conditional_p"] = float(neu.P.min()); S[f"{tag}_neuronal_min_holm_p"] = float(neu.holm_p_within_neuronal_primary.min())
            S[f"{tag}_neuronal_nominal_hits"] = neu.loc[neu.P < 0.05, "VARIABLE"].tolist()
            S[f"{tag}_prediction_holds"] = bool(neu.holm_p_within_neuronal_primary.min() >= 0.05)
json.dump(S, open(f"{R}/h3_summary.json", "w"), indent=1)
print(json.dumps(S, indent=1))
