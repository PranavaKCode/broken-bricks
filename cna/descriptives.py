#!/usr/bin/env python
"""Day-1 descriptive tables from the primary CNA runs (EXPLORATORY: none of these is a preregistered test).

  abeta_ptau_map_agreement.csv   per subclass: cell-level Spearman between the ncorr maps of CPS_Local, ABeta, pTau; global P / r2 of each
  top_supertypes.csv             per subclass: the most expanded and most depleted supertype (mean ncorr, CPS_Local)
  nosad_vs_all_agreement.csv     per subclass: cell-level Spearman between the primary map and the map without the 11 severely affected donors
  donor_cell_counts.csv          per subclass: nuclei per donor (min / median / max) and donors under 50 nuclei
usage: python descriptives.py [--res results/cna]
"""
import argparse, json, os
import numpy as np, pandas as pd
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
ap = argparse.ArgumentParser(); ap.add_argument("--res", default=os.path.join(HERE, "..", "results", "cna"))
args = ap.parse_args()
LAB = {"microglia": "Microglia-PVM", "astrocyte": "Astrocyte", "oligodendrocyte": "Oligodendrocyte", "opc": "OPC", "sst": "Sst", "pvalb": "Pvalb",
       "vip": "Vip", "lamp5": "Lamp5", "l23it": "L2/3 IT", "l4it": "L4 IT", "l5it": "L5 IT"}
rows, tops, nos, cnt = [], [], [], []
for l, disp in LAB.items():
    j = json.load(open(f"{args.res}/{l}__all.json")); c = pd.read_parquet(f"{args.res}/{l}__all_cells.parquet")
    a, b, p = c["ncorr_CPS_Local"], c["ncorr_CPS_Local_ABeta"], c["ncorr_CPS_Local_pTau"]
    m = a.notna() & b.notna() & p.notna()
    r = j["results"]
    rows.append({"subclass": disp, "rho_CPS_vs_ABeta": spearmanr(a[m], b[m])[0], "rho_CPS_vs_pTau": spearmanr(a[m], p[m])[0], "rho_ABeta_vs_pTau": spearmanr(b[m], p[m])[0],
                 "P_ABeta": r["CPS_Local_ABeta"]["global_p"], "P_pTau": r["CPS_Local_pTau"]["global_p"], "r2_ABeta": r["CPS_Local_ABeta"]["r2"], "r2_pTau": r["CPS_Local_pTau"]["r2"],
                 "fdr10_cells_ABeta": r["CPS_Local_ABeta"]["n_cells_fdr10"], "fdr10_cells_pTau": r["CPS_Local_pTau"]["n_cells_fdr10"]})
    st = pd.DataFrame(r["CPS_Local"]["supertypes"]).sort_values("mean_ncorr")
    st = st[st["n_cells"] >= 200]
    for tag, rr in [("most_depleted", st.iloc[0]), ("most_expanded", st.iloc[-1])]:
        tops.append({"subclass": disp, "which": tag, "supertype": rr["supertype"], "n_cells": int(rr["n_cells"]), "mean_ncorr": rr["mean_ncorr"]})
    cd = pd.read_parquet(f"{args.res}/{l}__nosad_cells.parquet", columns=["obs_name", "ncorr_CPS_Local"]).set_index("obs_name")["ncorr_CPS_Local"]
    both = pd.concat([c.set_index("obs_name")["ncorr_CPS_Local"].rename("all"), cd.rename("nosad")], axis=1, join="inner").dropna()
    nos.append({"subclass": disp, "n_cells_shared": int(len(both)), "rho_all_vs_nosad": spearmanr(both["all"], both["nosad"])[0]})
    cpd = c.groupby("donor_id").size()
    cnt.append({"subclass": disp, "min": int(cpd.min()), "median": float(cpd.median()), "max": int(cpd.max()), "donors_lt50": int((cpd < 50).sum())})
pd.DataFrame(rows).round(4).to_csv(f"{args.res}/abeta_ptau_map_agreement.csv", index=False)
pd.DataFrame(tops).round(4).to_csv(f"{args.res}/top_supertypes.csv", index=False)
pd.DataFrame(nos).round(4).to_csv(f"{args.res}/nosad_vs_all_agreement.csv", index=False)
pd.DataFrame(cnt).to_csv(f"{args.res}/donor_cell_counts.csv", index=False)
pd.set_option("display.width", 220)
for f in ["abeta_ptau_map_agreement", "top_supertypes", "nosad_vs_all_agreement", "donor_cell_counts"]:
    print(f); print(pd.read_csv(f"{args.res}/{f}.csv").to_string(index=False))
