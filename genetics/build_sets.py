#!/usr/bin/env python
"""Day 3: turn the per-subclass program tables into gene sets (A7 items 1-2).

Reads <ROOT>/results/day3/programs/<label>__<variant>_genes.csv, maps genes to NCBI Gene IDs (genemap.py), and writes
  programs_mapped/<label>__<variant>_genes_mapped.csv       tables with entrez IDs
  sets/<label>.sets.txt                                    MAGMA long format (entrez <tab> set), all kinds of set for that subclass
  programs.gmt                                             symbols, all sets (deliverable)
  program_membership.csv, program_summary.csv              long format with rho; per set: sizes, mapped counts, overlap with the primary program
Set names: <Covar>_red / _blue (primary, donor-blocked), _red_wd/_blue_wd (within-donor), _red_pc/_blue_pc (protein-coding), _red_nosad/_blue_nosad.
"""
import os, sys, json
import numpy as np, pandas as pd
sys.path.insert(0, "/home/pranavak/orcd/pool/mantis-films/research/broken-bricks/genetics")
from genemap import GeneMap
R = "/home/pranavak/orcd/pool/mantis-films/research/seaad/results/day3"
COVAR = {"microglia": "Microglia", "astrocyte": "Astrocyte", "oligodendrocyte": "Oligodendrocyte", "opc": "OPC", "sst": "Sst", "pvalb": "Pvalb",
         "vip": "Vip", "lamp5": "Lamp5", "l23it": "L23IT", "l4it": "L4IT", "l5it": "L5IT"}
os.makedirs(f"{R}/sets", exist_ok=True); os.makedirs(f"{R}/programs_mapped", exist_ok=True)
gm = GeneMap()
long_rows, gmt, summ = [], [], []
for lab, cv in COVAR.items():
    sets = {}                                                     # name -> DataFrame of member rows
    for variant in ["all", "nosad"]:
        f = f"{R}/programs/{lab}__{variant}_genes.csv"
        if not os.path.exists(f):
            print("missing", f); continue
        df = gm.map_table(pd.read_csv(f))
        df.to_csv(f"{R}/programs_mapped/{lab}__{variant}_genes_mapped.csv", index=False)
        if variant == "all":
            for col, suf, rho in [("program", "", "rho_donor_blocked"), ("program_within_donor", "_wd", "rho_within_donor"), ("program_protein_coding", "_pc", "rho_donor_blocked")]:
                for color in ["red", "blue"]:
                    m = df[df[col] == color]
                    sets[f"{cv}_{color}{suf}"] = m.assign(rho=m[rho])
        else:
            for color in ["red", "blue"]:
                m = df[df["program"] == color]
                sets[f"{cv}_{color}_nosad"] = m.assign(rho=m["rho_donor_blocked"])
    lines = []
    for name, m in sets.items():
        ent = m.loc[m["in_magma_loc"], "entrez"].astype(int).drop_duplicates().tolist()
        lines += [f"{e}\t{name}" for e in ent]
        gmt.append("\t".join([name, f"{cv} program, {len(m)} genes ({len(ent)} in MAGMA gene.loc)"] + m["symbol"].tolist()))
        for r in m.itertuples():
            long_rows.append({"set": name, "subclass": cv, "ensembl_id": r.ensembl_id, "symbol": r.symbol, "entrez": r.entrez, "in_magma_loc": r.in_magma_loc, "rho": r.rho})
        base = sets.get(f"{cv}_{name.split('_')[1]}")
        summ.append({"set": name, "subclass": cv, "color": name.split("_")[1], "kind": "primary" if name.count("_") == 1 else name.split("_", 2)[2],
                     "n_genes": len(m), "n_entrez": int(m["entrez"].notna().sum()), "n_in_magma_loc": len(ent),
                     "jaccard_with_primary": float(len(set(m.ensembl_id) & set(base.ensembl_id)) / len(set(m.ensembl_id) | set(base.ensembl_id))) if base is not None else np.nan,
                     "n_overlap_with_primary": int(len(set(m.ensembl_id) & set(base.ensembl_id))) if base is not None else np.nan,
                     "rho_min": float(m.rho.min()) if len(m) else np.nan, "rho_max": float(m.rho.max()) if len(m) else np.nan})
    open(f"{R}/sets/{lab}.sets.txt", "w").write("\n".join(lines) + "\n")
    print(lab, {k: len(v) for k, v in sets.items()}, flush=True)
open(f"{R}/programs.gmt", "w").write("\n".join(gmt) + "\n")
pd.DataFrame(long_rows).to_csv(f"{R}/program_membership.csv", index=False)
pd.DataFrame(summ).to_csv(f"{R}/program_summary.csv", index=False)
apoe = sorted(gm.apoe_region_genes())
pd.DataFrame({"entrez": apoe, "symbol": [gm.loc.loc[a, "symbol"] for a in apoe]}).to_csv(f"{R}/apoe_region_genes_b37.csv", index=False)
print("APOE-region genes (chr19:44-46 Mb, b37):", len(apoe))
