#!/usr/bin/env python
"""Day 3, A7 item 3 and 6: extra MAGMA input files.

 - bellenguez_rsid_p_neff.txt : SNP, P, N_eff = 4 / (1/n_cas + 1/n_con)  (effective-N sensitivity), same SNP set as the primary file
 - bellenguez_rsid_p_n.noAPOE.txt : the primary file without SNPs in chr19:44,000,000-46,000,000 (build 37 positions of g1000_eur.bim)
The primary file `bellenguez_rsid_p_n.txt` (N = n_cas + n_con) was built on Day 0.
"""
import time, pandas as pd, numpy as np
R = "/home/pranavak/orcd/pool/mantis-films/research/seaad"
t = time.time()
prim = pd.read_csv(f"{R}/results/bellenguez_rsid_p_n.txt", sep="\t")
print("primary file", prim.shape, prim.columns.tolist(), flush=True)
g = pd.read_csv(f"{R}/data/gwas/35379992-GCST90027158-MONDO_0004975.h.tsv.gz", sep="\t", usecols=["hm_rsid", "n_cas", "n_con"], low_memory=False)
g = g.dropna(subset=["hm_rsid"]).drop_duplicates("hm_rsid").set_index("hm_rsid")
g["neff"] = 4.0 / (1.0 / g["n_cas"] + 1.0 / g["n_con"])
prim["Neff"] = g["neff"].reindex(prim["SNP"]).values
print("Neff missing:", int(prim["Neff"].isna().sum()), " N summary", prim["N"].describe().round(0).to_dict(), " Neff summary", prim["Neff"].describe().round(0).to_dict(), flush=True)
prim.loc[prim.Neff.notna(), ["SNP", "P", "Neff"]].rename(columns={"Neff": "N"}).to_csv(f"{R}/results/bellenguez_rsid_p_neff.txt", sep="\t", index=False)
bim = pd.read_csv(f"{R}/data/magma/g1000_eur/g1000_eur.bim", sep="\t", header=None, names=["chr", "rsid", "cm", "pos", "a1", "a2"], usecols=[0, 1, 3])
ap = bim[(bim["chr"].astype(str) == "19") & bim["pos"].between(44_000_000, 46_000_000)]
print("reference SNPs in chr19:44-46Mb (b37):", len(ap), flush=True)
m = prim["SNP"].isin(set(ap["rsid"]))
print("GWAS SNPs removed:", int(m.sum()), flush=True)
prim.loc[~m, ["SNP", "P", "N"]].to_csv(f"{R}/results/bellenguez_rsid_p_n.noAPOE.txt", sep="\t", index=False)
print("done", round(time.time() - t, 1), "s")
