#!/usr/bin/env python
"""Day 3, A7 item 7: overlap of the programs with Open Targets genes (genetic_association score >= 0.5 for MONDO_0004975).

For each program (primary donor-blocked red/blue, 200 genes) of each subclass: overlap count with the OT set among program genes that are in the
subclass's expressed universe; null = 10,000 draws of an equally sized gene set matched on mean log-normalised expression (20 quantile bins of the
subclass's expressed genes, sampled without replacement within each bin, program genes excluded); one-sided empirical P = (1 + n_null >= obs) / 10,001.
Repeated without the genes in the APOE region (chr19:44-46 Mb, b37).  The OT gene list comes from the literature study's results/gene_master.csv
(G_genetic >= 0.5; 173 genes after amendment A3) and is read from <ROOT>/results/day3/ot_genes_G05.csv (gene_id = NCBI Gene ID, ensembl_id).
"""
import os, sys, json
import numpy as np, pandas as pd
R = "/home/pranavak/orcd/pool/mantis-films/research/seaad/results/day3"
SEED = 20261011; NNULL = 10000; NBIN = 20
LABELS = ["microglia", "astrocyte", "oligodendrocyte", "opc", "sst", "pvalb", "vip", "lamp5", "l23it", "l4it", "l5it"]
ot = pd.read_csv(f"{R}/ot_genes_G05.csv")
apoe = set(pd.read_csv(f"{R}/apoe_region_genes_b37.csv")["entrez"].astype(int))
ot_ent = set(ot["gene_id"].astype(int))
print("OT genes G>=0.5:", len(ot_ent), flush=True)
rows = []
for lab in LABELS:
    df = pd.read_csv(f"{R}/programs_mapped/{lab}__all_genes_mapped.csv")
    df = df[df["entrez"].notna()].copy(); df["entrez"] = df["entrez"].astype(int)
    df = df.sort_values("detect_frac", ascending=False).drop_duplicates("entrez")
    for apoe_ex in [False, True]:
        d = df[~df["entrez"].isin(apoe)] if apoe_ex else df
        d = d.reset_index(drop=True)
        d["bin"] = pd.qcut(d["mean_logexpr"].rank(method="first"), NBIN, labels=False)
        isot = d["entrez"].isin(ot_ent).values
        for color in ["red", "blue"]:
            prog = (d["program"] == color).values
            n_prog = int(prog.sum()); obs = int((prog & isot).sum())
            rng = np.random.default_rng(SEED + LABELS.index(lab) * 10 + (3 if apoe_ex else 0) + (1 if color == "blue" else 0))
            cnt_bin = d.loc[prog, "bin"].value_counts().to_dict()
            pools = {b: np.where((d["bin"].values == b) & ~prog)[0] for b in cnt_bin}
            null = np.empty(NNULL, dtype=int)
            for i in range(NNULL):
                s = 0
                for b, k in cnt_bin.items():
                    p = pools[b]; s += isot[rng.choice(p, size=min(k, len(p)), replace=False)].sum()
                null[i] = s
            rows.append({"subclass": lab, "program": color, "exclude_APOE_region": apoe_ex, "n_program_genes_mapped": n_prog, "n_universe": len(d), "n_ot_in_universe": int(isot.sum()),
                         "overlap_obs": obs, "overlap_null_mean": float(null.mean()), "fold_enrichment": float(obs / max(null.mean(), 1e-9)),
                         "p_one_sided": float((1 + (null >= obs).sum()) / (NNULL + 1)),
                         "ot_genes_in_program": ",".join(d.loc[prog & isot, "symbol"].tolist())})
    print(lab, "done", flush=True)
pd.DataFrame(rows).to_csv(f"{R}/ot_overlap.csv", index=False)
print(pd.DataFrame(rows).drop(columns=["ot_genes_in_program"]).round(3).to_string())
