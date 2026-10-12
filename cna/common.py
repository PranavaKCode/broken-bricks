"""Shared constants for the Broken Bricks Day-1 (CNA) analysis. Runs on ORCD with the `sc` conda env."""
ROOT = "/home/pranavak/orcd/pool/mantis-films/research/seaad"
H5AD_DIR = f"{ROOT}/data/cellxgene"
DONOR_TABLE = f"{ROOT}/results/seaad_mtg_donor_table.csv"
OUT = f"{ROOT}/results/day1"          # raw job outputs on ORCD (mirrored to results/cna/ in the repo)
SCCODA_SUMMARY = (f"{ROOT}/data/s3/Multiregion_2026/model_outputs/pertpy_compositional_modeling/"
                  "pertpy_summary_CPS_Local.20260622.csv")
SCCODA_DIR = f"{ROOT}/data/s3/Multiregion_2026/model_outputs/pertpy_compositional_modeling/CPS_Local"

# short label -> (h5ad file, display name).  The 11 non-vascular subclasses (preregistration C1).
SUBCLASSES = {
    "microglia":       ("Microglia-PVM_MTG.h5ad",   "Microglia-PVM"),
    "astrocyte":       ("Astrocyte_MTG.h5ad",       "Astrocyte"),
    "oligodendrocyte": ("Oligodendrocyte_MTG.h5ad", "Oligodendrocyte"),
    "opc":             ("OPC_MTG.h5ad",             "OPC"),
    "sst":             ("Sst_MTG.h5ad",             "Sst"),
    "pvalb":           ("Pvalb_MTG.h5ad",           "Pvalb"),
    "vip":             ("Vip_MTG.h5ad",             "Vip"),
    "lamp5":           ("Lamp5_MTG.h5ad",           "Lamp5"),
    "l23it":           ("L23IT_MTG.h5ad",           "L2/3 IT"),
    "l4it":            ("L4IT_MTG.h5ad",            "L4 IT"),
    "l5it":            ("L5IT_MTG.h5ad",            "L5 IT"),
}
PHENOS = ["CPS_Local", "CPS_Local_ABeta", "CPS_Local_pTau"]
NEIGHBORS = 30          # preregistered
NNULL = 10000           # global-P permutations (cna uses only the first 1,000 for the local FDR)
SEED = 0
