#!/usr/bin/env python
"""Broken Bricks Day 1: assemble CNA results; decide PC1 and H1 exactly as defined in preregistration.md + amendment A4.

Reads  results/day1/<label>__<variant>.json  (+ _cells.parquet)  and the SEA-AD scCODA outputs.
Writes (to --outdir):
  h1_global_p.csv               one row per subclass: global P, k, r2, FDR cell counts, for all three phenotypes and variants
  cna_supertype_means.csv       subclass x supertype x phenotype (primary run): mean ncorr, FDR-flagged cell counts
  pc1_sccoda_concordance.csv    credible SEA-AD supertypes: sign of scCODA effect vs sign of CNA supertype-mean ncorr
  pc1_summary.json              pass/fail of every PC1 criterion and of H1, with the numbers
"""
import argparse, json, os, sys
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, binomtest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *

ap = argparse.ArgumentParser()
ap.add_argument("--indir", default=OUT); ap.add_argument("--outdir", default=OUT)
args = ap.parse_args()
os.makedirs(args.outdir, exist_ok=True)

def load(label, variant):
    p = f"{args.indir}/{label}__{variant}.json"
    if not os.path.exists(p):
        return None
    j = json.load(open(p))
    return None if "identical_to" in j else j

# ---------------- H1 table ----------------
rows, st_rows = [], []
for label, (h5, disp) in SUBCLASSES.items():
    prim = load(label, "all")
    if prim is None:
        print("MISSING primary run for", label); continue
    row = {"subclass": disp, "label": label, "n_donors": prim["n_donors"], "n_cells": prim["n_cells_used"]}
    for ph in PHENOS:
        r = prim["results"][ph]
        short = {"CPS_Local": "CPS", "CPS_Local_ABeta": "ABeta", "CPS_Local_pTau": "pTau"}[ph]
        row.update({f"P_{short}": r["global_p"], f"k_{short}": r["k_pcs"], f"r2_{short}": r["r2"],
                    f"cells_fdr10_{short}": r["n_cells_fdr10"], f"pos_fdr10_{short}": r["pos_fdr10"], f"neg_fdr10_{short}": r["neg_fdr10"],
                    f"cells_fdr5_{short}": r["n_cells_fdr5"], f"nan_cells_{short}": r["n_cells_ncorr_nan"],
                    f"mean_ncorr_{short}": r["ncorr_mean"]})
        for s in r["supertypes"]:
            st_rows.append({"subclass": disp, "label": label, "phenotype": ph, **s})
    for v in ["nosad", "ge50"]:
        j = load(label, v)
        if j is None:
            row[f"P_CPS_{v}"] = prim["results"]["CPS_Local"]["global_p"] if v == "ge50" else np.nan
            row[f"donors_{v}"] = prim["n_donors"] if v == "ge50" else np.nan
            continue
        r = j["results"]["CPS_Local"]
        row[f"P_CPS_{v}"] = r["global_p"]; row[f"donors_{v}"] = j["n_donors"]; row[f"k_CPS_{v}"] = r["k_pcs"]; row[f"r2_CPS_{v}"] = r["r2"]
        row[f"pos_fdr10_CPS_{v}"] = r["pos_fdr10"]; row[f"neg_fdr10_CPS_{v}"] = r["neg_fdr10"]
        if v == "nosad":
            for ph, short in [("CPS_Local_ABeta", "ABeta"), ("CPS_Local_pTau", "pTau")]:
                if ph in j["results"]:
                    row[f"P_{short}_nosad"] = j["results"][ph]["global_p"]
    rows.append(row)
h1 = pd.DataFrame(rows)
st = pd.DataFrame(st_rows)
h1.to_csv(f"{args.outdir}/h1_global_p.csv", index=False)
st.to_csv(f"{args.outdir}/cna_supertype_means.csv", index=False)

out = {"n_subclasses_analysed": int(len(h1))}
alpha = 0.05
out["H1"] = {
    "criterion": ">= 4 of the analysed subclasses with global CNA P < 0.05 for CPS_Local (primary, all 84 donors)",
    "n_analysed": int(len(h1)),
    "n_significant": int((h1["P_CPS"] < alpha).sum()),
    "significant": h1.loc[h1["P_CPS"] < alpha, "subclass"].tolist(),
    "not_significant": h1.loc[h1["P_CPS"] >= alpha, "subclass"].tolist(),
    "pass": bool((h1["P_CPS"] < alpha).sum() >= 4),
    "sensitivity_without_11_severely_affected": {
        "n_significant": int((h1["P_CPS_nosad"] < alpha).sum()) if "P_CPS_nosad" in h1 else None,
        "n_analysed": int(h1["P_CPS_nosad"].notna().sum()) if "P_CPS_nosad" in h1 else None,
        "significant": h1.loc[h1["P_CPS_nosad"] < alpha, "subclass"].tolist() if "P_CPS_nosad" in h1 else None},
    "sensitivity_donors_ge50": {"n_significant": int((h1["P_CPS_ge50"] < alpha).sum()) if "P_CPS_ge50" in h1 else None},
    "ABeta_n_significant": int((h1["P_ABeta"] < alpha).sum()), "pTau_n_significant": int((h1["P_pTau"] < alpha).sum()),
    "p_floor": 1.0 / (NNULL + 1),
}
# descriptive multiplicity: BH over the 11 CPS_Local P
from statsmodels.stats.multitest import multipletests
h1["q_BH_CPS"] = multipletests(h1["P_CPS"], method="fdr_bh")[1]
out["H1"]["n_significant_BH_q_lt_0.05_descriptive"] = int((h1["q_BH_CPS"] < 0.05).sum())
h1.to_csv(f"{args.outdir}/h1_global_p.csv", index=False)

# ---------------- PC1 ----------------
pc1 = {}
# criterion 1: global P for Microglia-PVM and Sst
def gp(label): return load(label, "all")["results"]["CPS_Local"]
mic, sst = gp("microglia"), gp("sst")
pc1["crit1_global_P"] = {"microglia_P": mic["global_p"], "sst_P": sst["global_p"],
                         "pass": bool(mic["global_p"] < 0.05 and sst["global_p"] < 0.05)}
# criterion 2: Sst depletion (A4.4b)
sst_cells = pd.read_parquet(f"{args.indir}/sst__all_cells.parquet", columns=["ncorr_CPS_Local", "fdr_CPS_Local", "Supertype"])
nc, fd = sst_cells["ncorr_CPS_Local"], sst_cells["fdr_CPS_Local"]
pos10 = int(((nc > 0) & (fd <= 0.10)).sum()); neg10 = int(((nc < 0) & (fd <= 0.10)).sum())
pos5 = int(((nc > 0) & (fd <= 0.05)).sum()); neg5 = int(((nc < 0) & (fd <= 0.05)).sum())
pc1["crit2_sst_depleted"] = {"mean_ncorr_all_sst_cells": float(nc.mean()), "fdr10_neg_cells": neg10, "fdr10_pos_cells": pos10,
                             "fdr5_neg_cells": neg5, "fdr5_pos_cells": pos5,
                             "share_negative_among_fdr10": neg10 / max(1, neg10 + pos10),
                             "pass": bool(nc.mean() < 0 and neg10 > pos10)}
# microglia expansion (descriptive, A4.4e)
mc = pd.read_parquet(f"{args.indir}/microglia__all_cells.parquet", columns=["ncorr_CPS_Local", "fdr_CPS_Local", "Supertype"])
g = mc[mc["Supertype"] == "Micro-PVM_3-SEAAD"]
p10 = int(((g.ncorr_CPS_Local > 0) & (g.fdr_CPS_Local <= 0.10)).sum()); n10 = int(((g.ncorr_CPS_Local < 0) & (g.fdr_CPS_Local <= 0.10)).sum())
pc1["microglia_expansion_descriptive"] = {"supertype": "Micro-PVM_3-SEAAD", "n_cells": int(len(g)), "mean_ncorr": float(g.ncorr_CPS_Local.mean()),
                                          "fdr10_pos": p10, "fdr10_neg": n10, "holds": bool(g.ncorr_CPS_Local.mean() > 0 and p10 > n10),
                                          "all_microglia_fdr10_pos": int(((mc.ncorr_CPS_Local > 0) & (mc.fdr_CPS_Local <= 0.10)).sum()),
                                          "all_microglia_fdr10_neg": int(((mc.ncorr_CPS_Local < 0) & (mc.fdr_CPS_Local <= 0.10)).sum())}

# criterion 3: sign agreement with SEA-AD scCODA
summ = pd.read_csv(SCCODA_SUMMARY)
summ = summ[summ["Region"] == "MTG"].set_index("Cell Type")["Local Model"]
cna_st = (st[st["phenotype"] == "CPS_Local"].set_index("supertype")[["subclass", "n_cells", "mean_ncorr", "n_pos_fdr10", "n_neg_fdr10"]])
assert cna_st.index.is_unique, "supertype names are not unique across subclasses"

# own aggregation of the per-run results (A4.4d)
def own_aggregate():
    parts = []
    for fn in ["Neuronal: Glutamatergic Neuronal: GABAergic_Supertype_results.csv", "Non-neuronal and Non-neural_Supertype_results.csv"]:
        parts.append(pd.read_csv(f"{SCCODA_DIR}/{fn}", low_memory=False))
    m = pd.concat(parts)
    m = m[(m["Region"] == "MTG") & (m["Covariate"] == "CPS_Local")].copy()
    m["is_ref"] = m["Cell Type"] == m["Reference Cell Type"]
    m["cred"] = False
    for ref, gg in m.groupby("Reference Cell Type"):
        d = gg[~gg.is_ref].sort_values("Inclusion probability", ascending=False)
        p = d["Inclusion probability"].values
        ok = (np.cumsum(1 - p) / np.arange(1, len(p) + 1)) <= 0.05
        n = 0 if not ok.any() else int(np.where(ok)[0].max()) + 1
        m.loc[d.index[:n], "cred"] = True
    a = m[~m.is_ref].groupby("Cell Type").agg(own_median=("Final Parameter", "median"), own_frac_cred=("cred", "mean"), own_runs=("cred", "size"))
    return a
own = own_aggregate()

conc = pd.DataFrame(index=cna_st.index)
conc = cna_st.join(summ.rename("sea_ad_local_model")).join(own)
conc["in_scCODA_summary"] = conc["sea_ad_local_model"].notna()
conc["credible_primary"] = conc["sea_ad_local_model"].fillna(0) != 0
conc["credible_own"] = conc["own_frac_cred"].fillna(0) >= 0.5
conc["sign_cna"] = np.sign(conc["mean_ncorr"])
conc["agree_primary"] = np.where(conc["credible_primary"], np.sign(conc["sea_ad_local_model"]) == conc["sign_cna"], np.nan)
conc["agree_own"] = np.where(conc["credible_own"], np.sign(conc["own_median"]) == conc["sign_cna"], np.nan)
conc.reset_index().rename(columns={"index": "supertype"}).to_csv(f"{args.outdir}/pc1_sccoda_concordance.csv", index=False)

def summarize(flag, agree, effect):
    c = conc[conc[flag]]
    k = int(c[agree].astype(float).sum()); n = int(len(c))
    d = {"n_credible_matched": n, "n_agree": k, "agreement": (k / n) if n else None,
         "pass": bool(n > 0 and k / n >= 0.70),
         "binom_p_vs_50pct": float(binomtest(k, n, 0.5).pvalue) if n else None,
         "per_subclass": {s: [int(gg[agree].astype(float).sum()), int(len(gg))] for s, gg in c.groupby("subclass")}}
    cc = c.dropna(subset=[effect, "mean_ncorr"])
    d["spearman_credible"] = float(spearmanr(cc["mean_ncorr"], cc[effect])[0]) if len(cc) > 3 else None
    allm = conc.dropna(subset=[effect, "mean_ncorr"])
    d["spearman_all_modelled"] = float(spearmanr(allm["mean_ncorr"], allm[effect])[0]) if len(allm) > 3 else None
    d["n_all_modelled"] = int(len(allm))
    d["disagree"] = c.index[~c[agree].astype(bool)].tolist()
    return d
pc1["crit3_sign_agreement_primary_summary_table"] = summarize("credible_primary", "agree_primary", "sea_ad_local_model")
pc1["crit3_sign_agreement_sensitivity_own_aggregation"] = summarize("credible_own", "agree_own", "own_median")
# credible SEA-AD supertypes that could not be matched to a CNA supertype of the 11 subclasses
unmatched = summ[(summ.notna()) & (summ != 0) & (~summ.index.isin(conc.index))]
pc1["credible_summary_supertypes_not_in_the_11_subclasses"] = int(len(unmatched))
pc1["credible_summary_supertypes_unmatched_names"] = unmatched.index.tolist()

pc1["PC1_pass_primary"] = bool(pc1["crit1_global_P"]["pass"] and pc1["crit2_sst_depleted"]["pass"] and
                               pc1["crit3_sign_agreement_primary_summary_table"]["pass"])
pc1["PC1_pass_sensitivity"] = bool(pc1["crit1_global_P"]["pass"] and pc1["crit2_sst_depleted"]["pass"] and
                                   pc1["crit3_sign_agreement_sensitivity_own_aggregation"]["pass"])
out["PC1"] = pc1
json.dump(out, open(f"{args.outdir}/pc1_summary.json", "w"), indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
print(json.dumps(out, indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))[:6000])
