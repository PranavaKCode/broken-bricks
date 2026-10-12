"""Final check 2 (summary): CNA reruns with technical covariates (run_cna_tech.py).

Reads results/final_checks/cna_tech/<label>__all__<model>.json (+ _cells.parquet for the maps) and the Day-1 primary results results/cna/<label>__all.json.
Writes results/final_checks/cna_tech_summary.csv, cna_tech_supertypes.csv, cna_tech_map_agreement.csv and prints the tables used in 16-verification.md section 13.
"""
import glob, json, os
import numpy as np, pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
R = os.path.join(HERE, "..", "results")
T = f"{R}/final_checks/cna_tech"
LABELS = {"microglia": "Microglia-PVM", "astrocyte": "Astrocyte", "oligodendrocyte": "Oligodendrocyte", "opc": "OPC", "sst": "Sst", "pvalb": "Pvalb",
          "vip": "Vip", "lamp5": "Lamp5", "l23it": "L2/3 IT", "l4it": "L4 IT", "l5it": "L5 IT"}
rows = []
for lab, disp in LABELS.items():
    d1 = json.load(open(f"{R}/cna/{lab}__all.json"))["results"]["CPS_Local"]
    rows.append(dict(subclass=disp, label=lab, model="day1_primary", global_p=d1["global_p"], k=d1["k_pcs"], r2=d1["r2"], pos10=d1["pos_fdr10"], neg10=d1["neg_fdr10"],
                     n_cells_fdr10=d1["n_cells_fdr10"], covariates="age,sex,pmi,method_frac"))
    for m in ["base", "sat", "tech3", "neun", "rin", "nuc", "satperm"]:
        f = f"{T}/{lab}__all__{m}.json"
        if not os.path.exists(f):
            continue
        j = json.load(open(f)); r = j["result"]
        rows.append(dict(subclass=disp, label=lab, model=m, global_p=r["global_p"], k=r["k_pcs"], r2=r["r2"], pos10=r["pos_fdr10"], neg10=r["neg_fdr10"],
                         n_cells_fdr10=r["n_cells_fdr10"], covariates=",".join(j["covariates"])))
    # no-SAD variant (73 donors): Day-1 primary nosad run, rerun (base) and with saturation
    d1n = f"{R}/cna/{lab}__nosad.json"
    if os.path.exists(d1n):
        x = json.load(open(d1n))["results"]["CPS_Local"]
        rows.append(dict(subclass=disp, label=lab, model="nosad_day1", global_p=x["global_p"], k=x["k_pcs"], r2=x["r2"], pos10=x["pos_fdr10"], neg10=x["neg_fdr10"],
                         n_cells_fdr10=x["n_cells_fdr10"], covariates="age,sex,pmi,method_frac"))
    for m in ["base", "sat"]:
        f = f"{T}/{lab}__nosad__{m}.json"
        if not os.path.exists(f):
            continue
        j = json.load(open(f)); r = j["result"]
        rows.append(dict(subclass=disp, label=lab, model="nosad_" + m, global_p=r["global_p"], k=r["k_pcs"], r2=r["r2"], pos10=r["pos_fdr10"], neg10=r["neg_fdr10"],
                         n_cells_fdr10=r["n_cells_fdr10"], covariates=",".join(j["covariates"])))
S = pd.DataFrame(rows)
S.to_csv(f"{R}/final_checks/cna_tech_summary.csv", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
print("== global permutation P (10,000 permutations; floor 9.999e-5) and FDR-10% cells (positive/negative ncorr)")
piv = S.pivot(index="subclass", columns="model", values="global_p")[[c for c in ["day1_primary", "base", "sat", "tech3", "neun", "rin", "nuc", "satperm", "nosad_day1", "nosad_base", "nosad_sat"] if c in S.model.unique()]]
print(piv.reindex(list(LABELS.values())).to_string(float_format=lambda x: f"{x:.4g}"))
print()
for m in ["day1_primary", "sat", "tech3"]:
    if m in piv:
        v = piv[m].dropna()
        print(f"H1 count at P < 0.05, model {m:12s}: {(v < 0.05).sum()} of {len(v)} subclasses; Bonferroni (x{len(v)}) P < 0.05: {(v * len(v) < 0.05).sum()}; "
              f"hits: {', '.join(v[v < 0.05].index)}")
print()
print("== pos/neg nuclei at FDR 10% (Sst and Microglia)")
print(S[S.label.isin(["sst", "microglia"])][["subclass", "model", "global_p", "k", "r2", "pos10", "neg10", "n_cells_fdr10", "covariates"]].to_string(index=False, float_format=lambda x: f"{x:.4g}"))

# maps: agreement of adjusted with base
agr = []; sup = []
for lab in ["sst", "microglia", "astrocyte", "oligodendrocyte", "opc", "pvalb", "vip", "lamp5", "l23it", "l4it", "l5it"]:
    fb = f"{T}/{lab}__all__base_cells.parquet"
    base = None
    if os.path.exists(fb):
        base = pd.read_parquet(fb)
    else:
        base = pd.read_parquet(f"{R}/cna/{lab}__all_cells.parquet", columns=["obs_name", "donor_id", "Supertype", "ncorr_CPS_Local", "fdr_CPS_Local"]).rename(
            columns={"ncorr_CPS_Local": "ncorr", "fdr_CPS_Local": "fdr"})
    for m in ["sat", "tech3", "neun"]:
        f = f"{T}/{lab}__all__{m}_cells.parquet"
        if not os.path.exists(f):
            continue
        a = pd.read_parquet(f)
        x = base.merge(a, on=["obs_name"], suffixes=("_b", "_a"))
        ok = x["ncorr_b"].notna() & x["ncorr_a"].notna()
        both10 = ok & (x["fdr_b"] <= 0.10) & (x["fdr_a"] <= 0.10)
        sig_b = ok & (x["fdr_b"] <= 0.10)
        sig_a = x["ncorr_a"].notna() & (x["fdr_a"] <= 0.10)
        agr.append(dict(subclass=LABELS[lab], model=m, n_cells=int(ok.sum()), spearman_ncorr=stats.spearmanr(x.loc[ok, "ncorr_b"], x.loc[ok, "ncorr_a"])[0],
                        pearson_ncorr=np.corrcoef(x.loc[ok, "ncorr_b"], x.loc[ok, "ncorr_a"])[0, 1],
                        base_fdr10_cells=int(sig_b.sum()), adj_fdr10_cells=int(sig_a.sum()), both_fdr10=int(both10.sum()),
                        sign_agree_among_both=float((np.sign(x.loc[both10, "ncorr_b"]) == np.sign(x.loc[both10, "ncorr_a"])).mean()) if both10.any() else np.nan))
        # supertype means
        gb = x[ok].groupby("Supertype_b")
        sm = pd.DataFrame({"mean_ncorr_base": gb["ncorr_b"].mean(), "mean_ncorr_adj": gb["ncorr_a"].mean(), "n": gb.size(),
                           "pos10_adj": gb.apply(lambda g: int(((g["ncorr_a"] > 0) & (g["fdr_a"] <= 0.10)).sum()), include_groups=False),
                           "neg10_adj": gb.apply(lambda g: int(((g["ncorr_a"] < 0) & (g["fdr_a"] <= 0.10)).sum()), include_groups=False),
                           "pos10_base": gb.apply(lambda g: int(((g["ncorr_b"] > 0) & (g["fdr_b"] <= 0.10)).sum()), include_groups=False),
                           "neg10_base": gb.apply(lambda g: int(((g["ncorr_b"] < 0) & (g["fdr_b"] <= 0.10)).sum()), include_groups=False)})
        sm.insert(0, "model", m); sm.insert(0, "subclass", LABELS[lab]); sm.index.name = "supertype"
        sup.append(sm.reset_index())
A = pd.DataFrame(agr)
A.to_csv(f"{R}/final_checks/cna_tech_map_agreement.csv", index=False)
if sup:
    pd.concat(sup).to_csv(f"{R}/final_checks/cna_tech_supertypes.csv", index=False)
print()
print("== map agreement of the covariate-adjusted ncorr with the Day-1/base ncorr (per nucleus)")
print(A.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
if sup:
    SU = pd.concat(sup)
    print()
    print("== Sst supertypes (>= 200 nuclei): sign of mean ncorr, base vs adjusted")
    for m in ["sat", "tech3"]:
        g = SU[(SU.subclass == "Sst") & (SU.model == m) & (SU.n >= 200)]
        if len(g):
            print(f"  model {m}: {len(g)} supertypes; sign agreement base/adjusted {(np.sign(g.mean_ncorr_base) == np.sign(g.mean_ncorr_adj)).mean():.2f}; "
                  f"negative mean ncorr base {(g.mean_ncorr_base < 0).sum()} -> adjusted {(g.mean_ncorr_adj < 0).sum()}; "
                  f"mean ncorr base {g.mean_ncorr_base.mean():+.3f} -> adjusted {g.mean_ncorr_adj.mean():+.3f}")
    g = SU[(SU.supertype == "Micro-PVM_3-SEAAD")]
    print()
    print("== Micro-PVM_3-SEAAD (disease-associated microglia supertype)")
    print(g[["model", "n", "mean_ncorr_base", "mean_ncorr_adj", "pos10_base", "neg10_base", "pos10_adj", "neg10_adj"]].to_string(index=False, float_format=lambda x: f"{x:.3f}"))
