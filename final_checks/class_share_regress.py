"""Final check 3 (part 2): the class-level test that CNA cannot give. Per-donor share of each SEA-AD subclass / class among nuclei, against CPS_Local,
with and without technical covariates. Runs locally (numpy, pandas, scipy only).

Input : results/final_checks/donor_subclass_counts.csv   (ORCD job subclass_counts.py: nuclei per donor x Class x Subclass, 84 aged donors, 1,178,694 nuclei)
        results/bricks/donor_bricks.csv                    (CPS_Local, age, sex, PMI, technical covariates, severe-donor flag)
Output: results/final_checks/class_share_regression.csv    (one row per outcome x denominator x model)
        results/final_checks/class_share_sst_summary.txt

Outcome: arcsine(sqrt(share)), standardised over donors. Predictor of interest: z(CPS_Local). Models (covariates beyond CPS_Local):
  A   age, sex, PMI, multiome fraction                           (the Day-1 CNA covariates)
  B   A + NeuN-positive fraction on FANS                          (a sorting parameter; matters for shares of all nuclei)
  C   A + sequencing saturation
  D   A + NeuN fraction + sequencing saturation
  E   A + NeuN fraction + RIN + log10 nuclei + sequencing saturation
Denominators: 'all' (all nuclei), 'neuron' (Glutamatergic + GABAergic nuclei), 'class' (nuclei of the subclass's own class: the within-class share)
Variants: all 84 donors, and without the 11 Severely Affected Donors (73).
Also: centred log-ratio (CLR) of the 24-subclass composition (pseudocount 0.5), which does not depend on a chosen denominator.
"""
import os, sys, json
import numpy as np, pandas as pd
from scipy import stats

R = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")
cnt = pd.read_csv(f"{R}/final_checks/donor_subclass_counts.csv")
d = pd.read_csv(f"{R}/bricks/donor_bricks.csv").set_index("donor_id")
d["sex_male"] = (d["sex"] == "Male").astype(float)
cnt = cnt[cnt["donor_id"].isin(d.index)]
mat = cnt.pivot_table(index="donor_id", columns="Subclass", values="n", aggfunc="sum", fill_value=0).reindex(d.index).fillna(0)
cls_of = cnt.drop_duplicates("Subclass").set_index("Subclass")["Class"].to_dict()
classes = sorted(set(cls_of.values()))
cmat = pd.DataFrame({c: mat[[s for s in mat.columns if cls_of[s] == c]].sum(axis=1) for c in classes})
neur = cmat[[c for c in classes if c.startswith("Neuronal")]].sum(axis=1)
tot = mat.sum(axis=1)
assert tot.sum() == 1178694, tot.sum()

base = ["age", "sex_male", "pmi", "tech_multiome_frac"]
MODELS = {"A": base, "B": base + ["tech_neun_frac"], "C": base + ["tech_seq_saturation"],
          "D": base + ["tech_neun_frac", "tech_seq_saturation"],
          "E": base + ["tech_neun_frac", "tech_rin", "tech_log10_n_nuclei", "tech_seq_saturation"]}


def z(x):
    x = np.asarray(x, float)
    return (x - x.mean()) / x.std(ddof=1)


def ols_t(y, X):
    X1 = np.column_stack([np.ones(len(y)), X]); b, *_ = np.linalg.lstsq(X1, y, rcond=None)
    r = y - X1 @ b; df = len(y) - X1.shape[1]; s2 = (r ** 2).sum() / df
    se = np.sqrt(np.diag(s2 * np.linalg.inv(X1.T @ X1))); t = b / se
    return b, se, 2 * stats.t.sf(np.abs(t), df), df


def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p); q = np.empty(n); prev = 1.0
    for i, idx in enumerate(o[::-1]):
        prev = min(prev, p[idx] * n / (n - i)); q[idx] = prev
    return q


outcomes = {}   # (name, denominator) -> Series over donors
for s in mat.columns:
    outcomes[(s, "all")] = mat[s] / tot
    if cls_of[s].startswith("Neuronal"):
        outcomes[(s, "neuron")] = mat[s] / neur
    outcomes[(s, "class")] = mat[s] / cmat[cls_of[s]]
for c in classes:
    outcomes[(c, "all")] = cmat[c] / tot
    if c.startswith("Neuronal"):
        outcomes[(c, "neuron")] = cmat[c] / neur
clr = np.log(mat + 0.5); clr = clr.sub(clr.mean(axis=1), axis=0)
for s in mat.columns:
    outcomes[(s, "clr")] = clr[s]

rows = []
for variant, keep in [("all84", np.ones(len(d), bool)), ("nosad73", (d["severe_donor"] == 0).values)]:
    dd = d.loc[keep]
    for (name, den), ser in outcomes.items():
        v = ser.loc[dd.index]
        if den == "clr":
            y = z(v.values)
        else:
            if v.isna().any() or (v.values == 0).all():
                continue
            y = z(np.arcsin(np.sqrt(v.values)))
        sp = stats.spearmanr(dd["CPS_Local"], v.values)
        for m, covs in MODELS.items():
            X = np.column_stack([z(dd["CPS_Local"])] + [z(dd[c]) for c in covs])
            b, se, p, df = ols_t(y, X)
            rows.append(dict(variant=variant, outcome=name, denominator=den, model=m, n=len(dd), beta_cps=b[1], se=se[1], p=p[1],
                             spearman_raw=sp[0], spearman_p=sp[1], mean_share=float(np.mean(v.values)) if den != "clr" else np.nan))
res = pd.DataFrame(rows)
for (var, den, m), g in res.groupby(["variant", "denominator", "model"]):
    sub = g[~g["outcome"].str.startswith("Neuronal") & ~g["outcome"].str.startswith("Non-neuronal")] if den != "all" else g
    res.loc[sub.index, "q_bh_within_family"] = bh(sub["p"].values)
res.to_csv(f"{R}/final_checks/class_share_regression.csv", index=False)

out = []
def P(*a):
    s = " ".join(str(x) for x in a); print(s); out.append(s)

P("Donors", len(d), "nuclei", int(tot.sum()), "| SAD donors", int(d["severe_donor"].sum()))
P("Mean Sst share of all nuclei %.4f (SD %.4f); of neurons %.4f; of GABAergic nuclei %.4f" % (
    (mat["Sst"] / tot).mean(), (mat["Sst"] / tot).std(), (mat["Sst"] / neur).mean(), (mat["Sst"] / cmat["Neuronal: GABAergic"]).mean()))
sev = d["severe_donor"] == 1
for den, ser in [("all nuclei", mat["Sst"] / tot), ("neurons", mat["Sst"] / neur), ("GABAergic", mat["Sst"] / cmat["Neuronal: GABAergic"])]:
    P("  Sst share of %s: severe donors %.4f vs others %.4f (Mann-Whitney p = %.2g)" % (den, ser[sev].mean(), ser[~sev].mean(), stats.mannwhitneyu(ser[sev], ser[~sev]).pvalue))
P()
P("CPS_Local coefficient (per SD of CPS, outcome = SD units of arcsine-sqrt share). Models: A base | B +NeuN | C +seq saturation | D +NeuN +sat | E +NeuN +RIN +nuclei +sat")
for variant in ["all84", "nosad73"]:
    for den in ["all", "neuron", "class", "clr"]:
        for name in ["Sst", "Pvalb", "Vip", "Lamp5", "Immune", "Astrocyte", "Oligodendrocyte", "OPC", "L2/3 IT", "L4 IT", "L5 IT"]:
            g = res[(res.variant == variant) & (res.denominator == den) & (res.outcome == name)]
            if len(g) == 0:
                continue
            P(f"{variant:8s} {den:7s} {name:16s} n={g.n.iloc[0]}  rho_raw={g.spearman_raw.iloc[0]:+.2f} | " + " | ".join(
                f"{m} {r.beta_cps:+.2f} (p={r.p:.2g})" for m, r in g.set_index('model').iterrows()))
    P()
open(f"{R}/final_checks/class_share_sst_summary.txt", "w").write("\n".join(out))


# ---- cross-check with the CELLxGENE per-nucleus labels used by the CNA (results/cna/<label>__all_cells.parquet, 11 subclasses) -------------------------
# The CELLxGENE objects hold 1,071,242 aged nuclei in 13 subclass files; the final-nuclei metadata CSV holds fewer in the same subclasses (see
# cellxgene_vs_final_nuclei_counts.csv). Share of each subclass among the nuclei of the 11 analysed subclasses, CELLxGENE labels:
LAB = {"microglia": "Immune", "astrocyte": "Astrocyte", "oligodendrocyte": "Oligodendrocyte", "opc": "OPC", "sst": "Sst", "pvalb": "Pvalb", "vip": "Vip",
       "lamp5": "Lamp5", "l23it": "L2/3 IT", "l4it": "L4 IT", "l5it": "L5 IT"}
cx = {}
for k, v in LAB.items():
    c = pd.read_parquet(f"{R}/cna/{k}__all_cells.parquet", columns=["donor_id"])["donor_id"].value_counts()
    cx[v] = c.reindex(d.index).fillna(0)
cx = pd.DataFrame(cx); cxt = cx.sum(axis=1)
rows2 = []
for name in ["Sst", "Vip", "Pvalb", "Lamp5", "L2/3 IT", "L4 IT", "L5 IT", "Immune"]:
    y = z(np.arcsin(np.sqrt((cx[name] / cxt).values)))
    for m in ["A", "B", "C", "D", "E"]:
        X = np.column_stack([z(d["CPS_Local"])] + [z(d[c]) for c in MODELS[m]])
        b, se, p, df = ols_t(y, X)
        rows2.append(dict(source="CELLxGENE labels, share among the 11 analysed subclasses", outcome=name, model=m, n=len(d), beta_cps=b[1], se=se[1], p=p[1]))
cxr = pd.DataFrame(rows2); cxr.to_csv(f"{R}/final_checks/class_share_regression_cellxgene.csv", index=False)
P("")
P("Cross-check with the CELLxGENE labels (share among the 11 analysed subclasses), CPS_Local coefficient, models A to E:")
for name in ["Sst", "Vip", "L2/3 IT", "L4 IT"]:
    g = cxr[cxr.outcome == name].set_index("model")
    P(f"  {name:8s} " + " | ".join(f"{m} {r.beta_cps:+.2f} (p={r.p:.2g})" for m, r in g.iterrows()))
open(f"{R}/final_checks/class_share_sst_summary.txt", "w").write("\n".join(out))
