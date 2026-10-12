#!/usr/bin/env python
"""Broken Bricks Day 3 figures (local; venv-fig python: pandas, numpy, matplotlib).

fig4_genetic_fingerprints_{light,dark}.png     MAGMA: A) cell-type specificity gene-property test (positive control), B) red/blue program tests conditional on specificity
fig5_apoe4_oligo_cholesterol_{light,dark}.png  H4: module score by APOE4 status, per-gene effects
Inputs: results/genetics/h3_celltype_property.csv, h3_program_sets.csv, results/genetics/h4/*.csv
"""
import argparse, json, os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ap = argparse.ArgumentParser()
ap.add_argument("--res", default=os.path.join(HERE, "..", "results", "genetics"))
ap.add_argument("--out", default=os.path.join(HERE, "..", "figures", "genetics"))
args = ap.parse_args()
os.makedirs(args.out, exist_ok=True)
THEMES = {
    "light": dict(bg="#fcfcfb", fg="#0b0b0b", fg2="#52514e", grid="#e6e5e1", blue="#1d5db5", red="#c5302c", grey="#b9b8b2", c1="#2a78d6", c2="#eb6834", mark="#0b0b0b"),
    "dark": dict(bg="#1a1a19", fg="#ffffff", fg2="#c3c2b7", grid="#2e2e2c", blue="#4b98f2", red="#f06c6a", grey="#5c5b57", c1="#3987e5", c2="#d95926", mark="#ffffff"),
}
def style(theme):
    t = THEMES[theme]
    plt.rcParams.update({"figure.facecolor": t["bg"], "axes.facecolor": t["bg"], "savefig.facecolor": t["bg"], "text.color": t["fg"], "axes.labelcolor": t["fg2"],
                         "xtick.color": t["fg2"], "ytick.color": t["fg2"], "axes.edgecolor": t["grid"], "font.family": "DejaVu Sans", "font.size": 10,
                         "axes.titlesize": 11, "axes.titleweight": "regular"})
    return t

ORDER = ["Microglia", "Astrocyte", "Oligodendrocyte", "OPC", "Sst", "Pvalb", "Vip", "Lamp5", "L23IT", "L4IT", "L5IT"]
DISP = {"L23IT": "L2/3 IT", "L4IT": "L4 IT", "L5IT": "L5 IT", "Microglia": "Microglia-PVM"}
prop = pd.read_csv(os.path.join(args.res, "h3_celltype_property.csv"))
sets = pd.read_csv(os.path.join(args.res, "h3_program_sets.csv"))
summ = json.load(open(os.path.join(args.res, "h3_summary.json")))
h4 = pd.read_csv(os.path.join(args.res, "h4", "h4_ols.csv"))
h4g = pd.read_csv(os.path.join(args.res, "h4", "h4_per_gene_ols.csv"))
h4d = pd.read_csv(os.path.join(args.res, "h4", "h4_donor_scores.csv"), index_col=0)
h4ds = json.load(open(os.path.join(args.res, "h4", "h4_pydeseq2_summary.json")))

def fig4(theme):
    t = style(theme)
    fig = plt.figure(figsize=(13.5, 6.6))
    gs = fig.add_gridspec(1, 2, width_ratios=[1, 1.15], wspace=0.42, left=0.11, right=0.97, top=0.82, bottom=0.2)
    a1 = fig.add_subplot(gs[0]); a2 = fig.add_subplot(gs[1])
    # A: gene-property
    p = prop[prop.tag == "primary"].set_index("VARIABLE")
    types = ORDER + ["Endothelial", "VLMC"]
    y = np.arange(len(types))[::-1]
    vals = [(-np.log10(max(p.loc[k, "P"], 1e-300))) if p.loc[k, "BETA"] > 0 else 0.0 for k in types]
    cols = [t["c1"] if k == "Microglia" else (t["fg2"] if p.loc[k, "BETA"] > 0 else t["grey"]) for k in types]
    a1.barh(y, vals, color=cols, height=0.62)
    a1.axvline(4, color=t["red"], lw=1.2, ls=(0, (4, 3))); a1.text(4.15, y[5], "P = 1e-4\n(preregistered\nthreshold)", color=t["red"], fontsize=8.5, va="center")
    a1.set_yticks(y); a1.set_yticklabels([DISP.get(k, k) for k in types])
    for yy, k, v in zip(y, types, vals):
        if v > 1.3: a1.text(v + 0.15, yy, f"P = {p.loc[k, 'P']:.1e}", va="center", fontsize=8.5, color=t["fg"])
    a1.set_xlabel("-log10 P (one-sided, gene-property test of cell-type specificity)")
    a1.set_title("A. Which cell types carry the AD GWAS signal?", loc="left", fontsize=11)
    for s in ["top", "right"]: a1.spines[s].set_visible(False)
    a1.set_xlim(0, max(vals) * 1.35)
    # B: programs conditional on specificity
    q = sets[(sets.tag == "primary") & (sets.kind == "primary") & (sets.analysis == "conditional")].set_index("VARIABLE")
    y2 = np.arange(len(ORDER))[::-1]
    a2.axvline(-np.log10(0.05), color=t["grid"], lw=1, ls=(0, (4, 3))); a2.text(-np.log10(0.05) + 0.03, y2[-1] - 0.45, "P = 0.05", color=t["fg2"], fontsize=8, va="bottom", rotation=90)
    a2.axvline(-np.log10(0.05 / 22), color=t["red"], lw=1, ls=(0, (4, 3))); a2.text(-np.log10(0.05 / 22) + 0.03, y2[-1] - 0.45, "Bonferroni, 22 tests", color=t["red"], fontsize=8, va="bottom", rotation=90)
    for yy, k in zip(y2, ORDER):
        for color, mk, c, off in [("red", "o", t["red"], 0.14), ("blue", "o", t["blue"], -0.14)]:
            nm = f"{k}_{color}"
            if nm in q.index:
                v = -np.log10(q.loc[nm, "P"])
                a2.scatter([v], [yy + off], s=60, color=c, edgecolors=t["bg"], linewidths=1, zorder=3)
    a2.set_yticks(y2); a2.set_yticklabels([DISP.get(k, k) for k in ORDER])
    a2.axhline(y2[3] - 0.5, color=t["grid"], lw=1)
    a2.text(0.995, y2[3] + 0.42, "glia", ha="right", va="bottom", color=t["fg2"], fontsize=8.5, transform=a2.get_yaxis_transform())
    a2.text(0.995, y2[3] - 0.42, "neurons", ha="right", va="top", color=t["fg2"], fontsize=8.5, transform=a2.get_yaxis_transform())
    a2.set_xlabel("-log10 P (set test conditional on own-cell-type specificity)")
    a2.set_title("B. Is the pathology-red or -blue program enriched for AD risk genes?", loc="left", fontsize=11)
    a2.scatter([], [], s=60, color=t["red"], label="red program (neighborhoods expanding with pathology)"); a2.scatter([], [], s=60, color=t["blue"], label="blue program (shrinking)")
    a2.legend(frameon=False, fontsize=8.5, loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=2)
    for s in ["top", "right"]: a2.spines[s].set_visible(False)
    mc = summ["primary_microglia_red_conditional"]
    fig.text(0.04, 0.955, "Damage and genes: microglia carry the AD genetic signal; the pathology-linked programs add little beyond cell identity", fontsize=13, ha="left")
    fig.text(0.04, 0.925, f"MAGMA on Bellenguez et al. 2022 (GCST90027158), 1000G EUR, window 10 kb. Microglia red program, conditional on microglial specificity: beta = {mc['beta']:.3f}, one-sided P = {mc['p']:.2f} (not enriched beyond specificity).\n"
             f"Neuronal programs: none with P < 0.05, so the prediction of no genetic fingerprints on the neuronal damage programs holds.", fontsize=8.5, ha="left", va="top", color=t["fg2"])
    fig.text(0.04, 0.02, "Astrocyte-red and OPC-red programs reach nominal P = 0.005 and 0.004 but do not survive correction over 22 tests (exploratory). Sensitivities (effective N, APOE-region exclusion) in results/genetics/h3_*.csv.",
             fontsize=8, ha="left", color=t["fg2"])
    fig.savefig(os.path.join(args.out, f"fig4_genetic_fingerprints_{theme}.png"), dpi=170)
    plt.close(fig)

def fig5(theme):
    t = style(theme)
    prim = h4.iloc[0]
    fig = plt.figure(figsize=(13.5, 6.4))
    gs = fig.add_gridspec(1, 2, width_ratios=[0.8, 1.25], wspace=0.35, left=0.07, right=0.97, top=0.8, bottom=0.17)
    a1 = fig.add_subplot(gs[0]); a2 = fig.add_subplot(gs[1])
    rng = np.random.default_rng(20261011)
    for g, (lab, c) in enumerate([("non-carriers", t["fg2"]), ("APOE4 carriers", t["c2"])]):
        v = h4d.loc[h4d.apoe4 == g, "module_score_adjusted_for_covariates"].values
        x = g + rng.uniform(-0.17, 0.17, len(v))
        a1.scatter(x, v, s=24, color=c, alpha=0.8, edgecolors=t["bg"], linewidths=0.5, zorder=3)
        m = v.mean(); se = v.std(ddof=1) / np.sqrt(len(v))
        a1.plot([g - 0.3, g + 0.3], [m, m], color=t["mark"], lw=2, zorder=4)
        a1.plot([g, g], [m - 1.96 * se, m + 1.96 * se], color=t["mark"], lw=1.2, zorder=4)
        a1.text(g, -2.3, f"n = {len(v)}", ha="center", va="bottom", fontsize=9, color=t["fg2"])
    a1.set_xticks([0, 1]); a1.set_xticklabels(["non-carriers", "APOE4 carriers"]); a1.set_xlim(-0.6, 1.6)
    a1.set_ylabel("module score, adjusted for CPS, age, sex, PMI", fontsize=9.5); a1.set_ylim(-2.4, 2.4)
    a1.axhline(0, color=t["grid"], lw=1)
    for s in ["top", "right"]: a1.spines[s].set_visible(False)
    a1.set_title("A. Donor-level module score", loc="left", fontsize=11)
    g = h4g.sort_values("coef_apoe4_z_units")
    yy = np.arange(len(g))
    a2.axvline(0, color=t["fg2"], lw=1)
    a2.hlines(yy, g.ci95_low, g.ci95_high, color=t["grey"], lw=1.5)
    a2.scatter(g.coef_apoe4_z_units, yy, s=34, color=[t["c2"] if v > 0 else t["c1"] for v in g.coef_apoe4_z_units], zorder=3)
    a2.set_yticks(yy); a2.set_yticklabels(g.gene, fontsize=8.5)
    a2.set_xlabel("APOE4 effect on gene (z-scored log2 CPM units), 95% CI")
    a2.set_title("B. The 20 module genes, one at a time", loc="left", fontsize=11)
    for s in ["top", "right"]: a2.spines[s].set_visible(False)
    fig.text(0.04, 0.955, "APOE4 and oligodendrocyte cholesterol biosynthesis: not replicated in 84 open donors", fontsize=13, ha="left")
    fig.text(0.04, 0.925, f"Module score: APOE4 effect {prim.coef_apoe4:+.2f} z units (95% CI {prim.ci95_low:+.2f} to {prim.ci95_high:+.2f}), one-sided P = {prim.p_one_sided_greater:.2f} (25 carriers vs 59 non-carriers, adjusted for CPS, age, sex, PMI).\n"
             f"pydeseq2 on the 20 genes: {h4ds['genes_up']} up, 0 at one-sided P < 0.05; summed module log2FC {h4ds['module_sum_log2fc']:+.2f} (95% CI {h4ds['module_sum_ci'][0]:+.2f} to {h4ds['module_sum_ci'][1]:+.2f}).",
             fontsize=8.5, ha="left", va="top", color=t["fg2"])
    fig.text(0.04, 0.03, "Replication attempt of the APOE4-cholesterol-myelination finding (Blanchard et al. 2022, different cohort); direction positive, effect small and not distinguishable from zero here. Mean and 95% CI bars shown.",
             fontsize=8, ha="left", color=t["fg2"])
    fig.savefig(os.path.join(args.out, f"fig5_apoe4_oligo_cholesterol_{theme}.png"), dpi=170)
    plt.close(fig)

for th in ["light", "dark"]:
    fig4(th); fig5(th)
print("done")
