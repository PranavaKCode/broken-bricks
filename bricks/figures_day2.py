#!/usr/bin/env python
"""Broken Bricks Day 2 figures (local; venv-fig python: pandas, numpy, matplotlib).

fig2_brick_matrix_{light,dark}.png      donors x subclasses, ordered by CPS_Local: brick score B (left) and residual Z = B - f(CPS) (right)
fig3_discordant_pairs_{light,dark}.png  zG vs zN scatter (discordant donors marked) and two donors matched on CPS_Local with opposite profiles
Pair rule (amendment A7 item 9, fixed before running): among pairs of donors with opposite discordance sign, the smallest |CPS difference|,
ties by the larger sum of |zG - zN|.  If no donor is discordant in both directions the rule falls back to the two donors with the largest
opposite |delta| = |zG - zN| among pairs with |CPS difference| < 0.05; this is stated in the figure footer.
usage: python figures_day2.py --res results/bricks/day2/all --tier allcols --out figures/bricks
"""
import argparse, json, os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize

HERE = os.path.dirname(os.path.abspath(__file__))
ap = argparse.ArgumentParser()
ap.add_argument("--res", default=os.path.join(HERE, "..", "results", "bricks", "all"))
ap.add_argument("--out", default=os.path.join(HERE, "..", "figures", "bricks"))
ap.add_argument("--tier", default="allcols", help="tier whose composites / discordance labels are drawn (prereg | allcols | techNAM)")
ap.add_argument("--zsrc", default="Z", help="column prefix for the residual panel: Z (primary NAM) or Ztech")
args = ap.parse_args()
os.makedirs(args.out, exist_ok=True)

LAB = {"microglia": "Microglia", "astrocyte": "Astrocyte", "oligodendrocyte": "Oligodendrocyte", "opc": "OPC", "sst": "Sst", "pvalb": "Pvalb", "vip": "Vip",
       "lamp5": "Lamp5", "l23it": "L2/3 IT", "l4it": "L4 IT", "l5it": "L5 IT"}
GLIA = ["microglia", "astrocyte", "oligodendrocyte", "opc"]; NEUR = ["sst", "pvalb", "vip", "lamp5", "l23it", "l4it", "l5it"]
ORDER = GLIA + NEUR
THEMES = {
    "light": dict(bg="#fcfcfb", fg="#0b0b0b", fg2="#52514e", grid="#e6e5e1", mid="#f0efec", blue="#1d5db5", red="#c5302c", seq=["#f3efe6", "#7a5195"], grey="#d9d8d3", amber="#b8860b", cat=["#eb6834", "#1baf7a"]),
    "dark": dict(bg="#1a1a19", fg="#ffffff", fg2="#c3c2b7", grid="#2e2e2c", mid="#383835", blue="#4b98f2", red="#f06c6a", seq=["#2b2733", "#b79ad6"], grey="#3f3f3c", amber="#e0b84a", cat=["#d95926", "#199e70"]),
}
def style(theme):
    t = THEMES[theme]
    plt.rcParams.update({"figure.facecolor": t["bg"], "axes.facecolor": t["bg"], "savefig.facecolor": t["bg"], "text.color": t["fg"], "axes.labelcolor": t["fg2"],
                         "xtick.color": t["fg2"], "ytick.color": t["fg2"], "axes.edgecolor": t["grid"], "font.family": "DejaVu Sans", "font.size": 10,
                         "axes.titlesize": 11, "axes.titleweight": "regular"})
    return t

d = pd.read_csv(os.path.join(args.res, "donor_bricks.csv"), index_col="donor_id")
rel = pd.read_csv(os.path.join(args.res, "brick_reliability.csv")).set_index("label")
summ = json.load(open(os.path.join(args.res, "summary.json")))
d = d[d["in_variant"] == 1]
retained = [l for l in ORDER if bool(rel.loc[l, "retained"])]
d = d.sort_values("CPS_Local", ascending=False)

def fig_matrix(theme):
    t = style(theme)
    n = len(d)
    fig = plt.figure(figsize=(14.5, 11.5))
    gs = fig.add_gridspec(1, 5, width_ratios=[0.5, 0.34, 5.6, 5.6, 0.12], wspace=0.06, left=0.045, right=0.95, top=0.83, bottom=0.12)
    aC = fig.add_subplot(gs[0]); aF = fig.add_subplot(gs[1]); aB = fig.add_subplot(gs[2]); aZ = fig.add_subplot(gs[3]); aK = fig.add_subplot(gs[4])
    y = np.arange(n)
    # CPS bar
    aC.barh(y, d["CPS_Local"].values, color=t["fg2"], height=0.8)
    aC.set_ylim(n - 0.5, -0.5); aC.invert_xaxis(); aC.set_xlim(0.85, 0)
    aC.set_yticks([]); aC.set_xticks([0, 0.8]); aC.tick_params(labelsize=8); aC.set_xlabel("CPS_Local", fontsize=9)
    for s in ["top", "right", "left"]: aC.spines[s].set_visible(False)
    # flag strip: severe, resilient, demented
    aF.set_xlim(-0.5, 2.5); aF.set_ylim(n - 0.5, -0.5)
    for j, (col, c, lab) in enumerate([("severe_donor", t["red"], "SAD"), ("resilient_high_adnc", t["blue"], "Res."), ("discordant__" + args.tier, t["amber"], "Disc.")]):
        m = d[col].values == 1
        aF.scatter(np.full(m.sum(), j), y[m], s=14, color=c, marker="s", linewidths=0)
        aF.text(j, -2.2, lab, rotation=90, ha="center", va="bottom", fontsize=8, color=t["fg2"])
    aF.set_xticks([]); aF.set_yticks([])
    for s in aF.spines.values(): s.set_visible(False)
    seq = LinearSegmentedColormap.from_list("seq", t["seq"], N=256)
    div = LinearSegmentedColormap.from_list("div", [t["blue"], t["mid"], t["red"]], N=256)
    Bm = d[[f"B_{l}" for l in ORDER]].values; Zm = d[[f"{args.zsrc}_{l}" for l in ORDER]].values
    im1 = aB.imshow(np.ma.masked_invalid(Bm), aspect="auto", cmap=seq, vmin=0, vmax=0.8, interpolation="nearest")
    lim = 2.5
    im2 = aZ.imshow(np.ma.masked_invalid(Zm), aspect="auto", cmap=div, norm=Normalize(-lim, lim), interpolation="nearest")
    for a_, ttl in [(aB, "Brick score B: how sick this donor's cells of each type look (out-of-fold, in CPS units)"),
                    (aZ, "Residual Z: sicker (red) or healthier (blue) than the donor's overall CPS predicts")]:
        a_.set_xticks(range(len(ORDER)))
        a_.set_xticklabels([LAB[l] + (" *" if l in retained else "") for l in ORDER], rotation=50, ha="left", fontsize=9)
        a_.xaxis.tick_top(); a_.set_yticks([]); a_.tick_params(length=0)
        a_.axvline(len(GLIA) - 0.5, color=t["bg"], lw=3)
        a_.set_title(ttl, fontsize=9, loc="left", color=t["fg2"], pad=70)
        for s in a_.spines.values(): s.set_visible(False)
    # colorbars
    cax = fig.add_axes([0.07, 0.065, 0.2, 0.012]); cb = fig.colorbar(im1, cax=cax, orientation="horizontal"); cb.set_label("B (predicted CPS)", fontsize=8); cb.outline.set_visible(False)
    cax2 = fig.add_axes([0.52, 0.065, 0.2, 0.012]); cb2 = fig.colorbar(im2, cax=cax2, orientation="horizontal"); cb2.set_label("residual z (blue healthier, red sicker than CPS predicts)", fontsize=8); cb2.outline.set_visible(False)
    for c_ in (cb, cb2): c_.ax.tick_params(labelsize=8, length=0)
    aK.axis("off")
    fig.text(0.045, 0.985, "Broken Bricks: 84 donors by 11 cell types, ordered by overall pathology (CPS_Local, highest at top)", fontsize=13, ha="left", va="top", color=t["fg"])
    fig.text(0.045, 0.958, f"Each cell is one donor's score in one cell type. Blank = donor has < 50 nuclei of that type. * = column passes the preregistered filter "
             f"(split-half r >= 0.6, out-of-fold rho > 0.4, P < 0.01).\n"
             f"Residual panel: tier '{args.tier}', {'techNAM' if args.zsrc == 'Ztech' else 'primary NAM'}. The preregistered test of H2 did not pass: the residuals are mostly one shared axis, not separate bricks.", fontsize=8.5, ha="left", va="top", color=t["fg2"])
    fig.text(0.045, 0.015, "Squares: SAD = SEA-AD Severely Affected Donor; Res. = High-ADNC donor without dementia; Disc. = discordant (glia and neuron residual of opposite sign, both |z| > 1.5).\n"
             "Residual columns track sequencing saturation (see h2_technical_regression.csv): exploratory, not a claim.", fontsize=7.5, ha="left", va="bottom", color=t["fg2"])
    fig.savefig(os.path.join(args.out, f"fig2_brick_matrix_{theme}.png"), dpi=170)
    plt.close(fig)

def pick_pair(tier):
    zg = d[f"zG__{tier}"]; zn = d[f"zN__{tier}"]; fl = d[f"discordant__{tier}"] == 1
    a = d.index[fl & (zg > 0)]; b = d.index[fl & (zn > 0)]
    rule = "opposite discordant donors, smallest |CPS difference|"
    best = None
    if len(a) and len(b):
        for i in a:
            for j in b:
                dc = abs(d.loc[i, "CPS_Local"] - d.loc[j, "CPS_Local"]); sc = abs(zg[i] - zn[i]) + abs(zg[j] - zn[j])
                key = (round(dc, 6), -sc)
                if best is None or key < best[0]: best = (key, i, j)
    if best is None:
        rule = "fallback: no donor is discordant in both directions; largest opposite |zG - zN| among pairs with |CPS difference| < 0.05"
        dl = (zg - zn).dropna()
        for i in dl.index:
            for j in dl.index:
                if dl[i] > 0 > dl[j] and abs(d.loc[i, "CPS_Local"] - d.loc[j, "CPS_Local"]) < 0.05:
                    key = (-(dl[i] - dl[j]), 0)
                    if best is None or key < best[0]: best = (key, i, j)
    return (best[1], best[2], rule) if best else (None, None, "no pair available")

def fig_pairs(theme):
    t = style(theme)
    tier = args.tier
    zg = d[f"zG__{tier}"]; zn = d[f"zN__{tier}"]; fl = d[f"discordant__{tier}"] == 1
    i, j, rule = pick_pair(tier)
    fig = plt.figure(figsize=(13.5, 6.2))
    gs = fig.add_gridspec(1, 2, width_ratios=[1, 1.35], wspace=0.22, left=0.06, right=0.97, top=0.82, bottom=0.2)
    a1 = fig.add_subplot(gs[0]); a2 = fig.add_subplot(gs[1])
    ok = zg.notna() & zn.notna()
    a1.axhline(0, color=t["grid"], lw=1); a1.axvline(0, color=t["grid"], lw=1)
    for v in (-1.5, 1.5):
        a1.axhline(v, color=t["grid"], lw=0.8, ls=(0, (3, 3))); a1.axvline(v, color=t["grid"], lw=0.8, ls=(0, (3, 3)))
    a1.scatter(zg[ok & ~fl], zn[ok & ~fl], s=22, color=t["grey"], edgecolors=t["fg2"], linewidths=0.4, label="other donors", zorder=3)
    m1 = ok & fl & (zg > 0); m2 = ok & fl & (zn > 0)
    a1.scatter(zg[m1], zn[m1], s=46, color=t["red"], edgecolors=t["bg"], linewidths=1, label="glia sicker, neurons healthier", zorder=4)
    a1.scatter(zg[m2], zn[m2], s=46, color=t["blue"], edgecolors=t["bg"], linewidths=1, label="neurons sicker, glia healthier", zorder=4)
    for k, nm in [(i, "A"), (j, "B")]:
        if k is not None:
            a1.scatter([zg[k]], [zn[k]], s=190, facecolors="none", edgecolors=t["fg"], linewidths=1.6, zorder=5)
            a1.annotate(nm, (zg[k], zn[k]), xytext=(11, 6), textcoords="offset points", fontsize=11, color=t["fg"])
    a1.set_xlabel("glia residual (z)  - sicker than CPS predicts →"); a1.set_ylabel("neuron residual (z)  - sicker than CPS predicts →")
    a1.legend(frameon=False, fontsize=8, loc="upper left"); a1.set_title(f"{int(fl.sum())} of {int(ok.sum())} donors discordant (tier '{tier}')", loc="left", fontsize=10.5)
    for s in ["top", "right"]: a1.spines[s].set_visible(False)
    if i is not None:
        x = np.arange(len(ORDER)); w = 0.38
        za = d.loc[i, [f"{args.zsrc}_{l}" for l in ORDER]].values.astype(float); zb = d.loc[j, [f"{args.zsrc}_{l}" for l in ORDER]].values.astype(float)
        a2.axhline(0, color=t["fg2"], lw=0.8)
        a2.bar(x - w / 2, za, w, color=t["cat"][0], label=f"A  {i}  (CPS {d.loc[i, 'CPS_Local']:.2f}, {d.loc[i, 'cognitive_status']})")
        a2.bar(x + w / 2, zb, w, color=t["cat"][1], label=f"B  {j}  (CPS {d.loc[j, 'CPS_Local']:.2f}, {d.loc[j, 'cognitive_status']})")
        a2.set_xticks(x); a2.set_xticklabels([LAB[l] for l in ORDER], rotation=40, ha="right", fontsize=9)
        a2.axvline(len(GLIA) - 0.5, color=t["grid"], lw=1)
        a2.text(1.5, a2.get_ylim()[1] * 0.92, "glia", ha="center", color=t["fg2"], fontsize=9); a2.text(7, a2.get_ylim()[1] * 0.92, "neurons", ha="center", color=t["fg2"], fontsize=9)
        a2.set_ylabel("residual z (sicker → / ← healthier than CPS predicts)")
        a2.legend(frameon=False, fontsize=8.5, loc="lower left")
        a2.set_title("Same CPS, largest glia-vs-neuron contrast (not discordance)", loc="left", fontsize=10)
        for s in ["top", "right"]: a2.spines[s].set_visible(False)
    fig.text(0.06, 0.95, "Same pathology score, different broken bricks?", fontsize=13, ha="left")
    fig.text(0.06, 0.915, "Descriptive only. Preregistered H2 was not supported: no donor meets the discordance rule (opposite signs, both |z| > 1.5). Donor IDs only; no individual is identified.", fontsize=8.5, color=t["fg2"], ha="left")
    fig.text(0.06, 0.012, f"Pair rule: {rule}.\n"
             f"Residual z-scores: {'techNAM' if args.zsrc == 'Ztech' else 'primary NAM'}. Residual columns track sequencing saturation (h2_technical_regression.csv).", fontsize=8, color=t["fg2"], ha="left", va="bottom")
    fig.savefig(os.path.join(args.out, f"fig3_discordant_pairs_{theme}.png"), dpi=170)
    plt.close(fig)
    return i, j, rule

for th in ["light", "dark"]:
    fig_matrix(th)
    i, j, rule = fig_pairs(th)
print("pair:", i, j, rule)
json.dump({"pair_A": i, "pair_B": j, "rule": rule, "tier": args.tier}, open(os.path.join(args.out, "fig3_pair_choice.json"), "w"))
