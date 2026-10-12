#!/usr/bin/env python
"""Broken Bricks Day 1 figures (local; python with pandas, pyarrow, matplotlib, scipy).

fig1_redblue_<label>_{light,dark}.png       CNA with CPS_Local: all neighborhoods | FDR<=10% only
fig1b_redblue_abeta_ptau_<label>_*.png      CPS_Local_ABeta | CPS_Local_pTau (FDR<=10% only)
fig1_summary_global_p_*.png                 global P per subclass and phenotype; share of red/blue cells
fig2_pc1_*.png                              PC1: Sst depletion, microglia supertypes, sign agreement with SEA-AD scCODA
fig3_milo_vs_cna_microglia_*.png            Milo cross-check
fig4_floor_cells_grid_*.png                 the cells floor on its grid, coloured by ncorr

Colour: diverging blue <-> red with a neutral grey midpoint (dataviz skill: two hues + neutral midpoint, equal arms); symmetric limits.
usage: python figures.py [--res results/cna] [--out figures/cna] [--only fig1,summary,pc1,milo,floor]
"""
import argparse, json, os, sys, warnings
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable

warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
ap = argparse.ArgumentParser()
ap.add_argument("--res", default=os.path.join(HERE, "..", "results", "cna"))
ap.add_argument("--out", default=os.path.join(HERE, "..", "figures", "cna"))
ap.add_argument("--only", default="fig1,summary,pc1,milo,floor")
ap.add_argument("--themes", default="light,dark")
args = ap.parse_args()
os.makedirs(args.out, exist_ok=True)

SUBCLASSES = {"microglia": "Microglia-PVM", "astrocyte": "Astrocyte", "oligodendrocyte": "Oligodendrocyte", "opc": "OPC", "sst": "Sst",
              "pvalb": "Pvalb", "vip": "Vip", "lamp5": "Lamp5", "l23it": "L2/3 IT", "l4it": "L4 IT", "l5it": "L5 IT"}
THEMES = {
    "light": dict(bg="#fcfcfb", fg="#0b0b0b", fg2="#52514e", grid="#e6e5e1", mid="#f0efec", blue="#1d5db5", red="#c5302c", grey="#d9d8d3",
                  cat=["#2a78d6", "#eb6834", "#1baf7a"]),
    "dark": dict(bg="#1a1a19", fg="#ffffff", fg2="#c3c2b7", grid="#2e2e2c", mid="#383835", blue="#4b98f2", red="#f06c6a", grey="#3f3f3c",
                 cat=["#3987e5", "#d95926", "#199e70"]),
}
PH = [("CPS_Local", "CPS_Local"), ("CPS_Local_ABeta", "CPS_Local_ABeta (amyloid)"), ("CPS_Local_pTau", "CPS_Local_pTau (tau)")]

def style(theme):
    t = THEMES[theme]
    plt.rcParams.update({"figure.facecolor": t["bg"], "axes.facecolor": t["bg"], "savefig.facecolor": t["bg"], "text.color": t["fg"],
                         "axes.labelcolor": t["fg2"], "xtick.color": t["fg2"], "ytick.color": t["fg2"], "axes.edgecolor": t["grid"],
                         "font.family": "DejaVu Sans", "font.size": 10, "axes.titlesize": 11, "axes.titleweight": "regular"})
    return t

def cmap_for(t):
    return LinearSegmentedColormap.from_list("bwr_neutral", [t["blue"], t["mid"], t["red"]], N=256)

def load_json(label, variant="all"):
    p = os.path.join(args.res, f"{label}__{variant}.json")
    return json.load(open(p)) if os.path.exists(p) else None

def cells(label, variant="all"):
    return pd.read_parquet(os.path.join(args.res, f"{label}__{variant}_cells.parquet"))

def clean_axes(ax):
    ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal", adjustable="datalim")
    for s in ax.spines.values(): s.set_visible(False)

def pfmt(p):
    return "P < 1e-4" if p <= 1.0001e-4 else f"P = {p:.4f}".rstrip("0") if p < 0.01 else f"P = {p:.3f}"

def umap_panel(ax, df, ph, t, cm, masked, size):
    nc = df[f"ncorr_{ph}"].values; fd = df[f"fdr_{ph}"].values
    lim = float(np.nanmax(np.abs(nc))) or 0.1
    sig = (fd <= 0.10) & ~np.isnan(nc)
    x, y = df["umap_1"].values, df["umap_2"].values
    if masked:
        bg = ~sig
        ax.scatter(x[bg], y[bg], s=size, c=t["grey"], linewidths=0, rasterized=True)
        o = np.argsort(np.abs(nc[sig]))
        ax.scatter(x[sig][o], y[sig][o], s=size, c=nc[sig][o], cmap=cm, vmin=-lim, vmax=lim, linewidths=0, rasterized=True)
    else:
        v = np.nan_to_num(nc, nan=0.0); o = np.argsort(np.abs(v))
        ax.scatter(x[o], y[o], s=size, c=v[o], cmap=cm, vmin=-lim, vmax=lim, linewidths=0, rasterized=True)
    clean_axes(ax)
    return lim, int(((nc > 0) & sig).sum()), int(((nc < 0) & sig).sum())

def colorbar(fig, ax, lim, cm, t, label):
    cb = fig.colorbar(ScalarMappable(norm=Normalize(-lim, lim), cmap=cm), ax=ax, fraction=0.04, pad=0.02, shrink=0.8)
    cb.outline.set_visible(False); cb.ax.tick_params(labelsize=8, colors=t["fg2"]); cb.set_label(label, fontsize=9, color=t["fg2"])
    cb.set_ticks([-lim, 0, lim]); cb.set_ticklabels([f"-{lim:.2f}", "0", f"+{lim:.2f}"])

def two_panels():
    fig = plt.figure(figsize=(13, 6.3), dpi=170)
    gs = fig.add_gridspec(1, 2, left=0.01, right=0.99, top=0.83, bottom=0.10, wspace=0.04)
    return fig, fig.add_subplot(gs[0]), fig.add_subplot(gs[1])

def cbar_below(fig, ax, lim, cm, t, label):
    fig.canvas.draw()
    p = ax.get_position()
    cax = fig.add_axes([p.x0 + 0.01, p.y0 - 0.055, 0.22 * p.width, 0.016])
    cb = fig.colorbar(ScalarMappable(norm=Normalize(-lim, lim), cmap=cm), cax=cax, orientation="horizontal")
    cb.outline.set_visible(False); cb.ax.tick_params(labelsize=8, colors=t["fg2"], length=2)
    cb.set_ticks([-lim, 0, lim]); cb.set_ticklabels([f"-{lim:.2f}", "0", f"+{lim:.2f}"])
    fig.text(p.x0 + 0.01 + 0.22 * p.width + 0.008, p.y0 - 0.047, label, fontsize=8.5, color=t["fg2"], va="center")

def fig_subclass(label, theme):
    t = style(theme); cm = cmap_for(t)
    j = load_json(label); df = cells(label)
    r = j["results"]["CPS_Local"]
    size = 0.6 if len(df) > 150000 else (1.0 if len(df) > 60000 else 2.0)
    fig, a1, a2 = two_panels()
    for ax, masked in [(a1, False), (a2, True)]:
        lim, npos, nneg = umap_panel(ax, df, "CPS_Local", t, cm, masked, size)
        ax.set_title("All neighborhoods" if not masked else f"FDR <= 10% only  (red {npos:,} | blue {nneg:,} cells)", loc="left", color=t["fg2"], fontsize=10.5)
        cbar_below(fig, ax, lim, cm, t, "ncorr: neighborhood correlation with CPS_Local")
    fig.text(0.015, 0.965, f"{j['display']}: neighborhoods that expand (red) or shrink (blue) with CPS_Local",
             fontsize=14, color=t["fg"], va="top")
    fig.text(0.015, 0.915, f"CNA, {j['n_donors']} aged donors, {j['n_cells_used']:,} nuclei, global {pfmt(r['global_p'])} (10,000 permutations), "
             f"age/sex/PMI/method adjusted. SEA-AD native UMAP.", fontsize=9.5, color=t["fg2"], va="top")
    fig.savefig(os.path.join(args.out, f"fig1_redblue_{label}_{theme}.png")); plt.close(fig)

def fig_subclass_bricks(label, theme):
    t = style(theme); cm = cmap_for(t)
    j = load_json(label); df = cells(label)
    size = 0.6 if len(df) > 150000 else (1.0 if len(df) > 60000 else 2.0)
    fig, a1, a2 = two_panels()
    for ax, (ph, nm) in [(a1, PH[1]), (a2, PH[2])]:
        r = j["results"][ph]
        lim, npos, nneg = umap_panel(ax, df, ph, t, cm, True, size)
        ax.set_title(f"{nm}: {pfmt(r['global_p'])}; FDR <= 10% red {npos:,} | blue {nneg:,}", loc="left", color=t["fg2"], fontsize=10)
        cbar_below(fig, ax, lim, cm, t, "ncorr")
    fig.text(0.015, 0.965, f"{j['display']}: the amyloid and tau colorings", fontsize=14, color=t["fg"], va="top")
    fig.text(0.015, 0.915, "CNA on the two components of the pathology score; grey = not significant at FDR 10%. Same covariates and cells as the CPS_Local map.",
             fontsize=9.5, color=t["fg2"], va="top")
    fig.savefig(os.path.join(args.out, f"fig1b_redblue_abeta_ptau_{label}_{theme}.png")); plt.close(fig)

def fig_summary(theme):
    t = style(theme)
    h1 = pd.read_csv(os.path.join(args.res, "h1_global_p.csv"))
    h1["order"] = h1["label"].map({k: i for i, k in enumerate(SUBCLASSES)})
    h1 = h1.sort_values("order")
    fig, axs = plt.subplots(1, 2, figsize=(12.5, 5.6), dpi=170, gridspec_kw={"width_ratios": [1.15, 1]})
    ax = axs[0]; y = np.arange(len(h1))[::-1]
    for (ph, nm), short, c, off in zip(PH, ["CPS", "ABeta", "pTau"], t["cat"], [0.22, 0, -0.22]):
        ax.scatter(-np.log10(h1[f"P_{short}"]), y + off, s=46, color=c, label=nm, zorder=3, edgecolors=t["bg"], linewidths=1.2)
    ax.axvline(-np.log10(0.05), color=t["fg2"], lw=1, ls=(0, (4, 3))); ax.axvline(4.0, color=t["fg2"], lw=1, ls=":")
    ax.text(-np.log10(0.05) + 0.05, len(h1) - 0.55, "P = 0.05", fontsize=8.5, color=t["fg2"])
    ax.text(4.0 - 0.05, len(h1) - 0.55, "permutation floor", fontsize=8.5, color=t["fg2"], ha="right")
    ax.set_yticks(y); ax.set_yticklabels(h1["subclass"], color=t["fg"]); ax.set_xlim(0, 4.3)
    ax.set_xlabel("-log10 global P  (CNA, 84 donors, 10,000 permutations)")
    ax.grid(axis="x", color=t["grid"], lw=0.8); ax.set_axisbelow(True)
    for s in ["top", "right", "left"]: ax.spines[s].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=3, frameon=False, fontsize=8.5, labelcolor=t["fg2"])
    ax.set_title("Global association with pathology", loc="left", color=t["fg2"])
    ax = axs[1]
    tot = h1["n_cells"].values
    pos = h1["pos_fdr10_CPS"].values / tot * 100; neg = h1["neg_fdr10_CPS"].values / tot * 100
    ax.barh(y, pos, color=t["red"], height=0.62, label="expand with CPS (red)")
    ax.barh(y, -neg, color=t["blue"], height=0.62, label="shrink with CPS (blue)")
    ax.axvline(0, color=t["fg2"], lw=0.8)
    for yi, p, n in zip(y, pos, neg):
        ax.text(p + 0.4, yi, f"{p:.1f}", va="center", fontsize=8, color=t["fg2"]); ax.text(-n - 0.4, yi, f"{n:.1f}", va="center", ha="right", fontsize=8, color=t["fg2"])
    m = max(pos.max(), neg.max()) * 1.25 + 1
    ax.set_xlim(-m, m); ax.set_yticks(y); ax.set_yticklabels([]); ax.tick_params(axis="y", length=0)
    ax.set_xlabel("% of nuclei in neighborhoods with FDR <= 10%"); ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{abs(v):.0f}"))
    ax.grid(axis="x", color=t["grid"], lw=0.8); ax.set_axisbelow(True)
    for s in ["top", "right", "left"]: ax.spines[s].set_visible(False)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=2, frameon=False, fontsize=8.5, labelcolor=t["fg2"])
    ax.set_title("Which way the significant neighborhoods go (CPS_Local)", loc="left", color=t["fg2"])
    nsig = int((h1["P_CPS"] < 0.05).sum())
    fig.text(0.02, 0.965, f"H1: {nsig} of {len(h1)} subclasses show a global association with CPS_Local (P < 0.05)", fontsize=14, va="top")
    fig.text(0.02, 0.917, "Covarying neighborhood analysis on SEA-AD MTG, 84 aged donors; covariates age, sex, PMI, sequencing method.", fontsize=9.5, color=t["fg2"], va="top")
    fig.subplots_adjust(left=0.11, right=0.98, top=0.80, bottom=0.17, wspace=0.07)
    fig.savefig(os.path.join(args.out, f"fig1_summary_global_p_{theme}.png")); plt.close(fig)

def fig_pc1(theme):
    t = style(theme); cm = cmap_for(t)
    conc = pd.read_csv(os.path.join(args.res, "pc1_sccoda_concordance.csv"))
    pcs = json.load(open(os.path.join(args.res, "pc1_summary.json")))["PC1"]
    fig, axs = plt.subplots(1, 3, figsize=(15, 5.4), dpi=170, gridspec_kw={"width_ratios": [1, 0.9, 1.15]})
    # (a) Sst ncorr histogram
    ax = axs[0]; s = cells("sst"); nc = s["ncorr_CPS_Local"].dropna().values; fd = s.loc[s["ncorr_CPS_Local"].notna(), "fdr_CPS_Local"].values
    lim = np.abs(nc).max(); bins = np.linspace(-lim, lim, 61)
    ax.hist(nc[(fd > 0.10)], bins=bins, color=t["grey"], label="FDR > 10%")
    ax.hist(nc[(fd <= 0.10) & (nc < 0)], bins=bins, color=t["blue"], label="FDR <= 10%, depleted")
    ax.hist(nc[(fd <= 0.10) & (nc > 0)], bins=bins, color=t["red"], label="FDR <= 10%, enriched")
    ax.axvline(0, color=t["fg2"], lw=0.8)
    c2 = pcs["crit2_sst_depleted"]
    ax.set_title(f"Sst: depleted {c2['fdr10_neg_cells']:,} vs enriched {c2['fdr10_pos_cells']:,} cells at FDR 10%", loc="left", color=t["fg2"], fontsize=10)
    ax.set_xlabel("ncorr with CPS_Local"); ax.set_ylabel("nuclei"); ax.legend(frameon=False, fontsize=8, labelcolor=t["fg2"], loc="upper left")
    for sp_ in ["top", "right"]: ax.spines[sp_].set_visible(False)
    # (b) microglia supertypes
    ax = axs[1]; m = conc[conc["subclass"] == "Microglia-PVM"].sort_values("mean_ncorr")
    cols = [t["red"] if v > 0 else t["blue"] for v in m["mean_ncorr"]]
    yy = np.arange(len(m)); ax.barh(yy, m["mean_ncorr"], color=cols, height=0.62)
    ax.set_yticks(yy); ax.set_yticklabels(m["supertype"], fontsize=8.5, color=t["fg"]); ax.axvline(0, color=t["fg2"], lw=0.8)
    for yi, v, n in zip(yy, m["mean_ncorr"], m["n_cells"]):
        ax.text(v + (0.004 if v >= 0 else -0.004), yi, f"n={int(n):,}", fontsize=7.5, color=t["fg2"], va="center", ha="left" if v >= 0 else "right")
    ax.set_xlim(m["mean_ncorr"].min() * 1.9 - 0.02, m["mean_ncorr"].max() * 1.9 + 0.02)
    ax.set_xlabel("mean ncorr over the supertype's nuclei"); ax.set_title("Microglia supertypes", loc="left", color=t["fg2"], fontsize=10)
    for sp_ in ["top", "right", "left"]: ax.spines[sp_].set_visible(False)
    ax.tick_params(axis="y", length=0)
    # (c) concordance scatter
    ax = axs[2]; c = conc[conc["credible_primary"]].copy()
    ok = c["agree_primary"] == 1
    ax.axhline(0, color=t["fg2"], lw=0.8); ax.axvline(0, color=t["fg2"], lw=0.8)
    ax.scatter(c.loc[ok, "mean_ncorr"], c.loc[ok, "sea_ad_local_model"], s=34, color=t["cat"][2], label=f"sign agrees ({int(ok.sum())})", edgecolors=t["bg"], linewidths=0.8, zorder=3)
    ax.scatter(c.loc[~ok, "mean_ncorr"], c.loc[~ok, "sea_ad_local_model"], s=46, marker="D", color=t["cat"][1], label=f"sign differs ({int((~ok).sum())})", edgecolors=t["bg"], linewidths=0.8, zorder=4)
    offs = [(6, 5), (-34, 7), (7, -11), (6, 5)]
    for (_, r), off in zip(c.loc[~ok].iterrows(), offs):
        ax.annotate(r["supertype"], (r["mean_ncorr"], r["sea_ad_local_model"]), fontsize=7.5, color=t["fg2"], xytext=off, textcoords="offset points")
    k = pcs["crit3_sign_agreement_primary_summary_table"]
    ax.set_xlabel("CNA: supertype-mean ncorr (CPS_Local)"); ax.set_ylabel("SEA-AD scCODA effect (summary table, credible only)")
    ax.set_title(f"Sign agreement {k['n_agree']}/{k['n_credible_matched']} = {100*k['agreement']:.0f}%  (criterion: >= 70%)", loc="left", color=t["fg2"], fontsize=10)
    ax.legend(frameon=False, fontsize=8, labelcolor=t["fg2"], loc="lower right")
    for sp_ in ["top", "right"]: ax.spines[sp_].set_visible(False)
    verdict = "PASS" if pcs["PC1_pass_primary"] else "FAIL"
    fig.text(0.015, 0.965, f"PC1 positive controls: {verdict}", fontsize=14, va="top")
    c1 = pcs["crit1_global_P"]
    fig.text(0.015, 0.915, f"Microglia {pfmt(c1['microglia_P'])}, Sst {pfmt(c1['sst_P'])} (criterion P < 0.05); Sst depleted: {'yes' if c2['pass'] else 'no'}; "
             f"sign agreement with SEA-AD scCODA: {100*k['agreement']:.0f}% of {k['n_credible_matched']} credible supertypes.", fontsize=9.5, color=t["fg2"], va="top")
    fig.subplots_adjust(left=0.05, right=0.985, top=0.80, bottom=0.12, wspace=0.62)
    fig.savefig(os.path.join(args.out, f"fig2_pc1_{theme}.png")); plt.close(fig)

def fig_milo(theme):
    t = style(theme); cm = cmap_for(t)
    nh = pd.read_csv(os.path.join(args.res, "milo_microglia_nhoods.csv")); S = json.load(open(os.path.join(args.res, "milo_microglia_summary.json")))
    mc = pd.read_parquet(os.path.join(args.res, "milo_microglia_cells.parquet")).set_index("obs_name")
    df = cells("microglia").set_index("obs_name").join(mc[["milo_logFC_mean"]])
    fig, axs = plt.subplots(1, 3, figsize=(15.5, 5.3), dpi=170, gridspec_kw={"width_ratios": [1, 1, 1]})
    ax = axs[0]; ok = nh["logFC"].notna()
    sg = (nh["SpatialFDR"] < 0.10) & ok
    ax.scatter(nh.loc[ok & ~sg, "cna_mean_ncorr"], nh.loc[ok & ~sg, "logFC"], s=7, color=t["grey"], linewidths=0, label="Milo SpatialFDR >= 10%")
    ax.scatter(nh.loc[sg, "cna_mean_ncorr"], nh.loc[sg, "logFC"], s=9, color=t["cat"][2], linewidths=0, label="Milo SpatialFDR < 10%")
    ax.axhline(0, color=t["fg2"], lw=0.8); ax.axvline(0, color=t["fg2"], lw=0.8)
    ax.set_xlabel("CNA: mean ncorr of the neighborhood's nuclei"); ax.set_ylabel("Milo logFC per SD of CPS_Local")
    ax.set_title(f"{len(nh):,} neighborhoods: Spearman {S['spearman_nhood_logFC_vs_mean_ncorr'][0]:.2f}", loc="left", color=t["fg2"], fontsize=10)
    ax.legend(frameon=False, fontsize=8, labelcolor=t["fg2"], loc="upper left")
    for sp_ in ["top", "right"]: ax.spines[sp_].set_visible(False)
    for ax, col, ttl in [(axs[1], "ncorr_CPS_Local", "CNA ncorr"), (axs[2], "milo_logFC_mean", "Milo logFC (mean over neighborhoods)")]:
        v = df[col].values; lim = np.nanpercentile(np.abs(v), 99.5)
        o = np.argsort(np.abs(np.nan_to_num(v)))
        ax.scatter(df["umap_1"].values[o], df["umap_2"].values[o], s=1.0, c=np.nan_to_num(v)[o], cmap=cm, vmin=-lim, vmax=lim, linewidths=0, rasterized=True)
        clean_axes(ax); ax.set_title(ttl, loc="left", color=t["fg2"], fontsize=10); colorbar(fig, ax, lim, cm, t, "")
    fig.text(0.015, 0.965, "Milo cross-check on Microglia-PVM", fontsize=14, va="top")
    fig.text(0.015, 0.915, f"Cell-level Spearman between the two maps {S['spearman_cell_milo_vs_ncorr'][0]:.2f}; Milo (pydeseq2 solver) finds {S['n_nhoods_spatialFDR_lt_0.10']} "
             f"neighborhoods at SpatialFDR < 10% ({S['n_pos']} up, {S['n_neg']} down).", fontsize=9.5, color=t["fg2"], va="top")
    fig.subplots_adjust(left=0.05, right=0.95, top=0.80, bottom=0.12, wspace=0.22)
    fig.savefig(os.path.join(args.out, f"fig3_milo_vs_cna_microglia_{theme}.png")); plt.close(fig)

def fig_floor(theme):
    import anndata as ad
    t = style(theme); cm = cmap_for(t)
    a = ad.read_h5ad(os.path.join(args.res, "floor1_cells.h5ad"))
    U = a.obsm["X_umap"]; o = a.obs
    grid = {"Microglia-PVM": (0, 0), "Astrocyte": (1, 0), "Oligodendrocyte": (2, 0), "OPC": (3, 0), "Sst": (0, 1), "Pvalb": (1, 1), "L2/3 IT": (2, 1)}
    fig, axs = plt.subplots(2, 1, figsize=(11.5, 8.6), dpi=170)
    for ax, col in zip(axs, ["ncorr_CPS", "redblue_fdr10_CPS"]):
        if col == "ncorr_CPS":
            v = o[col].fillna(0).values; lim = 0.5; od = np.argsort(np.abs(v))
            ax.scatter(U[od, 0], U[od, 1], s=0.35, c=v[od], cmap=cm, vmin=-lim, vmax=lim, linewidths=0, rasterized=True)
        else:
            v = o[col].values
            ax.scatter(U[v == 0, 0], U[v == 0, 1], s=0.35, c=t["grey"], linewidths=0, rasterized=True)
            ax.scatter(U[v == 1, 0], U[v == 1, 1], s=0.35, c=t["red"], linewidths=0, rasterized=True)
            ax.scatter(U[v == -1, 0], U[v == -1, 1], s=0.35, c=t["blue"], linewidths=0, rasterized=True)
        for name, (gx, gy) in grid.items():
            n = int((o["subclass"] == name).sum())
            ax.text(gx * 26 + 0.5, -gy * 26 + 21.2, f"{name} ({n:,})", fontsize=8.5, color=t["fg2"], va="bottom")
        clean_axes(ax)
        ax.set_xlim(U[:, 0].min() - 2, U[:, 0].max() + 2); ax.set_ylim(U[:, 1].min() - 2, U[:, 1].max() + 6)
        ax.set_aspect("equal", adjustable="box", anchor="W")
    axs[0].set_title("ncorr with CPS_Local, all nuclei", loc="left", color=t["fg2"], fontsize=10.5)
    axs[1].set_title("FDR <= 10% flag: red = expands with CPS, blue = shrinks, grey = not significant", loc="left", color=t["fg2"], fontsize=10.5)
    fig.text(0.015, 0.985, f"Mantis floor 1 (cells): {len(o):,} nuclei, {o['donor_id'].nunique()} donors, 7 subclasses", fontsize=14, va="top")
    fig.text(0.015, 0.952, "Each subclass is drawn on its own SEA-AD native UMAP, placed on a grid; the scVI space (X_pca) is per subclass.", fontsize=9.5, color=t["fg2"], va="top")
    fig.subplots_adjust(left=0.01, right=0.99, top=0.89, bottom=0.08, hspace=0.10)
    cbar_below(fig, axs[1], 0.5, cm, t, "ncorr (panel above); red / blue flags (this panel)")
    fig.savefig(os.path.join(args.out, f"fig4_floor_cells_grid_{theme}.png")); plt.close(fig)

only = set(args.only.split(","))
for theme in args.themes.split(","):
    if "fig1" in only:
        for label in SUBCLASSES:
            if load_json(label):
                fig_subclass(label, theme); fig_subclass_bricks(label, theme); print("fig1", label, theme, flush=True)
    if "summary" in only: fig_summary(theme); print("summary", theme, flush=True)
    if "pc1" in only: fig_pc1(theme); print("pc1", theme, flush=True)
    if "milo" in only: fig_milo(theme); print("milo", theme, flush=True)
    if "floor" in only: fig_floor(theme); print("floor", theme, flush=True)
