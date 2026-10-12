#!/usr/bin/env python
"""Stage E: figures from out/results/*.csv (and work/ parquet for the map). Light and dark variants, 16:9."""
import os, sys
from pathlib import Path
import numpy as np
import polars as pl
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

LIT = Path(os.path.expanduser("~/orcd/pool/mantis-films/research/lit"))
WORK, RES = LIT / "work", LIT / "out" / "results"
FIG = LIT / "out" / "figures"; FIG.mkdir(parents=True, exist_ok=True)

CAT = {"light": ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"],
       "dark": ["#3987e5", "#d95926", "#199e70", "#c98500", "#d55181", "#008300", "#9085e9", "#e66767"]}
THEME = {"light": dict(surface="#fcfcfb", ink="#0b0b0b", ink2="#52514e", muted="#898781", grid="#e1e0d9", axis="#c3c2b7"),
         "dark": dict(surface="#1a1a19", ink="#ffffff", ink2="#c3c2b7", muted="#898781", grid="#2c2c2a", axis="#383835")}
SEQ = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]


def setup(mode):
    t = THEME[mode]
    plt.rcParams.update({"figure.facecolor": t["surface"], "axes.facecolor": t["surface"], "savefig.facecolor": t["surface"], "text.color": t["ink"],
                         "axes.labelcolor": t["ink2"], "xtick.color": t["muted"], "ytick.color": t["muted"], "axes.edgecolor": t["axis"], "axes.linewidth": 0.8,
                         "font.family": "sans-serif", "font.sans-serif": ["Segoe UI", "DejaVu Sans", "Arial"], "font.size": 11, "axes.titlesize": 15, "axes.titleweight": "bold",
                         "axes.titlelocation": "left", "axes.spines.top": False, "axes.spines.right": False, "axes.grid": False, "legend.frameon": False})
    return t, CAT[mode]


def save(fig, name, mode):
    fig.savefig(FIG / f"{name}_{mode}.png", dpi=200)
    plt.close(fig)


def grid(ax, t, axis="x"):
    ax.grid(axis=axis, color=t["grid"], linewidth=0.7)
    ax.set_axisbelow(True)


def fig_dumbbell(mode):
    t, C = setup(mode)
    d = pl.read_csv(RES / "biodomain_shares_dumbbell.csv").sort("gap_pp")
    n = d.height
    fig, ax = plt.subplots(figsize=(9.6, 5.4 + 0.12 * max(0, n - 12)))
    y = np.arange(n)
    a, g = d["attention_share"].to_numpy() * 100, d["genetic_share"].to_numpy() * 100
    for i in range(n):
        ax.plot([a[i], g[i]], [y[i], y[i]], color=t["axis"], lw=2, zorder=1, solid_capstyle="round")
    ax.scatter(a, y, s=70, color=C[0], zorder=3, edgecolor=t["surface"], linewidth=1.5, label="Share of AD abstracts (attention)")
    ax.scatter(g, y, s=70, color=C[1], zorder=3, edgecolor=t["surface"], linewidth=1.5, label="Share of genetic evidence (Open Targets)")
    ax.set_yticks(y); ax.set_yticklabels(d["domain"].to_list(), color=t["ink"], fontsize=10)
    xmax = max(a.max(), g.max())
    for i, r in enumerate(d.iter_rows(named=True)):
        sig = not r["gap_call"].startswith("within")
        txt = f"{r['gap_pp']:+.1f} pp" + (" *" if sig else "")
        ax.text(xmax * 1.04, y[i], txt, va="center", ha="left", fontsize=9, color=t["ink"] if sig else t["muted"], fontweight="bold" if sig else "normal")
    ax.set_xlim(0, xmax * 1.2); ax.set_xlabel("% of total (fractional biodomain membership)")
    grid(ax, t, "x")
    ax.set_title("Where the field looks vs. where the genes point", pad=24)
    ax.text(0, 1.02, "Alzheimer's biodomains: attention share vs. genetic-evidence share. Gap = genetic minus attention; * bootstrap 95% CI excludes 0 and |gap| >= 2 pp",
            transform=ax.transAxes, fontsize=9, color=t["ink2"], va="bottom")
    ax.legend(loc="lower right", fontsize=9)
    fig.tight_layout(); save(fig, "fig1_dumbbell_attention_vs_genetics", mode)


def fig_stream(mode):
    t, C = setup(mode)
    d = pl.read_csv(RES / "macro_topic_share_by_year.csv")
    mac = pl.read_csv(RES / "macro_topics.csv").sort("n_papers", descending=True)
    years = np.arange(1975, 2027)
    Y = np.zeros((len(mac), len(years)))
    for k, m in enumerate(mac["macro"].to_list()):
        s = d.filter(pl.col("macro") == m)
        for r in s.iter_rows(named=True):
            Y[k, r["year"] - 1975] = r["n"]
    tot = Y.sum(axis=0); tot[tot == 0] = 1
    S = Y / tot
    # 3-year centred smoothing for stability in the early sparse years (shares, not counts)
    Sm = np.copy(S)
    for j in range(len(years)):
        lo, hi = max(0, j - 1), min(len(years), j + 2)
        Sm[:, j] = Y[:, lo:hi].sum(axis=1) / Y[:, lo:hi].sum()
    fig, ax = plt.subplots(figsize=(9.6, 5.4))
    ax.stackplot(years, Sm, baseline="weighted_wiggle", colors=C[:len(mac)], edgecolor=t["surface"], linewidth=0.8,
                 labels=[f"{r['label']}  ({r['n_papers']:,})" for r in mac.iter_rows(named=True)])
    ax.axvspan(2025.5, 2026.5, color=t["surface"], alpha=0.6, zorder=3)
    ax.text(2026, ax.get_ylim()[1], "2026\npartial", ha="center", va="top", fontsize=7.5, color=t["muted"], zorder=4)
    ax.set_xlim(1975, 2026.5); ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_title("Alzheimer's literature, 1975-2026: topic share (macro-topics of the abstract map)", pad=8)
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), fontsize=8, title="Macro-topic (TF-IDF keywords; n papers)", title_fontsize=8)
    fig.tight_layout(); save(fig, "fig2_topic_stream", mode)


def fig_topics(mode):
    t, C = setup(mode)
    d = pl.read_csv(RES / "topic_share_by_year.csv").filter(pl.col("year") <= 2025)
    topics = d["topic"].unique(maintain_order=True).to_list()
    fig, ax = plt.subplots(figsize=(9.6, 5.4))
    ends = []
    for k, tp in enumerate(topics):
        s = d.filter(pl.col("topic") == tp).sort("year").filter(pl.col("year") >= 1980)
        ax.plot(s["year"], s["share"] * 100, color=C[k % 8], lw=2.2, label=tp)
        ends.append((s["share"][-1] * 100, tp, C[k % 8]))
    ends.sort()
    last_y = -99
    for v, tp, c in ends:
        yy = max(v, last_y + 1.6); last_y = yy
        ax.text(2025.6, yy, tp, fontsize=9, va="center", color=t["ink"])
        ax.plot([2025.1, 2025.5], [v, yy], color=c, lw=1)
    ax.set_xlim(1980, 2033); ax.set_ylabel("% of that year's AD abstracts mentioning the topic")
    grid(ax, t, "y"); ax.set_title("Share of AD abstracts by keyword topic (title or abstract)")
    ax.legend(ncol=4, loc="upper left", fontsize=8)
    fig.tight_layout(); save(fig, "fig2b_keyword_topics", mode)


def fig_lorenz(mode):
    t, C = setup(mode)
    L = pl.read_csv(RES / "lorenz_points.csv"); cg = pl.read_csv(RES / "concentration_gini.csv")
    fig, ax = plt.subplots(figsize=(7.2, 5.4))
    ax.plot([0, 1], [0, 1], color=t["axis"], lw=1.2, ls="--", label="Perfect equality")
    a = L.filter((pl.col("window") == "all_years") & (pl.col("series") == "attention"))
    g = L.filter((pl.col("window") == "all_years") & (pl.col("series") == "genetic_score"))
    b = L.filter((pl.col("window") == "1990-1999") & (pl.col("series") == "attention"))
    gi = lambda w, col="gini_attention_pc_with_zeros": cg.filter(pl.col("window") == w)[col][0]
    ax.plot(a["gene_frac"], a["cum_share"], color=C[0], lw=2.4, label=f"AD abstracts per gene, all years (Gini {gi('all_years'):.2f})")
    ax.plot(b["gene_frac"], b["cum_share"], color=C[2], lw=1.8, label=f"AD abstracts per gene, 1990-1999 (Gini {gi('1990-1999'):.2f})")
    ax.plot(g["gene_frac"], g["cum_share"], color=C[1], lw=2.4, label=f"Open Targets genetic score per gene (Gini {gi('all_years', 'gini_genetic_same_population'):.2f})")
    ax.set_xlabel("Cumulative share of protein-coding genes (least to most)"); ax.set_ylabel("Cumulative share of total")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); grid(ax, t, "both")
    ax.set_title("Attention is concentrated in few genes", pad=8)
    ax.legend(loc="upper left", fontsize=8.5)
    fig.tight_layout(); save(fig, "fig3_lorenz", mode)


def fig_under(mode):
    t, C = setup(mode)
    d = pl.read_csv(RES / "under_explored_top25.csv").sort("rank", descending=True)
    n = d.height
    fig, ax = plt.subplots(figsize=(9.6, 6.6))
    y = np.arange(n)
    A = d["A_total"].to_numpy() + 1
    cm = LinearSegmentedColormap.from_list("seq", SEQ[1:] if mode == "light" else SEQ[2:])
    gv = d["G_genetic"].to_numpy()
    for i in range(n):
        ax.plot([1, A[i]], [y[i], y[i]], color=t["axis"], lw=1.5, zorder=1)
    sc = ax.scatter(A, y, s=75, c=gv, cmap=cm, vmin=0.7, vmax=1.0, zorder=3, edgecolor=t["surface"], linewidth=1.2)
    ax.set_xscale("log"); ax.set_yticks(y); ax.set_yticklabels(d["symbol"].to_list(), color=t["ink"], fontsize=10)
    ax.axvline(100, color=C[1], lw=1.2, ls="--"); ax.text(105, n - 0.3, "100 AD abstracts", fontsize=8, color=t["ink2"], va="top")
    for i, r in enumerate(d.iter_rows(named=True)):
        ax.text(A[i] * 1.12, y[i], f"{r['A_total']:,}  |  G={r['G_genetic']:.2f}  |  first GWAS {r['y0'] if r['y0_source'] == 'gwas_catalog' else 'n/a'}", fontsize=8, va="center", color=t["ink2"])
    ax.set_xlim(1, A.max() * 60)
    ax.set_xlabel("AD abstracts mentioning the gene (PubTator3, log scale)")
    ax.set_title("Genetically supported genes with the fewest AD abstracts per year since first report", pad=8, fontsize=13)
    cb = fig.colorbar(sc, ax=ax, pad=0.01, fraction=0.025); cb.set_label("Open Targets genetic score", color=t["ink2"]); cb.outline.set_visible(False)
    grid(ax, t, "x")
    fig.tight_layout(); save(fig, "fig4_under_explored_top25", mode)


def fig_event(mode):
    t, C = setup(mode)
    d = pl.read_csv(RES / "event_study_counts.csv"); did = pl.read_csv(RES / "event_study_did.csv").filter(pl.col("event_year_source") == "protocol")
    genes = ["TREM2", "PLCG2", "ABI3"]
    fig, axes = plt.subplots(1, 3, figsize=(12.8, 5.0), sharey=False)
    for ax, gn in zip(axes, genes):
        s = d.filter(pl.col("gene") == gn).sort("t")
        if s.height == 0: continue
        ax.fill_between(s["t"], s["control_p25"], s["control_p75"], color=C[3], alpha=0.25, lw=0, label="20 matched genes without genetic evidence (IQR)")
        ax.plot(s["t"], s["control_mean"], color=C[3], lw=1.6)
        ax.plot(s["t"], s["count"], color=C[0], lw=2.6, marker="o", ms=4.5, label=f"{gn} AD abstracts")
        ax.axvline(0, color=t["ink2"], lw=1, ls="--")
        r = did.filter(pl.col("gene") == gn)
        ey = int(r["event_year"][0])
        ax.text(0.2, 0.97, f"first genetic report {ey}", transform=ax.get_xaxis_transform(), fontsize=8.5, color=t["ink2"], va="top")
        ax.set_title(gn, fontsize=13); ax.set_xlabel("years since first genetic report"); grid(ax, t, "y")
        if gn == genes[0]: ax.set_ylabel("AD abstracts mentioning the gene per year"); ax.legend(loc="center left", fontsize=7.5)
        ax.text(0.03, 0.55 if gn != "TREM2" else 0.45, f"post/pre x{r['post_over_pre_ratio'][0]:.1f}\nDiD {r['did_log1p'][0]:+.2f} log pts\n[{r['did_ci95_lo'][0]:+.2f}, {r['did_ci95_hi'][0]:+.2f}]", transform=ax.transAxes, fontsize=8.5, color=t["ink"], va="top")
    fig.suptitle("Attention after the genes were found: TREM2, PLCG2, ABI3", x=0.01, ha="left", fontsize=14, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.94)); save(fig, "fig5_event_study", mode)


def fig_pilot(mode):
    t, C = setup(mode)
    by = pl.read_csv(RES / "v1_corpus_by_year.csv").filter((pl.col("year") >= 1975) & (pl.col("year") <= 2025))
    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    ax.plot(by["year"], by["amyloid_share"] * 100, color=C[0], lw=2.4, label="Corpus: abstracts mentioning amyloid (regex, title/abstract)")
    pil = {1990: 13.8, 2000: 24.5, 2010: 33.4, 2020: 32.9, 2025: 31.0}
    ax.scatter(list(pil), list(pil.values()), s=60, color=C[1], zorder=3, edgecolor=t["surface"], linewidth=1.2, label="Pilot (PubMed field search)")
    ax.set_ylabel("% of AD abstracts"); ax.set_ylim(0, 45); grid(ax, t, "y")
    ax.set_title("Reproducing the pilot: about a third of AD abstracts mention amyloid", fontsize=13)
    ax.legend(loc="lower right", fontsize=8.5); fig.tight_layout(); save(fig, "fig7_pilot_amyloid", mode)


def fig_map(mode):
    t, C = setup(mode)
    p = pl.read_parquet(WORK / "papers_meta.parquet", columns=["year", "amyloid", "coordinate1", "coordinate2", "ot_genetic_max"])
    g = pl.read_parquet(WORK / "genes_meta.parquet", columns=["symbol", "coordinate1", "coordinate2", "G_genetic", "A_g", "strong_G07", "sparse_G07_A100"])
    fig, ax = plt.subplots(figsize=(9.6, 7.2))
    ax.scatter(p["coordinate1"], p["coordinate2"], s=0.25, color=t["axis"] if mode == "light" else "#5a5a55", alpha=0.35, linewidths=0, rasterized=True, label="AD abstracts (all)")
    am = p.filter(pl.col("amyloid"))
    ax.scatter(am["coordinate1"], am["coordinate2"], s=0.3, color=C[0], alpha=0.35, linewidths=0, rasterized=True, label="Abstracts mentioning amyloid")
    gs = g.filter(pl.col("G_genetic") >= 0.7)
    sp = gs.filter(pl.col("sparse_G07_A100")); st = gs.filter(~pl.col("sparse_G07_A100"))
    ax.scatter(st["coordinate1"], st["coordinate2"], s=26, color=C[1], edgecolor=t["surface"], linewidth=0.8, label="Genes, genetic score >= 0.7 (>= 100 abstracts)", zorder=4)
    ax.scatter(sp["coordinate1"], sp["coordinate2"], s=40, color=C[7], edgecolor=t["surface"], linewidth=0.8, marker="D", label="Genes, genetic score >= 0.7 and < 100 abstracts", zorder=5)
    for r in gs.sort("G_genetic", descending=True).head(14).iter_rows(named=True):
        ax.text(r["coordinate1"], r["coordinate2"], " " + r["symbol"], fontsize=7.5, color=t["ink"], zorder=6)
    ax.set_xticks([]); ax.set_yticks([]); [s.set_visible(False) for s in ax.spines.values()]
    ax.set_title("Map of 237k Alzheimer's abstracts with genes in the same embedding space (exploratory)", fontsize=12, pad=6)
    ax.legend(loc="lower left", fontsize=8, markerscale=1.3)
    fig.tight_layout(); save(fig, "fig6_map_literature_and_genes", mode)


if __name__ == "__main__":
    which = sys.argv[1:] or ["dumbbell", "stream", "topics", "lorenz", "under", "event", "pilot", "map"]
    fns = dict(dumbbell=fig_dumbbell, stream=fig_stream, topics=fig_topics, lorenz=fig_lorenz, under=fig_under, event=fig_event, pilot=fig_pilot, map=fig_map)
    for w in which:
        for mode in ("light", "dark"):
            try:
                fns[w](mode); print("ok", w, mode, flush=True)
            except Exception as e:
                print("FAIL", w, mode, type(e).__name__, e, flush=True)
