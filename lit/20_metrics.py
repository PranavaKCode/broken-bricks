#!/usr/bin/env python
"""Stage B: pre-registered literature-vs-genetics metrics (preregistration D1-D6, L1-L5, S1-S7).
Every number lands in out/results/*.csv.  Seed 20261011."""
import os, json, re, datetime, math
from pathlib import Path
import numpy as np
import polars as pl
import scipy.sparse as sp
from scipy.stats import spearmanr

SEED = 20261011
rng = np.random.default_rng(SEED)
LIT = Path(os.path.expanduser("~/orcd/pool/mantis-films/research/lit"))
RAW, WORK, RES = LIT / "raw", LIT / "work", LIT / "out" / "results"
RES.mkdir(parents=True, exist_ok=True)
LAST = 2025
NB = 1000
def log(*a): print(*a, flush=True)
def wcsv(df, name):
    df.write_csv(RES / name); log(f"  wrote {name} {df.shape}")

# ---------------------------------------------------------------------------------- load
hum = pl.read_parquet(WORK / "human_genes.parquet")
ga = pl.read_parquet(WORK / "gene_attention.parquet")
corp = pl.read_parquet(WORK / "corpus.parquet")
pairs = pl.read_parquet(WORK / "pairs_human.parquet").filter(pl.col("abs_verified"))   # primary attention = abstract-verified pairs (amendment A2)
Ay = pl.read_parquet(WORK / "A_gene_year.parquet")
Ay_s4 = pl.read_parquet(WORK / "A_gene_year_s4.parquet")
Ay_g2p = pl.read_parquet(WORK / "A_gene_year_gene2pubmed.parquet")
Ay_tiab = pl.read_parquet(WORK / "A_gene_year_tiab.parquet")
Ay_literal = pl.read_parquet(WORK / "A_gene_year_literal.parquet")
Ay_strict = pl.read_parquet(WORK / "A_gene_year_strict.parquet")
n_corpus = corp.height
gids = ga["gene_id"].to_list()
gpos = {g: i for i, g in enumerate(gids)}
NG = len(gids)
sym = ga["symbol"].to_list()
s2g = {s: g for s, g in zip(sym, gids)}
e2g = {e: g for e, g in zip(ga["ensembl_id"].to_list(), gids) if e}
is_pc = (ga["type_of_gene"] == "protein-coding").to_numpy()
log("genes", NG, "protein-coding", int(is_pc.sum()), "corpus", n_corpus)

YEARS = np.arange(1975, 2027)
def year_matrix(df):
    M = np.zeros((NG, len(YEARS)), dtype=np.float64)
    d = df.filter((pl.col("year") >= 1975) & (pl.col("year") <= 2026))
    gi = np.array([gpos[g] for g in d["gene_id"].to_list()], dtype=np.int64)
    M[gi, d["year"].to_numpy() - 1975] = d["n"].to_numpy()
    return M
M_lit = year_matrix(Ay)
M_literal = year_matrix(Ay_literal)
M_strict = year_matrix(Ay_strict)
M_s4 = year_matrix(Ay_s4)
M_g2p = year_matrix(Ay_g2p)
M_tiab = year_matrix(Ay_tiab)
yidx = lambda y0, y1: slice(int(y0 - 1975), int(y1 - 1975 + 1))
A_tot = ga["A_g"].to_numpy().astype(float)   # whole period, all years incl. pre-1975 and 2026 (literal D1)

# ---------------------------------------------------------------------------------- genetics (D2)
ot = pl.read_parquet(RAW / "ot" / "ot_assoc_indirect.parquet")
otd = pl.read_parquet(RAW / "ot" / "ot_assoc_direct.parquet")
def ot_to_G(df, col="dt_genetic_association"):
    G = np.zeros(NG)
    unm = 0
    for e, s, v in zip(df["ensembl_id"].to_list(), df["symbol"].to_list(), df[col].fill_null(0.0).to_list()):
        g = e2g.get(e, s2g.get(s))
        if g is None:
            unm += 1 if v > 0 else 0
            continue
        G[gpos[g]] = max(G[gpos[g]], v)
    return G, unm
G, unm_ind = ot_to_G(ot)
Gd, unm_dir = ot_to_G(otd)
log("G>0", int((G > 0).sum()), ">=0.5", int((G >= 0.5).sum()), ">=0.7", int((G >= 0.7).sum()), "unmapped with G>0:", unm_ind)

# GWAS first-report years
gw = pl.read_parquet(RAW / "gwas" / "gwas_ad_associations.parquet")
py = pl.read_csv(RAW / "gwas" / "gwas_pmid_years.csv", infer_schema_length=0)
def yr_of(row):
    for k in ("sortpubdate", "pubdate", "epubdate"):
        m = re.search(r"(1[89]\d\d|20\d\d)", row.get(k) or "")
        if m: return int(m.group(1))
    return None
pmid_year = {int(r["pmid"]): yr_of(r) for r in py.to_dicts()}
syn = {}
for g, s, sy in zip(gids, sym, ga["Synonyms"].fill_null("").to_list()):
    for x in sy.split("|"):
        if x and x != "-": syn.setdefault(x, []).append(g)
def sym_to_gid(s):
    if s in s2g: return s2g[s]
    c = syn.get(s)
    return c[0] if c and len(c) == 1 else None
first = {}; first_pmid = {}; n_assoc_g = {}
unmapped_syms = set()
for r in gw.to_dicts():
    y = pmid_year.get(int(r["pubmed_id"])) if r["pubmed_id"] and str(r["pubmed_id"]).isdigit() else None
    for s in (r["mapped_genes"] or "").split(";"):
        if not s: continue
        g = sym_to_gid(s)
        if g is None: unmapped_syms.add(s); continue
        n_assoc_g[g] = n_assoc_g.get(g, 0) + 1
        if y is not None and (g not in first or y < first[g]):
            first[g] = y; first_pmid[g] = r["pubmed_id"]
log("GWAS genes with first-report year", len(first), "unmapped mapped_gene symbols", len(unmapped_syms))
y0 = np.full(NG, 1975.0); y0_src = np.array(["no_gwas"] * NG, dtype=object)
y0_gwas = np.full(NG, np.nan)
for g, y in first.items():
    y0_gwas[gpos[g]] = y
    y0[gpos[g]] = min(max(y, 1975), LAST); y0_src[gpos[g]] = "gwas_catalog"
Yg = np.maximum(1, LAST - y0 + 1)

def Asince_of(M):
    out = np.zeros(NG)
    for i in np.where(G > 0)[0]:
        out[i] = M[i, yidx(max(y0[i], 1975), LAST)].sum()
    return out
Asince = Asince_of(M_lit)
with np.errstate(divide="ignore", invalid="ignore"):
    E = np.where(G > 0, (Asince / Yg) / G, np.nan)

# ---------------------------------------------------------------------------------- Agora (D3)
ag = {}
agp = RAW / "agora" / "agora_genes.jsonl"
if agp.exists():
    for line in open(agp):
        d = json.loads(line)
        if d.get("status") == 200:
            ag[d["ensembl_id"]] = d
log("Agora genes with data", len(ag))
dom_list = [[] for _ in range(NG)]
treatad = np.full(NG, np.nan)
for e, d in ag.items():
    g = e2g.get(e)
    if g is None: continue
    dom_list[gpos[g]] = sorted(set(d.get("biodomains") or []))
    if d.get("genetics_score") is not None: treatad[gpos[g]] = d["genetics_score"]
for e, s in zip(ot["ensembl_id"].to_list(), ot["symbol"].to_list()):
    g = e2g.get(e, s2g.get(s))   # OT genes without NCBI xref but present in Agora
    if g is None and e in ag:
        pass
n_dom = np.array([len(x) for x in dom_list])
in_U = (G > 0) | (A_tot >= 20)
assigned = in_U & (n_dom > 0)
domains = sorted({d for i in np.where(assigned)[0] for d in dom_list[i]})
log("universe", int(in_U.sum()), "assigned", int(assigned.sum()), "domains", len(domains))

# ---------------------------------------------------------------------------------- gene master table
gm = pl.DataFrame(dict(
    gene_id=gids, symbol=sym, ensembl_id=ga["ensembl_id"].to_list(), type_of_gene=ga["type_of_gene"].to_list(),
    A_g=A_tot.astype(int), A_2010_2025=ga["A_2010_2025"].to_list(), A_g_literal_D1=ga["A_g_literal"].to_list(), A_g_strict_A1=ga["A_g_strict"].to_list(), A_g_ortholog=ga["A_g_s4"].to_list(),
    A_g_gene2pubmed=ga["A_g_gene2pubmed"].to_list(), A_g_tiab=ga["A_g_tiab"].to_list(),
    G_genetic=G, G_genetic_direct=Gd, treatad_genetics_score=treatad,
    y0_gwas=y0_gwas, y0_used=y0.astype(int), y0_source=y0_src.tolist(), Y_g=Yg.astype(int), Asince_g=Asince.astype(int), E_g=E,
    n_biodomains=n_dom, biodomains="|".join([]) if False else ["|".join(x) for x in dom_list],
    in_universe=in_U, strong_G07=(G >= 0.7), moderate_G05=(G >= 0.5),
    sparse_G07_A100=((G >= 0.7) & (A_tot < 100)),
))
gm = gm.sort("A_g", descending=True)
wcsv(gm.filter((pl.col("A_g") > 0) | (pl.col("G_genetic") > 0)), "gene_master.csv")

# ---------------------------------------------------------------------------------- D4 concentration
def gini(x):
    x = np.sort(np.asarray(x, dtype=float)); n = len(x); c = np.cumsum(x)
    return float((n + 1 - 2 * c.sum() / c[-1]) / n) if c[-1] > 0 else float("nan")
def lorenz_pts(x, n_pts=2000):
    x = np.sort(np.asarray(x, dtype=float)); c = np.cumsum(x) / x.sum(); n = len(x)
    idx = np.unique(np.r_[np.linspace(0, n - 1, n_pts).astype(int), np.arange(max(0, n - 300), n)])
    return (idx + 1) / n, c[idx]
WIN = {"all_years": (None, None), "1975-1989": (1975, 1989), "1990-1999": (1990, 1999), "2000-2009": (2000, 2009), "2010-2019": (2010, 2019), "2020-2025": (2020, 2025), "2010-2025": (2010, 2025)}
pc_idx = np.where(is_pc)[0]
# paper x gene incidence for Poisson bootstrap (protein-coding only)
p_idx = {p: i for i, p in enumerate(corp["pmid"].to_list())}
pp = pairs.filter(pl.col("gene_id").is_in([gids[i] for i in pc_idx]))
pc_pos = {gids[i]: k for k, i in enumerate(pc_idx)}
rows_i = np.array([p_idx[p] for p in pp["pmid"].to_list()]); cols_i = np.array([pc_pos[g] for g in pp["gene_id"].to_list()])
P = sp.csr_matrix((np.ones(len(rows_i)), (rows_i, cols_i)), shape=(n_corpus, len(pc_idx)))
corp_year = corp["year"].to_numpy()
conc, lor = [], []
G_pc = G[pc_idx]
for w, (a, b) in WIN.items():
    mask = np.ones(n_corpus, bool) if a is None else ((corp_year >= a) & (corp_year <= b))
    cnt = np.asarray(P.T @ mask.astype(float)).ravel()
    gA, gA1 = gini(cnt), gini(cnt[cnt >= 1]) if (cnt >= 1).sum() > 1 else float("nan")
    boots = []
    for _ in range(NB):
        wts = rng.poisson(1.0, n_corpus) * mask
        boots.append(gini(np.asarray(P.T @ wts.astype(float)).ravel()))
    lo, hi = np.nanpercentile(boots, [2.5, 97.5])
    conc.append(dict(window=w, n_papers=int(mask.sum()), n_pc_genes=len(pc_idx), n_genes_A_ge1=int((cnt >= 1).sum()), gini_attention_pc_with_zeros=gA,
                     ci95_lo=lo, ci95_hi=hi, gini_attention_genes_ge1=gA1, gini_genetic_same_population=gini(G_pc) if w == "all_years" else None,
                     top1_gene=sym[pc_idx[int(np.argmax(cnt))]], top1_share=float(cnt.max() / cnt.sum()),
                     top10_share=float(np.sort(cnt)[-10:].sum() / cnt.sum()), top100_share=float(np.sort(cnt)[-100:].sum() / cnt.sum())))
    x, c = lorenz_pts(cnt)
    lor += [dict(window=w, series="attention", gene_frac=a_, cum_share=b_) for a_, b_ in zip(x, c)]
x, c = lorenz_pts(G_pc); lor += [dict(window="all_years", series="genetic_score", gene_frac=a_, cum_share=b_) for a_, b_ in zip(x, c)]
wcsv(pl.DataFrame(conc), "concentration_gini.csv"); wcsv(pl.DataFrame(lor), "lorenz_points.csv")
gA_all = [r for r in conc if r["window"] == "all_years"][0]
L1_hold = gA_all["gini_attention_pc_with_zeros"] > gA_all["gini_genetic_same_population"]
log("Gini attention", gA_all["gini_attention_pc_with_zeros"], "genetic", gA_all["gini_genetic_same_population"], "L1 holds:", L1_hold)

# ---------------------------------------------------------------------------------- L2 coupling
def boot_spearman(x, y):
    r = spearmanr(x, y).statistic
    n = len(x); bs = []
    for _ in range(NB):
        i = rng.integers(0, n, n); bs.append(spearmanr(x[i], y[i]).statistic)
    return r, np.nanpercentile(bs, [2.5, 97.5])
pos = np.where(G > 0)[0]
r_all, ci_all = boot_spearman(G[pos], A_tot[pos])
top32 = pos[np.argsort(-G[pos])[:32]]
r32, ci32 = boot_spearman(G[top32], A_tot[top32])
r07, ci07 = boot_spearman(G[G >= 0.5], A_tot[G >= 0.5]) if (G >= 0.5).sum() > 5 else (np.nan, (np.nan, np.nan))
coup = pl.DataFrame([dict(set="genes with G>0", n=len(pos), spearman=r_all, ci95_lo=ci_all[0], ci95_hi=ci_all[1]),
                     dict(set="top 32 by G (pilot comparator, pilot rho=0.41)", n=32, spearman=r32, ci95_lo=ci32[0], ci95_hi=ci32[1]),
                     dict(set="genes with G>=0.5", n=int((G >= 0.5).sum()), spearman=r07, ci95_lo=ci07[0], ci95_hi=ci07[1])])
wcsv(coup, "coupling_spearman.csv")
L2_hold = 0 < r_all < 0.5

# ---------------------------------------------------------------------------------- D3 biodomain shares
dom_pos = {d: k for k, d in enumerate(domains)}
rr, cc, vv = [], [], []
for i in np.where(assigned)[0]:
    for d in dom_list[i]:
        rr.append(i); cc.append(dom_pos[d]); vv.append(1.0 / n_dom[i])
W = sp.csr_matrix((vv, (rr, cc)), shape=(NG, len(domains)))
aidx = np.where(assigned)[0]
def shares(v, mask_idx=aidx, weights=None):
    vv_ = v.copy().astype(float)
    if weights is not None: vv_ = vv_ * weights
    tot = vv_[mask_idx].sum()
    return (W.T @ np.where(np.isin(np.arange(NG), mask_idx), vv_, 0.0)) / tot if tot > 0 else np.full(len(domains), np.nan)
G_for_shares = G.copy()
sA = shares(A_tot); sG = shares(G_for_shares)
# bootstrap over genes
boot_gap = np.zeros((NB, len(domains)))
for b in range(NB):
    cnt = np.bincount(rng.integers(0, len(aidx), len(aidx)), minlength=len(aidx)).astype(float)
    wv = np.zeros(NG); wv[aidx] = cnt
    boot_gap[b] = shares(G_for_shares, weights=wv) - shares(A_tot, weights=wv)
gap = sG - sA
gl, gh = np.percentile(boot_gap, [2.5, 97.5], axis=0)
# paper-level companion
pa = pairs.filter(pl.col("gene_id").is_in([gids[i] for i in aidx])).select(["pmid", "gene_id"])
dom_rows = [dict(gene_id=gids[i], domain=d) for i in aidx for d in dom_list[i]]
pa = pa.join(pl.DataFrame(dom_rows), on="gene_id").select(["pmid", "domain"]).unique()
paper_frac = {r["domain"]: r["n"] / n_corpus for r in pa.group_by("domain").agg(pl.len().alias("n")).to_dicts()}
unassigned = in_U & (n_dom == 0)
dtab = pl.DataFrame(dict(
    domain=domains, n_genes=[int(W[:, k].nnz) for k in range(len(domains))],
    attention_share=sA, genetic_share=sG, gap_pp=100 * gap, gap_ci95_lo_pp=100 * gl, gap_ci95_hi_pp=100 * gh,
    paper_level_fraction_mentioning=[paper_frac.get(d, 0.0) for d in domains],
)).with_columns(
    pl.when((pl.col("gap_pp") >= 2) & (pl.col("gap_ci95_lo_pp") > 0)).then(pl.lit("under-attended (genetic > attention)"))
      .when((pl.col("gap_pp") <= -2) & (pl.col("gap_ci95_hi_pp") < 0)).then(pl.lit("over-attended (attention > genetic)")).otherwise(pl.lit("within +/-2pp or CI crosses 0")).alias("gap_call"))
dtab = dtab.sort("gap_pp", descending=True)
wcsv(dtab, "biodomain_shares_dumbbell.csv")
dec_rows = []
for w, (a, b) in WIN.items():
    if w in ("all_years", "2010-2025"): continue
    Aw = M_lit[:, yidx(a, b)].sum(axis=1)
    sw = shares(Aw)
    for k, d in enumerate(domains):
        dec_rows.append(dict(window=w, domain=d, attention_share=sw[k], genetic_share=sG[k], gap_pp=100 * (sG[k] - sw[k])))
wcsv(pl.DataFrame(dec_rows), "biodomain_shares_by_decade.csv")
ua = dict(item=["universe genes (G>0 or A>=20)", "assigned to >=1 biodomain", "unassigned (no Agora biodomain)", "attention mass in unassigned (share of universe A)", "genetic mass in unassigned (share of universe G)"],
          value=[int(in_U.sum()), int(assigned.sum()), int(unassigned.sum()), float(A_tot[unassigned].sum() / A_tot[in_U].sum()), float(G[unassigned].sum() / G[in_U].sum())])
wcsv(pl.DataFrame(ua), "biodomain_unassigned.csv")
L3_under = dtab.filter(pl.col("gap_call").str.starts_with("under")).height
L3_over = dtab.filter(pl.col("gap_call").str.starts_with("over")).height
log("L3: under-attended domains", L3_under, "over-attended", L3_over)
# S3: universe sensitivity (G>0 only)
in_U3 = G > 0; assigned3 = in_U3 & (n_dom > 0); a3 = np.where(assigned3)[0]
sA3 = shares(A_tot, mask_idx=a3); sG3 = shares(G, mask_idx=a3)
wcsv(pl.DataFrame(dict(domain=domains, gap_pp_primary=100 * gap, gap_pp_universe_G_gt0=100 * (sG3 - sA3))).sort("gap_pp_primary", descending=True), "s3_universe_sensitivity.csv")

# ---------------------------------------------------------------------------------- D5 under-explored list
def under_list(Gv, Mv, thr=0.7, n=25, Asin_cached=None):
    A_tot_v = Mv.sum(axis=1)
    asin = np.zeros(NG)
    for i in np.where(Gv > 0)[0]:
        asin[i] = Mv[i, yidx(max(y0[i], 1975), LAST)].sum()
    with np.errstate(divide="ignore", invalid="ignore"):
        e = np.where(Gv > 0, (asin / Yg) / Gv, np.nan)
    cand = np.where(Gv >= thr)[0]
    return cand[np.argsort(e[cand], kind="stable")][:n], e
top25_idx = np.where(G >= 0.7)[0]
order = top25_idx[np.argsort(E[top25_idx], kind="stable")]
under = gm.filter(pl.col("G_genetic") >= 0.5).with_columns(pl.col("E_g").rank("ordinal").alias("rank_all_G05")).sort("E_g")
ue = pl.DataFrame(dict(rank=np.arange(1, len(order) + 1), gene_id=[gids[i] for i in order], symbol=[sym[i] for i in order], G_genetic=G[order], A_total=A_tot[order].astype(int),
                       A_since_y0=Asince[order].astype(int), y0=y0[order].astype(int), y0_source=y0_src[order].tolist(), Y_g=Yg[order].astype(int), E_g=E[order],
                       biodomains=["|".join(dom_list[i]) for i in order], sparse_A_lt_100=(A_tot[order] < 100)))
wcsv(ue, "under_explored_G07_ranked.csv")
wcsv(ue.head(25), "under_explored_top25.csv")
sparse_n = int(((G >= 0.7) & (A_tot < 100)).sum())
# contrast ratio
amax_i = int(np.argmax(A_tot)); studied = np.where((G >= 0.7) & (A_tot >= 1))[0]
least_i = studied[np.argmin(A_tot[studied])]
zero_i = np.where((G >= 0.7) & (A_tot == 0))[0]
contrast = pl.DataFrame([dict(item="most-mentioned gene (any)", gene=sym[amax_i], A=int(A_tot[amax_i])),
                         dict(item="least-mentioned gene with G>=0.7 and A>=1", gene=sym[least_i], A=int(A_tot[least_i])),
                         dict(item="ratio most-mentioned / least-mentioned strongly supported", gene=f"{sym[amax_i]}/{sym[least_i]}", A=float(A_tot[amax_i] / A_tot[least_i])),
                         dict(item="genes with G>=0.7 and A=0", gene=";".join(sym[i] for i in zero_i), A=len(zero_i)),
                         dict(item="genes with G>=0.7 and A<100 (sparsely studied)", gene="", A=sparse_n),
                         dict(item="genes with G>=0.7", gene="", A=int((G >= 0.7).sum())),
                         dict(item="genes with G>=0.5 and A<100", gene="", A=int(((G >= 0.5) & (A_tot < 100)).sum()))])
wcsv(contrast, "contrast_and_headline_counts.csv")
L4_hold = sparse_n >= 10
# amyloid-pathway comparison (S7)
AMY_SET = ["APP", "PSEN1", "PSEN2", "BACE1", "BACE2", "ADAM10", "APH1A", "APH1B", "NCSTN", "PSENEN"]
am = np.array([gpos[s2g[s]] for s in AMY_SET if s in s2g])
app_dom = np.array([i for i in range(NG) if "APP Metabolism" in dom_list[i]])
strong = np.where(G >= 0.7)[0]
strong_nonamy = np.setdiff1d(strong, am)
s7 = pl.DataFrame([
    dict(group="genes with G>=0.7 (all)", n=len(strong), mean_A=float(A_tot[strong].mean()), median_A=float(np.median(A_tot[strong]))),
    dict(group="genes with G>=0.7 excluding amyloid-pathway set", n=len(strong_nonamy), mean_A=float(A_tot[strong_nonamy].mean()), median_A=float(np.median(A_tot[strong_nonamy]))),
    dict(group="amyloid-pathway set (APP,PSEN1,PSEN2,BACE1,BACE2,ADAM10,APH1A,APH1B,NCSTN,PSENEN)", n=len(am), mean_A=float(A_tot[am].mean()), median_A=float(np.median(A_tot[am]))),
    dict(group="Agora 'APP Metabolism' biodomain genes in universe", n=int(np.intersect1d(app_dom, np.where(in_U)[0]).size), mean_A=float(A_tot[np.intersect1d(app_dom, np.where(in_U)[0])].mean()), median_A=float(np.median(A_tot[np.intersect1d(app_dom, np.where(in_U)[0])]))),
])
wcsv(s7, "s7_amyloid_pathway_comparison.csv")

# ---------------------------------------------------------------------------------- D6 event study
EV = {"TREM2": 2013, "PLCG2": 2017, "ABI3": 2017}
pre_sl, post_sl = (-3, -1), (1, 5)
Mlog = np.log1p(M_lit)
ctrl_pool = np.where(is_pc & (G == 0))[0]
ev_counts, ev_did = [], []
def run_event(gname, ey, tag):
    i = gpos[s2g[gname]]
    yrs = np.arange(max(1975, ey - 10), min(LAST, ey + 10) + 1)
    pre_y = np.arange(ey + pre_sl[0], ey + pre_sl[1] + 1); post_y = np.arange(ey + post_sl[0], min(LAST, ey + post_sl[1]) + 1)
    if post_y.size == 0 or pre_y.min() < 1975: return
    pre_cnt = M_lit[i, pre_y - 1975].mean(); post_cnt = M_lit[i, post_y - 1975].mean()
    t_pre = Mlog[i, pre_y - 1975].mean(); t_post = Mlog[i, post_y - 1975].mean()
    c_pre_mean = Mlog[ctrl_pool][:, pre_y - 1975].mean(axis=1)
    order_ = np.lexsort((rng.permutation(len(ctrl_pool)), np.abs(c_pre_mean - t_pre)))[:20]
    ctrl = ctrl_pool[order_]
    c_pre = Mlog[ctrl][:, pre_y - 1975].mean(axis=1); c_post = Mlog[ctrl][:, post_y - 1975].mean(axis=1)
    did = (t_post - t_pre) - (c_post - c_pre).mean()
    bs = [(t_post - t_pre) - (c_post - c_pre)[rng.integers(0, 20, 20)].mean() for _ in range(NB)]
    lo, hi = np.percentile(bs, [2.5, 97.5])
    ev_did.append(
                  dict(gene=gname, event_year=ey, event_year_source=tag, G_genetic=float(G[i]), pre_mean_count=float(pre_cnt), post_mean_count=float(post_cnt),
                       post_over_pre_ratio=float(post_cnt / pre_cnt) if pre_cnt > 0 else None, treated_delta_log1p=float(t_post - t_pre), control_delta_log1p=float((c_post - c_pre).mean()),
                       did_log1p=float(did), did_ci95_lo=float(lo), did_ci95_hi=float(hi), n_controls=20, control_pre_mean_count=float(M_lit[ctrl][:, pre_y - 1975].mean()),
                       n_post_years=int(post_y.size)))
    if tag == "protocol":
        for y in yrs:
            cm = M_lit[ctrl][:, y - 1975]
            ev_counts.append(dict(gene=gname, t=int(y - ey), year=int(y), count=float(M_lit[i, y - 1975]), control_mean=float(cm.mean()), control_p25=float(np.percentile(cm, 25)), control_p75=float(np.percentile(cm, 75))))
for gname, ey in EV.items():
    run_event(gname, ey, "protocol")
    gy = y0_gwas[gpos[s2g[gname]]]
    if not np.isnan(gy): run_event(gname, int(gy), "gwas_catalog_first_report")
did_df = pl.DataFrame(ev_did, infer_schema_length=None)
wcsv(did_df, "event_study_did.csv"); wcsv(pl.DataFrame(ev_counts), "event_study_counts.csv")
tr = did_df.filter((pl.col("gene") == "TREM2") & (pl.col("event_year_source") == "protocol")).to_dicts()[0]
L5_hold = (tr["post_over_pre_ratio"] is not None and tr["post_over_pre_ratio"] >= 5) and tr["did_ci95_lo"] > 0
log("TREM2 event:", tr)

# ---------------------------------------------------------------------------------- S1, S5, S6 sensitivity
defs = {"pubtator_abstract_verified_PRIMARY": M_lit, "pubtator_literal_D1": M_literal, "pubtator_strict_verified_A1": M_strict, "pubtator_plus_ortholog": M_s4, "gene2pubmed": M_g2p, "tiab_symbol": M_tiab}
lists, srows = {}, []
for k, Mv in defs.items():
    idx_, ev = under_list(G, Mv)
    lists[k] = [sym[i] for i in idx_]
base = set(lists["pubtator_abstract_verified_PRIMARY"])
anyc = {k: Mv.sum(axis=1) for k, Mv in defs.items()}
anyc["pubtator_abstract_verified_PRIMARY"] = A_tot
for k in defs:
    a, b = anyc["pubtator_abstract_verified_PRIMARY"], anyc[k]
    msk = (a > 0) | (b > 0)
    srows.append(dict(analysis="S1 attention definition", variant=k, spearman_A_vs_primary=float(spearmanr(a[msk], b[msk]).statistic), spearman_A_vs_primary_G_gt0=float(spearmanr(a[G > 0], b[G > 0]).statistic),
                      top25_jaccard_vs_primary=len(base & set(lists[k])) / len(base | set(lists[k])), top25_list=";".join(lists[k])))
# S5 genetic source variants
for name, Gv in (("OT_direct_only", Gd),):
    idx_, _ = under_list(Gv, M_lit); l = [sym[i] for i in idx_]
    srows.append(dict(analysis="S5 genetic source", variant=name, spearman_A_vs_primary=None, spearman_A_vs_primary_G_gt0=float(spearmanr(Gd[Gd > 0], A_tot[Gd > 0]).statistic),
                      top25_jaccard_vs_primary=len(base & set(l)) / len(base | set(l)), top25_list=";".join(l)))
ta = np.nan_to_num(treatad, nan=0.0)
nstrong = int((G >= 0.7).sum())
thr_ta = np.sort(ta)[-nstrong] if nstrong > 0 else np.inf
Gt = np.where(ta >= thr_ta, 1.0, 0.0) * np.clip(ta / (ta.max() or 1), 0.01, 1)
idx_, _ = under_list(np.where(ta >= thr_ta, 0.7 + 0.3 * ta / (ta.max() or 1), 0.0), M_lit)
l = [sym[i] for i in idx_]
srows.append(dict(analysis="S5 genetic source", variant=f"TREAT-AD genetics_score top-{nstrong} (as strong set)", spearman_A_vs_primary=None,
                  spearman_A_vs_primary_G_gt0=float(spearmanr(ta[ta > 0], A_tot[ta > 0]).statistic) if (ta > 0).sum() > 5 else None,
                  top25_jaccard_vs_primary=len(base & set(l)) / len(base | set(l)), top25_list=";".join(l)))
srows.append(dict(analysis="S5 genetic source", variant="Spearman(G primary, TREAT-AD genetics) among genes with both", spearman_A_vs_primary=float(spearmanr(G[(G > 0) & (ta > 0)], ta[(G > 0) & (ta > 0)]).statistic) if ((G > 0) & (ta > 0)).sum() > 5 else None,
                  spearman_A_vs_primary_G_gt0=None, top25_jaccard_vs_primary=None, top25_list=""))
wcsv(pl.DataFrame(srows, infer_schema_length=None), "s1_s5_sensitivity.csv")
sw = []
for t in (0.5, 0.6, 0.7, 0.8):
    for a_ in (50, 100, 200):
        sw.append(dict(G_threshold=t, A_threshold=a_, n_genes_G_ge_t=int((G >= t).sum()), n_genes_G_ge_t_and_A_lt_a=int(((G >= t) & (A_tot < a_)).sum())))
wcsv(pl.DataFrame(sw), "s6_threshold_sweep.csv")
# S2 earliest year definition: amyloid shares and TREM2
c2 = corp.group_by("year_earliest").agg(pl.len().alias("n"), pl.col("amyloid").sum().alias("na"), pl.col("trem2_flag").sum().alias("nt")).sort("year_earliest")
wcsv(c2.with_columns((pl.col("na") / pl.col("n")).alias("amyloid_share")).rename({"year_earliest": "year"}), "s2_earliest_year_corpus.csv")

# ---------------------------------------------------------------------------------- topic shares (stream graph data)
TOPICS = {"amyloid": r"\bamyloid", "tau": r"\btau\b|\bMAPT\b|neurofibrillary", "microglia/neuroinflammation": r"microglia|neuroinflamm", "lipid/cholesterol": r"cholesterol|\blipid|lipoprotein",
          "endosome/lysosome": r"endosom|lysosom|endocyt", "myelin/oligodendrocyte": r"myelin|oligodendro", "vascular": r"vascular|cerebrovascular|stroke|hypertension"}
txt = (corp["title"].fill_null("") + " " + corp["abstract"].fill_null("")).to_list()
tr_rows = []
flags = {k: np.array([bool(re.search(p, t, re.I)) for t in txt]) for k, p in TOPICS.items()}
for k, f in flags.items():
    s = corp.with_columns(pl.Series("f", f)).group_by("year").agg(pl.len().alias("n"), pl.col("f").sum().alias("k")).sort("year")
    tr_rows += [dict(topic=k, year=r["year"], n_papers=r["n"], n_topic=r["k"], share=r["k"] / r["n"]) for r in s.to_dicts() if 1975 <= r["year"] <= 2026]
wcsv(pl.DataFrame(tr_rows), "topic_share_by_year.csv")
pilot_counts = {k: int(f.sum()) for k, f in flags.items()}
wcsv(pl.DataFrame([dict(topic=k, n_papers=v, share=v / n_corpus) for k, v in pilot_counts.items()]), "topic_totals.csv")

# ---------------------------------------------------------------------------------- scorecard
card = [dict(id="L1", statement="Gini of attention > Gini of genetic evidence (protein-coding genes, all years)", observed=f"{gA_all['gini_attention_pc_with_zeros']:.3f} vs {gA_all['gini_genetic_same_population']:.3f}", held=bool(L1_hold)),
        dict(id="L2", statement="Spearman(G, A) among genes with G>0 positive and < 0.5", observed=f"{r_all:.3f} (95% CI {ci_all[0]:.3f} to {ci_all[1]:.3f}), n={len(pos)}", held=bool(L2_hold)),
        dict(id="L3", statement=">=1 domain with gap >= +2pp (CI>0) and >=1 with gap <= -2pp (CI<0)", observed=f"{L3_under} under-attended, {L3_over} over-attended", held=bool(L3_under >= 1 and L3_over >= 1)),
        dict(id="L4", statement=">=10 genes with G>=0.7 and A<100", observed=str(sparse_n), held=bool(L4_hold)),
        dict(id="L5", statement="TREM2 post/pre >= 5 and DiD CI excludes 0", observed=f"ratio {tr['post_over_pre_ratio']:.1f}; DiD {tr['did_log1p']:.2f} [{tr['did_ci95_lo']:.2f},{tr['did_ci95_hi']:.2f}]", held=bool(L5_hold))]
wcsv(pl.DataFrame(card), "scorecard_expectations.csv")
json.dump(dict(finished=datetime.datetime.utcnow().isoformat() + "Z", seed=SEED, n_boot=NB), open(RES / "stageB_meta.json", "w"))
log("STAGE_B_DONE")
