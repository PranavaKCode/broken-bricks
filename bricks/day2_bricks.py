#!/usr/bin/env python
"""Broken Bricks Day 2 (analysis part): brick scores, reliability, residuals, parallel analysis, discordance, H2b, H2c.

Follows amendment A6 (amendments.md) item by item. Runs on ORCD with the `sc` conda env (numpy, pandas, scipy, scikit-learn,
statsmodels, umap-learn). Input: <ROOT>/results/day2/<label>__<variant>.npz from day2_nam.py, the donor table, donor_tech.csv.
Output: <ROOT>/results/day2/analysis/<variant>/*.csv and summary.json.

usage: day2_bricks.py VARIANT [--quick]     VARIANT in all | nosad
"""
import argparse, json, os, sys, time, warnings
import numpy as np, pandas as pd, scipy.stats as st
from sklearn.linear_model import RidgeCV
from sklearn.covariance import MinCovDet
from sklearn.metrics import roc_auc_score
import statsmodels.api as sm
from statsmodels.stats.anova import anova_lm
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore")
ROOT = "/home/pranavak/orcd/pool/mantis-films/research/seaad"
IN = f"{ROOT}/results/day2"
SEED = 20261011
K = 15
ALPHAS = np.logspace(-4, 4, 33)
SUBS = {"microglia": ("Microglia-PVM", "glia"), "astrocyte": ("Astrocyte", "glia"), "oligodendrocyte": ("Oligodendrocyte", "glia"),
        "opc": ("OPC", "glia"), "sst": ("Sst", "neuron"), "pvalb": ("Pvalb", "neuron"), "vip": ("Vip", "neuron"),
        "lamp5": ("Lamp5", "neuron"), "l23it": ("L2/3 IT", "neuron"), "l4it": ("L4 IT", "neuron"), "l5it": ("L5 IT", "neuron")}
LABELS = list(SUBS)
GLIA = [l for l in SUBS if SUBS[l][1] == "glia"]; NEURON = [l for l in SUBS if SUBS[l][1] == "neuron"]

ap = argparse.ArgumentParser()
ap.add_argument("variant", choices=["all", "nosad"])
ap.add_argument("--quick", action="store_true", help="fewer permutations (debug)")
args = ap.parse_args()
V = args.variant
OUT = f"{IN}/analysis/{V}"; os.makedirs(OUT, exist_ok=True)
NPERM_RHO = 100 if args.quick else 1000
NPERM_PA = 100 if args.quick else 1000
NPERM_DISC = 200 if args.quick else 2000
NSIM_REL = 500 if args.quick else 10000
NPERM_LAB = 1000 if args.quick else 10000

donors = pd.read_csv(f"{ROOT}/results/seaad_mtg_donor_table.csv", index_col="donor_id")
tech = pd.read_csv(f"{IN}/donor_tech.csv", index_col="donor_id")
ALL84 = sorted(donors.index)
cps_all = donors["CPS_Local"].reindex(ALL84).values
sad_all = (donors["severe"].reindex(ALL84).values == 1)
SUMMARY = {"variant": V, "seed": SEED, "k": K}

# ------------------------------------------------------------------ helpers
def loess(x, y, span=0.75, degree=2):
    """local-quadratic LOESS, tricube weights, Gaussian family, no interpolation (R loess defaults)."""
    n = len(x); q = int(np.floor(span * n + 1e-5)); fit = np.empty(n)
    for i in range(n):
        d = np.abs(x - x[i]); idx = np.argsort(d, kind="stable")[:q]; h = d[idx].max()
        if h <= 0:
            fit[i] = y[idx].mean(); continue
        w = (1 - (d[idx] / h) ** 3) ** 3
        Xd = np.vander(x[idx] - x[i], degree + 1, increasing=True)
        sw = np.sqrt(w)
        beta = np.linalg.lstsq(Xd * sw[:, None], y[idx] * sw, rcond=None)[0]
        fit[i] = beta[0]
    return fit

def ridge_oof(X, y, fold, k=K, per_fold_scale=False):
    X = np.asarray(X[:, :k], dtype=float)
    pred = np.zeros(len(y)); g = X[:, 0].std(ddof=1) or 1.0
    for f in np.unique(fold):
        te = fold == f; tr = ~te
        s = (X[tr, 0].std(ddof=1) or 1.0) if per_fold_scale else g
        m = RidgeCV(alphas=ALPHAS).fit(X[tr] / s, y[tr])
        pred[te] = m.predict(X[te] / s)
    return pred

def ridge_oof_fold_pcs(FP, y, fold, k=K):
    """FP: (5, D, kmax) PCA refitted inside each training fold; held-out donors projected"""
    pred = np.zeros(len(y))
    for f in range(5):
        te = fold == f; tr = ~te
        X = FP[f][:, :k]; s = X[tr, 0].std(ddof=1) or 1.0
        m = RidgeCV(alphas=ALPHAS).fit(X[tr] / s, y[tr]); pred[te] = m.predict(X[te] / s)
    return pred

def pr(a, b):
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 3 or np.std(a[ok]) == 0 or np.std(b[ok]) == 0: return np.nan
    return float(np.corrcoef(a[ok], b[ok])[0, 1])

def resid_of(B, cps):
    return B - loess(cps, B)

def hedges_g(a, b):
    n1, n2 = len(a), len(b)
    sp_ = np.sqrt(((n1 - 1) * a.var(ddof=1) + (n2 - 1) * b.var(ddof=1)) / (n1 + n2 - 2))
    d = (a.mean() - b.mean()) / sp_
    J = 1 - 3 / (4 * (n1 + n2 - 2) - 1)
    g = d * J
    se = np.sqrt((n1 + n2) / (n1 * n2) + g ** 2 / (2 * (n1 + n2)))
    return g, se

# ------------------------------------------------------------------ 1. per-subclass brick scores
t0 = time.perf_counter()
cols = {}
rows = []
rng_fold = np.random.default_rng(SEED + 7)
extra_folds = [np.argsort(np.argsort(p)) % 5 for p in [rng_fold.permutation(84) for _ in range(20)]]   # 20 random fold assignments over the 84 donors
ALLpos = {d: i for i, d in enumerate(ALL84)}
supertype_ok = True
for lab in LABELS:
    f = f"{IN}/{lab}__{V}.npz"
    if not os.path.exists(f):
        print("MISSING", f); continue
    z = np.load(f, allow_pickle=True)
    ids = [str(d) for d in z["donors"]]; D = len(ids)
    cps = donors.loc[ids, "CPS_Local"].values.astype(float)
    fold = z["fold"].astype(int)
    P = z["pcs_full"]
    B = ridge_oof(P, cps, fold, K)
    B10 = ridge_oof(P, cps, fold, 10); B20 = ridge_oof(P, cps, fold, 20)
    Bcov = ridge_oof(z["pcs_cov"], cps, fold, K)
    Bfold = ridge_oof_fold_pcs(z["foldpcs"], cps, fold, K)
    rho = st.spearmanr(B, cps).statistic
    rng = np.random.default_rng(SEED + 1000 + LABELS.index(lab))
    null = np.empty(NPERM_RHO)
    for b in range(NPERM_RHO):
        yp = rng.permutation(cps)
        null[b] = st.spearmanr(ridge_oof(P, yp, fold, K), yp).statistic
    rho_p = (1 + (null >= rho).sum()) / (1 + NPERM_RHO)
    # fold-assignment stability
    pos = np.array([ALLpos[d] for d in ids])
    Bs = np.array([ridge_oof(P, cps, ef[pos], K) for ef in extra_folds])
    cm = np.corrcoef(Bs); stab = float((cm.sum() - len(cm)) / (len(cm) * (len(cm) - 1)))
    stab_prim = float(np.mean([pr(B, b) for b in Bs]))
    # split halves (A6 item 4/5)
    HP = z["halfpcs"]; r_half = []; rr_half = []
    for s in range(HP.shape[0]):
        bA = ridge_oof(HP[s, 0], cps, fold, K); bB = ridge_oof(HP[s, 1], cps, fold, K)
        r_half.append(pr(bA, bB))
        rr_half.append(pr(resid_of(bA, cps), resid_of(bB, cps)))
    # supertype-composition triangulation (protocol 4.4 item 2)
    bs = np.full(D, np.nan); rho_s = np.nan; r_ns = np.nan
    try:
        pq = pd.read_parquet(f"{ROOT}/results/day1/{lab}__all_cells.parquet", columns=["donor_id", "Supertype"])
        pq = pq[pq["donor_id"].isin(ids)]
        comp = pd.crosstab(pq["donor_id"], pq["Supertype"]).reindex(ids).fillna(0)
        comp = comp.div(comp.sum(axis=1), axis=0)
        comp = comp.loc[:, (comp > 0).sum(axis=0) >= 5]
        Xs = ((comp - comp.mean()) / comp.std(ddof=1).replace(0, 1)).values
        bs = ridge_oof(Xs, cps, fold, Xs.shape[1])
        rho_s = st.spearmanr(bs, cps).statistic; r_ns = pr(B, bs)
    except Exception as e:
        print("supertype triangulation failed", lab, repr(e)[:120])
    rr_med = float(np.nanmedian(rr_half))
    c = dict(label=lab, ids=ids, cps=cps, fold=fold, B=B, B10=B10, B20=B20, Bcov=Bcov, Bfold=Bfold, Bsup=bs,
             n_nuclei=z["n_nuclei"].astype(int), halves_r=r_half, halves_rr=rr_half)
    cols[lab] = c
    rows.append(dict(variant=V, label=lab, subclass=SUBS[lab][0], klass=SUBS[lab][1], n_donors=D, nsteps=int(z["nsteps"]),
                     split_half_r_median=float(np.nanmedian(r_half)), split_half_r_q25=float(np.nanpercentile(r_half, 25)),
                     split_half_r_q75=float(np.nanpercentile(r_half, 75)), split_half_r_min=float(np.nanmin(r_half)),
                     resid_split_half_r_median=rr_med, resid_split_half_r_q25=float(np.nanpercentile(rr_half, 25)),
                     resid_split_half_r_q75=float(np.nanpercentile(rr_half, 75)),
                     oof_rho=float(rho), oof_rho_perm_p=float(rho_p), oof_pearson=pr(B, cps),
                     passes_reliability=bool(np.nanmedian(r_half) >= 0.6),
                     retained=bool(np.nanmedian(r_half) >= 0.6 and rho > 0.4 and rho_p < 0.01),
                     fold_stability_pairwise_r=stab, fold_stability_r_with_primary=stab_prim,
                     rho_k10=float(st.spearmanr(B10, cps).statistic), rho_k20=float(st.spearmanr(B20, cps).statistic),
                     rho_covNAM=float(st.spearmanr(Bcov, cps).statistic), rho_foldPCA=float(st.spearmanr(Bfold, cps).statistic),
                     r_B_k10_vs_k15=pr(B, B10), r_B_k20_vs_k15=pr(B, B20), r_B_covNAM_vs_primary=pr(B, Bcov), r_B_foldPCA_vs_primary=pr(B, Bfold),
                     supertype_ridge_rho=float(rho_s) if rho_s == rho_s else np.nan, r_B_vs_supertype_ridge=r_ns,
                     varexp_pc1=float(z["sv_full"][0] / z["sv_full"].sum()), varexp_first15=float(z["sv_full"][:15].sum() / z["sv_full"].sum())))
    print(f"[{lab}] D={D} r_half={np.nanmedian(r_half):.3f} rho={rho:.3f} p={rho_p:.4f} rr={rr_med:.3f} stab={stab:.3f} "
          f"rho_sup={rho_s:.3f} r(B,Bsup)={r_ns:.3f}  t={time.perf_counter() - t0:.0f}s", flush=True)

rel = pd.DataFrame(rows)
rel.to_csv(f"{OUT}/brick_reliability.csv", index=False)
retained = [r["label"] for r in rows if r["retained"]]
SUMMARY["retained"] = retained
SUMMARY["n_pass_reliability"] = int(sum(r["passes_reliability"] for r in rows))
SUMMARY["pass_reliability_labels"] = [r["label"] for r in rows if r["passes_reliability"]]
print("retained:", retained, " reliability>=0.6:", SUMMARY["pass_reliability_labels"], flush=True)

# ------------------------------------------------------------------ 2. matrices on the 84 canonical donors
n84 = len(ALL84)
def to84(vec, ids):
    out = np.full(n84, np.nan)
    for d, v in zip(ids, vec): out[ALLpos[d]] = v
    return out
Bm = pd.DataFrame({lab: to84(c["B"], c["ids"]) for lab, c in cols.items()}, index=ALL84)
Rm = pd.DataFrame(index=ALL84, columns=list(cols), dtype=float)
Zm = Rm.copy()
for lab, c in cols.items():
    R = resid_of(c["B"], c["cps"])
    Rm[lab] = to84(R, c["ids"])
    Zm[lab] = to84(R / R.std(ddof=1), c["ids"])
Nn = pd.DataFrame({lab: to84(c["n_nuclei"], c["ids"]) for lab, c in cols.items()}, index=ALL84)
in_scope = np.ones(n84, bool) if V == "all" else ~sad_all           # donors that exist in this variant
use = [l for l in retained]
glia = [l for l in use if SUBS[l][1] == "glia"]; neur = [l for l in use if SUBS[l][1] == "neuron"]
SUMMARY.update({"glia_cols": glia, "neuron_cols": neur})

def composite(Z, cl):
    if not cl: return np.full(len(Z), np.nan)
    M = Z[cl].values
    ok = np.isfinite(M).any(axis=1)
    out = np.full(len(Z), np.nan); out[ok] = np.nanmean(M[ok], axis=1)
    return out
def zs(x):
    m = np.nanmean(x); s = np.nanstd(x, ddof=1)
    return (x - m) / s

def discord_flags(zg, zn, thr=1.5):
    return (np.sign(zg) != np.sign(zn)) & (np.abs(zg) > thr) & (np.abs(zn) > thr)

def disc_rate(Zmat, gcols_idx, ncols_idx):
    G = np.nanmean(Zmat[:, gcols_idx], axis=1); N = np.nanmean(Zmat[:, ncols_idx], axis=1)
    ok = np.isfinite(G) & np.isfinite(N)
    zg = zs(np.where(ok, G, np.nan)); zn = zs(np.where(ok, N, np.nan))
    fl = discord_flags(zg[ok], zn[ok])
    return fl.mean() if ok.sum() else np.nan

# ------------------------------------------------------------------ 3. parallel analysis, discordance
if len(use) >= 2:
    cc = Zm[use].dropna()
    ev = np.sort(np.linalg.eigvalsh(np.corrcoef(cc.values.T)))[::-1]
    rng = np.random.default_rng(SEED + 31)
    nulls = np.empty((NPERM_PA, len(use)))
    X0 = cc.values.copy()
    for b in range(NPERM_PA):
        Xp = np.column_stack([rng.permutation(X0[:, j]) for j in range(X0.shape[1])])
        nulls[b] = np.sort(np.linalg.eigvalsh(np.corrcoef(Xp.T)))[::-1]
    p95 = np.percentile(nulls, 95, axis=0); pmean = nulls.mean(axis=0)
    nfac = 0
    for j in range(len(ev)):
        if ev[j] > p95[j]: nfac += 1
        else: break
    pa = pd.DataFrame({"component": np.arange(1, len(ev) + 1), "eigenvalue": ev, "null_mean": pmean, "null_p95": p95,
                       "above_p95": ev > p95, "perm_p": [(1 + (nulls[:, j] >= ev[j]).sum()) / (1 + NPERM_PA) for j in range(len(ev))]})
    pa.insert(0, "variant", V); pa["n_donors_complete"] = len(cc); pa["columns"] = ",".join(use)
    pa.to_csv(f"{OUT}/h2_parallel_analysis.csv", index=False)
    SUMMARY.update({"pa_n_donors": int(len(cc)), "pa_factors": int(nfac), "pa_eigenvalues": [float(x) for x in ev],
                    "pa_null_p95": [float(x) for x in p95], "pa_pass": bool(nfac >= 2)})
    # PCA scores for floor / technical regression
    Xc = (cc.values - cc.values.mean(0)) / cc.values.std(0, ddof=1)
    U_, S_, Vt_ = np.linalg.svd(Xc, full_matrices=False)
    fac_scores = pd.DataFrame(U_[:, :3] * S_[:3], index=cc.index, columns=["factor1", "factor2", "factor3"])
    fac_loadings = pd.DataFrame(Vt_[:3].T * (S_[:3] / np.sqrt(len(cc) - 1)), index=use, columns=["factor1", "factor2", "factor3"])
    fac_loadings.to_csv(f"{OUT}/h2_factor_loadings.csv")
else:
    SUMMARY.update({"pa_factors": 0, "pa_pass": False, "pa_note": "fewer than 2 retained columns"})
    fac_scores = pd.DataFrame(index=ALL84)

disc = {}
if glia and neur:
    G = composite(Zm, glia); N = composite(Zm, neur)
    ok = np.isfinite(G) & np.isfinite(N) & in_scope
    zG = np.where(ok, zs(np.where(ok, G, np.nan)), np.nan); zN = np.where(ok, zs(np.where(ok, N, np.nan)), np.nan)
    flag = np.zeros(n84, bool); flag[ok] = discord_flags(zG[ok], zN[ok])
    obs_rate = flag[ok].mean(); n_ok = int(ok.sum())
    # permutation null (within-column)
    gi = [use.index(l) for l in glia]; ni = [use.index(l) for l in neur]
    Zu = Zm.loc[in_scope, use].values
    rng = np.random.default_rng(SEED + 41)
    perm_rates = np.empty(NPERM_DISC)
    for b in range(NPERM_DISC):
        Zp = Zu.copy()
        for j in range(Zp.shape[1]):
            m = np.isfinite(Zp[:, j]); Zp[m, j] = rng.permutation(Zp[m, j])
        perm_rates[b] = disc_rate(Zp, gi, ni)
    # reliability-matched null
    rr = np.array([max(0.0, float(rel.loc[rel.label == l, "resid_split_half_r_median"].iloc[0])) for l in use])
    nD = int(in_scope.sum()); sims = np.empty(NSIM_REL)
    rng = np.random.default_rng(SEED + 43)
    for b in range(NSIM_REL):
        F = rng.standard_normal(nD); E = rng.standard_normal((nD, len(use)))
        Zs = np.sqrt(rr)[None, :] * F[:, None] + np.sqrt(1 - rr)[None, :] * E
        sims[b] = disc_rate(Zs, gi, ni)
    disc = {"n_donors_with_both": n_ok, "n_discordant": int(flag.sum()), "observed_rate": float(obs_rate),
            "perm_null_mean": float(np.mean(perm_rates)), "perm_null_p95": float(np.percentile(perm_rates, 95)),
            "perm_p": float((1 + (perm_rates >= obs_rate).sum()) / (1 + NPERM_DISC)),
            "relmatched_null_mean": float(np.mean(sims)), "relmatched_null_p95": float(np.percentile(sims, 95)),
            "relmatched_p": float((1 + (sims >= obs_rate).sum()) / (1 + NSIM_REL)),
            "resid_reliability_used": {l: float(r) for l, r in zip(use, rr)},
            "corr_zG_zN": pr(zG, zN),
            "n_glia_to_sick_neuron_healthy": int((flag & (zG > 0)).sum()), "n_neuron_sick_glia_healthy": int((flag & (zN > 0)).sum()),
            "criterion_B_pass": bool(obs_rate >= 0.10 and float(np.mean(sims)) <= 0.05)}
    SUMMARY["discordance"] = disc
    pd.DataFrame([{**{"variant": V}, **{k: v for k, v in disc.items() if not isinstance(v, dict)}}]).to_csv(f"{OUT}/h2_discordance.csv", index=False)
else:
    G = N = zG = zN = np.full(n84, np.nan); flag = np.zeros(n84, bool); ok = np.zeros(n84, bool)
    SUMMARY["discordance"] = {"note": "no retained glial or neuronal column"}
delta = zG - zN
SUMMARY["H2_reliability_criterion_pass"] = bool(SUMMARY["n_pass_reliability"] >= 4)
critB = bool(disc.get("criterion_B_pass", False)) if disc else False
SUMMARY["H2_structure_criterion_A_pass"] = bool(SUMMARY.get("pa_pass", False))
SUMMARY["H2_structure_criterion_B_pass"] = critB
SUMMARY["H2_pass"] = bool(SUMMARY["H2_reliability_criterion_pass"] and (SUMMARY["H2_structure_criterion_A_pass"] or critB))
print("H2:", {k: SUMMARY[k] for k in ["n_pass_reliability", "pa_factors", "H2_reliability_criterion_pass", "H2_structure_criterion_A_pass", "H2_structure_criterion_B_pass", "H2_pass"]}, flush=True)
print("discordance:", {k: v for k, v in SUMMARY["discordance"].items() if not isinstance(v, dict)}, flush=True)

# ------------------------------------------------------------------ 4. technical covariates (A6 item 9)
T = pd.DataFrame({"log10_n_nuclei": np.log10(tech["n_nuclei"]), "rin": tech["rin_mean"], "pmi": donors["pmi"],
                  "neun_frac": tech["neun_frac"], "multiome_frac": tech["frac_multiome_all"], "seq_saturation": tech["seq_saturation"]}).reindex(ALL84)
Tz = (T - T.mean()) / T.std(ddof=1)
axes = {f"Z_{l}": Zm[l].values for l in use}
axes.update({"zG_glia": zG, "zN_neuron": zN, "delta_zG_minus_zN": delta})
for c_ in fac_scores.columns: axes[f"PAfactor_{c_}"] = fac_scores[c_].reindex(ALL84).values
# also the unfiltered columns for transparency
for l in cols:
    if l not in use: axes[f"Z_{l}(not retained)"] = Zm[l].values
trows = []
for name, y in axes.items():
    m = np.isfinite(y) & in_scope
    if m.sum() < 20 or not np.isfinite(y[m]).any(): continue
    X = sm.add_constant(Tz.values[m]); fit = sm.OLS(y[m], X).fit()
    sp_ = [st.spearmanr(Tz.values[m][:, j], y[m]) for j in range(Tz.shape[1])]
    minp = min(s.pvalue for s in sp_); jmin = int(np.argmin([s.pvalue for s in sp_]))
    tracks = bool(fit.f_pvalue < 0.05 or minp < 0.05 / 6)
    trows.append({"variant": V, "axis": name, "n": int(m.sum()), "adj_r2": float(fit.rsquared_adj), "joint_F": float(fit.fvalue), "joint_F_p": float(fit.f_pvalue),
                  "strongest_covariate": Tz.columns[jmin], "strongest_spearman": float(sp_[jmin].statistic), "strongest_spearman_p": float(minp),
                  "tracks_technical": tracks, **{f"beta_{c}": float(b) for c, b in zip(Tz.columns, fit.params[1:])}})
tech_df = pd.DataFrame(trows); tech_df.to_csv(f"{OUT}/h2_technical_regression.csv", index=False)
SUMMARY["technical_axes_flagged"] = tech_df.loc[tech_df.tracks_technical, "axis"].tolist()

# ------------------------------------------------------------------ 5. SAD recovery (A6 item 8)
def recovery(M, name):
    M = pd.DataFrame(M)
    Zs_ = (M - M.mean()) / M.std(ddof=1)
    Zs_ = Zs_.fillna(Zs_.median())
    mask = in_scope
    if V == "nosad": return None
    mcd = MinCovDet(support_fraction=0.75, random_state=0).fit(Zs_.values[mask])
    dist = mcd.mahalanobis(Zs_.values[mask]); lab_ = sad_all[mask]
    auc = roc_auc_score(lab_, dist)
    top = np.argsort(-dist)[:int(lab_.sum())]; prec = lab_[top].mean()
    rng = np.random.default_rng(SEED + 51); null = np.empty(NPERM_LAB)
    for b in range(NPERM_LAB):
        null[b] = roc_auc_score(rng.permutation(lab_), dist)
    p = (1 + (null >= auc).sum()) / (1 + NPERM_LAB)
    meanz = Zs_.values[mask].mean(axis=1); auc_mean = roc_auc_score(lab_, meanz)
    return {"variant": V, "matrix": name, "n_columns": M.shape[1], "n_donors": int(mask.sum()), "n_sad": int(lab_.sum()),
            "auc_robust_distance": float(auc), "precision_at_11": float(prec), "perm_p": float(p), "auc_mean_z": float(auc_mean),
            "n_sad_in_top11": int(lab_[top].sum())}
srows = []
if V == "all":
    for name, M in [("B_retained", Bm[use]), ("Z_retained", Zm[use]), ("B_all11", Bm)]:
        if M.shape[1] >= 2:
            r_ = recovery(M, name)
            if r_: srows.append(r_)
if srows:
    pd.DataFrame(srows).to_csv(f"{OUT}/h2_sad_recovery.csv", index=False)
    SUMMARY["sad_recovery"] = srows
    SUMMARY["sad_recovery_pass"] = bool(srows[0]["auc_robust_distance"] >= 0.75)

# ------------------------------------------------------------------ 6. H2b (A6 item 10)
def ordmap(s, mp):
    return s.map(lambda v: mp.get(v, np.nan) if isinstance(v, str) else np.nan)
dd = donors.reindex(ALL84)
lewy = ordmap(dd["Highest Lewy Body Disease"], {"Not Identified (olfactory bulb not assessed)": 0, "Not Identified (olfactory bulb assessed)": 0,
                                                "Olfactory bulb only": 1, "Amygdala-predominant": 1, "Brainstem-predominant": 1,
                                                "Limbic (Transitional)": 2, "Neocortical (Diffuse)": 3})
late = ordmap(dd["LATE"], {"Not Identified": 0, "LATE Stage 1": 1, "LATE Stage 2": 2, "LATE Stage 3": 3})
arter = ordmap(dd["Arteriolosclerosis"], {"None": 0, "Mild": 1, "Moderate": 2, "Severe": 3})
caa = ordmap(dd["Overall CAA Score"], {"Not identified": 0, "Mild": 1, "Moderate": 2, "Severe": 3})
micro = np.log1p(dd["Total Microinfarcts (not observed grossly)"].astype(float))
abdiff = dd["CPS_Local_ABeta"] - dd["CPS_Local_pTau"]
preds = pd.DataFrame({"LATE": late, "Lewy": lewy, "arteriolosclerosis": arter, "microinfarcts_log1p": micro, "CAA": caa, "ABeta_minus_pTau": abdiff}, index=ALL84)
preds["ABeta_minus_pTau"] = (preds["ABeta_minus_pTau"] - preds["ABeta_minus_pTau"].mean()) / preds["ABeta_minus_pTau"].std(ddof=1)
cps_s = pd.Series(cps_all, index=ALL84, name="CPS_Local")
def h2b_model(y, name):
    d = pd.concat([pd.Series(y, index=ALL84, name="y"), cps_s, preds], axis=1)[in_scope].dropna()
    if len(d) < 30 or d["y"].std() == 0: return None
    Xp = (d[preds.columns] - d[preds.columns].mean()) / d[preds.columns].std(ddof=1)
    X1 = sm.add_constant(pd.concat([d[["CPS_Local"]], Xp], axis=1)); X0 = sm.add_constant(d[["CPS_Local"]])
    m1 = sm.OLS(d["y"], X1).fit(); m0 = sm.OLS(d["y"], X0).fit()
    an = anova_lm(m0, m1); F = float(an["F"].iloc[1]); pF = float(an["Pr(>F)"].iloc[1])
    pv = m1.pvalues[preds.columns]; holm = multipletests(pv.values, method="holm")[1]
    out = []
    for i, c in enumerate(preds.columns):
        ci = m1.conf_int().loc[c]
        out.append({"variant": V, "outcome": name, "predictor": c, "n": int(len(d)), "coef_per_sd": float(m1.params[c]), "ci_low": float(ci[0]), "ci_high": float(ci[1]),
                    "p": float(pv[c]), "p_holm": float(holm[i]), "omnibus_F": F, "omnibus_p": pF, "adj_r2_model": float(m1.rsquared_adj),
                    "adj_r2_cps_only": float(m0.rsquared_adj), "delta_adj_r2": float(m1.rsquared_adj - m0.rsquared_adj)})
    return out
h2b = []
outcomes = {"delta_zG_minus_zN": delta, "abs_delta": np.abs(delta), "zG_glia": zG, "zN_neuron": zN}
for l in use: outcomes[f"Z_{l}"] = Zm[l].values
for name, y in outcomes.items():
    r_ = h2b_model(y, name)
    if r_: h2b += r_
h2b_df = pd.DataFrame(h2b); h2b_df.to_csv(f"{OUT}/h2b_copathology_regressions.csv", index=False)
if len(h2b_df):
    pri = h2b_df[h2b_df.outcome == "delta_zG_minus_zN"]
    if len(pri):
        ab = pri[pri.predictor == "ABeta_minus_pTau"].iloc[0]
        SUMMARY["H2b"] = {"omnibus_F": float(pri.omnibus_F.iloc[0]), "omnibus_p": float(pri.omnibus_p.iloc[0]),
                          "delta_adj_r2": float(pri.delta_adj_r2.iloc[0]), "n": int(pri.n.iloc[0]),
                          "abeta_minus_ptau_coef": float(ab.coef_per_sd), "abeta_minus_ptau_p": float(ab.p),
                          "pass": bool(pri.omnibus_p.iloc[0] < 0.05 or ab.p < 0.05),
                          "predictors": {r.predictor: {"coef": float(r.coef_per_sd), "p": float(r.p), "p_holm": float(r.p_holm)} for r in pri.itertuples()}}
# discordant donors vs co-pathologies (descriptive Fisher / Mann-Whitney)
if flag.any():
    drows = []
    for c in preds.columns:
        x = preds[c].values; m = np.isfinite(x) & ok
        if m.sum() > 10 and flag[m].any():
            u = st.mannwhitneyu(x[m & flag], x[m & ~flag]); drows.append({"variant": V, "predictor": c, "mean_discordant": float(np.mean(x[m & flag])), "mean_other": float(np.mean(x[m & ~flag])), "mannwhitney_p": float(u.pvalue)})
    pd.DataFrame(drows).to_csv(f"{OUT}/h2b_discordant_vs_other.csv", index=False)
# secondary: APOE4, dementia, CASI, MMSE with CPS as covariate
sec = []
second = {"apoe4_carrier": dd["apoe4_carrier"].astype(float), "dementia": dd["dementia"].astype(float), "last_CASI": dd["Last CASI Score"].astype(float),
          "last_MMSE": dd["Last MMSE Score"].astype(float)}
for name, y in {"delta_zG_minus_zN": delta, "zG_glia": zG, "zN_neuron": zN, **{f"Z_{l}": Zm[l].values for l in use}}.items():
    for pn, pv_ in second.items():
        d = pd.concat([pd.Series(y, index=ALL84, name="y"), cps_s, pv_.rename("x")], axis=1)[in_scope].dropna()
        if len(d) < 30 or d["y"].std() == 0: continue
        m_ = sm.OLS(d["y"], sm.add_constant(d[["x", "CPS_Local"]])).fit()
        ci = m_.conf_int().loc["x"]
        sec.append({"variant": V, "outcome": name, "predictor": pn, "n": int(len(d)), "coef": float(m_.params["x"]), "ci_low": float(ci[0]), "ci_high": float(ci[1]), "p": float(m_.pvalues["x"])})
pd.DataFrame(sec).to_csv(f"{OUT}/h2b_secondary_cognition_apoe.csv", index=False)

# ------------------------------------------------------------------ 7. H2c (A6 item 11)
hi = (dd["Overall AD neuropathological Change"] == "High").values & in_scope
res_ = (dd["Cognitive Status"] == "No dementia").values; dem_ = (dd["Cognitive Status"] == "Dementia").values
c_res = hi & res_; c_dem = hi & dem_
def h2c(vec, name, nboot=NPERM_LAB):
    a = vec[c_res]; b = vec[c_dem]; a = a[np.isfinite(a)]; b = b[np.isfinite(b)]
    if len(a) < 3 or len(b) < 3: return None
    g, se = hedges_g(a, b)
    rng = np.random.default_rng(SEED + 61)
    bs = np.array([hedges_g(rng.choice(a, len(a)), rng.choice(b, len(b)))[0] for _ in range(nboot)])
    t = st.ttest_ind(a, b, equal_var=False, alternative="less")
    diff = a.mean() - b.mean(); allv = np.concatenate([a, b]); nullm = np.empty(nboot)
    for i in range(nboot):
        pm = rng.permutation(allv); nullm[i] = pm[:len(a)].mean() - pm[len(a):].mean()
    # CPS-adjusted
    return {"variant": V, "outcome": name, "n_resilient": len(a), "n_demented": len(b), "mean_resilient": float(a.mean()), "mean_demented": float(b.mean()),
            "hedges_g_resilient_minus_demented": float(g), "ci95_low_HO": float(g - 1.96 * se), "ci95_high_HO": float(g + 1.96 * se),
            "ci95_low_boot": float(np.percentile(bs, 2.5)), "ci95_high_boot": float(np.percentile(bs, 97.5)),
            "welch_t": float(t.statistic), "welch_p_one_sided_less": float(t.pvalue), "perm_p_one_sided_less": float((1 + (nullm <= diff).sum()) / (1 + nboot))}
h2c_rows = []
neuron_comp_raw = composite(Zm, neur) if neur else np.full(n84, np.nan)
glia_comp_raw = composite(Zm, glia) if glia else np.full(n84, np.nan)
all_neur = [l for l in NEURON if l in cols]; all_glia = [l for l in GLIA if l in cols]
neur_all_raw = composite(Zm, all_neur)
for name, vec in [("neuronal_composite_retained(C5)", neuron_comp_raw), ("neuronal_composite_all_columns", neur_all_raw),
                  ("glia_composite_retained", glia_comp_raw), ("glia_composite_all_columns", composite(Zm, all_glia))] + \
                 [(f"Z_{l}", Zm[l].values) for l in cols]:
    r_ = h2c(np.asarray(vec, float), name)
    if r_: h2c_rows.append(r_)
# CPS-adjusted regression of the primary outcome
for name, vec in [("neuronal_composite_retained(C5)", neuron_comp_raw), ("neuronal_composite_all_columns", neur_all_raw)]:
    d = pd.DataFrame({"y": vec, "resilient": res_.astype(float), "CPS": cps_all}, index=ALL84)[c_res | c_dem].dropna()
    if len(d) > 10 and d["y"].std() > 0:
        m_ = sm.OLS((d["y"] - d["y"].mean()) / d["y"].std(ddof=1), sm.add_constant(d[["resilient", "CPS"]])).fit()
        ci = m_.conf_int().loc["resilient"]
        h2c_rows.append({"variant": V, "outcome": name + "|CPS-adjusted OLS (standardised outcome)", "n_resilient": int(d.resilient.sum()), "n_demented": int((1 - d.resilient).sum()),
                         "hedges_g_resilient_minus_demented": float(m_.params["resilient"]), "ci95_low_HO": float(ci[0]), "ci95_high_HO": float(ci[1]),
                         "welch_p_one_sided_less": float(m_.pvalues["resilient"] / 2 if m_.params["resilient"] < 0 else 1 - m_.pvalues["resilient"] / 2)})
h2c_df = pd.DataFrame(h2c_rows); h2c_df.to_csv(f"{OUT}/h2c_resilience.csv", index=False)
if len(h2c_df):
    p0 = h2c_df[h2c_df.outcome == "neuronal_composite_retained(C5)"]
    SUMMARY["H2c"] = (p0.iloc[0].to_dict() if len(p0) else {"note": "no retained neuronal column"})
    SUMMARY["H2c_all_columns"] = h2c_df[h2c_df.outcome == "neuronal_composite_all_columns"].iloc[0].to_dict() if (h2c_df.outcome == "neuronal_composite_all_columns").any() else None
SUMMARY["h2c_group_sizes"] = {"high_adnc": int(hi.sum()), "resilient": int(c_res.sum()), "demented": int(c_dem.sum())}

# ------------------------------------------------------------------ 8. donor_bricks table and Mantis donors floor
tab = pd.DataFrame(index=pd.Index(ALL84, name="donor_id"))
tab["in_variant"] = in_scope.astype(int)
tab["CPS_Local"] = cps_all; tab["CPS_Local_ABeta"] = dd["CPS_Local_ABeta"].values; tab["CPS_Local_pTau"] = dd["CPS_Local_pTau"].values
tab["severe_donor"] = sad_all.astype(int); tab["cognitive_status"] = dd["Cognitive Status"].values; tab["adnc"] = dd["Overall AD neuropathological Change"].values
tab["resilient_high_adnc"] = (hi & res_).astype(int); tab["demented_high_adnc"] = (hi & dem_).astype(int)
tab["apoe_genotype"] = dd["APOE Genotype"].values; tab["apoe4_carrier"] = dd["apoe4_carrier"].values
tab["age"] = dd["age"].values; tab["sex"] = dd["Sex"].values; tab["pmi"] = dd["pmi"].values
for c in preds.columns: tab[c] = preds[c].values
for l in LABELS:
    if l in cols:
        tab[f"n_nuclei_{l}"] = Nn[l].values; tab[f"B_{l}"] = Bm[l].values; tab[f"R_{l}"] = Rm[l].values; tab[f"Z_{l}"] = Zm[l].values
tab["glia_composite"] = G; tab["neuron_composite"] = N; tab["zG"] = zG; tab["zN"] = zN; tab["delta_zG_minus_zN"] = delta
tab["discordant"] = flag.astype(int)
tab["discordance_label"] = np.where(flag & (zG > 0), "glia sicker / neurons healthier than CPS predicts",
                             np.where(flag & (zN > 0), "neurons sicker / glia healthier than CPS predicts",
                             np.where(ok, "not discordant", "not scored")))
for c_ in fac_scores.columns: tab[c_] = fac_scores[c_].reindex(ALL84).values
tab.to_csv(f"{OUT}/donor_bricks.csv")
for name, key in [("B_k10", "B10"), ("B_k20", "B20"), ("B_covNAM", "Bcov"), ("B_foldPCA", "Bfold"), ("B_supertype_ridge", "Bsup")]:
    pd.DataFrame({lab: to84(c[key], c["ids"]) for lab, c in cols.items()}, index=pd.Index(ALL84, name="donor_id")).to_csv(f"{OUT}/brick_scores_{name}.csv")

# Mantis donors floor
try:
    import umap
    Zf = Zm[use].fillna(0).values if len(use) >= 2 else Zm.fillna(0).values
    Zf = Zf[in_scope] if V == "nosad" else Zf
    Xc = Zf - Zf.mean(0); U_, S_, Vt_ = np.linalg.svd(Xc, full_matrices=False)
    pc = U_[:, :2] * S_[:2]
    um = umap.UMAP(n_neighbors=10, min_dist=0.3, random_state=SEED, n_jobs=1).fit_transform(Xc)
    Bv = Bm.copy(); Bv = (Bv - Bv.mean()) / Bv.std(ddof=1); Bv = Bv.fillna(Bv.median())
    fl = pd.DataFrame(index=pd.Index(ALL84, name="donor_id"))
    sel = in_scope if V == "nosad" else np.ones(n84, bool)
    fl = tab[sel].copy()
    fl.insert(0, "floor_number", 2); fl.insert(1, "title", fl.index)
    fl.insert(2, "coordinate1", pc[:, 0]); fl.insert(3, "coordinate2", pc[:, 1])
    fl.insert(4, "umap1", um[:, 0]); fl.insert(5, "umap2", um[:, 1])
    for j, l in enumerate(LABELS):
        if l in Bv: fl[f"v{j}"] = Bv[l].values[sel]
    fl.to_csv(f"{OUT}/donors_floor.csv")
    SUMMARY["floor_pc_var_explained"] = [float(x) for x in (S_[:2] ** 2 / (S_ ** 2).sum())]
except Exception as e:
    print("floor export failed:", repr(e)[:200])

SUMMARY["seconds"] = round(time.perf_counter() - t0, 1)
json.dump(SUMMARY, open(f"{OUT}/summary.json", "w"), indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
print("DONE", V, SUMMARY["seconds"], "s", flush=True)
