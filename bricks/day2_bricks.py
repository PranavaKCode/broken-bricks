#!/usr/bin/env python
"""Broken Bricks Day 2 (analysis part): brick scores, reliability, residuals, parallel analysis, discordance, H2b, H2c.

Follows amendments A6 and A8 (amendments.md). Runs on ORCD with the `sc` conda env (numpy, pandas, scipy, scikit-learn, statsmodels).
Input : <ROOT>/results/day2/<label>__<variant>.npz from day2_nam.py, the donor table, donor_tech.csv.
Output: <ROOT>/results/day2/analysis/<variant>/*.csv and summary.json

Tiers (A8, A9): prereg    = preregistered retained columns (split-half r >= 0.6, rho > 0.4, perm P < 0.01), primary NAM   (verdicts)
            allcols   = all columns with split-half r >= 0.6 (no rho filter), primary NAM                           (sensitivity)
            techNAM   = NAM residualised for age, sex, PMI, method + 4 technical covariates, all columns            (sensitivity)
            techNAM_filtered = techNAM columns passing the same r / rho / P filter                                  (sensitivity)
            techResid = allcols residual columns regressed on the six technical covariates, re-standardised       (sensitivity, A9)

usage: day2_bricks.py VARIANT [--quick]     VARIANT in all | nosad
"""
import argparse, json, os, time, warnings
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
OUT = f"{IN}/analysis/{V}" + ("_quick" if args.quick else ""); os.makedirs(OUT, exist_ok=True)
NPERM_RHO = 100 if args.quick else 1000
NPERM_PA = 100 if args.quick else 1000
NPERM_DISC = 200 if args.quick else 2000
NSIM_REL = 500 if args.quick else 10000
NPERM_LAB = 1000 if args.quick else 10000

donors = pd.read_csv(f"{ROOT}/results/seaad_mtg_donor_table.csv", index_col="donor_id")
tech = pd.read_csv(f"{IN}/donor_tech.csv", index_col="donor_id")
ALL84 = sorted(donors.index); n84 = len(ALL84)
ALLpos = {d: i for i, d in enumerate(ALL84)}
cps_all = donors["CPS_Local"].reindex(ALL84).values
sad_all = (donors["severe"].reindex(ALL84).values == 1)
in_scope = np.ones(n84, bool) if V == "all" else ~sad_all
dd = donors.reindex(ALL84)
cps_s = pd.Series(cps_all, index=ALL84, name="CPS_Local")
SUMMARY = {"variant": V, "seed": SEED, "k": K, "tiers": {}}

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
        fit[i] = np.linalg.lstsq(Xd * sw[:, None], y[idx] * sw, rcond=None)[0][0]
    return fit

def ridge_oof(X, y, fold, k=K):
    X = np.asarray(X[:, :k], dtype=float)
    pred = np.zeros(len(y)); g = X[:, 0].std(ddof=1) or 1.0
    for f in np.unique(fold):
        te = fold == f; tr = ~te
        m = RidgeCV(alphas=ALPHAS).fit(X[tr] / g, y[tr])
        pred[te] = m.predict(X[te] / g)
    return pred

def ridge_oof_fold_pcs(FP, y, fold, k=K):
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

def resid_of(B, cps): return B - loess(cps, B)

def hedges_g(a, b):
    n1, n2 = len(a), len(b)
    sp_ = np.sqrt(((n1 - 1) * a.var(ddof=1) + (n2 - 1) * b.var(ddof=1)) / (n1 + n2 - 2))
    g = (a.mean() - b.mean()) / sp_ * (1 - 3 / (4 * (n1 + n2 - 2) - 1))
    return g, np.sqrt((n1 + n2) / (n1 * n2) + g ** 2 / (2 * (n1 + n2)))

def to84(vec, ids):
    out = np.full(n84, np.nan)
    for d, v in zip(ids, vec): out[ALLpos[d]] = v
    return out

def zs(x):
    return (x - np.nanmean(x)) / np.nanstd(x, ddof=1)

# ------------------------------------------------------------------ 1. per-subclass brick scores
t0 = time.perf_counter()
cols, rows = {}, []
rng_fold = np.random.default_rng(SEED + 7)
extra_folds = [np.argsort(np.argsort(p)) % 5 for p in [rng_fold.permutation(84) for _ in range(20)]]   # 20 random fold assignments

def perm_rho(P, cps, fold, seedoff):
    rng = np.random.default_rng(SEED + 1000 + seedoff)
    null = np.empty(NPERM_RHO)
    for b in range(NPERM_RHO):
        yp = rng.permutation(cps)
        null[b] = st.spearmanr(ridge_oof(P, yp, fold, K), yp).statistic
    return null

for li, lab in enumerate(LABELS):
    f = f"{IN}/{lab}__{V}.npz"
    if not os.path.exists(f):
        print("MISSING", f); continue
    z = np.load(f, allow_pickle=True)
    ids = [str(d) for d in z["donors"]]; D = len(ids)
    cps = donors.loc[ids, "CPS_Local"].values.astype(float)
    fold = z["fold"].astype(int)
    P = z["pcs_full"]; Pt = z["pcs_tech"]
    B = ridge_oof(P, cps, fold, K); Bt = ridge_oof(Pt, cps, fold, K)
    B10 = ridge_oof(P, cps, fold, 10); B20 = ridge_oof(P, cps, fold, 20)
    Bcov = ridge_oof(z["pcs_cov"], cps, fold, K); Bfold = ridge_oof_fold_pcs(z["foldpcs"], cps, fold, K)
    rho = st.spearmanr(B, cps).statistic; rho_t = st.spearmanr(Bt, cps).statistic
    null = perm_rho(P, cps, fold, li); rho_p = (1 + (null >= rho).sum()) / (1 + NPERM_RHO)
    null_t = perm_rho(Pt, cps, fold, 100 + li); rho_p_t = (1 + (null_t >= rho_t).sum()) / (1 + NPERM_RHO)
    pos = np.array([ALLpos[d] for d in ids])
    Bs = np.array([ridge_oof(P, cps, ef[pos], K) for ef in extra_folds])
    cm = np.corrcoef(Bs); stab = float((cm.sum() - len(cm)) / (len(cm) * (len(cm) - 1))); stab_prim = float(np.mean([pr(B, b) for b in Bs]))
    out_h = {}
    for nm, HPk, src in [("", "halfpcs", B), ("_tech", "halfpcs_tech", Bt)]:
        HP = z[HPk]; rh, rrh = [], []
        for s in range(HP.shape[0]):
            bA = ridge_oof(HP[s, 0], cps, fold, K); bB = ridge_oof(HP[s, 1], cps, fold, K)
            rh.append(pr(bA, bB)); rrh.append(pr(resid_of(bA, cps), resid_of(bB, cps)))
        out_h[nm] = (rh, rrh)
    rh, rrh = out_h[""]; rht, rrht = out_h["_tech"]
    # supertype-composition triangulation (protocol 4.4 item 2)
    bs = np.full(D, np.nan); rho_s = np.nan; r_ns = np.nan
    try:
        pq = pd.read_parquet(f"{ROOT}/results/day1/{lab}__all_cells.parquet", columns=["donor_id", "Supertype"])
        pq = pq[pq["donor_id"].isin(ids)]
        comp = pd.crosstab(pq["donor_id"], pq["Supertype"]).reindex(ids).fillna(0)
        comp = comp.div(comp.sum(axis=1), axis=0); comp = comp.loc[:, (comp > 0).sum(axis=0) >= 5]
        Xs = ((comp - comp.mean()) / comp.std(ddof=1).replace(0, 1)).values
        bs = ridge_oof(Xs, cps, fold, Xs.shape[1]); rho_s = st.spearmanr(bs, cps).statistic; r_ns = pr(B, bs)
    except Exception as e:
        print("supertype triangulation failed", lab, repr(e)[:120])
    cols[lab] = dict(ids=ids, cps=cps, fold=fold, B=B, Bt=Bt, B10=B10, B20=B20, Bcov=Bcov, Bfold=Bfold, Bsup=bs,
                     n_nuclei=z["n_nuclei"].astype(int), rr=float(np.nanmedian(rrh)), rr_t=float(np.nanmedian(rrht)))
    rows.append(dict(variant=V, label=lab, subclass=SUBS[lab][0], klass=SUBS[lab][1], n_donors=D, nsteps=int(z["nsteps"]),
                     split_half_r_median=float(np.nanmedian(rh)), split_half_r_q25=float(np.nanpercentile(rh, 25)), split_half_r_q75=float(np.nanpercentile(rh, 75)),
                     split_half_r_min=float(np.nanmin(rh)), resid_split_half_r_median=float(np.nanmedian(rrh)),
                     oof_rho=float(rho), oof_rho_perm_p=float(rho_p), oof_pearson=pr(B, cps),
                     passes_reliability=bool(np.nanmedian(rh) >= 0.6),
                     retained=bool(np.nanmedian(rh) >= 0.6 and rho > 0.4 and rho_p < 0.01),
                     fold_stability_pairwise_r=stab, fold_stability_r_with_primary=stab_prim,
                     rho_k10=float(st.spearmanr(B10, cps).statistic), rho_k20=float(st.spearmanr(B20, cps).statistic),
                     rho_covNAM=float(st.spearmanr(Bcov, cps).statistic), rho_foldPCA=float(st.spearmanr(Bfold, cps).statistic),
                     r_B_k10_vs_k15=pr(B, B10), r_B_k20_vs_k15=pr(B, B20), r_B_covNAM_vs_primary=pr(B, Bcov), r_B_foldPCA_vs_primary=pr(B, Bfold),
                     supertype_ridge_rho=float(rho_s) if rho_s == rho_s else np.nan, r_B_vs_supertype_ridge=r_ns,
                     varexp_pc1=float(z["sv_full"][0] / z["sv_full"].sum()), varexp_first15=float(z["sv_full"][:15].sum() / z["sv_full"].sum()),
                     tech_split_half_r_median=float(np.nanmedian(rht)), tech_resid_split_half_r_median=float(np.nanmedian(rrht)),
                     tech_oof_rho=float(rho_t), tech_oof_rho_perm_p=float(rho_p_t), r_B_tech_vs_primary=pr(B, Bt),
                     tech_passes_reliability=bool(np.nanmedian(rht) >= 0.6),
                     tech_retained=bool(np.nanmedian(rht) >= 0.6 and rho_t > 0.4 and rho_p_t < 0.01)))
    print(f"[{lab}] D={D} r_half={np.nanmedian(rh):.3f} rho={rho:.3f} p={rho_p:.4f} rr={np.nanmedian(rrh):.3f} | tech: r={np.nanmedian(rht):.3f} "
          f"rho={rho_t:.3f} p={rho_p_t:.4f} rr={np.nanmedian(rrht):.3f} | stab={stab:.3f} rho_sup={rho_s:.3f}  t={time.perf_counter() - t0:.0f}s", flush=True)

rel = pd.DataFrame(rows); rel.to_csv(f"{OUT}/brick_reliability.csv", index=False)
SUMMARY["n_pass_reliability"] = int(rel["passes_reliability"].sum())
SUMMARY["reliability_median_range"] = [float(rel.split_half_r_median.min()), float(rel.split_half_r_median.max())]
TIERS = {"prereg": ("B", [r.label for r in rel.itertuples() if r.retained]),
         "allcols": ("B", [r.label for r in rel.itertuples() if r.passes_reliability]),
         "techNAM": ("Bt", [r.label for r in rel.itertuples() if r.tech_passes_reliability]),
         "techNAM_filtered": ("Bt", [r.label for r in rel.itertuples() if r.tech_retained]),
         "techResid": ("B", [r.label for r in rel.itertuples() if r.passes_reliability])}
print("tiers:", TIERS, flush=True)

# ------------------------------------------------------------------ 2. co-pathology predictors (for H2b)
def ordmap(s, mp): return s.map(lambda v: mp.get(v, np.nan) if isinstance(v, str) else np.nan)
lewy = ordmap(dd["Highest Lewy Body Disease"], {"Not Identified (olfactory bulb not assessed)": 0, "Not Identified (olfactory bulb assessed)": 0, "Olfactory bulb only": 1,
                                                "Amygdala-predominant": 1, "Brainstem-predominant": 1, "Limbic (Transitional)": 2, "Neocortical (Diffuse)": 3})
late = ordmap(dd["LATE"], {"Not Identified": 0, "LATE Stage 1": 1, "LATE Stage 2": 2, "LATE Stage 3": 3})
arter = ordmap(dd["Arteriolosclerosis"], {"None": 0, "Mild": 1, "Moderate": 2, "Severe": 3})
caa = ordmap(dd["Overall CAA Score"], {"Not identified": 0, "Mild": 1, "Moderate": 2, "Severe": 3})
abdiff = dd["CPS_Local_ABeta"] - dd["CPS_Local_pTau"]
preds = pd.DataFrame({"LATE": late, "Lewy": lewy, "arteriolosclerosis": arter, "microinfarcts_log1p": np.log1p(dd["Total Microinfarcts (not observed grossly)"].astype(float)),
                      "CAA": caa, "ABeta_minus_pTau": (abdiff - abdiff.mean()) / abdiff.std(ddof=1)}, index=ALL84)
hi = (dd["Overall AD neuropathological Change"] == "High").values & in_scope
res_ = (dd["Cognitive Status"] == "No dementia").values; dem_ = (dd["Cognitive Status"] == "Dementia").values
c_res = hi & res_; c_dem = hi & dem_
SUMMARY["h2c_group_sizes"] = {"high_adnc": int(hi.sum()), "resilient": int(c_res.sum()), "demented": int(c_dem.sum())}
T = pd.DataFrame({"log10_n_nuclei": np.log10(tech["n_nuclei"]), "rin": tech["rin_mean"], "pmi": donors["pmi"], "neun_frac": tech["neun_frac"],
                  "multiome_frac": tech["frac_multiome_all"], "seq_saturation": tech["seq_saturation"]}).reindex(ALL84)
Tz = (T - T.mean()) / T.std(ddof=1)
SUMMARY["technical_covariate_correlations_with_CPS"] = {c: float(st.spearmanr(Tz[c].values, cps_all).statistic) for c in Tz.columns}

def composite(Zdf, cl):
    if not cl: return np.full(len(Zdf), np.nan)
    M = Zdf[cl].values; okr = np.isfinite(M).any(axis=1)
    o = np.full(len(Zdf), np.nan); o[okr] = np.nanmean(M[okr], axis=1); return o

def discord_flags(zg, zn, thr=1.5): return (np.sign(zg) != np.sign(zn)) & (np.abs(zg) > thr) & (np.abs(zn) > thr)

def disc_rate(Zmat, gi, ni):
    G = np.nanmean(Zmat[:, gi], axis=1); N = np.nanmean(Zmat[:, ni], axis=1)
    ok = np.isfinite(G) & np.isfinite(N)
    if ok.sum() < 5: return np.nan
    return discord_flags(zs(G[ok]), zs(N[ok])).mean()

def h2c_one(vec, name, tier, nboot=NPERM_LAB):
    a = vec[c_res]; b = vec[c_dem]; a = a[np.isfinite(a)]; b = b[np.isfinite(b)]
    if len(a) < 3 or len(b) < 3: return None
    g, se = hedges_g(a, b)
    rng = np.random.default_rng(SEED + 61)
    bs = np.array([hedges_g(rng.choice(a, len(a)), rng.choice(b, len(b)))[0] for _ in range(nboot)])
    tt = st.ttest_ind(a, b, equal_var=False, alternative="less")
    diff = a.mean() - b.mean(); allv = np.concatenate([a, b]); nullm = np.empty(nboot)
    for i in range(nboot):
        pm = rng.permutation(allv); nullm[i] = pm[:len(a)].mean() - pm[len(a):].mean()
    return {"variant": V, "tier": tier, "outcome": name, "n_resilient": len(a), "n_demented": len(b), "mean_resilient": float(a.mean()), "mean_demented": float(b.mean()),
            "hedges_g_resilient_minus_demented": float(g), "ci95_low_HO": float(g - 1.96 * se), "ci95_high_HO": float(g + 1.96 * se),
            "ci95_low_boot": float(np.percentile(bs, 2.5)), "ci95_high_boot": float(np.percentile(bs, 97.5)),
            "welch_t": float(tt.statistic), "welch_p_one_sided_less": float(tt.pvalue), "perm_p_one_sided_less": float((1 + (nullm <= diff).sum()) / (1 + nboot))}

# ------------------------------------------------------------------ 3. tiered cross-column analysis
ALLROWS = {k: [] for k in ["pa", "disc", "tech", "sad", "h2b", "h2b_disc", "sec", "h2c", "load"]}
TABLE = pd.DataFrame(index=pd.Index(ALL84, name="donor_id"))
MATS = {}

def run_tier(tier, src, use):
    S = {"columns": use, "n_columns": len(use)}
    Bm = pd.DataFrame({lab: to84(c[src], c["ids"]) for lab, c in cols.items()}, index=ALL84)
    Rm = pd.DataFrame({lab: to84(resid_of(c[src], c["cps"]), c["ids"]) for lab, c in cols.items()}, index=ALL84)
    Zm = Rm / Rm.std(ddof=1)
    if tier == "techResid":                      # A9: residual columns regressed on the six technical covariates (complete cases), re-standardised
        Za = Zm.copy(); okT = Tz.notna().all(axis=1).values
        for l in use:
            y = Zm[l].values; mm = np.isfinite(y) & in_scope & okT
            o_ = np.full(n84, np.nan); o_[mm] = sm.OLS(y[mm], sm.add_constant(Tz.values[mm])).fit().resid
            Za[l] = (o_ - np.nanmean(o_)) / np.nanstd(o_, ddof=1)
        Zm = Za
    MATS[tier] = (Bm, Rm, Zm)
    glia = [l for l in use if SUBS[l][1] == "glia"]; neur = [l for l in use if SUBS[l][1] == "neuron"]
    S.update({"glia_cols": glia, "neuron_cols": neur})
    fac_scores = pd.DataFrame(index=ALL84)
    # --- parallel analysis
    if len(use) >= 2:
        cc = Zm.loc[in_scope, use].dropna()
        ev = np.sort(np.linalg.eigvalsh(np.corrcoef(cc.values.T)))[::-1]
        rng = np.random.default_rng(SEED + 31); nulls = np.empty((NPERM_PA, len(use))); X0 = cc.values
        for b in range(NPERM_PA):
            Xp = np.column_stack([rng.permutation(X0[:, j]) for j in range(X0.shape[1])])
            nulls[b] = np.sort(np.linalg.eigvalsh(np.corrcoef(Xp.T)))[::-1]
        p95 = np.percentile(nulls, 95, axis=0); nfac = 0
        for j in range(len(ev)):
            if ev[j] > p95[j]: nfac += 1
            else: break
        for j in range(len(ev)):
            ALLROWS["pa"].append({"variant": V, "tier": tier, "component": j + 1, "eigenvalue": float(ev[j]), "null_mean": float(nulls[:, j].mean()), "null_p95": float(p95[j]),
                                  "above_p95": bool(ev[j] > p95[j]), "perm_p": float((1 + (nulls[:, j] >= ev[j]).sum()) / (1 + NPERM_PA)),
                                  "n_donors_complete": len(cc), "columns": ",".join(use)})
        S.update({"pa_n_donors": int(len(cc)), "pa_factors": int(nfac), "pa_eigenvalues": [float(x) for x in ev], "pa_null_p95": [float(x) for x in p95], "pa_pass": bool(nfac >= 2)})
        Xc = (cc.values - cc.values.mean(0)) / cc.values.std(0, ddof=1)
        U_, S_, Vt_ = np.linalg.svd(Xc, full_matrices=False)
        nk = min(3, len(use))
        fs = pd.DataFrame(U_[:, :nk] * S_[:nk], index=cc.index, columns=[f"factor{i + 1}" for i in range(nk)])
        sign = np.sign(Vt_[:nk].sum(axis=1)); sign[sign == 0] = 1
        fs = fs * sign
        fac_scores = fs.reindex(ALL84)
        for i in range(nk):
            for j, l in enumerate(use):
                ALLROWS["load"].append({"variant": V, "tier": tier, "factor": f"factor{i + 1}", "column": l, "loading": float(Vt_[i, j] * sign[i] * S_[i] / np.sqrt(len(cc) - 1))})
        S["pa_var_explained_first3"] = [float(x) for x in (S_[:nk] ** 2 / (S_ ** 2).sum())]
    else:
        S.update({"pa_factors": 0, "pa_pass": False, "pa_note": "fewer than 2 columns"})
    # --- composites and discordance
    G, N = composite(Zm, glia), composite(Zm, neur)
    okd = np.isfinite(G) & np.isfinite(N) & in_scope
    zG = np.full(n84, np.nan); zN = np.full(n84, np.nan)
    flag = np.zeros(n84, bool)
    if glia and neur and okd.sum() > 10:
        zG[okd] = zs(G[okd]); zN[okd] = zs(N[okd]); flag[okd] = discord_flags(zG[okd], zN[okd])
        obs_rate = flag[okd].mean()
        gi = [use.index(l) for l in glia]; ni = [use.index(l) for l in neur]
        Zu = Zm.loc[in_scope, use].values; rng = np.random.default_rng(SEED + 41)
        perm_rates = np.empty(NPERM_DISC)
        for b in range(NPERM_DISC):
            Zp = Zu.copy()
            for j in range(Zp.shape[1]):
                m = np.isfinite(Zp[:, j]); Zp[m, j] = rng.permutation(Zp[m, j])
            perm_rates[b] = disc_rate(Zp, gi, ni)
        rrv = np.array([max(0.0, cols[l]["rr" if src == "B" else "rr_t"]) for l in use])
        nD = int(in_scope.sum()); sims = np.empty(NSIM_REL); rng = np.random.default_rng(SEED + 43)
        for b in range(NSIM_REL):
            F = rng.standard_normal(nD); E = rng.standard_normal((nD, len(use)))
            sims[b] = disc_rate(np.sqrt(rrv)[None, :] * F[:, None] + np.sqrt(1 - rrv)[None, :] * E, gi, ni)
        d = {"variant": V, "tier": tier, "n_donors_with_both": int(okd.sum()), "n_discordant": int(flag.sum()), "observed_rate": float(obs_rate),
             "perm_null_mean": float(np.nanmean(perm_rates)), "perm_null_p95": float(np.nanpercentile(perm_rates, 95)), "perm_p": float((1 + (perm_rates >= obs_rate).sum()) / (1 + NPERM_DISC)),
             "relmatched_null_mean": float(np.nanmean(sims)), "relmatched_null_p95": float(np.nanpercentile(sims, 95)), "relmatched_p": float((1 + (sims >= obs_rate).sum()) / (1 + NSIM_REL)),
             "corr_zG_zN": pr(zG, zN), "n_glia_sicker_neuron_healthier": int((flag & (zG > 0)).sum()), "n_neuron_sicker_glia_healthier": int((flag & (zN > 0)).sum()),
             "criterion_B_pass": bool(obs_rate >= 0.10 and float(np.nanmean(sims)) <= 0.05), "mean_resid_reliability_used": float(rrv.mean())}
        ALLROWS["disc"].append(d); S["discordance"] = d
    else:
        S["discordance"] = {"note": "not evaluable: needs at least one glial and one neuronal column"}
    delta = zG - zN
    nrel = int(rel["tech_passes_reliability"].sum()) if src == "Bt" else int(rel["passes_reliability"].sum())
    S["n_pass_reliability_source"] = nrel
    S["H2_reliability_criterion_pass"] = bool(nrel >= 4)
    S["H2_structure_criterion_A_pass"] = bool(S.get("pa_pass", False))
    S["H2_structure_criterion_B_pass"] = bool(S["discordance"].get("criterion_B_pass", False))
    S["H2_pass"] = bool(S["H2_reliability_criterion_pass"] and (S["H2_structure_criterion_A_pass"] or S["H2_structure_criterion_B_pass"]))
    # --- technical covariates
    axes = {f"Z_{l}": Zm[l].values for l in use}
    axes.update({"zG_glia": zG, "zN_neuron": zN, "delta_zG_minus_zN": delta})
    for c_ in fac_scores.columns: axes[f"PAfactor_{c_}"] = fac_scores[c_].reindex(ALL84).values
    flagged = []
    for name, y in axes.items():
        m = np.isfinite(y) & in_scope
        if m.sum() < 20 or np.nanstd(y[m]) == 0: continue
        fit = sm.OLS(y[m], sm.add_constant(Tz.values[m])).fit()
        sps = [st.spearmanr(Tz.values[m][:, j], y[m]) for j in range(Tz.shape[1])]
        pmin = min(s.pvalue for s in sps); jm = int(np.argmin([s.pvalue for s in sps]))
        tr = bool(fit.f_pvalue < 0.05 or pmin < 0.05 / 6)
        if tr: flagged.append(name)
        ALLROWS["tech"].append({"variant": V, "tier": tier, "axis": name, "n": int(m.sum()), "adj_r2": float(fit.rsquared_adj), "joint_F": float(fit.fvalue), "joint_F_p": float(fit.f_pvalue),
                                "strongest_covariate": Tz.columns[jm], "strongest_spearman": float(sps[jm].statistic), "strongest_spearman_p": float(pmin), "tracks_technical": tr,
                                **{f"beta_{c}": float(b) for c, b in zip(Tz.columns, fit.params[1:])}})
    S["technical_axes_flagged"] = flagged
    # --- SAD recovery (only variant all)
    if V == "all" and len(use) >= 2:
        for name, M in [("B", Bm[use]), ("Z", Zm[use])]:
            Zs_ = (M - M.mean()) / M.std(ddof=1); Zs_ = Zs_.fillna(Zs_.median())
            mcd = MinCovDet(support_fraction=0.75, random_state=0).fit(Zs_.values)
            dist = mcd.mahalanobis(Zs_.values); auc = roc_auc_score(sad_all, dist)
            top = np.argsort(-dist)[:11]; rng = np.random.default_rng(SEED + 51)
            nl = np.array([roc_auc_score(rng.permutation(sad_all), dist) for _ in range(NPERM_LAB)])
            r_ = {"variant": V, "tier": tier, "matrix": name, "n_columns": len(use), "n_sad": int(sad_all.sum()), "auc_robust_distance": float(auc),
                  "precision_at_11": float(sad_all[top].mean()), "n_sad_in_top11": int(sad_all[top].sum()), "perm_p": float((1 + (nl >= auc).sum()) / (1 + NPERM_LAB)),
                  "auc_mean_z": float(roc_auc_score(sad_all, Zs_.values.mean(axis=1)))}
            ALLROWS["sad"].append(r_)
            if name == "B": S["sad_recovery"] = r_
    # --- H2b
    def h2b_model(y, name):
        d = pd.concat([pd.Series(y, index=ALL84, name="y"), cps_s, preds], axis=1)[in_scope].dropna()
        if len(d) < 30 or d["y"].std() == 0: return None
        Xp = (d[preds.columns] - d[preds.columns].mean()) / d[preds.columns].std(ddof=1)
        m1 = sm.OLS(d["y"], sm.add_constant(pd.concat([d[["CPS_Local"]], Xp], axis=1))).fit(); m0 = sm.OLS(d["y"], sm.add_constant(d[["CPS_Local"]])).fit()
        an = anova_lm(m0, m1); pv = m1.pvalues[preds.columns]; holm = multipletests(pv.values, method="holm")[1]
        o = []
        for i, c in enumerate(preds.columns):
            ci = m1.conf_int().loc[c]
            o.append({"variant": V, "tier": tier, "outcome": name, "predictor": c, "n": int(len(d)), "coef_per_sd": float(m1.params[c]), "ci_low": float(ci[0]), "ci_high": float(ci[1]),
                      "p": float(pv[c]), "p_holm": float(holm[i]), "omnibus_F": float(an["F"].iloc[1]), "omnibus_p": float(an["Pr(>F)"].iloc[1]),
                      "adj_r2_model": float(m1.rsquared_adj), "adj_r2_cps_only": float(m0.rsquared_adj), "delta_adj_r2": float(m1.rsquared_adj - m0.rsquared_adj)})
        return o
    outs = {"delta_zG_minus_zN": delta, "abs_delta": np.abs(delta), "zG_glia": zG, "zN_neuron": zN}
    for l in use: outs[f"Z_{l}"] = Zm[l].values
    for name, y in outs.items():
        r_ = h2b_model(y, name)
        if r_: ALLROWS["h2b"] += r_
    pri = [r for r in ALLROWS["h2b"] if r["tier"] == tier and r["variant"] == V and r["outcome"] == "delta_zG_minus_zN"]
    if pri:
        ab = [r for r in pri if r["predictor"] == "ABeta_minus_pTau"][0]
        S["H2b"] = {"omnibus_F": pri[0]["omnibus_F"], "omnibus_p": pri[0]["omnibus_p"], "delta_adj_r2": pri[0]["delta_adj_r2"], "n": pri[0]["n"],
                    "abeta_minus_ptau_coef": ab["coef_per_sd"], "abeta_minus_ptau_p": ab["p"], "pass": bool(pri[0]["omnibus_p"] < 0.05 or ab["p"] < 0.05),
                    "predictors": {r["predictor"]: {"coef": r["coef_per_sd"], "p": r["p"], "p_holm": r["p_holm"]} for r in pri}}
        if flag.any():
            for c in preds.columns:
                x = preds[c].values; m = np.isfinite(x) & okd
                if m.sum() > 10 and (m & flag).sum() >= 2:
                    u = st.mannwhitneyu(x[m & flag], x[m & ~flag])
                    ALLROWS["h2b_disc"].append({"variant": V, "tier": tier, "predictor": c, "n_discordant": int((m & flag).sum()), "mean_discordant": float(np.mean(x[m & flag])),
                                                "mean_other": float(np.mean(x[m & ~flag])), "mannwhitney_p": float(u.pvalue)})
    second = {"apoe4_carrier": dd["apoe4_carrier"].astype(float), "dementia": dd["dementia"].astype(float), "last_CASI": dd["Last CASI Score"].astype(float), "last_MMSE": dd["Last MMSE Score"].astype(float)}
    for name, y in {"delta_zG_minus_zN": delta, "zG_glia": zG, "zN_neuron": zN, **{f"Z_{l}": Zm[l].values for l in use}}.items():
        for pn, pv_ in second.items():
            d = pd.concat([pd.Series(y, index=ALL84, name="y"), cps_s, pv_.rename("x")], axis=1)[in_scope].dropna()
            if len(d) < 30 or d["y"].std() == 0: continue
            m_ = sm.OLS(d["y"], sm.add_constant(d[["x", "CPS_Local"]])).fit(); ci = m_.conf_int().loc["x"]
            ALLROWS["sec"].append({"variant": V, "tier": tier, "outcome": name, "predictor": pn, "n": int(len(d)), "coef": float(m_.params["x"]), "ci_low": float(ci[0]), "ci_high": float(ci[1]), "p": float(m_.pvalues["x"])})
    # --- H2c
    nc = composite(Zm, neur); gc = composite(Zm, glia)
    all_n = composite(Zm, [l for l in NEURON if l in cols]); all_g = composite(Zm, [l for l in GLIA if l in cols])
    for name, vec in [("neuronal_composite_tier_columns(C5)", nc), ("neuronal_composite_all_columns", all_n), ("glia_composite_tier_columns", gc), ("glia_composite_all_columns", all_g)] + \
                     [(f"Z_{l}", Zm[l].values) for l in cols]:
        r_ = h2c_one(np.asarray(vec, float), name, tier)
        if r_: ALLROWS["h2c"].append(r_)
    for name, vec in [("neuronal_composite_tier_columns(C5)", nc), ("neuronal_composite_all_columns", all_n)]:
        d = pd.DataFrame({"y": vec, "resilient": res_.astype(float), "CPS": cps_all}, index=ALL84)[c_res | c_dem].dropna()
        if len(d) > 10 and d["y"].std() > 0:
            m_ = sm.OLS((d["y"] - d["y"].mean()) / d["y"].std(ddof=1), sm.add_constant(d[["resilient", "CPS"]])).fit(); ci = m_.conf_int().loc["resilient"]
            ALLROWS["h2c"].append({"variant": V, "tier": tier, "outcome": name + "|CPS-adjusted OLS (standardised outcome)", "n_resilient": int(d.resilient.sum()), "n_demented": int((1 - d.resilient).sum()),
                                   "hedges_g_resilient_minus_demented": float(m_.params["resilient"]), "ci95_low_HO": float(ci[0]), "ci95_high_HO": float(ci[1]),
                                   "welch_p_one_sided_less": float(m_.pvalues["resilient"] / 2 if m_.params["resilient"] < 0 else 1 - m_.pvalues["resilient"] / 2)})
    c5 = [r for r in ALLROWS["h2c"] if r["tier"] == tier and r["variant"] == V and r["outcome"] == "neuronal_composite_tier_columns(C5)"]
    if c5:
        S["H2c"] = c5[0]
        nm = "neuronal_composite_tier_columns(C5)"
        S["H2c_neuron_axis_tracks_technical"] = bool(any(a in flagged for a in [f"Z_{l}" for l in neur]) or "zN_neuron" in flagged)
    # --- donor-level columns for the table
    TABLE[f"zG__{tier}"] = zG; TABLE[f"zN__{tier}"] = zN; TABLE[f"delta__{tier}"] = delta; TABLE[f"discordant__{tier}"] = flag.astype(int)
    TABLE[f"discordance_label__{tier}"] = np.where(flag & (zG > 0), "glia sicker / neurons healthier than CPS predicts",
                                           np.where(flag & (zN > 0), "neurons sicker / glia healthier than CPS predicts", np.where(okd, "not discordant", "not scored")))
    for c_ in fac_scores.columns: TABLE[f"{c_}__{tier}"] = fac_scores[c_].reindex(ALL84).values
    SUMMARY["tiers"][tier] = S
    print(f"[{tier}] cols={use} PA factors={S.get('pa_factors')} disc={ {k: S['discordance'].get(k) for k in ['observed_rate', 'perm_null_mean', 'relmatched_null_mean', 'perm_p', 'relmatched_p']} } "
          f"H2={S['H2_pass']}  flagged={flagged}", flush=True)

for tier, (src, use) in TIERS.items():
    run_tier(tier, src, use)

# ------------------------------------------------------------------ 4. write
for key, fn in [("pa", "h2_parallel_analysis"), ("disc", "h2_discordance"), ("tech", "h2_technical_regression"), ("sad", "h2_sad_recovery"), ("h2b", "h2b_copathology_regressions"),
                ("h2b_disc", "h2b_discordant_vs_other"), ("sec", "h2b_secondary_cognition_apoe"), ("h2c", "h2c_resilience"), ("load", "h2_factor_loadings")]:
    pd.DataFrame(ALLROWS[key]).to_csv(f"{OUT}/{fn}.csv", index=False)

base = pd.DataFrame(index=pd.Index(ALL84, name="donor_id"))
base["in_variant"] = in_scope.astype(int); base["CPS_Local"] = cps_all; base["CPS_Local_ABeta"] = dd["CPS_Local_ABeta"].values; base["CPS_Local_pTau"] = dd["CPS_Local_pTau"].values
base["severe_donor"] = sad_all.astype(int); base["cognitive_status"] = dd["Cognitive Status"].values; base["adnc"] = dd["Overall AD neuropathological Change"].values
base["resilient_high_adnc"] = (hi & res_).astype(int); base["demented_high_adnc"] = (hi & dem_).astype(int)
base["apoe_genotype"] = dd["APOE Genotype"].values; base["apoe4_carrier"] = dd["apoe4_carrier"].values
base["age"] = dd["age"].values; base["sex"] = dd["Sex"].values; base["pmi"] = dd["pmi"].values
for c in preds.columns: base[c] = preds[c].values
for c in Tz.columns: base[f"tech_{c}"] = T[c].values
Bm1, Rm1, Zm1 = MATS["allcols"]; Bm3, Rm3, Zm3 = MATS["techNAM"]; Bm4, Rm4, Zm4 = MATS["techResid"]
for l in cols:
    base[f"n_nuclei_{l}"] = to84(cols[l]["n_nuclei"], cols[l]["ids"])
    base[f"B_{l}"] = Bm1[l].values; base[f"R_{l}"] = Rm1[l].values; base[f"Z_{l}"] = Zm1[l].values
    base[f"Btech_{l}"] = Bm3[l].values; base[f"Ztech_{l}"] = Zm3[l].values; base[f"Zadj_{l}"] = Zm4[l].values
pd.concat([base, TABLE], axis=1).to_csv(f"{OUT}/donor_bricks.csv")
for name, key in [("B_k10", "B10"), ("B_k20", "B20"), ("B_covNAM", "Bcov"), ("B_foldPCA", "Bfold"), ("B_supertype_ridge", "Bsup"), ("B_techNAM", "Bt")]:
    pd.DataFrame({lab: to84(c[key], c["ids"]) for lab, c in cols.items()}, index=pd.Index(ALL84, name="donor_id")).to_csv(f"{OUT}/brick_scores_{name}.csv")

SUMMARY["seconds"] = round(time.perf_counter() - t0, 1)
json.dump(SUMMARY, open(f"{OUT}/summary.json", "w"), indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
print("DONE", V, SUMMARY["seconds"], "s", flush=True)
