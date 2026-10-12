#!/usr/bin/env python
"""H4: APOE4 carriers vs non-carriers on the oligodendrocyte cholesterol-biosynthesis module (preregistration 1.2; amendment A7 item 8).

Oligodendrocyte pseudobulk (SEA-AD, MTG rows of aged donors; rows = supertype x library prep; `.X` = raw counts) summed to donors.
Module score = mean over the 20 genes of z-scored (over donors) log2(CPM + 1), CPM = module counts / donor total counts x 1e6.
OLS: score ~ APOE4 carrier + CPS_Local + age + sex + PMI; one-sided test APOE4 > 0; 95% CI; standardised effect; Freedman-Lane permutation P.
pydeseq2 on the 20 raw module-gene counts (size factors from all genes), same design; plus the summed module count as one feature.
Sensitivities: donors with >= 50 oligodendrocyte nuclei; without the 11 Severely Affected Donors; without CPS; APOE e4 dose.
"""
import glob, json, os, sys, time, warnings
import numpy as np, pandas as pd, scipy.stats as st, scipy.sparse as sp, anndata as ad
import statsmodels.api as sm
warnings.filterwarnings("ignore")
R = "/home/pranavak/orcd/pool/mantis-films/research/seaad"
OUT = f"{R}/results/day3/h4"; os.makedirs(OUT, exist_ok=True)
SEED = 20261011
MODULE = ["HMGCR", "HMGCS1", "LDLR", "SREBF2", "MVD", "MVK", "FDFT1", "SQLE", "LSS", "CYP51A1", "DHCR7", "DHCR24", "INSIG1", "IDI1", "FDPS", "SC5D", "MSMO1", "NSDHL", "HSD17B7", "TM7SF2"]
t0 = time.perf_counter()

f = glob.glob(f"{R}/data/s3/Multiregion_2026/pseudobulk_objects/SEAAD_Oligodendrocyte_HIP_MEC_*pseudobulked.2026-06-22.h5ad")[0]
a = ad.read_h5ad(f, backed="r")
obs = a.obs[["Brain Region", "Neurotypical reference", "Donor ID", "Number of nuclei", "Supertype"]].copy()
m = ((obs["Brain Region"].astype(str) == "MTG") & (obs["Neurotypical reference"].astype(str) == "False")).values
idx = np.where(m)[0]
genes = pd.Index(a.var.index.astype(str))
blocks = []
for i in range(0, len(idx), 400):
    blk = idx[i:i + 400]
    xb = a.X[blk[0]:blk[-1] + 1]
    xb = xb.toarray() if sp.issparse(xb) else np.asarray(xb)
    blocks.append(xb[blk - blk[0]])
X = np.vstack(blocks); a.file.close()
o = obs.loc[m].reset_index(drop=True)
donor_counts = pd.DataFrame(X).groupby(o["Donor ID"].astype(str).values).sum()
donor_counts.columns = genes
nuc = o.groupby(o["Donor ID"].astype(str))["Number of nuclei"].sum()
print("MTG aged rows", len(o), "donors", donor_counts.shape[0], "genes", donor_counts.shape[1], "total counts median", donor_counts.sum(axis=1).median(), flush=True)
missing = [g for g in MODULE if g not in genes]
assert not missing, missing

donors = pd.read_csv(f"{R}/results/seaad_mtg_donor_table.csv", index_col="donor_id")
ids = sorted(donor_counts.index); assert len(ids) == 84 and set(ids) == set(donors.index)
donor_counts = donor_counts.loc[ids].astype(float); nuc = nuc.reindex(ids)
d = donors.loc[ids]
meta = pd.DataFrame({"apoe4": d["apoe4_carrier"].astype(int), "CPS": d["CPS_Local"], "age": d["age"], "sex_male": d["sex_male"].astype(int), "pmi": d["pmi"],
                     "severe": d["severe"].astype(int), "nuclei": nuc.values, "apoe": d["APOE Genotype"]}, index=ids)
print("APOE4 carriers", int(meta.apoe4.sum()), "non-carriers", int((1 - meta.apoe4).sum()), flush=True)

tot = donor_counts.sum(axis=1)
mod = donor_counts[MODULE]
lcpm = np.log2(mod.div(tot, axis=0) * 1e6 + 1)
zg = (lcpm - lcpm.mean()) / lcpm.std(ddof=1)
score = zg.mean(axis=1)

def ols_one(y, df, covars, label, nperm=10000, seed=SEED):
    X = sm.add_constant(df[["apoe4"] + covars].astype(float)); fit = sm.OLS(y, X).fit()
    b = fit.params["apoe4"]; se = fit.bse["apoe4"]; tval = fit.tvalues["apoe4"]
    p_one = fit.pvalues["apoe4"] / 2 if b > 0 else 1 - fit.pvalues["apoe4"] / 2
    ci = fit.conf_int().loc["apoe4"]
    # Freedman-Lane
    Z = sm.add_constant(df[covars].astype(float)) if covars else np.ones((len(df), 1))
    Zm = np.asarray(Z, float)
    H = Zm @ np.linalg.pinv(Zm); fitted0 = H @ y; res0 = y - fitted0
    rng = np.random.default_rng(seed); cnt = 0
    Xf = np.asarray(X, float); XtXi = np.linalg.inv(Xf.T @ Xf)
    for _ in range(nperm):
        ys = fitted0 + rng.permutation(res0)
        beta = XtXi @ Xf.T @ ys; r = ys - Xf @ beta
        s2 = (r @ r) / (len(ys) - Xf.shape[1]); tstar = beta[1] / np.sqrt(s2 * XtXi[1, 1])
        cnt += tstar >= tval
    pperm = (1 + cnt) / (1 + nperm)
    sd_res = np.sqrt(fit.mse_resid); sd_y = np.std(y, ddof=1)
    return {"analysis": label, "n": int(len(df)), "n_carriers": int(df.apoe4.sum()), "n_noncarriers": int((1 - df.apoe4).sum()),
            "coef_apoe4": float(b), "se": float(se), "ci95_low": float(ci[0]), "ci95_high": float(ci[1]), "t": float(tval), "p_one_sided_greater": float(p_one),
            "p_two_sided": float(fit.pvalues["apoe4"]), "perm_p_one_sided_FreedmanLane": float(pperm), "std_effect_coef_over_resid_sd": float(b / sd_res),
            "std_effect_coef_over_score_sd": float(b / sd_y), "mean_score_carriers": float(y[df.apoe4.values == 1].mean()), "mean_score_noncarriers": float(y[df.apoe4.values == 0].mean()),
            "adj_r2": float(fit.rsquared_adj)}

covs = ["CPS", "age", "sex_male", "pmi"]
rows = []
rows.append(ols_one(score.values, meta, covs, "primary (all 84 donors; CPS, age, sex, PMI)"))
ok50 = meta.nuclei >= 50
rows.append(ols_one(score[ok50].values, meta[ok50], covs, f"sensitivity: donors with >= 50 oligodendrocyte nuclei"))
nos = meta.severe == 0
rows.append(ols_one(score[nos].values, meta[nos], covs, "sensitivity: without the 11 Severely Affected Donors"))
rows.append(ols_one(score.values, meta, ["age", "sex_male", "pmi"], "sensitivity: without CPS"))
meta["log10_nuclei"] = np.log10(meta["nuclei"].clip(lower=1))
rows.append(ols_one(score.values, meta, covs + ["log10_nuclei"], "sensitivity: + log10 oligodendrocyte nuclei count"))
# unadjusted difference
rows.append(ols_one(score.values, meta, [], "unadjusted"))
# score recomputed with z-scoring on the ge50 set only is not needed; leave as is
res = pd.DataFrame(rows); res.to_csv(f"{OUT}/h4_ols.csv", index=False)
print(res[["analysis", "n", "coef_apoe4", "ci95_low", "ci95_high", "p_one_sided_greater", "perm_p_one_sided_FreedmanLane", "std_effect_coef_over_resid_sd"]].round(4).to_string(), flush=True)

# per-gene OLS (descriptive)
g_rows = []
for g in MODULE:
    y = zg[g].values
    fit = sm.OLS(y, sm.add_constant(meta[["apoe4"] + covs].astype(float))).fit()
    ci = fit.conf_int().loc["apoe4"]
    g_rows.append({"gene": g, "coef_apoe4_z_units": float(fit.params["apoe4"]), "ci95_low": float(ci[0]), "ci95_high": float(ci[1]), "p_two_sided": float(fit.pvalues["apoe4"]),
                   "mean_z_carriers": float(y[meta.apoe4.values == 1].mean()), "mean_z_noncarriers": float(y[meta.apoe4.values == 0].mean())})
pd.DataFrame(g_rows).to_csv(f"{OUT}/h4_per_gene_ols.csv", index=False)

# donor-level table for the figure
fig = meta.copy(); fig["module_score"] = score.values
fig["module_score_adjusted_for_covariates"] = sm.OLS(score.values, sm.add_constant(meta[covs].astype(float))).fit().resid
fig.to_csv(f"{OUT}/h4_donor_scores.csv")

# ---------------- pydeseq2 ----------------
try:
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    # size factors from all genes (median of ratios on genes with non-zero counts in every donor)
    allc = donor_counts.values
    nz = (allc > 0).all(axis=0)
    logc = np.log(allc[:, nz]); sf = np.exp(np.median(logc - logc.mean(axis=0), axis=1))
    sf = sf / np.exp(np.mean(np.log(sf)))
    print("size factors: genes used", int(nz.sum()), "range", sf.min().round(3), sf.max().round(3), flush=True)
    z = lambda v: (v - v.mean()) / v.std(ddof=1)
    md = pd.DataFrame({"apoe4": np.where(meta.apoe4 == 1, "carrier", "noncarrier"), "CPS": z(meta.CPS).values, "age": z(meta.age).values,
                       "sex_male": meta.sex_male.astype(str).values, "pmi": z(meta.pmi).values}, index=ids)

    def run_dds(counts_df, tag):
        dds = DeseqDataSet(counts=counts_df.round().astype(int), metadata=md, design="~ apoe4 + CPS + age + sex_male + pmi", refit_cooks=False, fit_type="mean", quiet=True)
        dds.obs["size_factors"] = sf
        dds.layers["normed_counts"] = dds.X / sf[:, None]
        dds.var["_normed_means"] = dds.layers["normed_counts"].mean(axis=0)
        with np.errstate(divide="ignore"):
            dds.logmeans = np.log(np.asarray(dds.X, float)).mean(axis=0)
        dds.fit_genewise_dispersions(); dds.fit_dispersion_trend(); dds.fit_dispersion_prior(); dds.fit_MAP_dispersions(); dds.fit_LFC(); dds.calculate_cooks()
        ds = DeseqStats(dds, contrast=["apoe4", "carrier", "noncarrier"], quiet=True, alt_hypothesis="greaterAbs")
        ds.summary()
        r = ds.results_df.copy(); r.insert(0, "feature", r.index); r["analysis"] = tag
        # one-sided P for log2FC > 0 from the Wald statistic
        r["p_one_sided_greater"] = np.where(r["stat"] > 0, r["pvalue"] / 2, 1 - r["pvalue"] / 2)
        r["ci95_low"] = r["log2FoldChange"] - 1.96 * r["lfcSE"]; r["ci95_high"] = r["log2FoldChange"] + 1.96 * r["lfcSE"]
        return r
    r1 = run_dds(donor_counts[MODULE], "pydeseq2_module_genes")
    r1.to_csv(f"{OUT}/h4_pydeseq2_genes.csv", index=False)
    print(r1[["feature", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "p_one_sided_greater"]].round(4).to_string(), flush=True)
    modsum = pd.DataFrame({"MODULE_SUM": donor_counts[MODULE].sum(axis=1)})
    r2 = run_dds(modsum, "pydeseq2_module_sum")
    r2.to_csv(f"{OUT}/h4_pydeseq2_module_sum.csv", index=False)
    print(r2.round(4).to_string(), flush=True)
    summ = {"genes_up": int((r1.log2FoldChange > 0).sum()), "genes_up_p_one_sided_lt_005": int(((r1.log2FoldChange > 0) & (r1.p_one_sided_greater < 0.05)).sum()),
            "mean_log2fc_genes": float(r1.log2FoldChange.mean()), "module_sum_log2fc": float(r2.log2FoldChange.iloc[0]),
            "module_sum_ci": [float(r2.ci95_low.iloc[0]), float(r2.ci95_high.iloc[0])], "module_sum_p_one_sided": float(r2.p_one_sided_greater.iloc[0])}
    json.dump(summ, open(f"{OUT}/h4_pydeseq2_summary.json", "w"), indent=1)
    print(summ, flush=True)
except Exception as e:
    import traceback; traceback.print_exc()
    json.dump({"error": repr(e)[:500]}, open(f"{OUT}/h4_pydeseq2_summary.json", "w"))
print("done", round(time.perf_counter() - t0, 1), "s")
