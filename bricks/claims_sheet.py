#!/usr/bin/env python
"""Claims sheet (protocol 06 section 7): every number is read from a cell of a results/*.csv file, never typed.

Writes results/claims_sheet.csv (claim, number, source file, row selector, column, status, note) and results/claims_sheet.md (the table pasted in the report).
usage: python bricks/claims_sheet.py     (python with pandas)
"""
import os, json
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
R = os.path.join(HERE, "..", "results")
def rd(f): return pd.read_csv(os.path.join(R, f))
def sci(x):
    return f"{x:.1e}" if x < 1e-3 else (f"{x:.4f}" if x < 0.01 else f"{x:.3f}".rstrip("0").rstrip("."))
C = []
def add(cid, claim, number, f, sel, col, status, note=""):
    C.append({"id": cid, "claim": claim, "number": number, "source_csv": f, "row_selector": sel, "column": col, "status": status, "note": note})

# ---- dataset facts
dt = rd("seaad_mtg_donor_table.csv")
add("D1", "84 donors, middle temporal gyrus (SEA-AD), all open data", f"{len(dt)}", "seaad_mtg_donor_table.csv", "count of rows", "donor_id", "allowed")
add("D2", "1,178,694 nuclei from SEA-AD in those 84 donors (final-nuclei metadata)", f"{int(dt.n_nuclei.sum()):,}", "seaad_mtg_donor_table.csv", "sum over the 84 rows", "n_nuclei", "allowed")
fs = rd("cna/floor1_sampling.csv")
add("D3", "149,995 nuclei on the Mantis cells floor (11 subclasses analysed; 7 on the floor)", f"{int(fs.n_selected.sum()):,}", "cna/floor1_sampling.csv", "sum over rows", "n_selected", "allowed")
sa = rd("stageA_summary.csv").set_index("item")["value"]
add("D4", "236,973 Alzheimer's abstracts with an abstract in the corpus (the protocol's pilot figure was 237,059)", f"{int(float(sa['n_corpus'])):,}", "stageA_summary.csv", "item = n_corpus", "value", "allowed", "use the measured value")
us = rd("universe_summary.csv").set_index("item")["value"]
add("D5", "13,006 genes linked to AD in Open Targets, 2,177 with genetic evidence (the protocol's pilot figure was 2,179)", f"{int(float(us['OT targets'])):,} / {int(float(us['OT targets with genetic_association > 0'])):,}",
    "universe_summary.csv", "items OT targets; OT targets with genetic_association > 0", "value", "allowed", "use the measured values")

# ---- reproductions (Day 1)
h1 = rd("cna/h1_global_p.csv").set_index("subclass")
add("R1", "Reproduced: early Sst loss, CNA global P (CPS_Local, 84 donors, 10,000 permutations)", sci(h1.loc["Sst", "P_CPS"]), "cna/h1_global_p.csv", "subclass = Sst", "P_CPS", "allowed")
add("R2", "Reproduced: microglial expansion, CNA global P (the permutation floor 1/10,001)", "< 1e-4", "cna/h1_global_p.csv", "subclass = Microglia-PVM", "P_CPS", "allowed", f"stored value {h1.loc['Microglia-PVM', 'P_CPS']:.2e}")
add("R3", "Subclasses with a global CNA association with CPS_Local at P < 0.05", f"{int((h1.P_CPS < 0.05).sum())} of {len(h1)}", "cna/h1_global_p.csv", "all rows", "P_CPS", "allowed")
add("R4", "... and without the 11 Severely Affected Donors", f"{int((h1.P_CPS_nosad < 0.05).sum())} of {len(h1)}", "cna/h1_global_p.csv", "all rows", "P_CPS_nosad", "allowed", "state it with R3")
cc = rd("cna/pc1_sccoda_concordance.csv")
cr = cc[cc.credible_primary == True]
add("R5", "Sign agreement with SEA-AD's scCODA calls among supertypes SEA-AD calls credible", f"{int(cr.agree_primary.sum())} of {len(cr)} ({cr.agree_primary.mean():.0%})", "cna/pc1_sccoda_concordance.csv", "rows with credible_primary = True", "agree_primary", "allowed")

# ---- measured gaps (literature)
v1 = rd("v1_pilot_reproduction.csv").set_index("check")
add("L1", "Since 2010, about a third of AD abstracts mention amyloid", f"{v1.loc['amyloid_share_pooled_2010_2025', 'observed']:.1%}", "v1_pilot_reproduction.csv", "check = amyloid_share_pooled_2010_2025", "observed", "allowed")
s7 = rd("s7_amyloid_pathway_comparison.csv").set_index("group")
g_ex = s7.loc["genes with G>=0.7 excluding amyloid-pathway set"]; g_am = s7.loc[[i for i in s7.index if i.startswith("amyloid-pathway set")][0]]
add("L2", "Genes with genetic score >= 0.7 (excluding the amyloid-pathway set) average this many AD abstracts, against the amyloid-pathway genes", f"{g_ex.mean_A:,.0f} (median {g_ex.median_A:,.0f}, n = {int(g_ex.n)}) vs {g_am.mean_A:,.0f} (median {g_am.median_A:,.0f}, n = {int(g_am.n)})",
    "s7_amyloid_pathway_comparison.csv", "groups 'genes with G>=0.7 excluding amyloid-pathway set' and 'amyloid-pathway set ...'", "mean_A, median_A", "allowed")
ch = rd("contrast_and_headline_counts.csv").set_index("item")
add("L3", "Genetically supported genes (score >= 0.7) in sparsely studied territory (fewer than 100 AD abstracts)", f"{int(ch.loc['genes with G>=0.7 and A<100 (sparsely studied)', 'A'])} of {int(ch.loc['genes with G>=0.7', 'A'])}",
    "contrast_and_headline_counts.csv", "items 'genes with G>=0.7 and A<100 (sparsely studied)' and 'genes with G>=0.7'", "A", "allowed")
add("L4", "Most-mentioned gene (APP) versus the least-mentioned strongly supported gene (SLC26A1)", f"{int(ch.loc['ratio most-mentioned / least-mentioned strongly supported', 'A']):,} to 1", "contrast_and_headline_counts.csv", "item ratio most-mentioned / least-mentioned strongly supported", "A", "allowed")

# ---- Day 2: bricks
rel = rd("bricks/all/brick_reliability.csv").set_index("label")
add("B1", "Split-half reliability of donor brick scores across the 11 subclasses (median over 20 splits)", f"{rel.split_half_r_median.min():.2f} to {rel.split_half_r_median.max():.2f}", "bricks/all/brick_reliability.csv", "all 11 rows", "split_half_r_median", "allowed", "scores are reproducible within a donor")
add("B2", "Out-of-fold correlation of the brick score with CPS_Local across the 11 subclasses (Spearman)", f"{rel.oof_rho.min():.2f} to {rel.oof_rho.max():.2f}; {int(rel.retained.sum())} subclasses pass rho > 0.4 (" + ", ".join(rel.index[rel.retained].tolist()) + ")",
    "bricks/all/brick_reliability.csv", "all 11 rows", "oof_rho, retained", "allowed")
pa = rd("bricks/all/h2_parallel_analysis.csv"); p1 = pa[(pa.tier == "allcols") & (pa.component == 1)].iloc[0]; p2 = pa[(pa.tier == "allcols") & (pa.component == 2)].iloc[0]
add("B3", "Residual structure after removing what CPS predicts: factors beyond a permutation null (parallel analysis, 11 subclasses)", f"1 (first eigenvalue {p1.eigenvalue:.2f} vs null 95th percentile {p1.null_p95:.2f}; second {p2.eigenvalue:.2f} vs {p2.null_p95:.2f})",
    "bricks/all/h2_parallel_analysis.csv", "tier = allcols, component 1 and 2", "eigenvalue, null_p95", "allowed", "one shared axis, not separate bricks")
dc = rd("bricks/all/h2_discordance.csv"); dca = dc[dc.tier == "allcols"].iloc[0]
add("B4", "Donors whose glia and neuron residuals have opposite signs, both |z| > 1.5 (preregistered discordance rule), tier allcols", f"{int(dca.n_discordant)} of {int(dca.n_donors_with_both)} (null expectation {dca.perm_null_mean:.1%})",
    "bricks/all/h2_discordance.csv", "tier = allcols", "n_discordant, n_donors_with_both, perm_null_mean", "NOT SUPPORTED", "H2 failed: do not write 'in N of 84 donors one cell type looks sicker while another looks healthier'")
tr = rd("bricks/all/h2_technical_regression.csv"); tra = tr[(tr.tier == "allcols") & tr.axis.str.startswith("Z_")]
add("B5", "Residual columns that track a technical covariate (joint F-test P < 0.05 or a single covariate at P < 0.05/6)", f"{int(tra.tracks_technical.sum())} of {len(tra)}; strongest: sequencing saturation (Spearman up to {tra.strongest_spearman.abs().max():.2f})",
    "bricks/all/h2_technical_regression.csv", "tier = allcols, axes Z_*", "tracks_technical, strongest_spearman", "allowed", "any residual-axis claim is dropped by the preregistered rule")
b2 = rd("bricks/all/h2b_copathology_regressions.csv"); b2a = b2[(b2.tier == "allcols") & (b2.outcome == "delta_zG_minus_zN")]
add("B6", "H2b: LATE, Lewy, arteriolosclerosis, microinfarcts, CAA and the Abeta-minus-pTau difference jointly explain the glia-minus-neuron residual (omnibus F-test)", f"P = {b2a.omnibus_p.iloc[0]:.2f}, adjusted R^2 gain {b2a.delta_adj_r2.iloc[0]:.3f}",
    "bricks/all/h2b_copathology_regressions.csv", "tier = allcols, outcome = delta_zG_minus_zN", "omnibus_p, delta_adj_r2", "NOT SUPPORTED")
h2c = rd("bricks/all/h2c_resilience.csv"); r = h2c[(h2c.tier == "prereg") & (h2c.outcome == "neuronal_composite_tier_columns(C5)")].iloc[0]
add("B7", "H2c, an association with n = 11: High-ADNC donors without dementia vs with dementia (n = 31), neuronal brick residual (Lamp5, L4 IT, L5 IT), Hedges' g (resilient minus demented)",
    f"g = {r.hedges_g_resilient_minus_demented:.2f} (95% CI {r.ci95_low_HO:.2f} to {r.ci95_high_HO:+.2f}; bootstrap {r.ci95_low_boot:.2f} to {r.ci95_high_boot:+.2f}); one-sided Welch P = {r.welch_p_one_sided_less:.3f}",
    "bricks/all/h2c_resilience.csv", "tier = prereg, outcome = neuronal_composite_tier_columns(C5)", "hedges_g_resilient_minus_demented, ci95_*, welch_p_one_sided_less", "CAVEATED",
    "the residual axis tracks sequencing saturation (preregistered rule: dropped from claims); see B8 for the analysis without the 11 severely affected donors")
n = rd("bricks/nosad/h2c_resilience.csv"); rn = n[(n.tier == "allcols") & (n.outcome == "neuronal_composite_all_columns")].iloc[0]
add("B8", "H2c without the 11 Severely Affected Donors (10 vs 21)", f"g = {rn.hedges_g_resilient_minus_demented:.2f} (95% CI {rn.ci95_low_HO:.2f} to {rn.ci95_high_HO:+.2f}); one-sided Welch P = {rn.welch_p_one_sided_less:.2f}",
    "bricks/nosad/h2c_resilience.csv", "tier = allcols, outcome = neuronal_composite_all_columns", "hedges_g_resilient_minus_demented, ci95_*, welch_p_one_sided_less", "allowed", "shows the effect depends on the extreme donors")
sd = rd("bricks/all/h2_sad_recovery.csv"); sdr = sd[(sd.tier == "prereg") & (sd.matrix == "B")].iloc[0]
add("B9", "Validation: SEA-AD Severely Affected Donors recovered as outliers (robust Mahalanobis distance of the brick scores, 3 retained subclasses), AUC", f"{sdr.auc_robust_distance:.2f} (permutation P = {sdr.perm_p:.3f}); mean-z AUC {sdr.auc_mean_z:.2f}",
    "bricks/all/h2_sad_recovery.csv", "tier = prereg, matrix = B", "auc_robust_distance, perm_p, auc_mean_z", "allowed")

# ---- Day 3: genetics
pr = rd("genetics/h3_celltype_property.csv"); m = pr[(pr.tag == "primary") & (pr.VARIABLE == "Microglia")].iloc[0]
add("G1", "Microglia-specific genes are enriched for AD risk (MAGMA gene-property test on Bellenguez 2022, one-sided)", f"P = {m.P:.1e} (beta {m.BETA:.3f}, SE {m.SE:.4f}, {int(m.NGENES):,} genes)", "genetics/h3_celltype_property.csv", "tag = primary, VARIABLE = Microglia", "P, BETA", "allowed", "remains P < 1e-4 with effective N and without chr19:44-46 Mb")
sp = rd("genetics/h3_program_sets.csv")
mr = sp[(sp.tag == "primary") & (sp.VARIABLE == "Microglia_red") & (sp.analysis == "conditional")].iloc[0]
add("G2", "The pathology-red microglial program is enriched for AD GWAS genes beyond microglial specificity (conditional MAGMA, one-sided)", f"NO: beta {mr.BETA:.3f}, P = {mr.P:.2f} ({int(mr.NGENES)} genes)", "genetics/h3_program_sets.csv", "tag = primary, VARIABLE = Microglia_red, analysis = conditional", "BETA, P", "NOT SUPPORTED")
neu = sp[(sp.tag == "primary") & (sp.kind == "primary") & (sp.analysis == "conditional") & (sp.neuronal == True)]
add("G3", "Neuronal red and blue programs carry no AD genetic fingerprint: smallest conditional P across the 14 programs", f"P = {neu.P.min():.2f}; none below 0.05 (Holm-adjusted minimum {neu.holm_p_within_neuronal_primary.min():.2f})", "genetics/h3_program_sets.csv",
    "tag = primary, kind = primary, analysis = conditional, neuronal subclasses", "P, holm_p_within_neuronal_primary", "allowed", "'damage without genetic fingerprints' is a prediction that held; absence of evidence at this power")
op = sp[(sp.tag == "primary") & (sp.VARIABLE == "OPC_red") & (sp.analysis == "conditional")].iloc[0]
add("G4", "Exploratory: the OPC red program is nominally enriched for AD GWAS genes beyond OPC specificity (not corrected for 22 tests; not preregistered as a claim)", f"P = {op.P:.4f} (Holm-adjusted over 22 would be {min(1, op.P * 22):.2f})", "genetics/h3_program_sets.csv", "tag = primary, VARIABLE = OPC_red, analysis = conditional", "P", "EXPLORATORY")
ot = rd("genetics/ot_overlap.csv"); o4 = ot[(ot.subclass == "l4it") & (ot.program == "red") & (~ot.exclude_APOE_region)].iloc[0]; om = ot[(ot.subclass == "microglia") & (ot.program == "red") & (~ot.exclude_APOE_region)].iloc[0]
add("G5", "Cross-check: overlap of red programs with the 173 Open Targets genes (score >= 0.5), versus expression-matched random gene sets", f"L4 IT red {int(o4.overlap_obs)} genes vs {o4.overlap_null_mean:.1f} expected (P = {o4.p_one_sided:.4f}); microglia red {int(om.overlap_obs)} vs {om.overlap_null_mean:.1f} (P = {om.p_one_sided:.3f})",
    "genetics/ot_overlap.csv", "subclass, program = red, exclude_APOE_region = False", "overlap_obs, overlap_null_mean, p_one_sided", "EXPLORATORY", "22 tests; this cross-check does not condition on cell-type specificity")

# ---- H4
h4 = rd("genetics/h4/h4_ols.csv").iloc[0]
add("H4a", "Replication attempt, APOE4 carriers (25) vs non-carriers (59), oligodendrocyte cholesterol-biosynthesis module score (adjusted for CPS, age, sex, PMI)", f"{h4.coef_apoe4:+.2f} z units (95% CI {h4.ci95_low:+.2f} to {h4.ci95_high:+.2f}), one-sided P = {h4.p_one_sided_greater:.2f}",
    "genetics/h4/h4_ols.csv", "analysis = primary", "coef_apoe4, ci95_low, ci95_high, p_one_sided_greater", "allowed", "not replicated: direction positive, effect not distinguishable from zero")
ms = rd("genetics/h4/h4_pydeseq2_module_sum.csv").iloc[0]; pg = rd("genetics/h4/h4_pydeseq2_genes.csv")
add("H4b", "pydeseq2 on the 20 module genes: summed-module log2 fold change (carriers vs non-carriers); genes with positive fold change; genes nominally significant", f"{ms.log2FoldChange:+.2f} (95% CI {ms.ci95_low:+.2f} to {ms.ci95_high:+.2f}), one-sided P = {ms.p_one_sided_greater:.2f}; {int((pg.log2FoldChange > 0).sum())} of 20 positive, {int(((pg.log2FoldChange > 0) & (pg.p_one_sided_greater < 0.05)).sum())} at P < 0.05",
    "genetics/h4/h4_pydeseq2_module_sum.csv; genetics/h4/h4_pydeseq2_genes.csv", "feature = MODULE_SUM; all 20 genes", "log2FoldChange, ci95_*, p_one_sided_greater", "allowed")

df = pd.DataFrame(C)
df.to_csv(os.path.join(R, "claims_sheet.csv"), index=False)
with open(os.path.join(R, "claims_sheet.md"), "w", encoding="utf-8") as f:
    f.write("| ID | Claim text | Number | results/*.csv cell | Status |\n|---|---|---|---|---|\n")
    for r_ in C:
        f.write(f"| {r_['id']} | {r_['claim']} | {r_['number']} | `{r_['source_csv']}` ({r_['row_selector']}; column {r_['column']}) | {r_['status']}{(': ' + r_['note']) if r_['note'] else ''} |\n")
print(df[["id", "number", "status"]].to_string())
