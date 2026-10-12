# Deviations log (Broken Bricks, Day 1 onward)

The preregistration (section 3) records deviations and implementation amendments in `amendments.md`, which remains the single source of truth (A1-A3: literature study; A4 onward: SEA-AD analysis). This file is the index the Day-1 task asked for: one line per entry with the reason. Anything that changes a threshold, covariate, exclusion or definition after the preregistration commit is marked **DEVIATION**; format-level and operational details are marked implementation. Entries A8 and A9 are additive deviations (extra analyses; see below); no entry changes a threshold, covariate, exclusion or decision rule.

| Entry | Kind | What | Reason |
|---|---|---|---|
| A4.1 | implementation | H1 denominator is 11 subclasses (C1 applied) | all 11 non-vascular subclasses have >= 81 donors with >= 50 nuclei |
| A4.2 | implementation | primary CNA uses all 84 donors; `nosad` and `ge50` variants as sensitivity | the >= 50 rule belongs to brick scores; "with and without severe donors" is in the protocol |
| A4.3 | implementation | Nnull = 10,000, seed 0, default ks, FDR flags from `cna` per-cell FDR | the protocol leaves them open; 1,000 permutations floors P at 0.001 |
| A4.4 | implementation | operational definitions of Sst depletion and of "credible" supertypes (SEA-AD summary table primary, own aggregation as sensitivity) | C3 left "credible" undefined; declared before the all-subclass concordance |
| A4.5 | implementation | Milo via pydeseq2 solver | R/edgeR not installed on the cluster |
| A4.7 | implementation | floor capped at 149,995 cells (18,515 per large type) | the protocol's "<= 150k" and "20k each" conflict; the cap is binding |
| A5.1-A5.5 | implementation / exploratory | post-result notes: outcome of the A4 definitions (PC1 pass on both, H1 11/11), FDR flag boundary, Milo worker cap, exploratory descriptive tables, what is committed | recorded for transparency; nothing changed |
| A6.1-A6.13 | implementation | Day-2 operational definitions (columns and donors, phenotype-free NAM with cna's routine, k = 15 PCs, ridge with inner LOO, folds from seed 20261011, split-half reliability, LOESS residuals, parallel analysis, discordance nulls, SAD recovery, technical regression, H2b model, H2c Hedges' g) | the protocol leaves them open; fixed before any brick score was computed |
| A7.1-A7.9 | implementation | Day-3 definitions (programs, donor-blocked Spearman, MAGMA settings, specificity covariates, set tests, APOE sensitivity, OT overlap, H4 model) | same; fixed before any program, MAGMA or H4 number |
| A8 | **DEVIATION (additive)** | three reporting tiers for H2: prereg (verdicts), allcols (no rho filter), techNAM (NAM adjusted for technical covariates) | the rho > 0.4 filter keeps no glial column, so the discordance arm was not evaluable; added after a debug run, before full results |
| A9 | **DEVIATION (additive)** | fourth tier techResid (residuals regressed on technical covariates) and the statement that techNAM is not a brick-score analysis | techNAM scores do not predict CPS (technical covariates and severity are confounded); written after the all-variant results |
| A10.1-A10.8 | implementation | MAGMA covariate tests derived one-sided; set-file column order; APOE sensitivity by SNP removal; pydeseq2 size-factor injection; donors-floor tier; list of what was not done | recorded for transparency |
