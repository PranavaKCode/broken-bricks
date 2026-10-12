# Deviations log (Broken Bricks, Day 1 onward)

The preregistration (section 3) records deviations and implementation amendments in `amendments.md`, which remains the single source of truth (A1-A3: literature study; A4 onward: SEA-AD analysis). This file is the index the Day-1 task asked for: one line per entry with the reason. Anything that changes a threshold, covariate, exclusion or definition after the preregistration commit is marked **DEVIATION**; format-level and operational details are marked implementation. No entry so far is a DEVIATION.

| Entry | Kind | What | Reason |
|---|---|---|---|
| A4.1 | implementation | H1 denominator is 11 subclasses (C1 applied) | all 11 non-vascular subclasses have >= 81 donors with >= 50 nuclei |
| A4.2 | implementation | primary CNA uses all 84 donors; `nosad` and `ge50` variants as sensitivity | the >= 50 rule belongs to brick scores; "with and without severe donors" is in the protocol |
| A4.3 | implementation | Nnull = 10,000, seed 0, default ks, FDR flags from `cna` per-cell FDR | the protocol leaves them open; 1,000 permutations floors P at 0.001 |
| A4.4 | implementation | operational definitions of Sst depletion and of "credible" supertypes (SEA-AD summary table primary, own aggregation as sensitivity) | C3 left "credible" undefined; declared before the all-subclass concordance |
| A4.5 | implementation | Milo via pydeseq2 solver | R/edgeR not installed on the cluster |
| A4.7 | implementation | floor capped at 149,995 cells (18,515 per large type) | the protocol's "<= 150k" and "20k each" conflict; the cap is binding |
| A5.1-A5.5 | implementation / exploratory | post-result notes: outcome of the A4 definitions (PC1 pass on both, H1 11/11), FDR flag boundary, Milo worker cap, exploratory descriptive tables, what is committed | recorded for transparency; nothing changed |
