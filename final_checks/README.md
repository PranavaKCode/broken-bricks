# final_checks/: verification follow-ups (2026-10-12)

Code for the open items of `mantis-demo-films/research/16-verification.md` section 12, run after the Day 1 to 3 results. Everything is exploratory or a recomputation; none of it changes a preregistered rule. Logged as amendment A11 in `amendments.md`. Results are in `results/final_checks/`. ORCD jobs ran from `~/orcd/pool/mantis-films/research/broken-bricks/final_checks/` with `cna/run_py.sbatch`; local scripts run with plain numpy, pandas, scipy.

| item | script | where | output |
|---|---|---|---|
| 1. `sacct` job times against the commit times | `sacct_vs_commits.py` (input `sacct_research_jobs.csv`, from `sacct -X -P`) | local | `results/final_checks/sacct_vs_commits.csv` |
| 2. CNA rerun with technical covariates (saturation, RIN, nuclei count; controls: NeuN fraction, permuted saturation; no-SAD variant) | `run_cna_tech.py`, `submit_cna_tech.sh`, `submit_cna_tech2.sh`, summary `cna_tech_summary.py` | ORCD, then local | `results/final_checks/cna_tech/` (JSON; per-cell parquet only for the Sst and Microglia runs), `cna_tech_summary.csv`, `cna_tech_map_agreement.csv`, `cna_tech_supertypes.csv` |
| 3. Share of each subclass and class among all nuclei, per donor, against CPS_Local with and without technical covariates | `subclass_counts.py` (reads the 1.15 GB final-nuclei metadata), `class_share_regress.py`, `count_source_compare.py` | ORCD, then local | `donor_subclass_counts.csv`, `class_share_regression.csv`, `class_share_regression_cellxgene.csv`, `class_share_sst_summary.txt`, `cellxgene_vs_final_nuclei_counts.csv` |
| 4. MAGMA gene-property test with a pooled "all neurons" specificity covariate | `specificity_pooled.py`, `magma_pooled.sbatch`, `parse_pooled.py` | ORCD, then local | `h3_pooled_neuron_property.csv`, `magma_pooled/<tag>/*.gsa.out`, `spec_pooled_meta.json` (the 5 MB covariate file `spec_covar_pooled.txt` is not committed) |
| 5. PubTator3 mention strings behind the 10 USP6NL abstracts | `usp6nl_mentions.py` | ORCD | `usp6nl_mentions.csv` |
| 6. Independent recomputation of split-half reliability (Sst, Microglia-PVM) | `reliability_check.py` | ORCD | `reliability_check/<label>.json` (`results/final_checks/reliability_check/sst.json`, `microglia.json`) |
