# cna/: Day-1 analysis code (SEA-AD MTG, CNA red/blue maps)

Everything here ran on MIT ORCD (conda env `~/orcd/pool/mantis-films/opt/conda/sc`: cna 0.2.6, scanpy 1.12.4, pertpy 1.4.0, anndata 0.13.4) except `figures.py` and `descriptives.py`, which run on the pulled results with plain pandas/pyarrow/matplotlib. Governing documents: `preregistration.md` and amendment A4 in `amendments.md` (index in `DEVIATIONS.md`).

| script | what | where it ran |
|---|---|---|
| `common.py` | paths, the 11 subclass labels, k = 30, Nnull = 10,000, seed 0 | imported |
| `run_cna.py LABEL VARIANT` | one CNA run: reads `obs` + `X_scVI` + `X_umap` (backed), drops the 5 reference donors, joins CPS from the donor table built from `Global_and_Local_CPS.20260105.csv`, kNN, `cna.tl.association` for CPS_Local, CPS_Local_ABeta, CPS_Local_pTau; writes `<label>__<variant>.json` and `_cells.parquet`. Variants: `all` (primary), `nosad` (without the 11 Severely Affected Donors), `ge50` (donors with >= 50 nuclei; CPS_Local only) | sbatch, `mit_quicktest`, 2 cores, 12-32 GB, 1-6 min each |
| `submit_day1.sh` | throttled submitter (the quicktest QOS allows about 8 queued jobs per user) | login node (sleep + sbatch) |
| `pc1_h1.py` | assembles `h1_global_p.csv`, `cna_supertype_means.csv`, `pc1_sccoda_concordance.csv`, `pc1_summary.json` | sbatch |
| `milo_microglia.py` | Milo (pertpy, pydeseq2 solver) cross-check on Microglia-PVM; needs `LOKY_MAX_CPU_COUNT` set or pydeseq2/joblib starts one worker per node core and is OOM-killed | sbatch, 4 cores, 32 GB |
| `export_floor.py` | the Mantis cells floor (`floor1_cells.h5ad`) | sbatch |
| `check_embeddings.py` | shows that `X_scVI` is a per-subclass latent | login node (reads obsm only) |
| `descriptives.py` | exploratory tables (Abeta vs pTau map agreement, top supertypes, no-SAD agreement) | local |
| `figures.py` | all figures in `figures/cna/` (light and dark variants) | local |

Outputs: `results/cna/` (JSON, CSV, per-cell parquet; only the 11 primary `*__all_cells.parquet` are committed), `figures/cna/`.

Reproduce one subclass: `sbatch -A mit_general -p mit_quicktest -t 00:15:00 --mem=16G -c 2 cna/run_py.sbatch cna/run_cna.py sst all`.
