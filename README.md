# broken-bricks

Open-data analyses for the Mantis research film in "The Seer's Atlas".

- `preregistration.md`: hypotheses PC1, H1, H2, H2b, H2c, H3, H4 (Broken Bricks, SEA-AD) and the pre-specified literature-versus-genetics metrics. Frozen by its git commit, made before any data were downloaded or analysed.
- `amendments.md`: every deviation or implementation amendment after the freeze.
- `protocol/06-research-project.md`: verbatim copy of the source protocol.
- `lit/`: code for the literature versus genetics study (runs on MIT ORCD).
- `results/`: every number reported, as CSV.
- `figures/`: figures.
- `floors/`, `data/`: not committed (large, regenerable); see `.gitignore`.

Private repository. Nothing here is a claim about individual people; all cohort data are public and de-identified.

Day 2-3 (donor brick profiles, genetics, APOE4):
- `bricks/`: Day-2 code (`day2_nam.py` and `day2_bricks.py` run on ORCD; `make_floor.py`, `figures_day2.py`, `claims_sheet.py` run locally). Results in `results/bricks/` (`all/` = preregistered variant, `nosad/` = without the Severely Affected Donors, `donor_bricks.csv`, `donors_floor.csv`), figures in `figures/bricks/`.
- `genetics/`: Day-3 code (programs, MAGMA, Open Targets overlap, H4). Results in `results/genetics/`, figures in `figures/genetics/`, `programs.gmt`.
- `results/claims_sheet.csv` / `.md`: the claims sheet (every number read from a results cell).
- Report: `mantis-demo-films/research/15-bricks-results.md`. Amendments A6-A10 in `amendments.md`.
