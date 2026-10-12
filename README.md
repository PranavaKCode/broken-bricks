# broken-bricks

Open-data analyses for the Mantis research film in "The Seer's Atlas".

- `preregistration.md`: hypotheses PC1, H1, H2, H2b, H2c, H3, H4 (Broken Bricks, SEA-AD) and the pre-specified literature-versus-genetics metrics. Frozen by its git commit, made before any data were downloaded or analysed.
- `amendments.md`: every deviation or implementation amendment after the freeze.
- `protocol/06-research-project.md`: verbatim copy of the source protocol.
- `lit/`: code for the literature versus genetics study (runs on MIT ORCD).
- `results/`: every number reported, as CSV.
- `figures/`: figures.
- `floors/`, `data/`: not committed (large, regenerable); see `.gitignore`.

Nothing here is a claim about individual people; all cohort data are public and de-identified.

Day 2-3 (donor brick profiles, genetics, APOE4):
- `bricks/`: Day-2 code (`day2_nam.py` and `day2_bricks.py` run on ORCD; `make_floor.py`, `figures_day2.py`, `claims_sheet.py` run locally). Results in `results/bricks/` (`all/` = preregistered variant, `nosad/` = without the Severely Affected Donors, `donor_bricks.csv`, `donors_floor.csv`), figures in `figures/bricks/`.
- `genetics/`: Day-3 code (programs, MAGMA, Open Targets overlap, H4). Results in `results/genetics/`, figures in `figures/genetics/`, `programs.gmt`.
- `results/claims_sheet.csv` / `.md`: the claims sheet (every number read from a results cell).
- Report: `mantis-demo-films/research/15-bricks-results.md`. Amendments A6-A10 in `amendments.md`.

Verification follow-ups (2026-10-12): `final_checks/` (code, see its README) and `results/final_checks/`; amendment A11. Report: `mantis-demo-films/research/16-verification.md` section 13.

## Data sources and attribution

All inputs are public. Derived results in `results/` and `figures/` carry the terms of their sources.

- **SEA-AD** (Seattle Alzheimer's Disease Brain Cell Atlas): Allen Institute for Brain Science, University of Washington and Kaiser Permanente Washington Health Research Institute. Gabitto et al., *Nature Neuroscience* 2024. Single-nucleus data via CZ CELLxGENE and the SEA-AD S3 bucket. Licensed CC BY 4.0.
- **AD GWAS**: Bellenguez et al., *Nature Genetics* 2022, harmonised summary statistics GCST90027158 from the NHGRI-EBI GWAS Catalog.
- **Open Targets Platform** release 26.09 (genetic association scores for MONDO_0004975).
- **PubMed** abstracts and **PubTator3** gene annotations (NCBI/NLM); **NCBI Gene** `gene_info`.
- **UniProt** reviewed human function text (CC BY 4.0).
- **Agora** (AD Knowledge Portal) biodomain annotations.
- **MAGMA** v1.10 (de Leeuw et al., *PLoS Comput Biol* 2015) with the 1000 Genomes EUR reference panel.

## Licence

Code in this repository is released under the MIT licence (see `LICENSE`). Data and derived results remain subject to the licences and terms of the sources above.
