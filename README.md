# broken-bricks

Open-data analyses for the Mantis research film in "The Seer's Atlas".

- `preregistration.md`: hypotheses PC1, H1, H2, H2b, H2c, H3, H4 (Broken Bricks, SEA-AD) and the pre-specified literature-versus-genetics metrics. Frozen by its git commit `ef5646b` (2026-10-11 22:39:30 EDT), made before the study's data downloads and analyses began. Exactly what the team had looked at by then is stated in "Timing of the preregistration" below.
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

Release tables (2026-10-14), each computed from the tables above by a script in this repository:
- `results/genetics/leads_microglia_red.csv`: genes with an Open Targets genetic score of at least 0.7, fewer than 100 Alzheimer's abstracts (abstract-level count) and membership in the 200-gene microglia red program. The filter is written in the file's header comment (read it with `pandas.read_csv(path, comment="#")`). It returns two genes, CTSH and PILRA, which are leads worth testing; no experiment on them is reported here. Script: `genetics/leads_microglia_red.py`.
- `results/provenance.csv`: the preregistration commit hash and time, the later commit times, the cluster job times around the preregistration commit, what had been looked at before it, and the data sources and versions. Script: `final_checks/make_provenance.py`.
- `results/agent_hypotheses.csv`: where the film's leads come from. They come from the Mantis map filter re-running the selection in `leads_microglia_red.csv`. No agent run is used in the film (the one recorded Composer attempt failed on a tool error and is not used).

## Timing of the preregistration

Correction, 2026-10-14. Earlier text in this README and in the header of `preregistration.md` said the preregistration was committed "before any data were downloaded or analysed" (the header: "before any data are downloaded or analysed"), and the message of commit `ef5646b` says "frozen before any data access". The record supports a narrower statement, given here. `preregistration.md` is frozen and is not edited, and the commit message cannot be edited; the preregistration's section 0.1 already discloses the pilot described below. The same correction is logged in `DEVIATIONS.md`.

What the record supports:
- **Committed before the study began.** The preregistration commit `ef5646b` was made at 2026-10-11 22:39:30 EDT (02:39:30 UTC on 2026-10-12). No study job had been submitted on the cluster before it. The first study job was submitted 91 seconds later (22:41:01), the first SEA-AD download job at 22:43:59, and the first single-nucleus analysis (CNA) job started at 22:53:34. No SEA-AD single-nucleus expression data had been downloaded, and none of the preregistered analyses had been run on the study data, when the preregistration was committed; the only numbers in hand were the pilot figures described next.
- **What was known first.** During planning on 2026-10-11 (`protocol/06-research-project.md`, written before the commit) the team had gathered pilot numbers from free API calls (PubMed counts, TREM2 counts, PubTator3 co-mention counts for a few genes and a rough correlation with genetic score among 32 top-genetic genes, Open Targets counts), listed the SEA-AD bucket and read its README and a CSV header, and made one streaming pass over the SEA-AD MTG metadata CSV to summarise the cohort per donor (84 donors, APOE4 carriers, co-pathology counts, the 11 Severely Affected Donors, nuclei per subclass), together with the median and range of the CPS scores. The pilot informed the expectations of the literature study, so those expectations confirm the pilot rather than being blind (`preregistration.md` section 0.1). The Broken Bricks hypotheses H1 to H4 had not been tested on any data.
- **Clocks.** The commit time is recorded by git. The independent evidence is GitHub (repository created 02:39:36 UTC, first push 02:39:42 UTC) and the cluster's Slurm records (`final_checks/sacct_research_jobs.csv`). The repository was private at the time and is public now.

All of these times, with their sources, are in `results/provenance.csv`.

## Notes on the named-entity audit

The 50-pair named-entity audit (`results/ner_audit_50.csv`, summary in `results/ner_audit_summary.csv`; called a hand-check in the protocol text) was judged by the analysis agent, an AI model (see `preregistration.md` V2 and `lit/26_ner_audit_judgments.py`). It has not been re-judged by a person. The CSV is left editable so that anyone can re-judge it.

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
