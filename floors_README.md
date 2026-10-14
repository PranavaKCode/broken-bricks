# Mantis floors (1 cells, 3 genes, 4 literature)

Files (not committed to git because of size; regenerable with `lit/35_assemble.py`). Local copies are in `floors/`; the originals are on ORCD at `~/orcd/pool/mantis-films/research/lit/out/floors/`.

| file | rows | size | md5 |
|---|---|---|---|
| `floor3_genes.csv` | 19,113 reviewed human genes | 23.7 MB | `0dad1bbe0c8700492bfaaf7f098db35f` |
| `floor4_literature_150k.csv` | 150,000 AD abstracts (stratified sample of the full set) | 205 MB | `f775fbe6b6604fb5064725fe3dfd8041` |
| `floor4_literature_full.csv` | 236,973 AD abstracts (all) | 324 MB | `085c04072e4ca82ee713645c95221a3a` |

## What "in register" means here

Both floors share one PCA-64 basis and one UMAP fit (cosine, n_neighbors 30, min_dist 0.05, seed 20261011) over the union of 236,973 abstracts and 19,113 gene texts, using `NeuML/pubmedbert-base-embeddings` (fp16, max_seq_length 256). So `coordinate1`/`coordinate2` and `v0` to `v63` are directly comparable between floors. The genes form their own island next to the paper cloud (UniProt function text reads differently from abstracts), so a gene sits near abstracts that describe similar biology, not inside the abstracts that mention it.

## Mantis ingestion settings (CSV route)

- `field_types['coordinate1_column'] = 'coordinate1'`, `field_types['coordinate2_column'] = 'coordinate2'`
- `field_types['vector_columns'] = ['v0', ..., 'v63']` (PCA-64 of the embedding; this keeps Mantis from calling the OpenAI embedding model)
- text column: `title` (floor 4 also has `abstract_snippet`, the first 500 characters; floor 3 has `text_snippet`)
- `floor_number` column: 3 for genes, 4 for literature
- Test with a 10k-row head first; ingestion speed at 150k rows is unverified.

## Columns

Floor 3 (genes): `floor_number, gene_id, title (= symbol), symbol, ensembl_id, uniprot, protein_name, text_snippet, text_source (uniprot_function or fallback_name), ot_genetic_score` (Open Targets 26.09 genetic_association for Alzheimer disease, 0 if none), `treatad_genetics_score` (Agora), `ad_paper_count` (abstract-verified PubTator3 mentions in AD abstracts, all years), `ad_paper_count_2010_2025, first_gwas_year, E_g` (abstracts per year since first report per unit genetic score; blank if no genetic evidence), `n_biodomains, biodomain` (pipe-separated Agora biodomains), `biodomain_primary` (alphabetically first, a convenience label only), `strong_G07` (genetic score >= 0.7), `sparse_G07_A100` (strong and fewer than 100 abstracts), `under_explored_rank` (rank among the 51 strong genes by E_g, 1 = least studied), `cluster, cluster_label, coordinate1, coordinate2, v0..v63`.

Floor 4 (literature): `floor_number, pmid, title, abstract_snippet, url, year` (PubMed issue-date year), `journal, genes` (human gene symbols found in the title/abstract, up to 30), `n_genes, biodomain` (dominant Agora biodomain of the genes mentioned by fractional membership, or Unassigned), `amyloid_flag` (1 if "amyloid" starts a word in title or abstract), `ot_genetic_score_max_gene` (highest genetic score among genes mentioned), `top_genetic_gene, cluster, cluster_label, macro_label, coordinate1, coordinate2, v0..v63`.

## Sampling of the 150k floor

Proportional allocation by publication year x amyloid flag, fixed seed 20261011, exactly 150,000 rows (`results/floor_sampling.csv`, `results/floor_sampling_check.csv`: amyloid share 30.964% in both the full set and the sample). All metrics in the report use the full corpus, never the sample. Use the full file when ingestion speed allows.

## Cautions

- Cluster labels are class-based TF-IDF keywords, not names given by a person or a model; some are odd. The agent should label them as keywords.
- A gene's `ad_paper_count` counts abstracts that mention it; 84% of 50 audited (PMID, gene) pairs were correct, lower for rarely mentioned genes (76% for genes with fewer than 100 abstracts). The audit was judged by an AI model and has not been re-judged by a person. See `results/ner_audit_summary.csv`.
- No donor or patient data are involved. Do not describe genes as cures or as "the answer"; see the "Never claim" list in the report.


---

# Floor 1 (cells): SEA-AD MTG nuclei with CNA red/blue maps

`floors/floor1_cells.h5ad` (21.5 MB, md5 `80179202e216069edb4e8cc92e4c8214`; not committed, regenerable with `cna/export_floor.py` from `results/cna/*__all_cells.parquet`). ORCD copy: `~/orcd/pool/mantis-films/research/seaad/results/day1/floor1_cells.h5ad`. Sampling table: `results/cna/floor1_sampling.csv` (nuclei available and selected per donor and subclass), summary `results/cna/floor1_cells_summary.json`.

- **149,995 nuclei, 84 aged donors, 7 subclasses, 62 supertypes.** All 38,905 aged Microglia-PVM nuclei plus 18,515 each of Astrocyte, Oligodendrocyte, OPC, Sst, Pvalb and L2/3 IT, allocated across donors in proportion to each donor's nuclei (largest remainder; seed 20261011). The protocol's "20k each" would give 158,905, so the 150k cap was treated as binding (amendment A4.7). The 5 young reference donors are not in it.
- `X` is a one-column all-zero placeholder (`placeholder`) until the Day-3 program genes exist; re-run the export with program genes then.
- `obsm['X_umap']`: the SEA-AD native per-subclass UMAPs, each scaled into a 20 x 20 box and placed on a 4 x 2 grid (Microglia, Astrocyte, Oligodendrocyte, OPC on the top row; Sst, Pvalb, L2/3 IT below); `obsm['X_umap_subclass']` keeps the native coordinates.
- `obsm['X_pca']` = `X_scVI`, 10 dimensions. **Caution:** the scVI latent in the CELLxGENE subclass files is fitted per subclass (each axis about N(0,1) in every file; a nearest-centroid classifier on it assigns cells to their subclass at 0.01-0.20 accuracy against 0.09 expected by chance, `results/cna/check_embeddings.json`). Distances in `X_pca` are therefore meaningful only within a subclass; do not run a joint kNN or UMAP over the whole floor and read it as biology.
- `obs`: `subclass, donor_id, supertype, leiden` (integer code per subclass+supertype pair, 0 to 61; mapping in `uns['leiden_codes']`), `CPS, CPS_ABeta, CPS_pTau` (donor-level CPS_Local and its two components, from the 2026 CSV), `ncorr_CPS, ncorr_ABeta, ncorr_pTau` (CNA neighborhood correlation; positive = neighborhood expands with pathology, negative = shrinks), `fdr_*` (per-cell empirical FDR), `redblue_fdr10_*` and `redblue_fdr5_*` (+1 red, -1 blue, 0 not significant), `apoe, apoe4_carrier, cognitive_status, dementia, adnc, braak, age_at_death, sex, severe_donor` (one of the 11 Severely Affected Donors).
- No NaN in any column. `ncorr_*` come from the primary all-donor run (all 84 donors, kNN 30 on X_scVI, covariates age, sex, PMI, method fraction).
- Mantis's h5ad adapter requirements have not been tested with this file (the Mantis codebase document does not specify them); run a 10k-cell head first.
- ncorr is computed over each subclass's own neighborhood graph, so values are comparable in sign and rough size across subclasses but not a common scale.


---

# Floor 2 (donors): SEA-AD MTG donors with brick scores

`floors/floor2_donors.csv` (84 rows, 198 KB; not committed; regenerated by `bricks/make_floor.py` from `results/bricks/all/donor_bricks.csv`, a copy of the same table is `results/bricks/donors_floor.csv`). One row per aged donor. Written by the Day 2-3 analysis (`research/15-bricks-results.md`).

- Ingestion (CSV route): `coordinate1`, `coordinate2` (first two principal components of the residual z matrix of the 11 subclass columns, 59.5% and 8.4% of its variance; positive PC1 = sicker than CPS predicts across cell types), `vector_columns` = `v0` to `v10` (the 11 brick scores B, z-scored over donors; 12 missing scores, from donors with < 50 nuclei of a type, replaced by the column median), `title` = donor ID, `floor_number` = 2. `umap1`, `umap2` are an alternative layout (n_neighbors 10, min_dist 0.3, seed 20261011).
- Columns: `CPS_Local`, `CPS_Local_ABeta`, `CPS_Local_pTau`, `cognitive_status`, `adnc`, `severe_donor`, `resilient_high_adnc`, `demented_high_adnc`, `apoe_genotype`, `apoe4_carrier`, age, sex, PMI, co-pathology scores (`LATE`, `Lewy`, `arteriolosclerosis`, `microinfarcts_log1p`, `CAA`, `ABeta_minus_pTau`), technical covariates (`tech_*`), per subclass `n_nuclei_<label>`, `B_<label>`, `R_<label>`, `Z_<label>`, `Ztech_<label>`, `Zadj_<label>`, tier composites and labels (`zG__*`, `zN__*`, `delta__*`, `discordant__*`, `discordance_label__*`, `factor1..3__*`), `profile_label`.
- **Cautions.** (1) The preregistered test of discordance found no donor with opposite glia and neuron residuals, so `discordance_label` is "not discordant" for every donor: do not build a "healthy microglia, sick neurons" story on this floor. (2) `profile_label` is an exploratory description of one shared residual axis, which tracks sequencing saturation (a technical covariate); it is not a claim about a donor. (3) Residuals are relative to the CPS-predicted value, so "healthier than CPS predicts" is not "healthy". (4) Use donor IDs only; never names or biographies.
