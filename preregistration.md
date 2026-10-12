# Preregistration: "Broken Bricks" (primary) and "Where the field looks vs. where the genes point" (literature study)

Project: The Seer's Atlas (animated demo films about Mantis). Repository: `PranavaKCode/broken-bricks` (private).
Protocol source: `protocol/06-research-project.md` (a verbatim copy of `mantis-demo-films/research/06-research-project.md`, SHA-256 `3E298DB34566B31D22700F81AEBAF4940A8810F87C59D4DBEFB070BF5989EAAE`, 516 lines).
Preregistration version: 1.0, written 2026-10-11, **frozen by the git commit that adds this file**. That commit is made before any data are downloaded or analysed. Its hash and timestamp are recorded in the report (`research/11-literature-results.md` in the mantis-demo-films folder) and are visible on GitHub.

## 0. How to read this document

1. Part A (section 1) fixes the hypotheses PC1, H1, H2, H2b, H2c, H3 and H4 of the primary study ("Broken Bricks", SEA-AD MTG). The hypothesis text is copied **verbatim** from `06-research-project.md` section 4.2 (lines 169-184). Frozen operational details come from sections 4.3, 4.4 and Appendix A of the same file; clarifications that the protocol leaves open are listed separately in section 1.3 and marked as clarifications, not as new hypotheses.
2. Part B (section 2) fixes the pre-specified metrics, definitions, expectations and decision rules of the literature-versus-genetics study. The metric text is copied **verbatim** from `06-research-project.md` section 5 (lines 358-368 and 372-385), followed by the exact operational definitions used in the code.
3. Section 3 says what happens when we need to deviate.
4. This document does not state results. Where it states an expectation, the expectation is falsifiable and the report will say whether it held.

### 0.1 What the analysts already know (disclosure)

The pilot numbers in `06-research-project.md` section 2 (and the Appendix A cohort summary) were gathered before this preregistration and informed the expectations below. In particular we already know: the amyloid share of PubMed AD papers by year (13.8% in 1990, 24.5% in 2000, 33.4% in 2010, 32.9% in 2020, 31.0% in 2025); the TREM2 PubMed counts (4 in 2012, 39 in 2013); a handful of PubTator3 co-mention counts (e.g. APP 99,555, PLCG2 717, USP6NL 75); the Open Targets 26.09 counts for MONDO_0004975 (13,006 associated targets, 2,179 with genetic evidence, 174 genes with genetic score >= 0.5, 51 with >= 0.7); and the SEA-AD MTG cohort summary (84 donors, 25 APOE4 carriers, 11 high-ADNC donors without dementia, 11 "Severely Affected Donors"). The literature study's expectations are therefore **confirmatory of the pilot, not blind**. The Broken Bricks hypotheses (H1 to H4) have not been tested on any data by anyone on the team; PC1 is a replication of published findings (Gabitto et al. 2024).

No data were downloaded and no analysis was run by the analysis agent before this commit, apart from reading the planning documents. File-format-dependent implementation details that cannot be known without opening the files (for example exact column order of `gene2pubtator3.gz`) are handled by the amendment procedure in section 3.

---

## 1. Part A. Broken Bricks (SEA-AD MTG): hypotheses

### 1.1 Hypotheses (verbatim from `06-research-project.md` section 4.2)

- **PC1, positive control.** CNA, with phenotype CPS_Local and covariates age, sex, PMI and method fraction:
  - finds a significant global association (P<0.05) for Microglia-PVM and for Sst;
  - shows Sst neighborhoods with negative `ncorr` (depleted);
  - shows supertype-level signs agreeing with SEA-AD's scCODA results for ≥70% of supertypes that SEA-AD calls credible.
- **H1, red/blue map.** ≥4 of the 10 analysed subclasses show a global CNA association with CPS_Local at P<0.05.
- **H2, broken bricks (headline).** After removing what overall pathology (CPS) predicts, donors' cell-type "sickness" residuals:
  - have reliability (split-half r) ≥0.6 in ≥4 subclasses; and
  - show structure beyond noise: ≥2 factors by parallel analysis against a permutation null, **or** ≥10% of donors discordant (glia residual and neuron residual of opposite sign, both |z|>1.5) versus ≤5% expected under the reliability-matched null.
- **H2b, co-pathologies.** Discordance is partly explained by independent "bricks" measured in the metadata (LATE/TDP-43, Lewy body, arteriolosclerosis/microinfarcts, CAA), or by Aβ-vs-pTau CPS components.
- **H2c, resilience.** Among High-ADNC donors, the 11 without dementia have lower (healthier) neuronal brick residuals than the 31 with dementia. Report Hedges' g with a 95% CI.
- **Brick-score rule.** A donor gets a brick score for a subclass only if that donor has ≥50 nuclei of it. Run all analyses with and without the 11 severely affected donors.
- **H3, genetics anchor.**
  - Positive control: microglia-specific genes are enriched for AD GWAS signal (MAGMA gene-property, Bellenguez 2022), P<1e-4.
  - Open question: is the pathology-red microglial program enriched for AD GWAS genes *beyond* microglial specificity (conditional MAGMA)?
  - Neuronal red/blue programs are predicted **not** to be GWAS-enriched: "damage without genetic fingerprints".
- **H4, APOE4–oligodendrocyte cholesterol.** One-sided test, APOE4 carriers > non-carriers on the cholesterol-biosynthesis module score, adjusted for CPS, age, sex, PMI. Report as a replication attempt regardless of outcome.

### 1.2 Frozen operational details (from sections 4.3, 4.4 and Appendix A of the protocol)

- **Cohort.** SEA-AD middle temporal gyrus (MTG), the 84 aged donors with quantitative neuropathology. The 5 young reference donors (29-68-year-old stages in the CELLxGENE metadata) are excluded from all CPS analyses. Raw data stay outside git.
- **Sources.** SEA-AD CELLxGENE subclass `.h5ad` files (collection `1ca90a2d-2943-483d-b678-b809bf464c30`), SEA-AD AWS Open Data (CPS: `Global_and_Local_CPS.20260105.csv`; metadata; pseudobulk objects; `pertpy_compositional_modeling` scCODA results), Bellenguez et al. 2022 GWAS (GCST90027158), Open Targets 26.09, Agora biodomains. The ROSMAP/Mathys data (Data Use Certificate) are **not** used in any form.
- **PC1, H1 (CNA).** Covarying Neighborhood Analysis, sample = donor, phenotype = `CPS_Local` (MTG), covariates = age, sex, PMI and the fraction of the donor's nuclei from the multiome method. kNN graph on `X_scVI`, `n_neighbors=30`. Outputs: per-cell `ncorr`, FDR flags at 5% and 10%, global P. Secondary phenotypes: `CPS_Local_ABeta` and `CPS_Local_pTau`. Cross-check: Milo (pertpy) on Microglia.
- **H2, H2b, H2c (brick scores).** Non-circular score: phenotype-free neighborhood-abundance-matrix PCs (k = 10-20), 5-fold-by-donor ridge regression of `CPS_Local`, out-of-fold prediction = brick score per donor per subclass (matrix B, 84 donors x subclasses). Residuals R = B - f(CPS) with a per-column LOESS on true CPS. Reliability by random split of each donor's nuclei into halves. Keep only subclasses with split-half r >= 0.6 and out-of-fold rho with true CPS > 0.4 (permutation P < 0.01). Parallel analysis against 1,000 within-column permutations. Discordance is evaluated against a reliability-matched null. Residual axes that track technical covariates (nuclei count, RIN, PMI, NeuN+ fraction, method, sequencing saturation) are dropped from the claims. Known-outlier check: SEA-AD "Severely Affected Donors" should emerge as outliers.
- **Brick-score rule and exclusions.** A donor gets a brick score for a subclass only if that donor has >= 50 nuclei of it. All analyses are run with and without the 11 Severely Affected Donors. The reference donors are excluded.
- **H3 (MAGMA).** MAGMA v1.10 with the 1000G EUR reference (build 37 positions; rsID matching to the GRCh38 summary statistics), window 10 kb up/down, N per SNP = n_cases + n_controls (effective-N sensitivity analysis). Positive control: cell-type specificity gene-property test. Program test conditional on microglia specificity. Sensitivity: exclude chr19:44-46 Mb (APOE). Cheaper cross-check: overlap with the Open Targets genes with genetic score >= 0.5 against a permutation background of expression-matched genes. Programs = top 200 positive and top 200 negative genes by donor-blocked Spearman correlation with `ncorr`.
- **H4.** Oligodendrocyte pseudobulk, module score of 20 cholesterol-biosynthesis genes (HMGCR, HMGCS1, LDLR, SREBF2, MVD, MVK, FDFT1, SQLE, LSS, CYP51A1, DHCR7, DHCR24, INSIG1, IDI1, FDPS, SC5D, MSMO1, NSDHL, HSD17B7, TM7SF2), OLS with covariates CPS, age, sex, PMI; one-sided test APOE4 carriers > non-carriers (25 vs 59 donors); pydeseq2 on raw-count pseudobulk restricted to the module genes as a second analysis. Reported as a replication attempt regardless of outcome.
- **Validation checks.** Section 4.4 of the protocol applies unchanged: positive controls must pass before any novel claim; CNA vs Milo triangulation; no label leakage (NAM and PCs phenotype-free, scores out-of-fold); reliability-matched nulls; technical-covariate regression; severe-donor recovery; conditional MAGMA and APOE-region sensitivity; this file's commit hash cited in the write-up with every deviation listed.

### 1.3 Clarifications of points the protocol leaves open (made before any data access)

- **C1, denominator of H1.** Section 4.2 says "the 10 analysed subclasses", but section 4.3 lists 13 subclass files, of which 11 are not vascular (Microglia-PVM, Astrocyte, Oligodendrocyte, OPC, Sst, Pvalb, Vip, Lamp5, L2/3 IT, L4 IT, L5 IT). Endothelial (2,069 nuclei) and VLMC (4,328) are excluded a priori as too small for donor-level inference. The analysed set for H1 is the non-vascular subclasses for which at least 40 donors have >= 50 nuclei. The threshold count "at least 4 subclasses with global P < 0.05" is unchanged. The realised number of analysed subclasses is reported next to the result.
- **C2, "method fraction".** The fraction of the donor's nuclei in that subclass that come from `10xMulti` (vs `10Xv3.1`).
- **C3, PC1 agreement with scCODA.** Sign concordance between the supertype-mean `ncorr` (CPS_Local) and the sign of SEA-AD's credible `CPS_Local` supertype effects in `pertpy_compositional_modeling`, among supertypes SEA-AD calls credible.
- **C4, split-half reliability.** Pearson r between the two half-sample brick-score vectors, uncorrected (no Spearman-Brown), median over 20 random splits.
- **C5, "neuronal brick residuals" (H2c).** The mean residual across the neuronal subclasses that pass the reliability filter (Sst, Pvalb, Vip, Lamp5, L2/3 IT, L4 IT, L5 IT).

---

## 2. Part B. Literature versus genetics ("Where the field looks vs. where the genes point")

### 2.1 Pre-specified metrics (verbatim from `06-research-project.md` section 5)

Data and attention metric (section 5, Day 1, items 1 to 3, summarised in its own words; the attention metric is **pre-registered** there as "abstract-level, normalized gene IDs, not raw string search"):

- Query: `("Alzheimer Disease"[MeSH Terms] OR alzheimer*[Title/Abstract]) AND hasabstract` (237,059 records in the pilot), split by publication month, `efetch` XML at <= 3 req/s with no API key.
- Attention metric: PubTator3 abstract-level normalized gene IDs in AD abstracts (`gene2pubtator3.gz`, filtered to the AD PMID set and to human genes via NCBI `Homo_sapiens.gene_info.gz`).
- Genetic metric: Open Targets 26.09 `genetic_association` score for `MONDO_0004975`; GWAS Catalog v2 AD associations for first-report years.

Metrics (verbatim, section 5, Day 2, item 3):

3. Metrics:
   - **Biodomain attention share vs. genetic share.** Use fractional membership for multi-domain genes (TREM2 sits in 10). Attention share = AD abstracts mentioning domain genes / total. Genetic share = sum of OT genetic scores of domain genes / total. Show per decade.
   - **Gini / Lorenz curve of attention over genes.**
   - **"Papers per unit genetic evidence"** per gene, normalized by years since first GWAS report, so recent discoveries are not penalized. This yields the under-explored list: genetic ≥0.7, top-decile low attention.
   - **Event study.** Attention trajectories before and after a gene's first genetic report (TREM2 2013, PLCG2 2017, ABI3 2017; confirm years from GWAS Catalog). Compare against genes without genetic evidence.
4. Figures:
   - dumbbell chart of attention vs. genetics by biodomain;
   - stream graph of topic share 1975–2026;
   - Lorenz curve;
   - top-25 under-explored genetically supported genes;
   - TREM2 event plot.

Validation, honest claims and never-claim lists (verbatim, section 5):

**Validation**
- Reproduce today's pilot shares (amyloid ≈31–33% of AD abstracts since 2010; TREM2 4→39 papers 2012→2013) with the bulk pipeline.
- Hand-check 50 random gene-mention assignments for NER errors.
- Sensitivity of rankings to the attention definition (gene2pubmed vs. PubTator3 vs. tiab) and to the genetic source (OT vs. TREAT-AD genetics_score).

**Honest claims**
- "Since 2010, about a third of Alzheimer's abstracts mention amyloid."
- "Genes with the strongest human-genetic support received a small fraction of the attention." Use the measured ratio, e.g. "APP has ~1,000× more AD co-mentions than USP6NL", per the chosen metric.
- "Map shows N genetically supported genes in sparsely studied territory."

**Never claim**
- "the field was wrong";
- "these genes are the cure";
- that attention should equal genetic score (genetics is one kind of evidence; drug tractability and biology matter too).

### 2.2 Operational definitions (frozen)

**Snapshot.** Every download records its retrieval date and (where available) the source release in `results/data_provenance.csv`. Open Targets is the 26.09 release. Seed for every random operation: `20261011`.

**Corpus (D0).**
- PubMed records matching the query above on the retrieval date. Records without an abstract are not in the corpus (`hasabstract`). Duplicate PMIDs across month bins are collapsed.
- Publication year = the year of the `PubDate` element of the PubMed XML (print/issue date; for `MedlineDate`, the first four-digit year), falling back to `ArticleDate`, then `DateCompleted`. This is the closest single-valued analogue of the `[dp]` year the pilot used.
- Because the retrieval date is in 2026, **2026 is a partial year**. It is shown in the figures but flagged, and it is excluded from every rate, share-per-decade ranking and event-study window (the last complete year is 2025).
- Amyloid flag (paper level): case-insensitive match of `\bamyloid` in title or abstract. Sensitivity flag: the same plus `\bA[bB]eta\b|\bA-?beta\b|\bAβ|\bbeta-?amyloid`. "Since 2010" means 2010-2025.
- TREM2 flag for validation: case-insensitive `\bTREM-?2\b` or "triggering receptor expressed on myeloid cells 2" in title or abstract.

**Attention (D1).** For each human gene g (NCBI Gene ID with taxonomy 9606 in `Homo_sapiens.gene_info.gz`):
- `A_g` = number of distinct AD-corpus PMIDs for which `gene2pubtator3.gz` lists the Gene ID of g (a Concept ID field holding several IDs is split and each ID counted once per PMID). `A_{g,y}` = the same restricted to year y.
- Non-human Gene IDs (mouse, rat, and so on) are **not** counted in the primary metric. Sensitivity S4 (below) maps mouse and rat IDs to human genes by symbol.
- "Human gene" for the Gini population means `type_of_gene == protein-coding`; the genetic and attention tables keep all gene types.

**Genetics (D2).**
- `G_g` = the Open Targets Platform 26.09 `genetic_association` datatype score of gene g for `MONDO_0004975` taken from the GraphQL `associatedTargets` object (the overall-association table that includes descendant diseases, as in the pilot). `G_g = 0` if the datatype has no score for g.
- Strong genetic support = `G_g >= 0.7`; moderate = `G_g >= 0.5` (pilot: 51 and 174 genes).
- Sensitivity S5: direct-only scores from the `association_by_datatype_direct` bulk parquet, and the TREAT-AD `genetics_score` from Agora.
- `y0_g` = first-report year of g = the earliest publication year of any GWAS Catalog v2 association for `MONDO_0004975` whose mapped genes include g (year from the PMID's PubMed record). If the GWAS Catalog has no such association for g, `y0_g` is set to 1975, the first year of the study window (flagged `y0_source = no_gwas` in the table).

**Biodomain shares (D3).**
- Biodomains come from Agora (`/api/v1/genes/<ENSG>`), requested for the gene universe U = genes with `G_g > 0` or `A_g >= 20`. Open Targets Ensembl IDs are mapped to NCBI Gene IDs through `Homo_sapiens.gene_info.gz` cross-references.
- A gene with n_g >= 1 biodomains contributes weight `1/n_g` to each of its domains (fractional membership; TREM2 has 10 domains). Genes in U without an Agora biodomain form the category "Unassigned", which is reported but excluded from the shares of named domains.
- Attention share of domain d = `sum_{g in d} A_g / n_g` divided by `sum_{g in U, assigned} A_g`. Genetic share of domain d = `sum_{g in d} G_g / n_g` divided by `sum_{g in U, assigned} G_g`. Gap_d = genetic share - attention share (percentage points).
- A paper-level companion statistic (the fraction of AD abstracts mentioning at least one gene of domain d; not fractional, so it does not sum to 1) is reported alongside.
- "Per decade": attention shares are recomputed for 1975-1989, 1990-1999, 2000-2009, 2010-2019 and 2020-2025. The genetic share is the current (26.09) share in every decade; genetic knowledge at the time is not available from this release.
- Uncertainty: 1,000 bootstrap resamples of genes (percentile 95% CI of Gap_d).

**Concentration (D4).**
- Primary: Gini coefficient and Lorenz curve of `A_g` over all human protein-coding genes (zeros included). Secondary: the same over genes with `A_g >= 1`. For comparison: Gini and Lorenz of `G_g` over the same population. Computed for the whole period and for each decade window above. Bootstrap CI by resampling PMIDs (1,000 resamples).

**Papers per unit genetic evidence (D5).**
- `Y_g = max(1, 2025 - y0_g + 1)` (years since first report, to the last complete year).
- `Asince_g` = number of AD abstracts (2025 or earlier) mentioning g published in year >= `y0_g`.
- `E_g = (Asince_g / Y_g) / G_g`, defined for `G_g > 0`: AD abstracts per year since first report, per unit of genetic score. Lower = more under-explored.
- **Under-explored list.** Candidate set S = {g : `G_g >= 0.5`}. Figure and table: the 25 genes of S with the lowest `E_g` among those with `G_g >= 0.7`, plus the full ranked list in the CSV. **Sparsely studied** = `G_g >= 0.7` and total `A_g < 100`. The headline count is the number of such genes (expected to be reported as "N genetically supported genes in sparsely studied territory").
- Contrast statistic: ratio of `A_g` of the most-studied gene to that of the least-studied strongly supported gene, under the primary metric, reported with the gene names (the pilot's "APP vs USP6NL" style).

**Event study (D6).**
- Genes: TREM2, PLCG2, ABI3. **Event year** (primary, taken from section 5 of the protocol): TREM2 2013, PLCG2 2017, ABI3 2017. The GWAS Catalog `y0_g` is reported as a cross-check; where it disagrees, the disagreement is reported and the primary analysis keeps the protocol years (the GWAS-year version is a sensitivity analysis).
- Outcome: `log(1 + A_{g,y})`, relative time t = year - event year, window t = -10 to +10 truncated at 2025.
- Controls: for each treated gene, the 20 genes with `G_g = 0` (or missing) whose pre-event mean `A` (t = -3 to -1) is closest in `log(1 + A)`, with replacement allowed across treated genes.
- Estimate: difference-in-differences of mean `log(1 + A)` over t = +1 to +5 minus t = -3 to -1, treated minus matched-control mean. 95% CI by bootstrapping the controls (1,000 resamples). Raw counts are also tabulated for each treated gene.

**Embedding and map (exploratory, descriptive; no hypothesis attached).**
- Model: `NeuML/pubmedbert-base-embeddings`, fp16, `max_seq_length = 256`, over `title + " " + abstract` (papers) and over UniProt "Function" text (reviewed human proteins), in one shared space. PCA-64 on the union, then UMAP 2-d (cosine, `n_neighbors = 30`, `min_dist = 0.05`), then Leiden clusters on the kNN graph of the PCA-64 space with TF-IDF keyword labels. No LLM calls.
- The literature floor for Mantis is capped at 150,000 points. If the corpus exceeds the cap, the floor is a sample stratified by publication year and amyloid flag, proportional allocation, fixed seed, documented in `results/floor_sampling.csv`. **All metrics use the full corpus**, never the sample.
- Stream graph: share of papers per year (1975-2026) of (a) Leiden clusters, and (b) the keyword topics amyloid, tau, microglia/neuroinflammation, lipid/cholesterol, endosome/lysosome, myelin/oligodendrocyte (the topics used in the pilot).

### 2.3 Pre-specified expectations and decision rules

These are the statements the report will be scored against. Each is reported as held / did not hold, with the number.

- **L-V1 (pilot reproduction, amyloid).** The amyloid share of corpus papers is within +/- 3.0 percentage points of the pilot value in each of 1990 (13.8%), 2000 (24.5%), 2010 (33.4%), 2020 (32.9%) and 2025 (31.0%), and the pooled 2010-2025 share is between 30% and 35%.
- **L-V2 (pilot reproduction, TREM2).** The TREM2-flag count is between 2 and 8 in 2012 and between 30 and 50 in 2013, with a 2013/2012 ratio of at least 5. (The pilot used PubMed field-tag searches over all AD records; our corpus is restricted to records with abstracts and uses a regular expression, so exact equality is not expected.)
- **L1 (concentration).** Gini of attention exceeds Gini of genetic evidence over the same gene population, for the whole period.
- **L2 (weak coupling).** Spearman correlation between `G_g` and `A_g` among genes with `G_g > 0` is positive and below 0.5 (the 32-gene pilot gave 0.41).
- **L3 (gap, direction by domain).** At least one named domain has Gap_d >= +2 percentage points with a bootstrap CI excluding 0 (genetic share exceeds attention share; expected among the Agora immune, endolysosome and lipid-metabolism domains), and at least one has Gap_d <= -2 points with a CI excluding 0 (attention exceeds genetic share; expected among APP metabolism and tau homeostasis). Domain names are the Agora names.
- **L4 (sparsely studied genes).** At least 10 genes have `G_g >= 0.7` and total `A_g < 100`.
- **L5 (event study).** For TREM2 the mean annual count over t = +1 to +5 is at least 5 times that over t = -3 to -1, and the difference-in-differences estimate against matched controls is positive with a 95% CI excluding 0. PLCG2 and ABI3: direction, size and CI are reported without a pass threshold.

### 2.4 Validation and sensitivity analyses (frozen list)

- **V1** Reproduce the pilot amyloid shares and TREM2 counts (L-V1, L-V2) with the bulk pipeline.
- **V2** Hand-check 50 random (PMID, gene) assignments for named-entity errors: 25 sampled uniformly from pairs of genes with total `A_g < 100` and 25 from genes with `A_g >= 100`. The audit reads the mention strings and the abstract text and records correct / incorrect / unclear in `results/ner_audit_50.csv`. The audit is performed by the analysis agent (not a human domain expert); precision is reported with a Wilson 95% CI, and the CSV is left in a form a human can re-judge.
- **S1** Attention definition: PubTator3 (primary) vs NCBI `gene2pubmed` (restricted to the AD PMIDs) vs title/abstract symbol search. Report Spearman correlation of `A_g` across definitions, and the Jaccard overlap of the top-25 under-explored lists.
- **S2** Exclude 2026 (already excluded from rates) and alternative year definition (earliest of print and electronic year).
- **S3** Gene universe: all genes with `G_g > 0` vs the primary universe U.
- **S4** Ortholog-by-symbol merge: add counts from mouse and rat Gene IDs whose upper-cased symbol equals a human symbol (only if the corresponding `gene_info` files are downloaded).
- **S5** Genetic source: Open Targets indirect (primary) vs direct-only vs TREAT-AD `genetics_score`.
- **S6** Under-explored threshold sweep: `G_g >= 0.5, 0.6, 0.7, 0.8` and `A_g < 50, 100, 200`.
- **S7** Amyloid-pathway reference set. Fixed list: APP, PSEN1, PSEN2, BACE1, BACE2, ADAM10, APH1A, APH1B, NCSTN, PSENEN. Sensitivity: genes in the Agora "APP Metabolism" biodomain.

### 2.5 Reporting rules

- Every number that appears in the report, a figure or the film comes from a cell in `results/*.csv`; the CSV name is quoted next to the number in the report.
- All pre-specified analyses are reported whether or not they support the expectations. A failed expectation is stated as failed.
- Exploratory analyses (the map, clusters, the stream graph, anything not listed above) are labelled exploratory.
- Claims language is limited to the "Honest claims" in section 2.1, with measured numbers. The "Never claim" list in section 2.1 applies to the report and to any film copy derived from it. Agent-written text is labelled "hypotheses generated by Mantis agents, not tested".

---

## 3. Deviations and amendments

- Any change to a threshold, covariate, exclusion, definition or the list of analyses after this commit is a **deviation**. Each is recorded in `amendments.md` with the date, the reason, what was known at the time, and the commit hash, and is listed in the report.
- File-format details that cannot be known before opening a file (column order, delimiter, field encoding, API pagination limits, renamed fields, an API that has been retired) are **amendments of kind "implementation"**: they must not change any threshold or decision rule. They are logged in `amendments.md`.
- If a pre-specified analysis cannot be run (a source is unavailable, a quota is hit), the report says so and does not substitute silently.
