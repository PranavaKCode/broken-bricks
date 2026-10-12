# 06 — The real research project behind the Mantis research-demo film

Planning-stage report. Investigated 2026-10-11. Everything marked **verified** was checked live today against the source (API call, FTP/S3 listing, or HTTP HEAD). Anything not checked is marked **UNVERIFIED**. No OpenRouter credit was spent. Nothing was published. No git state was touched. The only local writes were tiny JSON/py files in the session scratchpad (`...\scratchpad\rp\`).

---

## 0. Decision in one screen

| | Primary: **"Broken Bricks"** | Backup: **"Where the field looks vs. where the genes point"** |
|---|---|---|
| Question | Does Alzheimer's break the same "bricks" (cell-type programs) in every brain, or do patients with the same pathology carry different broken bricks? And which bricks carry the *genetic* (causal) signal and which only the *damage* (correlational)? | How far apart are 50 years of AD research attention and today's human-genetic evidence? Which genetically supported genes and modules have barely been studied? |
| Data (all open) | SEA-AD MTG snRNA-seq, 84 donors with quantitative neuropathology (CELLxGENE + AWS Open Data), Bellenguez 2022 GWAS, Open Targets 26.09, Agora biodomains | 237,059 PubMed AD abstracts, PubTator3 gene mentions, Open Targets 26.09 genetic scores, GWAS Catalog, Agora biodomains |
| Laptop cost | ~20–30 GB download, 3–5 working days, CPU-bound (minutes per step), GPU optional | ~2–3 GB download, 2–3 days, GPU used for one ~30–60 min embedding pass |
| Thesis fit (1–5) | **5** (modules/bricks, data-first, correlation→causation via genetics, representation learning, personalized medicine) | **4.5** (serendipity, "AI's clear-simple-wrong answer = amyloid", map vs. turn-by-turn) |
| Null-result risk | Medium for the headline (donor discordance); low for positive controls | Very low (pilot numbers below already show the gap) |
| Hero shot | Red/blue disease neighborhoods on 40k microglia; then two donors with the same pathology score whose "brick profiles" light up differently | Literature floor with an "amyloid continent"; gene floor above it where genetically supported genes glow over almost-empty seas of papers; time slider shows TREM2 erupting in 2013 |

**Recommendation:** run the backup's literature pipeline first (Days 1–2; it also becomes the literature floor of the primary), then the SEA-AD analysis (Days 3–6). The film ends up with two independent, honest results in one multi-floor Mantis space.

**Blocking issue (verified at 18:33 today):** drive `C:` (952.6 GB) had **0 bytes free**. A 33 MB test download failed with `curl: (23) Failure writing output`. Free space then moved between 32 KB and 150 MB, so something else on the machine is filling the disk (the top writers by bytes since boot were chrome.exe 30 GB, MBAMService 24 GB, SearchIndexer 16 GB, and claude.exe 11 GB). **Free at least 80 GB before Day 0.** Neither project can start until this is done.

---

## 1. Constraints discovered during investigation

1. **Disk is full.** See above. Budget: primary ≈ 25 GB raw + ≈ 15 GB intermediates; backup ≈ 3 GB.
2. **The Mathys 2023 data (the "2.4M cells / 430 donors" Kellis shows) is not usable openly.** Raw data and human metadata are behind a Synapse Data Use Certificate (syn52293417, syn3157322) per https://compbio.mit.edu/ad_aging_brain/. The UCSC Cell Browser copy is open: `ad430_matrix.mtx.gz` is 20.9 GB, `meta.tsv` is 273 MB, and it holds 2,327,742 cells (verified via https://cells.ucsc.edu/rosmap-ad-aging-brain/ad-aging-brain/dataset.json). But its metadata fields are only `cellName, Batch, n_counts, percent_ribo, percent_mito, n_genes, Dataset, Major_Cell_Type, Cell_Type, Individual`. There is **no pathology or diagnosis**, so red vs. blue cannot be computed from it. The `personal.broadinstitute.org/cboix/ad427_data/` links did not respond (HTTP 000). **Use SEA-AD instead.** It is the best open substitute and has *richer* donor phenotypes.
3. **`cellxgene-census` will not run natively on Windows.** Its dependency `tiledbsoma` 2.3.0 ships no Windows wheels (verified on PyPI JSON). Use direct `.h5ad` downloads from the CELLxGENE curation API (verified below), or WSL2.
4. **Python 3.14 is the system default.** Most of the stack is pure-Python, and numba 0.68 and torch 2.14.1 have cp314 Windows wheels (verified). Still, pin **Python 3.12 with uv** to avoid edge cases: pertpy 1.4.0 requires ≥3.12, scanpy 1.12.4 requires ≥3.12, and cna 0.2.6 requires numpy≥2.
5. **Mantis ingestion paths already exist (verified in `C:\Users\goodp\OneDrive\Mantis\MantisAPI\src`):**
   - **h5ad import:** `synthesis/h5ad_views.py` and `synthesis/tasks/h5ad_tasks.py`. It uses `obsm['X_umap']` as map coordinates and `obsm['X_pca']` as the raw embedding when they are provided and selected (`include_obsm_matrices`). It runs `sc.tl.leiden` **only if `obs['leiden']` is absent**, and casts `obs['leiden']` to int. So we can supply supertype integer codes as `leiden` and get clusters equal to SEA-AD supertypes. Note that `get_adata` loads the **whole file into RAM** (`sc.read_h5ad(BytesIO)`) and each cell becomes a Postgres `Idea` row. Keep each upload small: ≤150k cells, a tiny `X` (program genes only), and `X_scVI` copied into `X_pca`.
   - **CSV import with precomputed vectors or coordinates:** `synthesis/tasks/embedding.py` honors `field_types['coordinate1_column'/'coordinate2_column']`, which skips embedding, and `field_types['vector_columns']`, a list of numeric columns used as the raw embedding. This lets us embed locally (free) and avoid the default `text-embedding-3-small` OpenAI call (`synthesis/config.py`).
   - **Floors exist in the data model:** `floor_number` on Idea and Cluster (`db/models.py:945, 1137, 1298`), and `floor_number?: number` in the frontend `src/app/types.ts:1941`.
   - Time-slider support in this build: **UNVERIFIED**. I found no `TimeSlider` component. Fallback: year-bin filters, or a HyperFrames overlay.

---

## 2. Verified data inventory

| # | Source | What | Size | Access / license | Status |
|---|---|---|---|---|---|
| D1 | **SEA-AD on CELLxGENE**, collection `1ca90a2d-2943-483d-b678-b809bf464c30` (https://cellxgene.cziscience.com/collections/1ca90a2d-2943-483d-b678-b809bf464c30) | 50 h5ad files: 24 subclasses × {MTG, DLPFC} plus 2 whole-taxonomy files. Each has `X_scVI` and `X_umap`, raw counts in `raw.X`, and 35,483 features. MTG subsets: Microglia-PVM 40,000 cells / 0.51 GB; Astrocyte 70,009 / 1.22 GB; Oligodendrocyte 111,194 / 1.43 GB; OPC 32,493 / 0.58 GB; Sst 58,265 / 1.17 GB; Pvalb 90,804 / 2.31 GB; Vip 104,514 / 2.31 GB; Lamp5 42,921 / 0.87 GB; L2/3 IT 330,085 / 10.32 GB; L4 IT 168,860 / 4.78 GB; L5 IT 128,090 / 4.26 GB; Endothelial 2,069 / 0.03 GB; VLMC 4,328 / 0.06 GB. Whole MTG: 1,378,211 cells / 50.78 GB | Public, no login. CELLxGENE hosts data under CC BY 4.0 (per CELLxGENE submission policy; **confirm on the dataset page**). Cite Gabitto et al. 2024 | verified via `https://api.cellxgene.cziscience.com/curation/v1/collections/<id>` |
| D2 | **SEA-AD AWS Open Data** bucket `s3://sea-ad-single-cell-profiling` (https://sea-ad-single-cell-profiling.s3.amazonaws.com/README.md) | 2026-06-22 "Multiregion" release covering 10 regions (MTG, DFC, MEC in all 84 donors; HIP, LEC, ITG, STG, FI, AnG, V1C in 43). Key files: `MTG/RNAseq/SEAAD_MTG_RNAseq_final-nuclei.2026-06-22.h5ad` (32.98 GB); `..._final-nuclei_metadata.2026-06-22.csv` (1.15 GB, 131 columns incl. ADNC, Thal, Braak, CERAD, CAA, LATE, Lewy, microinfarcts, atherosclerosis, arteriolosclerosis, Cognitive Status, MMSE/MOCA/CASI, **APOE Genotype**, "Severely Affected Donor", Class/Subclass/Supertype); `Multiregion_2026/model_outputs/continuous_pseudo-progression_score/Global_and_Local_CPS.20260105.csv` (632 donor×region rows: CPS_Global, CPS_Local and their **ABeta / pTau components**; MTG has 84 donors); `Multiregion_2026/pseudobulk_objects/*_pseudobulked.2026-06-22.h5ad` (0.004–0.46 GB each); `Multiregion_2026/model_outputs/pertpy_compositional_modeling/` (SEA-AD's own scCODA results, small CSVs); `.../gpboost_differential_expression_modeling/CPS_Local/*_DE_results_corrected.csv`; `Multiregion_2026/subclass_objects/SEAAD_Immune_..._final-nuclei.2026-06-22.h5ad` (3.27 GB, microglia across 10 regions) | Allen Institute Terms of Use (https://alleninstitute.org/legal/terms-use/) via https://registry.opendata.aws/allen-sea-ad-atlas/. Public, no AWS account needed | verified (S3 XML listing, README, CSV header via HTTP range request) |
| D3 | SEA-AD donor metadata xlsx | ~44 KB, one row per donor | Allen terms | portal listing found via search; **file not downloaded** (disk) |
| D4 | **Bellenguez et al. 2022 AD/ADD GWAS** (GCST90027158) | `GCST90027158_buildGRCh38.tsv.gz` 720 MB; harmonised `35379992-GCST90027158-MONDO_0004975.h.tsv.gz` 1.1 GB (has rsIDs). European, N=487,511 | EBI GWAS Catalog FTP, free (https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90027001-GCST90028000/GCST90027158/) | verified |
| D5 | **Open Targets Platform 26.09** | GraphQL `https://api.platform.opentargets.org/api/v4/graphql`. For Alzheimer disease `MONDO_0004975`: **13,006 associated targets**; datatypes: literature 11,196, rna_expression 3,126, **genetic_association 2,179**, animal_model 832, affected_pathway 761, clinical 314. Genetic score ≥0.5: **174 genes**; ≥0.7: 51. Bulk parquet (CC0): `evidence_europepmc` 5.77 GB, `association_by_datatype_direct` 0.82 GB, `l2g_prediction` 0.49 GB, `credible_set` 2.85 GB, `literature_vector` 0.03 GB, `drug_molecule` 0.01 GB (https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/26.09/output/) | CC0 | verified (live queries) |
| D6 | **GWAS Catalog REST v2** | `https://www.ebi.ac.uk/gwas/rest/api/v2/associations?efo_id=MONDO_0004975`: 6,770 AD associations with mapped genes and PMIDs (gives "year a gene was first GWAS-linked") | EMBL-EBI terms | verified |
| D7 | **PubMed** (E-utilities) | `("Alzheimer Disease"[MeSH] OR alzheimer*[tiab])`: 257,274 records; **237,059 with abstracts**; MeSH-only 144,849 | NLM terms; no key needed at ≤3 req/s | verified |
| D8 | **PubTator3** | `gene2pubtator3.gz` 723 MB, `relation2pubtator3.gz` 286 MB, `disease2pubtator3.gz` 2.0 GB (https://ftp.ncbi.nlm.nih.gov/pub/lu/PubTator3/); search API returns entity co-mention counts (abstracts + PMC full text) | NCBI public data (**UNVERIFIED license wording**) | verified |
| D9 | NCBI `gene2pubmed.gz` | 289 MB (curated gene–PMID links) | public | verified |
| D10 | **UniProt** | REST: 20,431 reviewed human entries with `cc_function` text in one TSV stream. ProtT5 per-protein human embeddings `knowledgebase/embeddings/UP000005640_9606/per-protein.h5` 340 MB | CC BY 4.0 | verified |
| D11 | **Agora (AD Knowledge Portal) API** | `https://agora.adknowledgeportal.org/api/v1/genes/<ENSG>` returns `biodomains` (e.g. TREM2 has 10 domains), `overall_scores.{target_risk_score, genetics_score, multi_omics_score}` (TREAT-AD; Cary et al. 2024), RefSeq summary; `/api/v1/genes/nominated` returns **967 nominated targets** | public web API; cite AD Knowledge Portal / Cary 2024. Bulk Synapse copy: **UNVERIFIED** | verified |
| D12 | **ChEMBL 37** (2026-05-01) | 4,225 molecules with `max_phase=4`; 306 drug_indication rows for AD (MeSH D000544) | CC BY-SA 3.0 | verified |
| D13 | **scDRS data release** (figshare 19312583) | `gs_file.zip` 0.92 MB (GWAS gene sets for ~74 traits; AD set expected as `PASS_Alzheimers_Jansen2019`, **UNVERIFIED** until unzipped) | CC BY 4.0 | file list verified |
| D14 | **MAGMA v1.10** | Windows 64-bit binary; 1000G EUR reference 488 MB; gene locations build 37/38 (https://cncr.nl/research/magma/) | free to use; no redistribution | verified |
| D15 | Replication: **PsychAD RADC cohort** (CELLxGENE `84ce6837-...`) | 693,682 DLPFC nuclei, 152 donors (ROSMAP/RADC), 25 cell types, `X_umap` only, 6.28 GB, `https://datasets.cellxgene.cziscience.com/54293783-669c-410e-919d-474960f8761b.h5ad` | CC BY 4.0 (CELLxGENE) | verified. Pathology obs columns **UNVERIFIED** |
| D16 | Replication: **Human myeloid plasticity atlas** (Nat Genet 2026, doi 10.1038/s41588-026-02716-6; CELLxGENE `f6ea32cb-...`) | Frozen in-silico-sorted myeloid nuclei: 289,279 cells, **1,469 donors**, 3.01 GB (`.../44e42474-0297-463b-92b0-4b47decf699e.h5ad`). Fresh scRNA: 543,012 cells, 137 donors, 5.25 GB. The paper reports a GPNMB phagocytic subtype rising with AD burden and correlating with AD PRS | CC BY 4.0 (CELLxGENE) | verified |
| D17 | Cross-disease: **Rexach et al. 2024 Cell** cross-dementia (CELLxGENE `c53573b2-...`) | 432,555 nuclei, 40 donors, AD / Pick's / PSP / control, 4.99 GB (`.../aa87f914-07e0-48aa-b915-8de906c95baf.h5ad`) | CC BY 4.0 (CELLxGENE) | verified |
| D18 | Embedding models (Hugging Face API) | `NeuML/pubmedbert-base-embeddings` (Apache-2.0, sentence-transformers, 768-d); `allenai/specter2_base` (Apache-2.0); `ncbi/MedCPT-Article-Encoder` (license "other"); `BAAI/bge-small-en-v1.5` (MIT); `facebook/esm2_t33_650M_UR50D` (MIT); `ctheodoris/Geneformer` (Apache-2.0) | as listed | verified |

### Pilot numbers gathered today (free API calls; used to size the risk of a null result)

- **Amyloid share of AD literature (PubMed tiab):** 13.8% in 1990, 24.5% in 2000, 33.4% in 2010, 32.9% in 2020, 31.0% in 2025. Across all years: 82,038 of 257,274 AD papers (32%) mention amyloid, 35,281 mention tau, 12,847 microglia, 12,354 lipid/cholesterol, 4,429 endosome/lysosome/endocytosis, and **2,469 myelin/oligodendrocyte (1%)**.
- **TREM2 inflection:** 1, 2, 4, **39**, 46 AD+TREM2 papers in 2008, 2010, 2012, 2013, 2014, then 145 in 2024 and 223 in 2025. 2013 was the year of the TREM2 variant papers (Guerreiro et al. and Jonsson et al., NEJM 2013).
- **Attention vs. genetics pilot** (PubTator3 search API, AD co-mentions in abstracts plus full text; Open Targets genetic score in brackets):
  - Heavily studied: APP 99,555 [0.93]; MAPT 59,413 [none]; APOE 37,041 [0.91]; PSEN1 22,421 [0.95]; ACHE 15,376 [none]; GSK3B 13,473 [none]; BACE1 11,881 [none].
  - Barely studied despite strong genetics: PLCG2 717 [0.91]; ABI3 222 [0.87]; SHARPIN 147 [0.86]; NCK2 98 [0.86]; TSPAN14 91 [0.74]; SIGLEC9 90 [0.72]; USP6NL 75 [0.75].
  - Among 32 top-genetic genes, Spearman(genetic score, co-mentions) = 0.41.
  - Caveat: tiab counts and PubTator counts differ by ~10× (PLCG2 has 61 tiab papers vs. 717 PubTator co-mentions), so **the attention metric must be pre-registered**.

---

## 3. Candidate evaluation

Scoring: "Fit" = fit with Kellis's thesis (1–5). Compute estimates are for this laptop (Core Ultra 7 155H, 32 GB RAM, RTX 3050 4 GB) and are **estimates (UNVERIFIED until Day 0 timing runs)**.

### (a) Red vs. blue neighborhoods, reproduced and extended in SEA-AD (open) — **becomes the Primary**

- **Data.** D1 (MTG subclass h5ads) and D2 (CPS per donor, co-pathologies, APOE, cognition, pseudobulks, SEA-AD's own scCODA results for cross-checking).
- **Method.** Covarying Neighborhood Analysis (CNA; Reshef et al., *Nat Biotechnol* 2022; `pip install cna`, v0.2.6, AnnData-native since 0.2.0, https://github.com/immunogenomics/cna). It tests donor-level phenotypes against cell neighborhoods built on the provided `X_scVI` graph, with covariates. It outputs a per-cell neighborhood correlation (`ncorr`), which is exactly a red/blue coloring with FDR thresholds and a global P. As a cross-check, use Milo via `pertpy.tl.Milo` (pertpy 1.4.0) on one subclass.
- **Compute.** kNN on scVI latent plus CNA takes minutes per subclass, even for 330k cells (L2/3 IT), because CNA needs only `obs` and `obsm`. Read in backed mode, so RAM stays well under 32 GB. Gene-level steps run on a ≤50k-cell subsample per subclass.
- **Positive result looks like:**
  - Significant global association with `CPS_Local` (MTG) for microglia (expansion of disease-associated states) and for Sst interneurons (early depletion). Both are published SEA-AD findings: Gabitto et al., *Nat Neurosci* 2024, doi 10.1038/s41593-024-01774-5.
  - Supertype-level directions agree with SEA-AD's own scCODA output in `Multiregion_2026/model_outputs/pertpy_compositional_modeling/`.
- **Null risk.** Low for the reproduction. The extension (donor discordance) carries medium risk; see Primary.
- **Overclaiming guardrails.** Red means "associated with pathology", not "causal". Data are cross-sectional, one region (MTG), and NeuN-sorted (70% NeuN+ / 30% NeuN−), so compositional effects are interpretable only *within* a class.
- **Hero shot.** 40k microglia in Mantis. The map starts grey, then a red/blue wave rolls across it (`ncorr`). The human lassoes the red island; the agent names the state and lights up its top genes on the gene floor.
- **Fit: 5.**

### (b) Multi-floor AD "modular map" (genes / genetics / drugs / literature) — **becomes the Mantis presentation layer**

- **Data.** D5, D7, D8, D10, D11, D12. Genes floor: UniProt function text embedded with the same model as the abstracts, so genes and papers share one space. ProtT5 (D10) is an optional "protein-family" floor, but it clusters by fold and family, not by disease module. Drugs floor: ChEMBL approved (4,225), fingerprint UMAP.
- **Compute.** Trivial apart from the abstract embedding (see backup).
- **Positive result.** Biodomain modules form visible regions; genetics glow concentrates in immune / lipid / endolysosome regions while literature density concentrates in APP/amyloid.
- **Null risk.** Low as a visualization, but it is not a "result" by itself.
- **Fit: 4.** It is the vehicle, not the finding.

### (c) FTO / IRX3 / IRX5 obesity circuitry re-analysis

- **Data.**
  - GTEx v10 open eQTLs: bulk adipose famously shows little FTO-locus eQTL signal, which is part of the original story.
  - Roadmap Epigenomics chromatin states in E025 / E023 (adipose MSC / MSC-derived adipocytes).
  - JASPAR ARID5B motif scoring of rs1421085 T>C.
  - Optional AlphaGenome `score_variant`. It is free for non-commercial use but **needs the user to register for a key**, which we cannot do (https://deepmind.google/blog/alphagenome-ai-for-better-understanding-the-genome/).
- **Compute.** Small.
- **Positive result.** Reproduce the motif disruption and adipocyte-progenitor enhancer state of Claussnitzer et al., *NEJM* 2015.
- **Null risk.** Low (it is a reproduction), but there is **no novelty and weak map visuals**: it is one locus, not a landscape.
- **Fit: 2.5.** Better as a 15-second "flip one letter of 3.2 billion" animation than as the research project.

### (d) APOE4 → cholesterol → myelination signal check — **folded into the Primary as a pre-registered side test**

- **Prior.** Blanchard et al., *Nature* 2022 (Tsai/Kellis): APOE4 oligodendrocytes up-regulate cholesterol-biosynthesis genes, with a dose effect; snRNA-seq of 32 PFC samples.
- **Open test.** SEA-AD oligodendrocyte and OPC pseudobulks per donor (D2), split by APOE genotype in the metadata. Score a cholesterol-biosynthesis module (HMGCR, HMGCS1, LDLR, SREBF2, MVD, MVK, FDFT1, SQLE, LSS, CYP51A1, DHCR7, DHCR24, INSIG1, IDI1, FDPS, SC5D, MSMO1, NSDHL, HSD17B7, TM7SF2). Covariates: CPS, age, sex, PMI, method.
- **Null risk.** High: the SEA-AD APOE4 carrier count is modest (see §4.1) and the region differs (MTG vs. PFC). Report direction, effect size and CI whatever the outcome.
- **Fit: 4** as a "replication in an independent open cohort" beat.

### (e) Literature cartography: "where the field looks vs. where genetics points" — **becomes the Backup** (and the literature floor of the Primary)

- **Data.** D5, D6, D7, D8, D11.
- **Compute.** Embedding 237k abstracts with `NeuML/pubmedbert-base-embeddings` (max_seq_length 256, fp16, batch 64) on the RTX 3050 4 GB: an estimated **20–60 min** (UNVERIFIED; CPU fallback 2–4 h). UMAP on PCA-50 takes 10–20 min on CPU.
- **Positive result.** Already visible in the pilot: attention share ≫ genetic share for APP/amyloid; the reverse for immune / lipid / endolysosomal biodomains; dozens of genes with genetic score ≥0.7 and <100 AD abstracts.
- **Null risk.** Very low. The risk is "known qualitatively", so frame it as *quantification plus a navigable map*, not as a discovery.
- **Fit: 4.5.**

### (f) NEW — Cross-cohort microglia replication at 1,469 donors

- **Data.** D16.
- **Analysis.** Run CNA with AD vs. control on the frozen myeloid nuclei, then check whether the SEA-AD "red" microglial program and the paper's GPNMB/phagocytic subtype coincide.
- **Compute.** Needs PCA + kNN on 289k cells (no scVI in the file; `X_umap` only): about 30–60 min CPU and ~10 GB RAM (estimate).
- **Value.** Strong replication with 1,469 donors.
- **Fit: 4.** Use as an optional Day-6 extension.

### (g) NEW — Bricks reused across diseases

- **Data.** D17 (Rexach 2024: AD, Pick's, PSP, control; 40 donors).
- **Analysis.** Score SEA-AD-derived programs in each disease and ask which glial programs are shared ("reusable LEGO bricks") and which are disease-specific.
- **Fit: 5** for the LEGO metaphor, but only 40 donors and 4 groups: medium-high null risk. Stretch goal only.

### Score table

| Candidate | Thesis fit | Novelty | Null risk | Laptop cost | Visual payoff in Mantis | Verdict |
|---|---|---|---|---|---|---|
| (a) SEA-AD red/blue + bricks + genetics | 5 | medium (open-cohort replication + donor discordance + genetics overlay) | medium | 25 GB, 3–5 d | very high | **Primary** |
| (b) multi-floor modular map | 4 | low | low | small | very high | presentation layer |
| (c) FTO | 2.5 | none | low | small | low | animation only |
| (d) APOE4 cholesterol | 4 | low (replication) | high | small | medium | side test in Primary |
| (e) literature vs. genetics | 4.5 | low–medium | very low | 3 GB, 2–3 d | high | **Backup** + literature floor |
| (f) 1,469-donor microglia | 4 | medium | low–medium | 3 GB, 1 d | high | optional extension |
| (g) cross-disease bricks | 5 | medium | medium–high | 5 GB, 1–2 d | high | stretch |

---

## 4. PRIMARY: "Broken Bricks" — protocol

### 4.1 Cohort facts (SEA-AD MTG)

- 84 aged donors plus 5 young neurotypical reference donors. The reference donors appear as 29–68-year-old stages in the CELLxGENE metadata and must be excluded from CPS analyses. Mean age at death 88; 33 male, 51 female (Gabitto et al. 2024, via PMC11614693).
- CPS (MTG, from D2) covers 84 donors. CPS_Local median 0.523 (range 0.020–0.809); CPS_Global median 0.645 (0.018–0.990).
- The local Aβ and pTau components correlate at r = 0.68. They are related but separable, so they act as two candidate "bricks".
- Verified today by streaming the 2026 metadata (Appendix A):
  - **25 APOE4 carriers**;
  - 42 dementia / 42 no dementia;
  - **11 "resilient" donors** with High ADNC and no dementia;
  - 11 "Severely Affected Donors", who account for most neuronal dropout (e.g. all 5 donors with <100 Pvalb nuclei);
  - microglia median 380 nuclei per donor (min 67), which is adequate for donor-level scores.

### 4.2 Hypotheses (write these into `research/preregistration.md` and commit before running anything)

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

### 4.3 Day-by-day

All paths assume a project folder `C:\Users\goodp\OneDrive\Mantis\mantis-demo-films\research\broken-bricks\`. OneDrive syncing 25 GB is a bad idea, so put the raw data outside OneDrive, e.g. `C:\mantis-research\data\` (or another drive), and keep only code, results and figures in the repo folder.

**Day 0 (½ day): unblock and set up**

1. Free ≥80 GB on `C:`.
2. Create the environment:
   ```powershell
   uv venv C:\mantis-research\.venv --python 3.12
   uv pip install --python C:\mantis-research\.venv\Scripts\python.exe "scanpy>=1.12" anndata h5py cna==0.2.6 pertpy scikit-learn statsmodels pydeseq2 polars pyarrow openpyxl umap-learn pynndescent leidenalg igraph matplotlib seaborn tqdm requests
   uv pip install --python C:\mantis-research\.venv\Scripts\python.exe torch --index-url https://download.pytorch.org/whl/cu126   # CUDA tag UNVERIFIED for this driver; check nvidia-smi
   uv pip install --python C:\mantis-research\.venv\Scripts\python.exe sentence-transformers
   ```
3. Download SEA-AD (≈29 GB for the full list; ≈18 GB without L4/L5 IT). Re-resolve URLs first, because CELLxGENE version URLs change on revision. `GET https://api.cellxgene.cziscience.com/curation/v1/collections/1ca90a2d-2943-483d-b678-b809bf464c30`, then for each dataset take `assets[filetype==H5AD].url`. URLs as of today:
   ```text
   Microglia-PVM  https://datasets.cellxgene.cziscience.com/ab53c08f-03a3-4e0e-973a-44beabd78f60.h5ad   0.51 GB
   Astrocyte      https://datasets.cellxgene.cziscience.com/53f96bc7-3f23-444c-aa4f-4c69d62906bd.h5ad   1.22 GB
   Oligodendrocyte https://datasets.cellxgene.cziscience.com/021a2e09-19c6-42eb-9e52-3ed2cf3d6032.h5ad  1.43 GB
   OPC            https://datasets.cellxgene.cziscience.com/01a8996c-8817-43c5-ac45-d66342ddb97e.h5ad   0.58 GB
   Sst            https://datasets.cellxgene.cziscience.com/c3e6574c-c9ed-4495-9ccd-44146565c8f0.h5ad   1.17 GB
   Pvalb          https://datasets.cellxgene.cziscience.com/cf14bcbb-f691-4d39-98d1-59de53c666c4.h5ad   2.31 GB
   Vip            https://datasets.cellxgene.cziscience.com/d9d217bf-9a57-440c-9b0b-885185f94240.h5ad   2.31 GB
   Lamp5          https://datasets.cellxgene.cziscience.com/ffa13314-dc7c-4651-97c6-fc099e90ae64.h5ad   0.87 GB
   L2/3 IT        https://datasets.cellxgene.cziscience.com/18a1753c-4c65-4bc8-a16c-420cf84ad435.h5ad  10.32 GB
   L4 IT          https://datasets.cellxgene.cziscience.com/67385399-fe99-4f46-b8df-aa18fa4bd82d.h5ad   4.78 GB
   L5 IT          https://datasets.cellxgene.cziscience.com/4776ed90-3953-4c30-bc24-092d258a443b.h5ad   4.26 GB
   Endothelial    https://datasets.cellxgene.cziscience.com/69d18d59-a55c-43cf-9abf-5ed7bb31b39b.h5ad   0.03 GB
   VLMC           https://datasets.cellxgene.cziscience.com/60fdbb86-8c6e-48ec-9dd3-2bf2c9e883ea.h5ad   0.06 GB
   ```
   Use `curl -L --retry 5 -C - -o <file> <url>` (resumable). Observed throughput was ~4 MB/s, i.e. ~2 h for 29 GB.
4. Download from S3 (plain HTTPS, no AWS account):
   - `https://sea-ad-single-cell-profiling.s3.amazonaws.com/Multiregion_2026/model_outputs/continuous_pseudo-progression_score/Global_and_Local_CPS.20260105.csv`
   - The pseudobulk objects for Microglia-PVM ("Immune"), Astrocyte, Oligodendrocyte, OPC, Sst, Pvalb, Vip, L23-IT (≈1.8 GB total).
   - The `pertpy_compositional_modeling/CPS_Local/*_Supertype_results.csv` files (0.19 GB) and `pertpy_summary_CPS_Local.20260622.csv`.
   - The donor-metadata xlsx from the SEA-AD portal.
5. Download from EBI: Bellenguez harmonised sumstats (1.1 GB). From the MAGMA site: the Windows binary, `g1000_eur` (488 MB) and `NCBI37.3.gene.loc`.
6. **Timing runs** (go/no-go): load the Microglia h5ad in backed mode and confirm `obs` contains donor ID, Supertype and (ideally) CPS. Whether the CELLxGENE obs carries CPS is **UNVERIFIED**; the fallback is to join `Global_and_Local_CPS` on donor ID (ID format `H20.33.045` matches in both sources, verified). Then time CNA on Microglia.
7. Write and commit `preregistration.md` (§4.2), with the exact attention metric, thresholds and covariates.

**Day 1: red/blue maps (PC1, H1)**

For each subclass:

1. Read `obs` and `obsm['X_scVI']` only.
2. Drop the reference donors.
3. Build the kNN graph: `sc.pp.neighbors(adata, use_rep='X_scVI', n_neighbors=30)`.
4. Run CNA with sample = donor, phenotype = CPS_Local (MTG), covariates = age, sex, PMI and the fraction of the donor's nuclei from the multiome method. Using batch = library prep as the CNA batch argument is optional.
5. Save per-cell `ncorr`, FDR-significant flags (5% and 10%), and the global P.

Repeat with CPS_Local_ABeta and CPS_Local_pTau separately; these are the "amyloid brick" and "tau brick" colorings. Compare supertype-mean `ncorr` against SEA-AD's scCODA credible effects (agreement = sign concordance plus Spearman). Run Milo (pertpy) on Microglia as a method cross-check.

Deliverables: `fig1_redblue_<subclass>.png` (datamapplot / matplotlib, red-blue diverging, symmetric limits); `cna_results.parquet`.

**Day 2: donor brick profiles (H2, H2b)**

The non-circular score works as follows.

1. Compute the CNA neighborhood-abundance matrix (NAM, cells × donors) **without the phenotype** and take its top-k PCs (k = 10–20).
2. For each subclass, run 5-fold-by-donor ridge regression of CPS_Local on the NAM PCs.
3. The out-of-fold prediction is the donor's "how sick do *this donor's* microglia look, judged by models trained on other donors" score. It is computed per subclass, giving the brick matrix B (84 × 10).

Then:

- **Reliability.** Randomly split each donor's cells into halves and recompute B on each half: split-half r per subclass. Keep only subclasses with r ≥0.6 and an out-of-fold prediction ρ with true CPS >0.4 (permutation P<0.01).
- **Residuals.** R = B − f(CPS), using a per-column LOESS on true CPS. Then:
  - correlate residuals across subclasses and run parallel analysis against 1,000 within-column permutations;
  - cluster donors (k-means / Leiden on R);
  - list discordant donors ("healthy microglia, sick neurons" and the reverse).
- **Co-pathologies and outcomes.** Regress residuals on LATE, Lewy, arteriolosclerosis, microinfarcts, CAA, Aβ-vs-pTau CPS difference, APOE4, and cognitive status / last CASI / MMSE (with CPS as covariate).
- **Technical checks.** Regress the same residuals on technical covariates: nuclei count, RIN, PMI, NeuN+ fraction, method, sequencing saturation. Any residual axis that tracks technical covariates is dropped from the claims.
- **Known-outlier check.** Do SEA-AD's "Severely Affected Donors" (metadata flag) emerge as outliers? This is an external validation of the method.

Deliverables: `fig2_brick_matrix.png` (heatmap donors × subclasses, ordered by CPS); `fig3_discordant_pairs.png` (two donors matched on CPS with opposite profiles); `donor_bricks.csv`.

**Day 3: what the bricks are made of, and do they carry genetic fingerprints (H3, H4)**

1. **Programs.** On ≤50k cells per subclass, compute log-normalised expression and correlate each gene with `ncorr` (Spearman, donor-blocked: average within donor first to avoid pseudo-replication). Keep the top 200 positive genes (red program) and top 200 negative genes (blue program).
2. **Annotation.** Annotate programs with Agora biodomains (D11; per-gene API for program genes only, ~2k calls at ≤5/s) and GO enrichment. Agent-generated labels are "suggested names", not results.
3. **Genetics with MAGMA (Windows):**
   ```bash
   magma --annotate window=10,10 --snp-loc g1000_eur.bim --gene-loc NCBI37.3.gene.loc --out ad
   magma --bfile g1000_eur --pval bellenguez_rsid_p_n.txt use=rsid,p ncol=N --gene-annot ad.genes.annot --out ad_genes
   # cell-type specificity (positive control): covariate file = per-gene specificity per subclass from pseudobulk
   magma --gene-results ad_genes.genes.raw --gene-covar subclass_specificity.txt --out ad_celltype
   # programs, conditional on microglia specificity
   magma --gene-results ad_genes.genes.raw --set-annot programs.txt --gene-covar subclass_specificity.txt condition=Microglia --out ad_programs
   ```
   Matching is by rsID, so the build-37 reference works with the GRCh38 sumstats. N per SNP = n_cases + n_controls (effective-N sensitivity analysis). **Sensitivity: exclude chr19:44–46 Mb (APOE).**
4. **Cheaper cross-check.** Overlap of programs with the 174 Open Targets genes with genetic score ≥0.5, using a permutation background of expression-matched genes from the same subclass.
5. **Optional.** scDRS per-cell AD scores (D13 gene set) on Microglia and Astrocyte, giving a "genetic glow" layer.
6. **H4.** Use the oligodendrocyte pseudobulk (D2) and the cholesterol module score. Fit an OLS model with the covariates on the per-donor mean. Also run pydeseq2 on raw-count pseudobulk restricted to the module genes. Report the effect, 95% CI and one-sided P.

Deliverables:
- `fig4_genetic_fingerprints.png`: MAGMA −log10 P per subclass, plus red-program conditional P.
- `fig5_apoe4_oligo_cholesterol.png`
- `programs.gmt`

**Day 4: Mantis floors and agent rehearsal**

Build one Mantis space with these floors.

| Floor | Points | Coordinates | Raw embedding | Key features (obs / CSV columns) |
|---|---|---|---|---|
| 1 Cells | ≤150k cells: all 40k microglia + 20k each of astro, oligo, OPC, Sst, L2/3 IT, Pvalb (stratified by donor) | SEA-AD `X_umap` (per subclass, offset in a grid) or a joint UMAP of `X_scVI` | `X_scVI` → stored as `X_pca` | donor, supertype (→ `leiden` int codes), CPS, `ncorr_CPS`, `ncorr_ABeta`, `ncorr_pTau`, red/blue FDR flag, APOE, cognitive status |
| 2 Donors | 84 | UMAP of the residual brick matrix R (or PCA-2) | B (10-d) | CPS, cognition, co-pathologies, discordance label, per-subclass scores |
| 3 Genes | ~20k reviewed human genes | **shared** UMAP with floor 4 (same embedding space) | PubMedBERT of UniProt function text (PCA-64) | OT genetic score, TREAT-AD scores, biodomains, AD paper count, program membership |
| 4 Literature | 237k AD abstracts (or ≤150k if ingestion is slow) | shared UMAP with floor 3 | PubMedBERT (PCA-64) | year, journal, genes mentioned (PubTator3), biodomain, amyloid flag |

Ingestion:
- Floors 1–2 go in via the h5ad route. Write a slim h5ad: `X` = program genes only, plus `obsm['X_umap']`, `obsm['X_pca']`, and `obs` columns.
- Floors 3–4 go in via CSV with `coordinate1/2` columns, so the floors stay in register, plus `vector_columns` (64 PCs).
- Run a test with 10k points first and time it, because performance at 150–250k points on local Postgres is **UNVERIFIED**.

Rehearse the agent prompts (§6) and record clean screen captures at 4K/60 if possible.

**Day 5: write-up and figures.**
- `writeup.md`, 2 pages: question, data, methods, 5 figures, what replicated, what did not, limitations.
- Fill in the claims sheet (§7) with the actual numbers.
- Put an internal review checklist at the end: every on-screen number traces to a cell in `results/*.csv`.

**Day 6 (optional): replication.**
- D16, frozen myeloid nuclei, 1,469 donors: PCA on HVGs, kNN, CNA AD vs. control. Ask whether SEA-AD's red microglial program scores higher in AD-associated neighborhoods there, and whether it overlaps the GPNMB subtype.
- D15 (PsychAD RADC, 152 donors, all cell types): if its `obs` has pathology, run a discordance replication.

### 4.4 Validation checks (summary)

1. Positive controls PC1 (Sst down, microglia up; agreement with SEA-AD scCODA) must pass. If they fail, stop and debug before any novel claim.
2. Method triangulation: CNA vs. Milo on Microglia; supertype-composition ridge vs. NAM-PC ridge for B.
3. No label leakage: NAM and PCs are phenotype-free, and all donor scores are out-of-fold.
4. Reliability-matched nulls for discordance (split-half and permutation).
5. Technical-covariate regression of all residual axes.
6. SEA-AD "Severely Affected Donors" should be recovered as outliers.
7. The genetics test separates cell-type specificity from program membership (conditional MAGMA); APOE-region exclusion sensitivity.
8. A frozen `preregistration.md` commit hash appears in the write-up. Every deviation is listed.

### 4.5 What would be a positive vs. null result, and how the film handles either

| Outcome | What we say on screen |
|---|---|
| PC1 + H1 pass (expected) | "Open data, one laptop: Mantis reproduces the published early loss of Sst interneurons and the rise of disease-associated microglia in 84 brains." |
| H2 passes | "Same pathology score, different broken bricks: in N of 84 donors, one cell type looks far sicker than their overall pathology predicts, while another looks healthier." |
| H2 null | "At 84 donors, the bricks broke together, *within the noise we can measure*. To see individual bricks you need the 1,469-donor cohorts." This is still an honest, map-driven conclusion; switch emphasis to H3 and the literature gap. |
| H3 positive control passes (expected) | "The damage is everywhere, but the genetic fingerprints are on the microglia." |
| H3 program-level positive | "And the microglial state that expands with pathology is itself enriched for AD risk genes, beyond just being microglial." |
| H2c positive | "Eleven people had brains full of plaques and tangles and no dementia. Their neurons looked [g = x] healthier than those of people with the same pathology who had dementia." Phrase it as an association, n=11. |
| H4 either way | "We re-tested the APOE4–cholesterol–myelin finding in an independent open cohort (25 APOE4 carriers vs. 59 non-carriers): [direction, effect, CI]." |

---

## 5. BACKUP: "Where the field looks vs. where the genes point" — protocol

This also builds floors 3–4 of the Primary.

**Day 1: data (≈3 GB)**

1. **PubMed.**
   - Query: `("Alzheimer Disease"[MeSH Terms] OR alzheimer*[Title/Abstract]) AND hasabstract` (237,059 records today).
   - ESearch only returns the first 10,000 PMIDs per query, so split by publication month (`YYYY/MM[dp]`) with `usehistory=y`.
   - Then `efetch` XML in batches of ≤10,000 at ≤3 req/s; no API key is needed.
   - Store `pmid, year, journal, title, abstract, mesh` in `abstracts.parquet` (~400 MB).
2. **Gene mentions.** Download `gene2pubtator3.gz` (723 MB). Filter to the AD PMID set and to human genes (NCBI `Homo_sapiens.gene_info.gz` for the ID→symbol map). The result is per-gene, per-year counts of AD abstracts mentioning the gene. **Pre-register this as the primary attention metric:** abstract-level, normalized gene IDs, not raw string search, because symbols like APP, ACE, MME, HFE and MPO are ambiguous.
3. **Genetics.**
   - Open Targets GraphQL for `MONDO_0004975`: page through all 13,006 targets (14 calls of 1,000) and keep `genetic_association`.
   - GWAS Catalog v2: 6,770 associations → mapped genes → first-report year (via PMID).
   - Bellenguez 2022 75-locus gene list (paper's table; or OT credible sets for GCST90027158).
4. **Modules.** Agora biodomains for all genes with OT genetic >0 or ≥20 AD abstracts (~3–5k API calls at ≤5/s, ≈15 min). Alternative: rebuild from GO using the Cary 2024 supplementary GO-term table (medRxiv 10.1101/2022.12.15.22283478).

**Day 2: embeddings, map, metrics**

1. `sentence-transformers` with `NeuML/pubmedbert-base-embeddings` (fp16, `max_seq_length=256`, batch 64) on title+abstract (237k) and on UniProt function text (20,431 genes), in one shared space. Save as float16 `.npy`. Optional model comparison on a 10k subset: `allenai/specter2_base` or `ncbi/MedCPT-Article-Encoder`.
2. PCA-64 on the union, then UMAP 2-d (`n_neighbors=30, min_dist=0.05`, cosine), then Leiden / HDBSCAN clusters. Cluster names: TF-IDF keywords first; LLM labels only optionally (see budget).
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

**Day 3: Mantis floors 3–4, agent rehearsal, write-up.**

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

---

## 6. Mantis staging: hero shots and agent prompts

Each shot maps to a Kellis talk beat.

1. **"Chat answer vs. map"** (*clear, simple, and wrong*). The chatbot pane answers "What causes Alzheimer's?" with "amyloid plaques". Cut to the literature floor: the amyloid continent dominates. Pull up a floor to the genes: the genetic glow sits mostly *elsewhere* (immune / lipid / endolysosome).
2. **"Cave of Wonders → carpet → F-18"** (data-first). Open the 1.38M-cell SEA-AD universe, then fly into the 40k microglia floor. No hypothesis yet: color by CPS and watch the red/blue wave appear.
3. **"Lasso the red island"** (humans + agents). The human lassoes the red neighborhood. Agent prompt: *"What is this region and which genes define it? Light up the source cells and the genes."* The answer lights cells on floor 1 and program genes on floor 3, with abstracts on floor 4 as citations (traceability).
4. **"Correlation → causation"** (genetics). Toggle the "genetic fingerprint" layer (MAGMA / scDRS). Microglia glow; neurons, though red with damage, do not. Line: "Neurons carry the damage. Microglia carry the genes."
5. **"Personal LEGO bricks"** (personalized medicine). On the donor floor, pick two donors with near-identical CPS. Their brick bars differ: microglia vs. neurons. (Only if H2 passes; otherwise show the group-level brick decomposition: Aβ brick vs. tau brick coloring.)
6. **"Don't leave serendipity to chance."** Agent prompt: *"Find genes with strong human-genetic support for Alzheimer's that sit inside the expanding microglial program but have fewer than 100 Alzheimer's abstracts."* Points light up on the gene floor, with near-empty literature directly beneath. The agent drafts *hypotheses* with citations.
7. **Time slider 1975→2026.** The amyloid continent grows; TREM2 erupts in 2013 right after the genetics papers ("the genetics rerouted the field once; here are the next 25").

---

## 7. Claims sheet (fill numbers after the run; nothing goes on screen without a source cell)

**Allowed**, with numbers from `results/`:
- Dataset facts:
  - "84 donors, middle temporal gyrus, [N] nuclei from SEA-AD, all open."
  - "237,059 Alzheimer's abstracts; 13,006 genes linked to AD in Open Targets; 2,179 with genetic evidence."
- Reproductions: "Reproduced: early Sst loss and microglial expansion (CNA P = [x])."
- Measured gaps: "[x]% of AD abstracts mention amyloid; [genes] with genetic score ≥0.7 average [y] abstracts vs. [z] for amyloid-pathway genes."
- Conditional findings: "In N of 84 donors… [if H2 passes]"; "Microglia-specific genes are enriched for AD risk (MAGMA P = [x])."
- Agent output, always labelled: "Hypotheses generated by Mantis agents, not tested."

**Forbidden**
- "Mantis discovered a cure / new drug target."
- "Proved causation."
- "Personalized diagnosis."
- "Patients with X should…"
- Any claim about individual living people (all donors are deidentified and postmortem; use donor IDs only, never names or invented biographies).
- Implying that the SEA-AD or Kellis-lab results are ours.
- Using the ROSMAP/Mathys data (DUC) in any form.

**Attribution line in the credits:** SEA-AD (Allen Institute, Kaiser Permanente Washington Health Research Institute, UW), CZ CELLxGENE, Open Targets, GWAS Catalog / EADB (Bellenguez et al. 2022), NCBI PubMed / PubTator3, UniProt, AD Knowledge Portal / Agora (TREAT-AD).

---

## 8. Budget

- **OpenRouter: $0 required.** All embeddings are local, and all data APIs are free.
- Optional: LLM cluster labels for ~300 literature clusters and ~40 cell-state clusters, plus agent rehearsal calls inside Mantis. Use the cheapest capable model and cap at **$0.25**. Exact per-token prices are owned by the cost-planning report and are **UNVERIFIED here**. Make sure the local Mantis build does not silently call `text-embedding-3-small`: always pass `vector_columns` or coordinates.
- Time: Primary 5 working days (+1 optional replication day); Backup 3 days; shared literature work counted once, for ≈6–7 days total.

## 9. Risks and mitigations

| Risk | Mitigation |
|---|---|
| Disk full (current state) | Free ≥80 GB; store raw data outside OneDrive; delete whole-taxonomy files (never download the 50 GB file) |
| CELLxGENE obs lacks CPS / APOE fields | Join D2 CSV / donor xlsx on donor ID (format verified identical) |
| L2/3 IT (10.3 GB, 330k cells) strains RAM | Backed mode; CNA needs only obsm/obs; subsample for gene steps |
| H2 null at n=84 | Pre-registered fallback narrative; Day-6 replication at 1,469 donors (microglia) / 152 donors (all types) |
| Mantis ingestion slow at 150–250k points | Time a 10k test on Day 0; cap floors at 150k; literature floor can be 2000–2026 only |
| NER/symbol ambiguity in literature counts | PubTator3 normalized IDs; sensitivity analyses; manual audit of 50 |
| Overclaiming in a funny film | Claims sheet §7; every number traced; agent outputs labelled "hypothesis" |
| Name/brand issues | Use SEA-AD/CELLxGENE names only in credits/attribution, no logos implying endorsement |

## 10. Sources

- SEA-AD: https://registry.opendata.aws/allen-sea-ad-atlas/ ; https://sea-ad-single-cell-profiling.s3.amazonaws.com/README.md ; https://sea-ad-single-cell-profiling.s3.amazonaws.com/Multiregion_2026/README.md ; Gabitto et al. 2024, *Nat Neurosci*, https://doi.org/10.1038/s41593-024-01774-5 (PMC11614693); donor metadata via https://portal.brain-map.org/explore/seattle-alzheimers-disease/seattle-alzheimers-disease-brain-cell-atlas-download
- CELLxGENE curation API: https://api.cellxgene.cziscience.com/curation/v1/collections/1ca90a2d-2943-483d-b678-b809bf464c30 (and `84ce6837-548d-4a1f-919f-0bc0d9a3952f`, `f6ea32cb-a10e-46ba-ba39-38d469c9cade`, `c53573b2-eff4-4c5e-9ad0-b24d422dfd9b`)
- Mathys et al. 2023 *Cell*, https://doi.org/10.1016/j.cell.2023.08.039 ; access page https://compbio.mit.edu/ad_aging_brain/ ; UCSC browser https://cells.ucsc.edu/?ds=rosmap-ad-aging-brain+ad-aging-brain
- Mathys et al. 2024 *Nature* (multiregion, gene modules): https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11338834/
- Green et al. 2024 *Nature* 633:634–645 (cellular communities / BEYOND; two trajectories): https://doi.org/10.1038/s41586-024-07871-6
- Blanchard et al. 2022 *Nature* 611:769–779 (APOE4, cholesterol, myelination): https://compbio.mit.edu/publications/Blanchard_Nature_22.pdf
- Claussnitzer et al. 2015 *NEJM* 373:895–907 (FTO rs1421085)
- Bellenguez et al. 2022 *Nat Genet* (75 AD/ADD loci): GWAS Catalog GCST90027158, https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90027001-GCST90028000/GCST90027158/
- Open Targets Platform 26.09: https://api.platform.opentargets.org/api/v4/graphql ; https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/26.09/output/
- GWAS Catalog REST v2: https://www.ebi.ac.uk/gwas/rest/api/v2/associations?efo_id=MONDO_0004975
- PubMed E-utilities: https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi ; PubTator3: https://ftp.ncbi.nlm.nih.gov/pub/lu/PubTator3/ and https://www.ncbi.nlm.nih.gov/research/pubtator3-api/
- UniProt REST and embeddings: https://rest.uniprot.org/ ; https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/embeddings/
- Agora / TREAT-AD: https://agora.adknowledgeportal.org/api/v1/genes/ENSG00000095970 ; Cary et al. 2024 *Alzheimer's Dement (TRCI)* 10(2):e12461, https://doi.org/10.1002/trc2.12461
- ChEMBL 37: https://www.ebi.ac.uk/chembl/api/data/status.json
- Cummings et al. 2025 AD pipeline (138 drugs / 182 trials; amyloid ≈18% per UNLV summary): https://pubmed.ncbi.nlm.nih.gov/40463637/
- CNA: https://github.com/immunogenomics/cna (PyPI `cna` 0.2.6); Reshef et al. 2022 *Nat Biotechnol*
- scDRS: https://pypi.org/project/scdrs/ (1.0.2), data https://figshare.com/articles/dataset/scDRS_data_release_030122/19312583 ; Zhang et al. 2022 *Nat Genet*, https://www.nature.com/articles/s41588-022-01167-z
- MAGMA: https://cncr.nl/research/magma/
- AlphaGenome (optional, key required): https://deepmind.google/blog/alphagenome-ai-for-better-understanding-the-genome/
- PsychAD: https://doi.org/10.1101/2024.10.31.24316513 ; Myeloid plasticity: https://doi.org/10.1038/s41588-026-02716-6 ; Rexach et al. 2024 *Cell*: https://doi.org/10.1016/j.cell.2024.08.019
- Mantis code paths: `C:\Users\goodp\OneDrive\Mantis\MantisAPI\src\synthesis\tasks\h5ad_tasks.py`, `...\synthesis\tasks\embedding.py`, `...\synthesis\config.py`, `...\db\models.py`; frontend `C:\Users\goodp\OneDrive\Mantis\Mantis\src\app\types.ts`, `...\src\app\new-space\h5ad\`

---

## Appendix A. SEA-AD MTG per-donor summary

This was computed today by streaming `MTG/RNAseq/SEAAD_MTG_RNAseq_final-nuclei_metadata.2026-06-22.csv` (1.15 GB) through a parser, without saving it to disk. Summary JSON: `C:\Users\goodp\AppData\Local\Temp\claude\C--Users-goodp-OneDrive-Mantis\edb643d9-7b40-4f02-9b28-375542249b1d\scratchpad\rp\donor_summary.json`.

**Cohort size**
- 1,178,694 final nuclei, all flagged `Used in analysis = True`, from **84 donors**. The 5 young reference donors are *not* in this 2026 MTG file, although the CELLxGENE MTG files (1,378,211 nuclei) include them.
- Age at death: median 90 (range 65–102). Sex: 51 F / 33 M.
- Nuclei per donor: median 14,372 (range 1,409–27,431).
- 28 donors have some 10x Multiome nuclei (`method` = `10xMulti` vs. `10Xv3.1`).

**APOE genotype**

| 2/2 | 2/3 | 3/3 | 2/4 | 3/4 | 4/4 |
|---|---|---|---|---|---|
| 1 | 11 | 47 | 2 | 17 | 6 |

That is **25 APOE4 carriers** out of 84 (8 with ε4/ε4 or ε2/ε4). This is the power ceiling for H4: roughly 25 vs. 59 donors.

**Neuropathology and cognition**
- ADNC: Not AD 9, Low 12, Intermediate 21, High 42.
- Braak: 0 (2), II (4), III (6), IV (23), V (34), VI (15).
- LATE: none 27, stage 1 (15), stage 2 (36), stage 3 (5), unclassifiable 1.
- Lewy body disease: neocortical 14, limbic 9, amygdala-predominant 10, others/none 51.
- Cognitive status: 42 dementia / 42 no dementia.

Cognitive status × ADNC:

| | High | Intermediate | Low | Not AD |
|---|---|---|---|---|
| Dementia | 31 | 7 | 4 | 0 |
| No dementia | **11** | 14 | 8 | 9 |

**11 donors have High ADNC but no dementia.** They are a built-in "resilience" group, the open-data echo of the resilience analysis in Mathys et al. 2023. Add H2c to the pre-registration: *do resilient donors (high ADNC, no dementia) show healthier neuronal brick residuals than demented high-ADNC donors?* With 11 vs. 31 donors, report the effect size and CI, not just P.

**Severely Affected Donors and neuronal cell loss**
- 11 donors are flagged "Severely Affected Donor", all of them High ADNC.
- They account for most of the neuronal dropout:
  - Sst: 5 donors have <100 nuclei, and 4 of those 5 are severe.
  - Pvalb: 5 donors have <100 nuclei, and all 5 are severe.
  - L2/3 IT: 1 donor has <100 nuclei, and it is severe.
- Rule: require ≥50 nuclei per donor per subclass for that donor's brick score, and report results with and without the severe donors.

**Taxonomy note.** In the 2026 taxonomy, microglia sit in the subclass **"Immune"**: 36,686 nuclei, median 380 per donor (range 67–1,036). The CELLxGENE file still calls it "Microglia-PVM" (40,000 cells). Choose one source per analysis and do not mix taxonomies.

Nuclei per subclass in the 2026 file (top): L2/3 IT 289,940; L4 IT 147,987; L5 IT 105,259; Oligodendrocyte 90,578; Vip 89,661; Pvalb 78,746; Astrocyte 65,450; Sst 44,695; L6 IT 37,902; Lamp5 36,966; Immune 36,686; OPC 28,907.
