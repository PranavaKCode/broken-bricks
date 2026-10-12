# Mantis floors 3 and 4 (genes and literature, in register)

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
- A gene's `ad_paper_count` counts abstracts that mention it; 84% of 50 audited (PMID, gene) pairs were correct, lower for rarely mentioned genes (76% for genes with fewer than 100 abstracts). See `results/ner_audit_summary.csv`.
- No donor or patient data are involved. Do not describe genes as cures or as "the answer"; see the "Never claim" list in the report.
