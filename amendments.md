# Amendments and deviations

None recorded at the time of the preregistration commit. Each later entry: date, kind (deviation / implementation), what changed, why, what was known at the time, commit hash.

## A1 (2026-10-11, sensitivity added before any metric was computed; kind: implementation)

What was seen: only the first lines of `gene2pubtator3.gz` (five tab-separated columns `PMID, Type, ConceptID, Mentions, Resource`; no header; Resource = "PubTator3") and the `Homo_sapiens.gene_info.gz` and UniProt headers. No attention counts had been computed.
Issue: the file is keyed on PMID and does not say whether an annotation came from the abstract or from PMC full text.
Decision: the pre-registered literal definition D1 stays the primary attention metric. We add sensitivity **S8**: `A_g_abs` counts a (PMID, gene) pair only if at least one of its Mentions strings occurs (case-insensitive, whole-token) in the title or abstract of the corpus record. The report gives the fraction of pairs that pass, the Spearman correlation of `A_g` with `A_g_abs`, and the Jaccard overlap of the two top-25 under-explored lists, and states which headline claims hold under both. No result switches definition.
Thresholds and decision rules: unchanged.
