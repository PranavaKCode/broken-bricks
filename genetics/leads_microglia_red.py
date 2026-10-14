#!/usr/bin/env python
"""Leads in the microglia red program: recompute the selection from the repository's own tables.

Selection (all three must hold):
  1. genetic score >= 0.7        results/gene_master.csv, column G_genetic (Open Targets Platform 26.09, genetic_association
                                 score for MONDO_0004975, descendants included)
  2. fewer than 100 Alzheimer's  results/gene_master.csv, column A_g (abstract-level count, amendment A2: a PubTator3 (PMID, human
     abstracts                   gene) pair counts only if its mention string occurs in the title or abstract of the corpus record)
  3. in the 200-gene             results/genetics/program_membership.csv, rows with set == "Microglia_red" (top 200 genes by
     microglia_red program       donor-blocked Spearman correlation with the CNA ncorr of the CPS_Local map, primary `all` run)
Join key: Ensembl gene ID. Genes with neither a genetic score nor an abstract are not listed in gene_master and cannot pass 1 and 2.

Writes (nothing is typed by hand):
  results/genetics/leads_microglia_red.csv   the selected rows, with the filter written in a header comment (read with comment="#")
  results/agent_hypotheses.csv               plain statement of where the film's leads come from; the lead symbols are read from the
                                             selection above
Exits with an error, writing nothing, if the selection is not exactly {CTSH, PILRA} (the expectation recorded before this run).

Run from anywhere:  python genetics/leads_microglia_red.py
"""
import csv
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
GM = REPO / "results" / "gene_master.csv"
PM = REPO / "results" / "genetics" / "program_membership.csv"
OUT_LEADS = REPO / "results" / "genetics" / "leads_microglia_red.csv"
OUT_AGENT = REPO / "results" / "agent_hypotheses.csv"

G_MIN = 0.7
A_MAX_EXCL = 100
PROGRAM = "Microglia_red"
EXPECTED = {"CTSH", "PILRA"}

gm = pd.read_csv(GM)
pm = pd.read_csv(PM)

prog = pm[pm["set"] == PROGRAM].copy()
assert len(prog) == 200 and prog["ensembl_id"].is_unique, "program must be 200 unique genes"

m = prog.merge(gm, on="ensembl_id", how="left", suffixes=("", "_gm"))
# gene_master repeats some Ensembl IDs (shared by several NCBI records); none of them is a program gene
assert len(m) == len(prog), "an Ensembl ID of the program maps to more than one gene_master row"

n_program = len(prog)
n_in_gm = int(m["gene_id"].notna().sum())
n_g = int((m["G_genetic"] >= G_MIN).sum())
sel = m[(m["G_genetic"] >= G_MIN) & (m["A_g"] < A_MAX_EXCL)].copy()
n_pass = len(sel)

# cross-check against the flag written by lit/20_metrics.py (strong_G07 & A_g < 100)
flag = m[m["sparse_G07_A100"] == True]  # noqa: E712 (NaN = gene not in gene_master)
assert set(flag["ensembl_id"]) == set(sel["ensembl_id"]), "filter disagrees with gene_master.sparse_G07_A100"
# the program symbol and the gene_master symbol must agree for the selected rows
assert (sel["symbol"] == sel["symbol_gm"]).all(), "symbol mismatch between program and gene_master"

got = set(sel["symbol"])
if got != EXPECTED:
    sys.exit(f"STOP: the filter returned {sorted(got)}, expected {sorted(EXPECTED)}; nothing written")

sel = sel.sort_values("symbol")
cols = {
    "symbol": "symbol",
    "ensembl_id": "ensembl_id",
    "gene_id": "ncbi_gene_id",
    "G_genetic": "genetic_score_open_targets_26_09",
    "A_g": "alzheimers_abstracts_abstract_level",
    "rho": "program_rho_with_ncorr",
}
out = sel[list(cols)].rename(columns=cols)
out["ncbi_gene_id"] = out["ncbi_gene_id"].astype(int)

header = [
    "leads_microglia_red.csv: genes that pass all three conditions below. Computed by genetics/leads_microglia_red.py from the repository's own tables; do not edit by hand. Read with pandas.read_csv(path, comment='#').",
    "Filter: (1) genetic score >= 0.7 (results/gene_master.csv, G_genetic: Open Targets Platform 26.09 genetic_association score for MONDO_0004975);",
    "        (2) fewer than 100 Alzheimer's abstracts at the abstract-level count (gene_master.csv, A_g < 100; amendment A2: a PubTator3 gene annotation counts only if its mention string occurs in the title or abstract);",
    "        (3) member of the 200-gene microglia_red program (results/genetics/program_membership.csv, set == Microglia_red: top 200 genes by donor-blocked Spearman correlation with the CNA ncorr of the CPS_Local map, primary run).",
    "Join key: Ensembl gene ID.",
    f"Counts at each step: {n_program} program genes; {n_in_gm} have a row in gene_master; {n_g} of those have genetic score >= 0.7; {n_pass} also have fewer than 100 abstracts (the rows below).",
    "Status: leads to test. These genes were selected by a filter on public tables; no experiment on them was done in this repository.",
]
with open(OUT_LEADS, "w", encoding="utf-8", newline="") as fh:
    for line in header:
        fh.write("# " + line + "\n")
    out.to_csv(fh, index=False, lineterminator="\n", float_format="%.10g")

leads = ";".join(out["symbol"])
rows = [
    ("film_leads_source", "The leads shown in the film come from the Mantis map filter re-running the selection in results/genetics/leads_microglia_red.csv: genetic score >= 0.7, fewer than 100 Alzheimer's abstracts (abstract-level count), member of the 200-gene microglia_red program."),
    ("film_leads", leads),
    ("film_leads_count", str(len(out))),
    ("leads_table", "results/genetics/leads_microglia_red.csv (recomputed by genetics/leads_microglia_red.py from results/gene_master.csv and results/genetics/program_membership.csv)"),
    ("agent_run_used_in_film", "no"),
    ("agent_generated_hypotheses_in_this_repository", "0"),
    ("recorded_composer_attempt", "One recorded Mantis Composer attempt failed on a tool error and produced no usable output. It is not used in the film and nothing from it is in this repository."),
    ("status_of_leads", "Leads worth testing. Not discoveries and not drug targets: no experiment on them is reported here."),
]
with open(OUT_AGENT, "w", encoding="utf-8", newline="") as fh:
    w = csv.writer(fh, lineterminator="\n")
    w.writerow(["field", "value"])
    w.writerows(rows)

print(f"program genes {n_program}; in gene_master {n_in_gm}; G >= {G_MIN}: {n_g}; pass: {n_pass} -> {leads}")
print(out.to_string(index=False))
