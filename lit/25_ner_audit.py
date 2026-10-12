#!/usr/bin/env python
"""V2: draw the 50 (PMID, gene) pairs for the named-entity audit (preregistration 2.4): 25 from genes with total A_g < 100, 25 from A_g >= 100.
Writes out/results/ner_audit_50_sample.csv (judgment columns left empty; filled in by the audit)."""
import os, re
from pathlib import Path
import numpy as np
import polars as pl

SEED = 20261011
LIT = Path(os.path.expanduser("~/orcd/pool/mantis-films/research/lit"))
WORK, RES = LIT / "work", LIT / "out" / "results"
pairs = pl.read_parquet(WORK / "pairs_human.parquet")
ga = pl.read_parquet(WORK / "gene_attention.parquet").select(["gene_id", "symbol", "A_g"])
corp = pl.read_parquet(WORK / "corpus.parquet").select(["pmid", "title", "abstract"])
d = pairs.join(ga, on="gene_id").join(corp, on="pmid")
rng = np.random.default_rng(SEED)
out = []
for name, sel in (("A_lt_100", d.filter(pl.col("A_g") < 100)), ("A_ge_100", d.filter(pl.col("A_g") >= 100))):
    idx = rng.choice(sel.height, 25, replace=False)
    s = sel[idx.tolist()].with_columns(pl.lit(name).alias("stratum"))
    out.append(s)
s = pl.concat(out)
def ctx(title, abstract, mentions):
    t = f"{title} || {abstract}"
    for m in sorted(set(mentions.split("|")), key=len, reverse=True):
        m = m.strip()
        if not m: continue
        i = t.lower().find(m.lower())
        if i >= 0:
            return t[max(0, i - 170): i + len(m) + 170].replace("\n", " ")
    return "[mention string not found in title/abstract] " + t[:300]
s = s.with_columns(pl.struct(["title", "abstract", "mentions"]).map_elements(lambda r: ctx(r["title"], r["abstract"], r["mentions"]), return_dtype=pl.Utf8).alias("context"))
s = s.select(["stratum", "pmid", "gene_id", "symbol", "A_g", "mentions", "abs_verified", "year", "context"]).with_columns(
    pl.lit("").alias("judgment"), pl.lit("").alias("note"))
s.write_csv(RES / "ner_audit_50_sample.csv")
print(s.select(["stratum", "pmid", "symbol", "A_g", "mentions", "abs_verified"]))
