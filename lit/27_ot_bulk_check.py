#!/usr/bin/env python
"""Cross-check of the GraphQL direct-association genetic scores against the Open Targets 26.09 bulk parquet association_by_datatype_direct."""
import os, glob
from pathlib import Path
import polars as pl
import numpy as np
from scipy.stats import spearmanr
LIT = Path(os.path.expanduser("~/orcd/pool/mantis-films/research/lit"))
files = sorted(glob.glob(str(LIT / "raw" / "ot_direct" / "*.parquet")))
print(len(files), "parts")
lf = pl.scan_parquet(files)
print(lf.collect_schema())
b = lf.filter((pl.col("diseaseId") == "MONDO_0004975") & (pl.col("aggregationValue") == "genetic_association")).select(["targetId", pl.col("associationScore").alias("score"), "aggregationType"]).collect()
print(b["aggregationType"].value_counts())
g = pl.read_parquet(LIT / "raw" / "ot" / "ot_assoc_direct.parquet").select(["ensembl_id", "dt_genetic_association"]).drop_nulls("dt_genetic_association")
j = b.join(g, left_on="targetId", right_on="ensembl_id", how="full")
both = j.drop_nulls(["score", "dt_genetic_association"])
res = dict(bulk_rows=b.height, graphql_direct_rows=g.height, matched=both.height, only_bulk=int(j["dt_genetic_association"].is_null().sum()), only_graphql=int(j["score"].is_null().sum()),
           max_abs_diff=float((both["score"] - both["dt_genetic_association"]).abs().max()), spearman=float(spearmanr(both["score"], both["dt_genetic_association"]).statistic),
           bulk_ge_0_5=int((b["score"] >= 0.5).sum()), bulk_ge_0_7=int((b["score"] >= 0.7).sum()),
           graphql_direct_ge_0_5=int((g["dt_genetic_association"] >= 0.5).sum()), graphql_direct_ge_0_7=int((g["dt_genetic_association"] >= 0.7).sum()))
print(res)
out = LIT / "out" / "results" / "s5_ot_bulk_vs_graphql_direct.csv"
pl.DataFrame([dict(item=k, value=float(v)) for k, v in res.items()]).write_csv(out)
print("OT_BULK_CHECK_DONE")
