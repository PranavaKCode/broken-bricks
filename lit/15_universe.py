#!/usr/bin/env python
"""Gene universe U for the Agora query (preregistration D3): genes with G_g > 0 or A_g >= 20.
Writes work/agora_universe.csv (ensembl_id) and results/universe_summary.csv."""
import os
from pathlib import Path
import polars as pl

LIT = Path(os.path.expanduser("~/orcd/pool/mantis-films/research/lit"))
RAW, WORK, RES = LIT / "raw", LIT / "work", LIT / "out" / "results"
RES.mkdir(parents=True, exist_ok=True)
ga = pl.read_parquet(WORK / "gene_attention.parquet")
ot = pl.read_parquet(RAW / "ot" / "ot_assoc_indirect.parquet")
e2g = {e: g for e, g in zip(ga["ensembl_id"].to_list(), ga["gene_id"].to_list()) if e}
s2g = {s: g for s, g in zip(ga["symbol"].to_list(), ga["gene_id"].to_list())}
ot = ot.with_columns(pl.col("dt_genetic_association").fill_null(0.0).alias("G"))
ot = ot.with_columns(pl.col("ensembl_id").replace_strict(e2g, default=None, return_dtype=pl.Int64).alias("gene_id"))
ot = ot.with_columns(pl.when(pl.col("gene_id").is_null()).then(pl.col("symbol").replace_strict(s2g, default=None, return_dtype=pl.Int64)).otherwise(pl.col("gene_id")).alias("gene_id"))
unm = ot.filter(pl.col("gene_id").is_null())
print("OT targets", ot.height, "unmapped to NCBI gene", unm.height, "unmapped with G>0", unm.filter(pl.col("G") > 0).height, "G>=0.7", unm.filter(pl.col("G") >= 0.7).height)
ot_g = ot.filter(pl.col("G") > 0)
ids = set(ot_g["ensembl_id"].to_list())
big = ga.filter(pl.col("A_g") >= 20)
ids |= {e for e in big["ensembl_id"].to_list() if e}
n_noens = big.filter(pl.col("ensembl_id").is_null()).height
pl.DataFrame(dict(ensembl_id=sorted(ids))).write_csv(WORK / "agora_universe.csv")
pl.DataFrame([dict(item="OT targets", value=ot.height), dict(item="OT targets with genetic_association > 0", value=ot_g.height),
              dict(item="OT targets unmapped to NCBI Gene", value=unm.height),
              dict(item="OT unmapped with G >= 0.7", value=unm.filter(pl.col("G") >= 0.7).height),
              dict(item="genes with A_g >= 20", value=big.height), dict(item="  of which without Ensembl ID", value=n_noens),
              dict(item="universe size (Ensembl IDs queried at Agora)", value=len(ids))]).write_csv(RES / "universe_summary.csv")
print("universe", len(ids))
