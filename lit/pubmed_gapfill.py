#!/usr/bin/env python
"""Second pass over every bin: fetch PMIDs that the first pass returned from ESearch but did not parse
(NCBI Bookshelf documents, PubmedBookArticle), then rebuild abstracts.parquet. Same query, same rate limit."""
import os, sys, json, time, datetime
from pathlib import Path
import polars as pl
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pubmed_fetch as pf

OUT, CH = pf.OUT, pf.CH
def main():
    r0 = pf.call("esearch.fcgi", dict(db="pubmed", term=pf.TERM, retmax=0, retmode="json"))
    total_all = int(r0.json()["esearchresult"]["count"])
    still, added = [], 0
    for f in sorted(CH.glob("*.parquet")):
        label = f.stem
        if label.startswith("gap_"):
            continue
        a, b = (label, label) if "-" not in label else (label.replace("-", "/"), label.replace("-", "/"))
        count, ids = pf.esearch(a, b, retmax=10000)
        have = set(pl.read_parquet(f)["pmid"].to_list()) if f.stat().st_size > 0 and "pmid" in pl.read_parquet_schema(f) else set()
        gf = CH / f"gap_{label}.parquet"
        have_gap = set(pl.read_parquet(gf)["pmid"].to_list()) if gf.exists() else set()
        miss = [i for i in ids if int(i) not in have and int(i) not in have_gap]
        if not miss:
            continue
        rows = pf.fetch_ids(miss)
        if rows:
            df = pl.DataFrame(rows, schema_overrides={"year": pl.Int32, "year_pubdate": pl.Int32, "year_epub": pl.Int32}, infer_schema_length=None)
            df.write_parquet(gf, compression="zstd")
            added += len(rows)
        got = {r["pmid"] for r in rows}
        still += [i for i in miss if int(i) not in got]
        print(f"{label}: esearch {count}, chunk had {len(have)}, missing {len(miss)}, now parsed {len(rows)}", flush=True)
    parts = [pl.read_parquet(p) for p in sorted(CH.glob("*.parquet"))]
    parts = [p for p in parts if p.height and "abstract" in p.columns]
    df = pl.concat(parts, how="diagonal_relaxed")
    n_rows = df.height
    df = df.with_columns(pl.col("record_type").fill_null("article")).unique(subset="pmid", keep="first", maintain_order=True)
    n_noabs = df.filter(pl.col("abstract").str.len_chars() == 0).height
    df = df.filter(pl.col("abstract").str.len_chars() > 0)
    df.write_parquet(OUT / "abstracts.parquet", compression="zstd")
    (OUT / "missing_after_gapfill.txt").write_text("\n".join(still))
    prov = dict(retrieved_utc=datetime.datetime.utcnow().isoformat() + "Z", term=pf.TERM, esearch_total_all_dates=total_all, rows_before_dedupe=n_rows,
                unique_pmids_with_abstract=df.height, dropped_empty_abstract=n_noabs, book_records=int((df["record_type"] == "book").sum()),
                null_year=int(df["year"].null_count()), still_missing_after_gapfill=len(still), gap_rows_added=added)
    (OUT / "provenance.json").write_text(json.dumps(prov, indent=2))
    print(json.dumps(prov, indent=2), flush=True)
if __name__ == "__main__":
    main()
