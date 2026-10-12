#!/usr/bin/env python
"""Open Targets 26.09 AD associations (GraphQL) + GWAS Catalog v2 AD associations with first-report years.

Outputs (raw/ot, raw/gwas):
  ot_assoc_indirect.parquet   all associated targets of MONDO_0004975 incl. descendants (the pre-registered genetic metric source)
  ot_assoc_direct.parquet     direct-only (sensitivity S5)
  gwas_ad_associations.parquet, gwas_pmid_years.csv
"""
import os, sys, time, json, datetime
from pathlib import Path
import requests
import polars as pl

ROOT = Path(os.path.expanduser("~/orcd/pool/mantis-films/research/lit/raw"))
OT = ROOT / "ot"; OT.mkdir(parents=True, exist_ok=True)
GW = ROOT / "gwas"; GW.mkdir(parents=True, exist_ok=True)
EFO = "MONDO_0004975"
sess = requests.Session()
sess.headers["User-Agent"] = "mantis-films-lit/1.0"


def retry(fn, tries=8):
    for a in range(tries):
        try:
            return fn()
        except Exception as e:
            print(f"  retry {a}: {type(e).__name__} {e}", flush=True)
            time.sleep(min(90, 3 * 2 ** a))
    raise RuntimeError("giving up")


GQL = "https://api.platform.opentargets.org/api/v4/graphql"
Q = """query($id:String!,$i:Int!,$s:Int!,$ind:Boolean!){ disease(efoId:$id){ id name
 associatedTargets(page:{index:$i,size:$s}, enableIndirect:$ind){ count
  rows{ target{id approvedSymbol approvedName biotype} score datatypeScores{id score} } } } }"""


def ot_all(indirect):
    rows, i, s = [], 0, 1000
    while True:
        def go():
            r = sess.post(GQL, json=dict(query=Q, variables=dict(id=EFO, i=i, s=s, ind=indirect)), timeout=180)
            r.raise_for_status()
            j = r.json()
            if "errors" in j:
                raise RuntimeError(j["errors"])
            return j["data"]["disease"]["associatedTargets"]
        at = retry(go)
        for x in at["rows"]:
            d = dict(ensembl_id=x["target"]["id"], symbol=x["target"]["approvedSymbol"], name=x["target"]["approvedName"],
                     biotype=x["target"]["biotype"], overall_score=x["score"])
            for ds in x["datatypeScores"]:
                d["dt_" + ds["id"]] = ds["score"]
            rows.append(d)
        print(f"  OT indirect={indirect} page {i}: {len(at['rows'])} rows (total {at['count']})", flush=True)
        if len(rows) >= at["count"] or not at["rows"]:
            return rows, at["count"]
        i += 1
        time.sleep(0.3)


def ot():
    meta = retry(lambda: sess.post(GQL, json=dict(query="{ meta { apiVersion { x y z } dataVersion { year month iteration } } }"), timeout=60).json())
    for ind, name in ((True, "ot_assoc_indirect"), (False, "ot_assoc_direct")):
        rows, count = ot_all(ind)
        df = pl.DataFrame(rows, infer_schema_length=None)
        df.write_parquet(OT / f"{name}.parquet")
        print(name, df.shape, "count reported", count, flush=True)
    (OT / "meta.json").write_text(json.dumps(dict(retrieved_utc=datetime.datetime.utcnow().isoformat() + "Z", meta=meta["data"]["meta"], efo=EFO), indent=2))


def gwas():
    base = "https://www.ebi.ac.uk/gwas/rest/api/v2/associations"
    out, page, size = [], 0, 200
    while True:
        def go():
            r = sess.get(base, params=dict(efo_id=EFO, page=page, size=size), headers={"Accept": "application/json"}, timeout=180)
            r.raise_for_status()
            return r.json()
        j = retry(go)
        items = j.get("_embedded", {}).get("associations", [])
        for a in items:
            out.append(dict(association_id=a.get("association_id"), accession_id=a.get("accession_id"), pubmed_id=a.get("pubmed_id"),
                            first_author=a.get("first_author"), p_value=a.get("p_value"),
                            mapped_genes=";".join(a.get("mapped_genes") or []),
                            efo_traits=";".join(t.get("efo_id", "") for t in a.get("efo_traits", [])),
                            reported_trait=";".join(a.get("reported_trait") or []),
                            rsids=";".join(s.get("rs_id", "") for s in a.get("snp_allele", []))))
        tp = j["page"]["totalPages"]
        print(f"  GWAS page {page + 1}/{tp}: {len(items)}", flush=True)
        page += 1
        if page >= tp:
            break
        time.sleep(0.2)
    df = pl.DataFrame(out, infer_schema_length=None)
    df.write_parquet(GW / "gwas_ad_associations.parquet")
    print("GWAS associations", df.shape, "reported total", j["page"]["totalElements"], flush=True)
    # first-report years via PubMed esummary (<= 2.5 req/s)
    pmids = sorted({int(p) for p in df["pubmed_id"].drop_nulls().to_list() if str(p).isdigit()})
    rows = []
    for i in range(0, len(pmids), 150):
        b = pmids[i:i + 150]
        def go():
            r = sess.post("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
                          data=dict(db="pubmed", id=",".join(map(str, b)), retmode="json", tool="mantis-films-lit"), timeout=120)
            r.raise_for_status()
            return r.json()["result"]
        res = retry(go)
        for p in b:
            e = res.get(str(p), {})
            pd = e.get("pubdate", "") or ""
            sd = e.get("sortpubdate", "") or ""
            rows.append(dict(pmid=p, pubdate=pd, sortpubdate=sd, epubdate=e.get("epubdate", ""), title=e.get("title", ""),
                             source=e.get("source", "")))
        time.sleep(0.45)
    pl.DataFrame(rows).write_csv(GW / "gwas_pmid_years.csv")
    (GW / "meta.json").write_text(json.dumps(dict(retrieved_utc=datetime.datetime.utcnow().isoformat() + "Z", efo=EFO,
                                                  n_assoc=df.height, n_pmids=len(pmids)), indent=2))
    print("PMID years", len(rows), flush=True)


if __name__ == "__main__":
    which = sys.argv[1:] or ["ot", "gwas"]
    if "ot" in which:
        ot()
    if "gwas" in which:
        gwas()
    print("OT_GWAS_DONE", flush=True)
