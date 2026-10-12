#!/usr/bin/env python
"""Download all PubMed Alzheimer records with abstracts (preregistration D0).

Query: ("Alzheimer Disease"[MeSH Terms] OR alzheimer*[Title/Abstract]) AND hasabstract
Split by publication date (datetype=pdat): years, or months when a year holds > 9,000 records.
E-utilities, no API key, <= 2.5 requests/s (limit is 3/s).
Resumable: each bin is written to its own parquet chunk.
"""
import os, re, sys, time, json, datetime
from pathlib import Path
import requests
import polars as pl
from lxml import etree

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
TERM = '("Alzheimer Disease"[MeSH Terms] OR alzheimer*[Title/Abstract]) AND hasabstract'
OUT = Path(os.environ.get("PUBMED_OUT", os.path.expanduser("~/orcd/pool/mantis-films/research/lit/raw/pubmed")))
CH = OUT / "chunks"
CH.mkdir(parents=True, exist_ok=True)
MIN_INTERVAL = 0.4
FETCH_BATCH = 400
Y0, Y1 = 1800, datetime.date.today().year + 1

sess = requests.Session()
sess.headers["User-Agent"] = "mantis-films-lit/1.0 (research; tool=mantis-films-lit)"
_last = [0.0]


def call(endpoint, params, post=False):
    params = dict(params)
    params["tool"] = "mantis-films-lit"
    for attempt in range(10):
        wait = MIN_INTERVAL - (time.time() - _last[0])
        if wait > 0:
            time.sleep(wait)
        _last[0] = time.time()
        try:
            r = sess.post(BASE + endpoint, data=params, timeout=180) if post else sess.get(BASE + endpoint, params=params, timeout=180)
            if r.status_code == 200:
                return r
            if r.status_code in (429, 500, 502, 503, 504):
                print(f"  HTTP {r.status_code} on {endpoint}, retry {attempt}", flush=True)
                time.sleep(min(120, 3 * 2 ** attempt))
                continue
            r.raise_for_status()
        except requests.RequestException as e:
            print(f"  {type(e).__name__} on {endpoint}, retry {attempt}", flush=True)
            time.sleep(min(120, 3 * 2 ** attempt))
    raise RuntimeError(f"giving up on {endpoint} {params}")


def esearch(mindate, maxdate, retmax=0):
    r = call("esearch.fcgi", dict(db="pubmed", term=TERM, datetype="pdat", mindate=mindate, maxdate=maxdate,
                                  retmax=retmax, retmode="json"))
    j = r.json()["esearchresult"]
    return int(j["count"]), j.get("idlist", [])


def txt(el):
    return "".join(el.itertext()).strip() if el is not None else ""


def year_of(el):
    """Year from a PubDate-like element (Year or MedlineDate)."""
    if el is None:
        return None
    y = el.findtext("Year")
    if y and y.isdigit():
        return int(y)
    m = re.search(r"(1[5-9]\d\d|20\d\d)", el.findtext("MedlineDate") or "")
    return int(m.group(1)) if m else None


def parse(xml_bytes):
    root = etree.fromstring(xml_bytes)
    rows = []
    for pa in root.findall("PubmedArticle"):
        mc = pa.find("MedlineCitation")
        art = mc.find("Article")
        pmid = int(mc.findtext("PMID"))
        title = txt(art.find("ArticleTitle"))
        parts = []
        for at in art.findall("Abstract/AbstractText"):
            lab = at.get("Label")
            t = txt(at)
            parts.append(f"{lab}: {t}" if lab and t else t)
        abstract = " ".join(p for p in parts if p)
        pd_year = year_of(art.find("Journal/JournalIssue/PubDate"))
        ad_years = [int(a.findtext("Year")) for a in art.findall("ArticleDate") if (a.findtext("Year") or "").isdigit()]
        dc_year = year_of(mc.find("DateCompleted"))
        dr_year = year_of(mc.find("DateRevised"))
        year = pd_year or (min(ad_years) if ad_years else None) or dc_year or dr_year
        mesh = ";".join(txt(d) for d in mc.findall("MeshHeadingList/MeshHeading/DescriptorName"))
        rows.append(dict(
            pmid=pmid, record_type="article", year=year, year_pubdate=pd_year, year_epub=min(ad_years) if ad_years else None,
            journal=art.findtext("Journal/Title") or "", journal_iso=art.findtext("Journal/ISOAbbreviation") or "",
            title=title, abstract=abstract, mesh=mesh,
            pubtypes=";".join(txt(p) for p in art.findall("PublicationTypeList/PublicationType")),
            lang=";".join(txt(l) for l in art.findall("Language")),
        ))
    # NCBI Bookshelf documents indexed in PubMed (e.g. MICAD chapters) are part of the query result and kept (record_type = book)
    for pb in root.findall("PubmedBookArticle"):
        bd = pb.find("BookDocument")
        pmid = int(bd.findtext("PMID"))
        title = txt(bd.find("ArticleTitle")) or txt(bd.find("Book/BookTitle"))
        parts = []
        for at in bd.findall("Abstract/AbstractText"):
            lab = at.get("Label")
            t = txt(at)
            parts.append(f"{lab}: {t}" if lab and t else t)
        year = year_of(bd.find("Book/PubDate")) or year_of(bd.find("Book/BeginningDate"))
        rows.append(dict(
            pmid=pmid, record_type="book", year=year, year_pubdate=year, year_epub=None,
            journal=txt(bd.find("Book/BookTitle")), journal_iso="", title=title, abstract=" ".join(p for p in parts if p), mesh="",
            pubtypes=";".join(txt(p) for p in bd.findall("PublicationType")), lang=";".join(txt(l) for l in bd.findall("Language")),
        ))
    return rows


def fetch_ids(ids):
    rows = []
    for i in range(0, len(ids), FETCH_BATCH):
        b = ids[i:i + FETCH_BATCH]
        r = call("efetch.fcgi", dict(db="pubmed", id=",".join(b), retmode="xml", rettype="abstract"), post=True)
        rows.extend(parse(r.content))
    return rows


def do_bin(label, mindate, maxdate, count_hint=None):
    f = CH / f"{label}.parquet"
    if f.exists():
        return None
    count, ids = esearch(mindate, maxdate, retmax=10000)
    if count > 10000:
        raise RuntimeError(f"bin {label} has {count} > 10000")
    rows = fetch_ids(ids)
    df = pl.DataFrame(rows, schema_overrides={"year": pl.Int32, "year_pubdate": pl.Int32, "year_epub": pl.Int32}, infer_schema_length=None) if rows else pl.DataFrame()
    if rows:
        df.write_parquet(f, compression="zstd")
    else:
        pl.DataFrame({"pmid": []}).write_parquet(f)
    return dict(bin=label, esearch_count=count, ids_returned=len(ids), parsed=len(rows))


def main():
    t0 = time.time()
    r0 = call("esearch.fcgi", dict(db="pubmed", term=TERM, retmax=0, retmode="json"))
    total_all = int(r0.json()["esearchresult"]["count"])
    print(f"retrieval {datetime.datetime.utcnow().isoformat()}Z  total(all dates)={total_all}", flush=True)
    log = []
    for y in range(Y0, Y1 + 1):
        c, _ = esearch(str(y), str(y))
        if c == 0:
            continue
        if c <= 9000:
            r = do_bin(f"{y}", str(y), str(y))
            if r:
                log.append(r)
                print(f"{y}: {r}", flush=True)
        else:
            for m in range(1, 13):
                mm = f"{y}/{m:02d}"
                r = do_bin(f"{y}-{m:02d}", mm, mm)
                if r:
                    log.append(r)
                    print(f"{y}-{m:02d}: {r}", flush=True)
    if log:
        pl.DataFrame(log).write_csv(OUT / f"bins_log_{int(time.time())}.csv")
    # concatenate
    parts = [pl.read_parquet(p) for p in sorted(CH.glob("*.parquet"))]
    parts = [p for p in parts if p.height and "abstract" in p.columns]
    df = pl.concat(parts, how="diagonal_relaxed")
    n_rows = df.height
    df = df.unique(subset="pmid", keep="first", maintain_order=True)
    df = df.filter(pl.col("abstract").str.len_chars() > 0)
    df.write_parquet(OUT / "abstracts.parquet", compression="zstd")
    prov = dict(retrieved_utc=datetime.datetime.utcnow().isoformat() + "Z", term=TERM, esearch_total_all_dates=total_all,
                rows_before_dedupe=n_rows, unique_pmids_with_abstract=df.height, seconds=round(time.time() - t0))
    (OUT / "provenance.json").write_text(json.dumps(prov, indent=2))
    print(json.dumps(prov, indent=2), flush=True)


if __name__ == "__main__":
    main()
