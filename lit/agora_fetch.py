#!/usr/bin/env python
"""Agora (AD Knowledge Portal) per-gene biodomains and TREAT-AD scores for the gene universe U (preregistration D3).

Input : work/agora_universe.csv  (column ensembl_id)
Output: raw/agora/agora_genes.jsonl (one compact record per gene), raw/agora/meta.json
Rate  : <= 4 requests/s overall (limit is 5/s), 3 worker threads sharing one limiter. Resumable.
"""
import os, sys, json, time, threading, datetime, csv
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import requests

LIT = Path(os.path.expanduser("~/orcd/pool/mantis-films/research/lit"))
OUT = LIT / "raw" / "agora"; OUT.mkdir(parents=True, exist_ok=True)
UNIV = Path(os.environ.get("AGORA_UNIV", LIT / "work" / "agora_universe.csv"))
JL = OUT / os.environ.get("AGORA_JL", "agora_genes.jsonl")
SKIP = [OUT / x for x in os.environ.get("AGORA_SKIP", "").split(",") if x]
URL = "https://agora.adknowledgeportal.org/api/v1/genes/"
INTERVAL = float(os.environ.get("AGORA_INTERVAL", "0.25"))
lock = threading.Lock()
_last = [0.0]
tl = threading.local()


def sess():
    if not hasattr(tl, "s"):
        tl.s = requests.Session()
        tl.s.headers["User-Agent"] = "mantis-films-lit/1.0"
    return tl.s


def limiter():
    with lock:
        now = time.time()
        wait = _last[0] + INTERVAL - now
        _last[0] = max(now, _last[0] + INTERVAL)
    if wait > 0:
        time.sleep(wait)


def fetch(eid):
    for a in range(8):
        limiter()
        try:
            r = sess().get(URL + eid, timeout=(10, 40), headers={"Connection": "close"})
            if r.status_code == 404:
                return dict(ensembl_id=eid, status=404)
            if r.status_code == 200:
                if "json" not in r.headers.get("content-type", ""):
                    return dict(ensembl_id=eid, status=204)
                d = r.json()
                if d is None:
                    return dict(ensembl_id=eid, status=404, note="null body (gene not in Agora)")
                os_ = d.get("overall_scores") or {}
                return dict(ensembl_id=eid, status=200, symbol=d.get("hgnc_symbol"), name=d.get("name"),
                            biodomains=d.get("biodomains") or [], target_risk_score=os_.get("target_risk_score"),
                            genetics_score=os_.get("genetics_score"), multi_omics_score=os_.get("multi_omics_score"),
                            total_nominations=d.get("total_nominations"), is_igap=d.get("is_igap"), is_eqtl=d.get("is_eqtl"))
            time.sleep(min(30, 2 ** a))
        except Exception:
            time.sleep(min(30, 2 ** a))
    return dict(ensembl_id=eid, status=-1)


def main():
    ids = [r["ensembl_id"] for r in csv.DictReader(open(UNIV)) if r["ensembl_id"]]
    if os.environ.get("AGORA_REVERSE"):
        ids = ids[::-1]
    done = set()
    for f_ in [JL] + SKIP:
      if f_.exists():
        for line in open(f_):
            try:
                d = json.loads(line)
                if d.get("status") in (200, 404, 204):
                    done.add(d["ensembl_id"])
            except Exception:
                pass
    todo = [i for i in ids if i not in done]
    print(f"universe {len(ids)}; done {len(done)}; todo {len(todo)}", flush=True)
    t0 = time.time()
    with open(JL, "a") as fh, ThreadPoolExecutor(max_workers=3) as ex:
        for k, rec in enumerate(ex.map(fetch, todo), 1):
            fh.write(json.dumps(rec) + "\n")
            fh.flush()
            if k % 100 == 0:
                print(f"  {k}/{len(todo)}  {k / (time.time() - t0):.2f} req/s", flush=True)
    (OUT / "meta.json").write_text(json.dumps(dict(retrieved_utc=datetime.datetime.utcnow().isoformat() + "Z", n_requested=len(todo)), indent=2))
    print("AGORA_DONE", flush=True)


if __name__ == "__main__":
    main()
