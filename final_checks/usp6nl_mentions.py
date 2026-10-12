#!/usr/bin/env python
"""Final check 5: read the PubTator3 mention strings behind the 10 abstract-level USP6NL attention count (NCBI Gene 9712).

Reads  lit/work/pairs_human.parquet (pmid, gene_id, mentions, resource), lit/work/corpus.parquet (title, abstract, year),
       lit/raw/ncbi/Homo_sapiens.gene_info.gz (synonyms).
Applies the amendment-A2 verifier (A_g counts a pair if a mention string is found in the title/abstract after normalisation)
and prints, for each PMID, the mention strings, whether the verifier accepts the pair, and the sentence that contains the match.
Writes final_checks_out/usp6nl_mentions.csv
"""
import os, re, gzip, unicodedata, sys, json
import pandas as pd, pyarrow.parquet as pq

LIT = os.path.expanduser("~/orcd/pool/mantis-films/research/lit")
OUT = os.path.expanduser("~/orcd/pool/mantis-films/research/seaad/results/final_checks"); os.makedirs(OUT, exist_ok=True)
GID = 9712
pairs = pd.read_parquet(f"{LIT}/work/pairs_human.parquet")
print("pairs_human columns:", list(pairs.columns), pairs.shape, flush=True)
p = pairs[pairs["gene_id"] == GID].copy()
print("USP6NL rows in pairs_human:", len(p), flush=True)
pm = set(p["pmid"].astype(int))
cp = pq.read_table(f"{LIT}/work/corpus.parquet", columns=["pmid", "year", "title", "abstract"]).to_pandas()
cp = cp[cp["pmid"].isin(pm)].set_index("pmid")
# synonyms
with gzip.open(f"{LIT}/raw/ncbi/Homo_sapiens.gene_info.gz", "rt") as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if f[1] == str(GID):
            print("gene_info:", dict(zip(hdr, f))["Symbol"], "| synonyms:", dict(zip(hdr, f))["Synonyms"], "| description:", dict(zip(hdr, f))["description"])
            break

GREEK = str.maketrans({"α": "alpha", "β": "beta", "γ": "gamma", "δ": "delta", "ε": "epsilon", "κ": "kappa", "ω": "omega", "μ": "mu", "λ": "lambda",
                       "σ": "sigma", "θ": "theta", "ζ": "zeta", "η": "eta", "τ": "tau", "ρ": "rho", "π": "pi"})
NONALNUM = re.compile(r"[^0-9a-z]+")
def norm(s):
    return " " + NONALNUM.sub(" ", unicodedata.normalize("NFKC", s).lower().translate(GREEK)).strip() + " "

rows = []
for _, r in p.sort_values("pmid").iterrows():
    pmid = int(r["pmid"]); rec = cp.loc[pmid] if pmid in cp.index else None
    title = (rec["title"] or "") if rec is not None else ""; abstract = (rec["abstract"] or "") if rec is not None else ""
    text = (title + " " + abstract)
    nt = norm(text); nts = nt.replace(" ", "")
    ments = [m.strip() for m in set(str(r["mentions"]).split("|")) if m.strip()]
    hits = []
    for m in ments:
        nm = norm(m)
        if nm in nt or (len(m) >= 5 and nm.replace(" ", "") in nts):
            hits.append(m)
    # context: first sentence containing any hit (case-insensitive raw search)
    ctx = ""
    for m in hits or ments:
        i = text.lower().find(m.lower())
        if i >= 0:
            a = max(0, i - 160); b = min(len(text), i + len(m) + 160)
            ctx = ("..." if a else "") + text[a:b].replace("\n", " ") + ("..." if b < len(text) else "")
            break
    rows.append(dict(pmid=pmid, year=int(rec["year"]) if rec is not None and pd.notna(rec["year"]) else None, title=title[:140],
                     mentions="|".join(ments), verifier_accepts=bool(hits), matched=";".join(hits), resource=str(r["resource"])[:120], context=ctx))
df = pd.DataFrame(rows)
df.to_csv(f"{OUT}/usp6nl_mentions.csv", index=False)
print("\nverifier accepts", int(df["verifier_accepts"].sum()), "of", len(df), "USP6NL pairs\n")
for _, r in df.iterrows():
    print(f"PMID {r['pmid']} ({r['year']})  accepts={r['verifier_accepts']}  mentions={r['mentions']!r}")
    print(f"   title: {r['title']}")
    print(f"   context: {r['context']}\n")
