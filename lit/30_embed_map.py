#!/usr/bin/env python
"""Stage C: PubMedBERT embeddings of AD abstracts and UniProt function text in one space (GPU),
then PCA-64, kNN (cosine), UMAP-2d, Leiden clusters and TF-IDF cluster labels (CPU). No LLM calls. Exploratory (preregistration 2.2).

usage: 30_embed_map.py embed   (GPU node)
       30_embed_map.py map     (CPU node; needs the embeddings)
Outputs under work/embed/.
"""
import os, re, sys, json, time, gzip, io, datetime
from pathlib import Path
import numpy as np
import polars as pl

SEED = 20261011
STAGE = sys.argv[1] if len(sys.argv) > 1 else "all"
LIT = Path(os.path.expanduser("~/orcd/pool/mantis-films/research/lit"))
RAW, WORK = LIT / "raw", LIT / "work"
OUT = WORK / "embed"; OUT.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("HF_HOME", os.path.expanduser("~/orcd/pool/mantis-films/opt/hf"))
os.environ["TOKENIZERS_PARALLELISM"] = "true"
t00 = time.time()
def log(*a):
    print(f"[{time.time() - t00:7.0f}s]", *a, flush=True)


def clean(f):
    if not f: return ""
    f = re.sub(r"\{ECO:[^}]*\}", "", f)
    f = re.sub(r"\((?:PubMed|By similarity)[^)]*\)", "", f)
    f = re.sub(r"^\s*FUNCTION:\s*", "", f)
    return re.sub(r"\s+", " ", f).strip()


def load_texts():
    corp = pl.read_parquet(RAW / "pubmed" / "abstracts.parquet").sort("pmid")
    paper_text = (corp["title"].fill_null("") + " " + corp["abstract"].fill_null("")).str.strip_chars().to_list()
    up = pl.read_csv(gzip.open(RAW / "uniprot" / "uniprot_human_reviewed.tsv.gz").read(), separator="\t", quote_char=None, infer_schema_length=0)
    up = up.rename({"Entry": "uniprot", "Gene Names (primary)": "symbol_up", "Protein names": "protein_name", "Function [CC]": "function", "GeneID": "geneid_raw"})
    up = up.with_columns(pl.col("geneid_raw").fill_null("").str.extract(r"^(\d+)", 1).cast(pl.Int64, strict=False).alias("gene_id"))
    hum = pl.read_parquet(WORK / "human_genes.parquet")
    up = up.drop_nulls("gene_id").filter(pl.col("gene_id").is_in(hum["gene_id"])).unique(subset="gene_id", keep="first", maintain_order=True)
    ft = [clean(f) for f in up["function"].to_list()]
    gtext = [t if t else f"{(n or '')}. Gene {(s or '')}." for t, n, s in zip(ft, up["protein_name"].fill_null("").to_list(), up["symbol_up"].fill_null("").to_list())]
    genes = up.select(["gene_id", "uniprot", "symbol_up", "protein_name"]).with_columns(
        pl.Series("text_source", ["uniprot_function" if t else "fallback_name" for t in ft]), pl.Series("text", gtext))
    return corp, paper_text, genes, gtext


def embed():
    corp, paper_text, genes, gtext = load_texts()
    genes.write_parquet(OUT / "gene_texts.parquet")
    log("papers", corp.height, "genes", genes.height, "with function text", int((genes["text_source"] == "uniprot_function").sum()))
    import torch
    from sentence_transformers import SentenceTransformer
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    log("device", dev, torch.cuda.get_device_name(0) if dev == "cuda" else "")
    model = SentenceTransformer("NeuML/pubmedbert-base-embeddings", device=dev)
    model.max_seq_length = 256
    if dev == "cuda":
        model.half()
    def enc(texts, name):
        f = OUT / f"emb_{name}_fp16.npy"
        if f.exists():
            return np.load(f)
        t = time.time()
        e = model.encode(texts, batch_size=256, convert_to_numpy=True, normalize_embeddings=False, show_progress_bar=False).astype(np.float16)
        np.save(f, e)
        log(f"encoded {name} {e.shape} in {time.time() - t:.0f}s")
        return e
    enc(gtext, "genes")
    enc(paper_text, "papers")
    ids = pl.concat([pl.DataFrame(dict(kind="paper", id=corp["pmid"].to_list())), pl.DataFrame(dict(kind="gene", id=genes["gene_id"].to_list()))])
    ids.write_parquet(OUT / "ids_union.parquet")
    json.dump(dict(device=dev, gpu=torch.cuda.get_device_name(0) if dev == "cuda" else None, model="NeuML/pubmedbert-base-embeddings", max_seq_length=256, fp16=dev == "cuda",
                   n_papers=corp.height, n_genes=genes.height, finished=datetime.datetime.utcnow().isoformat() + "Z"), open(OUT / "embed_meta.json", "w"), indent=2)
    log("EMBED_DONE")


def mapstage():
    corp, paper_text, genes, gtext = load_texts()
    Ep = np.load(OUT / "emb_papers_fp16.npy"); Eg = np.load(OUT / "emb_genes_fp16.npy")
    assert Ep.shape[0] == corp.height and Eg.shape[0] == genes.height
    from sklearn.decomposition import PCA
    X = np.vstack([Ep, Eg]).astype(np.float32)
    pca = PCA(n_components=64, svd_solver="randomized", random_state=SEED)
    Z = pca.fit_transform(X).astype(np.float32)
    np.save(OUT / "pca64.npy", Z)
    np.savez(OUT / "pca_model.npz", components=pca.components_, mean=pca.mean_, evr=pca.explained_variance_ratio_)
    log("PCA done; explained variance (64 PCs)", float(pca.explained_variance_ratio_.sum()))
    del X
    import pynndescent
    t = time.time()
    nn = pynndescent.NNDescent(Z, n_neighbors=31, metric="cosine", random_state=SEED, n_jobs=-1, low_memory=True)
    idx, dist = nn.neighbor_graph
    log(f"kNN done {time.time() - t:.0f}s")
    np.save(OUT / "knn_idx.npy", idx.astype(np.int32)); np.save(OUT / "knn_dist.npy", dist.astype(np.float32))
    import umap
    t = time.time()
    reducer = umap.UMAP(n_neighbors=30, min_dist=0.05, metric="cosine", n_components=2, random_state=SEED, precomputed_knn=(idx, dist, nn), verbose=True)
    U = reducer.fit_transform(Z).astype(np.float32)
    np.save(OUT / "umap2.npy", U)
    log(f"UMAP done {time.time() - t:.0f}s")
    import igraph as ig, leidenalg
    n = Z.shape[0]; k = 30
    rows = np.repeat(np.arange(n), k); cols = idx[:, 1:k + 1].reshape(-1)
    wts = np.clip(1.0 - dist[:, 1:k + 1].reshape(-1), 1e-3, 1.0)
    a = np.minimum(rows, cols); b = np.maximum(rows, cols)
    key = a.astype(np.int64) * n + b
    order = np.argsort(key, kind="stable")
    key, a, b, wts = key[order], a[order], b[order], wts[order]
    first = np.r_[True, key[1:] != key[:-1]]
    grp = np.cumsum(first) - 1
    wmax = np.zeros(grp[-1] + 1); np.maximum.at(wmax, grp, wts)
    g = ig.Graph(n=n, edges=list(zip(a[first].tolist(), b[first].tolist())), directed=False)
    g.es["weight"] = wmax.tolist()
    log(f"graph {g.vcount()} nodes {g.ecount()} edges")
    RES = float(os.environ.get("LEIDEN_RES", "1.0"))
    part = leidenalg.find_partition(g, leidenalg.RBConfigurationVertexPartition, weights="weight", resolution_parameter=RES, n_iterations=4, seed=SEED)
    lab = np.array(part.membership, dtype=np.int32)
    log(f"Leiden resolution {RES}: {lab.max() + 1} clusters; modularity {part.modularity:.3f}")
    np.save(OUT / "leiden_labels.npy", lab)
    from sklearn.feature_extraction.text import TfidfVectorizer
    all_text = paper_text + gtext
    rng = np.random.default_rng(SEED)
    docs = []
    for c in range(lab.max() + 1):
        members = np.where(lab == c)[0]
        take = members if len(members) <= 400 else rng.choice(members, 400, replace=False)
        docs.append(" ".join(all_text[i] for i in take))
    vec = TfidfVectorizer(stop_words="english", ngram_range=(1, 2), min_df=2, max_df=0.3, sublinear_tf=True, max_features=300000,
                          token_pattern=r"(?u)\b[A-Za-z][A-Za-z0-9\-]{2,}\b")
    M = vec.fit_transform(docs)
    terms = np.array(vec.get_feature_names_out())
    npap = np.bincount(lab[:corp.height], minlength=lab.max() + 1); ngen = np.bincount(lab[corp.height:], minlength=lab.max() + 1)
    rows_ = []
    for c in range(lab.max() + 1):
        r = M.getrow(c).toarray().ravel(); top = terms[np.argsort(-r)[:8]]
        rows_.append(dict(cluster=c, n_papers=int(npap[c]), n_genes=int(ngen[c]), label=", ".join(top[:4]), keywords="; ".join(top)))
    pl.DataFrame(rows_).write_csv(OUT / "cluster_labels.csv")
    json.dump(dict(finished=datetime.datetime.utcnow().isoformat() + "Z", n=int(n), n_papers=int(corp.height), n_genes=int(genes.height), leiden_resolution=RES,
                   n_clusters=int(lab.max() + 1), seed=SEED, pca_evr=float(pca.explained_variance_ratio_.sum())), open(OUT / "run_meta.json", "w"), indent=2)
    log("MAP_DONE")


if __name__ == "__main__":
    if STAGE in ("embed", "all"):
        embed()
    if STAGE in ("map", "all"):
        mapstage()
