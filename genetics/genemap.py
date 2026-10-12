"""Shared gene-identifier helpers for Day 3 (amendment A7 item 2).

NCBI Gene IDs from `Homo_sapiens.gene_info.gz` with the A3 rule: a record whose symbol equals the source symbol first, then protein-coding,
then the lowest Gene ID. Runs on ORCD (reads the literature study's copy of gene_info).
"""
import gzip, os
import pandas as pd, numpy as np

ROOT = "/home/pranavak/orcd/pool/mantis-films/research"
GENE_INFO = f"{ROOT}/lit/raw/ncbi/Homo_sapiens.gene_info.gz"
GENELOC37 = f"{ROOT}/seaad/data/magma/geneloc37/NCBI37.3.gene.loc"

def load_gene_info():
    d = pd.read_csv(GENE_INFO, sep="\t", usecols=["GeneID", "Symbol", "Synonyms", "dbXrefs", "type_of_gene"], dtype=str)
    d["GeneID"] = d["GeneID"].astype(int)
    ens = d["dbXrefs"].fillna("").str.extract(r"Ensembl:(ENSG\d+)")[0]
    d["ensembl"] = ens
    return d

class GeneMap:
    def __init__(self):
        d = load_gene_info()
        self.info = d.set_index("GeneID")
        e = d.dropna(subset=["ensembl"])
        self.by_ens = {k: g for k, g in e.groupby("ensembl")}
        self.by_sym = {k: g for k, g in d.groupby("Symbol")}
        loc = pd.read_csv(GENELOC37, sep="\t", header=None, names=["GeneID", "chr", "start", "stop", "strand", "symbol"])
        self.loc = loc.set_index("GeneID")

    def map_one(self, ens, sym):
        cand = None
        if ens in self.by_ens: cand = self.by_ens[ens]
        elif sym in self.by_sym: cand = self.by_sym[sym]
        if cand is None: return (np.nan, "")
        c = cand.assign(symeq=(cand["Symbol"] == sym).astype(int), pc=(cand["type_of_gene"] == "protein-coding").astype(int))
        c = c.sort_values(["symeq", "pc", "GeneID"], ascending=[False, False, True])
        r = c.iloc[0]
        return (int(r["GeneID"]), r["type_of_gene"])

    def map_table(self, df, ens_col="ensembl_id", sym_col="symbol"):
        out = [self.map_one(e, s) for e, s in zip(df[ens_col].values, df[sym_col].values)]
        df = df.copy()
        df["entrez"] = [o[0] for o in out]; df["type_of_gene"] = [o[1] for o in out]
        df["in_magma_loc"] = df["entrez"].isin(self.loc.index)
        return df

    def apoe_region_genes(self, lo=44_000_000, hi=46_000_000):
        l = self.loc
        m = (l["chr"].astype(str) == "19") & (l["stop"] >= lo) & (l["start"] <= hi)
        return set(l.index[m])
