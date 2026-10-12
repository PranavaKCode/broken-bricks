#!/bin/bash
#SBATCH -J lit-stageA
#SBATCH -A mit_general
#SBATCH -p mit_normal
#SBATCH -c 4
#SBATCH --mem=48G
#SBATCH -t 02:00:00
#SBATCH -o /home/pranavak/orcd/pool/mantis-films/research/lit/logs/stageA-%j.log
set -euo pipefail
export LD_LIBRARY_PATH=$HOME/orcd/pool/mantis-films/opt/conda/lit/lib
PY=$HOME/orcd/pool/mantis-films/opt/conda/lit/bin/python
L=$HOME/orcd/pool/mantis-films/research/lit
date; hostname
# AD PMID list -> awk filter of gene2pubtator3 (no header; PMID is column 1)
$PY - <<'PYEOF'
import polars as pl, os
L=os.path.expanduser("~/orcd/pool/mantis-films/research/lit")
pl.read_parquet(f"{L}/raw/pubmed/abstracts.parquet", columns=["pmid"]).write_csv(f"{L}/work/ad_pmids.txt", include_header=False)
PYEOF
if [ ! -s $L/raw/pubtator/gene2pubtator3_AD.tsv ]; then
  zcat $L/raw/pubtator/gene2pubtator3.gz | awk -F'\t' 'NR==FNR{a[$1];next} ($1 in a)' $L/work/ad_pmids.txt - > $L/raw/pubtator/gene2pubtator3_AD.tsv
fi
wc -l $L/raw/pubtator/gene2pubtator3_AD.tsv; head -3 $L/raw/pubtator/gene2pubtator3_AD.tsv
$PY -u $L/code/10_attention.py
$PY -u $L/code/15_universe.py
date; echo STAGEA_JOB_DONE
