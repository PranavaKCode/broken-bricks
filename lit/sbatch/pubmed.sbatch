#!/bin/bash
#SBATCH -J lit-pubmed
#SBATCH -A mit_general
#SBATCH -p mit_data_transfer
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 12:00:00
#SBATCH -o /home/pranavak/orcd/pool/mantis-films/research/lit/logs/pubmed-%j.log
set -uo pipefail
export LD_LIBRARY_PATH=$HOME/orcd/pool/mantis-films/opt/conda/lit/lib
PY=$HOME/orcd/pool/mantis-films/opt/conda/lit/bin/python
C=$HOME/orcd/pool/mantis-films/research/lit/code
date; hostname
$PY -u $C/pubmed_fetch.py
# GWAS year lookups use the same NCBI rate budget, so run strictly after the PubMed download
$PY -u $C/ot_gwas_fetch.py gwas
date; echo PUBMED_JOB_DONE
