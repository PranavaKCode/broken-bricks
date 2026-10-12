#!/bin/bash
#SBATCH -J lit-agora-c
#SBATCH -A mit_general
#SBATCH -p mit_data_transfer
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 03:00:00
#SBATCH -o /home/pranavak/orcd/pool/mantis-films/research/lit/logs/agora-c-%j.log
set -uo pipefail
export LD_LIBRARY_PATH=$HOME/orcd/pool/mantis-films/opt/conda/lit/lib
PY=$HOME/orcd/pool/mantis-films/opt/conda/lit/bin/python
L=$HOME/orcd/pool/mantis-films/research/lit
export AGORA_UNIV=$L/work/agora_universe_primary.csv AGORA_JL=agora_genes_c.jsonl AGORA_SKIP=agora_genes.jsonl,agora_genes_b.jsonl AGORA_INTERVAL=0.3
date; hostname
$PY -u $L/code/agora_fetch.py
date; echo AGORA_C_DONE
