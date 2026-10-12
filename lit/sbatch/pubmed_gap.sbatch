#!/bin/bash
#SBATCH -J lit-pubmed-gap
#SBATCH -A mit_general
#SBATCH -p mit_data_transfer
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 04:00:00
#SBATCH -o /home/pranavak/orcd/pool/mantis-films/research/lit/logs/pubmed-gap-%j.log
set -uo pipefail
export LD_LIBRARY_PATH=$HOME/orcd/pool/mantis-films/opt/conda/lit/lib
PY=$HOME/orcd/pool/mantis-films/opt/conda/lit/bin/python
date; hostname
$PY -u $HOME/orcd/pool/mantis-films/research/lit/code/pubmed_gapfill.py
date; echo PUBMED_GAP_DONE
