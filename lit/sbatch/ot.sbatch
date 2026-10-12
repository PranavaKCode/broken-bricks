#!/bin/bash
#SBATCH -J lit-ot
#SBATCH -A mit_general
#SBATCH -p mit_data_transfer
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 02:00:00
#SBATCH -o /home/pranavak/orcd/pool/mantis-films/research/lit/logs/ot-%j.log
set -uo pipefail
export LD_LIBRARY_PATH=$HOME/orcd/pool/mantis-films/opt/conda/lit/lib
PY=$HOME/orcd/pool/mantis-films/opt/conda/lit/bin/python
C=$HOME/orcd/pool/mantis-films/research/lit/code
date; hostname
$PY -u $C/ot_gwas_fetch.py ot
date; echo OT_JOB_DONE
