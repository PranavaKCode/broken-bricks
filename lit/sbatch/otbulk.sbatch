#!/bin/bash
#SBATCH -J lit-otbulk
#SBATCH -A mit_general
#SBATCH -p mit_quicktest
#SBATCH -t 00:10:00
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -o /home/pranavak/orcd/pool/mantis-films/research/lit/logs/otbulk-%j.log
set -euo pipefail
export LD_LIBRARY_PATH=$HOME/orcd/pool/mantis-films/opt/conda/lit/lib
PY=$HOME/orcd/pool/mantis-films/opt/conda/lit/bin/python
$PY -u $HOME/orcd/pool/mantis-films/research/lit/code/27_ot_bulk_check.py
echo OTBULK_JOB_DONE
