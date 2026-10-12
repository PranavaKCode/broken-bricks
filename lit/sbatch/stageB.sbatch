#!/bin/bash
#SBATCH -J lit-stageB
#SBATCH -A mit_general
#SBATCH -p mit_normal
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 03:00:00
#SBATCH -o /home/pranavak/orcd/pool/mantis-films/research/lit/logs/stageB-%j.log
set -euo pipefail
export LD_LIBRARY_PATH=$HOME/orcd/pool/mantis-films/opt/conda/lit/lib
PY=$HOME/orcd/pool/mantis-films/opt/conda/lit/bin/python
date; hostname
$PY -u $HOME/orcd/pool/mantis-films/research/lit/code/20_metrics.py
date; echo STAGEB_JOB_DONE
