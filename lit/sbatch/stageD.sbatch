#!/bin/bash
#SBATCH -J lit-stageD
#SBATCH -A mit_general
#SBATCH -p mit_quicktest
#SBATCH -t 00:15:00
#SBATCH -c 4
#SBATCH --mem=48G
#SBATCH -o /home/pranavak/orcd/pool/mantis-films/research/lit/logs/stageD-%j.log
set -euo pipefail
export LD_LIBRARY_PATH=$HOME/orcd/pool/mantis-films/opt/conda/lit/lib
PY=$HOME/orcd/pool/mantis-films/opt/conda/lit/bin/python
C=$HOME/orcd/pool/mantis-films/research/lit/code
date; hostname
$PY -u $C/35_assemble.py
$PY -u $C/40_figures.py
date; echo STAGED_JOB_DONE
