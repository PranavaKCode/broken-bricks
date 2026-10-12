#!/bin/bash
#SBATCH -J lit-leiden
#SBATCH -A mit_general
#SBATCH -p mit_quicktest
#SBATCH -t 00:15:00
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -o /home/pranavak/orcd/pool/mantis-films/research/lit/logs/leiden-%j.log
set -euo pipefail
export LD_LIBRARY_PATH=$HOME/orcd/pool/mantis-films/opt/conda/lit/lib
PY=$HOME/orcd/pool/mantis-films/opt/conda/lit/bin/python
date; hostname
$PY -u $HOME/orcd/pool/mantis-films/research/lit/code/30_embed_map.py leiden
date; echo LEIDEN_JOB_DONE
