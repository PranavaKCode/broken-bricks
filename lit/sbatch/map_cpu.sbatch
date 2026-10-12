#!/bin/bash
#SBATCH -J lit-map
#SBATCH -A mit_general
#SBATCH -p mit_normal
#SBATCH -c 16
#SBATCH --mem=96G
#SBATCH -t 06:00:00
#SBATCH -o /home/pranavak/orcd/pool/mantis-films/research/lit/logs/map-%j.log
set -euo pipefail
export LD_LIBRARY_PATH=$HOME/orcd/pool/mantis-films/opt/conda/lit/lib
export NUMBA_NUM_THREADS=16 OMP_NUM_THREADS=16
PY=$HOME/orcd/pool/mantis-films/opt/conda/lit/bin/python
date; hostname
$PY -u $HOME/orcd/pool/mantis-films/research/lit/code/30_embed_map.py map
date; echo MAP_JOB_DONE
