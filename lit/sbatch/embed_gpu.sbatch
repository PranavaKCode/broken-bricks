#!/bin/bash
#SBATCH -J lit-embed
#SBATCH -A mit_general
#SBATCH -p mit_normal_gpu
#SBATCH --gres=gpu:1
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 02:00:00
#SBATCH -o /home/pranavak/orcd/pool/mantis-films/research/lit/logs/embed-%j.log
set -euo pipefail
export LD_LIBRARY_PATH=$HOME/orcd/pool/mantis-films/opt/conda/lit/lib
export HF_HOME=$HOME/orcd/pool/mantis-films/opt/hf HF_HUB_OFFLINE=1
PY=$HOME/orcd/pool/mantis-films/opt/conda/lit/bin/python
date; hostname; nvidia-smi -L
$PY -u $HOME/orcd/pool/mantis-films/research/lit/code/30_embed_map.py embed
date; echo EMBED_JOB_DONE
