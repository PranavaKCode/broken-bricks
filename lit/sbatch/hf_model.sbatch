#!/bin/bash
#SBATCH -J lit-hf
#SBATCH -A mit_general
#SBATCH -p mit_data_transfer
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 01:00:00
#SBATCH -o /home/pranavak/orcd/pool/mantis-films/research/lit/logs/hf-%j.log
set -uo pipefail
export LD_LIBRARY_PATH=$HOME/orcd/pool/mantis-films/opt/conda/lit/lib
B=$HOME/orcd/pool/mantis-films
export HF_HOME=$B/opt/hf
$B/opt/conda/lit/bin/python - <<'PY'
from huggingface_hub import snapshot_download
p = snapshot_download("NeuML/pubmedbert-base-embeddings")
print(p)
import os
for r,d,f in os.walk(p):
    for x in f: print(os.path.join(r,x), os.path.getsize(os.path.join(r,x)))
PY
echo HF_DONE
