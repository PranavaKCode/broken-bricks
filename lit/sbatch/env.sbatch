#!/bin/bash
#SBATCH -J lit-env
#SBATCH -A mit_general
#SBATCH -p mit_data_transfer
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 03:00:00
#SBATCH -o /home/pranavak/orcd/pool/mantis-films/research/lit/logs/env-%j.log
set -euo pipefail
module load miniforge/26.7.2-0
B=$HOME/orcd/pool/mantis-films
export CONDA_PKGS_DIRS=$B/opt/pkgs PIP_CACHE_DIR=$B/opt/pipcache
ENV=$B/opt/conda/lit
date; hostname
if [ ! -x $ENV/bin/python ]; then
  conda create -y -p $ENV -c conda-forge --override-channels python=3.12 pip numpy scipy pandas polars pyarrow requests lxml \
    matplotlib seaborn scikit-learn statsmodels leidenalg python-igraph tqdm numba pynndescent umap-learn
fi
$ENV/bin/pip install --no-cache-dir torch --index-url https://download.pytorch.org/whl/cu126
$ENV/bin/pip install --no-cache-dir sentence-transformers
$ENV/bin/python - <<'PY'
import sys, importlib
print(sys.version)
for m in ["torch","sentence_transformers","umap","pynndescent","leidenalg","igraph","polars","pyarrow","requests","matplotlib","seaborn","sklearn","statsmodels","lxml","scipy","numba"]:
    mod = importlib.import_module(m); print(m, getattr(mod,"__version__","?"))
import torch; print("torch cuda build", torch.version.cuda, "available", torch.cuda.is_available())
PY
rm -rf $B/opt/pkgs $B/opt/pipcache
date; echo ENV_DONE
