#!/bin/bash
#SBATCH -J lit-audit
#SBATCH -A mit_general
#SBATCH -p mit_quicktest
#SBATCH -t 00:10:00
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -o /home/pranavak/orcd/pool/mantis-films/research/lit/logs/audit-%j.log
set -euo pipefail
export LD_LIBRARY_PATH=$HOME/orcd/pool/mantis-films/opt/conda/lit/lib
PY=$HOME/orcd/pool/mantis-films/opt/conda/lit/bin/python
$PY -u $HOME/orcd/pool/mantis-films/research/lit/code/25_ner_audit.py
echo AUDIT_JOB_DONE
