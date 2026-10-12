#!/bin/bash
# generic wrapper (copied from the Day-0 setup): sbatch ... run_py.sbatch script.py [args]
set -uo pipefail
P=/home/pranavak/orcd/pool/mantis-films
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-2} NUMBA_NUM_THREADS=${SLURM_CPUS_PER_TASK:-2} MPLBACKEND=Agg
export NUMBA_CACHE_DIR=/tmp/numba-$SLURM_JOB_ID; mkdir -p $NUMBA_CACHE_DIR
date; hostname
echo "script: $@"
$P/opt/conda/sc/bin/python -I "$@"
rc=$?
date; echo "EXIT $rc"
exit $rc
