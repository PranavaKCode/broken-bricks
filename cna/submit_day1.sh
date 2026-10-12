#!/bin/bash
# Submit the Day-1 CNA jobs (one job per subclass x variant), throttled for the mit_quicktest per-user submit limit.
# Run on the ORCD login node (light: it only sleeps and calls sbatch):
#   nohup bash submit_day1.sh all nosad ge50 > submit.log 2>&1 &
# Jobs that already have a result JSON, or are queued/running under the same job name, are skipped.
# mit_quicktest (15 min) fits every job; if a job times out, rerun it with PART=mit_preemptable TL=02:00:00.
cd ~/orcd/pool/mantis-films/research/broken-bricks/cna || exit 1
L=~/orcd/pool/mantis-films/research/seaad/logs/day1
R=~/orcd/pool/mantis-films/research/seaad/results/day1
PART=${PART:-mit_quicktest}; TL=${TL:-00:15:00}
VARIANTS=${@:-all nosad ge50}
declare -A MEM=( [microglia]=12G [astrocyte]=12G [oligodendrocyte]=16G [opc]=12G [sst]=12G [pvalb]=16G [vip]=16G [lamp5]=12G [l23it]=32G [l4it]=20G [l5it]=20G )
# subclasses where every donor has >= 50 nuclei: ge50 is identical to all, skip
declare -A SKIPGE50=( [microglia]=1 [astrocyte]=1 [l23it]=1 [l4it]=1 )
for v in $VARIANTS; do
  for s in microglia astrocyte oligodendrocyte opc sst pvalb vip lamp5 l23it l4it l5it; do
    [ "$v" = ge50 ] && [ -n "${SKIPGE50[$s]:-}" ] && continue
    [ -s $R/${s}__${v}.json ] && { echo "skip $s $v (done)"; continue; }
    [ -n "$(squeue -u $USER -h -n d1-$s-$v)" ] && { echo "skip $s $v (queued)"; continue; }
    while true; do
      out=$(sbatch -A mit_general -p $PART -t $TL --mem=${MEM[$s]} -c 2 -J d1-$s-$v -o $L/${s}_${v}-%j.log \
            run_py.sbatch run_cna.py $s $v 2>&1)
      if echo "$out" | grep -q "QOSMaxSubmit"; then sleep 20; else echo "$s $v: $out"; break; fi
    done
  done
done
echo "ALL SUBMITTED"
