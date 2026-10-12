#!/bin/bash
# Submit the Day-2 NAM jobs (one per subclass x variant), throttled for the mit_quicktest per-user submit limit.
#   nohup bash submit_day2.sh all nosad > submit2.log 2>&1 &
cd ~/orcd/pool/mantis-films/research/broken-bricks/bricks || exit 1
L=~/orcd/pool/mantis-films/research/seaad/logs/day2
R=~/orcd/pool/mantis-films/research/seaad/results/day2
PART=${PART:-mit_quicktest}; TL=${TL:-00:15:00}
VARIANTS=${@:-all nosad}
declare -A MEM=( [microglia]=16G [astrocyte]=16G [oligodendrocyte]=24G [opc]=16G [sst]=16G [pvalb]=24G [vip]=24G [lamp5]=16G [l23it]=48G [l4it]=32G [l5it]=32G )
for v in $VARIANTS; do
  for s in microglia astrocyte oligodendrocyte opc sst pvalb vip lamp5 l23it l4it l5it; do
    [ -s $R/${s}__${v}.npz ] && { echo "skip $s $v (done)"; continue; }
    [ -n "$(squeue -u $USER -h -n d2-$s-$v)" ] && { echo "skip $s $v (queued)"; continue; }
    while true; do
      out=$(sbatch -A mit_general -p $PART -t $TL --mem=${MEM[$s]} -c 2 -J d2-$s-$v -o $L/${s}_${v}-%j.log \
            ../cna/run_py.sbatch day2_nam.py $s $v 2>&1)
      if echo "$out" | grep -q "QOSMaxSubmit"; then sleep 20; else echo "$s $v: $out"; break; fi
    done
  done
done
echo "ALL SUBMITTED"
