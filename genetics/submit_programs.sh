#!/bin/bash
# Submit Day-3 program jobs (one per subclass x variant) on mit_quicktest.   nohup bash submit_programs.sh all nosad > submit3.log 2>&1 &
cd ~/orcd/pool/mantis-films/research/broken-bricks/genetics || exit 1
L=~/orcd/pool/mantis-films/research/seaad/logs/day3
R=~/orcd/pool/mantis-films/research/seaad/results/day3/programs
mkdir -p $L $R
PART=${PART:-mit_quicktest}; TL=${TL:-00:15:00}
VARIANTS=${@:-all nosad}
declare -A MEM=( [microglia]=24G [astrocyte]=24G [oligodendrocyte]=32G [opc]=24G [sst]=24G [pvalb]=32G [vip]=32G [lamp5]=24G [l23it]=64G [l4it]=48G [l5it]=48G )
for v in $VARIANTS; do
  for s in microglia astrocyte oligodendrocyte opc sst pvalb vip lamp5 l23it l4it l5it; do
    [ -s $R/${s}__${v}_genes.csv ] && { echo "skip $s $v (done)"; continue; }
    [ -n "$(squeue -u $USER -h -n d3p-$s-$v)" ] && { echo "skip $s $v (queued)"; continue; }
    while true; do
      out=$(sbatch -A mit_general -p $PART -t $TL --mem=${MEM[$s]} -c 2 -J d3p-$s-$v -o $L/prog_${s}_${v}-%j.log \
            ../cna/run_py.sbatch programs.py $s $v 2>&1)
      if echo "$out" | grep -q -E "QOSMaxSubmit|QOSMax"; then sleep 20; else echo "$s $v: $out"; break; fi
    done
  done
done
echo "ALL SUBMITTED"
