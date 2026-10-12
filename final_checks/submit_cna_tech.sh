#!/bin/bash
# Final check 2: submit the CNA reruns with technical covariates (one job per subclass x model), throttled for the mit_quicktest submit limit.
# Run on the ORCD login node (light: sleeps and calls sbatch):   nohup bash submit_cna_tech.sh > submit_cna_tech.log 2>&1 &
# Order: Sst and Microglia first (base control, sat, tech3, neun), then the other nine subclasses (sat, tech3).
cd ~/orcd/pool/mantis-films/research/broken-bricks/final_checks || exit 1
L=~/orcd/pool/mantis-films/research/seaad/logs/final
R=~/orcd/pool/mantis-films/research/seaad/results/final_checks/cna_tech
mkdir -p $L $R
declare -A MEM=( [microglia]=12G [astrocyte]=12G [oligodendrocyte]=16G [opc]=12G [sst]=12G [pvalb]=16G [vip]=16G [lamp5]=12G [l23it]=32G [l4it]=20G [l5it]=20G )
sub() {  # label model
  s=$1; m=$2
  [ -s $R/${s}__all__${m}.json ] && { echo "skip $s $m (done)"; return; }
  [ -n "$(squeue -u $USER -h -n fc-$s-$m)" ] && { echo "skip $s $m (queued)"; return; }
  while true; do
    out=$(sbatch -A mit_general -p mit_quicktest -t 00:15:00 --mem=${MEM[$s]} -c 2 -J fc-$s-$m -o $L/cna_${s}_${m}-%j.log ../cna/run_py.sbatch ../final_checks/run_cna_tech.py $s $m 2>&1)
    if echo "$out" | grep -q "QOSMaxSubmit"; then sleep 20; else echo "$s $m: $out"; break; fi
  done
}
for s in sst microglia; do for m in base sat tech3 neun; do sub $s $m; done; done
for m in sat tech3; do
  for s in astrocyte oligodendrocyte opc pvalb vip lamp5 l23it l4it l5it; do sub $s $m; done
done
echo "ALL SUBMITTED"
