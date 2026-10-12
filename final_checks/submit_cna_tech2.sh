#!/bin/bash
# Final check 2, second batch: which technical covariate drives the change (rin, nuc), a negative control (satperm), and the no-SAD variant with saturation.
cd ~/orcd/pool/mantis-films/research/broken-bricks/final_checks || exit 1
L=~/orcd/pool/mantis-films/research/seaad/logs/final
R=~/orcd/pool/mantis-films/research/seaad/results/final_checks/cna_tech
declare -A MEM=( [microglia]=12G [sst]=12G )
sub() {  # label model variant
  s=$1; m=$2; v=$3
  [ -s $R/${s}__${v}__${m}.json ] && { echo "skip $s $v $m (done)"; return; }
  while true; do
    out=$(sbatch -A mit_general -p mit_quicktest -t 00:15:00 --mem=${MEM[$s]} -c 2 -J fc-$s-$m-$v -o $L/cna_${s}_${m}_${v}-%j.log ../cna/run_py.sbatch ../final_checks/run_cna_tech.py $s $m --variant $v 2>&1)
    if echo "$out" | grep -q "QOSMaxSubmit"; then sleep 20; else echo "$s $m $v: $out"; break; fi
  done
}
for s in sst microglia; do
  for m in rin nuc satperm; do sub $s $m all; done
  for m in base sat; do sub $s $m nosad; done
done
echo "ALL SUBMITTED"
