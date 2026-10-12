#!/bin/bash
# Submit remaining one-off jobs as submit slots free up (QOS limit on queued jobs): l23it nosad NAM, specificity, MAGMA neff / noapoe.
P=~/orcd/pool/mantis-films/research
cd $P/broken-bricks/genetics || exit 1
sub() {  # name, then sbatch args
  name=$1; shift
  while true; do
    out=$(sbatch "$@" 2>&1)
    if echo "$out" | grep -q -E "QOSMax|Submit"; then
      if echo "$out" | grep -q "Submitted"; then echo "$name: $out"; break; fi
      sleep 15
    else echo "$name: $out"; break; fi
  done
}
sub specificity -A mit_general -p mit_quicktest -t 00:15:00 --mem=32G -c 2 -J d3-spec -o $P/seaad/logs/day3/spec-%j.log ../cna/run_py.sbatch specificity.py
sub l23it_nosad -A mit_general -p mit_quicktest -t 00:15:00 --mem=48G -c 2 -J d2-l23it-nosad -o $P/seaad/logs/day2/l23it_nosad-%j.log ../cna/run_py.sbatch ../bricks/day2_nam.py l23it nosad
sub magma_neff -A mit_general -p mit_quicktest -t 00:15:00 --mem=48G -c 8 -J d3-magma-neff -o $P/seaad/logs/day3/magma_neff-%j.log magma_gene.sbatch neff bellenguez_rsid_p_neff.txt 8
sub magma_noapoe -A mit_general -p mit_quicktest -t 00:15:00 --mem=48G -c 8 -J d3-magma-noapoe -o $P/seaad/logs/day3/magma_noapoe-%j.log magma_gene.sbatch noapoe bellenguez_rsid_p_n.noAPOE.txt 8
echo ALL SUBMITTED
