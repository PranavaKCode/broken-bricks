#!/bin/bash
# Final check 4 (part 2): MAGMA gene-property test of the pooled-neuron specificity covariates.
#   usage: sbatch ... magma_pooled.sbatch
# Needs results/final_checks/spec_covar_pooled.txt (specificity_pooled.py) and the Day-3 gene-level results results/day3/magma/<tag>.genes.raw.
# Tags: primary | neff | noapoe.  Marginal test of every covariate, Neuron_all given Microglia, Microglia given Neuron_all.
set -uo pipefail
P=/home/pranavak/orcd/pool/mantis-films
D=$P/research/seaad
M=$P/opt/magma/magma
R=$D/results/day3
F=$D/results/final_checks
mkdir -p $F/magma_pooled
date; hostname
for TAG in primary neff noapoe; do
  O=$F/magma_pooled/$TAG; mkdir -p $O
  G=$R/magma/$TAG.genes.raw
  $M --gene-results $G --gene-covar $F/spec_covar_pooled.txt --out $O/pooled 2>&1 | tail -2
  $M --gene-results $G --gene-covar $F/spec_covar_pooled.txt --model condition=N7_Microglia --out $O/neuron_given_microglia 2>&1 | tail -2
  $M --gene-results $G --gene-covar $F/spec_covar_pooled.txt --model condition=N7_Neuron_all --out $O/microglia_given_neuron 2>&1 | tail -2
done
echo "== primary marginal"
grep -v "^#" $F/magma_pooled/primary/pooled.gsa.out | head -30
ls -la $F/magma_pooled/primary | head -20
date; echo EXIT 0
