#!/bin/bash
# MAGMA gene-property (cell-type specificity) and competitive set analyses for one gene-analysis tag.
#   usage: sbatch ... magma_sets.sbatch TAG        TAG = primary | neff | noapoe
# Outputs in results/day3/magma/<TAG>/ : celltype.gsa.out, <label>.marginal.gsa.out, <label>.conditional.gsa.out
set -uo pipefail
TAG=$1
P=/home/pranavak/orcd/pool/mantis-films
D=$P/research/seaad
M=$P/opt/magma/magma
R=$D/results/day3
O=$R/magma/$TAG; mkdir -p $O
G=$R/magma/$TAG.genes.raw
date; hostname; echo "tag=$TAG"
declare -A NAME=( [microglia]=Microglia [astrocyte]=Astrocyte [oligodendrocyte]=Oligodendrocyte [opc]=OPC [sst]=Sst [pvalb]=Pvalb [vip]=Vip [lamp5]=Lamp5 [l23it]=L23IT [l4it]=L4IT [l5it]=L5IT )
echo "== gene-property: cell-type specificity"
$M --gene-results $G --gene-covar $R/spec_covar.txt --out $O/celltype 2>&1 | tail -4
head -12 $O/celltype.gsa.out
for L in microglia astrocyte oligodendrocyte opc sst pvalb vip lamp5 l23it l4it l5it; do
  echo "== programs $L"
  $M --gene-results $G --set-annot $R/sets/$L.sets.txt --out $O/$L.marginal 2>&1 | grep -E "ERROR|WARNING|Error" | head -3
  $M --gene-results $G --set-annot $R/sets/$L.sets.txt --gene-covar $R/spec_covar.txt --model condition=${NAME[$L]} --out $O/$L.conditional 2>&1 | grep -E "ERROR|WARNING|Error" | head -3
done
ls $O | head -50
date; echo EXIT 0
