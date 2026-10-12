#!/bin/bash
# MAGMA gene analysis in NB parallel batches, then merge.   usage: sbatch ... magma_gene.sbatch TAG PFILE [NB]
#   TAG   output prefix (primary | neff | noapoe)
#   PFILE P/N file (SNP, P, N) in results/
# Writes results/day3/magma/<TAG>.genes.out and <TAG>.genes.raw
set -uo pipefail
TAG=$1; PF=$2; NB=${3:-6}
P=/home/pranavak/orcd/pool/mantis-films
D=$P/research/seaad
M=$P/opt/magma/magma
O=$D/results/day3/magma; mkdir -p $O
W=${TMPDIR:-/tmp}/magma-$SLURM_JOB_ID; mkdir -p $W
date; hostname; echo "tag=$TAG pfile=$PF batches=$NB"
ANN=$O/ad.genes.annot
if [ ! -s $ANN ]; then
  $M --annotate window=10,10 --snp-loc $D/data/magma/g1000_eur/g1000_eur.bim --gene-loc $D/data/magma/geneloc37/NCBI37.3.gene.loc --out $O/ad 2>&1 | tail -5
fi
for i in $(seq 1 $NB); do
  $M --bfile $D/data/magma/g1000_eur/g1000_eur --pval $D/results/$PF use=SNP,P ncol=N --gene-annot $ANN --gene-model snp-wise=mean \
     --batch $i $NB --out $W/$TAG > $W/$TAG.b$i.log 2>&1 &
done
wait
ls -la $W | head -40
tail -3 $W/$TAG.b1.log
$M --merge $W/$TAG --out $W/$TAG 2>&1 | tail -5
ls -la $W | head -40
cp $W/$TAG.genes.out $W/$TAG.genes.raw $O/ 2>&1
wc -l $O/$TAG.genes.out
date; echo EXIT 0
