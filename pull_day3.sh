#!/bin/bash
# pull Day-3 results from ORCD (run inside WSL Debian): bash pull_day3.sh
S=orcd:orcd/pool/mantis-films/research/seaad/results/day3
D=/mnt/c/mantis-films-work/research/broken-bricks/results/genetics
mkdir -p $D/magma
rsync -a --include="ot_overlap.csv" --include="h3_*.csv" --include="h3_summary.json" --include="program_summary.csv" --include="program_membership.csv" --include="programs.gmt" --include="specificity_table.csv" --include="specificity_meta.json" --include="apoe_region_genes_b37.csv" --exclude="*" $S/ $D/
rsync -a --include="*/" --include="*.gsa.out" --include="*.genes.out" --exclude="*" $S/magma/ $D/magma/
rsync -a $S/programs/ $D/programs/
rsync -a $S/sets/ $D/sets/
echo pulled
