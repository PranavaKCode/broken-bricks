#!/bin/bash
#SBATCH -J lit-dl-static
#SBATCH -A mit_general
#SBATCH -p mit_data_transfer
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 06:00:00
#SBATCH -o /home/pranavak/orcd/pool/mantis-films/research/lit/logs/dl-static-%j.log
set -uo pipefail
R=$HOME/orcd/pool/mantis-films/research/lit/raw
mkdir -p $R/ncbi $R/pubtator $R/ot_direct $R/uniprot; cd $R
date; hostname
dl() { # url out
  curl -fL --retry 8 --retry-delay 15 -C - -o "$2" "$1" && echo "OK $2 $(stat -c%s "$2") $(date -u +%FT%TZ)" || echo "FAIL $2"
}
dl https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz ncbi/Homo_sapiens.gene_info.gz
dl https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Mus_musculus.gene_info.gz ncbi/Mus_musculus.gene_info.gz
dl https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Rattus_norvegicus.gene_info.gz ncbi/Rattus_norvegicus.gene_info.gz
dl https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2pubmed.gz ncbi/gene2pubmed.gz
dl https://ftp.ncbi.nlm.nih.gov/pub/lu/PubTator3/gene2pubtator3.gz pubtator/gene2pubtator3.gz
# UniProt reviewed human, function text
UQ='https://rest.uniprot.org/uniprotkb/stream?compressed=true&fields=accession%2Cid%2Cgene_primary%2Cgene_names%2Cprotein_name%2Clength%2Ccc_function%2Cxref_geneid%2Cxref_ensembl%2Cxref_hgnc&format=tsv&query=%28reviewed%3Atrue+AND+organism_id%3A9606%29'
curl -fL --retry 8 --retry-delay 15 -o uniprot/uniprot_human_reviewed.tsv.gz "$UQ" && echo "OK uniprot $(stat -c%s uniprot/uniprot_human_reviewed.tsv.gz) $(date -u +%FT%TZ)" || echo "FAIL uniprot"
curl -fsI https://rest.uniprot.org/uniprotkb/search?query=organism_id:9606 | grep -i "x-uniprot-release" | tee uniprot/release.txt
# Open Targets bulk: association_by_datatype_direct (sensitivity S5)
OTB=https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/26.09/output/association_by_datatype_direct/
curl -fsL $OTB | grep -oE 'href="[^"]+\.parquet"' | sed 's/href="//;s/"$//' | sort -u > ot_direct/files.txt
wc -l ot_direct/files.txt
while read -r f; do dl "$OTB$f" "ot_direct/$f"; done < ot_direct/files.txt
sha256sum ncbi/* pubtator/* uniprot/*.gz > provenance_sha256.txt
ls -l ncbi pubtator uniprot ot_direct | head -50
date; echo DL_STATIC_DONE
