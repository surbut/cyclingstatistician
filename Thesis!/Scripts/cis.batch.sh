#!/bin/sh

rm -rf lists_snps 
mkdir lists_snps;
seq -w 1 $(ls lists_genes/* | wc -l) | \
./parallel 'i={};
bedtools window -w 100000 -a lists_genes/list_genes_${i}.bed.gz 
-b /mnt/lustre/data/external_private/gtex_pilot/Imputation_SNP_Locs_95cr.bed.gz| \ 
cut -f10 | sort | uniq | gzip > lists_snps/list_snps_${i}.txt.gz'

