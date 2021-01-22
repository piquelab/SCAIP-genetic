#!/bin/bash
lda=$1
cat ../datasets.txt | tr -d '\015' |  while read var; do
for j in $(seq 1 30); do
    sbatch -q primary -n 1 --mem=4G -t 1000 --wrap "fastQTL --vcf genotypes/SCAIP1-6_filtered_mock_reheaded.vcf.gz --bed ./normalized_data/${lda}/$var.bed.gz --out ./tests/${lda}/$var.nominals.chunk$j.txt --window 5e4 --chunk $j 30"
done;
done
