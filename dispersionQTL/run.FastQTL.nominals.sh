#!/bin/bash

cat ../datasets.txt | tr -d '\015' |  while read var; do
for i in $(seq 0 10); do
for j in $(seq 1 30); do
    sbatch -q express -N1-1 -n 1 --mem=2G -t 1000 --job-name=$j$var --wrap "cd disp-eQTL_output/; fastQTL --vcf ../../../vcf/SCAIP1-6_filtered_AF.vcf.gz --bed ../qnormed_dispersion/${var}_dispersion.bed.gz --out $var.GEPC${i}.nominals.chunk$j.txt --cov ../covariates/${var}_covs_${i}PCs.txt --window 5e4  --chunk $j 30"
sleep 0.5
done;
done;
done
