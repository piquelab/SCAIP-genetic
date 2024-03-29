#!/bin/bash

cat datasets.txt | tr -d '\015' |  while read var; do
for j in $(seq 1 30); do
    sbatch -q express -N1-1 -n 1 --mem=2G -t 1000 --job-name=$var --wrap "cd eQTL_output/; fastQTL --vcf ../../../vcf/SCAIP1-6_filtered_AF.vcf.gz --bed ../normalized_GE_residuals/$var.bed.gz --out $var.GEPC0.nominals.chunk$j.txt.gz --window 5e4  --chunk $j 30"
sleep 0.5
done;
done;


cat datasets.txt | tr -d '\015' |  while read var; do
for i in $(seq 1 10); do
for j in $(seq 1 30); do
    sbatch -q express -N1-1 -n 1 --mem=2G -t 1000 --job-name=$var --wrap "cd eQTL_output/; fastQTL --vcf ../../../vcf/SCAIP1-6_filtered_AF.vcf.gz --bed ../normalized_GE_residuals/$var.bed.gz --out $var.GEPC$i.nominals.chunk$j.txt.gz --window 5e4 --cov ../GE_PCs/${var}_${i}GEPCs.txt --chunk $j 30"
sleep 0.5
done;
done;
done
