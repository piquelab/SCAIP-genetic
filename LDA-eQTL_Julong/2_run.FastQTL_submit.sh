#!/bin/bash

cd $PWD

cat dataset_contrast.txt | #sed -n '1p' | 
while read -r cell lda treat; do
  for j in $(seq 1 30); do
  sbatch --export=cell=${cell},lda=${lda},treat=${treat},j=${j} --out=slurm.fastQTL.${cell}.${lda}.${treat}.${j}.output 2_run.FastQTL_one.sh
  done  
sleep 30s 
done
