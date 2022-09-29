#!/bin/bash/

cd $PWD

cat dataset_contrast.txt | 
while read -r cell lda treat; do
  sbatch --export=cell=${cell},lda=${lda},treat=${treat} --output=slurm.${cell}.${lda}.${treat}.output 1_normalize_one.sh 
  #echo ${cell} ${lda} ${treat}
  sleep 1
done
