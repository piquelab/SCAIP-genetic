#!/bin/bash
#SBATCH -q express
#SBATCH -p erprp
#SBATCH --mem=60G
#SBATCH --time=1-00:00:00
#SBATCH -N 1-1
#SBATCH -n 1

cd $PWD

module load R/4.0.3
module load samtools
module load bcftools

R CMD BATCH  "--args ${cell} ${lda} ${treat}" 3_get.dosages.R 3_get.dosages.${cell}.lda${lda}.trt${treat}.Rout
#normalize_${cell}.lda${lda}.trt${treat}.Rout --vanilla
