#!/bin/bash
#SBATCH -q express
#SBATCH --partition=erprp
#SBATCH --mem=10G
#SBATCH --time=1-00:00:00
#SBATCH -N 1-1
#SBATCH -n 1

cd $PWD

module load R/4.0.3

R CMD BATCH --vanilla "--args ${cell} ${lda} ${treat}" 1_normalize.R 
#normalize_${cell}.lda${lda}.trt${treat}.Rout
