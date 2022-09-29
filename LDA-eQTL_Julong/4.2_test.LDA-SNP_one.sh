#!/bin/bash
#SBATCH -q express
#SBATCH --partition=erprp
#SBATCH --mem=60G
#SBATCH --time=1-00:00:00
#SBATCH -N 1-1
#SBATCH -n 10

cd $PWD

module load R/4.0.3

R CMD BATCH  "--args ${cell} ${lda} ${treat}" 4.2_test.LDA-SNP.R  test_${cell}.${lda}.${treat}.Rout 
#normalize_${cell}.lda${lda}.trt${treat}.Rout --vanilla
