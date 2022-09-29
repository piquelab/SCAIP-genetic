#!/bin/bash
#SBATCH -q primary
#SBATCH --mem=4G
#SBATCH --time=1-00:00:00
#SBATCH -N 1-1
#SBATCH -n 1
##SBATCH --array=1-30 

#j=${SLURM_ARRAY_TASK_ID}
module load misc
option='DiagLDA2'

if [ ! -d "${option}/2_fastQTL.test/" ]; then
   mkdir -p "${option}/2_fastQTL.test/"
fi

fastQTL \
  --vcf SCAIP1-6_filtered_mock_reheaded.vcf.gz \
  --bed ${option}/1_normalized.data/${cell}_lda${lda}_trt${treat}.bed.gz \
  --out ${option}/2_fastQTL.test/${cell}_lda${lda}_trt${treat}.nominals.chunk${j}.txt.gz --window 5e4 --chunk ${j} 30
