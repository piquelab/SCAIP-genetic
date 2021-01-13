# SCAIP-genetic
1/6/2021 JR <br/>
This directory /wsu/home/groups/piquelab/SCAIP/SCAIP1-6_protein-coding contains all the scripts and results of SCAIP1-6 genetic analyses on protein-coding genes only <br/>

### eQTL mapping
./eQTL_mapping - all the scripts and results of eQTL mapping (and follow-up analyses) on SCAIP1-6 pseudo-bulk GE data<br/>
strategy: FastQTL on pseudo-bulk residuals <br/>
INPUT: /nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/6_DEG.CelltypeNew_output/Filter2/YtX_sel.comb.RData <br/>
OUTPUT: ./eQTL_mapping/eQTL_output/*.GEPC0.nominals.eQTL.txt.gz <br/>
###
Filters:<br/>
- inheritted from /nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/6_DEG.CelltypeNew.R<br/>
- eQTL mapping window: +/-50kb<br/>
- MAF: >=10% in cohort
- gene: 0.1 CPM in more than 20% of the samples<br/>
###
Steps:<br/>
1. Normalize GE data and extract residuals: ./eQTL/normalize-all.R
2. Calculate GE PCs (to correct for in eQTL mapping): ./eQTL_mapping/GE_PCA.sh
3. Run eQTL mapping: ./eQTL/run.FastQTL.nominals-all-covs.sh -> ./eQTL_mapping/concatenate.chunks.sh
4. Check removing how many GE PCs generates the most egenes: ./eQTL_mapping/process-nominals.sh
5. Get egenes and eQTL coordinates from optimal model results: ./eQTL_mapping/get-egenes-eQTLs.sh
### mashr_eQTL
./mashr_eQTL - all the scripts and results of running mashr (and follow-up analyses) on SCAIP1-6; relies on output in ./eQTL_mapping <br/>
strategy: mashr on eQTL results from all treatments and conditions at once <br/>
INPUT: ./eQTL_mapping/eQTL_output/*.GEPC0.nominals.eQTL.txt.gz  <br/>
OUTPUT:  <br/>
###
Filters:<br/>
- gene-variant pairs with estimates across all cell types and conditions
###
Steps:<br/>
1. Generate the files needed for mashr (SEs, pvalues, lfsr): ./mashr_eQTL/mashr_prep.sh
2. Run mashr on full eQTL results and save the model fit: ./mashr_eQTL/mashr-fit-model.R
3. Run mashr on chunked data using pre-saved model fit: ./mashr_eQTL/mashr_compute_posterior_chunks.sh
4. Make upset plot of mashr results across conditions: ./mashr_eQTL/plot_upset_mashr.R

### reQTL mapping
./reQTL_mapping - all the scripts and results of reQTL mapping (and follow-up analyses) on SCAIP1-6; relies on output in ./eQTL_mapping <br/>
strategy: lm testing dosage*treatment interaction using pair-wise trt-control models on union of significant eQTLs <br/>
INPUT:  <br/>
OUTPUT:  <br/>
###
Filters:<br/>
- tested only significant (10%FDR) eQTLs
###
Steps:<br/>
### eQTL mapinng on mean GE
./eQTL_mapinng-on-mean - all the scripts and results of eQTL mapping (and follow-up analyses) on SCAIP1-6 mean values from NB model
strategy: FastQTL on residuals extracted from mean values from NB model <br/>
INPUT: /nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/10_RNA.Variance_output/tmp9/1.2_Sel.Bx.RData <br/>
OUTPUT: ./eQTL_mapinng-on-mean/eQTL_output/*.GEPC0.nominals.eQTL.txt.gz <br/>
###
Filters:<br/>
- inheritted from /nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/10_RNA.variance.R<br/>
- min. 3 individuals per batch-condition<br/>
- eQTL mapping window: +/-50kb<br/>
- MAF: >=10% in cohort
### dispersion eQTL mapping
./dispersionQTL - all the scripts and results of running FastQTL (and follow-up analyses) on SCAIP1-6 dispersion data <br/>
strategy: FastQTL on dispersion residuals <br/>
INPUT: /nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/10_RNA.Variance_output/tmp9/1.2_Sel.PhxNew.RData <br/>
OUTPUT: dispersionQTL/disp-eQTL_output/*.GEPC0.nominals.eQTL.txt.gz <br/>
###
Filters:<br/>
- inheritted from /nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/10_RNA.variance.R<br/>
- min. 3 individuals per batch-condition<br/>
- eQTL mapping window: +/-50kb<br/>
- MAF: >=10% in cohort
###
Steps:<br/>
1. Quantile-normalize dispersion measure and extract residuals: dispersionQTL/get_residuals.R
2. Run dispersion-eQTL mapping: ./dispersionQTL/run.FastQTL.nominals.sh
3. Get the number of eQTLs, egenes per condition and save the egenes: ./dispersionQTL/process-nominals.sh

### LDA eQTL mapping
strategy: lm testing dosage*treatment interaction using GE bulked along 3 bins <br/>
INPUT: /nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/9_RNA.dynamic2_output/Filter2_DEG6571/Old/LDA{1,2}Bin/YtX.*.ave.RData <br/>
OUTPUT:  <br/>
###
Filters:<br/>
- inheritted from /nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/9_RNA.dynamic2.R
- eQTL mapping window: +/-50kb<br/>
- MAF: >=10% in cohort
- gene: 0.1 CPM in more than 20% of the samples<br/>
###
Steps:<br/>
