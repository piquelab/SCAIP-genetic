# this file tests for genotype by LDA bin interaction using a linear model
# based on /wsu/home/groups/piquelab/SCAIP/LDA/bin_model/test_LDA-SNP.R
# 1/20/2021 JR
# 5/28/2021, JW
# last motified 11/08/2021, JW

library(tidyverse)
library(data.table)
library(parallel)
library(qvalue)
library(lme4)
library(lmtest)
library(qqman)

rm(list=ls())

# passing arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)>0){
  cell <- args[1]
  lda <- args[2]
  treat <- args[3]
}else{
  cell <- "Bcell"
  lda <- "LPS"
  treat <- "LPS-EtOH"
} 

option <- "DiagLDA2"
# set FDR threshold;
FDR <- 0.1

##############################
### load the needed files: ###
##############################

## normalized GE files:
fn <- paste0(option, "/1_normalized.data/", cell, "_lda", lda, "_trt", treat, ".bed.gz")
data <- fread(fn, header=T, sep="\t", stringsAsFactors=FALSE, data.table=F)

## txt file with dosages:
fn <- paste0(option, "/3_genotypes/dosages/1_eQTL_signif_dosages_", cell, "_lda", lda, "_trt", treat, ".txt.gz")
dosages <- fread(fn, sep="\t", stringsAsFactors=F, data.table=F)
rownames(dosages) <- dosages$varID
dosages <- dosages[,4:ncol(dosages)]


## txt file with list of testable pairs:
fn <- paste0(option, "/3_genotypes/eQTL_coordinates/1_", cell, "_lda", lda, "_trt", treat, "_pairs.txt")
pairs <- read.table(fn, sep="\t", header=T, stringsAsFactors=F)
# repeat the dosages as in the pairs file:
## dosages <- dosages[pairs$varID, 4:ncol(dosages)]

## keep only the samples present in GE data:
keep <- intersect(colnames(data), colnames(dosages))
## keep <- colnames(data)[colnames(data) %in% colnames(dosages)]
## dos <- dosages[,keep]
## data <- data[,c(colnames(data)[1:4],keep)]

# add the bin info:

### Final data list used for testing:
GE <- data[,keep]
rownames(GE) <- data$ID
## expression <- GE[pairs$ENSG,]
## colnames(expression) <- gsub("[.]","-",colnames(expression))
dosages <- dosages[, keep]
bin <- as.numeric(gsub(".*_","", colnames(GE)))

# try for just ctrl for now:
results <- mclapply(1:nrow(pairs), function(i){    
     # build a data frame for this gene*SNP:
     ## df <- matrix(nrow=length(c(rep(ctrl,ncol(dos)),rep(trt,ncol(dos)))), ncol=4, dimnames=list(rep("",length(c(rep(ctrl,ncol(dos)),rep(trt,ncol(dos))))),c("GE","dbgap","dosage","bin")))
     snp <- pairs[i,2]
     gene <- pairs[i,1]
    
     ## df <- cbind(data.frame(t(dosages[snp,])),data.frame(t(expression[i,])))
     df <- data.frame(dosage=as.numeric(dosages[snp,]),
                      GE=as.numeric(GE[gene,]),
                      bin=bin)

     ### testing
     m0 <- try( lm(GE~ bin+dosage, data=df), silent=TRUE)
     m1 <- try( lm(GE~ bin*dosage, data=df), silent=TRUE)
     ##
     if (  (class(m0)!="try-error")&(class(m1)!="try-error")){
       m1_summ <-  summary(m1)
     # report the results:
     # grab the tvalue:
       ## res <- data.frame(ENSG=pairs[i,1], varID=pairs[i,2],
       ##   "bin_effect"=m1_summ$coefficients[2,1],
       ##   "bin_SE"=m1_summ$coefficients[2,2],
       ##   "bin_tvalue"=m1_summ$coefficients[2,3],
       ##   "bin_pvalue"=m1_summ$coefficients[2,4],
       ##   "dosage_effect"=m1_summ$coefficients[3,1],
       ##   "dosage_SE"=m1_summ$coefficients[3,2],
       ##   "dosage_tvalue"=m1_summ$coefficients[3,3],
       ##   "dosage_pvalue"=m1_summ$coefficients[3,4],
       ##   "interaction_effect"=m1_summ$coefficients[4,1],
       ##   "interaction_SE"=m1_summ$coefficients[4,2],
       ##   "interaction_tvalue"=m1_summ$coefficients[4,3],
       ##   "interaction_pvalue"=m1_summ$coefficients[4,4],
       ##   "interaction_ANOVA.pvalue"=anova(m0,m1, test="LRT")[2,5],
       ##   "interaction_lrtest.pvalue"=lrtest(m0,m1)[2,5])
          res <- data.frame(ENSG=gene, varID=snp,
             interaction_ANOVA.pvalue=anova(m0, m1, test="LRT")[2,5]) 
       }else{
         res <- NA
       }
 }, mc.cores=10)

results <- results[!is.na(results)]
results <- do.call(rbind, results)

# multiple test correct:
results$interaction_ANOVA.qvalue <- qvalue(results$interaction_ANOVA.pvalue)$qvalues
## results$interaction_lrtest.qvalue <- qvalue(results$interaction_lrtest.pvalue)$qvalues

# save the table with the results:
outdir <- paste(option, "/4.2_lm.results/", sep="")
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
###
opfn <- paste0(outdir, "1_", cell, "_lda", lda, "_trt", treat, ".results.txt") 
write.table(results, file=opfn, sep="\t", col.names=T, row.names=F, quote=F)


### multiple test correct by stratified FDR
# load the eQTLs:
eqtls <- read.table("/wsu/home/groups/piquelab/SCAIP/SCAIP-genetic/eQTL/eQTL_coordinates/all_signif_SNP-gene_pairs.txt", sep="\t", header=T,stringsAsFactors=F)

# split by eQTL or not:
geneSNP <- eqtls$geneSNP
res <- results%>%
   mutate(gene_SNP=paste(ENSG, varID, sep="_"),
          geneSNP_is_eQTL=gene_SNP%in%geneSNP,
          geneSNP_is_eQTL2=ifelse(geneSNP_is_eQTL, "yes", "no"))

res2 <- res%>%group_by(geneSNP_is_eQTL)%>%
   mutate(interaction_ANOVA.qvalue_stratified=qvalue(interaction_ANOVA.pvalue)$qvalues)%>%ungroup()

###
opfn <- paste0(outdir, "2_", cell, "_lda", lda, "_trt", treat, "_stratified.FDR.txt")
write.table(res2, file=opfn, sep="\t", col.names=T, row.names=F, quote=F)




# report number of significant interactions:
#signif <- sum(results$interaction_ANOVA.qvalue<FDR,na.rm=T)
#egenes <- length(unique(results[results$interaction_ANOVA.qvalue<FDR & !is.na(results$interaction_ANOVA.qvalue),"ENSG"]))
#tab <- t(c(dataset, signif, egenes, nrow(df)))

#write.table(tab, file=paste0("./ANOVA_signif_interactions_lm_",lda,".txt"), sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE, append=TRUE)
## END 1/21/2021
