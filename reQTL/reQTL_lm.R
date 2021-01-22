# this file tests for genotype by treatment interaction using a linear model
# relies on ../eQTL_mapping
# based on  ../../eQTL/FastQTL/nominals/reQTL/lm_SCAIP1-6/reQTL_lm.R
# 1/19/2021 JR

library(tidyr)
library(qvalue)
library(data.table)
library(qvalue)
library(lme4)
library(lmtest)
library(qqman)

cl <- "Bcell"
trt <- "_LPS-EtOH"
ctrl <- "_CTRL"

args = commandArgs(trailingOnly=TRUE)

# get the PC number:
if (length(args)>0){
      ## pc <- args[1]
      cl <- args[1]
      trt <- args[2]
      ctrl <- args[3]
      ## folder <- args[3]
      }

contrast <- paste0(cl, trt, ctrl)

# set FDR threshold;
FDR <- 0.1

# load the needed files:
# 1. normalized GE files:
# treatment:
treat <- fread(paste0("../eQTL_mapping/normalized_GE_residuals/",cl, trt,".bed.gz"), header=T, sep="\t", stringsAsFactors=FALSE)
# control:
contr <- fread(paste0("../eQTL_mapping/normalized_GE_residuals/",cl, ctrl,".bed.gz"), header=T, sep="\t", stringsAsFactors=FALSE)
# 2. txt file with dosages:
dosages <- read.table(paste0("../eQTL_mapping//eQTL_coordinates/all_uniq_eQTL_dosages.txt"), sep="\t", header=T,stringsAsFactors=F)
colnames(dosages) <- gsub("[.]","-",colnames(dosages))
# remove repeated dosages (where do they come from, though??):
dosages <- dosages[!duplicated(dosages$varID),]
rownames(dosages) <- dosages$varID
# 3. txt file with list of testable pairs:
pairs <- read.table(paste0("../eQTL_mapping/eQTL_coordinates/all_signif_SNP-gene_pairs.txt"),sep="\t", header=T,stringsAsFactors=F)
# remove the genes absent from eith ctrl or trt GE files:
test_pairs <- pairs[pairs$ENSG %in% treat$ID & pairs$ENSG %in% contr$ID,]
# repeat the dosages as in the pairs file:
dos <- dosages[test_pairs$varID,5:ncol(dosages)]

### TEST:
GE_c <- data.frame(contr[,5:ncol(contr)])
colnames(GE_c) <- gsub("[.]","-",colnames(GE_c))
rownames(GE_c) <- contr$ID
expression_c <- GE_c[test_pairs$ENSG,]
GE_t <- data.frame(treat[,5:ncol(treat)])
colnames(GE_t) <- gsub("[.]","-",colnames(GE_t))
rownames(GE_t) <- treat$ID
expression_t <- GE_t[test_pairs$ENSG,]

# subset to common individuals
common_dbgaps <- colnames(expression_c)[colnames(expression_c) %in% colnames(expression_t)]
expression_c <- expression_c[,common_dbgaps]
expression_t <- expression_t[,common_dbgaps]
dos <- dos[,common_dbgaps]

# make a data frame that will take the results:
results <- test_pairs[,1:2]
results$treatment_effect <- NA
results$treatment_SE <- NA
results$treatment_tvalue <- NA
results$dosage_effect <- NA
results$dosage_SE <- NA
results$dosage_tvalue <- NA
results$interaction_effect <- NA
results$interaction_SE <- NA
results$interaction_tvalue <- NA
results$interaction_ANOVA.pvalue <- NA
results$interaction_lrtest.pvalue <- NA

# try for just ctrl for now:
 for (i in 1:nrow(test_pairs)){
     # build a data frame for this gene*SNP:
     ## df <- matrix(nrow=length(c(rep(ctrl,ncol(dos)),rep(trt,ncol(dos)))), ncol=4, dimnames=list(rep("",length(c(rep(ctrl,ncol(dos)),rep(trt,ncol(dos))))),c("GE","dbgap","dosage","treatment")))
     df <- data.frame(matrix(nrow=length(c(rep(ctrl,ncol(dos)),rep(trt,ncol(dos)))), ncol=0))
     df <- cbind(df,t(cbind(expression_c[i,],expression_t[i,])))
     df <- cbind(df,rep(colnames(dos),2))
     df <- cbind(df,unlist(rep(dos[i,],2)))
     df <- cbind(df,c(rep(ctrl,ncol(dos)),rep(trt,ncol(dos))))
     colnames(df) <- c("GE","dbgap","dosage","treatment")   
     m0 <- lm(as.numeric(GE)~ dosage+treatment, data=df)
     m1 <- lm(as.numeric(GE)~ dosage*treatment, data=df)
     m1_summ <-  summary(m1)
     # report the results:
     # grab the tvalue:
     results[i,"treatment_effect"] <- m1_summ$coefficients[2,1]
     results[i,"treatment_SE"] <- m1_summ$coefficients[2,2]
     results[i,"treatment_tvalue"] <- m1_summ$coefficients[2,3]
     results[i,"dosage_effect"] <- m1_summ$coefficients[3,1]
     results[i,"dosage_SE"] <- m1_summ$coefficients[3,2]
     results[i,"dosage_tvalue"] <- m1_summ$coefficients[3,3]
     results[i,"interaction_effect"] <- m1_summ$coefficients[4,1]
     results[i,"interaction_SE"] <- m1_summ$coefficients[4,2]
     results[i,"interaction_tvalue"] <- m1_summ$coefficients[4,3]
# run the ANOVA and grab the pvalue:
     results[i,"interaction_ANOVA.pvalue"] <- anova(m0,m1, test="LRT")[2,5]
     results[i,"interaction_lrtest.pvalue"] <- lrtest(m0,m1)[2,5]
      }

# multiple test correct:
results$interaction_ANOVA.qvalue <- qvalue(results$interaction_ANOVA.pvalue)$qvalues
## results$interaction_lrtest.qvalue <- qvalue(results$interaction_lrtest.pvalue)$qvalues

# save the table with the results:
write.table(results,file=paste0("reQTL_lm_results/results_",cl, trt, "_vs", ctrl, ".txt"), sep="\t", col.names=T, row.names=F, quote=F)

# report number of significant interactions:
signif <- sum(results$interaction_ANOVA.qvalue<FDR)
egenes <- length(unique(results[results$interaction_ANOVA.qvalue<FDR,"ENSG"]))
tab <- t(c(cl,ctrl, trt, signif, egenes, length(common_dbgaps)))

write.table(tab, file="./ANOVA_signif_interactions_lm.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE, append=TRUE)

## make a qqplot and p-value histogram:
png(paste0("./plots/QQplots/QQplot_ANOVA_",cl,"_",ctrl, "_vs", trt, ".png"))
qq(results$interaction_ANOVA.pvalue)
dev.off()
png(paste0("./plots/QQplots/QQplot_lrtest_",cl,"_",ctrl, "_vs", trt, ".png"))
qq(results$interaction_lrtest.pvalue)
dev.off()
png(paste0("./plots/t.test/tvalue_histogram_",cl,"_",ctrl,"_vs",trt,".png"))
hist(results$interaction_tvalue)
dev.off()

sessionInfo()

## END 1/19/2021
