# this file runs permutations for genotype by treatment interaction test using a linear model
# based on ../../eQTL/FastQTL/nominals/reQTL/lm_SCAIP1-6/reQTL_lm_permutations_100.R
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

perms_number = 100
contrast <- paste0(cl, trt, ctrl)

# set FDR threshold;
FDR <- 0.1
## ncell=10

# load the needed files:
# 1. normalized GE files:
# treatment:
treat <- fread(paste0("../eQTL_mapping/normalized_GE_residuals/",cl, trt,".bed.gz"), header=T, sep="\t", stringsAsFactors=FALSE)
# control:
contr <- fread(paste0("../eQTL_mapping/normalized_GE_residuals/",cl, ctrl,".bed.gz"), header=T, sep="\t", stringsAsFactors=FALSE)
# 2. txt file with dosages:
dosages <- read.table(paste0("../eQTL_mapping/eQTL_coordinates/all_uniq_eQTL_dosages.txt"), sep="\t", header=T,stringsAsFactors=F)
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

# subset to common dbagps:
common_dbgaps <- colnames(expression_c)[colnames(expression_c) %in% colnames(expression_t)]
expression_c <- expression_c[,common_dbgaps]
expression_t <- expression_t[,common_dbgaps]
dos <- dos[,common_dbgaps]

# make a data frame that will take the results:
lm_ANOVA_pvals <- vector()

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
     # permute the interaction term using a design matrix:
     for(j in 1:perms_number){
X = model.matrix(~ dosage* treatment,data=df)
df$gxe <- X[,4]
df$permuted_gxe <- sample(df$gxe)
     lm0 <- lm(as.numeric(GE)~treatment + dosage, data=df)
     lm1 <- lm(as.numeric(GE)~treatment+dosage+permuted_gxe, data=df)
     ## m1_summ <-  summary(m1)
     # report the results:
     ## # grab the tvalue:
     ## results[i,"interaction_tvalue_lmer"] <- m1_summ$coefficients[4,3]
# run the ANOVA and grab the pvalue:
     lm_ANOVA_pvals <- append(lm_ANOVA_pvals,anova(lm0,lm1, test="LRT")[2,5])
     ## results[i,"interaction_lrtest.pvalue"] <- lrtest(m0,m1)[2,5]
      }
}


# save the permuted pvalues:
write.table(lm_ANOVA_pvals,file=paste0("reQTL_lm_results/permutations/",perms_number,"/",cl,"_", trt, "_vs", ctrl, ".txt"), sep="\n", col.names=F, row.names=F, quote=F)

## make a qqplot and p-value histogram:
png(paste0("./plots/QQplots/QQplot_ANOVA_lm_permutations_100_",cl,"_",ctrl, "_vs", trt, ".png"))
par(mfrow=c(1,2))
qq(lm_ANOVA_pvals)
dev.off()
sessionInfo()

## END 1/19/2021
