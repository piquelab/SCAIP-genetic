# this file permutation-correct reQTL p-values
# based on  ../../eQTL/FastQTL/nominals/reQTL/lm_SCAIP1-6/reQTL_lm_permutation-correct.100
# 1/19/2021 JR

library(tidyr)
library(qvalue)
library(data.table)
library(qvalue)
library(lme4)
library(lmtest)
library(qqman)
library(ggplot2)

cl <- "Bcell"
trt <- "_PHA-DEX"
ctrl <- "_PHA-EtOH"

args = commandArgs(trailingOnly=TRUE)

# get the PC number:
if (length(args)>0){
      ## pc <- args[1]
      cl <- args[1]
      trt <- args[2]
      ctrl <- args[3]
      ## folder <- args[3]
      }

perms_number <- 100

contrast <- paste0(cl, trt, ctrl)

# set FDR threshold;
FDR <- 0.1

# load the needed files:
# results:
res <- read.table(paste0("reQTL_lm_results/results_",cl, trt, "_vs", ctrl, ".txt"), sep="\t", header=T, stringsAsFactors=FALSE)
# permutation pvalues:
pvals <- read.table(paste0("reQTL_lm_results/permutations/100/",cl,"_", trt, "_vs", ctrl, ".txt"), sep="\t", stringsAsFactors=FALSE)

# get all pvalues:
pvalues <- pvals$V1

# correct p-values using Roger's method:
obspval <- res$interaction_ANOVA.pvalue

permpval = pvalues;
mesh = data.frame(pval=c(obspval,permpval),
cInd=c(rep(FALSE,length(obspval)),
rep(TRUE,length(permpval)))
                                    )

o <- order(mesh$pval)
cInd <- mesh$cInd
fp <- cumsum(cInd[o])
pncc <- (1:length(fp))-fp
fdr <-  (fp/sum(cInd)*sum(!cInd))/pncc
corr.pval <- fp
corr.pval[fp==0] <- 1;
corr.pval <- corr.pval/sum(cInd)
mesh$pcorr <- 1
mesh$pcorr[o] <- corr.pval
mesh$qval <- 1
qvcorr <- qvalue(mesh$pcorr[!cInd])
mesh$qval[!cInd] <- qvcorr$qvalues

mesh2 <- mesh[!cInd,]
mesh3 <- cbind(res,mesh2[,"pcorr"])
colnames(mesh3)[colnames(mesh3)=="mesh2[, \"pcorr\"]"] <- "interaction_A_lm_pcorr"

mesh3$interaction_A_lm_pcorr.qval <- qvalue(mesh3$interaction_A_lm_pcorr)$qvalues

eqtls <- sum(mesh3$interaction_A_lm_pcorr.qval<0.1)
egenes <- length(unique(mesh3[mesh3$interaction_A_lm_pcorr.qval<0.1,"ENSG"]))

tab <- t(c(contrast, eqtls,egenes))
write.table(tab, paste0("ANOVA_signif_interactions_corrected_lm_",perms_number,".txt"),sep="\t",append=T, col.names=F, row.names=F, quote=F)

# save the results:
write.table(mesh3, paste0("reQTL_lm_results/",contrast,"_corrected_",perms_number,".txt"),sep="\t", row.names=F, quote=F)

## make a qqplot and p-value histogram:
## qq(mesh3$interaction_ANOVA.pvalue)
## qq(mesh3$interaction_A_lmer_pcorr)
## dev.off()

permpoints <- -log10(sort(pvalues))
expected_permpoints <- -log10(ppoints(length(permpoints)))
perms <- data.frame(cbind(permpoints,expected_permpoints))

ci <- 0.95
## QQ plot:
     N  <- length(mesh3$interaction_ANOVA.pvalue)
      df <- data.frame(
                   observed = -log10(sort(mesh3$interaction_ANOVA.pvalue)),
                   corrected = -log10(sort(mesh3$interaction_A_lm_pcorr)),
                   expected = -log10(ppoints(N)),
          ## expected_perm = -log10(ppoints(length(permpoints))),
                  clower   = -log10(qbeta((1 - ci)/2, 1:N, N - 1:N+1)),
                   cupper   = -log10(qbeta((1 + ci)/2, 1:N, N - 1:N+1))
                              )

nas <- matrix(nrow=(nrow(perms)-nrow(df)), ncol=ncol(df))
colnames(nas) <- colnames(df)
df1 <- rbind(df, nas)
df <- cbind(df1,  perms)
      log10Pe <- expression(paste("Expected -log"[10], plain(P)))
      log10Po <- expression(paste("Observed -log"[10], plain(P)))

png(paste0("./plots/QQplots/QQplot_A_lm_pcorr_",cl,"_",ctrl, "_vs", trt, "_",perms_number,".png"))
ggplot(df) +
           theme_bw(base_size = 11, base_family = "") +
           ggrastr::geom_point_rast(aes(expected, observed,color='a'), shape=1, size=1)+
           ggrastr::geom_point_rast(aes(expected_permpoints, permpoints,color="b"), shape=1,size=1, alpha=0.4) +
           ggrastr::geom_point_rast(aes(expected, corrected, color="d"), shape=1, size=1)+
             geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
             geom_line(aes(expected, cupper), linetype = 2) +
             geom_line(aes(expected, clower), linetype = 2) +
                                       ggtitle(paste0(contrast)) +
                                       xlab(log10Pe) +
                                       theme(legend.position='right')+
                                       ylab(log10Po)+
      scale_colour_manual(name = 'color',
      values =c('a'='black','b'='slategrey',"d"="lightblue"), labels = c('observed pvalues', 'permuted pvalues', 'corrected pvalues'),guide = 'legend')
dev.off()


sessionInfo()

## END 1/19/2021
