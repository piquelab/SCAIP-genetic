# this script compares effect sizes between sc (SCAIP) and bulk data (ALOFT) and makes scatterplots
# based on /wsu/home/groups/piquelab/SCAIP/eQTL/FastQTL/nominals/analysis_SCAIP1-6/compare_to_bulk.R 
# 1/22/2021 JR

library(data.table)
library(qqman)
library(qvalue)
library(tidyverse)
library(ggplot2)

cluster <- "Bcell_CTRL"

args = commandArgs(trailingOnly=TRUE)

# get the PC number:
if (length(args)>0){
      pc <- args[1]
      cluster <- args[2]
      }

pc <- 3

trt <- gsub(".*_","",cluster)
# rename trts to their proper names:
trt <- gsub("-EtOH","",trt)
trt <- gsub("-","+",trt)

cell <- gsub("_.*","",cluster)
if(cell %in% c("Bcell", "Tcell","NKcell")){
cell <- paste(gsub('[a-z]','', cell), gsub('[A-Z]','', cell))
}
cell <- paste0(cell,"s")

# set FDR threshold;
FDR <- 0.1
paste0("zcat ./eQTL_output/", cluster, ".GEPC", pc, ".nominals.eQTL.txt.gz")
pcx <- fread(paste0("zcat ./eQTL_output/", cluster, ".GEPC", pc, ".nominals.eQTL.txt.gz"), sep=" ")

# colnames:
colnames(pcx) <- c("ENSG", "varID", "position", "pvalue", "slope_sc")

pcx$qvalue <- qvalue(pcx$pvalue)$qvalues
sum(pcx$V6<FDR)
## length(unique(pcx$V1[pcx$qvalues<FDR]))
pcx[pcx$qvalue<0.1,"significant_sc"] <- "yes"
pcx[pcx$qvalue>0.1,"significant_sc"] <- "no"
pcx$pair <- paste0(gsub("[.].*","",pcx$ENSG), "_", pcx$varID)

# read in the bulk ALOFT results:
bulk <- fread("/nfs/rprdata/ALOFT/AL1-6_ln/FastQTL/permutations/results/FastQTL_results_best_18GEPCs.txt", sep="\t", header=T)
bulk[bulk$bqval<0.1, "significant_bulk"] <- "yes"
bulk[bulk$bqval>0.1, "significant_bulk"] <- "no"
bulk$pair <- paste0(bulk$pid, "_", bulk$sid)

# inner-merge the data:
merged <- inner_join(bulk,pcx)
merged[merged$qvalue<0.1 & merged$bqval>0.1,"significance"] <- "single-cell only"
merged[merged$qvalue>0.1 & merged$bqval<0.1,"significance"] <- "bulk only"
merged[merged$qvalue>0.1 & merged$bqval>0.1,"significance"] <- "neither"
merged[merged$qvalue<0.1 & merged$bqval<0.1,"significance"] <- "both"

# sort to keep significant on top:
merged$significance <- factor(merged$significance, levels=c("neither", "both", "bulk only", "sc only"))

cor=cor.test(merged$slope, merged$slope_sc)$estimate
pval=cor.test(merged$slope, merged$slope_sc)$p.value

merged_sig <- filter(merged, !significance=="neither")
cor_sig=cor.test(merged_sig$slope, merged_sig$slope_sc)$estimate
pval_sig=cor.test(merged_sig$slope, merged_sig$slope_sc)$p.value

# save the correlation:
tab <- t(c(cluster, cor, pval, cor_sig, pval_sig))
write.table(tab, "bulk_correlations_Spearman.txt", sep="\t", append=T, col.names=F, row.names=F, quote=F)

# report eGenes at 10% FDR, eQTLs at 10%FDR, p<0.05:
egenes10 <- unique(pcx[pcx$qvalue <0.1, "ENSG"])
eqtls10 <- merged[merged$qvalue <0.1, c("ENSG","varID")]
eqtls05 <- merged[merged$pvalue <0.05, c("ENSG","varID")]

## write.table(egenes10, "SCAIP_eQTLs_egenes10.txt", sep="\t", append=T, col.names=F, row.names=F, quote=F)
write.table(eqtls10, "SCAIP_eQTLs_qval10.txt", sep="\t", append=T, col.names=F, row.names=F, quote=F)
## write.table(eqtls05, "SCAIP_eQTLs_pval05.txt", sep="\t", append=T, col.names=F, row.names=F, quote=F)

# change order of plotting (neither, bulk, both, sc):
## merged$significance <- factor(merged$significance, levels=c("neither","bulk only","single-cell only", "both"))
# sort the df:
mergedplot <- merged[merged$significance=="neither",]
mergedplot <- rbind(mergedplot, merged[merged$significance=="bulk only",])
mergedplot <- rbind(mergedplot, merged[merged$significance=="single-cell only",])
mergedplot <- rbind(mergedplot, merged[merged$significance=="both",])

# make scatterplots:
pdf(paste0("./plots/compare-bulk/bulk-vs-", cluster, ".pdf"))
gg0 <- ggplot(mergedplot, aes(slope,slope_sc, color=significance)) + ggtitle(paste0("rho=", round(cor,2), " pvalue=", round(pval,4))) + labs(x="Genetic effect size in bulk RNA-seq", y=paste0("Genetic effect size in ",trt, " condition ",cell, " from scRNA-seq")) +
 geom_vline(xintercept=0) + geom_hline(yintercept=0) +
## , cor_sig=", round(cor_sig,2), " pval_sig=", round(pval_sig,2)))
 geom_point()+ scale_color_manual(values=c("neither"="grey", "both"="black","single-cell only"="green", "bulk only"="blue")) + theme_bw()
gg0
dev.off()
ggsave(paste0("./plots/compare-bulk/bulk-vs-", cluster, ".png"),gg0)

### END 1/22/2021
