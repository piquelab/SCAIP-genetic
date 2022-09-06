# this script processed fastQTL output data to extract list of eQTL coordinates and egenes for each cell type*condition combination with pre-defined number of GEPCs removed
# based on ../../eQTL/FastQTL/nominals/analysis_SCAIP1-6/geteQTLs.R
# 1/13/2021 JR

library(data.table)
library(qqman)
library(qvalue)
library(tidyr)

args <- commandArgs(trailingOnly=TRUE)

fdr_thresh <- 0.1

dataset <- "Bcell_LPS-EtOH"

if (length(args)>0){
    dataset <- args[1]
   PCs <- args[2]
 }
PCs <- "4"

results <- fread(paste0("zcat ./eQTL_output/", dataset, ".GEPC", PCs,".nominals.eQTL.txt.gz"), sep=" ")
results <- data.frame(results)

colnames(results) <- c("ENSG", "varID", "position", "pvalue", "slope_sc")

## results$bpadj <- p.adjust(results$bpval, method="fdr")
results$qvalue <- (qvalue(results$pvalue))$qvalues

signif <- results[results$qvalue<fdr_thresh,]
save <- separate(signif,2,into=c("chr", "position",NA,NA),sep=":", remove=FALSE)[,c("chr", "position")]

# save the egenes:
write.table(unique(signif$ENSG),  paste0("egenes/", dataset,"_egenes.txt"),sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

# save the SNP-gene pairs to be tested:
write.table(signif, paste0("eQTL_coordinates/signif_", dataset,".txt"),sep="\t", quote=FALSE, row.names=FALSE)

write.table(save, paste0("eQTL_coordinates/eQTL_coordinates_", dataset,".txt"),sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)

system(paste0("uniq eQTL_coordinates/eQTL_coordinates_", dataset,".txt eQTL_coordinates/eQTL_coordinates_uniq_", dataset,".txt"))

### END 1/13/2021
