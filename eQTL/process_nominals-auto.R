# this script processes FastQTL results ran with different models (varying number of GE PCs removed) to find one maximizing the number of eGenes
# based on ../../eQTL/FastQTL/nominals/analysis_SCAIP1-6/process_nominals-auto.R
# 1/13/2021 JR

library(data.table)
library(qqman)
library(qvalue)

cluster <- "Tcell_CTRL"
pc = 0

args = commandArgs(trailingOnly=TRUE)

# get the PC number:
if (length(args)>0){
      pc <- args[1]
      cluster <- args[2]
      }


# set FDR threshold;
FDR <- 0.1
paste0("zcat ./eQTL_output/", cluster, ".GEPC", pc, ".nominals.eQTL.txt.gz")
pcx <- fread(paste0("zcat ./eQTL_output/", cluster, ".GEPC", pc, ".nominals.eQTL.txt.gz"), sep=" ")

pcx$V6 <- qvalue(pcx$V4)$qvalues
sum(pcx$V6<FDR)

## how many eQTLS:
eqtls <- sum(pcx$V6<FDR)
# how many egenes:
egenes <- length(unique(pcx$V1[pcx$V6<FDR]))

# write the number of eQTLs, egenes:
tab <- t(c(cluster, pc, eqtls, egenes))
write.table(tab, file=paste0("./nominals_eqtls-egenes_summary.txt"), sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE, append=TRUE)

# make QQ plots:
png(paste0("./plots/QQ_", cluster, "_PC1-", pc ,".png"))
temp <- paste0("pc", pc, "$V4")
qq(pcx$V4, main=paste0("QQ plot FastQTL analysis including ",pc, " PCs, ", cluster ))
dev.off()

### END 1/13/2021
