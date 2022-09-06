# this script fits mash on all the preprocessed SCAIP FastQTL eQTL mapping on mean (from NB) results and generates the input for mashr (pvalues,SEs)
# based on  ../mashr_eQTL/mashr_prep.R
# 2/15/2021 JR

library(ashr)
library(data.table)

dataset <- "Bcell_CTRL"

args = commandArgs(trailingOnly=TRUE)
if (length(args)>0){
      dataset <- args[1]
    }

pcs = 8

eQTL_dir <- "../eQTL-on-mean/mean-eQTL_output/"

# read in, grab effect estimate and its standard error
m <- fread(paste0(eQTL_dir, dataset,".GEPC",pcs,".nominals.eQTL.txt.gz"),sep=" ")
colnames(m) <- c("ENSG","varID", "distance", "pvalue", "slope")
# add unique SNP identifier:
m$uniqID <- paste0(m$ENSG,"_",m$varID)
# calculate SE:
m$z.abs=abs(qnorm(m$pvalue/2))
m$SE = abs(m$slope)/m$z.abs
# calculate lfsr:
ashr=ash(m$slope,m$SE)
m$lfsr <- ashr$result$lfsr
# save the slope and the SE in separate files:
slope <- m[,c("uniqID","slope")]
colnames(slope)[2] <- dataset
write.table(slope, paste0("input/", dataset, "_slope.txt"), sep="\t", col.names=T, row.names=F, quote=F)
SE <- m[,c("uniqID","SE")]
colnames(SE)[2] <- dataset
write.table(SE, paste0("input/", dataset, "_SE.txt"), sep="\t", col.names=T, row.names=F, quote=F)
# also save pvalues to select "strong" signals for estimating covariances step:
pvals <- m[,c("uniqID","pvalue")]
colnames(pvals)[2] <- dataset
# save lfsr to use instead of pvalues
write.table(pvals, paste0("input/", dataset, "_pvalue.txt"), sep="\t", col.names=T, row.names=F, quote=F)
lfsr <- m[,c("uniqID","lfsr")]
colnames(lfsr)[2] <- dataset
write.table(lfsr, paste0("input/", dataset, "_lfsr.txt"), sep="\t", col.names=T, row.names=F, quote=F)

# gzip all:
system(paste0("for file in input/", dataset,"*txt ; do gzip $file; done"))

### END 2/15/2021 JR
