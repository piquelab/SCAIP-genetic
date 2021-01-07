# this script takes in a normalized expression bed file, runs PCA on it and saves the PCs
# based on /wsu/home/groups/piquelab/SCAIP/eQTL/FastQTL/GE_PCA/GE_PCA.R
# 1/7/2021 JR

library(dplyr)
library(data.table)

data_dir <- "./normalized_GE_residuals/"

dataset <- 'Tcell_PHA-DEX'

args = commandArgs(trailingOnly=TRUE)

if (length(args)>0){
      dataset <- args[1]
    }

bed <- read.table(gzfile(paste0(data_dir, dataset, ".bed.gz")), comment.char="", header=TRUE)
colnames(bed)[1:4] <- c("#Chr", "start", "end", "ID")
colnames(bed) <- gsub("[.]", "-", colnames(bed))

# get the count matrix:
data <- bed[,5:ncol(bed)]
rownames(data) <- bed[,4]

#PCA
library(irlba)
PCs <- prcomp_irlba(t(data), n=20)
summary(PCs)

mypcs <- as.data.frame(t(PCs$x))
colnames(mypcs) <- colnames(data)

#save the residuals:
mypcs <- rbind(as.character(colnames(mypcs)), mypcs)
rownames(mypcs)[1] <- "id"
mypcs <- cbind(as.character(rownames(mypcs)), mypcs)
write.table(mypcs, paste0("GE_PCs/",dataset,"_20GE_PCs.txt"), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

# save the covariate files subsetted to 1-20 GEPCs to use as covariate by FastQTL:
for(i in 1:20){
write.table(mypcs[1:i,], file=paste0("./GE_PCs/", dataset, "_",i,"GEPCs.txt"), col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
}

### END 1/7/2021
