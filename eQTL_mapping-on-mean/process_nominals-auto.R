# this script processes FastQTL output to summarize numbers of mean GE eQTLs in SCAIP1-6
# based on ../meanersionQTL/process_nominals-auto.R
# 1/13/2021

library(data.table)
library(qqman)
library(qvalue)

dataset <- "Tcell_LPS-EtOH"
args = commandArgs(trailingOnly=TRUE)

# get the PC number:
if (length(args)>0){
      ## pc <- args[1]
      dataset <- args[1]
      }

pc <- 0

# set FDR threshold;
FDR <- 0.1

# concatenate the chunked FastQTL output:
# do only once:
if(!file.exists(paste0("mean-eQTL_output/",dataset,".GEPC0.nominals.eQTL.txt.gz"))){
system(paste0("for j in $(seq 1 30); do
     cat mean-eQTL_output/",dataset,".GEPC0.nominals.chunk$j.txt
done | gzip -c > mean-eQTL_output/",dataset,".GEPC0.nominals.eQTL.txt.gz;
done"
))
# remove the chunks:
system(paste0("rm mean-eQTL_output/",dataset,".GEPC0.nominals.chunk*"))
}


mean.qq <- fread(paste0("zcat ./mean-eQTL_output/", dataset, ".GEPC", pc, ".nominals.eQTL.txt.gz"), sep=" ",head=F,col.names=c("pid","sid","distace","pvalue","estimate"))
mean.qq$qvalue <- qvalue(mean.qq$pvalue)$qvalues

tab <- t(c(dataset,sum(mean.qq$qvalue<FDR), length(unique(mean.qq$pid[mean.qq$qvalue<FDR]))))

# write the number of eQTLs, egenes:
write.table(tab, file=paste0("mean_nominals_eqtls-egenes_summary.txt"), sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE, append=TRUE)

# save the mean genes:
if(length(unique(mean.qq$pid[mean.qq$qvalue<FDR]))>0){
write.table(unique(mean.qq$pid[mean.qq$qvalue<FDR]),paste0("mean_egenes/",dataset,"_genes.txt"),quote=F, row.names=F, col.names=F)
}

# make QQ plots:
png(paste0("./plots/QQ_", dataset, "_mean.png"))
qq(mean.qq$pvalue, main=paste0("QQ plot of mean FastQTL analysis for ", dataset ))
dev.off()


### END 1/13/2021
