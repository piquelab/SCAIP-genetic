# this script processes the uniqed eQTL coordinates for each condition to get unique eQTL coordinates overall and then dosages
# based on /wsu/home/groups/piquelab/SCAIP/eQTL/FastQTL/nominals/analysis_SCAIP1-6/get-all-uniq-eQTLs.R
# 1/19/2021 JR

library(data.table)
library(qqman)
library(qvalue)

datasets <- read.table("../datasets.txt",sep="\t",stringsAsFactors=F)[,1]

# read in each of the datasets and merge:
ll <- data.frame()

# run through all cell types and contrasts:
for(i in datasets){
  # load the data:
  tmp <- fread(paste0("eQTL_coordinates/eQTL_coordinates_uniq_", i, ".txt"))
  ll <- rbind(ll,tmp)
      }

# keep unique:
# add unique idenifier:
ll$uniq <- paste0(ll$V1,ll$V2)
df <- ll[!duplicated(ll$uniq),]
df <- df[order(as.numeric(df$V1),as.numeric(df$V2)),]

# extract the dosages from the vcf file:
write.table(df[,1:2],"eQTL_coordinates/all_eQTL_coordinates_uniq.txt",sep="\t",col.names=F,row.names=F, quote=F)
system("bcftools query -R eQTL_coordinates/all_eQTL_coordinates_uniq.txt -f '%CHROM\\t%POS\\t%ID\\t%AF[\\t%DS]\\n' ../../vcf/SCAIP1-6_filtered_AF.vcf.gz > eQTL_coordinates/all_uniq_eQTL_dosages.txt")

# add colnames:
dos <- fread("eQTL_coordinates/all_uniq_eQTL_dosages.txt")
dbgaps <- read.table("../SCAIP1-6_dbgaps.txt",stringsAsFactors=F)[,1]
colnames(dos) <- c("chr","pos","varID","AF",dbgaps)
write.table(dos,"eQTL_coordinates/all_uniq_eQTL_dosages.txt",sep="\t",col.names=T, row.names=F, quote=F)

# now run through all SNP-gene pairs:
# read in each of the datasets and merge:
ff <- data.frame()

# run through all cell types and contrasts:
for(i in datasets){
  # load the data:
 temp <- fread(paste0("eQTL_coordinates/signif_", i, ".txt"))[,1:2]
 ff <- rbind(ff,temp)
      }

# add gene-SNP column to keep only unique:
ff$geneSNP <- paste0(ff$ENSG,"_",ff$varID)
# uniq:
ff <- ff[!duplicated(ff$geneSNP),]
# save:
write.table(ff,"eQTL_coordinates/all_signif_SNP-gene_pairs.txt",sep="\t",col.names=T, row.names=F, quote=F)

sessionInfo()

### END 11/4/2021
