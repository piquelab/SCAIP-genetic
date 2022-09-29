# Thisre script extracts dosages from the vcf file for all the SNPs tested by FastQTL for each condition
# based on /wsu/home/groups/piquelab/SCAIP/LDA/bin_model/genotypes/get_dosages.R
# 1/20/2021 JR
# last modified 11/08/2021, JW

library(data.table)
library(tidyverse)

rm(list=ls())

# passing arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)>0){
  cell <- args[1]
  lda <- args[2]
  treat <- args[3]
}else{
  cell <- "Bcell"
  lda <- "LPS"
  treat <- "LPS-EtOH"
} 

##################################################
### 1. concatenate the chunked FastQTL output: ###
##################################################
# do only once:
option <- "DiagLDA2"
outdir <- paste(option, "/2_fastQTL.test/", sep="")
if ( !file.exists(outdir) ) dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

fn <- paste0(outdir, cell, "_lda", lda, "_trt", treat, ".nominals.eQTL.txt.gz")

if( !file.exists(fn)){
##     system(paste0("for j in $(seq 1 30); do
##      cat tests/", lda, "/",dataset,".nominals.chunk$j.txt
## done | gzip -c > tests/", lda, "/",dataset,".nominals.eQTL.txt.gz;
## done"
   system(paste0("zcat ", outdir, cell, "_lda", lda, "_trt", treat, ".nominals.chunk*.txt.gz",
          " | gzip -c > ", fn))

# remove the chunks:
system(paste0("rm ", outdir,  cell, "_lda", lda, "_trt", treat, ".nominals.chunk*.txt.gz" ))

}


#######################################
### 2. eQTL results and coordinates ###
#######################################

if(file.exists(fn)){
    
   res <- fread(fn,  data.table=F)
   colnames(res) <- c("ENSG", "varID", "position", "pvalue", "slope_sc")
   coords <- separate(res, 2, into=c("chr", "position", NA, NA), sep=":", remove=FALSE)[,c("chr", "position")]
# sort:
   coords$chr <- as.numeric(coords$chr)
   coords$position <- as.numeric(coords$position)
   corrds <-  setorder(coords, chr, position)

# save the SNP-gene pairs to be tested:
   outdir <- paste(option, "/3_genotypes/eQTL_coordinates/", sep="")
   if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=TRUE)

   opfn1 <- paste0(outdir, "1_", cell, "_lda", lda, "_trt", treat, "_pairs.txt")  
   write.table(res, opfn1, sep="\t", quote=FALSE, row.names=FALSE)
   ## 
   opfn2 <- paste0(outdir, "2_eQTL_coordinates_", cell, "_lda", lda, "_trt", treat, ".txt") 
   write.table(coords, opfn2, sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
    
   opfn3 <- paste0(outdir, "3_eQTL_coordinates_uniq_", cell, "_lda", lda, "_trt", treat, ".txt") 
   system(paste0("uniq ", opfn2, ">", opfn3))
}

    
####################################
### 3. extract dosages from vcf: ###
####################################
    
outdir <- paste(option, "/3_genotypes/dosages/", sep="")
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=TRUE)

#system("module load bcftools")
opfn3 <- paste0(option, "/3_genotypes/eQTL_coordinates/3_eQTL_coordinates_uniq_", cell, "_lda", lda, "_trt", treat, ".txt")     
opfn4 <- paste0(option, "/3_genotypes/dosages/1_eQTL_signif_dosages_", cell, "_lda", lda, "_trt", treat, ".txt")
system(paste0("bcftools view -R ", opfn3, " SCAIP1-6_filtered_mock_reheaded.vcf.gz | bcftools query -f '%CHROM\t%POS\t%ID[\t%DS]\n' > ", opfn4))          
### add a header:
dosages <- fread(opfn4, sep="\t", stringsAsFactors=F, data.table=F)
# read in the file with dbgap IDs:
IDs <- read.table("SCAIP1-6_sample_mock.txt", stringsAsFactors=F)[,1]
colnames(dosages) <- c("chr", "pos", "varID", IDs)
# remove repeated dosages (where do they come from, though??):
dosages <- dosages[!duplicated(dosages$varID),]
# save with colnames:
write.table(dosages, opfn4, sep="\t", col.names=T, row.names=F, quote=F)
system(paste0("bgzip ", opfn4))
### END 1/20/2021
