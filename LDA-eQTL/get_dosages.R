# Thisre script extracts dosages from the vcf file for all the SNPs tested by FastQTL for each condition
# based on /wsu/home/groups/piquelab/SCAIP/LDA/bin_model/genotypes/get_dosages.R
# 1/20/2021 JR

library(data.table)
library(tidyr)

args <- commandArgs(trailingOnly=TRUE)
lda <- "LDA1"
dataset <- "Bcell_LPS-DEX"

if (length(args)>0){
        dataset <- args[1]
        lda <- args[2]
         }


# concatenate the chunked FastQTL output:
# do only once:
if(!file.exists(paste0("tests/", lda, "/",dataset,".nominals.eQTL.txt.gz"))){
    system(paste0("for j in $(seq 1 30); do
     cat tests/", lda, "/",dataset,".nominals.chunk$j.txt
done | gzip -c > tests/", lda, "/",dataset,".nominals.eQTL.txt.gz;
done"
))
    # remove the chunks:
    system(paste0("rm tests/", lda, "/",dataset,".nominals.chunk*"))
    }

results <- fread(paste0("zcat tests/", lda, "/",dataset, ".nominals.eQTL.txt.gz"), sep=" ",data.table=F)

colnames(results) <- c("ENSG", "varID", "position", "pvalue", "slope_sc")

coords <- separate(results,2,into=c("chr", "position",NA,NA),sep=":", remove=FALSE)[,c("chr", "position")]
# sort:
coords$chr <- as.numeric(coords$chr)
coords$position <- as.numeric(coords$position)
corrds <-  setorder(coords,chr,position)

# save the SNP-gene pairs to be tested:
write.table(results, paste0("genotypes/eQTL_coordinates/",lda,"/", dataset,"_pairs.txt"),sep="\t", quote=FALSE, row.names=FALSE)
write.table(coords, paste0("genotypes/eQTL_coordinates/",lda,"/eQTL_coordinates_", dataset,".txt"),sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
system(paste0("uniq genotypes/eQTL_coordinates/",lda,"/eQTL_coordinates_", dataset,".txt  > genotypes/eQTL_coordinates/",lda,"/eQTL_coordinates_uniq_", dataset,".txt"))

# extract dosages from vcf:
system(paste0("bcftools view -R genotypes/eQTL_coordinates/",lda,"/eQTL_coordinates_uniq_", dataset,".txt genotypes/SCAIP1-6_filtered_mock_reheaded.vcf.gz | bcftools query -f '%CHROM\t%POS\t%ID[\t%DS]\n' > genotypes/dosages/",lda,"/eQTL_signif_dosages.", dataset,".txt"))

# add a header:
dosages <- read.table(paste0("genotypes/dosages/",lda,"/eQTL_signif_dosages.",dataset,".txt"), sep="\t", stringsAsFactors=F)
# read in the file with dbgap IDs:
IDs <- read.table(paste0("./genotypes/samples_SCAIP1-6_mock.txt"), stringsAsFactors=F)[,1]
colnames(dosages) <- c("chr","pos","varID",IDs)
# remove repeated dosages (where do they come from, though??):
dosages <- dosages[!duplicated(dosages$varID),]
# save with colnames:
write.table(dosages, paste0("genotypes/dosages/",lda,"/eQTL_signif_dosages.",dataset,".txt"), sep="\t", col.names=T,row.names=F,quote=F)
system(paste0("bgzip genotypes/dosages/",lda,"/eQTL_signif_dosages.",dataset,".txt"))
### END 1/20/2021
