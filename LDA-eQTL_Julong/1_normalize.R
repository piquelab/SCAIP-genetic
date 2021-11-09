# this script takes in expression bed files summed across individuals in each LDA bin, quantile-voom-normalizes it, extracts residuals and saves it in bed format
# based on /wsu/home/groups/piquelab/SCAIP/LDA/bin_model/normalized_data/normalize-all.R
# 1/20/2021 JR
# modified 5/28/2021 JW
# Last modified 11/8/2021, JW

library(edgeR)
library(tidyverse)
library(data.table)

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

# threshold of minimum cells/bin and minimum bins with ncellsbin:
ncellbin=5
minbin=2
option <- "DiagLDA2"

outdir <- paste(option, "/1_normalized.data/", sep="") 
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=FALSE, recursive=TRUE) 

cat(cell, lda, treat, "\n")

########################
### 1. load GE data: ###
########################

prefix <- paste("/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/9_RNA.dynamic2_output/Filter2_DEG6571/", option, "/DLDA_Bin/", sep="")
fn <- paste(prefix, "2_", cell, ".", lda, "_YtX.sum.rds", sep="")
YtX <- read_rds(fn)
fn <- paste(prefix, "2_", cell, ".", lda, "_ncell.rds", sep="")
ncell <- read_rds(fn)

# drop the samples with fewer than ncells per bin and :
ncell <- ncell%>%
  mutate(dbgap_cell_trt=gsub("_[123]", "", x), bin=gsub(".*_", "", x))
##
bin_count <- ncell%>%
   filter(ncell>=ncellbin)%>% 
   group_by(dbgap_cell_trt)%>%
   summarize(nbin=n_distinct(bin),.groups="drop")%>%as.data.frame()

## cell_counts <- data.frame(ncell %>% group_by(dbgap_cell_trt) %>% summarize(cells=sum(ncell)))
ncell <- ncell%>%left_join(bin_count, by="dbgap_cell_trt")
keep <- ncell%>%filter(ncell>=ncellbin, nbin>=minbin)%>%dplyr::pull(x)

YtX <- YtX[, keep]


######################################
### 2. Add coordinates into counts ###
######################################

# load the annotation file: gencode grch37 v31
anno <- read.table("/wsu/home/groups/piquelab/data/gencode/Gencode_human/release_31/GRCh37_mapping/gencode.v31lift37.annotation.gff3.gz", header=F, stringsAsFactors=F)

## grch38
## anno <- read.table("/wsu/home/groups/piquelab/data/gencode/Gencode_human/release_31/gencode.v31.annotation.gff3.gz", header=F, stringsAsFactors=F)

anno <- anno%>%mutate(V1=gsub("chr", "", V1),
   gene_id=gsub("ID=|;.*", "", V9))%>%
   filter(!grepl("PAR_Y", gene_id), V3=="gene") #grepl("protein_coding", anno$V9))

anno <- anno%>%dplyr::filter(V1%in%c(1:22))

rownames(anno) <- anno$gene_id

#
bed <- anno%>%
   dplyr::rename("Chr"="V1", "min"="V4", "max"="V5",
                 "ID"="gene_id", "strand"="V7")%>%
   mutate(start=ifelse(strand=="+", min, max),
          end=ifelse(strand=="+", max, min)) 
#
bed <- bed%>%dplyr::select(Chr, start, end, ID)

## anno$gene_id <- gsub("ID=","",anno$V9)
## anno$gene_id <- gsub(";.*","",anno$gene_id)
## # remove the _PAR_Y versions:
## anno <- anno[-grep("PAR_Y",anno$gene_id),]
## anno <- filter(anno,V3=="gene")
## rownames(anno) <- anno$gene_id

# cycle through all the trts:
# keep only the columns for the current trt:

counts <- YtX[, grepl(treat, colnames(YtX))]
colnames(counts) <- gsub("_.*_", "_", colnames(counts))
counts <- as.data.frame(counts)
counts$ID <- rownames(counts)   
# add coordinates::
counts <- bed%>%inner_join(counts, by="ID")
rownames(counts) <- counts$ID

# get the count matrix:
data <- counts[,5:ncol(counts)]
data <- mutate_all(data, function(x) as.numeric(as.character(x)))
data <- as.matrix(data)
rownames(data) <- counts$ID


### load covariate file:
cv <- read.table("/wsu/home/groups/piquelab/SCAIP/covariates/SCAIP1-6_ALOFT_cv.txt", sep="\t", header=T, comment="", quote='"', stringsAsFactors=F)%>%
    dplyr::select(dbgap.ID, Batch, Sex, cage1,
                  SCAIP1_6_genPC1, SCAIP1_6_genPC2, SCAIP1_6_genPC3)
cv2 <- str_split(colnames(data), "_", simplify=TRUE)%>%as.data.frame()
names(cv2) <- c("dbgap.ID", "Bin")
cv2 <- cv2%>%left_join(cv, by="dbgap.ID")
rownames(cv2) <- colnames(data)

#########################
### 3. Normalize data ###
#########################

## Normalization of data
# make edgeR object
dge <- DGEList(data)
sum(gsub("_.*","",colnames(dge$counts))==cv2$dbgap.ID)
#Transform counts to counts per million
cpm <- cpm(dge)
samples <- ncol(data)
#Remove genes that are lowly expressed
table(rowSums(dge$counts==0)==samples) #Shows how many transcripts have 0 count across all samples
# use same thresholds as GTEx:
keep.exprs <- rowSums(cpm>=0.1)>=(0.2*samples)
#& rowSums(data>=6)>=(0.2*samples) #Only keep transcripts that have cpm>0.1 and count of 6 in at least 20% of the samples
dge <- dge[keep.exprs,, keep.lib.sizes=FALSE]
#Normalize data
dge <- calcNormFactors(dge, method = "TMM")
head(dge$samples$norm.factors)

# baseline model:
design <- model.matrix(~0+ cv2$Batch + cv2$Sex + cv2$cage1 + cv2$SCAIP1_6_genPC1 + cv2$SCAIP1_6_genPC2 + cv2$SCAIP1_6_genPC3)
# it usually doesn't matter which model you use - same result:
v_e <- voom(dge, design, plot=FALSE, normalize.method='quantile')
genes_normed_baseline <- data.frame(v_e$E)

## ## extract the residuals:
#get residuals:
X <- design
H <- X %*% solve(t(X) %*% X) %*% t(X)
dim(H)
He <- (diag(rep(1,ncol(H)))-H)
Res <- as.matrix(genes_normed_baseline) %*% He
sum(abs(t(He)-He)) #Should be almost 0

#save the residuals:
Res <- data.frame(Res)
names(Res) <- gsub("\\.", "-", names(genes_normed_baseline))
bed <- counts[as.character(rownames(Res)), 1:4]
Resbed <- cbind(bed, Res)

# sort the bed file:
Resbed <- Resbed[order(as.numeric(Resbed[,1]), Resbed$start),]

# rename first line according to bed standard:
names(Resbed)[1] <- "#Chr"
nn <- as.character(names(Resbed))
# save:
opfn <- paste0(outdir, cell, "_lda", lda, "_trt", treat, ".bed")
write.table(Resbed, file=opfn, sep="\t", row.names=FALSE, col.names=nn, quote=FALSE)
# bgzip and index
system(paste0("bgzip ", opfn))
system(paste0("tabix -p bed ", opfn, ".gz"))
    
# save file with gene-regions to be tested pairs:
## distance <- 50000
# not worth it; easier to run FastQTL and get output


### END 1/20/2021
