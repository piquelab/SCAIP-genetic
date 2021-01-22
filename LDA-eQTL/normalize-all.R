# this script takes in expression bed files summed across individuals in each LDA bin, quantile-voom-normalizes it, extracts residuals and saves it in bed format
# based on /wsu/home/groups/piquelab/SCAIP/LDA/bin_model/normalized_data/normalize-all.R
# 1/20/2021 JR

library(edgeR)
library(dplyr)
library(data.table)

args = commandArgs(trailingOnly=TRUE)

cell <- 'Bcell'
lda= "LDA2"

if (length(args)>0){
    cell <- args[1]
    lda <- args[2]
  }

treats <- c("_CTRL","_LPS-DEX","_LPS-EtOH","_PHA-DEX","_PHA-EtOH")

# threshold of minimum cells/bin and minimum bins with ncellsbin:
ncellsbin=5
minbin=2

# load GE data:
load(paste0("/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/9_RNA.dynamic2_output/Filter2_DEG6571/Old/",lda,"Bin/YtX.",cell,".sum.RData"))
# load the files with cells/individual:
load(paste0("/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/9_RNA.dynamic2_output/Filter2_DEG6571/Old/",lda,"Bin/0_ncell.", cell,".ave.RData"))
ncell <- data.frame(ncell)
# drop the samples with fewer than ncells per bin and :
ncell$dbgap_cell_trt <- gsub("_[1-9]","",ncell[,1])
ncell$bin <- gsub(".*_","",ncell$x)
bin_counts <- data.frame(ncell %>% group_by(dbgap_cell_trt) %>% summarize(nbins=n_distinct(bin)))
## cell_counts <- data.frame(ncell %>% group_by(dbgap_cell_trt) %>% summarize(cells=sum(ncell)))
ncell <- left_join(ncell, bin_counts)
keep <- ncell[ncell$ncell>=ncellsbin & ncell$nbins>=minbin,"x"]
YtX <- YtX[,colnames(YtX) %in% keep]
    
# add coordinates to make the bed files:
# load the annotation file:
anno <- read.table("/wsu/home/groups/piquelab/data/gencode/Gencode_human/release_31/gencode.v31.annotation.gff3.gz",header=F,stringsAsFactors=F)
anno$gene_id <- gsub("ID=","",anno$V9)
anno$gene_id <- gsub(";.*","",anno$gene_id)
# remove the _PAR_Y versions:
anno <- anno[-grep("PAR_Y",anno$gene_id),]
anno <- filter(anno,V3=="gene")
rownames(anno) <- anno$gene_id

# cycle through all the trts:
for(trt in treats){
# keep only the columns for the current trt:
counts <- YtX[,grep(trt,colnames(YtX))]

# remove everything but dbgap and bin number from colnames:
colnames(counts) <- gsub(paste0("_",cell,trt,"_SCAIP[1-6]"), "", colnames(counts))
    
# add coordinates::
counts <- counts[rownames(counts) %in% rownames(anno),]
bed <- anno[rownames(counts),c(1,4,5,10,7)]
colnames(bed) <- c("Chr", "min", "max", "ID", "strand")
bed[,1] <- gsub("chr", "", bed[,1])
# make the TSS the start:
bed[bed$strand=="+","start"] <- bed[bed$strand=="+","min"]
bed[bed$strand=="-","start"] <- bed[bed$strand=="-","max"]
bed[bed$strand=="-","end"] <- bed[bed$strand=="-","min"]
bed[bed$strand=="+","end"] <- bed[bed$strand=="+","max"]
# keep only relevant comumns:
bed <- bed[,c("Chr","start","end","ID")]

counts <- cbind(bed,counts)

# remove sex and mt chr (done already but make sure):
counts <- counts[counts[,1] %in% c(1:22),]

# get the count matrix:
data <- counts[,5:ncol(counts)]
data <- mutate_all(data, function(x) as.numeric(as.character(x)))
data <- as.matrix(data)
rownames(data) <- counts[,4]
class(data)# <- numeric

# load covariate file:
cv <- read.table("../../covariates/SCAIP1-6_ALOFT_cv.txt", sep="\t", header=T, comment="", quote='"', stringsAsFactors=F)
# important! order rows in cv and colnames in data the same!
rownames(cv) <- cv$dbgap.ID
# and triplicate the rows (for the three bins):
cv <- cv[gsub("_.*","",colnames(data)),]

## Normalization of data
# make edgeR object
dge <- DGEList(data)
sum(gsub("_.*","",colnames(dge$counts))==cv$dbgap.ID)
#Transform counts to counts per million
cpm <- cpm(dge)
samples <- dim(data)[2]
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
design <- model.matrix(~0+ cv$Batch + cv$Sex + cv$cage1 + cv$SCAIP1_6_genPC1 + cv$SCAIP1_6_genPC2 + cv$SCAIP1_6_genPC3)
# it usually doesn't matter which model you use - same result:
v_e <- voom(dge, design, plot=FALSE,normalize.method='quantile')
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
colnames(Res) <- gsub("[.]","-",colnames(genes_normed_baseline))
bed <- counts[as.character(rownames(Res)),1:4]
Resbed <- cbind(bed, Res)

# sort the bed file:
Resbed <- Resbed[order(as.numeric(Resbed[,1]), Resbed$start),]

# rename first line according to bed standard:
colnames(Resbed)[1] <- "#Chr"

# save:
write.table(Resbed, paste0("normalized_data/",lda,"/", cell, trt, ".bed"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
# bgzip and index
system(paste0("bgzip normalized_data/",lda,"/", cell, trt, ".bed"))
system(paste0("tabix -p bed normalized_data/",lda,"/", cell, trt, ".bed.gz"))
    
    }

# save file with gene-regions to be tested pairs:
## distance <- 50000
# not worth it; easier to run FastQTL and get output


### END 1/20/2021
