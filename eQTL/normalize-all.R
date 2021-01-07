># this script takes in an R object matrix with summed GE counts for each individual-condition, voom-normalizes it, extracts residuals and saves it in bed format
# based on /wsu/home/groups/piquelab/SCAIP/eQTL/FastQTL/GE_files/normalize-all_JWcounts.R
# JR 1/7/2021

library(edgeR)
library(dplyr)

treats <- c("_CTRL","_LPS-DEX","_LPS-EtOH","_PHA-DEX","_PHA-EtOH")

# load GE data:
load("/nfs/rprdata/julong/SCAIP/analyses/SCAIP-ALL-2020.03.23/6_DEG.CelltypeNew_output/Filter2/YtX_sel.comb.RData")
YtX <- YtX_sel

# remove the batch from colnames:
colnames(YtX) <- gsub("_SCAIP[1-6]", "", colnames(YtX))

# grab the unique cell types:
cells <- unique(gsub("_.*","",colnames(YtX)))

# add coordinates to make the bed files:
# load the annotation file gencode gchr37 v31:
anno <- read.table("/wsu/home/groups/piquelab/data/gencode/Gencode_human/release_31/gencode.v31.annotation.gff3.gz",header=F,stringsAsFactors=F)
anno <- filter(anno,V3=="gene")
anno$gene_id <- gsub("ID=","",anno$V9)
anno$gene_id <- gsub(";.*","",anno$gene_id)
# remove the _PAR_Y versions:
anno <- anno[-grep("PAR_Y",anno$gene_id),]
anno$g.id <- gsub("[.].*","",anno$gene_id)
rownames(anno) <- anno$gene_id

# add annotation info to GE data:
# remove the genes missing from the annotation:
YtX <- YtX[rownames(YtX) %in% rownames(anno),] # 0 missing
bed <- anno[rownames(YtX),c(1,4,5,10,7)]
colnames(bed) <- c("Chr", "min", "max", "ID", "strand")
bed[,1] <- gsub("chr", "", bed[,1])
# make the TSS the start:
bed[bed$strand=="+","start"] <- bed[bed$strand=="+","min"]
bed[bed$strand=="-","start"] <- bed[bed$strand=="-","max"]
bed[bed$strand=="-","end"] <- bed[bed$strand=="-","min"]
bed[bed$strand=="+","end"] <- bed[bed$strand=="+","max"]
# keep only relevant columns:
bed <- bed[,c("Chr","start","end","ID")]
stopifnot(identical(rownames(YtX),rownames(bed)))
YtX <- cbind(bed,YtX)

# remove sex and mt chr (done already for this data, but make sure)
YtX <- YtX[YtX[,1] %in% c(1:22),]

# cycle through all the trts:
for(trt in treats){
# keep only the columns for the current trt:
countstrt <- YtX[,c(1:4,grep(trt,colnames(YtX)))]
# and all the cells:
for(cell in cells){
counts <- countstrt[,c(1:4,grep(cell,colnames(countstrt)))]
# remove everything but dbgap and bin number from colnames:
colnames(counts) <- gsub(paste0(cell,trt,"_"), "", colnames(counts))

# get the count matrix:
data <- counts[,5:ncol(counts)]
rownames(data) <- counts[,4]

# load covariate file:
cv <- read.table("../../covariates/SCAIP1-6_ALOFT_cv.txt", sep="\t", header=T, comment="", quote='"', stringsAsFactors=F)
# important! order rows in cv and colnames in data the same!
rownames(cv) <- cv$dbgap.ID
cv <- cv[colnames(data),]

## Normalization of data
# make edgeR object
dge <- DGEList(data)
sum(colnames(dge$counts)==cv$dbgap.ID)
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
v_e <- voom(dge, design, plot=FALSE, normalize.method="quantile")
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
bed <- counts[as.character(counts[,4]) %in% as.character(rownames(Res)),1:4]
Resbed <- cbind(bed, Res)

# sort the bed file:
Resbed <- Resbed[order(Resbed[,1], Resbed$start),]

# add the #Chr back:
colnames(Resbed)[1] <- "#Chr"

# save
write.table(Resbed, paste0("normalized_GE_residuals/",cell,trt, ".bed"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
# bgzip and index
system(paste0("bgzip normalized_GE_residuals/", cell, trt, ".bed"))
system(paste0("tabix -p bed normalized_GE_residuals/", cell, trt, ".bed.gz"))

}
    }

sessionInfo()
### END 1/7/2021
