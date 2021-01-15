# this script transforms mean-corrected dispersion data and extracts residuals for QTL mapping
# based on ../../varianceQTL/get_residuals.R but with different input and filters
# 1/11/2021 JR

library(preprocessCore)
library(dplyr)

# set min number of individuals per batch:
indivBatchFilter <- 3

# load the dispersion data from Julong:
load("/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/10_RNA.Variance_output/tmp9/1.2_Sel.PhxNew.RData")

# load the annotation to add chromosomal coordinates:
anno <- read.table("/wsu/home/groups/piquelab/data/gencode/Gencode_human/release_31/gencode.v31.annotation.gff3.gz",header=F,stringsAsFactors=F)
anno$gene_id <- gsub("ID=","",anno$V9)
anno$gene_id <- gsub(";.*","",anno$gene_id)
# remove the _PAR_Y versions:
anno <- anno[-grep("PAR_Y",anno$gene_id),]
## anno$g.id <- gsub("[.].*","",anno$gene_id) # NOT for SCAIP1-6!
anno <- filter(anno,V3=="gene")
rownames(anno) <- anno$gene_id
# remove non-autosomal chromosomes:
annoau <- anno[gsub("chr","",anno[,1]) %in% c(1:22),]
## # subset to protein-coding only:
## annoaupr <- annoau[grep("gene_type=protein_coding",annoau$V9),] # done by Julong already

# add gene coordinates:
# remove genes not present in annotaion (all are present, though)
Phx <- PhxNew2[rownames(PhxNew2) %in% rownames(annoau),]
# sort by position:
Phx <- Phx[rownames(annoau)[rownames(annoau) %in% rownames(Phx)],]
annobed <- annoau[,c(1,4,5,7,10)]
colnames(annobed) <- c("Chr","min","max","strand","ID")
annobed$Chr <- gsub("chr","",annobed$Chr)
# make the TSS the start:
annobed[annobed$strand=="+","start"] <- annobed[annobed$strand=="+","min"]
annobed[annobed$strand=="-","start"] <- annobed[annobed$strand=="-","max"]
annobed[annobed$strand=="-","end"] <- annobed[annobed$strand=="-","min"]
annobed[annobed$strand=="+","end"] <- annobed[annobed$strand=="+","max"]
# keep only relevant comumns:
annobed <- annobed[,c("Chr","start","end","ID")]
colnames(annobed)[1] <-"#Chr"


# load covariate file:
cvall <- read.table("../../covariates/SCAIP1-6_ALOFT_cv.txt", sep="\t", header=T, comment="", quote='"', stringsAsFactors=F)
# important! order rows in cv and colnames in data the same!
rownames(cvall) <- cvall$dbgap.ID

# go through the cell types:
clusters <- unique(gsub("_.*","",colnames(Phx)))
# and treatments:
trts <- c("_CTRL", "_PHA-DEX","_PHA-EtOH","_LPS-DEX","_LPS-EtOH")

#loop through all the trts:
for(t in trts){
        treat <- Phx[,grep(t,colnames(Phx))]
            # remove the "_CTRL":
            colnames(treat) <- gsub(t,"",colnames(treat))
            for(cn in clusters){
                   cellorig <- treat[,grep(cn,colnames(treat))]
                   #remove cell name from colnames:
                   colnames(cellorig) <- gsub(paste0(cn,"_"),"",colnames(cellorig))
# split by batch and then NA rows which have dispersion measure in <3 individuals:
cellsub <- cellorig[,-1:-ncol(cellorig)]
batches <- unique(gsub(".*_","",colnames(cellorig)))
for (b in batches){
# take out the columns for that batch:
current <- cellorig[,grep(b,colnames(cellorig))]
# set to NA if fewer than indivBatchFilter non-NAs:
current[rowSums(!is.na(current))<indivBatchFilter,] <- NA
#merge back:
cellsub <- cbind(cellsub,current)
}
# drop batch info from colnames:
colnames(cellsub) <- gsub("_SCAIP.*","",colnames(cellsub))
# order columns alphabetically as they are in vcf file:
                   cellsub <- cellsub[,sort(colnames(cellsub))]
                   # drop the genes with >80% missing data:
                   cell <- cellsub[rowSums(is.na(cellsub))<0.8*ncol(cellsub),]                   
                   # save qqnormed data:
                   qnorm <- normalize.quantiles(cell)
                   colnames(qnorm) <- colnames(cell)
                   rownames(qnorm) <- rownames(cell)
                   # add coordinates:
                   qnormbed <- cbind(annobed[rownames(qnorm),],qnorm)
                   # extract residuals
                   cv <- cvall[colnames(cell),]
                   # baseline model:
                   design <- model.matrix(~0+ cv$Batch + cv$Sex + cv$cage1 + cv$SCAIP1_6_genPC1 + cv$SCAIP1_6_genPC2 + cv$SCAIP1_6_genPC3)
                   #get residuals:
                   X <- design
                   H <- X %*% solve(t(X) %*% X) %*% t(X)
                   dim(H)
                   He <- (diag(rep(1,ncol(H)))-H)
                   Resqnorm <- as.matrix(qnorm) %*% He
                   sum(abs(t(He)-He)) #Should be almost 0
                   Resqnorm <- data.frame(Resqnorm)
                   colnames(Resqnorm) <- gsub("[.]","-",colnames(qnorm))
                   # calculate PCs of residuals and save as covariates:
                   # drop NAs:
                   Resqnormclean <- Resqnorm[rowSums(is.na(Resqnorm))==0,]
                   #PCA
PCs <- prcomp_irlba(t(Resqnormclean), n=20)
summary(PCs)
mypcs <- as.data.frame(PCs$x)
mypcs$dbgap.ID <- colnames(Resqnormclean)
covs <- left_join(cv[,c("dbgap.ID","Batch","Sex","cage1","SCAIP1_6_genPC1","SCAIP1_6_genPC2","SCAIP1_6_genPC3")],mypcs)
tcovs <- data.frame(t(covs))
# generate covariates for FastQTL:
for(i in 0:20){
write.table(tcovs[1:(7+i),],paste0("covariates/",cn,t,"_covs_",i,"PCs.txt"), sep="\t", row.names=TRUE, col.names=FALSE, quote=FALSE)
    }
                   ## # sort in order: # do not do this in R! doesn't work as it should for some reason
                   ## qnormbed[,1] <- as.numeric(qnormbed[,1])
                   ## qnormbed[,2] <- as.numeric(qnormbed[,2])
                   ## qnormbed <-qnormbed[order(qnormbed[,1],qnormbed[,2]),]
                   write.table(qnormbed,paste0("./normalized_dispersion_residuals/",cn, t,"_dispersion_unsorted.bed"), sep="\t", col.names=T, quote=F,row.names=F)
                   # sort the bed file:
                   system(paste0("(head -n 1 ./qnormed_dispersion/",cn, t,"_dispersion_unsorted.bed && tail -n +2 ./qnormed_dispersion/",cn, t,"_dispersion_unsorted.bed | sort -k 1,1 -k2,2n)> ./qnormed_dispersion/",cn, t,"_dispersion.bed"))
                   system(paste0("bgzip -f qnormed_dispersion/",cn, t,"_dispersion.bed"))
                   system(paste0("tabix -p bed qnormed_dispersion/",cn, t,"_dispersion.bed.gz"))
                        }
        }

### END 1/15/2021 JR
