# this script compares number of significant results from FastQTL and mashr on FastQTL results
# 2/1/2021 JR

library(data.table)
library(dplyr)

thresh <- 0.1

# load the mashr results:
mlfsr <- data.frame(fread("output/lfsr/lfsr_mashr_all.txt.gz"))
# add colnames:
cols <- read.table("FastQTL_all_conditions-colnames.txt",stringsAsFactors=F)[,1]
colnames(mlfsr) <- cols
# add rownames:
rows <- read.table("rownames_mash_updated.txt",stringsAsFactors=F)[,1]
## # remove the removed rows:
## rows <- rows[-grep("ENSG00000166710",rows)]
## rows <- rows[-grep("ENSG00000167996",rows)]
rownames(mlfsr) <- rows
## # save the updated rownames:
## write.table(rows,"rownames_mash_updated.txt",row.names=F,col.names=F,quote=F)

## sum(mlfsr<0.1)
## sum(mlfsr<0.1)/nrow(mlfsr)/ncol(mlfsr)
colSums(mlfsr<0.1)
## unique(rowSums(mlfsr<0.1))
## summary(as.factor(rowSums(mlfsr<0.1)))

# summarize the data per-gene:
# change all significant values to 1, and insignificant to 0:
mlfsr[mlfsr<thresh] <- -1
mlfsr[mlfsr>thresh] <- 0
mlfsr[mlfsr==-1] <- 1

mlfsr$gene <- gsub("_.*","",rownames(mlfsr))
mlfsr_gene <- data.frame(mlfsr %>% group_by(gene) %>% summarize_each(funs(sum)))
# drop gene column:
rownames(mlfsr_gene) <- mlfsr_gene$gene
mlfsr_gene <- mlfsr_gene[,-1]
# fix column names:
colnames(mlfsr_gene) <- gsub("[.]","+",colnames(mlfsr_gene))
colSums(mlfsr_gene>0)
sum(rowSums(mlfsr_gene>0)>0)

cells <- unique(gsub("_.*","",colnames(mlfsr)))
trts <- unique(gsub(".*_","",colnames(mlfsr)))
l=list()
for(cell in cells[1:4]){
    df <- mlfsr[,grep(cell,colnames(mlfsr))]
    print(cell)
       sdf <- data.frame(summary(as.factor(rowSums(df<0.1))))
    sdf$count <- rownames(sdf)
    l[[cell]] <- sdf
}
for(trt in trts[1:5]){
    df <- mlfsr[,grep(trt,colnames(mlfsr))]
    print(trt)
   sdf <- data.frame(summary(as.factor(rowSums(df<0.1))))
    sdf$count <- rownames(sdf)
    l[[trt]] <- sdf
    }


gl=list()
for(cell in cells[1:4]){
    df <- mlfsr_gene[,grep(cell,colnames(mlfsr_gene))]
    print(cell)
       sdf <- data.frame(summary(as.factor(rowSums(df==1))))
    sdf$count <- rownames(sdf)
    gl[[cell]] <- sdf
}
for(trt in trts[1:5]){
    df <- mlfsr_gene[,grep(trt,colnames(mlfsr_gene))]
    print(trt)
   sdf <- data.frame(summary(as.factor(rowSums(df==1))))
    sdf$count <- rownames(sdf)
    gl[[trt]] <- sdf
    }


# get number of unique to each column:
mlfsr_o <- mlfsr_gene[,1:20]

for(i in 1:ncol(mlfsr_o)){
    print(colnames(mlfsr_o)[i])
print(    sum(mlfsr_o[,i]==1 & rowSums(mlfsr_o[,-i])==0 ))
    }


# load the ashr results:
load(paste0("./mashr_input.Rd"))
lfsr <- lfsrs

# change all significant values to 1, and insignificant to 0:
lfsr[lfsr<thresh] <- -1
lfsr[lfsr>thresh] <- 0
lfsr[lfsr==-1] <- 1

# drop all 0 rows for faster processing:
lfsr <- lfsr[rowSums(lfsr,na.rm=T)>0,]
colSums(lfsr==1,na.rm=T)

# summarize the data per-gene:
lfsr$gene <- gsub("_.*","",rownames(lfsr))
lfsr_gene <- data.frame(lfsr %>% group_by(gene) %>% summarize_each(funs(sum)))
# drop gene column:
rownames(lfsr_gene) <- lfsr_gene$gene
lfsr_gene <- lfsr_gene[,-1]
# fix column names:
colnames(lfsr_gene) <- gsub("[.]","+",colnames(lfsr_gene))
colSums(lfsr_gene>0, na.rm = T)

sum(rowSums(lfsr_gene, na.rm = T)>0)
length(unique(gsub("_.*","",rownames(mlfsr[rowSums(mlfsr[,1:20], na.rm = T)>0,]))))

## # plot matrix of correlations:
## # set up colors:
## ann_colors = list(cell=c("Tcell"="#e41a1c", "NKcell"="#377eb8", "Bcell"="#4daf4a", "Monocyte"="#984ea3"),trt=c("CTRL"="#828282", "LPS"="#fb9a99", "LPS.EtOH"="#fb9a99", "LPS-EtOH"="#fb9a99", "LPS.DEX"="#e31a1c","LPS-DEX"="#e31a1c","LPS+DEX"="#e31a1c","PHA"="#a6cee3", "PHA.EtOH"="#a6cee3", "PHA-EtOH"="#a6cee3", "PHA.DEX"="#1f78b4", "PHA-DEX"="#1f78b4", "PHA+DEX"="#1f78b4"))
## annotation_col <- data.frame(matrix(nrow=20,ncol=2))
## rownames(annotation_col) <- colnames(lfsr)
## colnames(annotation_col) <- c("cell","trt")
## annotation_col$trt <- gsub(".*_","",rownames(annotation_col))
## annotation_col$cell <- gsub("_.*","",rownames(annotation_col))
## annotation_col$cell <- gsub("*_.*","",annotation_col$cell)
## library(pheatmap)
## pdf("./plots/mashr_pheatmap_corr_matrix_all.pdf",width=11, height=10)
## # calculate correlations:
## corr <- cor(lfsr, use='pairwise.complete.obs')
## library(psych)
## corr.psych <- corr.test(lfsr[,], adjust="none")
## # blank out non-significant correlations:
## corr <- corr*(corr.psych$p<0.05)
## pheatmap(corr, annotation_col = annotation_col, annotation_colors = ann_colors,cellheight=21,cellwidth=21,cluster_rows =F,cluster_cols =F)
## dev.off()




### END 11/30/2020

sessionInfo()

