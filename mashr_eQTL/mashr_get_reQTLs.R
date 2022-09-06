# this script gets the pair-wise unique eQTLs ~reQTLs
# 2/3/2021 JR

library(data.table)
library(dplyr)

thresh <- 0.1

# load the mashr results:
mlfsr <- data.frame(fread("output/lfsr/lfsr_mashr_all.txt.gz"))
rownames(mlfsr) <- mlfsr[,1]
mlfsr <- mlfsr[,-1]
# add colnames:
cols <- read.table("FastQTL_all_conditions-colnames.txt",stringsAsFactors=F)[,1]
colnames(mlfsr) <- gsub("[.]","-",cols)
# add rownames:
## rows <- read.table("rownames_mash.txt",stringsAsFactors=F)[,1]
## fail1 <- read.table("failed_rows_344.txt",stringsAsFactors=F)[,1]
## fail2 <- read.table("failed_rows_353.txt",stringsAsFactors=F)[,1]
## rows <- rows[!rows %in% c(fail1, fail2)]
## write.table(rows, "rownames_mash_updated.txt", col.names=F, row.names=F, quote=F)
## rows <- read.table("rownames_mash_updated.txt",stringsAsFactors=F)[,1]
## rownames(mlfsr) <- rows
mean <- fread("output/posterior_mean/posterior_mean_mashr_all.txt.gz",data.table=F)
rownames(mean) <- mean[,1]
mean <- mean[,-1]
colnames(mean) <- gsub("[.]","-",cols)

# get pair-wise unique eQTLs:
pairs <- read.table("contrasts.txt",sep="\n",stringsAsFactors=F)[,1]

# get sharing and specificity:
uum=list()
top =list()
eqtl <- list()
df <- data.frame(row.names=gsub("\t",":",pairs),trt=gsub("\t.*","",pairs), ctrl =gsub(".*\t","",pairs), union_shared_sign=as.numeric(""), union_unshared_sign=as.numeric(""), trtsignif_unshared_sign=as.numeric(""), union_shared_magnitude=as.numeric(""), union_unshared_magnitude=as.numeric(""), trtsignif_unshared_magnitude=as.numeric(""), union_unshared_magnitude_0.1=as.numeric(""))
for(pair in pairs){
trt = gsub("\t.*","",pair)
ctrl =gsub(".*\t","",pair)
pair = gsub("\t",":",pair)
# get mashr significant eQTLs:
eqtl[ctrl] <- list(rownames(mlfsr)[mlfsr[,ctrl]<0.1])
eqtl[trt] <- list(rownames(mlfsr)[mlfsr[,trt]<0.1])
# by sign:
df[pair,"union_shared_sign"] <- sum((mlfsr[,trt]<0.1 | mlfsr[,ctrl]<0.1) & mean[,trt]*mean[,ctrl]>0)
df[pair,"union_unshared_sign"] <-sum((mlfsr[,trt]<0.1 | mlfsr[,ctrl]<0.1) & mean[,trt]*mean[,ctrl]<0)
df[pair,"trtsignif_unshared_sign"] <- sum((mlfsr[,trt]<0.1) & mean[,trt]*mean[,ctrl]<0)
# by magnitude:
df[pair,"union_shared_magnitude"] <- sum((mlfsr[,trt]<0.1 | mlfsr[,ctrl]<0.1) & (mean[,trt]/mean[,ctrl]<=2 & mean[,trt]/mean[,ctrl]>=0.5))
df[pair,"union_unshared_magnitude"] <- sum((mlfsr[,trt]<0.1 | mlfsr[,ctrl]<0.1) & (mean[,trt]/mean[,ctrl]>2 | mean[,trt]/mean[,ctrl]<0.5))
df[pair,"trtsignif_unshared_magnitude"] <- sum((mlfsr[,trt]<0.1) & (mean[,trt]/mean[,ctrl]>2 | mean[,trt]/mean[,ctrl]<0.5))
uum[pair]= list(rownames(mlfsr)[(mlfsr[,trt]<0.1 | mlfsr[,ctrl]<0.1) & (mean[,trt]/mean[,ctrl]>2 | mean[,trt]/mean[,ctrl]<0.5)])
# introduce threshold on min effect size of 0.1 in either condition 3/9/2021
df[pair,"union_unshared_magnitude_0.1"] <- sum((mlfsr[,trt]<0.1 | mlfsr[,ctrl]<0.1) & (mean[,trt]/mean[,ctrl]>2 | mean[,trt]/mean[,ctrl]<0.5) & (abs(mean[,trt]) >0.1 | abs(mean[,ctrl] >0.1)))
top[pair]= list(rownames(mlfsr)[(mlfsr[,trt]<0.1 | mlfsr[,ctrl]<0.1) & (mean[,trt]/mean[,ctrl]>2 | mean[,trt]/mean[,ctrl]<0.5) & (mean[,trt]>0.5 | mean[,ctrl]>0.5)])
}

# save results:
save(top, file="mashr-reQTLs_union-unshared-magnitude2-mlfsr0.1-mean0.5.Rd") 
save(eqtl, file="mashr-eQTLs-mlfsr0.1.Rd")
write.table(df, "mashr_sharing_eQTL.txt", col.names=T, row.names=F, quote=F, sep="\t")
save(uum, file="mashr-reQTLs_union-unshared-magnitude2-mlfsr0.1.Rd")

## sum(mlfsr<0.1)
## sum(mlfsr<0.1)/nrow(mlfsr)/ncol(mlfsr)
colSums(mlfsr<0.1)
## unique(rowSums(mlfsr<0.1))
## summary(as.factor(rowSums(mlfsr<0.1)))
egenes <- lapply(eqtl, function(x) unique(gsub("[.].*","",x)))

# summarize the data per-gene:


# get sharing and specificity:
l5 <- character()
l6 <- character()
df <- data.frame(row.names=gsub("\t",":",pairs),trt=gsub("\t.*","",pairs), ctrl =gsub(".*\t","",pairs), union_shared_sign=as.numeric(""), union_unshared_sign=as.numeric(""), trtsignif_unshared_sign=as.numeric(""), union_shared_magnitude=as.numeric(""), union_unshared_magnitude=as.numeric(""), trtsignif_unshared_magnitude=as.numeric(""),union_unshared_magnitude_0.1=as.numeric(""))
for(pair in pairs){
trt = gsub("\t.*","",pair)
ctrl =gsub(".*\t","",pair)
pair = gsub("\t",":",pair)
# by sign:
## l[trt]= sum((mlfsr[,trt]<0.1 | mlfsr[,ctrl]<0.1) & mean[,trt]*mean[,ctrl]>0)
df[pair,"union_shared_sign"] <- length(unique(gsub("_.*","",rownames(mlfsr[(mlfsr[,trt]<0.1 | mlfsr[,ctrl]<0.1) & mean[,trt]*mean[,ctrl]>0,]))))
df[pair,"union_unshared_sign"] <-length(unique(gsub("_.*","",rownames(mlfsr[(mlfsr[,trt]<0.1 | mlfsr[,ctrl]<0.1) & mean[,trt]*mean[,ctrl]<0,]))))
df[pair,"trtsignif_unshared_sign"] <- length(unique(gsub("_.*","",rownames(mlfsr[(mlfsr[,trt]<0.1) & mean[,trt]*mean[,ctrl]<0,]))))
# by magnitude:
df[pair,"union_shared_magnitude"] <- length(unique(gsub("_.*","",rownames(mlfsr[((mlfsr[,trt]<0.1 | mlfsr[,ctrl]<0.1) & (mean[,trt]/mean[,ctrl]<=2 & mean[,trt]/mean[,ctrl]>=0.5)),]))))
df[pair,"trtsignif_unshared_magnitude"] <- length(unique(gsub("_.*","",rownames(mlfsr[(mlfsr[,trt]<0.1) & (mean[,trt]/mean[,ctrl]>2 | mean[,trt]/mean[,ctrl]<0.5),]))))
df[pair,"union_unshared_magnitude"] <- length(unique(gsub("_.*","",rownames(mlfsr[((mlfsr[,trt]<0.1 | mlfsr[,ctrl]<0.1) & (mean[,trt]/mean[,ctrl]>2 | mean[,trt]/mean[,ctrl]<0.5)),]))))
l5 <- c(l5, unique(gsub("_.*","",rownames(mlfsr[((mlfsr[,trt]<0.1 | mlfsr[,ctrl]<0.1) & (mean[,trt]/mean[,ctrl]>2 | mean[,trt]/mean[,ctrl]<0.5)),]))))
l6 <- c(l6, unique(gsub("_.*","",rownames(mlfsr[(mlfsr[,trt]<0.1) & (mean[,trt]/mean[,ctrl]>2 | mean[,trt]/mean[,ctrl]<0.5),]))))
df[pair,"union_unshared_magnitude_0.1"] <- length(unique(gsub("_.*","",rownames(mlfsr[((mlfsr[,trt]<0.1 | mlfsr[,ctrl]<0.1) & (mean[,trt]/mean[,ctrl]>2 | mean[,trt]/mean[,ctrl]<0.5)) & (abs(mean[,trt]) >0.1 | abs(mean[,ctrl] >0.1)),]))))
}

# save results:
write.table(df, "mashr_sharing_egene.txt", col.names=T, row.names=F, quote=F, sep="\t")

#####

## # change all significant values to 1, and insignificant to 0:
## mlfsr[mlfsr<thresh] <- -1
## mlfsr[mlfsr>thresh] <- 0
## mlfsr[mlfsr==-1] <- 1

## mlfsr$gene <- gsub("_.*","",rownames(mlfsr))
## mlfsr_gene <- data.frame(mlfsr %>% group_by(gene) %>% summarize_each(funs(sum)))
## # drop gene column:
## rownames(mlfsr_gene) <- mlfsr_gene$gene
## mlfsr_gene <- mlfsr_gene[,-1]
## # fix column names:
## colnames(mlfsr_gene) <- gsub("[.]","+",colnames(mlfsr_gene))
## colSums(mlfsr_gene>0)
## sum(rowSums(mlfsr_gene>0)>0)

## cells <- unique(gsub("_.*","",colnames(mlfsr)))
## trts <- unique(gsub(".*_","",colnames(mlfsr)))
## l=list()
## for(cell in cells[1:4]){
##     df <- mlfsr[,grep(cell,colnames(mlfsr))]
##     print(cell)
##        sdf <- data.frame(summary(as.factor(rowSums(df<0.1))))
##     sdf$count <- rownames(sdf)
##     l[[cell]] <- sdf
## }
## for(trt in trts[1:5]){
##     df <- mlfsr[,grep(trt,colnames(mlfsr))]
##     print(trt)
##    sdf <- data.frame(summary(as.factor(rowSums(df<0.1))))
##     sdf$count <- rownames(sdf)
##     l[[trt]] <- sdf
##     }


## gl=list()
## for(cell in cells[1:4]){
##     df <- mlfsr_gene[,grep(cell,colnames(mlfsr_gene))]
##     print(cell)
##        sdf <- data.frame(summary(as.factor(rowSums(df==1))))
##     sdf$count <- rownames(sdf)
##     gl[[cell]] <- sdf
## }
## for(trt in trts[1:5]){
##     df <- mlfsr_gene[,grep(trt,colnames(mlfsr_gene))]
##     print(trt)
##    sdf <- data.frame(summary(as.factor(rowSums(df==1))))
##     sdf$count <- rownames(sdf)
##     gl[[trt]] <- sdf
##     }


## # get number of unique to each column:
## mlfsr_o <- mlfsr_gene[,1:20]

## for(i in 1:ncol(mlfsr_o)){
##     print(colnames(mlfsr_o)[i])
## print(    sum(mlfsr_o[,i]==1 & rowSums(mlfsr_o[,-i])==0 ))
##     }


## # load the ashr results:
## load(paste0("./mashr_input.Rd"))
## lfsr <- lfsrs

## # change all significant values to 1, and insignificant to 0:
## lfsr[lfsr<thresh] <- -1
## lfsr[lfsr>thresh] <- 0
## lfsr[lfsr==-1] <- 1

## # drop all 0 rows for faster processing:
## lfsr <- lfsr[rowSums(lfsr,na.rm=T)>0,]
## colSums(lfsr==1,na.rm=T)

## # summarize the data per-gene:
## lfsr$gene <- gsub("_.*","",rownames(lfsr))
## lfsr_gene <- data.frame(lfsr %>% group_by(gene) %>% summarize_each(funs(sum)))
## # drop gene column:
## rownames(lfsr_gene) <- lfsr_gene$gene
## lfsr_gene <- lfsr_gene[,-1]
## # fix column names:
## colnames(lfsr_gene) <- gsub("[.]","+",colnames(lfsr_gene))
## colSums(lfsr_gene>0, na.rm = T)

## sum(rowSums(lfsr_gene, na.rm = T)>0)
## length(unique(gsub("_.*","",rownames(mlfsr[rowSums(mlfsr[,1:20], na.rm = T)>0,]))))

## ## # plot matrix of correlations:
## ## # set up colors:
## ## ann_colors = list(cell=c("Tcell"="#e41a1c", "NKcell"="#377eb8", "Bcell"="#4daf4a", "Monocyte"="#984ea3"),trt=c("CTRL"="#828282", "LPS"="#fb9a99", "LPS.EtOH"="#fb9a99", "LPS-EtOH"="#fb9a99", "LPS.DEX"="#e31a1c","LPS-DEX"="#e31a1c","LPS+DEX"="#e31a1c","PHA"="#a6cee3", "PHA.EtOH"="#a6cee3", "PHA-EtOH"="#a6cee3", "PHA.DEX"="#1f78b4", "PHA-DEX"="#1f78b4", "PHA+DEX"="#1f78b4"))
## ## annotation_col <- data.frame(matrix(nrow=20,ncol=2))
## ## rownames(annotation_col) <- colnames(lfsr)
## ## colnames(annotation_col) <- c("cell","trt")
## ## annotation_col$trt <- gsub(".*_","",rownames(annotation_col))
## ## annotation_col$cell <- gsub("_.*","",rownames(annotation_col))
## ## annotation_col$cell <- gsub("*_.*","",annotation_col$cell)
## ## library(pheatmap)
## ## pdf("./plots/mashr_pheatmap_corr_matrix_all.pdf",width=11, height=10)
## ## # calculate correlations:
## ## corr <- cor(lfsr, use='pairwise.complete.obs')
## ## library(psych)
## ## corr.psych <- corr.test(lfsr[,], adjust="none")
## ## # blank out non-significant correlations:
## ## corr <- corr*(corr.psych$p<0.05)
## ## pheatmap(corr, annotation_col = annotation_col, annotation_colors = ann_colors,cellheight=21,cellwidth=21,cluster_rows =F,cluster_cols =F)
## ## dev.off()




## ### END 11/30/2020

## sessionInfo()

