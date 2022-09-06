# this script gets the pair-wise unique dispersion eQTLs ~reQTLs
# based on ../mashr_eQTL/mashr_get_reQTLs.R
# 2/10/2021 JR


library(data.table)
library(dplyr)

thresh <- 0.1

# load the mashr results:
mlfsr <- fread("output/lfsr/lfsr_mashr_all.txt.gz", sep="\t",header=F, stringsAsFactors=F,data.table=F)
rownames(mlfsr) <- mlfsr[,1]
mlfsr <- mlfsr[,-1]
# add colnames:
cols <- read.table("FastQTL_all_conditions-colnames.txt",stringsAsFactors=F)[,1]
colnames(mlfsr) <- gsub("[.]","-",cols)
mean <- fread("output/posterior_mean/posterior_mean_mashr_all.txt.gz",data.table=F)
rownames(mean) <- mean[,1]
mean <- mean[,-1]
colnames(mean) <- gsub("[.]","-",cols)

# get pair-wise unique eQTLs:
pairs <- read.table("../mashr_eQTL/contrasts.txt",sep="\n",stringsAsFactors=F)[,1]

# get sharing and specificity:
uum=list()
top <- list()
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
df <- data.frame(row.names=gsub("\t","_",pairs),trt=gsub("\t.*","",pairs), ctrl =gsub(".*\t","",pairs), union_shared_sign=as.numeric(""), union_unshared_sign=as.numeric(""), trtsignif_unshared_sign=as.numeric(""), union_shared_magnitude=as.numeric(""), union_unshared_magnitude=as.numeric(""), trtsignif_unshared_magnitude=as.numeric(""))
for(pair in pairs){
trt = gsub("\t.*","",pair)
ctrl =gsub(".*\t","",pair)
pair = gsub("\t","_",pair)
# by sign:
## l[trt]= sum((mlfsr[,trt]<0.1 | mlfsr[,ctrl]<0.1) & mean[,trt]*mean[,ctrl]>0)
df[pair,"union_shared_sign"] <- length(unique(gsub("_.*","",rownames(mlfsr[(mlfsr[,trt]<0.1 | mlfsr[,ctrl]<0.1) & mean[,trt]*mean[,ctrl]>0,]))))
df[pair,"union_unshared_sign"] <-length(unique(gsub("_.*","",rownames(mlfsr[(mlfsr[,trt]<0.1 | mlfsr[,ctrl]<0.1) & mean[,trt]*mean[,ctrl]<0,]))))
df[pair,"trtsignif_unshared_sign"] <- length(unique(gsub("_.*","",rownames(mlfsr[(mlfsr[,trt]<0.1) & mean[,trt]*mean[,ctrl]<0,]))))
# by magnitude:
df[pair,"union_shared_magnitude"] <- length(unique(gsub("_.*","",rownames(mlfsr[((mlfsr[,trt]<0.1 | mlfsr[,ctrl]<0.1) & (mean[,trt]/mean[,ctrl]<=2 & mean[,trt]/mean[,ctrl]>=0.5)),]))))
df[pair,"union_unshared_magnitude"] <- length(unique(gsub("_.*","",rownames(mlfsr[((mlfsr[,trt]<0.1 | mlfsr[,ctrl]<0.1) & (mean[,trt]/mean[,ctrl]>2 | mean[,trt]/mean[,ctrl]<0.5)),]))))
df[pair,"trtsignif_unshared_magnitude"] <- length(unique(gsub("_.*","",rownames(mlfsr[(mlfsr[,trt]<0.1) & (mean[,trt]/mean[,ctrl]>2 | mean[,trt]/mean[,ctrl]<0.5),]))))
l5 <- c(l5, unique(gsub("_.*","",rownames(mlfsr[((mlfsr[,trt]<0.1 | mlfsr[,ctrl]<0.1) & (mean[,trt]/mean[,ctrl]>2 | mean[,trt]/mean[,ctrl]<0.5)),]))))
l6 <- c(l6, unique(gsub("_.*","",rownames(mlfsr[(mlfsr[,trt]<0.1) & (mean[,trt]/mean[,ctrl]>2 | mean[,trt]/mean[,ctrl]<0.5),]))))
}

# save results:
write.table(df, "mashr_sharing_egene.txt", col.names=T, row.names=F, quote=F, sep="\t")

##### END 2/10/2021
