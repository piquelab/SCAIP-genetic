# this script fits mash on all the preprocessed SCAIP FastQTL eQTL mapping results and saves the model to be run on chunked data
# 10/9/2020 JR
# based on ./mashr.R
# last edited and ran 11/24/2020 JR

library(ashr)
library(mashr)
library(data.table)
library(dplyr)
library(doParallel)
cores <- as.integer(Sys.getenv("SLURM_STEP_TASKS_PER_NODE"))
registerDoParallel(cores = cores)
library(RhpcBLASctl)
blas_set_num_threads(12)

folder = "SCAIP1-6"
pcs = 4

## # 1. read in all the data and convert into a mashr object:
## slope_files <-list.files(path=paste0("input/",folder),pattern=".*_slope.txt",full.name=T)
## ## slope_files <- grep("*narrow*",slope_files, invert=TRUE, value=T)
## datasets <- gsub("_slope.txt","",list.files(path=paste0("input/",folder),pattern=".*_slope.txt"))
## SE_files <- list.files(path=paste0("input/",folder),pattern=".*_SE.txt", full.name=T)
## ## SE_files <- grep("*narrow*",SE_files, invert=TRUE, value=T)
## p_files <- list.files(path=paste0("input/",folder),pattern=".*_pvalue.txt", full.name=T)
## lfsr_files <- list.files(path=paste0("input/",folder),pattern=".*lfsr.txt",full.name=T)
## ## p_files <- grep("*narrow*",p_files, invert=TRUE, value=T)
## # read in all the slopes:
## ldf_slope <- lapply(slope_files, read.table, sep="\t", header=T)
## ldf_SE <- lapply(SE_files, read.table, sep="\t", header=T)
## ldf_p <- lapply(p_files, read.table, sep="\t", header=T)
## ldf_lfsr <- lapply(lfsr_files, read.table, sep="\t", header=T)
## # merge them into one data frame:
## # reduced_slope <- Reduce(merge(ldf_slope)) # WHY won't this work??!
## slopes <- Reduce(function(x,y) merge(x = x, y = y, all=TRUE), ldf_slope)
## # remove the duplicated row for now:
## slopes <- slopes[!duplicated(slopes$uniqID),]
## rownames(slopes) <- slopes[,1]
## slopes <- slopes[,-1]
## SEs <- Reduce(function(x,y) merge(x = x, y = y, all=TRUE), ldf_SE)
## # remove the duplicated row for now:
## SEs <- SEs[!duplicated(SEs$uniqID),]
## rownames(SEs) <- SEs[,1]
## SEs <- SEs[,-1]
## pvalues <- Reduce(function(x,y) merge(x = x, y = y, all=TRUE), ldf_p)
## # remove the duplicated row for now:
## pvalues <- pvalues[!duplicated(pvalues$uniqID),]
## rownames(pvalues) <- pvalues[,1]
## pvalues <- pvalues[,-1]
## lfsrs <- Reduce(function(x,y) merge(x = x, y = y, all=TRUE), ldf_lfsr)
## # remove the duplicated row for now:
## lfsrs <- lfsrs[!duplicated(lfsrs$uniqID),]
## rownames(lfsrs) <- lfsrs[,1]
## lfsrs <- lfsrs[,-1]
## # save it for fututre use since this step takes forever:
## save(slopes, SEs, pvalues, lfsrs,file=paste0("./mashr_input_",folder,".Rd"))

load(paste0("./mashr_input_",folder,".Rd"))
# 48% of rows have some NAs (51% for SCAIP1-6), but I think it doesn't work with NAs
# so let's remove all NAs as well as Inf from SEs:
SEs <- SEs[!rowSums(is.finite(as.matrix(SEs)))<ncol(SEs),]
# remove the same rows from pvalue and slope matrixes:
slopes <- slopes[rownames(SEs),]
pvalues <- pvalues[rownames(SEs),]
# check:
sum(rowSums(is.na(as.matrix(slopes))>0))
sum(rowSums(is.na(as.matrix(pvalues))>0))
stopifnot(identical(rownames(slopes),rownames(SEs)))
stopifnot(identical(rownames(slopes),rownames(pvalues)))
# dim(SEs)
# 2601493      20
# let's try removing rows with >15 NAs:
## slopes <- slopes[!rowSums(is.na(slopes))>5,]
## SEs <- SEs[!rowSums(is.na(SEs))>5,]
## pvalues <- pvalues[!rowSums(is.na(pvalues))>5,]

# save the rownames:
write.table(rownames(pvalues),paste0("rownames_mash_",folder,".txt"),row.names=F,col.names=F,quote=F)

data.temp = mash_set_data(as.matrix(slopes), as.matrix(SEs))
# calculate Vhat:
Vhat = estimate_null_correlation_simple(data.temp)
data = mash_set_data(as.matrix(slopes), as.matrix(SEs),V=Vhat)

# 2. Set up the DATA-DRIVEN covariance matrices # NEXT TRY USING BOTH OR JUST CANONICAL
# 2.1. select strong signals
m.1by1 = mash_1by1(data)
strong = get_significant_results(m.1by1,0.05)
# 2.2. Obtain initial data-driven covariance matrices
# this step estimates the covariances of the observed data
U.pca = cov_pca(data,5,subset=strong)
print(names(U.pca))
# canonical matrix:
U.c = cov_canonical(data) 
# 2.3. Apply Extreme Deconvolution
# this step estimates the covariances of the actual underlying effects
U.ed = cov_ed(data, U.pca, subset=strong)

# Step 3: fit the model
# this fits a mixture model to the data, estimating the mixture proportions
## Sys.time()
## m.ed = mash(data, U.ed,outputlevel=1)

Sys.time()
m   = mash(data, c(U.c,U.ed),outputlevel=1)
Sys.time()
save(data,m,Vhat,file=paste0("mash-model-fit_",folder,".Rd"))

# save the colnames:
write.table(colnames(pvalues),"FastQTL_all_conditions.txt",sep="\n",col.names=F, row.names=F, quote=F)

### END 11/25/2020


## # select a good number of chunks:
## tests <- nrow(data$Bhat)
## for(i in 1:round(tests/2)){
## div <- tests/i
## if(floor(div)==div){
##     print(paste0(i,"  ",div))
##         }
## }

### END 10/9/2020
