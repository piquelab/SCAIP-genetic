# this script fits mash on all the preprocessed SCAIP FastQTL eQTL mapping results on mean (from NB) and saves the model to be run on chunked data
# based on ../mashr_eQTL/mashr-fit-model.R
# 2/15/2021 JR

library(ashr)
library(mashr)
library(data.table)
library(dplyr)
library(doParallel)
cores <- as.integer(Sys.getenv("SLURM_STEP_TASKS_PER_NODE"))
registerDoParallel(cores = cores)
library(RhpcBLASctl)
blas_set_num_threads(12)

pcs = 8

## # done once and for all 2/15/2021:
## # 1. read in all the data and convert into a mashr object:
## slope_files <-list.files(path=paste0("input/"),pattern=".*_slope.txt",full.name=T)
## datasets <- gsub("_slope.txt","",list.files(path=paste0("input/"),pattern=".*_slope.txt"))
## SE_files <- list.files(path=paste0("input/"),pattern=".*_SE.txt", full.name=T)
## p_files <- list.files(path=paste0("input/"),pattern=".*_pvalue.txt", full.name=T)
## lfsr_files <- list.files(path=paste0("input/"),pattern=".*lfsr.txt",full.name=T)
## # read in all the files:
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
## save(slopes, SEs, pvalues, lfsrs,file=paste0("./mashr_input.Rd"))

load(paste0("./mashr_input.Rd"))
# 31.5% of rows have some NAs, but I think it doesn't work with NAs
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

# save the rownames:
write.table(rownames(pvalues),paste0("rownames_mash.txt"),row.names=F,col.names=F,quote=F)

# sample 20k random tests:
set.seed(12)
sset <- sample(c(1:nrow(slopes)),20000)
data.temp = mash_set_data(as.matrix(slopes[sset,]), as.matrix(SEs[sset,]))
# calculate Vhat:
Vhat = estimate_null_correlation_simple(data.temp)
rm(data.temp)
data = mash_set_data(as.matrix(slopes), as.matrix(SEs),V=Vhat)
data.random <- mash_set_data(as.matrix(slopes[sset,]), as.matrix(SEs[sset,]),V=Vhat)


# 2. Set up the DATA-DRIVEN covariance matrices # NEXT TRY USING BOTH OR JUST CANONICAL
# 2.1. select strong signals
m.1by1 = mash_1by1(data)
strong = get_significant_results(m.1by1,0.05)
# 2.2. Obtain initial data-driven covariance matrices
# this step estimates the covariances of the observed data
U.pca = cov_pca(data,5,subset=strong)
print(names(U.pca))
# canonical matrix:
U.c = cov_canonical(data.random) 
# 2.3. Apply Extreme Deconvolution
# this step estimates the covariances of the actual underlying effects
U.ed = cov_ed(data, U.pca, subset=strong)

# save all objects created up to now:
save(data,U.c,U.ed,Vhat, file="mashr-fit-model_objects.Rd")

# Step 3: fit the model
# this fits a mixture model to the data, estimating the mixture proportions
## Sys.time()
## m.ed = mash(data, U.ed,outputlevel=1)

Sys.time()
m   = mash(data.random, c(U.c,U.ed),outputlevel=1)
Sys.time()
save(data,data.random,m,Vhat,file=paste0("mash-model-fit.Rd"))
# 8578 * 'error: chol(): decomposition failed'
# Warning message:
#In calc_lik_matrix(data, Ulist, log = TRUE, algorithm.version = algorithm.version) :
#  Some mixture components result in non-finite likelihoods, either
# due to numerical underflow/overflow, or due to invalid covariance matrices# 959, 990, 1021, 1052 


# save the colnames:
write.table(colnames(pvalues),"FastQTL_all_conditions-colnames.txt",sep="\n",col.names=F, row.names=F, quote=F)

### END 2/15/2021


## # select a good number of chunks:
## tests <- nrow(data$Bhat)
## for(i in 1:round(tests/2)){
## div <- tests/i
## if(floor(div)==div){
##     print(paste0(i,"  ",div))
##         }
## }

