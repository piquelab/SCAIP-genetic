# this script runs mashr on prepped eQTL results to find response and cell-type-specific eQTLs
# based on ../mashr_eQTL/mashr_compute_posterior_chunks.R
# prereq: ./mashr-fit-model.R
# 2/5/2021 JR

library(ashr)
library(mashr)
library(data.table)
library(dplyr)
library(doParallel)
cores <- as.integer(Sys.getenv("SLURM_STEP_TASKS_PER_NODE"))
registerDoParallel(cores = cores)
library(RhpcBLASctl)
blas_set_num_threads(cores)
library("profmem")

nchunks <- 500
chunk <- 5

args = commandArgs(trailingOnly=TRUE)
# get the PC number:
if (length(args)>0){
    chunk=as.numeric(args[1])
    nchunks <- as.numeric(args[2])
          }

chunk

# load the mashr prefit model:
load(paste0("mash-model-fit.Rd"))

# subset the data to the current chunk:
chunk.size = round(nrow(data$Bhat)/nchunks)

if(chunk==nchunks){
subset <- ((chunk-1)*chunk.size+1):nrow(data$Bhat)    
    }else{
subset <- ((chunk-1)*chunk.size+1):(chunk*chunk.size)
        }

# subset the data:
data.chunk <- mash_set_data(data$Bhat[subset,],data$Shat[subset,],V=Vhat)

# compute posterior matrices:
Sys.time()
m.chunk <- mash_compute_posterior_matrices(g=m,data=data.chunk)
Sys.time()

# save the lfsr results:
write.table(m.chunk$lfsr,paste0("output/lfsr/all_chunk", chunk,".txt"),sep="\t",quote=F, row.names=T,col.names=T)
# save the posterior means results:
write.table(m.chunk$PosteriorMean,paste0("output/posterior_mean/all_chunk", chunk,".txt"),sep="\t",quote=F, row.names=T,col.names=T)
# save the posterior SDs results:
write.table(m.chunk$PosteriorSD,paste0("output/posterior_SD/all_chunk", chunk,".txt"),sep="\t",quote=F, row.names=T,col.names=T)
write.table(m.chunk$lfdr,paste0("output/lfdr/all_chunk", chunk,".txt"),sep="\t",quote=F, row.names=T,col.names=T)
write.table(m.chunk$NegativeProb,paste0("output/NegativeProb/all_chunk", chunk,".txt"),sep="\t",quote=F, row.names=T,col.names=T)


### END 2/5/2021 JR

sessionInfo()

