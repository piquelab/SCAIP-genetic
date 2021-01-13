# this script runs mashr on prepped eQTL results to find response and cell-type-specific eQTLs
# 1/28/2020 JR
# based on mashr vignette https://stephenslab.github.io/mashr/articles/intro_mash.html
# last editted and ran 12/8/2020 JR: now saving all the matrices of mash output

# prereq: ./mashr-fit-model.R

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

folder = "SCAIP1-6"
nchunks <- 500
chunk <- 5

args = commandArgs(trailingOnly=TRUE)
# get the PC number:
if (length(args)>0){
    folder=args[1]
    chunk=as.numeric(args[2])
    nchunks <- as.numeric(args[3])
          }

# load the mashr prefit model:
## load(paste0("mash-model-fit_",folder,".Rd"))
load(paste0("mash-model-fit_",folder,".Rd"))

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
write.table(m.chunk$lfsr,paste0("output/lfsr/all_chunk", chunk,".txt"),sep="\t",quote=F, row.names=F,col.names=F)
# save the posterior means results:
write.table(m.chunk$PosteriorMean,paste0("output/posterior_mean/all_chunk", chunk,".txt"),sep="\t",quote=F, row.names=F,col.names=F)
# save the posterior SDs results:
write.table(m.chunk$PosteriorSD,paste0("output/posterior_SD/all_chunk", chunk,".txt"),sep="\t",quote=F, row.names=F,col.names=F)
write.table(m.chunk$lfdr,paste0("output/lfdr/all_chunk", chunk,".txt"),sep="\t",quote=F, row.names=F,col.names=F)
write.table(m.chunk$NegativeProb,paste0("output/NegativeProb/all_chunk", chunk,".txt"),sep="\t",quote=F, row.names=F,col.names=F)

# get pairwise sign sharing:
# 1. just the sign:
## share_sign <- get_pairwise_sharing(m.chunk, factor=0)
# doesn't work. I wonder why

### END 12/8/2020 JR

## # Step 3: fit the model
## # this fits a mixture model to the data, estimating the mixture proportions
## m.ed = mash(data, U.ed) #ISSUE: not done after two weeks


## #### END 10/7/2020 (never got past this step)
## m   = mash(data, c(U.c,U.ed))
## save(m,m.ed,file="m_m.ed_Tcell.Rd")
## # Step 4: Extract Posterior Summaries:
## # posterior mean:
## head(get_lfsr(m.ed))
## # posterior standard deviation:
## head(get_pm(m.ed))
## # local false sign rate (measure of significance)
## head(get_psd(m.ed))


sessionInfo()

