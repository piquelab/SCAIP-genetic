# this script takes mashr results and makes upset plot of genes in common
# 12/8/2020 JR
# last edited and ran 12/9/2020 JR: added subsetting to 1) shared sign 2) shared sign and magnitude when assessing sharing

library(data.table)
library(tidyr)
library(UpSetR)
## library(ggupset)

thresh <- 0.1

# get the colors for plotting:
anno_colors <- list(cell=c("Bcell"="#4daf4a", "Monocyte"="#984ea3","NKcell"="#AA4B56", "Tcell"="#FFAA00"), condition=c("CTRL"="#828282","LPS"="#fb9a99", "LPS+DEX"="#e31a1c","PHA"="#a6cee3", "PHA+DEX"="#1f78b4"), Batch=c("SCAIP1"="#e41a1c", "SCAIP2"="#377eb8","SCAIP3"="#4daf4a", "SCAIP4"="#984ea3", "SCAIP5"="#ff7f00", "SCAIP6"="#ffff33"))

# load the mashr results:
lfsr <- fread("output/lfsr/lfsr_mashr_all.txt.gz",data.table=F)
mean <- fread("output/posterior_mean/posterior_mean_mashr_all.txt.gz",data.table=F)

# add colnames:
cols <- read.table("FastQTL_all_conditions_proper.txt",stringsAsFactors=F)[,1]
colnames(lfsr) <- cols
colnames(mean) <- cols
# add rownames:
rows <- read.table("rownames_mash_SCAIP1-6.txt",stringsAsFactors=F)[,1]
rownames(lfsr) <- rows
rownames(mean) <- rows

# modify the data for upset function:
# change all significant values to 1, and insignificant to 0:
lfsr[lfsr<thresh] <- -1
lfsr[lfsr>thresh] <- 0
lfsr[lfsr==-1] <- 1
# drop the all all insignificant rows for faster processing:
lfsr <- lfsr[rowSums(lfsr==0)!=ncol(lfsr),]
# and from mean:
mean <- mean[rownames(lfsr),]

# change insignificant "means" to 0s:
mean <- mean*lfsr
# get the sign of the mean:
sign <- data.frame(apply(mean, c(1,2), sign))
# fix the colnames:
colnames(sign) <- gsub("[.]","+",colnames(sign))

# make the annotation:
anno <- data.frame(sets=colnames(lfsr))
anno$condition <- gsub(".*_","",anno$sets)
anno$cell <- gsub("_.*","",anno$sets)
rownames(anno) <- anno$sets

# make a per-eQTL plot:
pdf("plots/upset_eQTL_mashr-all.pdf",height=12, width=9)
upset(lfsr,nsets=20,order.by = "freq", number.angles = 90, show.number="no", text.scale=c(1.8, 1.4, 1.4, 1.2, 1.45, 1.5),
set.metadata =
    list(
        data = anno,
        plots =
            list(
                ## list(type = "hist",
                ##      column = "trt",
                ##      assign = 20
                ##      ),
                list(type = "matrix_rows",
                     column = "trt", colors = anno_colors$condition,alpha=1
                     )))
      )
## ggplot(lfsr,aes(x=))
dev.off()

pdf("plots/upset_eQTL_mashr-all_unordered.pdf",height=12, width=9)
upset(lfsr,nsets=20,order.by = "freq", number.angles = 90, show.number="no", text.scale=c(1.8, 1.4, 1.4, 1.2, 1.45, 1.5), keep.order = T,sets=rev(sort(rownames(anno))),
set.metadata =
    list(
        data = anno,
        plots =
            list(
                ## list(type = "hist",
                ##      column = "trt",
                ##      assign = 20
                ##      ),
                list(type = "matrix_rows",
                     column = "trt", colors = anno_colors$condition,alpha=1
                     )))
      )
## ggplot(lfsr,aes(x=))
dev.off()



# summarize the data per-gene:
lfsr$gene <- gsub("_.*","",rownames(lfsr))
lfsr_gene <- data.frame(lfsr %>% group_by(gene) %>% summarize_each(funs(sum)))
# drop gene column:
rownames(lfsr_gene) <- lfsr_gene$gene
lfsr_gene <- lfsr_gene[,-1]
# fix column names:
colnames(lfsr_gene) <- gsub("[.]","+",colnames(lfsr_gene))

# set all >0 to 1:
lfsr_gene[lfsr_gene>0] <- 1


pdf("plots/upset_eGene_mashr-all_unordered.pdf",height=12, width=9)
upsetplot=upset(lfsr_gene,nsets=20,order.by = "freq", number.angles = 90, show.number="no", text.scale=c(1.8, 1.4, 1.4, 1.2, 1.45, 1.5), keep.order = T,sets=rev(sort(rownames(anno))),
set.metadata =
    list(
        data = anno,
        plots =
            list(
                ## list(type = "hist",
                ##      column = "trt",
                ##      assign = 20
                ##      ),
                list(type = "matrix_rows",
                     column = "trt", colors = anno_colors$condition,alpha=1
                     )))
      )
## ggplot(lfsr,aes(x=))
dev.off()

#### added sign/magnitude sharing filters (12/9/2020):

# to get over the fact that upset() can only handle binary data, let's binarize the sign matrix as follows:
# extract all -1 as one matrix:
signneg <- sign
signneg[signneg>-1] <- 0
# convert -1 to 1:
signneg[signneg==-1] <- 1
# extract all 1 as one matrix:
signpos <- sign
signpos[signpos<1] <- 0
# concatenate the two:
signupset <- rbind(signneg,signpos)

# make a per-eQTL plot:
pdf("plots/upset_eQTL_mashr-all_unordered_sign_numbers.pdf",height=12, width=18)
upset(signupset,nsets=20,order.by = "freq", number.angles = 0, show.number="yes", text.scale=c(1.8, 1.4, 1.4, 1.2, 1.45, 1.5),keep.order = T,sets=rev(sort(rownames(anno))),
set.metadata =
    list(
        data = anno,
        plots =
            list(
                ## list(type = "hist",
                ##      column = "trt",
                ##      assign = 20
                ##      ),
                list(type = "matrix_rows",
                     column = "trt", colors = anno_colors$condition,alpha=1
                     )))
      )
## ggplot(lfsr,aes(x=))
dev.off()



## now make a heatmap for pairwise sharing of sign and magnitue(0.5<FC<2):
shared <- apply(mean,2,function(x)colSums(abs(x)<2*abs(mean) & abs(x)>0.5*abs(mean) & x!=0 & mean !=0 & sign(x)==sign(mean)))
# this is how you would calculate Jaccard index if shared was just an intersection
union <- apply(mean,2,function(x)colSums(x!=0 | mean!=0))
fakejacc <- apply(mean,2,function(x)colSums(abs(x)<2*abs(mean) & abs(x)>0.5*abs(mean) & x!=0 & mean !=0 & sign(x)==sign(mean)) / (colSums(x!=0 | mean!=0  )))
              

library(RColorBrewer)
# heatmap:
library(pheatmap)

pdf("plots/pheatmap_shared-by-magnitude2_eQTLs_unclustered.pdf")
fig1 <- pheatmap(shared,
               col=colorRampPalette(c("#FAE5D8","#B2182B"))(40),
                               scale="none",
                               border_color="NA",
                               cluster_rows=F, cluster_cols=F,
                               legend=T,
                               annotation_col=anno[,c("condition","cell")],
                               annotation_colors=anno_colors,
                               show_colnames=T, show_rownames=F,
                               na_col="white",
                 fontsize = 7, fontsize_row = 10, fontsize_col = 10,
                                  )
                 ## cell_width=10,
                 ## cellheight=10)
fig1
dev.off()


pdf("plots/pheatmap_shared-by-magnitude2_eQTLs_fakeJaccard_unclustered.pdf")
fig1 <- pheatmap(fakejacc,
               col=colorRampPalette(c("#FAE5D8","#B2182B"))(40),
                               scale="none",
                               border_color="NA",
                               cluster_rows=F, cluster_cols=F,
                               legend=T,
                               annotation_col=anno[,c("condition","cell")],
                               annotation_colors=anno_colors,
                               show_colnames=T, show_rownames=F,
                               na_col="white",
                 fontsize = 7, fontsize_row = 10, fontsize_col = 10,
                                  )
                 ## cell_width=10,
                 ## cellheight=10)
fig1
dev.off()

pdf("plots/pheatmap_shared-by-magnitude2_eQTLs_fakeJaccard_clustered.pdf")
fig1 <- pheatmap(fakejacc,
               col=colorRampPalette(c("#FAE5D8","#B2182B"))(40),
                               scale="none",
                               border_color="NA",
                               cluster_rows=T, cluster_cols=T,
                               legend=T,
                               annotation_col=anno[,c("condition","cell")],
                               annotation_colors=anno_colors,
                               show_colnames=T, show_rownames=F,
                               na_col="white",
                 fontsize = 7, fontsize_row = 10, fontsize_col = 10,
                                  )
                 ## cell_width=10,
                 ## cellheight=10)
fig1
dev.off()
