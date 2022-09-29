###
library(tidyverse)
library(data.table)
library(parallel)
library(qvalue)
library(qqman)
library(ggplot2)
library(grid)
library(gridExtra)
library(ggExtra)
#library(Seurat)
library(annotables)
library(ggrastr)

option <- "DiagLDA2"
outdir <- paste(option, "/summary/", sep="")
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)


#########################################################
### qq plots for LDA-interaction eQTL mapping results ###
#########################################################

contrast_ls <- list("LPS"=c("CTRL", "LPS-EtOH"),
   "LPS-DEX"=c("LPS-EtOH", "LPS-DEX"),
   "PHA"=c("CTRL", "PHA-EtOH"),
   "PHA-DEX"=c("PHA-EtOH", "PHA-DEX"))


dataset <- read.table("dataset_contrast.txt", header=F)
names(dataset) <- c("MCls", "LDA", "treats")

plotdata <- map_dfr(1:nrow(dataset), function(i){
###
   cell <- dataset[i,1]
   lda <- dataset[i,2]
   treat <- dataset[i,3] 
   fn <- paste0(option, "/4.2_lm.results/2_", cell, "_lda", lda, "_trt", treat, "_stratified.FDR.txt")
   res <- fread(fn, data.table=F)
    
   nsnp <- nrow(res)
   res2 <- res%>%dplyr::select(ENSG, varID, gene_SNP, interaction_ANOVA.pvalue)%>%
     arrange(interaction_ANOVA.pvalue)%>%
     mutate(observed=-log10(interaction_ANOVA.pvalue),
          expected=-log10(ppoints(nsnp)),
          MCls=cell, LDA=lda, treats=treat)
   res2
})


###
### QQ plots
plotdata <- plotdata%>%mutate(treat2=gsub("-EtOH","",treats))

p <- ggplot(plotdata, aes(x=expected, y=observed))+
   rasterize(geom_point(colour="grey30", size=0.2),dpi=300)+
   facet_grid(MCls~treat2, scales="free",
      labeller=labeller(treat2=c("LPS"="LPS", "LPS-DEX"="LPS+DEX",
                                 "PHA"="PHA", "PHA-DEX"="PHA+DEX")))+ 
   geom_abline(colour="red")+
   labs(x=bquote("Expected"~ -log[10]~"("~plain(P)~")"),
        y=bquote("Observed"~-log[10]~"("~plain(P)~")")) +
   theme_bw()+
   theme(axis.title=element_text(size=10),
         strip.text=element_text(size=12))

figfn <- paste(outdir, "Figure1_qq.anova.png", sep="")
png(figfn, width=720, height=750, res=120)
print(p)
dev.off()

###
figfn <- paste(outdir, "Figure1_qq.anova.pdf", sep="")
pdf(figfn, width=7.5, height=8)
print(p)
dev.off()



##########################################################
### QQ plots for LDA-interaction eQTL mapping, colored ###
##########################################################

## rm(list=ls())

contrast_ls <- list("LPS"=c("CTRL", "LPS-EtOH"),
   "LPS-DEX"=c("LPS-EtOH", "LPS-DEX"),
   "PHA"=c("CTRL", "PHA-EtOH"),
   "PHA-DEX"=c("PHA-EtOH", "PHA-DEX"))

option <- "DiagLDA2"

dataset <- read.table("dataset_contrast.txt", header=F)
names(dataset) <- c("MCls", "LDA", "treats")

### calculate qq
qq.fun <- function(df){
   nsnp <- nrow(df)
   df <- df%>%arrange(interaction_ANOVA.pvalue)%>%
     mutate(observed=-log10(interaction_ANOVA.pvalue),
            expected=-log10(ppoints(nsnp)))
   df         
}



##########################################
### qq plots show if enriched in eQTLs ###
##########################################

### read all LDA-interaction results
plotdata <- map_dfr(1:nrow(dataset), function(i){
###
   cell <- dataset[i,1]
   lda <- dataset[i,2]
   treat <- dataset[i,3] 
   fn <- paste0(option, "/4_lm.results/2_", cell, "_lda", lda, "_trt", treat, "_stratified.FDR.txt")
   res <- fread(fn, data.table=F)
   res2 <- res%>%group_by(geneSNP_is_eQTL2)%>%nest()%>%
      mutate(outlist=map(data,~qq.fun(.x)))%>%
      dplyr::select(-data)%>%
      unnest(outlist)%>%as.data.frame()%>%
      mutate(MCls=cell, LDA=lda, treats=treat)       
   res2
})


### QQ plots
plotdata2 <- plotdata%>%mutate(treat2=gsub("-EtOH", "", treats))

lab.treats <- c("LPS"="LPS", "LPS-DEX"="LPS+DEX",
                "PHA"="PHA", "PHA-DEX"="PHA+DEX")
p <- ggplot(plotdata2)+
   rasterize(geom_point(aes(x=expected, y=observed, colour=factor(geneSNP_is_eQTL2)), size=0.3), dpi=300)+
   scale_colour_manual("", values=c("yes"="green", "no"="grey30"),
      labels=c("yes"="eQTL", "no"="not eQTL"),
      guide=guide_legend(override.aes=list(size=2)))+
  geom_abline(colour="red")+
  labs(x=bquote("Expected"~ -log[10]~"("~plain(P)~")"),
       y=bquote("Observed"~-log[10]~"("~plain(P)~")"))+
  facet_grid(MCls~treat2, labeller=labeller(treat2=lab.treats))+
  theme_bw()

figfn <- paste(option, "/summary/Figure2.1_qq.anova_eQTLs.png", sep="")
png(figfn, width=750, height=750, res=120)
print(p)
dev.off()

figfn <- paste(option, "/summary/Figure2.1_qq.anova_eQTLs.pdf", sep="")
pdf(figfn, width=8, height=8)
print(p)
dev.off()


####################################
### if enriched in response eQTL ###
####################################

gene_SNP <- map(1:nrow(dataset), function(i){
  cell <- dataset[i,1]
  lda <- dataset[i,2]
  treat <- dataset[i,3]
  oneX <- contrast_ls[[lda]]
  treat1 <- oneX[2]
  treat0 <- oneX[1]
  fn <- paste("../reQTL_entire-cis/reQTL_lm_results/results_", cell, "_", treat1, "_vs_", treat0, "_stratified_FDR.txt", sep="")
  res <- fread(fn, data.table=F, stringsAsFactors=F)
  gene_SNP <- res%>%drop_na(interaction_ANOVA.qvalue_stratified)%>%
    filter(interaction_ANOVA.qvalue_stratified<0.1)%>%dplyr::pull(gene_SNP)
})
reQTL <- unique(unlist(gene_SNP))

load("../mashr_eQTL/mashr-reQTLs_union-unshared-magnitude2-mlfsr0.1.Rd")
reQTL <- unlist(uum)

### read all LDA-interaction results
plotdata <- map_dfr(1:nrow(dataset), function(i){
###
   cell <- dataset[i,1]
   lda <- dataset[i,2]
   treat <- dataset[i,3] 
   fn <- paste0(option, "/4_lm.results/2_", cell, "_lda", lda, "_trt", treat, "_stratified.FDR.txt")
   res <- fread(fn, data.table=F)
   res <- res%>%mutate(is_reQTL=ifelse(gene_SNP%in%reQTL, "yes", "no"))

   res2 <- res%>%group_by(is_reQTL)%>%nest()%>%
      mutate(outlist=map(data,~qq.fun(.x)))%>%
      dplyr::select(-data)%>%
      unnest(outlist)%>%as.data.frame()%>%
      mutate(MCls=cell, LDA=lda, treats=treat)       
   res2
})


### QQ plots
plotdata2 <- plotdata%>%mutate(treat2=gsub("-EtOH", "", treats))

lab.treats <- c("LPS"="LPS", "LPS-DEX"="LPS+DEX",
                "PHA"="PHA", "PHA-DEX"="PHA+DEX")
p <- ggplot(plotdata2)+
   rasterize(geom_point(aes(x=expected, y=observed, colour=factor(is_reQTL)),size=0.5),dpi=300)+
   scale_colour_manual("", values=c("yes"="green", "no"="grey30"),
      labels=c("yes"="reQTL", "no"="not reQTL"),
      guide=guide_legend(override.aes=list(size=2)))+
  geom_abline(colour="red")+
  labs(x=bquote("Expected"~ -log[10]~"("~plain(P)~")"),
       y=bquote("Observed"~-log[10]~"("~plain(P)~")"))+
  facet_grid(MCls~treat2, labeller=labeller(treat2=lab.treats))+
  theme_bw()

figfn <- paste(option, "/summary/Figure2.2_qq.anova_reGene.png", sep="")
png(figfn, width=800, height=750, res=120)
print(p)
dev.off()

figfn <- paste(option, "/summary/Figure2.2_qq.anova_reGene.pdf", sep="")
pdf(figfn, width=8, height=7.5)
print(p)
dev.off()



###########################################
### qq plots to test if enriched in DEG ###
###########################################

res_DEG <- read_rds("/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/6_DEG.CelltypeNew_output/Filter2/2_meta.rds")

### read all LDA-interaction results
plotdata <- map_dfr(1:nrow(dataset), function(i){
###
   cell <- dataset[i,1]
   lda <- dataset[i,2]
   treat <- dataset[i,3] 
   fn <- paste0(option, "/4_lm.results/2_", cell, "_lda", lda, "_trt", treat, "_stratified.FDR.txt")
   res <- fread(fn, data.table=F)
   DEG <- res_DEG%>%dplyr::filter(abs(beta)>0.5, qval<0.1)%>%
      mutate(MCls==cell, contrast==lda)%>%dplyr::pull(gene)
   res <- res%>%
      mutate(ENSG2=gsub("\\..*", "", ENSG),
             is_DEG=ifelse(ENSG2%in%DEG, "yes", "no"))
   ### 
   res2 <- res%>%group_by(is_DEG)%>%nest()%>%
      mutate(outlist=map(data,~qq.fun(.x)))%>%
      dplyr::select(-data)%>%
      unnest(outlist)%>%as.data.frame()%>%
      mutate(MCls=cell, LDA=lda, treats=treat)       
   res2
})

###
###
plotdata2 <- plotdata%>%mutate(treat2=gsub("-EtOH", "", treats))

lab.treats <- c("LPS"="LPS", "LPS-DEX"="LPS+DEX",
                "PHA"="PHA", "PHA-DEX"="PHA+DEX")
p <- ggplot(plotdata2)+
   rasterize(geom_point(aes(x=expected, y=observed, colour=factor(is_DEG)),size=0.3), dpi=300)+
   scale_colour_manual("", values=c("yes"="green", "no"="grey30"),
      labels=c("yes"="DEG", "no"="not DEG"),
      guide=guide_legend(override.aes=list(size=2)))+
  geom_abline(colour="red")+
  labs(x=bquote("Expected"~ -log[10]~"("~plain(P)~")"),
       y=bquote("Observed"~-log[10]~"("~plain(P)~")"))+
  facet_grid(MCls~treat2, labeller=labeller(treat2=lab.treats))+
  theme_bw()


figfn <- paste(option, "/summary/Figure2.3.1_qq.anova_DEGs.png", sep="")
png(figfn, width=750, height=750, res=120)
print(p)
dev.off()

figfn <- paste(option, "/summary/Figure2.3.1_qq.anova_DEGs.pdf", sep="")
pdf(figfn, width=7.5, height=7.5)
print(p)
dev.off()


####################
### union of DEG ###
####################
res_DEG <- read_rds("/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/6_DEG.CelltypeNew_output/Filter2/2_meta.rds")
DEG <- res_DEG%>%filter(abs(beta)>0.5, qval<0.1)%>%dplyr::pull(gene)
DEG <- unique(DEG)

### read all LDA-interaction results
plotdata <- map_dfr(1:nrow(dataset), function(i){
###
   cell <- dataset[i,1]
   lda <- dataset[i,2]
   treat <- dataset[i,3] 
   fn <- paste0(option, "/4_lm.results/2_", cell, "_lda", lda, "_trt", treat, "_stratified.FDR.txt")
   res <- fread(fn, data.table=F)
   res <- res%>%
      mutate(ENSG2=gsub("\\..*", "", ENSG),
             is_DEG=ifelse(ENSG2%in%DEG, "yes", "no"))
   ### 
   res2 <- res%>%group_by(is_DEG)%>%nest()%>%
      mutate(outlist=map(data,~qq.fun(.x)))%>%
      dplyr::select(-data)%>%
      unnest(outlist)%>%as.data.frame()%>%
      mutate(MCls=cell, LDA=lda, treats=treat)       
   res2
})

###
###
plotdata2 <- plotdata%>%mutate(treat2=gsub("-EtOH", "", treats))

lab.treats <- c("LPS"="LPS", "LPS-DEX"="LPS+DEX",
                "PHA"="PHA", "PHA-DEX"="PHA+DEX")
p <- ggplot(plotdata2)+
   rasterize(geom_point(aes(x=expected, y=observed, colour=factor(is_DEG)),size=0.3), dpi=300)+
   scale_colour_manual("", values=c("yes"="green", "no"="grey30"),
      labels=c("yes"="DEG", "no"="not DEG"),
      guide=guide_legend(override.aes=list(size=2)))+
  geom_abline(colour="red")+
  labs(x=bquote("Expected"~ -log[10]~"("~plain(P)~")"),
       y=bquote("Observed"~-log[10]~"("~plain(P)~")"))+
  facet_grid(MCls~treat2, labeller=labeller(treat2=lab.treats))+
  theme_bw()

figfn <- paste(option, "/summary/Figure2.3.2_qq.anova_DEGunion.png", sep="")
png(figfn, width=750, height=750, res=120)
print(p)
dev.off()


figfn <- paste(option, "/summary/Figure2.3.2_qq.anova_DEGunion.pdf", sep="")
pdf(figfn, width=7.5, height=7.5)
print(p)
dev.off()



#####################################################
### test if dynamical eQTLs are enriched in vQTLs ###
#####################################################

### read vQTLs
cells <- rep(c("Bcell", "Monocyte", "NKcell", "Tcell"), each=5)
treats <- rep(c("CTRL", "LPS-EtOH", "LPS-DEX", "PHA-EtOH", "PHA-DEX"), times=4)
xx <- data.frame(MCls=cells, treats=treats)
QTLs <- map(1:nrow(xx), function(i){
   cell <- xx[i,1]
   treat <- xx[i,2]
   cat(i, " cell:", cell, " treat:", treat, "\n")
   ## fn <- paste("../dispersionQTL/dispersion_egenes/5PCs_", cell, "_", treat, "_genes.txt", sep="")
   ## res <- read.table(fn)
   fn <- paste("../dispersionQTL/disp-eQTL_output/", cell, "_", treat, ".GEPC5.nominals.eQTL.txt.gz", sep="")
   res <- fread(fn, data.table=F, stringsAsFactors=F)%>%
       mutate(qvalue=qvalue(V4)$qvalues, gene_SNP=paste(V1, V2, sep="_"))
   res2 <- res%>%drop_na(V4)%>%dplyr::filter(qvalue<0.1)
   ## gene_SNP <- res2$gene_SNP
   ## gene_SNP
   res2$gene_SNP
})

QTLs <- unique(unlist(QTLs))

### read all LDA-interaction results
plotdata <- map_dfr(1:nrow(dataset), function(i){
###
   cell <- dataset[i,1]
   lda <- dataset[i,2]
   treat <- dataset[i,3]
   cat(i, " ", lda, " ", treat, "\n")  
   fn <- paste0(option, "/4_lm.results/2_", cell, "_lda", lda, "_trt", treat, "_stratified.FDR.txt")
   res <- fread(fn, data.table=F)
   res <- res%>%mutate(is_QTL=ifelse(gene_SNP%in%QTLs, "yes", "no"))

   res2 <- res%>%group_by(is_QTL)%>%nest()%>%
      mutate(outlist=map(data,~qq.fun(.x)))%>%
      dplyr::select(-data)%>%
      unnest(outlist)%>%as.data.frame()%>%
      mutate(MCls=cell, LDA=lda, treats=treat)       
   res2
})


### QQ plots
plotdata2 <- plotdata%>%mutate(treat2=gsub("-EtOH", "", treats))

lab.treats <- c("LPS"="LPS", "LPS-DEX"="LPS+DEX",
                "PHA"="PHA", "PHA-DEX"="PHA+DEX")
p <- ggplot(plotdata2)+
   rasterize(geom_point(aes(x=expected, y=observed, colour=factor(is_QTL)),size=0.5),dpi=300)+
   scale_colour_manual("", values=c("yes"="green", "no"="grey30"),
      labels=c("yes"="vGene", "no"="not vGene"),
      guide=guide_legend(override.aes=list(size=2)))+
  geom_abline(colour="red")+
  labs(x=bquote("Expected"~ -log[10]~"("~plain(P)~")"),
       y=bquote("Observed"~-log[10]~"("~plain(P)~")"))+
  facet_grid(MCls~treat2, labeller=labeller(treat2=lab.treats))+
  theme_bw()

figfn <- paste(option, "/summary/Figure2.4_qq.anova_vQTL.png", sep="")
png(figfn, width=750, height=750, res=120)
print(p)
dev.off()

figfn <- paste(option, "/summary/Figure2.4_qq.anova_vQTL.pdf", sep="")
pdf(figfn, width=7.5, height=7.5)
print(p)
dev.off()


####
#### test if enriched in eGene (mean expression)
### read mQTLs
cells <- rep(c("Bcell", "Monocyte", "NKcell", "Tcell"), each=5)
treats <- rep(c("CTRL", "LPS-EtOH", "LPS-DEX", "PHA-EtOH", "PHA-DEX"), times=4)
xx <- data.frame(MCls=cells, treats=treats)
QTLs <- map(1:nrow(xx), function(i){
   cell <- xx[i,1]
   treat <- xx[i,2]
   cat(i, " cell:", cell, " treat:", treat, "\n")
   ## fn <- paste("../dispersionQTL/dispersion_egenes/5PCs_", cell, "_", treat, "_genes.txt", sep="")
   ## res <- read.table(fn)
   fn <- paste("../eQTL-on-mean/mean-eQTL_output/", cell, "_", treat, ".GEPC7.nominals.eQTL.txt.gz", sep="")
   res <- fread(fn, data.table=F, stringsAsFactors=F)%>%
       mutate(qvalue=qvalue(V4)$qvalues, gene_SNP=paste(V1, V2, sep="_"))
   res2 <- res%>%drop_na(V4)%>%dplyr::filter(qvalue<0.1)
   ## gene_SNP <- res2$gene_SNP
   ## gene_SNP
   res2$gene_SNP
})

QTLs <- unique(unlist(QTLs))

### read all LDA-interaction results
plotdata <- map_dfr(1:nrow(dataset), function(i){
###
   cell <- dataset[i,1]
   lda <- dataset[i,2]
   treat <- dataset[i,3]
   cat(i, " ", lda, " ", treat, "\n")  
   fn <- paste0(option, "/4_lm.results/2_", cell, "_lda", lda, "_trt", treat, "_stratified.FDR.txt")
   res <- fread(fn, data.table=F)
   res <- res%>%mutate(is_meGene=ifelse(gene_SNP%in%QTLs, "yes", "no"))

   res2 <- res%>%group_by(is_meGene)%>%nest()%>%
      mutate(outlist=map(data,~qq.fun(.x)))%>%
      dplyr::select(-data)%>%
      unnest(outlist)%>%as.data.frame()%>%
      mutate(MCls=cell, LDA=lda, treats=treat)       
   res2
})


### QQ plots
plotdata2 <- plotdata%>%mutate(treat2=gsub("-EtOH", "", treats))

lab.treats <- c("LPS"="LPS", "LPS-DEX"="LPS+DEX",
                "PHA"="PHA", "PHA-DEX"="PHA+DEX")
p <- ggplot(plotdata2)+
   rasterize(geom_point(aes(x=expected, y=observed, colour=factor(is_meGene)),size=0.5),dpi=300)+
   scale_colour_manual("", values=c("yes"="green", "no"="grey30"),
      labels=c("yes"="eGenes", "no"="not eGenes"),
      guide=guide_legend(override.aes=list(size=2)))+
  geom_abline(colour="red")+
  labs(x=bquote("Expected"~ -log[10]~"("~plain(P)~")"),
       y=bquote("Observed"~-log[10]~"("~plain(P)~")"))+
  facet_grid(MCls~treat2, labeller=labeller(treat2=lab.treats))+
  theme_bw()

figfn <- paste(option, "/summary/Figure2.5_qq.anova_meGene.png", sep="")
png(figfn, width=750, height=750, res=120)
print(p)
dev.off()

figfn <- paste(option, "/summary/Figure2.5_qq.anova_meGene.pdf", sep="")
pdf(figfn, width=7.5, height=7.5)
print(p)
dev.off()



###
### eQTL 
cells <- rep(c("Bcell", "Monocyte", "NKcell", "Tcell"), each=5)
treats <- rep(c("CTRL", "LPS-EtOH", "LPS-DEX", "PHA-EtOH", "PHA-DEX"), times=4)
xx <- data.frame(MCls=cells, treats=treats)
gene_SNP <- map(1:nrow(xx), function(i){
   cell <- xx[i,1]
   treat <- xx[i,2]
   cat(i, " cell:", cell, " treat:", treat, "\n")
   ## fn <- paste("../dispersionQTL/dispersion_egenes/5PCs_", cell, "_", treat, "_genes.txt", sep="")
   ## res <- read.table(fn)
   fn <- paste("../eQTL-on-mean/mean-eQTL_output/", cell, "_", treat, ".GEPC7.nominals.eQTL.txt.gz", sep="")
   res <- fread(fn, data.table=F, stringsAsFactors=F)%>%
       mutate(qvalue=qvalue(V4)$qvalues, gene_SNP=paste(V1, V2, sep="_"))
   res2 <- res%>%drop_na(V4)%>%dplyr::filter(qvalue<0.1)
   gene_SNP <- res2$gene_SNP
   gene_SNP
   ## gene <- res2$V1
})

gene_SNP2 <- unique(unlist(gene_SNP))

### read all LDA-interaction results
plotdata <- map_dfr(1:nrow(dataset), function(i){
###
   cell <- dataset[i,1]
   lda <- dataset[i,2]
   treat <- dataset[i,3]
   cat(i, " ", lda, " ", treat, "\n")  
   fn <- paste0(option, "/4_lm.results/2_", cell, "_lda", lda, "_trt", treat, "_stratified.FDR.txt")
   res <- fread(fn, data.table=F)
   res <- res%>%mutate(is_meQTL=ifelse(gene_SNP%in%gene_SNP2, "yes", "no"))

   res2 <- res%>%group_by(is_meQTL)%>%nest()%>%
      mutate(outlist=map(data,~qq.fun(.x)))%>%
      dplyr::select(-data)%>%
      unnest(outlist)%>%as.data.frame()%>%
      mutate(MCls=cell, LDA=lda, treats=treat)       
   res2
})


### QQ plots
plotdata2 <- plotdata%>%mutate(treat2=gsub("-EtOH", "", treats))

lab.treats <- c("LPS"="LPS", "LPS-DEX"="LPS+DEX",
                "PHA"="PHA", "PHA-DEX"="PHA+DEX")
p <- ggplot(plotdata2)+
   rasterize(geom_point(aes(x=expected, y=observed, colour=factor(is_meQTL)),size=0.5),dpi=300)+
   scale_colour_manual("", values=c("yes"="green", "no"="grey30"),
      labels=c("yes"="eQTL", "no"="not eQTL"),
      guide=guide_legend(override.aes=list(size=2)))+
  geom_abline(colour="red")+
  labs(x=bquote("Expected"~ -log[10]~"("~plain(P)~")"),
       y=bquote("Observed"~-log[10]~"("~plain(P)~")"))+
  facet_grid(MCls~treat2, labeller=labeller(treat2=lab.treats))+
  theme_bw()

figfn <- paste(option, "/summary/Figure2.5.2_qq.anova_meQTL.png", sep="")
png(figfn, width=750, height=750, res=120)
print(p)
dev.off()
