##
library(tidyverse)
library(data.table)
library(parallel)
library(qvalue)
library(qqman)
library(KEGGREST)
library(clusterProfiler)
library(org.Hs.eg.db)
library(cowplot)
theme_set(theme_grey())


###############################################
### 1, summary LDA-interacting eQTL results ###
###############################################

rm(list=ls())

dataset <- read.table("dataset_contrast.txt",header=F)
names(dataset) <- c("MCls", "LDA", "treat")

option <- "DiagLDA2"

###
### summary LDA-interaction eQTL results
FDR <- 0.1
summ <- map_dfr(1:nrow(dataset), function(i){
   x <- dataset[i,]
   cell <- x[1]
   lda <- x[2]
   treat <- x[3]
   cat(paste(cell, lda, treat, sep="_"), "\n")
   fn <- paste0(option, "/4.2_lm.results/2_", cell, "_lda", lda, "_trt", treat, "_stratified.FDR.txt")
   res <- fread(fn, data.table=F)%>%
     drop_na(interaction_ANOVA.qvalue)
   res1 <- res%>%filter(interaction_ANOVA.qvalue<FDR)
   res2 <- res%>%filter(interaction_ANOVA.qvalue_stratified<FDR)%>%
     dplyr::select(ENSG, varID)%>%
     mutate(gene_SNP=paste(ENSG, "_", varID, sep=""))  
   df <- data.frame(MCls=cell, LDA=lda, treat=treat,
      signif=nrow(res1), egene=length(unique(res1$ENSG)),
      signif_stratified=nrow(res2), egene_stratified=length(unique(res2$ENSG)))
   df
})

opfn <- paste(option, "/summary/Table1_summary.ieGene.csv", sep="")
write.csv(summ, opfn, row.names=FALSE)


###
###
FDR <- 0.1
## prefix <- "/wsu/home/groups/piquelab/SCAIP/SCAIP1-6_protein-coding/LDA-eQTL_Julong/"
res <- map_dfr(1:nrow(dataset), function(i){
   x <- dataset[i,]
   cell <- x[1]
   lda <- x[2]
   treat <- x[3]
   cat(paste(cell, lda, treat, sep="_"), "\n")
   # new
   fn <- paste0(option, "/4.2_lm.results/2_", cell, "_lda", lda, "_trt", treat, "_stratified.FDR.txt")
   # old
   ## fn <- paste0(prefix, option, "/lm_results/2_", cell, "_lda", lda, "_trt", treat, "_stratified.FDR.txt")
   res <- fread(fn, data.table=F)%>%
     drop_na(interaction_ANOVA.qvalue)
   ## res1 <- res%>%filter(interaction_ANOVA.qvalue<FDR)
   res2 <- res%>%filter(interaction_ANOVA.qvalue_stratified<FDR)%>%
     ## dplyr::select(ENSG, varID)%>%
     mutate(gene_SNP=paste(ENSG, "_", varID, sep=""), condition=paste(cell, lda, treat, sep="_"))  
   res2
})

res <- res%>%as.data.frame()
#new
opfn <- paste(option, "/4.2_lm.results/zzz_dynamical.eQTL.rds", sep="")
write_rds(res, opfn)
### old
## opfn <- paste(prefix, option, "/lm_results/zzz_dynamical.eQTL.rds", sep="")
## write_rds(res, opfn)


## ### total old, LDA-ieGene and iQTL
## summ <- map_dfr(1:nrow(dataset), function(i){
##    x <- dataset[i,]
##    cell <- x[1]
##    lda <- x[2]
##    treat <- x[3]
##    cat(paste(cell, lda, treat, sep="_"), "\n")
##    prefix <- "/wsu/home/groups/piquelab/SCAIP/SCAIP1-6_protein-coding/LDA-eQTL_Julong/DiagLDA2/" 
##    fn <- paste(prefix, "lm_results/2_", cell, "_lda", lda, "_trt", treat, "_stratified.FDR.txt", sep="")
##    res <- fread(fn, data.table=F, stringsAsFactors=F)%>%
##       drop_na(interaction_ANOVA.qvalue)
##    res2 <- res%>%dplyr::filter(interaction_ANOVA.qvalue_stratified<FDR)%>%
##      dplyr::select(ENSG, varID)%>%   
##      mutate(gene_SNP=paste(ENSG, "_", varID, sep=""))
##    res2
## })
## length(unique(summ$ENSG))
## length(unique(summ$varID))


###
### compare with eGene, response eGenes, dGene and DEG
rm(list=ls())

option <- "DiagLDA2"
contrast_ls <- list("LPS"=c("CTRL", "LPS-EtOH"),
                    "LPS-DEX"=c("LPS-EtOH", "LPS-DEX"),
                    "PHA"=c("CTRL", "PHA-EtOH"),
                    "PHA-DEX"=c("PHA-EtOH", "PHA-DEX"))
dataset <- read.table("dataset_contrast.txt", header=F)
names(dataset) <- c("MCls", "LDA", "treats")
###
compared  <- mclapply(1:nrow(dataset), function(i){
###
   x <- dataset[i,]
   cell <- as.character(x[1])
   lda <- as.character(x[2])
   treat <- as.character(x[3])
   
   cat(paste(cell, lda, treat, sep="_"), "\n")

   ###New LDA-i eGenes
   fn <- paste0(option, "/4_lm.results/2_", cell, "_lda", lda, "_trt", treat, "_stratified.FDR.txt") 
   res <- fread(fn, data.table=F)
   gene1 <- res%>%
       drop_na(interaction_ANOVA.qvalue_stratified)%>%
       filter(interaction_ANOVA.qvalue_stratified<0.1)%>%
       dplyr::pull(ENSG)
   gene1 <- unique(gene1)

   ###
   ### response
   oneX <- contrast_ls[[lda]]
   treat1 <- oneX[2]
   treat0 <- oneX[1]
   ## fn <- paste("../reQTL_entire-cis/reQTL_lm_results/results_", cell, "_", treat1, "_vs_", treat0, "_stratified_FDR.txt", sep="")
   ## res <- fread(fn, data.table=F, stringsAsFactors=F)
   ## gene2 <- res%>%
   ##   drop_na(interaction_ANOVA.qvalue_stratified)%>%
   ##   filter(interaction_ANOVA.qvalue_stratified<0.1)%>%dplyr::pull(ENSG)
   ## gene2 <- unique(gene2)

   ###
   ### mash response
   fn <- "../mashr_eQTL/mashr-reQTLs_union-unshared-magnitude2-mlfsr0.1.Rd"
   load(fn)
   uum_name <- paste(cell, "_", treat1, ":", cell, "_", treat0, sep="")
   gene3 <- unique(gsub("_.*", "", uum[[uum_name]]))

   ###
   ### eGenes
   fn <- paste("/wsu/home/groups/piquelab/SCAIP/SCAIP-genetic/eQTL/egenes/", cell, "_", treat, "_egenes.txt", sep="")
   gene4 <- read.table(fn)%>%dplyr::pull(V1)

   ### 
   ### dGenes 
   fn <- paste("/wsu/home/groups/piquelab/SCAIP/SCAIP1-6_protein-coding/dispersionQTL/disp-eQTL_output/",
      cell, "_", treat, ".GEPC3.nominals.eQTL.txt.gz", sep="")
   res <- fread(fn, data.table=F, col.names=c("ENSG", "varID", "distance", "pvalue", "estimate"))%>%
     mutate(qvalue=qvalue(pvalue)$qvalues, gene_SNP=paste(ENSG, varID, sep="_"))
   gene5 <- res%>%drop_na(qvalue)%>%filter(qvalue<0.1)%>%
      dplyr::pull(ENSG)
   
   ###DEG
   res <- read_rds("/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/6_DEG.CelltypeNew_output/Filter2/2_meta.rds")
  gene6 <- res%>%filter(abs(beta)>0.5, qval<0.1)%>%
     mutate(MCls==cell, contrast==lda)%>%dplyr::pull(gene)
  gene1.1 <- gsub("\\..*", "", gene1) 
    
   tmp <- data.frame(MCls=cell, contrast=lda,
     ## shared_reGene=length(intersect(gene1, gene2)),
     "shared_reGene-mash"=length(intersect(gene1, gene3)),
     shared_eGene=length(intersect(gene1, gene4)),
     "shared_dGene"=length(intersect(gene1, gene5)),
     shared_DEG=length(intersect(gene1.1, gene6)) )
 tmp
}, mc.cores=1)

compared <- do.call(rbind, compared)

opfn <- paste(option, "/summary/table2_summary.overlap.csv", sep="")
write.csv(compared, opfn, row.names=FALSE)



###
###

option <- "DiagLDA2"
contrast_ls <- list("LPS"=c("CTRL", "LPS-EtOH"),
                    "LPS-DEX"=c("LPS-EtOH", "LPS-DEX"),
                    "PHA"=c("CTRL", "PHA-EtOH"),
                    "PHA-DEX"=c("PHA-EtOH", "PHA-DEX"))
dataset <- read.table("dataset_contrast.txt", header=F)
names(dataset) <- c("MCls", "LDA", "treats")
x
fn <- paste(option, "/4.2_lm.results/zzz_dynamical.eQTL.rds", sep="")
lda <- read_rds(fn)

###
load("../mashr_eQTL/mashr-eQTLs-mlfsr0.1.Rd")
mash_eqtl <- map_dfr(names(eqtl), function(ii){
   df <- data.frame(eqtl=eqtl[[ii]], ENSG2=gsub("\\..{1,2}_.*", "", eqtl[[ii]]), condition=ii)
   df
})

###
load("../mashr_dispersionQTL/mashr-eQTLs-mlfsr0.1.Rd")
mash_dqtl <- map_dfr(names(eqtl), function(ii){
   df <- data.frame(eqtl=eqtl[[ii]], ENSG2=gsub("\\..{1,2}_.*", "", eqtl[[ii]]), condition=ii)
   df
})

res <- map_dfr(1:16, function(i){
  ##
  cell <- dataset[i,1]
  LDA <- dataset[i,2]
  treat <- dataset[i,3]
  ###
  condition2 <- paste(cell, "_", LDA, "_", treat, sep="")
  condition3 <- paste(cell, "_", treat, sep="")
  ##
  lda2 <- lda%>%filter(condition==condition2)%>%
      dplyr::select(ENSG, varID, condition)%>%
      mutate(ENSG2=gsub("\\..*", "", ENSG))
  ## ##
  ## egene <- mash_eqtl%>%filter(condition==condition3)%>%dplyr::pull(ENSG2)
  ## egene <- unique(egene)
   ### eGenes
   fn <- paste("/wsu/home/groups/piquelab/SCAIP/SCAIP-genetic/eQTL/egenes/", condition3,  "_egenes.txt", sep="")
   egene <- read.table(fn)%>%dplyr::pull(V1)
   egene <- unique(gsub("\\..*", "", egene)) 
    
  ##
  ## dgene <- mash_dqtl%>%filter(condition==condition3)%>%dplyr::pull(ENSG2)
  fn <- paste("/wsu/home/groups/piquelab/SCAIP/SCAIP-genetic/dispersionQTL/disp-eQTL_output/",
        condition3,  ".GEPC5.nominals.eQTL.txt.gz", sep="")
   resD <- fread(fn, data.table=F, col.names=c("ENSG", "varID", "distance", "pvalue", "estimate"))%>%
     mutate(qvalue=qvalue(pvalue)$qvalues, gene_SNP=paste(ENSG, varID, sep="_"),
            ENSG2=gsub("\\..*", "", ENSG))
   dgene <- resD%>%drop_na(qvalue)%>%filter(qvalue<0.1)%>%
      dplyr::pull(ENSG)
   dgene <- unique(dgene) 
    
  ## dgene <- read.table(fn)%>%dplyr::pull(V1)  
  ## dgene <- unique(dgene)
  ##
  lda2 <- lda2%>%mutate(is_egene=ifelse(ENSG2%in%egene, 1, 0),
                        is_dgene=ifelse(ENSG2%in%dgene, 1,0))
  lda2
})    
    
    
    
##
x1 <- res%>%filter(is_egene==1)
length(unique(x1$ENSG2))

##
x2 <- res%>%filter(is_dgene==1)
length(unique(x2$ENSG2))



##################################################
### compare the corrected and previous results ###
##################################################

## option <- "DiagLDA2"
## contrast_ls <- list("LPS"=c("CTRL", "LPS-EtOH"),
##                     "LPS-DEX"=c("LPS-EtOH", "LPS-DEX"),
##                     "PHA"=c("CTRL", "PHA-EtOH"),
##                     "PHA-DEX"=c("PHA-EtOH", "PHA-DEX"))
## dataset <- read.table("dataset_contrast.txt", header=F)
## names(dataset) <- c("MCls", "LDA", "treats")

## ###
## compared  <- map_dfr(1:nrow(dataset), function(i){
## ###
##    x <- dataset[i,]
##    cell <- as.character(x[1])
##    lda <- as.character(x[2])
##    treat <- as.character(x[3])
   
##    cat(paste(cell, lda, treat, sep="_"), "\n")

##    ###New LDA-i eGenes
##    fn <- paste0(option, "/4_lm.results/2_", cell, "_lda", lda, "_trt", treat, "_stratified.FDR.txt") 
##    new <- fread(fn, data.table=F)
##    new2 <- new%>%
##        drop_na(interaction_ANOVA.qvalue_stratified)%>%
##        dplyr::filter(interaction_ANOVA.qvalue_stratified<0.1)%>%
##        dplyr::select(ENSG, varID)%>%
##        mutate(gene_SNP=paste(ENSG, "_", varID, sep="")) 
    
##    geneNew <- unique(new2$ENSG) 
    
##    ###
##    ### old
##    prefix <- "/wsu/home/groups/piquelab/SCAIP/SCAIP1-6_protein-coding/LDA-eQTL_Julong/DiagLDA2/" 
##    fn <- paste(prefix, "lm_results/2_", cell, "_lda", lda, "_trt", treat, "_stratified.FDR.txt", sep="")
##    old <- fread(fn, data.table=F, stringsAsFactors=F)
##    old2 <- old%>%
##      drop_na(interaction_ANOVA.qvalue_stratified)%>%
##      dplyr::filter(interaction_ANOVA.qvalue_stratified<0.1)%>%
##      dplyr::select(ENSG, varID)%>%
##      mutate(gene_SNP=paste(ENSG, "_", varID, sep="")) 

##    geneOld <- unique(old2$ENSG)

##    ###
##    ### 
##    tmp <- data.frame(MCls=cell, contrast=lda,
##       new_ieQTL=length(unique(new2$gene_SNP)), new_ieGene=length(geneNew),
##       old_ieQTL=length(unique(old2$gene_SNP)), old_ieGene=length(geneOld),
##       shared=length(intersect(geneNew, geneOld)) )
##    tmp
## })

## opfn <- paste(option, "/summary/table3_compare.New-old.csv", sep="")
## write.csv(compared, opfn, row.names=FALSE)



############################################################
### 3. test if LDA-eGenes are enriched in eGenes and DEG ###
############################################################

rm(list=ls())

option <- "DiagLDA2"
contrast_ls <- list("LPS"=c("CTRL", "LPS-EtOH"),
                    "LPS-DEX"=c("LPS-EtOH", "LPS-DEX"),
                    "PHA"=c("CTRL", "PHA-EtOH"),
                    "PHA-DEX"=c("PHA-EtOH", "PHA-DEX"))
dataset <- read.table("dataset_contrast.txt", header=F)
names(dataset) <- c("MCls", "LDA", "treats")

generateData <- function(celltype){
###
   ###
   dataset2 <- dataset%>%dplyr::filter(MCls==celltype)
   nr <- nrow(dataset2)

   ### mashr-eGenes
   load("../mashr_eQTL/mashr-eQTLs-mlfsr0.1.Rd")
   mash_eqtl <- unlist(eqtl)
   ### mashr-reGenes
   load("../mashr_eQTL/reQTLs_mashr_magnitude2-mlfsr0.1.Rd") 
   mash_reqtl <- unlist(uum) 

    
## Generate data for enrichment analysis of eGene, DEG, TWAS genes
   dfcomb<- map_dfr(1:nr, function(i){
###    
      x <- dataset2[i,]
      cell <- as.character(x[1])
      lda <- as.character(x[2])
      treat <- as.character(x[3])
   
      cat(paste(cell, lda, treat, sep="_"), "\n")

   ###New LDA-i eGenes
      fn <- paste0(option, "/4.2_lm.results/2_", cell, "_lda", lda, "_trt", treat, "_stratified.FDR.txt") 
      res <- fread(fn, data.table=F)
      res <- res%>%mutate(ENSG2=gsub("\\..*", "", ENSG)) 
      gene1 <- res%>%
         drop_na(interaction_ANOVA.qvalue_stratified)%>%
         filter(interaction_ANOVA.qvalue_stratified<0.1)%>%
         dplyr::pull(ENSG2)
      gene1 <- unique(gene1)
     
   ### fast eGenes
      fn <- paste("../eQTL/egenes/", cell, "_", treat, "_egenes.txt", sep="")
      eqtl <- read.table(fn)
      geneCategory <- unique(gsub("\\..*", "", eqtl$V1))
       
      ## fn <- paste("../eQTL/eQTL_output/", cell, "_", treat, ".GEPC4.nominals.eQTL.txt.gz", sep="") 
      ## TMP <- fread(fn, data.table=F, sep=" ")%>%mutate(ENSG2=gsub("\\..*", "", V1))       
      ## geneBg <- unique(intersect(res$ENSG2, TMP$ENSG2))
       
      df <- data.frame(LDA_eGene=length(gene1), 
                       ## LDA_eGene2=length(intersect(gene1, geneBg)),
                       overlap=length(intersect(gene1, geneCategory)),
                       interested.gene=length(unique(intersect(geneCategory, res$ENSG2))),
                       Bg.gene=length(unique(res$ENSG2)), category="eGene")

   ### mash eGenes
       condition <- paste(cell, treat, sep="_")
       geneCategory <- unique(gsub("\\..{1,2}_.*", "", mash_eqtl))
       df2 <- data.frame(LDA_eGene=length(gene1),
                         overlap=length(intersect(gene1, geneCategory)),
                         interested.gene=length(unique(intersect(geneCategory, res$ENSG2))),
                         Bg.gene=length(unique(res$ENSG2)), category="eGene-mash")
                         
   ### response eGenes
      contrast <- contrast_ls[[lda]] 
      treat1 <- paste(cell, contrast[2], sep="_")
      treat0 <- paste(cell, contrast[1], sep="_")
      condition <- paste(treat1, treat0, sep=":")
      geneCategory <- unique(gsub("\\..{1,2}_.*", "", mash_reqtl))
      df3 <- data.frame(LDA_eGene=length(gene1),
                         overlap=length(intersect(gene1, geneCategory)),
                         interested.gene=length(unique(intersect(geneCategory, res$ENSG2))),
                         Bg.gene=length(unique(res$ENSG2)), category="reGene-mash")
       
      df <- rbind(df, df2, df3)

      df$condition <- paste(cell, lda, treat, sep="_")
      df 
   })##
   
   dfcomb

 }
       
       
       
   ### DEGs
   ## res_DEG <- read_rds("/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/6_DEG.CelltypeNew_output/Filter2/2_meta.rds")
   ## res_DEG <- res_DEG%>%filter(MCls==cell, contrast==lda)
   ## geneBg <- unique(intersect(gsub("\\..*", "", res$ENSG), res_DEG$gene)) 
   ## gene3 <- res_DEG%>%filter(abs(beta)>0.5, qval<0.1)%>%dplyr::pull(gene)   
   ## gene1.1 <- gsub("\\..*", "", gene1)
   ## df0 <- data.frame(LDA_eGene=length(gene1),
   ##    LDA_eGene2=length(intersect(gene1.1, geneBg)),
   ##    LDA_eGene3.overlap=length(intersect(gene1.1, gene3)),
   ##    interested.gene=length(intersect(gene3, geneBg)),
   ##    Bg.gene=length(geneBg), category="DEG")
   ## df <- rbind(df, df0)
   ## df$condition <- paste(cell, lda, treat, sep="_")
   ## df
   ###
   ### response eGenes
    
## })
    
## opfn <- paste(option, "/summary/Table4_summary.Tcells.enriched.csv", sep="")
## write.csv(LDA_eGene, opfn, row.names=FALSE)



### Fisher.test
my.fisher <- function(df){
###    
   resFish <- map_dfr(1:nrow(df), function(i){   
      LDA.genes <- c(df[i,2], df[i,1]-df[i,2])
      Bg.genes <- c(df[i,3], df[i,4]-df[i,3])         
      dat <- data.frame(LDA.genes=LDA.genes, Bg.genes=Bg.genes)
      rownames(dat) <- c("in.category", "not.category")
      fish <- fisher.test(dat)
      res0 <- data.frame(odds=log2(as.numeric(fish$estimate)),
         CI.low=log2(fish$conf.int[1]), CI.high=log2(fish$conf.int[2]),
         pvalue=fish$p.value)
      res0
   })
   df <- cbind(df, resFish)
   df
}

###
dfcomb <- generateData("Tcell")
dfenriched <- my.fisher(dfcomb)

dfenriched <-  dfenriched%>%mutate(treats=gsub(".*_", "", condition),
      treat2=gsub("-EtOH", "", treats),
      is_sig=ifelse(pvalue<0.05, 1, 2))

write.csv(dfenriched, "./DiagLDA2/summary/Table4_summary.Tcells.enriched.csv", row.names=F)

###
col2 <- c("LPS"="#fb9a99", "LPS-DEX"="#e31a1c", "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")

###
### Figures of eGene enrichment
df2 <- dfenriched%>%dplyr::filter(category=="eGene")
p <- ggplot(df2, aes(x=odds, y=condition))+
   geom_vline(aes(xintercept=0), size=0.25, linetype="dashed")+
   geom_errorbarh(aes(xmax=CI.high, xmin=CI.low, colour=treat2),
                  size=0.5, height=0.2)+
   geom_point(aes(colour=treat2), shape=19, size=2.5)+
   scale_colour_manual(values=col2)+
   scale_y_discrete(labels=c("Tcell_LPS_LPS-EtOH"="LPS",
      "Tcell_LPS-DEX_LPS-DEX"="LPS-DEX",
      "Tcell_PHA_PHA-EtOH"="PHA",
      "Tcell_PHA-DEX_PHA-DEX"="PHA-DEX"))+
   ggtitle("T cells")+ 
   annotate("text", x=2.5, y=4.3,
      label=bquote(~"Fisher's Exact Test"~ " "~italic(p)~"<0.001"),
      size=3)+ 
   xlab(bquote(~log[2]~"(Odds ratio)"))+
   theme_bw()+
   theme(plot.title=element_text(hjust=0.5, size=12),
         axis.title.y=element_blank(),
         legend.position="none")

figfn <- "./DiagLDA2/summary/Table4_Figure1.enrich.eGene.png"
png(figfn, width=420, height=480, res=120)
print(p)
dev.off()


###
### figure enrichment analysis for eGene from mash

dfenriched <- read.csv("./DiagLDA2/summary/Table4_summary.Tcells.enriched.csv")

df2 <- dfenriched%>%dplyr::filter(category=="eGene-mash")%>%
    mutate(symbol=ifelse(pvalue<0.001, "***", "ns"))
###
p2 <- ggplot(df2, aes(x=odds, y=condition))+
   geom_vline(aes(xintercept=0), size=0.25, linetype="dashed")+
   geom_errorbarh(aes(xmax=CI.high, xmin=CI.low, colour=treat2),
                  size=0.5, height=0.2)+
   geom_point(aes(colour=treat2), shape=19, size=2.5)+
   scale_colour_manual(values=col2)+
   scale_y_discrete(labels=c("Tcell_LPS_LPS-EtOH"="LPS",
      "Tcell_LPS-DEX_LPS-DEX"="LPS+DEX",
      "Tcell_PHA_PHA-EtOH"="PHA",
      "Tcell_PHA-DEX_PHA-DEX"="PHA+DEX"))+
   ggtitle("eGenes (T cells)")+ 
   ## annotate("text",x=1.5, y=4.3,
   ##    label=bquote(atop("Fisher's Exact Test"~ " "~italic(p)~"< 0.001",
   ##                   "(except LPS-DEX ns)")),
   ##    size=3)+
   geom_text(aes(y=condition, x=CI.high+0.1, label=symbol), angle=-90)+ 
   xlab(bquote(~log[2]~"(Odds ratio)"))+
   theme_bw()+
   theme(plot.title=element_text(hjust=0.5, size=12),
         axis.title.y=element_blank(),
         legend.position="none")

figfn <- "./DiagLDA2/summary/Table4_Figure2.enrich.eGene-mash.png"
png(figfn, width=420, height=480, res=120)
print(p)
dev.off()


###
### figure of enrichment analysis for reGene
df3 <- dfenriched%>%dplyr::filter(category=="reGene-mash")%>%
   mutate(symbol=ifelse(pvalue<0.001, "***", "ns")) 
p3 <- ggplot(df3, aes(x=odds, y=condition))+
   geom_vline(aes(xintercept=0), size=0.25, linetype="dashed")+
   geom_errorbarh(aes(xmax=CI.high, xmin=CI.low, colour=treat2),
                  size=0.5, height=0.2)+
   geom_point(aes(colour=treat2), shape=19, size=2.5)+
   scale_colour_manual(values=col2)+
   scale_y_discrete(labels=c("Tcell_LPS_LPS-EtOH"="LPS",
      "Tcell_LPS-DEX_LPS-DEX"="LPS+DEX",
      "Tcell_PHA_PHA-EtOH"="PHA",
      "Tcell_PHA-DEX_PHA-DEX"="PHA+DEX"))+
   ggtitle("reGenes (T cells)")+ 
   ## annotate("text",x=1.5, y=4.3,
   ##    label=bquote(atop("Fisher's Exact Test"~ " "~italic(p)~"<0.001",
   ##                   "(except LPS-DEX ns)")),
   ##    size=3)+
   geom_text(aes(y=condition, x=CI.high+0.2, label=symbol), angle=-90)+ 
   xlab(bquote(~log[2]~"(Odds ratio)"))+
   theme_bw()+
   theme(plot.title=element_text(hjust=0.5, size=12),
         axis.title.y=element_blank(),
         legend.position="none")

figfn <- "./DiagLDA2/summary/Table4_Figure3.enrich.reGene-mash.png"
png(figfn, width=420, height=480, res=120)
print(p)
dev.off()

##
figfn <- "./DiagLDA2/summary/Table4_Figure4_enrich_eGene-reGene.pdf"
pdf(figfn, width=7, height=4)
plot_grid(p2, p3, ncol=2,
    labels="AUTO", label_size=14, label_x=0.1, label_fontface="plain", align="hv", axis="lr")
dev.off()
###








#####################
### 4. TWAS genes ###
#####################

option <- "DiagLDA2"
dataset <- read.table("dataset_contrast.txt", header=F)
names(dataset) <- c("MCls", "LDA", "treats")

###asthma genes
path <- keggGet("hsa05310") ##asthma
genes <- path[[1]]$GENE
entrez <- genes[(seq_len(length(genes))%%2)==1]
asthma <- bitr(entrez, fromType="ENTREZID", toType=c("ENSEMBL","SYMBOL"),
              OrgDb="org.Hs.eg.db") 
gene_asthma <- asthma$ENSEMBL

### twas genes
immue.traits <- read.table("../TWAS_overlap/PTWAS-immune-traits.txt")
TWAS <- read.table("../TWAS_overlap/media-4_bioRxiv_TableS2.txt", header=T)
TWAS2 <- TWAS%>%filter(grepl("Asthma|asthma", Trait))

##
twas0 <- gsub("\\..*", "", unique(TWAS$Gene))
twas2 <- gsub("\\..*", "",  unique(TWAS2$Gene))

##
fn <- paste(option, "/4.2_lm.results/zzz_dynamical.eQTL.rds", sep="")
res <- read_rds(fn)

## prefix <- "/wsu/home/groups/piquelab/SCAIP/SCAIP1-6_protein-coding/LDA-eQTL_Julong/"
## fn <- paste(prefix, option, "/lm_results/zzz_dynamical.eQTL.rds", sep="")
## res <- read_rds(fn)

LDA_eGene <- map_dfr(1:16, function(i){
###    
   x <- dataset[i,]
   cell <- as.character(x[1])
   lda <- as.character(x[2])
   treat <- as.character(x[3])
   ii <- paste(cell, lda, treat, sep="_")    
   cat(paste(cell, lda, treat, sep="_"), "\n")

   res2 <- res%>%dplyr::filter(condition==ii)  

###New LDA-i eGenes
   ## fn <- paste0(option, "/4.2_lm.results/2_", cell, "_lda", lda, "_trt", treat, "_stratified.FDR.txt") 
   ## res <- fread(fn, data.table=F)
   ## gene1 <- res%>%
   ##     drop_na(interaction_ANOVA.qvalue_stratified)%>%
   ##     filter(interaction_ANOVA.qvalue_stratified<0.1)%>%
   ##     dplyr::pull(ENSG)
    
   ENSG <- gsub("\\..*", "",  unique(res2$ENSG))
    
   ###   
   df <- data.frame(eGene=length(ENSG),
              overlap.twasImmune=length(intersect(ENSG, twas2)),
              overlap.twasAll=length(intersect(ENSG,twas0)),
              ## overlap.asthma=length(intersect(ENSG, gene_asthma)),
              condition=ii)
   df     
})
##
LDA_eGene <- as.data.frame(LDA_eGene)

##
opfn <- paste(option, "/summary/table2.1_summary.overlap.twas.csv", sep="")
write.csv(LDA_eGene, opfn, row.names=F)

###
### summary total twas gene in Tcells
fn <- paste(option, "/4.2_lm.results/zzz_dynamical.eQTL.rds", sep="")
res <- read_rds(fn)

LDA_eGene <- map_dfr(13:16, function(i){
###    
   x <- dataset[i,]
   cell <- as.character(x[1])
   lda <- as.character(x[2])
   treat <- as.character(x[3])
   ii <- paste(cell, lda, treat, sep="_")   
   cat(paste(cell, lda, treat, sep="_"), "\n")

   ###New LDA-i eGenes
   ## fn <- paste0(option, "/4_lm.results/2_", cell, "_lda", lda, "_trt", treat, "_stratified.FDR.txt") 
   ## res <- fread(fn, data.table=F)
   ## gene1 <- res%>%
   ##     drop_na(interaction_ANOVA.qvalue_stratified)%>%
   ##     filter(interaction_ANOVA.qvalue_stratified<0.1)%>%
   ##     dplyr::pull(ENSG)
   ## gene1 <- unique(gene1)

   res2 <- res%>%dplyr::filter(condition==ii) 

###    
   immue.traits <- read.table("../TWAS_overlap/PTWAS-immune-traits.txt")
   TWAS <- read.table("../TWAS_overlap/media-4_bioRxiv_TableS2.txt", header=T)
   TWAS2 <- TWAS%>%filter(Trait%in%immue.traits$V1) 
   twas0 <- unique(TWAS$Gene)
   twas2 <- unique(TWAS2$Gene)
   ###
   xx <- intersect(gene1, twas0) 
   df <- data.frame(x=xx)
})

###
immue.traits <- read.table("../TWAS_overlap/PTWAS-immune-traits.txt")
TWAS <- read.table("../TWAS_overlap/media-4_bioRxiv_TableS2.txt", header=T)
twas0 <- TWAS%>%filter(Trait%in%immue.traits$V1)
### op 1
opfn <- "./DiagLDA2/TWAS_immune-related-diseases.txt"
write.table(twas0, opfn, row.names=FALSE, quote=FALSE, sep="\t")
### op 2
opfn <- "./DiagLDA2/TWAS_all.txt"
write.table(TWAS, opfn, row.names=FALSE, quote=FALSE, sep="\t")

####################
### Asthma genes ###
####################

rm(list=ls())

option <- "DiagLDA2"
contrast_ls <- list("LPS"=c("CTRL", "LPS-EtOH"),
                    "LPS-DEX"=c("LPS-EtOH", "LPS-DEX"),
                    "PHA"=c("CTRL", "PHA-EtOH"),
                    "PHA-DEX"=c("PHA-EtOH", "PHA-DEX"))
dataset <- read.table("dataset_contrast.txt", header=F)
names(dataset) <- c("MCls", "LDA", "treats")

path <- keggGet("hsa05310") ##asthma
genes <- path[[1]]$GENE
entrez <- genes[(seq_len(length(genes))%%2)==1]
asthma <- bitr(entrez, fromType="ENTREZID", toType=c("ENSEMBL","SYMBOL"),
              OrgDb="org.Hs.eg.db") 
gene_asthma <- asthma$ENSEMBL
LDA_eGene <- map_dfr(1:16, function(i){
###    
   x <- dataset[i,]
   cell <- as.character(x[1])
   lda <- as.character(x[2])
   treat <- as.character(x[3])
   
   cat(paste(cell, lda, treat, sep="_"), "\n")

   ###New LDA-i eGenes
   fn <- paste0(option, "/4_lm.results/2_", cell, "_lda", lda, "_trt", treat, "_stratified.FDR.txt") 
   res <- fread(fn, data.table=F)
   gene1 <- res%>%
       drop_na(interaction_ANOVA.qvalue_stratified)%>%
       filter(interaction_ANOVA.qvalue_stratified<0.1)%>%
       dplyr::pull(ENSG)
   gene1 <- gsub("\\..*", "", unique(gene1))
   ## x <- intersect(gene1, gene_asthma)
   ## asthma%>%filter(ENSEMBL%in%x)
    
   df <- data.frame(LDA_eGene=length(gene1),
                    overlap=length(intersect(gene1,gene_asthma)))
})




############################
### compress LDA-results ###
############################
## dataset <- read.table("dataset_contrast.txt",header=F)
## names(dataset) <- c("MCls", "LDA", "treat")

## option <- "DiagLDA2"
## outdir <- paste(option, "/dynamic_eQTL_results/", sep="")
## if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)

## ###
## ### summary LDA-interaction eQTL results
## tmp <- lapply(1:nrow(dataset), function(i){
##    x <- dataset[i,]
##    cell <- x[1]
##    lda <- x[2]
##    treat <- x[3]
##    cat(paste(cell, lda, treat, sep="_"), "\n")
##    fn <- paste0(option, "/torm/4_lm.results/2_", cell, "_lda", lda, "_trt", treat, "_stratified.FDR.txt")
##    res <- fread(fn, data.table=F)%>%
##      drop_na(interaction_ANOVA.qvalue)
##    ###
##    res2 <- res%>%dplyr::select(ENSG, varID, interaction_effect, interaction_SE, interaction_ANOVA.pvalue, interaction_ANOVA.qvalue_stratified)
##    ###
##    opfn <- paste(outdir, cell, "_", lda, ".txt", sep="")
##    write.table(res2, opfn, row.names=F, col.names=T, quote=F)

##    system(paste("gzip ", opfn, sep=""))
   
## })
