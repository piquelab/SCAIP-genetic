### Example egene showing dynamical eQTLs
### 5/28/2021 By Julong wei
### last modified 11/10/2021, JW

library(tidyverse)
library(data.table)
library(parallel)
library(clusterProfiler)
library(KEGGREST)
library(qvalue)
library(qqman)
library(ggplot2)
library(grid)
library(gridExtra)
library(ggExtra)
library(plotrix)
#library(Seurat)
library(annotables)
library(org.Hs.eg.db)
    
### check
rm(list=ls())

outdir <- paste("../Figures_pub/Figure1/", sep="")
if (!file.exists(outdir))  dir.create(outdir, showWarnings=F, recursive=T)


####
####
option <- "DiagLDA2"
dataset <- read.table("dataset_contrast.txt", header=F)
names(dataset) <- c("MCls", "LDA", "treat")

col1 <- c("CTRL"="#828282", "LPS"="#fb9a99",
   "LPS-DEX"="#e31a1c", "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")

i <- 16
cell <- dataset[i,1]
lda <- dataset[i,2]
treat <- dataset[i,3]
col0 <- col1[[lda]]

treat2 <- gsub("-", "+", gsub("-EtOH", "", treat))


###
## normalized GE files:
fn <- paste0(option, "/1_normalized.data/", cell, "_lda", lda, "_trt", treat, ".bed.gz")
data <- fread(fn, header=T, sep="\t", stringsAsFactors=FALSE, data.table=F)
GE <- data[,5:ncol(data)]
rownames(GE) <- data$ID


## txt file with dosages:
fn <- paste0(option, "/3_genotypes/dosages/1_eQTL_signif_dosages_", cell, "_lda", lda, "_trt", treat, ".txt.gz")
dosages <- fread(fn, sep="\t", stringsAsFactors=F, data.table=F)
rownames(dosages) <- dosages$varID


### LDA-eGene
fn <- paste(option, "/examples/eGene_SNP.", i, "_", cell, ".lda",
   lda, ".trt", treat, ".txt", sep="")
res <- read.table(fn, header=TRUE)

## res2 <- res%>%dplyr::filter(symbol=="MRPL48", grepl("rs5792633", varID))
res2 <- res%>%dplyr::filter(symbol=="MRI1", grepl("rs35122230", varID))

## overlap with TWAS
## twas <- read.table("./DiagLDA2/TWAS_immune-related-diseases.txt", header=TRUE)%>%
##    mutate(ENSG2=gsub("\\..*", "", Gene))
## anno <- bitr(unique(twas$ENSG2), fromType="ENSEMBL", toType="SYMBOL", OrgDb=org.Hs.eg.db)
## twas <- twas%>%left_join(anno, by=c("ENSG2"="ENSEMBL"))
## gene <- unique(twas$ENSG2)

## res2 <- res%>%dplyr::filter(!ENSG2%in%gene) #symbol=="HLA-DRB5")


## ###
## reqtl <- read_rds("../mashr_eQTL/2_reQTL.output/reQTLs_All.infor.rds")%>%
##    dplyr::filter(condition=="Tcell_PHA-DEX:Tcell_PHA-EtOH", SYMBOL=="HLA-DRB5")
## ##
## res2 <- res%>%filter(symbol=="HLA-DRB5", varID%in%reqtl$varID)
## res2 <- res2%>%arrange(interaction_ANOVA.pvalue)

### athma genes
## path <- keggGet("hsa05310") ##asthma
## genes <- path[[1]]$GENE
## entrez <- genes[(seq_len(length(genes))%%2)==1]
## asthma <- bitr(entrez, fromType="ENTREZID", toType=c("ENSEMBL","SYMBOL"),
##               OrgDb="org.Hs.eg.db")
## asthma%>%filter(SYMBOL%in%res$symbol)

## res2 <- res%>%dplyr::filter(symbol=="ABO")
## res2 <- res%>%arrange(interaction_ANOVA.pvalue)

## res2 <- res%>%filter(symbol%in%gene_TWAS[-c(1,32,13,2,3,4,8,10,17,18)])#, grepl("rs71536576", varID))

### boxplot for bin separately
for (k in 1:nrow(res2)){
###
  ENSG <- res2$ENSG[k]
  varID <- res2$varID[k]
  symbol <- res2$symbol[k]  
  cat(ENSG, varID, "\n")  
  y <- GE[ENSG,]
  x <- round(dosages[varID, names(y)])
  df <- data.frame(rn=names(y), y=as.numeric(y), x=as.numeric(x))%>%
     mutate(Bin=gsub(".*_", "", rn))    
  ###

  # add genotype:
  ref <- strsplit(varID, ":")[[1]][3]
  alt <- gsub(";.*", "", strsplit(varID, ":")[[1]][4])
  rs <- gsub(".*;", "", varID)
  ###
  fig <- ggplot(data=df, aes(x=as.factor(x), y=y))+
     geom_boxplot(aes(alpha=Bin), fill=col0, outlier.shape=NA)+
   ## scale_fill_manual(values= c("0"="#4daf4a", "1"="#984ea3","2"="#ff7f00"))+
     scale_alpha_manual(values=c("1"=0.2, "2"=0.6, "3"=1), guide="none")+
     scale_x_discrete(rs,
         labels=c("0"=paste(ref, "/", ref, sep=""),
                  "1"=paste(ref, "/", alt, sep=""),
                  "2"=paste(alt, "/", alt, sep="")))+
     geom_smooth(method='lm', se=F, color="black", size=0.6, aes(group=Bin))+
     geom_jitter(width=0.25, size=0.8)+
     facet_grid(.~Bin, labeller=labeller(Bin=c("1"="1st tertile",
           "2"="2nd tertile", "3"="3rd tertile"))) +
     ylab("Normalized gene expresssion")+
     ggtitle(bquote(~.(cell)~.(treat2)~"("~italic(.(symbol))~")"))+
     theme_bw()+
     theme(axis.title=element_text(size=10),
           axis.text.x=element_text(size=8),
           axis.text.y=element_text(size=10),
           strip.text.x=element_text(size=13),
           plot.title=element_text(hjust=0.5, size=12))

     
  outdir2 <- paste(outdir, "/", symbol, sep="")
  if ( !file.exists(outdir2)) dir.create(outdir2,showWarnings=FALSE, recursive=TRUE)
    
  figfn <- paste(outdir2, "/", ENSG, "_", symbol, "_", varID, ".png", sep="")
  png(figfn, width=480, height=420, res=120)
  print(fig)
  dev.off()
}


## ### boxplots for bin together
## for (k in 1:20){
## ###    
##    ENSG <- res2$ENSG[k]
##    varID <- res2$varID[k]
##    symbol <- res2$symbol[k]  
##    cat(ENSG, varID, "\n")  
##    y <- GE[ENSG,]
##    x <- round(dosages[varID, names(y)])
##    df <- data.frame(rn=names(y), y=as.numeric(y), x=as.numeric(x))%>%
##       mutate(Bin=gsub(".*_", "", rn))     
##   m0 <- lm(y~0+as.numeric(x),data=df)
##   msum <- summary(m0)$coefficients  
##   xpos <- max(df$x)+0.2
##   ypos <- max(df$y)-0.1*max(df$y) 
##   pval <- round(as.numeric(msum[1,4]), digits=3)
##   beta <- round(as.numeric(msum[1,1]), digits=3)    
##   fig2 <- boxplot2(df, varID, symbol, col0, divided=FALSE)+
##      ggtitle(paste(cell,"_lda", lda, "_trt", treat, sep=""))+
##      annotate("text", x=xpos, y=ypos,
##         label=bquote(beta==.(beta)~", "~italic(pval)==.(pval)), size=3)+  
##      theme_bw()+
##      theme(plot.title=element_text(hjust=0.5, size=12))

##   outdir2 <- paste(outdir, "/", symbol, sep="")
##   if ( !file.exists(outdir2)) dir.create(outdir2,showWarnings=FALSE, recursive=TRUE)    
##   figfn <- paste(outdir2, "/", ENSG, "_", symbol, "_", varID, ".2.png", sep="")
##   png(figfn, width=420, height=600, res=120)
##   print(fig2)
##   dev.off()
## }
## ###

### eGenes
## fn <- paste(option, "/plots/eGene_SNP.", i, "_", cell, ".lda",
##    lda, ".trt", treat, ".txt", sep="")
## res <- read.table(fn, header=TRUE)
## fn <- paste("../eQTL/egenes/", cell, "_", treat, "_egenes.txt", sep="")
## gene4 <- read.table(fn)%>%dplyr::pull(V1)


###
### example genes use for paper


## ######################################################################
## ### 3. boxplots showing example genes with eQTL for each contidion ###
## ######################################################################

## rm(list=ls())
## option <- "DiagLDA2"
## outdir <- paste(option, "/plots/boxplot_eQTL/", sep="")
## if ( !file.exists(outdir)) dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

## ##
## contrast_ls <- list("LPS"=c("CTRL", "LPS-EtOH"),
##    "LPS-DEX"=c("LPS-EtOH", "LPS-DEX"),
##    "PHA"=c("CTRL", "PHA-EtOH"),
##    "PHA-DEX"=c("PHA-EtOH", "PHA-DEX"))

## col1 <- c("CTRL"="#828282", "LPS-EtOH"="#fb9a99",
##    "LPS-DEX"="#e31a1c", "PHA-EtOH"="#a6cee3", "PHA-DEX"="#1f78b4")

## dataset <- read.table("dataset_contrast.txt", header=F)
## names(dataset) <- c("MCls", "LDA", "treat")

## i <- 16
## cell <- dataset[i,1]
## lda <- dataset[i,2]
## treat <- dataset[i,3]
## contrast <- contrast_ls[[lda]]

## fn <- paste(option, "/plots/eGene_SNP.", i, "_", cell,
##    ".lda", lda, ".trt", treat, ".txt", sep="")
## res <- read.table(fn, header=TRUE)


## ##
## # load the needed files:
## # 1. normalized GE files:
## # treatment:
## treat1 <- contrast[2]
## fn <- paste("../eQTL/normalized_GE_residuals/",cell,"_", treat1, ".bed.gz", sep="")
## X <- fread(fn, header=T, sep="\t", data.table=F, stringsAsFactors=FALSE)
## GE_treat1 <- X[,5:ncol(X)]
## rownames(GE_treat1) <- X$ID
## fn <- paste("../eQTL/eQTL_output/", cell, "_", treat1, ".GEPC3.nominals.eQTL.txt.gz", sep="")
## eQTL_treat1 <- fread(fn, data.table=F,
##    col.names=c("ENSG", "varID", "position", "pvalue", "slope_sc"),
##    stringsAsFactors=FALSE)%>%
##    mutate(qvalue=qvalue(pvalue)$qvalues, pairXY=paste(ENSG, varID, sep="_"))
## rownames(eQTL_treat1) <- eQTL_treat1$pairXY

## ### contrast
## treat0 <- contrast[1]
## fn <- paste("../eQTL/normalized_GE_residuals/", cell, "_", treat0, ".bed.gz", sep="")
## X <- fread(fn, header=T, sep="\t", data.table=F, stringsAsFactors=FALSE)
## GE_treat0 <- X[,5:ncol(X)]
## rownames(GE_treat0) <- X$ID
## ##
## fn <- paste("../eQTL/eQTL_output/", cell, "_", treat0, ".GEPC3.nominals.eQTL.txt.gz", sep="")
## eQTL_treat0 <- fread(fn, data.table=F,
##    col.names=c("ENSG", "varID", "position", "pvalue", "slope_sc"),
##    stringsAsFactors=FALSE)%>%
##    mutate(qvalue=qvalue(pvalue)$qvalues, pairXY=paste(ENSG, varID, sep="_"))
## rownames(eQTL_treat0) <- eQTL_treat0$pairXY

## ### CTRL
## fn <- paste("../eQTL/normalized_GE_residuals/", cell, "_CTRL.bed.gz", sep="")
## X <- fread(fn, header=T, sep="\t", data.table=F, stringsAsFactors=FALSE)
## GE_ctrl <- X[,5:ncol(X)]
## rownames(GE_ctrl) <- X$ID
## fn <- paste("../eQTL/eQTL_output/", cell, "_CTRL.GEPC3.nominals.eQTL.txt.gz", sep="")
## eQTL_ctrl <- fread(fn, data.table=FALSE,
##    col.names=c("ENSG", "varID", "position", "pvalue", "slope_sc"),
##    stringsAsFactors=FALSE)%>%
##    mutate(qvalue=qvalue(pvalue)$qvalues, pairXY=paste(ENSG, varID, sep="_"))   
## rownames(eQTL_ctrl) <- eQTL_ctrl$pairXY


## ## fn <- "/wsu/home/groups/piquelab/SCAIP/eQTL/FastQTL/SCAIP1-6_filtered.vcf.gz"
## ## system(paste0("bcftools query -f '%CHROM\t%POS\t%ID[\t%DS]\n' ", fn, " > ./check_output/SCAIP1-6_filtered.vcf.txt") )
## ## system("bgzip ./check_output/SCAIP1-6_filtered.vcf.txt")       
## ## system(paste0("bcftools query -l ", fn, " > ./check_output/sample.txt" ))
## # 2. txt file with dosages
## dosages <- fread("./check_output/SCAIP1-6_filtered.vcf.txt.gz", header=F, data.table=F, stringsAsFactors=FALSE)
## IDs <- read.table("./check_output/sample.txt", stringsAsFactors=FALSE)[,1]
## colnames(dosages) <- c("chr", "pos", "varID", IDs)
## dosages <- dosages[!duplicated(dosages$varID),]
## rownames(dosages) <- dosages$varID

## ## fn <- paste0("./genotypes/dosages/1_eQTL_signif_dosages_", cell, "_lda", lda, "_trt", treat, ".txt.gz")
## ## dosages2 <- fread(fn, sep="\t", stringsAsFactors=F, data.table=F)

## ## dosages <- read.table(paste0("../eQTL/eQTL_coordinates/all_uniq_eQTL_dosages.txt"), sep="\t", header=T,stringsAsFactors=F)
## ## colnames(dosages) <- gsub("[.]","-",colnames(dosages))
## # remove repeated dosages (where do they come from, though??):

## ### plot function
## boxplot2 <- function(df, varID, symbol, col0){
## ###    
##   ref <- strsplit(varID, ":")[[1]][3]
##   alt <- gsub(";.*","",strsplit(varID, ":")[[1]][4])

##   ## col1 <- c("CTRL"="#828282",
##   ##    "LPS-EtOH"="#fb9a99", "LPS-DEX"="#e31a1c",
##   ##    "PHA-EtOH"="#a6cee3", "PHA-DEX"="#1f78b4")
##   ## col0 <- "#828282"
##   fig <- ggplot(data=df, aes(x=as.factor(x), y=y))+
##      geom_boxplot(alpha=0.5, fill=col0, outlier.shape=NA)+
##      ## scale_alpha_manual(values=c("0"=0.35, "1"=0.7, "2"=1))+   
##      scale_x_discrete("", labels=c("0"=paste(ref, "/", ref, sep=""),
##                                  "1"=paste(ref, "/", alt, sep=""),
##                                  "2"=paste(alt, "/", alt, sep="")))+
##    geom_smooth(method='lm')+
##    geom_jitter(width=0.25, size=1)+
##    ylab(bquote(~italic(.(symbol))~" normalized gene expresssion"))+
##    theme_bw()+
##    theme(legend.position="none",
##          axis.title.x = element_blank(),
##          axis.title.y=element_text(size=10))
##    fig
## } ###

## res2 <- res%>%filter(symbol=="CDKAL1", grepl("rs9348432",  varID))

## for (k in 1:nrow(res2)){

##     ENSG <- res2$ENSG[k]
##     varID <- res2$varID[k]
##     symbol <- res2$symbol[k]

##     pair <- paste(ENSG, "_", varID, sep="")

## ### treatment
##    y <- GE_treat1[ENSG,]
##    x <- dosages[varID, names(y)]
##    df <- data.frame(y=as.numeric(y), x=round(as.numeric(x)))%>%drop_na(x,y)
##    xpos <- max(df$x)+0.2
##    ypos <- max(df$y)-0.1*max(df$y) 
##    pval <- round(as.numeric(eQTL_treat1[pair,"pvalue"]), digits=3)
##    beta <- round(as.numeric(eQTL_treat1[pair,"slope_sc"]), digits=3)
##    fdr <- round(as.numeric(eQTL_treat1[pair,"qvalue"]), digits=3) 
##    col0 <- col1[[treat1]] 
##    fig1 <- boxplot2(df, varID, symbol, col0)+
##       ggtitle(paste(cell,"_", treat1, sep=""))+
##       annotate("text", x=xpos, y=ypos,
##         label=bquote(beta==.(beta)~", "~italic(pval)==.(pval)~
##           ", FDR"==.(fdr)), size=3)+ 
##       theme(plot.title=element_text(hjust=0.5, size=12))
    
##    figfn <- paste(option, "/plots/boxplot_eQTL/", cell, i, "_", ENSG, "_",
##       symbol, "_", varID, "_", treat1, ".png", sep="")
##    png(figfn, width=420, height=600, res=120)
##    print(fig1)
##    dev.off()

## ### contrast
##    y <- GE_treat0[ENSG,]
##    x <- dosages[varID, names(y)]    
##    df <- data.frame(y=as.numeric(y), x=round(as.numeric(x)))%>%drop_na(x,y)
##    xpos <- max(df$x)+0.2
##    ypos <- max(df$y)-0.1*max(df$y) 
##    pval <- round(as.numeric(eQTL_treat0[pair,"pvalue"]), digits=3)
##    beta <- round(as.numeric(eQTL_treat0[pair,"slope_sc"]), digits=3)
##    fdr <- round(as.numeric(eQTL_treat0[pair,"qvalue"]), digits=3) 
##    col0 <- col1[[treat0]] 
##    fig2 <- boxplot2(df, varID, symbol, col0)+
##       ggtitle(paste(cell,"_", treat0, sep=""))+
##       annotate("text", x=xpos, y=ypos,
##         label=bquote(beta==.(beta)~", "~italic(pval)==.(pval)~
##         ", FDR"==.(fdr)), size=3)+
##       theme(plot.title=element_text(hjust=0.5, size=12))
    
##   figfn <- paste(option, "/plots/boxplot_eQTL/", cell, i, "_", ENSG, "_",
##      symbol, "_", varID, "_", treat0, ".png", sep="")
##   png(figfn, width=420, height=600, res=120)
##   print(fig2)
##   dev.off()

## ### CTRL
##    if ( treat0!="CTRL"){ 
##       y <- GE_ctrl[ENSG,]
##       x <- dosages[varID, names(y)]   
##       df <- data.frame(y=as.numeric(y), x=round(as.numeric(x)))%>%drop_na(x,y)
##       xpos <- max(df$x)+0.2
##       ypos <- max(df$y)-0.1*max(df$y) 
##       pval <- round(as.numeric(eQTL_ctrl[pair, "pvalue"]), digits=3)
##       beta <- round(as.numeric(eQTL_ctrl[pair, "slope_sc"]), digits=3) 
##       fdr <- round(as.numeric(eQTL_ctrl[pair, "qvalue"]), digits=3)
##       col0 <- col1[["CTRL"]] 
##       fig3 <- boxplot2(df, varID, symbol, col0)+
##         ggtitle(paste(cell,"_CTRL", sep=""))+
##         annotate("text", x=xpos, y=ypos,
##           label=bquote(beta==.(beta)~", "~italic(pval)==.(pval)~
##           ", FDR"==.(fdr)), size=3)+        
##         theme(plot.title=element_text(hjust=0.5, size=12))
    
##      figfn <- paste(option, "/plots/boxplot_eQTL/", cell, i, "_", ENSG, "_",
##         symbol, "_", varID, "_CTRL.png", sep="")
##      png(figfn, width=420, height=600, res=120)
##      print(fig3)
##      dev.off()
##   }
    
## }
