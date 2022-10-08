##
library(tidyverse)
library(parallel)
library(purrr)
library(reshape)
library(qqman)
library(qvalue)
##
library(annotables)
library(biobroom)
library(clusterProfiler)
library(org.Hs.eg.db)
###
library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)
library(ggExtra)
library(gtable)
library(ggsignif)
library(pheatmap)
library(corrplot)
library(RColorBrewer)
library(patchwork)
library(viridis)
## theme_set(theme_grey())

rm(list=ls())

###
outdir <- "./5_twas_example.outs/"
if ( !(file.exists(outdir))) dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

 


#####################
### example genes ###
#####################


########################
### defined function ###
########################

adjGene <- function(cvt, center=T){
   cvt <- cvt%>%mutate(comb=paste(MCls, Batch, sep="_"))
   if(center){
      cvt <- cvt%>%group_by(comb)%>%mutate(yscale=y-mean(y,na.rm=T))%>%ungroup()
   }else{
      cvt <- cvt%>%mutate(yscale=y)
   }
}
###
adj2Gene <- function(cvt){
   cvt0 <- cvt%>%filter(treats=="CTRL")
   y0 <- cvt0$yscale
   names(y0) <- cvt0$sampleID
   y0[is.na(y0)] <- 0     
   cvt$y0 <- y0[cvt$sampleID]
   cvt <- cvt%>%mutate(yscale2=yscale-y0)
   cvt
}

###bulk, NB.mu, NB.phi
getData <- function(gene, datatype="bulk"){

### bulk data
if (datatype=="bulk"){
   fn <- "/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/6_DEG.CelltypeNew_output/Filter2/YtX_sel.comb.RData"
   load(fn)
   rn <- gsub("\\..*", "", rownames(YtX_sel))
   rownames(YtX_sel) <- rn
   X <- YtX_sel[rn%in%gene,]+1

   bti <- colnames(YtX_sel)
   cvt <- str_split(bti, "_", simplify=T)
   cvt <- data.frame(rn=bti, MCls=cvt[,1], treats=gsub("-EtOH", "", cvt[,2]), sampleID=cvt[,3], Batch=cvt[,4])

   fn <- "/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/6_DEG.CelltypeNew_output/Filter2/YtX.comb.RData"
   load(fn)
   counts <- colSums(YtX)
   counts <- counts[colnames(YtX_sel)]

   X <- (X/counts)*1e+06
   cvt$y <- log2(X)
   cvt <- adjGene(cvt, center=T)
} ###

### NB.mu
if(datatype=="NB.mu"){
###
   load("/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/10_RNA.Variance_output/tmp10/1.2_Sel.Bx.RData")
   rn <- gsub("\\..*", "", rownames(Bx))
   rownames(Bx) <- rn

   X <- Bx[rn%in%gene,]

   bti <- colnames(Bx)
   cvt <- str_split(bti, "_", simplify=T)
   cvt <- data.frame(rn=bti, MCls=cvt[,1], treats=gsub("-EtOH", "", cvt[,2]), sampleID=cvt[,3], Batch=cvt[,4])
   cvt$y <- log2(X)
   cvt <- adjGene(cvt, center=T) 
}

### NB.phi
if(datatype=="NB.phi"){
###
   load("/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/10_RNA.Variance_output/tmp10/1.2_Sel.PhxNew.RData")
   rn <- gsub("\\..*", "", rownames(PhxNew2))
   rownames(PhxNew2) <- rn
   X <- PhxNew2[rn%in%gene,]

   bti <- colnames(PhxNew2)
   cvt <- str_split(bti, "_", simplify=T)
   cvt <- data.frame(rn=bti, MCls=cvt[,1], treats=gsub("-EtOH", "", cvt[,2]), sampleID=cvt[,3], Batch=cvt[,4])
   X <- log2(X)
   cvt$y <- X
   cvt <- adjGene(cvt, center=T) 
}
cvt
}





####
###
                  
###
lab1 <- c("CTRL"="CTRL", 
          "LPS"="LPS", "LPS-DEX"="LPS+DEX",
          "PHA"="PHA", "PHA-DEX"="PHA+DEX")
col1 <- c("CTRL"="#828282", 
           "LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
           "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")

###
fn <- "/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/6_DEG.CelltypeNew_output/Filter2/2_meta.rds"
df1 <- read_rds(fn)
pval <- df1$pval
symbol <- rep("ns", nrow(df1))
symbol[pval<=0.05] <- "*"
symbol[pval<=0.01] <- "**"
symbol[pval<=0.001] <- "***"
df1$symbol <- symbol

anno <- bitr(unique(df1$gene), fromType="ENSEMBL", toType="SYMBOL", OrgDb=org.Hs.eg.db)


geneList <- c("IL1R2", "LTB", "TAP2", "CD52", "PDCD1")

###
### by cell-type
for (i in 1:length(geneList)){
###
   geneName <- geneList[i] 
   anno2 <- anno%>%filter(SYMBOL==geneName)
   symbol <- anno2[1,"SYMBOL"]
   ens <- anno2[1, "ENSEMBL"]

   cvt <- getData(gene=ens, datatype="bulk") 
    
   cat(symbol, "\n") 
   ### 
   MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell") 
   fig_ls <- lapply(1:4, function(j){
       ##
       oneMCl <- MCls[j]

   ### fig 1, gene mean
      cvt <- cvt%>%filter(MCls==oneMCl)
      cvt <- adj2Gene(cvt)
      cvt2 <- cvt%>%drop_na(y)%>%filter(treats!="CTRL")

### annotation
      sig_df <- df1%>%dplyr::filter(MCls==oneMCl, gene==ens, abs(beta)>0.5, qval<0.1)%>%
          dplyr::select(MCls, contrast, symbol)
      pos_df <- cvt2%>%group_by(treats)%>%summarise(ymax=max(yscale2, na.rm=T), .groups="drop")
      sig_df <- sig_df%>%
         left_join(pos_df, by=c("contrast"="treats"))%>%
         dplyr::rename("treats"="contrast") 

      p1 <- ggplot(cvt2,aes(x=factor(treats), y=yscale2, fill=treats))+
          geom_violin(width=0.8)+
          geom_boxplot(width=0.2, color="grey", outlier.shape=NA)+    
          geom_text(data=sig_df,
             aes(x=treats, y=ymax+0.2, label=symbol))+
          ylab(bquote(~log[2]~"(Expression)"))+
          scale_y_continuous(expand=expansion(mult=c(0.2, 0.2)))+
          scale_fill_manual("", values=col1, labels=lab1)+
          scale_x_discrete("", labels=lab1)+
          ggtitle(bquote(~.(oneMCl)~"("~italic(.(symbol))~")"))+
          theme_bw()+
          theme(axis.text.x=element_text(angle=-90, hjust=0, size=10),
          ## axis.text.x=element_blank(),
             axis.title.x=element_blank(),
             axis.ticks.x=element_blank(),
             axis.title.y=element_text(size=12),
             plot.title=element_text(hjust=0.5, size=14),
             legend.position="none")
       p1
  })
  ###  
  figfn <- paste(outdir, "Figure1.", i, "_", symbol, ".png", sep="")  
  png(figfn, width=800, height=800, res=120)
  print(plot_grid(plotlist=fig_ls, ncol=2))
  dev.off()
}###


col_MCls <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3", "NKcell"="#aa4b56", "Tcell"="#ffaa00")

###
###

geneName <- "IL1R2" 
anno2 <- anno%>%filter(SYMBOL==geneName)
symbol <- anno2[1,"SYMBOL"]
ens <- anno2[1, "ENSEMBL"]

oneMCl <- "Tcell"

cvt <- getData(gene=ens, datatype="bulk")%>%dplyr::filter(MCls==oneMCl) 
cvt <- adj2Gene(cvt)
cvt2 <- cvt%>%drop_na(y)%>%filter(treats!="CTRL")    


### annotation
sig_df <- df1%>%dplyr::filter(MCls==oneMCl, gene==ens, abs(beta)>0.5, qval<0.1)%>%
   dplyr::select(MCls, contrast, symbol)
pos_df <- cvt2%>%group_by(treats)%>%summarise(ymax=max(yscale2, na.rm=T), .groups="drop")
sig_df <- sig_df%>%
   left_join(pos_df, by=c("contrast"="treats"))%>%
   dplyr::rename("treats"="contrast")



p1 <- ggplot(cvt2,aes(x=factor(treats), y=yscale2, fill=treats))+
   geom_violin(width=0.8)+
   geom_boxplot(width=0.2, color="grey", outlier.shape=NA)+
   geom_signif(comparison=list(c("LPS", "LPS-DEX"), c("PHA", "PHA-DEX")),
       annotation=c(sig_df%>%filter(treats=="LPS-DEX")%>%pull(symbol),
                    sig_df%>%filter(treats=="PHA-DEX")%>%pull(symbol)),
       y_position=c(max(sig_df%>%filter(treats=="LPS")%>%pull(ymax)%>%as.numeric(),
                        sig_df%>%filter(treats=="LPS-DEX")%>%pull(ymax)%>%as.numeric())+0.5,
                    max(sig_df%>%filter(treats=="PHA")%>%pull(ymax)%>%as.numeric(),
                        sig_df%>%filter(treats=="PHA-DEX")%>%pull(ymax)%>%as.numeric())+0.5 ),
       tip_length=0.05, vjust=0, textsize=3)+   
   ## geom_text(data=sig_df, aes(x=treats, y=ymax+0.2, label=symbol))+
   ylab(bquote(~log[2]~"(Expression)"))+
   scale_y_continuous(expand=expansion(mult=c(0.2, 0.2)))+
   scale_fill_manual("", values=col1, labels=lab1)+
   scale_x_discrete("", labels=lab1)+
   ggtitle(bquote(~italic(.(symbol))~"("~.(oneMCl)~")"))+
   theme_bw()+
   theme(## axis.text.x=element_text(angle=-90, hjust=0, size=10),
         axis.text.x=element_text(size=8),
         axis.title.x=element_blank(),
         axis.ticks.x=element_blank(),
         axis.title.y=element_text(size=10),
         plot.title=element_text(hjust=0.5, size=12),
         panel.border=element_rect(color=col_MCls[oneMCl], size=1.5),
         legend.position="none")

###
figfn <- paste(outdir, "Figure2.1_IL1R2_coherent.png", sep="")
png(figfn, width=280, height=400, res=120)
p1
dev.off()


###
###
geneName <- "LTB" 
anno2 <- anno%>%filter(SYMBOL==geneName)
symbol <- anno2[1,"SYMBOL"]
ens <- anno2[1, "ENSEMBL"]

oneMCl <- "Tcell"

cvt <- getData(gene=ens, datatype="bulk")%>%dplyr::filter(MCls==oneMCl) 
cvt <- adj2Gene(cvt)
cvt2 <- cvt%>%drop_na(y)%>%filter(treats!="CTRL")    


### annotation
sig_df <- df1%>%dplyr::filter(MCls==oneMCl, gene==ens)%>%
   dplyr::select(MCls, contrast, symbol)
pos_df <- cvt2%>%group_by(treats)%>%summarise(ymax=max(yscale2, na.rm=T), .groups="drop")
sig_df <- sig_df%>%
   left_join(pos_df, by=c("contrast"="treats"))%>%
   dplyr::rename("treats"="contrast") 

p2 <- ggplot(cvt2,aes(x=factor(treats), y=yscale2, fill=treats))+
   geom_violin(width=0.8)+
   geom_boxplot(width=0.2, color="grey", outlier.shape=NA)+
   geom_signif(comparison=list(c("LPS", "LPS-DEX"), c("PHA", "PHA-DEX")),
       annotation=c(sig_df%>%filter(treats=="LPS-DEX")%>%pull(symbol),
                    sig_df%>%filter(treats=="PHA-DEX")%>%pull(symbol)),
       y_position=c(max(sig_df%>%filter(treats=="LPS")%>%pull(ymax)%>%as.numeric(),
                        sig_df%>%filter(treats=="LPS-DEX")%>%pull(ymax)%>%as.numeric())+0.2,
                    max(sig_df%>%filter(treats=="PHA")%>%pull(ymax)%>%as.numeric(),
                        sig_df%>%filter(treats=="PHA-DEX")%>%pull(ymax)%>%as.numeric())+0.2 ),
       tip_length=0.05, vjust=0, textsize=3)+   
   ## geom_text(data=sig_df, aes(x=treats, y=ymax+0.2, label=symbol))+
   ylab(bquote(~log[2]~"(Expression)"))+
   scale_y_continuous(expand=expansion(mult=c(0.2, 0.2)))+
   scale_fill_manual("", values=col1, labels=lab1)+
   scale_x_discrete("", labels=lab1)+
   ggtitle(bquote(~italic(.(symbol))~"("~.(oneMCl)~")"))+
   theme_bw()+
   theme(## axis.text.x=element_text(angle=-90, hjust=0, size=10),
         axis.text.x=element_text(size=8),
         axis.title.x=element_blank(),
         axis.ticks.x=element_blank(),
         axis.title.y=element_text(size=10),
         plot.title=element_text(hjust=0.5, size=12),
         panel.border=element_rect(color=col_MCls[oneMCl], size=1.5),
         legend.position="none")
 
figfn <- paste(outdir, "Figure2.2_LTB_coherent.png", sep="")
png(figfn, width=280, height=400, res=120)
p2
dev.off()


###
###
geneName <- "CD52"

anno2 <- anno%>%filter(SYMBOL==geneName)
symbol <- anno2[1,"SYMBOL"]
ens <- anno2[1, "ENSEMBL"]


cvt <- getData(gene=ens, datatype="bulk") 
cvt <- adj2Gene(cvt)
cvt2 <- cvt%>%drop_na(y)%>%filter(treats!="CTRL")%>%
    mutate(comb=paste(MCls, treats, sep="_"))    


### annotation
sig_df <- df1%>%dplyr::filter(gene==ens)%>%
   dplyr::select(MCls, treats=contrast, symbol)%>%mutate(comb=paste(MCls, treats, sep="_"))
###
pos_df <- cvt2%>%group_by(comb)%>%summarise(ymax=max(yscale2, na.rm=T), .groups="drop")
sig_df <- sig_df%>%left_join(pos_df, by="comb")

 
plotDF <- cvt2%>%filter(MCls=="Bcell", grepl("PHA", treats))
sig_df2 <- sig_df%>%filter(MCls=="Bcell", grepl("PHA", treats)) 
    
p3 <- ggplot(plotDF, aes(x=factor(treats), y=yscale2, fill=treats))+
   geom_violin(width=0.8)+
   geom_boxplot(width=0.2, color="grey", outlier.shape=NA)+
   geom_signif(comparison=list(c("PHA", "PHA-DEX")),
       annotation=sig_df2%>%filter(treats=="PHA-DEX")%>%pull(symbol),
       y_position=max(sig_df2%>%filter(treats=="PHA")%>%pull(ymax)%>%as.numeric(),
                      sig_df2%>%filter(treats=="PHA-DEX")%>%pull(ymax)%>%as.numeric())+0.2,
       tip_length=0.05, vjust=0, textsize=3)+      
   ## geom_text(data=sig_df2,aes(x=treats, y=ymax+0.2, label=symbol))+
   ylab(bquote(~log[2]~"(Expression)"))+
   scale_y_continuous(expand=expansion(mult=c(0.2, 0.2)))+
   scale_fill_manual("", values=col1, labels=lab1)+
   scale_x_discrete("", labels=lab1)+
   ggtitle(bquote(~italic(.(symbol))~"(Bcell)"))+
   theme_bw()+
   theme(## axis.text.x=element_text(angle=-90, hjust=0, size=10),
         axis.text.x=element_text(size=10),
         axis.title.x=element_blank(),
         axis.ticks.x=element_blank(),
         axis.title.y=element_text(size=10),
         plot.title=element_text(hjust=0.5, size=12),
         panel.border=element_rect(color=col_MCls["Bcell"], size=1.5),
         legend.position="none")

###
###
figfn <- paste(outdir, "Figure2.3_CD52_Bcell_inCoherent.png", sep="")
png(figfn, width=300, height=400, res=120)
p3
dev.off()



###
###
plotDF <- cvt2%>%filter(MCls=="NKcell", grepl("LPS", treats))
sig_df2 <- sig_df%>%filter(MCls=="NKcell", grepl("LPS", treats)) 
    
p4 <- ggplot(plotDF, aes(x=factor(treats), y=yscale2, fill=treats))+
   geom_violin(width=0.8)+
   geom_boxplot(width=0.2, color="grey", outlier.shape=NA)+
   geom_signif(comparison=list(c("LPS", "LPS-DEX")),
       annotation=sig_df2%>%filter(treats=="LPS-DEX")%>%pull(symbol),
       y_position=max(sig_df2%>%filter(treats=="LPS")%>%pull(ymax)%>%as.numeric(),
                      sig_df2%>%filter(treats=="LPS-DEX")%>%pull(ymax)%>%as.numeric())+0.2,
       tip_length=0.05, vjust=0, textsize=3)+          
   ## geom_text(data=sig_df2, aes(x=treats, y=ymax+0.2, label=symbol))+
   ylab(bquote(~log[2]~"(Expression)"))+
   scale_y_continuous(expand=expansion(mult=c(0.2, 0.2)))+
   scale_fill_manual("", values=col1, labels=lab1)+
   scale_x_discrete("", labels=lab1)+
   ggtitle(bquote(~italic(.(symbol))~"(NKcell)"))+
   theme_bw()+
   theme(## axis.text.x=element_text(angle=-90, hjust=0, size=10),
         axis.text.x=element_text(size=10),
         axis.title.x=element_blank(),
         axis.ticks.x=element_blank(),
         axis.title.y=element_text(size=10),
         plot.title=element_text(hjust=0.5, size=12),
         panel.border=element_rect(color=col_MCls["NKcell"], size=1.5),
         legend.position="none")

###
###
figfn <- paste(outdir, "Figure2.4_CD52_NKcell_Coherent.png", sep="")
png(figfn, width=300, height=400, res=120)
p4
dev.off()









## fig1 <- plot_grid(p1, p2, nrow=2,
##    align="v", axis="lr", rel_heights=c(1,1.2))


## fig3 <- plot_grid(p1, p2, nrow=2,
##     align="v", axis="lr", rel_heights=c(1,1.2))


###
## figfn <- "./10_RNA.Variance_output/tmp9_pub/Figure3.4_comb.png"
## png(figfn, width=1200, height=500, res=120)
## print(plot_grid(fig1, fig2, fig3, ncol=3,
##    labels=c("D", "E", "F"), label_x=c(0.05, 0.05, 0.1),
##    label_size=16, label_fontface="plain"))
## dev.off()


###
