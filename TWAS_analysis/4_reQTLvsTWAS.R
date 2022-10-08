##
library(tidyverse)
library(data.table)
library(clusterProfiler)
library(org.Hs.eg.db)

library(ComplexHeatmap)
library(RColorBrewer)
library(viridis)
library(circlize)
library(cowplot)
library(ggrepel)

rm(list=ls())

outdir <- "./4_reQTLvsTWAS.outs/"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)


###
### TWAS overlap:
# load the list of the traits of interest:
## immune.traits <- read.table("../TWAS_overlap/PTWAS-immune-traits.txt", sep="\t", header=FALSE, comment="=")[,1]
## # read in the TWAS file from William:
## TWAS <- fread("../TWAS_overlap/media-4_bioRxiv_TableS2.txt", data.table=F, stringsAsFactors=F)
## TWAS <- TWAS%>%mutate(ENSG2=gsub("\\..*", "", Gene))
## TWAS2 <- TWAS%>%dplyr::filter(Trait%in%immune.traits)
    

###
### PTWAS results-estimated effects
prefix <- "/wsu/home/groups/piquelab/gtex_v8/PTWAS_analysis_new/"
fn <- paste(prefix, "PTWAS_gtex_v8.effect_estimates.txt.gz", sep="")
res <- fread(fn, header=F, data.table=F, quote="")
colnames(res)<-c("Trait_Gene_Tissue_Combination", "Number_of_eQTL_Signal_Clusters",
    "Posterior_Expected_Number_of_eQTLs", "Strong_eQTL_Signal_Clusters_(_PIP_>_0.50_)",
    "Estimated_Effect", "Standard_Error", "I2")

###
res2 <- res%>%drop_na(Trait_Gene_Tissue_Combination)%>%
    mutate(Trait=gsub("_ENSG.*", "", Trait_Gene_Tissue_Combination),
           Gene=gsub("\\..*", "", gsub(".*_ENSG", "ENSG", Trait_Gene_Tissue_Combination)),
           Tissue=gsub(".*\\..{1,2}_", "", gsub(".*_ENSG", "ENSG", Trait_Gene_Tissue_Combination)))
###

### PTWAS-identified genes 
fn <- paste(prefix, "PTWAS.gtex_v8.all_traits.tgenes.fdr_05.txt", sep="")
x <- fread(fn, header=F, data.table=F)%>%
    mutate(Trait_Gene_Tissue_Combination=paste(V1, V2, V3, sep="_"),
           Gene=gsub("\\..*", "", V2))

## x2 <- x%>%filter(grepl("Depressive_Symptoms", V1))


x2 <- x%>%filter(V3=="Whole_Blood", grepl("Asthma|asthma", V1))## %>%
     ## dplyr::select(Trait_Gene_Tissue_Combination, FDR=V4, Gene, V2)

###
res3 <- res2%>%filter(Tissue=="Whole_Blood", grepl("Asthma|asthma", Trait), Gene%in%x2$Gene)%>%
    mutate(zscore=Estimated_Effect/Standard_Error)%>%
    dplyr::select(Trait_Gene_Tissue_Combination, Gene, Trait,
                  beta_TWAS=Estimated_Effect, zscore_TWAS=zscore)

res3 <- res3%>%
    group_by(Gene)%>%
    slice_max(order_by=abs(beta_TWAS), n=1)%>%
    distinct(Gene, .keep_all=T)%>%as.data.frame() ### 78 TWAS genes related to Asthma

geneID <- bitr(res3$Gene, fromType="ENSEMBL", toType="SYMBOL", OrgDb=org.Hs.eg.db)

res3 <- res3%>%left_join(geneID, by=c("Gene"="ENSEMBL"))


## x <- res4%>%group_by(Gene)%>%summarise(nt=n(), nsign=abs(sum(direction)), .groups="drop")%>%as.data.frame()

## res2$Trait <- gsub("_ENSG.*", "", res2$Trait_Gene_Tissue_Combination)
## res2$Gene <- gsub("\\..*", "", gsub(".*_ENSG", "ENSG", res2$Trait_Gene_Tissue_Combination))
## res2$Tissue <- gsub(".*\\..{1,2}_", "", gsub(".*_ENSG", "ENSG", res2$Trait_Gene_Tissue_Combination))
##


## dap <- read.table(fn, fill=T, row.names=NULL, header=F)

## dap <- map_dfr(16:20, function(i){
##    ##
##    condition <- conditions[i]
##    fn <- paste(prefix, "5_summary.outs/dap-g_pct_0.02/summary/", condition, "_gene_FDR.txt", sep="")
##    res2 <- read.table(fn, header=T)%>%inner_join(geneSel, by=c("gene"="Gene"))%>%filter(FDR<0.3)
##    res2
## })    



outdir <- "./4_reQTLvsTWAS.outs/"

############################
### differential results ###
############################

fn <- "/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/6_DEG.CelltypeNew_output/Filter2/2_meta.rds"
DE <- read_rds(fn)
DE2 <- DE%>%filter(gene%in%unique(res3$Gene))%>%
    mutate(gr2=ifelse(qval<0.1&abs(beta)>0.5, contrast, "0"), zscore_DE=beta/stderr)%>% 
    dplyr::select(MCls, contrast, gr2, Gene=gene, rn, beta_DE=beta, zscore_DE)
### 
plotDF <- DE2%>%inner_join(res3, by="Gene")

x <- plotDF%>%filter(gr2!="0")
length(unique(x$Gene))


## tmp <- plotDF%>%filter(gr2>0, SYMBOL!="IL1R1")%>%
##     dplyr::select(MCls, contrast, Gene, SYMBOL, zscore_DE, zscore_TWAS)%>%
##     mutate(coherent=sign(zscore_DE)*sign(zscore_TWAS))

## tmp2 <- tmp%>%group_by(MCls, contrast, coherent)%>%
##     summarize(ny=n(),.groups="drop")%>%as.data.frame()

## tmp2%>%pivot_wider(id_cols=MCls, names_from=contrast, values_from=ny)

## tmp%>%filter(MCls=="Tcell")

## tmp <- plotDF%>%filter(gr2!=0)%>%
##     group_by(MCls, contrast)%>%summarise(ngene=n(), .groups="drop")%>%
##     pivot_wider(id_cols="MCls", names_from="contrast", values_from="ngene") 
## length(unique(tmp$Gene))


### DEX 
## MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
## fig_ls <- lapply(MCls, function(oneMCl){
## ###    
##    plotDF2 <- plotDF%>%filter(MCls==oneMCl)
##    anno2 <- plotDF2%>%filter(gr2>0)
##    ## 
##    p0 <- ggplot(plotDF2, aes(x=zscore_DE, y=zscore_TWAS, color=factor(gr2)))+
##       geom_point(aes(shape=gr2, size=gr2))+
##       geom_text_repel(data=anno2, mapping=aes(x=zscore_DE, y=zscore_TWAS, label=SYMBOL, color=gr2),
##                       size=3, max.overlaps=10)+  
##       scale_color_manual(values=c("0"="grey", "LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
##                                   "PHA"="#a6cee3", "PHA-DEX"="#1f78b4"))+
##       scale_shape_manual(values=c("0"=4, "LPS"=2, "LPS-DEX"=2, "PHA"=2, "PHA-DEX"=2))+
##       scale_size_manual(values=c("0"=1, "LPS"=2.5, "LPS-DEX"=2.5, "PHA"=2.5, "PHA-DEX"=2.5))+ 
##       xlab("zscore of DE")+
##       ylab("zscore of TWAS")+
##       geom_hline(yintercept=0, linetype="dashed", color="grey30")+
##       geom_vline(xintercept=0, linetype="dashed", color="grey30")+
##       ggtitle(oneMCl)+ 
##       theme_bw()+
##       theme(legend.position="none", plot.title=element_text(hjust=0.5)) 
##     ###
##     p0
## })

## figfn <- paste(outdir, "Figure1.1_twas_DEG.png", sep="")
## png(figfn, width=720, height=720, res=100)
## plot_grid(plotlist=fig_ls)
## dev.off()


###
###
plotDF2 <- plotDF%>%dplyr::filter(gr2!="0", grepl("DEX", contrast), zscore_TWAS!=0)%>%
   mutate(gr_fill=ifelse(contrast=="PHA-DEX", MCls, "0")) 

p <- ggplot(plotDF2, aes(x=zscore_DE, y=zscore_TWAS, color=factor(MCls), fill=factor(gr_fill)))+
      geom_point(shape=24)+
      geom_text_repel(data=plotDF2,
          mapping=aes(x=zscore_DE, y=zscore_TWAS, label=SYMBOL, color=factor(MCls)), 
          size=3, max.overlaps=10)+  
      scale_color_manual(values=c("Bcell"="#4daf4a", "Monocyte"="#984ea3",
         "NKcell"="#aa4b56", "Tcell"="#ffaa00"), guide="none")+
      scale_fill_manual(values=c("0"="NA", "Bcell"="#4daf4a", "Monocyte"="#984ea3",
         "NKcell"="#aa4b56", "Tcell"="#ffaa00"),
         breaks=c("0", "Bcell"),
         labels=c("0"="LPS+DEX", "Bcell"="PHA+DEX"),
         guide=guide_legend(override.aes=list(color="#4daf4a")))+
      xlab("zscore of DE")+
      ylab("zscore of asthma TWAS")+
      xlim(-50, 20)+ylim(-8,5)+
      geom_hline(yintercept=0, linetype="dashed", color="grey30")+
      geom_vline(xintercept=0, linetype="dashed", color="grey30")+
      theme_bw()+
      theme(legend.position=c(0.25, 0.85),
            legend.title=element_blank(),
            legend.background=element_blank(),
            legend.box.background=element_blank(),
            legend.key=element_blank(),
            plot.title=element_text(hjust=0.5),
            axis.title=element_text(size=10))
###
figfn <- paste(outdir, "Figure1.1_twas_DEG.png", sep="")
png(figfn, width=480, height=420, res=120)
print(p)
dev.off()
####



###
p2 <- ggplot(plotDF2, aes(x=zscore_DE, y=zscore_TWAS, color=factor(MCls)))+
      geom_point(shape=24)+
      geom_text_repel(data=plotDF2,
          mapping=aes(x=zscore_DE, y=zscore_TWAS, label=SYMBOL, color=factor(MCls)), 
          size=2.5, max.overlaps=10)+  
      scale_color_manual(values=c("Bcell"="#4daf4a", "Monocyte"="#984ea3",
         "NKcell"="#aa4b56", "Tcell"="#ffaa00"), guide="none")+
      xlab("zscore of DE")+
      ylab("zscore of asthma TWAS")+
      xlim(-50, 20)+ylim(-8,5)+
      geom_hline(yintercept=0, linetype="dashed", color="grey30")+
      geom_vline(xintercept=0, linetype="dashed", color="grey30")+
      facet_grid(MCls~contrast, labeller=labeller(contrast=c("LPS-DEX"="LPS+DEX", "PHA-DEX"="PHA+DEX")))+ 
      theme_bw()+
      theme(## legend.position=c(0.25, 0.85),
            ## legend.title=element_blank(),
            ## legend.background=element_blank(),
            ## legend.box.background=element_blank(),
            ## legend.key=element_blank(),
            ## plot.title=element_text(hjust=0.5),
            strip.text=element_text(size=12),
            axis.title=element_text(size=10))
### 
figfn <- paste(outdir, "Figure1.2_twas_DEG.png", sep="")
png(figfn, width=420, height=600, res=120)
print(p2)
dev.off()



###
### LPS vs CTRL, PHA vs CTRL

plotDF2 <- plotDF%>%dplyr::filter(gr2!="0", !grepl("DEX", contrast), zscore_TWAS!=0)


p3 <- ggplot(plotDF2, aes(x=zscore_DE, y=zscore_TWAS, color=factor(MCls)))+
      geom_point(shape=24)+
      geom_text_repel(data=plotDF2,
          mapping=aes(x=zscore_DE, y=zscore_TWAS, label=SYMBOL, color=factor(MCls)), 
          size=2.5, max.overlaps=10)+  
      scale_color_manual(values=c("Bcell"="#4daf4a", "Monocyte"="#984ea3",
         "NKcell"="#aa4b56", "Tcell"="#ffaa00"), guide="none")+
      xlab("zscore of DE")+
      ylab("zscore of asthma TWAS")+
      xlim(-20, 20)+ylim(-8,5)+
      geom_hline(yintercept=0, linetype="dashed", color="grey30")+
      geom_vline(xintercept=0, linetype="dashed", color="grey30")+
      facet_grid(MCls~contrast)+ ##, labeller=labeller(contrast=c("LPS-DEX"="LPS+DEX", "PHA-DEX"="PHA+DEX")))+ 
      theme_bw()+
      theme(## legend.position=c(0.25, 0.85),
            ## legend.title=element_blank(),
            ## legend.background=element_blank(),
            ## legend.box.background=element_blank(),
            ## legend.key=element_blank(),
            ## plot.title=element_text(hjust=0.5),
            strip.text=element_text(size=12),
            axis.title=element_text(size=10))
### 
figfn <- paste(outdir, "Figure1.3_twas_DEG.png", sep="")
png(figfn, width=420, height=600, res=120)
print(p3)
dev.off()




###
###
## tmp2 <- plotDF2%>% dplyr::select(MCls, contrast, Gene, SYMBOL, zscore_DE, zscore_TWAS)%>%
##     mutate(coherent=sign(zscore_DE)*sign(zscore_TWAS))

## tmp2%>%filter(coherent==-1)%>%dplyr::pull(Gene)%>%unique()%>%length()

## tmp2 <- tmp%>%group_by(MCls, contrast, coherent)%>%
##     summarize(ny=n(),.groups="drop")%>%as.data.frame()



###############
### Heatmap ###
###############

## geneSel <- DE2%>%filter(gr2!="0")%>%dplyr::pull(Gene)%>%unique()
## ###  
## mat2 <- DE2%>%filter(Gene%in%geneSel)%>%mutate(conditions=paste(MCls, contrast, sep="_"))%>%
##     pivot_wider(id_cols=Gene, names_from=conditions, values_from=zscore_DE)%>%as.data.frame()
## rownames(mat2) <- mat2$Gene
## mat2 <- as.matrix(mat2[,-1])

## geneAnnot <- res3%>%
##     filter(Gene%in%geneSel)%>%
##     arrange(desc(zscore_TWAS))%>%dplyr::select(Gene, SYMBOL, zscore_TWAS)
## ### sort
## mat2 <- mat2[geneAnnot$Gene,]

## New_order <- c("Bcell_LPS", "Bcell_PHA", "Monocyte_LPS", "Monocyte_PHA", "NKcell_LPS", "NKcell_PHA",
##                "Tcell_LPS", "Tcell_PHA",
##                "Bcell_LPS-DEX", "Bcell_PHA-DEX", "Monocyte_LPS-DEX", "Monocyte_PHA-DEX",
##                "NKcell_LPS-DEX", "NKcell_PHA-DEX", "Tcell_LPS-DEX", "Tcell_PHA-DEX")
## mat2 <- mat2[, New_order]

## ###
## ### column annotation
## anno_df <- data.frame("celltype"=gsub("_.*", "", colnames(mat2)),
##     "Contrast"=gsub(".*_", "", colnames(mat2))) 
## ###
## col_ha <- HeatmapAnnotation(df=anno_df,                         
##    col=list("celltype"=c("Bcell"="#4daf4a", "Monocyte"="#984ea3", "NKcell"="#aa4b56", "Tcell"="#ffaa00"),
##    "Contrast"=c("LPS"="#fb9a99", "LPS-DEX"="#e31a1c", "PHA"="#a6cee3", "PHA-DEX"="#1f78b4")),
##    show_legend=T,
##    annotation_name_gp=gpar(fontsize=8))

## ###
## ### row annotation, z-score of TWAS
## z2 <- geneAnnot$zscore_TWAS
## breaks <- quantile(z2, probs=seq(0,1,length.out=10))
## col2 <- colorRamp2(breaks, colorRampPalette(rev(brewer.pal(n=4, name="RdBu")))(10) )
 
## row_ha <- rowAnnotation(zscore_TWAS=geneAnnot$zscore_TWAS,
##    col=list(zscore_TWAS=col2),
##    annotation_legend_param=list(zscore_TWAS=list(grid_width=unit(0.3, "cm"),
##        labels_gp=gpar(fontsize=8), title_gp=gpar(fontsize=10), title="zscore_TWAS")),
##    width=unit(0.5, "cm"), annotation_name_gp=gpar(fontsize=8))


## ###
## z3 <- as.vector(mat2)
## breaks <- quantile(z2, probs=seq(0,1,length.out=20))
## col_mat <- colorRamp2(breaks, colorRampPalette(rev(brewer.pal(n=4, name="RdBu")))(20))

## fig1 <- Heatmap(mat2, col=col_mat,
##    cluster_rows=F, cluster_columns=F,
##    show_row_names=T, show_column_names=F,
##    top_annotation=col_ha,
##    right_annotation=row_ha,
##    heatmap_legend_param=list(title="zscore_DE",
##       title_gp=gpar(fontsize=8),
##       grid_width=unit(0.3, "cm")),
##    raster_device="png")
## ###

## figfn <- paste(outdir, "Figure1.3_heatmap.png", sep="")
## png(figfn, height=720, width=720, res=120)
## fig1 <- draw(fig1)
## dev.off()




#############################
### DVG compared to TWAS  ###
#############################

fn <- "/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/10_RNA.Variance_output/tmp9/3_phiNew.meta"
DE <- read.table(fn, header=T)
DE2 <- DE%>%filter(gene%in%unique(res3$Gene))%>%
    mutate(gr2=ifelse(qval<0.1&abs(beta)>0.5, contrast, "0"), zscore_DE=beta/stderr)%>% 
    dplyr::select(MCls, contrast, gr2, Gene=gene, rn, beta_DE=beta, zscore_DE)
### 
plotDF <- DE2%>%inner_join(res3, by="Gene")

tmp <- plotDF%>%filter(gr2>0)%>%dplyr::select(MCls, contrast, Gene, SYMBOL, zscore_DE, zscore_TWAS)
## tmp%>%filter(MCls=="Tcell")

## tmp <- plotDF%>%filter(gr2!=0)%>%
##     group_by(MCls, contrast)%>%summarise(ngene=n(), .groups="drop")%>%
##     pivot_wider(id_cols="MCls", names_from="contrast", values_from="ngene") 
## length(unique(tmp$Gene))


### DEX
 
MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
fig_ls <- lapply(MCls, function(oneMCl){
###    
   plotDF2 <- plotDF%>%filter(MCls==oneMCl)
   anno2 <- plotDF2%>%filter(gr2>0)
   ## 
   p0 <- ggplot(plotDF2, aes(x=zscore_DE, y=zscore_TWAS, color=factor(gr2)))+
      geom_point(aes(shape=gr2, size=gr2))+
      geom_text_repel(data=anno2, mapping=aes(x=zscore_DE, y=zscore_TWAS, label=SYMBOL, color=gr2),
                      size=3, max.overlaps=10)+  
      scale_color_manual(values=c("0"="grey", "LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
                                  "PHA"="#a6cee3", "PHA-DEX"="#1f78b4"))+
      scale_shape_manual(values=c("0"=4, "LPS"=2, "LPS-DEX"=2, "PHA"=2, "PHA-DEX"=2))+
      scale_size_manual(values=c("0"=1, "LPS"=2.5, "LPS-DEX"=2.5, "PHA"=2.5, "PHA-DEX"=2.5))+ 
      xlab("zscore of DV")+
      ylab("zscore of TWAS")+
      geom_hline(yintercept=0, linetype="dashed", color="grey30")+
      geom_vline(xintercept=0, linetype="dashed", color="grey30")+
      ggtitle(oneMCl)+ 
      theme_bw()+
      theme(legend.position="none", plot.title=element_text(hjust=0.5)) 
    ###
    p0
})

figfn <- paste(outdir, "Figure1.2_twas_DVG.png", sep="")
png(figfn, width=720, height=720, res=100)
plot_grid(plotlist=fig_ls)
dev.off()




#############
### eQTLs ###
#############

###
### mashr egenes
load("mashr-eQTLs-mlfsr0.1.Rd")
conditions <- names(eqtl)
egene_DF <- map_dfr(conditions, function(ii){
   ## 
   egene <- unique(gsub("\\..{1,2}_.*", "", eqtl[[ii]]))
   tmp <- data.frame(egene=egene, conditions=gsub("-EtOH", "", ii))
   tmp
})


###
### FastQTL egenes
conditions <- names(eqtl)
egene_DF <- map_dfr(conditions, function(ii){
   ###
   fn <- paste("../eQTL/egenes/", ii, "_egenes.txt", sep="")
   egene <- read.table(fn)$V1
   egene <- unique(gsub("\\..*", "", egene))
   tmp <- data.frame(egene=egene, conditions=gsub("-EtOH", "", ii))
   tmp
})    


load("mashr_input.Rd")
eQTL <- slopes/SEs
colnames(eQTL) <- gsub("-EtOH", "", gsub("\\.", "-",  colnames(eQTL)))

eQTL2 <- eQTL%>%as.data.frame()%>%
   mutate(rn=rownames(eQTL), Gene=gsub("\\..{1,2}_.*", "", rn))%>%
   filter(Gene%in%unique(res3$Gene))

eQTL2_DF <- eQTL2%>%dplyr::select(-Gene)%>%
    pivot_longer(!rn, names_to="conditions", values_to="zscore_eQTL")%>%as.data.frame()%>%
    mutate(Gene=gsub("\\..{1,2}_.*", "", rn))
##
conditions <- sort(unique(eQTL2_DF$conditions))
eQTL2_DF <- map_dfr(conditions, function(ii){
   ##
   egene <- egene_DF%>%filter(conditions==ii)%>%pull(egene)%>%unique()
   ##
   tmp <- eQTL2_DF%>%filter(conditions==ii)%>%mutate(is_egene=ifelse(Gene%in%egene, 1, 0))
   tmp
})    

eQTL2_DF <- eQTL2_DF%>%group_by(conditions, Gene)%>%slice_max(order_by=abs(zscore_eQTL), n=1)%>%as.data.frame()
eQTL2_DF <- eQTL2_DF%>%
    mutate(MCls=gsub("_.*", "", conditions), treats=gsub(".*_", "", conditions))


plotDF <- eQTL2_DF%>%inner_join(res3, by="Gene")%>%mutate(gr2=ifelse(is_egene==0, "0", treats))


tmp <- plotDF%>%filter(gr2>0)%>%dplyr::select(conditions, MCls, treats, SYMBOL, zscore_TWAS, zscore_eQTL)
tmp%>%filter(MCls=="Bcell")

## x <- plotDF%>%group_by(conditions)%>%
##     summarise(n=sum(is_egene), .groups="drop")%>%
##     mutate(condition2=gsub("-EtOH", "", conditions),
##            MCls=gsub("_.*", "", condition2), treats=gsub(".*_", "", condition2))%>%
##     arrange(condition2)%>% 
##     pivot_wider(id_cols="MCls", names_from="treats", values_from="n")
 
MCls <- sort(unique(plotDF$MCls))
fig_ls <- lapply(MCls, function(oneMCl){
   ###
   plotDF2 <- plotDF%>%filter(MCls==oneMCl)
   anno2 <- plotDF2%>%filter(gr2>0)%>%
       mutate(sign_eqtl=sign(zscore_eQTL))%>%
       dplyr::distinct(Gene, sign_eqtl, .keep_all=T)
   ### 
   p0 <- ggplot(plotDF2, aes(x=zscore_eQTL, y=zscore_TWAS))+
      geom_point(aes(color=factor(gr2), shape=factor(gr2), size=factor(gr2)))+
      geom_text_repel(data=anno2, mapping=aes(x=zscore_eQTL, y=zscore_TWAS, label=SYMBOL, color=gr2),
                      size=3, max.overlaps=Inf)+         
      scale_color_manual(values=c("0"="grey80", "CTRL"="grey10", "LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
                                  "PHA"="#a6cee3", "PHA-DEX"="#1f78b4"))+
      scale_shape_manual(values=c("0"=4, "CTRL"=2, "LPS"=2, "LPS-DEX"=2, "PHA"=2, "PHA-DEX"=2))+
      scale_size_manual(values=c("0"=1, "CTRL"=2.5, "LPS"=2.5, "LPS-DEX"=2.5, "PHA"=2.5, "PHA-DEX"=2.5))+
      xlab("z-score of eQTL")+
      ylab("z score of TWAS")+
      ggtitle(oneMCl)+ 
      geom_hline(yintercept=0, linetype="dashed", color="grey30")+
      geom_vline(xintercept=0, linetype="dashed", color="grey30")+ 
      theme_bw()+
      theme(legend.position="none", plot.title=element_text(hjust=0.5))
   ##
   p0
})    

figfn <- paste(outdir, "Figure2.1_twas_eQTL_mashR.png", sep="")
png(figfn, width=720, height=720, res=100)
plot_grid(plotlist=fig_ls, ncol=2)
dev.off()


###
### eQTLs in T-cell

plotDF2 <- eQTL2_DF%>%inner_join(res3, by="Gene")%>%
    filter(is_egene>0,  MCls=="Tcell", grepl("DEX", treats))%>%
    mutate(gr_fill=ifelse(treats=="PHA-DEX", MCls, "0"))

##
p <- ggplot(plotDF2, aes(x=zscore_eQTL, y=zscore_TWAS, color=factor(treats)))+
      geom_point(shape=17, size=2.5)+
      geom_text_repel(data=plotDF2, mapping=aes(x=zscore_eQTL, y=zscore_TWAS, label=SYMBOL, color=factor(treats)),
                      size=2.5, max.overlaps=10)+
     scale_color_manual(values=c("LPS-DEX"="#e31a1c", "PHA-DEX"="#1f78b4"),
                        labels=c("LPS-DEX"="LPS+DEX", "PHA-DEX"="PHA+DEX"),
                        guide=guide_legend(override.aes=list(shape=17, size=3)))+   
      ## scale_color_manual(values=c("Bcell"="#4daf4a", "Monocyte"="#984ea3",
      ##    "NKcell"="#aa4b56", "Tcell"="#ffaa00"), guide="none")+
      ## scale_fill_manual(values=c("0"="NA", "Bcell"="#4daf4a", "Monocyte"="#984ea3",
      ##    "NKcell"="#aa4b56", "Tcell"="#ffaa00"),
      ##    breaks=c("0", "Tcell"),
      ##    labels=c("0"="LPS+DEX", "Bcell"="PHA+DEX"),
      ##    guide=guide_legend(override.aes=list(color="#ffaa00")))+
      xlab("zscore of eQTL")+
      ylab("zscore of TWAS")+
      geom_hline(yintercept=0, linetype="dashed", color="grey30")+
      geom_vline(xintercept=0, linetype="dashed", color="grey30")+
      ggtitle("T-cells")+
      theme_bw()+
      theme(## legend.position=c(0.25, 0.85),
            legend.title=element_blank(),
            legend.background=element_blank(),
            legend.box.background=element_blank(),
            plot.title=element_text(hjust=0.5))
###
figfn <- paste(outdir, "Figure2.3_twas_egene_Tcell.png", sep="")
png(figfn, width=520, height=420, res=120)
print(p)
dev.off()










######################
### response eQTLs ###
######################


###
### response eGenes
reQTL <- read_rds("./2_reQTL.output/reQTLs_All.infor.rds")%>%
    mutate(condition2=gsub("-EtOH", "", gsub(":.*", "", condition)))

## reGene <- unique(reQTL$ENSG2)


###
### mashr reGene
zdiff <- read_rds("./2_reQTL.output/reQTLs_mashr_Zdiff.rds")
colnames(zdiff) <- gsub("-EtOH", "", gsub(":.*", "", colnames(zdiff)))
###

zdiff2 <- zdiff%>%as.data.frame()%>%
   mutate(rn=rownames(zdiff), Gene=gsub("\\..{1,2}_.*", "", rn))%>%
   filter(Gene%in%unique(res3$Gene))

          ## varID=gsub(".*\\..{1,2}_", "", rn))
zdf <- zdiff2%>%dplyr::select(-Gene)%>%
    pivot_longer(!rn, names_to="conditions", values_to="zscore_reQTL")%>%as.data.frame()%>%
    mutate(Gene=gsub("\\..{1,2}_.*", "", rn))


##
conditions <- sort(unique(zdf$conditions))
zdf <- map_dfr(conditions, function(ii){
   ##
   reQTL2 <- reQTL%>%filter(condition2==ii)
   regene <- unique(reQTL2$ENSG2)
   ##
   tmp <- zdf%>%filter(conditions==ii)%>%mutate(is_reGene=ifelse(Gene%in%regene, 1, 0))
   tmp
})    

zdf2 <- zdf%>%group_by(conditions, Gene)%>%slice_max(order_by=abs(zscore_reQTL), n=1)%>%as.data.frame()
zdf2 <- zdf2%>%mutate(MCls=gsub("_.*", "", conditions), contrast=gsub(".*_", "", conditions)) 
plotDF <- zdf2%>%inner_join(res3, by="Gene")%>%
    mutate(gr2=ifelse(is_reGene==1, contrast, 0))

x <- plotDF%>%group_by(conditions)%>%
    summarise(n=sum(is_reGene), .groups="drop")%>%
    mutate(MCls=gsub("_.*", "", conditions), treats=gsub(".*_", "", conditions))%>%
    arrange(conditions)%>%
    pivot_wider(id_cols="MCls", names_from="treats", values_from="n")
 
MCls <- sort(unique(plotDF$MCls))
fig_ls <- lapply(MCls, function(oneMCl){
   ###
   plotDF2 <- plotDF%>%filter(MCls==oneMCl)
   anno2 <- plotDF2%>%filter(gr2>0)
   ### 
   p0 <- ggplot(plotDF2, aes(x=zscore_reQTL, y=zscore_TWAS, color=factor(gr2)))+
      geom_point(aes(shape=gr2, size=gr2))+
      geom_text_repel(data=anno2, mapping=aes(x=zscore_reQTL, y=zscore_TWAS, label=SYMBOL, color=gr2),
                      size=3, max.overlaps=Inf)+         
      scale_color_manual(values=c("0"="grey", "LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
                                  "PHA"="#a6cee3", "PHA-DEX"="#1f78b4"))+
      scale_shape_manual(values=c("0"=4, "LPS"=2, "LPS-DEX"=2, "PHA"=2, "PHA-DEX"=2))+
      scale_size_manual(values=c("0"=1, "LPS"=2.5, "LPS-DEX"=2.5, "PHA"=2.5, "PHA-DEX"=2.5))+
      xlab("z-score of reQTL")+
      ylab("z score of TWAS")+
      ggtitle(oneMCl)+ 
      geom_hline(yintercept=0, linetype="dashed", color="grey30")+
      geom_vline(xintercept=0, linetype="dashed", color="grey30")+ 
      theme_bw()+
      theme(legend.position="none", plot.title=element_text(hjust=0.5))
   ##
   p0
})    

figfn <- paste(outdir, "Figure3.1_twas_reQTL_mashr.png", sep="")
png(figfn, width=720, height=720, res=100)
plot_grid(plotlist=fig_ls, ncol=2)
dev.off()


###
### FastQTL
zdiff <- read_rds("./2_reQTL.output/reQTLs_FastQTL_Zdiff.rds")
colnames(zdiff) <- gsub("-EtOH", "", gsub(":.*", "", colnames(zdiff)))
###

zdiff2 <- zdiff%>%as.data.frame()%>%
   mutate(rn=rownames(zdiff), Gene=gsub("\\..{1,2}_.*", "", rn))%>%
   filter(Gene%in%unique(res3$Gene))

          ## varID=gsub(".*\\..{1,2}_", "", rn))
zdf <- zdiff2%>%dplyr::select(-Gene)%>%
    pivot_longer(!rn, names_to="conditions", values_to="zscore_reQTL")%>%as.data.frame()
zdf <- zdf%>%mutate(Gene=gsub("\\..{1,2}_.*", "", rn))    
zdf2 <- zdf%>%group_by(conditions, Gene)%>%slice_max(order_by=abs(zscore_reQTL), n=1)%>%as.data.frame()
zdf2 <- zdf2%>%mutate(MCls=gsub("_.*", "", conditions), contrast=gsub(".*_", "", conditions)) 

plotDF <- zdf2%>%inner_join(res3, by="Gene")%>%mutate(gr2=ifelse(Gene%in%reGene, contrast, 0))

MCls <- sort(unique(plotDF$MCls))
fig_ls <- lapply(MCls, function(oneMCl){
   ###
   plotDF2 <- plotDF%>%filter(MCls==oneMCl)
   anno2 <- plotDF2%>%filter(gr2>0)
   ## 
   p0 <- ggplot(plotDF, aes(x=zscore_reQTL, y=zscore_TWAS, color=factor(gr2)))+
      geom_point(aes(shape=gr2, size=gr2))+
      geom_text_repel(data=anno2, mapping=aes(x=zscore_reQTL, y=zscore_TWAS, label=SYMBOL, color=gr2),
                      size=3, max.overlaps=Inf)+   
      scale_color_manual(values=c("0"="grey", "LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
                                  "PHA"="#a6cee3", "PHA-DEX"="#1f78b4"))+
      scale_shape_manual(values=c("0"=4, "LPS"=2, "LPS-DEX"=2, "PHA"=2, "PHA-DEX"=2))+
      scale_size_manual(values=c("0"=1, "LPS"=2.5, "LPS-DEX"=2.5, "PHA"=2.5, "PHA-DEX"=2.5))+
      xlab("z-score of reQTL")+
      ylab("z score of TWAS")+
      ggtitle(oneMCl)+ 
      geom_hline(yintercept=0, linetype="dashed", color="grey30")+
      geom_vline(xintercept=0, linetype="dashed", color="grey30")+ 
      theme_bw()+
      theme(legend.position="none", plot.title=element_text(hjust=0.5))
   ##
   p0
})    



figfn <- paste(outdir, "Figure3.2_twas_reQTL_FastQTL.png", sep="")
png(figfn, width=720, height=720, res=100)
plot_grid(plotlist=fig_ls, ncol=2)
dev.off()



#############
#### vQTL ###
#############

### mashr
load("../mashr_dispersionQTL/mashr-eQTLs-mlfsr0.1.Rd")
conditions <- names(eqtl)
egene_DF <- map_dfr(conditions, function(ii){
   ## 
   egene <- unique(gsub("\\..{1,2}_.*", "", eqtl[[ii]]))
   tmp <- data.frame(egene=egene, conditions=gsub("-EtOH", "", ii))
   tmp
})


load("../mashr_dispersionQTL/mashr_input.Rd")

eQTL <- slopes/SEs
colnames(eQTL) <- gsub("-EtOH", "", gsub("\\.", "-", colnames(eQTL)))

eQTL2 <- eQTL%>%as.data.frame()%>%
   mutate(rn=rownames(eQTL), Gene=gsub("\\..{1,2}_.*", "", rn))

gene_overlap <- intersect(eQTL2$Gene, res3$Gene)
eQTL2 <- eQTL2%>%filter(Gene%in%gene_overlap)

eQTL2_DF <- eQTL2%>%dplyr::select(-Gene)%>%
    pivot_longer(!rn, names_to="conditions", values_to="zscore_eQTL")%>%as.data.frame()%>%
    mutate(Gene=gsub("\\..{1,2}_.*", "", rn))%>%drop_na(zscore_eQTL)

conditions <- sort(unique(eQTL2_DF$conditions))
eQTL2_DF <- map_dfr(conditions, function(ii){
    ##
    egene <- egene_DF%>%filter(conditions==ii)%>%pull(egene)
    if ( length(egene)>0){
       ##
       tmp <- eQTL2_DF%>%filter(conditions==ii)%>%mutate(is_egene=ifelse(Gene%in%egene, 1, 0))
    }else{
       tmp <- eQTL2_DF%>%filter(conditions==ii)
       tmp$is_egene <- 0
    }
    tmp
})    
###
###

eQTL2_DF <- eQTL2_DF%>%group_by(conditions, Gene)%>%slice_max(order_by=abs(zscore_eQTL), n=1)%>%as.data.frame()
  
eQTL2_DF <- eQTL2_DF%>%
    mutate(MCls=gsub("_.*", "", conditions), treats=gsub(".*_", "", conditions),
           gr2=ifelse(is_egene==0, "0", treats) )
 
plotDF <- eQTL2_DF%>%inner_join(res3, by="Gene")

## x <- plotDF%>%group_by(conditions)%>%
##     summarise(n=sum(is_egene), .groups="drop")%>%
##     mutate(MCls=gsub("_.*", "", conditions), treats=gsub(".*_", "", conditions))%>%
##     arrange(conditions)%>%
##     pivot_wider(id_cols="MCls", names_from="treats", values_from="n")

MCls <- sort(unique(plotDF$MCls))
fig_ls <- lapply(MCls, function(oneMCl){
   ###
   plotDF2 <- plotDF%>%filter(MCls==oneMCl)
   anno2 <- plotDF2%>%filter(gr2>0)
   ### 
   p0 <- ggplot(plotDF2, aes(x=zscore_eQTL, y=zscore_TWAS))+
      geom_point(aes(color=factor(gr2), shape=factor(gr2), size=factor(gr2)))+
      geom_text_repel(data=anno2, mapping=aes(x=zscore_eQTL, y=zscore_TWAS, label=SYMBOL, color=gr2),
                      size=3, max.overlaps=Inf)+    
      scale_color_manual(values=c("0"="grey80", "CTRL"="grey10", "LPS"="#fb9a99", "LPS-DEX"="#e31a1c",
                                  "PHA"="#a6cee3", "PHA-DEX"="#1f78b4"))+
      scale_shape_manual(values=c("0"=4, "CTRL"=2, "LPS"=2, "LPS-DEX"=2, "PHA"=2, "PHA-DEX"=2))+
      scale_size_manual(values=c("0"=1, "CTRL"=2.5, "LPS"=2.5, "LPS-DEX"=2.5, "PHA"=2.5, "PHA-DEX"=2.5))+
      xlab("z-score of vQTL")+
      ylab("z score of TWAS")+
      ggtitle(oneMCl)+ 
      geom_hline(yintercept=0, linetype="dashed", color="grey30")+
      geom_vline(xintercept=0, linetype="dashed", color="grey30")+ 
      theme_bw()+
      theme(legend.position="none", plot.title=element_text(hjust=0.5))
   ##
   p0
})    

figfn <- paste(outdir, "Figure4.1_twas_vQTL.png", sep="")
png(figfn, width=720, height=720, res=100)
plot_grid(plotlist=fig_ls, ncol=2)
dev.off()



###
### dynamical 

## fn <- "../LDA-eQTL_Julong/DiagLDA2/4.2_lm.results/zzz_dynamical.eQTL.rds"
## dqtl <- read_rds(fn)
## ##
## conditions <- sort(unique(dqtl$condition))
## ##
## summ2 <- map_dfr(conditions, function(ii){
##     ##
##     ENSG <- dqtl%>%filter(condition==ii)%>%dplyr::pull(ENSG)
##     ENSG <- unique(gsub("\\..*", "", ENSG))
##     tmp <- data.frame(conditions=ii, noverlap=length(intersect(ENSG, res3$Gene)))
##     tmp
## })

## ENSG <- unique(gsub("\\..*", "", dqtl$ENSG))
## res3%>%filter(Gene%in%ENSG)

