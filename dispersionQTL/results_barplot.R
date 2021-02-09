# this script makes publication-ready barplots of numbers of significant eQTLs
# based on /wsu/home/groups/piquelab/SCAIP/varianceQTL/FastQTL_on_residuals/nominals/analysis/results_barplot.R
# 1/25/2021 JR

library(ggplot2)

# colors:
scaip_colors <- list(cell=c("Bcell"="#4daf4a", "Monocyte"="#984ea3","NKcell"="#AA4B56", "Tcell"="#FFAA00"), treatment=c("CTRL"="#828282","LPS"="#fb9a99", "LPS+DEX"="#e31a1c","PHA"="#a6cee3", "PHA+DEX"="#1f78b4"), Batch=c("SCAIP1"="#e41a1c", "SCAIP2"="#377eb8","SCAIP3"="#4daf4a", "SCAIP4"="#984ea3", "SCAIP5"="#ff7f00", "SCAIP6"="#ffff33")) #brewer.pal(4,"Set1")

# read in the results:
res <- read.table(paste0("dispersion_nominals_eqtls-egenes_summary.txt"), sep="\t", stringsAsFactors=F)
colnames(res) <- c("condition","PCs","veQTLs","veGenes")

# subset to best number of PCs - 3:
res <- res[res$PCs==3,]

# subset to only the conditions you want to plot:
res <- res[grep("Tcell",res$condition),]

# rename treatments:
res$condition <- gsub("-DEX","+DEX",res$condition)
res$condition <- gsub("-EtOH","",res$condition)

# add cell type and trt for coloring and faceting:
res$cell <- gsub("_.*","",res$condition)
res$ctrl <- gsub(".*_","",res$condition)
res$trt <- gsub(".*cell_","",res$condition)
res$trt <- gsub("_.*","",res$trt)

# drop the empties:
res <- res[!res$veQTLs==0,]

# plot:
pdf("plots/dispersion_eQTL_results_barplot.pdf")#,res=120)# width=800, height=400, pointsize=12)
fig0 <- ggplot(res, aes(x=trt, y=veQTLs, fill=trt))+
            geom_bar(stat="identity")+
            scale_fill_manual(values=scaip_colors$treatment,labels="")+
    ## ylab(paste0())+
    ##         geom_text(data=ann2, aes(x=MCls, label=ngene, y=ngene+100, fill=NULL), size=3)+
            ## facet_grid(~trt)+
            theme_bw()+
    ylab(paste0("dispersion eQTLs"))+
            theme(legend.position="none",
               axis.title.x=element_blank(),
               axis.text.x=element_text(angle=-90,hjust=0, vjust=0.5))
fig0
dev.off()
ggsave("plots/dispersion_eQTL_results_barplot.png",fig0)

pdf("plots/dispersion_eGene_results_barplot.pdf")#, res=120)
fig0 <- ggplot(res, aes(x=trt, y=veGenes, fill=trt))+
            geom_bar(stat="identity")+
            scale_fill_manual(values=scaip_colors$treatment,labels="")+
    ## ylab(paste0())+
    ##         geom_text(data=ann2, aes(x=MCls, label=ngene, y=ngene+100, fill=NULL), size=3)+
            ## facet_grid(~trt)+
    ylab(paste0("dispersion eGenes"))+
            theme_bw()+
            theme(legend.position="none",
               axis.title.x=element_blank(),
               axis.text.x=element_text(angle=-90,hjust=0, vjust=0.5))
fig0
dev.off()
ggsave("plots/dispersion_eGene_results_barplot.png",fig0)
 
### END 1/25/2021
