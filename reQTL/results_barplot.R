# this script makes publication-ready barplots of numbers of significant eQTLs
# based on /wsu/home/groups/piquelab/SCAIP/eQTL/FastQTL/nominals/reQTL/lm_SCAIP1-6/results_barplot.R 
# 1/22/2021 JR

library(ggplot2)

# colors:
MCls <- c("Bcell", "Monocyte", "NKcell", "Tcell")
col2 <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3", "NKcell"="#aa4b56", "Tcell"="#ffaa00")  #T color "#ff9400" #NK color, "#a63728"
contrast <- c("LPS", "LPS-DEX", "PHA", "PHA-DEX")
col1 <- c("LPS"="#fb9a99", "LPS+DEX"="#e31a1c", "PHA"="#a6cee3", "PHA+DEX"="#1f78b4")
colors <- c(col1,col2)

# read in the results:
res <- read.table(paste0("ANOVA_signif_interactions_corrected_lm_100.txt"), sep="\t", stringsAsFactors=F)
colnames(res) <- c("condition","reQTLs","reGenes")

# rename treatments:
res$condition <- gsub("-DEX","+DEX",res$condition)
res$condition <- gsub("-EtOH","",res$condition)

# add cell type and trt for coloring and faceting:
res$cell <- gsub("_.*","",res$condition)
res$ctrl <- gsub(".*_","",res$condition)
res$trt <- gsub(".*cell_|.*cyte_","",res$condition)
res$trt <- gsub("_.*","",res$trt)

# drop the empties:
res <- res[!res$reQTLs==0,]

# plot:
## png(paste0("plots/eQTL_results_barplot.png"),res=120)# width=800, height=400, pointsize=12)
fig0 <- ggplot(res, aes(x=cell, y=reQTLs, fill=cell))+
            geom_bar(stat="identity")+
            scale_fill_manual(values=colors,labels="")+
    ## ylab(paste0())+
    ##         geom_text(data=ann2, aes(x=MCls, label=ngene, y=ngene+100, fill=NULL), size=3)+
            facet_grid(~trt)+
            theme_bw()+
    ylab(paste0("reQTLs"))+
            theme(legend.position="none",
               axis.title.x=element_blank(),
               axis.text.x=element_text(angle=-90,hjust=0, vjust=0.5))
ggsave(paste0("plots/reQTL_results_barplot.png"),fig0)

## png(paste0("plots/eGene_results_barplot.png"), res=120)
fig1 <- ggplot(res, aes(x=cell, y=reGenes, fill=cell))+
            geom_bar(stat="identity")+
            scale_fill_manual(values=colors,labels="")+
    ## ylab(paste0())+
    ##         geom_text(data=ann2, aes(x=MCls, label=ngene, y=ngene+100, fill=NULL), size=3)+
            facet_grid(~trt)+
    ylab(paste0("reGenes"))+
            theme_bw()+
            theme(legend.position="none",
               axis.title.x=element_blank(),
               axis.text.x=element_text(angle=-90,hjust=0, vjust=0.5))
ggsave(paste0("plots/reGene_results_barplot.png"),fig1)

### END 1/22/2021
