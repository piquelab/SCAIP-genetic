# this script makes QQplots for all interaction corrected p-values by cell type
# based on /wsu/home/groups/piquelab/SCAIP/eQTL/FastQTL/nominals/reQTL/lm_SCAIP1-6/QQplot_all_by-cell-type.R 
# 1/22/2021 JR

library(data.table)
library(dplyr)
library(qvalue)
library(ggplot2)
library(ggrastr)

## pthresh <- 0.1

# colors for plotting:
scaip_colors <- list(cell=c("Bcell"="#4daf4a", "Monocyte"="#984ea3","NKcell"="#AA4B56", "Tcell"="#FFAA00"), condition=c("CTRL"="#828282","LPS"="#fb9a99", "LPS+DEX"="#e31a1c","PHA"="#a6cee3", "PHA+DEX"="#1f78b4"), Batch=c("SCAIP1"="#e41a1c", "SCAIP2"="#377eb8","SCAIP3"="#4daf4a", "SCAIP4"="#984ea3", "SCAIP5"="#ff7f00", "SCAIP6"="#ffff33")) #brewer.pal(4,"Set1")

# load all the corrected pvalues:
files <- list.files(path=paste0("reQTL_lm_results"),pattern=".*_corrected_100.txt",full.names=T)
conditions <- gsub("_corrected_100.txt","",list.files(path=paste0("reQTL_lm_results"),pattern=".*_corrected_100.txt"))
conditions <- gsub("results","",conditions)
ldf <- lapply(files, fread)

## cors <- data.frame(geneSNP=character())
## for(i in 1: length(ldf)){
## ## ldf[[i]]$condition <- as.character(conditions[i])
## ldf[[i]]$geneSNP <- paste0(gsub("[.].*","",ldf[[i]]$ENSG),"_",ldf[[i]]$varID)
## colnames(ldf[[i]])[15] <- paste0(strsplit(conditions[i],"_")[[1]][1],"_",strsplit(conditions[i],"_")[[1]][2])
## cors <- full_join(cors,ldf[[i]][,c(15,18)])
## }

# merge :
corr <- data.frame()
for(i in 1: length(ldf)){
ldf[[i]]$condition <- as.character(conditions[i])
ldf[[i]]$treatment <- strsplit(ldf[[i]]$condition,"_")[[i]][2]
corr <- rbind(corr,ldf[[i]][,c(15,17:18)])
}

# add cell and trt columns:
corr$cell <- gsub("_.*","",corr$condition)
# rename pcorr:
colnames(corr)[1] <- "pcorr"

# fix the trt names:
corr$treatment <- gsub("-DEX","+DEX",corr$treatment)
corr$treatment <- gsub("-EtOH","",corr$treatment)

# calculate expected etc:
res2 <- corr %>% filter(!is.na(pcorr)) %>%
        arrange(pcorr) %>%
        group_by(treatment,cell) %>%
        mutate(r=rank(pcorr, ties.method = "random"),pexp=r/length(pcorr))

# QQplot:
p1 <- res2 %>%
        ggplot(aes(x=-log10(pexp),y=-log10(pcorr),color=cell)) +
        geom_point() +
        scale_color_manual(values=scaip_colors$cell) +
        guides(colour = guide_legend(override.aes = list(size=5),title="Cell type")) +
        geom_abline(slope=1,intercept=0) +
        facet_grid(~treatment) +
        xlab(expression(Expected -log[10](p))) +
        ylab(expression(Observed -log[10](p))) +
        theme_bw()
ggsave("plots/QQplots/QQ_trt_cell-color.png",p1)
pdf("plots/QQplots/QQ_trt_cell-color.pdf",height=6, width=9)
p1
dev.off()

p2 <- res2 %>%
        ggplot(aes(x=-log10(pexp),y=-log10(pcorr),color=treatment)) +
        geom_point() +
        scale_color_manual(values=scaip_colors$condition) +
        guides(colour = guide_legend(override.aes = list(size=5),title="Treatment")) +
        geom_abline(slope=1,intercept=0) +
        facet_grid(~cell) +
        xlab(expression(Expected -log[10](p))) +
        ylab(expression(Observed -log[10](p))) +
        theme_bw()
ggsave("plots/QQplots/QQ_cell_trt-color.png",p2)
pdf("plots/QQplots/QQ_cell_trt-color.pdf",height=6, width=9)
p2
dev.off()


### END 1/22/2021
