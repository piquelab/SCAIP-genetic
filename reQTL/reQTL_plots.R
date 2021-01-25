# this script makes reQTL plots for reQTLs based on lm model-ANOVA
# based on ../../eQTL/FastQTL/nominals/reQTL/lm_SCAIP1-6//reQTL_plots.R
# 1/25/2021 JR

library(dplyr)
library(data.table)
library(ggplot2)

cl <- "Tcell"
ctrl <- "_CTRL"
trt <- "_LPS-EtOH"

args = commandArgs(trailingOnly=TRUE)

# get the PC number:
if (length(args)>0){
      ## pc <- args[1]
      cl <- args[1]
      trt <- args[2]
      ctrl <- args[3]
      }

contrast <- paste0(cl, trt, ctrl)

# set FDR threshold;
FDR <- 0.1

# load the needed files:
# results:
res <- read.table(paste0("reQTL_lm_results/",cl, trt, ctrl, "_corrected_100.txt"), sep="\t", header=T, stringsAsFactors=FALSE)
# treatment GE:
treat <- fread(paste0("zcat ../eQTL/normalized_GE_residuals/",cl, trt,".bed.gz"), header=T, sep="\t", stringsAsFactors=FALSE)
# control GE:
contr <- fread(paste0("zcat ../eQTL/normalized_GE_residuals/",cl, ctrl,".bed.gz"), header=T, sep="\t", stringsAsFactors=FALSE)
# 2. txt file with dosages:
dosages <- read.table(paste0("../eQTL/eQTL_coordinates/all_uniq_eQTL_dosages.txt"), sep="\t", header=T,stringsAsFactors=F)
colnames(dosages) <- gsub("[.]","-",colnames(dosages))
# remove repeated dosages (where do they come from, though??):
dosages <- dosages[!duplicated(dosages$varID),]
rownames(dosages) <- dosages$varID
# 3. txt file with list of testable pairs:
pairs <- read.table(paste0("../eQTL/eQTL_coordinates/all_signif_SNP-gene_pairs.txt"),sep="\t", header=T,stringsAsFactors=F)
# remove the genes absent from eith ctrl or trt GE files:
test_pairs <- pairs[pairs$ENSG %in% treat$ID & pairs$ENSG %in% contr$ID,]
# repeat the dosages as in the pairs file:
dos <- dosages[test_pairs$varID,5:ncol(dosages)]

# ensg-to-gene-name key:
library(annotables)
gene_key <- data.frame(grch38[!duplicated(grch38$ensgene),c(1,3)])
rownames(gene_key) <- gene_key$ensgene

### prepare GE files:
GE_c <- data.frame(contr[,5:ncol(contr)])
colnames(GE_c) <- gsub("[.]","-",colnames(GE_c))
rownames(GE_c) <- contr$ID
GE_t <- data.frame(treat[,5:ncol(treat)])
colnames(GE_t) <- gsub("[.]","-",colnames(GE_t))
rownames(GE_t) <- treat$ID

# subset to common individuals only:
common_dbgaps <- colnames(GE_c)[colnames(GE_c) %in% colnames(GE_t)]
expression_c <- GE_c[,common_dbgaps]
expression_t <- GE_t[,common_dbgaps]
dos <- dosages[,common_dbgaps]

# subset to reQTLs you want to plot:
sig <- filter(res,interaction_A_lm_pcorr.qval<FDR)
expression_c <- expression_c[sig$ENSG,]
expression_t <- expression_t[sig$ENSG,]
dos <- dos[sig$varID,]


#Plot
# one plot per significant pair:
for (i in 1:nrow(dos)){

# make everything into a data frame:
df1 <- data.frame(cbind(t(dos)[,i,drop=F], t(expression_c[i,])))
df1$facet <- gsub("_","",ctrl)
df2 <- data.frame(cbind(t(dos)[,i,drop=F], t(expression_t[i,])))
df2$facet <- gsub("_","",trt)
df <- rbind(df1,df2)
colnames(df) <- c("Dosage", "Expression", "facet")
df$Dosage <- as.numeric(df$Dosage)
df$Dosage_categorical <- as.factor(round(df$Dosage))
# grab gene and variant names for labels:
gene <- gsub("^(.*?[.].*?)[.].*", "\\1",rownames(expression_c)[i])
variant <- gsub("[.].*","",rownames(dos)[i])
# get the gene name:
gene_name <- gene_key[gsub("[.].*","",gene),"symbol"]

# fix the trt names:
df$facet <- gsub("-EtOH","",df$facet)
df$facet <- gsub("-DEX","+DEX",df$facet)
    
# add genotype:
ref <- strsplit(variant,":")[[1]][3]
alt <- gsub(";.*","",strsplit(variant,":")[[1]][4])
df[df$Dosage_categorical=='0',"genotype"] <- paste0(ref,"/",ref)
df[df$Dosage_categorical=='1',"genotype"] <- paste0(ref,"/",alt)
df[df$Dosage_categorical=='2',"genotype"] <- paste0(alt,"/",alt)
    
# relevel the genotype, so it always goes: refref, refalt, altalt:
df$genotype <- factor(df$genotype,levels=c(paste0(ref,"/",ref),paste0(ref,"/",alt),paste0(alt,"/",alt)))

## # make a scatterplot:
## # get the regression coefficients for the gene-SNP pair:
##         model <- lm(df$Expression~df$Dosage*df$facet)
##         Intercept <- summary(model)$coefficients[1,1]
##         dosage_beta <- summary(model)$coefficients[2,1]
##         trt_beta <- summary(model)$coefficients[3,1]
##         interaction_beta  <- summary(model)$coefficients[4,1]
##         in0 <- Intercept
##         in1 <- Intercept+dosage_beta
##         in2 <- Intercept+2*dosage_beta
##         slop0 <- trt_beta
##         slop1 <- trt_beta+interaction_beta
##         slop2 <- trt_beta+2*interaction_beta
##         df[df$Dosage_categorical==1,"alpha_dosage"] <- 0.35
##         df[df$Dosage_categorical==2,"alpha_dosage"] <- 0.7
##         df[df$Dosage_categorical==3,"alpha_dosage"] <- 1

##         pdf(paste0("./plots/eQTLs/scatterplots/", cl, ctrl,"_vs",trt,"_scatterplot_", variant,"_", gene, ".pdf"))
##          p <- ggplot(df, aes(genotype, Expression, color=Dosage_categorical,alpha=alpha_dosage))
##          print(p + geom_point() + labs(x="", y=paste0(gene_name, " gene expresssion")) +
##                             # fit a smoothing spline to each genotype class:
##                             geom_smooth(alpha=0.7)+
##                             ## geom_abline(intercept=in1,slope=slop1,color="#984ea3", size=1)+
##                             ## geom_abline(intercept=in2,slope=slop2,color="#ff7f00", size=1)+
##                             theme_bw(base_size = 11, base_family = "") +
##                             scale_colour_manual(values= c("0"="#4daf4a", "1"="#984ea3","2"="#ff7f00")) + theme(axis.title.y = element_text(size = 18),axis.title.x = element_text(size = 18),legend.position="none"))
##         dev.off()
## ggsave(paste0("./plots/eQTLs/scatterplots/", cl, ctrl,"_vs",trt,"_scatterplot_", variant,"_", gene, ".png"),p)
    
# make a boxplot:
pdf(paste0("./plots/eQTLs/boxplots/boxplot_", cl, ctrl,"_vs",trt,"_",gene,"_",variant,"_boxplot_notched_colored.pdf"))
p<-ggplot(data=df, aes(x=as.factor(genotype), fill=facet, y=Expression))  +xlab("")+ ylab(paste0(gene_name, " normalized gene expression")) + geom_boxplot(outlier.shape=NA) + geom_smooth(method='lm')+ geom_jitter(width = 0.25)+ theme_bw()+facet_grid(~facet) +  scale_fill_manual(values= c("CTRL"="#828282", "LPS"="#fb9a99", "LPS+DEX"="#e31a1c","PHA"="#a6cee3", "PHA+DEX"="#1f78b4"))+ theme(legend.position="none") + theme(axis.title.y = element_text(size = 14),strip.text.x = element_text(size = 16))
print(p)
dev.off()
ggsave(paste0("./plots/eQTLs/boxplots/boxplot_", cl, ctrl,"_vs",trt,"_",gene,"_",variant,"_boxplot_notched_colored.png"),p)    
}

sessionInfo()

## END 1/25/2021
