###########################
# Kali Pruss
# Pruss and Sonnenburg, Nature 2021

# Code for in vitro transcriptional profiling (Extended Data Fig. 8c)

###########################

## Input: variance stabilizing-transformed counts (vsd), see `in vitro RNA-seq.R`

library(ggplot2); library(DESeq2)

pcaData <- plotPCA(vsd, intgroup=c("group"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=group)) +
  geom_point(size=4, alpha = 0.8) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + theme_bw() +
  scale_color_manual(values= c("#9F9D95","#F5DA1C", "#8481E7" , "#04A90B"))
