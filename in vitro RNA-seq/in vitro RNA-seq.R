###########################
# Kali Pruss
# Pruss and Sonnenburg, Nature 2021

# Code for in vitro transcriptional profiling differential gene expression (Extended Data Fig. 8, Supplementary Table 3)

###########################

## Input: takes raw counts (exp.all) table with rRNA reads removed where:
# rows are genes, with unique row IDs (gene IDs) and
# columns are samples, with unique column names that match sample names in metadata (library_trmts)

## Input: genome annotation dataframe with gene IDs that match rownames of raw counts table

library(tidyverse); library(DESeq2); library(data.table) # required

# create metadata table
sample <- c('KP1', 'KP2', 'KP3', 'KP4', 'KP5', 'KP6', 'KP7', 'KP8', 'KP9', 'KP10', 'KP11', 'KP12') # sample IDs, also column headers of counts table (exp.all)
# three biological replicates of base media (control), sorbitol (sorb), mannitol (mann), or glucose (gluc) supplementation:
group <- c('control', 'control', 'control', 'sorb', 'sorb', 'sorb', 'mann', 'mann', 'mann', 'gluc', 'gluc', 'gluc') # treatment groups (media conditions)
library_trmts <- data.frame(sample,group)

# DESeq 
exp.all
IDs = row.names(exp.all)
colData <- library_trmts
ds <- DESeqDataSetFromMatrix(countData = as.data.frame(exp.all),
                             colData = colData, design = ~group)

dds <- estimateSizeFactors(ds)
dds <- dds[ rowSums(counts(dds)) > 3, ] # min n=3/group
dds <- DESeq(dds) # dds is key

vsd <- vst(dds, blind=FALSE) # vst-tranformed counts

# pairwise differential gene expression
b.sorb_ctrl <- results(dds, pAdjustMethod = "bonferroni", alpha = 0.01, contrast = c("group", "sorb", "control"))
b.sorb_mann <- results(dds, pAdjustMethod = "bonferroni", alpha = 0.01, contrast = c("group", "sorb", "mann"))
b.sorb_gluc <- results(dds, pAdjustMethod = "bonferroni", alpha = 0.01, contrast = c("group", "sorb", "gluc"))
b.mann_ctrl <- results(dds, pAdjustMethod = "bonferroni", alpha = 0.01, contrast = c("group", "mann", "control"))
b.gluc_ctrl <- results(dds, pAdjustMethod = "bonferroni", alpha = 0.01, contrast = c("group", "gluc", "control"))
b.mann_gluc <- results(dds, pAdjustMethod = "bonferroni", alpha = 0.01, contrast = c("group", "mann", "gluc"))

# store in list
list.res.b <- list(b.sorb_ctrl, b.sorb_mann, b.sorb_gluc, b.mann_ctrl, b.gluc_ctrl, b.mann_gluc)

# subset to adjusted P-value < 0.01
b.subset_p01 <- lapply(list.res.b, function(res){
  sub <- subset(res, padj < 0.01)
  df <- as.data.frame(sub)
})

# optional: annotate the differential gene expression results with locus tags and gene descriptions
b.subanno <- lapply(b.subset_p01, function(x){
  mat <- match(row.names(x), anno$GeneID) # Name not locus tag here
  genes <- anno[mat,] 
  x$name <- genes$Name
  x$locus_tag <- genes$locus_tag
  x$description <- genes$product
  x
})


btransform <- lapply(b.subanno, function(df){
  df$padj[is.na(df$padj)] <- 1
  df$neg.log.padj <- -log10(df$padj)
  df$log2FoldChange[is.na(df$log2FoldChange)] <- 0
  df$flip = 0
  idxs = df$log2FoldChange > 0
  df$flip[idxs] = df$neg.log.padj[idxs]
  df$flip[!idxs] = -1*df$neg.log.padj[!idxs]
  df$geneid <- row.names(df)
  df$sig <- cut(df$flip, c(-100,-2,2,100), labels = c("down", "n.s.", "up"))
  df
})

btransform[[1]]$comparison <- "sorb_ctrl" # comparison sorbitol vs. base medium
btransform[[2]]$comparison <- "sorb_mann" # comparison sorbitol vs. mannitol
btransform[[3]]$comparison <- "sorb_gluc" # comparison sorbitol  vs. glucose
btransform[[4]]$comparison <- "mann_ctrl" # comparison mannitol vs. base medium
btransform[[5]]$comparison <- "gluc_ctrl" # comparison glucose vs. base medium
btransform[[6]]$comparison <- "mann_gluc" # comparison mannitol vs. glucose
 
bon.comb <- rbindlist(btransform, use.names = TRUE)

# Output: combined dataframe with 6 pairwise-comparisons for differential gene expression 




