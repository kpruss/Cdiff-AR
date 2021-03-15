###########################
# Kali Pruss
# Pruss and Sonnenburg, Nature 2021

# Code for in vivo transcriptional profiling differential gene expression (Fig. 1b, Supplementary Table 2)

###########################

## Input: takes raw counts (exp.all) table with rRNA reads removed where:
# rows are genes, with unique row IDs and
# columns are samples, with unique column names that match sample names in metadata provided

# Input: takes metadata dataframe (library_trmts) where first column (sample IDs) matches column headers in counts table and
# additional metadata is provided in subsequent columns, in this case Treatment is the comparison of interest

library(DESEq2) # required

exp.all # raw counts dataframe
library_trmts # metadata dataframe

geneIDs = row.names(exp.all)

colData <- DataFrame(library_trmts)
ds <- DESeqDataSetFromMatrix(countData = as.data.frame(exp.all),
                             colData = colData, design = ~Treatment)
ds

dds <- estimateSizeFactors(ds)
dds <- dds[ rowSums(counts(dds)) > 3, ] # min group n=3

dds <- DESeq(dds)
res <- results(dds)

resOrdered <- res[order(res$pvalue),]
resSig05 <- subset(resOrdered, padj < 0.05)
resSig05 <- as.data.frame(resSig05) # p value cut-off 0.05

l <- list(resSig05) # alter this for following function

transform <- lapply(l, function(df){
  df$padj[is.na(df$padj)] <- 1
  df$neg.log.padj <- -log10(df$padj)
  df$log2FoldChange[is.na(df$log2FoldChange)] <- 0
  df$flip = 0
  idxs = df$log2FoldChange > 0
  df$flip[idxs] = df$neg.log.padj[idxs]
  df$flip[!idxs] = -1*df$neg.log.padj[!idxs]
  df$geneid <- row.names(df); df
})

df <-as.data.frame(transform[[1]])

# Output: dataframe with positive log2-fold change for significant in group of interest (here, WT), 
# negative log2-fold change for significance in other group (here, Tox-)
