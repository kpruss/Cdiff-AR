###########################
# Kali Pruss
# Pruss and Sonnenburg, Nature 2021

# Code for transcripts-per million calculation (Fig. 2a)

###########################

## Input: takes raw counts table (x) with rRNA reads removed where:
# rows are genes, with unique row IDs and
# columns are samples

# Input: genome annotation (anno.630) with gene IDs that match rownames for counts table and contains
# gene lengths

library(dplyr); library(tidyr)

# add gene lengths to counts table 
mat = match(row.names(x), anno.630$GeneID) # appropriate matching is crucial
add =  anno.630[mat, ]
x$length = add$Length # Length in annotation file specifies gene length in bases 

no.na = x %>% drop_na() # remove loci with NAs (if not all loci in counts table have annotations)

counts = no.na %>% select(-length, -Name) # keep raw counts columns only
genes =  no.na %>% select(Name, length) # keep gene names and lengths only

# function to calculate rpkm
rpkm <- function(counts, lengths) {
  rate <- counts / lengths 
  rate / sum(counts) * 1e6
}

# function to calculate tpm
tpm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}

# apply functions to gene counts
rpkms <- apply(counts, 2, function(x) rpkm(x, genes$length))
tpms <- apply(counts, 2, function(x) tpm(x, genes$length))

# Output: matrix of rpkm- or tpm-normalized expression values

