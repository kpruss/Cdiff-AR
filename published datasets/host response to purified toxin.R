###########################
# Kali Pruss
# Pruss and Sonnenburg, Nature 2021

# Code for differential correlation analysis for host aldose reductase in response to C. difficile TcdA (Fig. 3g)

###########################


## Input: microarray expression data from D'auria et al 2013 (exp) where rows are microarray probes and columns are values for each mouse
## Input: metadata from the same study, (metadata), where rows are samples (mice) and columns include treatment group
## data table for affymetrix mouse microarray (ref)

# filter genes

library(DGCA); library(dplyr) # required

# filter microarray data 
exp_fil = filterGenes(exp, 
                      filterTypes = "central", filterCentralType = "median", filterCentralPercentile = 0.3)

# remove the single sample that was injected with both TcdA and B
exp_fil <- exp_fil %>% dplyr::select(-TcdAB_6)

# set the experimental design
groups <- factor(metadata$Injection)
des_inj <- model.matrix( ~ 0 + groups ) # contrasts for each treatment: sham injection, TcdA, and TcdB

set.seed(151)
res = ddcorAll(inputMat = exp_fil, design = des_inj,
                compare = c("groupsinjection: TcdA", "groupsinjection: Sham"), # change comparison here
                adjust = "perm", nPerm = 5, splitSet = c("1437133_x_at", "1456590_x_at","1448319_at")) # change gene here - switching leads to same

mat <- match(res$Gene1, ref$Affy_ID)
add <- ref[mat, ]
res$desc <- add$Gene.Title # add the names of genes
res$name <- add$Gene.Symbol # add gene symbols
res$entrez <- add$ENTREZ_GENE_ID # add Entrez IDs
res[order(res$pValDiff_adj), ] 

# output: ordered list of genes that are significantly correlated (positively or negatively) with three Akr1b3 probes