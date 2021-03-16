###########################
# Kali Pruss
# Pruss and Sonnenburg, Nature 2021

# Code for single cell aldose reductase expression (Extended Data Figs. 6d,e)

###########################


## Input: dataframe (cell.table1, cell.table2) with gene expression for gene(s) of interest in columns for each single cell in rows. Associated metadata (cell type and condition) for each cell in columns
## Here, data was downloaded as a gct file from the Broad Institute Single Cell Portal

library(tidyverse); library(rstatix) # required


## ED Fig 6d
no01 = cell.table1 %>% filter(AKR1B1 > 0) # filter for cells with positive expression

# pairwise wilcox text
no01 %>%
  group_by(Cluster) %>% 
  pairwise_wilcox_test(
    AKR1B1 ~ Health,
    comparisons = list(c("Inflamed", "Uninflamed"), c("Inflamed", "Healthy"), c("Healthy", "Uninflamed")),
    p.adjust.method = "bonferroni"
  )
no01 %>% group_by(Cluster, Health) %>% summarise(mean(AKR1B1)) # mean expression

## ED Fig 6e
no02 = cell.table2 %>% filter(Akr1b3 > 0) # filter for cells with positive expression

# pairwise wilcox text
no02 %>%
  group_by(Cluster) %>% 
  pairwise_wilcox_test(
    Akr1b3 ~ treatment,
    p.adjust.method = "bonferroni"
  )

no02 %>% group_by(Cluster, treatment) %>% summarise(mean(Akr1b3)) # mean expression

