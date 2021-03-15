###########################
# Kali Pruss
# Pruss and Sonnenburg, Nature 2021

# Code for generation of in vitro transcriptional profiling heatmap (Extended Data Fig. 8d)

###########################


## Input: vst-normalized counts table (vsd). See `in vitro RNA-seq.R`
## Input: DESeq2 object for same counts table (dds). See`in vitro RNA-seq.R`
## Input: combined results of differential expression analysis (bon.comb), see `in vitro RNA-seq.R`

library(tidyverse); library(DESeq2); library(ComplexHeatmap); library(RColorBrewer) # required

m = assay(vsd) # make matrix
scaled_mat = t(scale(t(m))) # row-normalized 
mat_df = as.data.frame(scaled_mat)

df <- as.data.frame(colData(dds)[,c("sample", "group")])
df = df %>% select(-sample) 
ha = HeatmapAnnotation(df = df, 
                       col = list(group=c("control" = "#9F9D95",
                                          "gluc" = "#F5DA1C",
                                          "mann" = "#8481E7",
                                          "sorb" = "#04A90B")))

# Take significant differences between sorbitol supplementation vs. base medium and vs. mannitol
sorb_only_bon <- bon.comb %>% filter(comparison == "sorb_ctrl" | comparison == "sorb_mann") # 

mat_sorb_bon = subset(mat_df, rownames(mat_df) %in% sorb_only_bon$geneid) 
matrix_sorb_bon = as.matrix(mat_sorb_bon)

rows = match(rownames(matrix_sorb_bon), sorb_only_bon$geneid) 
sorb_only_bon_df = as.data.frame(sorb_only_bon) 
sorb_only_bon_df = sorb_only_bon_df[rows, ] 

rows_df = data.frame(rownames(matrix_sorb_bon), sorb_only_bon_df$locus_tag, sorb_only_bon_df$description) # change plot rownames to description instead of locus tag
rows_df = as.data.frame(sapply(rows_df, function(x) gsub("\n", "", x))) # remove line break signals for labeling

mat2 = matrix_sorb_bon # to test
rownames(mat2) = rows_df$sorb_only_bon_df.description

h = Heatmap(mat2, cluster_columns=FALSE, 
            show_row_names = TRUE, 
            col = colorRampPalette(rev(brewer.pal(n=11, name = "RdBu")))(200),
            top_annotation=ha,
            name = "normalized expression", show_column_names = FALSE, row_names_max_width = unit(2, "in"),
            row_names_gp = gpar(fontsize = 5))

draw(h, heatmap_legend_side = "left", 
     annotation_legend_side = "left")

# Output: heatmap of row-normalized vst-transformed counts
