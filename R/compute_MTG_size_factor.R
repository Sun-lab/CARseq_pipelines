# Compute MTG size factor
#
# Calculate the cell size factors using Paul's "MTG/all_genes_MTG.rds" as the starting point.
# Using the intersection of genes in MTG single cell, the SCZ dataset and the ASD dataset,
# and then calculate the cell size factors.
# The calculated size factor would work for transforming CIBERSORT or ICeDT
# estimates to cell fractions.

library(SummarizedExperiment)
library(tidyverse)

all_genes_MTG = readRDS("../MTG/all_genes_MTG.rds")
rse_ASD = readRDS("../data/ASD_rse_filtered_with_SVs.rds")
rse_SCZ = readRDS("../data/rse_filtered_SV.rds")
gene_names_in_common = intersect(intersect(rowData(rse_ASD)$gene_name,
                                           rowData(rse_SCZ)$gene_name),
                                 all_genes_MTG$anno$gene)
length(gene_names_in_common)
# check if any of these are duplicated:  X-Y homologous genes
length(rowData(rse_ASD)$gene_name[rowData(rse_ASD)$gene_name %in% gene_names_in_common])
length(rowData(rse_SCZ)$gene_name[rowData(rse_SCZ)$gene_name %in% gene_names_in_common])
length(all_genes_MTG$anno$gene[all_genes_MTG$anno$gene %in% gene_names_in_common])
# remove X-Y homologous genes
gene_names = gene_names_in_common %>%
    setdiff(rowData(rse_ASD)$gene_name[duplicated(rowData(rse_ASD)$gene_name)]) %>%
    setdiff(rowData(rse_SCZ)$gene_name[duplicated(rowData(rse_SCZ)$gene_name)])
length(gene_names)

# calculate size factors
cell_sizes = colSums(all_genes_MTG$SIG / all_genes_MTG$anno$gene_length)
cell_sizes = cell_sizes / max(cell_sizes)
cell_sizes
saveRDS(file = "../MTG/cell_sizes_MTG.rds",
        list(cell_sizes = cell_sizes,
             gene_names = gene_names))

# further validate that the gene lengths are generally similar
cor(as.matrix(data.frame(
  ASD_gene_length = rowData(rse_ASD)$gene_length[match(gene_names, rowData(rse_ASD)$gene_name)],
  SCZ_gene_length = rowData(rse_SCZ)$gene_length[match(gene_names, rowData(rse_SCZ)$gene_name)],
  MTG_gene_length = all_genes_MTG$anno$gene_length[match(gene_names, all_genes_MTG$anno$gene)])))
