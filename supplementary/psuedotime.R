library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(RColorBrewer)
library(monocle)
library(gprofiler2)
library(stringi)

## Importing Data
setwd("/wynton/home/brack/smahin/ISR_R_08_09_24/09_03_24/adult_aged/")

plot_save_path <- "/wynton/home/brack/smahin/ISR_R_08_09_24/09_03_24/adult_aged/figures_pseudo_5k/"

pp_data <- readRDS("ad_ag_musc_5k.rds")

expression_matrix <- GetAssayData(pp_data, assay = "RNA", slot = "counts")
expression_matrix <- as(expression_matrix, "sparseMatrix")
cell_metadata <- pp_data@meta.data
genes <- rownames(expression_matrix)
feature_data <- data.frame(gene_short_name = genes, row.names = genes)
cds <- newCellDataSet(expression_matrix,
                      phenoData = new("AnnotatedDataFrame", data = cell_metadata),
                      featureData = new("AnnotatedDataFrame", data = feature_data),
                      lowerDetectionLimit = 0.1)
pData(cds)$cell_id <- rownames(cds@phenoData@data)
fData(cds)$gene_id <- rownames(cds@featureData@data)

theme_aes <- theme(axis.text.x = element_blank(),
                   axis.ticks.x = element_blank(),  
                   axis.text.y = element_blank(),
                   axis.ticks.y = element_blank(),  
                   axis.title.x = element_text(size = 20),
                   axis.title.y = element_text(size = 20))
set.seed(123)

desired_genes <- c("Ccnd1", "Stoml2", "Ndufb7", "Ndufab1", "Cox7c", 
                   "Ldha", "Ndufa12", "Ndufc2", "Uqcrc2", "Atp5g2")

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

##ALL GENES
cds <- reduceDimension(cds, reduction_method='DDRTree', max_components = 2)
cds <- orderCells(cds)

AG_ct <- plot_cell_trajectory(cds, color_by="seurat_clusters", 
                              cell_size = 0.5, show_branch_points = FALSE) + theme_aes
ggsave(paste0(plot_save_path, "AG_cell_trajectory1.png"), plot = AG_ct, width = 7, height = 6,  dpi = 600)

gene_indices <- which(rownames(fData(cds)) %in% desired_genes)
cds_subset <- cds[gene_indices,]
AG_gp <- plot_genes_in_pseudotime(cds_subset, color_by='seurat_clusters',
                                  horizontal_jitter = TRUE)
ggsave(paste0(plot_save_path, "AG_genes_pseudotime1.png"), plot = AG_gp, dpi = 600, width = 7, height = 21)

