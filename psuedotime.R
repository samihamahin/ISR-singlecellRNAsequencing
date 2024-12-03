BiocManager::install(version = "3.19")

BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'HDF5Array',
                       'terra', 'ggrastr'))

install.packages("devtools")
devtools::install_github('cole-trapnell-lab/monocle3')

library(monocle3)
library(Seurat)
library(ggplot2)

setwd("/Users/SamihaMahin/Desktop/brack_lab/ISR_R_08_09_24/09_03_24/")

plot_save_path <- "/Users/SamihaMahin/Desktop/brack_lab/ISR_R_08_09_24/09_03_24/figures/"

#musc_data <- readRDS("musc_data.rds")
musc_data <- readRDS("adult_aged/ad_ag_musc_5k.rds")

cds <- SeuratWrappers::as.cell_data_set(musc_data)


fData(cds)$gene_short_name <- rownames(fData(cds))
recreate.partition <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partition) <- cds@colData@rownames
recreate.partition <- as.factor(recreate.partition)

cds@clusters$UMAP$partitions <- recreate.partition

list_cluster <- musc_data@active.ident
cds@clusters$UMAP$clusters <- list_cluster

cds@int_colData@listData$reducedDims$UMAP <- musc_data@reductions$umap@cell.embeddings

cluster_before_trajectory <- plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = FALSE, 
                                        group_label_size = 5) + theme(legend.position = "right")

cds <- learn_graph(cds, use_partition = FALSE)

cds <- order_cells(cds, reduction_method = 'UMAP', root_cells = colnames(cds[,clusters(cds) == 'Cluster 0']))

plot_cells(cds,
           color_cells_by = 'pseudotime',
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE,
           show_trajectory_graph = FALSE) + theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
                                         plot.title = element_text(size= 25),axis.text.y = element_text(size = 14),
                                         axis.text.x = element_text(size = 14),legend.text = element_text(size = 14),
                                         legend.title = element_text(size = 14)) + guides(color = guide_colorbar(barwidth = 1.5, 
                                                                                                         barheight = 10)) + labs(color = "Psuedotime") 
pseudo_plot <- plot_cells(cds,
           color_cells_by = 'pseudotime',
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE,
           show_trajectory_graph = FALSE,
           cell_size = 0.2) +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        plot.title = element_text(size= 25),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14)) +
  guides(color = guide_colorbar(barwidth = 1.5, barheight = 10)) +
  scale_color_viridis_c(name = "Pseudotime", option = "turbo")

pseudo_plot

pseudo_path <- paste0(plot_save_path, "pseudo_plot_ADAG.png")

ggsave(pseudo_path, plot = pseudo_plot, width= 5, height= 3, dpi = 600)


pseudotime(cds)
cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))

ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(seurat_clusters, monocle3_pseudotime, median), fill = 'seurat_clusters')) +
  geom_boxplot()
