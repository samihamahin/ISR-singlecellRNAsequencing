library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(RColorBrewer)
library(monocle)
library(gprofiler2)
library(stringi)
library(openxlsx)
library(dittoSeq)

## Importing Data

setwd("/Users/SamihaMahin/Desktop/brack_lab/ISR_R_08_09_24/09_03_24/")

plot_save_path <- "/Users/SamihaMahin/Desktop/brack_lab/ISR_R_08_09_24/09_03_24/figures/"
data_dir_path <- "/wynton/home/brack/smahin/isr_scRNA_data_dir/"

theme_aes <- theme(plot.title = element_blank(), axis.title.x = element_blank(), 
                   axis.title.y = element_blank(), axis.text.x = element_text(size = 14),
                   axis.text.y = element_text(size = 14), legend.text = element_text(size = 14))
set.seed(12345)

data_dir_1A <- paste0(data_dir_path, "1A/filtered_feature_bc_matrix")
data_dir_1B <- paste0(data_dir_path, "1B/filtered_feature_bc_matrix")

data_dir_2A <- paste0(data_dir_path, "2A/filtered_feature_bc_matrix")
data_dir_2B <- paste0(data_dir_path, "2B/filtered_feature_bc_matrix")

data_dir_3A <- paste0(data_dir_path, "3A/filtered_feature_bc_matrix")
data_dir_3B <- paste0(data_dir_path, "3B/filtered_feature_bc_matrix")

data_dir_1A_1 <- paste0(data_dir_path, "1A_1/filtered_feature_bc_matrix")
data_dir_1B_1 <- paste0(data_dir_path, "1B_1/filtered_feature_bc_matrix")

data_dir_2A_1 <- paste0(data_dir_path, "2A_1/filtered_feature_bc_matrix")
data_dir_2B_1 <- paste0(data_dir_path, "2B_1/filtered_feature_bc_matrix")

data_dir_3A_1 <- paste0(data_dir_path, "3A_1/filtered_feature_bc_matrix")
data_dir_3B_1 <- paste0(data_dir_path, "3B_1/filtered_feature_bc_matrix")

create_seurat_from_10x <- function(data_dir, sample, group) {
  cell_data <- Read10X(data.dir = data_dir)
  seurat_object <- CreateSeuratObject(counts = cell_data)
  seurat_object <- AddMetaData(object = seurat_object, metadata = sample, col.name = "Sample")
  seurat_object <- AddMetaData(object = seurat_object, metadata = group, col.name = "Group")
  return(seurat_object)}



## Create Seurat Object

ad_1A <- create_seurat_from_10x(data_dir_1A, "Adult-N1", "Adult")
ad_1B <- create_seurat_from_10x(data_dir_1B, "Adult-N2", "Adult")

ad_1A_1 <- create_seurat_from_10x(data_dir_1A_1, "Adult-N3", "Adult")
ad_1B_1 <- create_seurat_from_10x(data_dir_1B_1, "Adult-N4", "Adult")

ag_2A <- create_seurat_from_10x(data_dir_2A, "Aged-N1", "Aged")
ag_2B <- create_seurat_from_10x(data_dir_2B, "Aged-N2", "Aged")

ag_2A_1 <- create_seurat_from_10x(data_dir_2A_1, "Aged-N3", "Aged")
ag_2B_1 <- create_seurat_from_10x(data_dir_2B_1, "Aged-N4", "Aged")

agS_3A <- create_seurat_from_10x(data_dir_3A, "Aged Sal-N1", "Aged Sal" )
agS_3B <- create_seurat_from_10x(data_dir_3B, "Aged Sal-N2", "Aged Sal")

agS_3A_1 <- create_seurat_from_10x(data_dir_3A_1, "Aged Sal-N3", "Aged Sal" )
agS_3B_1 <- create_seurat_from_10x(data_dir_3B_1, "Aged Sal-N4", "Aged Sal")

all_data <- merge(ad_1A, y = c(ad_1B, ag_2A, ag_2B, agS_3A, agS_3B, 
                               ad_1A_1, ad_1B_1, ag_2A_1, ag_2B_1, agS_3A_1, agS_3B_1), 
                  add.cell.ids = c("AD-1A", "AD-1B", "AG-2A", "AG-2B","AGS-3A", "AGS-3B", 
                                   "AD-1A-1", "AD-1B-1", "AG-2A-1", "AG-2B-1","AGS-3A-1", "AGS-3B-1"), 
                  project = "ALLDATA")
all_data <- JoinLayers(all_data)

## Data Preprocesssing
all_data[["percent.mt"]] <- PercentageFeatureSet(all_data, pattern = "^mt-")

SaveSeuratRds(all_data, file = "all_data.rds")

all_data <- readRDS("all_data.rds")

plot1 <- FeatureScatter(all_data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all_data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

VlnPlot(all_data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

pp_data <- subset(all_data, subset = percent.mt < 10 & nFeature_RNA < 5000)

bad_genes <- c("Gm42418", "AY036118")
pp_data <- subset(pp_data,features=setdiff(rownames(pp_data),bad_genes))

pp_data <- NormalizeData(pp_data, normalization.method = "LogNormalize", scale.factor = 10000)
pp_data <- FindVariableFeatures(pp_data, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(pp_data), 10)
plot1 <- VariableFeaturePlot(pp_data)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
ggsave(paste0(plot_save_path, "var_feat_plot1.png"), plot = plot1, width = 6, height = 4, dpi = 600)
ggsave(paste0(plot_save_path, "var_feat_plot2.png"), plot = plot2, width = 6, height = 4, dpi = 600)

all.genes <- rownames(pp_data)
pp_data <- ScaleData(pp_data, features = all.genes)

pp_data <- RunPCA(pp_data, npcs = 10, ndims.print = 1:10, features = VariableFeatures(object = pp_data))
print(pp_data[["pca"]], dims = 1:10, nfeatures = 5)
VizDimLoadings(pp_data, dims = 1:2, reduction = "pca")  
DimPlot(pp_data, reduction = "pca")
ElbowPlot(pp_data)


## Clustering to find MuSCs
pp_data <- FindNeighbors(pp_data, dims = 1:10)
pp_data <- FindClusters(pp_data, resolution = 0.25) 

pp_data <- RunUMAP(pp_data, dims = 1:10)
pp_cluster <- DimPlot(pp_data, reduction = "umap") + theme_aes
pp_cluster
ggsave(paste0(plot_save_path, "pp_cluster_umap.png"), plot = pp_cluster, width = 6, height = 4, dpi = 600)

pp_expr <- FeaturePlot(pp_data, features= c("Pax7", "Myod1"))
pp_expr
ggsave(paste0(plot_save_path, "pp_expr.png"), plot = pp_expr, width = 12, height = 4, dpi = 600)

### Subsetting Clusters that contian MuSC markers (Pax7 and Myod1)
musc_data <- subset(pp_data, idents=c(0, 1))
DimPlot(musc_data, reduction = "umap")

## ISR Cluster analysis

### Re-clustering to find three clusters
musc_data <- FindNeighbors(musc_data)
musc_data <- FindClusters(musc_data, resolution = 0.15)
DimPlot(musc_data, reduction = "umap")

### Refactoring Data according to ISR gene levels 
init_vln = VlnPlot(musc_data, features = c("Atf4", "Atf6","Ppp1r15b", 
                                           "Ddit3", "Nars"), pt.size = 0, group.by = "seurat_clusters")
init_vln
ggsave(paste0(plot_save_path, "init_vln.png"), plot = init_vln, width = 6, height = 4, dpi = 600)

### Swap cluster labels (0->2, 1->0, 2->1)
DimPlot(musc_data, reduction = "umap")
musc_data$seurat_clusters <- ifelse(musc_data$seurat_clusters == 0, 2,
                                    ifelse(musc_data$seurat_clusters == 1, 0,
                                           ifelse(musc_data$seurat_clusters == 2, 1, NA)))
VlnPlot(musc_data, features = c("Atf4", "Atf6","Ppp1r15b", 
                                "Ddit3", "Nars"), pt.size = 0, group.by = "seurat_clusters")
musc_data$seurat_clusters <- factor(musc_data$seurat_clusters, levels = c("0", "1", "2"),
                                    labels = c("Cluster 0", "Cluster 1", "Cluster 2"))

refactor_vln = VlnPlot(musc_data, features = c("Atf4", "Atf6","Ppp1r15b", 
                                           "Ddit3", "Nars"), pt.size = 0, group.by = "seurat_clusters")
refactor_vln
ggsave(paste0(plot_save_path, "refactor_vln.png"), plot = refactor_vln, width = 6, height = 4, dpi = 600)
### UMAP CLUSTER 
cluster_umap <- DimPlot(musc_data, reduction = "umap", group.by = "seurat_clusters") + theme_aes
cluster_umap
ggsave(paste0(plot_save_path, "redo_cluster_umap.png"), plot = cluster_umap, width = 6, height = 4, dpi = 600)

### UMAP GROUP
group_colors <- c("Adult" = "#1b9e77", "Aged" = "#d95f02", "Aged Sal" = "#7570b3") 
group_umap <- DimPlot(musc_data, reduction = "umap", group.by = "Group", cols = group_colors, pt.size=0.25) + theme_aes
group_umap
ggsave(paste0(plot_save_path, "group_umap.png"), plot = group_umap, width = 6, height = 4, dpi = 600)

### VIOLIN GENE LIST
isr_genes_vln <- c("Atf4", "Atf6","Ppp1r15b", "Ddit3", "Nars")

violin_plot_gene_lst <- function(data, genes, plot_save_path) {
  for (gene in genes) {
    gene_vln <- VlnPlot(data, features = c(gene), pt.size = 0, group.by = "seurat_clusters") + 
      theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
            axis.text.y = element_text(size = 14), axis.text.x = element_blank(),
            legend.position = "none", plot.title = element_text(size= 25))
    file_name <- paste0(plot_save_path, gene, "_vln.png")
    ggsave(file_name, plot = gene_vln, width= 4, height = 3, dpi = 600)
  }
}

violin_plot_gene_lst(musc_data, isr_genes_vln, "/wynton/home/brack/smahin/ISR_R_08_09_24/09_03_24/figures/vln_genes/")

### UMAP GENE LIST
isr_genes_umap <- c("Pax7", "Ccnd1", "Cdkn1a", "Atf4", "Ppp1r15a", "Ddit3")

umap_plot_gene_lst <- function(data, genes, plot_save_path) {
  for (gene in genes) {
    gene_umap <- FeaturePlot(musc_data, features= c(gene)) + 
      theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
            plot.title = element_text(size= 25),axis.text.y = element_text(size = 14),
            axis.text.x = element_text(size = 14),legend.text = element_text(size = 14)) +
      guides(color = guide_colorbar(barwidth = 1.5, barheight = 10))
    
    file_name <- paste0(plot_save_path, gene, "_umap.png")
    ggsave(file_name, plot = gene_umap, width= 4, height= 3, dpi = 600)
  }
}
umap_plot_gene_lst(musc_data, isr_genes_umap, "/wynton/home/brack/smahin/ISR_R_08_09_24/09_03_24/figures/umap_genes/")

## HISTOGRAMS 
counts <- musc_data@meta.data %>%
  group_by(seurat_clusters, Group) %>%
  summarise(Count = n(), .groups = 'drop')
total_counts_per_cluster <- counts %>%
  group_by(seurat_clusters) %>%
  summarise(TotalClusterCount = sum(Count), .groups = 'drop')
total_counts_per_group <- counts %>%
  group_by(Group) %>%
  summarise(TotalGroupCount = sum(Count), .groups = 'drop')
counts <- counts %>%
  left_join(total_counts_per_cluster, by = "seurat_clusters") %>%
  left_join(total_counts_per_group, by = "Group")
counts <- counts %>%
  mutate(NormalizedPercentagePerCluster = (Count / TotalClusterCount) * 100,
         NormalizedPercentagePerGroup = (Count / TotalGroupCount) * 100)

theme_hist <- theme(axis.title.x= element_blank(), axis.text.x = element_text(size = 30, color="black"),
                    axis.text.y = element_text(size = 30, color= "black"), axis.ticks = element_blank(),
                    axis.title.y = element_text(face = "bold", size = 50, color="black"), legend.text = element_text(size=30),
                    legend.title = element_blank())

hist_group <- ggplot(counts, aes(x = factor(seurat_clusters), y = NormalizedPercentagePerCluster, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") + labs(y = "Percentage") + 
  scale_fill_manual(values = group_colors) + theme_hist
hist_group
ggsave(paste0(plot_save_path, "group_hist.png"), plot = hist_group, width = 10, height = 6,  dpi = 600)

hist_cluster <- ggplot(counts, aes(x = Group, y = NormalizedPercentagePerGroup, fill = factor(seurat_clusters))) +
  geom_bar(stat = "identity", position = "dodge") + labs(y = "Percentage") + theme_hist
hist_cluster
ggsave(paste0(plot_save_path, "cluster_hist.png"), plot = hist_cluster, width = 10, height = 6,  dpi = 600)

## DIFFERENTIAL EXPRESSION ANALYSIS
Idents(object = musc_data) <- "seurat_clusters"

musc_data <- readRDS("musc_data.rds")

#logfc_threshold = 0.25, min_pct = 0.10 
make_GO_text <- function(data, cluster1, cluster2, comparison,
                         logfc_threshold = 0.25, min_pct = 0.10, p_val_thresh = 0.05, dir) {
  markers <- FindMarkers(data, ident.1 = cluster1, ident.2 = cluster2,
                         logfc.threshold = logfc_threshold, min.pct = min_pct, 
                         test.use = "wilcox", return.thresh = p_val_thresh)
  
  markers <- markers[markers$p_val_adj < p_val_thresh, ]
  
  cluster1_markers <- markers[markers$avg_log2FC > 0, ]
  cluster2_markers <- markers[markers$avg_log2FC < 0, ]
  
  cluster1_name = paste0(stri_sub(cluster1, 1, 1), stri_sub(cluster1, -1, -1))
  cluster2_name = paste0(stri_sub(cluster2, 1, 1), stri_sub(cluster2, -1, -1))
  
  write.table(rownames(cluster1_markers), file = paste0(dir, cluster1_name, comparison, "markers_genes.txt"), 
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(rownames(cluster2_markers), file = paste0(dir, cluster2_name, comparison, "markers_genes.txt"), 
              row.names = FALSE, col.names = FALSE, quote = FALSE)}

make_GO_text(musc_data, "Cluster 0", "Cluster 1", "_c0vsc1_", dir= "GO/")
make_GO_text(musc_data, "Cluster 0", "Cluster 2", "_c0vsc2_", dir= "GO/")
make_GO_text(musc_data, "Cluster 1", "Cluster 2", "_c1vsc2_", dir= "GO/")

## HEATMAP

musc_data$Cluster <- musc_data$seurat_clusters

ordered_cells <- rownames(musc_data@meta.data[order(musc_data@meta.data$Cluster, musc_data@meta.data$Group), ])

custom_colors <- colorRampPalette(c("magenta", "black", "yellow"))

dittoHeatmap(musc_data,
             genes = c("Ccnd1", "Stoml2", "Ndufb7", "Ndufab1", "Cox7c",
                       "Ldha", "Ndufa12", "Ndufc2", "Uqcrc2", "Atp5g2"),
             annot.by = c("Cluster", "Group"),
             cells.use = ordered_cells,  # Order cells by Cluster and then Group
             cluster_rows = FALSE,
             heatmap.colors = custom_colors(1000),
             annotation_colors = list(Group = c("Adult" = "#1b9e77", "Aged" = "#d95f02", "Aged Sal" = "#7570b3"), 
                                      Cluster = c("Cluster 0" = "#f8766d", "Cluster 1" = "#00ba38", "Cluster 2" = "#609cff")))
make_GO_excel <- function(url, savename, dir){
  david_data <- read.table(url, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  david_data$Term <- sub(".*~", "", david_data$Term)
  david_data$log10_PValue <- -log10(david_data$PValue)
  first_15_rows <- head(david_data, 15)
  write.xlsx(first_15_rows, file=paste0(dir, savename, ".xlsx"))
}

c0_c0vsc1_url <- "https://david.ncifcrf.gov/data/download/chart_9C43A83792021726770533350.txt"
c1_c0vsc1_url <- "https://david.ncifcrf.gov/data/download/chart_9C43A83792021726770729634.txt"

c0_c0vsc2_url <- "https://david.ncifcrf.gov/data/download/chart_9C43A83792021726770887653.txt"
c2_c0vsc2_url <- "https://david.ncifcrf.gov/data/download/chart_18A554CA153F1726771416899.txt"

c1_c1vsc2_url <- "https://david.ncifcrf.gov/data/download/chart_18A554CA153F1726771548782.txt"
c2_c1vsc2_url <- "https://david.ncifcrf.gov/data/download/chart_18A554CA153F1726771636390.txt"

make_GO_excel(c0_c0vsc1_url, "C0_c0vsc1", "GO/")
make_GO_excel(c1_c0vsc1_url, "C1_c0vsc1", "GO/")

make_GO_excel(c0_c0vsc2_url, "C0_c0vsc2", "GO/")
make_GO_excel(c2_c0vsc2_url, "C2_c0vsc2", "GO/")

make_GO_excel(c1_c1vsc2_url, "C1_c1vsc2", "GO/")
make_GO_excel(c2_c1vsc2_url, "C2_c1vsc2", "GO/")

## Save Seurat object
SaveSeuratRds(musc_data, file = "musc_data.rds")
