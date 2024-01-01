# scRNA-seq downstream analysis
# Part III: Clustering

# Major steps in this section of the analysis:
# Step 1: Significant PC identification
# Step 2: Clustering cells based on top PCs 
# Step 3: Exploration of quality control metrics
# Step 4: Searching for expected cell types using known cell type-specific gene markers

# Load packages ----------------------------------------------------------------
# install.packages(c("Seurat", "tidyverse", "RCurl", "cowplot", "dplyr"))
library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)
library(dplyr) 

# Set dir ----------------------------------------------------------------------
data_dir <- "/proj/SN_data/"
plots_dir <- "/proj/SN_plots/"

# Load data --------------------------------------------------------------------
setwd(data_dir)
load("integrated_seurat.RData")

# Step 1: Significant PCs ------------------------------------------------------
# PCA was performed on integrated_seurat.RData in 02_Integration.R

# Option 1: Use heatmap
# figure 27_PCs_heatmap_top9 ---------------------------------------------------
setwd(plots_dir)
DimHeatmap(integrated_seurat, dims = 1:9, cells = 500, balanced = TRUE)
# figure 28_PCs_heatmap_top21 --------------------------------------------------
# DimHeatmap(integrated_seurat, dims = 1:21, cells = 500, balanced = TRUE)

# Printing out the most variable genes driving PCs
feature_list <- list(x = integrated_seurat[["pca"]]) 
for (i in 1:length(x = integrated_seurat[["pca"]])) {
  feature_list[[i]] <- TopFeatures(object = integrated_seurat[["pca"]], dim = i, nfeatures = 20, balanced = T)
}

set(data_dir)
write.csv(feature_list, file ='Most_variable_genes_driving_PCs.csv', row.names = FALSE)

# Option 2: Use elbow plot
# figure 29_PCs_elbow_plot -----------------------------------------------------
ElbowPlot(object = integrated_seurat, ndims = 40)

# Option 3: Use jackstraw plot
setwd(data_dir)
load("filtered_seurat_rmDoublet.RData")
# JackStraw cannot be run on SCTransform-normalized data.
integrated_seurat_jackstraw <- JackStraw(filtered_seurat_rmDoublet, num.replicate = 100)
integrated_seurat_jackstraw <- ScoreJackStraw(integrated_seurat_jackstraw, dims = 1:20)

# Step 2: Clustering -----------------------------------------------------------
# Determine the K-nearest neighbor graph using the first 40 PCs to generate clusters
integrated_seurat <- FindNeighbors(object = integrated_seurat, dims = 1:40)
# Determine the clusters for various resolutions                                
integrated_seurat <- FindClusters(object = integrated_seurat, 
                                  resolution = c(0.4, 0.6, 0.8, 1.0, 1.4))
# Check different resolution outcomes 
integrated_seurat@meta.data %>% 
  View()

# Start from res0.8
# figure 31_UMAP_res0.8 --------------------------------------------------------
Idents(object = integrated_seurat) <- "integrated_snn_res.0.8"
DimPlot(integrated_seurat, reduction = "umap", label = TRUE, label.size = 6)

# figure 32_UMAP_res0.4 --------------------------------------------------------
Idents(object = integrated_seurat) <- "integrated_snn_res.0.4"
DimPlot(integrated_seurat,reduction = "umap", label = TRUE, label.size = 6)

# figure 33_UMAP_res0.2 --------------------------------------------------------
# figure 34_UMAP_res0.3 --------------------------------------------------------

# Step 3: QC -------------------------------------------------------------------
# Extract identity and sample information from seurat object to determine the number of cells per cluster per sample
n_cells <- FetchData(integrated_seurat, 
                     vars = c("ident", "orig.ident")) %>%
  dplyr::count(ident, orig.ident) %>%
  tidyr::spread(ident, n)

# View table
View(n_cells)

# figure 34_UMAP_res0.3_by_sample ----------------------------------------------
DimPlot(integrated_seurat, label = TRUE, split.by = "sample") + NoLegend()

# figure 35_UMAP_res0.3_by_cell_cycle ------------------------------------------
DimPlot(integrated_seurat, label = TRUE,  split.by = "Phase")  + NoLegend()

# Determine metrics to plot present in integrated_seurat@meta.data
metrics <-  c("nCount_RNA", "nFeature_RNA", "S.Score", "G2M.Score", "percent.mt")
# figure 36_UMAP_metrics -------------------------------------------------------
FeaturePlot(integrated_seurat, 
            reduction = "umap", 
            features = metrics,
            pt.size = 0.4, 
            sort.cell = TRUE,
            min.cutoff = 'q10',
            label = TRUE)

# Define the information in the seurat object of interest
columns <- c(paste0("PC_", 1:16), "ident",  "UMAP_1", "UMAP_2")
# Extract PCs from the seurat object
# Use the column names (PC_1, PC_2, PC_3, etc.) to pull out the coordinates or PC scores corresponding to each cell for each of the PCs.
pc_data <- FetchData(integrated_seurat, vars = columns)
View(integrated_seurat@reductions)

# Extract the UMAP coordinates for the first 10 cells
integrated_seurat@reductions$umap@cell.embeddings[1:10, 1:2]

# Visualize top 16 PCs
#   Adding cluster label to center of cluster on UMAP
umap_label <- FetchData(integrated_seurat, 
                        vars = c("ident", "UMAP_1", "UMAP_2"))  %>%
  group_by(ident) %>%
  summarise(x=mean(UMAP_1), y=mean(UMAP_2))

# Plot a UMAP plot for each of the PCs
# figure 37_UMAP_PCs -----------------------------------------------------------
map(paste0("PC_", 1:16), function(pc){
  ggplot(pc_data, 
         aes(UMAP_1, UMAP_2)) +
    geom_point(aes_string(color=pc), 
               alpha = 0.7) +
    scale_color_gradient(guide = FALSE, 
                         low = "grey90", 
                         high = "blue")  +
    geom_text(data=umap_label, 
              aes(label=ident, x, y)) +
    ggtitle(pc)
}) %>% 
  plot_grid(plotlist = .)

# Examine PCA results 
print(integrated_seurat[["pca"]], dims = 1:16, nfeatures = 5)

# Step 4: Marker gene exploration ----------------------------------------------
# Explore known cell type markers

# Select the RNA counts slot to be the default assay
DefaultAssay(integrated_seurat) <- "RNA"

# Normalize RNA data for visualization purposes
integrated_seurat <- NormalizeData(integrated_seurat, verbose = FALSE)

# CD14+ monocyte markers
FeaturePlot(integrated_seurat, 
            reduction = "umap", 
            features = c("CD14", "LYZ"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)

# KRAS+ markers 
FeaturePlot(integrated_seurat, 
            reduction = "umap", 
            features = c("KRAS","AURKA","EGFR","MAPK3","FGFR2","MET","MYC","IDH1","IDH2"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)

FeaturePlot(integrated_seurat, 
            reduction = "umap", 
            features = c("TP53","STK11","KEAP1","PIK3CA","SMAD4","PTCH1","APC","PTEN","BRAF"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)

FeaturePlot(integrated_seurat, 
            reduction = "umap", 
            features = c("SMARCA4","RB1","NF1","CTNNB1","RICTOR","NRAS"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)

# Epithelial-mesenchymal transition (EMT) markers
# In cancer cells: 
#  gain of expresison in VIM
#  loss of E-cadherin (CDH1) is typically associated with a gain of expression of the mesenchymal marker, N-cadherin (CDH1)
#  CLDN7 executes an oncogenic function so is typically upregulated in cancer cells
#  SNAI1 provokes the loss of epithelial markers and is typically observed to be overexpressed in cancers
#  ZEB1 and ZEB2 are typically upregulated 
FeaturePlot(integrated_seurat, 
            reduction = "umap", 
            features = c("CDH1","CDH2","VIM", "SNAI1","ZEB1","ZEB2","CLDN7"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)

# Check NRAS by samples
FeaturePlot(integrated_seurat, 
            reduction = "umap", 
            features = c("NRAS"), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE,
            split.by = 'sample')

# [OPTIONAL] Identification of all markers for each cluster 
# Find markers for every cluster compared to all remaining cells, report only the positive ones
# markers <- FindAllMarkers(object = integrated_seurat,
#                           only.pos = TRUE,
#                           logfc.threshold = 0.25)


# Identification of conserved markers in all conditions
# Use the original counts and not the integrated data
#  raw and normalized counts are stored in "RNA"
DefaultAssay(integrated_seurat) <- "RNA"
cluster0_conserved_markers <- FindConservedMarkers(integrated_seurat,
                                                   ident.1 = 0,
                                                   grouping.var = "sample",
                                                   only.pos = TRUE,
                                                   logfc.threshold = 0.25)

View(cluster0_conserved_markers)
# The output is a matrix containing a ranked list of putative markers listed by gene ID for the cluster 
#     gene: gene symbol
#     condition_p_val: p-value not adjusted for multiple test correction for condition
#     condition_avg_logFC: average log2 fold change for condition. Positive values indicate that the gene is more highly expressed in the cluster.
#     condition_pct.1: percentage of cells where the gene is detected in the cluster for condition
#     condition_pct.2: percentage of cells where the gene is detected on average in the other clusters for condition
#     condition_p_val_adj: adjusted p-value for condition, based on bonferroni correction using all genes in the dataset, used to determine significance
#     max_pval: largest p value of p value calculated by each group/condition
#     minimump_p_val: combined p value

# Note: Since each cell is being treated as a replicate this will result in inflated p-values within each group
#      A gene may have an incredibly low p-value < 1e-50 but that doesn't translate as a highly reliable marker gene.


sessionInfo()