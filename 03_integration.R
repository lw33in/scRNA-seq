# scRNA-seq downstream analysis
# Part II: Integration

# Major steps in this section of the analysis:
# Step 1: Normalization, highly variable gene (HVG) identification, variance stabilization, and regression of unwanted variation
# Step 2: Cell cycle scoring 
# Step 3: SCTransform, a modeling framework for the normalization and variance stabilization of molecular count data 
# Step 4: [OPTIONAL] Integration of the samples using shared HGVs 


# Load packages ----------------------------------------------------------------
# install.packages("Seurat", "igraph", "tidyverse", "dplyr")
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
load("filtered_seurat_rmDoublet.RData")


# Step 1: Normalization --------------------------------------------------------
# LogNormalize
normalized_seurat <- NormalizeData(filtered_seurat_rmDoublet, normalization.method = "LogNormalize")

# HVG identification
# Use top 2000 HVGs to reflect high expressions
hvg_seurat <- FindVariableFeatures(normalized_seurat, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
# Label the top 15 HGVs 
top15 <- head(VariableFeatures(hvg_seurat), 15)
# figure 17_ctrl_HVGs_top2000-----------------------------------------------
# split by sample to plot separately below
split <- SplitObject(hvg_seurat, split.by = "sample") 

plot1 <- VariableFeaturePlot(split$ctrl) + 
  theme(legend.position="top")
plot2 <- LabelPoints(plot = plot1, points = top15, repel = TRUE) + 
  theme(legend.position="none")
plot2
rm(plot2)

# figure 18_exp_HVGs_top2000----------------------------------------------------
plot1 <- VariableFeaturePlot(split$exp) + 
  theme(legend.position="top")
plot2 <- LabelPoints(plot = plot1, points = top15, repel = TRUE) + 
  theme(legend.position="none")
plot2

# Scale the data to scale variation with expression level
#   1) adjusting the expression of each gene to give a mean expression across cells to be 0
#   2) scaling expression of each gene to give a variance across cells to be 1
hvg_seurat <- ScaleData(hvg_seurat)

# Export .RData object 
setwd(data_dir)
save(hvg_seurat, file="hvg_seurat.RData")

# Step 2: Cell cycle scoring ---------------------------------------------------
# https://satijalab.org/seurat/articles/cell_cycle_vignette.html
# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  
# Segregate this list into markers of G2/M phase and markers of S phase.
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
# Scoring
cell_cycle_seurat <- CellCycleScoring(hvg_seurat, g2m.features = g2m.genes, s.features = s.genes, set.ident = TRUE)

# figure 19_PCA_by_cell_cycle---------------------------------------------------
# Perform PCA
cell_cycle_seurat <- RunPCA(cell_cycle_seurat)
# Plot by cell cycle phase
DimPlot(cell_cycle_seurat,
        reduction = "pca",
        group.by= "Phase",
        split.by = "Phase")

# Note: use figure_19 to assess the necessity for cell cycle regression
#  If there are no large differences due to cell cycle phase, regression can be skipped.
#  Othersie, unwanted cell cycle variance needs to be regressed out. 

# ONLY perform the following chunk of code if cell cycle regression is needed.
# If not, skip to Step 3: SCTransform

# View cell cycle scores and phase assignments
head(cell_cycle_seurat[[]])
View(cell_cycle_seurat@meta.data) 
# figure 20_cell_cycle_markers -------------------------------------------------
RidgePlot(cell_cycle_seurat, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)

# Run a PCA on cell cycle genes reveals to check if cells separate entirely by phase
cell_cycle_seurat <- RunPCA(cell_cycle_seurat, features = c(s.genes, g2m.genes))
# figure 21_PCA_cell_cycle_genes -----------------------------------------------
DimPlot(cell_cycle_seurat)

# Cell cycle regression 
regress_cell_cycle_seurat <- ScaleData(cell_cycle_seurat, vars.to.regress = c("S.Score", "G2M.Score"), 
                                       features = rownames(cell_cycle_seurat))

# Re-run a PCA on the variable genes 
#  it should no longer return components associated with cell cycle in PCs
regress_cell_cycle_seurat <- RunPCA(regress_cell_cycle_seurat, features = VariableFeatures(regress_cell_cycle_seurat), nfeatures.print = 10)

# Run a PCA on only cell cycle genes
#  cells should no longer separate by cell-cycle phase
check_regression <- RunPCA(regress_cell_cycle_seurat, features = c(s.genes, g2m.genes))
# figure 22_PCA_post_cell_cycle_regress ----------------------------------------
DimPlot(check_regression)

# Step 3: SCTransform ----------------------------------------------------------
# sctransform accounts for cellular sequencing depth, or nUMIs.

# Note: Oftentimes, it is useful to regress out variation due to mitochondrial expression. 
#   However, if the differences in mitochondrial gene expression represent a biological phenomenon 
#   that may help to distinguish cell clusters, then do not regress the mitochondrial expression.

# Mitochondrial expression is an unwanted variance in this analysis.

# Adjust the limit for allowable object sizes within R. Otherwise, memory exhausted error may appear.
options(future.globals.maxSize = 4000 * 1024^2)

# Section 1: If cell cycle regression is performed, run this chunk. Otherwise, skip to Section 2
# Input: regress_cell_cycle_seurat
# Split seurat object by condition to perform cell cycle scoring and SCT on all samples
# SCTransform will rank the genes by residual variance and output the 3000 most variant genes.
split_seurat <- SplitObject(regress_cell_cycle_seurat, split.by = "sample") %>% 
  split_seurat[c("ctrl", "exp")]
for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- SCTransform(split_seurat[[i]], vars.to.regress = c("percent.mt"))
}

# Section 2: If cell cycle regression is NOT performed, run this chunk.
# Input: hvg_seurat
# Split seurat object by condition to perform cell cycle scoring and SCT on all samples
# SCTransform will rank the genes by residual variance and output the 3000 most variant genes.
split_seurat <- SplitObject(hvg_seurat, split.by = "sample") %>% 
  split_seurat[c("ctrl", "exp")]
for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- NormalizeData(split_seurat[[i]], verbose = TRUE)
  split_seurat[[i]] <- CellCycleScoring(split_seurat[[i]], g2m.features=g2m.genes, s.features=s.genes)
  split_seurat[[i]] <- SCTransform(split_seurat[[i]], vars.to.regress = c("percent.mt"))
}

# View the different assays stored in objects
# The most variable features will be the only genes stored inside the SCT assay.
split_seurat$ctrl@assays

# Export .RData object 
setwd(data_dir)
save(split_seurat, file="split_seurat.RData")

# Step 4: [OPTIONAL] Integration of the samples using shared HGVs  -------------
# If condition-specific clustering of the cells is seen, integration of the cells across conditions is needed.
# Input: split_seurat
integration_features <- SelectIntegrationFeatures(object.list = split_seurat, 
                                                  nfeatures = 3000) 
# Prepare the SCT list object for integration
split_seurat <- PrepSCTIntegration(object.list = split_seurat, 
                                   anchor.features = integration_features)
# Perform canonical correlation analysis (CCA) to find the best buddies or anchors and filter incorrect anchors
#   CCA identifies shared sources of variation between the conditions/groups. 
#   It is a form of PCA, in that it identifies the greatest sources of variation in the data, 
#   but only if it is shared or conserved across the conditions/groups. 
# This step roughly aligns the cells using the greatest shared sources of variation.
#   The shared highly variable genes are used because they are the most likely to represent those genes distinguishing the different cell types present.
integration_anchors <- FindIntegrationAnchors(object.list = split_seurat, 
                                              normalization.method = "SCT", 
                                              anchor.features = integration_features)
# Integrate across conditions
integrated_seurat <- IntegrateData(anchorset = integration_anchors, normalization.method = "SCT")

# PCA
integrated_seurat <- RunPCA(object = integrated_seurat)
# figure 23_PCA_integrated -----------------------------------------------------
setwd(plots_dir)
PCAPlot(integrated_seurat) 
# figure 24_PCA_integrated_by_sample -------------------------------------------
PCAPlot(integrated_seurat, split.by = "sample") 

#  UMAP
integrated_seurat <- RunUMAP(integrated_seurat, dims = 1:40, reduction = "pca")
# figure 25_UMAP_integrated ----------------------------------------------------
DimPlot(integrated_seurat, label = FALSE, group.by = 'sample')   
# figure 26_UMAP_integrated_by_sample -------------------------------------------
DimPlot(integrated_seurat, split.by = "sample", )     

# Export .RData object 
setwd(data_dir)
save(integrated_seurat, file="integrated_seurat.RData")


sessionInfo()
