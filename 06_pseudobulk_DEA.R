# scRNA-seq downstream analysis
# Part V: Pseudobulk differential expression analysis (DEA)

# Rationale for Pseudobulk DEA:
# Seurat has the ability to perform DEA, however, the p-values are often inflated as 
#  each cell is treated as an independent sample. In this case, variances detected represent
#  the differences among individual cells, not across populations / clusters. 
#  The best practice to perform DEA in scRNA-seq is the pseudobulk method.

# Major steps in this section of the analysis:
# Step 1: Extract raw counts and metadata
# Step 2: Transform data into one sample-level dataset per cluster / cell type 
# Step 3: Aggregate the counts and metadata to the sample level 
# Step 4: Splitting the counts matrix by cell type
# Step 5: Generating matching metadata at the sample-level
# Step 6: DEA using DESeq2 
  # Step 6.1: Normalization
  # Step 6.2: Unsupervised clustering analysis 
  # Step 6.3: Modeling and shrinking
  # Step 6.4: DEA

# Note: At least two biological repolicates are required for DEA.

# Load packages ----------------------------------------------------------------
# install.packages(c("Matrix", "Seurat", "dplyr", "RCurl", "cowplot"))
# BiocManager::install("SingleCellExperiment")
# BiocManager::install(c("DESeq2", "edgeR"))
library(Matrix)
library(Seurat)
library(dplyr)
library(DESeq2)
library(SingleCellExperiment)
library(tidyverse)
library(cowplot)
library(Matrix.utils)
library(edgeR)
library(Matrix)
library(reshape2)
library(S4Vectors)
library(pheatmap)
library(apeglm)
library(png)
library(RColorBrewer)
library(data.table)

# Set dir ----------------------------------------------------------------------
data_dir <- "/proj/SN_data/"
plots_dir <- "/proj/SN_plots/"

# Load data --------------------------------------------------------------------
setwd(data_dir)
load("integrated_seurat.RData")

# Step 1: Extract raw counts and metadata --------------------------------------
counts <- integrated_seurat@assays$RNA@counts 
metadata <- integrated_seurat@meta.data
metadata$cluster_id <- factor(integrated_seurat@active.ident)
# Create SingleCellExperiment (sce) object
#  This object is a type of specialized list. See details in PMID: 31792435
pseudobulk_sce <- SingleCellExperiment(assays = list(counts = counts), 
                            colData = metadata)

# Export data 
save(pseudobulk_sce, file="pseudobulk_sce.RData")

# Check the assays
assays(pseudobulk_sce)

# Explore the raw counts: a cell by gene sparse matrix
#  rows = genes, columns = cells 
# `head()` is not helpful here as it would show all the columns
dim(counts(pseudobulk_sce))  
counts(pseudobulk_sce)[1:6, 1:6]

# Explore the cellular metadata
#   `cluster_id` shows the cell types assigned
dim(colData(pseudobulk_sce))
head(colData(pseudobulk_sce))

# Step 2: Transform data into one sample-level dataset per cluster / cell type -
# Determine the number of clusters and the cluster names (cell types)
# Extract unique names of clusters (= levels of cluster_id factor variable)
cluster_names <- levels(colData(pseudobulk_sce)$cluster_id)
cluster_names
length(cluster_names)

# Extract unique names of samples (= levels of sample_id factor variable)
sample_names <- levels(colData(pseudobulk_sce)$sample_id)
sample_names

# Total number of samples
length(sample_names)

# Step 3: Aggregate the counts and metadata to the sample level ----------------
# Subset metadata to include only the variables you want to aggregate across 
groups <- colData(pseudobulk_sce)[, c("cluster_id", "sample_id")] # a df
head(groups)

# Aggregate across cluster-sample groups
# transposing row/columns to have cell_ids as row names matching those of groups
aggr_counts <- aggregate.Matrix(t(counts(pseudobulk_sce)), 
                                groupings = groups, fun = "sum") # a sparse matrix 

# Explore output matrix
class(aggr_counts)
dim(aggr_counts)
aggr_counts[1:6, 1:6]

# Step 4: Splitting the counts matrix by cell type -----------------------------
#   1) Transform the matrix back, so that the genes are listed in rows and the samples are in columns
#   2) Split our matrix by cell type
# Transpose aggregated matrix to have genes as rows and samples as columns
aggr_counts <- t(aggr_counts)
aggr_counts[1:6, 1:6]

## Exploring structure of function output (list)
tstrsplit(colnames(aggr_counts), "_") %>% str()

## Comparing the first 10 elements of our input and output strings
head(colnames(aggr_counts), n = 10)
head(tstrsplit(colnames(aggr_counts), "_")[[1]], n = 10)

# Using which() to look up tstrsplit() output
macrophage_idx <- which(tstrsplit(colnames(aggr_counts), "_")[[1]] == "Macrophage")
macrophage_idx
colnames(aggr_counts)[macrophage_idx]
aggr_counts[1:10, macrophage_idx]

# Cell types in a vector called cluster_names
# cluster_names

# Loop over all cell types to extract corresponding counts, and store information in a list
# Initiate empty list
counts_ls <- list()
for (i in 1:length(cluster_names)) {
  # Extract indexes of columns in the global matrix that match a given cluster
  column_idx <- which(tstrsplit(colnames(aggr_counts), "_")[[1]] == cluster_names[i])
  # Store corresponding sub-matrix as one element of a list
  counts_ls[[i]] <- aggr_counts[, column_idx]
  names(counts_ls)[i] <- cluster_names[i]
  
}

# Explore the different components of the list
str(counts_ls)

# Step 5: Generating matching metadata at the sample-level----------------------
head(colData(sce))

# Extract sample-level variables
metadata <- colData(pseudobulk_sce) %>% 
  as.data.frame() %>% 
  dplyr::select(group_id, patient_id, sample_id)
dim(metadata)
head(metadata)

# Exclude duplicated rows
metadata <- metadata[!duplicated(metadata), ]
# Rename rows
rownames(metadata) <- metadata$sample_id
dim(metadata)
head(metadata)

# Check the number of cells per sample and cluster
t <- table(colData(pseudobulk_sce)$sample_id,
           colData(pseudobulk_sce)$cluster_id)
t[1:6, 1:6]

# Append this cell count information to our generic metadata table
#   and generate one metadata data frame specific of each cell type
#   make sure that for each cell type (cluster), the row names of the metadata data frame 
#   match the column names of the corresponding counts matrix.
## Initiate empty list
metadata_ls <- list()
for (i in 1:length(counts_ls)) {
  # Initiate a data frame for cluster i with one row per sample (matching column names in the counts matrix)
  df <- data.frame(cluster_sample_id = colnames(counts_ls[[i]]))
  # Use tstrsplit() to separate cluster (cell type) and sample IDs
  df$cluster_id <- tstrsplit(df$cluster_sample_id, "_")[[1]]
  df$sample_id  <- tstrsplit(df$cluster_sample_id, "_")[[2]]
  # Retrieve cell count information for this cluster from global cell count table
  idx <- which(colnames(t) == unique(df$cluster_id))
  cell_counts <- t[, idx]
  # Remove samples with zero cell contributing to the cluster
  cell_counts <- cell_counts[cell_counts > 0]
  # Match order of cell_counts and sample_ids
  sample_order <- match(df$sample_id, names(cell_counts))
  cell_counts <- cell_counts[sample_order]
  # Append cell_counts to data frame
  df$cell_count <- cell_counts
  # Join data frame (capturing metadata specific to cluster) to generic metadata
  df <- plyr::join(df, metadata, 
                   by = intersect(names(df), names(metadata)))
  # Update rownames of metadata to match colnames of count matrix, as needed later for DE
  rownames(df) <- df$cluster_sample_id
  # Store complete metadata for cluster i in list
  metadata_ls[[i]] <- df
  names(metadata_ls)[i] <- unique(df$cluster_id)
  
}

# Explore the different components of the list
str(metadata_ls)

# Summary
# 1) Cell types are stored in a vector called `cluster_names`
# 2) Cell types are also stored in names(`counts_ls`) and names(`metdata_ls`), which should match each other

# Step 6: DEA using DESeq2 -----------------------------------------------------
# Create a DESeq2 object
cluster_names
# Double-check that both lists have same names
all(names(counts_ls) == names(metadata_ls))

# Use Macrophage as an example 
idx <- which(names(counts_ls) == "Macrophage")
cluster_counts <- counts_ls[[idx]]
cluster_metadata <- metadata_ls[[idx]]
# Check contents of extracted objects
cluster_counts[1:6, 1:6]
head(cluster_metadata)
# Check matching of matrix columns and metadata rows
all(colnames(cluster_counts) == rownames(cluster_metadata))

# Create object
pseudobulk_dds <- DESeqDataSetFromMatrix(cluster_counts, 
                              colData = cluster_metadata, 
                              design = ~ group_id)

# Export data 
save(pseudobulk_dds, file="pseudobulk_dds.RData")


sessionInfo()
