# scRNA-seq downstream analysis
# Part I: QC

# Load packages ----------------------------------------------------------------
library(magrittr) 
library(dplyr) 
library(DoubletFinder)
library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(dplyr)
# library(DESeq2)

# Set dir ----------------------------------------------------------------------
count_dir <- "/proj/"
# dir.create("SN_data/", recursive = TRUE)
# dir.create("SN_plots/", recursive = TRUE)
data_dir <- "/proj/SN_data/"
plots_dir <- "/proj/SN_plots/"

# Load data --------------------------------------------------------------------
# Set sample info
setwd(count_dir)
sample_info <- data.frame(ctrl = "ctrl", exp = "exp")
# Read in 10X data for a single sample (output is a sparse matrix)
ctrl_counts <- Read10X(data.dir = "./sample-ctrl/outs/raw_feature_bc_matrix/")
exp_counts <- Read10X(data.dir = "./sample-exp/outs/raw_feature_bc_matrix/")
# Turn count matrix into a Seurat object (output is a Seurat object)
#   cells with less than 100 genes detected are not considered for analysis
ctrl <- CreateSeuratObject(counts = ctrl_counts, min.features = 100, project = sample_info$ctrl)
exp <- CreateSeuratObject(counts = exp_counts, min.features = 100, project = sample_info$exp)
# Check 
ctrl 
head(ctrl@meta.data)
exp
head(exp@meta.data)
# 1) orig.ident: this often contains the sample identity if known, but will default to "SeuratProject"
# 2) nCount_RNA: number of UMIs per cell
# 3) nFeature_RNA: number of genes detected per cell

# Read in multiple samples with a for loop
# # Create each individual Seurat object for every sample
# for (experiment in c(sample.info$ctrl, sample.info$exp)){
#  seurat_count =  Read10X(data.dir = paste0(data_dir, experiment, "/raw_feature_bc_matrix/"))
#  seurat_obj = CreateSeuratObject(counts = seurat_count, min.features = 100, project = experiment)
#  assign(experiment,seurat_obj)
# }

# Create merged seurat object --------------------------------------------------
merged_seurat <- merge(x = ctrl, y = exp, add.cell.id = c('ctrl', 'exp'))
View(merged_seurat@meta.data)

# Calculate log10GenesPerUMI ---------------------------------------------------
#   log10 transform the result for better comparison between samples
#   add number of genes per UMI for each cell to metadata
merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA) 
head(merged_seurat)

# Compute uninterested gene % --------------------------------------------------
merged_seurat$percent.mt <- PercentageFeatureSet(object = merged_seurat, pattern = "^MT-")
merged_seurat$percent.rb <- PercentageFeatureSet(object = merged_seurat, pattern = "^RP[SL]")
merged_seurat$percent.hb <- PercentageFeatureSet(object = merged_seurat, pattern = "^HB[^(P)]")

# Create metadata df -----------------------------------------------------------
metadata <-  merged_seurat@meta.data
# Add a new column with cell IDs 
metadata$cells <- rownames(metadata)
# Add a new column with sample names for each of the cells based on the cell prefix
metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "_ctrl"))] <- "ctrl"
metadata$sample[which(str_detect(metadata$cells, "_exp"))] <- "exp"
View(metadata)

# Re-attach the updated metadata back to Seurat object
merged_seurat@meta.data <- metadata

# Export .RData object 
setwd(data_dir)
save(merged_seurat, file="./merged_seurat.RData")

# QC plots ---------------------------------------------------------------------
setwd(plots_dir)
# figure 1_cells_per_sample ----------------------------------------------------
# Visualize the number of cell counts per sample 
metadata %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("nCells") 

# figure 2_nCount_RNA ----------------------------------------------------------
# visualize the number UMIs/transcripts per cell 
metadata %>% 
  ggplot(aes(color=sample, x=nCount_RNA, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)

# figure 3_genes_per_cell ------------------------------------------------------
# visualize the distribution of genes detected per cell via histogram
metadata %>% 
  ggplot(aes(color=sample, x=nFeature_RNA, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300) +
  xlim (100, 15000)

# figure 4_cells_vs_genes ------------------------------------------------------
# visualize the distribution of genes detected per cell via boxplot
metadata %>% 
  ggplot(aes(x=sample, y=log10(nFeature_RNA), fill=sample)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("nCells vs nFeature_RNAs") 

# figure 5_corr_genes_UMIs_mt --------------------------------------------------
# Visualize the correlation between genes detected and number of UMIs and determine 
#   whether strong presence of cells with low numbers of genes/UMIs
metadata %>% 
  ggplot(aes(x=nFeature_RNA, y=nCount_RNA, color=percent.mt)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~sample)

# figure 6_corr_genes_UMIs_rb --------------------------------------------------
metadata %>% 
  ggplot(aes(x=nFeature_RNA, y=nCount_RNA, color=percent.rb)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~sample)

# figure 7_QC_violin_plot ------------------------------------------------------
VlnPlot(merged_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt",
                                    "percent.rb", "percent.hb"), ncol = 5)

# figure 8_raw_feature_scatter -------------------------------------------------
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(merged_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(merged_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1+plot2

# figure 9_raw_complexity ------------------------------------------------------
# another way of visualizing complexity
# visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)

# figure 10_counts_each_gene ---------------------------------------------------
par(mar = c(4, 8, 2, 1))
C <- merged_seurat@assays$RNA@counts
C <- Matrix::t(Matrix::t(C)/Matrix::colSums(C)) * 100
most_expressed <- order(Matrix::rowSums(C), decreasing = T)[20:1]
boxplot(as.matrix(t(C[most_expressed, ])), cex = 0.1, las = 1, xlab = "% total count per cell",
        col = (scales::hue_pal())(20)[20:1], horizontal = TRUE)

rm(C)

# Cell-level filtering ---------------------------------------------------------
# setwd(data_dir)
# load("merged_seurat.RData")
filtered_seurat <- subset(x = merged_seurat, 
                         subset= 
                           (nFeature_RNA >= 500) &
                           (nCount_RNA >= 250) & 
                           (log10GenesPerUMI > 0.80) & 
                           (percent.mt < 20) &
                           (percent.rb > 5))
# Check stats
# Pre-filtering cell numbers
table(merged_seurat@meta.data$sample) 
# Post-filtering cell numbers
table(filtered_seurat@meta.data$sample) 

# Gene-level filtering ---------------------------------------------------------
filtered_seurat <- filtered_seurat[!grepl(c("MALAT1", "^MT-", "^RP[SL]", "^HB[^(P)]"), 
                                          rownames(filtered_seurat)), ]
# Check stats
# Pre-filtering gene numbers
dim(merged_seurat) 
# Post-filtering gene numbers
dim(filtered_seurat) 

# Output a logical vector for every gene on whether the more than zero counts per cell
counts <- GetAssayData(object = filtered_seurat, assay = "RNA", slot = "counts")
nonzero <- counts > 0
# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10
filtered_counts <- counts[keep_genes, ]

# Reassign to filtered Seurat object
filtered_seurat <- CreateSeuratObject(counts = filtered_counts, meta.data = filtered_seurat@meta.data)
filtered_seurat@active.ident <- merged_seurat@active.ident

# Check stats
# Pre-filtering gene numbers
dim(merged_seurat) 
# Post-filtering gene numbers
dim(filtered_seurat) 

# Update metadata after filtering
metadata_clean <- filtered_seurat@meta.data

# Export .RData object 
setwd(data_dir)
save(filtered_seurat, file="filtered_seurat.RData")

# Post-filtering QC plots ------------------------------------------------------
setwd(plots_dir)
# figure 11_filtered_counts_each_gene ------------------------------------------
par(mar = c(4, 8, 2, 1))
C <- filtered_seurat@assays$RNA@counts
C <- Matrix::t(Matrix::t(C)/Matrix::colSums(C)) * 100
most_expressed <- order(Matrix::rowSums(C), decreasing = T)[20:1]
boxplot(as.matrix(t(C[most_expressed, ])), cex = 0.1, las = 1, xlab = "% total count per cell",
        col = (scales::hue_pal())(20)[20:1], horizontal = TRUE)

rm(C)
# figure 12_filtered_violin_plot -----------------------------------------------
VlnPlot(filtered_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt",
                                      "percent.rb", "percent.hb"), ncol = 5)

# figure 13_filtered_feature_scatter--------------------------------------------
plot1 <- FeatureScatter(filtered_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(filtered_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1+plot2

# figure 14_filtered_complexity-------------------------------------------------
metadata_clean %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)

# DoubletFinder ----------------------------------------------------------------
suppressMessages(require(DoubletFinder))
filtered_seurat_Doublet = NormalizeData(filtered_seurat, verbose = F)
filtered_seurat_Doublet = FindVariableFeatures(filtered_seurat_Doublet, selection.method = "vst", nfeatures = 2000, verbose = F)
filtered_seurat_Doublet = ScaleData(filtered_seurat_Doublet, verbose = F)
filtered_seurat_Doublet = RunPCA(filtered_seurat_Doublet, verbose = F, npcs = 20)
filtered_seurat_Doublet = RunUMAP(filtered_seurat_Doublet, dims = 1:10, verbose = F)
# Define the expected number of doublet cellscells.
nExp <- round(ncol(filtered_seurat_Doublet) * 0.076)  # expect 5% doublets 
nExp 
# Use the table here to estimate: https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_01_qc.html

# doubleFinder function code 
# fixed bug from source code (https://github.com/chris-mcginnis-ucsf/DoubletFinder)
doubletFinder_sn <- function(seu, PCs, pN = 0.25, pK, nExp, reuse.pANN = FALSE, sct = FALSE, annotations = NULL) {
  require(Seurat); require(fields); require(KernSmooth)
  
  ## Generate new list of doublet classificatons from existing pANN vector to save time
  if (reuse.pANN != FALSE ) {
    pANN.old <- seu@meta.data[ , reuse.pANN]
    classifications <- rep("Singlet", length(pANN.old))
    classifications[order(pANN.old, decreasing=TRUE)[1:nExp]] <- "Doublet"
    seu@meta.data[, paste("DF.classifications",pN,pK,nExp,sep="_")] <- classifications
    return(seu)
  }
  
  if (reuse.pANN == FALSE) {
    ## Make merged real-artifical data
    real.cells <- rownames(seu@meta.data)
    data <- seu@assays$RNA@counts[, real.cells]
    n_real.cells <- length(real.cells)
    n_doublets <- round(n_real.cells/(1 - pN) - n_real.cells)
    print(paste("Creating",n_doublets,"artificial doublets...",sep=" "))
    real.cells1 <- sample(real.cells, n_doublets, replace = TRUE)
    real.cells2 <- sample(real.cells, n_doublets, replace = TRUE)
    doublets <- (data[, real.cells1] + data[, real.cells2])/2
    colnames(doublets) <- paste("X", 1:n_doublets, sep = "")
    data_wdoublets <- cbind(data, doublets)
    # Keep track of the types of the simulated doublets
    if(!is.null(annotations)){
      stopifnot(typeof(annotations)=="character")
      stopifnot(length(annotations)==length(Cells(seu)))
      stopifnot(!any(is.na(annotations)))
      annotations <- factor(annotations)
      names(annotations) <- Cells(seu)
      doublet_types1 <- annotations[real.cells1]
      doublet_types2 <- annotations[real.cells2]
    }
    ## Store important pre-processing information
    orig.commands <- seu@commands
    
    ## Pre-process Seurat object
    if (sct == FALSE) {
      print("Creating Seurat object...")
      seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)
      
      print("Normalizing Seurat object...")
      seu_wdoublets <- NormalizeData(seu_wdoublets,
                                     normalization.method = orig.commands$NormalizeData.RNA@params$normalization.method,
                                     scale.factor = orig.commands$NormalizeData.RNA@params$scale.factor,
                                     margin = orig.commands$NormalizeData.RNA@params$margin)
      
      print("Finding variable genes...")
      seu_wdoublets <- FindVariableFeatures(seu_wdoublets,
                                            selection.method = orig.commands$FindVariableFeatures.RNA$selection.method,
                                            loess.span = orig.commands$FindVariableFeatures.RNA$loess.span,
                                            clip.max = orig.commands$FindVariableFeatures.RNA$clip.max,
                                            mean.function = orig.commands$FindVariableFeatures.RNA$mean.function,
                                            dispersion.function = orig.commands$FindVariableFeatures.RNA$dispersion.function,
                                            num.bin = orig.commands$FindVariableFeatures.RNA$num.bin,
                                            binning.method = orig.commands$FindVariableFeatures.RNA$binning.method,
                                            nfeatures = orig.commands$FindVariableFeatures.RNA$nfeatures,
                                            mean.cutoff = orig.commands$FindVariableFeatures.RNA$mean.cutoff,
                                            dispersion.cutoff = orig.commands$FindVariableFeatures.RNA$dispersion.cutoff)
      
      print("Scaling data...")
      seu_wdoublets <- ScaleData(seu_wdoublets,
                                 features = orig.commands$ScaleData.RNA$features,
                                 model.use = orig.commands$ScaleData.RNA$model.use,
                                 do.scale = orig.commands$ScaleData.RNA$do.scale,
                                 do.center = orig.commands$ScaleData.RNA$do.center,
                                 scale.max = orig.commands$ScaleData.RNA$scale.max,
                                 block.size = orig.commands$ScaleData.RNA$block.size,
                                 min.cells.to.block = orig.commands$ScaleData.RNA$min.cells.to.block)
      
      print("Running PCA...")
      seu_wdoublets <- RunPCA(seu_wdoublets,
                              features = orig.commands$ScaleData.RNA$features,
                              npcs = length(PCs),
                              rev.pca =  orig.commands$RunPCA.RNA$rev.pca,
                              weight.by.var = orig.commands$RunPCA.RNA$weight.by.var,
                              verbose=FALSE)
      pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[ , PCs]
      cell.names <- rownames(seu_wdoublets@meta.data)
      nCells <- length(cell.names)
      rm(seu_wdoublets); gc() # Free up memory
    }
    
    if (sct == TRUE) {
      require(sctransform)
      print("Creating Seurat object...")
      seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)
      
      print("Running SCTransform...")
      seu_wdoublets <- SCTransform(seu_wdoublets)
      
      print("Running PCA...")
      seu_wdoublets <- RunPCA(seu_wdoublets, npcs = length(PCs))
      pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[ , PCs]
      cell.names <- rownames(seu_wdoublets@meta.data)
      nCells <- length(cell.names)
      rm(seu_wdoublets); gc()
    }
    
    ## Compute PC distance matrix
    print("Calculating PC distance matrix...")
    dist.mat <- fields::rdist(pca.coord)
    
    ## Compute pANN
    print("Computing pANN...")
    pANN <- as.data.frame(matrix(0L, nrow = n_real.cells, ncol = 1))
    if(!is.null(annotations)){
      neighbor_types <- as.data.frame(matrix(0L, nrow = n_real.cells, ncol = length(levels(doublet_types1))))
    }
    rownames(pANN) <- real.cells
    colnames(pANN) <- "pANN"
    k <- round(nCells * pK)
    for (i in 1:n_real.cells) {
      neighbors <- order(dist.mat[, i])
      neighbors <- neighbors[2:(k + 1)]
      pANN$pANN[i] <- length(which(neighbors > n_real.cells))/k
      if(!is.null(annotations)){
        for(ct in unique(annotations)){
          neighbors_that_are_doublets = neighbors[neighbors>n_real.cells]
          if(length(neighbors_that_are_doublets) > 0){
            neighbor_types[i,] <-
              table( doublet_types1[neighbors_that_are_doublets - n_real.cells] ) +
              table( doublet_types2[neighbors_that_are_doublets - n_real.cells] )
            neighbor_types[i,] <- neighbor_types[i,] / sum( neighbor_types[i,] )
          } else {
            neighbor_types[i,] <- NA
          }
        }
      }
    }
    print("Classifying doublets..")
    classifications <- rep("Singlet",n_real.cells)
    classifications[order(pANN$pANN[1:n_real.cells], decreasing=TRUE)[1:nExp]] <- "Doublet"
    seu@meta.data[, paste("pANN",pN,pK,nExp,sep="_")] <- pANN[rownames(seu@meta.data), 1]
    seu@meta.data[, paste("DF.classifications",pN,pK,nExp,sep="_")] <- classifications
    if(!is.null(annotations)){
      colnames(neighbor_types) = levels(doublet_types1)
      for(ct in levels(doublet_types1)){
        seu@meta.data[, paste("DF.doublet.contributors",pN,pK,nExp,ct,sep="_")] <- neighbor_types[,ct]
      }
    }
    return(seu)
  }
}
# Identify doublets
filtered_seurat_rmDoublet <- doubletFinder_sn(filtered_seurat_Doublet, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10,
                                              reuse.pANN = FALSE, sct = FALSE)
# Name of the DF prediction can change, so extract the correct column name.
DF.name = colnames(filtered_seurat_rmDoublet@meta.data)[grepl("DF.classification", colnames(filtered_seurat_rmDoublet@meta.data))]

# figure 15_doubletFinder_0.76 -------------------------------------------------
cowplot::plot_grid(ncol = 2, DimPlot(filtered_seurat_rmDoublet, group.by = "orig.ident" ) + NoAxes(),
                   DimPlot(filtered_seurat_rmDoublet, group.by = DF.name, cols = c("red", "lightgreen")) + NoAxes())

# figure 16_doublet_features_0.76 ----------------------------------------------
VlnPlot(filtered_seurat_rmDoublet, features = "nFeature_RNA", group.by = DF.name, pt.size = 0.1, cols = c("red", "lightgreen"))

# Remove all predicted doublets
filtered_seurat_rmDoublet = filtered_seurat_rmDoublet[, filtered_seurat_rmDoublet@meta.data[, DF.name] == "Singlet"]
nExp 

# Export .RData object 
setwd(data_dir)
save(filtered_seurat_rmDoublet, file="filtered_seurat_rmDoublet.RData")


sessionInfo()