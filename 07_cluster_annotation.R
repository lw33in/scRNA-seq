# scRNA-seq downstream analysis

# Use clustering output in excel format to identify origin of cells and calculate percentages in each cluster

# Load packages ----------------------------------------------------------------
# install.packages(c("readxl", "dplyr"))
library(readxl)
# library(dplyr)

# Load data --------------------------------------------------------------------
setwd("/scRNA/output/")
exp1 <- read_excel("exp1.clusters_info.xlsx")

# Calculation ------------------------------------------------------------------
head(exp1.clusters_info.xlsx)
table <- table(exp1$orig.ident, exp1$seurat_clusters) 
cluster_stats <- data.frame(exp = table[1,], ctrl = table[2,])

# Calculate cell origin % in each cluster
cluster_stats$cluster <- row.names(cluster_stats)
cluster_stats$exp_perc  <- paste0(round(with(cluster_stats, exp  / (exp + ctrl) * 100), 2), "%" )
cluster_stats$ctrl_perc <- paste0(round(with(cluster_stats, ctrl / (exp + ctrl) * 100), 2), "%" )
cluster_stats <- cluster_stats[, c(ncol(cluster_stats), 1:(ncol(cluster_stats)-1))]

write.table(cluster_stats, 
            file = "cluster_stats.txt", sep = "\t", row.names = FALSE, quote = FALSE)


sessionInfo()
