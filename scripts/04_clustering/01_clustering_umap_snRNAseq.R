# 01_clustering_umap_snRNAseq.R
# ------------------------------------------------------------
# Clustering and UMAP visualization for snRNA-seq
#
# This step operates on an integrated Seurat object
# (e.g. Harmony-corrected) and performs:
#   - neighbor graph construction
#   - clustering
#   - UMAP embedding
#
# Author: Naghmeh Rezaei
# Repository: snRNAseq-workflows
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
})

# ---------------------------
# 1. Load integrated object
# ---------------------------

input_path <- "results/seurat_harmony_integrated.rds"
seu <- readRDS(input_path)

# ---------------------------
# 2. Dimensionality choice
# ---------------------------

# NOTE:
# Number of dimensions should be guided by:
# - ElbowPlot
# - cumulative variance
# - dataset complexity
#
# Typical range for snRNA-seq: 20–40 PCs

dims_use <- 1:30

# ---------------------------
# 3. Neighbor graph
# ---------------------------

seu <- FindNeighbors(
  seu,
  reduction = "harmony",
  dims = dims_use
)

# ---------------------------
# 4. Clustering
# ---------------------------

# NOTE:
# Resolution controls cluster granularity.
# Typical starting range: 0.3–1.0

seu <- FindClusters(
  seu,
  resolution = 0.6
)

# ---------------------------
# 5. UMAP
# ---------------------------

seu <- RunUMAP(
  seu,
  reduction = "harmony",
  dims = dims_use
)

# ---------------------------
# 6. Save clustered object
# ---------------------------

output_path <- "results/seurat_harmony_clustered.rds"
saveRDS(seu, file = output_path)

cat("Clustered Seurat object saved to:", output_path, "\n")
