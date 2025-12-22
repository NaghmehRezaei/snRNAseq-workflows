# 02_doublet_handling_snRNAseq.R
# ------------------------------------------------------------
# Doublet detection strategies for snRNA-seq
#
# This workflow supports TWO commonly used approaches:
#
# 1) CellBender (preferred when raw Cell Ranger outputs are available)
#    - Removes ambient RNA
#    - Models technical background
#    - Often sufficient without additional doublet filtering
#
# 2) scds (used when starting from processed Seurat objects)
#    - cxds / bcds / hybrid scoring
#    - Applied post-QC
#
# Choice of method depends on dataset origin and preprocessing history.
#
# Author: Naghmeh Rezaei
# Repository: snRNAseq-workflows
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

# ---------------------------
# 1. Load QC-filtered object
# ---------------------------

qc_path <- "results/seurat_qc_filtered.rds"
seu_qc <- readRDS(qc_path)

# ---------------------------
# 2. Strategy A: CellBender
# ---------------------------

# NOTE:
# If CellBender was run upstream (recommended when raw matrices exist),
# the resulting filtered matrix is assumed to be largely free of doublets.
#
# In this case, this step serves as documentation only,
# and no additional filtering is applied here.

# Example (documentation only):
# seu_qc <- Read10X_h5("cellbender_filtered.h5")

cat("CellBender strategy: no additional doublet filtering applied\n")

# ---------------------------
# 3. Strategy B: scds (R-based)
# ---------------------------

# NOTE:
# Use scds when starting from an existing Seurat object
# and CellBender was not applied.

# Example (not executed by default):
#
# library(scds)
# sce <- as.SingleCellExperiment(seu_qc)
# sce <- cxds(sce)
# sce <- bcds(sce)
# sce <- cxds_bcds_hybrid(sce)
#
# seu_qc$doublet_score <- colData(sce)$hybrid_score
#
# seu_qc <- subset(
#   seu_qc,
#   subset = doublet_score < quantile(doublet_score, 0.95)
# )

cat("scds strategy documented (execution optional)\n")

# ---------------------------
# 4. Save post-doublet object
# ---------------------------

output_path <- "results/seurat_post_doublet.rds"
saveRDS(seu_qc, file = output_path)

cat("Post-doublet Seurat object saved to:", output_path, "\n")
