# 01_normalization_snRNAseq.R
# ------------------------------------------------------------
# Normalization strategies for snRNA-seq
#
# This workflow supports TWO commonly used normalization approaches:
#
# 1) LogNormalize
#    - Lightweight
#    - Preserves raw count structure
#    - Often used for exploratory analyses or small datasets
#
# 2) SCTransform
#    - Variance-stabilizing normalization
#    - Regresses technical effects
#    - Preferred for multi-sample integration
#
# Choice depends on dataset size, batch structure, and downstream goals.
#
# Author: Naghmeh Rezaei
# Repository: snRNAseq-workflows
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
})

# ---------------------------
# 1. Load post-doublet object
# ---------------------------

input_path <- "results/seurat_post_doublet.rds"
seu <- readRDS(input_path)

# ---------------------------
# 2. Strategy A: LogNormalize
# ---------------------------

# NOTE:
# LogNormalize is often sufficient for:
# - single-sample analyses
# - exploratory clustering
# - datasets with minimal batch effects

seu_log <- NormalizeData(
  seu,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

seu_log <- FindVariableFeatures(
  seu_log,
  selection.method = "vst",
  nfeatures = 2000
)

# ---------------------------
# 3. Strategy B: SCTransform
# ---------------------------

# NOTE:
# SCTransform is preferred for:
# - multi-sample datasets
# - integration workflows
# - robust downstream comparisons
#
# percent.mt regression is optional for snRNA-seq

# Example (not run by default):
#
# seu_sct <- SCTransform(
#   seu,
#   vars.to.regress = NULL,
#   verbose = FALSE
# )

# ---------------------------
# 4. Save outputs
# ---------------------------

saveRDS(seu_log, file = "results/seurat_log_normalized.rds")
cat("Log-normalized object saved\n")

# Uncomment if SCTransform is used
# saveRDS(seu_sct, file = "results/seurat_sct_normalized.rds")
