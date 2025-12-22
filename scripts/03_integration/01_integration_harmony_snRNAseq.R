# 01_integration_harmony_snRNAseq.R
# ------------------------------------------------------------
# Dataset integration for snRNA-seq using Harmony
#
# This script integrates multiple snRNA-seq samples while
# preserving biological variation and correcting technical
# batch effects.
#
# Typical use cases:
# - Multiple samples / mice
# - Multiple experimental conditions
# - snRNA-seq with prep- or batch-related variability
#
# Author: Naghmeh Rezaei
# Repository: snRNAseq-workflows
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
})

# ---------------------------
# 1. Load normalized object
# ---------------------------

# Use LogNormalize output by default
# (SCTransform-based integration can be added if needed)
input_path <- "results/seurat_log_normalized.rds"
seu <- readRDS(input_path)

# ---------------------------
# 2. Metadata sanity checks
# ---------------------------

# REQUIRED metadata columns (typical)
# These should already exist in your Seurat object
required_meta <- c("sample", "condition")

missing_meta <- setdiff(required_meta, colnames(seu@meta.data))
if (length(missing_meta) > 0) {
  stop(
    "Missing required metadata columns: ",
    paste(missing_meta, collapse = ", ")
  )
}

# ---------------------------
# 3. Scaling & PCA
# ---------------------------

# Scale only variable features
seu <- ScaleData(
  seu,
  features = VariableFeatures(seu),
  verbose = FALSE
)

seu <- RunPCA(
  seu,
  features = VariableFeatures(seu),
  npcs = 50,
  verbose = FALSE
)

# ---------------------------
# 4. Harmony integration
# ---------------------------

# batch variable typically used in real projects:
# - sample
# - mouse_id
# - sequencing_run
#
# Here we use "sample" as a generic batch covariate

seu <- RunHarmony(
  object = seu,
  group.by.vars = "sample",
  reduction = "pca",
  assay.use = "RNA",
  verbose = TRUE
)

# ---------------------------
# 5. Save integrated object
# ---------------------------

output_path <- "results/seurat_harmony_integrated.rds"
saveRDS(seu, file = output_path)

cat("Harmony-integrated Seurat object saved to:", output_path, "\n")
