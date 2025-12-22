# 01_qc_seurat_snRNAseq.R
# ------------------------------------------------------------
# Single-nucleus RNA-seq (snRNA-seq) quality control workflow
#
# Steps:
#   1. Load raw Seurat object
#   2. Compute QC metrics
#   3. Filter low-quality nuclei
#   4. (Optional) Doublet detection placeholder
#
# NOTE:
# - Paths and sample names are intentionally generic
# - This script is adapted from real research analyses
#
# Author: Naghmeh Rezaei
# Repository: snRNAseq-workflows
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

# ---------------------------
# 1. Load data
# ---------------------------

# Example placeholder path
# Replace with your own Seurat object
seurat_path <- "data/seurat_raw_example.rds"

# This file is NOT included in the repo
# It demonstrates expected input format
seu <- readRDS(seurat_path)

# ---------------------------
# 2. Compute QC metrics
# ---------------------------

# Percentage of mitochondrial genes
# (snRNA-seq often has lower mitochondrial content than scRNA-seq)
seu[["percent.mt"]] <- PercentageFeatureSet(
  seu,
  pattern = "^mt-"
)

# Basic QC inspection
qc_metrics <- seu@meta.data %>%
  select(nFeature_RNA, nCount_RNA, percent.mt)

print(summary(qc_metrics))

# ---------------------------
# 3. Filter nuclei
# ---------------------------

# Thresholds should be dataset-specific
# These are typical starting values for snRNA-seq
seu_qc <- subset(
  seu,
  subset =
    nFeature_RNA > 200 &
    nFeature_RNA < 6000 &
    percent.mt < 5
)

cat(
  "Nuclei before QC:", ncol(seu), "\n",
  "Nuclei after QC :", ncol(seu_qc), "\n"
)

# ---------------------------
# 4. Doublet detection (placeholder)
# ---------------------------

# In real analyses, this step may use:
# - Scrublet (Python)
# - scds (R)
#
# This repo includes a placeholder to indicate where
# doublet detection would be integrated.

# Example:
# seu_qc$doublet_score <- NA

# ---------------------------
# 5. Save QC-filtered object
# ---------------------------

output_path <- "results/seurat_qc_filtered.rds"
saveRDS(seu_qc, file = output_path)

cat("QC-filtered Seurat object saved to:", output_path, "\n")
