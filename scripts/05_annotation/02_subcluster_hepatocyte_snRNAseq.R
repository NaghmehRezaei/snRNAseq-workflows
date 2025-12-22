# 02_subcluster_hepatocyte_snRNAseq.R
# ------------------------------------------------------------
# Hepatocyte subclustering workflow for snRNA-seq data
#
# Description:
# This script performs cell-type–specific subclustering of
# hepatocytes from an already QC-filtered, integrated, and
# annotated snRNA-seq Seurat object.
#
# Rationale:
# - Broad clustering masks hepatocyte heterogeneity
# - Subclustering reveals zonation, metabolic states,
#   stress responses, and disease-associated programs
#
# Assumptions:
# - QC and doublet removal already performed
# - Batch correction / integration already completed
# - Metadata column `cell_type` exists
#
# Notes:
# - Re-normalization is required for valid subclustering
# - Conservative clustering resolution is recommended
#
# Author: Naghmeh Rezaei
# Repository: snRNAseq-workflows
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

# ---------------------------
# 1. Load annotated Seurat object
# ---------------------------

input_path <- "results/seurat_annotated.rds"
seu <- readRDS(input_path)

stopifnot("cell_type" %in% colnames(seu@meta.data))

# ---------------------------
# 2. Subset hepatocytes
# ---------------------------

hep <- subset(seu, subset = cell_type == "Hepatocyte")

cat("Number of hepatocyte nuclei:", ncol(hep), "\n")

# ---------------------------
# 3. Re-normalize within hepatocytes
# ---------------------------
# Required to avoid dominance of global variance structure

hep <- NormalizeData(hep)

hep <- FindVariableFeatures(
  hep,
  selection.method = "vst",
  nfeatures = 3000
)

hep <- ScaleData(
  hep,
  vars.to.regress = c("nCount_RNA", "percent.mt")
)

hep <- RunPCA(
  hep,
  npcs = 50,
  verbose = FALSE
)

# ---------------------------
# 4. Determine dimensionality
# ---------------------------
# Inspect elbow plot manually when running interactively

ElbowPlot(hep, ndims = 50)

# Typical range for hepatocyte subclustering
pcs_use <- 1:25

# ---------------------------
# 5. Subclustering
# ---------------------------

hep <- FindNeighbors(hep, dims = pcs_use)

hep <- FindClusters(
  hep,
  resolution = 0.3  # conservative to avoid overclustering
)

hep <- RunUMAP(
  hep,
  dims = pcs_use,
  reduction.name = "umap_hepatocyte"
)

# ---------------------------
# 6. Marker discovery (QC + interpretation)
# ---------------------------

hep_markers <- FindAllMarkers(
  hep,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

write.csv(
  hep_markers,
  file = "results/hepatocyte_subcluster_markers.csv",
  row.names = FALSE
)

# ---------------------------
# 7. Save outputs
# ---------------------------

saveRDS(
  hep,
  file = "results/hepatocyte_subclustered.rds"
)

cat("Hepatocyte subclustering completed successfully\n")
