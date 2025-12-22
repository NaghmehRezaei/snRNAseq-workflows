# 01_annotation_marker_based_snRNAseq.R
# ------------------------------------------------------------
# Marker-based cell type annotation for snRNA-seq
#
# This script assigns biological cell type labels to clusters
# using canonical marker genes and expert-driven inspection.
#
# Annotation is an iterative process and should be validated
# with known biology and visualization.
#
# Author: Naghmeh Rezaei
# Repository: snRNAseq-workflows
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

# ---------------------------
# 1. Load clustered object
# ---------------------------

input_path <- "results/seurat_harmony_clustered.rds"
seu <- readRDS(input_path)

# ---------------------------
# 2. Identify cluster markers
# ---------------------------

markers <- FindAllMarkers(
  seu,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

# Save marker table for manual inspection
write.csv(
  markers,
  file = "results/cluster_markers_all.csv",
  row.names = FALSE
)

# ---------------------------
# 3. Define canonical markers
# ---------------------------

# Example marker lists (to be adapted per tissue)
marker_reference <- list(
  Hepatocyte = c("Alb", "Ttr", "Apoa1"),
  Kupffer = c("Lyz2", "Clec4f", "Adgre1"),
  LSEC = c("Kdr", "Pecam1", "Rspo3"),
  Stellate = c("Col1a1", "Des", "Lrat"),
  Cholangiocyte = c("Krt19", "Krt7"),
  T_cell = c("Cd3d", "Cd3e"),
  B_cell = c("Cd79a", "Ms4a1")
)

# ---------------------------
# 4. Visual inspection
# ---------------------------

# DotPlot for marker validation
DotPlot(
  seu,
  features = marker_reference
) + RotatedAxis()

# ---------------------------
# 5. Assign cluster annotations
# ---------------------------

# NOTE:
# Cluster-to-cell type mapping is dataset-specific.
# This example mapping must be updated after inspection.

cluster_to_celltype <- c(
  "0" = "Hepatocyte",
  "1" = "Kupffer",
  "2" = "LSEC",
  "3" = "Stellate",
  "4" = "T_cell",
  "5" = "B_cell"
)

seu$cell_type <- plyr::mapvalues(
  seu$seurat_clusters,
  from = names(cluster_to_celltype),
  to = unname(cluster_to_celltype)
)

# ---------------------------
# 6. Save annotated object
# ---------------------------

output_path <- "results/seurat_annotated.rds"
saveRDS(seu, file = output_path)

cat("Annotated Seurat object saved to:", output_path, "\n")
