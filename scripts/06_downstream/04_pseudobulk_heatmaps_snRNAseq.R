# 04_pseudobulk_heatmaps_snRNAseq.R
# ------------------------------------------------------------
# Pseudobulk heatmap visualization for snRNA-seq
#
# This script generates heatmaps of top DE genes
# from pseudobulk cell-type–specific analyses.
#
# Heatmaps are used for:
# - visual validation of DE
# - comparison across samples / conditions
# - publication-quality figures
#
# Author: Naghmeh Rezaei
# Repository: snRNAseq-workflows
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ComplexHeatmap)
  library(circlize)
})

# ---------------------------
# 1. Load annotated object
# ---------------------------

seu <- readRDS("results/seurat_annotated.rds")

# ---------------------------
# 2. Select cell type to plot
# ---------------------------

# NOTE:
# Heatmaps are usually generated per cell type.
# Change this value as needed.

cell_type_of_interest <- "Hepatocyte"

# ---------------------------
# 3. Load DE results
# ---------------------------

de_file <- paste0(
  "results/celltype_DE/DE_",
  gsub(" ", "_", cell_type_of_interest),
  ".csv"
)

de <- read.csv(de_file)

# Select top DE genes
top_genes <- de %>%
  filter(FDR < 0.05) %>%
  arrange(desc(abs(logFC))) %>%
  slice_head(n = 50) %>%
  pull(gene)

# ---------------------------
# 4. Build pseudobulk matrix
# ---------------------------

seu_ct <- subset(seu, subset = cell_type == cell_type_of_interest)

pb_mat <- AggregateExpression(
  seu_ct,
  group.by = "sample",
  assays = "RNA",
  slot = "counts",
  return.seurat = FALSE
)$RNA

# Log-transform
pb_log <- log2(pb_mat + 1)

# Subset to DE genes
pb_log <- pb_log[top_genes, , drop = FALSE]

# ---------------------------
# 5. Sample annotation
# ---------------------------

sample_meta <- seu_ct@meta.data %>%
  distinct(sample, condition) %>%
  arrange(match(sample, colnames(pb_log)))

ha_col <- HeatmapAnnotation(
  condition = sample_meta$condition,
  col = list(
    condition = c(
      Control = "#4DAF4A",
      Treatment = "#E41A1C"
    )
  )
)

# ---------------------------
# 6. Heatmap
# ---------------------------

ht <- Heatmap(
  pb_log,
  name = "log2(counts + 1)",
  top_annotation = ha_col,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_title = paste(cell_type_of_interest, "pseudobulk expression"),
  use_raster = TRUE
)

# ---------------------------
# 7. Save figure
# ---------------------------

pdf(
  file = paste0(
    "figures/pseudobulk_heatmap_",
    gsub(" ", "_", cell_type_of_interest),
    ".pdf"
  ),
  width = 8,
  height = 10
)

draw(ht)
dev.off()

cat("Pseudobulk heatmap saved\n")
