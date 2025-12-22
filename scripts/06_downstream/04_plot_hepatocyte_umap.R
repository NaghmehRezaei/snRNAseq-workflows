# 04_plot_hepatocyte_umap.R
# ------------------------------------------------------------
# UMAP visualization of hepatocyte subclusters
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

hep <- readRDS("results/hepatocyte_subcluster_annotated.rds")

p <- DimPlot(
  hep,
  reduction = "umap_hepatocyte",
  group.by = "hep_subtype",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("Hepatocyte subclusters (snRNA-seq)") +
  theme_minimal(base_size = 14)

ggsave(
  filename = "figures/hepatocyte_subclusters_umap.png",
  plot = p,
  width = 7,
  height = 6,
  dpi = 300
)

cat("UMAP figure saved to figures/\n")
