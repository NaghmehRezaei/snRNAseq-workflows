# 05_cell_composition_analysis_snRNAseq.R
# ------------------------------------------------------------
# Cell type composition analysis for snRNA-seq
#
# This script quantifies and visualizes differences in
# cell-type proportions across samples and conditions.
#
# Composition analysis is complementary to DE and helps
# distinguish abundance changes from transcriptional changes.
#
# Author: Naghmeh Rezaei
# Repository: snRNAseq-workflows
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

# ---------------------------
# 1. Load annotated object
# ---------------------------

seu <- readRDS("results/seurat_annotated.rds")

# ---------------------------
# 2. Metadata checks
# ---------------------------

required_meta <- c("sample", "condition", "cell_type")
missing_meta <- setdiff(required_meta, colnames(seu@meta.data))
if (length(missing_meta) > 0) {
  stop("Missing required metadata: ", paste(missing_meta, collapse = ", "))
}

# ---------------------------
# 3. Compute cell-type proportions
# ---------------------------

comp_table <- seu@meta.data %>%
  count(sample, condition, cell_type) %>%
  group_by(sample) %>%
  mutate(
    proportion = n / sum(n)
  ) %>%
  ungroup()

# Save composition table
write.csv(
  comp_table,
  file = "results/cell_composition_table.csv",
  row.names = FALSE
)

# ---------------------------
# 4. Visualization
# ---------------------------

p <- ggplot(
  comp_table,
  aes(
    x = sample,
    y = proportion,
    fill = cell_type
  )
) +
  geom_bar(stat = "identity") +
  facet_wrap(~ condition, scales = "free_x") +
  theme_bw() +
  ylab("Cell-type proportion") +
  xlab("Sample") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

# ---------------------------
# 5. Save figure
# ---------------------------

ggsave(
  filename = "figures/cell_type_composition_barplot.pdf",
  plot = p,
  width = 10,
  height = 6
)

cat("Cell composition analysis completed\n")
