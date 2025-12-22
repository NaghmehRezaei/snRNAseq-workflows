# 01_pseudobulk_DE_snRNAseq.R
# ------------------------------------------------------------
# Pseudobulk differential expression analysis for snRNA-seq
#
# This script aggregates expression at the sample × cell type
# level and performs DE analysis using bulk RNA-seq methods.
#
# Pseudobulk is preferred over single-cell DE for:
# - robust statistics
# - proper replication
# - interpretability
#
# Author: Naghmeh Rezaei
# Repository: snRNAseq-workflows
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(edgeR)
})

# ---------------------------
# 1. Load annotated object
# ---------------------------

input_path <- "results/seurat_annotated.rds"
seu <- readRDS(input_path)

# ---------------------------
# 2. Metadata requirements
# ---------------------------

required_meta <- c("sample", "condition", "cell_type")
missing_meta <- setdiff(required_meta, colnames(seu@meta.data))
if (length(missing_meta) > 0) {
  stop("Missing required metadata: ", paste(missing_meta, collapse = ", "))
}

# ---------------------------
# 3. Create pseudobulk counts
# ---------------------------

pb_counts <- AggregateExpression(
  seu,
  group.by = c("sample", "cell_type"),
  assays = "RNA",
  slot = "counts",
  return.seurat = FALSE
)$RNA

# Create sample metadata
pb_meta <- expand.grid(
  sample = unique(seu$sample),
  cell_type = unique(seu$cell_type)
) %>%
  filter(paste(sample, cell_type, sep = "_") %in% colnames(pb_counts))

pb_meta$condition <- seu@meta.data %>%
  distinct(sample, condition) %>%
  right_join(pb_meta, by = "sample") %>%
  pull(condition)

# ---------------------------
# 4. Differential expression (edgeR)
# ---------------------------

dge <- DGEList(counts = pb_counts)
dge <- calcNormFactors(dge)

design <- model.matrix(~ condition)
dge <- estimateDisp(dge, design)
fit <- glmQLFit(dge, design)

qlf <- glmQLFTest(fit, coef = 2)

de_results <- topTags(qlf, n = Inf)$table
de_results$gene <- rownames(de_results)

# ---------------------------
# 5. Save DE results
# ---------------------------

write.csv(
  de_results,
  file = "results/pseudobulk_DE_results.csv",
  row.names = FALSE
)

cat("Pseudobulk DE results saved\n")
