# 03_celltype_specific_DE_snRNAseq.R
# ------------------------------------------------------------
# Cell-type–specific pseudobulk differential expression
#
# This script performs DE analysis separately for each
# annotated cell type using pseudobulk aggregation and edgeR.
#
# This approach:
# - respects biological replication
# - avoids single-cell DE pitfalls
# - is suitable for publication-quality analyses
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
# 2. Metadata checks
# ---------------------------

required_meta <- c("sample", "condition", "cell_type")
missing_meta <- setdiff(required_meta, colnames(seu@meta.data))
if (length(missing_meta) > 0) {
  stop("Missing required metadata: ", paste(missing_meta, collapse = ", "))
}

cell_types <- sort(unique(seu$cell_type))

# Create output directory
dir.create("results/celltype_DE", showWarnings = FALSE)

# ---------------------------
# 3. Loop over cell types
# ---------------------------

for (ct in cell_types) {

  message("Processing cell type: ", ct)

  # Subset to one cell type
  seu_ct <- subset(seu, subset = cell_type == ct)

  # Require at least 2 samples per condition
  ct_meta <- seu_ct@meta.data %>%
    distinct(sample, condition)

  if (length(unique(ct_meta$condition)) < 2 ||
      min(table(ct_meta$condition)) < 2) {
    message("  Skipping ", ct, " (insufficient replication)")
    next
  }

  # ---------------------------
  # 4. Pseudobulk aggregation
  # ---------------------------

  pb_counts <- AggregateExpression(
    seu_ct,
    group.by = "sample",
    assays = "RNA",
    slot = "counts",
    return.seurat = FALSE
  )$RNA

  pb_meta <- ct_meta %>%
    distinct(sample, condition) %>%
    arrange(match(sample, colnames(pb_counts)))

  # ---------------------------
  # 5. edgeR DE
  # ---------------------------

  dge <- DGEList(counts = pb_counts)
  dge <- calcNormFactors(dge)

  design <- model.matrix(~ condition, data = pb_meta)
  dge <- estimateDisp(dge, design)
  fit <- glmQLFit(dge, design)

  qlf <- glmQLFTest(fit, coef = 2)

  de_res <- topTags(qlf, n = Inf)$table
  de_res$gene <- rownames(de_res)
  de_res$cell_type <- ct

  # ---------------------------
  # 6. Save results
  # ---------------------------

  out_file <- paste0(
    "results/celltype_DE/DE_",
    gsub(" ", "_", ct),
    ".csv"
  )

  write.csv(de_res, out_file, row.names = FALSE)
}

cat("Cell-type–specific pseudobulk DE completed\n")
