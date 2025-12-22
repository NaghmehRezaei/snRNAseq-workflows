# 02_pathway_GSEA_snRNAseq.R
# ------------------------------------------------------------
# Gene Set Enrichment Analysis (GSEA) for snRNA-seq pseudobulk DE
#
# This script performs pathway enrichment using ranked
# pseudobulk differential expression results.
#
# GSEA is preferred over over-representation analysis because:
# - it uses all genes
# - it avoids arbitrary cutoffs
# - it is robust for transcriptomic data
#
# Author: Naghmeh Rezaei
# Repository: snRNAseq-workflows
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(clusterProfiler)
  library(msigdbr)
  library(enrichplot)
})

# ---------------------------
# 1. Load pseudobulk DE results
# ---------------------------

de_path <- "results/pseudobulk_DE_results.csv"
de <- read.csv(de_path)

# ---------------------------
# 2. Prepare ranked gene list
# ---------------------------

# Use logFC as ranking metric
# Gene symbols assumed
gene_ranks <- de %>%
  filter(!is.na(logFC)) %>%
  arrange(desc(logFC)) %>%
  distinct(gene, .keep_all = TRUE) %>%
  select(gene, logFC)

ranked_genes <- gene_ranks$logFC
names(ranked_genes) <- gene_ranks$gene

# ---------------------------
# 3. Load pathway gene sets
# ---------------------------

# Hallmark gene sets (recommended starting point)
msig_hallmark <- msigdbr(
  species = "Mus musculus",
  category = "H"
)

hallmark_list <- split(
  msig_hallmark$gene_symbol,
  msig_hallmark$gs_name
)

# ---------------------------
# 4. Run GSEA
# ---------------------------

gsea_res <- GSEA(
  geneList = ranked_genes,
  TERM2GENE = hallmark_list,
  pvalueCutoff = 0.05,
  verbose = FALSE
)

# ---------------------------
# 5. Save results
# ---------------------------

gsea_table <- as.data.frame(gsea_res)

write.csv(
  gsea_table,
  file = "results/GSEA_hallmark_results.csv",
  row.names = FALSE
)

# ---------------------------
# 6. Basic visualization
# ---------------------------

pdf("figures/GSEA_hallmark_dotplot.pdf", width = 8, height = 6)
dotplot(gsea_res, showCategory = 20)
dev.off()

cat("GSEA analysis completed and results saved\n")
