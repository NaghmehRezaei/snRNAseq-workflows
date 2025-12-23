# 06A_ORA_pathway_enrichment.R
# ------------------------------------------------------------
# Over-representation analysis (ORA) for snRNA-seq
#
# Input:
#   - Cell-type–specific DE results (CSV)
# Output:
#   - ORA tables and dotplots
#
# Author: Naghmeh Rezaei
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(dplyr)
  library(ggplot2)
})

# ---------------------------
# 1. Load DE results
# ---------------------------

de <- read.csv(
  "results/hepatocyte_DE_results.csv",
  stringsAsFactors = FALSE
)

# ---------------------------
# 2. Filter significant genes
# ---------------------------

sig_genes <- de %>%
  filter(
    padj < 0.05,
    abs(log2FoldChange) > 0.5
  )

gene_symbols <- unique(sig_genes$gene)

# ---------------------------
# 3. Convert to Entrez IDs
# ---------------------------

entrez <- bitr(
  gene_symbols,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

# ---------------------------
# 4. ORA using GO
# ---------------------------

ego <- enrichGO(
  gene          = entrez$ENTREZID,
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

# ---------------------------
# 5. Save results
# ---------------------------

write.csv(
  as.data.frame(ego),
  file = "results/ORA_GO_hepatocyte.csv",
  row.names = FALSE
)

p <- dotplot(ego, showCategory = 20) +
  ggtitle("ORA – Hepatocyte (GO Biological Process)")

ggsave(
  filename = "figures/ORA_GO_hepatocyte.png",
  plot = p,
  width = 8,
  height = 6,
  dpi = 300
)

cat("ORA analysis completed\n")
