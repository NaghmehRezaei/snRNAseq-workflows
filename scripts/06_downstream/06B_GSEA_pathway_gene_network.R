# 06B_GSEA_pathway_gene_network.R
# ------------------------------------------------------------
# GSEA-based pathway–gene network construction
#
# Input:
#   - Ranked gene list (log2FC)
# Output:
#   - Pathway–gene network figure
#
# Author: Naghmeh Rezaei
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(dplyr)
  library(igraph)
  library(ggraph)
  library(ggplot2)
})

# ---------------------------
# 1. Load ranked gene list
# ---------------------------

de <- read.csv(
  "results/hepatocyte_DE_results.csv",
  stringsAsFactors = FALSE
)

ranked_genes <- de$log2FoldChange
names(ranked_genes) <- de$gene
ranked_genes <- sort(ranked_genes, decreasing = TRUE)

# ---------------------------
# 2. Run GSEA (Hallmark)
# ---------------------------

hallmark <- msigdbr::msigdbr(
  species = "Mus musculus",
  category = "H"
)

gsea <- GSEA(
  ranked_genes,
  TERM2GENE = hallmark[, c("gs_name", "gene_symbol")],
  pvalueCutoff = 0.05
)

# ---------------------------
# 3. Extract leading-edge genes
# ---------------------------

le <- gsea@result %>%
  filter(NES > 0, p.adjust < 0.05) %>%
  select(ID, leading_edge)

edges <- le %>%
  tidyr::separate_rows(leading_edge, sep = ",") %>%
  rename(pathway = ID, gene = leading_edge)

# ---------------------------
# 4. Build network
# ---------------------------

g <- graph_from_data_frame(edges, directed = FALSE)

# ---------------------------
# 5. Plot network
# ---------------------------

p <- ggraph(g, layout = "fr") +
  geom_edge_link(alpha = 0.3) +
  geom_node_point(size = 3, color = "red") +
  geom_node_text(aes(label = name), repel = TRUE, size = 2) +
  ggtitle("GSEA Pathway–Gene Network (Hepatocyte)")

ggsave(
  "figures/GSEA_pathway_gene_network_hepatocyte.png",
  p,
  width = 10,
  height = 8,
  dpi = 300
)

cat("GSEA network completed\n")
