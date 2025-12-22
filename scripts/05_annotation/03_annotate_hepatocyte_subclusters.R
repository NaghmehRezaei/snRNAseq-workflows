# 03_annotate_hepatocyte_subclusters.R
# ------------------------------------------------------------
# Annotation of hepatocyte subclusters (snRNA-seq)
#
# Description:
# Assigns biological labels to hepatocyte subclusters
# based on known zonation and functional marker genes.
#
# Author: Naghmeh Rezaei
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

# ---------------------------
# 1. Load subclustered hepatocytes
# ---------------------------

hep <- readRDS("results/hepatocyte_subclustered.rds")

# ---------------------------
# 2. Inspect cluster markers (manual step)
# ---------------------------
# Use marker file generated previously:
# results/hepatocyte_subcluster_markers.csv

# Example zonation markers:
# Periportal: CPS1, HAL, ASS1
# Pericentral: CYP2E1, GLUL, CYP1A2
# Stress/inflammatory: JUN, FOS, ATF3

# ---------------------------
# 3. Define cluster annotations
# ---------------------------

hep$hep_subtype <- dplyr::case_when(
  hep$seurat_clusters %in% c("0", "2") ~ "Periportal-like hepatocytes",
  hep$seurat_clusters %in% c("1")      ~ "Mid-zonal hepatocytes",
  hep$seurat_clusters %in% c("3")      ~ "Pericentral-like hepatocytes",
  hep$seurat_clusters %in% c("4")      ~ "Stress-associated hepatocytes",
  TRUE ~ "Unassigned"
)

hep$hep_subtype <- factor(
  hep$hep_subtype,
  levels = c(
    "Periportal-like hepatocytes",
    "Mid-zonal hepatocytes",
    "Pericentral-like hepatocytes",
    "Stress-associated hepatocytes",
    "Unassigned"
  )
)

# ---------------------------
# 4. Save annotated object
# ---------------------------

saveRDS(
  hep,
  file = "results/hepatocyte_subcluster_annotated.rds"
)

cat("Hepatocyte subclusters annotated successfully\n")
