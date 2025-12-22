# 06_scVI_latent_space_snRNAseq.py
# ------------------------------------------------------------
# scVI latent space modeling for snRNA-seq
#
# This script trains an scVI model to learn a low-dimensional,
# denoised latent representation of snRNA-seq data.
#
# scVI is used here for:
# - robustness checks
# - latent-space exploration
# - downstream machine learning analyses
#
# It is NOT a replacement for Seurat/Harmony integration.
#
# Author: Naghmeh Rezaei
# Repository: snRNAseq-workflows
# ------------------------------------------------------------

import scvi
import scanpy as sc
import anndata as ad

# ---------------------------
# 1. Load data
# ---------------------------

# Input: counts matrix exported from Seurat
# (e.g. via SeuratDisk::SaveH5Seurat + Convert)
adata = sc.read_h5ad("results/seurat_for_scvi.h5ad")

# ---------------------------
# 2. Setup scVI model
# ---------------------------

# Register AnnData with scVI
scvi.model.SCVI.setup_anndata(
    adata,
    batch_key="sample"
)

# ---------------------------
# 3. Train model
# ---------------------------

model = scvi.model.SCVI(
    adata,
    n_latent=30
)

model.train(
    max_epochs=100,
    early_stopping=True
)

# ---------------------------
# 4. Extract latent space
# ---------------------------

adata.obsm["X_scVI"] = model.get_latent_representation()

# ---------------------------
# 5. Save outputs
# ---------------------------

adata.write("results/seurat_scvi_latent.h5ad")

print("scVI latent space modeling completed")
