# snRNAseq-workflows

This repository contains reproducible **single-nucleus RNA-seq (snRNA-seq) analysis workflows**
used in real research projects, covering quality control, normalization, integration,
clustering, annotation, and downstream analyses.

The goal of this repository is to demonstrate **practical, research-grade snRNA-seq pipelines**
with an emphasis on:
- biological interpretability
- reproducibility
- clean, modular code

---

## 📁 Repository structure
snRNAseq-workflows/
├── data/
│   └── example_metadata/        # Example metadata tables (no raw sequencing data)
├── scripts/
│   ├── 01_qc/                   # QC and filtering
│   ├── 02_normalization/        # Normalization and scaling
│   ├── 03_integration/          # Dataset integration (Harmony, Seurat)
│   ├── 04_clustering/           # Dimensionality reduction and clustering
│   ├── 05_annotation/           # Cell type annotation
│   └── 06_downstream/           # DE, pathway analysis, ML
├── envs/                        # Conda / R environment specifications
├── figures/                     # Generated figures (PNG/PDF)
├── results/                     # Result tables and summaries
└── README.md


---

## 🧬 Typical workflow

1. **Quality control**
   - Filtering nuclei by counts, features, and mitochondrial content
   - Doublet detection (Scrublet / SCDS)

2. **Normalization**
   - Log-normalization / SCTransform
   - Variable feature selection

3. **Integration**
   - Batch correction using Harmony or Seurat integration
   - Support for multi-sample and multi-condition designs

4. **Clustering & visualization**
   - PCA, UMAP
   - Cluster resolution tuning

5. **Annotation**
   - Marker-based cell type annotation
   - Reference-guided approaches (when applicable)

6. **Downstream analysis**
   - Differential expression
   - Pathway enrichment
   - Machine-learning–based representations (e.g. scVI)

---

## ⚠️ Notes

- No raw sequencing data are included in this repository.
- Scripts are adapted from real analyses and cleaned for public sharing.
- Paths and sample identifiers are generalized.

---

## 👩‍🔬 Author

**Naghmeh Rezaei**  
Computational Biology · Single-cell & Spatial Omics · Machine Learning
