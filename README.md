# bulk-rnaseq-shiny-app

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18905852.svg)](https://doi.org/10.5281/zenodo.18905852)

A suite to analyze bulk RNA-seq data, including differential gene expression, overrepresentation analysis, gene set enrichment analysis, and custom plots.

---

## Features

- **Multi-batch support** — upload up to 3 batches of samples with optional limma batch correction
- **Differential gene expression** — powered by edgeR (quasi-likelihood F-test), with volcano, MA, and heatmap plots
- **Overrepresentation analysis (ORA)** — for up- and downregulated genes using clusterProfiler
- **Gene set enrichment analysis (GSEA)** — using fgsea with MSigDB collections (Hallmark, GO, KEGG, Reactome)
- **Pathway Explorer** — compare enrichment across conditions, view enrichment traces, and gene-level heatmaps
- **Custom gene plots** — heatmap, dotplot, or violin plots for genes of interest or top N DE genes
- **Quality control** — PCA, library size, expression distribution, and mitochondrial fraction plots
- **Gene filtering** — filter by biotype (protein-coding, lncRNA, miRNA), remove mitochondrial and Rik genes
- **Multi-species support** — mouse, human, rat, zebrafish, fly, and worm

---

## Supported Species

| Species | Genome |
|---|---|
| *Mus musculus* | GRCm38/mm10 |
| *Homo sapiens* | GRCh38/hg38 |
| *Rattus norvegicus* | Rnor6 |
| *Danio rerio* | GRCz11 |
| *Drosophila melanogaster* | BDGP6 |
| *Caenorhabditis elegans* | WBcel235 |

---

## Input Format

Each sample should be a tab-separated (`.tsv` or `.txt`) count file with:
- Gene IDs (Ensembl or gene symbols) as row names
- A single column of raw counts

Example:
```
ENSMUSG00000000001    145
ENSMUSG00000000003    0
ENSMUSG00000000028    892
```

---

## Installation

### Requirements

- R ≥ 4.1
- The following R packages:

```r
# CRAN
install.packages(c("shiny", "ggplot2", "dplyr", "tidyr", "ggrepel",
                   "writexl", "openxlsx", "DT", "pheatmap", "gridExtra"))

# Bioconductor
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c(
  "edgeR", "limma", "biomaRt", "clusterProfiler", "enrichplot",
  "fgsea", "msigdbr", "BiocParallel", "AnnotationDbi",
  "org.Mm.eg.db", "org.Hs.eg.db", "org.Rn.eg.db",
  "org.Dr.eg.db", "org.Dm.eg.db", "org.Ce.eg.db"
))
```

### Running the App

```r
# Clone the repo
git clone https://github.com/eabikh/bulk-rnaseq-shiny-app.git
cd bulk-rnaseq-shiny-app

# Launch in R
shiny::runApp("bbulk.r")
```

---

## Usage

1. **Upload** your count files (one `.tsv` per sample) under Batch 1
2. **Enter group labels** matching the file order, comma-separated (e.g. `Control,Control,KO,KO`)
3. Optionally **add more batches** and enable limma batch correction
4. Select your **reference genome** and configure filtering/threshold options
5. Click **Run Analysis**
6. Explore results across the tabs: QC → Differential Expression → ORA → GSEA → Pathway Explorer → Custom Genes
7. **Download** results as `.xlsx`, `.csv`, or `.png` files

---

## Batch Correction

When multiple batches are provided, the app uses `limma::removeBatchEffect` on logCPM values for visualization (PCA, heatmaps). The differential expression model includes batch as a covariate in the design matrix.

---

## Output

| Download | Format | Contents |
|---|---|---|
| Raw Counts | `.csv` | Raw count matrix |
| DE Results | `.xlsx` | All comparisons, one sheet each |
| Sig Genes | `.csv` | Significant genes for selected comparison |
| GSEA Results | `.xlsx` | fgsea results with summary sheet |
| All plots | `.png` | Publication-ready figures (300 dpi) |

---

## Citation

If you use this app in your research, please cite the underlying tools:

- **edgeR**: Robinson et al., *Bioinformatics* (2010)
- **limma**: Ritchie et al., *Nucleic Acids Research* (2015)
- **clusterProfiler**: Wu et al., *The Innovation* (2021)
- **fgsea**: Korotkevich et al., *bioRxiv* (2021)
- **MSigDB**: Liberzon et al., *Cell Systems* (2015)
