# Rice RNA-seq Transcriptome Analysis Pipeline

This repository stores reproducible R scripts for RNA-seq analysis of rice cold-stress experiments, developed for screening cold-responsive splicing factors and functional characterization of candidate genes (OsSR32/OsSR33).

> 
> Core highlight: Resolved the lack of native MSU7 OrgDb annotation package for rice, implemented two independent ID-mapping strategies for gene annotation.

## Project Structure

```
├── annotation_enrichment/       # Gene ID mapping & functional enrichment
│   ├── AnnotationHub-mapping.R   # Dual ID annotation strategy: biomaRt (primary) + httr2 REST API (fallback)
│   ├── AnnotationHub.sh          # Shell helper for HPC
│   ├── GO_batch.R                # Batch GO enrichment, table export & dotplot visualization
│   └── KEGG_batch.R              # Batch KEGG enrichment using offline gson database
├── diff_expression/              # Differential expression analysis
│   ├── DESeq2.R                  # DESeq2 model fitting, contrast setup, DEG extraction
│   └── count_table_clean.R       # Clean featureCounts raw count matrix, sample name preprocessing
├── visualization/                # Plotting scripts
│   ├── pheatmap.R               # Heatmap for DEG expression (row-wise Z-score normalization)
│   └── volcano_plot.R           # Volcano plot for differentially expressed genes
├── DESeq2.R                      # Main workflow entry
└── README.md
```

## Workflow Overview

1. **Count matrix preprocessing**: Clean featureCounts output, remove redundant columns, standardize sample names.
2. **Differential expression analysis**: Construct interaction model `~ genotype + time + genotype:time` with DESeq2, extract DEGs (padj<0.05, |log2FC| ≥1).
3. **Gene ID mapping (key innovation)**:
   - Challenge: Rice MSU7 gene IDs do not have an official OrgDb annotation library.
   - Strategy 1: Use `biomaRt` to convert RAP IDs to ENTREZ IDs (preferred).
   - Strategy 2: Fallback to Biomart REST API via `httr2` when biomaRt connection fails.
4. **Functional enrichment**: Batch GO / KEGG enrichment for all comparison groups, simplify GO terms, export tables and PDF dotplots.
5. **Visualization**: Volcano plots for DEGs, heatmap of normalized expression of significant genes.

## Dependencies

```
library(DESeq2)
library(tidyverse)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(clusterProfiler)
library(AnnotationHub)
library(biomaRt)
library(httr2)
library(openxlsx)
```

## Usage

1. Prepare raw count matrix from featureCounts.
2. Run `count_table_clean.R` to clean count data.
3. Run `DESeq2.R` to fit model and extract DEGs.
4. Run `AnnotationHub-mapping.R` to perform two-pathway gene ID conversion.
5. Run `GO_batch.R` and `KEGG_batch.R` for batch functional enrichment.
6. Run scripts under `visualization/` to generate publication-ready figures.

## Notes

- Relative file paths are used, compatible with local machine and HPC environment.
- Offline KEGG `gson` object is saved locally to avoid unstable online database access.
- All loops and functions are modularized for reuse on additional comparison groups.
