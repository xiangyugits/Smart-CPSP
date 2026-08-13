# Smart-CPSP

Single-cell RNA-seq analysis pipeline for the Smart-CPSP mouse study.

## Overview

This repository contains the R analysis code used to process and analyze Smart-seq3 single-cell RNA-seq data from the Smart-CPSP mouse model.

## Experimental methods

- **Cell extraction**: Single cells were sorted into 384-well plates containing lysis buffer. Total RNA was extracted from individual cells following the Smart-seq3 protocol.
- **Library construction**: cDNA was synthesized and amplified from single-cell lysates using the Smart-seq3 protocol. Libraries were constructed using tagmentation with Tn5 transposase. The resulting libraries were purified, pooled, and sequenced on the MGISEQ-2000RS platform with paired-end 150 bp reads.

## Data processing

- **Read alignment and quantification**: The mouse reference genome GRCm38.p6 with GENCODE vM25 annotation was supplemented with the WPRE sequence. A kallisto transcriptome index was built from the combined reference. Raw FASTQ files were processed using kb-python with kallisto as the pseudoalignment engine. Reads were pseudoaligned to the custom index, and gene-level count matrices were generated from the pseudoalignment output.
- **Single-cell downstream analysis**: Gene expression count matrices were analyzed using the Seurat standard workflow. Cells were filtered, normalized, and scaled. Dimensionality reduction was performed using PCA followed by UMAP visualization. Cell clustering and annotation were performed to identify cell types. Processed data files include the annotated expression matrix and cell metadata.

## Repository structure

| File | Description |
|------|-------------|
| `analysis_code.md` | Consolidated R analysis pipeline, organized into 9 sequential sections (one code block per analysis step). |

## Sections in `analysis_code.md`

1. SingleR Annotation (reference-based)
2. Cell Clustering (SCTransform + PCA + UMAP + resolution)
3. Main Cell Type Annotation
4. Neuron Sub-clustering
5. Neuron DEA (MAST)
6. Neuron DEA (Wilcoxon)
7. MG Clustering
8. MG CPSP DEG
9. CellChat

## Requirements

- R (>= 4.2)
- Key R packages: `Seurat`, `SingleR`, `celldex`, `tidyverse`, `MAST`, `clusterProfiler`, `org.Mm.eg.db`, `CellChat`, `patchwork`

## Citation

If you use this code, please cite the Smart-CPSP study.

## License

For research use. Contact the authors for licensing and data access.
