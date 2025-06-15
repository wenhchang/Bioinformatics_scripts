# Bioinformatics Scripts

Project-Root/
├── RNA-seq/
│   ├── QC_STAR_featureCounts/
│   │   ├── FastQC.sh
│   │   ├── STAR_FC.sh
│   │   └── STAR_index.sh
│   ├── 0.1_biomart.R
│   ├── 0.2_combined_count_table.R
│   ├── 0.3_TPM_conversion.R
│   ├── 1.1_DESeq2_LFCS.R
│   ├── 1.2_Limma_blocking_term.R
│   ├── 1.5_Volcano_signatures.R
│   ├── 2.0_GSEA_ORA.R
│   ├── 2.5_visualize_Hallmark.R
│   ├── 3.0_resistance_analysis.R
│   ├── 4.0_Euler_diagrams.R
│   ├── 5.0_TCGA_analysis.R
│   └── 6.0_DepMap_analysis.R
│
├── RPPA/
│   ├── 1.0_preprocessing_QC.R
│   ├── 2.0_Differential_analysis_Limma.R
│   └── 3.0_Heatmap_RPPA.R
│
├── README.md
└── LICENSE

# RNA-seq Analysis Pipeline

This repository contains a RNA-seq data analysis pipeline designed for bulk mRNA transcriptomic profiling. It integrates shell-based preprocessing (on HPC clusters) with R-based downstream analysis to support a wide range of experimental designs and biological questions.

This repository is designed as templates for RNA-seq analysis and does NOT contain the reference data tables used in each script. The reference data table and an example R markdown summary file are available upon reasonable request.

### Workflow Overview

**The following steps are carried out in an HPC environment:**

1. **Quality control (QC)** of raw FASTQ files using `FastQC` to assess sequence quality, GC content, adapter contamination, etc.
2. **Read alignment** to the reference genome using `STAR`. This includes a script to build the STAR genome index and to run alignment in two-pass mode.
3. **Read quantification** using `featureCounts` to generate a gene-level raw count matrix for all samples.

**Following alignment, downstream analysis is performed in R:**

1. **Differential expression analysis** using either `DESeq2` or `limma`, depending on the experimental design (e.g., with or without batch/blocking terms).
2. **Pathway enrichment analysis**, including Gene Set Enrichment Analysis (GSEA) and Overrepresentation Analysis (ORA), to interpret gene expression changes in the context of known biological pathways.
3. **Functional profiling and integration** with external datasets such as TCGA (for clinical correlation) and DepMap (for gene dependency analysis), as well as gene overlap visualization using Euler plots.
4. **Resistance-related analysis**, including identification of pathways or genes associated with resistance phenotypes.


### NOTE ###
####### An example R Markdown summary report and reference files are available upon request. 

## Directory Overview

This folder includes the following components:

### Shell Scripts
- `FastQC.sh`  
  Perform quality control on raw FASTQ files using **FastQC**.

- `STAR_index.sh`  
  Generate genome index files for STAR alignment using reference genome and GTF annotation.

- `STAR_FC.sh`  
  Full pipeline for:
  - Adapter trimming with **cutadapt**
  - Read alignment with **STAR**
  - Gene-level quantification with **featureCounts**

### R Scripts

#### [0.X] Preprocessing and Normalization
- `0.1_biomart.R`: Retrieve gene annotations using `biomaRt`.
- `0.2_combine_count_table.R`: Merge multiple count matrices into a single expression table.
- `0.3_TPM_conversion.R`: Convert raw counts to TPM values.

#### [1.] Differential Expression Analysis
- `1.1_DEseq2_LFCS.R`: Perform DE analysis using DESeq2.
- `1.2_Limma_blocking_term.R`: Limma-based DE analysis with blocking terms.
- `1.5_Volcanoplots_with_other_signatures.R`: Visualization of differential genes with external gene sets.

#### [2.] Pathway and Signature Analysis
- `2.0_GSEA_ORA.R`: Gene set enrichment analysis (GSEA) and overrepresentation analysis.
- `2.5_visualization_GSEA_Hallmark.R`: Custom visualization of hallmark GSEA results.

#### [3.–6.] Functional, TCGA, and DepMap Analysis
- `3.0_resistance analysis.R`: Analysis related to drug resistance (e.g., gene scores, clustering).
- `4.0_Euler_diagrams.R`: Generate Euler plots for gene set overlaps.
- `5.0_TCGA_analysis.R`: Use TCGA datasets for clinical or expression-based correlation.
- `6.0_DepMap_analysis.R`: Analyze gene dependency data from the DepMap portal.

---

## Requirements

- R (≥ 4.1.0)
- STAR aligner
- Subread (for featureCounts)
- Cutadapt
- FastQC

---

# RPPA (Reverse Phase Protein Array) analysis pipeline

This repository contains a pipeline for processing and analysing Reverse Phase Protein Array (RPPA) data.

### Workflow overview
All the following analysese are done in R envirnoment.
1. **Normalization** of the raw RPPA intensity file is performed using median centering.
2. **Quality control** of the raw RPPA data is performed by PCA analysis and correlation matrix. This ensures the reproducibility between replicates within the same conditons.
3. **Differential Expression (DE)** is performed using Limma.
4. **Visualization** of the DE result by heatmap.

### R scripts
- `1.0_preprocessing_QC.R`: Median centering normalization, PCA analysis, and correlation coeficients matrix analysis to normalize the raw intensities and ensure the quality of the data.
- `2.0_Differential_analysis_Limma.R`: Conduct Limma-based differential expression analysis.
- `3.0_Heatmap_RPPA.R`: Visualize the DE result for RPPA.


---

## License

This project is licensed under the MIT License. See the [LICENSE](../LICENSE) file for details.

---

## Author

Wen-Hsuan Chang  
wenhsuanc@icloud.com
