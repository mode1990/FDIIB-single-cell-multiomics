# FOUNDIN-PD: Multiomics Analysis Pipeline

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R Version](https://img.shields.io/badge/R-%3E%3D4.3.0-blue)](https://www.r-project.org/)

> Integrative analysis of scRNA-seq and scATAC-seq data from iPSC-derived dopaminergic neurons in Parkinson's disease

## 🔬 Overview

This pipeline analyzes **10x Genomics Multiome data** to compare cells with Parkinson's disease mutations (SNCA, LRRK2, GBA1) to healthy controls.

### Key Features

- ✅ **Quality Control**: Automated filtering for RNA and ATAC data
- ✅ **Doublet Removal**: scDblFinder integration
- ✅ **Batch Correction**: Harmony-based integration
- ✅ **Cell Type Annotation**: 7 neuronal populations (ENP, LNP, EPN, NEP, PFPP, IDN, DAN)
- ✅ **Differential Analysis**: Gene expression and chromatin accessibility
- ✅ **Motif Enrichment**: JASPAR-based TF binding site analysis
- ✅ **Multi-omics Integration**: MOFA2 factor analysis

## 📊 Data

- **Assays**: Gene expression (scRNA-seq) + Chromatin accessibility (scATAC-seq)
- **Platform**: 10x Genomics Multiome
- **Cell Types**: 7 neuronal populations
- **Conditions**: Healthy control (HC), SNCA, LRRK2, GBA1 mutations

## 🚀 Quick Start

### Prerequisites

- R ≥ 4.3.0
- Python ≥ 3.8 (for MACS2)
- 32GB RAM minimum

### Installation
```bash
# Clone repository
git clone https://github.com/mode1990/FOUNDIN-PD-analysis.git
cd FOUNDIN-PD-analysis

# Install R packages (in R console)
source("scripts/install_packages.R")

# Set up Python environment
conda env create -f environment.yml
conda activate foundin-pd
```

### Running the Pipeline

Execute scripts in numerical order:
```bash
# 1. RNA preprocessing
Rscript scripts/01_preprocessing/01_qc_rna.R
Rscript scripts/01_preprocessing/02_merge_rna_samples.R

# 2. Integration
Rscript scripts/02_integration/01_integrate_rna.R
Rscript scripts/02_integration/02_integrate_atac.R

# 3. Analysis
Rscript scripts/03_analysis/01_annotate_cell_types.R
Rscript scripts/03_analysis/02_differential_expression.R
Rscript scripts/03_analysis/03_differential_accessibility.R

# 4. Visualization
Rscript scripts/04_visualization/01_make_figures.R
```

## 📂 Repository Structure
## 📖 Documentation

- **[Installation Guide](docs/installation.md)** - Detailed setup instructions
- **[Tutorial](docs/tutorial.md)** - Step-by-step walkthrough
- **[Troubleshooting](docs/troubleshooting.md)** - Common issues and solutions

## 📦 Dependencies

**Core R Packages:**
- Seurat (5.0.3) - Single-cell RNA-seq analysis
- Signac (1.13.0) - Single-cell ATAC-seq analysis
- Harmony - Batch correction
- MOFA2 - Multi-omics integration
- scDblFinder - Doublet detection

**Python:**
- MACS2 (2.2.7) - Peak calling

See [requirements.txt](requirements.txt) for complete list.



## 📄 License

This project is licensed under the MIT License - see [LICENSE](LICENSE) file.

---

**Author**: Mo Dehestani  
**Last Updated**: Nov 2025
