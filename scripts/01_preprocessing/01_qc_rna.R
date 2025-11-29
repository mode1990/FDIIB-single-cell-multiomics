{\rtf1\ansi\ansicpg1252\cocoartf2761
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\paperw11900\paperh16840\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs24 \cf0 ################################################################################\
# STEP 1: scRNA-seq Quality Control\
# Description: Load, QC, and filter single-cell RNA-seq data\
# Input: Raw 10X output\
# Output: QC-filtered RDS objects\
################################################################################\
\
# Load libraries\
library(Seurat)\
library(dplyr)\
library(ggplot2)\
\
# Set paths (CHANGE THESE TO YOUR PATHS)\
sample_path <- "/path/to/sample1/filtered_feature_bc_matrix/"\
output_path <- "results/rds_objects/sample1_rna_qc.rds"\
\
# Load data\
rna_counts <- Read10X(sample_path)$`Gene Expression`\
pbmc_rna <- CreateSeuratObject(counts = rna_counts)\
\
# Calculate QC metrics\
pbmc_rna[["percent.mt"]] <- PercentageFeatureSet(pbmc_rna, pattern = "^MT-")\
\
# Pre-QC plot\
print("Creating pre-QC plots...")\
p1 <- VlnPlot(pbmc_rna, \
              features = c("nCount_RNA", "percent.mt"), \
              ncol = 2, \
              log = TRUE, \
              pt.size = 0) + NoLegend()\
ggsave("results/figures/01_rna_pre_qc.png", p1, width = 10, height = 6)\
\
# Filter cells\
print("Filtering cells...")\
pbmc_rna <- subset(pbmc_rna, \
                   subset = nCount_RNA < 25000 & \
                            nCount_RNA > 1000 & \
                            percent.mt < 20)\
\
# Post-QC plot\
print("Creating post-QC plots...")\
p2 <- VlnPlot(pbmc_rna, \
              features = c("nCount_RNA", "percent.mt"), \
              ncol = 2, \
              log = TRUE, \
              pt.size = 0) + NoLegend()\
ggsave("results/figures/01_rna_post_qc.png", p2, width = 10, height = 6)\
\
# Save\
print("Saving QC'd object...")\
saveRDS(pbmc_rna, file = output_path)\
\
print("\uc0\u10003  RNA QC complete!")\
print(paste("Cells retained:", ncol(pbmc_rna)))}