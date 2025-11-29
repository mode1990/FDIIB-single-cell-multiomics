{\rtf1\ansi\ansicpg1252\cocoartf2761
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\paperw11900\paperh16840\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs24 \cf0 ################################################################################\
# STEP 2: Merge RNA Samples & Remove Doublets\
# Description: Combine multiple samples and filter doublets\
# Input: Individual QC'd RDS files\
# Output: Merged RDS object\
################################################################################\
\
library(Seurat)\
library(scDblFinder)\
library(SingleCellExperiment)\
\
# Set paths\
rds_directory <- "results/rds_objects/"\
output_file <- "results/rds_objects/merged_rna_qc_doublet_removed.rds"\
\
# Load all samples\
print("Loading samples...")\
rna_samples <- list.files(path = rds_directory, \
                          pattern = "sample.*_rna_qc.rds", \
                          full.names = TRUE)\
rna_objects <- lapply(rna_samples, readRDS)\
\
# Add sample metadata\
sample_names <- paste0("Sample", seq_along(rna_objects))\
for (i in seq_along(rna_objects)) \{\
  rna_objects[[i]]@meta.data$Sample <- sample_names[i]\
\}\
\
# Remove doublets\
print("Removing doublets...")\
filtered_rna_objects <- list()\
for (i in seq_along(rna_objects)) \{\
  print(paste("Processing", sample_names[i]))\
  sce <- as.SingleCellExperiment(rna_objects[[i]])\
  sce <- scDblFinder(sce)\
  rna_objects[[i]]$scDblFinder_class <- sce$scDblFinder.class\
  filtered_rna_objects[[i]] <- subset(rna_objects[[i]], \
                                      subset = scDblFinder_class == "singlet")\
\}\
\
# Merge\
print("Merging samples...")\
merged_rna <- do.call(merge, filtered_rna_objects)\
\
# Save\
print("Saving merged object...")\
saveRDS(merged_rna, file = output_file)\
\
print("\uc0\u10003  Merging complete!")\
print(paste("Total cells:", ncol(merged_rna)))}