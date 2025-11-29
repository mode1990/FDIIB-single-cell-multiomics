{\rtf1\ansi\ansicpg1252\cocoartf2761
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\paperw11900\paperh16840\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs24 \cf0 ################################################################################\
# STEP 3: RNA Integration with Harmony\
# Description: Integrate samples using Harmony batch correction\
# Input: Merged RNA object\
# Output: Integrated RNA object with UMAP\
################################################################################\
\
library(Seurat)\
library(harmony)\
library(dplyr)\
\
# Load data\
print("Loading merged RNA data...")\
rna <- readRDS("results/rds_objects/merged_rna_qc_doublet_removed.rds")\
\
# Filter to specific mutations (optional)\
rna_filtered <- subset(rna, subset = Mutation %in% c("GBA1", "HC", "LRRK2", "SNCA"))\
\
# Split by sample\
print("Splitting by sample...")\
rna_list <- SplitObject(rna_filtered, split.by = "Sample")\
\
# SCTransform normalization\
print("Running SCTransform on each sample...")\
options(future.globals.maxSize = Inf)\
for (i in seq_along(rna_list)) \{\
  print(paste("SCTransform:", names(rna_list)[i]))\
  rna_list[[i]] <- SCTransform(rna_list[[i]], \
                                vars.to.regress = "percent.mt", \
                                verbose = FALSE)\
\}\
\
# Select integration features\
print("Selecting integration features...")\
features <- SelectIntegrationFeatures(object.list = rna_list, nfeatures = 3000)\
\
# Prepare for integration\
print("Preparing for integration...")\
rna_list <- PrepSCTIntegration(object.list = rna_list, \
                               anchor.features = features)\
\
# Merge\
print("Merging objects...")\
rna_merged <- merge(rna_list[[1]], y = rna_list[-1], merge.data = TRUE)\
DefaultAssay(rna_merged) <- "SCT"\
\
# Run PCA\
print("Running PCA...")\
rna_merged <- RunPCA(rna_merged, npcs = 50, verbose = FALSE, features = features)\
\
# Harmony integration\
print("Running Harmony integration...")\
rna_merged <- IntegrateLayers(\
  object = rna_merged, \
  method = HarmonyIntegration,\
  orig.reduction = "pca", \
  new.reduction = "harmony",\
  group.by.vars = "Sample",\
  verbose = FALSE\
)\
\
# Clustering and UMAP\
print("Clustering...")\
rna_merged <- FindNeighbors(rna_merged, reduction = "harmony", dims = 1:20)\
rna_merged <- FindClusters(rna_merged, resolution = 0.2)\
\
print("Running UMAP...")\
rna_merged <- RunUMAP(rna_merged, reduction = "harmony", dims = 1:20)\
\
# Save\
print("Saving integrated object...")\
saveRDS(rna_merged, "results/rds_objects/rna_integrated_harmony.rds")\
\
# Create UMAP plot\
print("Creating visualization...")\
p <- DimPlot(rna_merged, group.by = "seurat_clusters", label = TRUE)\
ggsave("results/figures/03_rna_umap_clusters.png", p, width = 10, height = 8)\
\
print("\uc0\u10003  Integration complete!")}