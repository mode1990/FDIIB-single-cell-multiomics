{\rtf1\ansi\ansicpg1252\cocoartf2761
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\paperw11900\paperh16840\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs24 \cf0 ################################################################################\
# STEP 4: Annotate Cell Types\
# Description: Annotate clusters based on marker genes\
# Input: Integrated RNA object\
# Output: Annotated RNA object\
################################################################################\
\
library(Seurat)\
library(ggplot2)\
\
# Load data\
print("Loading integrated data...")\
rna <- readRDS("results/rds_objects/rna_integrated_harmony.rds")\
\
# Define marker genes\
marker_genes <- list(\
  ENP = c("RFX4", "HES1", "SLIT2"),\
  LNP = c("DLK1", "LGALS1", "VCAN"),\
  DAN = c("TH", "ATP1A3", "ZCCHC12", "MAP2", "SYT1"),\
  IDN = c("TPH1", "SLC18A1", "SLC18A2", "SNAP25"),\
  PFPP = c("HMGB2", "TOP2A", "MKI67"),\
  NEP = c("KRT19", "KRT8", "COL17A1"),\
  EPN = c("MLF1", "STOML3", "FOXJ1")\
)\
\
# Visualize markers\
DefaultAssay(rna) <- "SCT"\
print("Creating marker plots...")\
for (cell_type in names(marker_genes)) \{\
  print(paste("Plotting:", cell_type))\
  p <- FeaturePlot(rna, \
                   features = marker_genes[[cell_type]], \
                   reduction = "umap", \
                   ncol = 2)\
  ggsave(paste0("results/figures/04_markers_", cell_type, ".png"), \
         p, width = 12, height = 8)\
\}\
\
# Manual annotation (ADJUST THESE BASED ON YOUR MARKER PLOTS)\
Idents(rna) <- "seurat_clusters"\
\
# Remove low-quality cluster if needed\
# rna <- subset(rna, seurat_clusters != 6)\
\
# Define cluster to cell type mapping\
# IMPORTANT: Adjust these numbers based on your FeaturePlot results!\
new_cluster_ids <- c(\
  "0" = "IDN",\
  "1" = "LNP",\
  "2" = "LNP",\
  "3" = "ENP",\
  "4" = "DAN",\
  "5" = "NEP",\
  "7" = "PFPP",\
  "8" = "LNP",\
  "9" = "DAN",\
  "10" = "EPN",\
  "11" = "NEP"\
)\
\
# Apply annotation\
rna <- RenameIdents(rna, new_cluster_ids)\
rna@meta.data$BroadCellType <- Idents(rna)\
\
# Define color palette\
celltype_cols <- c(\
  LNP = "#4CA3DD",\
  ENP = "#A0D8F1",\
  EPN = "#7CB67C",\
  NEP = "#FFD700",\
  PFPP = "#C8B5D0",\
  DAN = "#FF6F61",\
  IDN = "#F7CAC9"\
)\
\
# Save color palette\
saveRDS(celltype_cols, "config/celltype_colors.rds")\
\
# Create annotated UMAP\
print("Creating annotated UMAP...")\
p <- DimPlot(rna, \
             reduction = "umap", \
             cols = celltype_cols,\
             label = TRUE,\
             pt.size = 0.5) + \
  ggtitle("Cell Type Annotations")\
ggsave("results/figures/04_umap_annotated.png", p, width = 10, height = 8)\
\
# Save\
print("Saving annotated object...")\
saveRDS(rna, "results/rds_objects/rna_annotated_final.rds")\
\
print("\uc0\u10003  Annotation complete!")\
print(table(rna$BroadCellType))}