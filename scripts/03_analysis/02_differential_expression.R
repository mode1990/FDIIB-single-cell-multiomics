{\rtf1\ansi\ansicpg1252\cocoartf2761
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\paperw11900\paperh16840\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs24 \cf0 ################################################################################\
# STEP 5: Differential Expression Analysis\
# Description: Find DEGs between mutations in dopaminergic neurons\
# Input: Annotated RNA object\
# Output: DEG tables\
################################################################################\
\
library(Seurat)\
library(dplyr)\
\
# Load data\
print("Loading annotated data...")\
rna <- readRDS("results/rds_objects/rna_annotated_final.rds")\
\
# Subset to dopaminergic neurons\
print("Subsetting to DAN...")\
Idents(rna) <- "BroadCellType"\
DAN <- subset(rna, idents = "DAN")\
\
# Set identity to mutation\
Idents(DAN) <- "Mutation"\
\
# Prepare for DE testing\
print("Preparing for DE testing...")\
DAN <- PrepSCTFindMarkers(DAN, verbose = FALSE)\
\
# Run DE analysis for each mutation vs HC\
mutations <- c("LRRK2", "SNCA", "GBA1")\
\
for (mut in mutations) \{\
  print(paste("Running DE for:", mut, "vs HC"))\
  \
  markers <- FindMarkers(\
    DAN,\
    recorrect_umi = FALSE,\
    ident.1 = mut,\
    ident.2 = "HC",\
    min.pct = 0.25,\
    logfc.threshold = 0.1,\
    verbose = FALSE\
  ) %>%\
    subset(p_val_adj < 0.05) %>%\
    transform(Gene = rownames(.))\
  \
  # Save results\
  output_file <- paste0("results/tables/DEGs_DAN_", mut, "_vs_HC.csv")\
  write.csv(markers, output_file, row.names = FALSE)\
  \
  print(paste("  Found", nrow(markers), "DEGs"))\
  print(paste("  Saved to:", output_file))\
\}\
\
print("\uc0\u10003  Differential expression complete!")}