################################################################################
# Install Required R Packages
# Run this once: Rscript scripts/install_packages.R
################################################################################

cat("Installing R packages for FOUNDIN-PD pipeline...\n\n")

# Install BiocManager if needed
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Read requirements
requirements <- readLines("requirements.txt")
requirements <- requirements[!grepl("^#", requirements) & nzchar(requirements)]

# Extract package names
packages <- sub(">=.*", "", requirements)
packages <- trimws(packages)

# Bioconductor packages (start with capital letter or contain ".")
bioc_packages <- c(
  "BSgenome.Hsapiens.UCSC.hg38", "EnsDb.Hsapiens.v86",
  "GenomicFeatures", "GenomicRanges", "rtracklayer", "Biostrings",
  "scDblFinder", "SingleCellExperiment", "JASPAR2020", "TFBSTools",
  "ComplexHeatmap"
)

# CRAN packages
cran_packages <- setdiff(packages, bioc_packages)

# Install CRAN packages
cat("Installing CRAN packages...\n")
for (pkg in cran_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(paste("Installing", pkg, "\n"))
    install.packages(pkg, dependencies = TRUE)
  }
}

# Install Bioconductor packages
cat("\nInstalling Bioconductor packages...\n")
for (pkg in bioc_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(paste("Installing", pkg, "\n"))
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  }
}

# Install MOFA2 from GitHub
cat("\nInstalling MOFA2 from GitHub...\n")
if (!require("MOFA2", quietly = TRUE)) {
  remotes::install_github("bioFAM/MOFA2", subdir = "MOFA2")
}

cat("\n✓ All packages installed successfully!\n")
cat("\nVerifying key packages:\n")
key_packages <- c("Seurat", "Signac", "harmony", "MOFA2", "scDblFinder")
for (pkg in key_packages) {
  version <- packageVersion(pkg)
  cat(paste("  ", pkg, "version:", version, "\n"))
}
