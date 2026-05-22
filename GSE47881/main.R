# ============================================================
# GSE47881 Master Execution Script
# Run this script with working directory set to GSE47881/
# ============================================================

# 1. Setup Environment & Directories
message("=== Initializing GSE47881 Project ===")
source("src/setup.R")

message("\n=== Starting GSE47881 Pipeline ===")
source("src/01_data_prep_and_qc.R")
source("src/02_limma_analysis.R")
source("src/03_visualization.R")
source("src/04_pathway_analysis.R")

message(paste("\n[DONE] Pipeline execution complete. Check the directory:", output_dir))