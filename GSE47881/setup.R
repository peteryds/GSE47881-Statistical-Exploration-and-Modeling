# ============================================================
# Setup Environment & Directories
# ============================================================

# 1. Load Environment & Functions
source("R/load_packages.R")
source("R/load_data.R")
source("R/munging.R")  
source("R/models.R")   
source("R/eda.R")      

# Setup libraries
setup_environment()

# 2. Setup Directories (Aligned with README structure)
output_dir <- "results"
if (!dir.exists(output_dir)) dir.create(output_dir)

dir_raw <- file.path(output_dir, "QC_1_Raw")
dir_clean <- file.path(output_dir, "QC_2_Cleaned")
dir_plots <- file.path(output_dir, "Plots_Interaction")

for (d in c(dir_raw, dir_clean, dir_plots)) {
  if (!dir.exists(d)) dir.create(d)
}