# ============================================================
# Gene Expression Analysis
# Main Execution Script
# ============================================================

# ============================================================
# 1. Load Environment & Functions
# ============================================================

source("R/load_packages.R")
source("R/load_data.R")
source("R/munging.R")
source("R/diff_models.R")

# Setup libraries
setup_environment()

# ============================================================
# 2. Data loading and munging 
# ============================================================

# Download/Load Data
eset <- get_geo_data(gse_id = "GSE47881")

# Process Data (Long -> Wide -> Calculate Diffs)
# This creates a master dataframe with diff scores for all the genes
final_df <- process_gene_data(eset)

# Preview data
print(head(final_df[, 1:5]))

# ============================================================
# 3. Analysis: Targeted Gene List
# ============================================================

# Define the 8 Genes of Interest
target_genes <- c(
  "211980_at", "204114_at", "212013_at", "204008_at",
  "218429_s_at", "203477_at", "205656_at", "204115_at"
)

# --- Model A: Is there a significant change? (Intercept only) ---
results_intercept <- run_intercept_models(final_df, target_genes)

cat("\n===== Model A Results: Gene Change (Intercept) =====\n")
print(results_intercept)


# --- Model B: Does Age affect the change? (Diff ~ Age) ---
results_age <- run_age_models(final_df, target_genes)

cat("\n===== Model B Results: Age Effect =====\n")
print(results_age)


# ============================================================
# 4. Export Results
# ============================================================

# Define output directory
output_dir <- "output"

# Create directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
  message(paste("📂 Created output directory:", output_dir))
}

# Construct file paths
file_path_intercept <- file.path(output_dir, "results_intercept_model.csv")
file_path_age <- file.path(output_dir, "results_age_model.csv")

# Write CSVs
write.csv(results_intercept, file_path_intercept, row.names = FALSE)
write.csv(results_age, file_path_age, row.names = FALSE)

message(paste("Results saved to", output_dir))