# ============================================================
# Gene Expression Analysis: Resistance Training & Age
# Main Execution Script (Final Updated Version)
# ============================================================

# 1. Load Environment & Functions
source("R/load_packages.R")
source("R/load_data.R")
source("R/munging.R")  # Must contain: process_gene_data, clean_and_normalize_data
source("R/models.R")   # Must contain: run_limma_screening, run_limma_interaction
source("R/eda.R")      # Must contain: plotting functions
source("R/Descriptive_analysis.R")

# Setup libraries
setup_environment()

# 2. Setup Directories
output_dir <- "output"
if (!dir.exists(output_dir)) dir.create(output_dir)

dir_raw <- file.path(output_dir, "QC_1_Raw")
if (!dir.exists(dir_raw)) dir.create(dir_raw)

dir_clean <- file.path(output_dir, "QC_2_Cleaned")
if (!dir.exists(dir_clean)) dir.create(dir_clean)

dir_plots <- file.path(output_dir, "Plots_Interaction")
if (!dir.exists(dir_plots)) dir.create(dir_plots)

# ============================================================
# 3. Data Loading & Phase 1 QC (Raw Data)
# ============================================================
message("\n=== STEP 1: Loading Raw Data ===")
eset_raw <- get_geo_data(gse_id = "GSE47881")

message("Generating Phase 1 QC Plots (Raw Data)...")

# 3.1 PCA (Raw)
# This will likely show the outlier and lack of clustering
p_pca_raw <- plot_pca(eset_raw, title = "PCA: Raw Data (Before QC)")
ggsave(file.path(dir_raw, "QC_Raw_PCA.png"), plot = p_pca_raw, width = 8, height = 6)

# 3.2 Heatmap (Raw)
# Use png() device for pheatmap
png(file.path(dir_raw, "QC_Raw_Heatmap.png"), width = 800, height = 800)
plot_sample_heatmap(eset_raw)
dev.off()

# 3.3 Density (Raw)
# This will likely show the linear scale (not Log2) and the "spike" at 0
p_dens_raw <- plot_density(eset_raw) + 
  labs(subtitle = "Raw Data: Note Linear Scale & Outlier Spike")
ggsave(file.path(dir_raw, "QC_Raw_Density.png"), plot = p_dens_raw, width = 8, height = 6)


# ============================================================
# 4. Data Munging (Cleaning & Normalization)
# ============================================================
message("\n=== STEP 2: Data Cleaning & Normalization ===")

# This function performs:
# 1. Log2 transformation (if data is raw)
# 2. Removal of specific outliers (optional)
# 3. Removal of unpaired 'orphan' subjects (Crucial for Paired Analysis)
# NOTE: If you identified specific GSM IDs to remove from the Raw Heatmap, add them here.
# e.g., outliers = c("GSM1161833")

eset_clean <- clean_and_normalize_data(eset_raw, outliers_to_remove = NULL)


# ============================================================
# 5. Phase 2 QC (Cleaned Data)
# ============================================================
message("Generating Phase 2 QC Plots (Cleaned Data)...")

# 5.1 PCA (Clean)
# Outlier should be gone, scale should be normalized
p_pca_clean <- plot_pca(eset_clean, title = "PCA: Cleaned & Log2 Transformed")
ggsave(file.path(dir_clean, "QC_Clean_PCA.png"), plot = p_pca_clean, width = 8, height = 6)

# 5.2 Heatmap (Clean)
png(file.path(dir_clean, "QC_Clean_Heatmap.png"), width = 800, height = 800)
plot_sample_heatmap(eset_clean)
dev.off()

# 5.3 Density (Clean)
# Curves should overlap (Bell shape)
p_dens_clean <- plot_density(eset_clean) + 
  labs(subtitle = "Cleaned Data: Log2 Scale & Normalized")
ggsave(file.path(dir_clean, "QC_Clean_Density.png"), plot = p_dens_clean, width = 8, height = 6)


# ============================================================
# 6. Statistical Analysis (Limma)
# ============================================================
message("\n=== STEP 3: Running Limma Models ===")

# A. Standard Paired Analysis (Main Effect: Post vs Pre)
# Tests: Does exercise change gene expression (on average)?
limma_res <- run_limma_screening(eset_clean, p_cutoff = 0.05)

# B. Interaction Analysis (Age Effect)
# Tests: Does Age affect the MAGNITUDE of change? (Timepoint * Age)
limma_int <- run_limma_interaction(eset_clean, p_cutoff = 0.05)

# Save Results
write.csv(limma_res$full_results, file.path(output_dir, "Results_Main_Effect.csv"))
write.csv(limma_res$sig_genes_df, file.path(output_dir, "Results_Main_Effect_Sig.csv"))
write.csv(limma_int$full_results, file.path(output_dir, "Results_Interaction_Age.csv"))
write.csv(limma_int$sig_genes_df, file.path(output_dir, "Results_Interaction_Age_Sig.csv"))


# ============================================================
# 7. Visualization of Results
# ============================================================
message("\n=== STEP 4: Visualizing Significant Findings ===")

# 7.1 Volcano Plot (Main Effect)
volcano_plot <- plot_volcano(limma_res$full_results, p_cutoff = 0.05)
ggsave(file.path(output_dir, "Volcano_Main_Effect.png"), plot = volcano_plot, width = 8, height = 6)

# 7.2 Detailed Plots for Top Interaction Genes
# We need the wide-format dataframe (diff) for plotting the scatter plots
# We generate this strictly for visualization purposes
final_df_viz <- process_gene_data(eset_clean) 

# Select top genes: 
# Priority 1: Significant Interaction Genes (FDR < 0.05)
# Priority 2: Top 10 genes by P-value (if no sig genes)
top_genes <- rownames(limma_int$sig_genes_df)
if (length(top_genes) == 0) {
  message("No significant interaction genes found (FDR < 0.05). Plotting top 10 by P-value.")
  top_genes <- rownames(head(limma_int$full_results, 10))
} else {
  # Limit to top 20 to avoid too many files
  top_genes <- head(top_genes, 20)
}

message(paste("Generating plots for", length(top_genes), "interaction candidates..."))

for (gene in top_genes) {
  # A. Scatter Plot (Age vs Change) - Proves the Interaction
  p_scatter <- plot_gene_age_scatter(final_df_viz, gene)
  if (!is.null(p_scatter)) {
    # Add subtitle to indicate this was found via Limma
    p_scatter <- p_scatter + labs(subtitle = "Selected via Limma Interaction (Age * Time)")
    ggsave(file.path(dir_plots, paste0("Scatter_", gene, ".png")), p_scatter, width = 5, height = 4)
  }
  
  # B. Violin Plot (Pre vs Post) - Shows the raw change
  p_violin <- plot_gene_violin(eset_clean, gene)
  if (!is.null(p_violin)) {
    ggsave(file.path(dir_plots, paste0("Violin_", gene, ".png")), p_violin, width = 5, height = 4)
  }
}

message("\n[DONE] Pipeline Finished Successfully!")
message(paste("Check output directory:", output_dir))

# ============================================================
# 3. Raw Data Summary & QC
# ============================================================

message("\n===== RAW DATA SUMMARY & QC =====")
raw_summary <- summarize_raw_data(eset, output_dir = "output/raw_qc")

# Print summary statistics
print(raw_summary$summary_stats)
cat("Global Mean: ", raw_summary$global_mean, "\n")
cat("Global SD: ", raw_summary$global_sd, "\n")
cat("Pre Mean: ", raw_summary$pre_mean, "\n")
cat("Pre SD: ", raw_summary$pre_sd, "\n")
cat("Post Mean: ", raw_summary$post_mean, "\n")
cat("Post SD: ", raw_summary$post_sd, "\n")