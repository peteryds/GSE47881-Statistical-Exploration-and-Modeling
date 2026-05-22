# ============================================================
# Data Loading & QC
# ============================================================

message("\n=== STEP 1: Loading Raw Data ===")
eset_raw <- get_geo_data(gse_id = "GSE47881")

message("Generating Phase 1 QC Plots (Raw Data)...")

# 1.1 PCA (Raw)
p_pca_raw <- plot_pca(eset_raw, title = "PCA: Raw Data (Before QC)")
ggsave(file.path(dir_raw, "QC_Raw_PCA.png"), plot = p_pca_raw, width = 8, height = 6)

# 1.2 Heatmap (Raw)
png(file.path(dir_raw, "QC_Raw_Heatmap.png"), width = 800, height = 800)
plot_sample_heatmap(eset_raw)
dev.off()

# 1.3 Density (Raw)
p_dens_raw <- plot_density(eset_raw) + 
  labs(subtitle = "Raw Data: Note Linear Scale & Outlier Spike")
ggsave(file.path(dir_raw, "QC_Raw_Density.png"), plot = p_dens_raw, width = 8, height = 6)

# ============================================================
# Data Munging (Cleaning & Normalization)
# ============================================================
message("\n=== STEP 2: Data Cleaning & Normalization ===")
eset_clean <- clean_and_normalize_data(eset_raw, outliers_to_remove = NULL)

message("Generating Phase 2 QC Plots (Cleaned Data)...")

# 2.1 PCA (Clean)
p_pca_clean <- plot_pca(eset_clean, title = "PCA: Cleaned & Log2 Transformed")
ggsave(file.path(dir_clean, "QC_Clean_PCA.png"), plot = p_pca_clean, width = 8, height = 6)

# 2.2 Heatmap (Clean)
png(file.path(dir_clean, "QC_Clean_Heatmap.png"), width = 800, height = 800)
plot_sample_heatmap(eset_clean)
dev.off()

# 2.3 Density (Clean)
p_dens_clean <- plot_density(eset_clean) + 
  labs(subtitle = "Cleaned Data: Log2 Scale & Normalized")
ggsave(file.path(dir_clean, "QC_Clean_Density.png"), plot = p_dens_clean, width = 8, height = 6)