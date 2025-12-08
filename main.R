# ============================================================
# Gene Expression Analysis: Resistance Training & Age
# Main Execution Script
# ============================================================

# 1. Load Environment & Functions
source("R/load_packages.R")
source("R/load_data.R")
source("R/munging.R")
source("R/diff_models.R")
source("R/eda.R") 

# Setup libraries
setup_environment()

output_dir <- "output"
if (!dir.exists(output_dir)) dir.create(output_dir)
eda_dir <- file.path(output_dir, "EDA_Global")
if (!dir.exists(eda_dir)) dir.create(eda_dir)

# ============================================================
# 2. Data Loading & Global QC (Part 1 EDA)
# ============================================================

eset <- get_geo_data(gse_id = "GSE47881")

message("Generating Global QC Plots...")

# 2.1 PCA Plot
pca_plot <- plot_pca(eset, title = "PCA: All Samples")
ggsave(file.path(eda_dir, "QC_Global_PCA.png"), plot = pca_plot, width = 8, height = 6)

# 2.2 Sample Distance Heatmap
# Note: pheatmap draws directly to file using png() device
png(file.path(eda_dir, "QC_Sample_Distance_Heatmap.png"), width = 800, height = 800)
plot_sample_heatmap(eset)
dev.off() # Close device

# 2.3 Density Plot
density_plot <- plot_density(eset)
ggsave(file.path(eda_dir, "QC_Density_Plot.png"), plot = density_plot, width = 8, height = 6)


# ============================================================
# 3. Genome-Wide Screening (LIMMA) & Results EDA (Part 2 EDA)
# ============================================================

# Run limma (Paired Analysis)
limma_res <- run_limma_screening(eset, p_cutoff = 0.05)

# Save Results
write.csv(limma_res$full_results, file.path(output_dir, "limma_all_genes.csv"))
write.csv(limma_res$sig_genes_df, file.path(output_dir, "limma_significant_genes.csv"))

message("Generating Significant Results EDA...")

# 3.1 Volcano Plot
volcano_plot <- plot_volcano(limma_res$full_results, p_cutoff = 0.05)
ggsave(file.path(output_dir, "EDA_Volcano_Plot.png"), plot = volcano_plot, width = 8, height = 6)

# 3.2 Heatmap of Top Significant Genes
# Get top 50 significant genes (or fewer if not enough) for heatmap
top_n_heatmap <- 50
sig_genes_all <- rownames(limma_res$sig_genes_df)

if (length(sig_genes_all) > 0) {
  genes_for_heatmap <- head(sig_genes_all, top_n_heatmap)
  
  png(file.path(output_dir, "EDA_Sig_Genes_Heatmap.png"), width = 800, height = 1000)
  plot_sig_heatmap(eset, genes_for_heatmap)
  dev.off()
} else {
  warning("[WARNING] No significant genes for heatmap.")
}


# ============================================================
# 4. Data Munging (Prepare for Age Analysis)
# ============================================================

final_df <- process_gene_data(eset)


# ============================================================
# 5. Detailed Analysis: Targeted Genes
# ============================================================

# Dynamic Selection (Tier 1: FDR < 0.05, Tier 2: P < 0.001)
target_genes <- NULL
if (length(sig_genes_all) > 0) {
  target_genes <- sig_genes_all
} else {
  loose_p <- 0.001
  target_genes <- rownames(limma_res$full_results %>% filter(P.Value < loose_p))
  if (length(target_genes) == 0) target_genes <- rownames(head(limma_res$full_results, 5))
}

cat("\nAnalyzing", length(target_genes), "Target Genes.\n")

# Run Models
results_intercept <- run_intercept_models(final_df, target_genes)
results_age <- run_age_models(final_df, target_genes)

# ============================================================
# 6. Specific Gene Visualization (Scatter & Violin)
# ============================================================

gene_plot_dir <- file.path(output_dir, "gene_plots")
if (!dir.exists(gene_plot_dir)) dir.create(gene_plot_dir)

# Cap at top 20 for plotting
genes_to_plot <- head(target_genes, 5)
message("Generating detailed plots for top ", length(genes_to_plot), " genes...")

for (gene in genes_to_plot) {
  
  # A. Age vs Change Scatter Plot
  p_scatter <- plot_gene_age_scatter(final_df, gene)
  if (!is.null(p_scatter)) {
    ggsave(file.path(gene_plot_dir, paste0(gene, "_Age_Scatter.png")), p_scatter, width = 5, height = 4)
  }
  
  # B. Pre vs Post Violin Plot
  # Note: This uses 'eset' directly, not 'final_df'
  p_violin <- plot_gene_violin(eset, gene)
  if (!is.null(p_violin)) {
    ggsave(file.path(gene_plot_dir, paste0(gene, "_PrePost_Violin.png")), p_violin, width = 5, height = 4)
  }
}

# Export Results
write.csv(results_intercept, file.path(output_dir, "results_intercept_model.csv"), row.names = FALSE)
write.csv(results_age, file.path(output_dir, "results_age_model.csv"), row.names = FALSE)

message(paste("[DONE] Analysis & EDA complete! Check", output_dir))