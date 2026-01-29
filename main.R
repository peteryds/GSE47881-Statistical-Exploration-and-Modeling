# ============================================================
# Gene Expression Analysis: Resistance Training & Age
# Main Execution Script
# ============================================================

# 1. Load Environment & Functions
source("R/load_packages.R")
source("R/load_data.R")
source("R/munging.R")  # Must contain: process_gene_data, clean_and_normalize_data
source("R/models.R")   # Must contain: run_limma_screening, run_limma_interaction
source("R/eda.R")      # Must contain: plotting functions

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

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

BiocManager::install(c("hgu133plus2.db", "AnnotationDbi", "enrichR"))

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


# ============================================================
# 6.5 Annotation (ID Mapping) - NEW SECTION
# ============================================================
message("\n=== STEP 3.5: Mapping Probe IDs to Gene Symbols ===")
library(hgu133plus2.db)
library(AnnotationDbi)

# Define a helper function to perform ID mapping
annotate_results <- function(df) {
  # Safety Check: If dataframe is empty, return it immediately
  if (nrow(df) == 0) {
    message("  [Alert] Dataframe is empty. Skipping annotation.")
    return(df)
  }
  
  # Get Probe IDs (Row names)
  probes <- rownames(df)
  
  # Use mapIds to convert Probes to Gene Symbols
  symbols <- mapIds(hgu133plus2.db,
                    keys = probes,
                    column = "SYMBOL",
                    keytype = "PROBEID",
                    multiVals = "first")
  
  # Add Gene Symbol back to the dataframe
  df$Gene_Symbol <- symbols
  
  # Optionally add full Gene Name
  genenames <- mapIds(hgu133plus2.db,
                      keys = probes,
                      column = "GENENAME",
                      keytype = "PROBEID",
                      multiVals = "first")
  df$Gene_Name <- genenames
  
  # Reorder columns
  df <- df[, c("Gene_Symbol", setdiff(names(df), "Gene_Symbol"))]
  
  return(df)
}

# Execute Annotation
message("Annotating Main Effect Results...")
limma_res$full_results <- annotate_results(limma_res$full_results)
limma_res$sig_genes_df <- annotate_results(limma_res$sig_genes_df)

message("Annotating Interaction Results...")
limma_int$full_results <- annotate_results(limma_int$full_results)
limma_int$sig_genes_df <- annotate_results(limma_int$sig_genes_df)


# Save Annotated Results
# The output CSVs will now contain human-readable Gene Symbols
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
# Use the processed data for visualization
final_df_viz <- process_gene_data(eset_clean) 

# Select top genes: Prioritize significant interaction genes
top_genes_probes <- rownames(limma_int$sig_genes_df)

if (length(top_genes_probes) == 0) {
  message("No significant interaction genes found (FDR < 0.05). Plotting top 10 by P-value.")
  top_genes_probes <- rownames(head(limma_int$full_results, 10))
} else {
  # Limit to top 20 to avoid generating too many files
  top_genes_probes <- head(top_genes_probes, 20)
}

message(paste("Generating plots for", length(top_genes_probes), "interaction candidates..."))

for (probe in top_genes_probes) {
  # Retrieve Gene Symbol for plot titles
  gene_symbol <- limma_int$full_results[probe, "Gene_Symbol"]
  if (is.na(gene_symbol)) gene_symbol <- probe # Fallback to Probe ID if NA
  
  # A. Scatter Plot (Age vs Change)
  p_scatter <- plot_gene_age_scatter(final_df_viz, probe)
  if (!is.null(p_scatter)) {
    p_scatter <- p_scatter + 
      labs(title = paste("Interaction:", gene_symbol),
           subtitle = paste("Probe:", probe, "| Selected via Limma (Age * Time)"))
    ggsave(file.path(dir_plots, paste0("Scatter_", gene_symbol, "_", probe, ".png")), p_scatter, width = 5, height = 4)
  }
  
  # B. Violin Plot (Pre vs Post)
  p_violin <- plot_gene_violin(eset_clean, probe)
  if (!is.null(p_violin)) {
    p_violin <- p_violin + labs(title = paste("Change:", gene_symbol))
    ggsave(file.path(dir_plots, paste0("Violin_", gene_symbol, "_", probe, ".png")), p_violin, width = 5, height = 4)
  }
}


# ============================================================
# 8. Pathway Analysis (Split by Direction: UP vs DOWN)
# ============================================================
message("\n=== STEP 5: Pathway Enrichment Analysis (Directional) ===")
library(enrichR)

# 1. Setup Databases and Thresholds
dbs <- c("KEGG_2021_Human", "GO_Biological_Process_2021")
p_threshold <- 0.05  # Standard cutoff for saving files (you can filter to 0.01 later)

# 2. Prepare Gene Lists (Split by logFC)
# We assume 'limma_res$sig_genes_df' contains a 'logFC' column and 'Gene_Symbol'
sig_df <- limma_res$sig_genes_df

# Check if logFC exists (Safety check)
if (!"logFC" %in% colnames(sig_df)) {
  stop("Error: 'logFC' column not found in results. Cannot split by direction.")
}

# Extract UP-regulated genes (logFC > 0)
genes_up <- unique(na.omit(sig_df$Gene_Symbol[sig_df$logFC > 0]))
genes_up

# Extract DOWN-regulated genes (logFC < 0)
genes_down <- unique(na.omit(sig_df$Gene_Symbol[sig_df$logFC < 0]))

message(paste0("Found ", length(genes_up), " Upregulated genes."))
message(paste0("Found ", length(genes_down), " Downregulated genes."))


# 3. Define Analysis Function
run_directional_enrichment <- function(gene_list, direction_label) {
  
  if (length(gene_list) < 5) {
    message(paste0("\n[Skipping] Not enough genes in ", direction_label, " list (<5)."))
    return(NULL)
  }
  
  message(paste0("\nRunning Enrichment for: ", direction_label, "..."))
  enriched <- enrichr(gene_list, dbs)
  
  # --- Process Each Database ---
  for (db_name in dbs) {
    if (!is.null(enriched[[db_name]])) {
      # 1. Sort by P-value
      res <- enriched[[db_name]]
      res <- res[order(res$P.value), ]
      
      # 2. Filter Significant (Standard P < 0.05)
      res_sig <- res[res$P.value < p_threshold, ]
      
      # 3. Save to CSV (e.g., Pathways_KEGG_UP.csv)
      # Clean DB name for filename (remove year if desired, keeping it simple here)
      clean_db_name <- strsplit(db_name, "_")[[1]][1] # e.g., "KEGG"
      fname <- file.path(output_dir, paste0("Pathways_", clean_db_name, "_", direction_label, ".csv"))
      
      # Only write CSV when there are significant pathways to avoid empty files
      if (nrow(res_sig) > 0) {
        write.csv(res_sig, fname)
      }
      
      # 4. Print Top Results to Console
      message(paste0("  [", clean_db_name, "] Top 3 Significant Pathways:"))
      if (nrow(res_sig) > 0) {
        print(head(res_sig[, c("Term", "P.value", "Overlap")], 3))
      } else {
        message("    No significant pathways found (P < 0.05).")
      }
    }
  }
}

# 4. Execute Analysis
run_directional_enrichment(genes_up, "UP")
run_directional_enrichment(genes_down, "DOWN")

message("\n[DONE] Directional Pathway Analysis Finished!")
message(paste("Check output directory:", output_dir))
