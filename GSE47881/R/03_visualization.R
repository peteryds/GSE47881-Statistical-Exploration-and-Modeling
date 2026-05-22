# ============================================================
# Visualization of Results
# ============================================================

message("\n=== STEP 4: Visualizing Significant Findings ===")

# 4.1 Volcano Plot (Main Effect)
volcano_plot <- plot_volcano(limma_res$full_results, p_cutoff = 0.05)
ggsave(file.path(output_dir, "Volcano_Main_Effect.png"), plot = volcano_plot, width = 8, height = 6)

# 4.2 Detailed Plots for Top Interaction Genes
final_df_viz <- process_gene_data(eset_clean) 
top_genes_probes <- rownames(limma_int$sig_genes_df)

if (length(top_genes_probes) == 0) {
  message("No significant interaction genes found (FDR < 0.05). Plotting top 10 by P-value.")
  top_genes_probes <- rownames(head(limma_int$full_results, 10))
} else {
  top_genes_probes <- head(top_genes_probes, 20)
}

message(paste("Generating plots for", length(top_genes_probes), "interaction candidates..."))

for (probe in top_genes_probes) {
  gene_symbol <- limma_int$full_results[probe, "Gene_Symbol"]
  if (is.na(gene_symbol)) gene_symbol <- probe 
  
  p_scatter <- plot_gene_age_scatter(final_df_viz, probe)
  if (!is.null(p_scatter)) {
    p_scatter <- p_scatter + labs(title = paste("Interaction:", gene_symbol), subtitle = paste("Probe:", probe, "| Selected via Limma (Age * Time)"))
    ggsave(file.path(dir_plots, paste0("Scatter_", gene_symbol, "_", probe, ".png")), p_scatter, width = 5, height = 4)
  }
  
  p_violin <- plot_gene_violin(eset_clean, probe)
  if (!is.null(p_violin)) {
    p_violin <- p_violin + labs(title = paste("Change:", gene_symbol))
    ggsave(file.path(dir_plots, paste0("Violin_", gene_symbol, "_", probe, ".png")), p_violin, width = 5, height = 4)
  }
}