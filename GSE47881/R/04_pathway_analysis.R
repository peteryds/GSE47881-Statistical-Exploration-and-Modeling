# ============================================================
# Pathway Analysis (Split by Direction: UP vs DOWN)
# ============================================================

message("\n=== STEP 5: Pathway Enrichment Analysis (Directional) ===")

dbs <- c("KEGG_2021_Human", "GO_Biological_Process_2021")
p_threshold <- 0.05 
sig_df <- limma_res$sig_genes_df

if (is.null(sig_df) || nrow(sig_df) == 0) {
  message("\n[Skipping] No significant genes available for directional pathway analysis.")
} else {
  if (!"logFC" %in% colnames(sig_df)) {
    stop("Error: 'logFC' column not found in results. Cannot split by direction.")
  }
  
  genes_up <- unique(na.omit(sig_df$Gene_Symbol[sig_df$logFC > 0]))
  genes_down <- unique(na.omit(sig_df$Gene_Symbol[sig_df$logFC < 0]))
  
  message(paste0("Found ", length(genes_up), " Upregulated genes."))
  message(paste0("Found ", length(genes_down), " Downregulated genes."))
  
  run_directional_enrichment <- function(gene_list, direction_label) {
    if (length(gene_list) < 5) {
      message(paste0("\n[Skipping] Not enough genes in ", direction_label, " list (<5)."))
      return(NULL)
    }
    
    message(paste0("\nRunning Enrichment for: ", direction_label, "..."))
    enriched <- enrichr(gene_list, dbs)
    
    for (db_name in dbs) {
      if (!is.null(enriched[[db_name]])) {
        res <- enriched[[db_name]]
        res <- res[order(res$P.value), ]
        res_sig <- res[res$P.value < p_threshold, ]
        
        clean_db_name <- strsplit(db_name, "_")[[1]][1] 
        fname <- file.path(output_dir, paste0("Pathways_", clean_db_name, "_", direction_label, ".csv"))
        
        write.csv(res_sig, fname)
        message(paste0("  [", clean_db_name, "] Top Significant Pathways Saved."))
      }
    }
  }
  
  run_directional_enrichment(genes_up, "UP_Activated")
  run_directional_enrichment(genes_down, "DOWN_Inhibited")
  
  message("\n[DONE] Directional Pathway Analysis Finished!")
}