# ============================================================
# Statistical Analysis (Limma) & Annotation
# ============================================================

message("\n=== STEP 3: Running Limma Models ===")

# A. Standard Paired Analysis (Main Effect: Post vs Pre)
limma_res <- run_limma_screening(eset_clean, p_cutoff = 0.05)

# B. Interaction Analysis (Age Effect)
limma_int <- run_limma_interaction(eset_clean, p_cutoff = 0.05)

# ============================================================
# Annotation (ID Mapping)
# ============================================================
message("\n=== Section 3.5: Mapping Probe IDs to Gene Symbols ===")

annotate_results <- function(df) {
  if (nrow(df) == 0) {
    message("  [Alert] Dataframe is empty. Skipping annotation.")
    return(df)
  }
  
  probes <- rownames(df)
  symbols <- mapIds(hgu133plus2.db, keys = probes, column = "SYMBOL", keytype = "PROBEID", multiVals = "first")
  df$Gene_Symbol <- symbols
  
  genenames <- mapIds(hgu133plus2.db, keys = probes, column = "GENENAME", keytype = "PROBEID", multiVals = "first")
  df$Gene_Name <- genenames
  
  df <- df[, c("Gene_Symbol", "Gene_Name", setdiff(names(df), c("Gene_Symbol", "Gene_Name")))]
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
message("Saving analysis results to CSV...")
write.csv(limma_res$full_results, file.path(output_dir, "Results_Main_Effect.csv"))
write.csv(limma_res$sig_genes_df, file.path(output_dir, "Results_Main_Effect_Sig.csv"))
write.csv(limma_int$full_results, file.path(output_dir, "Results_Interaction_Age.csv"))
write.csv(limma_int$sig_genes_df, file.path(output_dir, "Results_Interaction_Age_Sig.csv"))