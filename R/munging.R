#' Prepare Analysis DataFrame
#'
#' Extracts expression data and phenotype data, merges them, 
#' and reshapes from Long to Wide to calculate (Post - Pre).
#'
#' @param eset An ExpressionSet object from GEOquery.
#'
#' @return A data frame with Subject_ID, Age, and columns for gene differences.
process_gene_data <- function(eset) {
  
  message("Processing data... (This may take a moment for all genes)")
  
  # 1. Extract Components
  expr_matrix <- exprs(eset)
  pdata <- pData(eset)
  
  # 2. Clean Phenotype Data
  # Extract patient ID and numeric Age.
  # Note: Adjust column names based on actual pData structure inspect
  clean_pdata <- pdata %>%
    select(geo_accession, 
           subject_id = `patientid:ch1`, 
           timepoint = `time:ch1`, 
           age = `age:ch1`) %>%
    mutate(
      age = as.numeric(age),
      # Ensure subject_id is a factor or character to prevent issues
      subject_id = as.character(subject_id) 
    )
  
  # 3. Transpose Expression Matrix (Samples as Rows)
  # Result: Rows = Samples, Cols = Probes
  expr_t <- t(expr_matrix) %>% as.data.frame()
  expr_t$geo_accession <- rownames(expr_t)
  
  # 4. Merge Metadata with Expression
  full_df <- left_join(clean_pdata, expr_t, by = "geo_accession")
  
  # 5. Pivot Wider: Subject as ID, split by Timepoint
  # We exclude geo_accession here as we merge by subject
  gene_cols <- setdiff(colnames(full_df), c("geo_accession", "subject_id", "timepoint", "age"))
  
  message("   ... Pivoting to wide format (Pre vs Post)")
  
  wide_df <- full_df %>%
    select(-geo_accession) %>%
    pivot_wider(
      id_cols = c(subject_id, age),
      names_from = timepoint,
      values_from = all_of(gene_cols),
      names_glue = "{.value}_{timepoint}" # e.g., 211980_at_pre-training
    )
  
  # 6. Calculate Differences (Post - Pre)
  # We identify columns ending in "_pre-training" and match them to "_post-training"
  
  # Identify gene probe names from the original list
  probes <- gene_cols
  
  # Create a matrix for fast subtraction
  # Find indices
  pre_suffix <- "_pre-training"
  post_suffix <- "_post-training"
  
  message("   ... Calculating (Post - Pre) differences")
  
  # Initialize a list to store diff vectors
  diff_list <- list()
  
  # Base R matrix subtraction for speed.
  
  pre_cols <- paste0(probes, pre_suffix)
  post_cols <- paste0(probes, post_suffix)
  
  # Ensure columns exist
  valid_indices <- pre_cols %in% colnames(wide_df) & post_cols %in% colnames(wide_df)
  valid_probes <- probes[valid_indices]
  
  # Create matrices
  mat_pre <- as.matrix(wide_df[, paste0(valid_probes, pre_suffix)])
  mat_post <- as.matrix(wide_df[, paste0(valid_probes, post_suffix)])
  
  mat_diff <- mat_post - mat_pre
  colnames(mat_diff) <- paste0(valid_probes, "_diff")
  
  # Combine back with metadata
  final_df <- cbind(
    wide_df[, c("subject_id", "age")],
    as.data.frame(mat_diff)
  )
  
  # Remove rows with NA (e.g. if a subject missed a post-test)
  final_df <- na.omit(final_df)
  
  message(paste("Data Munging Complete. Subjects:", nrow(final_df), "Genes:", length(valid_probes)))
  
  return(final_df)
}