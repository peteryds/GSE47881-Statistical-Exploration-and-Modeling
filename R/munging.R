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
  # FIX: Use dplyr::select explicitly to avoid conflict with AnnotationDbi::select
  clean_pdata <- pdata %>%
    dplyr::select(geo_accession, 
                  subject_id = `patientid:ch1`, 
                  timepoint = `time:ch1`, 
                  age = `age:ch1`) %>%
    mutate(
      age = as.numeric(age),
      subject_id = as.character(subject_id) 
    )
  
  # 3. Transpose Expression Matrix (Samples as Rows)
  expr_t <- t(expr_matrix) %>% as.data.frame()
  expr_t$geo_accession <- rownames(expr_t)
  
  # 4. Merge Metadata with Expression
  full_df <- left_join(clean_pdata, expr_t, by = "geo_accession")
  
  # 5. Pivot Wider: Subject as ID, split by Timepoint
  gene_cols <- setdiff(colnames(full_df), c("geo_accession", "subject_id", "timepoint", "age"))
  
  message("   ... Pivoting to wide format (Pre vs Post)")
  
  # FIX: Use dplyr::select explicitly here as well
  wide_df <- full_df %>%
    dplyr::select(-geo_accession) %>%
    pivot_wider(
      id_cols = c(subject_id, age),
      names_from = timepoint,
      values_from = all_of(gene_cols),
      names_glue = "{.value}_{timepoint}"
    )
  
  # 6. Calculate Differences (Post - Pre)
  probes <- gene_cols
  pre_suffix <- "_pre-training"
  post_suffix <- "_post-training"
  
  message("   ... Calculating (Post - Pre) differences")
  
  pre_cols <- paste0(probes, pre_suffix)
  post_cols <- paste0(probes, post_suffix)
  
  valid_indices <- pre_cols %in% colnames(wide_df) & post_cols %in% colnames(wide_df)
  valid_probes <- probes[valid_indices]
  
  mat_pre <- as.matrix(wide_df[, paste0(valid_probes, pre_suffix)])
  mat_post <- as.matrix(wide_df[, paste0(valid_probes, post_suffix)])
  
  mat_diff <- mat_post - mat_pre
  colnames(mat_diff) <- paste0(valid_probes, "_diff")
  
  final_df <- cbind(
    wide_df[, c("subject_id", "age")],
    as.data.frame(mat_diff)
  )
  
  final_df <- na.omit(final_df)
  
  message(paste("[DONE] Data Munging Complete. Subjects:", nrow(final_df), "Genes:", length(valid_probes)))
  
  return(final_df)
}