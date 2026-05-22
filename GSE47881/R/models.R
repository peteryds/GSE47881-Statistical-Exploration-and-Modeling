#' Run Intercept Model (Gene Change ~ 1)
#'
#' Fits a linear model to test if the mean difference is significantly different from 0.
#' (Intercept model).
#'
#' @param df The clean dataframe containing `_diff` columns.
#' @param genes Vector of probe IDs to test (e.g., c("211980_at", ...)).
#'
#' @return A data frame of results.
run_intercept_models <- function(df, genes) {
  
  message("Running Intercept Models (Mean Change)...")
  
  results_list <- lapply(genes, function(g) {
    # 1. Determine column name
    # The names produced by our munging step are typically "211980_at_diff"
    col_name <- paste0(g, "_diff") 
    
    # Safety check: If not found, try with 'X' prefix (R sometimes adds X to numeric names)
    if (!col_name %in% colnames(df)) {
      col_name <- paste0("X", g, "_diff")
    }
    
    if (!col_name %in% colnames(df)) {
      warning(paste("Column not found for:", g))
      return(NULL)
    }
    
    # 2. Build formula
    # FIX: Use backticks around column name to prevent errors with numeric starts
    # (e.g. use `211980_at_diff` ~ 1 instead of 211980_at_diff ~ 1)
    f <- as.formula(paste0("`", col_name, "` ~ 1"))
    
    model <- lm(f, data = df)
    s <- summary(model)
    
    data.frame(
      Gene = g,
      Model = "Intercept_Only",
      Estimate_Mean_Diff = coef(s)[1, 1], # Intercept estimate
      SE = coef(s)[1, 2],
      T_Value = coef(s)[1, 3],
      P_Value = coef(s)[1, 4]
    )
  })
  
  do.call(rbind, results_list)
}

#' Run Age Interaction Model (Gene Change ~ Age)
#'
#' Fits a linear model to test if Age affects the change in expression.
#'
#' @param df The clean dataframe containing `_diff` columns.
#' @param genes Vector of probe IDs to test.
#'
#' @return A data frame of results.
run_age_models <- function(df, genes) {
  
  message("Running Age Models (Change ~ Age)...")
  
  results_list <- lapply(genes, function(g) {
    col_name <- paste0(g, "_diff")
    
    if (!col_name %in% colnames(df)) {
      col_name <- paste0("X", g, "_diff")
    }
    
    if (!col_name %in% colnames(df)) return(NULL)
    
    # Formula: Diff ~ Age
    # FIX: Use backticks here as well
    f <- as.formula(paste0("`", col_name, "` ~ age"))
    
    model <- lm(f, data = df)
    s <- summary(model)
    
    # Extract Age term
    data.frame(
      Gene = g,
      Model = "Age_Effect",
      Intercept = coef(s)[1, 1],
      Slope_Age = coef(s)[2, 1], # Age slope
      SE_Age = coef(s)[2, 2],
      P_Value_Age = coef(s)[2, 4]
    )
  })
  
  do.call(rbind, results_list)
}

#' Run Limma Screening for Paired Data
#'
#' Uses limma to perform a genome-wide paired analysis (Post vs Pre).
#' This is much more efficient and statistically robust than looping lm().
#'
#' @param eset The ExpressionSet object containing raw expression and phenotype data.
#' @param p_cutoff Adjusted P-value cutoff (FDR) for significance.
#' @param lfc_cutoff Log Fold Change cutoff (absolute value).
#'
#' @return A list containing:
#'   - top_table: A dataframe of all genes sorted by significance.
#'   - sig_genes: A vector of probe IDs that passed the threshold.
run_limma_screening <- function(eset, p_cutoff = 0.05, lfc_cutoff = 0) {
  
  message("Running limma genome-wide screening...")
  
  # 1. Prepare Data
  # Extract Phenotype Data
  pdata <- pData(eset)
  subjects <- factor(pdata$`patientid:ch1`)
  timepoints <- factor(pdata$`time:ch1`, levels = c("pre-training", "post-training"))
  
  # 2. Design Matrix (Paired Design)
  # ~ subjects + timepoints
  # This accounts for baseline differences between subjects (pairing)
  # and tests the effect of time (Post vs Pre).
  design <- model.matrix(~ subjects + timepoints)
  
  # 3. Fit Linear Model
  fit <- lmFit(exprs(eset), design)
  
  # 4. Empirical Bayes Smoothing
  fit <- eBayes(fit)
  
  # 5. Extract Results (Top Table)
  # The coefficient for timepoints usually is the last column if we used the formula above.
  # Let's verify the column name for "post-training".
  coef_name <- grep("timepoints", colnames(design), value = TRUE)
  
  all_results <- topTable(fit, coef = coef_name, number = Inf, adjust.method = "BH")
  
  # 6. Filter Significant Genes
  sig_results <- all_results %>%
    filter(adj.P.Val < p_cutoff & abs(logFC) > lfc_cutoff)
  
  message(paste("Limma Complete."))
  message(paste("Total Genes Tested:", nrow(all_results)))
  message(paste("Significant Genes (FDR <", p_cutoff, "):", nrow(sig_results)))
  
  return(list(
    full_results = all_results,
    sig_genes_df = sig_results,
    sig_probes = rownames(sig_results)
  ))
}

#' Run Limma Interaction Analysis (Timepoint x Age)
#'
#' Fixes "Coefficients not estimable" by manually creating the interaction term.
#' Model: ~ subjects + timepoints + age_interaction_term (where age_interaction_term = is_post * age_centered)
#' This tests if the *Change* (Post - Pre) is associated with Age.
#'
#' @param eset ExpressionSet object.
#' @param p_cutoff FDR cutoff.
#' @return List of results.
run_limma_interaction <- function(eset, p_cutoff = 0.05) {
  message("Running Limma Interaction Model (Timepoint * Age)...")
  
  pdata <- pData(eset)
  
  # --- 1. Column Detection ---
  col_time <- grep("time", colnames(pdata), value = TRUE, ignore.case = TRUE)[1]
  col_subj <- grep("patient|subject", colnames(pdata), value = TRUE, ignore.case = TRUE)[1]
  col_age  <- grep("age", colnames(pdata), value = TRUE, ignore.case = TRUE)[1]
  
  # Validation for required columns
  if (is.na(col_time)) stop("Error: Time column not found.")
  if (is.na(col_age)) stop("Error: Age column not found.")
  
  # --- 2. Check for Unpaired Samples ---
  subject_counts <- table(pdata[[col_subj]])
  orphans <- names(subject_counts[subject_counts != 2])
  
  if (length(orphans) > 0) {
    stop(paste("Error: Unpaired subjects detected (should have been removed upstream):", paste(orphans, collapse=", ")))
  }
  
  # All subjects are paired; proceed
  pdata_clean <- pdata
  eset_clean <- eset
  
  # --- 3. Prepare Variables ---
  subjects <- factor(pdata_clean[[col_subj]])
  
  # Timepoints: Ensure Pre is Reference (0)
  timepoints <- factor(pdata_clean[[col_time]])
  ref_level <- grep("pre|base|0", levels(timepoints), value = TRUE, ignore.case = TRUE)[1]
  if (!is.na(ref_level)) {
    timepoints <- relevel(timepoints, ref = ref_level)
  }
  
  # Age: Clean & Center
  age_numeric <- as.numeric(gsub("[^0-9.]", "", pdata_clean[[col_age]]))
  if (any(is.na(age_numeric))) {
    stop("Error: Missing or invalid age values detected in the data.")
  }
  age_centered <- age_numeric - mean(age_numeric)
  
  # --- 4. Manual Interaction Term (THE FIX) ---
  # We create a numeric vector that is:
  # - 0 for Pre-training samples
  # - Age_Centered for Post-training samples
  # This avoids collinearity with Subjects.
  
  # Create numeric 0/1 for Time (0=Ref/Pre, 1=Treat/Post)
  is_post <- as.numeric(timepoints) - 1 
  
  # Interaction Vector
  age_interaction_term <- is_post * age_centered
  
  # --- 5. Design Matrix ---
  # Formula: ~ Subject + Time + Manual_Interaction
  design <- model.matrix(~ subjects + timepoints + age_interaction_term)
  
  # Rename the interaction column to something readable
  # usually it gets named "age_interaction_term"
  # First, make all column names safe
  colnames(design) <- make.names(colnames(design))
  # Then, explicitly set the last column name to "Interaction_Age"
  colnames(design)[ncol(design)] <- "Interaction_Age"
  
  # --- 6. Run Limma ---
  fit <- lmFit(exprs(eset_clean), design)
  fit <- eBayes(fit)
  
  # --- 7. Extract Results ---
  # Use the actual name of the last column for the coefficient
  coef_name <- colnames(design)[ncol(design)]
  
  if (!coef_name %in% colnames(design)) {
    # Fallback search if renaming failed
    coef_name <- grep("age", colnames(design), value = TRUE, ignore.case = TRUE)
    coef_name <- tail(coef_name, 1) # usually the last one
  }
  
  message(paste("Targeting Coefficient:", coef_name))
  
  full_results <- topTable(fit, coef = coef_name, number = Inf, adjust.method = "BH")
  sig_genes <- full_results %>% filter(adj.P.Val < p_cutoff)
  
  message(paste("Analysis Complete. Significant Interaction Genes:", nrow(sig_genes)))
  
  return(list(
    full_results = full_results,
    sig_genes_df = sig_genes
  ))
}