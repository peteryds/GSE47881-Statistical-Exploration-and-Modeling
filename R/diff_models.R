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