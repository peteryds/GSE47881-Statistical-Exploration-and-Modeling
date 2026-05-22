# =========================================================================
# Script: R/02_lasso_predictive_modeling.R
# Objective: Phase 2 - Matrix Transposition & LASSO Logistic Regression
# =========================================================================
library(glmnet)
library(dplyr)

cat("Step 1: Constructing unique feature names...\n")
# Combine assay, tissue, and feature_ID to create unique feature names
feature_names <- paste(
  multi_omics_matrix$assay,
  multi_omics_matrix$tissue,
  multi_omics_matrix$feature_ID,
  sep = "_"
)

cat("Step 2: Extracting numeric matrix and transposing...\n")
# Remove the first 4 metadata columns to keep the pure numeric matrix
numeric_matrix <- as.matrix(multi_omics_matrix[, -(1:4)])
rownames(numeric_matrix) <- feature_names

# 轉置矩陣：讓 Rows = 樣本 (pid), Columns = 特徵 (81565)
X_transposed <- t(numeric_matrix) 

cat("Step 3: Aligning Features (X) with Target (Y)...\n")
# Ensure the sample rows in the feature matrix align with our target_y vector
shared_pids <- intersect(rownames(X_transposed), names(target_y))
X_model <- X_transposed[shared_pids, ]
Y_model <- target_y[shared_pids]

cat("Final modeling cohort:", length(shared_pids), "rats.\n")

cat("Step 4: Handling missing values (Imputation)...\n")
# Multi-omics data often contains NAs. To prevent errors in glmnet,
# we impute NAs with 0. This can be changed to median/mean if needed.
X_model[is.na(X_model)] <- 0

cat("Step 5: Training LASSO Logistic Regression Model...\n")
# Set a random seed to ensure reproducibility
set.seed(2026) 

# Run 5-fold cross-validated LASSO logistic regression.
# alpha = 1 specifies LASSO; family = "binomial" for 0/1 classification.
cv_lasso <- cv.glmnet(
  x = X_model,
  y = Y_model,
  family = "binomial",
  alpha = 1,
  nfolds = 5
)

# 畫出交叉驗證曲線 (這張圖強烈建議匯出放上 SMH 海報！)
plot(cv_lasso)
title("LASSO Cross-Validation Curve", line = 2.5)

cat("Step 6: Extracting Key Systemic Molecular Predictors...\n")
# 取得最佳的 L1 正規化懲罰值 (lambda)
best_lambda <- cv_lasso$lambda.min
cat("Optimal Lambda:", best_lambda, "\n")

# Extract features with non-zero coefficients at the optimal lambda
lasso_coefs <- coef(cv_lasso, s = "lambda.min")
non_zero_indices <- which(lasso_coefs != 0)

non_zero_features <- rownames(lasso_coefs)[non_zero_indices]
non_zero_values <- lasso_coefs[non_zero_indices]

# 整理成 Dataframe 並依權重絕對值排序
important_features <- data.frame(
  Feature = non_zero_features,
  Coefficient = non_zero_values
) %>% 
  filter(Feature != "(Intercept)") %>% 
  arrange(desc(abs(Coefficient)))

cat("\n=== Top Systemic Predictors for VO2 Max Response ===\n")
print(head(important_features, 15))


# Save the model and key features for subsequent biological interpretation
save(
  cv_lasso,
  important_features,
  file = "data/LASSO_Model_Results.RData"
)
cat("\nLASSO modeling complete. Results saved.\n")
