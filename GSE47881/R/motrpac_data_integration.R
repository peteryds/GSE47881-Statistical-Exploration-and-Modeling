# =========================================================================
# Script: R/01_motrpac_data_integration.R
# Project: Exploring Omic Responses to Exercise: An Initial Analysis
# Author: Yu-Hsin Yeh & Maria Teresa Lara Moran
# Objective: Phase 2 - MoTrPAC Multi-Omics Data Acquisition &
#            Pre-processing
# =========================================================================

# --- Step 1: Install and load required packages ---
# Use build_vignettes = FALSE to avoid installation errors
if (!requireNamespace("MotrpacRatTraining6mo", quietly = TRUE)) {
  devtools::install_github("MoTrPAC/MotrpacRatTraining6mo",
    build_vignettes = FALSE
  )
}
if (!requireNamespace("MotrpacRatTraining6moData", quietly = TRUE)) {
  devtools::install_github("MoTrPAC/MotrpacRatTraining6moData")
}

library(MotrpacRatTraining6mo)
library(MotrpacRatTraining6moData)
library(dplyr)
library(tidyr)
library(glmnet)

# --- Step 2: Define prediction target (Y: VO2 max High/Low Responders) ---
cat("Extracting Phenotypic Data and Calculating VO2 max Responses...\n")

# Extract data from PHENO, using official columns and deduplicating by "pid"
animal_target_data <- PHENO %>%
  filter(group == "8w") %>% # Filter for the 8-week training group
  mutate(
    # Use the official pre-calculated VO2 max percentage change column
    vo2_change_pct = calculated.variables.vo2_max_change,
    # Split into High/Low Responders based on the median
    responder_status = ifelse(vo2_change_pct >
      median(vo2_change_pct, na.rm = TRUE), "High", "Low"),
    # Convert to 0/1 required for LASSO
    responder_label = ifelse(responder_status == "High", 1, 0)
  ) %>%
  # Ensure one row per rat (pid) to avoid duplication from multiple tissues
  distinct(pid, .keep_all = TRUE) %>%
  select(pid, sex, vo2_change_pct, responder_status, responder_label)

# Create the Target Y vector with pid as the index names
target_y <- animal_target_data$responder_label
names(target_y) <- animal_target_data$pid

# Check the distribution of High/Low Responder rats
cat("VO2 Max Responders Distribution:\n")
print(table(animal_target_data$responder_status))

# --- Step 3: Integrate multi-tissue, multi-omics feature matrix (X) ---
cat("\nIntegrating Multi-Omics Feature Matrix...\n")
cat("This might take a few minutes depending on your RAM...\n")

target_tissues <- c("BLOOD", "SKM-GN", "LIVER", "HEART")
target_assays <- c("TRNSCRPT", "PROT")

# Extract all aligned features for the target tissues and assays
multi_omics_matrix <- combine_normalized_data(
  tissues = target_tissues,
  assays = target_assays
)

cat("\nDimensions of Initial Feature Matrix:\n")
print(dim(multi_omics_matrix))
head(multi_omics_matrix)

# --- Modified Step 4 & Step 5 ---
cat("\nData Acquired. Transposition and alignment handled in next script.\n")

# Create a 'data' directory and save the raw merged matrix and Y labels
if (!dir.exists("data")) dir.create("data")
save(
  multi_omics_matrix,
  target_y,
  animal_target_data,
  file = "data/Phase2_MoTrPAC_Cleaned.RData"
)

cat("\nSuccess! Saved to data/Phase2_MoTrPAC_Cleaned.RData\n")