# ==============================================================================
# MoTrPAC Rat Long-Term Training Project: Finding Plasma Metabolite Markers for Exercise Performance at Week 8 (Concurrent Biomarkers)
# Final Complete Debugged Version
# ==============================================================================

# Load necessary packages
if (!require("devtools", quietly = TRUE)) install.packages("devtools")
# devtools::install_github("MoTrPAC/MotrpacRatTraining6moData") # Uncomment if not installed
library(MotrpacRatTraining6mo)
library(MotrpacRatTraining6moData)
library(tidyverse)

# ------------------------------------------------------------------------------
# Task 1: Organize basic data and VO2 max changes for the 8-week training group rats
# ------------------------------------------------------------------------------
cat("\n[1/5] Extracting 8-week rat physical performance data...\n")

PHENO$study_group_timepoint

column_names <- colnames(PHENO)
column_names

# data_8w <- PHENO %>%s
#   filter(group == "8w")  %>%
#   select(pid, sex, vo2_max_delta = `calculated.variables.vo2_max_change`, group,
#          vo2.max.test.vo2_max_1, vo2.max.test.vo2_max_2,
#          vo2.max.test.days_vo2_1, vo2.max.test.days_vo2_2, vo2.max.test.t_complete_1, 
#          vo2.max.test.t_complete_2, vo2.max.test.speed_max_1, vo2.max.test.blactate_begin_2) %>%
#   distinct()

vo2_8w_data <- PHENO %>%
  filter(group == "8w") %>%
  select(
    pid, 
    sex, 
    vo2_max_delta = `calculated.variables.vo2_max_change`
  ) %>%
  distinct() %>%                  # Remove duplicate sample rows, revert to individual level
  drop_na(vo2_max_delta) %>%      # Remove individuals missing physical performance data
  mutate(pid = as.character(pid)) # Ensure pid is a string type for easier merging later

vo2_8w_data

cat(sprintf("✅ Successfully extracted data for %d 8-week rats.\n", nrow(vo2_8w_data)))

# ------------------------------------------------------------------------------
# Task 2: Find plasma metabolites with significant changes at week 8 (8w)
# ------------------------------------------------------------------------------
cat("\n[2/5] Searching for significant metabolites at week 8...\n")

# Obtain plasma metabolome timewise differential analysis results
plasma_metab_timewise <- metab_timewise_da(tissue = "PLASMA")

column_names_timewise <- colnames(plasma_metab_timewise)
column_names_timewise

# Filter conditions: 8w and FDR < 0.05
sig_8w_features <- plasma_metab_timewise %>%
  filter(comparison_group == "8w") %>% 
  group_by(comparison_group) %>%
  mutate(fdr = p.adjust(p_value, method = "fdr")) %>%
  ungroup() %>%
  filter(fdr < 0.05) %>%
  pull(feature_ID) %>%
  unique()

sig_8w_features

cat(sprintf("✅ Found %d metabolite features with significant changes at week 8.\n", length(sig_8w_features)))

# ------------------------------------------------------------------------------
# Task 3: Get the large expression matrix and convert to long format
# ------------------------------------------------------------------------------
cat("\n[3/5] Reading the large expression matrix and organizing data...\n")

plasma_metab_expr <- load_sample_data(
  tissue = "PLASMA", 
  assay = "METAB", 
  normalized = TRUE, 
  training_regulated_only = FALSE
)
head(plasma_metab_expr)
colnames(plasma_metab_expr)
row.names(plasma_metab_expr)

# Wide to long format, ensure column names are directly pid
plasma_long <- plasma_metab_expr %>%
  pivot_longer(
    cols = -c(feature, feature_ID, tissue, assay, dataset), 
    names_to = "pid", 
    values_to = "intensity"
  ) %>%
  mutate(pid = as.character(pid))

head(plasma_long)
plasma_long

cat("✅ Expression matrix organization complete.\n")

# ------------------------------------------------------------------------------
# Task 4: Spearman correlation analysis by sex (Bulletproof revised version)
# ------------------------------------------------------------------------------
cat("\n[4/5] Calculating Spearman correlation between metabolites and VO2 max changes...\n")

# a <- plasma_long %>%
#   filter(feature_ID %in% sig_8w_features)

# Merge data (Perfectly joined using pid)
analysis_8w_correlation <- plasma_long %>%
  filter(feature_ID %in% sig_8w_features) %>%
  inner_join(vo2_8w_data, by = "pid")
analysis_8w_correlation


# sig_1w_features
# data_1w

# analysis_1w_correlation <- plasma_long %>%
#   filter(feature_ID %in% sig_1w_features) %>%
#   inner_join(data_1w, by = "pid")
# analysis_1w_correlation


# ------------------------------------------------------------------------------
# Task 4.5: Data health diagnosis (Inventory of missing values and sample sizes)
# ------------------------------------------------------------------------------
cat("\n📊 Conducting data health diagnosis...\n")

# 1. Count the total number of missing values (NA) in each column
cat("\n[Metric 1] Number of NAs in the merged data:\n")
na_counts <- colSums(is.na(analysis_8w_correlation))
print(na_counts)

# Find the 2 records where intensity expression is NA
missing_intensity_rows <- analysis_8w_correlation %>%
  filter(is.na(intensity))

cat("\n⚠️ The following are the data rows where intensity is NA:\n")
print(missing_intensity_rows)

# 2. Calculate the "total number of rats" and "valid number of rats" for each group by sex and metabolite
cat("\n[Metric 2] Sample size distribution for each group (Sex + Metabolite):\n")
sample_size_summary <- analysis_8w_correlation %>%
  group_by(sex, feature_ID) %>%
  summarize(
    total_rats = n(), # Total number of rats matched in this group
    valid_rats = sum(!is.na(intensity) & !is.na(vo2_max_delta)), # Number of rats truly available for correlation calculation after deducting NAs
    .groups = "drop"
  )

sample_size_summary

# Display a summary of the sample size distribution (e.g., see the min/max number of rats)
print(summary(sample_size_summary[, c("total_rats", "valid_rats")]))

# Check how many groups will be filtered out by your original bulletproof mechanism (n >= 4)
dropped_groups <- sample_size_summary %>% filter(valid_rats < 4)
cat(sprintf("⚠️ Warning: A total of %d groups will be excluded from subsequent calculations due to having less than 4 valid samples.\n", nrow(dropped_groups)))

cat("\n✅ Data diagnosis complete, preparing to enter correlation calculations.\n")

# ------------------------------------------------------------------------------

# Calculate correlations by sex
correlation_results <- analysis_8w_correlation %>%
  # Strictly filter missing values
  drop_na(intensity, vo2_max_delta) %>%
  group_by(sex, feature_ID) %>%
  # Bulletproof check: Ensure valid pairs >= 4, and both variables have variance
  filter(n() >= 4, sd(intensity) > 0, sd(vo2_max_delta) > 0) %>%
  summarize(
    cor_coef = cor(intensity, vo2_max_delta, method = "spearman"),
    p_val = cor.test(intensity, vo2_max_delta, method = "spearman", exact = FALSE)$p.value,
    .groups = "drop"
  ) %>%
  group_by(sex) %>%
  mutate(fdr_cor = p.adjust(p_val, method = "BH")) %>% 
  ungroup()

correlation_results

# Extract the top 10 for males and females separately (Fix sex strings to "male" and "female")
top_male_raw <- correlation_results %>% 
  filter(sex == "male") %>% 
  arrange(desc(abs(cor_coef))) %>% 
  head(10)

top_female_raw <- correlation_results %>% 
  filter(sex == "female") %>% 
  arrange(desc(abs(cor_coef))) %>% 
  head(10)

top_male_raw
top_female_raw

cat("✅ Correlation calculation and sorting complete.\n")

# ------------------------------------------------------------------------------
# Task 5: Use a dictionary to translate Feature IDs into human-readable chemical names (Fully joined version)
# ------------------------------------------------------------------------------
cat("\n[5/5] Translating metabolite names using the MoTrPAC dictionary...\n")

# Load the real metabolite dictionary from the package
data("METAB_FEATURE_ID_MAP", package = "MotrpacRatTraining6moData")

# Organize dictionary: Rename feature_ID_da to feature_ID for joining
dict_clean <- METAB_FEATURE_ID_MAP %>%
  select(feature_ID = feature_ID_da,          # This is the key to joining!
         metabolite_name = metabolite_refmet) %>%
  # Ensure only one name remains for the same feature to avoid explosion during merge
  distinct(feature_ID, .keep_all = TRUE)

# Translate male and female leaderboards
top_male_final <- top_male_raw %>%
  left_join(dict_clean, by = "feature_ID") %>%
  # If the name is not found in the dictionary, keep the original ID
  mutate(metabolite_name = coalesce(metabolite_name, feature_ID)) %>%
  select(sex, cor_coef, p_val, metabolite_name, feature_ID)

top_female_final <- top_female_raw %>%
  left_join(dict_clean, by = "feature_ID") %>%
  mutate(metabolite_name = coalesce(metabolite_name, feature_ID)) %>%
  select(sex, cor_coef, p_val, metabolite_name, feature_ID)

cat("\n=======================================================\n")
cat("🏆 Ultimate Version: Male Rat TOP 10 Metabolites \n")
cat("=======================================================\n")
print(top_male_final, width = Inf)

cat("\n=======================================================\n")
cat("🏆 Ultimate Version: Female Rat TOP 10 Metabolites \n")
cat("=======================================================\n")
print(top_female_final, width = Inf)



# ------------------------------------------------------------------------------
# Task 6: Simpson's Paradox Real Data Visualization
# ------------------------------------------------------------------------------
cat("\n[6/6] Plotting Simpson's Paradox chart in 'Storytelling with Data' style...\n")

library(ggplot2)

# 1. Select 4 representative metabolites with strong sex differences at 8w (Real existing IDs)
target_paradox_features <- c(
  "PC(35:3)_feature1",                               
  "Cystathionine 1",                             
  "kynurenine",                                  
  "LPC(18:1/0:0)_feature2"                       
)

# 2. Prepare plotting data: Filter target features and create readable labels
plot_data_paradox <- analysis_8w_correlation %>%
  filter(feature_ID %in% target_paradox_features) %>%
  # For poster communication effectiveness, attach biologically meaningful labels
  mutate(Biomarker = case_when(
    feature_ID == "PC(35:3)_feature1" ~ "PC(35:3) [Lipid Clearance]",
    feature_ID == "Cystathionine 1" ~ "Cystathionine [Antioxidant]",
    feature_ID == "kynurenine" ~ "Kynurenine [Neuro-Immune]",
    feature_ID == "LPC(18:1/0:0)_feature2" ~ "LPC(18:1/0:0) [Lyso-lipid]"
  )) %>%
  # Fix the arrangement order of the facet plots
  mutate(Biomarker = factor(Biomarker, levels = c(
    "PC(35:3) [Lipid Clearance]", "Cystathionine [Antioxidant]",
    "Kynurenine [Neuro-Immune]", "LPC(18:1/0:0) [Lyso-lipid]"
  )))

# 3. Draw high-communication chart based on Storytelling with Data (SWD) principles
p_simpson <- ggplot(plot_data_paradox, aes(x = intensity, y = vo2_max_delta)) +
  
  # [Background Layer: Misleading Illusion] Mixed-sex trend line (Thick dashed line, grey)
  # Use aes(group = 1) to force a single regression line for all points
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, 
              color = "grey70", linetype = "dashed", size = 1.2, alpha = 0.8) +
  
  # [Foreground Layer: Real Data] Data points colored by sex
  geom_point(aes(color = sex), size = 3, alpha = 0.8) +
  
  # [Focus Layer: Real Signal] Sex-stratified independent trend lines
  geom_smooth(aes(color = sex), method = "lm", se = FALSE, size = 1.5) +
  
  # Cut into 2x2 Small Multiples
  facet_wrap(~ Biomarker, scales = "free_x", ncol = 2) +
  
  # Set high contrast and colorblind-friendly colors (Male = Dark Blue, Female = Bright Orange)
  scale_color_manual(
    values = c("male" = "#0072B2", "female" = "#D55E00"),
    labels = c("male" = "Male (n=5)", "female" = "Female (n=5)")
  ) +
  
  # SWD Title Strategy: State the conclusion directly
  labs(
    title = "Simpson's Paradox in Exercise Adaptation (Week 8)",
    subtitle = "Combining sexes (grey dashed line) obscures the true, divergent physiological \nstrategies utilized by Males (blue) and Females (orange).",
    x = "Log2 Normalized Plasma Intensity",
    y = "Change in VO2 Max (ml/kg/min)",
    color = "Biological Sex"
  ) +
  
  # SWD Theme settings (Eliminate noise, emphasize contrast)
  theme_minimal(base_size = 14) +
  theme(
    # Remove unnecessary background grid lines
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey90"),
    
    # Keep facet titles clean and sharp
    strip.text = element_text(face = "bold", size = 12, hjust = 0),
    strip.background = element_blank(),
    
    # Strengthen title visual hierarchy
    plot.title = element_text(face = "bold", size = 18, color = "black"),
    plot.subtitle = element_text(size = 13, color = "grey30", margin = margin(b = 15)),
    
    # Move legend to the top left to align with reading flow
    legend.position = "top",
    legend.justification = "left",
    legend.text = element_text(size = 12, face = "bold"),
    legend.title = element_blank(), 
    
    # Add a slight border to hold the four plots together
    panel.border = element_rect(color = "black", fill = NA, size = 0.5)
  )

# Display the chart
print(p_simpson)

# Output as a high-resolution PDF for posters (Uncomment if needed)
ggsave("Poster_Simpsons_Paradox.pdf", p_simpson, width = 11, height = 8, dpi = 300)

cat("✅ Chart generation complete! Please check the Viewer window.\n")
