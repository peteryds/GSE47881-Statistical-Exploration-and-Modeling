# ==============================================================================
# MoTrPAC Rat Long-Term Training Project: Finding Plasma Metabolite Markers 
# for Exercise Performance at Week 8 (Concurrent Biomarkers)
# Ultimate Pooled Interaction Edition (% Change + N=10 LM + FDR)
# ==============================================================================

# Load required packages
if (!require("devtools", quietly = TRUE)) install.packages("devtools")
# devtools::install_github("MoTrPAC/MotrpacRatTraining6moData") # Uncomment if not installed
if (!require("broom", quietly = TRUE)) install.packages("broom")

library(MotrpacRatTraining6mo)
library(MotrpacRatTraining6moData)
library(tidyverse)
library(broom)

# ------------------------------------------------------------------------------
# Task 1: Organize basic data and calculate VO2 Max percent improvement (Dimensionality reduction to mitigate ceiling effect)
# ------------------------------------------------------------------------------
cat("\n[1/6] Extracting 8-week rat physical performance data and calculating percent improvement...\n")

vo2_8w_data <- PHENO %>%
  filter(group == "8w") %>%
  select(
    pid, 
    sex, 
    vo2_max_delta = `calculated.variables.vo2_max_change`,
    baseline_vo2 = `vo2.max.test.vo2_max_1` # Extract baseline physical performance
  ) %>%
  distinct() %>%                  
  drop_na(vo2_max_delta, baseline_vo2) %>%      
  mutate(
    # Core dimensionality reduction: Calculate percent improvement, blending baseline into the Y-axis
    vo2_percent_change = (vo2_max_delta / baseline_vo2) * 100,
    pid = as.character(pid)
  )

cat(sprintf("✅ Successfully extracted data for %d 8-week rats and converted to percent improvement.\n", nrow(vo2_8w_data)))

# ------------------------------------------------------------------------------
# Task 2: Identify plasma metabolites with significant changes at week 8 (8w) (approx. 68 candidates)
# ------------------------------------------------------------------------------
cat("\n[2/6] Searching for significant metabolites at week 8...\n")

plasma_metab_timewise <- metab_timewise_da(tissue = "PLASMA")

sig_8w_features <- plasma_metab_timewise %>%
  filter(comparison_group == "8w") %>% 
  group_by(comparison_group) %>%
  mutate(fdr = p.adjust(p_value, method = "fdr")) %>%
  ungroup() %>%
  filter(fdr < 0.05) %>%
  pull(feature_ID) %>%
  unique()

cat(sprintf("✅ Found %d metabolite features with significant changes at week 8.\n", length(sig_8w_features)))

# ------------------------------------------------------------------------------
# Task 3: Obtain the large expression matrix and convert to long format
# ------------------------------------------------------------------------------
cat("\n[3/6] Reading large expression matrix and organizing data...\n")

plasma_metab_expr <- load_sample_data(
  tissue = "PLASMA", 
  assay = "METAB", 
  normalized = TRUE, 
  training_regulated_only = FALSE
)

plasma_long <- plasma_metab_expr %>%
  pivot_longer(
    cols = -c(feature, feature_ID, tissue, assay, dataset), 
    names_to = "pid", 
    values_to = "intensity"
  ) %>%
  mutate(pid = as.character(pid))

cat("✅ Expression matrix organization complete.\n")

# ------------------------------------------------------------------------------
# Task 4: Execute N=10 Pooled Interaction Model
# ------------------------------------------------------------------------------
cat("\n[4/6] Executing N=10 mixed interaction regression model...\n")

# Merge data (mixed sex, N=10)
analysis_8w_pooled <- plasma_long %>%
  filter(feature_ID %in% sig_8w_features) %>%
  inner_join(vo2_8w_data, by = "pid") %>%
  drop_na(intensity, vo2_percent_change)

# --- Data health bulletproof check ---
sample_size_check <- analysis_8w_pooled %>%
  group_by(feature_ID) %>%
  summarize(valid_rats = n(), .groups = "drop")
dropped_features <- sample_size_check %>% filter(valid_rats < 8)
cat(sprintf("⚠️ Warning: Removed %d features with fewer than 8 valid samples.\n", nrow(dropped_features)))
# --------------------------

# Batch execute regression model for surviving features
lm_results <- analysis_8w_pooled %>%
  group_by(feature_ID) %>%
  filter(n() >= 8, sd(intensity) > 0, sd(vo2_percent_change) > 0) %>% 
  do({
    # 🌟 Core model: Percent improvement ~ expression * sex
    fit <- lm(vo2_percent_change ~ intensity * sex, data = .)
    tidy(fit)
  }) %>%
  ungroup()

# ------------------------------------------------------------------------------
# Task 5: Extract main effects and interactions, and perform FDR multiple testing correction
# ------------------------------------------------------------------------------
cat("\n[5/6] Extracting main effects and interactions, and performing FDR correction...\n")

# Dictionary translation preparation
data("METAB_FEATURE_ID_MAP", package = "MotrpacRatTraining6moData")
dict_clean <- METAB_FEATURE_ID_MAP %>%
  select(feature_ID = feature_ID_da, metabolite_name = metabolite_refmet) %>%
  distinct(feature_ID, .keep_all = TRUE)

# 1. Universal markers (Main Effect: intensity)
# Significance: Regardless of sex, whether metabolite intensity can predict physical improvement
main_effect_markers <- lm_results %>%
  filter(term == "intensity") %>%
  mutate(fdr_main = p.adjust(p.value, method = "BH")) %>%
  arrange(p.value) %>%
  left_join(dict_clean, by = "feature_ID") %>%
  mutate(metabolite_name = coalesce(metabolite_name, feature_ID)) %>%
  select(metabolite_name, estimate, p.value, fdr_main, feature_ID)
main_effect_markers

# 2. Sex-specific markers (Interaction Effect: intensity:sex)
# Significance: Whether the predictive power of the metabolite differs significantly between male and female rats
interaction_markers <- lm_results %>%
  filter(str_detect(term, ":sex")) %>%
  mutate(fdr_interaction = p.adjust(p.value, method = "BH")) %>%
  arrange(p.value) %>%
  left_join(dict_clean, by = "feature_ID") %>%
  mutate(metabolite_name = coalesce(metabolite_name, feature_ID)) %>%
  select(metabolite_name, estimate, p.value, fdr_interaction, feature_ID)
interaction_markers

cat("\n=======================================================\n")
cat("🏆 Universal Biomarkers (Main Effect Top 10) \n")
cat("=======================================================\n")
print(head(main_effect_markers, 10), width = Inf)

cat("\n=======================================================\n")
cat("🏆 Sex-specific Biomarkers (Interaction Effect Top 10) \n")
cat("=======================================================\n")
print(head(interaction_markers, 10), width = Inf)

# ------------------------------------------------------------------------------
# Task 6: Sex-specific marker visualization (Simpson's Paradox / Interaction Plot)
# ------------------------------------------------------------------------------
cat("\n[6/6] Plotting Week 8 sex-specific markers (Interaction Effect) chart...\n")

library(ggplot2)

# Automatically extract the top 4 feature IDs with the smallest Interaction P-values
top_4_interaction_ids <- interaction_markers$feature_ID[1:4]
top_4_interaction_names <- interaction_markers$metabolite_name[1:4]

plot_data_interaction <- analysis_8w_pooled %>%
  filter(feature_ID %in% top_4_interaction_ids) %>%
  # Apply human-readable labels
  mutate(Biomarker = factor(feature_ID, 
                            levels = top_4_interaction_ids, 
                            labels = top_4_interaction_names))

# 🚨 Foolproof checkpoint (Ensure data was successfully generated)
if(nrow(plot_data_interaction) == 0) {
  stop("❌ Error: Plotting data not found! Please ensure analysis_8w_pooled was generated correctly.")
}

fig_3 <- ggplot(plot_data_interaction, aes(x = intensity, y = vo2_percent_change)) +
  
  # Background line: Overall trend regardless of sex (Demonstrates misinterpretation caused by ignoring sex)
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, 
              color = "grey75", linetype = "dashed", linewidth = 1.2) +
  
  # Foreground: Data points
  geom_point(aes(color = sex), size = 4, alpha = 0.85) +
  
  # Highlight: Independent regression lines by sex (Demonstrates divergence / interaction)
  geom_smooth(aes(color = sex), method = "lm", se = FALSE, linewidth = 2) +
  
  # Facet into a 2x2 grid
  facet_wrap(~ Biomarker, scales = "free_x", ncol = 2) +
  
  # Poster-specific color scheme (Colorblind-friendly and high contrast)
  scale_color_manual(
    values = c("male" = "#0072B2", "female" = "#D55E00"),
    labels = c("male" = "Male (n=5)", "female" = "Female (n=5)")
  ) +
  
  # Title and axis labels (Unified format with W4, removed aggressive wording)
  labs(
    title = "Figure 2: Sexually Dimorphic Exercise Adaptation at Week 8",
    subtitle = "Top biomarkers exhibit strong Sex × Metabolite interactions.",
    x = "Log2 Normalized Plasma Intensity",
    y = "VO2 Max % Change from Baseline (%)",
    color = "Biological Sex"
  ) +
  
  # Poster-level theme settings (Enlarged fonts, removed noise)
  theme_minimal(base_size = 16) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey90"),
    strip.text = element_text(face = "bold", size = 14, hjust = 0),
    strip.background = element_rect(fill = "grey95", color = NA),
    plot.title = element_text(face = "bold", size = 20, color = "black"),
    plot.subtitle = element_text(size = 14, color = "grey40", margin = margin(b = 20)),
    legend.position = "top",
    legend.justification = "left",
    legend.text = element_text(size = 14, face = "bold"),
    legend.title = element_blank(), 
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )

print(fig_3)

# Export high-resolution PDF for poster use
ggsave("Figure_3_W8_Interaction.pdf", plot = fig_3, width = 12, height = 9, dpi = 300)
cat("✅ Figure 3 has been saved as Figure_3_W8_Interaction.pdf\n")