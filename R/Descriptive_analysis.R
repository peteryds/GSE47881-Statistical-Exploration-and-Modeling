#' Summarize Raw Gene Expression Data
#'
#' Produces summary statistics and QC plots (raw + log scale),
#' including PRE vs POST comparisons. These plots justify
#' log-transformation, normalization, and Empirical Bayes shrinkage.
#'
#' @param eset ExpressionSet object
#' @param output_dir Directory where QC plots will be saved
#'
#' @return A list containing global and group-wise summary statistics
#' @export

summarize_raw_data <- function(eset, output_dir = "output/raw_qc") {
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  message("🔍 Running Raw QC summarization...")
  
  # ============================================================
  # Extract Expression Data
  # ============================================================
  mat <- Biobase::exprs(eset)
  n_genes <- nrow(mat)
  
  time_vector <- ifelse(
    pData(eset)$`time:ch1` == "pre-training", "Pre", "Post"
  )
  
subject_ids <- pData(eset)$`patientid:ch1`
  
  df_all <- data.frame(
    expr     = as.numeric(mat),
    expr_log = log10(as.numeric(mat)),
    Time     = rep(time_vector, each = n_genes),
    Subject  = rep(subject_ids, each = n_genes)
  )
  
  df_pre  <- subset(df_all, Time == "Pre")
  df_post <- subset(df_all, Time == "Post")
  
  # Summary stats
  summary_stats <- summary(as.vector(mat))
  global_mean   <- mean(mat)
  global_sd     <- sd(as.vector(mat))
  pre_mean      <- mean(df_pre$expr)
  pre_sd        <- sd(df_pre$expr)
  post_mean     <- mean(df_post$expr)
  post_sd       <- sd(df_post$expr)
  
  # ============================================================
  # Identify paired subjects (subjects with BOTH PRE and POST)
  # ============================================================
  subjects_pre  <- unique(df_pre$Subject)
  subjects_post <- unique(df_post$Subject)
  paired_subjects <- intersect(subjects_pre, subjects_post)
  
  df_pre_paired  <- subset(df_pre,  Subject %in% paired_subjects)
  df_post_paired <- subset(df_post, Subject %in% paired_subjects)
  
  # ============================================================
  # Build long data for PRE–POST paired lines
  # ============================================================
  library(dplyr)
  library(tidyr)
  
  df_pairs_long <- full_join(
    df_pre_paired  %>% group_by(Subject) %>% summarise(expr_log_pre  = mean(expr_log)),
    df_post_paired %>% group_by(Subject) %>% summarise(expr_log_post = mean(expr_log)),
    by = "Subject"
  ) %>%
    pivot_longer(
      cols = c(expr_log_pre, expr_log_post),
      names_to = "Time",
      values_to = "Expr_log"
    ) %>%
    mutate(Time = ifelse(Time == "expr_log_pre", "Pre", "Post"),
           Time = factor(Time, levels = c("Pre","Post")))
  
  # ============================================================
  # 1. RAW BOXPLOT (PRE vs POST)
  # ============================================================
  p_box_raw <- ggplot(df_all, aes(Time, expr_log, fill = Time)) +
    geom_boxplot(outlier.alpha = 0.05) +
    scale_fill_manual(values = c(Pre="#00A6D6", Post="#FF6F61")) +
    labs(title = "Raw Intensities: PRE vs POST",
         y = "Raw expression") +
    theme_minimal(14)
  
  ggsave("box_raw_pre_post.png", p_box_raw, path = output_dir,
         width = 8, height = 6, dpi = 300)
  
  # ============================================================
  # 2. LOG VIOLIN + BOX + MEAN + CI 
  # ============================================================
  p_violin_log <- ggplot(df_all, aes(Time, expr_log, fill = Time)) +
    geom_violin(alpha = 0.3) +
    geom_boxplot(width = 0.2, outlier.alpha = 0.05) +
    stat_summary(fun = mean, geom = "point", size = 3, color = "black") +
    stat_summary(fun.data = mean_cl_normal, geom = "errorbar",
                 width = 0.15, color = "black") +
    scale_fill_manual(values = c(Pre="#00A6D6", Post="#FF6F61")) +
    labs(title = "Log10 Intensities: Violin + Boxplot + Mean + 95% CI",
         y = "Log10 expression") +
    theme_minimal(14)
  
  ggsave("box_log_violin_pre_post.png", p_violin_log, path = output_dir,
         width = 8, height = 6, dpi = 300)
  
  # ============================================================
  # 3A. PRE: RAW vs LOG HISTOGRAM (facet)
  # ============================================================
  df_pre_long <- rbind(
    data.frame(value=df_pre$expr,     Scale="Raw"),
    data.frame(value=df_pre$expr_log, Scale="Log10")
  )
  df_pre_long$Scale <- factor(df_pre_long$Scale, levels=c("Raw","Log10"))
  
  p_hist_pre <- ggplot(df_pre_long, aes(value)) +
    geom_histogram(bins=80, fill="#00A6D6", color="white") +
    facet_wrap(~Scale, ncol=2, scales="free_x") +
    labs(title="PRE Samples: RAW vs LOG", x="Expression", y="Frequency") +
    theme_minimal(14)
  
  ggsave("hist_pre_raw_vs_log.png", p_hist_pre,
         path = output_dir, width=12, height=6, dpi=300)
  
  # ============================================================
  # 3B. POST: RAW vs LOG HISTOGRAM (facet)
  # ============================================================
  df_post_long <- rbind(
    data.frame(value=df_post$expr,     Scale="Raw"),
    data.frame(value=df_post$expr_log, Scale="Log10")
  )
  df_post_long$Scale <- factor(df_post_long$Scale, levels=c("Raw","Log10"))
  
  p_hist_post <- ggplot(df_post_long, aes(value)) +
    geom_histogram(bins=80, fill="#FF6F61", color="white") +
    facet_wrap(~Scale, ncol=2, scales="free_x") +
    labs(title="POST Samples: RAW vs LOG", x="Expression", y="Frequency") +
    theme_minimal(14)
  
  ggsave("hist_post_raw_vs_log.png", p_hist_post,
         path = output_dir, width=12, height=6, dpi=300)
  
  # ============================================================
  # 4. LOG HISTOGRAM + NORMAL CURVE (PRE vs POST)
  # ============================================================
  mu_pre  <- mean(df_pre$expr_log)
  sd_pre  <- sd(df_pre$expr_log)
  mu_post <- mean(df_post$expr_log)
  sd_post <- sd(df_post$expr_log)
  
  x_pre  <- seq(min(df_pre$expr_log),  max(df_pre$expr_log),  length.out=400)
  x_post <- seq(min(df_post$expr_log), max(df_post$expr_log), length.out=400)
  
  curve_df <- rbind(
    data.frame(x=x_pre,  Time="Pre",  y=dnorm(x_pre,  mu_pre,  sd_pre)),
    data.frame(x=x_post, Time="Post", y=dnorm(x_post, mu_post, sd_post))
  )
  
  p_log_norm <- ggplot(df_all, aes(expr_log, fill=Time)) +
    geom_histogram(aes(y=after_stat(density)),
                   bins=80, color="white", alpha=0.8) +
    geom_line(data=curve_df, aes(x,y), color="black", linewidth=1) +
    facet_wrap(~Time, ncol=2, scales="free_y") +
    scale_fill_manual(values=c(Pre="#00A6D6", Post="#FF6F61")) +
    labs(title="Log10 Histogram + Normal Curve: PRE vs POST",
         x="Log10 expression", y="Density") +
    theme_minimal(14)
  
  ggsave("hist_log_pre_post_norm.png", p_log_norm,
         path = output_dir, width=12, height=6, dpi=300)
  
  # ============================================================
  # 5. RAW + LOG HISTOGRAM (ALL samples, orange) — FACET
  # ============================================================
  
  df_all_long <- rbind(
    data.frame(value = df_all$expr,     Scale = "Raw"),
    data.frame(value = df_all$expr_log, Scale = "Log10")
  )
  
  df_all_long$Scale <- factor(df_all_long$Scale, levels = c("Raw", "Log10"))
  
  p_hist_all_facet <- ggplot(df_all_long, aes(value)) +
    geom_histogram(bins = 80, fill = "orange", color = "white") +
    facet_wrap(~Scale, ncol = 2, scales = "free_x") +
    labs(title = "ALL Samples: RAW vs LOG10 Histograms",
         x = "Expression", y = "Frequency") +
    theme_minimal(14)
  
  ggsave("hist_all_raw_vs_log_orange_facet.png", 
         p_hist_all_facet,
         path = output_dir,
         width = 12, height = 6, dpi = 300)
  
  
  # ============================================================
  # RETURN SUMMARY
  # ============================================================
  return(list(
    global_mean = global_mean,
    global_sd   = global_sd,
    pre_mean    = pre_mean,
    post_mean   = post_mean,
    pre_sd      = pre_sd,
    post_sd     = post_sd,
    summary_stats = summary_stats
  ))
}
