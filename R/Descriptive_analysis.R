#' Summarize Raw Gene Expression Data
#'
#' Produces summary stats + QC plots for Empirical Bayes justification.
#'
#' @param eset ExpressionSet
#' @param output_dir Directory for plots
#' @return List with mean, sd, summary stats
summarize_raw_data <- function(eset, output_dir = "output/raw_qc") {
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  message("🔍 Summarizing raw data...")
  
  # Extract expression matrix
  mat <- Biobase::exprs(eset)
  
  # Summary
  summary_stats <- summary(as.vector(mat))
  global_mean <- mean(mat)
  global_sd   <- sd(as.vector(mat))
  
  message(paste("Global Mean:", round(global_mean, 3)))
  message(paste("Global SD:",   round(global_sd, 3)))
  
  # --------------------------------------------------
  # 1. Boxplot per sample
  # --------------------------------------------------
  png(file.path(output_dir, "boxplot_raw.png"), width = 1200, height = 800)
  boxplot(mat,
          main = "Raw Expression Values: Boxplot per Sample",
          xlab = "Sample",
          ylab = "Log Expression",
          outline = FALSE,
          col = "lightblue")
  dev.off()
  
  # --------------------------------------------------
  # 2. Histogram (log10 scale)
  # --------------------------------------------------
  df_plot <- data.frame(expr = as.vector(mat))
  
  p_hist <- ggplot(df_plot, aes(x = expr)) +
    geom_histogram(bins = 80, fill = "orange", color = "white") +
    scale_x_log10() +
    labs(
      title = "Histogram of Raw Expression Values (Log10 Scale)",
      x = "Expression (log10)",
      y = "Frequency"
    ) +
    theme_minimal(base_size = 14)
  
  ggsave("histogram_raw_log10.png", plot = p_hist, path = output_dir,
         width = 10, height = 6, dpi = 300)
  
  # --------------------------------------------------
  # 3. Density
  # --------------------------------------------------
  png(file.path(output_dir, "density_raw.png"), width = 1200, height = 800)
  plot(density(as.vector(mat)),
       main = "Density of Raw Expression Values",
       xlab = "Log Expression",
       lwd = 3, col = "darkblue")
  dev.off()
  
  # --------------------------------------------------
  # 4. Sample mean barplot
  # --------------------------------------------------
  sample_means <- colMeans(mat)
  png(file.path(output_dir, "sample_mean_barplot.png"), width = 1200, height = 800)
  barplot(sample_means,
          main = "Sample-Wise Mean Expression",
          xlab = "Sample",
          ylab = "Mean Log Expression",
          col = "orange")
  dev.off()
  
  list(
    global_mean = global_mean,
    global_sd = global_sd,
    summary_stats = summary_stats
  )
}


