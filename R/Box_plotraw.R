#' Global PRE vs POST Boxplot (Raw Intensities)
#'
#' @param eset ExpressionSet object
#' @param output_dir Directory to save plot
#'
#' @return Path to saved PNG plot
plot_raw_pre_post_boxplot <- function(eset, output_dir = "output/raw_qc") {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  mat <- Biobase::exprs(eset)
  
  # Extract the PRE/POST labels from phenotype
  time_vector <- pData(eset)$`time:ch1`
  time_vector <- ifelse(time_vector == "pre-training", "Pre", "Post")
  
  # Build long dataframe
  df <- data.frame(
    value = as.numeric(mat),
    sample = rep(colnames(mat), each = nrow(mat)),
    Time = rep(time_vector, each = nrow(mat))
  )
  
  # Boxplot
  p_box <- ggplot(df, aes(x = Time, y = value, fill = Time)) +
    geom_boxplot(outlier.alpha = 0.05) +
    scale_y_log10() +
    scale_fill_manual(values = c("Pre" = "#00A6D6", "Post" = "#FF6F61")) +
    labs(
      title = "PRE vs POST Raw Expression (All 54,675 Genes)",
      x = "",
      y = "Raw Log-Intensity (log10)"
    ) +
    theme_minimal(base_size = 14)
  
  out_path <- file.path(output_dir, "raw_pre_post_boxplot.png")
  
  ggsave(
    filename = out_path,
    plot = p_box,
    width = 10,
    height = 6,
    dpi = 300
  )
  
  message("Saved PRE vs POST raw boxplot → ", out_path)
  return(out_path)
}

