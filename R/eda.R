#' Plot PCA (Principal Component Analysis)
#'
#' Performs PCA on the expression matrix to visualize sample structure.
#'
#' @param eset ExpressionSet object.
#' @param title Plot title.
#'
#' @return A ggplot object.
plot_pca <- function(eset, title = "PCA: All Genes") {
  require(ggplot2)
  
  expr_t <- t(exprs(eset))
  pca <- prcomp(expr_t, scale. = TRUE) 
  pca_df <- as.data.frame(pca$x)
  
  pca_df$geo_accession <- rownames(pca_df)
  pca_df$Timepoint <- pData(eset)$`time:ch1`
  pca_df$Subject <- pData(eset)$`patientid:ch1`
  
  var_explained <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 1)
  
  ggplot(pca_df, aes(x = PC1, y = PC2, color = Timepoint)) +
    geom_point(size = 3, alpha = 0.7) +
    labs(
      title = title,
      x = paste0("PC1 (", var_explained[1], "%)"),
      y = paste0("PC2 (", var_explained[2], "%)")
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
}

#' Plot Sample Distance Heatmap
#'
#' Calculates Euclidean distance between samples to check for outliers/batch effects.
#'
#' @param eset ExpressionSet object.
#'
#' @return A pheatmap object (drawn on device).
plot_sample_heatmap <- function(eset) {
  require(pheatmap)
  require(RColorBrewer)
  
  # Calculate distance matrix (transpose so we calculate distance between samples)
  sample_dists <- dist(t(exprs(eset)))
  sample_dist_matrix <- as.matrix(sample_dists)
  
  # Annotation for columns/rows
  annotation_df <- data.frame(
    Timepoint = pData(eset)$`time:ch1`
  )
  rownames(annotation_df) <- colnames(exprs(eset))
  
  colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
  
  pheatmap(sample_dist_matrix,
           clustering_distance_rows = sample_dists,
           clustering_distance_cols = sample_dists,
           col = colors,
           annotation_col = annotation_df,
           show_colnames = FALSE,
           main = "Sample-to-Sample Euclidean Distance")
}

#' Plot Density Distribution
#'
#' visualizes the distribution of expression values for all samples.
#' Good for checking normalization.
#'
#' @param eset ExpressionSet object.
#'
#' @return A ggplot object.
plot_density <- function(eset) {
  require(ggplot2)
  require(reshape2)
  
  # Melt expression matrix to long format
  expr_df <- melt(exprs(eset))
  colnames(expr_df) <- c("Gene", "Sample", "Expression")
  
  # Add timepoint info
  pdata <- pData(eset)
  expr_df$Timepoint <- pdata$`time:ch1`[match(expr_df$Sample, rownames(pdata))]
  
  ggplot(expr_df, aes(x = Expression, group = Sample, color = Timepoint)) +
    geom_density(alpha = 0.3, size = 0.5) +
    labs(title = "Density Plot of Expression Values",
         subtitle = "Check for normalization (curves should align)",
         x = "Log2 Expression",
         y = "Density") +
    theme_minimal() +
    theme(legend.position = "right")
}

#' Plot Volcano
#'
#' Visualizes statistical significance vs magnitude of change.
#'
#' @param limma_results Dataframe from topTable.
#' @param p_cutoff Significance threshold.
#' @param lfc_cutoff LogFC threshold.
#'
#' @return A ggplot object.
plot_volcano <- function(limma_results, p_cutoff = 0.05, lfc_cutoff = 0) {
  require(ggplot2)
  require(dplyr)
  
  df <- limma_results %>%
    mutate(
      Significance = case_when(
        adj.P.Val < p_cutoff & logFC > lfc_cutoff ~ "Up",
        adj.P.Val < p_cutoff & logFC < -lfc_cutoff ~ "Down",
        TRUE ~ "Not Sig"
      )
    )
  
  ggplot(df, aes(x = logFC, y = -log10(adj.P.Val), color = Significance)) +
    geom_point(alpha = 0.6, size = 1) +
    scale_color_manual(values = c("Up" = "red", "Down" = "blue", "Not Sig" = "grey")) +
    geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed") +
    labs(
      title = "Volcano Plot: Resistance Training Effect",
      x = "Log2 Fold Change (Post - Pre)",
      y = "-Log10 Adjusted P-value"
    ) +
    theme_minimal()
}

#' Plot Heatmap of Significant Genes
#'
#' Plots Z-score scaled expression for top genes.
#'
#' @param eset ExpressionSet object.
#' @param sig_genes Vector of gene IDs (probes).
#'
#' @return A pheatmap object.
plot_sig_heatmap <- function(eset, sig_genes) {
  require(pheatmap)
  
  if (length(sig_genes) < 2) return(NULL)
  
  # Subset data
  mat <- exprs(eset)[sig_genes, ]
  
  # Annotation
  annotation_df <- data.frame(
    Timepoint = pData(eset)$`time:ch1`
  )
  rownames(annotation_df) <- colnames(mat)
  
  # Scale by row (gene) to see patterns (Z-score)
  pheatmap(mat,
           scale = "row",
           annotation_col = annotation_df,
           show_rownames = FALSE,
           show_colnames = FALSE,
           main = paste("Heatmap of Top", length(sig_genes), "Significant Genes"),
           color = colorRampPalette(c("navy", "white", "firebrick3"))(100))
}

#' Plot Gene Regression (Age vs Change)
#'
#' @param df Wide dataframe with `_diff` columns.
#' @param gene_id Probe ID.
#'
#' @return A ggplot object.
plot_gene_age_scatter <- function(df, gene_id) {
  col_name <- paste0(gene_id, "_diff")
  if (!col_name %in% colnames(df)) col_name <- paste0("X", gene_id, "_diff")
  if (!col_name %in% colnames(df)) return(NULL)
  
  plot_data <- data.frame(Age = df$age, Change = df[[col_name]])
  
  ggplot(plot_data, aes(x = Age, y = Change)) +
    geom_point(size = 3, color = "darkgreen", alpha = 0.7) +
    geom_smooth(method = "lm", color = "black", linetype = "dashed", se = TRUE) +
    labs(
      title = paste("Age vs Change:", gene_id),
      x = "Age (years)", y = "Log2 Change (Post - Pre)"
    ) +
    theme_light()
}

#' Plot Gene Violin (Pre vs Post)
#'
#' Visualizes expression distribution before and after training for a single gene.
#'
#' @param eset ExpressionSet object.
#' @param gene_id Probe ID.
#'
#' @return A ggplot object.
plot_gene_violin <- function(eset, gene_id) {
  require(ggplot2)
  
  # Extract data for single gene
  expr_vals <- exprs(eset)[gene_id, ]
  plot_df <- data.frame(
    Expression = expr_vals,
    Timepoint = pData(eset)$`time:ch1`,
    Subject = pData(eset)$`patientid:ch1`
  )
  
  # Order timepoints
  plot_df$Timepoint <- factor(plot_df$Timepoint, levels = c("pre-training", "post-training"))
  
  ggplot(plot_df, aes(x = Timepoint, y = Expression, fill = Timepoint)) +
    geom_violin(trim = FALSE, alpha = 0.5) +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    geom_point(position = position_jitter(width = 0.1), size = 1, alpha = 0.6) +
    # Optional: connect dots for paired samples (looks messy with 44 lines, but useful)
    geom_line(aes(group = Subject), alpha = 0.2, color = "grey50") +
    scale_fill_manual(values = c("pre-training" = "#A4A4A4", "post-training" = "#D55E00")) +
    labs(
      title = paste("Expression Change:", gene_id),
      subtitle = "Lines connect the same subject (Pre -> Post)",
      y = "Log2 Expression Level"
    ) +
    theme_minimal() +
    theme(legend.position = "none")
}