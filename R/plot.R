# Plotting helpers. Currently just the concordance scatter used by the
# cross-modality integration block.

# Concordance scatter plot for two log2FC vectors with Spearman + Pearson rho.
concordance_plot <- function(df, x_col, y_col, label_col = "gene",
                             x_lab = NULL, y_lab = NULL, title = NULL,
                             top_n = 15) {
  df <- df[is.finite(df[[x_col]]) & is.finite(df[[y_col]]), , drop = FALSE]
  rho_s <- suppressWarnings(cor(df[[x_col]], df[[y_col]], method = "spearman"))
  rho_p <- suppressWarnings(cor(df[[x_col]], df[[y_col]], method = "pearson"))
  df$score <- abs(df[[x_col]]) + abs(df[[y_col]])
  top <- df |> dplyr::arrange(dplyr::desc(score)) |> head(top_n)
  ggplot(df, aes(.data[[x_col]], .data[[y_col]])) +
    geom_hline(yintercept = 0, colour = "grey70") +
    geom_vline(xintercept = 0, colour = "grey70") +
    geom_point(alpha = 0.3, size = 0.6) +
    geom_smooth(method = "lm", se = FALSE, colour = "#1f77b4") +
    ggrepel::geom_text_repel(data = top, aes(label = .data[[label_col]]),
                             size = 3, max.overlaps = 50) +
    labs(
      x = x_lab %||% x_col,
      y = y_lab %||% y_col,
      title = title,
      subtitle = sprintf("Spearman rho = %.3f, Pearson r = %.3f, n = %d",
                         rho_s, rho_p, nrow(df))
    ) +
    theme_bw()
}
