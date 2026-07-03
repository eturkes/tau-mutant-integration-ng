# Shared ggplot2 plotting layer: a project base theme, muted colour/fill scales, and the
# cross-modality concordance scatter. Helpers consumed by P1+ analysis chapters. All
# ggplot2/ggrepel/stats calls are namespace-qualified so the file sources cleanly into any
# session (targets attaches only `quarto`; a Quarto chunk may attach ggplot2 for inline
# plotting -- these helpers need neither). Plot data colours use a deliberately plain,
# colourblind-aware project palette; the HTML chrome (IBM Plex + crimson) stays decoupled.

tau_discrete_colours <- c(
  "#4C78A8",  # blue
  "#F58518",  # orange
  "#54A24B",  # green
  "#C44E52",  # muted red
  "#72B7B2",  # teal
  "#B279A2",  # mauve
  "#9D755D",  # brown
  "#8C8C8C"   # grey
)
tau_discrete_scale_types <- lapply(seq_along(tau_discrete_colours), function(i) {
  tau_discrete_colours[seq_len(i)]
})

genotype_colours <- c(
  MAPTKI      = "#4C78A8",
  P301S       = "#54A24B",
  NLGF_MAPTKI = "#F58518",
  NLGF_P301S  = "#C44E52"
)

substate_colours <- c(
  Homeostatic   = "#4C78A8",
  DAM           = "#C44E52",
  IFN           = "#F58518",
  Proliferative = "#B279A2"
)

set_tau_plot_defaults <- function() {
  options(
    ggplot2.discrete.colour = tau_discrete_scale_types,
    ggplot2.discrete.fill   = tau_discrete_scale_types
  )
  invisible(tau_discrete_colours)
}

# Project base theme. base_family defaults to "" (device default sans) on purpose: a named
# family (e.g. "IBM Plex Sans") not registered with the graphics device warns at draw time and
# is non-portable across machines -> the Plex identity is carried by the HTML chrome (theme.scss)
# while the figures stay warning-free. Pass base_family explicitly once a phase registers the font.
theme_tau <- function(base_size = 11, base_family = "") {
  set_tau_plot_defaults()
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "grey92"),
      plot.title       = ggplot2::element_text(face = "bold"),
      plot.subtitle    = ggplot2::element_text(colour = "grey30"),
      strip.text       = ggplot2::element_text(face = "bold"),
      axis.title       = ggplot2::element_text(colour = "grey20")
    )
}

# Genotype colour/fill scales: map the canonical 4 genotypes to a stable muted palette
# while pinning the domain and legend order. limits/breaks = genotype_levels include all
# genotypes even when a subset is plotted (drop = FALSE alone does so only for a complete
# factor); `...` forwards name/labels/guide to the underlying manual scale. The data column
# must be a factor/character over genotype_levels.
scale_colour_genotype <- function(...) {
  ggplot2::scale_colour_manual(values = genotype_colours, limits = genotype_levels,
                               breaks = genotype_levels, drop = FALSE, ...)
}
scale_fill_genotype <- function(...) {
  ggplot2::scale_fill_manual(values = genotype_colours, limits = genotype_levels,
                             breaks = genotype_levels, drop = FALSE, ...)
}
scale_color_genotype <- scale_colour_genotype   # US-spelling alias

# Microglia-substate colour/fill scales: keep DAM red and IFN orange wherever
# coherent substates are mapped. `breaks=` can be narrowed to present_sub in a qmd.
scale_colour_substate <- function(..., breaks = microglia_substate_levels) {
  ggplot2::scale_colour_manual(values = substate_colours, limits = microglia_substate_levels,
                               breaks = breaks, drop = FALSE, ...)
}
scale_fill_substate <- function(..., breaks = microglia_substate_levels) {
  ggplot2::scale_fill_manual(values = substate_colours, limits = microglia_substate_levels,
                             breaks = breaks, drop = FALSE, ...)
}
scale_color_substate <- scale_colour_substate

# RWB heatmap scales: blue lows, white midpoint, red highs. `midpoint=NULL` maps
# the observed continuous range; signed effects pass `midpoint=0` so zero is white.
rwb_colours <- c(low = "#4C78A8", mid = "#F7F7F7", high = "#C44E52")
scale_fill_rwb <- function(..., midpoint = NULL, colours = rwb_colours) {
  stopifnot(is.character(colours), length(colours) == 3L)
  if (is.null(midpoint)) {
    ggplot2::scale_fill_gradientn(colours = unname(colours), ...)
  } else {
    ggplot2::scale_fill_gradient2(low = colours[[1]], mid = colours[[2]],
                                  high = colours[[3]], midpoint = midpoint, ...)
  }
}
scale_colour_rwb <- function(..., midpoint = NULL, colours = rwb_colours) {
  stopifnot(is.character(colours), length(colours) == 3L)
  if (is.null(midpoint)) {
    ggplot2::scale_colour_gradientn(colours = unname(colours), ...)
  } else {
    ggplot2::scale_colour_gradient2(low = colours[[1]], mid = colours[[2]],
                                    high = colours[[3]], midpoint = midpoint, ...)
  }
}
scale_color_rwb <- scale_colour_rwb

# Concordance scatter for two effect-size vectors (e.g. log2FC of the same features under two
# contrasts or modalities): faint points, an OLS trend line, zero crosshairs, and the top_n
# features by |x|+|y| labelled. Subtitle reports Spearman + Pearson correlation and n. `df` is
# filtered to rows finite in both columns first. Used by P4 cross-modality concordance.
concordance_plot <- function(df, x_col, y_col, label_col = "gene",
                             x_lab = NULL, y_lab = NULL, title = NULL, top_n = 15) {
  df    <- df[is.finite(df[[x_col]]) & is.finite(df[[y_col]]), , drop = FALSE]
  rho_s <- suppressWarnings(stats::cor(df[[x_col]], df[[y_col]], method = "spearman"))
  rho_p <- suppressWarnings(stats::cor(df[[x_col]], df[[y_col]], method = "pearson"))
  ord   <- order(abs(df[[x_col]]) + abs(df[[y_col]]), decreasing = TRUE)   # base order -> no temp col
  top   <- df[utils::head(ord, top_n), , drop = FALSE]
  ggplot2::ggplot(df, ggplot2::aes(.data[[x_col]], .data[[y_col]])) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey70") +
    ggplot2::geom_vline(xintercept = 0, colour = "grey70") +
    ggplot2::geom_point(alpha = 0.3, size = 0.6) +
    # formula spelt out -> silence geom_smooth()'s default-formula message (keeps render logs clean)
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
    ggrepel::geom_text_repel(data = top, ggplot2::aes(label = .data[[label_col]]),
                             size = 3, max.overlaps = 50) +
    ggplot2::labs(
      x = x_lab %||% x_col, y = y_lab %||% y_col, title = title,
      subtitle = sprintf("Spearman rho = %.3f, Pearson r = %.3f, n = %d",
                         rho_s, rho_p, nrow(df))
    ) +
    theme_tau()
}
