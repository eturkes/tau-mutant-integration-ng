# Shared ggplot2 plotting layer: a project base theme, restrained colour/fill scales,
# and the cross-modality concordance scatter. Helpers consumed by P1+ analysis chapters.
# All ggplot2/ggrepel/stats calls are namespace-qualified so the file sources cleanly into
# any session. Plot data colours use a quiet, colourblind-aware journal palette; HTML
# chrome lives in theme.scss and stays decoupled.

tau_discrete_colours <- c(
  "#56616D",  # graphite
  "#3F6F6A",  # deep teal
  "#9B7A3C",  # muted ochre
  "#7A4052",  # wine
  "#647C8A",  # blue-grey
  "#8A6F83",  # muted mauve
  "#7C745F",  # umber-grey
  "#8A8A84"   # warm grey
)
tau_discrete_scale_types <- lapply(seq_along(tau_discrete_colours), function(i) {
  tau_discrete_colours[seq_len(i)]
})

genotype_colours <- c(
  MAPTKI      = "#56616D",
  P301S       = "#3F6F6A",
  NLGF_MAPTKI = "#9B7A3C",
  NLGF_P301S  = "#7A4052"
)

substate_colours <- c(
  Homeostatic   = "#5E7483",
  DAM           = "#7A4052",
  IFN           = "#9B7A3C",
  Proliferative = "#8A6F83"
)

tau_binary_colours <- c(
  `FALSE` = "#8A8A84",
  `TRUE`  = "#315E6F"
)

tau_direction_colours <- c(
  down = "#5E7483",
  up   = "#7A4052"
)

tau_background_colours <- c(
  MAPTKI = "#56616D",
  P301S  = "#3F6F6A"
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
      panel.grid.major = ggplot2::element_line(colour = "#E6E3DD", linewidth = 0.25),
      axis.line        = ggplot2::element_line(colour = "#3B3B37", linewidth = 0.25),
      plot.title       = ggplot2::element_text(face = "bold", colour = "#262B2F"),
      plot.subtitle    = ggplot2::element_text(colour = "#53504A"),
      strip.background = ggplot2::element_rect(fill = "#F2F0EA", colour = NA),
      strip.text       = ggplot2::element_text(face = "bold", colour = "#262B2F"),
      axis.title       = ggplot2::element_text(colour = "#333333"),
      axis.text        = ggplot2::element_text(colour = "#333333"),
      legend.title     = ggplot2::element_text(colour = "#333333")
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

scale_colour_tau_background <- function(...) {
  ggplot2::scale_colour_manual(values = tau_background_colours,
                               limits = names(tau_background_colours),
                               breaks = names(tau_background_colours), drop = FALSE, ...)
}
scale_color_tau_background <- scale_colour_tau_background

# Microglia-substate colour/fill scales: keep coherent biology colours stable wherever
# substates are mapped. `breaks=` can be narrowed to present_sub in a qmd.
scale_colour_substate <- function(..., breaks = microglia_substate_levels) {
  ggplot2::scale_colour_manual(values = substate_colours, limits = microglia_substate_levels,
                               breaks = breaks, drop = FALSE, ...)
}
scale_fill_substate <- function(..., breaks = microglia_substate_levels) {
  ggplot2::scale_fill_manual(values = substate_colours, limits = microglia_substate_levels,
                             breaks = breaks, drop = FALSE, ...)
}
scale_color_substate <- scale_colour_substate

# Binary and signed-direction helpers keep TRUE/FALSE and up/down states out of ggplot's
# high-saturation defaults.
scale_colour_tau_binary <- function(..., values = tau_binary_colours) {
  ggplot2::scale_colour_manual(values = values, breaks = names(values), drop = FALSE, ...)
}
scale_fill_tau_binary <- function(..., values = tau_binary_colours) {
  ggplot2::scale_fill_manual(values = values, breaks = names(values), drop = FALSE, ...)
}
scale_color_tau_binary <- scale_colour_tau_binary

scale_fill_direction <- function(..., values = tau_direction_colours) {
  ggplot2::scale_fill_manual(values = values, breaks = names(values), drop = FALSE, ...)
}

# Diverging/sequential heatmap scales. Signed effects pass `midpoint=0` so zero is the
# paper tone; count-density panels pass midpoint=NULL and get a neutral sequential scale.
rwb_colours <- c(low = "#4F6D7A", mid = "#F7F5F0", high = "#7A4052")
tau_sequential_colours <- c(low = "#F0EEE8", high = "#3F5F6F")
scale_fill_rwb <- function(..., midpoint = NULL, colours = rwb_colours) {
  stopifnot(is.character(colours), length(colours) == 3L)
  if (is.null(midpoint)) {
    ggplot2::scale_fill_gradient(low = tau_sequential_colours[["low"]],
                                 high = tau_sequential_colours[["high"]], ...)
  } else {
    ggplot2::scale_fill_gradient2(low = colours[[1]], mid = colours[[2]],
                                  high = colours[[3]], midpoint = midpoint, ...)
  }
}
scale_colour_rwb <- function(..., midpoint = NULL, colours = rwb_colours) {
  stopifnot(is.character(colours), length(colours) == 3L)
  if (is.null(midpoint)) {
    ggplot2::scale_colour_gradient(low = tau_sequential_colours[["low"]],
                                   high = tau_sequential_colours[["high"]], ...)
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
