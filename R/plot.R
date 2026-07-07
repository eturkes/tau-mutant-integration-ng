# Shared ggplot2 plotting layer: a project base theme, restrained colour/fill scales,
# and the cross-modality concordance scatter. Helpers consumed by P1+ analysis chapters.
# All ggplot2/ggrepel/stats calls are namespace-qualified so the file sources cleanly into
# any session. Plot data colours use a quiet, colourblind-aware journal palette; HTML
# chrome lives in theme.scss and stays decoupled.

tau_discrete_colours <- c(
  "#3F5F7F",  # steel blue
  "#0B7A75",  # teal
  "#C8841C",  # amber
  "#A63A50",  # cranberry
  "#2F7EA8",  # azure
  "#7D5CB8",  # violet
  "#8A6A32",  # bronze
  "#6F7782"   # cool grey
)
tau_discrete_scale_types <- lapply(seq_along(tau_discrete_colours), function(i) {
  tau_discrete_colours[seq_len(i)]
})

genotype_colours <- c(
  MAPTKI      = "#3F5F7F",
  P301S       = "#0B7A75",
  NLGF_MAPTKI = "#C8841C",
  NLGF_P301S  = "#A63A50"
)

substate_colours <- c(
  Homeostatic   = "#2F78A0",
  DAM           = "#A63A50",
  IFN           = "#C8841C",
  Proliferative = "#7D5CB8"
)

tau_binary_colours <- c(
  `FALSE` = "#7C838A",
  `TRUE`  = "#0B6F7E"
)

tau_direction_colours <- c(
  down = "#2F78A0",
  up   = "#A63A50"
)

tau_background_colours <- c(
  MAPTKI = "#3F5F7F",
  P301S  = "#0B7A75"
)

figure7_score_fill_colours <- c(
  MAPTKI = "#245A9A",
  P301S  = "#C65A1E"
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

# Genotype colour/fill scales: map the canonical 4 genotypes to a stable saturated palette
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
rwb_colours <- c(low = "#2F78A0", mid = "#F8F5ED", high = "#A63A50")
tau_sequential_colours <- c(low = "#F1EEE5", high = "#1F6F8B")
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

# Amyloid-response interaction scatter (ONE modality panel) --------------------------------
# Per-feature amyloid effect on the tau-KO background (y = logFC NLGF_MAPTKI vs MAPTKI)
# against the mutant-tau background (x = logFC NLGF_P301S vs P301S). The dashed y=x identity
# line is the null of a tau-INDEPENDENT amyloid response; signed distance from it (y - x) is
# exactly the -interaction contrast, so features far off the diagonal are where mutant tau
# reshapes the amyloid response. Faint points, zero crosshairs, an OLS trend (tilt vs the
# diagonal = systematic interaction), and empirical off-diagonal outliers labelled. coord_equal
# on a symmetric square so the diagonal reads at 45 deg. `df` needs numeric x/y + a label column.
modality_interaction_scatter <- function(df, title = NULL, n_label = NULL,
                                         label_col = "label",
                                         label_tail_quantile = 0.99,
                                         x_lab = "log2FC  NLGF_P301S vs P301S",
                                         y_lab = "log2FC  NLGF_MAPTKI vs MAPTKI",
                                         point_colour = "#6F7782", label_colour = "#A63A50") {
  stopifnot(is.data.frame(df), all(c("x", "y", label_col) %in% names(df)))
  label_cutoff <- attr(df, "offdiag_cutoff", exact = TRUE)
  label_cutoff_source <- attr(df, "offdiag_cutoff_source", exact = TRUE)
  df    <- df[is.finite(df$x) & is.finite(df$y), , drop = FALSE]
  stopifnot(nrow(df) > 0L)
  rho_s <- suppressWarnings(stats::cor(df$x, df$y, method = "spearman"))
  rho_p <- suppressWarnings(stats::cor(df$x, df$y, method = "pearson"))
  lim   <- max(abs(c(df$x, df$y)), na.rm = TRUE)
  lim   <- if (is.finite(lim) && lim > 0) lim else 1
  # one label per feature: bulk assays reuse a site_id / gene across measured rows, so keep each label's
  # most-divergent threshold-passing instance (order is |y-x| desc) -> no duplicate repel labels.
  top   <- modality_scatter_label_rows(df, n_label = n_label, label_col = label_col,
                                       tail_quantile = label_tail_quantile,
                                       cutoff = label_cutoff,
                                       cutoff_source = label_cutoff_source)
  label_n <- nrow(top)
  if (label_n >= 80L) lim <- lim * 1.08
  label_size <- if (label_n >= 150L) {
    1.55
  } else if (label_n >= 80L) {
    1.85
  } else if (label_n >= 40L) {
    2.15
  } else {
    2.6
  }
  label_box_padding <- if (label_n >= 80L) 0.08 else 0.25
  label_point_padding <- if (label_n >= 80L) 0.03 else 0.10
  cutoff <- attr(top, "offdiag_cutoff")
  cutoff_source <- attr(top, "offdiag_cutoff_source") %||% ""
  cutoff_name <- if (grepl("^pooled", cutoff_source)) {
    "pooled cutoff"
  } else if (grepl("^within-method", cutoff_source) || grepl("^panel", cutoff_source)) {
    "within-method cutoff"
  } else {
    "cutoff"
  }
  label_subtitle <- if (is.finite(cutoff)) {
    sprintf("labels = %s, %s |x-y| >= %.2f",
            format(nrow(top), big.mark = ","), cutoff_name, cutoff)
  } else "labels = 0"
  ggplot2::ggplot(df, ggplot2::aes(x, y)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.25) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.25) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                         colour = "grey55", linewidth = 0.4) +
    ggplot2::geom_point(alpha = 0.25, size = 0.5, colour = point_colour) +
    # formula spelt out -> silence geom_smooth()'s default-formula message (keeps render logs clean)
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                         colour = "#2F7EA8", linewidth = 0.6) +
    ggplot2::geom_point(data = top, ggplot2::aes(x, y),
                        size = if (label_n >= 80L) 0.85 else 1.1,
                        colour = label_colour) +
    # max.overlaps = Inf + fixed seed -> deterministic layout, no "unlabeled points" warning (warn=2)
    ggrepel::geom_text_repel(data = top, ggplot2::aes(x, y, label = .data[[label_col]]),
                             size = label_size, colour = "#20242A", max.overlaps = Inf,
                             box.padding = label_box_padding,
                             point.padding = label_point_padding,
                             seed = 42L, min.segment.length = 0,
                             max.iter = 20000L, max.time = 3,
                             segment.colour = "grey65") +
    ggplot2::coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
    ggplot2::labs(
      x = x_lab, y = y_lab, title = title,
      subtitle = sprintf("Spearman = %.2f, Pearson = %.2f, n = %s\n%s",
                         rho_s, rho_p, format(nrow(df), big.mark = ","), label_subtitle)
    ) +
    theme_tau()
}

# Functional-category aggregate scores for empirical off-diagonal genes/proteins in the
# four-method amyloid-response scatter; phosphoproteomics uses the displayed parent-protein
# means. Rows are modality-specific role categories, facets are modalities. Each
# segment connects the aggregate amyloid logFC under MAPTKI to the aggregate amyloid logFC under
# P301S; segment colour is the requested contrast, P301S minus MAPTKI.
functional_group_score_plot <- function(group_summary, title = NULL) {
  stopifnot(is.data.frame(group_summary))
  need <- c("modality", "group_label_plot", "n_feature", "score_maptki", "score_p301s", "delta")
  miss <- setdiff(need, names(group_summary))
  if (length(miss)) {
    stop("group_summary missing columns: ", paste(miss, collapse = ", "), call. = FALSE)
  }
  x <- group_summary[group_summary$n_feature > 0L &
                       is.finite(group_summary$score_maptki) &
                       is.finite(group_summary$score_p301s) &
                       is.finite(group_summary$delta), , drop = FALSE]
  if (!nrow(x)) stop("group_summary has no finite aggregate scores", call. = FALSE)
  lim <- max(abs(c(x$score_maptki, x$score_p301s)), na.rm = TRUE)
  lim <- if (is.finite(lim) && lim > 0) lim else 1
  delta_lim <- max(abs(x$delta), na.rm = TRUE)
  delta_lim <- if (is.finite(delta_lim) && delta_lim > 0) delta_lim else 1
  n_min <- min(x$n_feature)
  n_max <- max(x$n_feature)
  size_breaks <- if (n_max - n_min <= 8L) {
    as.integer(unique(round(seq(n_min, n_max, length.out = min(3L, n_max - n_min + 1L)))))
  } else {
    as.integer(unique(round(c(n_min, pretty(c(n_min, n_max), n = 3), n_max))))
  }
  size_breaks <- size_breaks[size_breaks >= n_min & size_breaks <= n_max]
  if (!length(size_breaks)) size_breaks <- sort(unique(x$n_feature))
  point_df <- rbind(
    data.frame(x[, c("modality", "group_label_plot", "n_feature", "delta"), drop = FALSE],
               background = "MAPTKI", score = x$score_maptki, stringsAsFactors = FALSE),
    data.frame(x[, c("modality", "group_label_plot", "n_feature", "delta"), drop = FALSE],
               background = "P301S", score = x$score_p301s, stringsAsFactors = FALSE)
  )
  point_df$background <- factor(point_df$background, levels = c("MAPTKI", "P301S"))
  bg_fill <- figure7_score_fill_colours

  ggplot2::ggplot(x, ggplot2::aes(y = group_label_plot)) +
    ggplot2::geom_vline(xintercept = 0, colour = "#BFB8AA", linewidth = 0.3) +
    ggplot2::geom_segment(
      ggplot2::aes(x = score_maptki, xend = score_p301s, yend = group_label_plot,
                   colour = delta),
      linewidth = 0.75, lineend = "round") +
    ggplot2::geom_point(
      data = point_df,
      ggplot2::aes(x = score, fill = background, size = n_feature),
      shape = 21, colour = "#2B2A27", stroke = 0.25) +
    scale_colour_rwb(midpoint = 0, limits = c(-delta_lim, delta_lim), oob = scales::squish,
                     name = "P301S - MAPTKI") +
    ggplot2::scale_fill_manual(values = bg_fill, breaks = names(bg_fill),
                               labels = c(MAPTKI = "NLGF_MAPTKI", P301S = "NLGF_P301S"),
                               name = "score") +
    ggplot2::scale_size_area(max_size = 5.8, breaks = size_breaks, name = "scored items") +
    ggplot2::scale_x_continuous(limits = c(-lim, lim), oob = scales::squish) +
    ggplot2::scale_y_discrete(drop = TRUE) +
    ggplot2::facet_wrap(ggplot2::vars(modality), ncol = 2, scales = "free_y") +
    ggplot2::guides(
      colour = ggplot2::guide_colourbar(order = 1, barheight = grid::unit(0.45, "lines"),
                                        barwidth = grid::unit(4.5, "lines")),
      fill = ggplot2::guide_legend(order = 2,
                                   override.aes = list(size = 3.5, shape = 21,
                                                       colour = "#2B2A27")),
      size = ggplot2::guide_legend(order = 3)
    ) +
    ggplot2::labs(
      x = "Aggregate amyloid log2FC", y = NULL, title = title,
      subtitle = "Rows categorize within-method Q99 amyloid-effect scatter outliers; colour is P301S - MAPTKI"
    ) +
    theme_tau(base_size = 10) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 8, lineheight = 0.92),
      panel.grid.major.y = ggplot2::element_line(colour = "#ECE8DF", linewidth = 0.25),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.spacing.y = grid::unit(0.1, "lines")
    )
}

# Single-modality landscape plate ---------------------------------------------------------
# One figure per non-snRNAseq modality: left = the same per-feature amyloid-response
# interaction scatter used in the four-method panel, with a tighter default label cap for
# standalone readability; right = that modality's functional categories over the empirical
# off-diagonal tail. This keeps the modality-specific readout visual, target-derived, and
# free of duplicated prose.
modality_group_score_plot <- function(group_summary, modality, title = NULL, top_n = 6L) {
  stopifnot(is.data.frame(group_summary), is.character(modality), length(modality) == 1L,
            is.numeric(top_n), length(top_n) == 1L, top_n >= 1L)
  need <- c("modality", "group_label_plot", "n_feature", "score_maptki", "score_p301s",
            "delta")
  miss <- setdiff(need, names(group_summary))
  if (length(miss)) {
    stop("group_summary missing columns: ", paste(miss, collapse = ", "), call. = FALSE)
  }
  x <- group_summary[as.character(group_summary$modality) == modality &
                       group_summary$n_feature > 0L &
                       is.finite(group_summary$score_maptki) &
                       is.finite(group_summary$score_p301s) &
                       is.finite(group_summary$delta), , drop = FALSE]
  if (!nrow(x)) stop("group_summary has no finite rows for modality: ", modality, call. = FALSE)
  rank_score <- if ("rank_score" %in% names(x)) x$rank_score else abs(x$delta) * log1p(x$n_feature)
  group_priority <- if ("group_priority" %in% names(x)) x$group_priority else seq_len(nrow(x))
  group_label <- if ("group_label" %in% names(x)) as.character(x$group_label) else as.character(x$group_label_plot)
  x <- x[order(-rank_score, group_priority, group_label, method = "radix"), , drop = FALSE]
  x <- utils::head(x, as.integer(top_n))
  x$group_label_plot <- factor(as.character(x$group_label_plot),
                               levels = rev(as.character(x$group_label_plot)))

  lim <- max(abs(c(x$score_maptki, x$score_p301s)), na.rm = TRUE)
  lim <- if (is.finite(lim) && lim > 0) lim else 1
  delta_lim <- max(abs(x$delta), na.rm = TRUE)
  delta_lim <- if (is.finite(delta_lim) && delta_lim > 0) delta_lim else 1
  n_min <- min(x$n_feature)
  n_max <- max(x$n_feature)
  size_breaks <- as.integer(unique(round(seq(n_min, n_max, length.out = min(3L, n_max - n_min + 1L)))))
  size_breaks <- size_breaks[size_breaks >= n_min & size_breaks <= n_max]
  if (!length(size_breaks)) size_breaks <- sort(unique(x$n_feature))
  point_df <- rbind(
    data.frame(x[, c("group_label_plot", "n_feature", "delta"), drop = FALSE],
               background = "MAPTKI", score = x$score_maptki, stringsAsFactors = FALSE),
    data.frame(x[, c("group_label_plot", "n_feature", "delta"), drop = FALSE],
               background = "P301S", score = x$score_p301s, stringsAsFactors = FALSE)
  )
  point_df$background <- factor(point_df$background, levels = c("MAPTKI", "P301S"))

  ggplot2::ggplot(x, ggplot2::aes(y = group_label_plot)) +
    ggplot2::geom_vline(xintercept = 0, colour = "#BFB8AA", linewidth = 0.3) +
    ggplot2::geom_segment(
      ggplot2::aes(x = score_maptki, xend = score_p301s, yend = group_label_plot,
                   colour = delta),
      linewidth = 0.75, lineend = "round") +
    ggplot2::geom_point(
      data = point_df,
      ggplot2::aes(x = score, fill = background, size = n_feature),
      shape = 21, colour = "#2B2A27", stroke = 0.25) +
    scale_colour_rwb(midpoint = 0, limits = c(-delta_lim, delta_lim), oob = scales::squish,
                     name = "P301S - MAPTKI") +
    ggplot2::scale_fill_manual(values = figure7_score_fill_colours,
                               breaks = names(figure7_score_fill_colours),
                               labels = c(MAPTKI = "NLGF_MAPTKI", P301S = "NLGF_P301S"),
                               name = "score") +
    ggplot2::scale_size_area(max_size = 5.2, breaks = size_breaks, name = "scored items") +
    ggplot2::scale_x_continuous(limits = c(-lim, lim), oob = scales::squish) +
    ggplot2::labs(
      x = "Aggregate amyloid log2FC", y = NULL, title = title %||% "Off-diagonal categories",
      subtitle = paste0(modality, " Q99 scatter outliers; colour is P301S - MAPTKI")
    ) +
    theme_tau(base_size = 9) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 7.2, lineheight = 0.9),
      panel.grid.major.y = ggplot2::element_line(colour = "#ECE8DF", linewidth = 0.25),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.spacing.y = grid::unit(0.08, "lines")
    )
}

modality_landscape_plot <- function(modality_scatter_figures, modality, title = NULL,
                                    top_groups = 6L, scatter_n_label = 30L) {
  stopifnot(is.list(modality_scatter_figures), is.list(modality_scatter_figures$panels),
            is.list(modality_scatter_figures$groups),
            is.data.frame(modality_scatter_figures$groups$summary),
            is.character(modality), length(modality) == 1L)
  if (!modality %in% names(modality_scatter_figures$panels)) {
    stop("unknown modality panel: ", modality, call. = FALSE)
  }
  panel <- modality_scatter_figures$panels[[modality]]
  stopifnot(is.list(panel), is.data.frame(panel$data))
  scatter <- modality_interaction_scatter(
    panel$data,
    title = title %||% panel$title %||% modality,
    n_label = scatter_n_label
  ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 10.5),
      plot.subtitle = ggplot2::element_text(size = 8.2),
      axis.title = ggplot2::element_text(size = 8.5)
    )
  groups <- modality_group_score_plot(
    modality_scatter_figures$groups$summary,
    modality = modality,
    title = "Functional categories",
    top_n = top_groups
  ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 10.5),
      plot.subtitle = ggplot2::element_text(size = 8.2),
      axis.title = ggplot2::element_text(size = 8.5)
    )
  patchwork::wrap_plots(list(scatter, groups), ncol = 2, widths = c(1.05, 0.95),
                        guides = "collect")
}

# Cross-modality support matrix ------------------------------------------------------------
# One bubble per (feature x modality) cell of a cross-modality effect plate. Modality sits on
# the x-axis -- never colour -- so a reader scans a single feature row across assays for
# directional agreement, the standard multi-omics concordance read. Fill = signed effect on a
# diverging, symmetric scale (zero = paper); a dark ring flags a supported/earned call versus a
# faint ring for measured-but-unsupported; supported bubbles are drawn slightly larger. Faint
# x / open-triangle glyphs mark not-observed / blocked cells so coverage gaps stay explicit.
# `measured` and `missing` are the two row groups the axis-effect spine splits out
# (measured_state_chr == "measured" vs the rest). Facets are added by the caller.
plate_support_matrix <- function(measured, missing = NULL,
                                 effect_name = "signed effect",
                                 title = NULL, x_lab = NULL, y_lab = NULL,
                                 point_size = 3.7) {
  stopifnot(is.data.frame(measured), nrow(measured) > 0L,
            all(c("modality_label", "feature_label_plot", "plot_effect", "plot_status") %in%
                  names(measured)))
  measured <- measured[is.finite(measured$plot_effect), , drop = FALSE]
  stopifnot(nrow(measured) > 0L)
  measured$support_ring_lvl <- factor(
    ifelse(as.character(measured$plot_status) == "supported/earned", "supported", "measured"),
    levels = c("supported", "measured"))
  eff_max <- max(abs(measured$plot_effect), na.rm = TRUE)
  eff_max <- if (is.finite(eff_max) && eff_max > 0) eff_max else 1
  support_ring <- c(supported = "#20242A", measured = "#C7C1B4")

  p <- ggplot2::ggplot(measured,
                       ggplot2::aes(modality_label, feature_label_plot)) +
    ggplot2::geom_point(ggplot2::aes(fill = plot_effect, colour = support_ring_lvl,
                                     size = support_ring_lvl),
                        shape = 21, stroke = 0.7)
  # Draw the not-observed / blocked glyphs (and their shape scale) only when such cells exist --
  # a shape scale with no matching data trips ggplot's "no shared levels" warning (fatal at warn=2).
  has_missing <- FALSE
  if (!is.null(missing) && nrow(missing) > 0L) {
    missing <- missing[as.character(missing$plot_status) %in% c("not observed", "blocked"), ,
                       drop = FALSE]
    if (nrow(missing) > 0L) {
      has_missing <- TRUE
      missing$plot_status <- factor(as.character(missing$plot_status),
                                    levels = c("not observed", "blocked"))
      p <- p + ggplot2::geom_point(
        data = missing,
        ggplot2::aes(modality_label, feature_label_plot, shape = plot_status),
        inherit.aes = FALSE, colour = "#B4AFA3", size = 1.7, stroke = 0.5) +
        ggplot2::scale_shape_manual(values = c(`not observed` = 4, blocked = 2), name = NULL)
    }
  }
  p <- p +
    scale_fill_rwb(midpoint = 0, name = effect_name,
                   limits = c(-eff_max, eff_max), oob = scales::squish) +
    ggplot2::scale_colour_manual(values = support_ring, breaks = c("supported", "measured"),
                                 name = NULL) +
    ggplot2::scale_size_manual(values = c(supported = point_size,
                                          measured = point_size * 0.8), guide = "none") +
    ggplot2::labs(title = title, x = x_lab, y = y_lab) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(colour = "#EFEBE2", linewidth = 0.2),
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
      legend.position = "bottom")
  guide_args <- list(
    fill = ggplot2::guide_colourbar(order = 1, barheight = grid::unit(0.45, "lines"),
                                    barwidth = grid::unit(4.5, "lines")),
    colour = ggplot2::guide_legend(order = 2,
      override.aes = list(size = 3.2, fill = "#ECE7DC", shape = 21)))
  if (has_missing) {
    guide_args$shape <- ggplot2::guide_legend(order = 3,
      override.aes = list(colour = "#8C877B", size = 2.4))
  }
  p + do.call(ggplot2::guides, guide_args)
}

# Pair-support count matrix ----------------------------------------------------------------
# Compact matrix for the CCC-lite clearance pairs: one cell per (pair x contrast). Fill and the
# printed number both give how many modalities coherently support the pair; an "earned" pair
# gets a dark ring. Replaces the earlier per-slot dot strip, which read as cryptic. `pairs` needs
# columns x (contrast), y (pair), count (integer), earned (logical).
plate_pair_matrix <- function(pairs, title = NULL, count_name = "supported modalities") {
  stopifnot(is.data.frame(pairs), nrow(pairs) > 0L,
            all(c("x", "y", "count", "earned") %in% names(pairs)))
  pairs$count <- as.integer(round(pairs$count))
  pairs$earned <- factor(as.logical(pairs$earned), levels = c(TRUE, FALSE))
  count_max <- max(2L, max(pairs$count, na.rm = TRUE))
  ggplot2::ggplot(pairs, ggplot2::aes(x, y)) +
    ggplot2::geom_point(ggplot2::aes(fill = count, colour = earned),
                        shape = 21, size = 8, stroke = 0.9) +
    ggplot2::geom_text(ggplot2::aes(label = count), size = 2.9, colour = "#20242A") +
    ggplot2::scale_fill_gradient(low = "#F1EEE5", high = "#1F6F8B",
                                 limits = c(0, count_max), breaks = 0:count_max,
                                 name = count_name) +
    ggplot2::scale_colour_manual(values = c(`TRUE` = "#20242A", `FALSE` = "#D8D3C8"),
                                 breaks = c("TRUE", "FALSE"),
                                 labels = c(`TRUE` = "earned", `FALSE` = "not earned"),
                                 name = NULL) +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    ggplot2::guides(
      fill = ggplot2::guide_colourbar(order = 1, barheight = grid::unit(0.45, "lines"),
                                      barwidth = grid::unit(4, "lines")),
      colour = ggplot2::guide_legend(order = 2,
        override.aes = list(size = 3.5, fill = "#ECE7DC"))) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(colour = "#EFEBE2", linewidth = 0.2),
      axis.text.x = ggplot2::element_text(angle = 20, hjust = 1),
      legend.position = "bottom")
}
