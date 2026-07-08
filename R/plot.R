# Shared ggplot2 plotting layer: a project base theme, restrained colour/fill scales,
# and helpers consumed by the rendered report.
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

# Modality-native descriptive plates -------------------------------------------------------
# GeoMx gets a spatial AOI score plate; bulk proteome gets sample PCA + protein volcano; phospho
# gets phosphosite volcano + top-site abundance heatmap. These are intentionally separate from
# the integrated amyloid-effect scatter below.
modality_volcano_plot <- function(volcano, title = NULL,
                                  x_lab = "log2FC  NLGF_P301S vs P301S",
                                  alpha = 0.10) {
  stopifnot(is.data.frame(volcano),
            all(c("effect", "neg_log10_fdr", "direction", "label", "label_show") %in%
                  names(volcano)),
            is.numeric(alpha), length(alpha) == 1L, alpha > 0, alpha < 1)
  x <- volcano[is.finite(volcano$effect) & is.finite(volcano$neg_log10_fdr), ,
               drop = FALSE]
  if (!nrow(x)) stop("volcano has no finite rows", call. = FALSE)
  labels <- x[x$label_show %in% TRUE, , drop = FALSE]
  lim <- max(abs(x$effect), na.rm = TRUE)
  lim <- if (is.finite(lim) && lim > 0) lim else 1
  y_cut <- -log10(alpha)
  direction_colours <- c(down = "#2F78A0", `not significant` = "#AFA89C", up = "#A63A50")

  ggplot2::ggplot(x, ggplot2::aes(effect, neg_log10_fdr)) +
    ggplot2::geom_hline(yintercept = y_cut, colour = "#BFB8AA", linewidth = 0.3,
                        linetype = "dotted") +
    ggplot2::geom_vline(xintercept = 0, colour = "#BFB8AA", linewidth = 0.3) +
    ggplot2::geom_point(ggplot2::aes(colour = direction), alpha = 0.45, size = 0.75) +
    ggrepel::geom_text_repel(
      data = labels,
      ggplot2::aes(label = label),
      size = 2.35, colour = "#20242A", max.overlaps = Inf, seed = 42L,
      min.segment.length = 0, segment.colour = "grey65", box.padding = 0.16,
      point.padding = 0.06, max.iter = 20000L, max.time = 3) +
    ggplot2::scale_colour_manual(values = direction_colours, drop = FALSE, name = NULL) +
    ggplot2::scale_x_continuous(limits = c(-lim, lim), oob = scales::squish) +
    ggplot2::labs(
      x = x_lab, y = "-log10 FDR", title = title,
      subtitle = sprintf("FDR threshold = %.2f; labelled points are top-ranked by FDR and effect",
                         alpha)
    ) +
    theme_tau(base_size = 9.5) +
    ggplot2::theme(legend.position = "bottom")
}

geomx_qc_atlas_plot <- function(qc, title = "GeoMx AOI QC atlas") {
  stopifnot(is.list(qc), is.data.frame(qc$metrics), is.data.frame(qc$flag_counts))
  metrics <- qc$metrics
  need <- c("slide", "segment", "genotype", "metric_label", "value", "flag_metric")
  miss <- setdiff(need, names(metrics))
  if (length(miss)) stop("GeoMx QC metrics missing columns: ", paste(miss, collapse = ", "),
                         call. = FALSE)
  metrics <- metrics[is.finite(metrics$value), , drop = FALSE]
  if (!nrow(metrics)) stop("GeoMx QC metrics have no finite values", call. = FALSE)
  metrics$slide <- factor(metrics$slide)
  metrics$segment <- factor(metrics$segment)
  metrics$genotype <- factor(as.character(metrics$genotype), levels = genotype_levels)
  stopifnot(!anyNA(metrics$genotype), !anyNA(metrics$slide), !anyNA(metrics$segment))

  segment_levels <- levels(droplevels(metrics$segment))
  shape_values <- stats::setNames(
    rep(c(16, 17, 15, 18, 3, 4, 7, 8, 0, 1, 2, 5, 6), length.out = length(segment_levels)),
    segment_levels
  )
  flagged <- metrics[metrics$flag_metric %in% TRUE, , drop = FALSE]

  metric_plot <- ggplot2::ggplot(metrics, ggplot2::aes(slide, value)) +
    ggplot2::geom_boxplot(
      ggplot2::aes(group = slide),
      width = 0.54, outlier.shape = NA, fill = "#F2F0EA", colour = "#8A8174",
      linewidth = 0.28) +
    ggplot2::geom_point(
      ggplot2::aes(colour = genotype, shape = segment),
      position = ggplot2::position_jitter(width = 0.14, height = 0, seed = 42L),
      size = 1.05, alpha = 0.68) +
    ggplot2::geom_point(
      data = flagged,
      position = ggplot2::position_jitter(width = 0.14, height = 0, seed = 42L),
      shape = 21, fill = NA, colour = "#20242A", stroke = 0.44, size = 1.9) +
    scale_colour_genotype(name = "genotype") +
    ggplot2::scale_shape_manual(values = shape_values, name = "segment", drop = FALSE) +
    ggplot2::facet_wrap(ggplot2::vars(metric_label), scales = "free_y", ncol = 3) +
    ggplot2::scale_y_continuous(labels = scales::label_number(big.mark = ",")) +
    ggplot2::labs(
      x = "slide", y = NULL, title = title,
      subtitle = "AOIs by slide; colour = genotype, shape = segment, black ring = metric-specific QC flag"
    ) +
    theme_tau(base_size = 8.8) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
      panel.spacing = grid::unit(0.65, "lines")
    )

  flags <- qc$flag_counts
  need_flag <- c("slide", "segment", "flag", "n")
  miss_flag <- setdiff(need_flag, names(flags))
  if (length(miss_flag)) {
    stop("GeoMx QC flag counts missing columns: ", paste(miss_flag, collapse = ", "),
         call. = FALSE)
  }
  flags <- flags[is.finite(flags$n), , drop = FALSE]
  if (!nrow(flags)) stop("GeoMx QC flag counts have no finite rows", call. = FALSE)
  label_flags <- flags[flags$n > 0, , drop = FALSE]
  max_n <- max(flags$n, na.rm = TRUE)
  max_n <- if (is.finite(max_n) && max_n > 0) max_n else 1

  flag_plot <- ggplot2::ggplot(flags, ggplot2::aes(slide, flag)) +
    ggplot2::geom_point(
      ggplot2::aes(size = n, fill = n),
      shape = 21, colour = "#554F47", stroke = 0.25, alpha = 0.92) +
    ggplot2::geom_text(
      data = label_flags,
      ggplot2::aes(label = n),
      size = 2.25, colour = "#20242A", fontface = "bold") +
    scale_fill_rwb(midpoint = NULL, limits = c(0, max_n), name = "flagged AOIs") +
    ggplot2::scale_size_area(max_size = 8, limits = c(0, max_n), guide = "none") +
    ggplot2::facet_wrap(ggplot2::vars(segment), nrow = 1) +
    ggplot2::labs(
      x = "slide", y = NULL,
      title = "QC flag counts",
      subtitle = qc$provenance$flag_rule %||% "descriptive QC flags only; no AOIs excluded"
    ) +
    theme_tau(base_size = 8.6) +
    ggplot2::theme(
      legend.position = "bottom",
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 7.4),
      panel.spacing.x = grid::unit(0.7, "lines")
    )

  patchwork::wrap_plots(list(metric_plot, flag_plot), ncol = 1, heights = c(1.75, 1))
}

geomx_spatial_modality_plot <- function(spatial, title = "GeoMx spatial AOIs") {
  stopifnot(is.list(spatial), is.data.frame(spatial$aoi), is.data.frame(spatial$genes))
  aoi <- spatial$aoi
  need <- c("slide", "roi", "genotype", "x_coord", "y_coord", "aoi_area",
            "signed_response_score", "score_abs")
  miss <- setdiff(need, names(aoi))
  if (length(miss)) stop("GeoMx spatial data missing columns: ", paste(miss, collapse = ", "),
                         call. = FALSE)
  aoi <- aoi[is.finite(aoi$x_coord) & is.finite(aoi$y_coord) &
               is.finite(aoi$aoi_area) & aoi$aoi_area > 0 &
               is.finite(aoi$signed_response_score) & is.finite(aoi$score_abs), ,
             drop = FALSE]
  if (!nrow(aoi)) stop("GeoMx spatial data has no finite AOIs", call. = FALSE)
  score_lim <- max(abs(aoi$signed_response_score), na.rm = TRUE)
  score_lim <- if (is.finite(score_lim) && score_lim > 0) score_lim else 1
  area_lab <- if (max(aoi$aoi_area, na.rm = TRUE) >= 1000) {
    function(z) paste0(format(round(z / 1000, 1), trim = TRUE), "k")
  } else {
    function(z) format(round(z, 1), trim = TRUE)
  }
  label_aoi <- aoi[0, , drop = FALSE]
  for (sl in levels(droplevels(aoi$slide))) {
    d <- aoi[aoi$slide == sl, , drop = FALSE]
    keep <- utils::head(order(-d$score_abs, d$roi, method = "radix"), 2L)
    label_aoi <- rbind(label_aoi, d[keep, , drop = FALSE])
  }

  spatial_map <- ggplot2::ggplot(aoi, ggplot2::aes(x_coord, y_coord)) +
    ggplot2::geom_point(
      ggplot2::aes(fill = signed_response_score, colour = genotype, size = aoi_area),
      shape = 21, stroke = 0.42, alpha = 0.95) +
    ggrepel::geom_text_repel(
      data = label_aoi,
      ggplot2::aes(label = roi),
      size = 2.25, colour = "#20242A", max.overlaps = Inf, seed = 42L,
      min.segment.length = 0, segment.colour = "grey65", box.padding = 0.12,
      point.padding = 0.08, max.iter = 20000L, max.time = 3) +
    scale_fill_rwb(midpoint = 0, limits = c(-score_lim, score_lim), oob = scales::squish,
                   name = "signed score") +
    scale_colour_genotype(name = "genotype") +
    ggplot2::scale_size_continuous(range = c(1.4, 5.0), labels = area_lab,
                                   name = "AOI area") +
    ggplot2::facet_wrap(ggplot2::vars(slide), ncol = 4) +
    ggplot2::coord_equal() +
    ggplot2::scale_y_reverse() +
    ggplot2::labs(x = "ROI coordinate X", y = "ROI coordinate Y",
                  title = title,
                  subtitle = "AOI coordinates by slide; fill = score, ring = genotype, size = AOI area") +
    theme_tau(base_size = 9) +
    ggplot2::guides(
      fill = ggplot2::guide_colourbar(order = 1),
      colour = ggplot2::guide_legend(order = 2, override.aes = list(fill = "white", size = 3)),
      size = ggplot2::guide_legend(order = 3)
    ) +
    ggplot2::theme(legend.position = "bottom", legend.box = "horizontal")

  score_dist <- ggplot2::ggplot(aoi, ggplot2::aes(genotype, signed_response_score)) +
    ggplot2::geom_hline(yintercept = 0, colour = "#BFB8AA", linewidth = 0.3) +
    ggplot2::geom_violin(ggplot2::aes(fill = genotype), width = 0.82, alpha = 0.26,
                         colour = NA, trim = FALSE) +
    ggplot2::geom_boxplot(ggplot2::aes(colour = genotype), width = 0.18,
                          fill = "white", outlier.shape = NA, linewidth = 0.35) +
    ggplot2::geom_point(ggplot2::aes(colour = genotype),
                        position = ggplot2::position_jitter(width = 0.10, height = 0, seed = 42L),
                        size = 1.15, alpha = 0.68) +
    scale_fill_genotype(guide = "none") +
    scale_colour_genotype(guide = "none") +
    ggplot2::labs(x = NULL, y = "signed score",
                  title = "AOI score by genotype",
                  subtitle = sprintf("n = %s AOIs; descriptive, not an independent-replicate test",
                                     format(nrow(aoi), big.mark = ","))) +
    theme_tau(base_size = 9) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))

  genes <- spatial$genes
  need_gene <- c("symbol", "y", "x", "mean_effect", "rank_score")
  miss_gene <- setdiff(need_gene, names(genes))
  if (length(miss_gene)) {
    stop("GeoMx score-gene data missing columns: ", paste(miss_gene, collapse = ", "),
         call. = FALSE)
  }
  genes <- genes[is.finite(genes$y) & is.finite(genes$x) &
                   is.finite(genes$mean_effect) & is.finite(genes$rank_score), ,
                 drop = FALSE]
  if (!nrow(genes)) stop("GeoMx score-gene data has no finite rows", call. = FALSE)
  genes <- utils::head(genes[order(-genes$rank_score, -abs(genes$mean_effect),
                                  genes$symbol, method = "radix"), , drop = FALSE], 12L)
  genes$symbol_plot <- factor(genes$symbol, levels = rev(genes$symbol))
  gene_long <- rbind(
    data.frame(symbol_plot = genes$symbol_plot, background = "MAPTKI",
               effect = genes$y, stringsAsFactors = FALSE),
    data.frame(symbol_plot = genes$symbol_plot, background = "P301S",
               effect = genes$x, stringsAsFactors = FALSE)
  )
  gene_long$background <- factor(gene_long$background, levels = c("MAPTKI", "P301S"))
  driver_plot <- ggplot2::ggplot(genes, ggplot2::aes(y = symbol_plot)) +
    ggplot2::geom_vline(xintercept = 0, colour = "#BFB8AA", linewidth = 0.3) +
    ggplot2::geom_segment(ggplot2::aes(x = y, xend = x, yend = symbol_plot),
                          linewidth = 0.55, colour = "#8A8174", alpha = 0.8) +
    ggplot2::geom_point(data = gene_long, ggplot2::aes(x = effect, colour = background),
                        size = 1.9, alpha = 0.95) +
    scale_colour_tau_background(name = "amyloid effect in") +
    ggplot2::labs(x = "GeoMx log2FC", y = NULL,
                  title = "Score-gene drivers",
                  subtitle = "Top genes ranked by amyloid-effect size and FDR") +
    theme_tau(base_size = 9) +
    ggplot2::theme(legend.position = "bottom",
                   axis.text.y = ggplot2::element_text(size = 7.4))

  bottom <- patchwork::wrap_plots(list(score_dist, driver_plot), ncol = 2,
                                  widths = c(0.9, 1.1))
  patchwork::wrap_plots(list(spatial_map, bottom), ncol = 1,
                        heights = c(1.12, 0.88))
}

proteome_modality_plot <- function(proteome, title = "Bulk proteome") {
  stopifnot(is.list(proteome), is.data.frame(proteome$pca), is.data.frame(proteome$volcano))
  pca <- proteome$pca
  need <- c("pc1", "pc2", "pc1_var", "pc2_var", "genotype", "run_index")
  miss <- setdiff(need, names(pca))
  if (length(miss)) stop("proteome PCA data missing columns: ", paste(miss, collapse = ", "),
                         call. = FALSE)
  pca_plot <- ggplot2::ggplot(pca, ggplot2::aes(pc1, pc2, colour = genotype)) +
    ggplot2::geom_hline(yintercept = 0, colour = "#D5D0C4", linewidth = 0.25) +
    ggplot2::geom_vline(xintercept = 0, colour = "#D5D0C4", linewidth = 0.25) +
    ggplot2::geom_point(size = 2.4, alpha = 0.95) +
    ggrepel::geom_text_repel(ggplot2::aes(label = run_index), size = 2.3,
                             max.overlaps = Inf, seed = 42L, min.segment.length = 0,
                             segment.colour = "grey65", box.padding = 0.12,
                             point.padding = 0.08, max.iter = 20000L, max.time = 3) +
    scale_colour_genotype(name = "genotype") +
    ggplot2::labs(
      x = sprintf("PC1 (%.1f%%)", 100 * unique(pca$pc1_var)[[1]]),
      y = sprintf("PC2 (%.1f%%)", 100 * unique(pca$pc2_var)[[1]]),
      title = title, subtitle = "Sample PCA of protein-group intensities") +
    theme_tau(base_size = 9.5) +
    ggplot2::theme(legend.position = "bottom")

  volcano <- modality_volcano_plot(
    proteome$volcano,
    title = "Protein differential abundance",
    alpha = proteome$provenance$alpha %||% 0.10
  )
  patchwork::wrap_plots(list(pca_plot, volcano), ncol = 2, widths = c(0.9, 1.1))
}

phospho_site_heatmap_plot <- function(heatmap, title = "Top phosphosite abundance") {
  stopifnot(is.data.frame(heatmap),
            all(c("sample_label", "site_label_plot", "genotype", "z") %in% names(heatmap)))
  x <- heatmap[is.finite(heatmap$z), , drop = FALSE]
  if (!nrow(x)) stop("phosphosite heatmap has no finite rows", call. = FALSE)
  lim <- max(abs(x$z), na.rm = TRUE)
  lim <- if (is.finite(lim) && lim > 0) lim else 1

  ggplot2::ggplot(x, ggplot2::aes(sample_label, site_label_plot, fill = z)) +
    ggplot2::geom_tile(colour = "#F8F5ED", linewidth = 0.18) +
    ggplot2::facet_grid(. ~ genotype, scales = "free_x", space = "free_x") +
    scale_fill_rwb(midpoint = 0, limits = c(-lim, lim), oob = scales::squish,
                   name = "row z") +
    ggplot2::scale_x_discrete(drop = TRUE) +
    ggplot2::labs(x = "24M run index", y = NULL, title = title,
                  subtitle = "Rows are top-ranked sites in the mutant-tau amyloid contrast") +
    theme_tau(base_size = 8.5) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6.5),
      axis.text.y = ggplot2::element_text(size = 6.7, lineheight = 0.9),
      panel.spacing.x = grid::unit(0.16, "lines"),
      legend.position = "bottom"
    )
}

phospho_modality_plot <- function(phospho, title = "Bulk phosphoproteome") {
  stopifnot(is.list(phospho), is.data.frame(phospho$volcano), is.data.frame(phospho$heatmap))
  volcano <- modality_volcano_plot(
    phospho$volcano,
    title = "Phosphosite differential abundance",
    x_lab = "site log2FC  NLGF_P301S vs P301S",
    alpha = phospho$provenance$alpha %||% 0.10
  )
  heatmap <- phospho_site_heatmap_plot(phospho$heatmap, title = "Top phosphosite heatmap")
  patchwork::wrap_plots(list(volcano, heatmap), ncol = 2, widths = c(0.95, 1.05)) +
    patchwork::plot_annotation(title = title)
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
