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

subpopulation_colours <- c(
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

amyloid_status_colours <- c(
  `NLGF-` = "#6F7782",
  `NLGF+` = "#C8841C"
)

figure7_score_fill_colours <- c(
  MAPTKI = "#245A9A",
  P301S  = "#C65A1E"
)

tau_report_base_size <- 22
tau_report_axis_size <- 18.8
tau_report_dense_axis_size <- 16.8
tau_report_label_size <- 5.3
tau_report_dense_label_size <- 4.6

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
theme_tau <- function(base_size = tau_report_base_size, base_family = "") {
  set_tau_plot_defaults()
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "#E6E3DD", linewidth = 0.25),
      axis.line        = ggplot2::element_line(colour = "#3B3B37", linewidth = 0.25),
      plot.title       = ggplot2::element_text(face = "bold", colour = "#262B2F",
                                               size = ggplot2::rel(0.9)),
      plot.subtitle    = ggplot2::element_text(colour = "#53504A",
                                               size = ggplot2::rel(0.76)),
      strip.background = ggplot2::element_rect(fill = "#F2F0EA", colour = NA),
      strip.text       = ggplot2::element_text(face = "bold", colour = "#262B2F"),
      axis.title       = ggplot2::element_text(colour = "#333333"),
      axis.text        = ggplot2::element_text(colour = "#333333"),
      legend.title     = ggplot2::element_text(colour = "#333333"),
      legend.text      = ggplot2::element_text(colour = "#333333")
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

# Microglia-subpopulation colour/fill scales: keep coherent biology colours stable wherever
# subpopulations are mapped. `breaks=` can be narrowed to present_sub in a qmd.
scale_colour_subpopulation <- function(..., breaks = microglia_subpopulation_levels) {
  ggplot2::scale_colour_manual(values = subpopulation_colours, limits = microglia_subpopulation_levels,
                               breaks = breaks, drop = FALSE, ...)
}
scale_fill_subpopulation <- function(..., breaks = microglia_subpopulation_levels) {
  ggplot2::scale_fill_manual(values = subpopulation_colours, limits = microglia_subpopulation_levels,
                             breaks = breaks, drop = FALSE, ...)
}
scale_color_subpopulation <- scale_colour_subpopulation

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
                                         point_colour = "#6F7782", label_colour = "#A63A50",
                                         xy_lim = NULL) {
  stopifnot(is.data.frame(df), all(c("x", "y", label_col) %in% names(df)),
            is.null(xy_lim) ||
              (is.numeric(xy_lim) && length(xy_lim) == 1L && is.finite(xy_lim) && xy_lim > 0))
  label_cutoff <- attr(df, "offdiag_cutoff", exact = TRUE)
  label_cutoff_source <- attr(df, "offdiag_cutoff_source", exact = TRUE)
  threshold_cutoff <- if (!is.null(label_cutoff) && length(label_cutoff) == 1L &&
                            is.finite(label_cutoff)) {
    as.numeric(label_cutoff)
  } else {
    NA_real_
  }
  df    <- df[is.finite(df$x) & is.finite(df$y), , drop = FALSE]
  stopifnot(nrow(df) > 0L)
  lim   <- max(abs(c(df$x, df$y)), na.rm = TRUE)
  lim   <- if (is.finite(lim) && lim > 0) lim else 1
  if (!is.null(xy_lim)) lim <- as.numeric(xy_lim)
  # one label per feature: bulk assays reuse a site_id / gene across measured rows, so keep each label's
  # most-divergent threshold-passing instance (order is |y-x| desc) -> no duplicate repel labels.
  top   <- modality_scatter_label_rows(df, n_label = n_label, label_col = label_col,
                                       tail_quantile = label_tail_quantile,
                                       cutoff = label_cutoff,
                                       cutoff_source = label_cutoff_source)
  label_n <- nrow(top)
  if (label_n >= 80L && is.null(xy_lim)) lim <- lim * 1.08
  label_size <- if (label_n >= 150L) {
    tau_report_dense_label_size
  } else if (label_n >= 80L) {
    tau_report_dense_label_size
  } else if (label_n >= 40L) {
    4.8
  } else {
    tau_report_label_size
  }
  label_box_padding <- if (label_n >= 80L) 0.08 else 0.25
  label_point_padding <- if (label_n >= 80L) 0.03 else 0.10
  p <- ggplot2::ggplot(df, ggplot2::aes(x, y)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.25) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.25) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                         colour = "grey55", linewidth = 0.4)
  if (is.finite(threshold_cutoff)) {
    p <- p +
      ggplot2::geom_abline(slope = 1, intercept = c(-threshold_cutoff, threshold_cutoff),
                           linetype = "dotted", colour = "#B7AA97", linewidth = 0.4)
  }
  p +
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
    ggplot2::labs(x = x_lab, y = y_lab, title = title) +
    theme_tau()
}

# Functional-category aggregate scores for empirical off-diagonal genes/proteins in the
# amyloid-response scatter; phosphoproteomics uses the displayed parent-protein means. Rows
# are modality-specific role categories, and facets are the modalities with retained summary
# rows. Each segment connects the aggregate amyloid logFC under MAPTKI to the aggregate
# amyloid logFC under P301S; segment colour is the log2FC difference between tau backgrounds.
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
  modality_levels <- levels(droplevels(x$modality))
  if (is.null(modality_levels)) modality_levels <- unique(as.character(x$modality))
  modality_plot_levels <- modality_levels
  modality_plot_levels[modality_plot_levels == "snRNAseq"] <- "snRNA"
  x$modality_plot <- factor(
    modality_plot_levels[match(as.character(x$modality), modality_levels)],
    levels = modality_plot_levels
  )
  facet_ncol <- max(1L, length(modality_plot_levels))
  lim <- max(abs(c(x$score_maptki, x$score_p301s)), na.rm = TRUE)
  lim <- if (is.finite(lim) && lim > 0) lim else 1
  lim <- lim * 1.12
  delta_breaks <- -3:3
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
    data.frame(x[, c("modality_plot", "group_label_plot", "n_feature", "delta"), drop = FALSE],
               background = "MAPTKI", score = x$score_maptki, stringsAsFactors = FALSE),
    data.frame(x[, c("modality_plot", "group_label_plot", "n_feature", "delta"), drop = FALSE],
               background = "P301S", score = x$score_p301s, stringsAsFactors = FALSE)
  )
  point_df$background <- factor(point_df$background, levels = c("MAPTKI", "P301S"))
  bg_fill <- figure7_score_fill_colours
  point_outline <- "#262B2F"
  visual_scale <- 1.08
  base_size <- tau_report_base_size * visual_scale
  axis_y_size <- tau_report_axis_size * visual_scale
  title_size <- base_size * 1.06

  ggplot2::ggplot(x, ggplot2::aes(y = group_label_plot)) +
    ggplot2::geom_vline(xintercept = 0, colour = "#BFB8AA", linewidth = 0.35) +
    ggplot2::geom_segment(
      ggplot2::aes(x = score_maptki, xend = score_p301s, yend = group_label_plot,
                   colour = delta),
      linewidth = 0.85, lineend = "round") +
    ggplot2::geom_point(
      data = point_df,
      ggplot2::aes(x = score, size = n_feature),
      shape = 21, colour = "white", fill = "white", stroke = 1.85,
      show.legend = FALSE) +
    ggplot2::geom_point(
      data = point_df,
      ggplot2::aes(x = score, fill = background, size = n_feature),
      shape = 21, colour = point_outline, stroke = 0.82, alpha = 0.98) +
    scale_colour_rwb(midpoint = 0, limits = c(-3, 3), breaks = delta_breaks,
                     labels = function(x) sprintf("%d", x),
                     oob = scales::squish, name = "log2FC difference") +
    ggplot2::scale_fill_manual(values = bg_fill, breaks = names(bg_fill),
                               labels = c(MAPTKI = "NLGF_MAPTKI", P301S = "NLGF_P301S"),
                               name = "genotype") +
    ggplot2::scale_size_continuous(range = c(6.2, 18.5), breaks = size_breaks,
                                   name = "scored items") +
    ggplot2::scale_x_continuous(limits = c(-lim, lim), breaks = scales::breaks_width(2),
                                oob = scales::squish) +
    ggplot2::scale_y_discrete(drop = TRUE) +
    ggplot2::facet_wrap(ggplot2::vars(modality_plot), ncol = facet_ncol, scales = "free_y") +
    ggplot2::guides(
      colour = ggplot2::guide_colourbar(order = 1, barheight = grid::unit(0.52, "lines"),
                                        barwidth = grid::unit(9.2, "lines"),
                                        theme = ggplot2::theme(
                                          legend.spacing.x = grid::unit(0.65, "lines"),
                                          legend.title = ggplot2::element_text(
                                            margin = ggplot2::margin(r = 30)
                                          )
                                        )),
      fill = ggplot2::guide_legend(order = 2,
                                   override.aes = list(size = 6.4, shape = 21,
                                                       colour = point_outline)),
      size = ggplot2::guide_legend(order = 3)
    ) +
    ggplot2::labs(x = NULL, y = NULL, title = title) +
    theme_tau(base_size = base_size) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = axis_y_size, lineheight = 0.92),
      panel.grid.major.y = ggplot2::element_line(colour = "#ECE8DF", linewidth = 0.35),
      legend.title = ggplot2::element_text(size = base_size),
      legend.text = ggplot2::element_text(size = base_size * 0.8),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.box.just = "center",
      legend.spacing.y = grid::unit(0.1, "lines"),
      plot.title = ggplot2::element_text(face = "bold", colour = "#262B2F",
                                         size = title_size, hjust = 0.64,
                                         margin = ggplot2::margin(t = 10, b = 26)),
      plot.title.position = "plot"
    )
}

# Modality-native descriptive plates -------------------------------------------------------
# Live report uses the GeoMx DAM-gene AOI heatmap plus one bulk context plate:
# proteome PCA beside the phosphoproteome abundance heatmap.

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

geomx_normalization_rle_plot <- function(normalization,
                                         title = "GeoMx normalization and RLE") {
  stopifnot(is.list(normalization), is.data.frame(normalization$distribution),
            is.data.frame(normalization$rle), is.data.frame(normalization$background),
            is.list(normalization$voom), is.data.frame(normalization$voom$points),
            is.data.frame(normalization$voom$line))
  distribution <- normalization$distribution
  need_dist <- c("slide", "method", "q05", "q25", "q50", "q75", "q95")
  miss_dist <- setdiff(need_dist, names(distribution))
  if (length(miss_dist)) {
    stop("GeoMx normalization distribution missing columns: ",
         paste(miss_dist, collapse = ", "), call. = FALSE)
  }
  distribution <- distribution[stats::complete.cases(distribution[, need_dist]), ,
                               drop = FALSE]
  if (!nrow(distribution)) stop("GeoMx normalization distribution has no finite rows",
                                call. = FALSE)
  distribution$method <- factor(as.character(distribution$method),
                                levels = c("Raw logCPM", "TMM logCPM"))

  method_colours <- c(`Raw logCPM` = "#6F7782", `TMM logCPM` = "#0B7A75")
  dist_pos <- ggplot2::position_dodge(width = 0.55)
  distribution_plot <- ggplot2::ggplot(
    distribution, ggplot2::aes(slide, q50, colour = method)
  ) +
    ggplot2::geom_linerange(ggplot2::aes(ymin = q05, ymax = q95),
                            position = dist_pos, alpha = 0.28, linewidth = 0.35) +
    ggplot2::geom_pointrange(ggplot2::aes(ymin = q25, ymax = q75),
                             position = dist_pos, alpha = 0.82, linewidth = 0.42,
                             size = 0.55) +
    ggplot2::scale_colour_manual(values = method_colours, drop = FALSE,
                                 name = "scale") +
    ggplot2::labs(
      x = "slide", y = "logCPM quantile",
      title = "Raw versus TMM logCPM",
      subtitle = "Per-AOI gene-value median with IQR and 5-95% range"
    ) +
    theme_tau(base_size = 8.8) +
    ggplot2::theme(legend.position = "bottom",
                   axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))

  rle <- normalization$rle
  need_rle <- c("slide", "segment", "genotype", "q10", "q25", "q50", "q75", "q90")
  miss_rle <- setdiff(need_rle, names(rle))
  if (length(miss_rle)) {
    stop("GeoMx RLE data missing columns: ", paste(miss_rle, collapse = ", "),
         call. = FALSE)
  }
  rle <- rle[stats::complete.cases(rle[, need_rle]), , drop = FALSE]
  if (!nrow(rle)) stop("GeoMx RLE data has no finite rows", call. = FALSE)
  rle$genotype <- factor(as.character(rle$genotype), levels = genotype_levels)
  rle$segment <- factor(rle$segment)
  segment_levels <- levels(droplevels(rle$segment))
  shape_values <- stats::setNames(
    rep(c(16, 17, 15, 18, 3, 4, 7, 8, 0, 1, 2, 5, 6), length.out = length(segment_levels)),
    segment_levels
  )
  rle_pos <- ggplot2::position_jitter(width = 0.13, height = 0, seed = 42L)
  rle_plot <- ggplot2::ggplot(rle, ggplot2::aes(slide, q50)) +
    ggplot2::geom_hline(yintercept = 0, colour = "#BFB8AA", linewidth = 0.3) +
    ggplot2::geom_linerange(ggplot2::aes(ymin = q10, ymax = q90, colour = genotype),
                            position = rle_pos, alpha = 0.32, linewidth = 0.38) +
    ggplot2::geom_pointrange(
      ggplot2::aes(ymin = q25, ymax = q75, colour = genotype, shape = segment),
      position = rle_pos, alpha = 0.88, linewidth = 0.42, size = 0.78
    ) +
    scale_colour_genotype(name = "genotype") +
    ggplot2::scale_shape_manual(values = shape_values, name = "segment", drop = FALSE) +
    ggplot2::labs(
      x = "slide", y = "relative log expression",
      title = "TMM RLE by AOI",
      subtitle = "Median with IQR and 10-90% range after subtracting each gene median"
    ) +
    theme_tau(base_size = 8.8) +
    ggplot2::theme(legend.position = "none",
                   axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))

  background <- normalization$background
  need_bg <- c("slide", "segment", "genotype", "q3_factor", "neg_background")
  miss_bg <- setdiff(need_bg, names(background))
  if (length(miss_bg)) {
    stop("GeoMx background data missing columns: ", paste(miss_bg, collapse = ", "),
         call. = FALSE)
  }
  background <- background[is.finite(background$q3_factor) &
                             background$q3_factor > 0 &
                             is.finite(background$neg_background) &
                             background$neg_background > 0, , drop = FALSE]
  if (!nrow(background)) stop("GeoMx background data has no positive finite rows",
                              call. = FALSE)
  background$genotype <- factor(as.character(background$genotype), levels = genotype_levels)
  background$segment <- factor(background$segment)
  rho <- normalization$provenance$q3_neg_background_spearman %||% NA_real_
  rho_label <- if (is.finite(rho)) {
    sprintf("Spearman rho on log10 scale = %.2f", rho)
  } else {
    "Spearman rho unavailable"
  }
  background_plot <- ggplot2::ggplot(
    background, ggplot2::aes(q3_factor, neg_background, colour = genotype, shape = segment)
  ) +
    ggplot2::geom_point(size = 1.65, alpha = 0.82) +
    ggplot2::scale_x_log10(labels = scales::label_number()) +
    ggplot2::scale_y_log10(labels = scales::label_number()) +
    scale_colour_genotype(name = "genotype") +
    ggplot2::scale_shape_manual(values = shape_values, name = "segment", drop = FALSE) +
    ggplot2::facet_wrap(ggplot2::vars(slide), nrow = 1) +
    ggplot2::labs(
      x = "Q3 normalization factor", y = "negative-control background",
      title = "Q3 factor versus background",
      subtitle = rho_label
    ) +
    theme_tau(base_size = 8.6) +
    ggplot2::theme(legend.position = "bottom", legend.box = "horizontal")

  voom_points <- normalization$voom$points
  voom_line <- normalization$voom$line
  need_voom <- c("mean_log_count", "sqrt_sd")
  miss_voom <- setdiff(need_voom, names(voom_points))
  if (length(miss_voom) || !all(need_voom %in% names(voom_line))) {
    stop("GeoMx voom trend missing columns", call. = FALSE)
  }
  voom_points <- voom_points[stats::complete.cases(voom_points[, need_voom]), ,
                             drop = FALSE]
  voom_line <- voom_line[stats::complete.cases(voom_line[, need_voom]), , drop = FALSE]
  if (!nrow(voom_points) || !nrow(voom_line)) {
    stop("GeoMx voom trend has no finite rows", call. = FALSE)
  }
  voom_plot <- ggplot2::ggplot(voom_points, ggplot2::aes(mean_log_count, sqrt_sd)) +
    ggplot2::geom_point(colour = "#6F7782", alpha = 0.28, size = 0.48) +
    ggplot2::geom_line(data = voom_line, colour = "#A63A50", linewidth = 0.8) +
    ggplot2::labs(
      x = normalization$voom$labels$x %||% "log2 count size",
      y = normalization$voom$labels$y %||% "sqrt standard deviation",
      title = "voom mean-variance trend",
      subtitle = "Trend from the primary slide-adjusted, bio-unit-blocked GeoMx fit"
    ) +
    theme_tau(base_size = 8.8)

  subtitle <- sprintf(
    "%s filter-passing genes across %s AOIs; primary model unchanged: limma-voom + slide fixed effect + bio-unit duplicateCorrelation",
    format(normalization$provenance$n_kept_features %||% nrow(voom_points), big.mark = ","),
    format(normalization$provenance$n_aoi %||% nrow(background), big.mark = ",")
  )
  patchwork::wrap_plots(
    list(distribution_plot, rle_plot, background_plot, voom_plot),
    ncol = 2
  ) + patchwork::plot_annotation(title = title, subtitle = subtitle)
}

geomx_ordination_plot <- function(ordination, title = "GeoMx ordination") {
  stopifnot(is.list(ordination), is.data.frame(ordination$sample),
            is.data.frame(ordination$scree), is.data.frame(ordination$loadings))
  sample <- ordination$sample
  need_sample <- c("slide", "segment", "genotype", "pc1", "pc2", "pc1_var", "pc2_var",
                   "mds1", "mds2", "mds1_var", "mds2_var")
  miss_sample <- setdiff(need_sample, names(sample))
  if (length(miss_sample)) {
    stop("GeoMx ordination sample data missing columns: ",
         paste(miss_sample, collapse = ", "), call. = FALSE)
  }
  sample <- sample[is.finite(sample$pc1) & is.finite(sample$pc2) &
                     is.finite(sample$mds1) & is.finite(sample$mds2), , drop = FALSE]
  if (!nrow(sample)) stop("GeoMx ordination sample data has no finite AOIs", call. = FALSE)
  sample$genotype <- factor(as.character(sample$genotype), levels = genotype_levels)
  sample$slide <- factor(sample$slide)
  sample$segment <- factor(sample$segment)
  stopifnot(!anyNA(sample$genotype), !anyNA(sample$slide), !anyNA(sample$segment))

  segment_levels <- levels(droplevels(sample$segment))
  shape_values <- stats::setNames(
    rep(c(16, 17, 15, 18, 3, 4, 7, 8, 0, 1, 2, 5, 6), length.out = length(segment_levels)),
    segment_levels
  )
  axis_label <- function(axis, value) {
    if (is.finite(value)) sprintf("%s (%.1f%%)", axis, 100 * value) else axis
  }
  pc1_var <- unique(sample$pc1_var[is.finite(sample$pc1_var)])
  pc2_var <- unique(sample$pc2_var[is.finite(sample$pc2_var)])
  mds1_var <- unique(sample$mds1_var[is.finite(sample$mds1_var)])
  mds2_var <- unique(sample$mds2_var[is.finite(sample$mds2_var)])
  pc1_var <- if (length(pc1_var)) pc1_var[[1L]] else NA_real_
  pc2_var <- if (length(pc2_var)) pc2_var[[1L]] else NA_real_
  mds1_var <- if (length(mds1_var)) mds1_var[[1L]] else NA_real_
  mds2_var <- if (length(mds2_var)) mds2_var[[1L]] else NA_real_

  ord_theme <- theme_tau(base_size = 8.4) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      panel.spacing.x = grid::unit(0.5, "lines"),
      axis.text = ggplot2::element_text(size = 6.8)
    )
  scatter_guides <- ggplot2::guides(
    colour = ggplot2::guide_legend(order = 1, override.aes = list(size = 2.6)),
    shape = if (length(segment_levels) > 1L) {
      ggplot2::guide_legend(order = 2, override.aes = list(size = 2.6))
    } else {
      "none"
    }
  )

  pca_plot <- ggplot2::ggplot(sample, ggplot2::aes(pc1, pc2)) +
    ggplot2::geom_hline(yintercept = 0, colour = "#D5D0C4", linewidth = 0.25) +
    ggplot2::geom_vline(xintercept = 0, colour = "#D5D0C4", linewidth = 0.25) +
    ggplot2::geom_point(ggplot2::aes(colour = genotype, shape = segment),
                        size = 1.35, alpha = 0.88) +
    scale_colour_genotype(name = "genotype") +
    ggplot2::scale_shape_manual(values = shape_values, name = "segment", drop = FALSE) +
    ggplot2::facet_wrap(ggplot2::vars(slide), nrow = 1) +
    ggplot2::labs(
      x = axis_label("PC1", pc1_var),
      y = axis_label("PC2", pc2_var),
      title = "PCA of AOI expression",
      subtitle = "TMM logCPM, top variable filter-passing genes; facet = slide"
    ) +
    ord_theme + scatter_guides

  mds_plot <- ggplot2::ggplot(sample, ggplot2::aes(mds1, mds2)) +
    ggplot2::geom_hline(yintercept = 0, colour = "#D5D0C4", linewidth = 0.25) +
    ggplot2::geom_vline(xintercept = 0, colour = "#D5D0C4", linewidth = 0.25) +
    ggplot2::geom_point(ggplot2::aes(colour = genotype, shape = segment),
                        size = 1.35, alpha = 0.88) +
    scale_colour_genotype(name = "genotype") +
    ggplot2::scale_shape_manual(values = shape_values, name = "segment", drop = FALSE) +
    ggplot2::facet_wrap(ggplot2::vars(slide), nrow = 1) +
    ggplot2::labs(
      x = axis_label("MDS1", mds1_var),
      y = axis_label("MDS2", mds2_var),
      title = "Classical MDS of AOI distances",
      subtitle = "Euclidean distances on the same scaled expression matrix"
    ) +
    ord_theme + scatter_guides

  scree <- ordination$scree
  need_scree <- c("pc", "pc_num", "variance_percent")
  miss_scree <- setdiff(need_scree, names(scree))
  if (length(miss_scree)) {
    stop("GeoMx ordination scree data missing columns: ",
         paste(miss_scree, collapse = ", "), call. = FALSE)
  }
  scree <- scree[is.finite(scree$pc_num) & is.finite(scree$variance_percent), ,
                 drop = FALSE]
  if (nrow(scree) < 2L) stop("GeoMx ordination scree data has too few PCs",
                             call. = FALSE)
  scree$pc <- factor(as.character(scree$pc), levels = as.character(scree$pc))
  scree_plot <- ggplot2::ggplot(scree, ggplot2::aes(pc_num, variance_percent)) +
    ggplot2::geom_line(colour = "#3F5F7F", linewidth = 0.55) +
    ggplot2::geom_point(colour = "#3F5F7F", fill = "#F8F5ED", shape = 21,
                        size = 2.0, stroke = 0.45) +
    ggplot2::scale_x_continuous(breaks = scree$pc_num, labels = scree$pc) +
    ggplot2::labs(
      x = "component", y = "variance explained (%)",
      title = "PCA variance"
    ) +
    theme_tau(base_size = 8.8)

  loadings <- ordination$loadings
  need_loading <- c("symbol", "pc", "pc_num", "loading", "abs_loading", "direction")
  miss_loading <- setdiff(need_loading, names(loadings))
  if (length(miss_loading)) {
    stop("GeoMx ordination loading data missing columns: ",
         paste(miss_loading, collapse = ", "), call. = FALSE)
  }
  loadings <- loadings[is.finite(loadings$loading) & is.finite(loadings$abs_loading), ,
                       drop = FALSE]
  if (!nrow(loadings)) stop("GeoMx ordination loading data has no finite rows",
                            call. = FALSE)
  loadings <- loadings[order(loadings$pc_num, loadings$loading,
                             loadings$symbol, method = "radix"), , drop = FALSE]
  loadings$symbol_pc <- paste(loadings$symbol, loadings$pc, sep = " | ")
  loadings$symbol_plot <- factor(loadings$symbol_pc, levels = rev(loadings$symbol_pc))
  loadings$direction <- factor(as.character(loadings$direction),
                               levels = c("negative", "positive"))
  loading_plot <- ggplot2::ggplot(loadings, ggplot2::aes(loading, symbol_plot)) +
    ggplot2::geom_vline(xintercept = 0, colour = "#BFB8AA", linewidth = 0.3) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = loading, yend = symbol_plot,
                                       colour = direction),
                          linewidth = 0.46, alpha = 0.78) +
    ggplot2::geom_point(ggplot2::aes(colour = direction), size = 1.55, alpha = 0.95) +
    ggplot2::scale_colour_manual(values = c(negative = "#2F78A0",
                                            positive = "#A63A50"),
                                 drop = FALSE, name = "loading sign") +
    ggplot2::scale_y_discrete(labels = function(x) sub(" \\| PC[0-9]+$", "", x)) +
    ggplot2::facet_wrap(ggplot2::vars(pc), scales = "free_y", nrow = 1) +
    ggplot2::labs(
      x = "signed loading", y = NULL,
      title = "Top loading genes"
    ) +
    theme_tau(base_size = 8.3) +
    ggplot2::theme(
      legend.position = "bottom",
      axis.text.y = ggplot2::element_text(size = 6.6, lineheight = 0.9),
      panel.spacing.x = grid::unit(0.65, "lines")
    )

  subtitle <- sprintf(
    "%s variable genes from %s filter-passing genes across %s AOIs; ordination is descriptive and excludes no AOIs",
    format(ordination$provenance$n_variable_features %||% nrow(loadings), big.mark = ","),
    format(ordination$provenance$n_kept_features %||% NA_integer_, big.mark = ","),
    format(ordination$provenance$n_aoi %||% nrow(sample), big.mark = ",")
  )
  top <- patchwork::wrap_plots(list(pca_plot, mds_plot), ncol = 1, guides = "collect") &
    ggplot2::theme(legend.position = "bottom")
  bottom <- patchwork::wrap_plots(list(scree_plot, loading_plot), ncol = 2,
                                  widths = c(0.62, 1.38))
  patchwork::wrap_plots(list(top, bottom), ncol = 1, heights = c(1.42, 1)) +
    patchwork::plot_annotation(title = title, subtitle = subtitle)
}

geomx_gene_detection_plot <- function(detection, title = "GeoMx gene detectability") {
  stopifnot(is.list(detection), is.data.frame(detection$genes),
            is.data.frame(detection$filter_bins), is.data.frame(detection$marker_top),
            is.data.frame(detection$marker_labels), is.data.frame(detection$top_detected))
  genes <- detection$genes
  need_gene <- c("symbol", "mean_logcpm", "detect_fraction_min_count", "filter_status",
                 "filter_pass")
  miss_gene <- setdiff(need_gene, names(genes))
  if (length(miss_gene)) {
    stop("GeoMx gene detection data missing columns: ",
         paste(miss_gene, collapse = ", "), call. = FALSE)
  }
  genes <- genes[is.finite(genes$mean_logcpm) &
                   is.finite(genes$detect_fraction_min_count), , drop = FALSE]
  if (!nrow(genes)) stop("GeoMx gene detection data has no finite genes", call. = FALSE)
  genes$filter_status <- factor(as.character(genes$filter_status),
                                levels = c("low coverage", "filter passing"))

  marker_labels <- detection$marker_labels
  need_label <- c("symbol", "mean_logcpm", "detect_fraction_min_count", "marker_class")
  miss_label <- setdiff(need_label, names(marker_labels))
  if (length(miss_label)) {
    stop("GeoMx marker-label data missing columns: ",
         paste(miss_label, collapse = ", "), call. = FALSE)
  }
  marker_labels <- marker_labels[is.finite(marker_labels$mean_logcpm) &
                                   is.finite(marker_labels$detect_fraction_min_count), ,
                                 drop = FALSE]
  if (!nrow(marker_labels)) stop("GeoMx marker-label data has no finite rows",
                                 call. = FALSE)
  marker_labels$marker_class <- factor(as.character(marker_labels$marker_class),
                                       levels = c("Microglia", "Homeostatic", "DAM"))

  marker_colours <- c(Microglia = "#0B7A75",
                      Homeostatic = subpopulation_colours[["Homeostatic"]],
                      DAM = subpopulation_colours[["DAM"]])
  filter_colours <- c(`low coverage` = "#AFA89C", `filter passing` = "#3F5F7F")

  detect_plot <- ggplot2::ggplot(
    genes, ggplot2::aes(mean_logcpm, detect_fraction_min_count)
  ) +
    ggplot2::geom_point(ggplot2::aes(colour = filter_status),
                        size = 0.42, alpha = 0.34) +
    ggplot2::geom_point(
      data = marker_labels,
      ggplot2::aes(fill = marker_class),
      shape = 21, colour = "#20242A", stroke = 0.28, size = 1.8,
      inherit.aes = TRUE) +
    ggrepel::geom_text_repel(
      data = marker_labels,
      ggplot2::aes(label = symbol),
      size = 2.2, colour = "#20242A", max.overlaps = Inf, seed = 42L,
      min.segment.length = 0, segment.colour = "grey65", box.padding = 0.12,
      point.padding = 0.06, max.iter = 20000L, max.time = 3) +
    ggplot2::scale_colour_manual(values = filter_colours, drop = FALSE,
                                 name = "filter") +
    ggplot2::scale_fill_manual(values = marker_colours, drop = FALSE,
                               name = "labelled marker") +
    ggplot2::scale_y_continuous(labels = scales::label_percent(accuracy = 1),
                                limits = c(0, 1)) +
    ggplot2::labs(
      x = "mean TMM logCPM",
      y = sprintf("AOIs with count >= %s", detection$provenance$min_count %||% 5),
      title = "Expression versus detection",
      subtitle = "All GeoMx WTA genes; labelled points are high-detection microglia markers"
    ) +
    theme_tau(base_size = 8.7) +
    ggplot2::theme(legend.position = "bottom", legend.box = "horizontal")

  bins <- detection$filter_bins
  need_bin <- c("filter_status", "bin_mid", "n")
  miss_bin <- setdiff(need_bin, names(bins))
  if (length(miss_bin)) {
    stop("GeoMx gene-detection filter bins missing columns: ",
         paste(miss_bin, collapse = ", "), call. = FALSE)
  }
  bins <- bins[is.finite(bins$bin_mid) & is.finite(bins$n) & bins$n > 0, , drop = FALSE]
  if (!nrow(bins)) stop("GeoMx gene-detection filter bins have no positive rows",
                        call. = FALSE)
  bins$filter_status <- factor(as.character(bins$filter_status),
                               levels = c("low coverage", "filter passing"))
  filter_plot <- ggplot2::ggplot(bins, ggplot2::aes(bin_mid, n, colour = filter_status)) +
    ggplot2::geom_line(linewidth = 0.56, alpha = 0.85) +
    ggplot2::geom_point(size = 1.25, alpha = 0.9) +
    ggplot2::scale_colour_manual(values = filter_colours, drop = FALSE,
                                 name = "filter") +
    ggplot2::scale_x_continuous(labels = scales::label_percent(accuracy = 1),
                                limits = c(0, 1)) +
    ggplot2::scale_y_continuous(labels = scales::label_number(big.mark = ",")) +
    ggplot2::labs(
      x = sprintf("AOI fraction with count >= %s", detection$provenance$min_count %||% 5),
      y = "genes",
      title = "Low-coverage filter boundary",
      subtitle = "edgeR filterByExpr decision over the primary GeoMx design"
    ) +
    theme_tau(base_size = 8.7) +
    ggplot2::theme(legend.position = "bottom")

  marker_top <- detection$marker_top
  need_marker <- c("symbol", "marker_class", "mean_logcpm", "detect_fraction_min_count",
                   "filter_pass")
  miss_marker <- setdiff(need_marker, names(marker_top))
  if (length(miss_marker)) {
    stop("GeoMx marker-top data missing columns: ",
         paste(miss_marker, collapse = ", "), call. = FALSE)
  }
  marker_top <- marker_top[is.finite(marker_top$mean_logcpm) &
                             is.finite(marker_top$detect_fraction_min_count), ,
                           drop = FALSE]
  if (!nrow(marker_top)) stop("GeoMx marker-top data has no finite rows", call. = FALSE)
  marker_top$marker_class <- factor(as.character(marker_top$marker_class),
                                    levels = c("Microglia", "Homeostatic", "DAM"))
  marker_top <- marker_top[order(marker_top$marker_class,
                                 marker_top$detect_fraction_min_count,
                                 marker_top$mean_logcpm,
                                 marker_top$symbol,
                                 method = "radix"), , drop = FALSE]
  marker_top$symbol_key <- paste(marker_top$marker_class, marker_top$symbol, sep = " | ")
  marker_top$symbol_plot <- factor(marker_top$symbol_key, levels = unique(marker_top$symbol_key))
  marker_plot <- ggplot2::ggplot(marker_top, ggplot2::aes(detect_fraction_min_count,
                                                          symbol_plot)) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = detect_fraction_min_count,
                                       yend = symbol_plot),
                          colour = "#8A8174", linewidth = 0.36, alpha = 0.62) +
    ggplot2::geom_point(ggplot2::aes(fill = mean_logcpm),
                        shape = 21, colour = "#20242A", stroke = 0.24,
                        size = 1.9, alpha = 0.95) +
    ggplot2::scale_fill_gradient(low = "#F8F5ED", high = "#A63A50",
                                 name = "mean logCPM") +
    ggplot2::scale_x_continuous(labels = scales::label_percent(accuracy = 1),
                                limits = c(0, 1)) +
    ggplot2::scale_y_discrete(labels = function(x) sub("^.* \\| ", "", x)) +
    ggplot2::facet_grid(ggplot2::vars(marker_class), scales = "free_y", space = "free_y") +
    ggplot2::labs(
      x = sprintf("AOI fraction with count >= %s", detection$provenance$min_count %||% 5),
      y = NULL,
      title = "Marker-gene measurability",
      subtitle = "Top retained markers by detection within each signature"
    ) +
    theme_tau(base_size = 8.4) +
    ggplot2::theme(
      legend.position = "bottom",
      axis.text.y = ggplot2::element_text(size = 6.9, lineheight = 0.9),
      panel.spacing.y = grid::unit(0.45, "lines")
    )

  top_detected <- detection$top_detected
  need_top <- c("symbol", "mean_logcpm", "detect_fraction_min_count", "marker_class")
  miss_top <- setdiff(need_top, names(top_detected))
  if (length(miss_top)) {
    stop("GeoMx top-detected data missing columns: ",
         paste(miss_top, collapse = ", "), call. = FALSE)
  }
  top_detected <- top_detected[is.finite(top_detected$mean_logcpm) &
                                 is.finite(top_detected$detect_fraction_min_count), ,
                               drop = FALSE]
  if (!nrow(top_detected)) stop("GeoMx top-detected data has no finite rows",
                                call. = FALSE)
  top_detected <- top_detected[order(top_detected$mean_logcpm,
                                     top_detected$symbol,
                                     method = "radix"), , drop = FALSE]
  top_detected$symbol_plot <- factor(top_detected$symbol, levels = top_detected$symbol)
  top_detected$marker_hit <- ifelse(top_detected$marker_class == "other",
                                    "other", "marker")
  top_plot <- ggplot2::ggplot(top_detected, ggplot2::aes(mean_logcpm, symbol_plot)) +
    ggplot2::geom_segment(ggplot2::aes(x = min(mean_logcpm), xend = mean_logcpm,
                                       yend = symbol_plot),
                          colour = "#8A8174", linewidth = 0.38, alpha = 0.65) +
    ggplot2::geom_point(ggplot2::aes(fill = marker_hit),
                        shape = 21, colour = "#20242A", stroke = 0.26,
                        size = 2.0, alpha = 0.94) +
    ggplot2::scale_fill_manual(values = c(marker = "#0B7A75", other = "#C8841C"),
                               name = NULL) +
    ggplot2::labs(
      x = "mean TMM logCPM", y = NULL,
      title = "Highest-detected genes",
      subtitle = "Filter-passing genes ranked by detection and expression"
    ) +
    theme_tau(base_size = 8.5) +
    ggplot2::theme(legend.position = "bottom",
                   axis.text.y = ggplot2::element_text(size = 6.9))

  subtitle <- sprintf(
    "%s/%s genes pass filterByExpr(min.count=%s); marker genes present/pass: %s",
    format(detection$provenance$n_kept_features %||% sum(genes$filter_pass), big.mark = ","),
    format(detection$provenance$n_input_features %||% nrow(genes), big.mark = ","),
    detection$provenance$min_count %||% 5,
    paste(sprintf("%s %s/%s",
                  names(detection$provenance$n_marker_filter_passing),
                  detection$provenance$n_marker_filter_passing,
                  detection$provenance$n_marker_present),
          collapse = "; ")
  )
  left <- patchwork::wrap_plots(list(detect_plot, marker_plot), ncol = 1,
                                heights = c(1.0, 1.05))
  right <- patchwork::wrap_plots(list(filter_plot, top_plot), ncol = 1,
                                 heights = c(0.88, 1.12))
  patchwork::wrap_plots(list(left, right), ncol = 2, widths = c(1.25, 1)) +
    patchwork::plot_annotation(title = title, subtitle = subtitle)
}

geomx_sample_heatmap_plot <- function(sample_heatmap,
                                      title = NULL) {
  stopifnot(is.list(sample_heatmap), is.data.frame(sample_heatmap$heatmap),
            is.data.frame(sample_heatmap$sample), is.data.frame(sample_heatmap$features))
  heatmap <- sample_heatmap$heatmap
  sample <- sample_heatmap$sample
  need_heat <- c("symbol", "sample_rank", "symbol_plot", "z")
  miss_heat <- setdiff(need_heat, names(heatmap))
  if (length(miss_heat)) {
    stop("GeoMx sample heatmap data missing columns: ",
         paste(miss_heat, collapse = ", "), call. = FALSE)
  }
  heatmap <- heatmap[is.finite(heatmap$sample_rank) & is.finite(heatmap$z), ,
                     drop = FALSE]
  if (!nrow(heatmap)) stop("GeoMx sample heatmap has no finite rows", call. = FALSE)
  need_sample <- c("sample_rank", "genotype", "slide", "bio_unit", "roi",
                   "signed_response_score")
  miss_sample <- setdiff(need_sample, names(sample))
  if (length(miss_sample)) {
    stop("GeoMx sample heatmap sample data missing columns: ",
         paste(miss_sample, collapse = ", "), call. = FALSE)
  }
  sample <- sample[is.finite(sample$sample_rank) &
                     is.finite(sample$signed_response_score), , drop = FALSE]
  if (!nrow(sample)) stop("GeoMx sample heatmap has no finite samples", call. = FALSE)
  sample <- sample[order(sample$sample_rank, method = "radix"), , drop = FALSE]
  if (!"tau_background" %in% names(sample)) {
    sample$tau_background <- ifelse(
      as.character(sample$genotype) %in% c("P301S", "NLGF_P301S"),
      "P301S", "MAPTKI"
    )
  }
  if (!"amyloid_status" %in% names(sample)) {
    sample$amyloid_status <- ifelse(startsWith(as.character(sample$genotype), "NLGF"),
                                    "NLGF+", "NLGF-")
  }
  sample$tau_background <- factor(as.character(sample$tau_background),
                                  levels = names(tau_background_colours))
  sample$amyloid_status <- factor(as.character(sample$amyloid_status),
                                  levels = names(amyloid_status_colours))
  stopifnot(!anyNA(sample$tau_background), !anyNA(sample$amyloid_status))
  n_sample <- max(sample$sample_rank)
  x_limits <- c(0.5, n_sample + 0.5)

  id_palette <- function(values, saturation = 0.48, value = 0.70) {
    levels <- sort(unique(as.character(values)))
    if (!length(levels)) return(character())
    hue <- seq(0, 0.86, length.out = length(levels))
    stats::setNames(grDevices::hsv(hue, saturation, value), levels)
  }
  tracks <- sample_heatmap$metadata_tracks
  if (!is.data.frame(tracks)) {
    track_one <- function(track_id, track_label, track_group, track_type, values,
                          numeric_value = rep(NA_real_, length(values))) {
      data.frame(
        sample_rank = sample$sample_rank,
        track_id = track_id,
        track_label = track_label,
        track_group = track_group,
        track_type = track_type,
        source = "fallback",
        value = as.character(values),
        numeric_value = numeric_value,
        stringsAsFactors = FALSE
      )
    }
    tracks <- rbind(
      track_one("genotype", "genotype", "design / derived", "categorical", sample$genotype),
      track_one("tau_background", "tau", "design / derived", "categorical",
                sample$tau_background),
      track_one("amyloid_status", "amyloid", "design / derived", "categorical",
                sample$amyloid_status),
      track_one("slide", "slide", "design / derived", "categorical", sample$slide),
      track_one("bio_unit", "bio-unit", "design / derived", "categorical",
                sample$bio_unit),
      track_one("roi", "ROI", "ROI / spatial", "categorical", sample$roi),
      track_one("signed_response_score", "amyloid score", "design / derived", "numeric",
                sample$signed_response_score, sample$signed_response_score)
    )
  }
  need_tracks <- c("sample_rank", "track_id", "track_label", "track_group",
                   "track_type", "value", "numeric_value")
  miss_tracks <- setdiff(need_tracks, names(tracks))
  if (length(miss_tracks)) {
    stop("GeoMx metadata tracks missing columns: ",
         paste(miss_tracks, collapse = ", "), call. = FALSE)
  }
  tracks <- tracks[is.finite(tracks$sample_rank), , drop = FALSE]
  if (!nrow(tracks)) stop("GeoMx metadata track frame has no finite rows", call. = FALSE)
  tracks$track_id <- as.character(tracks$track_id)
  tracks$track_label <- as.character(tracks$track_label)
  tracks$track_group <- factor(as.character(tracks$track_group),
                               levels = unique(as.character(tracks$track_group)))
  track_levels <- unique(tracks$track_id)
  label_map <- stats::setNames(
    tracks$track_label[match(track_levels, tracks$track_id)],
    track_levels
  )
  tracks$track_id <- factor(tracks$track_id, levels = rev(track_levels))

  sequential <- grDevices::colorRampPalette(c("#F1EEE5", "#8BAFC0", "#1F6F8B"))(101)
  diverging <- grDevices::colorRampPalette(c("#2F78A0", "#F8F5ED", "#A63A50"))(101)
  tau_amyloid_colours <- c(no = "#6E6AA8", yes = "#B66A7A")
  tau_track_colours <- c(MAPTKI = tau_amyloid_colours[["no"]],
                         P301S = tau_amyloid_colours[["yes"]])
  amyloid_track_colours <- c(`NLGF-` = tau_amyloid_colours[["no"]],
                             `NLGF+` = tau_amyloid_colours[["yes"]])
  z_limit <- as.numeric(sample_heatmap$provenance$z_limit %||% 2.5)
  if (length(z_limit) != 1L || !is.finite(z_limit) || z_limit <= 0) z_limit <- 2.5
  colour_numeric <- function(x, diverge = FALSE, limit = NULL) {
    finite <- is.finite(x)
    fill <- rep("#D8D2C8", length(x))
    if (!any(finite)) return(fill)
    if (isTRUE(diverge)) {
      lim <- if (is.null(limit)) max(abs(x[finite]), na.rm = TRUE) else as.numeric(limit)
      if (!is.finite(lim) || lim <= 0) return(fill)
      scaled <- (pmax(pmin(x, lim), -lim) + lim) / (2 * lim)
      fill[finite] <- diverging[pmax(1L, pmin(101L, floor(scaled[finite] * 100) + 1L))]
      return(fill)
    }
    lim <- stats::quantile(x[finite], c(0.02, 0.98), names = FALSE, na.rm = TRUE)
    if (!all(is.finite(lim)) || lim[[1L]] == lim[[2L]]) lim <- range(x[finite])
    if (!all(is.finite(lim)) || lim[[1L]] == lim[[2L]]) return(fill)
    scaled <- (pmax(pmin(x, lim[[2L]]), lim[[1L]]) - lim[[1L]]) /
      (lim[[2L]] - lim[[1L]])
    fill[finite] <- sequential[pmax(1L, pmin(101L, floor(scaled[finite] * 100) + 1L))]
    fill
  }
  colour_categorical <- function(id, values) {
    values <- as.character(values)
    palette <- switch(
      id,
      genotype = genotype_colours,
      tau_background = tau_track_colours,
      amyloid_status = amyloid_track_colours,
      id_palette(values, saturation = 0.38, value = 0.72)
    )
    fill <- unname(palette[values])
    fill[is.na(fill)] <- "#D8D2C8"
    fill
  }
  track_id_chr <- as.character(tracks$track_id)
  tracks$fill <- "#D8D2C8"
  for (id in unique(track_id_chr)) {
    idx <- which(track_id_chr == id)
    type <- unique(as.character(tracks$track_type[idx]))
    type <- type[[1L]]
    if (identical(type, "numeric")) {
      diverge <- identical(id, "signed_response_score") ||
        startsWith(id, "signature_") || startsWith(id, "gene_")
      tracks$fill[idx] <- colour_numeric(
        as.numeric(tracks$numeric_value[idx]),
        diverge = diverge,
        limit = if (startsWith(id, "gene_")) z_limit else NULL
      )
    } else {
      tracks$fill[idx] <- colour_categorical(id, tracks$value[idx])
    }
  }
  if (any(is.na(tracks$fill))) stop("GeoMx metadata track palette has missing fills",
                                    call. = FALSE)
  tau_amyloid_key <- rep(NA_character_, nrow(tracks))
  tau_amyloid_key[track_id_chr == "tau_background"] <-
    ifelse(as.character(tracks$value[track_id_chr == "tau_background"]) == "P301S",
           "tau_amyloid_yes", "tau_amyloid_no")
  tau_amyloid_key[track_id_chr == "amyloid_status"] <-
    ifelse(as.character(tracks$value[track_id_chr == "amyloid_status"]) == "NLGF+",
           "tau_amyloid_yes", "tau_amyloid_no")
  tracks$fill_key <- tracks$fill
  tracks$fill_key[!is.na(tau_amyloid_key)] <- tau_amyloid_key[!is.na(tau_amyloid_key)]
  non_legend_fills <- unique(tracks$fill[is.na(tau_amyloid_key)])
  fill_values <- stats::setNames(non_legend_fills, non_legend_fills)
  fill_values <- c(
    fill_values,
    tau_amyloid_no = tau_amyloid_colours[["no"]],
    tau_amyloid_yes = tau_amyloid_colours[["yes"]]
  )

  genotype_values <- unique(as.character(tracks$value[track_id_chr == "genotype"]))
  genotype_levels_present <- genotype_levels[genotype_levels %in% genotype_values]
  legend_square_size <- 5.4
  genotype_legend_fills <- unname(genotype_colours[genotype_levels_present])
  if (!length(genotype_legend_fills)) genotype_legend_fills <- "#333333"
  genotype_key <- data.frame()
  if (length(genotype_levels_present)) {
    genotype_key <- data.frame(
      genotype = factor(genotype_levels_present, levels = genotype_levels_present),
      sample_rank = x_limits[[1L]],
      track_id = factor(levels(tracks$track_id)[[1L]], levels = levels(tracks$track_id)),
      stringsAsFactors = FALSE
    )
  }

  main_plot <- ggplot2::ggplot(tracks, ggplot2::aes(sample_rank, track_id, fill = fill_key)) +
    ggplot2::geom_tile(width = 1, height = 0.9, colour = "#F8F5ED", linewidth = 0.05,
                       key_glyph = ggplot2::draw_key_point) +
    ggplot2::geom_point(
      data = genotype_key,
      mapping = ggplot2::aes(sample_rank, track_id, colour = genotype),
      inherit.aes = FALSE,
      shape = 15,
      size = 0,
      alpha = 0,
      show.legend = TRUE
    ) +
    ggplot2::scale_fill_manual(
      values = fill_values,
      breaks = c("tau_amyloid_no", "tau_amyloid_yes"),
      labels = c("no", "yes"),
      name = "tau/amyloid",
      guide = ggplot2::guide_legend(
        order = 2,
        nrow = 1,
        byrow = TRUE,
        title.position = "left",
        label.position = "right",
        override.aes = list(shape = 22, size = legend_square_size, alpha = 1,
                            colour = "#F8F5ED", stroke = 0.12)
      )
    ) +
    ggplot2::scale_colour_manual(
      values = genotype_colours[genotype_levels_present],
      breaks = genotype_levels_present,
      name = "genotype",
      guide = ggplot2::guide_legend(
        order = 1,
        nrow = 1,
        byrow = TRUE,
        title.position = "left",
        label.position = "right",
        override.aes = list(shape = 22, size = legend_square_size, alpha = 1,
                            fill = genotype_legend_fills, colour = "#F8F5ED",
                            stroke = 0.12)
      )
    ) +
    ggplot2::scale_y_discrete(labels = label_map,
                              expand = ggplot2::expansion(add = 0.05)) +
    ggplot2::scale_x_continuous(limits = x_limits, expand = c(0, 0), breaks = NULL) +
    ggplot2::labs(
      x = NULL,
      y = NULL,
      title = title
    ) +
    theme_tau() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = tau_report_axis_size, lineheight = 0.9),
      panel.grid = ggplot2::element_blank(),
      legend.position = if (length(genotype_levels_present)) "top" else "none",
      legend.justification = "left",
      legend.box = "vertical",
      legend.box.just = "left",
      legend.direction = "horizontal",
      legend.title = ggplot2::element_text(face = "bold", size = tau_report_axis_size,
                                           colour = "#333333"),
      legend.text = ggplot2::element_text(size = tau_report_axis_size, colour = "#333333",
                                          margin = ggplot2::margin(l = 6, r = 10)),
      legend.key.size = grid::unit(0.29, "in"),
      legend.spacing.x = grid::unit(0.2, "in"),
      legend.spacing.y = grid::unit(0.05, "in"),
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      legend.box.margin = ggplot2::margin(0, 0, -2, 0),
      plot.margin = ggplot2::margin(1, 8, 6, 6)
    )
  main_plot
}

geomx_spatial_program_overlay_plot <- function(program_overlay,
                                               title = "GeoMx spatial program overlays") {
  stopifnot(is.list(program_overlay), is.data.frame(program_overlay$aoi),
            is.data.frame(program_overlay$programs),
            is.data.frame(program_overlay$genotype_summary))
  aoi <- program_overlay$aoi
  need <- c("slide", "genotype", "program_label", "x_coord", "y_coord", "score")
  miss <- setdiff(need, names(aoi))
  if (length(miss)) {
    stop("GeoMx spatial-program AOI data missing columns: ",
         paste(miss, collapse = ", "), call. = FALSE)
  }
  aoi <- aoi[is.finite(aoi$x_coord) & is.finite(aoi$y_coord) &
               is.finite(aoi$score), , drop = FALSE]
  if (!nrow(aoi)) stop("GeoMx spatial-program AOI data has no finite rows",
                       call. = FALSE)
  programs <- program_overlay$programs
  need_program <- c("program_label", "program_type", "n_scored_features")
  miss_program <- setdiff(need_program, names(programs))
  if (length(miss_program)) {
    stop("GeoMx spatial-program catalog missing columns: ",
         paste(miss_program, collapse = ", "), call. = FALSE)
  }
  program_levels <- as.character(programs$program_label)
  aoi$program_label <- factor(as.character(aoi$program_label), levels = program_levels)
  aoi$genotype <- factor(as.character(aoi$genotype), levels = genotype_levels)
  aoi$slide <- factor(aoi$slide)
  stopifnot(!anyNA(aoi$program_label), !anyNA(aoi$genotype), !anyNA(aoi$slide))
  score_lim <- program_overlay$provenance$z_limit %||% max(abs(aoi$score), na.rm = TRUE)
  score_lim <- if (is.finite(score_lim) && score_lim > 0) score_lim else 1

  spatial_map <- ggplot2::ggplot(aoi, ggplot2::aes(x_coord, y_coord)) +
    ggplot2::geom_point(
      ggplot2::aes(fill = score, colour = genotype),
      shape = 21, size = 1.55, stroke = 0.34, alpha = 0.94) +
    scale_fill_rwb(midpoint = 0, limits = c(-score_lim, score_lim),
                   oob = scales::squish, name = "row z") +
    scale_colour_genotype(name = "genotype") +
    ggplot2::facet_grid(ggplot2::vars(program_label), ggplot2::vars(slide)) +
    ggplot2::coord_equal() +
    ggplot2::scale_y_reverse() +
    ggplot2::labs(
      x = "ROI coordinate X", y = "ROI coordinate Y",
      title = title,
      subtitle = program_overlay$provenance$coordinate_status %||%
        "coordinate-only GeoMx overlay"
    ) +
    theme_tau(base_size = 7.8) +
    ggplot2::guides(
      fill = ggplot2::guide_colourbar(order = 1, barheight = grid::unit(0.42, "lines"),
                                      barwidth = grid::unit(4.0, "lines")),
      colour = ggplot2::guide_legend(order = 2, override.aes = list(fill = "white",
                                                                    size = 2.5))
    ) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      axis.text = ggplot2::element_text(size = 5.6),
      axis.title = ggplot2::element_text(size = 7.4),
      strip.text.x = ggplot2::element_text(size = 7.0),
      strip.text.y = ggplot2::element_text(size = 7.0, angle = 0),
      panel.grid = ggplot2::element_line(colour = "#ECE8DF", linewidth = 0.18),
      panel.spacing = grid::unit(0.32, "lines")
    )

  summary <- program_overlay$genotype_summary
  need_summary <- c("program_label", "genotype", "n_aoi", "median_score",
                    "q25_score", "q75_score")
  miss_summary <- setdiff(need_summary, names(summary))
  if (length(miss_summary)) {
    stop("GeoMx spatial-program summary missing columns: ",
         paste(miss_summary, collapse = ", "), call. = FALSE)
  }
  summary <- summary[is.finite(summary$n_aoi) & is.finite(summary$median_score) &
                       is.finite(summary$q25_score) & is.finite(summary$q75_score), ,
                     drop = FALSE]
  if (!nrow(summary)) stop("GeoMx spatial-program summary has no finite rows",
                           call. = FALSE)
  summary$program_label <- factor(as.character(summary$program_label),
                                  levels = program_levels)
  summary$genotype <- factor(as.character(summary$genotype), levels = genotype_levels)
  stopifnot(!anyNA(summary$program_label), !anyNA(summary$genotype))
  summary_plot <- ggplot2::ggplot(
    summary,
    ggplot2::aes(genotype, median_score, colour = genotype)
  ) +
    ggplot2::geom_hline(yintercept = 0, colour = "#BFB8AA", linewidth = 0.25) +
    ggplot2::geom_pointrange(
      ggplot2::aes(ymin = q25_score, ymax = q75_score),
      linewidth = 0.32, size = 0.82, alpha = 0.95) +
    scale_colour_genotype(guide = "none") +
    ggplot2::facet_wrap(ggplot2::vars(program_label), nrow = 1) +
    ggplot2::scale_y_continuous(limits = c(-score_lim, score_lim),
                                oob = scales::squish) +
    ggplot2::labs(
      x = NULL, y = "median row z",
      title = "Program score by genotype",
      subtitle = "Points show AOI median; ranges show interquartile interval"
    ) +
    theme_tau(base_size = 8.0) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1, size = 5.8),
      strip.text = ggplot2::element_text(size = 6.6),
      panel.spacing.x = grid::unit(0.28, "lines")
    )

  patchwork::wrap_plots(list(spatial_map, summary_plot), ncol = 1,
                        heights = c(3.2, 1.0))
}

geomx_contrast_diagnostic_plot <- function(contrast_diagnostics,
                                           title = "GeoMx contrast diagnostics") {
  stopifnot(is.list(contrast_diagnostics),
            is.data.frame(contrast_diagnostics$volcano),
            is.data.frame(contrast_diagnostics$labels),
            is.data.frame(contrast_diagnostics$support_counts),
            is.data.frame(contrast_diagnostics$interaction_top),
            is.data.frame(contrast_diagnostics$summary))
  volcano <- contrast_diagnostics$volcano
  need <- c("symbol", "contrast_label", "effect", "ave_expr", "neg_log10_fdr",
            "direction", "label", "label_show", "fdr", "ci_low", "ci_high",
            "supported")
  miss <- setdiff(need, names(volcano))
  if (length(miss)) {
    stop("GeoMx contrast diagnostic data missing columns: ",
         paste(miss, collapse = ", "), call. = FALSE)
  }
  volcano <- volcano[is.finite(volcano$effect) &
                       is.finite(volcano$ave_expr) &
                       is.finite(volcano$neg_log10_fdr) &
                       is.finite(volcano$fdr), , drop = FALSE]
  if (!nrow(volcano)) stop("GeoMx contrast diagnostic has no finite rows",
                           call. = FALSE)
  label_rows <- volcano[volcano$label_show %in% TRUE, , drop = FALSE]
  if (!nrow(label_rows)) stop("GeoMx contrast diagnostic has no label rows",
                              call. = FALSE)
  alpha <- contrast_diagnostics$provenance$alpha %||% 0.10
  effect_lim <- max(abs(c(volcano$effect, volcano$ci_low, volcano$ci_high)),
                    na.rm = TRUE)
  effect_lim <- if (is.finite(effect_lim) && effect_lim > 0) effect_lim else 1
  y_cut <- -log10(alpha)
  direction_colours <- c(down = "#2F78A0",
                         `not supported` = "#AFA89C",
                         up = "#A63A50")

  volcano_plot <- ggplot2::ggplot(
    volcano,
    ggplot2::aes(effect, neg_log10_fdr)
  ) +
    ggplot2::geom_hline(yintercept = y_cut, colour = "#BFB8AA",
                        linewidth = 0.25, linetype = "dotted") +
    ggplot2::geom_vline(xintercept = 0, colour = "#BFB8AA", linewidth = 0.25) +
    ggplot2::geom_point(ggplot2::aes(colour = direction),
                        size = 0.30, alpha = 0.28) +
    ggrepel::geom_text_repel(
      data = label_rows,
      ggplot2::aes(label = label),
      size = 1.75, colour = "#20242A", max.overlaps = Inf, seed = 42L,
      min.segment.length = 0, segment.colour = "grey65", box.padding = 0.08,
      point.padding = 0.03, max.iter = 20000L, max.time = 3) +
    ggplot2::scale_colour_manual(values = direction_colours, drop = FALSE,
                                 name = NULL) +
    ggplot2::scale_x_continuous(limits = c(-effect_lim, effect_lim),
                                oob = scales::squish) +
    ggplot2::facet_wrap(ggplot2::vars(contrast_label), nrow = 1) +
    ggplot2::labs(
      x = "GeoMx log2FC", y = "-log10 FDR",
      title = "Volcano diagnostics",
      subtitle = sprintf("Dotted line = FDR %.2f; labels are top-ranked per contrast",
                         alpha)
    ) +
    theme_tau(base_size = 7.7) +
    ggplot2::theme(
      legend.position = "bottom",
      panel.spacing.x = grid::unit(0.18, "lines"),
      strip.text = ggplot2::element_text(size = 6.1, lineheight = 0.88),
      axis.text = ggplot2::element_text(size = 5.9),
      axis.title = ggplot2::element_text(size = 7.0)
    )

  ma_plot <- ggplot2::ggplot(
    volcano,
    ggplot2::aes(ave_expr, effect)
  ) +
    ggplot2::geom_hline(yintercept = 0, colour = "#BFB8AA", linewidth = 0.25) +
    ggplot2::geom_point(ggplot2::aes(colour = direction),
                        size = 0.30, alpha = 0.28) +
    ggrepel::geom_text_repel(
      data = label_rows,
      ggplot2::aes(label = label),
      size = 1.65, colour = "#20242A", max.overlaps = Inf, seed = 42L,
      min.segment.length = 0, segment.colour = "grey65", box.padding = 0.08,
      point.padding = 0.03, max.iter = 20000L, max.time = 3) +
    ggplot2::scale_colour_manual(values = direction_colours, drop = FALSE,
                                 name = NULL) +
    ggplot2::scale_y_continuous(limits = c(-effect_lim, effect_lim),
                                oob = scales::squish) +
    ggplot2::facet_wrap(ggplot2::vars(contrast_label), nrow = 1) +
    ggplot2::labs(
      x = "mean log2 expression", y = "GeoMx log2FC",
      title = "MA diagnostics",
      subtitle = "Same primary limma-voom fit; colour marks FDR-supported direction"
    ) +
    theme_tau(base_size = 7.7) +
    ggplot2::theme(
      legend.position = "bottom",
      panel.spacing.x = grid::unit(0.18, "lines"),
      strip.text = ggplot2::element_text(size = 6.1, lineheight = 0.88),
      axis.text = ggplot2::element_text(size = 5.9),
      axis.title = ggplot2::element_text(size = 7.0)
    )

  counts <- contrast_diagnostics$support_counts
  need_counts <- c("contrast_label", "direction", "n", "signed_n")
  miss_counts <- setdiff(need_counts, names(counts))
  if (length(miss_counts)) {
    stop("GeoMx support-count data missing columns: ",
         paste(miss_counts, collapse = ", "), call. = FALSE)
  }
  counts <- counts[is.finite(counts$n) & is.finite(counts$signed_n), ,
                   drop = FALSE]
  if (!nrow(counts)) stop("GeoMx support-count data has no finite rows",
                          call. = FALSE)
  count_lim <- max(abs(counts$signed_n), na.rm = TRUE)
  count_lim <- if (is.finite(count_lim) && count_lim > 0) count_lim else 1
  count_plot <- ggplot2::ggplot(
    counts,
    ggplot2::aes(signed_n, contrast_label, colour = direction)
  ) +
    ggplot2::geom_vline(xintercept = 0, colour = "#BFB8AA", linewidth = 0.28) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = signed_n,
                                       yend = contrast_label),
                          linewidth = 0.65, alpha = 0.82) +
    ggplot2::geom_point(size = 2.0, alpha = 0.94) +
    ggplot2::scale_colour_manual(values = direction_colours[c("down", "up")],
                                 drop = FALSE, name = "direction") +
    ggplot2::scale_x_continuous(
      limits = c(-count_lim, count_lim),
      labels = function(x) abs(x),
      oob = scales::squish
    ) +
    ggplot2::labs(
      x = sprintf("FDR <= %.2f genes", alpha), y = NULL,
      title = "Signed support counts",
      subtitle = "Left = lower in numerator; right = higher in numerator"
    ) +
    theme_tau(base_size = 8.2) +
    ggplot2::theme(
      legend.position = "bottom",
      axis.text.y = ggplot2::element_text(size = 6.3, lineheight = 0.88)
    )

  interaction_top <- contrast_diagnostics$interaction_top
  need_ix <- c("symbol_plot", "effect", "ci_low", "ci_high", "direction", "fdr",
               "supported")
  miss_ix <- setdiff(need_ix, names(interaction_top))
  if (length(miss_ix)) {
    stop("GeoMx interaction-top data missing columns: ",
         paste(miss_ix, collapse = ", "), call. = FALSE)
  }
  interaction_top <- interaction_top[is.finite(interaction_top$effect) &
                                       is.finite(interaction_top$ci_low) &
                                       is.finite(interaction_top$ci_high) &
                                       is.finite(interaction_top$fdr), ,
                                     drop = FALSE]
  if (!nrow(interaction_top)) {
    stop("GeoMx interaction-top data has no finite rows", call. = FALSE)
  }
  ix_support_n <- sum(interaction_top$supported)
  interaction_plot <- ggplot2::ggplot(
    interaction_top,
    ggplot2::aes(effect, symbol_plot)
  ) +
    ggplot2::geom_vline(xintercept = 0, colour = "#BFB8AA", linewidth = 0.28) +
    ggplot2::geom_segment(ggplot2::aes(x = ci_low, xend = ci_high,
                                       yend = symbol_plot, colour = direction),
                          linewidth = 0.62, alpha = 0.86) +
    ggplot2::geom_point(ggplot2::aes(fill = direction),
                        shape = 21, colour = "#20242A", stroke = 0.24,
                        size = 2.0, alpha = 0.96) +
    ggplot2::scale_colour_manual(values = direction_colours, drop = FALSE,
                                 guide = "none") +
    ggplot2::scale_fill_manual(values = direction_colours, drop = FALSE,
                               name = "direction") +
    ggplot2::scale_x_continuous(limits = c(-effect_lim, effect_lim),
                                oob = scales::squish) +
    ggplot2::labs(
      x = "interaction log2FC with 95% CI", y = NULL,
      title = "Interaction emphasis",
      subtitle = sprintf("Top %s genes by interaction rank; %s FDR-supported",
                         nrow(interaction_top), ix_support_n)
    ) +
    theme_tau(base_size = 8.2) +
    ggplot2::theme(
      legend.position = "bottom",
      axis.text.y = ggplot2::element_text(size = 6.5, lineheight = 0.90)
    )

  subtitle <- contrast_diagnostics$provenance$display %||%
    "primary GeoMx topTables only; no new inference model"
  bottom <- patchwork::wrap_plots(list(count_plot, interaction_plot),
                                  ncol = 2, widths = c(0.95, 1.05))
  patchwork::wrap_plots(list(volcano_plot, ma_plot, bottom),
                        ncol = 1, heights = c(1.05, 1.05, 0.9)) +
    patchwork::plot_annotation(title = title, subtitle = subtitle)
}

geomx_roi_replicate_plot <- function(roi_replicates,
                                     title = "GeoMx ROI and block audit") {
  stopifnot(is.list(roi_replicates),
            is.data.frame(roi_replicates$support),
            is.data.frame(roi_replicates$block),
            is.data.frame(roi_replicates$pair_correlation),
            is.data.frame(roi_replicates$pair_summary))

  support <- roi_replicates$support
  need_support <- c("genotype", "bio_rep", "n_aoi", "slide", "present", "label")
  miss_support <- setdiff(need_support, names(support))
  if (length(miss_support)) {
    stop("GeoMx ROI-replicate support data missing columns: ",
         paste(miss_support, collapse = ", "), call. = FALSE)
  }
  support <- support[is.finite(support$n_aoi), , drop = FALSE]
  if (!nrow(support)) stop("GeoMx ROI-replicate support data has no finite rows",
                           call. = FALSE)
  support$genotype <- factor(as.character(support$genotype), levels = genotype_levels)
  support$bio_rep <- factor(as.character(support$bio_rep),
                            levels = sort(unique(as.character(support$bio_rep))))
  support$present <- support$present %in% TRUE
  stopifnot(!anyNA(support$genotype), !anyNA(support$bio_rep))
  support_present <- support[support$present & support$n_aoi > 0, , drop = FALSE]
  support_missing <- support[!support$present | support$n_aoi <= 0, , drop = FALSE]
  if (!nrow(support_present)) {
    stop("GeoMx ROI-replicate support data has no present bio-units", call. = FALSE)
  }
  slide_levels <- sort(unique(as.character(support$slide)))
  slide_palette <- stats::setNames(rep(tau_discrete_colours, length.out = length(slide_levels)),
                                   slide_levels)
  if ("absent" %in% slide_levels) slide_palette[["absent"]] <- "#D8D2C8"
  max_aoi <- max(support$n_aoi, na.rm = TRUE)
  max_aoi <- if (is.finite(max_aoi) && max_aoi > 0) max_aoi else 1

  support_plot <- ggplot2::ggplot(support, ggplot2::aes(bio_rep, genotype)) +
    ggplot2::geom_point(
      data = support_present,
      ggplot2::aes(size = n_aoi, fill = slide),
      shape = 21, colour = "#20242A", stroke = 0.34, alpha = 0.94) +
    ggplot2::geom_text(
      data = support_present,
      ggplot2::aes(label = label),
      size = 2.45, colour = "#20242A", fontface = "bold") +
    ggplot2::geom_text(
      data = support_missing,
      ggplot2::aes(label = label),
      size = 2.35, colour = "#8B857B", fontface = "bold") +
    ggplot2::scale_fill_manual(values = slide_palette, name = "slide", drop = FALSE) +
    ggplot2::scale_size_area(max_size = 9, limits = c(0, max_aoi),
                             name = "AOIs") +
    ggplot2::labs(
      x = "bio-replicate", y = NULL,
      title = "Bio-unit support",
      subtitle = "Point area and label = AOIs; grey 0 marks an absent genotype-by-replicate unit"
    ) +
    theme_tau(base_size = 8.4) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      panel.grid.minor = ggplot2::element_blank()
    )

  block <- roi_replicates$block
  need_block <- c("bio_unit_plot", "bio_unit", "genotype", "slide", "n_aoi",
                  "score_min", "score_max", "score_mean")
  miss_block <- setdiff(need_block, names(block))
  if (length(miss_block)) {
    stop("GeoMx ROI-replicate block data missing columns: ",
         paste(miss_block, collapse = ", "), call. = FALSE)
  }
  block <- block[is.finite(block$n_aoi) & is.finite(block$score_min) &
                   is.finite(block$score_max) & is.finite(block$score_mean), ,
                 drop = FALSE]
  if (!nrow(block)) stop("GeoMx ROI-replicate block data has no finite rows",
                         call. = FALSE)
  block$genotype <- factor(as.character(block$genotype), levels = genotype_levels)
  block$slide <- factor(block$slide)
  block$bio_unit_plot <- factor(block$bio_unit_plot,
                                levels = levels(block$bio_unit_plot))
  stopifnot(!anyNA(block$genotype), !anyNA(block$slide), !anyNA(block$bio_unit_plot))
  slide_block_levels <- levels(droplevels(block$slide))
  shape_values <- stats::setNames(
    rep(c(21, 22, 24, 25, 23, 21, 22, 24), length.out = length(slide_block_levels)),
    slide_block_levels
  )

  block_plot <- ggplot2::ggplot(block, ggplot2::aes(n_aoi, bio_unit_plot)) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, xend = n_aoi, yend = bio_unit_plot, colour = genotype),
      linewidth = 0.55, alpha = 0.78) +
    ggplot2::geom_point(
      ggplot2::aes(fill = genotype, shape = slide),
      colour = "#20242A", stroke = 0.24, size = 2.15, alpha = 0.95) +
    ggplot2::geom_text(ggplot2::aes(label = n_aoi),
                       nudge_x = 0.34, size = 2.15, colour = "#20242A") +
    scale_colour_genotype(guide = "none") +
    scale_fill_genotype(guide = "none") +
    ggplot2::scale_shape_manual(values = shape_values, guide = "none", drop = FALSE) +
    ggplot2::scale_x_continuous(limits = c(0, max(block$n_aoi) + 1),
                                breaks = seq(0, max(block$n_aoi) + 1, by = 2),
                                expand = ggplot2::expansion(mult = c(0.02, 0.03))) +
    ggplot2::labs(
      x = "AOIs in duplicateCorrelation block", y = NULL,
      title = "AOIs per bio-unit block",
      subtitle = "Every present bio-unit is a block; all AOIs are retained"
    ) +
    theme_tau(base_size = 8.2) +
    ggplot2::theme(
      legend.position = "bottom",
      axis.text.y = ggplot2::element_text(size = 6.3, lineheight = 0.88),
      panel.grid.minor = ggplot2::element_blank()
    )

  pairs <- roi_replicates$pair_correlation
  summary <- roi_replicates$pair_summary
  need_pairs <- c("pair_type", "correlation")
  miss_pairs <- setdiff(need_pairs, names(pairs))
  if (length(miss_pairs)) {
    stop("GeoMx ROI-replicate pair-correlation data missing columns: ",
         paste(miss_pairs, collapse = ", "), call. = FALSE)
  }
  pairs <- pairs[is.finite(pairs$correlation), , drop = FALSE]
  if (!nrow(pairs)) stop("GeoMx ROI-replicate pair-correlation data has no finite rows",
                         call. = FALSE)
  pair_levels <- levels(pairs$pair_type)
  if (is.null(pair_levels)) pair_levels <- unique(as.character(pairs$pair_type))
  pairs$pair_type <- factor(as.character(pairs$pair_type), levels = pair_levels)
  need_summary <- c("pair_type", "n_pairs", "median_correlation",
                    "q25_correlation", "q75_correlation")
  miss_summary <- setdiff(need_summary, names(summary))
  if (length(miss_summary)) {
    stop("GeoMx ROI-replicate pair-summary data missing columns: ",
         paste(miss_summary, collapse = ", "), call. = FALSE)
  }
  summary <- summary[is.finite(summary$n_pairs) &
                       is.finite(summary$median_correlation) &
                       is.finite(summary$q25_correlation) &
                       is.finite(summary$q75_correlation), , drop = FALSE]
  summary$pair_type <- factor(as.character(summary$pair_type), levels = pair_levels)
  stopifnot(!anyNA(pairs$pair_type), !anyNA(summary$pair_type))
  pair_palette <- stats::setNames(c("#3F5F7F", "#C8841C", "#7C838A")[seq_along(pair_levels)],
                                  pair_levels)
  summary$label <- paste0("n=", scales::comma(summary$n_pairs))
  summary$label_y <- pmin(1, summary$q75_correlation + 0.045)

  cor_plot <- ggplot2::ggplot(pairs, ggplot2::aes(pair_type, correlation)) +
    ggplot2::geom_hline(yintercept = 0, colour = "#BFB8AA", linewidth = 0.25) +
    ggplot2::geom_violin(ggplot2::aes(fill = pair_type),
                         width = 0.82, alpha = 0.34, colour = "#5D5851",
                         linewidth = 0.24, trim = TRUE) +
    ggplot2::geom_boxplot(width = 0.14, outlier.shape = NA, colour = "#20242A",
                          fill = "#F8F5ED", linewidth = 0.26) +
    ggplot2::geom_pointrange(
      data = summary,
      ggplot2::aes(x = pair_type, y = median_correlation, ymin = q25_correlation,
                   ymax = q75_correlation),
      inherit.aes = FALSE, linewidth = 0.34, size = 0.95, colour = "#20242A") +
    ggplot2::geom_text(
      data = summary,
      ggplot2::aes(x = pair_type, y = label_y, label = label),
      inherit.aes = FALSE, size = 2.45, colour = "#20242A") +
    ggplot2::scale_fill_manual(values = pair_palette, guide = "none") +
    ggplot2::coord_cartesian(ylim = c(-1, 1), clip = "off") +
    ggplot2::labs(
      x = NULL, y = "AOI-pair expression correlation",
      title = "DuplicateCorrelation block audit",
      subtitle = roi_replicates$provenance$correlation_basis %||%
        "Pearson correlations over top-variable TMM-logCPM genes"
    ) +
    theme_tau(base_size = 8.6) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 7.0, lineheight = 0.88),
      legend.position = "none",
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(4, 10, 6, 6)
    )

  prov <- roi_replicates$provenance
  dc <- prov$duplicate_correlation %||% NA_real_
  dc_label <- if (is.finite(dc)) sprintf("%.3f", dc) else "NA"
  segment_label <- paste(prov$segments %||% character(), collapse = ", ")
  if (!nzchar(segment_label)) segment_label <- "not recorded"
  block_range <- prov$block_size_range
  block_min <- if (!is.null(block_range) && length(block_range) >= 1L) {
    block_range[[1L]]
  } else {
    min(block$n_aoi)
  }
  block_max <- if (!is.null(block_range) && length(block_range) >= 2L) {
    block_range[[2L]]
  } else {
    max(block$n_aoi)
  }
  subtitle <- sprintf(
    "%s/%s expected bio-units; %s AOIs in %s duplicateCorrelation blocks (%s-%s AOIs); consensus = %s; segments = %s.",
    prov$n_present_bio_units %||% nrow(block),
    prov$n_expected_bio_units %||% NA_integer_,
    prov$n_aoi %||% nrow(support_present),
    prov$n_bio_units %||% nrow(block),
    block_min,
    block_max,
    dc_label,
    segment_label
  )
  if (isTRUE(prov$all_segments_single_level)) {
    subtitle <- paste(subtitle, "Paired segment differences are unavailable because all AOIs use one segment level.")
  }
  subtitle <- paste(strwrap(subtitle, width = 128), collapse = "\n")

  top <- patchwork::wrap_plots(list(support_plot, block_plot),
                               ncol = 2, widths = c(0.95, 1.05))
  patchwork::wrap_plots(list(top, cor_plot), ncol = 1, heights = c(1.05, 1.0)) +
    patchwork::plot_annotation(title = title, subtitle = subtitle)
}

geomx_decon_feasibility_plot <- function(decon_feasibility,
                                         title = "GeoMx decon feasibility") {
  stopifnot(is.list(decon_feasibility),
            is.data.frame(decon_feasibility$coverage),
            is.data.frame(decon_feasibility$aoi),
            is.data.frame(decon_feasibility$status_counts))

  coverage <- decon_feasibility$coverage
  need_coverage <- c("component", "n_signature", "n_present", "n_filter_passing",
                     "coverage_fraction", "filter_fraction", "median_detect_fraction",
                     "status")
  miss_coverage <- setdiff(need_coverage, names(coverage))
  if (length(miss_coverage)) {
    stop("GeoMx decon-feasibility coverage data missing columns: ",
         paste(miss_coverage, collapse = ", "), call. = FALSE)
  }
  coverage <- coverage[is.finite(coverage$n_signature) &
                         is.finite(coverage$n_present) &
                         is.finite(coverage$n_filter_passing) &
                         is.finite(coverage$coverage_fraction) &
                         is.finite(coverage$filter_fraction) &
                         is.finite(coverage$median_detect_fraction), ,
                       drop = FALSE]
  if (!nrow(coverage)) stop("GeoMx decon-feasibility coverage data has no finite rows",
                            call. = FALSE)
  component_levels <- as.character(coverage$component)
  coverage$component_plot <- factor(as.character(coverage$component),
                                    levels = rev(component_levels))
  coverage$status <- factor(as.character(coverage$status),
                            levels = c("absent", "thin", "covered"))
  coverage$label <- paste0(coverage$n_filter_passing, "/", coverage$n_signature)
  stopifnot(!anyNA(coverage$component_plot), !anyNA(coverage$status))
  coverage_palette <- c(absent = "#B9B1A4", thin = "#C8841C", covered = "#0B7A75")

  coverage_plot <- ggplot2::ggplot(coverage, ggplot2::aes(y = component_plot)) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, xend = coverage_fraction, yend = component_plot),
      colour = "#BFB8AA", linewidth = 0.45, alpha = 0.75) +
    ggplot2::geom_point(
      ggplot2::aes(x = coverage_fraction),
      shape = 21, size = 2.1, fill = "#F8F5ED", colour = "#5E584F", stroke = 0.28) +
    ggplot2::geom_point(
      ggplot2::aes(x = filter_fraction, fill = status, size = n_filter_passing),
      shape = 21, colour = "#20242A", stroke = 0.28, alpha = 0.96) +
    ggplot2::geom_text(
      ggplot2::aes(x = pmin(1, pmax(coverage_fraction, filter_fraction) + 0.055),
                   label = label),
      size = 2.25, colour = "#20242A", hjust = 0) +
    ggplot2::scale_fill_manual(values = coverage_palette, drop = FALSE,
                               name = "coverage") +
    ggplot2::scale_size_area(max_size = 5.2, name = "kept genes") +
    ggplot2::scale_x_continuous(
      limits = c(0, 1.13),
      breaks = seq(0, 1, by = 0.25),
      labels = scales::percent_format(accuracy = 1),
      expand = ggplot2::expansion(mult = c(0.02, 0.01))) +
    ggplot2::labs(
      x = "Marker-set fraction", y = NULL,
      title = "Reference marker coverage",
      subtitle = "Open point = present in WTA; filled point = filterByExpr-kept"
    ) +
    theme_tau(base_size = 8.3) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      axis.text.y = ggplot2::element_text(size = 7.2),
      panel.grid.minor = ggplot2::element_blank()
    )

  aoi <- decon_feasibility$aoi
  need_aoi <- c("slide", "genotype", "x_coord", "y_coord", "library_size",
                "q3_scaled_background", "input_status", "marker_residual_rmse")
  miss_aoi <- setdiff(need_aoi, names(aoi))
  if (length(miss_aoi)) {
    stop("GeoMx decon-feasibility AOI data missing columns: ",
         paste(miss_aoi, collapse = ", "), call. = FALSE)
  }
  aoi <- aoi[is.finite(aoi$x_coord) & is.finite(aoi$y_coord) &
               is.finite(aoi$library_size) &
               is.finite(aoi$q3_scaled_background) &
               is.finite(aoi$marker_residual_rmse), ,
             drop = FALSE]
  if (!nrow(aoi)) stop("GeoMx decon-feasibility AOI data has no finite rows",
                       call. = FALSE)
  aoi$slide <- factor(aoi$slide)
  aoi$genotype <- factor(as.character(aoi$genotype), levels = genotype_levels)
  status_levels <- c("no local blocker", "low-input tail", "background/Q3 tail",
                     "absolute-count blocked")
  aoi$input_status <- factor(as.character(aoi$input_status), levels = status_levels)
  stopifnot(!anyNA(aoi$slide), !anyNA(aoi$genotype), !anyNA(aoi$input_status))
  status_palette <- c(
    `no local blocker` = "#0B7A75",
    `low-input tail` = "#C8841C",
    `background/Q3 tail` = "#7D5CB8",
    `absolute-count blocked` = "#A63A50"
  )
  lib_breaks <- pretty(log10(pmax(aoi$library_size, 1)), n = 3)
  lib_breaks <- lib_breaks[is.finite(lib_breaks)]

  blocker_map <- ggplot2::ggplot(aoi, ggplot2::aes(x_coord, y_coord)) +
    ggplot2::geom_point(
      ggplot2::aes(fill = input_status, colour = genotype,
                   size = log10(pmax(library_size, 1))),
      shape = 21, stroke = 0.36, alpha = 0.95) +
    ggplot2::scale_fill_manual(values = status_palette, drop = FALSE,
                               name = "input status") +
    scale_colour_genotype(name = "genotype") +
    ggplot2::scale_size_continuous(breaks = lib_breaks, name = "log10 library") +
    ggplot2::facet_wrap(ggplot2::vars(slide), ncol = 4) +
    ggplot2::coord_equal() +
    ggplot2::scale_y_reverse() +
    ggplot2::labs(
      x = "ROI coordinate X", y = "ROI coordinate Y",
      title = "AOI blocker map",
      subtitle = "Fill marks decon preconditions; ring marks genotype"
    ) +
    theme_tau(base_size = 8.0) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(order = 1, override.aes = list(size = 3.0)),
      colour = ggplot2::guide_legend(order = 2, override.aes = list(fill = "white",
                                                                    size = 2.8)),
      size = ggplot2::guide_legend(order = 3)
    ) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.box = "vertical",
      axis.text = ggplot2::element_text(size = 5.7),
      axis.title = ggplot2::element_text(size = 7.2),
      strip.text = ggplot2::element_text(size = 7.0),
      panel.spacing = grid::unit(0.32, "lines")
    )

  status_counts <- decon_feasibility$status_counts
  need_counts <- c("genotype", "input_status", "n_aoi")
  miss_counts <- setdiff(need_counts, names(status_counts))
  if (length(miss_counts)) {
    stop("GeoMx decon-feasibility status-count data missing columns: ",
         paste(miss_counts, collapse = ", "), call. = FALSE)
  }
  status_counts <- status_counts[is.finite(status_counts$n_aoi), , drop = FALSE]
  if (!nrow(status_counts)) stop("GeoMx decon-feasibility status counts are empty",
                                call. = FALSE)
  status_counts$genotype <- factor(as.character(status_counts$genotype),
                                   levels = genotype_levels)
  status_counts$input_status <- factor(as.character(status_counts$input_status),
                                       levels = status_levels)
  stopifnot(!anyNA(status_counts$genotype), !anyNA(status_counts$input_status))
  status_counts$label <- ifelse(status_counts$n_aoi > 0,
                                as.character(status_counts$n_aoi), "")
  count_plot <- ggplot2::ggplot(
    status_counts,
    ggplot2::aes(genotype, n_aoi, fill = input_status)
  ) +
    ggplot2::geom_col(width = 0.72, colour = "#F8F5ED", linewidth = 0.2) +
    ggplot2::geom_text(
      ggplot2::aes(label = label),
      position = ggplot2::position_stack(vjust = 0.5),
      size = 2.3, colour = "#20242A") +
    ggplot2::scale_fill_manual(values = status_palette, drop = FALSE,
                               name = "input status") +
    ggplot2::labs(
      x = NULL, y = "AOIs",
      title = "Precondition counts",
      subtitle = "Nuclei sentinels block absolute-count scaling; tails flag input-risk AOIs"
    ) +
    theme_tau(base_size = 8.2) +
    ggplot2::theme(
      legend.position = "bottom",
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1)
    )

  residual_plot <- ggplot2::ggplot(
    aoi,
    ggplot2::aes(input_status, marker_residual_rmse, colour = genotype)
  ) +
    ggplot2::geom_boxplot(
      data = aoi,
      mapping = ggplot2::aes(x = input_status, y = marker_residual_rmse,
                             group = input_status),
      inherit.aes = FALSE,
      width = 0.24, outlier.shape = NA, fill = "#F8F5ED",
      colour = "#20242A", linewidth = 0.26) +
    ggplot2::geom_point(
      position = ggplot2::position_jitter(width = 0.13, height = 0, seed = 42L),
      size = 1.15, alpha = 0.72) +
    scale_colour_genotype(name = "genotype") +
    ggplot2::labs(
      x = NULL, y = "marker residual RMSE",
      title = "Proxy fit residual",
      subtitle = "Marker-coherence RMSE; proxy QC, not a SpatialDecon residual"
    ) +
    theme_tau(base_size = 8.2) +
    ggplot2::theme(
      legend.position = "bottom",
      axis.text.x = ggplot2::element_text(angle = 25, hjust = 1, size = 6.8),
      plot.margin = ggplot2::margin(4, 8, 4, 4)
    )

  prov <- decon_feasibility$provenance
  subtitle <- sprintf(
    "%s/%s AOIs carry nuclei sentinels; %s/%s candidate components have >=2 filter-passing marker genes. %s",
    prov$n_nuclei_sentinel %||% sum(!aoi$nuclei_usable),
    prov$n_aoi %||% nrow(aoi),
    prov$n_score_components %||% sum(coverage$n_filter_passing >= 2L),
    prov$n_components %||% nrow(coverage),
    prov$live_status %||% "blocked diagnostic; no cell-abundance claim"
  )
  subtitle <- paste(strwrap(subtitle, width = 126), collapse = "\n")
  top <- patchwork::wrap_plots(list(coverage_plot, blocker_map),
                               ncol = 2, widths = c(0.92, 1.08))
  bottom <- patchwork::wrap_plots(list(count_plot, residual_plot),
                                  ncol = 2, widths = c(0.86, 1.14))
  patchwork::wrap_plots(list(top, bottom), ncol = 1, heights = c(1.12, 0.88)) +
    patchwork::plot_annotation(title = title, subtitle = subtitle)
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

proteome_pca_plot <- function(proteome, title = "Bulk proteome PCA") {
  stopifnot(is.list(proteome), is.data.frame(proteome$pca))
  pca <- proteome$pca
  need <- c("pc1", "pc2", "pc1_var", "pc2_var", "genotype", "run_index")
  miss <- setdiff(need, names(pca))
  if (length(miss)) stop("proteome PCA data missing columns: ", paste(miss, collapse = ", "),
                         call. = FALSE)
  ggplot2::ggplot(pca, ggplot2::aes(pc1, pc2, colour = genotype)) +
    ggplot2::geom_hline(yintercept = 0, colour = "#D5D0C4", linewidth = 0.25) +
    ggplot2::geom_vline(xintercept = 0, colour = "#D5D0C4", linewidth = 0.25) +
    ggplot2::geom_point(size = 2.4, alpha = 0.95) +
    ggrepel::geom_text_repel(ggplot2::aes(label = run_index), size = tau_report_label_size,
                             max.overlaps = Inf, seed = 42L, min.segment.length = 0,
                             segment.colour = "grey65", box.padding = 0.12,
                             point.padding = 0.08, max.iter = 20000L, max.time = 3,
                             show.legend = FALSE) +
    scale_colour_genotype(name = "genotype") +
    ggplot2::guides(colour = ggplot2::guide_legend(nrow = 2, byrow = TRUE)) +
    ggplot2::labs(
      x = sprintf("PC1 (%.1f%%)", 100 * unique(pca$pc1_var)[[1]]),
      y = sprintf("PC2 (%.1f%%)", 100 * unique(pca$pc2_var)[[1]]),
      title = title) +
    theme_tau() +
    ggplot2::theme(legend.position = "bottom",
                   legend.box = "vertical")
}

phospho_site_heatmap_plot <- function(heatmap, title = "Top phosphosite abundance") {
  stopifnot(is.data.frame(heatmap),
            all(c("sample_label", "site_label_plot", "genotype", "z") %in% names(heatmap)))
  x <- heatmap[is.finite(heatmap$z), , drop = FALSE]
  if (!nrow(x)) stop("phosphosite heatmap has no finite rows", call. = FALSE)
  lim <- max(abs(x$z), na.rm = TRUE)
  lim <- if (is.finite(lim) && lim > 0) lim else 1
  genotype_strip_labels <- c(
    MAPTKI = "MAPTKI",
    P301S = "P301S",
    NLGF_MAPTKI = "NLGF\nMAPTKI",
    NLGF_P301S = "NLGF\nP301S"
  )

  ggplot2::ggplot(x, ggplot2::aes(sample_label, site_label_plot, fill = z)) +
    ggplot2::geom_tile(colour = "#F8F5ED", linewidth = 0.18) +
    ggplot2::facet_grid(. ~ genotype, scales = "free_x", space = "free_x",
                         labeller = ggplot2::labeller(genotype = genotype_strip_labels)) +
    scale_fill_rwb(midpoint = 0, limits = c(-lim, lim), oob = scales::squish,
                   name = "z-score") +
    ggplot2::scale_x_discrete(drop = TRUE) +
    ggplot2::labs(x = NULL, y = NULL, title = title) +
    theme_tau() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1,
                                          size = tau_report_dense_axis_size),
      axis.text.y = ggplot2::element_text(size = tau_report_dense_axis_size,
                                          lineheight = 0.9),
      strip.text.x = ggplot2::element_text(lineheight = 0.9),
      panel.spacing.x = grid::unit(0.16, "lines"),
      legend.position = "bottom",
      legend.box.spacing = grid::unit(0, "pt"),
      legend.margin = ggplot2::margin(t = -24, r = 0, b = -8, l = 0)
    )
}

bulk_modality_context_plot <- function(proteome, phospho,
                                       proteome_title = "Bulk proteome PCA",
                                       phospho_title = "Top phosphosite z-score") {
  stopifnot(is.list(proteome), is.data.frame(proteome$pca),
            is.list(phospho), is.data.frame(phospho$heatmap))
  pca <- proteome_pca_plot(proteome, title = proteome_title)
  heatmap <- phospho_site_heatmap_plot(phospho$heatmap, title = phospho_title)
  widths <- c(pca = 0.78, heatmap = 1.22)
  total_width <- sum(widths)
  top <- patchwork::wrap_plots(
    list(patchwork::plot_spacer(), pca, patchwork::plot_spacer()),
    ncol = 3,
    widths = c((total_width - widths[["pca"]]) / 2, widths[["pca"]],
               (total_width - widths[["pca"]]) / 2)
  )
  bottom <- patchwork::wrap_plots(
    list(patchwork::plot_spacer(), heatmap, patchwork::plot_spacer()),
    ncol = 3,
    widths = c((total_width - widths[["heatmap"]]) / 2, widths[["heatmap"]],
               (total_width - widths[["heatmap"]]) / 2)
  )
  patchwork::wrap_plots(list(top, bottom), ncol = 1, heights = c(1, 1))
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
