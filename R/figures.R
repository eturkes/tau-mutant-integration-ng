# Inline figure-data contracts for the post-report visual-density pass. These
# builders return compact, qmd-ready data frames/lists; they do not draw plots.

.fig_contrasts <- function() {
  c("tau_alone", "nlgf_in_maptki", "nlgf_in_p301s", "tau_in_nlgf", "interaction")
}

.fig_focus_contrasts <- function() {
  c("interaction", "nlgf_in_maptki", "nlgf_in_p301s", "tau_in_nlgf")
}

figure_manifest <- function(chapter = NULL) {
  out <- data.frame(
    figure_id = c(
      "fig-synthesis-evidence-map",
      "fig-microglia-umap-substate",
      "fig-microglia-score-triptych",
      "fig-microglia-unit-composition",
      "fig-microglia-score-distribution",
      "fig-microglia-composition-concordance",
      "fig-microglia-whole-volcano",
      "fig-microglia-substate-audit",
      "fig-microglia-substate-volcano",
      "fig-trajectory-pt-density",
      "fig-trajectory-unit-pt-dam",
      "fig-trajectory-kitagawa-forest",
      "fig-trajectory-audit",
      "fig-mechanism-project-pathway",
      "fig-mechanism-go-dotplot",
      "fig-mechanism-tf-lollipop",
      "fig-mechanism-nfkb-discordance",
      "fig-mechanism-kinase-heatmap",
      "fig-crossmodality-geomx-volcano",
      "fig-crossmodality-geomx-sensitivity",
      "fig-crossmodality-bulk-run-index",
      "fig-crossmodality-phospho-correction",
      "fig-crossmodality-anchor-heatmap",
      "fig-crossmodality-clearance-grid",
      "fig-crossmodality-symbol-matrix",
      "fig-crossmodality-pathway-heatmap"
    ),
    chapter = c(
      "synthesis",
      rep("microglia", 8),
      rep("trajectory", 4),
      rep("mechanism", 5),
      rep("crossmodality", 8)
    ),
    target = c(
      "synthesis_report",
      rep("microglia_figures", 8),
      rep("trajectory_figures", 4),
      rep("mechanism_figures", 5),
      rep("crossmodality_figures", 8)
    ),
    slot = c(
      "evidence_map",
      "umap_by_substate",
      "score_triptych",
      "unit_composition",
      "score_distribution",
      "composition_concordance",
      "whole_de_volcano",
      "substate_fit_audit",
      "substate_de_volcano",
      "pt_density",
      "unit_pt_vs_dam",
      "kitagawa_forest",
      "trajectory_audit",
      "project_pathway_heatmap",
      "go_top_dotplot",
      "tf_lollipop",
      "nfkb_discordance",
      "kinase_heatmap",
      "geomx_volcano",
      "geomx_sensitivity",
      "bulk_run_index",
      "phospho_raw_corrected",
      "anchor_effect_heatmap",
      "clearance_pair_grid",
      "symbol_modality_matrix",
      "pathway_axis_heatmap"
    ),
    stringsAsFactors = FALSE
  )
  stopifnot(!anyDuplicated(out$figure_id), !any(grepl("_", out$figure_id, fixed = TRUE)))
  if (!is.null(chapter)) out <- out[out$chapter %in% chapter, , drop = FALSE]
  rownames(out) <- NULL
  out
}

.fig_require_cols <- function(x, cols, label) {
  if (!is.data.frame(x)) stop(label, " must be a data.frame", call. = FALSE)
  miss <- setdiff(cols, names(x))
  if (length(miss)) {
    stop(label, " missing columns: ", paste(miss, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

.fig_bind <- function(xs) {
  xs <- Filter(function(x) is.data.frame(x) && nrow(x) > 0L, xs)
  if (!length(xs)) return(data.frame())
  out <- do.call(rbind, xs)
  rownames(out) <- NULL
  out
}

.fig_assert_nonempty <- function(x, label) {
  if (!is.data.frame(x) || nrow(x) < 1L) {
    stop(label, " is empty", call. = FALSE)
  }
  invisible(TRUE)
}

.fig_assert_finite <- function(x, cols, label) {
  .fig_require_cols(x, cols, label)
  bad <- cols[!vapply(cols, function(cl) all(is.finite(x[[cl]])), logical(1))]
  if (length(bad)) {
    stop(label, " has non-finite numeric columns: ", paste(bad, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

.fig_midpoints <- function(breaks) {
  (utils::head(breaks, -1L) + utils::tail(breaks, -1L)) / 2
}

.fig_breaks <- function(x, n = 50L, symmetric = FALSE, lower = NULL, upper_q = 0.995) {
  x <- x[is.finite(x)]
  if (!length(x)) stop("cannot bin an all-non-finite vector", call. = FALSE)
  if (symmetric) {
    hi <- as.numeric(stats::quantile(abs(x), upper_q, names = FALSE, na.rm = TRUE))
    hi <- max(hi, 1e-6)
    return(seq(-hi, hi, length.out = n + 1L))
  }
  lo <- if (is.null(lower)) min(x, na.rm = TRUE) else lower
  hi <- as.numeric(stats::quantile(x, upper_q, names = FALSE, na.rm = TRUE))
  hi <- max(hi, lo + 1e-6)
  seq(lo, hi, length.out = n + 1L)
}

.fig_clip <- function(x, breaks) {
  pmin(pmax(x, min(breaks)), max(breaks))
}

.fig_grid_bins <- function(x, y, group, x_bins = 60L, y_bins = 50L,
                           symmetric_x = TRUE, lower_y = 0) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]; group <- group[ok]
  xb <- .fig_breaks(x, n = x_bins, symmetric = symmetric_x)
  yb <- .fig_breaks(y, n = y_bins, symmetric = FALSE, lower = lower_y)
  bx <- cut(.fig_clip(x, xb), breaks = xb, include.lowest = TRUE, labels = FALSE)
  by <- cut(.fig_clip(y, yb), breaks = yb, include.lowest = TRUE, labels = FALSE)
  raw <- data.frame(group = group, bin_x = bx, bin_y = by, n = 1L,
                    stringsAsFactors = FALSE)
  out <- stats::aggregate(n ~ group + bin_x + bin_y, data = raw, FUN = sum)
  xm <- .fig_midpoints(xb); ym <- .fig_midpoints(yb)
  out$x_mid <- xm[out$bin_x]
  out$y_mid <- ym[out$bin_y]
  out$x_min <- xb[out$bin_x]
  out$x_max <- xb[out$bin_x + 1L]
  out$y_min <- yb[out$bin_y]
  out$y_max <- yb[out$bin_y + 1L]
  out[order(out$group, out$bin_x, out$bin_y, method = "radix"), , drop = FALSE]
}

.fig_histogram <- function(df, value_col, group_cols, bins = 50L, lower = NULL) {
  .fig_require_cols(df, c(value_col, group_cols), paste0("histogram(", value_col, ")"))
  x <- df[[value_col]]
  ok <- is.finite(x)
  df <- df[ok, , drop = FALSE]
  x <- x[ok]
  br <- .fig_breaks(x, n = bins, lower = lower, symmetric = FALSE)
  b <- cut(.fig_clip(x, br), breaks = br, include.lowest = TRUE, labels = FALSE)
  raw <- cbind(df[group_cols], data.frame(bin = b, n = 1L, stringsAsFactors = FALSE))
  out <- stats::aggregate(n ~ ., data = raw, FUN = sum)
  mids <- .fig_midpoints(br)
  out$x_mid <- mids[out$bin]
  out$x_min <- br[out$bin]
  out$x_max <- br[out$bin + 1L]
  rownames(out) <- NULL
  out
}

.fig_symbol_labels <- function(tt, symbol_map = NULL) {
  if ("symbol" %in% names(tt)) {
    lab <- as.character(tt$symbol)
  } else if ("gene" %in% names(tt) && is.data.frame(symbol_map)) {
    .fig_require_cols(symbol_map, c("ensembl", "symbol"), "symbol_map")
    map <- stats::setNames(as.character(symbol_map$symbol), as.character(symbol_map$ensembl))
    lab <- unname(map[as.character(tt$gene)])
    lab[is.na(lab) | lab == ""] <- as.character(tt$gene)[is.na(lab) | lab == ""]
  } else if ("site_id" %in% names(tt)) {
    lab <- as.character(tt$site_id)
  } else if ("feature" %in% names(tt)) {
    lab <- as.character(tt$feature)
  } else {
    lab <- paste0("row", seq_len(nrow(tt)))
  }
  lab[is.na(lab) | lab == ""] <- paste0("row", which(is.na(lab) | lab == ""))
  lab
}

.fig_volcano_data <- function(top_list, contrasts = .fig_contrasts(), symbol_map = NULL,
                              alpha = 0.10, n_label = 10L, x_bins = 70L, y_bins = 55L) {
  stopifnot(is.list(top_list), all(contrasts %in% names(top_list)))
  pieces <- lapply(contrasts, function(cn) {
    tt <- top_list[[cn]]
    .fig_require_cols(tt, c("logFC", "P.Value", "adj.P.Val"), paste0("top$", cn))
    lab <- .fig_symbol_labels(tt, symbol_map)
    x <- as.numeric(tt$logFC)
    p <- pmax(as.numeric(tt$P.Value), 1e-300)
    f <- as.numeric(tt$adj.P.Val)
    ok <- is.finite(x) & is.finite(p) & is.finite(f)
    y <- -log10(p[ok])
    sig <- f[ok] < alpha
    bins <- .fig_grid_bins(x[ok], y, group = rep(cn, sum(ok)),
                           x_bins = x_bins, y_bins = y_bins,
                           symmetric_x = TRUE, lower_y = 0)
    names(bins)[names(bins) == "group"] <- "contrast"
    raw_lab <- data.frame(
      contrast = cn,
      label = lab[ok],
      effect = x[ok],
      neg_log10_p = y,
      p_value = p[ok],
      fdr = f[ok],
      significant = sig,
      stringsAsFactors = FALSE
    )
    ord <- order(raw_lab$fdr, raw_lab$p_value, -abs(raw_lab$effect),
                 raw_lab$label, method = "radix", na.last = TRUE)
    labels <- raw_lab[utils::head(ord, n_label), , drop = FALSE]
    counts <- data.frame(
      contrast = cn,
      n_features = nrow(raw_lab),
      n_fdr_0_10 = sum(raw_lab$fdr < alpha),
      n_up_fdr_0_10 = sum(raw_lab$fdr < alpha & raw_lab$effect > 0),
      n_down_fdr_0_10 = sum(raw_lab$fdr < alpha & raw_lab$effect < 0),
      stringsAsFactors = FALSE
    )
    list(bins = bins, labels = labels, counts = counts)
  })
  out <- list(
    bins = .fig_bind(lapply(pieces, `[[`, "bins")),
    labels = .fig_bind(lapply(pieces, `[[`, "labels")),
    counts = .fig_bind(lapply(pieces, `[[`, "counts"))
  )
  .fig_assert_nonempty(out$bins, "volcano bins")
  .fig_assert_nonempty(out$labels, "volcano labels")
  .fig_assert_finite(out$bins, c("x_mid", "y_mid", "n"), "volcano bins")
  .fig_assert_finite(out$labels, c("effect", "neg_log10_p", "p_value", "fdr"),
                     "volcano labels")
  out
}

.fig_unit_composition <- function(counts) {
  mat <- as.matrix(counts)
  stopifnot(nrow(mat) > 0L, ncol(mat) > 0L, all(is.finite(mat)), all(mat >= 0))
  rn <- rownames(mat) %||% paste0("unit", seq_len(nrow(mat)))
  cn <- colnames(mat) %||% paste0("substate", seq_len(ncol(mat)))
  dimnames(mat) <- list(rn, cn)
  out <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
  names(out) <- c("genotype_batch", "substate", "n_cells")
  out$n_cells <- as.numeric(out$n_cells)
  out$genotype <- vapply(out$genotype_batch, function(u) {
    hit <- genotype_levels[startsWith(u, paste0(genotype_levels, "_"))]
    if (length(hit)) hit[1L] else NA_character_
  }, character(1))
  out$batch <- mapply(function(u, g) {
    if (is.na(g)) NA_character_ else substring(u, nchar(g) + 2L)
  }, out$genotype_batch, out$genotype, USE.NAMES = FALSE)
  totals <- stats::aggregate(n_cells ~ genotype_batch, data = out, FUN = sum)
  names(totals)[2L] <- "unit_total"
  out <- merge(out, totals, by = "genotype_batch", all.x = TRUE, sort = FALSE)
  out$proportion <- ifelse(out$unit_total > 0, out$n_cells / out$unit_total, NA_real_)
  out$genotype <- factor(out$genotype, levels = genotype_levels)
  out <- out[order(match(out$genotype, genotype_levels), out$genotype_batch, out$substate,
                   method = "radix", na.last = TRUE), , drop = FALSE]
  rownames(out) <- NULL
  .fig_assert_finite(out, c("n_cells", "unit_total", "proportion"), "unit composition")
  stopifnot(!anyNA(out$genotype))
  out
}

.fig_score_triptych <- function(cf, score_cols) {
  .fig_require_cols(cf, c("umap_1", "umap_2", score_cols), "microglia cell_frame")
  pieces <- lapply(score_cols, function(sc) {
    df <- data.frame(umap_1 = cf$umap_1, umap_2 = cf$umap_2, score = cf[[sc]],
                     score_name = sub("_UCell_z$", "", sc), stringsAsFactors = FALSE)
    ok <- is.finite(df$umap_1) & is.finite(df$umap_2) & is.finite(df$score)
    df <- df[ok, , drop = FALSE]
    xb <- .fig_breaks(df$umap_1, n = 80L, symmetric = FALSE)
    yb <- .fig_breaks(df$umap_2, n = 80L, symmetric = FALSE)
    bx <- cut(.fig_clip(df$umap_1, xb), xb, include.lowest = TRUE, labels = FALSE)
    by <- cut(.fig_clip(df$umap_2, yb), yb, include.lowest = TRUE, labels = FALSE)
    raw <- data.frame(score_name = df$score_name, bin_x = bx, bin_y = by,
                      score_sum = df$score, n = 1L, stringsAsFactors = FALSE)
    ag <- stats::aggregate(cbind(score_sum, n) ~ score_name + bin_x + bin_y,
                           data = raw, FUN = sum)
    ag$score_mean <- ag$score_sum / ag$n
    xm <- .fig_midpoints(xb); ym <- .fig_midpoints(yb)
    ag$umap_1 <- xm[ag$bin_x]
    ag$umap_2 <- ym[ag$bin_y]
    ag[c("score_name", "bin_x", "bin_y", "umap_1", "umap_2", "score_mean", "n")]
  })
  out <- .fig_bind(pieces)
  .fig_assert_finite(out, c("umap_1", "umap_2", "score_mean", "n"), "score triptych")
  out
}

microglia_figure_data <- function(microglia_report, composition_results,
                                  pb_de_microglia, pb_de_substate, symbol_map,
                                  alpha = 0.10) {
  stopifnot(is.list(microglia_report), is.list(composition_results),
            is.list(pb_de_microglia), is.list(pb_de_substate), is.data.frame(symbol_map))
  cf <- microglia_report$cell_frame
  score_cols <- c("Homeostatic_UCell_z", "DAM_UCell_z", "MHC_APC_UCell_z")
  .fig_require_cols(cf, c("umap_1", "umap_2", "genotype", "substate", score_cols),
                    "microglia_report$cell_frame")
  umap <- cf[c("umap_1", "umap_2", "genotype", "substate")]
  score_long <- .fig_bind(lapply(score_cols, function(sc) {
    data.frame(genotype = cf$genotype, substate = cf$substate,
               score_name = sub("_UCell_z$", "", sc), score = cf[[sc]],
               stringsAsFactors = FALSE)
  }))
  score_dist <- .fig_histogram(score_long, "score", c("genotype", "substate", "score_name"),
                               bins = 50L)

  comp <- .fig_unit_composition(composition_results$counts)
  .fig_require_cols(composition_results$concordance,
                    c("contrast", "substate", "dir_logit", "sig_logit",
                      "dir_asin", "sig_asin", "dir_concordant", "sig_concordant", "flag"),
                    "composition_results$concordance")
  concordance <- composition_results$concordance

  whole_volcano <- .fig_volcano_data(pb_de_microglia$top, .fig_contrasts(), symbol_map,
                                     alpha = alpha, n_label = 8L)
  substate_volcano_pieces <- lapply(names(pb_de_substate$per_substate), function(st) {
    res <- pb_de_substate$per_substate[[st]]
    if (!is.list(res) || !identical(res$status, "fit") || !"top" %in% names(res)) return(NULL)
    v <- .fig_volcano_data(res$top, .fig_contrasts(), symbol_map,
                           alpha = alpha, n_label = 5L)
    v$bins$substate <- st
    v$labels$substate <- st
    v$counts$substate <- st
    v
  })
  substate_volcano_pieces <- Filter(Negate(is.null), substate_volcano_pieces)
  substate_volcano <- list(
    bins = .fig_bind(lapply(substate_volcano_pieces, `[[`, "bins")),
    labels = .fig_bind(lapply(substate_volcano_pieces, `[[`, "labels")),
    counts = .fig_bind(lapply(substate_volcano_pieces, `[[`, "counts"))
  )

  cc <- as.matrix(pb_de_substate$cell_counts)
  cell_counts <- as.data.frame(as.table(cc), stringsAsFactors = FALSE)
  names(cell_counts) <- c("substate", "genotype_batch", "n_cells")
  cell_counts$n_cells <- as.numeric(cell_counts$n_cells)
  fit_status <- .fig_bind(lapply(names(pb_de_substate$per_substate), function(st) {
    res <- pb_de_substate$per_substate[[st]]
    data.frame(
      substate = st,
      status = as.character(res$status %||% NA_character_),
      n_cells = as.numeric(res$n_cells %||% NA_real_),
      units = as.numeric(res$units %||% NA_real_),
      reason = as.character(res$reason %||% NA_character_),
      stringsAsFactors = FALSE
    )
  }))

  out <- list(
    manifest = figure_manifest("microglia"),
    umap_by_substate = umap,
    score_triptych = .fig_score_triptych(cf, score_cols),
    unit_composition = comp,
    score_distribution = score_dist,
    composition_concordance = concordance,
    whole_de_volcano = whole_volcano,
    substate_fit_audit = list(cell_counts = cell_counts, fit_status = fit_status),
    substate_de_volcano = substate_volcano,
    provenance = list(
      source_targets = c("microglia_report", "composition_results",
                         "pb_de_microglia", "pb_de_substate", "symbol_map"),
      alpha = alpha,
      contract = "compact inline microglia figure data; no heavy Seurat object"
    )
  )

  .fig_assert_finite(out$umap_by_substate, c("umap_1", "umap_2"), "umap_by_substate")
  .fig_assert_nonempty(out$score_triptych, "score_triptych")
  .fig_assert_nonempty(out$score_distribution, "score_distribution")
  .fig_assert_nonempty(out$composition_concordance, "composition_concordance")
  .fig_assert_nonempty(out$substate_de_volcano$bins, "substate_de_volcano bins")
  out
}

trajectory_figure_data <- function(trajectory_report, composition_results, alpha = 0.10) {
  stopifnot(is.list(trajectory_report), is.list(composition_results))
  cf <- trajectory_report$cell_frame
  .fig_require_cols(cf, c("genotype", "substate", "on_lineage", "pt_raw"), "trajectory cell_frame")
  pt <- cf[cf$on_lineage %in% TRUE & is.finite(cf$pt_raw), , drop = FALSE]
  density <- .fig_histogram(pt, "pt_raw", c("genotype", "substate"), bins = 55L, lower = 0)

  comp <- .fig_unit_composition(composition_results$counts)
  dam <- comp[comp$substate == "DAM", c("genotype_batch", "proportion"), drop = FALSE]
  names(dam)[2L] <- "dam_fraction"
  unit <- merge(trajectory_report$per_unit, dam, by = "genotype_batch", all.x = TRUE, sort = FALSE)
  .fig_require_cols(unit, c("genotype_batch", "genotype", "batch", "mean_pt",
                            "within_homeostatic", "within_dam", "dam_fraction"),
                    "trajectory unit_pt_vs_dam")

  measures <- c("mean_pt", "comp_cf", "progression_cf", "cross",
                "within_homeostatic", "within_dam")
  forest <- .fig_bind(lapply(names(trajectory_report$weighted_top), function(cn) {
    w <- trajectory_report$weighted_top[[cn]]
    .fig_require_cols(w, c("measure", "coef", "ci_l", "ci_r", "p_value"), paste0("weighted_top$", cn))
    w <- w[w$measure %in% measures, , drop = FALSE]
    data.frame(contrast = cn, w, stringsAsFactors = FALSE)
  }))
  forest$significant <- is.finite(forest$p_value) & forest$p_value < alpha

  audit <- list(
    sensitivity = trajectory_report$sensitivity,
    lineage_per_unit = trajectory_report$lineage_per_unit,
    glmm = as.data.frame(trajectory_report$glmm, stringsAsFactors = FALSE),
    provenance = as.data.frame(trajectory_report$provenance[c("concordance_rho",
                                                              "omitted_frac_overall",
                                                              "composition_loading",
                                                              "progression_loading",
                                                              "cross_loading")],
                               stringsAsFactors = FALSE)
  )

  out <- list(
    manifest = figure_manifest("trajectory"),
    pt_density = density,
    unit_pt_vs_dam = unit,
    kitagawa_forest = forest,
    trajectory_audit = audit,
    provenance = list(
      source_targets = c("trajectory_report", "composition_results"),
      alpha = alpha,
      contract = "compact trajectory figure data; composition joined only from compact counts"
    )
  )
  .fig_assert_finite(out$pt_density, c("x_mid", "n"), "pt_density")
  .fig_assert_finite(out$unit_pt_vs_dam, c("mean_pt", "dam_fraction"), "unit_pt_vs_dam")
  .fig_assert_finite(out$kitagawa_forest, c("coef", "ci_l", "ci_r", "p_value"), "kitagawa_forest")
  out
}

mechanism_figure_data <- function(mechanism_report, alpha = 0.10) {
  stopifnot(is.list(mechanism_report))
  project <- mechanism_report$pathway_project
  go_top <- mechanism_report$pathway_go_top
  tf <- mechanism_report$tf_highlights
  nfkb <- mechanism_report$nfkb$table
  kinase <- mechanism_report$kinase$table
  .fig_require_cols(project, c("pathway", "population", "contrast", "NES", "fdr"),
                    "mechanism pathway_project")
  .fig_require_cols(go_top, c("pathway", "population", "contrast", "NES", "fdr", "size"),
                    "mechanism pathway_go_top")
  .fig_require_cols(tf, c("population", "source", "contrast", "score", "fdr", "selection"),
                    "mechanism tf_highlights")
  .fig_require_cols(nfkb, c("population", "contrast", "test", "score", "p_value",
                            "primary_test", "supportive_only", "primary_family_fdr"),
                    "mechanism nfkb")
  .fig_require_cols(kinase, c("source", "contrast", "score", "fdr", "significant",
                              "run_index_score", "run_index_fdr", "run_index_supports",
                              "include_reason"),
                    "mechanism kinase")
  out <- list(
    manifest = figure_manifest("mechanism"),
    project_pathway_heatmap = project,
    go_top_dotplot = go_top,
    tf_lollipop = tf,
    nfkb_discordance = list(table = nfkb, verdict = mechanism_report$nfkb$verdict),
    kinase_heatmap = kinase,
    provenance = list(
      source_targets = "mechanism_report",
      alpha = alpha,
      contract = "compact mechanism figure data from the guarded report bundle"
    )
  )
  .fig_assert_finite(out$project_pathway_heatmap, c("NES", "fdr"), "project_pathway_heatmap")
  .fig_assert_finite(out$go_top_dotplot, c("NES", "fdr", "size"), "go_top_dotplot")
  .fig_assert_finite(out$tf_lollipop, c("score", "fdr"), "tf_lollipop")
  .fig_assert_finite(out$kinase_heatmap, c("score", "fdr", "run_index_score"), "kinase_heatmap")
  out
}

.fig_sensitivity_counts <- function(primary_top, sensitivity, contrasts, alpha = 0.10) {
  stopifnot(is.list(primary_top), is.list(sensitivity), all(contrasts %in% names(primary_top)))
  fits <- names(sensitivity)
  .fig_bind(lapply(fits, function(fit) {
    sens <- sensitivity[[fit]]
    .fig_bind(lapply(contrasts, function(cn) {
      p <- primary_top[[cn]]
      .fig_require_cols(p, c("symbol", "logFC", "adj.P.Val"), paste0("primary$", cn))
      if (!is.list(sens) || !identical(sens$status %||% NA_character_, "fit") ||
          !"top" %in% names(sens) || !cn %in% names(sens$top)) {
        return(data.frame(fit = fit, contrast = cn, status = as.character(sens$status %||% NA_character_),
                          n_primary_sig = sum(p$adj.P.Val < alpha, na.rm = TRUE),
                          n_sensitivity_sig = NA_real_, n_lost = NA_real_,
                          n_gained = NA_real_, n_sign_flip = NA_real_,
                          stringsAsFactors = FALSE))
      }
      s <- sens$top[[cn]]
      .fig_require_cols(s, c("symbol", "logFC", "adj.P.Val"), paste0(fit, "$", cn))
      m <- merge(p[c("symbol", "logFC", "adj.P.Val")],
                 s[c("symbol", "logFC", "adj.P.Val")],
                 by = "symbol", suffixes = c("_primary", "_sensitivity"))
      ps <- is.finite(m$adj.P.Val_primary) & m$adj.P.Val_primary < alpha
      ss <- is.finite(m$adj.P.Val_sensitivity) & m$adj.P.Val_sensitivity < alpha
      data.frame(
        fit = fit, contrast = cn, status = as.character(sens$status %||% "fit"),
        n_primary_sig = sum(ps),
        n_sensitivity_sig = sum(ss),
        n_lost = sum(ps & !ss),
        n_gained = sum(!ps & ss),
        n_sign_flip = sum(ps & ss & sign(m$logFC_primary) != sign(m$logFC_sensitivity)),
        stringsAsFactors = FALSE
      )
    }))
  }))
}

.fig_best_by_id <- function(x, id_col, effect_col = "logFC", fdr_col = "adj.P.Val") {
  .fig_require_cols(x, c(id_col, effect_col, fdr_col), "best_by_id input")
  x <- x[!is.na(x[[id_col]]) & x[[id_col]] != "" &
           is.finite(x[[effect_col]]) & is.finite(x[[fdr_col]]), , drop = FALSE]
  x <- x[order(x[[id_col]], x[[fdr_col]], -abs(x[[effect_col]]),
               method = "radix", na.last = TRUE), , drop = FALSE]
  x[!duplicated(x[[id_col]]), , drop = FALSE]
}

.fig_phospho_correction <- function(raw, corrected, contrasts, alpha = 0.10) {
  stopifnot(is.list(raw), is.list(corrected), all(contrasts %in% names(raw)),
            all(contrasts %in% names(corrected)))
  pieces <- lapply(contrasts, function(cn) {
    r <- .fig_best_by_id(raw[[cn]], "site_id")
    c <- .fig_best_by_id(corrected[[cn]], "site_id")
    m <- merge(r[c("site_id", "logFC", "adj.P.Val")],
               c[c("site_id", "logFC", "adj.P.Val")],
               by = "site_id", suffixes = c("_raw", "_corrected"))
    sig <- (m$adj.P.Val_raw < alpha) | (m$adj.P.Val_corrected < alpha)
    bins <- .fig_grid_bins(m$logFC_raw, m$logFC_corrected, group = rep(cn, nrow(m)),
                           x_bins = 70L, y_bins = 70L, symmetric_x = TRUE, lower_y = NULL)
    names(bins)[names(bins) == "group"] <- "contrast"
    ord <- order(pmin(m$adj.P.Val_raw, m$adj.P.Val_corrected),
                 -abs(m$logFC_raw) - abs(m$logFC_corrected),
                 m$site_id, method = "radix", na.last = TRUE)
    labels <- m[utils::head(ord, 12L), , drop = FALSE]
    data <- data.frame(
      contrast = cn,
      site_id = labels$site_id,
      raw_logFC = labels$logFC_raw,
      corrected_logFC = labels$logFC_corrected,
      raw_fdr = labels$adj.P.Val_raw,
      corrected_fdr = labels$adj.P.Val_corrected,
      stringsAsFactors = FALSE
    )
    counts <- data.frame(
      contrast = cn, n_sites = nrow(m), n_any_fdr_0_10 = sum(sig),
      n_raw_only = sum(m$adj.P.Val_raw < alpha & !(m$adj.P.Val_corrected < alpha)),
      n_corrected_only = sum(!(m$adj.P.Val_raw < alpha) & m$adj.P.Val_corrected < alpha),
      n_both = sum(m$adj.P.Val_raw < alpha & m$adj.P.Val_corrected < alpha),
      stringsAsFactors = FALSE
    )
    list(bins = bins, labels = data, counts = counts)
  })
  out <- list(
    bins = .fig_bind(lapply(pieces, `[[`, "bins")),
    labels = .fig_bind(lapply(pieces, `[[`, "labels")),
    counts = .fig_bind(lapply(pieces, `[[`, "counts"))
  )
  .fig_assert_finite(out$bins, c("x_mid", "y_mid", "n"), "phospho correction bins")
  .fig_assert_finite(out$labels, c("raw_logFC", "corrected_logFC", "raw_fdr", "corrected_fdr"),
                     "phospho correction labels")
  out
}

crossmodality_figure_data <- function(crossmodality_report, geomx_de, bulk_omics_summary,
                                      phospho_de_24m, phospho_corrected_24m,
                                      alpha = 0.10) {
  stopifnot(is.list(crossmodality_report), is.list(geomx_de),
            is.list(bulk_omics_summary), is.list(phospho_de_24m),
            is.list(phospho_corrected_24m))
  contrasts <- .fig_focus_contrasts()
  geomx_volcano <- .fig_volcano_data(geomx_de$primary$top, contrasts,
                                     alpha = alpha, n_label = 8L)
  geomx_sens <- .fig_sensitivity_counts(geomx_de$primary$top, geomx_de$sensitivity,
                                        contrasts, alpha = alpha)
  run_index <- crossmodality_report$bulk$run_index
  run_index$loss_fraction <- ifelse(run_index$n_primary_sig > 0,
                                    run_index$n_lost_or_flipped / run_index$n_primary_sig,
                                    0)
  phospho_correction <- .fig_phospho_correction(phospho_de_24m$top,
                                                phospho_corrected_24m$top,
                                                contrasts, alpha = alpha)

  out <- list(
    manifest = figure_manifest("crossmodality"),
    geomx_volcano = geomx_volcano,
    geomx_sensitivity = geomx_sens,
    bulk_run_index = run_index,
    phospho_raw_corrected = phospho_correction,
    anchor_effect_heatmap = crossmodality_report$bulk$anchor_rows,
    clearance_pair_grid = crossmodality_report$clearance$pair_support,
    symbol_modality_matrix = crossmodality_report$divergence$axis_symbols,
    pathway_axis_heatmap = crossmodality_report$pathway$axis_summary,
    provenance = list(
      source_targets = c("crossmodality_report", "geomx_de",
                         "bulk_omics_summary", "phospho_de_24m",
                         "phospho_corrected_24m"),
      alpha = alpha,
      contract = "compact cross-modality figure data; heavy top tables reduced to bins"
    )
  )
  .fig_assert_finite(out$geomx_sensitivity,
                     c("n_primary_sig", "n_sensitivity_sig", "n_lost", "n_gained", "n_sign_flip"),
                     "geomx_sensitivity")
  .fig_assert_finite(out$bulk_run_index, c("n_primary_sig", "n_lost_or_flipped", "loss_fraction"),
                     "bulk_run_index")
  .fig_assert_finite(out$anchor_effect_heatmap, c("logFC", "t", "P.Value", "adj.P.Val"),
                     "anchor_effect_heatmap")
  .fig_assert_finite(out$symbol_modality_matrix,
                     c("n_modalities_present", "n_modalities_sig", "min_fdr", "rank_score"),
                     "symbol_modality_matrix")
  .fig_assert_finite(out$pathway_axis_heatmap,
                     c("n_modalities_present", "n_modalities_sig", "rank_score"),
                     "pathway_axis_heatmap")
  out
}
