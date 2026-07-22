# Inline figure-data contracts for the post-report visual-density pass. These
# builders return compact, qmd-ready data frames/lists; they do not draw plots.

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

microglia_figure_data <- function(microglia_report) {
  stopifnot(is.list(microglia_report))
  cf <- microglia_report$cell_frame
  score_cols <- c("Homeostatic_UCell_z", "DAM_UCell_z", "MHC_APC_UCell_z")
  .fig_require_cols(cf, c("umap_1", "umap_2", "genotype", "subpopulation", score_cols),
                    "microglia_report$cell_frame")
  .fig_require_cols(microglia_report$unit_composition,
                    c("genotype_batch", "subpopulation", "n_cells", "genotype", "batch",
                      "unit_total", "proportion"),
                    "microglia_report$unit_composition")
  umap <- cf[c("umap_1", "umap_2", "genotype", "subpopulation")]
  comp <- microglia_report$unit_composition

  out <- list(
    umap_by_subpopulation = umap,
    unit_composition = comp,
    provenance = list(
      source_targets = "microglia_report",
      contract = "compact data for rendered microglia figures only"
    )
  )

  .fig_assert_finite(out$umap_by_subpopulation, c("umap_1", "umap_2"), "umap_by_subpopulation")
  .fig_assert_finite(out$unit_composition, c("n_cells", "unit_total", "proportion"),
                     "unit_composition")
  out
}

trajectory_figure_data <- function(trajectory_report) {
  stopifnot(is.list(trajectory_report))
  cf <- trajectory_report$cell_frame
  .fig_require_cols(cf, c("genotype", "subpopulation", "on_lineage", "pt_raw"),
                    "trajectory cell_frame")
  pt <- cf[cf$on_lineage %in% TRUE & is.finite(cf$pt_raw), , drop = FALSE]
  density <- .fig_histogram(pt, "pt_raw", c("genotype", "subpopulation"), bins = 55L, lower = 0)

  out <- list(
    pt_density = density,
    provenance = list(
      source_targets = "trajectory_report",
      contract = "compact data for rendered trajectory figures only"
    )
  )
  .fig_assert_finite(out$pt_density, c("x_mid", "n"), "pt_density")
  out
}

.fig_nonblank <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x)] <- ""
  x
}

.fig_first_nonblank <- function(...) {
  xs <- list(...)
  n <- max(vapply(xs, length, integer(1)))
  out <- rep("", n)
  for (x in xs) {
    z <- rep(.fig_nonblank(x), length.out = n)
    take <- out == "" & z != ""
    out[take] <- z[take]
  }
  out[out == ""] <- paste0("feature_", which(out == ""))
  out
}

.fig_contrast_rank_points <- function(tt, labels, feature_col = "feature") {
  .fig_require_cols(tt, c(feature_col, "logFC", "P.Value", "adj.P.Val"),
                    "modality contrast points")
  stopifnot(length(labels) == nrow(tt))
  out <- data.frame(
    feature = as.character(tt[[feature_col]]),
    label = .fig_first_nonblank(labels, tt[[feature_col]]),
    effect = as.numeric(tt$logFC),
    p_value = as.numeric(tt$P.Value),
    fdr = as.numeric(tt$adj.P.Val),
    stringsAsFactors = FALSE
  )
  out <- out[is.finite(out$effect) & is.finite(out$p_value) & is.finite(out$fdr), ,
             drop = FALSE]
  .fig_assert_nonempty(out, "modality contrast points")
  out$neg_log10_p <- -log10(pmax(out$p_value, 1e-300))
  out$rank_score <- abs(out$effect) * out$neg_log10_p
  .fig_assert_finite(out, c("effect", "p_value", "fdr", "neg_log10_p", "rank_score"),
                     "modality contrast points")
  rownames(out) <- NULL
  out
}

.fig_bulk_pca_data <- function(mat, meta, label = "bulk matrix") {
  mat <- as.matrix(mat)
  stopifnot(is.numeric(mat), nrow(mat) >= 2L, ncol(mat) >= 3L, is.data.frame(meta),
            identical(colnames(mat), rownames(meta)))
  .fig_require_cols(meta, c("genotype", "run_index"), paste0(label, " metadata"))
  keep <- rowSums(is.finite(mat)) >= 2L
  mat <- mat[keep, , drop = FALSE]
  if (nrow(mat) < 2L) stop(label, " has too few rows for PCA", call. = FALSE)
  med <- apply(mat, 1L, stats::median, na.rm = TRUE)
  miss <- which(!is.finite(mat), arr.ind = TRUE)
  if (nrow(miss)) mat[miss] <- med[miss[, "row"]]
  row_sd <- apply(mat, 1L, stats::sd)
  mat <- mat[is.finite(row_sd) & row_sd > 0, , drop = FALSE]
  if (nrow(mat) < 2L) stop(label, " has too few variable rows for PCA", call. = FALSE)
  pc <- stats::prcomp(t(mat), center = TRUE, scale. = TRUE)
  var <- pc$sdev^2 / sum(pc$sdev^2)
  out <- data.frame(
    sample_id = rownames(meta),
    genotype = factor(as.character(meta$genotype), levels = genotype_levels),
    run_index = as.integer(meta$run_index),
    pc1 = as.numeric(pc$x[, 1L]),
    pc2 = as.numeric(pc$x[, 2L]),
    pc1_var = var[[1L]],
    pc2_var = var[[2L]],
    stringsAsFactors = FALSE
  )
  .fig_assert_finite(out, c("run_index", "pc1", "pc2", "pc1_var", "pc2_var"),
                     paste0(label, " PCA"))
  stopifnot(!anyNA(out$genotype))
  rownames(out) <- NULL
  out
}

.fig_matrix_heatmap_data <- function(mat, meta, feature_ids, labels) {
  mat <- as.matrix(mat)
  stopifnot(is.numeric(mat), is.data.frame(meta), identical(colnames(mat), rownames(meta)),
            length(feature_ids) == length(labels), all(feature_ids %in% rownames(mat)))
  .fig_require_cols(meta, c("genotype", "run_index"), "heatmap metadata")
  sample_order <- rownames(meta)[order(match(as.character(meta$genotype), genotype_levels),
                                       as.integer(meta$run_index), rownames(meta),
                                       method = "radix")]
  z <- mat[feature_ids, sample_order, drop = FALSE]
  med <- apply(z, 1L, stats::median, na.rm = TRUE)
  miss <- which(!is.finite(z), arr.ind = TRUE)
  if (nrow(miss)) z[miss] <- med[miss[, "row"]]
  row_mean <- rowMeans(z)
  row_sd <- apply(z, 1L, stats::sd)
  row_sd[!is.finite(row_sd) | row_sd <= 0] <- 1
  z <- sweep(sweep(z, 1L, row_mean, "-"), 1L, row_sd, "/")
  label_plot <- make.unique(.fig_first_nonblank(labels, feature_ids), sep = " ")
  names(label_plot) <- feature_ids
  long <- as.data.frame(as.table(z), stringsAsFactors = FALSE)
  names(long) <- c("feature", "sample_id", "z")
  sm <- meta[as.character(long$sample_id), , drop = FALSE]
  long$genotype <- factor(as.character(sm$genotype), levels = genotype_levels)
  long$run_index <- as.integer(sm$run_index)
  long$sample_label <- factor(as.character(long$run_index),
                              levels = as.character(meta$run_index[match(sample_order, rownames(meta))]))
  long$site_label_plot <- factor(label_plot[as.character(long$feature)],
                                 levels = rev(label_plot[feature_ids]))
  .fig_assert_finite(long, c("z", "run_index"), "modality heatmap")
  stopifnot(!anyNA(long$genotype), !anyNA(long$site_label_plot))
  rownames(long) <- NULL
  long
}

.fig_matrix_profile_key <- function(mat) {
  mat <- as.matrix(mat)
  stopifnot(is.numeric(mat))
  vapply(seq_len(nrow(mat)), function(i) {
    x <- mat[i, ]
    paste(ifelse(is.na(x), "<NA>", sprintf("%.17g", x)), collapse = "\r")
  }, character(1))
}

# Historical helper token retained; it produces feature labels only and never an assay-name label.
.fig_proteome_labels <- function(tt) {
  .fig_first_nonblank(if ("gene_first" %in% names(tt)) tt$gene_first else "",
                      if ("gene_symbols" %in% names(tt)) tt$gene_symbols else "",
                      tt$feature)
}

.fig_phosphosite_labels <- function(tt) {
  .fig_first_nonblank(if ("site_id" %in% names(tt)) tt$site_id else "",
                      if ("gene" %in% names(tt)) tt$gene else "",
                      tt$feature)
}

# Historical `proteome` symbols/keys stay stable, but this descriptor is the TiO2 phospho-PTM
# report summed to protein groups rather than an independent global-proteome assay.
proteome_modality_descriptor <- function(proteome_de_24m) {
  stopifnot(is.list(proteome_de_24m),
            is.matrix(proteome_de_24m$matrix), is.data.frame(proteome_de_24m$meta))
  list(
    # `proteome` is an internal diagnostic tag retained for compatibility, not rendered text.
    pca = .fig_bulk_pca_data(proteome_de_24m$matrix, proteome_de_24m$meta,
                             label = "proteome"),
    provenance = list(
      n_features = nrow(proteome_de_24m$matrix),
      n_samples = ncol(proteome_de_24m$matrix),
      display = "sample PCA of median-normalised TiO2 phospho-PTM protein-group-sum intensities"
    )
  )
}

phospho_modality_descriptor <- function(phospho_de_24m,
                                        contrast = "nlgf_in_p301s",
                                        n_heatmap = 20L,
                                        heatmap_exclude_genes = c("Plcb1", "Arhgef7")) {
  stopifnot(is.list(phospho_de_24m), contrast %in% names(phospho_de_24m$top),
            is.matrix(phospho_de_24m$matrix), is.data.frame(phospho_de_24m$meta))
  tt <- phospho_de_24m$top[[contrast]]
  labels <- .fig_phosphosite_labels(tt)
  contrast_points <- .fig_contrast_rank_points(tt, labels, feature_col = "feature")
  ranked <- contrast_points[order(contrast_points$fdr, -contrast_points$rank_score,
                                  contrast_points$label, contrast_points$feature,
                                  method = "radix"), , drop = FALSE]
  heatmap_exclude_genes <- unique(as.character(heatmap_exclude_genes))
  heatmap_exclude_genes <- heatmap_exclude_genes[!is.na(heatmap_exclude_genes) &
                                                   nzchar(heatmap_exclude_genes)]
  if (length(heatmap_exclude_genes)) {
    .fig_require_cols(tt, c("feature", "gene"), "phospho top table")
    ranked_gene <- as.character(tt$gene[match(ranked$feature, tt$feature)])
    ranked <- ranked[is.na(ranked_gene) | !(ranked_gene %in% heatmap_exclude_genes), ,
                     drop = FALSE]
  }
  ranked <- ranked[ranked$feature %in% rownames(phospho_de_24m$matrix), , drop = FALSE]
  n_heatmap_candidates <- nrow(ranked)
  heatmap_effect_sign <- sign(ranked$effect[which(ranked$effect != 0)[1L]])
  if (!is.finite(heatmap_effect_sign) || heatmap_effect_sign == 0) {
    stop("phosphosite heatmap candidates have no non-zero effect direction", call. = FALSE)
  }
  ranked <- ranked[sign(ranked$effect) == heatmap_effect_sign, , drop = FALSE]
  if (!nrow(ranked)) stop("no phosphosite heatmap candidates remain after direction filter", call. = FALSE)
  n_heatmap_opposite_direction_candidates <- n_heatmap_candidates - nrow(ranked)
  n_heatmap_direction_candidates <- nrow(ranked)
  profile_key <- .fig_matrix_profile_key(phospho_de_24m$matrix[ranked$feature, , drop = FALSE])
  ranked <- ranked[!duplicated(profile_key), , drop = FALSE]
  heat_features <- utils::head(ranked$feature, as.integer(n_heatmap))
  if (!length(heat_features)) stop("no phosphosite heatmap features overlap matrix", call. = FALSE)
  heat_labels <- ranked$label[match(heat_features, ranked$feature)]
  list(
    heatmap = .fig_matrix_heatmap_data(phospho_de_24m$matrix, phospho_de_24m$meta,
                                       heat_features, heat_labels),
    provenance = list(
      contrast = contrast,
      n_features = nrow(phospho_de_24m$matrix),
      n_samples = ncol(phospho_de_24m$matrix),
      n_heatmap = length(heat_features),
      heatmap_exclude_genes = heatmap_exclude_genes,
      n_heatmap_duplicate_candidates = n_heatmap_direction_candidates - nrow(ranked),
      n_heatmap_opposite_direction_candidates = n_heatmap_opposite_direction_candidates,
      heatmap_effect_direction = if (heatmap_effect_sign > 0) "positive" else "negative",
      heatmap_deduplicate = "exact log2 median-normalised input profiles; first ranked representative retained",
      display = "z-scored top-site abundance heatmap for the mutant-tau amyloid contrast"
    )
  )
}

# Per-modality amyloid-response logFC pairs for fig-modality-amyloid-effect (one scatter per
# method), plus compact functional-category scores for shared-cutoff off-diagonal features. y = logFC of
# `nlgf_in_maptki` (amyloid effect on the WT-humanized-tau (MAPTKI) background), x = logFC of
# `nlgf_in_p301s` (amyloid effect on the mutant-tau / P301S background). Both per-contrast
# topTables come from ONE fit per modality (identical feature rows), aligned by the modality's
# feature key. Compact per-modality frames {feature, label, gene_symbols, x, y, interaction}
# -> the qmd reads this small target, never a heavy DE object. Feature keys / display labels
# differ by assay: snRNAseq = Ensembl gene (mapped to symbol), GeoMx = gene symbol, the
# historical Proteome key = TiO2 phospho-PTM protein-group sum (gene_first label, all group
# symbols for scoring), and Phospho = parent-protein mean of TiO2 phosphosite rows.
modality_logfc_scatter_data <- function(pb_de_microglia, symbol_map, geomx_de,
                                         proteome_de_24m, phospho_de_24m,
                                         y_contrast = "nlgf_in_maptki",
                                         x_contrast = "nlgf_in_p301s",
                                         group_gene_sets = NULL,
                                         offdiag_cutoff = 3.5,
                                         group_min_genes = 1L,
                                         group_max_groups = 10L,
                                         group_min_abs_delta = 0.5) {
  stopifnot(is.list(pb_de_microglia), is.data.frame(symbol_map), is.list(geomx_de),
            is.list(proteome_de_24m), is.list(phospho_de_24m),
            is.numeric(offdiag_cutoff), length(offdiag_cutoff) == 1L,
            is.finite(offdiag_cutoff), offdiag_cutoff > 0,
            is.numeric(group_min_abs_delta), length(group_min_abs_delta) == 1L,
            is.finite(group_min_abs_delta), group_min_abs_delta >= 0)
  offdiag_cutoff <- as.numeric(offdiag_cutoff)

  pair <- function(top_list, key_col, label_fun, gene_fun, modality) {
    stopifnot(is.list(top_list), all(c(y_contrast, x_contrast) %in% names(top_list)))
    ty <- top_list[[y_contrast]]; tx <- top_list[[x_contrast]]
    .fig_require_cols(ty, c(key_col, "logFC"), paste0(modality, " top$", y_contrast))
    .fig_require_cols(tx, c(key_col, "logFC"), paste0(modality, " top$", x_contrast))
    ky <- as.character(ty[[key_col]]); kx <- as.character(tx[[key_col]])
    # one fit -> both contrasts share an identical, unique feature set. Assert BOTH keys unique + equal
    # length before match(): setequal() ignores multiplicity, so a duplicated kx would first-match silently.
    stopifnot(anyDuplicated(ky) == 0L, anyDuplicated(kx) == 0L,
              length(ky) == length(kx), setequal(ky, kx))
    idx <- match(ky, kx)
    df <- data.frame(
      feature = ky,
      label = as.character(label_fun(ty)),
      gene_symbols = vapply(gene_fun(ty), function(z) {
        paste(.fig_gene_tokens(z), collapse = ";")
      }, character(1), USE.NAMES = FALSE),
      y = as.numeric(ty$logFC),                                # nlgf_in_maptki (amyloid | MAPTKI)
      x = as.numeric(tx$logFC)[idx],                           # nlgf_in_p301s  (amyloid | P301S)
      stringsAsFactors = FALSE
    )
    df <- df[is.finite(df$x) & is.finite(df$y), , drop = FALSE]
    blank <- is.na(df$label) | df$label == ""
    df$label[blank] <- df$feature[blank]
    df$interaction <- df$x - df$y
    df$abs_interaction <- abs(df$interaction)
    .fig_assert_nonempty(df, paste0(modality, " logFC pairs"))
    .fig_assert_finite(df, c("x", "y", "interaction", "abs_interaction"),
                       paste0(modality, " logFC pairs"))
    rownames(df) <- NULL
    df
  }

  symbol_label <- function(tt) .fig_symbol_labels(tt, symbol_map)
  symbol_gene <- function(tt) symbol_label(tt)
  gene_first_label <- function(tt) {
    lab <- if ("gene_first" %in% names(tt)) as.character(tt$gene_first) else rep(NA_character_, nrow(tt))
    bad <- is.na(lab) | lab == ""
    lab[bad] <- as.character(tt$feature)[bad]
    lab
  }
  protein_group_genes <- function(tt) {
    lab <- if ("gene_symbols" %in% names(tt)) as.character(tt$gene_symbols) else gene_first_label(tt)
    bad <- is.na(lab) | lab == ""
    lab[bad] <- gene_first_label(tt)[bad]
    lab
  }
  site_id_label <- function(tt) {
    lab <- if ("site_id" %in% names(tt)) as.character(tt$site_id) else rep(NA_character_, nrow(tt))
    bad <- is.na(lab) | lab == ""
    lab[bad] <- as.character(tt$feature)[bad]
    lab
  }
  phosphosite_parent_gene <- function(x) {
    lab <- site_id_label(x)
    lab <- sub("_[A-Za-z][0-9].*$", "", lab, perl = TRUE)
    lab
  }
  phospho_gene <- function(tt) {
    lab <- if ("gene" %in% names(tt)) as.character(tt$gene) else rep(NA_character_, nrow(tt))
    fallback <- phosphosite_parent_gene(tt)
    bad <- is.na(lab) | lab == "" | grepl("_[A-Za-z][0-9].*$", lab, perl = TRUE)
    lab[bad] <- fallback[bad]
    lab
  }
  phospho_parent_label <- function(gene_symbols, labels, features) {
    vapply(seq_along(gene_symbols), function(i) {
      for (z in list(gene_symbols[i], labels[i], features[i])) {
        tok <- .fig_gene_tokens(z)
        if (length(tok)) return(tok[[1]])
      }
      fallback <- as.character(features[i])
      if (!is.na(fallback) && fallback != "") fallback else paste0("phospho_feature_", i)
    }, character(1), USE.NAMES = FALSE)
  }
  collapse_phospho_by_protein <- function(df) {
    .fig_require_cols(df, c("feature", "label", "gene_symbols", "x", "y"),
                      "phospho logFC pairs")
    parent <- phospho_parent_label(df$gene_symbols, df$label, df$feature)
    parent <- trimws(parent)
    parent[is.na(parent) | parent == ""] <- paste0("phospho_feature_",
                                                   which(is.na(parent) | parent == ""))
    idx <- split(seq_len(nrow(df)), parent, drop = TRUE)
    parents <- sort(names(idx))
    out <- .fig_bind(lapply(parents, function(p) {
      z <- df[idx[[p]], , drop = FALSE]
      y <- mean(z$y)
      x <- mean(z$x)
      data.frame(
        feature = paste0("phospho_protein:", p),
        label = p,
        gene_symbols = p,
        y = y,
        x = x,
        interaction = x - y,
        abs_interaction = abs(x - y),
        n_phosphosite = nrow(z),
        stringsAsFactors = FALSE
      )
    }))
    .fig_assert_nonempty(out, "phospho parent-protein logFC pairs")
    .fig_assert_finite(out, c("x", "y", "interaction", "abs_interaction", "n_phosphosite"),
                       "phospho parent-protein logFC pairs")
    attr(out, "phospho_site_n") <- nrow(df)
    attr(out, "phospho_parent_collapse") <- "mean x/y by best-fit parent gene"
    rownames(out) <- NULL
    out
  }

  panels <- list(
    snRNAseq = list(
      title = "snRNAseq microglia",
      data  = pair(pb_de_microglia$top, "gene",
                   symbol_label, symbol_gene, "snRNAseq")),
    GeoMx = list(
      title = "GeoMx microglia",
      data  = pair(geomx_de$primary$top, "symbol",
                   function(tt) as.character(tt$symbol),
                   function(tt) as.character(tt$symbol), "GeoMx")),
    # Historical list keys stay stable; displayed titles distinguish two levels of one TiO2 assay.
    Proteome = list(
      title = "Bulk phospho (protein-group)",
      data  = pair(proteome_de_24m$top, "feature", gene_first_label,
                   protein_group_genes, "proteome")),
    Phospho = list(
      title = "Bulk phospho (site)",
      data  = collapse_phospho_by_protein(pair(phospho_de_24m$top, "feature", site_id_label,
                                               phospho_gene, "phospho")))
  )
  order <- c("snRNAseq", "GeoMx", "Proteome", "Phospho")
  offdiag_cutoff_by_modality <- stats::setNames(rep(offdiag_cutoff, length(order)), order)
  offdiag_cutoff_source <- sprintf("shared |x-y| >= %.1f", offdiag_cutoff)
  offdiag_thresholds <- .fig_bind(lapply(order, function(m) {
    d <- panels[[m]]$data
    dist <- abs(d$y - d$x)
    dist <- dist[is.finite(dist) & dist > 0]
    data.frame(
      modality = m,
      cutoff = offdiag_cutoff,
      max_distance = if (length(dist)) max(dist) else NA_real_,
      n_at_cutoff = sum(dist >= offdiag_cutoff),
      n_features = length(dist),
      stringsAsFactors = FALSE
    )
  }))
  for (m in order) {
    d <- panels[[m]]$data
    attr(d, "offdiag_cutoff") <- offdiag_cutoff
    attr(d, "offdiag_cutoff_source") <- offdiag_cutoff_source
    panels[[m]]$data <- d
  }
  groups <- modality_offdiag_group_score_data(
    list(panels = panels, order = order),
    group_sets = group_gene_sets,
    min_genes = group_min_genes,
    max_groups = group_max_groups,
    min_abs_delta = group_min_abs_delta
  )
  descriptive <- list(
    GeoMx = list(sample_heatmap = geomx_de$sample_heatmap),
    # Historical payload key retained; this is the TiO2 protein-group-sum view.
    Proteome = proteome_modality_descriptor(proteome_de_24m),
    Phospho = phospho_modality_descriptor(phospho_de_24m)
  )
  stopifnot(is.list(descriptive$GeoMx),
            is.list(descriptive$GeoMx$sample_heatmap),
            is.data.frame(descriptive$GeoMx$sample_heatmap$heatmap),
            is.list(descriptive$Proteome), is.data.frame(descriptive$Proteome$pca),
            is.list(descriptive$Phospho), is.data.frame(descriptive$Phospho$heatmap))

  list(
    panels = panels,
    order = order,
    groups = groups,
    descriptive = descriptive,
    provenance = list(
      y_contrast = y_contrast,
      x_contrast = x_contrast,
      y_meaning = "amyloid effect on the WT-humanized-tau (MAPTKI) background",
      x_meaning = "amyloid effect on the mutant-tau (P301S) background",
      interaction = "x - y is the tau-by-amyloid interaction contrast per feature; |x - y| ranks off-diagonal distance",
      offdiag_cutoff = offdiag_cutoff_by_modality,
      offdiag_thresholds = offdiag_thresholds,
      offdiag_cutoff_source = offdiag_cutoff_source,
      scatter_label_rule = "all points passing the stored |x-y| cutoff",
      group_min_abs_delta = group_min_abs_delta,
      n_features = vapply(order, function(m) nrow(panels[[m]]$data), integer(1)),
      phospho_site_features = attr(panels$Phospho$data, "phospho_site_n", exact = TRUE),
      phospho_parent_proteins = nrow(panels$Phospho$data),
      phospho_scatter = "phosphosite rows collapsed to best-fit parent proteins; x/y are arithmetic means across finite sites per protein",
      feature_key = c(snRNAseq = "Ensembl gene (symbol label)", GeoMx = "gene symbol",
                      Proteome = "TiO2 phospho-PTM protein-group sum (gene_first label)",
                      Phospho = "TiO2 parent-protein mean of phosphosite rows (best-fit gene label)"),
      source_targets = c("pb_de_microglia", "symbol_map", "geomx_de",
                         "proteome_de_24m", "phospho_de_24m"),
      contract = "compact per-modality amyloid-response logFC pairs + shared-cutoff off-diagonal functional-category aggregate scores + modality-native descriptive figure data: GeoMx sample heatmap, TiO2 phospho protein-group PCA, TiO2 phosphosite heatmap; no heavy DE object"
    )
  )
}

.fig_gene_tokens <- function(x) {
  if (!length(x)) return(character())
  raw <- unlist(strsplit(paste(as.character(x), collapse = ";"), "[;,]", perl = TRUE),
                use.names = FALSE)
  raw <- trimws(raw)
  raw <- raw[!is.na(raw) & raw != ""]
  raw[raw %in% c("hMapt", "hMAPT", "MAPT")] <- "Mapt"
  unique(raw)
}

.fig_pathway_label <- function(x) {
  lab <- gsub("^(HALLMARK|GOBP|REACTOME)_", "", as.character(x))
  lab <- gsub("_", " ", lab, fixed = TRUE)
  lab <- tools::toTitleCase(tolower(lab))
  lab <- gsub("\\bDna\\b", "DNA", lab)
  lab <- gsub("\\bRna\\b", "RNA", lab)
  lab <- gsub("\\bUv\\b", "UV", lab)
  lab <- gsub("\\bGtpase\\b", "GTPase", lab)
  lab <- gsub("\\bIl([0-9]+)\\b", "IL\\1", lab)
  lab <- gsub("\\bTnfa\\b", "TNFA", lab)
  lab <- gsub("\\bNfkb\\b", "NF-kB", lab)
  lab <- gsub("Signaling", "signalling", lab, fixed = TRUE)
  lab
}

.fig_functional_group_label <- function(x) {
  lab <- as.character(x)
  complement <- grepl("^Complement\\s*/\\s*Antigen\\s+Presentation\\b", lab)
  lab[complement] <- "Complement / MHC"
  phagocytosis <- lab == "Phagocytosis"
  chemotaxis <- lab == "Chemotaxis"
  leukocyte_adhesion <- lab == "Leukocyte Adhesion"
  lab[leukocyte_adhesion] <- "Cell-Cell Adhesion"
  leukocyte_migration <- lab == "Leukocyte Migration"
  matrix <- lab == "Matrix"
  lab[matrix] <- "Extracellular Matrix"
  immune_inflammatory <- grepl("^Immune\\s*/\\s*Inflammatory\\b", lab)
  lab[immune_inflammatory] <- "Immune / Inflammatory"
  keep_slash <- complement | phagocytosis | chemotaxis | leukocyte_adhesion |
    leukocyte_migration | matrix | immune_inflammatory
  slash <- grepl("/", lab, fixed = TRUE) & !keep_slash
  lab[slash] <- trimws(sub("\\s*/.*$", "", lab[slash]))
  lab[lab == "Synapse"] <- "Synaptic"
  lab
}

.fig_default_pathway_sets <- function(collection = "M5", subcollection = "GO:BP",
                                      min_size = 10L, max_size = 500L) {
  x <- suppressMessages(msigdbr::msigdbr(db_species = "MM", species = "Mus musculus",
                                         collection = collection,
                                         subcollection = subcollection))
  .fig_require_cols(x, c("gene_symbol", "gs_name"), "msigdbr pathway gene sets")
  sets <- split(as.character(x$gene_symbol), as.character(x$gs_name))
  sets <- lapply(sets, .fig_gene_tokens)
  n <- vapply(sets, length, integer(1))
  sets[n >= min_size & n <= max_size]
}

.fig_default_functional_groups <- function(pathway_sets = NULL) {
  pathway_sets <- pathway_sets %||% .fig_default_pathway_sets()
  stopifnot(is.list(pathway_sets), length(pathway_sets) >= 1L, !is.null(names(pathway_sets)),
            !any(names(pathway_sets) == ""))
  pathway_sets <- lapply(pathway_sets, .fig_gene_tokens)
  pathway_sets <- pathway_sets[vapply(pathway_sets, length, integer(1)) > 0L]
  stopifnot(length(pathway_sets) >= 1L)

  role_patterns <- list(
    `Complement / antigen presentation` =
      c("COMPLEMENT", "OPSONIZATION", "ANTIGEN", "MHC"),
    Phagocytosis =
      c("PHAGOCYTOSIS"),
    Chemotaxis =
      c("CHEMOTAXIS", "GRANULOCYTE_MIGRATION", "NEUTROPHIL_MIGRATION"),
    `Lipid handling / sterol biology` =
      c("LIPID", "STEROL", "CHOLESTEROL", "LIPOPROTEIN", "FATTY_ACID"),
    `Endolysosome / vesicle traffic` =
      c("LYSOSOME", "LYSOSOMAL", "ENDOSOME", "ENDOCYTOSIS", "VESICLE",
        "VACUOLE", "AUTOPHAGY", "PHAGOSOME"),
    `Synapse / neuronal signalling` =
      c("SYNAP", "NEURON", "AXON", "DENDRITE", "NEUROTRANSMITTER",
        "ACTION_POTENTIAL", "MEMBRANE_POTENTIAL"),
    `Leukocyte adhesion` =
      c("LEUKOCYTE_CELL_CELL_ADHESION", "LEUKOCYTE_ADHESION"),
    `Leukocyte migration` =
      c("LEUKOCYTE_MIGRATION", "DENDRITIC_CELL_MIGRATION",
        "MONONUCLEAR_CELL_MIGRATION", "MYELOID_LEUKOCYTE_MIGRATION"),
    Matrix =
      c("EXTRACELLULAR_MATRIX", "CELL_MATRIX_ADHESION", "CELL_SUBSTRATE_ADHESION",
        "FOCAL_ADHESION"),
    `Cell motility / cytoskeleton` =
      c("MOTILITY", "AMEBOIDAL_TYPE_CELL_MIGRATION", "TISSUE_MIGRATION",
        "CYTOSKELETON", "ACTIN", "MICROTUBULE", "TUBULIN", "SARCOMERE",
        "MYOFIBRIL", "CONTRACTILE"),
    `Proteostasis / RNA translation` =
      c("TRANSLATION", "RIBOSOM", "PROTEASOM", "UBIQUITIN",
        "PROTEIN_FOLDING", "RNA_PROCESSING", "MRNA"),
    `Mitochondrial metabolism / oxidative stress` =
      c("MITOCHONDR", "OXIDATIVE", "RESPIRATORY_CHAIN", "ATP",
        "ELECTRON_TRANSPORT", "REACTIVE_OXYGEN"),
    `Immune / inflammatory signalling` =
      c("IMMUNE", "INFLAMMATORY", "CYTOKINE", "INTERFERON", "LEUKOCYTE",
        "MYELOID", "MACROPHAGE", "MICROGLIA", "TOLL_LIKE")
  )
  pathway_names <- toupper(names(pathway_sets))
  out <- lapply(role_patterns, function(patterns) {
    hit <- vapply(pathway_names, function(nm) {
      any(vapply(patterns, grepl, logical(1), x = nm, fixed = TRUE))
    }, logical(1), USE.NAMES = FALSE)
    .fig_gene_tokens(unlist(pathway_sets[hit], use.names = FALSE))
  })
  out[vapply(out, length, integer(1)) > 0L]
}

.fig_offdiag_gene_rows <- function(panels, order) {
  .fig_bind(lapply(order, function(m) {
    stopifnot(m %in% names(panels), is.data.frame(panels[[m]]$data))
    d <- panels[[m]]$data
    .fig_require_cols(d, c("feature", "label", "gene_symbols", "x", "y", "interaction",
                           "abs_interaction"),
                      paste0("modality panel ", m))
    rows <- lapply(seq_len(nrow(d)), function(i) {
      g <- .fig_gene_tokens(d$gene_symbols[i])
      if (!length(g)) return(NULL)
      label_rank <- if ("scatter_label_rank" %in% names(d)) d$scatter_label_rank[i] else NA_integer_
      data.frame(
        modality = m,
        feature = d$feature[i],
        label = d$label[i],
        score_feature = d$feature[i],
        score_label = d$label[i],
        gene_symbol = g,
        scatter_label_rank = label_rank,
        x = d$x[i],
        y = d$y[i],
        interaction = d$interaction[i],
        abs_interaction = d$abs_interaction[i],
        offdiag_distance = if ("offdiag_distance" %in% names(d)) d$offdiag_distance[i] else d$abs_interaction[i],
        stringsAsFactors = FALSE
      )
    })
    out <- .fig_bind(rows)
    if (!nrow(out)) return(out)
    out$modality <- factor(out$modality, levels = order)
    out
  }))
}

.fig_predicted_unannotated_role <- "Predicted / unannotated loci"
.fig_olfactory_role <- "Olfactory receptor"
.fig_other_annotated_role <- "Other annotated / no role-set hit"
.fig_uncategorized_roles <- c(.fig_predicted_unannotated_role, .fig_other_annotated_role)

# Functional-category score summary for shared-cutoff off-diagonal features in
# fig-modality-amyloid-effect. Default role categories are priority-ordered GO-BP term-family
# unions: complement/MHC, phagocytosis, chemotaxis, lipid/endolysosome/synapse, then broader
# cell-cell-adhesion/extracellular-matrix/motility buckets, with broad immune/inflammatory last. Explicit classified
# fallbacks are retained, but uncategorized fallback buckets are omitted from the visible summary.
# The phosphoproteomics panel is already parent-protein-collapsed upstream, so category scores
# use the same averaged protein points displayed in the amyloid-effect scatter.
.fig_fallback_role <- function(gene_symbol, label) {
  token <- .fig_gene_tokens(c(gene_symbol, label))
  if (!length(token)) token <- as.character(label)
  if (any(grepl("^Olfr[0-9]", token, perl = TRUE))) return(.fig_olfactory_role)
  if (any(grepl("^(Gm[0-9]+|LOC[0-9]+)$", token, perl = TRUE) |
          grepl("Rik$", token, perl = TRUE))) {
    return(.fig_predicted_unannotated_role)
  }
  .fig_other_annotated_role
}

.fig_wrapped_feature_list <- function(x, per_line = 3L, max_lines = 2L) {
  stopifnot(is.numeric(per_line), length(per_line) == 1L, per_line >= 1L,
            is.numeric(max_lines), length(max_lines) == 1L, max_lines >= 1L)
  x <- unique(trimws(as.character(x)))
  x <- x[!is.na(x) & x != ""]
  if (!length(x)) return("")
  per_line <- max(as.integer(per_line), ceiling(length(x) / as.integer(max_lines)))
  line <- if (as.integer(max_lines) == 2L && length(x) > per_line) {
    c(rep(1L, floor(length(x) / 2)), rep(2L, ceiling(length(x) / 2)))
  } else {
    ceiling(seq_along(x) / as.integer(per_line))
  }
  paste(vapply(split(x, line), paste, character(1), collapse = ", "),
        collapse = "\n")
}

.fig_primary_role <- function(gene_symbol, label, group_sets, group_labels, fallback_priority) {
  g <- .fig_gene_tokens(gene_symbol)
  hit <- which(vapply(group_sets, function(set) any(g %in% set), logical(1), USE.NAMES = FALSE))
  if (length(hit)) {
    i <- hit[[1]]
    return(data.frame(group = names(group_sets)[i], group_label = group_labels[i],
                      group_priority = i, stringsAsFactors = FALSE))
  }
  fallback_order <- c(.fig_predicted_unannotated_role, .fig_olfactory_role,
                      .fig_other_annotated_role)
  lab <- .fig_fallback_role(gene_symbol, label)
  data.frame(group = lab, group_label = lab,
             group_priority = fallback_priority + match(lab, fallback_order),
             stringsAsFactors = FALSE)
}

modality_offdiag_group_score_data <- function(modality_scatter_figures,
                                              group_sets = NULL,
                                              tail_quantile = 0.99,
                                              offdiag_cutoff = NULL,
                                              cutoff_source = NULL,
                                              min_genes = 1L,
                                              max_groups = 10L,
                                              min_abs_delta = 0) {
  stopifnot(is.list(modality_scatter_figures), is.list(modality_scatter_figures$panels),
            is.character(modality_scatter_figures$order),
            length(modality_scatter_figures$order) >= 1L,
            is.numeric(tail_quantile), length(tail_quantile) == 1L,
            tail_quantile > 0, tail_quantile < 1,
            is.null(offdiag_cutoff) ||
              (is.numeric(offdiag_cutoff) && length(offdiag_cutoff) == 1L),
            is.numeric(min_genes), length(min_genes) == 1L, min_genes >= 1L,
            is.numeric(max_groups), length(max_groups) == 1L, max_groups >= 1L,
            is.numeric(min_abs_delta), length(min_abs_delta) == 1L,
            is.finite(min_abs_delta), min_abs_delta >= 0)
  order <- modality_scatter_figures$order
  group_set_source <- if (is.null(group_sets)) {
    "priority-ordered MSigDB mouse GO Biological Process keyword unions via msigdbr"
  } else {
    "custom functional groups"
  }
  group_sets <- group_sets %||% .fig_default_functional_groups()
  stopifnot(is.list(group_sets), length(group_sets) >= 1L, !is.null(names(group_sets)),
            !any(names(group_sets) == ""))
  group_sets <- lapply(group_sets, .fig_gene_tokens)
  group_sets <- group_sets[vapply(group_sets, length, integer(1)) > 0L]
  stopifnot(length(group_sets) >= 1L)
  group_labels <- .fig_functional_group_label(.fig_pathway_label(names(group_sets)))
  fallback_priority <- length(group_sets)

  selected <- .fig_bind(lapply(order, function(m) {
    d <- modality_scatter_figures$panels[[m]]$data
    d <- modality_scatter_label_rows(d, label_col = "label",
                                     tail_quantile = tail_quantile,
                                     cutoff = offdiag_cutoff,
                                     cutoff_source = cutoff_source)
    d <- .fig_offdiag_gene_rows(list(.panel = list(data = d)), ".panel")
    if (nrow(d)) d$modality <- factor(m, levels = order)
    if (!nrow(d)) return(d)
    d <- d[order(d$scatter_label_rank, d$gene_symbol, d$feature, method = "radix"), , drop = FALSE]
    if (identical(m, "Phospho")) {
      d$score_feature <- paste0("phospho_gene:", d$gene_symbol)
      d$score_label <- d$gene_symbol
    }
    role <- .fig_bind(lapply(seq_len(nrow(d)), function(i) {
      .fig_primary_role(d$gene_symbol[i], d$score_label[i], group_sets, group_labels,
                        fallback_priority)
    }))
    d <- data.frame(d, role, stringsAsFactors = FALSE)
    d <- d[order(d$group_priority, d$scatter_label_rank, d$gene_symbol, d$feature,
                 method = "radix"), , drop = FALSE]
    d <- d[!duplicated(d$score_feature), , drop = FALSE]
    d$scatter_offdiag <- TRUE
    d
  }))
  .fig_assert_nonempty(selected, "empirical off-diagonal genes/proteins")

  rows <- .fig_bind(lapply(order, function(m) {
    sel <- selected[as.character(selected$modality) == m, , drop = FALSE]
    if (!nrow(sel)) return(data.frame())
    group_order <- unique(sel$group[order(sel$group_priority, sel$group, method = "radix")])
    .fig_bind(lapply(group_order, function(grp) {
      feature_hits <- sel[sel$group == grp, , drop = FALSE]
      feature_hits <- feature_hits[order(feature_hits$scatter_label_rank,
                                         feature_hits$score_label,
                                         feature_hits$feature, method = "radix"), , drop = FALSE]
      k_feature <- length(unique(feature_hits$score_feature))
      if (k_feature < min_genes) return(data.frame())
      k_gene <- length(unique(feature_hits$gene_symbol))
      score_maptki <- mean(feature_hits$y)
      score_p301s <- mean(feature_hits$x)
      delta <- score_p301s - score_maptki
      data.frame(
        modality = m,
        group = grp,
        group_label = feature_hits$group_label[[1]],
        group_priority = feature_hits$group_priority[[1]],
        n_gene = k_gene,
        n_feature = k_feature,
        n_selected = length(unique(sel$score_feature)),
        n_labeled_feature = length(unique(sel$score_feature)),
        n_p301s_higher = sum(feature_hits$interaction > 0),
        n_maptki_higher = sum(feature_hits$interaction < 0),
        score_maptki = score_maptki,
        score_p301s = score_p301s,
        delta = delta,
        abs_delta = abs(delta),
        mean_abs_feature_delta = mean(feature_hits$abs_interaction),
        direction = if (delta >= 0) "P301S higher" else "MAPTKI higher",
        top_genes = paste(unique(feature_hits$gene_symbol), collapse = ", "),
        top_features = paste(unique(feature_hits$score_label), collapse = ", "),
        top_features_plot = .fig_wrapped_feature_list(feature_hits$score_label,
                                                      per_line = 3L),
        stringsAsFactors = FALSE
      )
    }))
  }))
  .fig_assert_nonempty(rows, "off-diagonal functional-category score summary")
  rows <- rows[is.finite(rows$score_maptki) & is.finite(rows$score_p301s) &
                 is.finite(rows$delta), , drop = FALSE]
  .fig_assert_nonempty(rows, "finite off-diagonal functional-category score summary")
  rows <- rows[!rows$group %in% .fig_uncategorized_roles, , drop = FALSE]
  .fig_assert_nonempty(rows, "categorized off-diagonal functional-category score summary")
  if (min_abs_delta > 0) {
    rows <- rows[abs(rows$delta) >= min_abs_delta, , drop = FALSE]
    .fig_assert_nonempty(rows, "delta-filtered off-diagonal functional-category score summary")
  }
  rows$rank_score <- rows$abs_delta * log1p(rows$n_feature)
  rows <- .fig_bind(lapply(order, function(m) {
    z <- rows[as.character(rows$modality) == m, , drop = FALSE]
    z <- z[order(-z$rank_score, z$group_priority, z$group_label, method = "radix"), , drop = FALSE]
    utils::head(z, as.integer(max_groups))
  }))
  rownames(rows) <- NULL

  rows$group_label_plot <- paste0(rows$group_label, "\n", rows$top_features_plot)
  rows <- rows[order(match(as.character(rows$modality), order), -rows$rank_score,
                     rows$group_priority, rows$group_label, method = "radix"), , drop = FALSE]
  rows$group_label_plot <- factor(rows$group_label_plot, levels = rev(unique(rows$group_label_plot)))
  rows$group_label <- factor(rows$group_label, levels = rev(unique(rows$group_label)))
  rows$direction <- factor(rows$direction, levels = c("MAPTKI higher", "P301S higher"))
  rows$modality <- factor(rows$modality, levels = order)
  .fig_assert_finite(rows, c("n_gene", "n_feature", "n_selected", "n_labeled_feature",
                            "n_p301s_higher", "n_maptki_higher", "score_maptki",
                            "score_p301s", "delta", "abs_delta", "mean_abs_feature_delta",
                            "rank_score"),
                     "off-diagonal functional-category score summary")
  .fig_assert_finite(selected, c("x", "y", "interaction", "abs_interaction",
                                "offdiag_distance", "scatter_label_rank"),
                     "empirical off-diagonal genes/proteins")
  panel_attr_cutoffs <- stats::setNames(vapply(order, function(m) {
    x <- attr(modality_scatter_figures$panels[[m]]$data, "offdiag_cutoff", exact = TRUE)
    if (is.null(x)) NA_real_ else as.numeric(x[[1]])
  }, numeric(1)), order)
  if (is.null(offdiag_cutoff)) {
    panel_missing <- !is.finite(panel_attr_cutoffs)
    if (any(panel_missing)) {
      panel_attr_cutoffs[panel_missing] <-
        modality_scatter_panel_cutoffs(modality_scatter_figures$panels, order,
                                       tail_quantile = tail_quantile)[panel_missing]
    }
    provenance_cutoff <- panel_attr_cutoffs
  } else {
    provenance_cutoff <- stats::setNames(rep(as.numeric(offdiag_cutoff), length(order)), order)
  }
  panel_sources <- unique(vapply(order, function(m) {
    attr(modality_scatter_figures$panels[[m]]$data, "offdiag_cutoff_source", exact = TRUE) %||%
      "panel_tail_quantile"
  }, character(1), USE.NAMES = FALSE))
  provenance_cutoff_source <- cutoff_source %||%
    if (length(panel_sources) == 1L) panel_sources else paste(panel_sources, collapse = "; ")
  list(
    summary = rows[order(as.integer(rows$modality), -rows$rank_score, rows$group_priority,
                         method = "radix"), , drop = FALSE],
    selected_genes = selected,
    selected_outliers = selected,
    provenance = list(
      group_set_source = group_set_source,
      offdiag_tail_quantile = tail_quantile,
      offdiag_cutoff = provenance_cutoff,
      offdiag_cutoff_source = provenance_cutoff_source,
      min_abs_delta = min_abs_delta,
      min_genes = as.integer(min_genes),
      max_groups = as.integer(max_groups),
      n_group_sets = length(group_sets),
      selection = "same shared off-diagonal rule as fig-modality-amyloid-effect: features pass the stored |x-y| cutoff; duplicate display labels collapsed after thresholding",
      category_assignment = "one primary role per scored item: first matching priority-ordered GO-BP term-family union, with complement/MHC, phagocytosis, and chemotaxis separated before broader cell-cell-adhesion/extracellular-matrix/motility and immune/inflammatory residual buckets; otherwise predicted/unannotated, olfactory receptor, or other annotated fallback; visible summary excludes predicted/unannotated and other annotated/no role-set buckets; each visible category label lists every retained scored feature, wrapped into at most two lines",
      phosphoproteomics_scoring = "phosphoproteomics points are parent-protein aggregates of finite phosphosite logFC pairs; category scores use those displayed protein points",
      n_labeled_features = stats::setNames(
        vapply(order, function(m) length(unique(selected$score_feature[as.character(selected$modality) == m])),
               integer(1)), order),
      n_offdiag_features = stats::setNames(
        vapply(order, function(m) length(unique(selected$score_feature[as.character(selected$modality) == m])),
               integer(1)), order)
    )
  )
}

# P6 Figure 10 -------------------------------------------------------------------------------
# Compact, deterministic leaf for the gene-resolved state plate. The accepted occupancy panel is
# retained; the other panels expose direct paired state differences and conventional factorial
# response plots. No programme score, fitted object, count matrix, or outcome-selected gene set
# reaches the report.
state_decomposition_figure_data <- function(state_response, gene_atlas, alpha = 0.05,
                                            max_target_bytes = 6 * 1024^2) {
  stopifnot(
    is.list(state_response), identical(state_response$schema, "p6_state_response_v1"),
    is.list(state_response$occupancy),
    is.list(gene_atlas), identical(gene_atlas$schema, "p6_state_gene_atlas_v2"),
    is.data.frame(gene_atlas$features), is.list(gene_atlas$effects),
    is.data.frame(gene_atlas$omnibus), is.data.frame(gene_atlas$markers),
    length(alpha) == 1L, is.finite(alpha), alpha == 0.05,
    length(max_target_bytes) == 1L, is.finite(max_target_bytes), max_target_bytes > 0
  )

  # A: raw replicate occupancy plus beta-binomial standardized genotype means.
  unit <- state_response$occupancy$unit
  means <- state_response$occupancy$probability_means
  contrasts <- state_response$occupancy$probability_contrasts
  .fig_require_cols(unit, c("genotype_batch", "genotype", "batch", "n_Homeostatic",
                            "n_DAM", "n_primary", "coverage", "DAM_fraction"),
                    "state decomposition occupancy units")
  .fig_require_cols(means, c("genotype", "estimate", "se", "ci_l", "ci_r"),
                    "state decomposition occupancy means")
  .fig_require_cols(contrasts, c("contrast", "estimate", "se", "ci_l", "ci_r", "margin",
                                 "p_zero", "fdr_zero", "p_minimum", "fdr_minimum"),
                    "state decomposition occupancy contrasts")
  unit <- unit[order(match(as.character(unit$genotype), genotype_levels),
                     as.character(unit$batch), method = "radix"), , drop = FALSE]
  means <- means[match(genotype_levels, as.character(means$genotype)), , drop = FALSE]
  occupancy_interaction <- contrasts[contrasts$contrast == "interaction", , drop = FALSE]
  stopifnot(
    nrow(unit) == 16L, !anyDuplicated(unit$genotype_batch),
    identical(as.integer(table(factor(as.character(unit$genotype), levels = genotype_levels))),
              rep.int(4L, length(genotype_levels))),
    nrow(means) == 4L, identical(as.character(means$genotype), genotype_levels),
    nrow(occupancy_interaction) == 1L, occupancy_interaction$margin == 0.10
  )
  unit$genotype <- factor(as.character(unit$genotype), levels = genotype_levels)
  means$genotype <- factor(as.character(means$genotype), levels = genotype_levels)
  .fig_assert_finite(unit, c("n_Homeostatic", "n_DAM", "n_primary", "coverage",
                             "DAM_fraction"), "state decomposition occupancy units")
  .fig_assert_finite(means, c("estimate", "se", "ci_l", "ci_r"),
                     "state decomposition occupancy means")
  stopifnot(all(unit$DAM_fraction >= 0 & unit$DAM_fraction <= 1),
            all(means$estimate >= 0 & means$estimate <= 1),
            all(means$ci_l <= means$estimate & means$estimate <= means$ci_r))

  # B/C share the paired atlas. C consumes the four state/background amyloid responses; B uses
  # their two direct DAM-minus-Homeostatic contrasts, including paired-model CI and BH inference.
  response_meta <- data.frame(
    contrast = c("homeostatic_maptki", "dam_maptki",
                 "homeostatic_p301s", "dam_p301s"),
    background = c("MAPTKI", "MAPTKI", "P301S", "P301S"),
    state = c("Homeostatic", "DAM", "Homeostatic", "DAM"),
    stringsAsFactors = FALSE
  )
  difference_meta <- data.frame(
    background = c("MAPTKI", "P301S"),
    homeostatic = c("homeostatic_maptki", "homeostatic_p301s"),
    dam = c("dam_maptki", "dam_p301s"),
    difference = c("dam_minus_homeostatic_maptki",
                   "dam_minus_homeostatic_p301s"),
    stringsAsFactors = FALSE
  )
  contrast_order <- c(
    "homeostatic_maptki", "homeostatic_p301s", "dam_maptki", "dam_p301s",
    "dam_minus_homeostatic_maptki", "dam_minus_homeostatic_p301s",
    "homeostatic_interaction", "dam_interaction", "dam_minus_homeostatic_interaction"
  )
  matrix_names <- c("estimate", "se", "ci95_l", "ci95_r",
                    "p_value", "fdr", "treat_p", "treat_fdr")
  stopifnot(
    all(matrix_names %in% names(gene_atlas$effects)),
    all(vapply(gene_atlas$effects[matrix_names], is.matrix, logical(1))),
    identical(colnames(gene_atlas$effects$estimate), contrast_order),
    all(vapply(gene_atlas$effects[matrix_names], function(x)
      identical(dimnames(x), dimnames(gene_atlas$effects$estimate)), logical(1)))
  )
  features <- gene_atlas$features
  omnibus <- gene_atlas$omnibus
  .fig_require_cols(features, c("gene", "symbol", "ave_expr"), "state gene atlas features")
  .fig_require_cols(
    omnibus,
    c("gene", "response_F", "response_p", "response_fdr",
      "interaction_F", "interaction_p", "interaction_fdr"),
    "state gene atlas omnibus inference"
  )
  genes <- rownames(gene_atlas$effects$estimate)
  stopifnot(identical(features$gene, genes), identical(omnibus$gene, genes),
            !anyDuplicated(genes), all(nzchar(genes)))

  programs <- names(canonical_microglia_markers)
  markers <- gene_atlas$markers
  .fig_require_cols(markers, c("program", "symbol", "gene", "present", "detected"),
                    "state gene atlas markers")
  declared <- do.call(rbind, lapply(programs, function(program)
    data.frame(program = program, symbol = canonical_microglia_markers[[program]],
               stringsAsFactors = FALSE)))
  rownames(declared) <- NULL
  stopifnot(
    nrow(markers) == nrow(declared), nrow(markers) == sum(lengths(canonical_microglia_markers)),
    identical(markers[c("program", "symbol")], declared),
    !anyDuplicated(markers[c("program", "symbol")]),
    all(markers$detected <= markers$present)
  )
  marker_genes <- unique(markers[c("symbol", "gene", "present", "detected")])
  stopifnot(nrow(marker_genes) == length(unique(markers$gene)),
            all(marker_genes$present), !anyNA(marker_genes$gene),
            !anyDuplicated(marker_genes$gene), !anyDuplicated(marker_genes$symbol))
  marker_i <- match(marker_genes$gene, genes)
  detected <- marker_genes$detected & !is.na(marker_i)
  stopifnot(identical(detected, marker_genes$detected))
  marker_rows <- lapply(seq_len(nrow(response_meta)), function(j) {
    contrast <- response_meta$contrast[[j]]
    estimate <- ci_l <- ci_r <- fdr <- treat_fdr <- rep(NA_real_, nrow(marker_genes))
    estimate[detected] <- gene_atlas$effects$estimate[marker_i[detected], contrast]
    ci_l[detected] <- gene_atlas$effects$ci95_l[marker_i[detected], contrast]
    ci_r[detected] <- gene_atlas$effects$ci95_r[marker_i[detected], contrast]
    fdr[detected] <- gene_atlas$effects$fdr[marker_i[detected], contrast]
    treat_fdr[detected] <- gene_atlas$effects$treat_fdr[marker_i[detected], contrast]
    data.frame(
      marker_genes[c("symbol", "gene")],
      detected = detected,
      contrast = contrast,
      background = response_meta$background[[j]],
      state = response_meta$state[[j]],
      estimate = estimate, ci_l = ci_l, ci_r = ci_r,
      fdr = fdr, treat_fdr = treat_fdr,
      row.names = NULL, stringsAsFactors = FALSE
    )
  })
  marker_effects <- do.call(rbind, marker_rows)
  rownames(marker_effects) <- NULL
  marker_effects$background <- factor(
    marker_effects$background, levels = difference_meta$background
  )
  marker_effects$state <- factor(
    marker_effects$state, levels = c("Homeostatic", "DAM")
  )
  marker_effects$support <- factor(
    ifelse(!marker_effects$detected, "below count filter",
           ifelse(marker_effects$treat_fdr <= alpha, "minimum-effect FDR <= 0.05",
                  ifelse(marker_effects$fdr <= alpha, "nonzero FDR <= 0.05",
                         "not supported"))),
    levels = c("minimum-effect FDR <= 0.05", "nonzero FDR <= 0.05",
               "not supported", "below count filter")
  )
  detected_effects <- marker_effects[marker_effects$detected, , drop = FALSE]
  .fig_assert_finite(detected_effects, c("estimate", "ci_l", "ci_r", "fdr", "treat_fdr"),
                     "state gene atlas detected marker effects")
  stopifnot(
    nrow(marker_effects) == nrow(marker_genes) * nrow(response_meta),
    all(detected_effects$ci_l <= detected_effects$estimate &
          detected_effects$estimate <= detected_effects$ci_r),
    all(is.na(marker_effects$estimate[!marker_effects$detected])),
    all(is.na(marker_effects$ci_l[!marker_effects$detected])),
    all(is.na(marker_effects$ci_r[!marker_effects$detected])),
    all(is.na(marker_effects$fdr[!marker_effects$detected])),
    all(is.na(marker_effects$treat_fdr[!marker_effects$detected])),
    !anyNA(marker_effects$support)
  )

  # Evidence-first ordering changes facet position only; all 52 declared unique genes remain fixed.
  gene_summary <- do.call(rbind, lapply(
    split(detected_effects, detected_effects$gene),
    function(x) {
      state_difference <- vapply(difference_meta$background, function(background) {
        z <- x[as.character(x$background) == background, , drop = FALSE]
        stopifnot(nrow(z) == 2L, identical(as.character(z$state), c("Homeostatic", "DAM")))
        z$estimate[[2L]] - z$estimate[[1L]]
      }, numeric(1))
      data.frame(
        gene = x$gene[[1L]], symbol = x$symbol[[1L]],
        any_minimum = any(x$treat_fdr <= alpha), any_nonzero = any(x$fdr <= alpha),
        max_state_difference = max(abs(state_difference)),
        max_abs_effect = max(abs(x$estimate)),
        row.names = NULL, stringsAsFactors = FALSE
      )
    }
  ))
  gene_summary <- gene_summary[order(
    !gene_summary$any_minimum, !gene_summary$any_nonzero,
    -gene_summary$max_state_difference, -gene_summary$max_abs_effect,
    gene_summary$symbol, gene_summary$gene, method = "radix"
  ), , drop = FALSE]
  marker_effects$gene_label <- factor(marker_effects$symbol, levels = gene_summary$symbol)
  stopifnot(nrow(gene_summary) == sum(marker_genes$detected),
            !anyDuplicated(gene_summary$gene), !anyDuplicated(gene_summary$symbol),
            !anyNA(marker_effects$gene_label[marker_effects$detected]))

  # Every filter-passing gene appears once per tau background. The ordinate and CI are direct
  # paired contrasts; point labels are restricted to predeclared markers clearing treat BH.
  marker_gene_ids <- marker_genes$gene[marker_genes$detected]
  state_difference <- do.call(rbind, lapply(seq_len(nrow(difference_meta)), function(j) {
    meta <- difference_meta[j, , drop = FALSE]
    homeostatic_effect <- gene_atlas$effects$estimate[, meta$homeostatic]
    dam_effect <- gene_atlas$effects$estimate[, meta$dam]
    difference_effect <- gene_atlas$effects$estimate[, meta$difference]
    stopifnot(max(abs(difference_effect - (dam_effect - homeostatic_effect))) <= 1e-10)
    data.frame(
      gene = genes, symbol = features$symbol,
      background = meta$background,
      mean_effect = (homeostatic_effect + dam_effect) / 2,
      state_difference = difference_effect,
      ci_l = gene_atlas$effects$ci95_l[, meta$difference],
      ci_r = gene_atlas$effects$ci95_r[, meta$difference],
      fdr = gene_atlas$effects$fdr[, meta$difference],
      treat_fdr = gene_atlas$effects$treat_fdr[, meta$difference],
      is_declared_marker = genes %in% marker_gene_ids,
      row.names = NULL, stringsAsFactors = FALSE
    )
  }))
  state_difference$background <- factor(
    state_difference$background, levels = difference_meta$background
  )
  state_difference$support <- factor(
    ifelse(state_difference$treat_fdr <= alpha, "minimum-effect FDR <= 0.05",
           ifelse(state_difference$fdr <= alpha, "nonzero FDR <= 0.05",
                  "not supported")),
    levels = c("minimum-effect FDR <= 0.05", "nonzero FDR <= 0.05", "not supported")
  )
  state_difference$label <- ifelse(
    state_difference$is_declared_marker & state_difference$treat_fdr <= alpha,
    state_difference$symbol, ""
  )
  .fig_assert_finite(
    state_difference,
    c("mean_effect", "state_difference", "ci_l", "ci_r", "fdr", "treat_fdr"),
    "state gene direct-difference map"
  )
  stopifnot(
    nrow(state_difference) == length(genes) * nrow(difference_meta),
    all(state_difference$ci_l <= state_difference$state_difference &
          state_difference$state_difference <= state_difference$ci_r),
    !anyNA(state_difference$symbol), all(nzchar(state_difference$symbol)),
    !anyNA(state_difference$support)
  )

  difference_fdr_counts <- tapply(
    state_difference$fdr <= alpha, state_difference$background, sum
  )
  difference_minimum_counts <- tapply(
    state_difference$treat_fdr <= alpha, state_difference$background, sum
  )
  stopifnot(
    identical(as.integer(difference_fdr_counts),
              as.integer(gene_atlas$audit$state_difference_fdr_supported)),
    identical(as.integer(difference_minimum_counts),
              as.integer(gene_atlas$audit$state_difference_minimum_supported))
  )

  out <- list(
    schema = "p6_state_decomposition_figures_v4",
    occupancy = list(unit = unit, means = means),
    gene_atlas = list(marker_effects = marker_effects, state_difference = state_difference),
    provenance = list(
      source_targets = c("microglia_state_response", "microglia_state_gene_atlas"),
      defining_contrasts = difference_meta$difference,
      method = gene_atlas$audit$method,
      marker_program_order = programs,
      response_contrast_order = response_meta$contrast,
      selection = paste(
        "all 52 declared unique marker genes are fixed before outcomes;",
        "B contains every filter-passing gene in both backgrounds;",
        "labels require declared-marker status plus direct minimum-effect BH support"
      ),
      marker_ordering = "minimum-effect, then nonzero evidence, state difference, magnitude, symbol",
      alpha = alpha,
      gene_lfc = gene_atlas$audit$thresholds$gene_lfc,
      occupancy_margin = occupancy_interaction$margin
    ),
    audit = list(
      row_counts = c(
        occupancy_unit = nrow(unit), occupancy_means = nrow(means),
        marker_effects = nrow(marker_effects), state_difference = nrow(state_difference)
      ),
      n_marker_memberships = nrow(markers),
      n_marker_genes = nrow(marker_genes),
      n_marker_genes_detected = sum(marker_genes$detected),
      n_marker_response_minimum_genes = length(unique(
        detected_effects$gene[detected_effects$treat_fdr <= alpha]
      )),
      n_state_difference_fdr_supported = difference_fdr_counts,
      n_state_difference_minimum_supported = difference_minimum_counts,
      n_marker_state_difference_minimum_rows = sum(nzchar(state_difference$label)),
      n_interaction_fdr_supported = sum(omnibus$interaction_fdr <= alpha),
      below_filter_symbols = marker_genes$symbol[!marker_genes$detected],
      all_declared_markers_represented = TRUE,
      fixed_marker_rows = TRUE,
      deterministic_order = TRUE,
      parent_isolated = NA,
      in_memory_bytes = NA_real_, serialized_bytes = NA_real_,
      max_target_bytes = max_target_bytes
    )
  )
  stopifnot(!state_substrate_contains_parent(out))
  out$audit$parent_isolated <- TRUE
  out$audit$in_memory_bytes <- as.numeric(object.size(out))
  out$audit$serialized_bytes <- as.numeric(length(qs2::qs_serialize(out)))
  stopifnot(out$audit$in_memory_bytes <= max_target_bytes,
            out$audit$serialized_bytes <= max_target_bytes)
  out
}


# P8 Figures 11-13: compact report-only summaries of the three integration leaves. The
# pathway payload is deliberately reduced to counts and a fixed top-set display; no full
# score or consensus table crosses the report boundary.
integration_figure_data <- function(integration_decomposition,
                                    integration_concordance,
                                    integration_pathway,
                                    top_pathways_per_contrast = 5L,
                                    max_target_bytes = 1024^2) {
  stopifnot(
    is.list(integration_decomposition), is.list(integration_concordance),
    is.list(integration_pathway),
    length(top_pathways_per_contrast) == 1L,
    is.finite(top_pathways_per_contrast), top_pathways_per_contrast >= 1,
    top_pathways_per_contrast == as.integer(top_pathways_per_contrast),
    length(max_target_bytes) == 1L, is.finite(max_target_bytes), max_target_bytes > 0
  )
  top_pathways_per_contrast <- as.integer(top_pathways_per_contrast)
  modalities <- c("snRNAseq", "GeoMx", "bulk")
  modality_labels <- c(snRNAseq = "snRNA-seq", GeoMx = "GeoMx", bulk = "bulk TiO2")
  contrasts <- c("tau_alone", "nlgf_in_maptki", "nlgf_in_p301s",
                 "tau_in_nlgf", "interaction")
  contrast_labels <- c(
    tau_alone = "tau | no amyloid",
    nlgf_in_maptki = "amyloid | MAPTKI",
    nlgf_in_p301s = "amyloid | P301S",
    tau_in_nlgf = "tau | amyloid",
    interaction = "interaction"
  )
  amyloid_contrasts <- c("nlgf_in_maptki", "nlgf_in_p301s")
  shared_universe <- integration_decomposition$provenance$universe
  shared_n_genes <- integration_decomposition$provenance$n_genes
  stopifnot(
    identical(shared_universe, "complete_case"),
    identical(shared_n_genes, 3109L),
    identical(integration_concordance$provenance$universe, shared_universe),
    identical(integration_concordance$provenance$n_genes, shared_n_genes)
  )

  # Figure 11: the selected joint rank is zero, so retain the variance partition and
  # the top unselected shared-candidate diagnostic instead of empty joint loadings.
  decomposition <- integration_decomposition$primary
  stopifnot(
    is.list(decomposition), identical(decomposition$statistic, "standardized_logFC"),
    identical(integration_decomposition$provenance$modalities, modalities),
    identical(integration_decomposition$provenance$contrasts, contrasts),
    is.matrix(decomposition$joint_gene_scores), ncol(decomposition$joint_gene_scores) == 0L,
    is.list(decomposition$joint_contrast_loadings),
    all(vapply(decomposition$joint_contrast_loadings, ncol, integer(1)) == 0L)
  )
  variance <- decomposition$variance
  .fig_require_cols(variance, c("modality", "component", "sum_squares", "fraction"),
                    "integration decomposition variance")
  variance <- variance[
    variance$modality %in% modalities &
      variance$component %in% c("joint", "individual", "residual") &
      is.finite(variance$sum_squares) & is.finite(variance$fraction),
    c("modality", "component", "sum_squares", "fraction"), drop = FALSE
  ]
  variance <- variance[order(match(variance$modality, modalities),
                             match(variance$component,
                                   c("joint", "individual", "residual"))), , drop = FALSE]
  rownames(variance) <- NULL
  stopifnot(nrow(variance) == length(modalities) * 3L,
            !anyDuplicated(variance[c("modality", "component")]),
            all(variance$fraction >= 0 & variance$fraction <= 1))
  fraction_sums <- vapply(modalities, function(modality) {
    sum(variance$fraction[variance$modality == modality])
  }, numeric(1))
  stopifnot(max(abs(fraction_sums - 1)) <= 1e-8)
  variance$modality_label <- unname(modality_labels[variance$modality])
  variance$component_label <- unname(c(
    joint = "joint", individual = "block-individual", residual = "residual"
  )[variance$component])

  ranks <- decomposition$ranks
  .fig_require_cols(ranks, c("modality", "signal_rank", "joint_rank",
                             "individual_rank", "residual_rank_budget"),
                    "integration decomposition ranks")
  joint_rank <- unique(as.integer(ranks$joint_rank))
  stopifnot(length(joint_rank) == 1L, !is.na(joint_rank), joint_rank == 0L)
  ranks <- ranks[match(modalities, as.character(ranks$modality)),
                 c("modality", "signal_rank", "joint_rank", "individual_rank",
                   "residual_rank_budget"), drop = FALSE]
  stopifnot(!anyNA(ranks$modality))
  rownames(ranks) <- NULL

  joint_diagnostic <- decomposition$diagnostics$joint
  stopifnot(is.list(joint_diagnostic), is.matrix(joint_diagnostic$block_alignment),
            is.data.frame(joint_diagnostic$candidates),
            is.list(joint_diagnostic$thresholds))
  candidates <- joint_diagnostic$candidates
  .fig_require_cols(candidates, c("axis", "stacked_singular_value", "selected"),
                    "integration joint candidates")
  candidate_axis <- as.character(candidates$axis[[which.max(candidates$stacked_singular_value)]])
  candidate_row <- match(candidate_axis, rownames(joint_diagnostic$block_alignment))
  candidate_meta <- candidates[match(candidate_axis, candidates$axis), , drop = FALSE]
  stopifnot(!is.na(candidate_row), nrow(candidate_meta) == 1L,
            identical(candidate_meta$selected, FALSE),
            all(modalities %in% colnames(joint_diagnostic$block_alignment)))
  squared_cosine <- as.numeric(
    joint_diagnostic$block_alignment[candidate_row, modalities, drop = TRUE]
  )
  squared_cosine <- pmin(1, pmax(0, squared_cosine))
  alignment_threshold <- joint_diagnostic$thresholds$block_alignment
  stopifnot(is.numeric(alignment_threshold), all(modalities %in% names(alignment_threshold)))
  alignment <- data.frame(
    modality = modalities,
    modality_label = unname(modality_labels[modalities]),
    candidate_axis = candidate_axis,
    squared_cosine = squared_cosine,
    principal_angle_degrees = acos(sqrt(squared_cosine)) * 180 / pi,
    alignment_threshold = as.numeric(alignment_threshold[modalities]),
    selected = rep(FALSE, length(modalities)),
    stringsAsFactors = FALSE
  )
  .fig_assert_finite(alignment,
                     c("squared_cosine", "principal_angle_degrees", "alignment_threshold"),
                     "integration shared-candidate alignment")

  # Figure 12: raw-logFC Spearman rho plus exact same-sign fractions on the common universe.
  rho_matrix <- integration_concordance$primary$spearman_logFC
  stopifnot(is.matrix(rho_matrix), identical(colnames(rho_matrix), contrasts),
            identical(rownames(rho_matrix),
                      c("snRNAseq_GeoMx", "snRNAseq_bulk", "GeoMx_bulk")),
            identical(integration_concordance$provenance$n_genes, 3109L))
  pair_labels <- c(
    snRNAseq_GeoMx = "snRNA-seq - GeoMx",
    snRNAseq_bulk = "snRNA-seq - bulk",
    GeoMx_bulk = "GeoMx - bulk"
  )
  rho <- expand.grid(
    pair = rownames(rho_matrix), contrast = colnames(rho_matrix),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  rho$rho <- as.numeric(rho_matrix)
  rho$pair_label <- unname(pair_labels[rho$pair])
  rho$contrast_label <- unname(contrast_labels[rho$contrast])
  rho$n <- as.integer(integration_concordance$provenance$n_genes)
  .fig_assert_finite(rho, c("rho", "n"), "integration concordance heatmap")

  directional <- integration_concordance$directional_overlap$sign_concordance
  .fig_require_cols(
    directional,
    c("pair", "contrast", "concordant", "discordant", "concordant_fraction", "n"),
    "integration directional sign concordance"
  )
  directional <- directional[
    directional$pair %in% names(pair_labels) & directional$contrast %in% contrasts &
      is.finite(directional$concordant) & is.finite(directional$discordant) &
      is.finite(directional$concordant_fraction) & is.finite(directional$n),
    c("pair", "contrast", "concordant", "discordant", "concordant_fraction", "n"),
    drop = FALSE
  ]
  directional_key <- paste(directional$pair, directional$contrast, sep = "\r")
  rho_key <- paste(rho$pair, rho$contrast, sep = "\r")
  directional <- directional[match(rho_key, directional_key), , drop = FALSE]
  stopifnot(nrow(directional) == nrow(rho), !anyNA(directional$pair),
            all(directional$concordant + directional$discordant == directional$n),
            all(directional$n == rho$n))
  directional$pair_label <- unname(pair_labels[directional$pair])
  directional$contrast_label <- unname(contrast_labels[directional$contrast])
  rownames(directional) <- NULL

  # Figure 13: count every scored set for the directional-consensus bars, then retain only
  # all-three-covered top GO-BP rows for the small amyloid-contrast score heatmap.
  consensus <- integration_pathway$consensus
  scores <- integration_pathway$scores
  .fig_require_cols(
    consensus,
    c("collection", "set", "contrast", "set_size", "n_modalities_covered",
      "n_up", "n_down", "consensus_direction"),
    "integration pathway consensus"
  )
  .fig_require_cols(
    scores,
    c("collection", "set", "contrast", "modality", "set_size", "coverage",
      "score_logFC"),
    "integration pathway scores"
  )
  stopifnot(
    identical(integration_pathway$provenance$gene_sets$total_sets, 7540L),
    identical(integration_pathway$provenance$scoring$min_coverage, 5L),
    identical(integration_pathway$provenance$scoring$score_threshold, 0.5)
  )
  min_consensus_modalities <- 2L
  directions <- c("down", "none", "up")
  consensus <- consensus[
    consensus$contrast %in% contrasts & consensus$consensus_direction %in% directions,
    , drop = FALSE
  ]
  consensus_table <- table(
    factor(consensus$contrast, levels = contrasts),
    factor(consensus$consensus_direction, levels = directions)
  )
  consensus_counts <- as.data.frame(consensus_table, stringsAsFactors = FALSE)
  names(consensus_counts) <- c("contrast", "direction", "n")
  consensus_counts$contrast <- as.character(consensus_counts$contrast)
  consensus_counts$direction <- as.character(consensus_counts$direction)
  consensus_counts$n <- as.integer(consensus_counts$n)
  consensus_counts$contrast_label <- unname(contrast_labels[consensus_counts$contrast])
  consensus_counts$signed_count <- ifelse(
    consensus_counts$direction == "up", consensus_counts$n,
    ifelse(consensus_counts$direction == "down", -consensus_counts$n, 0L)
  )
  stopifnot(sum(consensus_counts$n) == nrow(consensus),
            all(consensus$n_modalities_covered %in% 0:3))
  coverage_summary <- data.frame(
    contrast = contrasts,
    contrast_label = unname(contrast_labels[contrasts]),
    undercovered = vapply(contrasts, function(contrast) {
      as.integer(sum(consensus$contrast == contrast & consensus$n_modalities_covered < min_consensus_modalities))
    }, integer(1)),
    eligible = vapply(contrasts, function(contrast) {
      as.integer(sum(consensus$contrast == contrast & consensus$n_modalities_covered >= min_consensus_modalities))
    }, integer(1)),
    stringsAsFactors = FALSE
  )
  stopifnot(
    all(coverage_summary$undercovered + coverage_summary$eligible ==
          integration_pathway$provenance$gene_sets$total_sets),
    all(coverage_summary$undercovered > 0L), all(coverage_summary$eligible > 0L)
  )

  eligible <- consensus[
    consensus$collection == "GO:BP" & consensus$contrast %in% amyloid_contrasts &
      consensus$consensus_direction %in% c("up", "down"),
    c("collection", "set", "contrast", "set_size", "n_modalities_covered",
      "consensus_direction"), drop = FALSE
  ]
  stopifnot(nrow(eligible) > 0L,
            !anyDuplicated(eligible[c("collection", "set", "contrast")]))
  eligible_key <- paste(eligible$collection, eligible$set, eligible$contrast, sep = "\r")
  finite_scores <- scores[
    scores$collection == "GO:BP" & scores$contrast %in% amyloid_contrasts &
      scores$modality %in% modalities & is.finite(scores$set_size) &
      is.finite(scores$coverage) & is.finite(scores$score_logFC),
    c("collection", "set", "contrast", "modality", "set_size", "coverage",
      "score_logFC"), drop = FALSE
  ]
  finite_key <- paste(finite_scores$collection, finite_scores$set,
                      finite_scores$contrast, sep = "\r")
  finite_scores <- finite_scores[finite_key %in% eligible_key, , drop = FALSE]
  stopifnot(nrow(finite_scores) > 0L,
            !anyDuplicated(finite_scores[c("collection", "set", "contrast", "modality")]))

  score_summary <- stats::aggregate(
    score_logFC ~ collection + set + contrast,
    data = finite_scores, FUN = mean
  )
  names(score_summary)[names(score_summary) == "score_logFC"] <- "mean_score_logFC"
  modality_summary <- stats::aggregate(
    modality ~ collection + set + contrast,
    data = finite_scores, FUN = function(x) length(unique(x))
  )
  names(modality_summary)[names(modality_summary) == "modality"] <- "n_modalities"
  score_summary <- merge(
    score_summary, modality_summary,
    by = c("collection", "set", "contrast"), sort = FALSE
  )
  score_key <- paste(score_summary$collection, score_summary$set,
                     score_summary$contrast, sep = "\r")
  eligible_index <- match(score_key, eligible_key)
  keep_summary <- score_summary$n_modalities == length(modalities) & !is.na(eligible_index)
  score_summary <- score_summary[keep_summary, , drop = FALSE]
  eligible_index <- eligible_index[keep_summary]
  score_summary$consensus_direction <- eligible$consensus_direction[eligible_index]
  score_summary$n_modalities_covered <- eligible$n_modalities_covered[eligible_index]
  .fig_assert_finite(score_summary, c("mean_score_logFC", "n_modalities",
                                      "n_modalities_covered"),
                     "integration top pathway summary")

  selected_summary <- .fig_bind(lapply(amyloid_contrasts, function(contrast) {
    z <- score_summary[score_summary$contrast == contrast, , drop = FALSE]
    z <- z[order(-abs(z$mean_score_logFC), z$set, method = "radix"), , drop = FALSE]
    z <- utils::head(z, top_pathways_per_contrast)
    z$rank <- seq_len(nrow(z))
    z
  }))
  stopifnot(nrow(selected_summary) ==
              length(amyloid_contrasts) * top_pathways_per_contrast)
  selected_key <- paste(selected_summary$collection, selected_summary$set,
                        selected_summary$contrast, sep = "\r")
  top_consensus <- finite_scores[
    paste(finite_scores$collection, finite_scores$set,
          finite_scores$contrast, sep = "\r") %in% selected_key,
    , drop = FALSE
  ]
  top_index <- match(
    paste(top_consensus$collection, top_consensus$set,
          top_consensus$contrast, sep = "\r"),
    selected_key
  )
  stopifnot(!anyNA(top_index),
            nrow(top_consensus) == nrow(selected_summary) * length(modalities))
  top_consensus$rank <- selected_summary$rank[top_index]
  top_consensus$mean_score_logFC <- selected_summary$mean_score_logFC[top_index]
  top_consensus$consensus_direction <- selected_summary$consensus_direction[top_index]
  top_consensus$n_modalities_covered <-
    selected_summary$n_modalities_covered[top_index]
  top_consensus$modality_label <- unname(modality_labels[top_consensus$modality])
  top_consensus$contrast_label <- unname(contrast_labels[top_consensus$contrast])
  top_consensus$set_label <- .fig_pathway_label(top_consensus$set)
  top_consensus$set_label <- vapply(top_consensus$set_label, function(label) {
    paste(strwrap(label, width = 46L), collapse = "\n")
  }, character(1), USE.NAMES = FALSE)
  top_consensus$pathway_key <- paste(top_consensus$contrast, top_consensus$set, sep = "::")
  top_consensus <- top_consensus[order(
    match(top_consensus$contrast, amyloid_contrasts), top_consensus$rank,
    match(top_consensus$modality, modalities), method = "radix"
  ), , drop = FALSE]
  rownames(top_consensus) <- NULL
  .fig_assert_finite(
    top_consensus,
    c("set_size", "coverage", "score_logFC", "rank", "mean_score_logFC",
      "n_modalities_covered"),
    "integration top pathway scores"
  )

  format_integer <- function(x) {
    format(as.integer(x), big.mark = ",", scientific = FALSE, trim = TRUE)
  }
  format_signed_word <- function(x) {
    paste0(ifelse(x < 0, "minus ", "plus "), sprintf("%.3f", abs(x)))
  }
  joint_fraction <- variance$fraction[variance$component == "joint"]
  individual_fraction <- variance$fraction[variance$component == "individual"]
  residual_fraction <- variance$fraction[variance$component == "residual"]
  joint_rank_word <- if (joint_rank == 0L) "zero" else format_integer(joint_rank)
  alignment_index <- match(modalities, alignment$modality)
  pair_order <- names(pair_labels)
  p301s_index <- match(
    paste(pair_order, "nlgf_in_p301s", sep = "\r"), rho_key
  )
  interaction_index <- match(
    paste(pair_order, "interaction", sep = "\r"), rho_key
  )
  p301s_sign_match_percent <- round(
    100 * directional$concordant_fraction[p301s_index]
  )
  undercovered_per_contrast <- unique(coverage_summary$undercovered)
  eligible_per_contrast <- unique(coverage_summary$eligible)
  stopifnot(
    length(joint_fraction) == length(modalities), all(joint_fraction == 0),
    length(individual_fraction) == length(modalities),
    length(residual_fraction) == length(modalities),
    !anyNA(alignment_index), !anyNA(p301s_index), !anyNA(interaction_index),
    length(p301s_sign_match_percent) == length(pair_order),
    all(is.finite(p301s_sign_match_percent)),
    length(undercovered_per_contrast) == 1L,
    length(eligible_per_contrast) == 1L
  )
  consensus_value <- function(contrast, direction) {
    value <- consensus_counts$n[
      consensus_counts$contrast == contrast & consensus_counts$direction == direction
    ]
    stopifnot(length(value) == 1L)
    as.integer(value)
  }
  alt_scalars <- list(
    figure_11 = c(
      shared_universe_genes = paste0(format_integer(shared_n_genes), " genes"),
      zero_joint_variance = paste0(joint_rank_word, " joint variance"),
      selected_joint_rank = paste0("r_J equals ", joint_rank_word),
      block_individual_range = paste0(
        sprintf("%.1f", 100 * min(individual_fraction)), " to ",
        sprintf("%.1f", 100 * max(individual_fraction)), " percent"
      ),
      residual_range = paste0(
        sprintf("%.1f", 100 * min(residual_fraction)), " to ",
        sprintf("%.1f", 100 * max(residual_fraction)), " percent"
      ),
      snRNAseq_squared_cosine = paste0(
        sprintf("%.3f", alignment$squared_cosine[alignment_index[[1L]]]),
        " for snRNA-seq"
      ),
      GeoMx_squared_cosine = paste0(
        sprintf("%.3f", alignment$squared_cosine[alignment_index[[2L]]]),
        " for GeoMx"
      ),
      bulk_squared_cosine = paste0(
        sprintf("%.3f", alignment$squared_cosine[alignment_index[[3L]]]),
        " for bulk"
      ),
      principal_angles = paste0(
        sprintf("%.1f", alignment$principal_angle_degrees[alignment_index[[1L]]]), ", ",
        sprintf("%.1f", alignment$principal_angle_degrees[alignment_index[[2L]]]),
        ", and ",
        sprintf("%.1f", alignment$principal_angle_degrees[alignment_index[[3L]]]),
        " degrees"
      )
    ),
    figure_12 = c(
      shared_universe_genes = paste0(format_integer(shared_n_genes), " common genes"),
      p301s_snRNAseq_GeoMx_rho = paste0(
        "rho ", sprintf("%.3f", rho$rho[p301s_index[[1L]]]),
        " for snRNA-seq versus GeoMx"
      ),
      p301s_snRNAseq_bulk_rho = paste0(
        sprintf("%.3f", rho$rho[p301s_index[[2L]]]),
        " for snRNA-seq versus bulk"
      ),
      p301s_GeoMx_bulk_rho = paste0(
        sprintf("%.3f", rho$rho[p301s_index[[3L]]]),
        " for GeoMx versus bulk"
      ),
      interaction_rho = paste0(
        format_signed_word(rho$rho[interaction_index[[1L]]]), ", ",
        format_signed_word(rho$rho[interaction_index[[2L]]]), ", and ",
        format_signed_word(rho$rho[interaction_index[[3L]]])
      ),
      p301s_sign_match_peak = paste0(
        min(p301s_sign_match_percent), " to ",
        max(p301s_sign_match_percent), " percent"
      )
    ),
    figure_13 = c(
      total_pathway_sets = paste0(
        format_integer(integration_pathway$provenance$gene_sets$total_sets),
        " GO biological-process"
      ),
      eligible_sets_per_contrast = paste0(
        format_integer(eligible_per_contrast), " sets per contrast"
      ),
      undercovered_sets_per_contrast = paste0(
        format_integer(undercovered_per_contrast), " sets per contrast"
      ),
      min_modalities = paste0(">=", min_consensus_modalities, " modalities"),
      min_covered_symbols = paste0(
        ">=", integration_pathway$provenance$scoring$min_coverage,
        " covered symbols"
      ),
      score_threshold = paste0(
        "|score| >=",
        sprintf("%.1f", integration_pathway$provenance$scoring$score_threshold)
      ),
      p301s_up = paste0(format_integer(consensus_value("nlgf_in_p301s", "up")), " up"),
      p301s_down = paste0(
        format_integer(consensus_value("nlgf_in_p301s", "down")), " down"
      ),
      maptki_up = paste0(format_integer(consensus_value("nlgf_in_maptki", "up")), " up"),
      maptki_down = paste0(
        format_integer(consensus_value("nlgf_in_maptki", "down")), " down"
      ),
      interaction_up = paste0(
        "only ", format_integer(consensus_value("interaction", "up")), " up"
      ),
      interaction_down = paste0(
        format_integer(consensus_value("interaction", "down")), " down"
      )
    )
  )
  stopifnot(
    identical(names(alt_scalars), c("figure_11", "figure_12", "figure_13")),
    all(vapply(alt_scalars, function(values) {
      is.character(values) && length(values) > 0L && !is.null(names(values)) &&
        all(nzchar(names(values))) && all(nzchar(values)) && !anyDuplicated(names(values))
    }, logical(1)))
  )

  out <- list(
    schema = "p8_integration_figures_v1",
    alt_scalars = alt_scalars,
    decomposition = list(
      variance = variance,
      alignment = alignment,
      ranks = ranks,
      joint_rank = joint_rank,
      candidate_axis = candidate_axis
    ),
    concordance = list(
      rho = rho,
      directional = directional
    ),
    pathway = list(
      consensus_counts = consensus_counts,
      coverage_summary = coverage_summary,
      top_consensus = top_consensus
    ),
    provenance = list(
      source_targets = c("integration_decomposition", "integration_concordance",
                         "integration_pathway"),
      modalities = modalities,
      modality_labels = modality_labels,
      contrasts = contrasts,
      contrast_labels = contrast_labels,
      amyloid_contrasts = amyloid_contrasts,
      universe = shared_universe,
      n_genes = shared_n_genes,
      total_pathway_sets = integration_pathway$provenance$gene_sets$total_sets,
      pathway_min_coverage = integration_pathway$provenance$scoring$min_coverage,
      pathway_min_modalities = min_consensus_modalities,
      pathway_score_threshold = integration_pathway$provenance$scoring$score_threshold,
      top_pathway_rule = paste0(
        "top ", top_pathways_per_contrast,
        " GO:BP directional-consensus sets per amyloid contrast by absolute mean finite ",
        "standardized-logFC score; all three modalities coverage-gated"
      ),
      alignment_field =
        "primary$diagnostics$joint$block_alignment['joint_candidate_1', ]",
      interpretation = paste(
        "DESCRIPTIVE ONLY: point-estimate variance, raw-logFC rank concordance,",
        "directional overlap, and pathway scores; no calibrated cross-modality p-value"
      ),
      contract = paste(
        "compact report-only variance/alignment, 3 x 5 concordance summaries,",
        "consensus counts, and fixed top-pathway rows; no full pathway score/consensus table"
      )
    ),
    audit = list(
      row_counts = c(
        variance = nrow(variance), alignment = nrow(alignment),
        rho = nrow(rho), directional = nrow(directional),
        consensus_counts = nrow(consensus_counts),
        coverage_summary = nrow(coverage_summary),
        top_consensus = nrow(top_consensus)
      ),
      consensus_direction_totals = stats::setNames(
        vapply(directions, function(direction) {
          sum(consensus_counts$n[consensus_counts$direction == direction])
        }, integer(1)), directions
      ),
      joint_rank_zero = TRUE,
      empty_joint_payload_respected = TRUE,
      deterministic_order = TRUE,
      parent_isolated = NA,
      in_memory_bytes = NA_real_, serialized_bytes = NA_real_,
      max_target_bytes = max_target_bytes
    )
  )
  stopifnot(!integration_contains_parent(out))
  out$audit$parent_isolated <- TRUE
  out$audit$in_memory_bytes <- as.numeric(object.size(out))
  out$audit$serialized_bytes <- as.numeric(length(qs2::qs_serialize(out)))
  stopifnot(out$audit$in_memory_bytes <= max_target_bytes,
            out$audit$serialized_bytes <= max_target_bytes)
  out
}
