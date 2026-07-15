# Differential expression for the three non-snRNAseq modalities (GeoMx spatial, 24M bulk
# proteome, 24M bulk phosphosite), restored lean from the pre-teardown P4 crossmodality /
# mechanism modules. SCOPE = the PRIMARY per-contrast topTables only: each producer returns
# `$top[[contrast]]` (or `$primary$top[[contrast]]` for GeoMx) with a `logFC` column keyed by
# the 5 canonical contrasts. The amyloid-response logFC scatter (R/report/figures.R ->
# modality_logfc_scatter_data) reads `nlgf_in_maptki` (y) vs `nlgf_in_p301s` (x) from these.
# The torn-down auxiliary arms (GeoMx unblocked / bio-unit-collapsed sensitivities +
# SpatialDecon abundance fitting; proteome / phospho additive run-index sensitivity) stay
# absent. The only GeoMx modality-native descriptor in the live report is the compact
# DAM-gene AOI sample heatmap. Shared machinery reused from HEAD: fit_limma_log /
# median_normalise / prevalence_filter (R/analysis/de_pb.R), factorial_design / make_contrast_matrix
# (R/core/design.R), match_intensity_columns / load_geomx / read_spectronaut_tsv /
# proteomics_sample_meta (R/core/io.R). All calls namespace-qualified so the file sources cleanly
# into any session.

# --- GeoMx spatial DE (RNA counts; slide fixed effect + bio-unit duplicateCorrelation) -------

geomx_required_meta_cols <- function() {
  c("genotype", "slide_rep", "bio_rep", "roi", "SampleID",
    "segment", "area", "ROI Coordinate X", "ROI Coordinate Y",
    "q_norm_qFactors", "NegGeoMean_Mm_R_NGS_WTA_v1.0", "nuclei")
}

# Pull the explicit RNA counts layer (the GeoMx object default assay is SCT), drop all-zero
# genes, and coerce to integer only when the layer is already integer-valued (records residues
# otherwise). Provenance rides on attr() for the report.
geomx_count_matrix <- function(geomx, assay = "RNA", layer = "counts", integer_tol = 1e-6) {
  stopifnot(assay %in% names(geomx@assays), nzchar(layer), is.numeric(integer_tol), integer_tol >= 0)
  cnt <- SeuratObject::GetAssayData(geomx, assay = assay, layer = layer)
  cnt <- as.matrix(cnt)
  stopifnot(!is.null(rownames(cnt)), !is.null(colnames(cnt)),
            all(is.finite(cnt)), all(cnt >= 0))
  empty <- Matrix::rowSums(cnt) <= 0
  cnt <- cnt[!empty, , drop = FALSE]
  residue <- abs(cnt - round(cnt))
  max_residue <- if (length(residue)) max(residue) else 0
  n_non_integer <- sum(residue > integer_tol)
  coerced_integer <- identical(n_non_integer, 0L)
  if (coerced_integer) {
    cnt <- round(cnt)
    storage.mode(cnt) <- "integer"
  } else {
    storage.mode(cnt) <- "double"
  }
  attr(cnt, "geomx_count_provenance") <- list(
    assay = assay,
    layer = layer,
    n_genes_input = nrow(cnt) + sum(empty),
    n_genes_dropped_empty = sum(empty),
    n_aoi = ncol(cnt),
    integer_tol = integer_tol,
    max_abs_integer_residue = max_residue,
    n_non_integer = n_non_integer,
    coerced_integer = coerced_integer
  )
  cnt
}

# AOI metadata -> the replicate-aware design covariates. bio_unit = genotype:bio_rep is the
# duplicateCorrelation block (repeated AOIs within a biological unit). Fails loud on any missing
# column or non-finite covariate (a dropped AOI would corrupt the design row shape).
geomx_meta <- function(geomx) {
  md <- geomx@meta.data
  missing <- setdiff(geomx_required_meta_cols(), names(md))
  if (length(missing)) {
    stop("GeoMx metadata missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  stopifnot(identical(rownames(md), colnames(geomx)))
  out <- data.frame(
    genotype = factor(as.character(md$genotype), levels = genotype_levels),
    slide = factor(as.character(md$slide_rep)),
    bio_rep = as.character(md$bio_rep),
    roi = as.character(md$roi),
    segment = as.character(md$segment),
    SampleID = as.character(md$SampleID),
    area = as.numeric(md$area),
    x = as.numeric(md[["ROI Coordinate X"]]),
    y = as.numeric(md[["ROI Coordinate Y"]]),
    q3_factor = as.numeric(md$q_norm_qFactors),
    neg_background = as.numeric(md[["NegGeoMean_Mm_R_NGS_WTA_v1.0"]]),
    nuclei = as.numeric(md$nuclei),
    row.names = rownames(md),
    stringsAsFactors = FALSE
  )
  out$bio_unit <- paste(out$genotype, out$bio_rep, sep = ":")
  stopifnot(!anyNA(out$genotype), !anyNA(out$slide), !anyNA(out$bio_unit),
            !anyNA(out$roi), !anyNA(out$segment), !anyNA(out$SampleID),
            all(is.finite(out$area)), all(out$area > 0),
            all(is.finite(out$x)), all(is.finite(out$y)),
            all(is.finite(out$q3_factor)), all(is.finite(out$neg_background)),
            all(is.finite(out$nuclei)))
  out$genotype <- factor(out$genotype, levels = genotype_levels)
  out$bio_unit <- factor(out$bio_unit)
  attr(out, "geomx_meta_provenance") <- list(
    n_aoi = nrow(out),
    n_bio_units = nlevels(out$bio_unit),
    n_slides = nlevels(out$slide),
    segments = sort(unique(out$segment)),
    area_range = range(out$area),
    nuclei_sentinel_count = sum(out$nuclei < 0)
  )
  attr(out, "geomx_raw_meta") <- md
  out
}

# Cell-means design (~ 0 + genotype [+ slide]); genotype columns renamed to bare levels so
# make_contrast_matrix's makeContrasts sees valid R names. Batch-free modality -> slide is the
# only nuisance term.
geomx_slide_design <- function(meta, include_slide = TRUE) {
  stopifnot(is.data.frame(meta), all(c("genotype", "slide") %in% names(meta)))
  df <- data.frame(
    genotype = factor(as.character(meta$genotype), levels = genotype_levels),
    slide = factor(as.character(meta$slide)),
    row.names = rownames(meta)
  )
  stopifnot(!anyNA(df$genotype), !anyNA(df$slide))
  if (include_slide) {
    if (nlevels(df$slide) < 2L) stop("GeoMx slide design needs >=2 slide levels", call. = FALSE)
    design <- stats::model.matrix(~ 0 + genotype + slide, data = df)
  } else {
    design <- stats::model.matrix(~ 0 + genotype, data = df)
  }
  colnames(design) <- sub("^genotype", "", colnames(design))
  if (qr(design)$rank != ncol(design)) {
    stop("GeoMx design is rank-deficient", call. = FALSE)
  }
  list(design = design, contrasts = make_contrast_matrix(design))
}

.geomx_top_tables <- function(fit, contrasts) {
  out <- lapply(colnames(contrasts), function(cn) {
    tt <- tibble::rownames_to_column(
      limma::topTable(fit, coef = cn, number = Inf, sort.by = "none", confint = TRUE),
      "symbol"
    )
    tt$contrast <- cn
    tt[, c("symbol", "contrast", "logFC", "P.Value", "adj.P.Val", "t", "CI.L", "CI.R",
           setdiff(names(tt), c("symbol", "contrast", "logFC", "P.Value", "adj.P.Val", "t", "CI.L", "CI.R"))),
       drop = FALSE]
  })
  names(out) <- colnames(contrasts)
  out
}

geomx_spatial_descriptor <- function(counts, meta, top,
                                     y_contrast = "nlgf_in_maptki",
                                     x_contrast = "nlgf_in_p301s",
                                     top_n = 24L) {
  stopifnot(is.matrix(counts), is.data.frame(meta), is.list(top),
            identical(colnames(counts), rownames(meta)),
            all(c(y_contrast, x_contrast) %in% names(top)),
            is.numeric(top_n), length(top_n) == 1L, top_n >= 1L)
  need_meta <- c("slide", "roi", "segment", "SampleID", "genotype", "area", "x", "y",
                 "q3_factor", "neg_background", "nuclei")
  missing_meta <- setdiff(need_meta, names(meta))
  if (length(missing_meta)) {
    stop("GeoMx spatial descriptor metadata missing columns: ",
         paste(missing_meta, collapse = ", "), call. = FALSE)
  }
  ty <- top[[y_contrast]]
  tx <- top[[x_contrast]]
  need <- c("symbol", "logFC", "adj.P.Val")
  missing_y <- setdiff(need, names(ty))
  missing_x <- setdiff(need, names(tx))
  if (length(missing_y) || length(missing_x)) {
    stop("GeoMx top tables missing columns for spatial descriptor", call. = FALSE)
  }
  sy <- as.character(ty$symbol)
  sx <- as.character(tx$symbol)
  stopifnot(anyDuplicated(sy) == 0L, anyDuplicated(sx) == 0L,
            length(sy) == length(sx), setequal(sy, sx))
  idx <- match(sy, sx)
  rank <- data.frame(
    symbol = sy,
    y = as.numeric(ty$logFC),
    x = as.numeric(tx$logFC)[idx],
    fdr_y = as.numeric(ty$adj.P.Val),
    fdr_x = as.numeric(tx$adj.P.Val)[idx],
    stringsAsFactors = FALSE
  )
  rank <- rank[rank$symbol %in% rownames(counts) &
                 is.finite(rank$y) & is.finite(rank$x), , drop = FALSE]
  if (!nrow(rank)) stop("no GeoMx spatial descriptor genes overlap counts", call. = FALSE)
  min_fdr <- pmin(rank$fdr_y, rank$fdr_x)
  min_fdr[!is.finite(min_fdr)] <- 1
  rank$mean_effect <- rowMeans(rank[, c("y", "x"), drop = FALSE])
  rank$rank_score <- pmax(abs(rank$y), abs(rank$x)) * -log10(pmax(min_fdr, 1e-300))
  rank <- rank[order(-rank$rank_score, -abs(rank$mean_effect), rank$symbol,
                     method = "radix"), , drop = FALSE]
  rank <- rank[!duplicated(rank$symbol), , drop = FALSE]
  genes <- utils::head(rank$symbol, as.integer(top_n))

  dge <- edgeR::normLibSizes(edgeR::DGEList(counts = counts))
  logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
  mat <- logcpm[genes, , drop = FALSE]
  row_mean <- rowMeans(mat)
  row_sd <- apply(mat, 1L, stats::sd)
  row_sd[!is.finite(row_sd) | row_sd <= 0] <- 1
  z <- sweep(sweep(mat, 1L, row_mean, "-"), 1L, row_sd, "/")
  gene_direction <- sign(rank$mean_effect[match(rownames(z), rank$symbol)])
  gene_direction[!is.finite(gene_direction) | gene_direction == 0] <- 1
  signed_z <- z * gene_direction
  score <- colMeans(signed_z)

  aoi <- data.frame(
    aoi = rownames(meta),
    slide = factor(as.character(meta$slide), levels = sort(unique(as.character(meta$slide)))),
    roi = as.character(meta$roi),
    segment = as.character(meta$segment),
    sample_id = as.character(meta$SampleID),
    genotype = factor(as.character(meta$genotype), levels = genotype_levels),
    aoi_area = as.numeric(meta$area),
    x_coord = as.numeric(meta$x),
    y_coord = as.numeric(meta$y),
    q3_factor = as.numeric(meta$q3_factor),
    neg_background = as.numeric(meta$neg_background),
    nuclei = as.numeric(meta$nuclei),
    signed_response_score = as.numeric(score[rownames(meta)]),
    stringsAsFactors = FALSE
  )
  aoi$score_abs <- abs(aoi$signed_response_score)
  aoi$amyloid <- factor(ifelse(startsWith(as.character(aoi$genotype), "NLGF"),
                               "NLGF+", "NLGF-"),
                        levels = c("NLGF-", "NLGF+"))
  stopifnot(!anyNA(aoi$genotype), all(is.finite(aoi$x_coord)), all(is.finite(aoi$y_coord)),
            !anyNA(aoi$segment), all(is.finite(aoi$aoi_area)), all(aoi$aoi_area > 0),
            all(is.finite(aoi$signed_response_score)), all(is.finite(aoi$score_abs)))

  selected <- rank[match(genes, rank$symbol), , drop = FALSE]
  rownames(selected) <- NULL
  list(
    aoi = aoi,
    genes = selected,
    provenance = list(
      y_contrast = y_contrast,
      x_contrast = x_contrast,
      score = "mean signed AOI z-score across top GeoMx amyloid-response genes; sign follows the mean of the two amyloid contrasts",
      n_score_genes = nrow(selected),
      top_genes = paste(selected$symbol, collapse = ", ")
    )
  )
}

geomx_contrast_diagnostic_descriptor <- function(top, alpha = 0.10,
                                                 labels_per_contrast = 4L,
                                                 interaction_top_n = 12L) {
  contrast_order <- c("tau_alone", "nlgf_in_maptki", "nlgf_in_p301s",
                      "tau_in_nlgf", "interaction")
  contrast_labels <- c(
    tau_alone = "Tau alone\nP301S - MAPTKI",
    nlgf_in_maptki = "Amyloid | MAPTKI\nNLGF_MAPTKI - MAPTKI",
    nlgf_in_p301s = "Amyloid | P301S\nNLGF_P301S - P301S",
    tau_in_nlgf = "Tau | amyloid\nNLGF_P301S - NLGF_MAPTKI",
    interaction = "Interaction\namyloid x tau"
  )
  stopifnot(is.list(top), all(contrast_order %in% names(top)),
            is.numeric(alpha), length(alpha) == 1L, alpha > 0, alpha < 1,
            is.numeric(labels_per_contrast), length(labels_per_contrast) == 1L,
            labels_per_contrast >= 1L,
            is.numeric(interaction_top_n), length(interaction_top_n) == 1L,
            interaction_top_n >= 1L)

  build_one <- function(cn) {
    tt <- top[[cn]]
    need <- c("symbol", "logFC", "P.Value", "adj.P.Val", "AveExpr", "CI.L", "CI.R")
    missing <- setdiff(need, names(tt))
    if (length(missing)) {
      stop("GeoMx top table ", cn, " missing columns: ",
           paste(missing, collapse = ", "), call. = FALSE)
    }
    symbol <- as.character(tt$symbol)
    if (anyDuplicated(symbol)) {
      stop("GeoMx top table ", cn, " has duplicated symbols", call. = FALSE)
    }
    out <- data.frame(
      symbol = symbol,
      contrast = cn,
      contrast_label = unname(contrast_labels[[cn]]),
      effect = as.numeric(tt$logFC),
      ave_expr = as.numeric(tt$AveExpr),
      p_value = as.numeric(tt$P.Value),
      fdr = as.numeric(tt$adj.P.Val),
      ci_low = as.numeric(tt$CI.L),
      ci_high = as.numeric(tt$CI.R),
      stringsAsFactors = FALSE
    )
    out <- out[is.finite(out$effect) & is.finite(out$ave_expr) &
                 is.finite(out$p_value) & is.finite(out$fdr) &
                 is.finite(out$ci_low) & is.finite(out$ci_high), ,
               drop = FALSE]
    if (!nrow(out)) stop("GeoMx contrast ", cn, " has no finite diagnostic rows",
                         call. = FALSE)
    out$neg_log10_fdr <- -log10(pmax(out$fdr, 1e-300))
    out$supported <- out$fdr <= alpha
    out$direction <- ifelse(
      out$supported,
      ifelse(out$effect < 0, "down", "up"),
      "not supported"
    )
    out$abs_effect <- abs(out$effect)
    out$rank_score <- out$neg_log10_fdr * sqrt(pmax(out$abs_effect,
                                                    .Machine$double.eps))
    out
  }

  volcano <- do.call(rbind, lapply(contrast_order, build_one))
  rownames(volcano) <- NULL
  volcano$contrast <- factor(volcano$contrast, levels = contrast_order)
  volcano$contrast_label <- factor(volcano$contrast_label,
                                   levels = unname(contrast_labels[contrast_order]))
  volcano$direction <- factor(volcano$direction,
                              levels = c("down", "not supported", "up"))
  stopifnot(all(is.finite(volcano$effect)), all(is.finite(volcano$ave_expr)),
            all(is.finite(volcano$neg_log10_fdr)), !anyNA(volcano$contrast),
            !anyNA(volcano$contrast_label), !anyNA(volcano$direction))

  label_rows <- do.call(rbind, lapply(contrast_order, function(cn) {
    x <- volcano[as.character(volcano$contrast) == cn, , drop = FALSE]
    x <- x[order(!x$supported, -x$rank_score, -x$neg_log10_fdr,
                 -x$abs_effect, x$symbol, method = "radix"), , drop = FALSE]
    utils::head(x, as.integer(labels_per_contrast))
  }))
  label_key <- paste(as.character(label_rows$contrast), label_rows$symbol, sep = "\r")
  volcano$label_show <- paste(as.character(volcano$contrast), volcano$symbol, sep = "\r") %in%
    label_key
  volcano$label <- ifelse(volcano$label_show, volcano$symbol, "")
  label_rows <- volcano[volcano$label_show, , drop = FALSE]
  rownames(label_rows) <- NULL

  support_grid <- expand.grid(
    contrast = contrast_order,
    direction = c("down", "up"),
    stringsAsFactors = FALSE
  )
  support <- volcano[volcano$supported & as.character(volcano$direction) %in% c("down", "up"),
                     c("contrast", "direction"), drop = FALSE]
  support$n <- 1L
  if (nrow(support)) {
    support_count <- stats::aggregate(n ~ contrast + direction, data = support, FUN = sum)
  } else {
    support_count <- data.frame(contrast = character(), direction = character(),
                                n = integer(), stringsAsFactors = FALSE)
  }
  support_counts <- merge(support_grid, support_count,
                          by = c("contrast", "direction"), all.x = TRUE, sort = FALSE)
  support_counts$n[is.na(support_counts$n)] <- 0L
  support_counts$contrast_label <- unname(contrast_labels[support_counts$contrast])
  support_counts$signed_n <- ifelse(support_counts$direction == "down",
                                    -support_counts$n, support_counts$n)
  support_counts$contrast <- factor(support_counts$contrast, levels = contrast_order)
  support_counts$contrast_label <- factor(support_counts$contrast_label,
                                          levels = unname(contrast_labels[contrast_order]))
  support_counts$direction <- factor(support_counts$direction, levels = c("down", "up"))
  support_counts <- support_counts[order(support_counts$contrast,
                                         support_counts$direction, method = "radix"), ,
                                   drop = FALSE]
  rownames(support_counts) <- NULL

  interaction_top <- volcano[as.character(volcano$contrast) == "interaction", ,
                             drop = FALSE]
  interaction_top <- interaction_top[order(!interaction_top$supported,
                                           -interaction_top$rank_score,
                                           -interaction_top$neg_log10_fdr,
                                           -interaction_top$abs_effect,
                                           interaction_top$symbol,
                                           method = "radix"), , drop = FALSE]
  interaction_top <- utils::head(interaction_top, as.integer(interaction_top_n))
  interaction_top <- interaction_top[order(interaction_top$effect,
                                           interaction_top$symbol,
                                           method = "radix"), , drop = FALSE]
  interaction_top$symbol_plot <- factor(interaction_top$symbol,
                                        levels = interaction_top$symbol)
  rownames(interaction_top) <- NULL

  summary <- do.call(rbind, lapply(contrast_order, function(cn) {
    x <- volcano[as.character(volcano$contrast) == cn, , drop = FALSE]
    data.frame(
      contrast = cn,
      contrast_label = unname(contrast_labels[[cn]]),
      n_features = nrow(x),
      n_supported = sum(x$supported),
      n_up = sum(x$supported & as.character(x$direction) == "up"),
      n_down = sum(x$supported & as.character(x$direction) == "down"),
      stringsAsFactors = FALSE
    )
  }))
  summary$contrast <- factor(summary$contrast, levels = contrast_order)
  summary$contrast_label <- factor(summary$contrast_label,
                                   levels = unname(contrast_labels[contrast_order]))
  rownames(summary) <- NULL

  list(
    volcano = volcano,
    labels = label_rows,
    support_counts = support_counts,
    interaction_top = interaction_top,
    summary = summary,
    provenance = list(
      alpha = alpha,
      contrast_order = contrast_order,
      contrast_labels = contrast_labels[contrast_order],
      labels_per_contrast = as.integer(labels_per_contrast),
      interaction_top_n = as.integer(interaction_top_n),
      n_features_per_contrast = stats::setNames(
        vapply(contrast_order, function(cn) nrow(top[[cn]]), integer(1)),
        contrast_order
      ),
      display = "volcano and MA diagnostics from existing GeoMx primary topTables; no new model or exclusion"
    )
  )
}

geomx_roi_replicate_descriptor <- function(counts, meta, design, duplicate_correlation,
                                           spatial_aoi, min_count = 5,
                                           n_variable_features = 2000L) {
  stopifnot(is.matrix(counts), is.data.frame(meta), is.matrix(design),
            is.list(duplicate_correlation), is.data.frame(spatial_aoi),
            identical(colnames(counts), rownames(meta)),
            identical(colnames(counts), rownames(design)),
            qr(design)$rank == ncol(design),
            is.numeric(min_count), length(min_count) == 1L, min_count >= 0,
            is.numeric(n_variable_features), length(n_variable_features) == 1L,
            n_variable_features >= 2L,
            anyDuplicated(rownames(counts)) == 0L)
  need <- c("slide", "roi", "segment", "genotype", "SampleID", "bio_unit", "bio_rep")
  missing <- setdiff(need, names(meta))
  if (length(missing)) {
    stop("GeoMx ROI-replicate metadata missing columns: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  need_spatial <- c("aoi", "signed_response_score")
  missing_spatial <- setdiff(need_spatial, names(spatial_aoi))
  if (length(missing_spatial)) {
    stop("GeoMx ROI-replicate spatial AOI data missing columns: ",
         paste(missing_spatial, collapse = ", "), call. = FALSE)
  }
  stopifnot(!is.null(duplicate_correlation$used),
            !is.null(duplicate_correlation$consensus_correlation))

  meta_frame <- data.frame(
    aoi = rownames(meta),
    genotype = factor(as.character(meta$genotype), levels = genotype_levels),
    bio_rep = as.character(meta$bio_rep),
    bio_unit = factor(as.character(meta$bio_unit)),
    slide = factor(as.character(meta$slide), levels = sort(unique(as.character(meta$slide)))),
    roi = as.character(meta$roi),
    segment = factor(as.character(meta$segment),
                     levels = sort(unique(as.character(meta$segment)))),
    sample_id = as.character(meta$SampleID),
    stringsAsFactors = FALSE
  )
  meta_frame$signed_response_score <- as.numeric(
    spatial_aoi$signed_response_score[match(meta_frame$aoi, spatial_aoi$aoi)]
  )
  stopifnot(!anyNA(meta_frame$genotype), !anyNA(meta_frame$bio_rep),
            !anyNA(meta_frame$bio_unit), !anyNA(meta_frame$slide),
            !anyNA(meta_frame$segment),
            all(is.finite(meta_frame$signed_response_score)))

  bio_reps <- sort(unique(meta_frame$bio_rep))
  support_grid <- expand.grid(
    genotype = genotype_levels,
    bio_rep = bio_reps,
    stringsAsFactors = FALSE
  )
  support_obs <- stats::aggregate(
    list(n_aoi = meta_frame$aoi, n_roi = meta_frame$roi,
         n_segments = meta_frame$segment),
    by = list(genotype = as.character(meta_frame$genotype),
              bio_rep = meta_frame$bio_rep),
    FUN = function(x) length(unique(x))
  )
  slide_by_unit <- stats::aggregate(
    list(slide = as.character(meta_frame$slide)),
    by = list(genotype = as.character(meta_frame$genotype),
              bio_rep = meta_frame$bio_rep),
    FUN = function(x) paste(sort(unique(x)), collapse = "+")
  )
  support_obs <- merge(support_obs, slide_by_unit,
                       by = c("genotype", "bio_rep"), all.x = TRUE,
                       sort = FALSE)
  support <- merge(support_grid, support_obs, by = c("genotype", "bio_rep"),
                   all.x = TRUE, sort = FALSE)
  support$n_aoi[is.na(support$n_aoi)] <- 0L
  support$n_roi[is.na(support$n_roi)] <- 0L
  support$n_segments[is.na(support$n_segments)] <- 0L
  support$slide[is.na(support$slide)] <- "absent"
  support$present <- support$n_aoi > 0
  support$label <- ifelse(support$present, as.character(support$n_aoi), "0")
  support$genotype <- factor(support$genotype, levels = genotype_levels)
  support$bio_rep <- factor(support$bio_rep, levels = bio_reps)
  rownames(support) <- NULL

  block <- do.call(rbind, lapply(split(meta_frame, meta_frame$bio_unit), function(x) {
    data.frame(
      bio_unit = as.character(x$bio_unit[[1L]]),
      genotype = x$genotype[[1L]],
      bio_rep = x$bio_rep[[1L]],
      slide = x$slide[[1L]],
      n_aoi = nrow(x),
      n_roi = length(unique(x$roi)),
      n_segments = length(unique(x$segment)),
      score_mean = mean(x$signed_response_score),
      score_sd = if (nrow(x) > 1L) stats::sd(x$signed_response_score) else 0,
      score_min = min(x$signed_response_score),
      score_max = max(x$signed_response_score),
      stringsAsFactors = FALSE
    )
  }))
  block$genotype <- factor(as.character(block$genotype), levels = genotype_levels)
  block$slide <- factor(as.character(block$slide), levels = levels(meta_frame$slide))
  block <- block[order(block$genotype, as.integer(block$bio_rep), block$bio_unit,
                       method = "radix"), , drop = FALSE]
  block$block_rank <- seq_len(nrow(block))
  block$bio_unit_plot <- factor(block$bio_unit, levels = rev(block$bio_unit))
  rownames(block) <- NULL
  stopifnot(!anyNA(block$genotype), !anyNA(block$slide),
            all(is.finite(block$n_aoi)), all(block$n_aoi >= 1L),
            all(is.finite(block$n_roi)), all(is.finite(block$n_segments)),
            all(is.finite(block$score_mean)), all(is.finite(block$score_sd)))

  dge0 <- edgeR::DGEList(counts = counts)
  keep <- edgeR::filterByExpr(dge0, design = design, min.count = min_count)
  if (!any(keep)) stop("no GeoMx features passed filterByExpr for ROI replicate audit",
                       call. = FALSE)
  dge <- edgeR::normLibSizes(dge0[keep, , keep.lib.sizes = FALSE], method = "TMM")
  logcpm <- edgeR::cpm(dge, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 1)
  stopifnot(identical(colnames(logcpm), rownames(meta)), all(is.finite(logcpm)))
  feature_var <- apply(logcpm, 1L, stats::var)
  var_rank <- data.frame(symbol = rownames(logcpm),
                         variance = as.numeric(feature_var),
                         stringsAsFactors = FALSE)
  var_rank <- var_rank[is.finite(var_rank$variance) & var_rank$variance > 0, ,
                       drop = FALSE]
  if (nrow(var_rank) < 2L) {
    stop("GeoMx ROI replicate audit has too few variable filter-passing genes",
         call. = FALSE)
  }
  var_rank <- var_rank[order(-var_rank$variance, var_rank$symbol, method = "radix"), ,
                       drop = FALSE]
  variable_genes <- utils::head(var_rank$symbol, as.integer(n_variable_features))
  mat <- logcpm[variable_genes, , drop = FALSE]
  sample_cor <- stats::cor(mat, method = "pearson")
  stopifnot(all(is.finite(sample_cor)))
  pairs <- utils::combn(colnames(sample_cor), 2L)
  sample_lookup <- meta_frame
  rownames(sample_lookup) <- sample_lookup$aoi
  pair_correlation <- data.frame(
    aoi_1 = pairs[1L, ],
    aoi_2 = pairs[2L, ],
    correlation = as.numeric(sample_cor[cbind(pairs[1L, ], pairs[2L, ])]),
    stringsAsFactors = FALSE
  )
  pair_correlation$bio_unit_1 <- as.character(sample_lookup[pair_correlation$aoi_1, "bio_unit"])
  pair_correlation$bio_unit_2 <- as.character(sample_lookup[pair_correlation$aoi_2, "bio_unit"])
  pair_correlation$genotype_1 <- as.character(sample_lookup[pair_correlation$aoi_1, "genotype"])
  pair_correlation$genotype_2 <- as.character(sample_lookup[pair_correlation$aoi_2, "genotype"])
  pair_correlation$pair_type <- ifelse(
    pair_correlation$bio_unit_1 == pair_correlation$bio_unit_2,
    "same bio-unit",
    ifelse(pair_correlation$genotype_1 == pair_correlation$genotype_2,
           "same genotype, different bio-unit", "different genotype")
  )
  pair_levels <- c("same bio-unit", "same genotype, different bio-unit",
                   "different genotype")
  pair_correlation$pair_type <- factor(pair_correlation$pair_type, levels = pair_levels)
  pair_correlation <- pair_correlation[is.finite(pair_correlation$correlation), ,
                                       drop = FALSE]
  if (!nrow(pair_correlation)) {
    stop("GeoMx ROI replicate audit has no finite AOI-pair correlations",
         call. = FALSE)
  }

  pair_summary <- do.call(rbind, lapply(split(pair_correlation, pair_correlation$pair_type,
                                              drop = TRUE), function(x) {
    data.frame(
      pair_type = x$pair_type[[1L]],
      n_pairs = nrow(x),
      median_correlation = stats::median(x$correlation),
      q25_correlation = as.numeric(stats::quantile(x$correlation, 0.25, names = FALSE)),
      q75_correlation = as.numeric(stats::quantile(x$correlation, 0.75, names = FALSE)),
      mean_correlation = mean(x$correlation),
      stringsAsFactors = FALSE
    )
  }))
  pair_summary$pair_type <- factor(as.character(pair_summary$pair_type),
                                   levels = pair_levels)
  rownames(pair_summary) <- NULL
  stopifnot(!anyNA(pair_summary$pair_type),
            all(is.finite(pair_summary$n_pairs)),
            all(is.finite(pair_summary$median_correlation)),
            all(is.finite(pair_summary$q25_correlation)),
            all(is.finite(pair_summary$q75_correlation)))

  list(
    support = support,
    block = block,
    pair_correlation = pair_correlation,
    pair_summary = pair_summary,
    provenance = list(
      n_aoi = ncol(counts),
      n_bio_units = nlevels(meta_frame$bio_unit),
      n_expected_bio_units = length(genotype_levels) * length(bio_reps),
      n_present_bio_units = sum(support$present),
      n_missing_bio_units = sum(!support$present),
      n_slides = nlevels(meta_frame$slide),
      segments = levels(meta_frame$segment),
      n_segment_levels = nlevels(meta_frame$segment),
      all_segments_single_level = nlevels(meta_frame$segment) == 1L,
      n_repeated_blocks = sum(block$n_aoi > 1L),
      block_size_range = range(block$n_aoi),
      min_count = min_count,
      n_input_features = nrow(counts),
      n_kept_features = sum(keep),
      n_variable_features = length(variable_genes),
      duplicate_correlation_used = isTRUE(duplicate_correlation$used),
      duplicate_correlation = unname(as.numeric(duplicate_correlation$consensus_correlation)),
      design_rank = qr(design)$rank,
      residual_df = nrow(design) - qr(design)$rank,
      model = "primary GeoMx DE uses slide fixed effect and bio-unit duplicateCorrelation; this audit changes no AOI exclusions or model terms",
      correlation_basis = "Pearson AOI-pair correlations over TMM-logCPM top-variable filterByExpr-kept genes"
    )
  )
}

geomx_decon_feasibility_descriptor <- function(counts, meta, design, spatial_aoi,
                                               min_count = 5,
                                               marker_sets = list(
                                                 Microglia = microglia_identity_markers,
                                                 Homeostatic = canonical_microglia_markers$Homeostatic,
                                                 DAM = canonical_microglia_markers$DAM,
                                                 IFN = canonical_microglia_markers$IFN,
                                                 MHC_APC = canonical_microglia_markers$MHC_APC,
                                                 Astrocyte = contam_signatures$Astro,
                                                 Oligodendrocyte = contam_signatures$Oligo,
                                                 Neuron = contam_signatures$Neuron
                                               )) {
  stopifnot(is.matrix(counts), is.data.frame(meta), is.matrix(design),
            is.data.frame(spatial_aoi),
            identical(colnames(counts), rownames(meta)),
            identical(colnames(counts), rownames(design)),
            qr(design)$rank == ncol(design),
            is.numeric(min_count), length(min_count) == 1L, min_count >= 0,
            is.list(marker_sets), length(marker_sets) >= 1L, !is.null(names(marker_sets)),
            !any(names(marker_sets) == ""),
            anyDuplicated(rownames(counts)) == 0L)
  need <- c("slide", "roi", "segment", "SampleID", "genotype", "area", "x", "y",
            "q3_factor", "neg_background", "nuclei")
  missing <- setdiff(need, names(meta))
  if (length(missing)) {
    stop("GeoMx decon-feasibility metadata missing columns: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  need_spatial <- c("aoi", "signed_response_score")
  missing_spatial <- setdiff(need_spatial, names(spatial_aoi))
  if (length(missing_spatial)) {
    stop("GeoMx decon-feasibility spatial data missing columns: ",
         paste(missing_spatial, collapse = ", "), call. = FALSE)
  }

  marker_sets <- lapply(marker_sets, function(x) unique(trimws(as.character(x))))
  marker_sets <- lapply(marker_sets, function(x) x[!is.na(x) & x != ""])
  if (any(vapply(marker_sets, length, integer(1)) == 0L)) {
    stop("GeoMx decon-feasibility marker set is empty", call. = FALSE)
  }

  dge0 <- edgeR::DGEList(counts = counts)
  keep <- edgeR::filterByExpr(dge0, design = design, min.count = min_count)
  if (!any(keep)) stop("no GeoMx features passed filterByExpr for decon feasibility",
                       call. = FALSE)
  names(keep) <- rownames(counts)
  dge <- edgeR::normLibSizes(dge0, method = "TMM")
  logcpm <- edgeR::cpm(dge, normalized.lib.sizes = TRUE, log = TRUE,
                       prior.count = 1)
  stopifnot(identical(rownames(logcpm), rownames(counts)),
            identical(colnames(logcpm), colnames(counts)),
            all(is.finite(logcpm)))

  genes <- data.frame(
    symbol = rownames(counts),
    mean_count = as.numeric(rowMeans(counts)),
    mean_logcpm = as.numeric(rowMeans(logcpm)),
    detect_fraction_min_count = as.numeric(rowMeans(counts >= min_count)),
    filter_pass = as.logical(keep),
    stringsAsFactors = FALSE
  )
  stopifnot(all(is.finite(genes$mean_count)), all(genes$mean_count >= 0),
            all(is.finite(genes$mean_logcpm)),
            all(is.finite(genes$detect_fraction_min_count)))

  component_levels <- names(marker_sets)
  coverage <- do.call(rbind, lapply(component_levels, function(component) {
    sig <- marker_sets[[component]]
    present <- intersect(sig, rownames(counts))
    passing <- present[keep[present]]
    metrics <- genes[match(present, genes$symbol), , drop = FALSE]
    data.frame(
      component = component,
      n_signature = length(sig),
      n_present = length(present),
      n_filter_passing = length(passing),
      coverage_fraction = length(present) / length(sig),
      filter_fraction = length(passing) / length(sig),
      median_detect_fraction = if (nrow(metrics)) {
        stats::median(metrics$detect_fraction_min_count)
      } else 0,
      median_mean_logcpm = if (nrow(metrics)) stats::median(metrics$mean_logcpm) else 0,
      stringsAsFactors = FALSE
    )
  }))
  coverage$status <- ifelse(coverage$n_filter_passing >= 2L, "covered",
                            ifelse(coverage$n_filter_passing == 1L, "thin", "absent"))
  coverage$component <- factor(coverage$component, levels = component_levels)
  coverage$status <- factor(coverage$status, levels = c("absent", "thin", "covered"))
  rownames(coverage) <- NULL
  stopifnot(all(is.finite(coverage$n_signature)), all(coverage$n_signature >= 1L),
            all(is.finite(coverage$n_present)), all(is.finite(coverage$n_filter_passing)),
            all(is.finite(coverage$coverage_fraction)),
            all(is.finite(coverage$filter_fraction)),
            all(is.finite(coverage$median_detect_fraction)),
            all(is.finite(coverage$median_mean_logcpm)),
            !anyNA(coverage$status))

  score_components <- as.character(coverage$component[coverage$n_filter_passing >= 2L])
  if (!length(score_components)) {
    stop("GeoMx decon-feasibility marker sets have no covered components", call. = FALSE)
  }
  logcpm_keep <- logcpm[keep, , drop = FALSE]
  row_mean <- rowMeans(logcpm_keep)
  row_sd <- apply(logcpm_keep, 1L, stats::sd)
  row_sd[!is.finite(row_sd) | row_sd <= 0] <- 1
  z <- sweep(sweep(logcpm_keep, 1L, row_mean, "-"), 1L, row_sd, "/")
  z[!is.finite(z)] <- 0

  score_list <- lapply(score_components, function(component) {
    present <- intersect(marker_sets[[component]], rownames(z))
    colMeans(z[present, , drop = FALSE])
  })
  names(score_list) <- score_components
  score_matrix <- do.call(cbind, score_list)
  rownames(score_matrix) <- colnames(counts)
  colnames(score_matrix) <- score_components
  stopifnot(all(is.finite(score_matrix)), identical(rownames(score_matrix), colnames(counts)))

  resid_ss <- stats::setNames(rep(0, ncol(counts)), colnames(counts))
  resid_n <- stats::setNames(rep(0L, ncol(counts)), colnames(counts))
  for (component in score_components) {
    present <- intersect(marker_sets[[component]], rownames(z))
    res <- sweep(z[present, , drop = FALSE], 2L, score_matrix[, component], "-")
    resid_ss <- resid_ss + colSums(res^2)
    resid_n <- resid_n + nrow(res)
  }
  marker_residual_rmse <- sqrt(resid_ss / pmax(resid_n, 1L))
  stopifnot(all(is.finite(marker_residual_rmse)))

  library_size <- Matrix::colSums(counts)
  q3_scaled_background <- as.numeric(meta$neg_background) / as.numeric(meta$q3_factor)
  qn <- function(x, q) as.numeric(stats::quantile(as.numeric(x), q, names = FALSE, type = 7))
  thresholds <- c(
    library_size_low = qn(library_size, 0.05),
    area_low = qn(meta$area, 0.05),
    q3_factor_high = qn(meta$q3_factor, 0.95),
    q3_scaled_background_high = qn(q3_scaled_background, 0.95)
  )
  aoi <- data.frame(
    aoi = rownames(meta),
    slide = factor(as.character(meta$slide), levels = sort(unique(as.character(meta$slide)))),
    roi = as.character(meta$roi),
    segment = factor(as.character(meta$segment),
                     levels = sort(unique(as.character(meta$segment)))),
    sample_id = as.character(meta$SampleID),
    genotype = factor(as.character(meta$genotype), levels = genotype_levels),
    aoi_area = as.numeric(meta$area),
    x_coord = as.numeric(meta$x),
    y_coord = as.numeric(meta$y),
    library_size = as.numeric(library_size[rownames(meta)]),
    q3_factor = as.numeric(meta$q3_factor),
    neg_background = as.numeric(meta$neg_background),
    q3_scaled_background = q3_scaled_background,
    nuclei = as.numeric(meta$nuclei),
    signed_response_score = as.numeric(spatial_aoi$signed_response_score[
      match(rownames(meta), spatial_aoi$aoi)
    ]),
    marker_residual_rmse = as.numeric(marker_residual_rmse[rownames(meta)]),
    max_abs_marker_score = as.numeric(apply(abs(score_matrix[rownames(meta), , drop = FALSE]),
                                            1L, max)),
    stringsAsFactors = FALSE
  )
  aoi$nuclei_usable <- aoi$nuclei >= 0
  aoi$low_library <- aoi$library_size <= thresholds[["library_size_low"]] &
    aoi$library_size < max(aoi$library_size)
  aoi$small_area <- aoi$aoi_area <= thresholds[["area_low"]] &
    aoi$aoi_area < max(aoi$aoi_area)
  aoi$high_q3 <- aoi$q3_factor >= thresholds[["q3_factor_high"]] &
    aoi$q3_factor > min(aoi$q3_factor)
  aoi$high_q3_scaled_background <- aoi$q3_scaled_background >=
    thresholds[["q3_scaled_background_high"]] &
    aoi$q3_scaled_background > min(aoi$q3_scaled_background)
  aoi$input_status <- ifelse(
    !aoi$nuclei_usable, "absolute-count blocked",
    ifelse(aoi$high_q3 | aoi$high_q3_scaled_background, "background/Q3 tail",
           ifelse(aoi$low_library | aoi$small_area, "low-input tail",
                  "no local blocker"))
  )
  status_levels <- c("no local blocker", "low-input tail", "background/Q3 tail",
                     "absolute-count blocked")
  aoi$input_status <- factor(aoi$input_status, levels = status_levels)
  rownames(aoi) <- NULL
  stopifnot(!anyNA(aoi$slide), !anyNA(aoi$segment), !anyNA(aoi$genotype),
            all(is.finite(aoi$aoi_area)), all(aoi$aoi_area > 0),
            all(is.finite(aoi$x_coord)), all(is.finite(aoi$y_coord)),
            all(is.finite(aoi$library_size)), all(aoi$library_size >= 0),
            all(is.finite(aoi$q3_factor)), all(aoi$q3_factor > 0),
            all(is.finite(aoi$neg_background)), all(is.finite(aoi$q3_scaled_background)),
            all(is.finite(aoi$nuclei)), all(is.finite(aoi$signed_response_score)),
            all(is.finite(aoi$marker_residual_rmse)),
            all(is.finite(aoi$max_abs_marker_score)), !anyNA(aoi$input_status))

  count_grid <- expand.grid(
    genotype = genotype_levels,
    input_status = status_levels,
    stringsAsFactors = FALSE
  )
  count_obs <- stats::aggregate(
    list(n_aoi = aoi$aoi),
    by = list(genotype = as.character(aoi$genotype),
              input_status = as.character(aoi$input_status)),
    FUN = length
  )
  status_counts <- merge(count_grid, count_obs, by = c("genotype", "input_status"),
                         all.x = TRUE, sort = FALSE)
  status_counts$n_aoi[is.na(status_counts$n_aoi)] <- 0L
  status_counts$genotype <- factor(status_counts$genotype, levels = genotype_levels)
  status_counts$input_status <- factor(status_counts$input_status, levels = status_levels)
  status_counts$n_aoi <- as.integer(status_counts$n_aoi)
  rownames(status_counts) <- NULL
  stopifnot(!anyNA(status_counts$genotype), !anyNA(status_counts$input_status),
            all(is.finite(status_counts$n_aoi)), all(status_counts$n_aoi >= 0))

  list(
    coverage = coverage,
    aoi = aoi,
    status_counts = status_counts,
    thresholds = thresholds,
    provenance = list(
      n_aoi = ncol(counts),
      n_input_features = nrow(counts),
      n_kept_features = sum(keep),
      min_count = min_count,
      n_components = length(component_levels),
      n_score_components = length(score_components),
      score_components = score_components,
      n_nuclei_sentinel = sum(!aoi$nuclei_usable),
      n_background_q3_tail = sum(aoi$input_status == "background/Q3 tail"),
      n_low_input_tail = sum(aoi$input_status == "low-input tail"),
      spatialdecon_dependency_declared = FALSE,
      spatialdecon_live_targets = FALSE,
      live_status = "blocked diagnostic: current lean DAG has no SpatialDecon dependency, reference-profile, beta, or abundance-DE target; this descriptor reports preconditions and marker-coherence residuals only",
      coverage_basis = "candidate decon marker-set overlap with GeoMx WTA genes and the existing filterByExpr-kept feature set",
      residual_basis = "per-AOI RMSE after each covered marker gene z-score is reconstructed by its component mean score; proxy fit QC only, not a SpatialDecon residual",
      input_status_rule = "absolute-count blocked = nuclei sentinel; background/Q3 tail = top 5% Q3 factor or Q3-scaled negative background; low-input tail = bottom 5% library size or AOI area"
    )
  )
}

geomx_qc_descriptor <- function(counts, meta, lower_q = 0.05, upper_q = 0.95) {
  stopifnot(is.matrix(counts), is.data.frame(meta),
            identical(colnames(counts), rownames(meta)),
            is.numeric(lower_q), length(lower_q) == 1L, lower_q > 0, lower_q < 0.5,
            is.numeric(upper_q), length(upper_q) == 1L, upper_q > 0.5, upper_q < 1)
  need <- c("slide", "roi", "segment", "SampleID", "genotype", "area",
            "q3_factor", "neg_background", "nuclei")
  missing <- setdiff(need, names(meta))
  if (length(missing)) {
    stop("GeoMx QC descriptor metadata missing columns: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  library_size <- Matrix::colSums(counts)
  detected_genes <- Matrix::colSums(counts > 0)
  qn <- function(x, q) {
    as.numeric(stats::quantile(as.numeric(x), probs = q, names = FALSE, type = 7))
  }
  aoi <- data.frame(
    aoi = rownames(meta),
    slide = factor(as.character(meta$slide), levels = sort(unique(as.character(meta$slide)))),
    roi = as.character(meta$roi),
    segment = factor(as.character(meta$segment), levels = sort(unique(as.character(meta$segment)))),
    sample_id = as.character(meta$SampleID),
    genotype = factor(as.character(meta$genotype), levels = genotype_levels),
    library_size = as.numeric(library_size[rownames(meta)]),
    detected_genes = as.numeric(detected_genes[rownames(meta)]),
    nuclei = as.numeric(meta$nuclei),
    aoi_area = as.numeric(meta$area),
    neg_background = as.numeric(meta$neg_background),
    q3_factor = as.numeric(meta$q3_factor),
    stringsAsFactors = FALSE
  )
  stopifnot(!anyNA(aoi$genotype), !anyNA(aoi$slide), !anyNA(aoi$segment),
            all(is.finite(aoi$library_size)), all(aoi$library_size >= 0),
            all(is.finite(aoi$detected_genes)), all(aoi$detected_genes >= 0),
            all(is.finite(aoi$nuclei)), all(is.finite(aoi$aoi_area)),
            all(aoi$aoi_area > 0), all(is.finite(aoi$neg_background)),
            all(is.finite(aoi$q3_factor)))

  thresholds <- c(
    library_size_low = qn(aoi$library_size, lower_q),
    detected_genes_low = qn(aoi$detected_genes, lower_q),
    nuclei_sentinel_lt = 0,
    aoi_area_low = qn(aoi$aoi_area, lower_q),
    neg_background_high = qn(aoi$neg_background, upper_q),
    q3_factor_high = qn(aoi$q3_factor, upper_q)
  )
  low_flag <- function(x, threshold) {
    x <- as.numeric(x)
    x <= threshold & x < max(x)
  }
  high_flag <- function(x, threshold) {
    x <- as.numeric(x)
    x >= threshold & x > min(x)
  }
  flags <- data.frame(
    aoi = aoi$aoi,
    low_library = low_flag(aoi$library_size, thresholds[["library_size_low"]]),
    low_detected_genes = low_flag(aoi$detected_genes, thresholds[["detected_genes_low"]]),
    nuclei_sentinel = aoi$nuclei < thresholds[["nuclei_sentinel_lt"]],
    small_area = low_flag(aoi$aoi_area, thresholds[["aoi_area_low"]]),
    high_background = high_flag(aoi$neg_background, thresholds[["neg_background_high"]]),
    high_q3_factor = high_flag(aoi$q3_factor, thresholds[["q3_factor_high"]]),
    stringsAsFactors = FALSE
  )
  flags$flag_any <- rowSums(flags[setdiff(names(flags), "aoi")]) > 0
  aoi$flag_any <- flags$flag_any[match(aoi$aoi, flags$aoi)]

  metric_map <- data.frame(
    metric = c("library_size", "detected_genes", "nuclei", "aoi_area",
               "neg_background", "q3_factor"),
    metric_label = c("Library size", "Detected genes", "Nuclei", "AOI area",
                     "Negative background", "Q3 factor"),
    flag = c("low_library", "low_detected_genes", "nuclei_sentinel", "small_area",
             "high_background", "high_q3_factor"),
    stringsAsFactors = FALSE
  )
  metrics <- do.call(rbind, lapply(seq_len(nrow(metric_map)), function(i) {
    m <- metric_map$metric[[i]]
    f <- metric_map$flag[[i]]
    data.frame(
      aoi = aoi$aoi,
      slide = aoi$slide,
      roi = aoi$roi,
      segment = aoi$segment,
      genotype = aoi$genotype,
      metric = m,
      metric_label = metric_map$metric_label[[i]],
      value = as.numeric(aoi[[m]]),
      flag_metric = flags[[f]],
      flag_any = aoi$flag_any,
      stringsAsFactors = FALSE
    )
  }))
  metrics$metric_label <- factor(metrics$metric_label, levels = metric_map$metric_label)
  rownames(metrics) <- NULL
  stopifnot(nrow(metrics) == nrow(aoi) * nrow(metric_map),
            all(is.finite(metrics$value)), !anyNA(metrics$metric_label))

  flag_labels <- c(
    low_library = "low library",
    low_detected_genes = "low genes",
    nuclei_sentinel = "nuclei sentinel",
    small_area = "small area",
    high_background = "high background",
    high_q3_factor = "high q3"
  )
  flag_long <- do.call(rbind, lapply(names(flag_labels), function(f) {
    data.frame(
      slide = aoi$slide,
      segment = aoi$segment,
      flag = flag_labels[[f]],
      n = as.integer(flags[[f]]),
      stringsAsFactors = FALSE
    )
  }))
  count <- stats::aggregate(n ~ slide + segment + flag, data = flag_long, FUN = sum)
  grid <- expand.grid(
    slide = levels(aoi$slide),
    segment = levels(aoi$segment),
    flag = unname(flag_labels),
    stringsAsFactors = FALSE
  )
  flag_counts <- merge(grid, count, by = c("slide", "segment", "flag"), all.x = TRUE,
                       sort = FALSE)
  flag_counts$n[is.na(flag_counts$n)] <- 0L
  flag_counts$slide <- factor(flag_counts$slide, levels = levels(aoi$slide))
  flag_counts$segment <- factor(flag_counts$segment, levels = levels(aoi$segment))
  flag_counts$flag <- factor(flag_counts$flag, levels = rev(unname(flag_labels)))
  flag_counts$n <- as.integer(flag_counts$n)
  rownames(flag_counts) <- NULL
  stopifnot(all(is.finite(flag_counts$n)), all(flag_counts$n >= 0))

  list(
    aoi = aoi,
    metrics = metrics,
    flag_counts = flag_counts,
    thresholds = thresholds,
    provenance = list(
      lower_q = lower_q,
      upper_q = upper_q,
      flag_rule = "descriptive QC flags only: varying bottom-tail library/detected/area, varying top-tail background/q3, and nuclei < 0 sentinel; no AOIs excluded",
      n_aoi = nrow(aoi),
      n_flagged_any = sum(aoi$flag_any)
    )
  )
}

.geomx_sample_quantiles <- function(mat, meta, label,
                                    probs = c(0.05, 0.25, 0.50, 0.75, 0.95)) {
  stopifnot(is.matrix(mat), is.data.frame(meta), identical(colnames(mat), rownames(meta)),
            length(label) == 1L, nzchar(label), is.numeric(probs),
            all(is.finite(probs)), all(probs > 0), all(probs < 1))
  qnames <- paste0("q", sprintf("%02d", round(probs * 100)))
  qq <- t(vapply(seq_len(ncol(mat)), function(j) {
    x <- as.numeric(mat[, j])
    stopifnot(all(is.finite(x)))
    as.numeric(stats::quantile(x, probs = probs, names = FALSE, type = 7))
  }, numeric(length(probs))))
  colnames(qq) <- qnames
  out <- data.frame(
    aoi = colnames(mat),
    slide = factor(as.character(meta$slide), levels = sort(unique(as.character(meta$slide)))),
    roi = as.character(meta$roi),
    segment = factor(as.character(meta$segment), levels = sort(unique(as.character(meta$segment)))),
    genotype = factor(as.character(meta$genotype), levels = genotype_levels),
    method = label,
    stringsAsFactors = FALSE
  )
  out <- cbind(out, as.data.frame(qq, check.names = FALSE))
  rownames(out) <- NULL
  stopifnot(!anyNA(out$slide), !anyNA(out$segment), !anyNA(out$genotype),
            all(vapply(qnames, function(nm) all(is.finite(out[[nm]])), logical(1))))
  out
}

.geomx_voom_trend_descriptor <- function(v, max_points = 1600L) {
  stopifnot(is.list(v), is.numeric(max_points), length(max_points) == 1L, max_points >= 100L)
  if (is.null(v$voom.xy) || is.null(v$voom.line)) {
    stop("GeoMx voom object lacks saved trend data", call. = FALSE)
  }
  points <- data.frame(
    mean_log_count = as.numeric(v$voom.xy$x),
    sqrt_sd = as.numeric(v$voom.xy$y),
    stringsAsFactors = FALSE
  )
  points <- points[is.finite(points$mean_log_count) & is.finite(points$sqrt_sd), ,
                   drop = FALSE]
  if (nrow(points) > max_points) {
    points <- points[order(points$mean_log_count, points$sqrt_sd, method = "radix"), ,
                     drop = FALSE]
    keep <- unique(round(seq(1, nrow(points), length.out = max_points)))
    points <- points[keep, , drop = FALSE]
  }
  line <- data.frame(
    mean_log_count = as.numeric(v$voom.line$x),
    sqrt_sd = as.numeric(v$voom.line$y),
    stringsAsFactors = FALSE
  )
  line <- line[is.finite(line$mean_log_count) & is.finite(line$sqrt_sd), , drop = FALSE]
  line <- line[order(line$mean_log_count, method = "radix"), , drop = FALSE]
  rownames(points) <- rownames(line) <- NULL
  if (!nrow(points) || !nrow(line)) stop("GeoMx voom trend has no finite rows", call. = FALSE)
  list(
    points = points,
    line = line,
    labels = list(
      x = v$voom.xy$xlab %||% "log2 count size",
      y = v$voom.xy$ylab %||% "sqrt standard deviation"
    ),
    provenance = list(n_points = nrow(points), max_points = max_points)
  )
}

geomx_normalization_descriptor <- function(counts, meta, design, voom_trend,
                                           min_count = 5) {
  stopifnot(is.matrix(counts), is.data.frame(meta), is.matrix(design),
            identical(colnames(counts), rownames(meta)),
            identical(colnames(counts), rownames(design)),
            qr(design)$rank == ncol(design),
            is.list(voom_trend),
            is.numeric(min_count), length(min_count) == 1L, min_count >= 0)
  need <- c("slide", "roi", "segment", "genotype", "q3_factor", "neg_background")
  missing <- setdiff(need, names(meta))
  if (length(missing)) {
    stop("GeoMx normalization descriptor metadata missing columns: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  dge0 <- edgeR::DGEList(counts = counts)
  keep <- edgeR::filterByExpr(dge0, design = design, min.count = min_count)
  if (!any(keep)) stop("no GeoMx features passed filterByExpr for normalization descriptor",
                       call. = FALSE)
  dge_raw <- dge0[keep, , keep.lib.sizes = FALSE]
  raw_logcpm <- edgeR::cpm(dge_raw, normalized.lib.sizes = FALSE, log = TRUE,
                           prior.count = 1)
  dge_norm <- edgeR::normLibSizes(dge_raw, method = "TMM")
  norm_logcpm <- edgeR::cpm(dge_norm, normalized.lib.sizes = TRUE, log = TRUE,
                            prior.count = 1)
  stopifnot(identical(colnames(raw_logcpm), rownames(meta)),
            identical(colnames(norm_logcpm), rownames(meta)))

  distribution <- rbind(
    .geomx_sample_quantiles(raw_logcpm, meta, "Raw logCPM"),
    .geomx_sample_quantiles(norm_logcpm, meta, "TMM logCPM")
  )
  distribution$method <- factor(distribution$method,
                                levels = c("Raw logCPM", "TMM logCPM"))

  gene_median <- apply(norm_logcpm, 1L, stats::median)
  rle_mat <- sweep(norm_logcpm, 1L, gene_median, "-")
  rle <- .geomx_sample_quantiles(rle_mat, meta, "TMM logCPM RLE",
                                 probs = c(0.10, 0.25, 0.50, 0.75, 0.90))

  library_size <- Matrix::colSums(counts)
  norm_factor <- dge_norm$samples$norm.factors[match(rownames(meta), rownames(dge_norm$samples))]
  background <- data.frame(
    aoi = rownames(meta),
    slide = factor(as.character(meta$slide), levels = sort(unique(as.character(meta$slide)))),
    roi = as.character(meta$roi),
    segment = factor(as.character(meta$segment), levels = sort(unique(as.character(meta$segment)))),
    genotype = factor(as.character(meta$genotype), levels = genotype_levels),
    q3_factor = as.numeric(meta$q3_factor),
    neg_background = as.numeric(meta$neg_background),
    library_size = as.numeric(library_size[rownames(meta)]),
    tmm_norm_factor = as.numeric(norm_factor),
    stringsAsFactors = FALSE
  )
  stopifnot(!anyNA(background$slide), !anyNA(background$segment),
            !anyNA(background$genotype), all(is.finite(background$q3_factor)),
            all(background$q3_factor > 0), all(is.finite(background$neg_background)),
            all(background$neg_background > 0), all(is.finite(background$library_size)),
            all(background$library_size >= 0), all(is.finite(background$tmm_norm_factor)),
            all(background$tmm_norm_factor > 0))
  q3_bg_rho <- suppressWarnings(stats::cor(log10(background$q3_factor),
                                           log10(background$neg_background),
                                           method = "spearman"))
  if (!is.finite(q3_bg_rho)) q3_bg_rho <- NA_real_

  list(
    distribution = distribution,
    rle = rle,
    background = background,
    voom = voom_trend,
    provenance = list(
      n_aoi = ncol(counts),
      n_input_features = nrow(counts),
      n_kept_features = sum(keep),
      min_count = min_count,
      distribution = "raw logCPM versus TMM-normalized logCPM on the same filterByExpr-kept GeoMx genes",
      rle = "relative log expression computed from TMM-normalized logCPM after subtracting each gene median across AOIs",
      q3_neg_background_spearman = q3_bg_rho,
      model = "primary GeoMx DE still uses limma-voom with slide fixed effect and bio-unit duplicateCorrelation"
    )
  )
}

geomx_ordination_descriptor <- function(counts, meta, design, min_count = 5,
                                        n_variable_features = 2000L,
                                        n_pc = 6L,
                                        n_loadings = 12L) {
  stopifnot(is.matrix(counts), is.data.frame(meta), is.matrix(design),
            identical(colnames(counts), rownames(meta)),
            identical(colnames(counts), rownames(design)),
            qr(design)$rank == ncol(design),
            is.numeric(min_count), length(min_count) == 1L, min_count >= 0,
            is.numeric(n_variable_features), length(n_variable_features) == 1L,
            n_variable_features >= 2L,
            is.numeric(n_pc), length(n_pc) == 1L, n_pc >= 2L,
            is.numeric(n_loadings), length(n_loadings) == 1L, n_loadings >= 1L,
            anyDuplicated(rownames(counts)) == 0L)
  need <- c("slide", "roi", "segment", "genotype", "SampleID", "bio_unit")
  missing <- setdiff(need, names(meta))
  if (length(missing)) {
    stop("GeoMx ordination descriptor metadata missing columns: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }

  dge0 <- edgeR::DGEList(counts = counts)
  keep <- edgeR::filterByExpr(dge0, design = design, min.count = min_count)
  if (!any(keep)) stop("no GeoMx features passed filterByExpr for ordination descriptor",
                       call. = FALSE)
  dge <- edgeR::normLibSizes(dge0[keep, , keep.lib.sizes = FALSE], method = "TMM")
  logcpm <- edgeR::cpm(dge, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 1)
  stopifnot(identical(colnames(logcpm), rownames(meta)), all(is.finite(logcpm)))

  feature_sd <- apply(logcpm, 1L, stats::sd)
  feature_var <- feature_sd^2
  var_rank <- data.frame(
    symbol = rownames(logcpm),
    variance = as.numeric(feature_var),
    stringsAsFactors = FALSE
  )
  var_rank <- var_rank[is.finite(var_rank$variance) & var_rank$variance > 0, ,
                       drop = FALSE]
  if (nrow(var_rank) < 2L) {
    stop("GeoMx ordination has too few variable filter-passing genes", call. = FALSE)
  }
  var_rank <- var_rank[order(-var_rank$variance, var_rank$symbol, method = "radix"), ,
                       drop = FALSE]
  variable_genes <- utils::head(var_rank$symbol, as.integer(n_variable_features))

  mat <- logcpm[variable_genes, , drop = FALSE]
  row_mean <- rowMeans(mat)
  row_sd <- apply(mat, 1L, stats::sd)
  row_sd[!is.finite(row_sd) | row_sd <= 0] <- 1
  z <- sweep(sweep(mat, 1L, row_mean, "-"), 1L, row_sd, "/")
  z[!is.finite(z)] <- 0

  pc <- stats::prcomp(t(z), center = FALSE, scale. = FALSE)
  if (ncol(pc$x) < 2L) stop("GeoMx PCA returned fewer than two components", call. = FALSE)
  for (j in seq_len(ncol(pc$rotation))) {
    k <- which.max(abs(pc$rotation[, j]))
    if (is.finite(pc$rotation[k, j]) && pc$rotation[k, j] < 0) {
      pc$rotation[, j] <- -pc$rotation[, j]
      pc$x[, j] <- -pc$x[, j]
    }
  }
  variance_fraction <- pc$sdev^2 / sum(pc$sdev^2)
  pc_keep <- seq_len(min(as.integer(n_pc), length(variance_fraction)))
  scree <- data.frame(
    pc = paste0("PC", pc_keep),
    pc_num = pc_keep,
    variance_fraction = as.numeric(variance_fraction[pc_keep]),
    variance_percent = 100 * as.numeric(variance_fraction[pc_keep]),
    stringsAsFactors = FALSE
  )

  loading_keep <- seq_len(min(2L, ncol(pc$rotation)))
  loadings <- do.call(rbind, lapply(loading_keep, function(j) {
    ld <- data.frame(
      symbol = rownames(pc$rotation),
      pc = paste0("PC", j),
      pc_num = j,
      loading = as.numeric(pc$rotation[, j]),
      variance_fraction = as.numeric(variance_fraction[j]),
      stringsAsFactors = FALSE
    )
    ld$abs_loading <- abs(ld$loading)
    ld <- ld[is.finite(ld$loading), , drop = FALSE]
    ld <- ld[order(-ld$abs_loading, ld$symbol, method = "radix"), , drop = FALSE]
    ld <- utils::head(ld, as.integer(n_loadings))
    ld$loading_rank <- seq_len(nrow(ld))
    ld
  }))
  rownames(loadings) <- NULL
  loadings$direction <- factor(ifelse(loadings$loading >= 0, "positive", "negative"),
                               levels = c("negative", "positive"))

  mds <- stats::cmdscale(stats::dist(t(z)), k = 2L, eig = TRUE)
  mds_points <- as.matrix(mds$points)
  if (ncol(mds_points) < 2L) stop("GeoMx MDS returned fewer than two dimensions",
                                  call. = FALSE)
  for (j in seq_len(2L)) {
    k <- which.max(abs(mds_points[, j]))
    if (is.finite(mds_points[k, j]) && mds_points[k, j] < 0) {
      mds_points[, j] <- -mds_points[, j]
    }
  }
  positive_eig <- mds$eig[is.finite(mds$eig) & mds$eig > 0]
  mds_fraction <- rep(NA_real_, 2L)
  if (length(positive_eig) >= 2L && sum(positive_eig) > 0) {
    mds_fraction <- positive_eig[seq_len(2L)] / sum(positive_eig)
  }

  sample <- data.frame(
    aoi = rownames(meta),
    slide = factor(as.character(meta$slide), levels = sort(unique(as.character(meta$slide)))),
    roi = as.character(meta$roi),
    segment = factor(as.character(meta$segment), levels = sort(unique(as.character(meta$segment)))),
    sample_id = as.character(meta$SampleID),
    genotype = factor(as.character(meta$genotype), levels = genotype_levels),
    bio_unit = factor(as.character(meta$bio_unit)),
    pc1 = as.numeric(pc$x[rownames(meta), 1L]),
    pc2 = as.numeric(pc$x[rownames(meta), 2L]),
    pc1_var = as.numeric(variance_fraction[1L]),
    pc2_var = as.numeric(variance_fraction[2L]),
    mds1 = as.numeric(mds_points[rownames(meta), 1L]),
    mds2 = as.numeric(mds_points[rownames(meta), 2L]),
    mds1_var = as.numeric(mds_fraction[1L]),
    mds2_var = as.numeric(mds_fraction[2L]),
    stringsAsFactors = FALSE
  )
  rownames(sample) <- NULL
  stopifnot(!anyNA(sample$slide), !anyNA(sample$segment), !anyNA(sample$genotype),
            all(is.finite(sample$pc1)), all(is.finite(sample$pc2)),
            all(is.finite(sample$pc1_var)), all(is.finite(sample$pc2_var)),
            all(is.finite(sample$mds1)), all(is.finite(sample$mds2)),
            all(is.finite(scree$variance_fraction)), all(is.finite(loadings$loading)))

  list(
    sample = sample,
    scree = scree,
    loadings = loadings,
    provenance = list(
      n_aoi = ncol(counts),
      n_input_features = nrow(counts),
      n_kept_features = sum(keep),
      n_variable_features = length(variable_genes),
      min_count = min_count,
      transform = "TMM-normalized logCPM, row-centered and row-scaled over top variable filterByExpr-kept genes",
      pca = "stats::prcomp on AOIs x selected genes; PC signs oriented so the largest absolute loading is positive",
      mds = "classical MDS on Euclidean AOI distances over the same scaled expression matrix",
      loading_rule = sprintf("top %s absolute loadings for PC1 and PC2", as.integer(n_loadings))
    )
  )
}

geomx_sample_heatmap_descriptor <- function(counts, meta, design, top,
                                            min_count = 5,
                                            n_variable_features = 40L,
                                            n_track_features = 5L,
                                            z_limit = 2.5) {
  stopifnot(is.matrix(counts), is.data.frame(meta), is.matrix(design),
            is.list(top),
            identical(colnames(counts), rownames(meta)),
            identical(colnames(counts), rownames(design)),
            qr(design)$rank == ncol(design),
            is.numeric(min_count), length(min_count) == 1L, min_count >= 0,
            is.numeric(n_variable_features), length(n_variable_features) == 1L,
            n_variable_features >= 2L,
            is.numeric(n_track_features), length(n_track_features) == 1L,
            n_track_features >= 2L,
            is.numeric(z_limit), length(z_limit) == 1L, is.finite(z_limit),
            z_limit > 0,
            anyDuplicated(rownames(counts)) == 0L)
  need <- c("slide", "roi", "genotype", "SampleID", "bio_unit")
  missing <- setdiff(need, names(meta))
  if (length(missing)) {
    stop("GeoMx sample heatmap metadata missing columns: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }

  dge0 <- edgeR::DGEList(counts = counts)
  keep <- edgeR::filterByExpr(dge0, design = design, min.count = min_count)
  if (!any(keep)) stop("no GeoMx features passed filterByExpr for sample heatmap",
                       call. = FALSE)
  dge <- edgeR::normLibSizes(dge0[keep, , keep.lib.sizes = FALSE], method = "TMM")
  logcpm <- edgeR::cpm(dge, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 1)
  stopifnot(identical(colnames(logcpm), rownames(meta)), all(is.finite(logcpm)))

  feature_var <- apply(logcpm, 1L, stats::var)
  var_rank <- data.frame(
    symbol = rownames(logcpm),
    variance = as.numeric(feature_var),
    stringsAsFactors = FALSE
  )
  var_rank <- var_rank[is.finite(var_rank$variance) & var_rank$variance > 0, ,
                       drop = FALSE]
  if (nrow(var_rank) < 2L) {
    stop("GeoMx sample heatmap has too few variable filter-passing genes",
         call. = FALSE)
  }
  var_rank <- var_rank[order(-var_rank$variance, var_rank$symbol, method = "radix"), ,
                       drop = FALSE]
  feature_ids <- utils::head(var_rank$symbol, as.integer(n_variable_features))
  dam_track_genes <- geomx_dam_track_genes(logcpm)

  mat <- logcpm[feature_ids, , drop = FALSE]
  row_mean <- rowMeans(mat)
  row_sd <- apply(mat, 1L, stats::sd)
  row_sd[!is.finite(row_sd) | row_sd <= 0] <- 1
  z <- sweep(sweep(mat, 1L, row_mean, "-"), 1L, row_sd, "/")
  z[!is.finite(z)] <- 0

  raw_meta <- attr(meta, "geomx_raw_meta")
  if (!is.data.frame(raw_meta)) raw_meta <- meta
  if (!identical(rownames(raw_meta), rownames(meta))) {
    stop("GeoMx raw metadata rows do not align with model metadata", call. = FALSE)
  }
  if (length(dam_track_genes) < 2L) {
    stop("GeoMx DAM-gene track clustering needs >=2 displayed DAM genes", call. = FALSE)
  }
  track_mat <- logcpm[dam_track_genes, , drop = FALSE]
  track_mean <- rowMeans(track_mat)
  track_sd <- apply(track_mat, 1L, stats::sd)
  track_sd[!is.finite(track_sd) | track_sd <= 0] <- 1
  track_z <- sweep(sweep(track_mat, 1L, track_mean, "-"), 1L, track_sd, "/")
  track_z[!is.finite(track_z)] <- 0
  feature_hclust <- stats::hclust(stats::dist(z), method = "average")
  track_gene_hclust <- stats::hclust(stats::dist(track_z), method = "average")
  dam_track_genes <- utils::head(rev(rownames(track_z)[track_gene_hclust$order]),
                                 as.integer(n_track_features))
  track_z <- track_z[dam_track_genes, , drop = FALSE]
  sample_hclust <- stats::hclust(stats::dist(t(track_z)), method = "average")
  sample_dend <- stats::reorder(stats::as.dendrogram(sample_hclust),
                                wts = colMeans(track_z),
                                agglo.FUN = mean)
  sample_order <- labels(sample_dend)
  feature_order <- rownames(z)[feature_hclust$order]
  z_plot <- pmax(pmin(z[feature_order, sample_order, drop = FALSE], z_limit), -z_limit)

  heatmap <- as.data.frame(as.table(z_plot), stringsAsFactors = FALSE)
  names(heatmap) <- c("symbol", "aoi", "z")
  heatmap$sample_rank <- match(heatmap$aoi, sample_order)
  heatmap$feature_rank <- match(heatmap$symbol, feature_order)
  heatmap$symbol_plot <- factor(heatmap$symbol, levels = rev(feature_order))
  rownames(heatmap) <- NULL

  spatial <- geomx_spatial_descriptor(counts, meta, top)
  sample <- meta[sample_order, , drop = FALSE]
  sample_frame <- data.frame(
    aoi = rownames(sample),
    sample_rank = seq_along(sample_order),
    slide = factor(as.character(sample$slide), levels = sort(unique(as.character(meta$slide)))),
    roi = factor(as.character(sample$roi)),
    sample_id = as.character(sample$SampleID),
    genotype = factor(as.character(sample$genotype), levels = genotype_levels),
    tau_background = factor(
      ifelse(as.character(sample$genotype) %in% c("P301S", "NLGF_P301S"),
             "P301S", "MAPTKI"),
      levels = c("MAPTKI", "P301S")
    ),
    amyloid_status = factor(
      ifelse(startsWith(as.character(sample$genotype), "NLGF"), "NLGF+", "NLGF-"),
      levels = c("NLGF-", "NLGF+")
    ),
    bio_unit = factor(as.character(sample$bio_unit)),
    signed_response_score = as.numeric(spatial$aoi$signed_response_score[
      match(sample_order, spatial$aoi$aoi)
    ]),
    stringsAsFactors = FALSE
  )
  feature_frame <- var_rank[match(feature_order, var_rank$symbol), , drop = FALSE]
  feature_frame$feature_rank <- seq_along(feature_order)
  rownames(feature_frame) <- NULL

  raw_sample <- raw_meta[sample_order, , drop = FALSE]
  biology_tracks <- geomx_sample_biology_tracks(logcpm, sample_order, sample_frame,
                                                marker_genes = dam_track_genes,
                                                track_group = "DAM genes",
                                                source = "DAM marker TMM logCPM row z-score")
  metadata_tracks <- geomx_sample_metadata_tracks(raw_sample, sample_frame,
                                                  biology_tracks = biology_tracks)

  stopifnot(all(is.finite(heatmap$z)), all(is.finite(heatmap$sample_rank)),
            all(is.finite(heatmap$feature_rank)), !anyNA(heatmap$symbol_plot),
            !anyNA(sample_frame$slide), !anyNA(sample_frame$genotype),
            !anyNA(sample_frame$tau_background), !anyNA(sample_frame$amyloid_status),
            !anyNA(sample_frame$bio_unit),
            all(is.finite(sample_frame$signed_response_score)),
            all(is.finite(feature_frame$variance)),
            all(c("sample_rank", "track_id", "track_label", "track_group",
                  "track_type", "value") %in% names(metadata_tracks)))

  list(
    heatmap = heatmap,
    sample = sample_frame,
    features = feature_frame,
    metadata_tracks = metadata_tracks,
    provenance = list(
      n_aoi = ncol(counts),
      n_input_features = nrow(counts),
      n_kept_features = sum(keep),
      n_variable_features = length(feature_ids),
      n_track_features = length(dam_track_genes),
      min_count = min_count,
      z_limit = z_limit,
      transform = "TMM-normalized logCPM, top variable filterByExpr-kept genes, row z-scored and clipped for display",
      clustering = "AOIs average-linkage clustered on the displayed five DAM-gene z-score tracks, then dendrogram-rotated by mean displayed DAM z-score; rows preserve the first five genes from the prior full DAM-gene row cluster",
      sample_tracks = unique(metadata_tracks$track_id),
      n_metadata_columns = ncol(raw_sample),
      n_visible_tracks = length(unique(metadata_tracks$track_id)),
      display = "focused AOI design and first-five DAM-gene tracks with clustered AOI columns; no AOIs or genes are excluded beyond the existing primary-model filter"
    )
  )
}

geomx_dam_track_genes <- function(logcpm,
                                  omit = c("Fth1", "Lyz2")) {
  stopifnot(is.matrix(logcpm), is.character(omit))
  genes <- setdiff(canonical_microglia_markers$DAM, omit)
  genes[genes %in% rownames(logcpm)]
}

geomx_sample_biology_tracks <- function(logcpm, sample_order, sample_frame,
                                        marker_genes,
                                        track_group,
                                        source,
                                        z_limit = 2.5) {
  stopifnot(is.matrix(logcpm), is.character(sample_order),
            identical(sample_order, sample_frame$aoi),
            all(sample_order %in% colnames(logcpm)),
            is.character(marker_genes), length(marker_genes) >= 1L,
            is.character(track_group), length(track_group) == 1L,
            !is.na(track_group), track_group != "",
            is.character(source), length(source) == 1L,
            !is.na(source), source != "",
            is.numeric(z_limit), length(z_limit) == 1L, is.finite(z_limit),
            z_limit > 0)

  z_rows <- function(symbols) {
    symbols <- unique(trimws(as.character(symbols)))
    symbols <- symbols[!is.na(symbols) & symbols != ""]
    symbols <- symbols[symbols %in% rownames(logcpm)]
    if (!length(symbols)) return(NULL)
    mat <- logcpm[symbols, sample_order, drop = FALSE]
    row_mean <- rowMeans(mat)
    row_sd <- apply(mat, 1L, stats::sd)
    row_sd[!is.finite(row_sd) | row_sd <= 0] <- 1
    z <- sweep(sweep(mat, 1L, row_mean, "-"), 1L, row_sd, "/")
    z[!is.finite(z)] <- 0
    list(symbols = symbols, z = z)
  }
  make_numeric_track <- function(track_id, track_label, track_group, values,
                                 source, n_features = NA_integer_) {
    values <- pmax(pmin(as.numeric(values), z_limit), -z_limit)
    data.frame(
      sample_rank = sample_frame$sample_rank,
      track_id = track_id,
      track_label = track_label,
      track_group = track_group,
      track_type = "numeric",
      source = source,
      value = format(signif(values, 5), trim = TRUE, scientific = FALSE),
      numeric_value = values,
      n_features = as.integer(n_features),
      stringsAsFactors = FALSE
    )
  }

  tracks <- list()
  for (gene in marker_genes) {
    zr <- z_rows(gene)
    if (is.null(zr)) next
    tracks[[length(tracks) + 1L]] <- make_numeric_track(
      paste0("gene_", gene),
      gene,
      track_group,
      as.numeric(zr$z[gene, ]),
      source = source,
      n_features = 1L
    )
  }
  out <- if (length(tracks)) do.call(rbind, tracks) else data.frame(
    sample_rank = integer(), track_id = character(), track_label = character(),
    track_group = character(), track_type = character(), source = character(),
    value = character(), numeric_value = numeric(), n_features = integer(),
    stringsAsFactors = FALSE
  )
  rownames(out) <- NULL
  out
}

geomx_sample_metadata_tracks <- function(raw_sample, sample_frame,
                                         biology_tracks = NULL) {
  stopifnot(is.data.frame(raw_sample), is.data.frame(sample_frame),
            identical(rownames(raw_sample), sample_frame$aoi),
            all(c("sample_rank", "genotype", "tau_background", "amyloid_status",
                  "bio_unit", "signed_response_score") %in% names(sample_frame)))

  label_map <- c(
    "ROI Coordinate X" = "ROI X",
    "ROI Coordinate Y" = "ROI Y",
    "NegGeoMean_Mm_R_NGS_WTA_v1.0" = "NegGeoMean",
    "q_norm_qFactors" = "Q3 factor",
    "nCount_RNA" = "library size"
  )
  categorical_numeric <- c("bio_rep", "tech_rep", "slide_rep")
  group_for <- function(id) {
    if (id %in% c("genotype", "tau_background", "amyloid_status")) {
      "design"
    } else {
      "other"
    }
  }
  label_for <- function(id) {
    if (id %in% names(label_map)) unname(label_map[[id]]) else id
  }
  value_string <- function(x) {
    if (is.numeric(x) || is.integer(x)) {
      out <- format(signif(as.numeric(x), 5), trim = TRUE, scientific = FALSE)
    } else {
      out <- as.character(x)
    }
    out[is.na(out)] <- "NA"
    out[out == ""] <- "<blank>"
    out
  }
  make_track <- function(track_id, values, track_label = label_for(track_id),
                         track_group = group_for(track_id), source = "metadata",
                         track_type = NULL) {
    if (is.null(track_type)) {
      track_type <- if ((is.numeric(values) || is.integer(values)) &&
                        !track_id %in% categorical_numeric) {
        "numeric"
      } else {
        "categorical"
      }
    }
    numeric_value <- if (identical(track_type, "numeric")) {
      suppressWarnings(as.numeric(values))
    } else {
      rep(NA_real_, length(values))
    }
    data.frame(
      sample_rank = sample_frame$sample_rank,
      track_id = track_id,
      track_label = track_label,
      track_group = track_group,
      track_type = track_type,
      source = source,
      value = value_string(values),
      numeric_value = numeric_value,
      n_features = NA_integer_,
      stringsAsFactors = FALSE
    )
  }

  tracks <- list()
  added_raw <- character()
  append_track <- function(x) {
    tracks[[length(tracks) + 1L]] <<- x
    invisible(NULL)
  }
  append_raw <- function(col) {
    if (col %in% names(raw_sample) && !col %in% added_raw) {
      append_track(make_track(col, raw_sample[[col]]))
      added_raw <<- c(added_raw, col)
    }
    invisible(NULL)
  }

  append_raw("genotype")
  append_track(make_track("tau_background", sample_frame$tau_background,
                          track_label = "tau", source = "derived",
                          track_type = "categorical"))
  append_track(make_track("amyloid_status", sample_frame$amyloid_status,
                          track_label = "amyloid", source = "derived",
                          track_type = "categorical"))
  if (is.data.frame(biology_tracks) && nrow(biology_tracks)) {
    append_track(biology_tracks)
  }
  out <- do.call(rbind, tracks)
  out$track_group <- factor(out$track_group,
                            levels = c("design", "DAM genes"))
  rownames(out) <- NULL
  stopifnot(!anyNA(out$sample_rank), !anyNA(out$track_id),
            !anyNA(out$track_label), !anyNA(out$track_group),
            !anyNA(out$track_type), !anyNA(out$value))
  out
}

geomx_spatial_program_descriptor <- function(counts, meta, design,
                                             min_count = 5,
                                             signature_sets = list(
                                               Homeostatic = canonical_microglia_markers$Homeostatic,
                                               DAM = canonical_microglia_markers$DAM,
                                               IFN = canonical_microglia_markers$IFN,
                                               MHC_APC = canonical_microglia_markers$MHC_APC
                                             ),
                                             gene_features = c("Apoe", "Trem2"),
                                             z_limit = 2.5) {
  stopifnot(is.matrix(counts), is.data.frame(meta), is.matrix(design),
            identical(colnames(counts), rownames(meta)),
            identical(colnames(counts), rownames(design)),
            qr(design)$rank == ncol(design),
            is.numeric(min_count), length(min_count) == 1L, min_count >= 0,
            is.list(signature_sets), length(signature_sets) >= 1L,
            !is.null(names(signature_sets)), !any(names(signature_sets) == ""),
            is.character(gene_features), length(gene_features) >= 1L,
            is.numeric(z_limit), length(z_limit) == 1L, is.finite(z_limit),
            z_limit > 0,
            anyDuplicated(rownames(counts)) == 0L)
  need <- c("slide", "roi", "segment", "genotype", "SampleID", "bio_unit",
            "area", "x", "y")
  missing <- setdiff(need, names(meta))
  if (length(missing)) {
    stop("GeoMx spatial-program metadata missing columns: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }

  dge0 <- edgeR::DGEList(counts = counts)
  keep <- edgeR::filterByExpr(dge0, design = design, min.count = min_count)
  if (!any(keep)) stop("no GeoMx features passed filterByExpr for spatial programs",
                       call. = FALSE)
  dge <- edgeR::normLibSizes(dge0[keep, , keep.lib.sizes = FALSE], method = "TMM")
  logcpm <- edgeR::cpm(dge, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 1)
  stopifnot(identical(colnames(logcpm), rownames(meta)), all(is.finite(logcpm)))

  score_rows <- function(symbols) {
    symbols <- unique(trimws(as.character(symbols)))
    symbols <- symbols[!is.na(symbols) & symbols != ""]
    symbols[symbols %in% rownames(logcpm)]
  }
  z_score <- function(symbols) {
    symbols <- score_rows(symbols)
    if (!length(symbols)) {
      stop("GeoMx spatial-program score has no filter-passing genes", call. = FALSE)
    }
    mat <- logcpm[symbols, , drop = FALSE]
    row_mean <- rowMeans(mat)
    row_sd <- apply(mat, 1L, stats::sd)
    row_sd[!is.finite(row_sd) | row_sd <= 0] <- 1
    z <- sweep(sweep(mat, 1L, row_mean, "-"), 1L, row_sd, "/")
    z[!is.finite(z)] <- 0
    colMeans(z)
  }

  sig_levels <- names(signature_sets)
  program_catalog <- do.call(rbind, c(
    lapply(sig_levels, function(nm) {
      used <- score_rows(signature_sets[[nm]])
      data.frame(
        program_id = nm,
        program_label = switch(nm,
                               Homeostatic = "Homeostatic signature",
                               DAM = "DAM signature",
                               IFN = "IFN signature",
                               MHC_APC = "MHC/APC signature",
                               nm),
        program_type = "signature",
        n_input_features = length(unique(signature_sets[[nm]])),
        n_scored_features = length(used),
        scored_features = paste(used, collapse = ", "),
        stringsAsFactors = FALSE
      )
    }),
    lapply(gene_features, function(gene) {
      used <- score_rows(gene)
      data.frame(
        program_id = gene,
        program_label = gene,
        program_type = "single gene",
        n_input_features = 1L,
        n_scored_features = length(used),
        scored_features = paste(used, collapse = ", "),
        stringsAsFactors = FALSE
      )
    })
  ))
  if (any(program_catalog$n_scored_features < 1L)) {
    miss <- program_catalog$program_label[program_catalog$n_scored_features < 1L]
    stop("GeoMx spatial-program score(s) lack filter-passing features: ",
         paste(miss, collapse = ", "), call. = FALSE)
  }
  program_catalog$program_label <- factor(program_catalog$program_label,
                                          levels = program_catalog$program_label)
  rownames(program_catalog) <- NULL

  score_list <- c(
    lapply(sig_levels, function(nm) z_score(signature_sets[[nm]])),
    lapply(gene_features, function(gene) z_score(gene))
  )
  names(score_list) <- program_catalog$program_id
  base <- data.frame(
    aoi = rownames(meta),
    slide = factor(as.character(meta$slide), levels = sort(unique(as.character(meta$slide)))),
    roi = as.character(meta$roi),
    segment = factor(as.character(meta$segment), levels = sort(unique(as.character(meta$segment)))),
    sample_id = as.character(meta$SampleID),
    genotype = factor(as.character(meta$genotype), levels = genotype_levels),
    bio_unit = factor(as.character(meta$bio_unit)),
    aoi_area = as.numeric(meta$area),
    x_coord = as.numeric(meta$x),
    y_coord = as.numeric(meta$y),
    stringsAsFactors = FALSE
  )
  stopifnot(!anyNA(base$slide), !anyNA(base$segment), !anyNA(base$genotype),
            !anyNA(base$bio_unit), all(is.finite(base$aoi_area)),
            all(base$aoi_area > 0), all(is.finite(base$x_coord)),
            all(is.finite(base$y_coord)))

  aoi <- do.call(rbind, lapply(seq_len(nrow(program_catalog)), function(i) {
    pr <- program_catalog[i, , drop = FALSE]
    score <- as.numeric(score_list[[as.character(pr$program_id)]][base$aoi])
    data.frame(
      base,
      program_id = as.character(pr$program_id),
      program_label = pr$program_label,
      program_type = as.character(pr$program_type),
      score = score,
      score_abs = abs(score),
      stringsAsFactors = FALSE
    )
  }))
  aoi$program_label <- factor(as.character(aoi$program_label),
                              levels = levels(program_catalog$program_label))
  rownames(aoi) <- NULL
  stopifnot(!anyNA(aoi$program_label), all(is.finite(aoi$score)),
            all(is.finite(aoi$score_abs)))

  genotype_summary <- do.call(rbind, lapply(split(aoi, list(aoi$program_label, aoi$genotype),
                                                drop = TRUE), function(x) {
    if (!nrow(x)) return(NULL)
    data.frame(
      program_label = x$program_label[[1L]],
      genotype = x$genotype[[1L]],
      n_aoi = nrow(x),
      median_score = stats::median(x$score),
      q25_score = as.numeric(stats::quantile(x$score, 0.25, names = FALSE)),
      q75_score = as.numeric(stats::quantile(x$score, 0.75, names = FALSE)),
      mean_score = mean(x$score),
      stringsAsFactors = FALSE
    )
  }))
  genotype_summary$program_label <- factor(as.character(genotype_summary$program_label),
                                           levels = levels(program_catalog$program_label))
  genotype_summary$genotype <- factor(as.character(genotype_summary$genotype),
                                      levels = genotype_levels)
  rownames(genotype_summary) <- NULL
  stopifnot(!anyNA(genotype_summary$program_label), !anyNA(genotype_summary$genotype),
            all(is.finite(genotype_summary$n_aoi)),
            all(is.finite(genotype_summary$median_score)),
            all(is.finite(genotype_summary$q25_score)),
            all(is.finite(genotype_summary$q75_score)),
            all(is.finite(genotype_summary$mean_score)))

  list(
    aoi = aoi,
    programs = program_catalog,
    genotype_summary = genotype_summary,
    provenance = list(
      n_aoi = ncol(counts),
      n_input_features = nrow(counts),
      n_kept_features = sum(keep),
      min_count = min_count,
      z_limit = z_limit,
      transform = "TMM-normalized logCPM over filterByExpr-kept GeoMx genes; each feature row is z-scored across AOIs",
      score_rule = "signature scores are the mean row z-score over present filter-passing marker genes; single-gene panels show that gene's row z-score",
      coordinate_status = "coordinate-only GeoMx overlay: source tissue images/OME-TIFFs are not in the live repo",
      display = "descriptive spatial program overlay; no AOIs or genes are excluded beyond the existing primary-model filter"
    )
  )
}

geomx_gene_detection_descriptor <- function(counts, meta, design, min_count = 5,
                                            marker_sets = list(
                                              Microglia = microglia_identity_markers,
                                              Homeostatic = canonical_microglia_markers$Homeostatic,
                                              DAM = canonical_microglia_markers$DAM
                                            ),
                                            n_marker_label = 4L,
                                            n_marker_top = 8L,
                                            n_top_detected = 14L,
                                            n_bins = 24L) {
  stopifnot(is.matrix(counts), is.data.frame(meta), is.matrix(design),
            identical(colnames(counts), rownames(meta)),
            identical(colnames(counts), rownames(design)),
            qr(design)$rank == ncol(design),
            is.numeric(min_count), length(min_count) == 1L, min_count >= 0,
            is.list(marker_sets), length(marker_sets) >= 1L, !is.null(names(marker_sets)),
            !any(names(marker_sets) == ""),
            is.numeric(n_marker_label), length(n_marker_label) == 1L, n_marker_label >= 1L,
            is.numeric(n_marker_top), length(n_marker_top) == 1L, n_marker_top >= 1L,
            is.numeric(n_top_detected), length(n_top_detected) == 1L, n_top_detected >= 1L,
            is.numeric(n_bins), length(n_bins) == 1L, n_bins >= 4L,
            anyDuplicated(rownames(counts)) == 0L)

  dge0 <- edgeR::DGEList(counts = counts)
  keep <- edgeR::filterByExpr(dge0, design = design, min.count = min_count)
  if (!any(keep)) stop("no GeoMx features passed filterByExpr for gene detection descriptor",
                       call. = FALSE)
  dge <- edgeR::normLibSizes(dge0, method = "TMM")
  logcpm <- edgeR::cpm(dge, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 1)
  stopifnot(identical(rownames(logcpm), rownames(counts)),
            identical(colnames(logcpm), colnames(counts)),
            all(is.finite(logcpm)))

  genes <- data.frame(
    symbol = rownames(counts),
    mean_count = as.numeric(rowMeans(counts)),
    mean_logcpm = as.numeric(rowMeans(logcpm)),
    detect_fraction_any = as.numeric(rowMeans(counts > 0)),
    detect_fraction_min_count = as.numeric(rowMeans(counts >= min_count)),
    detected_aoi_any = as.integer(rowSums(counts > 0)),
    detected_aoi_min_count = as.integer(rowSums(counts >= min_count)),
    filter_pass = as.logical(keep),
    stringsAsFactors = FALSE
  )
  genes$filter_status <- factor(ifelse(genes$filter_pass, "filter passing", "low coverage"),
                                levels = c("low coverage", "filter passing"))
  stopifnot(all(is.finite(genes$mean_count)), all(genes$mean_count >= 0),
            all(is.finite(genes$mean_logcpm)),
            all(is.finite(genes$detect_fraction_any)),
            all(genes$detect_fraction_any >= 0), all(genes$detect_fraction_any <= 1),
            all(is.finite(genes$detect_fraction_min_count)),
            all(genes$detect_fraction_min_count >= 0),
            all(genes$detect_fraction_min_count <= 1),
            all(is.finite(genes$detected_aoi_any)),
            all(is.finite(genes$detected_aoi_min_count)))

  marker_sets <- lapply(marker_sets, function(x) unique(trimws(as.character(x))))
  marker_sets <- lapply(marker_sets, function(x) x[!is.na(x) & x != ""])
  if (any(vapply(marker_sets, length, integer(1)) == 0L)) {
    stop("GeoMx gene detection marker set is empty", call. = FALSE)
  }
  class_levels <- names(marker_sets)
  marker_catalog <- do.call(rbind, lapply(class_levels, function(cls) {
    data.frame(
      marker_class = cls,
      symbol = marker_sets[[cls]],
      n_signature = length(marker_sets[[cls]]),
      stringsAsFactors = FALSE
    )
  }))
  marker_catalog$marker_class <- factor(marker_catalog$marker_class, levels = class_levels)
  marker_catalog$present <- marker_catalog$symbol %in% genes$symbol
  marker_metrics <- merge(
    marker_catalog[marker_catalog$present, c("marker_class", "symbol", "n_signature"),
                   drop = FALSE],
    genes,
    by = "symbol",
    all.x = TRUE,
    sort = FALSE
  )
  marker_metrics$marker_class <- factor(as.character(marker_metrics$marker_class),
                                        levels = class_levels)
  if (!nrow(marker_metrics)) stop("no GeoMx marker genes present in count matrix",
                                  call. = FALSE)
  missing_classes <- class_levels[!class_levels %in% as.character(marker_metrics$marker_class)]
  if (length(missing_classes)) {
    stop("no GeoMx marker genes present for class(es): ",
         paste(missing_classes, collapse = ", "), call. = FALSE)
  }

  marker_order <- order(marker_metrics$marker_class,
                        -as.integer(marker_metrics$filter_pass),
                        -marker_metrics$detect_fraction_min_count,
                        -marker_metrics$mean_logcpm,
                        marker_metrics$symbol,
                        method = "radix")
  marker_metrics <- marker_metrics[marker_order, , drop = FALSE]
  marker_top <- do.call(rbind, lapply(class_levels, function(cls) {
    x <- marker_metrics[as.character(marker_metrics$marker_class) == cls, , drop = FALSE]
    utils::head(x, as.integer(n_marker_top))
  }))
  rownames(marker_top) <- NULL

  marker_labels <- do.call(rbind, lapply(class_levels, function(cls) {
    x <- marker_metrics[as.character(marker_metrics$marker_class) == cls, , drop = FALSE]
    utils::head(x, as.integer(n_marker_label))
  }))
  marker_labels <- marker_labels[order(-marker_labels$detect_fraction_min_count,
                                       -marker_labels$mean_logcpm,
                                       marker_labels$symbol,
                                       method = "radix"), , drop = FALSE]
  marker_labels <- marker_labels[!duplicated(marker_labels$symbol), , drop = FALSE]
  marker_labels$marker_classes <- vapply(marker_labels$symbol, function(s) {
    paste(as.character(unique(marker_metrics$marker_class[marker_metrics$symbol == s])),
          collapse = "/")
  }, character(1), USE.NAMES = FALSE)
  rownames(marker_labels) <- NULL

  breaks <- seq(0, 1, length.out = as.integer(n_bins) + 1L)
  bin_id <- cut(pmin(pmax(genes$detect_fraction_min_count, 0), 1),
                breaks = breaks, include.lowest = TRUE, labels = FALSE)
  bin_raw <- data.frame(
    filter_status = genes$filter_status,
    bin = as.integer(bin_id),
    n = 1L,
    stringsAsFactors = FALSE
  )
  bin_count <- stats::aggregate(n ~ filter_status + bin, data = bin_raw, FUN = sum)
  bin_grid <- expand.grid(
    filter_status = levels(genes$filter_status),
    bin = seq_len(as.integer(n_bins)),
    stringsAsFactors = FALSE
  )
  filter_bins <- merge(bin_grid, bin_count, by = c("filter_status", "bin"),
                       all.x = TRUE, sort = FALSE)
  filter_bins$n[is.na(filter_bins$n)] <- 0L
  filter_bins$filter_status <- factor(filter_bins$filter_status,
                                      levels = levels(genes$filter_status))
  filter_bins$bin_mid <- (breaks[filter_bins$bin] + breaks[filter_bins$bin + 1L]) / 2
  filter_bins$bin_min <- breaks[filter_bins$bin]
  filter_bins$bin_max <- breaks[filter_bins$bin + 1L]
  filter_bins$n <- as.integer(filter_bins$n)
  rownames(filter_bins) <- NULL

  top_detected <- genes[genes$filter_pass %in% TRUE, , drop = FALSE]
  top_detected <- top_detected[order(-top_detected$detect_fraction_min_count,
                                     -top_detected$mean_logcpm,
                                     -top_detected$mean_count,
                                     top_detected$symbol,
                                     method = "radix"), , drop = FALSE]
  top_detected <- utils::head(top_detected, as.integer(n_top_detected))
  top_detected$marker_class <- vapply(top_detected$symbol, function(s) {
    cls <- as.character(unique(marker_metrics$marker_class[marker_metrics$symbol == s]))
    if (length(cls)) paste(cls, collapse = "/") else "other"
  }, character(1), USE.NAMES = FALSE)
  rownames(top_detected) <- NULL

  marker_summary <- do.call(rbind, lapply(class_levels, function(cls) {
    cat <- marker_catalog[as.character(marker_catalog$marker_class) == cls, , drop = FALSE]
    met <- marker_metrics[as.character(marker_metrics$marker_class) == cls, , drop = FALSE]
    data.frame(
      marker_class = cls,
      n_signature = nrow(cat),
      n_present = nrow(met),
      n_filter_passing = sum(met$filter_pass),
      n_detected_half = sum(met$detect_fraction_min_count >= 0.5),
      median_detect_fraction_min_count = stats::median(met$detect_fraction_min_count),
      median_mean_logcpm = stats::median(met$mean_logcpm),
      stringsAsFactors = FALSE
    )
  }))
  marker_summary$marker_class <- factor(marker_summary$marker_class, levels = class_levels)
  rownames(marker_summary) <- NULL

  list(
    genes = genes,
    filter_bins = filter_bins,
    marker_top = marker_top,
    marker_labels = marker_labels,
    top_detected = top_detected,
    marker_summary = marker_summary,
    provenance = list(
      n_aoi = ncol(counts),
      n_input_features = nrow(counts),
      n_kept_features = sum(keep),
      n_low_coverage_features = sum(!keep),
      min_count = min_count,
      marker_sets = vapply(marker_sets, length, integer(1)),
      n_marker_present = stats::setNames(marker_summary$n_present,
                                         as.character(marker_summary$marker_class)),
      n_marker_filter_passing = stats::setNames(marker_summary$n_filter_passing,
                                                as.character(marker_summary$marker_class)),
      detect_fraction = "fraction of AOIs with raw count >= min_count; expression is TMM-normalized mean logCPM",
      filter_rule = "edgeR::filterByExpr on the primary GeoMx slide-adjusted design with min.count",
      display = "descriptive gene detectability / marker measurability diagnostic; no AOIs or genes are removed beyond the existing primary-model filter"
    )
  )
}

# voom + TMM + limma, optionally blocking repeated AOIs by bio_unit via duplicateCorrelation
# (the primary GeoMx model). Fails loud if the design has no residual df or the consensus
# correlation is non-finite.
.fit_geomx_voom <- function(counts, design, contrasts, block = NULL, min_count = 5) {
  stopifnot(identical(colnames(counts), rownames(design)),
            identical(rownames(contrasts), colnames(design)),
            qr(design)$rank == ncol(design))
  if (nrow(design) <= ncol(design)) {
    stop("GeoMx design has no residual degrees of freedom", call. = FALSE)
  }
  dge <- edgeR::DGEList(counts = counts)
  keep <- edgeR::filterByExpr(dge, design = design, min.count = min_count)
  if (!any(keep)) stop("no GeoMx features passed filterByExpr", call. = FALSE)
  dge <- edgeR::normLibSizes(dge[keep, , keep.lib.sizes = FALSE], method = "TMM")

  duplicate <- list(used = FALSE, consensus_correlation = NA_real_)
  if (!is.null(block)) {
    stopifnot(length(block) == ncol(counts), !anyNA(block))
    block <- factor(block)
    v0 <- limma::voom(dge, design = design, plot = FALSE)
    corfit <- limma::duplicateCorrelation(v0, design = design, block = block)
    if (!is.finite(corfit$consensus.correlation)) {
      stop("GeoMx duplicateCorrelation returned non-finite consensus correlation", call. = FALSE)
    }
    duplicate <- list(used = TRUE, consensus_correlation = unname(corfit$consensus.correlation))
    v <- limma::voom(dge, design = design, plot = FALSE, save.plot = TRUE, block = block,
                     correlation = duplicate$consensus_correlation)
    fit0 <- limma::lmFit(v, design = design, block = block,
                         correlation = duplicate$consensus_correlation)
  } else {
    v <- limma::voom(dge, design = design, plot = FALSE, save.plot = TRUE)
    fit0 <- limma::lmFit(v, design = design)
  }
  fit <- limma::eBayes(limma::contrasts.fit(fit0, contrasts), robust = TRUE)
  list(
    kept = sum(keep),
    n_input_features = nrow(counts),
    n_samples = ncol(counts),
    design_cols = colnames(design),
    duplicate_correlation = duplicate,
    voom_trend = .geomx_voom_trend_descriptor(v),
    top = .geomx_top_tables(fit, contrasts)
  )
}

run_geomx_de <- function(geomx, min_count = 5) {
  counts <- geomx_count_matrix(geomx)
  meta <- geomx_meta(geomx)
  stopifnot(identical(colnames(counts), rownames(meta)))
  fd <- geomx_slide_design(meta, include_slide = TRUE)
  primary <- .fit_geomx_voom(counts, fd$design, fd$contrasts, block = meta$bio_unit,
                             min_count = min_count)
  list(
    n_aoi = ncol(counts),
    n_bio_units = length(unique(as.character(meta$bio_unit))),
    primary = c(list(status = "fit", model = "voom_tmm_slide_duplicateCorrelation"), primary),
    sample_heatmap = geomx_sample_heatmap_descriptor(counts, meta, fd$design, primary$top,
                                                     min_count = min_count),
    provenance = list(
      counts = attr(counts, "geomx_count_provenance"),
      meta = attr(meta, "geomx_meta_provenance"),
      thresholds = list(min_count = min_count)
    )
  )
}

# --- 24M bulk sample-key match (shared by proteome + phospho) --------------------------------

# Normalise a Spectronaut export's intensity columns to run stubs, join to the 24M sample key,
# and assert exactly n_keep matched columns balanced 4/genotype. Returns the ordered column
# vector + a per-sample meta frame (sample_id/column/stub/label/genotype/run_index).
match_24m_bulk_columns <- function(tbl, sample_key, n_keep = 16L, modality = "bulk") {
  stopifnot(is.data.frame(tbl), is.data.frame(sample_key),
            all(c("genotype", "col_stub", "label") %in% names(sample_key)),
            length(modality) == 1L, nzchar(modality))
  if (nrow(sample_key) != n_keep) {
    stop("sample_key must contain exactly ", n_keep, " 24M rows", call. = FALSE)
  }
  if (anyDuplicated(sample_key$col_stub)) stop("duplicate sample-key stubs", call. = FALSE)
  hits <- match_intensity_columns(names(tbl), sample_key)
  hits <- hits[!is.na(hits$key_idx), , drop = FALSE]
  if (nrow(hits) != n_keep || anyDuplicated(hits$key_idx) ||
      !setequal(hits$key_idx, seq_len(n_keep))) {
    stop("expected exactly ", n_keep, "/", n_keep, " matched 24M ", modality,
         " intensity columns", call. = FALSE)
  }
  if (anyDuplicated(hits$stub)) {
    stop("duplicate matched ", modality, " intensity stubs", call. = FALSE)
  }
  hits <- hits[order(hits$key_idx), , drop = FALSE]
  geno <- factor(as.character(hits$genotype), levels = genotype_levels)
  if (anyNA(geno) || !all(as.integer(table(geno)) == n_keep / 4L)) {
    stop("24M ", modality, " columns are not balanced 4/genotype", call. = FALSE)
  }
  meta <- data.frame(
    sample_id = hits$stub,
    column = hits$column,
    stub = hits$stub,
    label = hits$label,
    genotype = geno,
    run_index = seq_len(nrow(hits)),
    stringsAsFactors = FALSE
  )
  rownames(meta) <- meta$sample_id
  list(columns = hits$column, meta = meta,
       matched = hits, n_expected = n_keep, n_matched = nrow(hits))
}

# --- 24M bulk proteome DE (protein-group-summed, log2 median-normalised, no batch) ----------

.split_gene_symbols <- function(x) {
  x <- trimws(as.character(x))
  x <- x[!is.na(x) & x != ""]
  if (!length(x)) return(character())
  out <- unlist(strsplit(x, "[;,]", perl = TRUE), use.names = FALSE)
  out <- trimws(out)
  out[!is.na(out) & out != ""]
}

protein_group_features <- function(proteomics_tbl) {
  stopifnot(is.data.frame(proteomics_tbl))
  need <- c("PG.ProteinGroups", "PG.Genes")
  stopifnot(all(need %in% names(proteomics_tbl)))
  protein_group <- trimws(as.character(proteomics_tbl[["PG.ProteinGroups"]]))
  gene_raw <- trimws(as.character(proteomics_tbl[["PG.Genes"]]))
  missing_group <- is.na(protein_group) | protein_group == ""
  keep <- !missing_group
  if (!any(keep)) stop("no non-empty PG.ProteinGroups in proteomics table", call. = FALSE)

  row_info <- data.frame(
    original_row = seq_len(nrow(proteomics_tbl)),
    protein_group = protein_group,
    gene_raw = gene_raw,
    stringsAsFactors = FALSE
  )[keep, , drop = FALSE]
  by_group <- split(row_info, row_info$protein_group)
  out <- do.call(rbind, lapply(sort(names(by_group), method = "radix"), function(pg) {
    rows <- by_group[[pg]]
    syms_all <- .split_gene_symbols(rows$gene_raw)
    syms <- unique(syms_all[order(syms_all, method = "radix")])
    data.frame(
      protein_group = pg,
      gene_first = if (length(syms_all)) syms_all[1] else NA_character_,
      gene_symbols = if (length(syms)) paste(syms, collapse = ";") else NA_character_,
      n_gene_symbols = length(syms),
      n_raw_rows = nrow(rows),
      first_original_row = min(rows$original_row),
      stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- out$protein_group
  attr(out, "protein_group_counts") <- list(
    n_rows = nrow(proteomics_tbl),
    n_missing_protein_group = sum(missing_group),
    n_rows_with_protein_group = sum(keep),
    n_protein_groups = nrow(out),
    n_duplicate_rows_by_group = sum(duplicated(protein_group[keep]))
  )
  out
}

aggregate_proteome_raw <- function(proteomics_tbl, columns, features) {
  stopifnot(is.data.frame(proteomics_tbl), is.character(columns), length(columns) >= 1L,
            all(columns %in% names(proteomics_tbl)), is.data.frame(features),
            "protein_group" %in% names(features))
  raw <- as.matrix(proteomics_tbl[, columns, drop = FALSE])
  storage.mode(raw) <- "double"
  protein_group <- trimws(as.character(proteomics_tbl[["PG.ProteinGroups"]]))
  missing_group <- is.na(protein_group) | protein_group == ""
  raw <- raw[!missing_group, , drop = FALSE]
  protein_group <- protein_group[!missing_group]

  nonpositive <- !is.na(raw) & raw <= 0
  raw[nonpositive] <- NA_real_
  by_group <- split(seq_len(nrow(raw)), protein_group)
  agg <- vapply(sort(names(by_group), method = "radix"), function(pg) {
    block <- raw[by_group[[pg]], , drop = FALSE]
    present <- !is.na(block)
    vals <- colSums(block, na.rm = TRUE)
    vals[colSums(present) == 0L] <- NA_real_
    vals
  }, numeric(ncol(raw)))
  agg <- t(agg)
  colnames(agg) <- columns
  stopifnot(identical(rownames(agg), features$protein_group))
  attr(agg, "proteome_aggregate_counts") <- list(
    n_raw_values = length(raw),
    n_missing_input = sum(is.na(as.matrix(proteomics_tbl[!missing_group, columns, drop = FALSE]))),
    n_nonpositive_to_na = sum(nonpositive),
    n_protein_groups = nrow(agg)
  )
  agg
}

prepare_proteome_24m_matrix <- function(proteomics_tbl, sample_key,
                                        min_present = 2L, min_groups = 4L) {
  match <- match_24m_bulk_columns(proteomics_tbl, sample_key, modality = "proteome")
  features <- protein_group_features(proteomics_tbl)
  feature_counts <- attr(features, "protein_group_counts")
  raw_sum <- aggregate_proteome_raw(proteomics_tbl, match$columns, features)
  colnames(raw_sum) <- rownames(match$meta)
  log_mat <- log2(raw_sum)
  attr(log_mat, "log2_counts") <- list(
    n_values = length(raw_sum),
    n_missing_input = sum(is.na(raw_sum)),
    n_missing_output = sum(is.na(log_mat))
  )
  norm <- median_normalise(log_mat)
  filt <- prevalence_filter(norm, match$meta$genotype,
                            min_present = min_present, min_groups = min_groups)
  features <- features[match(rownames(filt), features$protein_group), , drop = FALSE]
  stopifnot(identical(features$protein_group, rownames(filt)),
            identical(colnames(filt), rownames(match$meta)))
  list(
    matrix = filt,
    meta = match$meta,
    features = features,
    matched = match$matched,
    counts = c(feature_counts,
               attr(raw_sum, "proteome_aggregate_counts"),
               attr(log_mat, "log2_counts"),
               list(n_features_raw = nrow(raw_sum),
                    n_features_filtered = nrow(filt),
                    min_present = min_present,
                    min_groups = min_groups))
  )
}

annotate_proteome_top_tables <- function(top, features) {
  stopifnot(is.list(top), is.data.frame(features), "protein_group" %in% names(features))
  lapply(top, function(tbl) {
    stopifnot(is.data.frame(tbl), "feature" %in% names(tbl))
    idx <- match(tbl$feature, features$protein_group)
    if (anyNA(idx)) stop("top table contains protein group absent from annotation", call. = FALSE)
    cbind(tbl, features[idx, setdiff(names(features), "protein_group"), drop = FALSE])
  })
}

.limma_log_de_from_matrix <- function(mat, meta) {
  fd <- factorial_design(meta, add_batch = FALSE)
  fit <- fit_limma_log(mat, fd$design, fd$contrasts)
  list(fit = fit, design = list(add_batch = FALSE, design_cols = colnames(fd$design),
                                contrast_names = colnames(fd$contrasts),
                                residual_df = nrow(fd$design) - qr(fd$design)$rank))
}

run_proteome_de_24m <- function(proteomics_tbl, sample_key,
                                min_present = 2L, min_groups = 4L) {
  prep <- prepare_proteome_24m_matrix(proteomics_tbl, sample_key,
                                      min_present = min_present, min_groups = min_groups)
  de <- .limma_log_de_from_matrix(prep$matrix, prep$meta)
  top <- annotate_proteome_top_tables(de$fit$top, prep$features)
  list(
    n_samples = ncol(prep$matrix),
    n_features = nrow(prep$matrix),
    matrix = prep$matrix,
    meta = prep$meta,
    features = prep$features,
    top = top,
    design = de$design,
    filters = prep$counts,
    provenance = list(
      feature_id = "PG.ProteinGroups",
      aggregation = "sum raw positive peptide/PTM-row intensities by protein group before log2",
      transform = "log2 summed positive intensities",
      normalisation = "sample-wise median shift on log2 scale",
      prevalence = list(min_present = min_present, min_groups = min_groups)
    )
  )
}

# --- 24M bulk phosphosite DE (log2 median-normalised, no batch) ------------------------------

phosphosite_probability_column <- function(phospho_tbl) {
  candidates <- c("Phosphosite probability", "Phopshosite probability")
  hit <- intersect(candidates, names(phospho_tbl))
  if (!length(hit)) {
    stop("missing phosphosite probability column; expected one of: ",
         paste(candidates, collapse = ", "), call. = FALSE)
  }
  hit[1]
}

# Per-row phosphosite feature frame. feature = "row<n>|<PTM.CollapseKey>" is the unique id;
# site_id = "<gene>_<AA><loc>" (e.g. Mapt_S404), NA for missing/multi-gene rows (multiple
# feature rows can share one site_id).
phospho_feature_frame <- function(phospho_tbl) {
  stopifnot(is.data.frame(phospho_tbl))
  need <- c("PG.Genes", "PTM.SiteAA", "PTM.SiteLocation", "PTM.CollapseKey")
  stopifnot(all(need %in% names(phospho_tbl)))
  prob_col <- phosphosite_probability_column(phospho_tbl)
  row_index <- seq_len(nrow(phospho_tbl))
  collapse_key <- trimws(as.character(phospho_tbl[["PTM.CollapseKey"]]))
  collapse_blank <- is.na(collapse_key) | collapse_key == ""
  collapse_trace <- ifelse(collapse_blank, "<blank>", collapse_key)
  feature <- paste0("row", row_index, "|", collapse_trace)
  stopifnot(anyDuplicated(feature) == 0L)

  gene <- trimws(as.character(phospho_tbl[["PG.Genes"]]))
  aa <- trimws(as.character(phospho_tbl[["PTM.SiteAA"]]))
  loc <- trimws(as.character(phospho_tbl[["PTM.SiteLocation"]]))
  missing_gene <- is.na(gene) | gene == ""
  multi_gene <- grepl("[;,]", gene)
  missing_site <- is.na(aa) | aa == "" | is.na(loc) | loc == "" | tolower(loc) %in% c("na", "nan")
  site_id <- ifelse(!missing_gene & !multi_gene & !missing_site, paste0(gene, "_", aa, loc), NA_character_)
  probability <- suppressWarnings(as.numeric(phospho_tbl[[prob_col]]))
  out <- data.frame(
    feature = feature,
    original_row = row_index,
    collapse_key = collapse_key,
    collapse_key_blank = collapse_blank,
    gene = gene,
    site_aa = aa,
    site_location = loc,
    site_id = site_id,
    phosphosite_probability = probability,
    stringsAsFactors = FALSE
  )
  attr(out, "phospho_feature_counts") <- list(
    probability_col = prob_col,
    n_rows = nrow(out),
    n_blank_collapse_key = sum(collapse_blank),
    n_duplicate_collapse_key = sum(duplicated(collapse_key[!collapse_blank])),
    n_missing_gene = sum(missing_gene),
    n_multi_gene = sum(multi_gene & !missing_gene),
    n_missing_site = sum(missing_site),
    n_site_rows = sum(!is.na(site_id)),
    n_unique_sites = length(unique(site_id[!is.na(site_id)]))
  )
  out
}

positive_log2_matrix <- function(mat) {
  stopifnot(is.matrix(mat), is.numeric(mat))
  nonpositive <- !is.na(mat) & mat <= 0
  out <- mat
  out[nonpositive] <- NA_real_
  out <- log2(out)
  attr(out, "log2_counts") <- list(
    n_values = length(mat),
    n_missing_input = sum(is.na(mat)),
    n_nonpositive_to_na = sum(nonpositive),
    n_missing_output = sum(is.na(out))
  )
  out
}

prepare_phospho_24m_matrix <- function(phospho_tbl, sample_key,
                                       min_present = 2L, min_groups = 4L) {
  match <- match_24m_bulk_columns(phospho_tbl, sample_key, modality = "phospho")
  feat <- phospho_feature_frame(phospho_tbl)
  feature_counts <- attr(feat, "phospho_feature_counts")
  raw <- as.matrix(phospho_tbl[, match$columns, drop = FALSE])
  storage.mode(raw) <- "double"
  rownames(raw) <- feat$feature
  colnames(raw) <- rownames(match$meta)
  log_mat <- positive_log2_matrix(raw)
  norm <- median_normalise(log_mat)
  filt <- prevalence_filter(norm, match$meta$genotype,
                            min_present = min_present, min_groups = min_groups)
  feat <- feat[match(rownames(filt), feat$feature), , drop = FALSE]
  stopifnot(identical(feat$feature, rownames(filt)),
            identical(colnames(filt), rownames(match$meta)))
  list(
    matrix = filt,
    meta = match$meta,
    features = feat,
    matched = match$matched,
    counts = c(feature_counts, attr(log_mat, "log2_counts"),
               list(n_features_raw = nrow(raw),
                    n_features_filtered = nrow(filt),
                    n_filtered_site_rows = sum(!is.na(feat$site_id) & feat$site_id != ""),
                    n_filtered_unique_sites = length(unique(feat$site_id[!is.na(feat$site_id) &
                                                                          feat$site_id != ""])),
                    min_present = min_present,
                    min_groups = min_groups))
  )
}

annotate_phospho_top_tables <- function(top, features) {
  stopifnot(is.list(top), is.data.frame(features), "feature" %in% names(features))
  lapply(top, function(tbl) {
    stopifnot(is.data.frame(tbl), "feature" %in% names(tbl))
    idx <- match(tbl$feature, features$feature)
    if (anyNA(idx)) stop("top table contains feature absent from annotation", call. = FALSE)
    cbind(tbl, features[idx, setdiff(names(features), "feature"), drop = FALSE])
  })
}

run_phospho_de_24m <- function(phospho_tbl, sample_key,
                               min_present = 2L, min_groups = 4L) {
  prep <- prepare_phospho_24m_matrix(phospho_tbl, sample_key,
                                     min_present = min_present, min_groups = min_groups)
  fd <- factorial_design(prep$meta, add_batch = FALSE)
  fit <- fit_limma_log(prep$matrix, fd$design, fd$contrasts)
  top <- annotate_phospho_top_tables(fit$top, prep$features)
  list(
    n_samples = ncol(prep$matrix),
    n_features = nrow(prep$matrix),
    matrix = prep$matrix,
    meta = prep$meta,
    features = prep$features,
    top = top,
    design = list(add_batch = FALSE, design_cols = colnames(fd$design),
                  contrast_names = colnames(fd$contrasts),
                  residual_df = nrow(fd$design) - qr(fd$design)$rank),
    filters = prep$counts,
    provenance = list(
      transform = "log2 positive intensities; nonpositive values set to NA",
      normalisation = "sample-wise median shift on log2 scale",
      prevalence = list(min_present = min_present, min_groups = min_groups),
      feature_id = "row<original_row>|<PTM.CollapseKey>"
    )
  )
}
