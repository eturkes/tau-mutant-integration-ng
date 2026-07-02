# P4 cross-modality helpers. S1 covers GeoMx DE plus a deconvolution preflight only;
# later steps add bulk omics, clearance-axis, integration, and report bundles.

geomx_required_meta_cols <- function() {
  c("genotype", "slide_rep", "bio_rep", "roi", "SampleID",
    "ROI Coordinate X", "ROI Coordinate Y",
    "q_norm_qFactors", "NegGeoMean_Mm_R_NGS_WTA_v1.0", "nuclei")
}

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
    SampleID = as.character(md$SampleID),
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
            !anyNA(out$roi), !anyNA(out$SampleID),
            all(is.finite(out$x)), all(is.finite(out$y)),
            all(is.finite(out$q3_factor)), all(is.finite(out$neg_background)),
            all(is.finite(out$nuclei)))
  out$genotype <- factor(out$genotype, levels = genotype_levels)
  out$bio_unit <- factor(out$bio_unit)
  attr(out, "geomx_meta_provenance") <- list(
    n_aoi = nrow(out),
    n_bio_units = nlevels(out$bio_unit),
    n_slides = nlevels(out$slide),
    nuclei_sentinel_count = sum(out$nuclei < 0)
  )
  out
}

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
    v <- limma::voom(dge, design = design, plot = FALSE, block = block,
                     correlation = duplicate$consensus_correlation)
    fit0 <- limma::lmFit(v, design = design, block = block,
                         correlation = duplicate$consensus_correlation)
  } else {
    v <- limma::voom(dge, design = design, plot = FALSE)
    fit0 <- limma::lmFit(v, design = design)
  }
  fit <- limma::eBayes(limma::contrasts.fit(fit0, contrasts), robust = TRUE)
  list(
    kept = sum(keep),
    n_input_features = nrow(counts),
    n_samples = ncol(counts),
    design_cols = colnames(design),
    duplicate_correlation = duplicate,
    top = .geomx_top_tables(fit, contrasts)
  )
}

collapse_geomx_by_bio_unit <- function(counts, meta) {
  stopifnot(identical(colnames(counts), rownames(meta)),
            all(c("bio_unit", "genotype") %in% names(meta)))
  unit <- as.character(meta$bio_unit)
  idx <- split(seq_along(unit), unit)
  collapsed <- vapply(names(idx), function(u) Matrix::rowSums(counts[, idx[[u]], drop = FALSE]),
                      numeric(nrow(counts)))
  collapsed <- matrix(collapsed, nrow = nrow(counts),
                      dimnames = list(rownames(counts), names(idx)))
  unit_meta <- lapply(names(idx), function(u) {
    rows <- meta[idx[[u]], , drop = FALSE]
    if (length(unique(as.character(rows$genotype))) != 1L) {
      stop("GeoMx bio_unit has non-constant genotype: ", u, call. = FALSE)
    }
    data.frame(genotype = as.character(rows$genotype[1]), bio_unit = u,
               n_aoi = nrow(rows), row.names = u, stringsAsFactors = FALSE)
  })
  unit_meta <- do.call(rbind, unit_meta)
  unit_meta$genotype <- factor(unit_meta$genotype, levels = genotype_levels)
  list(counts = collapsed, meta = unit_meta)
}

fit_geomx_collapsed_sensitivity <- function(counts, meta, min_count = 5) {
  collapsed <- collapse_geomx_by_bio_unit(counts, meta)
  fd <- try(geomx_slide_design(data.frame(genotype = collapsed$meta$genotype,
                                          slide = factor("collapsed"),
                                          row.names = rownames(collapsed$meta)),
                               include_slide = FALSE),
            silent = TRUE)
  if (inherits(fd, "try-error")) {
    return(list(status = "skipped", reason = conditionMessage(attr(fd, "condition")),
                n_bio_units = ncol(collapsed$counts),
                genotype_counts = table(collapsed$meta$genotype)))
  }
  if (nrow(fd$design) <= ncol(fd$design)) {
    return(list(status = "skipped", reason = "GeoMx collapsed design has no residual degrees of freedom",
                n_bio_units = ncol(collapsed$counts),
                genotype_counts = table(collapsed$meta$genotype)))
  }
  fit <- .fit_geomx_voom(collapsed$counts, fd$design, fd$contrasts, block = NULL,
                         min_count = min_count)
  c(list(status = "fit", reason = NA_character_, n_bio_units = ncol(collapsed$counts),
         genotype_counts = table(collapsed$meta$genotype)), fit)
}

repo_package_available <- function(package, repos = getOption("repos")) {
  warnings <- character()
  messages <- character()
  result <- tryCatch(
    withCallingHandlers({
      ap <- utils::available.packages(repos = repos)
      hit <- package %in% rownames(ap)
      list(package = package, available = hit,
           version = if (hit) unname(ap[package, "Version"]) else NA_character_,
           repos = unname(repos), error = NA_character_)
    }, warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }, message = function(m) {
      messages <<- c(messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    }),
    error = function(e) {
      list(package = package, available = NA, version = NA_character_,
           repos = unname(repos), error = conditionMessage(e))
    }
  )
  result$warnings <- warnings
  result$messages <- messages
  result
}

geomx_decon_preflight <- function(meta, counts, profile = NULL, profile_corr_threshold = 0.95,
                                  spatialdecon = NULL) {
  stopifnot(is.data.frame(meta), identical(colnames(counts), rownames(meta)),
            all(c("q3_factor", "neg_background", "nuclei") %in% names(meta)))
  if (is.null(spatialdecon)) spatialdecon <- repo_package_available("SpatialDecon")

  q3_ok <- all(is.finite(meta$q3_factor) & meta$q3_factor > 0)
  bg_ok <- all(is.finite(meta$neg_background) & meta$neg_background > 0)
  nuclei_sentinel <- sum(meta$nuclei < 0)
  profile_tested <- !is.null(profile)
  max_profile_cor <- NA_real_
  profile_ok <- FALSE
  if (profile_tested) {
    stopifnot(is.matrix(profile), ncol(profile) >= 2L, nrow(profile) >= 2L)
    cc <- stats::cor(profile, use = "pairwise.complete.obs")
    off <- abs(cc[upper.tri(cc)])
    max_profile_cor <- if (length(off)) max(off, na.rm = TRUE) else NA_real_
    profile_ok <- is.finite(max_profile_cor) && max_profile_cor < profile_corr_threshold
  }

  reasons <- character()
  if (!q3_ok) reasons <- c(reasons, "Q3 normalisation factors are missing, non-finite, or non-positive")
  if (!bg_ok) reasons <- c(reasons, "negative-probe background is missing, non-finite, or non-positive")
  if (nuclei_sentinel > 0L) {
    reasons <- c(reasons, sprintf("nuclei contains %d sentinel value(s); absolute nuclei rescaling stays disabled",
                                  nuclei_sentinel))
  }
  if (!isTRUE(spatialdecon$available)) {
    reasons <- c(reasons, "SpatialDecon is not confirmed available from the pinned repositories")
  }
  if (!profile_tested) {
    reasons <- c(reasons, "no compact reference profile was built in S1; deconvolution deferred to S3")
  } else if (!profile_ok) {
    reasons <- c(reasons, sprintf("reference profile collinearity %.3f exceeds threshold %.3f",
                                  max_profile_cor, profile_corr_threshold))
  }

  hard_block <- !q3_ok || !bg_ok || (profile_tested && !profile_ok)
  status <- if (hard_block) "blocked" else if (isTRUE(spatialdecon$available) && profile_ok) "earned" else "defer"
  list(
    status = status,
    reasons = reasons,
    background = list(q3_ok = q3_ok, negative_background_ok = bg_ok,
                      q3_factor_range = range(meta$q3_factor),
                      neg_background_range = range(meta$neg_background),
                      instruction = "scale negative-probe background onto the Q3-normalised expression scale before deconvolution"),
    nuclei = list(n_sentinel = nuclei_sentinel, absolute_rescaling_enabled = FALSE),
    reference = list(profile_tested = profile_tested, profile_ok = profile_ok,
                     max_abs_correlation = max_profile_cor,
                     threshold = profile_corr_threshold),
    spatialdecon = spatialdecon,
    memory = list(count_matrix_mb = as.numeric(utils::object.size(counts)) / 1024^2,
                  estimated_peak_mb = as.numeric(utils::object.size(counts)) / 1024^2 * 6)
  )
}

fit_geomx_de <- function(counts, meta, min_count = 5) {
  stopifnot(identical(colnames(counts), rownames(meta)),
            all(c("genotype", "slide", "bio_unit") %in% names(meta)))
  fd <- geomx_slide_design(meta, include_slide = TRUE)
  primary <- .fit_geomx_voom(counts, fd$design, fd$contrasts, block = meta$bio_unit,
                             min_count = min_count)
  unblocked <- .fit_geomx_voom(counts, fd$design, fd$contrasts, block = NULL,
                               min_count = min_count)
  collapsed <- fit_geomx_collapsed_sensitivity(counts, meta, min_count = min_count)
  list(
    n_aoi = ncol(counts),
    n_bio_units = length(unique(as.character(meta$bio_unit))),
    primary = c(list(status = "fit", model = "voom_tmm_slide_duplicateCorrelation"), primary),
    sensitivity = list(
      unblocked = c(list(status = "fit", model = "voom_tmm_slide_unblocked"), unblocked),
      collapsed_bio_unit = collapsed
    )
  )
}

run_geomx_de <- function(geomx, min_count = 5) {
  counts <- geomx_count_matrix(geomx)
  meta <- geomx_meta(geomx)
  stopifnot(identical(colnames(counts), rownames(meta)))
  de <- fit_geomx_de(counts, meta, min_count = min_count)
  de$decon_preflight <- geomx_decon_preflight(meta, counts)
  de$provenance <- list(
    counts = attr(counts, "geomx_count_provenance"),
    meta = attr(meta, "geomx_meta_provenance"),
    thresholds = list(min_count = min_count)
  )
  de
}
