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

geomx_q3_scaled_background <- function(meta) {
  stopifnot(is.data.frame(meta), all(c("q3_factor", "neg_background") %in% names(meta)))
  q3 <- as.numeric(meta$q3_factor)
  bg <- as.numeric(meta$neg_background)
  stopifnot(length(q3) == length(bg),
            all(is.finite(q3) & q3 > 0),
            all(is.finite(bg) & bg > 0))
  out <- bg / q3
  attr(out, "scale") <- "negative-probe background divided by q_norm_qFactors"
  out
}

profile_collinearity <- function(profile) {
  stopifnot(is.matrix(profile), nrow(profile) >= 2L, ncol(profile) >= 2L,
            all(is.finite(profile)))
  cc <- suppressWarnings(stats::cor(profile, use = "pairwise.complete.obs"))
  off <- abs(cc[upper.tri(cc)])
  max_cor <- if (length(off) && any(is.finite(off))) max(off, na.rm = TRUE) else NA_real_
  list(n_features = nrow(profile), n_profiles = ncol(profile),
       max_abs_correlation = max_cor)
}

geomx_decon_preflight <- function(meta, counts, profile = NULL, profile_corr_threshold = 0.95,
                                  spatialdecon = NULL) {
  stopifnot(is.data.frame(meta), identical(colnames(counts), rownames(meta)),
            all(c("q3_factor", "neg_background", "nuclei") %in% names(meta)))
  if (is.null(spatialdecon)) spatialdecon <- repo_package_available("SpatialDecon")

  bg_scaled <- try(geomx_q3_scaled_background(meta), silent = TRUE)
  q3_ok <- all(is.finite(meta$q3_factor) & meta$q3_factor > 0)
  bg_ok <- !inherits(bg_scaled, "try-error") && all(is.finite(bg_scaled) & bg_scaled > 0)
  nuclei_sentinel <- sum(meta$nuclei < 0)
  profile_tested <- !is.null(profile)
  max_profile_cor <- NA_real_
  profile_ok <- FALSE
  if (profile_tested) {
    pc <- profile_collinearity(profile)
    max_profile_cor <- pc$max_abs_correlation
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
                      q3_scaled_background_range = if (bg_ok) range(bg_scaled) else c(NA_real_, NA_real_),
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

.geomx_abundance_top_tables <- function(fit, contrasts) {
  out <- lapply(colnames(contrasts), function(cn) {
    tt <- tibble::rownames_to_column(
      limma::topTable(fit, coef = cn, number = Inf, sort.by = "none", confint = TRUE),
      "feature"
    )
    tt$contrast <- cn
    tt[, c("feature", "contrast", "logFC", "P.Value", "adj.P.Val", "t", "CI.L", "CI.R",
           setdiff(names(tt), c("feature", "contrast", "logFC", "P.Value", "adj.P.Val", "t", "CI.L", "CI.R"))),
       drop = FALSE]
  })
  names(out) <- colnames(contrasts)
  out
}

fit_geomx_abundance_de <- function(abundance, meta, offset = 1e-4) {
  stopifnot(is.matrix(abundance), !is.null(rownames(abundance)), !is.null(colnames(abundance)),
            is.data.frame(meta), identical(colnames(abundance), rownames(meta)),
            all(c("genotype", "slide", "bio_unit") %in% names(meta)),
            is.numeric(offset), length(offset) == 1L, is.finite(offset), offset > 0,
            all(is.finite(abundance)), all(abundance >= 0))
  fd <- geomx_slide_design(meta, include_slide = TRUE)
  if (nrow(fd$design) <= ncol(fd$design)) {
    stop("GeoMx abundance design has no residual degrees of freedom", call. = FALSE)
  }
  log_abund <- log(abundance + offset)
  corfit <- limma::duplicateCorrelation(log_abund, design = fd$design, block = meta$bio_unit)
  if (!is.finite(corfit$consensus.correlation)) {
    stop("GeoMx abundance duplicateCorrelation returned non-finite consensus correlation",
         call. = FALSE)
  }
  fit0 <- limma::lmFit(log_abund, design = fd$design, block = meta$bio_unit,
                       correlation = corfit$consensus.correlation)
  fit <- limma::eBayes(limma::contrasts.fit(fit0, fd$contrasts), robust = TRUE)
  fit_unblocked <- limma::eBayes(limma::contrasts.fit(
    limma::lmFit(log_abund, design = fd$design), fd$contrasts
  ), robust = TRUE)
  list(
    status = "fit",
    model = "log_beta_slide_duplicateCorrelation",
    n_aoi = ncol(abundance),
    n_features = nrow(abundance),
    offset = offset,
    duplicate_correlation = list(used = TRUE,
                                 consensus_correlation = unname(corfit$consensus.correlation)),
    top = .geomx_abundance_top_tables(fit, fd$contrasts),
    sensitivity = list(
      unblocked = list(status = "fit", model = "log_beta_slide_unblocked",
                       top = .geomx_abundance_top_tables(fit_unblocked, fd$contrasts))
    )
  )
}

# ---- P4-S2: bulk proteome + protein-corrected phospho -------------------------------

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

  run_index <- run_index_factorial_design(prep$meta)
  if (identical(run_index$status, "fit")) {
    fit_ri <- fit_limma_log(prep$matrix, run_index$design, run_index$contrasts)
    run_index$top <- annotate_proteome_top_tables(fit_ri$top, prep$features)
  }
  list(
    n_samples = ncol(prep$matrix),
    n_features = nrow(prep$matrix),
    matrix = prep$matrix,
    meta = prep$meta,
    features = prep$features,
    top = top,
    design = de$design,
    run_index = run_index,
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

prepare_phospho_corrected_24m_matrix <- function(phospho_tbl, sample_key, proteome_de_24m,
                                                 min_present = 2L, min_groups = 4L) {
  stopifnot(is.data.frame(phospho_tbl), is.list(proteome_de_24m),
            is.matrix(proteome_de_24m$matrix), is.data.frame(proteome_de_24m$meta),
            is.data.frame(proteome_de_24m$features))
  phospho <- prepare_phospho_24m_matrix(phospho_tbl, sample_key,
                                        min_present = min_present, min_groups = min_groups)
  if (!identical(rownames(phospho$meta), rownames(proteome_de_24m$meta)) ||
      !identical(colnames(phospho$matrix), colnames(proteome_de_24m$matrix))) {
    stop("phospho and proteome 24M sample order must be identical for protein correction",
         call. = FALSE)
  }
  stopifnot("PG.ProteinGroups" %in% names(phospho_tbl),
            "original_row" %in% names(phospho$features))
  parent <- trimws(as.character(phospho_tbl[["PG.ProteinGroups"]][phospho$features$original_row]))
  missing_parent_id <- is.na(parent) | parent == ""
  parent_idx <- match(parent, rownames(proteome_de_24m$matrix))
  matched <- !missing_parent_id & !is.na(parent_idx)
  if (!any(matched)) stop("no phosphosites have a matched parent protein group", call. = FALSE)

  corrected <- phospho$matrix[matched, , drop = FALSE] -
    proteome_de_24m$matrix[parent_idx[matched], , drop = FALSE]
  rownames(corrected) <- rownames(phospho$matrix)[matched]
  corrected <- prevalence_filter(corrected, phospho$meta$genotype,
                                 min_present = min_present, min_groups = min_groups)

  features <- phospho$features[match(rownames(corrected), phospho$features$feature), , drop = FALSE]
  parent_kept <- parent[match(features$feature, phospho$features$feature)]
  features$parent_protein_group <- parent_kept
  pfeat <- proteome_de_24m$features[match(parent_kept, proteome_de_24m$features$protein_group), ,
                                    drop = FALSE]
  features$parent_gene_symbols <- pfeat$gene_symbols
  features$parent_n_raw_rows <- pfeat$n_raw_rows
  stopifnot(identical(features$feature, rownames(corrected)),
            identical(colnames(corrected), rownames(phospho$meta)),
            !anyNA(features$parent_protein_group))

  list(
    matrix = corrected,
    meta = phospho$meta,
    features = features,
    matched = phospho$matched,
    counts = c(phospho$counts,
               list(n_phospho_features_for_correction = nrow(phospho$matrix),
                    n_missing_parent_id = sum(missing_parent_id),
                    n_parent_not_in_filtered_proteome = sum(!missing_parent_id & is.na(parent_idx)),
                    n_parent_matched = sum(matched),
                    n_unique_parent_matched = length(unique(parent[matched])),
                    n_values_corrected = length(corrected),
                    n_missing_corrected_output = sum(is.na(corrected)),
                    n_features_corrected_filtered = nrow(corrected),
                    correction_min_present = min_present,
                    correction_min_groups = min_groups))
  )
}

run_phospho_corrected_24m <- function(phospho_tbl, sample_key, proteome_de_24m,
                                      min_present = 2L, min_groups = 4L) {
  prep <- prepare_phospho_corrected_24m_matrix(phospho_tbl, sample_key, proteome_de_24m,
                                               min_present = min_present, min_groups = min_groups)
  de <- .limma_log_de_from_matrix(prep$matrix, prep$meta)
  top <- annotate_phospho_top_tables(de$fit$top, prep$features)

  run_index <- run_index_factorial_design(prep$meta)
  if (identical(run_index$status, "fit")) {
    fit_ri <- fit_limma_log(prep$matrix, run_index$design, run_index$contrasts)
    run_index$top <- annotate_phospho_top_tables(fit_ri$top, prep$features)
  }
  list(
    n_samples = ncol(prep$matrix),
    n_features = nrow(prep$matrix),
    matrix = prep$matrix,
    meta = prep$meta,
    features = prep$features,
    top = top,
    design = de$design,
    run_index = run_index,
    filters = prep$counts,
    provenance = list(
      feature_id = "row<original_row>|<PTM.CollapseKey>",
      correction = "median-normalised phosphosite log2 intensity minus matched parent-protein log2 intensity",
      parent_match = "exact PG.ProteinGroups match against filtered proteome_de_24m matrix",
      normalisation = "inherits phosphosite and parent protein sample-wise median shifts before subtraction",
      prevalence = list(min_present = min_present, min_groups = min_groups)
    )
  )
}

bulk_significant_counts <- function(top, layer) {
  stopifnot(is.list(top))
  out <- do.call(rbind, lapply(names(top), function(cn) {
    tt <- top[[cn]]
    stopifnot(all(c("adj.P.Val", "logFC") %in% names(tt)))
    data.frame(
      layer = layer,
      contrast = cn,
      n_features = nrow(tt),
      n_fdr_0_05 = sum(is.finite(tt$adj.P.Val) & tt$adj.P.Val < 0.05),
      n_fdr_0_10 = sum(is.finite(tt$adj.P.Val) & tt$adj.P.Val < 0.10),
      n_up_fdr_0_10 = sum(is.finite(tt$adj.P.Val) & tt$adj.P.Val < 0.10 & tt$logFC > 0),
      n_down_fdr_0_10 = sum(is.finite(tt$adj.P.Val) & tt$adj.P.Val < 0.10 & tt$logFC < 0),
      stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL
  out
}

bulk_run_index_summary <- function(primary_top, run_index, layer, alpha = 0.10) {
  stopifnot(is.list(primary_top), is.list(run_index))
  if (!identical(run_index$status, "fit") || !is.list(run_index$top)) {
    return(data.frame(layer = layer, contrast = names(primary_top),
                      status = run_index$status %||% "skipped",
                      reason = run_index$reason %||% NA_character_,
                      n_primary_sig = NA_integer_, n_lost_or_flipped = NA_integer_,
                      stringsAsFactors = FALSE))
  }
  out <- lapply(names(primary_top), function(cn) {
    p <- primary_top[[cn]]
    r <- run_index$top[[cn]]
    stopifnot(all(c("feature", "adj.P.Val", "logFC") %in% names(p)),
              all(c("feature", "adj.P.Val", "logFC") %in% names(r)))
    idx <- match(p$feature, r$feature)
    if (anyNA(idx)) stop("run-index top table missing primary feature(s)", call. = FALSE)
    p_sig <- is.finite(p$adj.P.Val) & p$adj.P.Val < alpha
    r_sig <- is.finite(r$adj.P.Val[idx]) & r$adj.P.Val[idx] < alpha
    flip <- is.finite(p$logFC) & is.finite(r$logFC[idx]) & sign(p$logFC) != sign(r$logFC[idx])
    data.frame(layer = layer, contrast = cn, status = "fit", reason = NA_character_,
               n_primary_sig = sum(p_sig),
               n_lost_or_flipped = sum(p_sig & (!r_sig | flip)),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

bulk_anchor_dictionary <- function() {
  data.frame(
    symbol = c("Gsk3b", "Mapt", "App", "Apoe", "Trem2", "Cd74", "Mertk", "Pros1",
               "C1qa", "C1qb", "C1qc", "C3", "Syn1", "Syp", "Snap25", "Dlg4", "Grin1"),
    anchor_class = c("gsk3b", "tau", rep("clearance", 6), rep("complement", 4),
                     rep("synaptic", 5)),
    stringsAsFactors = FALSE
  )
}

.row_symbols <- function(tbl) {
  if ("gene_symbols" %in% names(tbl)) {
    vapply(tbl$gene_symbols, function(x) paste(.split_gene_symbols(x), collapse = ";"), character(1))
  } else if ("gene" %in% names(tbl)) {
    as.character(tbl$gene)
  } else {
    rep(NA_character_, nrow(tbl))
  }
}

bulk_anchor_rows <- function(top, layer, anchors = bulk_anchor_dictionary()) {
  stopifnot(is.list(top), is.data.frame(anchors), all(c("symbol", "anchor_class") %in% names(anchors)))
  rows <- list()
  for (cn in names(top)) {
    tt <- top[[cn]]
    syms <- .row_symbols(tt)
    hit <- vapply(syms, function(x) {
      any(.split_gene_symbols(x) %in% anchors$symbol)
    }, logical(1))
    if (!any(hit)) next
    d <- tt[hit, , drop = FALSE]
    ds <- syms[hit]
    hit_symbols <- vapply(ds, function(x) {
      paste(intersect(.split_gene_symbols(x), anchors$symbol), collapse = ";")
    }, character(1))
    hit_class <- vapply(hit_symbols, function(x) {
      paste(unique(anchors$anchor_class[match(.split_gene_symbols(x), anchors$symbol)]),
            collapse = ";")
    }, character(1))
    feature <- if ("feature" %in% names(d)) d$feature else rep(NA_character_, nrow(d))
    site <- if ("site_id" %in% names(d)) d$site_id else rep(NA_character_, nrow(d))
    rows[[length(rows) + 1L]] <- data.frame(
      layer = layer,
      contrast = cn,
      feature = feature,
      site_id = site,
      symbols = ds,
      anchor_symbols = hit_symbols,
      anchor_class = hit_class,
      logFC = d$logFC,
      t = d$t,
      P.Value = d$P.Value,
      adj.P.Val = d$adj.P.Val,
      stringsAsFactors = FALSE
    )
  }
  if (!length(rows)) {
    return(data.frame(layer = character(), contrast = character(), feature = character(),
                      site_id = character(), symbols = character(), anchor_symbols = character(),
                      anchor_class = character(), logFC = numeric(), t = numeric(),
                      P.Value = numeric(), adj.P.Val = numeric(), stringsAsFactors = FALSE))
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

bulk_omics_summary_data <- function(proteome_de_24m, phospho_de_24m, phospho_corrected_24m,
                                    alpha = 0.10) {
  layers <- list(
    proteome = proteome_de_24m,
    phospho_raw = phospho_de_24m,
    phospho_corrected = phospho_corrected_24m
  )
  stopifnot(all(vapply(layers, is.list, logical(1))))
  feature_counts <- do.call(rbind, lapply(names(layers), function(layer) {
    x <- layers[[layer]]
    data.frame(
      layer = layer,
      n_samples = x$n_samples,
      n_features = x$n_features,
      n_missing_output = x$filters$n_missing_corrected_output %||% x$filters$n_missing_output %||% NA_integer_,
      n_nonpositive_to_na = x$filters$n_nonpositive_to_na %||% NA_integer_,
      n_parent_matched = x$filters$n_parent_matched %||% NA_integer_,
      n_parent_not_in_filtered_proteome = x$filters$n_parent_not_in_filtered_proteome %||% NA_integer_,
      stringsAsFactors = FALSE
    )
  }))
  significant_counts <- do.call(rbind, lapply(names(layers), function(layer) {
    bulk_significant_counts(layers[[layer]]$top, layer)
  }))
  run_index <- do.call(rbind, lapply(names(layers), function(layer) {
    bulk_run_index_summary(layers[[layer]]$top, layers[[layer]]$run_index, layer, alpha = alpha)
  }))
  anchors <- do.call(rbind, lapply(names(layers), function(layer) {
    bulk_anchor_rows(layers[[layer]]$top, layer)
  }))
  anchor_dict <- bulk_anchor_dictionary()
  anchor_coverage <- do.call(rbind, lapply(seq_len(nrow(anchor_dict)), function(i) {
    sym <- anchor_dict$symbol[i]
    data.frame(
      symbol = sym,
      anchor_class = anchor_dict$anchor_class[i],
      n_rows = if (nrow(anchors)) sum(vapply(strsplit(anchors$anchor_symbols, ";", fixed = TRUE),
                                             function(x) sym %in% x, logical(1))) else 0L,
      stringsAsFactors = FALSE
    )
  }))
  out <- list(
    feature_counts = feature_counts,
    significant_counts = significant_counts,
    run_index = run_index,
    anchors = anchors,
    anchor_coverage = anchor_coverage,
    provenance = list(
      alpha = alpha,
      layers = names(layers),
      bulk_context = "24M bulk hippocampus proteome/phosphoproteome; not microglia-sorted",
      summary = "compact S2 report/input summary; margins computed from source targets"
    )
  )
  stopifnot(nrow(out$feature_counts) == 3L,
            all(out$feature_counts$n_samples == 16L),
            all(c("Gsk3b", "Mapt") %in% out$anchor_coverage$symbol))
  out
}

# ---- P4-S3: spatial-composition decision + clearance-axis CCC-lite -------------------

clearance_axis_dictionary <- function() {
  data.frame(
    symbol = c("Apoe", "Trem2", "App", "Cd74", "Pros1", "Mertk",
               "C1qa", "C1qb", "C1qc", "C3",
               "Syn1", "Syp", "Snap25", "Dlg4", "Grin1"),
    axis = c(rep("clearance", 6), rep("complement", 4), rep("synaptic", 5)),
    role = c("ligand", "receptor", "ligand", "receptor", "ligand", "receptor",
             rep("component", 4), rep("marker", 5)),
    pair = c("Apoe_Trem2", "Apoe_Trem2", "App_Cd74", "App_Cd74",
             "Pros1_Mertk", "Pros1_Mertk", rep(NA_character_, 9)),
    stringsAsFactors = FALSE
  )
}

.empty_clearance_rows <- function() {
  data.frame(axis = character(), pair = character(), role = character(), symbol = character(),
             modality = character(), layer = character(), contrast = character(),
             feature = character(), site_id = character(), logFC = numeric(), t = numeric(),
             P.Value = numeric(), adj.P.Val = numeric(), significant = logical(),
             direction = integer(), stringsAsFactors = FALSE)
}

.clearance_rows_from_vectors <- function(tbl, symbols, anchors, modality, layer, contrast,
                                         feature, site_id = NULL, alpha = 0.10) {
  stopifnot(is.data.frame(tbl), length(symbols) == nrow(tbl), length(feature) == nrow(tbl),
            all(c("symbol", "axis", "role", "pair") %in% names(anchors)))
  if (is.null(site_id)) site_id <- rep(NA_character_, nrow(tbl))
  out <- list()
  k <- 0L
  for (i in seq_len(nrow(tbl))) {
    hit <- intersect(.split_gene_symbols(symbols[i]), anchors$symbol)
    if (!length(hit)) next
    for (sym in hit) {
      ai <- match(sym, anchors$symbol)
      k <- k + 1L
      out[[k]] <- data.frame(
        axis = anchors$axis[ai],
        pair = anchors$pair[ai],
        role = anchors$role[ai],
        symbol = sym,
        modality = modality,
        layer = layer,
        contrast = contrast,
        feature = as.character(feature[i]),
        site_id = as.character(site_id[i] %||% NA_character_),
        logFC = as.numeric(tbl$logFC[i]),
        t = as.numeric(tbl$t[i]),
        P.Value = as.numeric(tbl$P.Value[i]),
        adj.P.Val = as.numeric(tbl$adj.P.Val[i]),
        significant = is.finite(as.numeric(tbl$adj.P.Val[i])) &&
          as.numeric(tbl$adj.P.Val[i]) < alpha,
        direction = as.integer(sign(as.numeric(tbl$logFC[i]))),
        stringsAsFactors = FALSE
      )
    }
  }
  if (!length(out)) return(.empty_clearance_rows())
  ans <- do.call(rbind, out)
  rownames(ans) <- NULL
  ans
}

clearance_rows_from_rna_top <- function(top, symbol_map, layer, anchors,
                                        modality = "snRNAseq_microglia", alpha = 0.10) {
  stopifnot(is.list(top), is.data.frame(symbol_map),
            all(c("ensembl", "symbol") %in% names(symbol_map)))
  map <- symbol_map[!is.na(symbol_map$ensembl) & !is.na(symbol_map$symbol) &
                      symbol_map$ensembl != "" & symbol_map$symbol != "",
                    c("ensembl", "symbol")]
  map <- map[!duplicated(map$ensembl), , drop = FALSE]
  rows <- lapply(intersect(mechanism_contrasts(), names(top)), function(cn) {
    tt <- top[[cn]]
    stopifnot(is.data.frame(tt),
              all(c("gene", "logFC", "t", "P.Value", "adj.P.Val") %in% names(tt)))
    symbols <- map$symbol[match(as.character(tt$gene), map$ensembl)]
    .clearance_rows_from_vectors(tt, symbols, anchors, modality, layer, cn,
                                 feature = tt$gene, alpha = alpha)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

clearance_rows_from_geomx <- function(geomx_de, anchors, alpha = 0.10) {
  stopifnot(is.list(geomx_de), is.list(geomx_de$primary), is.list(geomx_de$primary$top))
  rows <- lapply(intersect(mechanism_contrasts(), names(geomx_de$primary$top)), function(cn) {
    tt <- geomx_de$primary$top[[cn]]
    stopifnot(is.data.frame(tt),
              all(c("symbol", "logFC", "t", "P.Value", "adj.P.Val") %in% names(tt)))
    .clearance_rows_from_vectors(tt, tt$symbol, anchors, "GeoMx_spatial", "primary", cn,
                                 feature = tt$symbol, alpha = alpha)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

clearance_rows_from_bulk_summary <- function(bulk_omics_summary, anchors, alpha = 0.10) {
  stopifnot(is.list(bulk_omics_summary), is.data.frame(bulk_omics_summary$anchors))
  b <- bulk_omics_summary$anchors
  if (!nrow(b)) return(.empty_clearance_rows())
  stopifnot(all(c("layer", "contrast", "feature", "site_id", "anchor_symbols",
                  "logFC", "t", "P.Value", "adj.P.Val") %in% names(b)))
  rows <- lapply(seq_len(nrow(b)), function(i) {
    .clearance_rows_from_vectors(b[i, , drop = FALSE], b$anchor_symbols[i], anchors,
                                 "bulk_hippocampus", b$layer[i], b$contrast[i],
                                 feature = b$feature[i], site_id = b$site_id[i],
                                 alpha = alpha)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.collapse_clearance_symbol <- function(d) {
  d <- d[is.finite(d$logFC) & is.finite(d$t), , drop = FALSE]
  if (!nrow(d)) return(d[FALSE, , drop = FALSE])
  ord <- order(d$adj.P.Val, -abs(d$t), -abs(d$logFC), d$layer, d$feature,
               method = "radix", na.last = TRUE)
  d[ord[1], , drop = FALSE]
}

clearance_pair_support <- function(measured, anchors = clearance_axis_dictionary(),
                                   alpha = 0.10) {
  pairs <- na.omit(unique(anchors$pair))
  contrasts <- mechanism_contrasts()
  out <- list()
  k <- 0L
  for (pair in pairs) {
    pair_symbols <- anchors$symbol[anchors$pair %in% pair]
    for (cn in contrasts) {
      d <- measured[measured$pair %in% pair & measured$contrast %in% cn, , drop = FALSE]
      modalities <- sort(unique(d$modality), method = "radix")
      coherent_supported <- character()
      coherent_measured <- character()
      for (mod in modalities) {
        dm <- d[d$modality == mod, , drop = FALSE]
        reps <- do.call(rbind, lapply(pair_symbols, function(sym) {
          .collapse_clearance_symbol(dm[dm$symbol == sym, , drop = FALSE])
        }))
        if (!is.data.frame(reps) || nrow(reps) != length(pair_symbols)) next
        signs <- sign(reps$logFC)
        coherent <- all(signs != 0) && length(unique(signs)) == 1L
        if (coherent) {
          coherent_measured <- c(coherent_measured, mod)
          if (any(reps$significant, na.rm = TRUE)) coherent_supported <- c(coherent_supported, mod)
        }
      }
      mg <- d[d$modality == "snRNAseq_microglia" & d$layer == "whole_microglia", , drop = FALSE]
      mg_reps <- do.call(rbind, lapply(pair_symbols, function(sym) {
        .collapse_clearance_symbol(mg[mg$symbol == sym, , drop = FALSE])
      }))
      microglia_strong <- is.data.frame(mg_reps) && nrow(mg_reps) == length(pair_symbols) &&
        all(sign(mg_reps$logFC) != 0) && length(unique(sign(mg_reps$logFC))) == 1L &&
        all(is.finite(mg_reps$adj.P.Val) & mg_reps$adj.P.Val < alpha)
      non_micro_supported <- setdiff(coherent_supported, "snRNAseq_microglia")
      earned <- length(unique(coherent_supported)) >= 2L ||
        (length(non_micro_supported) >= 1L && microglia_strong)
      k <- k + 1L
      out[[k]] <- data.frame(
        pair = pair,
        contrast = cn,
        n_sides_measured = length(unique(d$symbol)),
        modalities_measured = paste(sort(unique(d$modality), method = "radix"), collapse = ";"),
        coherent_measured_modalities = paste(sort(unique(coherent_measured), method = "radix"),
                                             collapse = ";"),
        coherent_supported_modalities = paste(sort(unique(coherent_supported), method = "radix"),
                                              collapse = ";"),
        n_coherent_supported_modalities = length(unique(coherent_supported)),
        microglia_strong = microglia_strong,
        status = if (earned) "earned" else "not_earned",
        rule = "earned requires two coherent supported modalities, or one non-microglia modality plus a strong whole-microglia anchor",
        stringsAsFactors = FALSE
      )
    }
  }
  ans <- do.call(rbind, out)
  rownames(ans) <- NULL
  ans
}

clearance_axis_coverage <- function(measured, anchors = clearance_axis_dictionary()) {
  rows <- lapply(seq_len(nrow(anchors)), function(i) {
    d <- measured[measured$symbol == anchors$symbol[i], , drop = FALSE]
    data.frame(
      symbol = anchors$symbol[i],
      axis = anchors$axis[i],
      role = anchors$role[i],
      pair = anchors$pair[i],
      measured = nrow(d) > 0L,
      n_rows = nrow(d),
      n_modalities = length(unique(d$modality)),
      modalities = paste(sort(unique(d$modality), method = "radix"), collapse = ";"),
      n_significant = sum(d$significant, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

clearance_axis_modality_support <- function(measured) {
  if (!nrow(measured)) {
    return(data.frame(axis = character(), contrast = character(), modality = character(),
                      n_symbols = integer(), n_significant = integer(),
                      direction_balance = integer(), stringsAsFactors = FALSE))
  }
  key <- interaction(measured$axis, measured$contrast, measured$modality,
                     drop = TRUE, lex.order = TRUE)
  rows <- lapply(split(measured, key), function(d) {
    sig <- d[d$significant, , drop = FALSE]
    data.frame(
      axis = d$axis[1],
      contrast = d$contrast[1],
      modality = d$modality[1],
      n_symbols = length(unique(d$symbol)),
      n_significant = length(unique(sig$symbol)),
      direction_balance = sum(sign(d$logFC), na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out <- out[order(out$axis, match(out$contrast, mechanism_contrasts()), out$modality,
                   method = "radix", na.last = TRUE), , drop = FALSE]
  rownames(out) <- NULL
  out
}

synaptic_gene_set_rows <- function(mechanism_gene_sets,
                                   anchors = clearance_axis_dictionary(),
                                   pattern = "SYNAP") {
  stopifnot(is.list(mechanism_gene_sets), is.list(mechanism_gene_sets$sets))
  syn <- anchors$symbol[anchors$axis == "synaptic"]
  rows <- list()
  k <- 0L
  for (collection in names(mechanism_gene_sets$sets)) {
    sets <- mechanism_gene_sets$sets[[collection]]
    hit <- grep(pattern, names(sets), value = TRUE, ignore.case = TRUE)
    if (!length(hit)) next
    for (set in hit) {
      genes <- unique(sets[[set]])
      overlap <- intersect(syn, genes)
      k <- k + 1L
      rows[[k]] <- data.frame(collection = collection, set = set,
                              size = length(genes),
                              n_synaptic_anchor_overlap = length(overlap),
                              synaptic_anchor_overlap = paste(overlap, collapse = ";"),
                              stringsAsFactors = FALSE)
    }
  }
  if (!length(rows)) {
    return(data.frame(collection = character(), set = character(), size = integer(),
                      n_synaptic_anchor_overlap = integer(),
                      synaptic_anchor_overlap = character(), stringsAsFactors = FALSE))
  }
  out <- do.call(rbind, rows)
  out <- out[order(out$collection, out$set, method = "radix"), , drop = FALSE]
  rownames(out) <- NULL
  out
}

clearance_axis_data <- function(pb_de_microglia, pb_de_substate, symbol_map, geomx_de,
                                bulk_omics_summary, mechanism_gene_sets,
                                alpha = 0.10) {
  stopifnot(is.list(pb_de_microglia), is.list(pb_de_substate), is.list(geomx_de),
            is.list(bulk_omics_summary), is.list(mechanism_gene_sets),
            is.numeric(alpha), length(alpha) == 1L, alpha > 0, alpha < 1)
  decon <- geomx_de$decon_preflight
  stopifnot(is.list(decon), "status" %in% names(decon))
  if (identical(decon$status, "earned")) {
    stop("GeoMx decon preflight is earned; add geomx_decon and geomx_abundance_de targets before clearance_axis",
         call. = FALSE)
  }

  anchors <- clearance_axis_dictionary()
  rna_whole <- clearance_rows_from_rna_top(pb_de_microglia$top, symbol_map, "whole_microglia",
                                           anchors, alpha = alpha)
  rna_sub <- lapply(pb_de_substate$per_substate, function(x) {
    if (!identical(x$status, "fit")) return(.empty_clearance_rows())
    clearance_rows_from_rna_top(x$top, symbol_map, x$substate, anchors, alpha = alpha)
  })
  geomx_rows <- clearance_rows_from_geomx(geomx_de, anchors, alpha = alpha)
  bulk_rows <- clearance_rows_from_bulk_summary(bulk_omics_summary, anchors, alpha = alpha)
  measured <- do.call(rbind, c(list(rna_whole), rna_sub, list(geomx_rows, bulk_rows)))
  rownames(measured) <- NULL
  measured <- measured[order(measured$axis, measured$pair, measured$symbol,
                             match(measured$contrast, mechanism_contrasts()),
                             measured$modality, measured$layer, measured$feature,
                             method = "radix", na.last = TRUE), , drop = FALSE]
  rownames(measured) <- NULL

  coverage <- clearance_axis_coverage(measured, anchors)
  pair_support <- clearance_pair_support(measured, anchors, alpha = alpha)
  modality_support <- clearance_axis_modality_support(measured)
  syn_sets <- synaptic_gene_set_rows(mechanism_gene_sets, anchors)
  earned_pairs <- unique(pair_support$pair[pair_support$status == "earned"])
  verdict_status <- if (length(earned_pairs)) "earned" else "not_earned"
  out <- list(
    dictionary = anchors,
    measured = measured,
    coverage = coverage,
    pair_support = pair_support,
    modality_support = modality_support,
    synaptic_gene_sets = syn_sets,
    spatial_decon = list(status = decon$status,
                         action = "skipped",
                         reasons = decon$reasons,
                         background = decon$background,
                         nuclei = decon$nuclei,
                         reference = decon$reference),
    verdict = list(status = verdict_status,
                   ccc_called = FALSE,
                   earned_pairs = earned_pairs,
                   rule = "CCC-lite only; no communication claim unless pair support is earned by measured modalities"),
    provenance = list(alpha = alpha,
                      modalities = sort(unique(measured$modality), method = "radix"),
                      contrasts = mechanism_contrasts(),
                      decon_policy = "SpatialDecon targets are added only when S1 preflight status is earned",
                      bulk_context = bulk_omics_summary$provenance$bulk_context %||%
                        "24M bulk hippocampus; not microglia-sorted")
  )
  stopifnot(nrow(out$coverage) == nrow(anchors),
            all(anchors$symbol %in% out$coverage$symbol),
            nrow(out$pair_support) == length(na.omit(unique(anchors$pair))) *
              length(mechanism_contrasts()),
            out$verdict$status %in% c("earned", "not_earned"),
            identical(out$verdict$ccc_called, FALSE),
            decon$status %in% c("defer", "blocked"),
            nrow(out$synaptic_gene_sets) >= 1L)
  out
}
