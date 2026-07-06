# Differential expression for the three non-snRNAseq modalities (GeoMx spatial, 24M bulk
# proteome, 24M bulk phosphosite), restored lean from the pre-teardown P4 crossmodality /
# mechanism modules. SCOPE = the PRIMARY per-contrast topTables only: each producer returns
# `$top[[contrast]]` (or `$primary$top[[contrast]]` for GeoMx) with a `logFC` column keyed by
# the 5 canonical contrasts. The amyloid-response logFC scatter (R/figures.R ->
# modality_logfc_scatter_data) reads `nlgf_in_maptki` (y) vs `nlgf_in_p301s` (x) from these.
# The torn-down auxiliary arms (GeoMx unblocked / bio-unit-collapsed sensitivities +
# SpatialDecon decon-preflight; proteome / phospho additive run-index sensitivity) are
# INTENTIONALLY not restored -- they served the deconvolution / bulk-run-order caveats, not the
# per-feature effect sizes this figure needs. Shared machinery reused from HEAD: fit_limma_log /
# median_normalise / prevalence_filter (R/de_pb.R), factorial_design / make_contrast_matrix
# (R/design.R), match_intensity_columns / load_geomx / read_spectronaut_tsv / proteomics_sample_meta
# (R/io.R). All calls namespace-qualified so the file sources cleanly into any session.

# --- GeoMx spatial DE (RNA counts; slide fixed effect + bio-unit duplicateCorrelation) -------

geomx_required_meta_cols <- function() {
  c("genotype", "slide_rep", "bio_rep", "roi", "SampleID",
    "ROI Coordinate X", "ROI Coordinate Y",
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
