# Pseudobulk + bulk DE machinery: limma-voom for aggregated counts (snRNAseq pseudobulk),
# limma-trend for log-intensities (proteomics, phospho), plus the matrix-prep helpers.
# Pure functions consumed by P1+ DE targets; S3 builds + unit-tests the machinery and runs
# no analysis. Single-cell DE (NEBULA/glmmTMB) lives in P1's de_sc.R; edgeR-QLF / DESeq2 /
# dream are deferred until a phase needs them (KISS). Designs + contrasts come from
# R/design.R (the 5 canonical contrast names flow through unchanged). All non-base calls are
# namespace-qualified -- targets attaches only `quarto`, so nothing else is on the search path.

# Aggregate single-cell counts into per-sample pseudobulk: sum raw counts across the cells of
# each group_col level. Returns a features x samples numeric matrix, columns in sorted level
# order (deterministic). Matrix::rowSums handles the sparse dgCMatrix counts; vapply forces a
# matrix result for any (>=1 feature, >=1 group) shape.
pseudobulk_counts <- function(seurat_obj, group_col, assay = "RNA", layer = "counts") {
  counts <- SeuratObject::GetAssayData(seurat_obj, assay = assay, layer = layer)
  groups <- as.character(seurat_obj@meta.data[[group_col]])
  # groups is positional; Seurat keeps meta.data rows cell-aligned to assay columns. Assert it so a
  # violated invariant fails loud instead of silently summing the wrong cells into each group.
  stopifnot(identical(rownames(seurat_obj@meta.data), colnames(counts)),
            length(groups) == ncol(counts), !anyNA(groups))
  uniq <- sort(unique(groups), method = "radix")   # locale-independent -> deterministic column order
  idx  <- split(seq_along(groups), groups)
  out  <- vapply(uniq, function(g) Matrix::rowSums(counts[, idx[[g]], drop = FALSE]),
                 numeric(nrow(counts)))
  # vapply collapses to a vector when nrow(counts)==1 -> re-wrap so a single-feature object still
  # returns a features x samples matrix (matches the docstring guarantee + prevalence_filter).
  matrix(out, nrow = nrow(counts), dimnames = list(rownames(counts), uniq))
}

# Build (counts, aligned meta) for a pseudobulk fit. Sums counts by sample_col, then attaches
# one metadata row per sample (covariate_cols). Fails loud unless each covariate is constant
# within a sample (else the one-row-per-sample pick would silently lose information) and the
# meta rows align 1:1 with the count columns.
build_pseudobulk <- function(seurat_obj, sample_col, covariate_cols, assay = "RNA") {
  cnt  <- pseudobulk_counts(seurat_obj, sample_col, assay = assay, layer = "counts")
  md   <- seurat_obj@meta.data[, c(sample_col, covariate_cols), drop = FALSE]
  stopifnot(!anyNA(md))   # an all-NA covariate reads as "constant" (length(unique(NA))==1); reject it
  samp <- as.character(md[[sample_col]])
  for (cc in covariate_cols) {
    n_unique <- tapply(as.character(md[[cc]]), samp, function(v) length(unique(v)))
    stopifnot(all(n_unique == 1L))
  }
  meta <- md[!duplicated(samp), , drop = FALSE]
  rownames(meta) <- as.character(meta[[sample_col]])
  meta <- meta[colnames(cnt), , drop = FALSE]
  stopifnot(identical(rownames(meta), colnames(cnt)))
  list(counts = cnt, meta = meta)
}

# limma-voom fit for pseudobulk counts: drop low-count features (filterByExpr, min.count),
# TMM-normalise, voom-weight, fit the supplied design, apply the contrast matrix, moderate
# (robust eBayes). Returns the DGEList, voom object, fit, per-contrast topTables (sort.by =
# "none" -> stable feature order), and the kept-feature count.
fit_limma_voom <- function(counts, design, contrasts, min_count = 5) {
  stopifnot(identical(colnames(counts), rownames(design)),      # limma fits by position, not by name
            identical(rownames(contrasts), colnames(design)),   # contrasts.fit only warns on a mismatch
            qr(design)$rank == ncol(design))                    # full rank -> contrasts estimable
  dge  <- edgeR::DGEList(counts = counts)
  keep <- edgeR::filterByExpr(dge, design = design, min.count = min_count)
  dge  <- dge[keep, , keep.lib.sizes = FALSE]
  dge  <- edgeR::normLibSizes(dge, method = "TMM")   # calcNormFactors renamed in edgeR 4.x
  v    <- limma::voom(dge, design = design, plot = FALSE)
  fit  <- limma::contrasts.fit(limma::lmFit(v, design = design), contrasts)
  fit  <- limma::eBayes(fit, robust = TRUE)
  top  <- lapply(colnames(contrasts), function(cn)
    tibble::rownames_to_column(
      limma::topTable(fit, coef = cn, number = Inf, sort.by = "none"), "gene"))
  names(top) <- colnames(contrasts)
  list(dge = dge, voom = v, fit = fit, top = top, kept = sum(keep))
}

# limma-trend fit for an already-log-transformed intensity matrix (proteomics, phospho): no
# voom (continuous, not counts); trend = TRUE models the mean-variance trend on the log scale.
fit_limma_log <- function(mat, design, contrasts) {
  stopifnot(identical(colnames(mat), rownames(design)),         # limma fits by position, not by name
            identical(rownames(contrasts), colnames(design)),   # contrasts.fit only warns on a mismatch
            qr(design)$rank == ncol(design))                    # full rank -> contrasts estimable
  fit <- limma::contrasts.fit(limma::lmFit(mat, design = design), contrasts)
  fit <- limma::eBayes(fit, robust = TRUE, trend = TRUE)
  top <- lapply(colnames(contrasts), function(cn)
    tibble::rownames_to_column(
      limma::topTable(fit, coef = cn, number = Inf, sort.by = "none"), "feature"))
  names(top) <- colnames(contrasts)
  list(fit = fit, top = top)
}

# Sample-wise median normalisation of a log-intensity matrix: shift each column so its median
# matches the global median of the per-column medians. NA-robust.
median_normalise <- function(mat) {
  med    <- apply(mat, 2, stats::median, na.rm = TRUE)
  global <- stats::median(med, na.rm = TRUE)
  sweep(mat, 2, med - global, "-")
}

# Prevalence filter: keep features present (non-NA) in >= min_present samples within each of
# >= min_groups groups. Drops features with missing/empty rownames first. Returns the
# row-subset matrix.
prevalence_filter <- function(mat, group, min_present = 3L, min_groups = 2L) {
  stopifnot(!is.null(rownames(mat)), length(group) == ncol(mat))
  mat   <- mat[!is.na(rownames(mat)) & rownames(mat) != "", , drop = FALSE]
  group <- as.character(group)
  stopifnot(!anyNA(group))   # NA group -> NA column index -> fabricated all-NA rows; fail loud
  present   <- !is.na(mat)
  per_group <- vapply(sort(unique(group), method = "radix"),
                      function(g) rowSums(present[, group == g, drop = FALSE]),
                      numeric(nrow(mat)))
  if (!is.matrix(per_group)) per_group <- matrix(per_group, nrow = nrow(mat))
  keep <- rowSums(per_group >= min_present) >= min_groups
  mat[keep, , drop = FALSE]
}
