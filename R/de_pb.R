# Pseudobulk + bulk DE machinery + the P1 microglia DE orchestration. limma-voom (quality-
# weighted) for aggregated snRNAseq pseudobulk counts, limma-trend for log-intensities
# (proteomics, phospho), a stageR family-screen, and the per-population / per-substate target
# builders. Single-cell DE was DROPPED at the P1 gate -- pseudobulk is the sole condition/
# population inference (Squair 2021 / Murphy-Skene 2022: cell-level DE is pseudoreplicated ->
# FDR-inflated); edgeR-QLF / DESeq2 / dream stay deferred until a phase needs them (KISS).
# Designs + contrasts come from R/design.R (the 5 canonical contrast names flow through
# unchanged). All non-base calls are namespace-qualified -- targets attaches only `quarto`, so
# nothing else is on the search path.

# Aggregate single-cell counts into per-sample pseudobulk: sum raw counts across the cells of
# each group_col level. Returns a features x samples numeric matrix, columns in sorted level
# order (deterministic). Matrix::rowSums handles the sparse dgCMatrix counts; vapply forces a
# matrix result for any (>=1 feature, >=1 group) shape.
pseudobulk_counts <- function(seurat_obj, group_col, assay = "RNA", layer = "counts", cells = NULL) {
  counts <- SeuratObject::GetAssayData(seurat_obj, assay = assay, layer = layer)
  groups <- as.character(seurat_obj@meta.data[[group_col]])
  # groups is positional; Seurat keeps meta.data rows cell-aligned to assay columns. Assert it so a
  # violated invariant fails loud instead of silently summing the wrong cells into each group.
  stopifnot(identical(rownames(seurat_obj@meta.data), colnames(counts)),
            length(groups) == ncol(counts), !anyNA(groups))
  if (!is.null(cells)) {   # restrict to a cell subset (e.g. one substate) BEFORE aggregating
    stopifnot(!anyDuplicated(cells), all(cells %in% colnames(counts)))   # bad/dup cell -> fail loud
    sel <- colnames(counts) %in% cells
    counts <- counts[, sel, drop = FALSE]
    groups <- groups[sel]
  }
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
build_pseudobulk <- function(seurat_obj, sample_col, covariate_cols, assay = "RNA", cells = NULL) {
  cnt  <- pseudobulk_counts(seurat_obj, sample_col, assay = assay, layer = "counts", cells = cells)
  md   <- seurat_obj@meta.data[, c(sample_col, covariate_cols), drop = FALSE]
  if (!is.null(cells)) md <- md[rownames(md) %in% cells, , drop = FALSE]   # match the count subset
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
# (robust eBayes). quality_weights (default TRUE) uses voomWithQualityWeights -> sample-level
# weights down-weight empirically-noisy pseudobulk units (the SOTA P1 default; correlates with,
# does not equal, low cell count). Returns the DGEList, voom object, fit, per-contrast topTables
# (sort.by = "none" -> stable feature order; confint = TRUE -> CI.L/CI.R columns for effect-size
# reporting), and the kept-feature count.
fit_limma_voom <- function(counts, design, contrasts, min_count = 5, quality_weights = TRUE) {
  stopifnot(identical(colnames(counts), rownames(design)),      # limma fits by position, not by name
            identical(rownames(contrasts), colnames(design)),   # contrasts.fit only warns on a mismatch
            qr(design)$rank == ncol(design))                    # full rank -> contrasts estimable
  dge  <- edgeR::DGEList(counts = counts)
  keep <- edgeR::filterByExpr(dge, design = design, min.count = min_count)
  dge  <- dge[keep, , keep.lib.sizes = FALSE]
  dge  <- edgeR::normLibSizes(dge, method = "TMM")   # calcNormFactors renamed in edgeR 4.x
  v    <- if (quality_weights) limma::voomWithQualityWeights(dge, design = design, plot = FALSE)
          else                 limma::voom(dge, design = design, plot = FALSE)
  fit  <- limma::contrasts.fit(limma::lmFit(v, design = design), contrasts)
  fit  <- limma::eBayes(fit, robust = TRUE)
  top  <- lapply(colnames(contrasts), function(cn)
    tibble::rownames_to_column(
      limma::topTable(fit, coef = cn, number = Inf, sort.by = "none", confint = TRUE), "gene"))
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

# ---- P1-S4: pseudobulk DE orchestration (voomWithQualityWeights + stageR family screen) ----

# stageR two-stage family test across the 5-contrast set. SCREEN on the moderated omnibus F
# (H0: every contrast == 0, i.e. no genotype effect on the gene) with BH at the target OFDR;
# CONFIRM per-contrast within each screened gene by stageR's modified post-screen Holm -- FWER-valid
# under arbitrary dependence (safe for the rank-deficient contrast family), and NOT plain p.adjust
# Holm: stageR folds the OFDR screen scaling into the per-gene step. Returns the raw screen p-values,
# the stage-wise adjusted matrix (genes x {padjScreen, the 5 contrasts}; confirmation NA where the
# gene fails the screen -> a contrast is significant at OFDR alpha iff its column <= alpha), and
# the screened-gene count. stageR is on-snapshot -> a locked, reproducible inference layer.
stage_wise_test <- function(fit, alpha = 0.05) {
  stopifnot(!is.null(fit$F.p.value), !is.null(fit$p.value),
            nrow(fit$p.value) == length(fit$F.p.value))
  pScreen <- fit$F.p.value; names(pScreen) <- rownames(fit)
  pConf   <- as.matrix(fit$p.value); rownames(pConf) <- rownames(fit)
  obj  <- stageR::stageR(pScreen = pScreen, pConfirmation = pConf, pScreenAdjusted = FALSE)
  obj  <- stageR::stageWiseAdjustment(obj, method = "holm", alpha = alpha)
  # getAdjustedPValues restates the fixed-OFDR caveat via message() (informational, deterministic,
  # text has no ^Warning anchor) -> muffle so heavy fresh builds keep a clean log
  padj <- suppressMessages(stageR::getAdjustedPValues(obj, onlySignificantGenes = FALSE, order = FALSE))
  list(alpha = alpha, screen_p = pScreen, stage_padj = padj,
       n_screened = sum(padj[, "padjScreen"] <= alpha, na.rm = TRUE))
}

# Honest power statement for the under-powered diff-of-differences `interaction` contrast: the
# NOMINAL minimum detectable log2FC at `power` for the median gene = median posterior SE x the
# two-sided-t multiplier. This is per-test nominal-t power -- NOT stageR/OFDR/BH discovery power,
# and NOT gene-specific. Read it as "the median gene needs an effect >= mde to be detectable at
# this alpha/power": context for the threshold count, never a bare "0 genes" null (absence of
# evidence != evidence of absence). ~9 residual df + eBayes shrinkage -> df.total is the reference.
interaction_power <- function(fit, contrast = "interaction", power = 0.80, alpha = 0.05) {
  stopifnot(contrast %in% colnames(fit$coefficients))
  se     <- sqrt(fit$s2.post) * fit$stdev.unscaled[, contrast]
  df     <- stats::median(fit$df.total, na.rm = TRUE)
  mult   <- stats::qt(1 - alpha / 2, df) + stats::qt(power, df)   # two-sided alpha, one-sided power
  med_se <- stats::median(se, na.rm = TRUE)
  list(contrast = contrast, df = df, median_se = med_se,
       mde = mult * med_se, power = power, alpha = alpha)
}

# Single-population pseudobulk DE: aggregate counts by sample_col -> factorial design + the 5
# contrasts -> voomWQW + robust eBayes -> per-contrast topTables (with CI) -> stageR family
# screen -> interaction power. `cells` (optional) restricts to a cell subset (per-substate)
# before aggregation. The caller guarantees an estimable design (all genotype levels present).
de_pseudobulk <- function(seurat, sample_col = "genotype_batch",
                          covariate_cols = c("genotype", "batch"),
                          assay = "RNA", min_count = 5, alpha = 0.05, cells = NULL) {
  pb <- build_pseudobulk(seurat, sample_col, covariate_cols, assay = assay, cells = cells)
  fd <- factorial_design(pb$meta)
  vf <- fit_limma_voom(pb$counts, fd$design, fd$contrasts, min_count = min_count)
  list(n_samples = ncol(pb$counts), lib_size = colSums(pb$counts), kept = vf$kept,
       top = vf$top, stageR = stage_wise_test(vf$fit, alpha = alpha),
       interaction = interaction_power(vf$fit), thresholds = list(fdr = alpha, lfc = 0.5))
}

# Report (never gate) the amyloid -> DAM activation direction: among DAM-signature genes present
# in the fit, the fraction up-regulated + mean logFC + sig-up count in the two amyloid contrasts.
# v1 prior = DAM up with amyloid -> a concordance check, reported whatever the rebuild shows.
dam_direction <- function(top, symbol_map, markers = canonical_microglia_markers$DAM,
                          contrasts = c("nlgf_in_maptki", "nlgf_in_p301s"),
                          lfc = 0.5, fdr = 0.05) {
  ens   <- unique(symbol_map$ensembl[match(markers, symbol_map$symbol)])   # first hit/symbol, dedup
  ens   <- ens[!is.na(ens)]
  n_req <- length(markers); n_mapped <- length(ens)   # surface map attrition, not just the in-fit count
  lapply(stats::setNames(contrasts, contrasts), function(cn) {
    sub <- top[[cn]][top[[cn]]$gene %in% ens, , drop = FALSE]
    n   <- nrow(sub)                                   # NA (not NaN) when no marker survives the fit
    list(n_markers_requested = n_req, n_markers_mapped = n_mapped, n_markers_in_fit = n,
         frac_up    = if (n) mean(sub$logFC > 0) else NA_real_,
         mean_logFC = if (n) mean(sub$logFC)     else NA_real_,
         n_sig_up   = sum(sub$logFC > lfc & sub$adj.P.Val < fdr))
  })
}

# Guarantee the replicate units form a COMPLETE crossing of the covariates (the balanced
# genotype x batch design). A unit absent from the object would silently shrink the design and fit
# on fewer units -> assert n_units == prod(covariate levels) so a broken object fails loud, not quiet.
assert_complete_crossing <- function(meta, sample_col, covariate_cols = c("genotype", "batch")) {
  n_units <- length(unique(as.character(meta[[sample_col]])))
  n_cross <- prod(vapply(covariate_cols,
                         function(cc) length(unique(as.character(meta[[cc]]))), integer(1)))
  stopifnot(n_units == n_cross)
}

# S4 target: whole-microglia pseudobulk DE (the headline amyloid -> DAM activation programme) +
# the DAM-direction concordance vs the v1 prior. RAW RNA counts -> pseudobulk by genotype_batch
# (the 16 replicate units).
run_pb_de_microglia <- function(seurat, symbol_map, assay = "RNA", min_count = 5, alpha = 0.05) {
  assert_complete_crossing(seurat@meta.data, "genotype_batch")   # all 16 units present, else fail loud
  res <- de_pseudobulk(seurat, assay = assay, min_count = min_count, alpha = alpha)
  res$level          <- "whole_microglia"
  res$n_cells        <- ncol(seurat)
  res$dam_concordance <- dam_direction(res$top, symbol_map)
  res
}

# S4 target: per-substate pseudobulk DE under a PRE-DECLARED min-cell floor. A substate is FIT
# only if EVERY genotype_batch unit carries >= min_cells of it (full estimable factorial design,
# thin units down-weighted by WQW); otherwise SKIPPED -> descriptive-only. The substate x unit
# cell-count table is ALWAYS stored (report the dropout asymmetry, never hide it). Real argmax
# substates only (ambiguous / unassigned are not biological states).
run_pb_de_substate <- function(seurat, substate_col = "microglia_substate",
                               sample_col = "genotype_batch",
                               substates = microglia_substate_levels,
                               min_cells = 10, assay = "RNA", min_count = 5, alpha = 0.05) {
  assert_complete_crossing(seurat@meta.data, sample_col)   # all 16 units present, else fail loud
  sub  <- as.character(seurat@meta.data[[substate_col]])
  unit <- droplevels(factor(as.character(seurat@meta.data[[sample_col]])))
  ct   <- table(substate = factor(sub, levels = substates), unit = unit)   # substate x unit
  per  <- lapply(stats::setNames(substates, substates), function(s) {
    u <- ct[s, ]
    if (any(u < min_cells)) {   # every unit must clear the floor (0-cell units fail it too)
      return(list(status = "skipped", substate = s, n_cells = sum(u), units = u,
                  reason = sprintf("min cells/unit = %d < floor %d (units >= floor: %d/%d)",
                                   min(u), min_cells, sum(u >= min_cells), length(u))))
    }
    cells <- colnames(seurat)[sub == s]
    res <- de_pseudobulk(seurat, sample_col = sample_col, assay = assay,
                         min_count = min_count, alpha = alpha, cells = cells)
    c(list(status = "fit", substate = s, n_cells = length(cells), units = u), res)
  })
  list(min_cells = min_cells, substates = substates, cell_counts = ct, per_substate = per)
}
