# Pseudobulk / bulk DE methods (limma-voom family + edgeR + DESeq2 + dream)
# plus the small matrix-prep helpers used by the proteomics/phospho path.
# Single-cell DE machinery (NEBULA, glmmTMB) is in de_sc.R.
#
# Each fitting helper returns a list of the form
#   list(top = named-list of 5 contrast tibbles, fit = method-specific raw
#        object, n_genes, n_samples_or_cells, design)
# so downstream consumers can iterate over methods uniformly. Contrasts
# follow the project-wide naming: tau_alone, nlgf_in_maptki, nlgf_in_p301s,
# tau_in_nlgf, interaction.

# Pseudobulk: aggregate counts of a Seurat object by a metadata column.
pseudobulk_counts <- function(seurat_obj, group_col, assay = "RNA", slot = "counts") {
  counts <- GetAssayData(seurat_obj, assay = assay, layer = slot)
  groups <- as.character(seurat_obj@meta.data[[group_col]])
  stopifnot(length(groups) == ncol(counts))
  uniq <- sort(unique(groups))
  idx  <- split(seq_along(groups), groups)
  out  <- sapply(uniq, function(g) {
    rowSums(counts[, idx[[g]], drop = FALSE])
  })
  storage.mode(out) <- "integer"
  out
}

# Pseudobulk: build (counts matrix, sample metadata) for limma/edgeR.
build_pseudobulk <- function(seurat_obj, sample_col, covariate_cols, assay = "RNA") {
  cnt  <- pseudobulk_counts(seurat_obj, sample_col, assay = assay, slot = "counts")
  meta <- seurat_obj@meta.data[, c(sample_col, covariate_cols), drop = FALSE] |>
    distinct(.data[[sample_col]], .keep_all = TRUE)
  rownames(meta) <- meta[[sample_col]]
  meta <- meta[colnames(cnt), , drop = FALSE]
  list(counts = cnt, meta = meta)
}

# Standard limma-voom fit returning fits and topTables per contrast.
fit_limma_voom <- function(counts, group, design, contrasts, min_count = 5, weights = NULL) {
  dge <- edgeR::DGEList(counts = counts, group = group)
  keep <- edgeR::filterByExpr(dge, design = design, min.count = min_count)
  dge  <- dge[keep, , keep.lib.sizes = FALSE]
  dge  <- edgeR::calcNormFactors(dge, method = "TMM")
  v    <- limma::voom(dge, design = design, plot = FALSE)
  fit  <- limma::lmFit(v, design = design, weights = weights)
  fit2 <- limma::contrasts.fit(fit, contrasts)
  fit2 <- limma::eBayes(fit2, robust = TRUE)
  tts <- lapply(colnames(contrasts), function(cn) {
    limma::topTable(fit2, coef = cn, number = Inf, sort.by = "none") |>
      tibble::rownames_to_column("gene")
  })
  names(tts) <- colnames(contrasts)
  list(dge = dge, voom = v, fit = fit2, top = tts, kept = sum(keep))
}

# Standard limma fit on already-log-transformed matrix (proteomics, phospho).
fit_limma_log <- function(mat, group, design, contrasts) {
  fit <- limma::lmFit(mat, design = design)
  fit <- limma::contrasts.fit(fit, contrasts)
  fit <- limma::eBayes(fit, robust = TRUE, trend = TRUE)
  tts <- lapply(colnames(contrasts), function(cn) {
    limma::topTable(fit, coef = cn, number = Inf, sort.by = "none") |>
      tibble::rownames_to_column("feature")
  })
  names(tts) <- colnames(contrasts)
  list(fit = fit, top = tts)
}

# Simple median normalisation (sample-wise) of a log-intensity matrix.
median_normalise <- function(mat) {
  med <- apply(mat, 2, median, na.rm = TRUE)
  global <- median(med, na.rm = TRUE)
  sweep(mat, 2, med - global, "-")
}

# Filter a matrix: features with >= min_present samples per group in >= min_groups groups.
prevalence_filter <- function(mat, group, min_present = 3, min_groups = 2) {
  stopifnot(!is.null(rownames(mat)))
  mat <- mat[!is.na(rownames(mat)) & rownames(mat) != "", , drop = FALSE]
  group <- as.character(group)
  stopifnot(length(group) == ncol(mat))
  present <- !is.na(mat)
  group_levels <- unique(group)
  per_group_counts <- vapply(group_levels, function(g) {
    rowSums(present[, group == g, drop = FALSE])
  }, numeric(nrow(mat)))
  if (length(group_levels) == 1) per_group_counts <- matrix(per_group_counts, ncol = 1)
  ok <- rowSums(per_group_counts >= min_present) >= min_groups
  mat[ok, , drop = FALSE]
}

# edgeR Quasi-likelihood F-test on pseudobulks (muscat::pbDS method="edgeR"
# is a wrapper around the same `glmQLFit`/`glmQLFTest` calls; we go direct).
fit_edger_qlf <- function(counts, design, contrasts, min_count = 10,
                          symbol_map = NULL) {
  dge  <- edgeR::DGEList(counts = counts)
  keep <- edgeR::filterByExpr(dge, design = design, min.count = min_count)
  dge  <- dge[keep, , keep.lib.sizes = FALSE]
  dge  <- edgeR::calcNormFactors(dge, method = "TMM")
  dge  <- edgeR::estimateDisp(dge, design = design, robust = TRUE)
  fit  <- edgeR::glmQLFit(dge, design = design, robust = TRUE)
  tts  <- lapply(colnames(contrasts), function(cn) {
    qlf <- edgeR::glmQLFTest(fit, contrast = contrasts[, cn])
    tt  <- edgeR::topTags(qlf, n = Inf, sort.by = "none")$table
    tt  <- tibble::rownames_to_column(tt, "gene")
    tt$P.Value   <- tt$PValue
    tt$adj.P.Val <- tt$FDR
    tt$t         <- tt$`F`
    tt$se        <- NA_real_
    tt$symbol    <- if (!is.null(symbol_map)) {
      symbol_map$symbol[match(tt$gene, symbol_map$ensembl)]
    } else NA_character_
    tibble::as_tibble(tt[, c("gene", "symbol", "logFC", "se", "t",
                             "P.Value", "adj.P.Val")])
  })
  names(tts) <- colnames(contrasts)
  list(fit = fit, top = tts, n_genes = sum(keep),
       n_samples = ncol(counts), design = design)
}

# DESeq2 Wald test on pseudobulks (muscat::pbDS method="DESeq2" equivalent).
# DESeq2 estimates its own size factors and dispersion shrinkage; we hand it
# the same factorial design so contrasts are aligned with the other methods.
fit_deseq2_pb <- function(counts, meta_df, design, contrasts,
                          min_count = 10, symbol_map = NULL) {
  stopifnot(requireNamespace("DESeq2", quietly = TRUE))
  dge_for_filt <- edgeR::DGEList(counts = counts)
  keep <- edgeR::filterByExpr(dge_for_filt, design = design,
                              min.count = min_count)
  counts <- counts[keep, , drop = FALSE]

  storage.mode(counts) <- "integer"
  dds <- DESeq2::DESeqDataSetFromMatrix(counts, colData = meta_df,
                                         design = design)
  dds <- DESeq2::DESeq(dds, fitType = "parametric", quiet = TRUE,
                       parallel = FALSE)

  tts <- lapply(colnames(contrasts), function(cn) {
    res <- DESeq2::results(dds, contrast = contrasts[, cn],
                           independentFiltering = TRUE, alpha = 0.05,
                           cooksCutoff = TRUE)
    tt <- as.data.frame(res) |> tibble::rownames_to_column("gene")
    tt$logFC     <- tt$log2FoldChange
    tt$se        <- tt$lfcSE
    tt$t         <- tt$stat
    tt$P.Value   <- tt$pvalue
    tt$adj.P.Val <- tt$padj
    tt$symbol    <- if (!is.null(symbol_map)) {
      symbol_map$symbol[match(tt$gene, symbol_map$ensembl)]
    } else NA_character_
    tibble::as_tibble(tt[, c("gene", "symbol", "logFC", "se", "t",
                             "P.Value", "adj.P.Val")])
  })
  names(tts) <- colnames(contrasts)
  list(fit = dds, top = tts, n_genes = sum(keep),
       n_samples = ncol(counts), design = design)
}

# dreamlet-style fit using variancePartition::dream directly. dreamlet is
# itself a thin wrapper around dream that handles multi-cluster aggregation;
# we already have single-cluster pseudobulks so calling dream directly is
# cleaner and avoids re-aggregating from cell-level SCE.
#
# Mixed-effects formula: ~ tau + nlgf + tau_nlgf + (1 | batch). batch enters
# as a random intercept (4 levels) instead of the fixed-effect coding used
# by limma-voom, giving a different empirical-Bayes shrinkage profile.
fit_dream_pb <- function(counts, meta_df, contrasts,
                         formula = ~ tau + nlgf + tau_nlgf + (1 | batch),
                         min_count = 10, symbol_map = NULL,
                         BPPARAM = NULL) {
  if (is.null(BPPARAM)) {
    nc <- max(1L, parallel::detectCores() - 1L)
    BPPARAM <- if (nc > 1L) BiocParallel::MulticoreParam(workers = nc)
               else BiocParallel::SerialParam()
  }
  stopifnot(requireNamespace("variancePartition", quietly = TRUE))
  # Pre-filter using a fixed-effect projection of the same design so the
  # gene set matches the other pseudobulk methods.
  design_fix <- model.matrix(~ tau + nlgf + tau_nlgf + batch, data = meta_df)
  dge  <- edgeR::DGEList(counts = counts)
  keep <- edgeR::filterByExpr(dge, design = design_fix, min.count = min_count)
  dge  <- dge[keep, , keep.lib.sizes = FALSE]
  dge  <- edgeR::calcNormFactors(dge, method = "TMM")

  v <- variancePartition::voomWithDreamWeights(dge, formula, meta_df,
                                               BPPARAM = BPPARAM)

  # Build dream contrast objects from the fixed-effect coefficient names.
  L <- variancePartition::makeContrastsDream(
    formula, meta_df,
    contrasts = c(
      tau_alone      = "tau",
      nlgf_in_maptki = "nlgf",
      interaction    = "tau_nlgf",
      nlgf_in_p301s  = "nlgf + tau_nlgf",
      tau_in_nlgf    = "tau + tau_nlgf"
    )
  )

  fit <- variancePartition::dream(v, formula, meta_df, L = L,
                                  BPPARAM = BPPARAM)
  fit <- variancePartition::eBayes(fit)

  tts <- lapply(colnames(L), function(cn) {
    tt <- variancePartition::topTable(fit, coef = cn, number = Inf,
                                      sort.by = "none")
    tt <- tibble::rownames_to_column(tt, "gene")
    tt$se <- tt$logFC / tt$t
    tt$symbol <- if (!is.null(symbol_map)) {
      symbol_map$symbol[match(tt$gene, symbol_map$ensembl)]
    } else NA_character_
    tibble::as_tibble(tt[, c("gene", "symbol", "logFC", "se", "t",
                             "P.Value", "adj.P.Val")])
  })
  names(tts) <- colnames(L)
  list(fit = fit, top = tts, n_genes = sum(keep),
       n_samples = ncol(counts), design = L)
}
