# --------------------------------------------------------------------
# GeoMx spatial deconvolution (plan arc L). Pure functions that turn the
# snRNAseq reference + the 91-ROI GeoMx-WTA object into per-ROI cell-type
# / microglial-substate ABUNDANCE estimates (SpatialDecon; Danaher 2022),
# push those abundances through the project's locked 5-contrast 2x2
# factorial (R/design.R + R/de_pb.R::fit_limma_log), and run a per-slide
# spatial-autocorrelation check (hand-rolled Moran's I + permutation).
# The I/O driver is scripts/build_spatial_deconvolution.R (step L3); the
# display chapter is rmd/21 (L4). NO function here writes to disk.
#
# Why a new readout: every prior arc (D pathway .. K SCENIC) reads
# EXPRESSION or ACTIVITY; none reads HOW MUCH of each cell type is present
# in tissue. The GeoMx is the only location-aware modality. Deconvolving
# its UNSORTED ROIs gives a tissue-composition estimate the sorting-biased
# snRNAseq cannot. Framed as orthogonal corroboration, NOT proof.
#
# Locked facts wired in (L1, verified at L2 against the live caches):
#  * Reference = seurat_full_processed.rds; rownames are MGI SYMBOLS (not
#    Ensembl) -> intersect directly with the GeoMx symbols, NO mapping.
#  * GeoMx `data` layer = vendor Q3-normalised = counts / q_norm_qFactors
#    (verified: data/counts == 1/qFactor, constant across genes per ROI,
#    qFactor 0.10-7.49). => the negative-probe background MUST be Q3-scaled
#    the SAME way: bg_roi = NegGeoMean_roi / qFactor_roi. This corrects the
#    L1 smoke, which broadcast the RAW NegGeoMean (a scale mismatch up to
#    ~7x); on a Q3-normalised `norm` the background has to live in Q3 space.
#  * `nuclei` carries a -1 missing sentinel in 42/91 ROIs -> absolute-count
#    scaling via spatialdecon's `cell_counts` would NA half the ROIs, so it
#    is NOT used; beta stays on its native (normalize=TRUE) abundance scale.
#  * Two-stage granularity (user gate): stage-1 6-level (5 broad + pooled
#    Microglia = robust total) -> stage-2 sub-deconvolve the 4 substates and
#    anchor within-microglia fractions to the stage-1 total.
#
# Caveats carried downstream (state at EVERY interpretation): the snRNAseq
# reference proportions are sorting/loading-biased; the 4 substates are
# transcriptionally close (collinear profile columns -> unstable substate
# betas -- the build reports the condition number + max column correlation
# so this is visible, never hidden); GeoMx ROIs are whole-tissue geometric
# samples with NO plaque/Iba1 targeting (shifts are regional-bulk, not
# plaque-niche); single-nucleus reference vs GeoMx-WTA platform differ;
# slides MIX genotypes (per-slide Moran's I is confounded by within-slide
# genotype composition -> optional `covar` residualisation).
# --------------------------------------------------------------------

# Build a gene x cell-type expression profile from the snRNAseq reference
# via SpatialDecon::create_profile_matrix. Cells are capped per type
# (default 3000) before the call because create_profile_matrix densifies
# its input internally; the profile is a per-type MEAN, so subsampling is
# ~lossless for the profile while bounding a multi-GB densification at the
# 286k-cell reference. NO symbol mapping -- `sc` rownames are already MGI
# symbols (L1 fact); the profile gene set is intersect(rownames(sc),
# geomx_symbols) (mtx is pre-subset to that intersection to halve memory).
#
# Returns list(profile = gene x cell-type matrix, qc = list(...)). qc
# reports gene coverage, cells used per type, the profile condition number
# (kappa; large => collinear columns / unstable betas) and the maximum
# off-diagonal Pearson correlation among the log1p cell-type columns (the
# substate-collinearity flag from the guardrails) plus the full cor matrix.
build_reference_profile <- function(sc, cell_type_col, geomx_symbols,
                                    cap = 3000L, seed = 1L,
                                    minCellNum = 15L, minGenes = 100L) {
  stopifnot(requireNamespace("SpatialDecon", quietly = TRUE),
            cell_type_col %in% colnames(sc@meta.data))
  ct <- as.character(sc@meta.data[[cell_type_col]])
  # Cap cells per type (reproducible) to bound the internal densification.
  set.seed(seed)
  idx_by <- split(seq_along(ct), ct)
  keep <- sort(unlist(lapply(idx_by, function(ix)
    if (length(ix) > cap) sample(ix, cap) else ix), use.names = FALSE))
  common <- intersect(rownames(sc), geomx_symbols)
  cnt <- GetAssayData(sc, assay = "RNA", layer = "counts")[common, keep, drop = FALSE]
  annots <- data.frame(CellID = colnames(cnt), ct = ct[keep],
                       stringsAsFactors = FALSE)
  names(annots)[2] <- cell_type_col
  prof <- suppressMessages(SpatialDecon::create_profile_matrix(
    mtx = cnt, cellAnnots = annots,
    cellTypeCol = cell_type_col, cellNameCol = "CellID",
    geneList = rownames(cnt), normalize = TRUE,
    minCellNum = as.integer(minCellNum), minGenes = as.integer(minGenes),
    discardCellTypes = FALSE, outDir = NULL))
  prof <- as.matrix(prof)
  cc <- suppressWarnings(stats::cor(log1p(prof)))
  diag(cc) <- NA_real_
  qc <- list(
    n_common_genes  = length(common),
    n_profile_genes = nrow(prof),
    n_cells_used    = length(keep),
    cells_per_type  = table(ct[keep]),
    cell_types      = colnames(prof),
    condition_number     = kappa(prof, exact = FALSE),
    max_abs_celltype_cor = max(abs(cc), na.rm = TRUE),
    celltype_cor    = cc)
  list(profile = prof, qc = qc)
}

# Derive the per-ROI negative-probe background as a gene x ROI matrix in the
# SAME scale as the `data` (Q3) layer fed to spatialdecon. The WTA object
# collapses the negative probes to a per-ROI geomean in metadata (no
# neg-probe ROWS for SpatialDecon::derive_GeoMx_background to read), so the
# background is built manually: take NegGeoMean per ROI, Q3-scale it by the
# SAME factor that maps counts -> data (data = counts / q_norm_qFactors,
# verified), then broadcast across all genes. If `qfactor_col` is absent the
# raw geomean is broadcast (with a warning) -- correct only if `norm` is on
# the raw-count scale.
derive_geomx_background <- function(geomx,
                                    neg_col = "NegGeoMean_Mm_R_NGS_WTA_v1.0",
                                    qfactor_col = "q_norm_qFactors",
                                    assay = "RNA", layer = "data") {
  norm <- GetAssayData(geomx, assay = assay, layer = layer)
  md <- geomx@meta.data
  stopifnot(neg_col %in% colnames(md))
  negmean <- as.numeric(md[[neg_col]])
  if (!is.null(qfactor_col) && qfactor_col %in% colnames(md)) {
    bg_roi <- negmean / as.numeric(md[[qfactor_col]])   # Q3-scale to `data`
  } else {
    warning("qfactor_col absent; broadcasting RAW NegGeoMean (raw-scale only)")
    bg_roi <- negmean
  }
  bg <- matrix(1, nrow(norm), ncol(norm),
               dimnames = list(rownames(norm), colnames(norm)))
  sweep(bg, 2, bg_roi, "*")
}

# Run SpatialDecon on the Q3 `data` layer against a reference profile.
# align_genes=TRUE reconciles the gene sets of norm / bg / profile. Returns
# list(beta = cell-type x ROI abundance, prop = column-normalised beta
# [self-computed; spatialdecon's prop_of_all returned NaN here, L1 to-do],
# fit_qc = per-ROI + overall residual RMSE on the log-fit scale + aligned
# gene count, sigmas). `cell_counts` (nuclei) is deliberately NOT passed
# (42/91 ROIs carry the -1 sentinel), so beta keeps its native scale.
run_spatialdecon <- function(geomx, profile, background,
                             assay = "RNA", layer = "data") {
  stopifnot(requireNamespace("SpatialDecon", quietly = TRUE))
  norm <- as.matrix(GetAssayData(geomx, assay = assay, layer = layer))
  res <- SpatialDecon::spatialdecon(norm = norm, bg = background,
                                    X = profile, align_genes = TRUE)
  beta <- res$beta
  cs <- colSums(beta)
  prop <- sweep(beta, 2, ifelse(cs > 0, cs, NA_real_), "/")
  resid <- res$resids %||% res$resid
  per_roi_rmse <- if (!is.null(resid)) sqrt(colMeans(resid^2, na.rm = TRUE)) else NA_real_
  list(beta = beta, prop = prop,
       fit_qc = list(
         per_roi_rmse    = per_roi_rmse,
         resid_rmse      = if (!is.null(resid)) sqrt(mean(resid^2, na.rm = TRUE)) else NA_real_,
         n_genes_aligned = length(intersect(rownames(norm), rownames(profile)))),
       sigmas = res$sigmas)
}

# Assemble the two-stage abundance layers. beta6 = 6-level deconvolution
# (5 broad + pooled `micro_label`); beta9 = 9-level (5 broad + the 4
# `substate_prefix` substates). Stage-1 keeps beta6 verbatim (the robust
# total-microglia + broad layer). Stage-2 takes the WITHIN-microglia
# substate fractions from beta9 and anchors them to the stage-1 total, so
# substate abundance is conditional on the stable total (the locked
# two-stage mechanic). On RESOLVED ROIs colSums(stage2) == microglia_total
# (fractions sum to 1) -> consistency_spearman == 1; reported as a check.
# An ROI whose 9-level substate betas are all zero (no resolvable microglia
# in the substate fit) is UNRESOLVED: its fraction is 0/0, so the substate
# columns are set NA (honest -- a zero cannot be split) and counted in
# n_unresolved; the consistency check then uses only complete columns.
combine_two_stage <- function(beta6, beta9, micro_label = "Microglia",
                              substate_prefix = "Microglia_") {
  stopifnot(micro_label %in% rownames(beta6),
            all(colnames(beta6) == colnames(beta9)))
  subs <- grep(paste0("^", substate_prefix), rownames(beta9), value = TRUE)
  stopifnot(length(subs) >= 2)
  micro_total <- beta6[micro_label, ]
  sub9 <- beta9[subs, , drop = FALSE]
  cs <- colSums(sub9)
  frac <- sweep(sub9, 2, ifelse(cs > 0, cs, NA_real_), "/")
  stage2 <- sweep(frac, 2, micro_total, "*")
  consistency <- suppressWarnings(stats::cor(colSums(stage2), micro_total,
                                  use = "complete.obs", method = "spearman"))
  list(stage1 = beta6, stage2 = stage2,
       microglia_total = micro_total, substate_fractions = frac,
       n_unresolved = sum(cs == 0), consistency_spearman = consistency)
}

# Fit the locked 2x2 factorial + slide-block design to a cell-type x ROI
# abundance matrix on the log scale and extract all 5 contrasts per cell
# type via limma-trend. The design mirrors the GeoMx gene DE (rmd/04,
# genotype + slide): factorial_design encodes the SAME genotype effects as
# `~ 0 + genotype + slide` + make_contrast_matrix but in the 2x2 + batch
# parameterisation, so the interaction is byte-comparable to every other
# modality. `slide_col` is passed as factorial_design's batch_col. Default
# offset = 1 is a one-cell pseudocount (normalize=TRUE betas are on a cell-
# abundance scale); record the realised value. Run separately on the
# stage-1 (6 cell types) and stage-2 (4 substates) layers.
#
# Returns tibble(cell_type, contrast, logFC, t, P.Value, adj.P.Val, sig,
#                offset); adj.P.Val is limma's BH across cell types within
# each contrast.
fit_abundance_contrasts <- function(abund, meta, offset = 1,
                                    genotype_col = "genotype",
                                    slide_col = "slide name",
                                    padj_cut = 0.10) {
  stopifnot(all(colnames(abund) %in% rownames(meta)))
  meta <- meta[colnames(abund), , drop = FALSE]
  log_abund <- log(abund + offset)
  fd <- factorial_design(meta, genotype_col = genotype_col,
                         batch_col = slide_col, add_batch = TRUE)
  mat <- log_abund[, rownames(fd$design), drop = FALSE]
  fit <- fit_limma_log(mat, group = meta[[genotype_col]],
                       design = fd$design, contrasts = fd$contrasts)
  rows <- lapply(names(fit$top), function(cn) {
    tt <- fit$top[[cn]]
    data.frame(cell_type = tt$feature, contrast = cn,
               logFC = tt$logFC, t = tt$t, P.Value = tt$P.Value,
               adj.P.Val = tt$adj.P.Val, stringsAsFactors = FALSE)
  })
  res <- dplyr::bind_rows(rows)
  res$sig <- res$adj.P.Val < padj_cut
  res$offset <- offset
  tibble::as_tibble(res)
}

# Moran's I for a vector x under a weights matrix W (diag 0). NA if x has
# no variance. Used per (slide, cell_type) by spatial_autocorrelation.
.morans_I <- function(x, W) {
  n <- length(x); xc <- x - mean(x); s2 <- sum(xc^2)
  if (s2 == 0) return(NA_real_)
  (n / sum(W)) * (sum(W * tcrossprod(xc)) / s2)
}

# Per-slide spatial autocorrelation of each cell-type abundance, with a
# permutation test. Slides are handled SEPARATELY (per-slide scan offsets
# differ; within a slide the raw ROI X/Y are comparable). Weights are
# inverse-distance (w_ij = 1/d_ij, diag 0) from the per-slide coords. The
# reported p is one-sided for POSITIVE autocorrelation (the spatial-gradient
# / clustering hypothesis): p = (1 + #{I_perm >= I_obs}) / (n_perm + 1).
# `covar` (optional, e.g. genotype) residualises each abundance vector
# within the slide before the test -- slides mix genotypes, so the raw
# Moran's I is confounded by within-slide genotype composition; pass
# genotype to isolate residual spatial structure.
#
# abund: cell-type x ROI. coords: ROI x 2 (X, Y), rows aligned to columns of
# abund. slide: per-ROI label (length = ncol(abund)). Returns tibble(slide,
# cell_type, n_roi, morans_I, expected_I, p_perm, residualised).
spatial_autocorrelation <- function(abund, coords, slide, covar = NULL,
                                    n_perm = 999L, seed = 1L, min_roi = 5L) {
  coords <- as.matrix(coords)
  stopifnot(ncol(abund) == nrow(coords), ncol(abund) == length(slide))
  slide <- as.character(slide)
  set.seed(seed)
  rows <- list()
  for (s in sort(unique(slide))) {
    ix <- which(slide == s)
    if (length(ix) < min_roi) next
    d <- as.matrix(stats::dist(coords[ix, , drop = FALSE]))
    W <- 1 / d; diag(W) <- 0; W[!is.finite(W)] <- 0
    exp_I <- -1 / (length(ix) - 1)
    cv <- if (!is.null(covar)) covar[ix] else NULL
    for (k in rownames(abund)) {
      x <- as.numeric(abund[k, ix])
      if (!is.null(cv) && length(unique(cv)) > 1)
        x <- as.numeric(residuals(stats::lm(x ~ factor(cv))))
      I_obs <- .morans_I(x, W)
      p <- if (is.na(I_obs)) NA_real_ else {
        I_perm <- replicate(n_perm, .morans_I(sample(x), W))
        (1 + sum(I_perm >= I_obs, na.rm = TRUE)) / (n_perm + 1)
      }
      rows[[length(rows) + 1L]] <- data.frame(
        slide = s, cell_type = k, n_roi = length(ix),
        morans_I = I_obs, expected_I = exp_I, p_perm = p,
        residualised = !is.null(cv), stringsAsFactors = FALSE)
    }
  }
  tibble::as_tibble(dplyr::bind_rows(rows))
}
