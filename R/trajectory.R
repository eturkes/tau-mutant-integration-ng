# --------------------------------------------------------------------
# Microglial activation-trajectory / dynamics layer (plan arc M). Pure
# functions that infer a homeostatic -> DAM activation pseudotime from the
# snRNAseq microglia (Slingshot on the harmony embedding, the 4 `state`s as
# clusters, root = homeostatic), validate the asserted root with a
# root-free potency proxy, summarise pseudotime PER REPLICATE
# (genotype_batch, the locked 16-id unit), push that summary through the
# project's locked 5-contrast 2x2 factorial (R/design.R + R/de_pb.R::
# fit_limma_log), and -- the load-bearing piece -- DECOMPOSE the per-
# replicate shift into a between-state COMPOSITION channel and a within-
# state PROGRESSION channel so a mean-pseudotime move that is really just
# cells reshuffling into the larger DAM cluster (the arc-L composition
# result) is never mistaken for genuine activation advancement. The I/O
# driver is scripts/build_trajectory.R (step M3); the display chapter is
# rmd/22 (M4). NO function here writes to disk.
#
# Why a new readout: every prior arc reads a STATIC quantity -- expression
# DE (02b), pathway/module/TF/kinase/CCC/causal/SCENIC activity (D..K), or
# discrete cell COMPOSITION (L). None reads the ORDERING of microglia along
# an activation continuum, i.e. how far along the homeostatic->DAM program
# each cell sits, including within-cluster advancement that neither DE nor
# discrete composition can see. The interaction
# (NLGF_P301S-P301S)-(NLGF_MAPTKI-MAPTKI) is here re-expressed as a
# PROGRESSION effect.
#
# Locked facts wired in (M1, verified against the live cache):
#  * Substrate = microglia_seurat_processed.rds (26,104 cells); embedding =
#    `harmony` (batch-corrected; NOT raw pca). The caller truncates harmony
#    to n_dims BEFORE calling build_microglia_trajectory -- M1 found the
#    full 30 dims let cell-cycle pull `proliferative` in as a spurious
#    intermediate, so 10-15 dims + treating proliferative as a cycling
#    side-state (omit_clusters) is the documented mitigation. The policy
#    (n_dims, which clusters to omit) lives in the build script, not here:
#    these functions stay un-opinionated and report whatever topology
#    Slingshot returns (guardrail #8: locks fixed before inspection).
#  * `state` (homeostatic/DAM/IFN/proliferative) is re-derived at build
#    time by label_microglia_states (R/microglia.R); NOT cached.
#  * Replication unit = genotype_batch (16 ids); FDR<0.10; contrast set =
#    the 5 canonical. Root = homeostatic, to be VALIDATED by potency here,
#    not asserted. Velocity is OUT (no spliced layers).
#
# Caveats carried downstream (state at EVERY interpretation): pseudotime is
# a single inferred ACTIVATION ORDERING, not chronological/developmental
# time (guardrail #3); the substates are transcriptionally close so a
# per-replicate mean shift is substantially composition-driven -- ALWAYS
# read the decomposition + cross-method concordance, never the raw mean
# alone (guardrails #1, #6); a null/additive interaction is a finding, not
# a failure (guardrail #4); all states / all lineages are reported, absence
# of an IFN or proliferative effect included (guardrail #7).
# --------------------------------------------------------------------

# Build a Slingshot trajectory on a low-dim embedding with the discrete
# states as clusters, rooted at `start_clus`. `embedding` is a cell x dim
# matrix (the caller truncates harmony to n_dims); `clusters` is the per-
# cell state label (length = nrow(embedding)). `end_clus` optionally fixes
# terminal cluster(s); `omit_clusters` drops cells of those states BEFORE
# lineage inference (the documented cell-cycle mitigation -- omitted cells
# get NA pseudotime but are kept in the returned `clusters` so downstream
# reporting still sees them, guardrail #7). Slingshot's topology is
# INFERRED, never asserted (only the root is fixed).
#
# Returns list(pseudotime = cell x lineage matrix [NA off-lineage / omitted],
# weights = cell x lineage curve weights, lineages = named list of cluster
# sequences, mst, sds = the PseudotimeOrdering [for embedCurves/plots],
# clusters, terminal_states = per-lineage last cluster, dam_lineage = name
# of the lineage terminating in `dam_state` (NA if none/ambiguous), params).
build_microglia_trajectory <- function(embedding, clusters,
                                       start_clus = "homeostatic",
                                       end_clus = NULL, omit_clusters = NULL,
                                       dam_state = "DAM", seed = 1L) {
  stopifnot(requireNamespace("slingshot", quietly = TRUE),
            is.matrix(embedding), nrow(embedding) == length(clusters))
  clusters <- as.character(clusters)
  cells <- rownames(embedding) %||% as.character(seq_len(nrow(embedding)))
  rownames(embedding) <- cells
  keep <- if (is.null(omit_clusters)) rep(TRUE, length(clusters)) else
    !clusters %in% omit_clusters
  stopifnot(start_clus %in% clusters[keep])

  set.seed(seed)
  sds <- slingshot::slingshot(embedding[keep, , drop = FALSE],
                              clusterLabels = clusters[keep],
                              start.clus = start_clus, end.clus = end_clus)
  pt_sub <- slingshot::slingPseudotime(sds)            # kept-cell x lineage
  w_sub  <- slingshot::slingCurveWeights(sds)
  lineages <- slingshot::slingLineages(sds)

  # Re-expand to ALL cells (omitted / off-lineage rows = NA).
  pt <- matrix(NA_real_, length(cells), ncol(pt_sub),
               dimnames = list(cells, colnames(pt_sub)))
  w  <- matrix(NA_real_, length(cells), ncol(w_sub),
               dimnames = list(cells, colnames(w_sub)))
  pt[cells[keep], ] <- pt_sub
  w[cells[keep], ]  <- w_sub

  terminal <- vapply(lineages, function(z) z[length(z)], character(1))
  dam_lin  <- names(terminal)[terminal == dam_state]
  dam_lineage <- if (length(dam_lin) == 1L) dam_lin else NA_character_

  list(pseudotime = pt, weights = w, lineages = lineages,
       mst = slingshot::slingMST(sds), sds = sds,
       clusters = stats::setNames(clusters, cells),
       terminal_states = terminal, dam_lineage = dam_lineage,
       params = list(start_clus = start_clus, end_clus = end_clus,
                     omit_clusters = omit_clusters, dam_state = dam_state,
                     n_dims = ncol(embedding), n_cells = length(cells),
                     n_cells_kept = sum(keep), seed = seed))
}

# Per-cell differentiation-potency proxy from a genes x cells COUNTS matrix
# -- the root-free cross-check / fallback for CytoTRACE2 (guardrail #2).
# Convention: HIGHER = more potent / less differentiated (so the validated
# root should carry the MAX mean potency).
#   "entropy"  Shannon entropy of the library-normalised expression
#              (H = log S - sum_i x_i log x_i / S; more uniform => higher).
#   "n_genes"  number of detected genes (classic CytoTRACE intuition).
# Returns a named numeric vector (one per cell).
cell_potency <- function(counts, method = c("entropy", "n_genes")) {
  method <- match.arg(method)
  counts <- methods::as(counts, "CsparseMatrix")
  if (method == "n_genes") {
    return(stats::setNames(diff(counts@p), colnames(counts)))
  }
  S <- Matrix::colSums(counts)
  xlogx <- counts
  xlogx@x <- counts@x * log(counts@x)                  # 0 log 0 := 0 (no zeros stored)
  H <- log(S) - Matrix::colSums(xlogx) / S
  H[!is.finite(H)] <- NA_real_
  stats::setNames(as.numeric(H), colnames(counts))
}

# Validate the asserted root: does `root` carry the HIGHEST mean potency
# (least differentiated)? `potency` = per-cell numeric (any proxy where
# higher = more potent); `states` = per-cell label. Returns list(
# per_state = tibble(state, mean_potency, sd_potency, n, rank [1 = most
# potent]), root, root_rank, root_is_most_potent, delta_to_next [root mean
# minus the next-highest state mean; positive => root cleanly most potent],
# n_states). Reports disagreement rather than enforcing it (guardrail #2).
validate_root_potency <- function(potency, states, root = "homeostatic") {
  stopifnot(length(potency) == length(states))
  df <- data.frame(state = as.character(states), potency = as.numeric(potency),
                   stringsAsFactors = FALSE)
  df <- df[is.finite(df$potency), , drop = FALSE]
  agg <- df |>
    dplyr::group_by(.data$state) |>
    dplyr::summarise(mean_potency = mean(.data$potency),
                     sd_potency = stats::sd(.data$potency),
                     n = dplyr::n(), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(.data$mean_potency)) |>
    dplyr::mutate(rank = dplyr::row_number())
  root_rank <- agg$rank[agg$state == root]
  if (length(root_rank) == 0L) root_rank <- NA_integer_
  means <- agg$mean_potency[order(-agg$mean_potency)]
  root_mean <- agg$mean_potency[agg$state == root]
  next_mean <- max(agg$mean_potency[agg$state != root], na.rm = TRUE)
  list(per_state = tibble::as_tibble(agg), root = root,
       root_rank = root_rank,
       root_is_most_potent = isTRUE(root_rank == 1L),
       delta_to_next = if (length(root_mean)) root_mean - next_mean else NA_real_,
       n_states = nrow(agg))
}

# Per-replicate (genotype_batch) pseudotime summary on ONE lineage.
# `pseudotime` = per-cell pseudotime for the chosen lineage (a column of
# build_microglia_trajectory()$pseudotime; NA = off-lineage); `states` and
# `replicate` are per-cell. Cells with NA pseudotime are dropped (off the
# lineage). For each replicate it returns: n_on_lineage, mean_pt (optionally
# curve-weighted), frac_past (fraction past `dam_threshold` -- a progression
# measure; if NULL the threshold defaults to the GLOBAL median pseudotime of
# `dam_state` cells, a fixed data-driven onset), and within_<state> (mean
# pseudotime within each state level = the composition-decoupled advancement
# signal). `replicate_meta` (rownames = replicate id, cols genotype/batch)
# is joined when supplied.
#
# Returns list(by_replicate = tibble [replicate (+genotype/batch), ...],
# dam_threshold, state_levels).
pseudotime_per_replicate <- function(pseudotime, states, replicate,
                                     weights = NULL, dam_state = "DAM",
                                     dam_threshold = NULL, state_levels = NULL,
                                     replicate_meta = NULL) {
  stopifnot(length(pseudotime) == length(states),
            length(pseudotime) == length(replicate))
  on <- is.finite(pseudotime)
  pt <- pseudotime[on]; st <- as.character(states)[on]
  rp <- as.character(replicate)[on]
  wt <- if (is.null(weights)) rep(1, length(pt)) else weights[on]
  if (is.null(state_levels)) state_levels <- sort(unique(as.character(states)))
  if (is.null(dam_threshold))
    dam_threshold <- stats::median(pt[st == dam_state], na.rm = TRUE)

  ids <- sort(unique(rp))
  rows <- lapply(ids, function(r) {
    i <- rp == r
    base <- data.frame(
      replicate    = r,
      n_on_lineage = sum(i),
      mean_pt      = stats::weighted.mean(pt[i], wt[i], na.rm = TRUE),
      frac_past    = sum(i & pt >= dam_threshold) / sum(i),
      stringsAsFactors = FALSE)
    ws <- vapply(state_levels, function(s) {
      j <- i & st == s
      if (any(j)) mean(pt[j]) else NA_real_
    }, numeric(1))
    cbind(base, as.data.frame(as.list(stats::setNames(ws,
            paste0("within_", state_levels)))))
  })
  by_rep <- dplyr::bind_rows(rows)
  if (!is.null(replicate_meta)) {
    m <- as.data.frame(replicate_meta)
    by_rep$genotype <- m$genotype[match(by_rep$replicate, rownames(m))]
    by_rep$batch    <- m$batch[match(by_rep$replicate, rownames(m))]
  }
  list(by_replicate = tibble::as_tibble(by_rep),
       dam_threshold = dam_threshold, state_levels = state_levels)
}

# Decompose the per-replicate pseudotime into a COMPOSITION channel and a
# PROGRESSION channel (guardrail #1 -- the load-bearing piece). With state
# fraction pi_{r,s} and within-state mean mu_{r,s} per replicate r:
#   observed_r       = sum_s pi_{r,s} mu_{r,s}   (= plain mean pseudotime)
#   composition_cf_r = sum_s pi_{r,s} mu_bar_s   (vary COMPOSITION; within-
#                                                 state means held at global)
#   progression_cf_r = sum_s pi_bar_s mu*_{r,s}  (vary PROGRESSION; fractions
#                                                 held at global; mu* falls
#                                                 back to mu_bar_s when a
#                                                 replicate lacks state s)
# All three are on the pseudotime scale, so the build script can fit the
# locked factorial on rbind(observed, composition_cf, progression_cf) and
# read which channel carries the interaction: composition-only => it reduces
# to the arc-L composition result; progression => genuine within-state
# advancement. `states`/`replicate` per cell; NA-pseudotime cells dropped.
#
# Returns list(by_replicate = tibble(replicate, observed, composition_cf,
# progression_cf), global = tibble(state, pi_bar, mu_bar), n_dropped_na,
# state_levels).
decompose_progression_vs_composition <- function(pseudotime, states, replicate,
                                                 state_levels = NULL) {
  stopifnot(length(pseudotime) == length(states),
            length(pseudotime) == length(replicate))
  on <- is.finite(pseudotime)
  n_drop <- sum(!on)
  pt <- pseudotime[on]; st <- as.character(states)[on]
  rp <- as.character(replicate)[on]
  # Only states actually PRESENT on the lineage enter the partition: a state
  # entirely absent on-lineage (e.g. a cycling side-state the caller omitted)
  # has pi_bar = 0 and undefined mu_bar, and would poison the sums (0 * NA);
  # its true contribution is zero, so it is dropped.
  observed_states <- sort(unique(st))
  state_levels <- if (is.null(state_levels)) observed_states else
    intersect(state_levels, observed_states)

  pi_bar <- vapply(state_levels, function(s) mean(st == s), numeric(1))
  mu_bar <- vapply(state_levels, function(s)
    if (any(st == s)) mean(pt[st == s]) else NA_real_, numeric(1))

  ids <- sort(unique(rp))
  rows <- lapply(ids, function(r) {
    i <- rp == r; n <- sum(i)
    pi_r <- vapply(state_levels, function(s) sum(i & st == s) / n, numeric(1))
    mu_r <- vapply(state_levels, function(s) {
      j <- i & st == s; if (any(j)) mean(pt[j]) else NA_real_
    }, numeric(1))
    mu_star <- ifelse(is.na(mu_r), mu_bar, mu_r)        # neutral fallback
    data.frame(replicate = r,
               observed       = sum(pi_r * mu_star),
               composition_cf = sum(pi_r * mu_bar),
               progression_cf = sum(pi_bar * mu_star),
               stringsAsFactors = FALSE)
  })
  list(by_replicate = tibble::as_tibble(dplyr::bind_rows(rows)),
       global = tibble::tibble(state = state_levels,
                               pi_bar = pi_bar, mu_bar = mu_bar),
       n_dropped_na = n_drop, state_levels = state_levels)
}

# Fit the LOCKED 2x2 factorial + batch design to a measure x replicate
# summary matrix and extract all 5 contrasts per measure via limma-trend
# (mirrors fit_aucell_contrasts / fit_abundance_contrasts). `summary_mat`
# rows = measures (lineage mean pseudotimes, decomposition channels, ...),
# cols = the 16 genotype_batch ids; `meta` rownames = those ids with
# genotype + batch. transform "logit" variance-stabilises bounded [0,1]
# measures (fractions); default "none" for raw pseudotime. `se` (= logFC/t)
# is carried so the §23 forest can draw CIs.
#
# Returns tibble(measure, contrast, logFC, se, t, P.Value, adj.P.Val, sig,
#                transform); adj.P.Val = limma BH across measures within
# each contrast.
fit_trajectory_contrasts <- function(summary_mat, meta,
                                     transform = c("none", "logit"),
                                     padj_cut = 0.10) {
  transform <- match.arg(transform)
  stopifnot(all(colnames(summary_mat) %in% rownames(meta)))
  meta <- meta[colnames(summary_mat), , drop = FALSE]
  # Drop measure rows with any non-finite cell (limma-trend's eBayes covariate
  # rejects NA/Inf); a measure undefined for some replicate is not fittable.
  finite_rows <- apply(summary_mat, 1, function(z) all(is.finite(z)))
  if (any(!finite_rows))
    warning(sprintf("fit_trajectory_contrasts: dropping non-finite measure(s): %s",
                    paste(rownames(summary_mat)[!finite_rows], collapse = ", ")))
  summary_mat <- summary_mat[finite_rows, , drop = FALSE]
  stopifnot(nrow(summary_mat) >= 1L)
  mat <- switch(transform,
                none  = summary_mat,
                logit = { stopifnot(all(summary_mat > 0 & summary_mat < 1))
                          qlogis(summary_mat) })
  fd <- factorial_design(meta, genotype_col = "genotype",
                         batch_col = "batch", add_batch = TRUE)
  fit <- fit_limma_log(mat[, rownames(fd$design), drop = FALSE],
                       group = meta$genotype,
                       design = fd$design, contrasts = fd$contrasts)
  rows <- lapply(names(fit$top), function(cn) {
    tt <- fit$top[[cn]]
    se <- ifelse(tt$t == 0, NA_real_, tt$logFC / tt$t)
    data.frame(measure = tt$feature, contrast = cn, logFC = tt$logFC,
               se = se, t = tt$t, P.Value = tt$P.Value,
               adj.P.Val = tt$adj.P.Val, stringsAsFactors = FALSE)
  })
  res <- dplyr::bind_rows(rows)
  res$sig <- res$adj.P.Val < padj_cut
  res$transform <- transform
  tibble::as_tibble(res)
}

# Cross-method pseudotime concordance (guardrail #6: discordance is a
# finding, not smoothed). `pt_df` = cells x methods (a data frame / matrix;
# columns = method labels, e.g. slingshot_DAM, dpt, cellrank_dam_fate).
# Returns list(pairwise = tibble(method_a, method_b, rho, n [complete
# pairs]), cor_mat = method x method Spearman matrix). Per-cell, on the
# intersection of finite values for each pair.
pseudotime_concordance <- function(pt_df) {
  m <- as.matrix(pt_df)
  meth <- colnames(m)
  stopifnot(length(meth) >= 2)
  cor_mat <- suppressWarnings(stats::cor(m, method = "spearman",
                                         use = "pairwise.complete.obs"))
  rows <- list()
  for (a in seq_along(meth)) for (b in seq_along(meth)) if (a < b) {
    ok <- is.finite(m[, a]) & is.finite(m[, b])
    rows[[length(rows) + 1L]] <- data.frame(
      method_a = meth[a], method_b = meth[b],
      rho = cor_mat[a, b], n = sum(ok), stringsAsFactors = FALSE)
  }
  list(pairwise = tibble::as_tibble(dplyr::bind_rows(rows)), cor_mat = cor_mat)
}

# ---- differential dynamics: 2x2 interaction-along-pseudotime (arc O) --------
# tradeSeq NB-GAM with conditions = genotype fits one smoother of pseudotime per
# genotype on a single lineage. All four genotypes share that lineage's spline
# basis (same knots), so a difference of two genotypes' smoothers equals the
# difference of their basis-coefficient blocks EXACTLY. The tau x amyloid
# interaction-along-pseudotime is therefore the difference-of-differences of the
# four coefficient blocks, knot by knot:
#   L = (NLGF_P301S - P301S) - (NLGF_MAPTKI - MAPTKI)            (one column / knot)
# tested per gene by the same rank-aware Wald tradeSeq's conditionTest uses
# (getEigenStatGAMFC, l2fc = 0): W = (L'b)' (L' S L)^-1 (L'b), df = rank(L' S L)
# via an eigen pseudo-inverse. SAME inference core as conditionTest, with a
# custom 2x2 contrast in place of the omnibus all-pairwise contrast. The
# coefficient-space contrast is exactly equal (identical W, df) to the K-point
# "fitted-curve difference-of-differences" via predictGAM, because the K-point
# basis map has full column rank and the Wald is invariant under it; the O2 smoke
# test verifies this. Pinned to tradeSeq 1.24.0 (beta/Sigma rowData layout +
# predictGAM/.getPredictRangeDf internals).

# Thin fitGAM wrapper for ONE lineage with conditions = genotype. `counts` =
# genes x cells raw matrix, `pseudotime` = per-cell numeric, `conditions` =
# per-cell genotype (coerced to factor at genotype_levels), `genes` = optional
# row subset. Cells with non-finite pseudotime are dropped (off-lineage). Returns
# the fitted SingleCellExperiment (per-gene beta/Sigma in rowData(sce)$tradeSeq).
# `parallel`/`BPPARAM` pass straight through to fitGAM for the genome-wide O3 run;
# defaults (SerialParam) reproduce the O2 smoke-test path byte-for-byte — the
# per-gene NB-GAM fit is deterministic, so MulticoreParam == SerialParam.
fit_lineage_gam <- function(counts, pseudotime, conditions, genes = NULL,
                            nknots = 6L, verbose = FALSE,
                            parallel = FALSE, BPPARAM = BiocParallel::SerialParam()) {
  if (!is.null(genes)) counts <- counts[genes, , drop = FALSE]
  fin <- is.finite(pseudotime)
  cond <- factor(as.character(conditions)[fin], levels = genotype_levels)
  tradeSeq::fitGAM(
    counts      = as.matrix(counts[, fin, drop = FALSE]),
    pseudotime  = matrix(pseudotime[fin], ncol = 1),
    cellWeights = matrix(1, nrow = sum(fin), ncol = 1),
    conditions  = cond,
    nknots      = as.integer(nknots), verbose = verbose,
    parallel    = parallel, BPPARAM = BPPARAM)
}

# Build the coefficient-space difference-of-differences contrast L
# (ncol(X) x nknots; rownames = coefficient names; column k = the 2x2 interaction
# across the four genotypes' knot-k smoother coefficients). Weights default to
# the locked interaction derived from contrast_definitions (NLGF_P301S +1,
# P301S -1, NLGF_MAPTKI -1, MAPTKI +1). Smooth basis coefficients are named
# `s(t{lin}):l{lin}_{cond}.{knot}`; the (cond, knot) index is parsed from the
# name so the layout is not assumed. Internal.
.interaction_contrast_L <- function(sce, weights = NULL) {
  X <- SummarizedExperiment::colData(sce)$tradeSeq$X
  conditions <- SummarizedExperiment::colData(sce)$tradeSeq$conditions
  cond_levels <- levels(conditions)
  if (is.null(weights)) {
    weights <- stats::setNames(numeric(length(genotype_levels)), genotype_levels)
    ip <- contrast_definitions$nlgf_in_p301s   # c(NLGF_P301S, P301S)
    im <- contrast_definitions$nlgf_in_maptki  # c(NLGF_MAPTKI, MAPTKI)
    weights[ip[1]] <- weights[ip[1]] + 1       # +NLGF_P301S
    weights[ip[2]] <- weights[ip[2]] - 1       # -P301S
    weights[im[1]] <- weights[im[1]] - 1       # -NLGF_MAPTKI
    weights[im[2]] <- weights[im[2]] + 1       # +MAPTKI
  }
  smooth_idx <- grep("^s\\(t", colnames(X))
  nk <- length(smooth_idx) / length(cond_levels)
  stopifnot(nk == as.integer(nk))
  nk <- as.integer(nk)
  after  <- sub(".*:l", "", colnames(X)[smooth_idx])             # "{lin}_{cond}.{knot}"
  cond_i <- as.integer(sub("^[0-9]+_([0-9]+)\\..*$", "\\1", after))
  knot_i <- as.integer(sub("^.*\\.", "", after))
  L <- matrix(0, nrow = ncol(X), ncol = nk,
              dimnames = list(colnames(X), paste0("knot", seq_len(nk))))
  for (j in seq_along(smooth_idx))
    L[smooth_idx[j], knot_i[j]] <- weights[[cond_levels[cond_i[j]]]]
  attr(L, "weights") <- weights
  L
}

# Rank-aware Wald statistic for one gene; mirrors tradeSeq::getEigenStatGAMFC at
# l2fc = 0. `beta` = p x 1, `Sigma` = p x p, `L` = p x m. Returns c(stat, df);
# NA when the contrast covariance is singular or collapses to rank 1. Internal.
.wald_eigen <- function(beta, Sigma, L, eigen_thresh = 0.01) {
  est <- t(L) %*% beta
  sigma <- t(L) %*% Sigma %*% L
  eS <- eigen(sigma, symmetric = TRUE)
  r <- try(sum(eS$values / eS$values[1] > eigen_thresh), silent = TRUE)
  if (inherits(r, "try-error") || is.na(r) || r < 2) return(c(NA_real_, NA_real_))
  half <- eS$vectors[, seq_len(r), drop = FALSE] %*%
    diag(1 / sqrt(eS$values[seq_len(r)]), nrow = r, ncol = r)
  halfStat <- t(est) %*% half
  c(as.numeric(crossprod(t(halfStat))), r)
}

# Per-gene 2x2 interaction-along-pseudotime Wald over a fitted-GAM SCE. Returns a
# tibble(gene, symbol, waldStat, df, pvalue, adj.P.Val, effect_peak, effect_l2)
# arranged by pvalue. effect_peak = signed largest-magnitude per-knot diff-of-
# differences of the (log-link) smoother coefficients; effect_l2 = its L2 norm
# (overall interaction-dynamics magnitude). NA stats propagate to NA p-values
# (genes tradeSeq could not fit). `weights` overrides the locked interaction.
interaction_dynamics_contrast <- function(sce, symbol_map = NULL,
                                          weights = NULL, eigen_thresh = 0.01) {
  L <- .interaction_contrast_L(sce, weights = weights)
  betaAll  <- SummarizedExperiment::rowData(sce)$tradeSeq$beta[[1]]
  sigmaAll <- SummarizedExperiment::rowData(sce)$tradeSeq$Sigma
  genes <- rownames(sce)
  mat <- vapply(seq_len(nrow(sce)), function(ii) {
    beta <- matrix(as.numeric(betaAll[ii, ]), ncol = 1)
    Sigma <- sigmaAll[[ii]]
    if (anyNA(beta) || anyNA(Sigma)) return(c(NA_real_, NA_real_, NA_real_, NA_real_))
    w  <- .wald_eigen(beta, Sigma, L, eigen_thresh)
    dd <- as.numeric(t(L) %*% beta)                       # per-knot diff-of-diffs (link)
    c(w, dd[which.max(abs(dd))], sqrt(sum(dd^2)))
  }, numeric(4))
  mat <- t(mat)
  sym <- if (!is.null(symbol_map))
    symbol_map$symbol[match(genes, symbol_map$ensembl)] else NA_character_
  out <- tibble::tibble(
    gene = genes, symbol = sym,
    waldStat = mat[, 1], df = mat[, 2],
    pvalue = 1 - stats::pchisq(mat[, 1], df = mat[, 2]),
    effect_peak = mat[, 3], effect_l2 = mat[, 4])
  out$adj.P.Val <- stats::p.adjust(out$pvalue, "BH")
  out[order(out$pvalue), ]
}

# Side-by-side companions to the targeted interaction (anti-cherry-pick guardrail):
# omnibus = "do the four genotype smoothers differ at all" (conditionTest global);
# association = "is the gene dynamic along pseudotime at all" (associationTest).
omnibus_dynamics <- function(sce, symbol_map = NULL) {
  ct <- tradeSeq::conditionTest(sce, global = TRUE, pairwise = FALSE)
  ct$gene <- rownames(ct)
  ct$adj.P.Val <- stats::p.adjust(ct$pvalue, "BH")
  if (!is.null(symbol_map))
    ct$symbol <- symbol_map$symbol[match(ct$gene, symbol_map$ensembl)]
  tibble::as_tibble(ct)
}
association_dynamics <- function(sce, symbol_map = NULL) {
  at <- tradeSeq::associationTest(sce, global = TRUE, lineages = FALSE)
  at$gene <- rownames(at)
  at$adj.P.Val <- stats::p.adjust(at$pvalue, "BH")
  if (!is.null(symbol_map))
    at$symbol <- symbol_map$symbol[match(at$gene, symbol_map$ensembl)]
  tibble::as_tibble(at)
}

# Per-genotype fitted smoother + 95% CI for one gene, all conditions evaluated on
# a COMMON pseudotime grid (so the difference-of-differences is visually
# comparable) over the lineage's overall range. Returns long
# tibble(pseudotime, genotype, fit, lo, hi, se, gene) on the response (count)
# scale. fit = exp(X_c b) at the mean offset; CI = exp(X_c b +/- 1.96 se). Uses
# the same predictGAM design route the smoke test validates against the contrast.
extract_condition_smoothers <- function(sce, gene, K = 100L) {
  X  <- SummarizedExperiment::colData(sce)$tradeSeq$X
  dm <- SummarizedExperiment::colData(sce)$tradeSeq$dm
  conditions  <- SummarizedExperiment::colData(sce)$tradeSeq$conditions
  cond_levels <- levels(conditions)
  crv <- SummarizedExperiment::colData(sce)$crv
  pt  <- crv[, grep("pseudotime", colnames(crv))]
  pt  <- if (is.null(dim(pt))) matrix(pt, ncol = 1) else as.matrix(pt)
  id  <- match(gene, rownames(sce))
  beta  <- matrix(as.numeric(SummarizedExperiment::rowData(sce)$tradeSeq$beta[[1]][id, ]), ncol = 1)
  Sigma <- SummarizedExperiment::rowData(sce)$tradeSeq$Sigma[[id]]
  rng  <- range(pt[is.finite(pt[, 1]), 1])
  grid <- seq(rng[1], rng[2], length.out = K)
  rows <- lapply(seq_along(cond_levels), function(ci) {
    df <- tradeSeq:::.getPredictRangeDf(dm, lineageId = 1, conditionId = ci, nPoints = K)
    df[, grep("^t[1-9]", colnames(df))] <- grid     # share the grid across genotypes
    Xc <- tradeSeq:::predictGAM(lpmatrix = X, df = df, pseudotime = pt, conditions = conditions)
    fit <- as.numeric(Xc %*% beta)
    se  <- sqrt(pmax(0, diag(Xc %*% Sigma %*% t(Xc))))
    data.frame(pseudotime = grid, genotype = cond_levels[ci],
               fit = exp(fit), lo = exp(fit - 1.96 * se), hi = exp(fit + 1.96 * se),
               se = se, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  out$genotype <- factor(out$genotype, levels = genotype_levels)
  out$gene <- gene
  tibble::as_tibble(out)
}

# ---- plot helpers (render INSIDE the shared knit session, per rmd/15,20) ----

# UMAP coloured by a discrete state factor or a continuous pseudotime, with
# optional Slingshot lineage curves overlaid. `coords` = cell x 2 (umap);
# `colour` = per-cell factor (discrete) or numeric (pseudotime); `curves` =
# optional list of 2-col matrices (embedCurves(sds, umap)$<lineage> $s) drawn
# as ordered paths.
plot_trajectory_umap <- function(coords, colour, curves = NULL, title = NULL,
                                 colour_lab = NULL, point_size = 0.3) {
  df <- data.frame(x = coords[, 1], y = coords[, 2], col = colour)
  p <- ggplot2::ggplot(df, ggplot2::aes(x, y, colour = col)) +
    ggplot2::geom_point(size = point_size, alpha = 0.6) +
    ggplot2::labs(x = "UMAP-1", y = "UMAP-2",
                  colour = colour_lab %||% "", title = title) +
    ggplot2::theme_classic()
  if (is.numeric(colour))
    p <- p + ggplot2::scale_colour_viridis_c()
  if (!is.null(curves)) for (cv in curves) {
    cd <- as.data.frame(cv); names(cd)[1:2] <- c("x", "y")
    p <- p + ggplot2::geom_path(data = cd, ggplot2::aes(x, y),
                                inherit.aes = FALSE, linewidth = 0.9)
  }
  p
}

# Per-cell pseudotime density split by genotype (the 2x2 read). `df` carries
# a pseudotime column + a genotype column; uses the project genotype palette.
plot_pseudotime_density <- function(df, pt_col = "pseudotime",
                                    genotype_col = "genotype", title = NULL) {
  d <- data.frame(pt = df[[pt_col]],
                  genotype = factor(df[[genotype_col]], levels = genotype_levels))
  d <- d[is.finite(d$pt), , drop = FALSE]
  ggplot2::ggplot(d, ggplot2::aes(pt, colour = genotype, fill = genotype)) +
    ggplot2::geom_density(alpha = 0.15, linewidth = 0.8) +
    ggplot2::scale_colour_manual(values = genotype_colours) +
    ggplot2::scale_fill_manual(values = genotype_colours) +
    ggplot2::labs(x = "pseudotime (homeostatic -> DAM)", y = "density",
                  title = title) +
    ggplot2::theme_classic()
}

# Forest of the 5 canonical contrasts (logFC +/- 95% CI) for one measure of
# a fit_trajectory_contrasts() table. CI = logFC +/- 1.96 * se; the vertical
# guide at 0 is the no-effect line; significant contrasts (adj.P<padj_cut)
# are filled. The `interaction` contrast is the divergence headline.
plot_interaction_forest <- function(contrast_tbl, measure,
                                    contrast_order = c("tau_alone",
                                      "nlgf_in_maptki", "nlgf_in_p301s",
                                      "tau_in_nlgf", "interaction"),
                                    padj_cut = 0.10, title = NULL) {
  d <- contrast_tbl[contrast_tbl$measure == measure, , drop = FALSE]
  d$contrast <- factor(d$contrast, levels = rev(contrast_order))
  d$lo <- d$logFC - 1.96 * d$se
  d$hi <- d$logFC + 1.96 * d$se
  d$signif <- !is.na(d$adj.P.Val) & d$adj.P.Val < padj_cut
  ggplot2::ggplot(d, ggplot2::aes(logFC, contrast)) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = lo, xmax = hi),
                           orientation = "y", width = 0.2) +
    ggplot2::geom_point(ggplot2::aes(fill = signif), shape = 21, size = 3) +
    ggplot2::scale_fill_manual(values = c(`TRUE` = "#B0344D", `FALSE` = "white"),
                               name = sprintf("adj.P<%.2f", padj_cut)) +
    ggplot2::labs(x = "logFC (pseudotime units)", y = NULL,
                  title = title %||% measure) +
    ggplot2::theme_classic()
}

# Per-genotype fitted NB-GAM expression smoothers along pseudotime -- the visual
# of the difference-of-differences (interaction) dynamics contrast. `df` is the
# long table produced by extract_condition_smoothers() (and pre-stored in
# trajectory_dynamics.rds$smoothers): columns pseudotime, genotype, fit, lo, hi.
# A facet column (default `symbol`) panels multiple genes with free y-scales
# (fitted expression scales differ across genes). The four genotype smoothers
# differing in shape IS the interaction; the project genotype palette keeps the
# 2x2 readable (cool = amyloid-free, warm = amyloid). Renders inside the shared
# knit session, per the rmd/15,20 convention.
plot_condition_smoothers <- function(df, panel_col = "symbol", title = NULL,
                                     ncol = 3L) {
  d <- df
  d$genotype <- factor(d$genotype, levels = genotype_levels)
  p <- ggplot2::ggplot(d, ggplot2::aes(pseudotime, fit,
                                       colour = genotype, fill = genotype)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi),
                         alpha = 0.12, colour = NA) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::scale_colour_manual(values = genotype_colours) +
    ggplot2::scale_fill_manual(values = genotype_colours) +
    ggplot2::labs(x = "pseudotime (homeostatic -> DAM)",
                  y = "fitted expression (NB-GAM)", title = title) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "bottom")
  if (!is.null(panel_col) && panel_col %in% names(d))
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste0("~ ", panel_col)),
                                 scales = "free_y", ncol = ncol)
  p
}
