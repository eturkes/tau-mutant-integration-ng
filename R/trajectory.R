# Activation-trajectory pseudotime (P2-S1). Order microglia along the homeostatic->DAM
# axis with slingshot on the BATCH-CORRECTED harmony embedding, seeded by the substate
# biology. PRIMARY lineage = a forced single clean Homeostatic->DAM curve (clusterLabels =
# the 2 substate super-clusters -> MST on 2 nodes = one edge = one lineage; principal curve
# still fits the full high-dim cell cloud). IFN/Proliferative cells are OMITTED from the
# lineage (off-lineage flag + NA pt), not deleted -> per-unit omitted fraction is reported,
# never hidden. Pseudotime is an ACTIVATION ORDERING (position/extent of advance), NOT
# developmental time or potency -- direction rests on slingshot rooting + DAM-marker
# rank association, validated post-hoc, never a potency claim. All non-base calls namespace-
# qualified (targets attaches only `quarto`). Pure helpers unit-tested; the heavy
# orchestrator (build_activation_trajectory) is smoke-tested live on microglia_annotated.

# Marker SCORE-AXIS pseudotime proxy = raw UCell DAM minus Homeostatic. Assumption-light
# (no population z-centring): both are UCell scores on the same rank-based [0,1] scale, so
# their difference is a direct per-cell activation contrast. Used ONLY as a CONCORDANCE
# anchor for the slingshot ordering (it shares the marker SYSTEM that defines the substate
# labels -> catches gross trajectory failure, NOT a statistically independent robustness
# check). Pure.
score_axis_pseudotime <- function(dam_score, homeo_score) {
  stopifnot(is.numeric(dam_score), is.numeric(homeo_score),
            length(dam_score) == length(homeo_score),
            all(is.finite(dam_score)), all(is.finite(homeo_score)))
  dam_score - homeo_score
}

# Smithson-Verkuilen squeeze of a pseudotime vector into the OPEN interval (0,1) for a beta
# GLMM (P2-S3). Two steps: (1) min-max scale finite values to [0,1]; (2) y' = (y*(n-1)+0.5)/n
# with n = number of finite obs -> nudges the 0/1 endpoints inward (~0.5/n) so beta's
# unbounded logit stays finite. NA preserved as NA (off-lineage cells). Degenerate constant
# input -> all 0.5. Pure.
squeeze_unit_interval <- function(x) {
  stopifnot(is.numeric(x))
  ok <- is.finite(x)
  n <- sum(ok)
  stopifnot(n >= 2L)
  rng <- range(x[ok])
  scaled <- if (diff(rng) > 0) (x - rng[1L]) / diff(rng) else rep(0.5, length(x))
  out <- (scaled * (n - 1L) + 0.5) / n
  out[!ok] <- NA_real_
  out
}

# Spearman concordance between the slingshot pseudotime and the marker score-axis, over the
# cells where BOTH are finite (off-lineage pt is NA). Returns list(rho, n) -- rho recorded in
# provenance (concordance is reported, not gated; a low rho is an honest flag, not a build
# failure). Pure.
trajectory_concordance <- function(pt, score_axis) {
  stopifnot(is.numeric(pt), is.numeric(score_axis), length(pt) == length(score_axis))
  ok <- is.finite(pt) & is.finite(score_axis)
  stopifnot(sum(ok) >= 2L,                             # >=2 finite pairs for a correlation
            length(unique(pt[ok])) >= 2L,              # constant pt or score -> Spearman undefined
            length(unique(score_axis[ok])) >= 2L)      # (fail loud, not a silent NA rho)
  list(rho = stats::cor(pt[ok], score_axis[ok], method = "spearman"), n = sum(ok))
}

# Fit slingshot on a reduced-dim embedding + cluster labels, return the pseudotime of the
# lineage TERMINATING in `terminal_clus` (Homeostatic->DAM). matrix in, list out -> pure +
# unit-testable on a synthetic embedding. 2-cluster labels guarantee a single lineage; >2
# (all-retained sensitivity, IFN present) may branch -> the DAM-terminal lineage is selected
# (longest if tied) and off-lineage cells keep slingshot's NA. Rooted at `start_clus` ->
# pseudotime increases away from the homeostatic root. RNG seed + all three kinds pinned for
# the fit, and the caller's RNG state restored on exit (pure); the curve is RNG-free given inputs.
run_slingshot_lineage <- function(embedding, cluster_labels, start_clus, terminal_clus,
                                   seed = 42L) {
  stopifnot(is.matrix(embedding), !is.null(rownames(embedding)), ncol(embedding) >= 2L,
            all(is.finite(embedding)),
            length(cluster_labels) == nrow(embedding),
            start_clus %in% cluster_labels, terminal_clus %in% cluster_labels,
            all(table(as.character(cluster_labels)) >= 2L)) # >=2 cells/cluster (covariance; singular-cov gotcha)
  # Pin seed + all three RNG kinds for the fit, then RESTORE the caller's RNG state on exit
  # (kind + .Random.seed) so the helper leaves no random-stream side effect downstream.
  old_kind <- RNGkind()
  has_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (has_seed) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    RNGkind(old_kind[1], old_kind[2], old_kind[3])
    if (has_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion",
           sample.kind = "Rejection")
  res <- slingshot::slingshot(embedding, clusterLabels = as.character(cluster_labels),
                              start.clus = start_clus)
  lins <- slingshot::slingLineages(res)
  term_ok <- vapply(lins, function(l) l[length(l)] == terminal_clus, logical(1))
  if (!any(term_ok)) {
    stop("run_slingshot_lineage: no lineage terminates in '", terminal_clus,
         "' (lineages: ", paste(vapply(lins, paste, character(1), collapse = "->"),
                                collapse = " | "), ")")
  }
  idx <- which(term_ok)
  if (length(idx) > 1L) idx <- idx[which.max(vapply(lins[idx], length, integer(1)))]
  lin_name <- names(lins)[idx]
  pt_mat <- slingshot::slingPseudotime(res)
  pt <- pt_mat[, lin_name]
  names(pt) <- rownames(pt_mat)
  list(pt = pt, lineage = unname(lins[[idx]]), lineage_name = lin_name,
       n_lineages = length(lins), all_lineages = lapply(lins, unname))
}

# Pinned-stack + run-environment provenance for the trajectory build (mirror of
# reprocess_provenance): pkg versions, seed/RNG, thread snapshot, R version. Pure.
trajectory_provenance <- function(seed) {
  pkgs <- c("slingshot", "princurve", "TrajectoryUtils", "SingleCellExperiment", "UCell")
  list(
    versions = stats::setNames(
      vapply(pkgs, function(p) as.character(utils::packageVersion(p)), character(1)), pkgs),
    seed = seed,
    rng_kind = RNGkind(),
    threads = reprocess_thread_env(),
    r_version = as.character(getRversion())
  )
}

# Fail-loud audit of per-unit metadata before the per-unit omitted-fraction summary: every
# cell needs a non-missing replicate unit + a genotype in the declared levels, and each unit
# must carry EXACTLY ONE genotype (the summary takes the first genotype per unit, so mixed
# metadata would silently corrupt the audit). Pure + data.frame-testable (no Seurat needed).
validate_trajectory_units <- function(unit_vec, geno_vec, levels = genotype_levels) {
  stopifnot(is.character(unit_vec), is.character(geno_vec),
            length(unit_vec) == length(geno_vec),
            !anyNA(unit_vec), all(nzchar(unit_vec)),               # no missing/blank unit id
            !anyNA(geno_vec), all(geno_vec %in% levels),          # genotype in declared levels
            all(tapply(geno_vec, unit_vec,                        # one genotype per unit
                       function(g) length(unique(g))) == 1L))
  invisible(TRUE)
}

# Orchestrate the activation-trajectory target from microglia_annotated. Reads the cached
# harmony embedding + per-cell UCell scores + substate labels (NO recompute), fits the
# PRIMARY Homeostatic->DAM slingshot pseudotime (dims 1:primary_dims), a FIXED sensitivity
# table (alt dims + an all-cells-retained re-fit), and the marker score-axis concordance.
# Returns a COMPACT target (per-cell frame + per-unit table + lineage/sensitivity/provenance
# lists) -- NEVER the 612MB Seurat (cheap-render invariant; only the report-side qmd must stay
# light, but the target itself is small too). FORK CHOICES (dims, lineage substates) are
# PRE-DECLARED here, not retuned after inspecting contrasts.
build_activation_trajectory <- function(seurat_obj,
                                         primary_dims = 15L,
                                         sensitivity_dims = c(10L, 20L),
                                         lineage_substates = c("Homeostatic", "DAM"),
                                         root_substate = "Homeostatic",
                                         terminal_substate = "DAM",
                                         reduction = "harmony",
                                         substate_col = "microglia_substate",
                                         dam_col = "DAM_UCell",
                                         homeo_col = "Homeostatic_UCell",
                                         unit_col = "genotype_batch",
                                         genotype_col = "genotype",
                                         concordance_floor = 0.5,
                                         seed = 42L) {
  stopifnot(
    inherits(seurat_obj, "Seurat"),
    reduction %in% SeuratObject::Reductions(seurat_obj),
    all(c(substate_col, dam_col, homeo_col, unit_col, genotype_col) %in%
          colnames(seurat_obj@meta.data)),
    is.numeric(primary_dims), length(primary_dims) == 1L,
    primary_dims >= 2L, primary_dims == round(primary_dims),       # positive whole-number dims
    is.numeric(sensitivity_dims), length(sensitivity_dims) >= 1L,
    all(is.finite(sensitivity_dims)), all(sensitivity_dims >= 2L),
    all(sensitivity_dims == round(sensitivity_dims)))
  md <- seurat_obj@meta.data
  emb_full <- SeuratObject::Embeddings(seurat_obj, reduction)
  stopifnot(identical(rownames(emb_full), rownames(md)),
            ncol(emb_full) >= max(c(primary_dims, sensitivity_dims)),
            is.numeric(md[[dam_col]]), is.numeric(md[[homeo_col]])) # score cols numeric, not factor

  sub <- as.character(md[[substate_col]])
  dam <- as.numeric(md[[dam_col]]); homeo <- as.numeric(md[[homeo_col]])
  on_lineage <- sub %in% lineage_substates
  stopifnot(sum(on_lineage) >= 2L,
            root_substate %in% sub[on_lineage], terminal_substate %in% sub[on_lineage])

  # PRIMARY: forced single Homeostatic->DAM lineage on the on-lineage cells, dims 1:primary.
  emb_on <- emb_full[on_lineage, seq_len(primary_dims), drop = FALSE]
  prim <- run_slingshot_lineage(emb_on, sub[on_lineage],
                                root_substate, terminal_substate, seed)
  pt_raw <- stats::setNames(rep(NA_real_, nrow(md)), rownames(md))
  pt_raw[names(prim$pt)] <- prim$pt

  # marker score-axis (all retained cells) + slingshot-vs-score-axis concordance (on-lineage).
  score_axis_pt <- score_axis_pseudotime(dam, homeo)
  conc <- trajectory_concordance(pt_raw, score_axis_pt)

  # pt01 squeeze (on-lineage cells only -> off-lineage stay NA).
  pt01 <- squeeze_unit_interval(pt_raw)

  # SENSITIVITY: alt dims re-fit (Spearman vs primary on shared on-lineage cells) + an
  # all-cells-retained re-fit (IFN/prolif present -> may branch; take the DAM-terminal
  # lineage) to gauge lineage-conditioning selection effects.
  sens_rows <- list()
  for (d in sensitivity_dims) {
    alt <- run_slingshot_lineage(emb_full[on_lineage, seq_len(d), drop = FALSE],
                                 sub[on_lineage], root_substate, terminal_substate, seed)
    ok <- is.finite(alt$pt) & is.finite(prim$pt[names(alt$pt)])
    sens_rows[[length(sens_rows) + 1L]] <- data.frame(
      variant = paste0("dims_", d),
      spearman_vs_primary = stats::cor(prim$pt[names(alt$pt)][ok], alt$pt[ok],
                                       method = "spearman"),
      n_lineages = alt$n_lineages, n_shared = sum(ok), stringsAsFactors = FALSE)
  }
  allret <- run_slingshot_lineage(emb_full[, seq_len(primary_dims), drop = FALSE],
                                  sub, root_substate, terminal_substate, seed)
  shared <- intersect(names(prim$pt), names(allret$pt)[is.finite(allret$pt)])
  sens_rows[[length(sens_rows) + 1L]] <- data.frame(
    variant = "all_retained", n_lineages = allret$n_lineages, n_shared = length(shared),
    spearman_vs_primary = stats::cor(prim$pt[shared], allret$pt[shared], method = "spearman"),
    stringsAsFactors = FALSE)
  sensitivity <- do.call(rbind, lapply(sens_rows,
    function(r) r[c("variant", "spearman_vs_primary", "n_lineages", "n_shared")]))

  # fail-loud per-unit metadata audit before the frame (mixed/missing unit or off-levels
  # genotype would silently corrupt the per-unit omitted-fraction summary below).
  unit_vec <- as.character(md[[unit_col]]); geno_vec <- as.character(md[[genotype_col]])
  validate_trajectory_units(unit_vec, geno_vec)

  # per-cell compact frame (all retained cells; off-lineage pt NA, score-axis always defined).
  cell_frame <- data.frame(
    cell = rownames(md),
    genotype_batch = unit_vec,
    genotype = factor(geno_vec, levels = genotype_levels),
    substate = factor(sub, levels = sort(unique(sub))),
    on_lineage = on_lineage,
    pt_raw = pt_raw, pt01 = pt01, score_axis_pt = score_axis_pt,
    DAM_UCell = dam, Homeostatic_UCell = homeo,
    row.names = NULL, stringsAsFactors = FALSE, check.names = FALSE)

  # per-unit on-lineage / omitted fraction (16 genotype_batch units; conditioning-selection
  # audit -- a genotype-skewed omitted fraction would bias the summaries).
  u <- factor(cell_frame$genotype_batch)
  per_unit <- data.frame(
    genotype_batch = levels(u),
    genotype = vapply(levels(u), function(g) as.character(cell_frame$genotype[u == g][1]),
                      character(1)),
    n_cells = as.integer(table(u)),
    n_on_lineage = as.integer(tapply(cell_frame$on_lineage, u, sum)),
    row.names = NULL, stringsAsFactors = FALSE)
  per_unit$omitted_frac <- 1 - per_unit$n_on_lineage / per_unit$n_cells

  prov <- trajectory_provenance(seed)
  prov$primary_dims <- primary_dims
  prov$sensitivity_dims <- sensitivity_dims
  prov$lineage_substates <- lineage_substates
  prov$root_substate <- root_substate
  prov$terminal_substate <- terminal_substate
  prov$reduction <- reduction
  prov$concordance_rho <- conc$rho
  prov$concordance_n <- conc$n
  # gross-failure floor (0.5 = conventional large-correlation boundary), RECORDED not gated:
  # rho below it (or negative) means the transcriptome ordering disagrees with marker biology
  # -> trajectory suspect. A moderate rho is EXPECTED (slingshot = transcriptome geometry vs
  # score-axis = marker contrast; they are related, not identical), so this is an honest flag.
  prov$concordance_floor <- concordance_floor
  prov$concordant <- isTRUE(conc$rho >= concordance_floor)
  # directional rooting (recorded here, GATED in the postconditions): DAM up, Homeostatic down.
  on0 <- is.finite(pt_raw)
  prov$dam_pt_rho <- stats::cor(pt_raw[on0], dam[on0], method = "spearman")
  prov$homeo_pt_rho <- stats::cor(pt_raw[on0], homeo[on0], method = "spearman")
  prov$n_on_lineage <- sum(on_lineage)
  prov$n_omitted <- sum(!on_lineage)
  prov$omitted_frac_overall <- mean(!on_lineage)

  lineage <- list(primary = prim$lineage, primary_n_lineages = prim$n_lineages,
                  all_retained = allret$all_lineages,
                  all_retained_dam_lineage = allret$lineage)

  # build-time postconditions (the gate's warn=2 unit tests skip this heavy body -> the S1
  # acceptance INVARIANTS are pinned here; a silent regression fails tar_make, gate red):
  on <- is.finite(pt_raw)
  stopifnot(
    prim$n_lineages == 1L,                                   # single clean primary lineage
    identical(prim$lineage, c(root_substate, terminal_substate)),
    sum(on) == sum(on_lineage), !anyNA(pt_raw[on_lineage]),  # every on-lineage cell ordered
    prov$dam_pt_rho > 0, prov$homeo_pt_rho < 0,              # DAM rises, Homeo falls along pt (rooting)
    all(pt01[on] > 0 & pt01[on] < 1),                        # squeeze opened the interval
    is.finite(conc$rho),                                     # concordance recorded
    nrow(cell_frame) == nrow(md), nrow(per_unit) == nlevels(u),
    sum(per_unit$n_cells) == nrow(cell_frame),              # per-unit partition covers all cells
    all(is.finite(sensitivity$spearman_vs_primary)))        # every sensitivity recorded

  list(cell_frame = cell_frame, per_unit = per_unit,
       lineage = lineage, sensitivity = sensitivity, provenance = prov)
}

# ============================================================================================
# P2-S2a: per-replicate summary + contrast fit + Kitagawa decomposition (estimation core).
# Collapse per-cell pseudotime to the 16 genotype_batch replicate summaries, then run the 2x2
# factorial interaction as the TRAJECTORY ANALOGUE of the pseudobulk DE -- reuse factorial_design
# + its 5 canonical contrasts (ordinary t at 16 - 7 = 9 residual df, NO eBayes: a handful of
# heterogeneous pseudotime endpoints, not many shrinkable features). The 3-channel Kitagawa /
# Oaxaca shift-share (composition / progression / cross) splits the interaction EXACTLY on the
# additive pt_raw scale. Pure + deterministic-fixture-testable; NO target here (S2b wires it).
# All non-base calls namespace-qualified.

# ONE within-state column sanitizer: S2a builds the within_<state> columns, S2b reads them back
# by the SAME rule -> never hand-case the state label (a verbatim-cased lookup silently misses).
within_state_col <- function(state) paste0("within_", tolower(state))

# Extract the batch label from a genotype_batch id, VECTORISED over its length. Genotypes
# themselves contain "_" (e.g. NLGF_P301S) -> a naive strsplit/regex on "_" mis-splits, and base
# sub() does NOT vectorise over a length>1 pattern (warns + uses pattern[1] -> wrong under
# warn=2). Strip the per-element literal "<genotype>_" prefix instead, then round-trip ASSERT.
derive_batch <- function(genotype_batch, genotype) {
  stopifnot(is.character(genotype_batch),
            length(genotype_batch) == length(genotype))
  geno   <- as.character(genotype)
  prefix <- paste0(geno, "_")                              # per-element literal prefix
  stopifnot(all(startsWith(genotype_batch, prefix)))
  batch <- substring(genotype_batch, nchar(prefix) + 1L)   # vectorised over text + start
  stopifnot(all(nzchar(batch)),
            identical(paste(geno, batch, sep = "_"), genotype_batch))  # round-trip, fail loud
  batch
}

# Collapse on-lineage per-cell pseudotime to per-replicate (genotype_batch) summaries + the
# state x unit composition (pi) and within-state mean (mu) matrices the Kitagawa decomposition
# consumes. Filters to on-lineage cells (finite pt_raw); asserts every cell is in a lineage state
# + the per-unit metadata is clean. dam_onset (pooled DAM median) is PRE-DECLARED -> frac_past =
# the only genuinely [0,1] measure. within_<state> per-unit means stay on the raw pt scale
# (UNtransformed) so they share additive units with the Kitagawa progression channel. within_skip
# flags a state whose smallest unit count is below the floor (S2b drops its within-state measure).
# Pure.
pseudotime_per_replicate <- function(cell_frame, lineage_states, dam_state = "DAM",
                                     min_within = 10L) {
  stopifnot(is.data.frame(cell_frame),
            all(c("genotype_batch", "genotype", "substate", "pt_raw") %in% names(cell_frame)),
            is.character(lineage_states), dam_state %in% lineage_states)
  cf   <- cell_frame[is.finite(cell_frame$pt_raw), , drop = FALSE]   # on-lineage cells
  sub  <- as.character(cf$substate)
  unit <- as.character(cf$genotype_batch)
  geno <- as.character(cf$genotype)
  stopifnot(all(sub %in% lineage_states))                  # on-lineage cells are all lineage states
  validate_trajectory_units(unit, geno)                    # one geno/unit, geno in levels, no blank unit

  states <- lineage_states
  units  <- sort(unique(unit), method = "radix")           # locale-independent reproducible order
  f_state <- factor(sub,  levels = states)
  f_unit  <- factor(unit, levels = units)

  cnt <- table(f_state, f_unit)                            # state x unit counts
  cnt <- matrix(as.integer(cnt), nrow = length(states), ncol = length(units),
                dimnames = list(states, units))
  stopifnot(all(cnt >= 1L))                                # state present in every unit (Kitagawa needs mu_su)

  pi_mat <- sweep(cnt, 2, colSums(cnt), "/")               # column-normalised composition
  mu_mat <- tapply(cf$pt_raw, list(f_state, f_unit), mean)
  sd_mat <- tapply(cf$pt_raw, list(f_state, f_unit), stats::sd)
  dimnames(mu_mat) <- list(states, units)                  # pin identical dimnames (Kitagawa asserts it)
  dimnames(sd_mat) <- list(states, units)
  pi_bar <- rowSums(cnt) / sum(cnt)                        # CELL-weighted pooled composition
  # per-state pooled mean (states order). as.numeric STRIPS tapply's 1-D `dim` attribute: a 1-D
  # array * matrix is "non-conformable" (both carry dims), whereas a plain vector recycles.
  mu_bar <- stats::setNames(as.numeric(tapply(cf$pt_raw, f_state, mean)[states]), states)
  within_skip <- apply(cnt, 1, function(r) any(r < min_within))

  dam_onset <- stats::median(cf$pt_raw[sub == dam_state])  # PRE-DECLARED progression landmark

  geno_by_unit  <- vapply(units, function(u) geno[unit == u][1L], character(1))
  per_unit <- data.frame(
    genotype_batch = units,
    genotype  = geno_by_unit,
    batch     = derive_batch(units, geno_by_unit),
    n_cells   = as.integer(colSums(cnt)),
    sd_pt     = as.numeric(tapply(cf$pt_raw, f_unit, stats::sd)),
    mean_pt   = as.numeric(tapply(cf$pt_raw, f_unit, mean)),
    median_pt = as.numeric(tapply(cf$pt_raw, f_unit, stats::median)),
    q90       = as.numeric(tapply(cf$pt_raw, f_unit,
                                  function(v) stats::quantile(v, 0.9, names = FALSE))),
    frac_past = as.numeric(tapply(cf$pt_raw, f_unit, function(v) mean(v > dam_onset))),
    row.names = NULL, stringsAsFactors = FALSE, check.names = FALSE)
  for (s in states) per_unit[[within_state_col(s)]] <- as.numeric(mu_mat[s, units])

  list(per_unit = per_unit, states = states, units = units, counts = cnt,
       pi = pi_mat, mu = mu_mat, sd = sd_mat, pi_bar = pi_bar, mu_bar = mu_bar,
       dam_onset = dam_onset, within_skip = within_skip, min_within = min_within)
}

# Ordinary (UN-moderated) t table for one contrast of a contrasts.fit'd limma fit. topTable
# REQUIRES eBayes -> compute the ordinary t by hand (se = sigma * stdev.unscaled, df = residual
# df). For a handful of heterogeneous pseudotime endpoints the moderated/shrunk variance is
# incoherent -> plain OLS t at 9 df. Returns one row per measure (fit feature). Pure.
ordinary_t_table <- function(fit, contrast_name, conf_level = 0.95) {
  stopifnot(contrast_name %in% colnames(fit$coefficients),
            length(conf_level) == 1L, conf_level > 0, conf_level < 1)
  coef <- fit$coefficients[, contrast_name]
  se   <- fit$sigma * fit$stdev.unscaled[, contrast_name]
  df   <- fit$df.residual
  t    <- coef / se
  q    <- stats::qt(1 - (1 - conf_level) / 2, df)
  data.frame(
    measure = rownames(fit$coefficients), contrast = contrast_name,
    coef = coef, se = se, t = t, df = df,
    p_value = 2 * stats::pt(-abs(t), df), ci_l = coef - q * se, ci_r = coef + q * se,
    row.names = NULL, stringsAsFactors = FALSE)
}

# Fit a measures x units matrix against the factorial design + 5 contrasts (ordinary t, NO
# eBayes). limma consumes weights BY POSITION -> assert dimnames match (not just dims) to catch
# row/unit drift. Returns list(fit = contrasts.fit object, top = NAMED-by-contrast list of
# ordinary_t_tables) so top$interaction AND fit$coefficients[, "interaction"] both resolve
# downstream. Do NOT route through fit_limma_log (limma-TREND + eBayes, incoherent here). Pure.
fit_trajectory_contrasts <- function(measure_mat, design, contrasts, weights = NULL,
                                     conf_level = 0.95) {
  stopifnot(is.matrix(measure_mat), !is.null(rownames(measure_mat)),
            identical(colnames(measure_mat), rownames(design)),
            identical(rownames(contrasts), colnames(design)),
            qr(design)$rank == ncol(design),
            all(is.finite(measure_mat)))
  if (!is.null(weights)) {
    stopifnot(is.matrix(weights),
              identical(dimnames(weights), dimnames(measure_mat)),  # weights apply BY POSITION
              all(is.finite(weights)), all(weights > 0))
  }
  cfit <- limma::contrasts.fit(
    limma::lmFit(measure_mat, design = design, weights = weights), contrasts)
  # limma drops the feature rowname on a single-row fit -> restore the measure labels (1:1 with
  # the input rows; lmFit preserves feature count + order) so ordinary_t_table keys them.
  rownames(cfit$coefficients)    <- rownames(measure_mat)
  rownames(cfit$stdev.unscaled)  <- rownames(measure_mat)
  top <- stats::setNames(
    lapply(colnames(contrasts), function(cn) ordinary_t_table(cfit, cn, conf_level)),
    colnames(contrasts))
  list(fit = cfit, top = top)
}

# EXACT 3-channel Kitagawa / Oaxaca shift-share of the per-unit mean pseudotime:
#   mean_pt_r = comp_cf_r + prog_cf_r + cross_r - const, where (state s, unit r)
#   comp_cf = sum_s pi_sr * mu_bar_s  (vary composition, hold within-state means pooled)
#   prog_cf = sum_s pi_bar_s * mu_sr  (vary within-state means, hold composition pooled)
#   cross   = sum_s (pi_sr - pi_bar_s)(mu_sr - mu_bar_s),   const = sum_s pi_bar_s * mu_bar_s.
# Every broadcast (pi_bar/mu_bar down columns) is POSITION-based -> assert the state order is
# identical across pi/mu/pi_bar/mu_bar FIRST (a silent order drift would mis-split). Holds on the
# RAW additive pt scale only (logit/asin break additivity). Returns a per-unit data.frame. Pure.
kitagawa_channels <- function(pi_mat, mu_mat, pi_bar, mu_bar, tol = 1e-8) {
  stopifnot(is.matrix(pi_mat), is.matrix(mu_mat),
            identical(dimnames(pi_mat), dimnames(mu_mat)),
            identical(rownames(pi_mat), names(pi_bar)),
            identical(names(pi_bar), names(mu_bar)),
            all(is.finite(pi_mat)), all(is.finite(mu_mat)),
            all(is.finite(pi_bar)), all(is.finite(mu_bar)))
  const   <- sum(pi_bar * mu_bar)
  mean_pt <- colSums(pi_mat * mu_mat)
  comp_cf <- colSums(pi_mat * mu_bar)
  prog_cf <- colSums(pi_bar * mu_mat)
  cross   <- colSums((pi_mat - pi_bar) * (mu_mat - mu_bar))
  stopifnot(max(abs(mean_pt - (comp_cf + prog_cf + cross - const))) < tol)  # reconstruction
  data.frame(
    genotype_batch = colnames(pi_mat),
    mean_pt = mean_pt, comp_cf = comp_cf, prog_cf = prog_cf, cross = cross, const = const,
    row.names = NULL, stringsAsFactors = FALSE)
}

# Decompose the 2x2 interaction on mean pseudotime into its composition / progression / cross
# channels. The interaction contrast L() is LINEAR + intercept-free -> L(mean_pt) = L(comp_cf) +
# L(prog_cf) + L(cross) EXACTLY (const is unit-constant, annihilated). The exactness needs ONE
# shared per-unit weight vector replicated across the 4 channel-rows (the SAME WLS operator hits
# each row); differing per-row weights would break reconstruction. loadings = the 3 channel
# interaction coefs / L(mean_pt) (NA if |L(mean_pt)| < tol). Cell-weighted pooled anchors =
# PRIMARY; replicate-balanced (rowMeans) anchors = a SENSITIVITY. Pure.
decompose_progression_vs_composition <- function(per_rep, design, contrasts,
                                                 weights = NULL, conf_level = 0.95) {
  tol <- 1e-8
  channel_rows <- c("mean_pt", "comp_cf", "prog_cf", "cross")
  fit_channels <- function(pi_bar, mu_bar) {
    ch <- kitagawa_channels(per_rep$pi, per_rep$mu, pi_bar, mu_bar)
    M4 <- t(as.matrix(ch[, channel_rows]))
    dimnames(M4) <- list(channel_rows, ch$genotype_batch)
    stopifnot(identical(colnames(M4), rownames(design)))
    W4 <- NULL
    if (!is.null(weights)) {
      stopifnot(is.numeric(weights), !is.null(names(weights)),
                all(colnames(M4) %in% names(weights)))
      wv <- weights[colnames(M4)]                          # index by unit -> M4 column order
      stopifnot(all(is.finite(wv)), all(wv > 0))
      W4 <- matrix(wv, nrow(M4), ncol(M4), byrow = TRUE, dimnames = dimnames(M4))
    }
    fit   <- fit_trajectory_contrasts(M4, design, contrasts, weights = W4, conf_level = conf_level)
    L_int <- fit$fit$coefficients[, "interaction"]         # named over channel-rows
    recon <- unname(abs(L_int["mean_pt"] -
                        (L_int["comp_cf"] + L_int["prog_cf"] + L_int["cross"])))
    loadings <- if (abs(L_int["mean_pt"]) < tol)
      stats::setNames(rep(NA_real_, 3L), c("comp_cf", "prog_cf", "cross"))
    else L_int[c("comp_cf", "prog_cf", "cross")] / unname(L_int["mean_pt"])
    list(channels = ch, fit = fit, L_int = L_int, loadings = loadings, recon_resid_max = recon)
  }
  primary  <- fit_channels(per_rep$pi_bar, per_rep$mu_bar)              # cell-weighted anchors
  balanced <- fit_channels(rowMeans(per_rep$pi), rowMeans(per_rep$mu)) # replicate-balanced anchors
  stopifnot(primary$recon_resid_max < tol, balanced$recon_resid_max < tol)
  list(channels = primary$channels, fit = primary$fit, L_int = primary$L_int,
       loadings = primary$loadings, interaction = primary$fit$top$interaction,
       recon_resid_max = primary$recon_resid_max, balanced = balanced)
}
