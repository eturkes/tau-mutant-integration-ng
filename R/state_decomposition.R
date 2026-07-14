# P6 state composition-versus-regulation. S1 reads the heavy annotated Seurat once and emits raw
# Homeostatic/DAM pseudobulks + exact unit/score/map audits; S2 fits the compact substrate. No
# Seurat/S4 parent, fitted model, or cell-level frame may cross either target boundary.

# Sum raw feature x cell counts into state x replicate-unit pseudobulks with one sparse
# membership multiply. Output = named list of dense feature x unit matrices, state order fixed
# by `states`, unit order fixed by `units`.
aggregate_state_unit_counts <- function(counts, state, unit, states, units) {
  stopifnot(
    length(state) == ncol(counts), length(unit) == ncol(counts),
    length(states) == 2L, !anyDuplicated(states), all(nzchar(states)),
    length(units) >= 1L, !anyDuplicated(units), all(nzchar(units)),
    !is.null(rownames(counts)), !anyNA(rownames(counts)),
    !anyDuplicated(rownames(counts)), all(nzchar(rownames(counts)))
  )
  nonzero <- if (inherits(counts, "sparseMatrix")) counts@x else as.numeric(counts)
  stopifnot(all(is.finite(nonzero)), all(nonzero >= 0),
            all(abs(nonzero - round(nonzero)) <= sqrt(.Machine$double.eps)))

  state_i <- match(state, states)
  unit_i <- match(unit, units)
  primary <- !is.na(state_i)
  stopifnot(any(primary), !anyNA(unit_i), all(tabulate(state_i[primary], length(states)) > 0L))
  group_i <- unit_i[primary] + length(units) * (state_i[primary] - 1L)
  membership <- Matrix::sparseMatrix(
    i = which(primary), j = group_i, x = 1,
    dims = c(ncol(counts), length(states) * length(units))
  )
  aggregated <- as.matrix(counts %*% membership)
  dimnames(aggregated) <- list(
    rownames(counts),
    unlist(lapply(states, function(s) paste(s, units, sep = "::")), use.names = FALSE)
  )
  cell_library <- Matrix::colSums(counts)
  expected_library <- rowsum(
    matrix(cell_library[primary], ncol = 1L),
    factor(group_i, levels = seq_len(length(states) * length(units))),
    reorder = TRUE
  )[, 1L]
  stopifnot(all(is.finite(aggregated)), all(aggregated >= 0),
            isTRUE(all.equal(unname(colSums(aggregated)), unname(expected_library),
                            tolerance = 0)))

  out <- stats::setNames(lapply(seq_along(states), function(i) {
    idx <- seq_len(length(units)) + length(units) * (i - 1L)
    x <- aggregated[, idx, drop = FALSE]
    colnames(x) <- units
    x
  }), states)
  stopifnot(all(vapply(out, function(x) identical(rownames(x), rownames(counts)), logical(1))),
            all(vapply(out, function(x) identical(colnames(x), units), logical(1))))
  out
}

# Preserve every declared marker row, including absent raw-count mappings. A zero-row programme
# is explicit `testable = FALSE`; downstream rotation inference must not fabricate coverage.
state_marker_mapping <- function(marker_sets, feature_map) {
  stopifnot(is.list(marker_sets), length(marker_sets) >= 1L,
            !is.null(names(marker_sets)), !anyDuplicated(names(marker_sets)),
            is.data.frame(feature_map),
            all(c("ensembl", "symbol") %in% names(feature_map)),
            !anyNA(feature_map$ensembl), !anyDuplicated(feature_map$ensembl),
            all(nzchar(feature_map$ensembl)),
            !anyNA(feature_map$symbol), !anyDuplicated(feature_map$symbol),
            all(nzchar(feature_map$symbol)))
  rows <- lapply(names(marker_sets), function(program) {
    symbols <- marker_sets[[program]]
    stopifnot(is.character(symbols), length(symbols) >= 1L,
              !anyNA(symbols), !anyDuplicated(symbols), all(nzchar(symbols)))
    idx <- match(symbols, feature_map$symbol)
    data.frame(
      program = program, symbol = symbols,
      ensembl = feature_map$ensembl[idx], present = !is.na(idx),
      stringsAsFactors = FALSE
    )
  })
  mapping <- do.call(rbind, rows)
  rownames(mapping) <- NULL
  coverage <- do.call(rbind, lapply(names(marker_sets), function(program) {
    hit <- mapping$program == program
    data.frame(
      program = program, n_declared = sum(hit), n_present = sum(mapping$present[hit]),
      coverage = mean(mapping$present[hit]), testable = any(mapping$present[hit]),
      stringsAsFactors = FALSE
    )
  }))
  rownames(coverage) <- NULL
  stopifnot(nrow(mapping) == sum(lengths(marker_sets)),
            identical(coverage$program, names(marker_sets)),
            all(is.finite(coverage$coverage)),
            all(coverage$coverage >= 0 & coverage$coverage <= 1))
  list(mapping = mapping, coverage = coverage)
}

# Recursive isolation guard: target payloads remain inert base vectors/matrices/data frames/lists.
# Any heavy-parent, fitted-model, function, S4/environment, or external-pointer reachability fails.
state_substrate_contains_parent <- function(x) {
  forbidden <- c("Seurat", "Assay", "Assay5", "DimReduc", "Neighbor",
                 "glmmTMB", "lm", "glm", "MArrayLM", "EList", "DGEList")
  if (inherits(x, forbidden) || is.function(x) || isS4(x) || is.environment(x) ||
      typeof(x) == "externalptr") return(TRUE)
  if (!is.list(x)) return(FALSE)
  any(vapply(unclass(x), state_substrate_contains_parent, logical(1)))
}

# Extract + gate the fixed Homeostatic/DAM universe from microglia_annotated. Runtime assertions
# encode the S1 preconditions; elapsed time/peak RSS and qs bytes are measured around the fresh
# target build because embedding nondeterministic runtime metrics in the target would stale it.
build_microglia_state_substrate <- function(
    seurat_obj, symbol_map,
    states = c("Homeostatic", "DAM"),
    state_col = "microglia_subpopulation",
    unit_col = "genotype_batch", genotype_col = "genotype", batch_col = "batch",
    score_sets = canonical_microglia_markers,
    assay = "RNA", layer = "counts",
    expected_units = 16L, min_cells = 20L,
    min_overall_coverage = 0.95, min_unit_coverage = 0.90,
    max_in_memory_bytes = 25 * 1024^2,
    max_serialized_bytes = 25 * 1024^2) {
  programs <- names(score_sets)
  score_cols <- paste0(programs, "_UCell")
  required_md <- c(state_col, unit_col, genotype_col, batch_col, score_cols)
  stopifnot(
    inherits(seurat_obj, "Seurat"), assay %in% SeuratObject::Assays(seurat_obj),
    identical(states, c("Homeostatic", "DAM")),
    identical(score_sets, canonical_microglia_markers),
    length(programs) == 5L, all(required_md %in% colnames(seurat_obj@meta.data)),
    is.data.frame(symbol_map), all(c("ensembl", "symbol") %in% names(symbol_map)),
    length(expected_units) == 1L, expected_units == 16L,
    length(min_cells) == 1L, is.finite(min_cells), min_cells >= 20L,
    length(min_overall_coverage) == 1L, is.finite(min_overall_coverage),
    min_overall_coverage >= 0, min_overall_coverage <= 1,
    length(min_unit_coverage) == 1L, is.finite(min_unit_coverage),
    min_unit_coverage >= 0, min_unit_coverage <= 1,
    length(max_in_memory_bytes) == 1L, is.finite(max_in_memory_bytes),
    max_in_memory_bytes > 0,
    length(max_serialized_bytes) == 1L, is.finite(max_serialized_bytes),
    max_serialized_bytes > 0
  )
  md <- seurat_obj@meta.data
  counts <- SeuratObject::GetAssayData(seurat_obj, assay = assay, layer = layer)
  stopifnot(identical(rownames(md), colnames(counts)),
            all(vapply(md[score_cols], is.numeric, logical(1))))

  state <- as.character(md[[state_col]])
  unit <- as.character(md[[unit_col]])
  genotype <- as.character(md[[genotype_col]])
  batch <- as.character(md[[batch_col]])
  stopifnot(!anyNA(state), all(nzchar(state)), !anyNA(unit), all(nzchar(unit)),
            !anyNA(genotype), all(genotype %in% genotype_levels),
            !anyNA(batch), all(nzchar(batch)))

  # One row/unit; exact genotype x batch bijection; canonical factorial design with 9 residual df.
  unit_rows <- unique(data.frame(
    genotype_batch = unit, genotype = genotype, batch = batch,
    stringsAsFactors = FALSE
  ))
  stopifnot(nrow(unit_rows) == expected_units,
            !anyDuplicated(unit_rows$genotype_batch),
            !anyDuplicated(unit_rows[c("genotype", "batch")]))
  units <- sort(unit_rows$genotype_batch, method = "radix")
  unit_meta <- unit_rows[match(units, unit_rows$genotype_batch), , drop = FALSE]
  unit_meta$genotype <- factor(unit_meta$genotype, levels = genotype_levels)
  unit_meta$batch <- factor(unit_meta$batch,
                            levels = sort(unique(as.character(unit_meta$batch)), method = "radix"))
  rownames(unit_meta) <- units
  stopifnot(!anyNA(unit_meta$genotype), !anyNA(unit_meta$batch),
            identical(paste(as.character(unit_meta$genotype), as.character(unit_meta$batch),
                            sep = "_"), units))
  assert_complete_crossing(md, unit_col, c(genotype_col, batch_col))
  fd <- factorial_design(unit_meta)
  design_rank <- qr(fd$design)$rank
  stopifnot(design_rank == ncol(fd$design),
            nrow(fd$design) - design_rank == 9L,
            identical(colnames(fd$contrasts),
                      c("tau_alone", "nlgf_in_maptki", "nlgf_in_p301s",
                        "tau_in_nlgf", "interaction")))

  # Fixed two-state coverage + per-unit cell-count gates.
  unit_i <- match(unit, units)
  state_i <- match(state, states)
  primary <- !is.na(state_i)
  group_i <- unit_i[primary] + expected_units * (state_i[primary] - 1L)
  state_counts <- matrix(
    tabulate(group_i, nbins = expected_units * length(states)),
    nrow = expected_units, ncol = length(states), dimnames = list(units, states)
  )
  unit_total <- tabulate(unit_i, nbins = expected_units)
  unit_primary <- rowSums(state_counts)
  unit_coverage_value <- unit_primary / unit_total
  overall_coverage <- sum(primary) / nrow(md)
  stopifnot(
    all(unit_total > 0L), all(state_counts >= min_cells),
    all(colSums(state_counts > 0L) == expected_units),
    overall_coverage >= min_overall_coverage,
    all(unit_coverage_value >= min_unit_coverage)
  )
  unit_coverage <- data.frame(
    genotype_batch = units,
    genotype = unit_meta$genotype, batch = unit_meta$batch,
    n_retained = unit_total, n_primary = unit_primary,
    coverage = unit_coverage_value,
    n_Homeostatic = state_counts[, "Homeostatic"],
    n_DAM = state_counts[, "DAM"],
    DAM_fraction = state_counts[, "DAM"] / unit_primary,
    row.names = NULL, stringsAsFactors = FALSE
  )
  declared_states <- if (is.factor(md[[state_col]])) levels(md[[state_col]]) else character()
  all_states <- unique(c(declared_states, sort(unique(state), method = "radix")))
  state_summary <- data.frame(
    state = all_states,
    n_cells = as.integer(table(factor(state, levels = all_states))),
    primary = all_states %in% states,
    stringsAsFactors = FALSE
  )

  # Raw UCell values: finite + nonconstant in the pooled primary universe; compact unit/state means.
  score_mat <- as.matrix(md[primary, score_cols, drop = FALSE])
  colnames(score_mat) <- programs
  pooled_sd <- apply(score_mat, 2L, stats::sd)
  stopifnot(all(is.finite(score_mat)), all(is.finite(pooled_sd)), all(pooled_sd > 0))
  score_sums <- rowsum(
    score_mat,
    factor(group_i, levels = seq_len(expected_units * length(states))),
    reorder = TRUE
  )
  stopifnot(identical(rownames(score_sums),
                      as.character(seq_len(expected_units * length(states)))))
  score_mean_mat <- sweep(score_sums, 1L, as.numeric(state_counts), "/")
  score_scale <- data.frame(
    program = programs, pooled_sd = unname(pooled_sd),
    pooled_min = unname(apply(score_mat, 2L, min)),
    pooled_max = unname(apply(score_mat, 2L, max)),
    row.names = NULL, stringsAsFactors = FALSE
  )

  # Exact feature map + all declared marker rows for later fixed-programme rotations.
  feature_ids <- rownames(counts)
  map_i <- match(feature_ids, symbol_map$ensembl)
  stopifnot(!anyNA(map_i), !anyDuplicated(symbol_map$ensembl),
            !anyNA(symbol_map$symbol), !anyDuplicated(symbol_map$symbol),
            all(nzchar(symbol_map$symbol)))
  feature_map <- symbol_map[map_i, c("ensembl", "symbol"), drop = FALSE]
  rownames(feature_map) <- NULL
  stopifnot(identical(feature_map$ensembl, feature_ids))
  marker <- state_marker_mapping(score_sets, feature_map)

  # One aggregation pass; align libraries, cell counts, metadata, and score means on 32 rows.
  pseudobulk <- aggregate_state_unit_counts(counts, state, unit, states, units)
  library_size <- unlist(lapply(pseudobulk, colSums), use.names = FALSE)
  unit_state <- data.frame(
    state = factor(rep(states, each = expected_units), levels = states),
    genotype_batch = rep(units, times = length(states)),
    genotype = factor(rep(as.character(unit_meta$genotype), times = length(states)),
                      levels = genotype_levels),
    batch = factor(rep(as.character(unit_meta$batch), times = length(states)),
                   levels = levels(unit_meta$batch)),
    n_cells = as.numeric(state_counts), library_size = library_size,
    row.names = NULL, stringsAsFactors = FALSE
  )
  score_means <- cbind(unit_state[c("state", "genotype_batch", "genotype", "batch", "n_cells")],
                       as.data.frame(score_mean_mat, row.names = NULL, check.names = FALSE))
  stopifnot(all(is.finite(unit_state$library_size)), all(unit_state$library_size > 0),
            all(is.finite(as.matrix(score_means[programs]))),
            identical(as.character(score_means$state), as.character(unit_state$state)),
            identical(score_means$genotype_batch, unit_state$genotype_batch))

  out <- list(
    schema = "p6_state_substrate_v1",
    states = states,
    counts = pseudobulk,
    unit_meta = unit_meta,
    unit_state = unit_state,
    unit_coverage = unit_coverage,
    state_summary = state_summary,
    score_means = score_means,
    score_scale = score_scale,
    feature_map = feature_map,
    marker_map = marker$mapping,
    marker_coverage = marker$coverage,
    audit = list(
      n_retained = nrow(md), n_primary = sum(primary),
      overall_coverage = overall_coverage,
      min_unit_coverage_observed = min(unit_coverage_value),
      n_units = expected_units, min_state_unit_cells = min(state_counts),
      n_features = nrow(counts), design_rank = design_rank,
      design_columns = ncol(fd$design), residual_df = nrow(fd$design) - design_rank,
      contrast_names = colnames(fd$contrasts), score_programs = programs,
      assay = assay, layer = layer,
      thresholds = list(
        min_cells = min_cells, min_overall_coverage = min_overall_coverage,
        min_unit_coverage = min_unit_coverage,
        max_in_memory_bytes = max_in_memory_bytes,
        max_serialized_bytes = max_serialized_bytes
      ),
      parent_isolated = NA, in_memory_bytes = NA_real_
    )
  )
  stopifnot(!state_substrate_contains_parent(out))
  out$audit$parent_isolated <- TRUE
  out$audit$in_memory_bytes <- as.numeric(object.size(out))
  stopifnot(out$audit$in_memory_bytes <= max_in_memory_bytes,
            length(qs2::qs_serialize(out)) <= max_serialized_bytes)
  out
}

# ============================================================================================
# P6-S2: retained-nuclei occupancy + state-conditional raw-count response. Every inferential
# endpoint uses the 16 replicate units. Model/EList/MArrayLM objects stay local; the target emits
# only finite tables, compact diagnostics, and bridge summaries for S3/S4.

state_contrast_names <- function() {
  c("tau_alone", "nlgf_in_maptki", "nlgf_in_p301s", "tau_in_nlgf", "interaction")
}

state_genotype_contrasts <- function() {
  d <- diag(length(genotype_levels))
  dimnames(d) <- list(genotype_levels, genotype_levels)
  out <- make_contrast_matrix(d)
  stopifnot(identical(rownames(out), genotype_levels),
            identical(colnames(out), state_contrast_names()))
  out
}

# Muffle + retain messages, but convert every warning into a labelled build failure. targets only
# records warnings by default; S2's contract makes warning cleanliness an immediate model gate.
state_capture_clean <- function(expr, label) {
  warns <- character(); messages <- character()
  value <- withCallingHandlers(
    force(expr),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning")
    },
    message = function(m) {
      messages <<- c(messages, conditionMessage(m)); invokeRestart("muffleMessage")
    }
  )
  if (length(warns)) {
    stop(label, " emitted warning(s): ", paste(unique(warns), collapse = " | "), call. = FALSE)
  }
  list(value = value, messages = unique(trimws(messages[nzchar(trimws(messages))])))
}

# RNG-pure fixed stream for rotation/permutation inference.
state_with_seed <- function(seed, expr) {
  stopifnot(length(seed) == 1L, is.finite(seed), seed == round(seed))
  old_kind <- RNGkind()
  has_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (has_seed) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    RNGkind(old_kind[1L], old_kind[2L], old_kind[3L])
    if (has_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion",
           sample.kind = "Rejection")
  force(expr)
}

state_wald_contrasts <- function(beta, vcov, contrast_matrix, family, conf_level = 0.95) {
  stopifnot(is.numeric(beta), !is.null(names(beta)), is.matrix(vcov),
            identical(rownames(vcov), names(beta)), identical(colnames(vcov), names(beta)),
            is.matrix(contrast_matrix), identical(rownames(contrast_matrix), names(beta)),
            identical(colnames(contrast_matrix), state_contrast_names()),
            all(is.finite(beta)), all(is.finite(vcov)), all(is.finite(contrast_matrix)),
            length(family) == 1L, nzchar(family), conf_level > 0, conf_level < 1)
  estimate <- drop(crossprod(contrast_matrix, beta))
  contrast_vcov <- crossprod(contrast_matrix, vcov %*% contrast_matrix)
  se <- sqrt(pmax(0, diag(contrast_vcov)))
  stopifnot(all(is.finite(estimate)), all(is.finite(se)), all(se > 0))
  z <- estimate / se
  p <- 2 * stats::pnorm(-abs(z))
  q <- stats::qnorm(1 - (1 - conf_level) / 2)
  data.frame(
    contrast = state_contrast_names(), estimate = estimate, se = se, z = z,
    p_value = p, fdr = stats::p.adjust(p, "BH"),
    ci_l = estimate - q * se, ci_r = estimate + q * se,
    family = family, row.names = NULL, stringsAsFactors = FALSE
  )
}

# Equal-batch standardization on the response scale. Analytic gradients propagate the fixed-effect
# covariance through genotype means + all five contrasts; central differences guard the calculus.
state_probability_standardization <- function(beta, vcov, x_new, genotype,
                                               margin = 0.10, conf_level = 0.95) {
  stopifnot(is.numeric(beta), !is.null(names(beta)), is.matrix(vcov),
            identical(rownames(vcov), names(beta)), identical(colnames(vcov), names(beta)),
            is.matrix(x_new), identical(colnames(x_new), names(beta)),
            length(genotype) == nrow(x_new), all(genotype %in% genotype_levels),
            all(is.finite(beta)), all(is.finite(vcov)), all(is.finite(x_new)),
            length(margin) == 1L, is.finite(margin), margin > 0, margin < 1,
            conf_level > 0, conf_level < 1)
  groups <- lapply(genotype_levels, function(g) which(genotype == g))
  stopifnot(all(lengths(groups) == lengths(groups)[1L]), lengths(groups)[1L] >= 2L)
  mean_at <- function(b) vapply(groups, function(i) mean(stats::plogis(x_new[i, ] %*% b)),
                                numeric(1))
  estimate <- mean_at(beta)
  eta <- drop(x_new %*% beta); p <- stats::plogis(eta)
  gradient <- t(vapply(groups, function(i)
    colMeans(x_new[i, , drop = FALSE] * (p[i] * (1 - p[i]))), numeric(length(beta))))
  dimnames(gradient) <- list(genotype_levels, names(beta))
  numeric_gradient <- matrix(NA_real_, nrow(gradient), ncol(gradient), dimnames = dimnames(gradient))
  for (j in seq_along(beta)) {
    h <- 1e-6 * (1 + abs(beta[j]))
    bp <- bm <- beta; bp[j] <- bp[j] + h; bm[j] <- bm[j] - h
    numeric_gradient[, j] <- (mean_at(bp) - mean_at(bm)) / (2 * h)
  }
  gradient_residual <- max(abs(gradient - numeric_gradient))
  stopifnot(is.finite(gradient_residual), gradient_residual < 1e-7)

  mean_vcov <- gradient %*% vcov %*% t(gradient)
  mean_se <- sqrt(pmax(0, diag(mean_vcov)))
  stopifnot(all(is.finite(estimate)), all(estimate > 0 & estimate < 1),
            all(is.finite(mean_vcov)), all(is.finite(mean_se)), all(mean_se > 0))
  q <- stats::qnorm(1 - (1 - conf_level) / 2)
  means <- data.frame(
    genotype = genotype_levels, estimate = estimate, se = mean_se,
    ci_l = pmax(0, estimate - q * mean_se), ci_r = pmin(1, estimate + q * mean_se),
    row.names = NULL, stringsAsFactors = FALSE
  )

  cg <- state_genotype_contrasts()
  contrast_estimate <- drop(crossprod(cg, estimate))
  contrast_gradient <- crossprod(cg, gradient)
  contrast_vcov <- contrast_gradient %*% vcov %*% t(contrast_gradient)
  contrast_se <- sqrt(pmax(0, diag(contrast_vcov)))
  stopifnot(all(is.finite(contrast_estimate)), all(is.finite(contrast_se)), all(contrast_se > 0))
  z_zero <- contrast_estimate / contrast_se
  p_zero <- 2 * stats::pnorm(-abs(z_zero))
  beyond <- pmax(abs(contrast_estimate) - margin, 0)
  z_minimum <- beyond / contrast_se
  p_minimum <- ifelse(abs(contrast_estimate) <= margin, 1,
                      2 * stats::pnorm(-z_minimum))
  contrasts <- data.frame(
    contrast = state_contrast_names(), estimate = contrast_estimate, se = contrast_se,
    ci_l = contrast_estimate - q * contrast_se,
    ci_r = contrast_estimate + q * contrast_se,
    z_zero = z_zero, p_zero = p_zero, fdr_zero = stats::p.adjust(p_zero, "BH"),
    margin = margin, z_minimum = z_minimum, p_minimum = p_minimum,
    fdr_minimum = stats::p.adjust(p_minimum, "BH"),
    family_zero = "occupancy_probability_zero_all_contrasts",
    family_minimum = "occupancy_probability_minimum_all_contrasts",
    row.names = NULL, stringsAsFactors = FALSE
  )
  list(means = means, mean_vcov = mean_vcov, contrasts = contrasts,
       gradient_residual = gradient_residual)
}

# Freedman-Lane interaction sensitivity with residuals shuffled only within batch strata. The
# primary beta-binomial Wald tests remain load-bearing; this OLS permutation is finite-sample audit.
freedman_lane_stratified_interaction <- function(y, design, strata, int_col = "tau_nlgf",
                                                 n_perm = 9999L, seed = 614L) {
  stopifnot(is.numeric(y), is.matrix(design), length(y) == nrow(design),
            length(strata) == length(y), all(is.finite(y)), all(is.finite(design)),
            !anyNA(strata), int_col %in% colnames(design), qr(design)$rank == ncol(design),
            length(n_perm) == 1L, is.finite(n_perm), n_perm >= 1L, n_perm == round(n_perm))
  strata_i <- split(seq_along(y), factor(strata))
  stopifnot(length(strata_i) >= 2L, all(lengths(strata_i) >= 2L))
  xtx_inv <- chol2inv(chol(crossprod(design)))
  j <- match(int_col, colnames(design)); df <- nrow(design) - ncol(design)
  t_stat <- function(yv) {
    f <- stats::lm.fit(design, yv)
    sigma2 <- sum(f$residuals^2) / df
    unname(f$coefficients[j] / sqrt(sigma2 * xtx_inv[j, j]))
  }
  t_obs <- t_stat(y)
  reduced <- design[, colnames(design) != int_col, drop = FALSE]
  f0 <- stats::lm.fit(reduced, y)
  t_star <- state_with_seed(seed, vapply(seq_len(n_perm), function(i) {
    perm <- seq_along(y)
    for (ii in strata_i) perm[ii] <- sample(ii, length(ii), replace = FALSE)
    t_stat(f0$fitted.values + f0$residuals[perm])
  }, numeric(1)))
  stopifnot(is.finite(t_obs), all(is.finite(t_star)))
  list(t_obs = t_obs, n_perm = as.integer(n_perm), seed = as.integer(seed),
       perm_p = (1 + sum(abs(t_star) >= abs(t_obs))) / (n_perm + 1))
}

fit_state_occupancy <- function(substrate, margin = 0.10, n_perm = 9999L, seed = 614L,
                                max_abs_coef = 20, max_vcov_condition = 1e10) {
  d <- substrate$unit_coverage
  required <- c("genotype_batch", "genotype", "batch", "n_Homeostatic", "n_DAM",
                "n_primary", "coverage", "DAM_fraction")
  stopifnot(all(required %in% names(d)), nrow(d) == 16L,
            identical(d$genotype_batch, rownames(substrate$unit_meta)),
            d$n_primary == d$n_Homeostatic + d$n_DAM,
            isTRUE(all.equal(d$DAM_fraction, d$n_DAM / d$n_primary, tolerance = 0)))
  d$genotype <- factor(as.character(d$genotype), levels = genotype_levels)
  d$batch <- factor(as.character(d$batch),
                    levels = levels(substrate$unit_meta$batch))
  rownames(d) <- d$genotype_batch
  fd <- factorial_design(d)
  cap <- state_capture_clean(
    glmmTMB::glmmTMB(cbind(n_DAM, n_Homeostatic) ~ 0 + genotype + batch,
                     data = d, family = glmmTMB::betabinomial(link = "logit")),
    "beta-binomial DAM occupancy"
  )
  fit <- cap$value
  beta <- glmmTMB::fixef(fit)$cond
  vc <- stats::vcov(fit)$cond
  se_beta <- sqrt(diag(vc))
  dispersion <- stats::sigma(fit)
  vcov_condition <- kappa(vc)
  stopifnot(
    fit$fit$convergence == 0L, isTRUE(fit$sdr$pdHess),
    all(is.finite(beta)), all(is.finite(se_beta)), all(se_beta > 0),
    max(abs(beta)) <= max_abs_coef,
    all(is.finite(vc)), min(eigen(vc, symmetric = TRUE, only.values = TRUE)$values) > 0,
    is.finite(vcov_condition), vcov_condition <= max_vcov_condition,
    is.finite(dispersion), dispersion > 0
  )
  model_x <- stats::model.matrix(~ 0 + genotype + batch, d)
  stopifnot(identical(colnames(model_x), names(beta)), qr(model_x)$rank == ncol(model_x),
            max(abs(drop(model_x %*% beta) - stats::predict(fit, type = "link"))) < 1e-8)

  cm <- matrix(0, nrow = length(beta), ncol = length(state_contrast_names()),
               dimnames = list(names(beta), state_contrast_names()))
  cm[paste0("genotype", genotype_levels), ] <- state_genotype_contrasts()
  log_odds <- state_wald_contrasts(beta, vc, cm,
                                   family = "occupancy_logodds_zero_all_contrasts")

  grid <- expand.grid(genotype = genotype_levels, batch = levels(d$batch),
                      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  grid$genotype <- factor(grid$genotype, levels = genotype_levels)
  grid$batch <- factor(grid$batch, levels = levels(d$batch))
  x_new <- stats::model.matrix(~ 0 + genotype + batch, grid)
  stopifnot(identical(colnames(x_new), names(beta)))
  probability <- state_probability_standardization(beta, vc, x_new,
                                                    as.character(grid$genotype), margin)

  empirical_logit <- log((d$n_DAM + 0.5) / (d$n_Homeostatic + 0.5))
  emat <- matrix(empirical_logit, nrow = 1L,
                 dimnames = list("empirical_logit", rownames(fd$design)))
  ols_cap <- state_capture_clean(
    fit_trajectory_contrasts(emat, fd$design, fd$contrasts),
    "empirical-logit occupancy OLS"
  )
  ols <- do.call(rbind, ols_cap$value$top)
  rownames(ols) <- NULL
  ols$fdr <- stats::p.adjust(ols$p_value, "BH")
  ols$family <- "occupancy_empirical_logit_zero_all_contrasts"
  fl <- freedman_lane_stratified_interaction(
    empirical_logit, fd$design, d$batch, n_perm = n_perm, seed = seed
  )
  stopifnot(abs(fl$t_obs - ols$t[ols$contrast == "interaction"]) < 1e-8)

  list(
    unit = d[required], log_odds = log_odds,
    probability_means = probability$means,
    probability_vcov = probability$mean_vcov,
    probability_contrasts = probability$contrasts,
    empirical_logit = ols,
    permutation = data.frame(
      contrast = "interaction", t_obs = fl$t_obs, n_perm = fl$n_perm,
      perm_p = fl$perm_p, seed = fl$seed,
      family = "occupancy_empirical_logit_interaction_permutation",
      row.names = NULL, stringsAsFactors = FALSE
    ),
    diagnostics = list(
      convergence = fit$fit$convergence, pdHess = fit$sdr$pdHess,
      dispersion = unname(dispersion), max_abs_coef = max(abs(beta)),
      max_abs_coef_gate = max_abs_coef, vcov_condition = unname(vcov_condition),
      vcov_condition_gate = max_vcov_condition,
      gradient_residual = probability$gradient_residual,
      n_units = nrow(d), design_rank = qr(fd$design)$rank,
      residual_df = nrow(fd$design) - ncol(fd$design), messages = cap$messages
    )
  )
}

state_de_table <- function(efit, tfit, endpoint, contrast_names = state_contrast_names()) {
  stopifnot(inherits(efit, "MArrayLM"), inherits(tfit, "MArrayLM"),
            identical(colnames(efit$coefficients), contrast_names),
            identical(colnames(tfit$coefficients), contrast_names),
            identical(rownames(efit$coefficients), rownames(tfit$coefficients)))
  rows <- lapply(contrast_names, function(cn) {
    raw <- limma::topTable(efit, coef = cn, number = Inf, sort.by = "none", confint = TRUE)
    threshold <- limma::topTreat(tfit, coef = cn, number = Inf, sort.by = "none")
    stopifnot(identical(rownames(raw), rownames(threshold)),
              isTRUE(all.equal(raw$logFC, threshold$logFC, tolerance = 1e-12)))
    p <- raw$P.Value; tp <- threshold$P.Value
    data.frame(
      endpoint = endpoint, contrast = cn, gene = rownames(raw),
      logFC = raw$logFC, ci_l = raw$CI.L, ci_r = raw$CI.R,
      ave_expr = raw$AveExpr, t = raw$t, p_value = p, fdr = stats::p.adjust(p, "BH"),
      treat_t = threshold$t, treat_p = tp, treat_fdr = stats::p.adjust(tp, "BH"),
      family = paste0("gene_zero_", endpoint, "_", cn),
      treat_family = paste0("gene_minimum_", endpoint, "_", cn),
      row.names = NULL, stringsAsFactors = FALSE, check.names = FALSE
    )
  })
  out <- do.call(rbind, rows); rownames(out) <- NULL
  numeric_cols <- c("logFC", "ci_l", "ci_r", "ave_expr", "t", "p_value", "fdr",
                    "treat_t", "treat_p", "treat_fdr")
  stopifnot(!anyDuplicated(out[c("endpoint", "contrast", "gene")]),
            all(vapply(out[numeric_cols], function(x) all(is.finite(x)), logical(1))),
            all(out$p_value >= 0 & out$p_value <= 1), all(out$fdr >= 0 & out$fdr <= 1),
            all(out$treat_p >= 0 & out$treat_p <= 1),
            all(out$treat_fdr >= 0 & out$treat_fdr <= 1))
  out
}

fit_state_expression <- function(expression, design, contrasts, endpoint, weights = NULL,
                                 lfc = 0.5) {
  stopifnot(is.matrix(expression), !is.null(rownames(expression)),
            !anyDuplicated(rownames(expression)), identical(colnames(expression), rownames(design)),
            identical(rownames(contrasts), colnames(design)), all(is.finite(expression)),
            qr(design)$rank == ncol(design), length(lfc) == 1L, is.finite(lfc), lfc > 0)
  if (!is.null(weights)) {
    stopifnot(is.matrix(weights), identical(dimnames(weights), dimnames(expression)),
              all(is.finite(weights)), all(weights > 0))
  }
  cap <- state_capture_clean({
    base <- limma::lmFit(expression, design = design, weights = weights)
    contrast_fit <- limma::contrasts.fit(base, contrasts)
    efit <- limma::eBayes(contrast_fit, robust = TRUE)
    tfit <- limma::treat(contrast_fit, lfc = lfc, robust = TRUE)
    list(contrast_fit = contrast_fit, efit = efit, tfit = tfit)
  }, paste0(endpoint, " weighted expression fit"))
  fit <- cap$value
  table <- state_de_table(fit$efit, fit$tfit, endpoint)
  stopifnot(all(fit$contrast_fit$df.residual == 9L))
  list(table = table, contrast_fit = fit$contrast_fit, messages = cap$messages)
}

fit_state_counts <- function(counts, fd, endpoint, min_count = 5, lfc = 0.5) {
  stopifnot(is.matrix(counts), identical(colnames(counts), rownames(fd$design)),
            identical(rownames(fd$contrasts), colnames(fd$design)),
            all(is.finite(counts)), all(counts >= 0))
  cap <- state_capture_clean({
    dge <- edgeR::DGEList(counts = counts)
    keep <- edgeR::filterByExpr(dge, design = fd$design, min.count = min_count)
    stopifnot(sum(keep) >= 100L)
    dge <- edgeR::normLibSizes(dge[keep, , keep.lib.sizes = FALSE], method = "TMM")
    voom <- limma::voomWithQualityWeights(dge, design = fd$design, plot = FALSE)
    list(voom = voom, keep = keep, norm_factors = dge$samples$norm.factors)
  }, paste0(endpoint, " voomWithQualityWeights"))
  v <- cap$value$voom
  # voom weights are position-aligned but intentionally un-dimnamed; pin the alignment before any
  # paired arithmetic/rotation so every later weight operation can assert exact feature/unit keys.
  dimnames(v$weights) <- dimnames(v$E)
  stopifnot(identical(colnames(v$E), rownames(fd$design)),
            identical(dimnames(v$weights), dimnames(v$E)),
            all(is.finite(v$E)), all(is.finite(v$weights)), all(v$weights > 0),
            all(is.finite(v$targets$sample.weights)), all(v$targets$sample.weights > 0))
  fitted <- fit_state_expression(v$E, fd$design, fd$contrasts, endpoint,
                                 weights = v$weights, lfc = lfc)
  list(
    voom = v, table = fitted$table, contrast_fit = fitted$contrast_fit,
    audit = list(
      endpoint = endpoint, n_input = nrow(counts), n_kept = nrow(v$E),
      min_count = min_count, lfc_margin = lfc,
      raw_library_size = colSums(counts), norm_factors = cap$value$norm_factors,
      sample_weights = v$targets$sample.weights,
      weight_range = range(v$weights), residual_df = 9L
    ),
    messages = unique(c(cap$messages, fitted$messages))
  )
}

# Paired state response: primary precision is the harmonic combination of the two voom weight
# matrices; an unweighted delta fit is retained as the fixed sensitivity.
fit_paired_state_difference <- function(homeostatic, dam, fd, lfc = 0.5) {
  hgenes <- rownames(homeostatic$voom$E); dgenes <- rownames(dam$voom$E)
  common <- hgenes[hgenes %in% dgenes]
  stopifnot(length(common) >= 100L, !anyDuplicated(common))
  hi <- match(common, hgenes); di <- match(common, dgenes)
  h <- homeostatic$voom$E[hi, , drop = FALSE]
  d <- dam$voom$E[di, , drop = FALSE]
  wh <- homeostatic$voom$weights[hi, , drop = FALSE]
  wd <- dam$voom$weights[di, , drop = FALSE]
  delta <- d - h
  weights <- 1 / (1 / wd + 1 / wh)
  dimnames(weights) <- dimnames(delta)
  primary <- fit_state_expression(delta, fd$design, fd$contrasts,
                                  "DAM_minus_Homeostatic", weights, lfc)
  sensitivity <- fit_state_expression(delta, fd$design, fd$contrasts,
                                      "DAM_minus_Homeostatic_unweighted", NULL, lfc)
  delta_voom <- dam$voom[di, ]
  delta_voom$E <- delta; delta_voom$weights <- weights
  stopifnot(identical(dimnames(delta_voom$E), dimnames(delta_voom$weights)))
  list(
    voom = delta_voom, table = primary$table, sensitivity_table = sensitivity$table,
    contrast_fit = primary$contrast_fit,
    audit = list(n_common = length(common), weight_range = range(weights),
                 residual_df = 9L, primary_weight = "1/(1/w_DAM+1/w_Homeostatic)",
                 sensitivity = "unweighted paired delta"),
    messages = unique(c(primary$messages, sensitivity$messages))
  )
}

state_rotation_table <- function(voom, contrast_fit, marker_map, endpoint, contrasts,
                                 nrot = 9999L, seed = 614L) {
  programs <- names(canonical_microglia_markers)
  stopifnot(inherits(voom, "EList"), inherits(contrast_fit, "MArrayLM"),
            is.data.frame(marker_map), all(c("program", "ensembl", "present") %in% names(marker_map)),
            identical(colnames(contrast_fit$coefficients), state_contrast_names()),
            identical(rownames(contrast_fit$coefficients), rownames(voom$E)),
            is.matrix(contrasts), identical(rownames(contrasts), colnames(voom$design)),
            identical(colnames(contrasts), state_contrast_names()),
            all(programs %in% marker_map$program), nrot == 9999L)
  index <- stats::setNames(lapply(programs, function(program) {
    ids <- marker_map$ensembl[marker_map$program == program & marker_map$present]
    i <- match(ids, rownames(voom$E), nomatch = 0L)
    unique(i[i > 0L])
  }), programs)
  rows <- lapply(state_contrast_names(), function(cn) {
    testable <- lengths(index) > 0L
    result <- NULL; messages <- character()
    if (any(testable)) {
      cap <- state_capture_clean(
        state_with_seed(seed, limma::mroast(
          voom, index = index[testable], design = voom$design,
          contrast = contrasts[, cn],
          set.statistic = "mean", nrot = nrot, approx.zscore = TRUE,
          adjust.method = "BH", midp = TRUE, sort = "none"
        )),
        paste0(endpoint, " mroast ", cn)
      )
      result <- cap$value[programs[testable], , drop = FALSE]
      messages <- cap$messages
    }
    direction <- rep("untestable", length(programs))
    if (any(testable)) direction[testable] <- result$Direction
    out <- data.frame(
      endpoint = endpoint, contrast = cn, program = programs,
      n_genes = lengths(index), testable = testable,
      mean_logFC = vapply(index, function(i) if (length(i))
        mean(contrast_fit$coefficients[i, cn]) else NA_real_, numeric(1)),
      direction = direction,
      p_value = NA_real_, fdr = NA_real_,
      nrot = as.integer(nrot), seed = as.integer(seed),
      family = paste0("rotation_", endpoint, "_", cn, "_programmes"),
      messages = paste(messages, collapse = " | "),
      row.names = NULL, stringsAsFactors = FALSE
    )
    if (any(testable)) {
      out$p_value[testable] <- result$PValue
      out$fdr[testable] <- stats::p.adjust(result$PValue, "BH")
    }
    out
  })
  out <- do.call(rbind, rows); rownames(out) <- NULL
  finite <- out$testable
  stopifnot(nrow(out) == length(programs) * length(state_contrast_names()),
            all(is.finite(out$mean_logFC[finite])),
            all(is.finite(out$p_value[finite])), all(is.finite(out$fdr[finite])),
            all(out$direction[finite] %in% c("Up", "Down")),
            all(out$p_value[finite] >= 0 & out$p_value[finite] <= 1),
            all(out$fdr[finite] >= 0 & out$fdr[finite] <= 1),
            all(is.na(out$p_value[!finite])), all(is.na(out$fdr[!finite])))
  out
}

state_safe_spearman <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  if (length(x) < 3L || length(unique(x)) < 2L || length(unique(y)) < 2L) return(NA_real_)
  unname(stats::cor(x, y, method = "spearman"))
}

# Compare pooled Homeostatic+DAM effects to the established whole-microglia target. This is a
# descriptive exclusion bridge: every contrast/program is carried, with no pass threshold.
state_pooled_bridge <- function(pooled_table, whole_pb, marker_map) {
  stopifnot(is.data.frame(pooled_table), is.list(whole_pb), is.list(whole_pb$top),
            identical(names(whole_pb$top), state_contrast_names()))
  programs <- names(canonical_microglia_markers)
  summaries <- lapply(state_contrast_names(), function(cn) {
    p <- pooled_table[pooled_table$contrast == cn, c("gene", "logFC"), drop = FALSE]
    w <- whole_pb$top[[cn]][, c("gene", "logFC"), drop = FALSE]
    common <- p$gene[p$gene %in% w$gene]
    pi <- match(common, p$gene); wi <- match(common, w$gene)
    x <- p$logFC[pi]; y <- w$logFC[wi]
    stopifnot(length(common) >= 100L, all(is.finite(x)), all(is.finite(y)))
    data.frame(
      contrast = cn, n_pooled = nrow(p), n_whole = nrow(w), n_common = length(common),
      spearman = state_safe_spearman(x, y),
      direction_concordance = mean(sign(x) == sign(y)),
      mean_abs_logFC_delta = mean(abs(x - y)),
      row.names = NULL, stringsAsFactors = FALSE
    )
  })
  summary <- do.call(rbind, summaries); rownames(summary) <- NULL
  stopifnot(all(is.finite(as.matrix(summary[c("n_pooled", "n_whole", "n_common", "spearman",
                                                     "direction_concordance", "mean_abs_logFC_delta")]))))

  programme_rows <- unlist(lapply(state_contrast_names(), function(cn) {
    p <- pooled_table[pooled_table$contrast == cn, c("gene", "logFC"), drop = FALSE]
    w <- whole_pb$top[[cn]][, c("gene", "logFC"), drop = FALSE]
    lapply(programs, function(program) {
      declared <- unique(marker_map$ensembl[marker_map$program == program & marker_map$present])
      common <- declared[declared %in% p$gene & declared %in% w$gene]
      x <- p$logFC[match(common, p$gene)]; y <- w$logFC[match(common, w$gene)]
      stopifnot(length(common) >= 1L, all(is.finite(x)), all(is.finite(y)))
      rho <- state_safe_spearman(x, y)
      data.frame(
        contrast = cn, program = program, n_genes = length(common),
        pooled_mean_logFC = mean(x), whole_mean_logFC = mean(y),
        mean_direction_concordant = sign(mean(x)) == sign(mean(y)),
        gene_spearman = rho, correlation_testable = is.finite(rho),
        mean_abs_logFC_delta = mean(abs(x - y)),
        row.names = NULL, stringsAsFactors = FALSE
      )
    })
  }), recursive = FALSE)
  programme <- do.call(rbind, programme_rows); rownames(programme) <- NULL
  stopifnot(nrow(programme) == length(programs) * length(state_contrast_names()),
            all(is.finite(programme$n_genes)),
            all(is.finite(programme$pooled_mean_logFC)),
            all(is.finite(programme$whole_mean_logFC)),
            all(is.finite(programme$mean_abs_logFC_delta)),
            all(is.finite(programme$gene_spearman[programme$correlation_testable])),
            all(is.na(programme$gene_spearman[!programme$correlation_testable])))
  list(summary = summary, programme = programme,
       note = "Pooled Homeostatic+DAM versus all retained microglia; descriptive, no agreement gate.")
}

# Deterministic algebra oracle: paired delta coefficients must equal DAM minus Homeostatic when
# all three fits share weights; canonical and probability interaction signs/gradients are checked.
state_response_synthetic_checks <- function(fd, unit_meta) {
  x <- fd$design; contrasts <- fd$contrasts
  stopifnot(nrow(x) == 16L, ncol(x) == 7L,
            identical(colnames(contrasts), state_contrast_names()))
  ng <- 3L
  bh <- matrix(0, ng, ncol(x), dimnames = list(paste0("g", seq_len(ng)), colnames(x)))
  bd <- bh
  bh[, c("(Intercept)", "tau", "nlgf", "tau_nlgf")] <- rbind(
    c(3, 0.1, 0.2, 0.1), c(4, -0.2, 0.1, -0.1), c(2, 0.3, -0.1, 0.2))
  bd[, c("(Intercept)", "tau", "nlgf", "tau_nlgf")] <- rbind(
    c(3.5, 0.2, 0.4, 2.0), c(3.7, -0.1, 0.3, -0.4), c(2.3, 0.1, 0.2, 0.7))
  noise_h <- outer(seq_len(ng), seq_len(nrow(x)), function(i, j) sin(i + 0.7 * j) / 20)
  noise_d <- outer(seq_len(ng), seq_len(nrow(x)), function(i, j) cos(0.4 * i + j) / 20)
  yh <- bh %*% t(x) + noise_h; yd <- bd %*% t(x) + noise_d
  weights <- 1 + abs(outer(seq_len(ng), seq_len(nrow(x)), function(i, j) sin(i * j)))
  dimnames(yh) <- dimnames(yd) <- dimnames(weights) <- list(rownames(bh), rownames(x))
  coef_for <- function(y) limma::contrasts.fit(
    limma::lmFit(y, design = x, weights = weights), contrasts)$coefficients
  ch <- coef_for(yh); cd <- coef_for(yd); delta <- coef_for(yd - yh)
  paired_residual <- max(abs(delta - (cd - ch)))
  stopifnot(paired_residual < 1e-10, delta["g1", "interaction"] > 0)

  grid <- expand.grid(genotype = genotype_levels, batch = levels(unit_meta$batch),
                      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  grid$genotype <- factor(grid$genotype, levels = genotype_levels)
  grid$batch <- factor(grid$batch, levels = levels(unit_meta$batch))
  x_new <- stats::model.matrix(~ 0 + genotype + batch, grid)
  beta <- c(-1.2, -0.8, -0.1, 1.0, 0.15, -0.10, 0.05)
  names(beta) <- colnames(x_new)
  vc <- diag(rep(0.02, length(beta)))
  dimnames(vc) <- list(names(beta), names(beta))
  probability <- state_probability_standardization(beta, vc, x_new,
                                                    as.character(grid$genotype), 0.10)
  pint <- probability$contrasts$estimate[probability$contrasts$contrast == "interaction"]
  stopifnot(pint > 0)
  canonical <- drop(crossprod(state_genotype_contrasts(), c(0.1, 0.2, 0.3, 0.7)))
  stopifnot(canonical["interaction"] > 0)
  list(paired_contrast_residual = paired_residual,
       probability_gradient_residual = probability$gradient_residual,
       paired_interaction = unname(delta["g1", "interaction"]),
       probability_interaction = unname(pint), canonical_interaction = canonical["interaction"])
}

run_microglia_state_response <- function(substrate, pb_de_microglia,
                                         min_count = 5, gene_lfc = 0.5,
                                         occupancy_margin = 0.10,
                                         nrot = 9999L, n_perm = 9999L, seed = 614L,
                                         max_target_bytes = 60 * 1024^2) {
  contrasts <- state_contrast_names(); states <- c("Homeostatic", "DAM")
  stopifnot(
    is.list(substrate), identical(substrate$schema, "p6_state_substrate_v1"),
    identical(substrate$states, states), identical(names(substrate$counts), states),
    identical(substrate$audit$contrast_names, contrasts),
    substrate$audit$n_units == 16L, substrate$audit$residual_df == 9L,
    isTRUE(substrate$audit$parent_isolated),
    all(vapply(substrate$counts, is.matrix, logical(1))),
    all(vapply(substrate$counts, function(x)
      identical(colnames(x), rownames(substrate$unit_meta)), logical(1))),
    is.list(pb_de_microglia), identical(names(pb_de_microglia$top), contrasts),
    nrot == 9999L, n_perm >= 999L, max_target_bytes > 0
  )
  fd <- factorial_design(substrate$unit_meta)
  stopifnot(identical(colnames(fd$contrasts), contrasts),
            nrow(fd$design) - ncol(fd$design) == 9L)
  synthetic <- state_response_synthetic_checks(fd, substrate$unit_meta)
  occupancy <- fit_state_occupancy(substrate, occupancy_margin, n_perm, seed)

  fits <- stats::setNames(lapply(states, function(state)
    fit_state_counts(substrate$counts[[state]], fd, state, min_count, gene_lfc)), states)
  paired <- fit_paired_state_difference(fits$Homeostatic, fits$DAM, fd, gene_lfc)
  rotations <- do.call(rbind, c(
    lapply(states, function(state)
      state_rotation_table(fits[[state]]$voom, fits[[state]]$contrast_fit,
                           substrate$marker_map, state, fd$contrasts, nrot, seed)),
    list(state_rotation_table(paired$voom, paired$contrast_fit,
                              substrate$marker_map, "DAM_minus_Homeostatic",
                              fd$contrasts, nrot, seed))
  ))
  rownames(rotations) <- NULL
  marker_coverage <- unique(rotations[c("endpoint", "program", "n_genes", "testable")])
  declared_i <- match(marker_coverage$program, substrate$marker_coverage$program)
  stopifnot(!anyNA(declared_i))
  marker_coverage$n_declared <- substrate$marker_coverage$n_declared[declared_i]
  marker_coverage$n_mapped <- substrate$marker_coverage$n_present[declared_i]
  rownames(marker_coverage) <- NULL

  pooled <- fit_state_counts(substrate$counts$Homeostatic + substrate$counts$DAM,
                             fd, "Pooled_primary_states", min_count, gene_lfc)
  bridge <- state_pooled_bridge(pooled$table, pb_de_microglia, substrate$marker_map)

  state_gene <- do.call(rbind, lapply(states, function(state) fits[[state]]$table))
  rownames(state_gene) <- NULL
  out <- list(
    schema = "p6_state_response_v1",
    occupancy = occupancy,
    gene = list(
      state = state_gene,
      state_difference = paired$table,
      state_difference_unweighted = paired$sensitivity_table
    ),
    rotations = rotations,
    marker_coverage = marker_coverage,
    bridge = bridge,
    audit = list(
      n_units = 16L, residual_df = 9L, contrast_names = contrasts,
      states = states, programmes = names(canonical_microglia_markers),
      state_fits = lapply(fits, `[[`, "audit"), paired = paired$audit,
      pooled = pooled$audit, synthetic = synthetic,
      messages = list(
        occupancy = occupancy$diagnostics$messages,
        state = lapply(fits, `[[`, "messages"), paired = paired$messages,
        pooled = pooled$messages
      ),
      thresholds = list(
        min_count = min_count, gene_lfc = gene_lfc,
        occupancy_margin = occupancy_margin, nrot = as.integer(nrot),
        n_perm = as.integer(n_perm), seed = as.integer(seed),
        max_target_bytes = max_target_bytes
      ),
      versions = list(
        r = as.character(getRversion()),
        glmmTMB = as.character(utils::packageVersion("glmmTMB")),
        edgeR = as.character(utils::packageVersion("edgeR")),
        limma = as.character(utils::packageVersion("limma"))
      ),
      parent_isolated = NA, in_memory_bytes = NA_real_
    )
  )

  finite_table <- function(x, cols) all(vapply(x[cols], function(z) all(is.finite(z)), logical(1)))
  stopifnot(
    nrow(out$occupancy$log_odds) == 5L,
    nrow(out$occupancy$probability_contrasts) == 5L,
    nrow(out$occupancy$empirical_logit) == 5L,
    identical(out$occupancy$log_odds$contrast, contrasts),
    identical(out$occupancy$probability_contrasts$contrast, contrasts),
    finite_table(out$occupancy$log_odds,
                 c("estimate", "se", "z", "p_value", "fdr", "ci_l", "ci_r")),
    finite_table(out$occupancy$probability_contrasts,
                 c("estimate", "se", "ci_l", "ci_r", "p_zero", "fdr_zero",
                   "p_minimum", "fdr_minimum")),
    finite_table(out$occupancy$probability_means,
                 c("estimate", "se", "ci_l", "ci_r")),
    all(is.finite(out$occupancy$probability_vcov)),
    finite_table(out$occupancy$empirical_logit,
                 c("coef", "se", "t", "df", "p_value", "fdr", "ci_l", "ci_r")),
    all(vapply(out$gene, function(x) finite_table(
      x, c("logFC", "ci_l", "ci_r", "p_value", "fdr", "treat_p", "treat_fdr")), logical(1))),
    all(vapply(out$gene, function(x) identical(unique(x$contrast), contrasts), logical(1))),
    identical(unique(out$gene$state$endpoint), states),
    identical(unique(out$gene$state_difference$endpoint), "DAM_minus_Homeostatic"),
    identical(unique(out$gene$state_difference_unweighted$endpoint),
              "DAM_minus_Homeostatic_unweighted"),
    all(vapply(out$gene, function(x)
      all(nzchar(x$family)) && all(nzchar(x$treat_family)), logical(1))),
    nrow(out$rotations) == 3L * 5L * 5L,
    all(table(out$rotations$endpoint, out$rotations$contrast) == 5L),
    nrow(out$marker_coverage) == 3L * 5L,
    !anyDuplicated(out$marker_coverage[c("endpoint", "program")]),
    all(out$marker_coverage$n_genes <= out$marker_coverage$n_mapped),
    nrow(out$bridge$summary) == 5L, nrow(out$bridge$programme) == 25L,
    !state_substrate_contains_parent(out)
  )
  out$audit$parent_isolated <- TRUE
  out$audit$in_memory_bytes <- as.numeric(object.size(out))
  stopifnot(out$audit$in_memory_bytes <= max_target_bytes,
            length(qs2::qs_serialize(out)) <= max_target_bytes)
  out
}
