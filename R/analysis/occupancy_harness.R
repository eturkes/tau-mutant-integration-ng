# P7.4 preregistration-only DAM-occupancy harness. This module maps a supplied membership
# labeling to the frozen occupancy family; it does not generate or select labeling variants.

occupancy_harness_empty_aggregation <- function(reason) {
  unit_meta <- data.frame(
    genotype_batch = character(),
    genotype = factor(character(), levels = genotype_levels),
    batch = factor(character()),
    stringsAsFactors = FALSE
  )
  unit_coverage <- data.frame(
    genotype_batch = character(),
    genotype = factor(character(), levels = genotype_levels),
    batch = factor(character()),
    n_retained = integer(), n_primary = numeric(), coverage = numeric(),
    n_Homeostatic = integer(), n_DAM = integer(), DAM_fraction = numeric(),
    stringsAsFactors = FALSE
  )
  list(
    unit_coverage = unit_coverage,
    unit_meta = unit_meta,
    diagnostics = list(
      status = "non_estimable", reason = as.character(reason), n_units = 0L,
      design_rank = NA_integer_, design_columns = NA_integer_, residual_df = NA_integer_,
      missing_state_units = character()
    )
  )
}

membership_to_unit_coverage <- function(
    membership, meta, states = c("Homeostatic", "DAM"),
    unit_col = "genotype_batch", genotype_col = "genotype", batch_col = "batch") {
  required <- c(unit_col, genotype_col, batch_col)
  if (!is.data.frame(meta)) {
    return(occupancy_harness_empty_aggregation("meta must be a data frame"))
  }
  if (length(membership) != nrow(meta)) {
    return(occupancy_harness_empty_aggregation(
      "membership length must equal the number of metadata rows"
    ))
  }
  if (!identical(states, c("Homeostatic", "DAM"))) {
    return(occupancy_harness_empty_aggregation(
      "states must be exactly c('Homeostatic', 'DAM')"
    ))
  }
  if (!all(required %in% names(meta))) {
    return(occupancy_harness_empty_aggregation(
      paste("metadata is missing required column(s):",
            paste(setdiff(required, names(meta)), collapse = ", "))
    ))
  }

  state <- as.character(membership)
  unit <- as.character(meta[[unit_col]])
  genotype <- as.character(meta[[genotype_col]])
  batch <- as.character(meta[[batch_col]])
  if (anyNA(state) || any(!nzchar(state)) || anyNA(unit) || any(!nzchar(unit)) ||
      anyNA(genotype) || any(!nzchar(genotype)) || anyNA(batch) || any(!nzchar(batch))) {
    return(occupancy_harness_empty_aggregation(
      "membership and unit/genotype/batch metadata must be non-missing, non-empty values"
    ))
  }
  if (!all(genotype %in% genotype_levels)) {
    return(occupancy_harness_empty_aggregation(
      "all genotypes must belong to the canonical genotype_levels"
    ))
  }

  unit_rows <- unique(data.frame(
    genotype_batch = unit, genotype = genotype, batch = batch,
    stringsAsFactors = FALSE
  ))
  if (anyDuplicated(unit_rows$genotype_batch)) {
    return(occupancy_harness_empty_aggregation(
      "each unit must map to exactly one genotype and batch"
    ))
  }
  if (anyDuplicated(unit_rows[c("genotype", "batch")])) {
    return(occupancy_harness_empty_aggregation(
      "each genotype-by-batch combination must map to at most one unit"
    ))
  }

  units <- sort(unit_rows$genotype_batch, method = "radix")
  n_units <- length(units)
  if (!n_units) {
    return(occupancy_harness_empty_aggregation("at least one unit is required"))
  }
  unit_meta <- unit_rows[match(units, unit_rows$genotype_batch), , drop = FALSE]
  unit_meta$genotype <- factor(unit_meta$genotype, levels = genotype_levels)
  batch_levels <- sort(unique(as.character(unit_meta$batch)), method = "radix")
  unit_meta$batch <- droplevels(factor(as.character(unit_meta$batch), levels = batch_levels))
  rownames(unit_meta) <- units

  unit_i <- match(unit, units)
  state_i <- match(state, states)
  primary <- !is.na(state_i)
  group_i <- unit_i[primary] + n_units * (state_i[primary] - 1L)
  state_counts <- matrix(
    tabulate(group_i, nbins = n_units * length(states)),
    nrow = n_units, ncol = length(states), dimnames = list(units, states)
  )
  unit_total <- tabulate(unit_i, nbins = n_units)
  unit_primary <- rowSums(state_counts)
  coverage <- unit_primary / unit_total
  unit_coverage <- data.frame(
    genotype_batch = units,
    genotype = unit_meta$genotype,
    batch = unit_meta$batch,
    n_retained = unit_total,
    n_primary = unit_primary,
    coverage = coverage,
    n_Homeostatic = state_counts[, "Homeostatic"],
    n_DAM = state_counts[, "DAM"],
    DAM_fraction = state_counts[, "DAM"] / unit_primary,
    row.names = NULL, stringsAsFactors = FALSE
  )

  reasons <- character()
  missing_state_units <- units[rowSums(state_counts > 0L) != length(states)]
  if (length(missing_state_units)) {
    reasons <- c(
      reasons,
      paste0("each retained unit must contain both primary states; failing unit(s): ",
             paste(missing_state_units, collapse = ", "))
    )
  }
  fd <- tryCatch(factorial_design(unit_meta), error = function(e) e)
  if (inherits(fd, "error")) {
    reasons <- c(reasons, paste0(
      "the present-unit 2x2 factorial design is not full-rank/estimable: ",
      conditionMessage(fd)
    ))
    design_rank <- design_columns <- residual_df <- NA_integer_
  } else {
    design_rank <- qr(fd$design)$rank
    design_columns <- ncol(fd$design)
    residual_df <- nrow(fd$design) - design_rank
  }

  list(
    unit_coverage = unit_coverage,
    unit_meta = unit_meta,
    diagnostics = list(
      status = if (length(reasons)) "non_estimable" else "estimable",
      reason = if (length(reasons)) paste(reasons, collapse = " | ") else NA_character_,
      n_units = n_units,
      design_rank = design_rank,
      design_columns = design_columns,
      residual_df = residual_df,
      missing_state_units = missing_state_units
    )
  )
}

occupancy_harness_capture <- function(expr) {
  tryCatch(
    list(ok = TRUE, value = force(expr), reason = NA_character_),
    error = function(e) list(ok = FALSE, value = NULL, reason = conditionMessage(e))
  )
}

occupancy_harness_status <- function(requested) {
  allowed <- c("beta_binomial", "empirical_logit", "simple_proportion")
  stats::setNames(lapply(allowed, function(estimator) list(
    status = if (estimator %in% requested) "pending" else "not_requested",
    reason = NA_character_, messages = character()
  )), allowed)
}

occupancy_family <- function(
    unit_coverage, unit_meta, margin = 0.10, n_perm = 9999L, seed = 614L,
    estimators = c("beta_binomial", "empirical_logit", "simple_proportion"),
    max_abs_coef = 20, max_vcov_condition = 1e10) {
  allowed_estimators <- c("beta_binomial", "empirical_logit", "simple_proportion")
  if (missing(estimators)) {
    estimators <- allowed_estimators
  } else {
    estimators <- unique(as.character(estimators))
  }
  stopifnot(
    length(estimators) >= 1L, all(estimators %in% allowed_estimators),
    length(margin) == 1L, is.finite(margin), margin > 0, margin < 1,
    length(n_perm) == 1L, is.finite(n_perm), n_perm >= 1L, n_perm == round(n_perm),
    length(seed) == 1L, is.finite(seed), seed == round(seed),
    length(max_abs_coef) == 1L, is.finite(max_abs_coef), max_abs_coef > 0,
    length(max_vcov_condition) == 1L, is.finite(max_vcov_condition),
    max_vcov_condition > 0
  )

  statuses <- occupancy_harness_status(estimators)
  n_units <- if (is.data.frame(unit_coverage)) nrow(unit_coverage) else 0L
  out <- list(
    unit = if (is.data.frame(unit_coverage)) unit_coverage else data.frame(),
    log_odds = NULL,
    probability_means = NULL,
    probability_vcov = NULL,
    probability_contrasts = NULL,
    empirical_logit = NULL,
    permutation = NULL,
    simple_proportion = NULL,
    diagnostics = list(
      convergence = NA_integer_, pdHess = NA, dispersion = NA_real_,
      max_abs_coef = NA_real_, max_abs_coef_gate = max_abs_coef,
      vcov_condition = NA_real_, vcov_condition_gate = max_vcov_condition,
      gradient_residual = NA_real_, n_units = n_units,
      design_rank = NA_integer_, design_columns = NA_integer_, residual_df = NA_integer_,
      messages = character(), estimators = statuses,
      status = "estimator_failed", reason = NA_character_, structural_reason = NA_character_
    )
  )

  prepared <- occupancy_harness_capture({
    if (!is.data.frame(unit_coverage) || !is.data.frame(unit_meta)) {
      stop("unit_coverage and unit_meta must be data frames", call. = FALSE)
    }
    required_schema <- c(
      "genotype_batch", "genotype", "batch", "n_retained", "n_Homeostatic", "n_DAM",
      "n_primary", "coverage", "DAM_fraction"
    )
    if (!all(required_schema %in% names(unit_coverage))) {
      stop("unit_coverage is missing required schema columns", call. = FALSE)
    }
    if (!all(c("genotype", "batch") %in% names(unit_meta))) {
      stop("unit_meta is missing genotype or batch", call. = FALSE)
    }
    if (nrow(unit_coverage) == 0L || nrow(unit_meta) != nrow(unit_coverage)) {
      stop("unit_coverage and unit_meta must contain the same positive number of units",
           call. = FALSE)
    }
    if (is.null(rownames(unit_meta)) ||
        !identical(as.character(unit_coverage$genotype_batch), rownames(unit_meta))) {
      stop("unit order must match rownames(unit_meta)", call. = FALSE)
    }
    if (anyDuplicated(unit_coverage$genotype_batch)) {
      stop("genotype_batch values must be unique", call. = FALSE)
    }
    count_cols <- c("n_retained", "n_Homeostatic", "n_DAM", "n_primary")
    if (!all(vapply(unit_coverage[count_cols], is.numeric, logical(1))) ||
        !all(is.finite(as.matrix(unit_coverage[count_cols]))) ||
        !all(as.matrix(unit_coverage[count_cols]) == round(as.matrix(unit_coverage[count_cols])))) {
      stop("unit counts must be finite integer-valued numerics", call. = FALSE)
    }
    if (any(unit_coverage$n_retained <= 0L) ||
        any(unit_coverage$n_Homeostatic <= 0L) || any(unit_coverage$n_DAM <= 0L)) {
      stop("each retained unit must contain positive counts for both primary states",
           call. = FALSE)
    }
    if (!all(unit_coverage$n_primary ==
             unit_coverage$n_Homeostatic + unit_coverage$n_DAM)) {
      stop("n_primary must equal n_Homeostatic + n_DAM", call. = FALSE)
    }
    if (any(unit_coverage$n_primary > unit_coverage$n_retained)) {
      stop("n_primary cannot exceed n_retained", call. = FALSE)
    }
    if (!isTRUE(all.equal(
      unit_coverage$coverage, unit_coverage$n_primary / unit_coverage$n_retained,
      tolerance = 0
    ))) {
      stop("coverage must equal n_primary / n_retained", call. = FALSE)
    }
    if (!isTRUE(all.equal(
      unit_coverage$DAM_fraction, unit_coverage$n_DAM / unit_coverage$n_primary,
      tolerance = 0
    ))) {
      stop("DAM_fraction must equal n_DAM / n_primary", call. = FALSE)
    }
    if (!identical(as.character(unit_coverage$genotype),
                   as.character(unit_meta$genotype)) ||
        !identical(as.character(unit_coverage$batch), as.character(unit_meta$batch))) {
      stop("unit_coverage genotype/batch metadata must match unit_meta", call. = FALSE)
    }

    d <- unit_coverage
    d$genotype <- factor(as.character(d$genotype), levels = genotype_levels)
    if (anyNA(d$genotype)) {
      stop("all units must use canonical genotype levels", call. = FALSE)
    }
    batch_levels <- sort(unique(as.character(d$batch)), method = "radix")
    d$batch <- droplevels(factor(as.character(d$batch), levels = batch_levels))
    if (anyNA(d$batch)) {
      stop("all units must have a present batch level", call. = FALSE)
    }
    rownames(d) <- d$genotype_batch
    fd <- tryCatch(
      factorial_design(d),
      error = function(e) stop(
        "present-unit factorial design is not full-rank/estimable: ",
        conditionMessage(e), call. = FALSE
      )
    )
    list(d = d, fd = fd)
  })

  if (!prepared$ok) {
    for (estimator in estimators) {
      statuses[[estimator]]$status <- "estimator_failed"
      statuses[[estimator]]$reason <- prepared$reason
    }
    out$diagnostics$estimators <- statuses
    out$diagnostics$reason <- prepared$reason
    out$diagnostics$structural_reason <- prepared$reason
    return(out)
  }

  d <- prepared$value$d
  fd <- prepared$value$fd
  unit_required <- c(
    "genotype_batch", "genotype", "batch", "n_Homeostatic", "n_DAM",
    "n_primary", "coverage", "DAM_fraction"
  )
  out$unit <- d[unit_required]
  out$diagnostics$n_units <- nrow(d)
  out$diagnostics$design_rank <- qr(fd$design)$rank
  out$diagnostics$design_columns <- ncol(fd$design)
  out$diagnostics$residual_df <- nrow(fd$design) - ncol(fd$design)

  if ("beta_binomial" %in% estimators) {
    beta_result <- occupancy_harness_capture({
      cap <- state_capture_clean(
        glmmTMB::glmmTMB(
          cbind(n_DAM, n_Homeostatic) ~ 0 + genotype + batch,
          data = d, family = glmmTMB::betabinomial(link = "logit")
        ),
        "beta-binomial DAM occupancy"
      )
      fit <- cap$value
      beta <- glmmTMB::fixef(fit)$cond
      vc <- stats::vcov(fit)$cond
      se_beta <- sqrt(diag(vc))
      dispersion <- stats::sigma(fit)
      vcov_condition <- kappa(vc)
      if (is.null(fit$fit$convergence) || fit$fit$convergence != 0L) {
        stop("beta-binomial did not converge", call. = FALSE)
      }
      if (!isTRUE(fit$sdr$pdHess)) {
        stop("beta-binomial Hessian is not positive definite", call. = FALSE)
      }
      if (!all(is.finite(beta)) || !all(is.finite(se_beta)) || any(se_beta <= 0)) {
        stop("beta-binomial fixed effects or standard errors are invalid", call. = FALSE)
      }
      if (max(abs(beta)) > max_abs_coef) {
        stop("beta-binomial fixed-effect magnitude exceeds its gate", call. = FALSE)
      }
      if (!all(is.finite(vc)) ||
          min(eigen(vc, symmetric = TRUE, only.values = TRUE)$values) <= 0) {
        stop("beta-binomial fixed-effect covariance is not positive definite",
             call. = FALSE)
      }
      if (!is.finite(vcov_condition) || vcov_condition > max_vcov_condition) {
        stop("beta-binomial covariance condition number exceeds its gate", call. = FALSE)
      }
      if (!is.finite(dispersion) || dispersion <= 0) {
        stop("beta-binomial dispersion is invalid", call. = FALSE)
      }

      model_x <- stats::model.matrix(~ 0 + genotype + batch, d)
      if (!identical(colnames(model_x), names(beta)) ||
          qr(model_x)$rank != ncol(model_x)) {
        stop("beta-binomial model matrix is rank-deficient or misaligned", call. = FALSE)
      }
      if (max(abs(drop(model_x %*% beta) - stats::predict(fit, type = "link"))) >= 1e-8) {
        stop("beta-binomial fixed-effect reconstruction failed", call. = FALSE)
      }

      cm <- matrix(
        0, nrow = length(beta), ncol = length(state_contrast_names()),
        dimnames = list(names(beta), state_contrast_names())
      )
      cm[paste0("genotype", genotype_levels), ] <- state_genotype_contrasts()
      log_odds <- state_wald_contrasts(
        beta, vc, cm, family = "occupancy_logodds_zero_all_contrasts"
      )

      grid <- expand.grid(
        genotype = genotype_levels, batch = levels(d$batch),
        KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
      )
      grid$genotype <- factor(grid$genotype, levels = genotype_levels)
      grid$batch <- factor(grid$batch, levels = levels(d$batch))
      x_new <- stats::model.matrix(~ 0 + genotype + batch, grid)
      if (!identical(colnames(x_new), names(beta))) {
        stop("standardization grid is misaligned with beta-binomial coefficients",
             call. = FALSE)
      }
      probability <- state_probability_standardization(
        beta, vc, x_new, as.character(grid$genotype), margin
      )
      list(
        log_odds = log_odds,
        probability = probability,
        convergence = fit$fit$convergence,
        pdHess = fit$sdr$pdHess,
        dispersion = unname(dispersion),
        max_abs_coef = max(abs(beta)),
        vcov_condition = unname(vcov_condition),
        messages = cap$messages
      )
    })
    if (beta_result$ok) {
      b <- beta_result$value
      out$log_odds <- b$log_odds
      out$probability_means <- b$probability$means
      out$probability_vcov <- b$probability$mean_vcov
      out$probability_contrasts <- b$probability$contrasts
      out$diagnostics$convergence <- b$convergence
      out$diagnostics$pdHess <- b$pdHess
      out$diagnostics$dispersion <- b$dispersion
      out$diagnostics$max_abs_coef <- b$max_abs_coef
      out$diagnostics$vcov_condition <- b$vcov_condition
      out$diagnostics$gradient_residual <- b$probability$gradient_residual
      out$diagnostics$messages <- b$messages
      statuses$beta_binomial <- list(
        status = "ok", reason = NA_character_, messages = b$messages
      )
    } else {
      statuses$beta_binomial$status <- "estimator_failed"
      statuses$beta_binomial$reason <- beta_result$reason
    }
  }

  if ("empirical_logit" %in% estimators) {
    empirical_result <- occupancy_harness_capture({
      empirical_logit <- log((d$n_DAM + 0.5) / (d$n_Homeostatic + 0.5))
      emat <- matrix(
        empirical_logit, nrow = 1L,
        dimnames = list("empirical_logit", rownames(fd$design))
      )
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
      if (abs(fl$t_obs - ols$t[ols$contrast == "interaction"]) >= 1e-8) {
        stop("empirical-logit OLS and permutation t statistics disagree", call. = FALSE)
      }
      list(
        table = ols,
        permutation = data.frame(
          contrast = "interaction", t_obs = fl$t_obs, n_perm = fl$n_perm,
          perm_p = fl$perm_p, seed = fl$seed,
          family = "occupancy_empirical_logit_interaction_permutation",
          row.names = NULL, stringsAsFactors = FALSE
        ),
        messages = ols_cap$messages
      )
    })
    if (empirical_result$ok) {
      out$empirical_logit <- empirical_result$value$table
      out$permutation <- empirical_result$value$permutation
      statuses$empirical_logit <- list(
        status = "ok", reason = NA_character_,
        messages = empirical_result$value$messages
      )
    } else {
      statuses$empirical_logit$status <- "estimator_failed"
      statuses$empirical_logit$reason <- empirical_result$reason
    }
  }

  if ("simple_proportion" %in% estimators) {
    simple_result <- occupancy_harness_capture({
      smat <- matrix(
        d$DAM_fraction, nrow = 1L,
        dimnames = list("simple_proportion", rownames(fd$design))
      )
      ols_cap <- state_capture_clean(
        fit_trajectory_contrasts(smat, fd$design, fd$contrasts),
        "simple-proportion occupancy OLS"
      )
      ols <- do.call(rbind, ols_cap$value$top)
      rownames(ols) <- NULL
      ols$fdr <- stats::p.adjust(ols$p_value, "BH")
      ols$family <- "occupancy_simple_proportion_zero_all_contrasts"
      list(table = ols, messages = ols_cap$messages)
    })
    if (simple_result$ok) {
      out$simple_proportion <- simple_result$value$table
      statuses$simple_proportion <- list(
        status = "ok", reason = NA_character_, messages = simple_result$value$messages
      )
    } else {
      statuses$simple_proportion$status <- "estimator_failed"
      statuses$simple_proportion$reason <- simple_result$reason
    }
  }

  failed <- names(statuses)[vapply(
    statuses, function(x) identical(x$status, "estimator_failed"), logical(1)
  )]
  out$diagnostics$estimators <- statuses
  out$diagnostics$status <- if (length(failed)) "estimator_failed" else "ok"
  out$diagnostics$reason <- if (length(failed)) {
    paste(vapply(failed, function(estimator) paste0(
      estimator, ": ", statuses[[estimator]]$reason
    ), character(1)), collapse = " | ")
  } else {
    NA_character_
  }
  out
}

occupancy_from_membership <- function(
    membership, meta, states = c("Homeostatic", "DAM"),
    unit_col = "genotype_batch", genotype_col = "genotype", batch_col = "batch",
    margin = 0.10, n_perm = 9999L, seed = 614L,
    estimators = c("beta_binomial", "empirical_logit", "simple_proportion"),
    max_abs_coef = 20, max_vcov_condition = 1e10) {
  aggregated <- membership_to_unit_coverage(
    membership = membership, meta = meta, states = states,
    unit_col = unit_col, genotype_col = genotype_col, batch_col = batch_col
  )
  out <- occupancy_family(
    unit_coverage = aggregated$unit_coverage,
    unit_meta = aggregated$unit_meta,
    margin = margin, n_perm = n_perm, seed = seed, estimators = estimators,
    max_abs_coef = max_abs_coef, max_vcov_condition = max_vcov_condition
  )
  out$membership_diagnostics <- aggregated$diagnostics
  out
}

# Deterministic arbitrary counts used only to exercise reduced-design plumbing. None of these
# counts, omissions, or unit rows is derived from the project data.
occupancy_harness_fabricated_unit_data <- function(batches, drop_unit = NULL) {
  grid <- expand.grid(
    genotype = genotype_levels, batch = as.character(batches),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  grid$genotype_batch <- paste(grid$genotype, grid$batch, sep = "_")
  if (!is.null(drop_unit)) {
    grid <- grid[grid$genotype_batch != drop_unit, , drop = FALSE]
  }
  grid <- grid[order(grid$genotype_batch, method = "radix"), , drop = FALSE]
  gi <- match(grid$genotype, genotype_levels)
  bi <- match(grid$batch, as.character(batches))
  residual <- ((gi * 11L + bi * 7L) %% 13L) - 6L
  n_dam <- as.integer(28L + 7L * gi + 3L * bi + 2L * (gi == 4L) * bi + residual)
  n_homeostatic <- as.integer(
    128L - 5L * gi + 4L * bi + ((gi * 3L + bi * 5L) %% 11L)
  )
  n_primary <- n_homeostatic + n_dam
  n_retained <- n_primary + as.integer((gi + bi) %% 3L)
  batch_levels <- sort(unique(grid$batch), method = "radix")
  unit_meta <- data.frame(
    genotype_batch = grid$genotype_batch,
    genotype = factor(grid$genotype, levels = genotype_levels),
    batch = droplevels(factor(grid$batch, levels = batch_levels)),
    row.names = grid$genotype_batch, stringsAsFactors = FALSE
  )
  unit_coverage <- data.frame(
    genotype_batch = grid$genotype_batch,
    genotype = unit_meta$genotype,
    batch = unit_meta$batch,
    n_retained = n_retained,
    n_primary = n_primary,
    coverage = n_primary / n_retained,
    n_Homeostatic = n_homeostatic,
    n_DAM = n_dam,
    DAM_fraction = n_dam / n_primary,
    row.names = NULL, stringsAsFactors = FALSE
  )
  list(unit_coverage = unit_coverage, unit_meta = unit_meta)
}

occupancy_harness_finite_table <- function(x) {
  if (!is.data.frame(x) || nrow(x) == 0L) return(FALSE)
  numeric_columns <- vapply(x, is.numeric, logical(1))
  any(numeric_columns) && all(is.finite(as.matrix(x[numeric_columns])))
}

occupancy_harness_assert_smoke <- function(x, n_units, n_batches) {
  stopifnot(
    is.list(x), is.data.frame(x$unit), nrow(x$unit) == n_units,
    x$diagnostics$n_units == n_units,
    x$diagnostics$design_rank == x$diagnostics$design_columns,
    x$diagnostics$residual_df == n_units - x$diagnostics$design_columns,
    length(unique(as.character(x$unit$batch))) == n_batches
  )
  check_failed <- function(estimator) {
    status <- x$diagnostics$estimators[[estimator]]
    stopifnot(
      identical(status$status, "estimator_failed"),
      length(status$reason) == 1L, !is.na(status$reason), nzchar(status$reason)
    )
  }

  if (identical(x$diagnostics$estimators$beta_binomial$status, "ok")) {
    stopifnot(
      nrow(x$probability_contrasts) == 5L,
      nrow(x$probability_means) == length(genotype_levels),
      identical(dim(x$probability_vcov), c(length(genotype_levels), length(genotype_levels))),
      nrow(x$log_odds) == 5L,
      occupancy_harness_finite_table(x$probability_contrasts),
      occupancy_harness_finite_table(x$probability_means),
      all(is.finite(x$probability_vcov)),
      occupancy_harness_finite_table(x$log_odds)
    )
  } else {
    check_failed("beta_binomial")
  }
  if (identical(x$diagnostics$estimators$empirical_logit$status, "ok")) {
    stopifnot(
      nrow(x$empirical_logit) == 5L, nrow(x$permutation) == 1L,
      occupancy_harness_finite_table(x$empirical_logit),
      occupancy_harness_finite_table(x$permutation)
    )
  } else {
    check_failed("empirical_logit")
  }
  if (identical(x$diagnostics$estimators$simple_proportion$status, "ok")) {
    stopifnot(
      nrow(x$simple_proportion) == 5L,
      occupancy_harness_finite_table(x$simple_proportion)
    )
  } else {
    check_failed("simple_proportion")
  }
  invisible(TRUE)
}

occupancy_harness_reduced_design_smoke <- function() {
  # LOU-shaped: one arbitrary fabricated unit omitted from a fabricated 4x4 layout.
  lou <- occupancy_harness_fabricated_unit_data(
    batches = as.character(seq_len(4L)), drop_unit = "MAPTKI_1"
  )
  lou_fit <- occupancy_family(
    lou$unit_coverage, lou$unit_meta, n_perm = 99L, seed = 614L
  )
  occupancy_harness_assert_smoke(lou_fit, n_units = 15L, n_batches = 4L)

  # LOBO-shaped: a separate fabricated complete 4x3 layout; no project batch is removed.
  lobo <- occupancy_harness_fabricated_unit_data(batches = as.character(seq_len(3L)))
  lobo_fit <- occupancy_family(
    lobo$unit_coverage, lobo$unit_meta, n_perm = 99L, seed = 614L
  )
  occupancy_harness_assert_smoke(lobo_fit, n_units = 12L, n_batches = 3L)
  "ok"
}

occupancy_harness_max_abs_diff <- function(x, y) {
  if (is.data.frame(x) && is.data.frame(y)) {
    numeric_columns <- names(x)[vapply(x, is.numeric, logical(1))]
    values <- unlist(lapply(numeric_columns, function(column) {
      as.numeric(x[[column]]) - as.numeric(y[[column]])
    }), use.names = FALSE)
    return(max(c(0, abs(values))))
  }
  if ((is.matrix(x) || is.numeric(x)) && (is.matrix(y) || is.numeric(y))) {
    return(max(c(0, abs(as.numeric(x) - as.numeric(y)))))
  }
  0
}

check_occupancy_harness <- function(substrate, response) {
  stopifnot(
    is.list(substrate), is.list(response), is.list(response$occupancy),
    is.data.frame(substrate$unit_coverage), is.data.frame(substrate$unit_meta)
  )

  # Layer B: current-label estimator core only. This is the sole real occupancy evaluation in P7.4.
  hc <- occupancy_family(substrate$unit_coverage, substrate$unit_meta)
  reference <- response$occupancy
  stopifnot(
    isTRUE(all.equal(
      hc$probability_contrasts, reference$probability_contrasts, tolerance = 1e-8
    )),
    isTRUE(all.equal(hc$probability_means, reference$probability_means, tolerance = 1e-8)),
    isTRUE(all.equal(hc$probability_vcov, reference$probability_vcov, tolerance = 1e-8)),
    isTRUE(all.equal(hc$log_odds, reference$log_odds, tolerance = 1e-8)),
    isTRUE(all.equal(hc$empirical_logit, reference$empirical_logit, tolerance = 1e-8)),
    isTRUE(all.equal(hc$permutation, reference$permutation, tolerance = 1e-8)),
    identical(hc$diagnostics$estimators$beta_binomial$status, "ok"),
    identical(hc$diagnostics$estimators$empirical_logit$status, "ok"),
    identical(hc$diagnostics$estimators$simple_proportion$status, "ok")
  )

  # Layer A: expand only the established aggregate counts into arithmetic membership rows. This
  # reconstructs no cell identity and creates no alternate labeling.
  d <- substrate$unit_coverage
  expanded <- do.call(rbind, lapply(seq_len(nrow(d)), function(i) {
    state <- c(
      rep("Homeostatic", as.integer(d$n_Homeostatic[i])),
      rep("DAM", as.integer(d$n_DAM[i]))
    )
    data.frame(
      membership = state,
      genotype_batch = rep(as.character(d$genotype_batch[i]), length(state)),
      genotype = rep(as.character(d$genotype[i]), length(state)),
      batch = rep(as.character(d$batch[i]), length(state)),
      stringsAsFactors = FALSE
    )
  }))
  aggregation <- membership_to_unit_coverage(
    expanded$membership,
    expanded[c("genotype_batch", "genotype", "batch")]
  )
  ai <- match(d$genotype_batch, aggregation$unit_coverage$genotype_batch)
  stopifnot(!anyNA(ai), identical(aggregation$diagnostics$status, "estimable"))
  actual <- aggregation$unit_coverage[ai, , drop = FALSE]
  stopifnot(
    identical(as.integer(actual$n_Homeostatic), as.integer(d$n_Homeostatic)),
    identical(as.integer(actual$n_DAM), as.integer(d$n_DAM)),
    isTRUE(all.equal(actual$n_primary, d$n_primary, tolerance = 0)),
    isTRUE(all.equal(actual$DAM_fraction, d$DAM_fraction, tolerance = 0))
  )

  reduced_design_smoke <- occupancy_harness_reduced_design_smoke()
  e1 <- hc$probability_contrasts[hc$probability_contrasts$contrast == "interaction", ]
  e2 <- hc$empirical_logit[hc$empirical_logit$contrast == "interaction", ]
  e3 <- hc$simple_proportion[hc$simple_proportion$contrast == "interaction", ]
  stopifnot(nrow(e1) == 1L, nrow(e2) == 1L, nrow(e3) == 1L)

  list(
    status = "reproduced",
    e1_interaction = data.frame(
      estimate = unname(e1$estimate), ci_l = unname(e1$ci_l), ci_r = unname(e1$ci_r),
      fdr_zero = unname(e1$fdr_zero), fdr_minimum = unname(e1$fdr_minimum),
      row.names = NULL
    ),
    e2_empirical_logit_interaction = data.frame(
      coef = unname(e2$coef), perm_p = unname(hc$permutation$perm_p),
      row.names = NULL
    ),
    e3_simple_proportion_interaction = data.frame(
      estimate = unname(e3$coef), fdr = unname(e3$fdr), row.names = NULL
    ),
    max_abs_diff = c(
      probability_contrasts = occupancy_harness_max_abs_diff(
        hc$probability_contrasts, reference$probability_contrasts
      ),
      probability_means = occupancy_harness_max_abs_diff(
        hc$probability_means, reference$probability_means
      ),
      probability_vcov = occupancy_harness_max_abs_diff(
        hc$probability_vcov, reference$probability_vcov
      ),
      log_odds = occupancy_harness_max_abs_diff(hc$log_odds, reference$log_odds),
      empirical_logit = occupancy_harness_max_abs_diff(
        hc$empirical_logit, reference$empirical_logit
      ),
      permutation = occupancy_harness_max_abs_diff(hc$permutation, reference$permutation)
    ),
    reduced_design_smoke = reduced_design_smoke
  )
}
