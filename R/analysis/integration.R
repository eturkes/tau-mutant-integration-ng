# P8 integration is restricted to the P7.3-compliant five-contrast x symbol
# effect-size space: no animal-level pairing is attempted. The protein-group
# view is the single bulk modality; phosphosite-to-parent collapse is retained
# only as a within-assay alternate and is never counted as a fourth modality.

integration_contrast_names <- function() {
  c("tau_alone", "nlgf_in_maptki", "nlgf_in_p301s", "tau_in_nlgf", "interaction")
}

integration_radix_sort <- function(x) {
  sort(as.character(x), method = "radix")
}

integration_validate_top_tables <- function(top, contrasts, feature_col, extra_cols = character()) {
  required <- unique(c(feature_col, "AveExpr", "logFC", "t", extra_cols))
  stopifnot(is.list(top), identical(names(top), contrasts))

  first <- top[[contrasts[[1L]]]]
  stopifnot(is.data.frame(first), all(required %in% names(first)), nrow(first) > 0L)
  features <- as.character(first[[feature_col]])
  ave_expr <- first$AveExpr
  stopifnot(
    !anyNA(features), all(nzchar(features)), anyDuplicated(features) == 0L,
    is.numeric(ave_expr), all(is.finite(ave_expr))
  )

  for (contrast in contrasts) {
    tab <- top[[contrast]]
    stopifnot(is.data.frame(tab), all(required %in% names(tab)), nrow(tab) == nrow(first))
    tab_features <- as.character(tab[[feature_col]])
    stopifnot(!anyNA(tab_features), all(nzchar(tab_features)), anyDuplicated(tab_features) == 0L)
    idx <- match(features, tab_features)
    stopifnot(
      !anyNA(idx), identical(tab_features[idx], features),
      is.numeric(tab$AveExpr),
      identical(as.vector(tab$AveExpr[idx]), as.vector(ave_expr)),
      is.numeric(tab$logFC), is.numeric(tab$t),
      all(is.finite(tab$logFC)), all(is.finite(tab$t))
    )
  }

  list(first = first, features = features, ave_expr = ave_expr)
}

integration_select_representatives <- function(features, symbols, ave_expr, valid) {
  stopifnot(
    length(features) == length(symbols), length(features) == length(ave_expr),
    length(features) == length(valid), is.logical(valid), !anyNA(valid), any(valid)
  )
  symbols <- as.character(symbols)
  candidates <- data.frame(
    symbol = symbols[valid],
    feature = as.character(features[valid]),
    AveExpr = ave_expr[valid],
    stringsAsFactors = FALSE
  )
  stopifnot(
    !anyNA(candidates$symbol), all(nzchar(candidates$symbol)),
    identical(candidates$symbol, trimws(candidates$symbol)),
    !anyNA(candidates$feature), all(nzchar(candidates$feature)),
    all(is.finite(candidates$AveExpr))
  )

  ord <- order(
    candidates$symbol, -candidates$AveExpr, candidates$feature,
    method = "radix"
  )
  candidates <- candidates[ord, , drop = FALSE]
  selected <- candidates[!duplicated(candidates$symbol), , drop = FALSE]
  selected <- selected[order(selected$symbol, method = "radix"), , drop = FALSE]
  rownames(selected) <- NULL

  # Oracle the representative rule itself: maximum AveExpr, then radix-min feature id.
  max_ave <- tapply(candidates$AveExpr, candidates$symbol, max)
  at_max <- candidates$AveExpr == unname(max_ave[candidates$symbol])
  min_feature <- tapply(
    candidates$feature[at_max], candidates$symbol[at_max],
    function(x) integration_radix_sort(x)[[1L]]
  )
  stopifnot(
    identical(as.vector(selected$AveExpr), as.vector(unname(max_ave[selected$symbol]))),
    identical(as.vector(selected$feature), as.vector(unname(min_feature[selected$symbol]))),
    identical(selected$symbol, integration_radix_sort(selected$symbol)),
    anyDuplicated(selected$symbol) == 0L
  )
  selected
}

integration_effect_matrix <- function(top, contrasts, feature_col, representatives, statistic) {
  stopifnot(statistic %in% c("logFC", "t"), is.data.frame(representatives))
  rows <- lapply(contrasts, function(contrast) {
    tab <- top[[contrast]]
    idx <- match(representatives$feature, as.character(tab[[feature_col]]))
    stopifnot(!anyNA(idx))
    values <- tab[[statistic]][idx]
    stopifnot(is.numeric(values), all(is.finite(values)))
    values
  })
  out <- do.call(rbind, rows)
  dimnames(out) <- list(contrasts, representatives$symbol)

  # Exact, tolerance-zero reconstruction from every representative source row.
  for (j in seq_along(contrasts)) {
    tab <- top[[contrasts[[j]]]]
    idx <- match(representatives$feature, as.character(tab[[feature_col]]))
    stopifnot(identical(as.vector(out[j, ]), as.vector(tab[[statistic]][idx])))
  }
  out
}

robust_z <- function(x) {
  stopifnot(
    is.matrix(x), is.numeric(x), nrow(x) > 0L, ncol(x) > 0L,
    !is.null(rownames(x)), !is.null(colnames(x)), all(is.finite(x))
  )
  center <- apply(x, 1L, stats::median)
  scale <- vapply(
    seq_len(nrow(x)),
    function(i) stats::mad(x[i, ], center = center[[i]]),
    numeric(1)
  )
  names(center) <- rownames(x)
  names(scale) <- rownames(x)
  stopifnot(all(is.finite(center)), all(is.finite(scale)), all(scale > 0))

  z <- sweep(sweep(x, 1L, center, "-"), 1L, scale, "/")
  dimnames(z) <- dimnames(x)
  reconstructed <- sweep(sweep(z, 1L, scale, "*"), 1L, center, "+")
  stopifnot(all(is.finite(z)), max(abs(reconstructed - x)) < 1e-9)

  list(z = z, center = center, scale = scale)
}

integration_contains_parent <- function(x) {
  forbidden <- c(
    "Seurat", "Assay", "Assay5", "DimReduc", "Neighbor",
    "DGEList", "EList", "MArrayLM", "lm", "glm", "glmmTMB"
  )
  if (inherits(x, forbidden) || isS4(x) || is.function(x) || is.environment(x) ||
      typeof(x) == "externalptr") return(TRUE)
  if (!is.list(x)) return(FALSE)
  any(vapply(unclass(x), integration_contains_parent, logical(1)))
}

build_integration_substrate <- function(
    pb_de_microglia, symbol_map, geomx_de, proteome_de_24m, phospho_de_24m) {
  modalities <- c("snRNAseq", "GeoMx", "bulk")
  contrasts <- integration_contrast_names()
  stopifnot(
    is.list(pb_de_microglia), is.list(geomx_de),
    is.list(proteome_de_24m), is.list(phospho_de_24m),
    is.data.frame(symbol_map), all(c("ensembl", "symbol") %in% names(symbol_map)),
    !anyNA(symbol_map$ensembl), !anyNA(symbol_map$symbol),
    all(nzchar(as.character(symbol_map$ensembl))),
    all(nzchar(as.character(symbol_map$symbol))),
    anyDuplicated(symbol_map$ensembl) == 0L,
    anyDuplicated(symbol_map$symbol) == 0L
  )

  sn_top <- pb_de_microglia$top
  geomx_top <- geomx_de$primary$top
  bulk_top <- proteome_de_24m$top
  phospho_top <- phospho_de_24m$top

  sn_info <- integration_validate_top_tables(sn_top, contrasts, "gene")
  geomx_info <- integration_validate_top_tables(geomx_top, contrasts, "symbol")
  bulk_info <- integration_validate_top_tables(bulk_top, contrasts, "feature", "gene_first")
  phospho_info <- integration_validate_top_tables(phospho_top, contrasts, "feature", "gene")

  sn_map_i <- match(sn_info$features, as.character(symbol_map$ensembl))
  stopifnot(!anyNA(sn_map_i))
  sn_symbol_all <- as.character(symbol_map$symbol[sn_map_i])
  sn_valid <- !is.na(sn_symbol_all) & nzchar(trimws(sn_symbol_all))

  geomx_symbol_all <- as.character(geomx_info$first$symbol)
  geomx_valid <- !is.na(geomx_symbol_all) & nzchar(trimws(geomx_symbol_all))

  bulk_symbol_all <- as.character(bulk_info$first$gene_first)
  bulk_valid <- !is.na(bulk_symbol_all) & nzchar(trimws(bulk_symbol_all))

  phospho_symbol_all <- as.character(phospho_info$first$gene)
  phospho_valid <- !is.na(phospho_symbol_all) &
    nzchar(trimws(phospho_symbol_all)) &
    !grepl("[;,]", phospho_symbol_all)

  sn_representatives <- integration_select_representatives(
    sn_info$features, sn_symbol_all, sn_info$ave_expr, sn_valid
  )
  geomx_representatives <- integration_select_representatives(
    geomx_info$features, geomx_symbol_all, geomx_info$ave_expr, geomx_valid
  )
  bulk_representatives <- integration_select_representatives(
    bulk_info$features, bulk_symbol_all, bulk_info$ave_expr, bulk_valid
  )
  phospho_representatives <- integration_select_representatives(
    phospho_info$features, phospho_symbol_all, phospho_info$ave_expr, phospho_valid
  )

  # The two transcriptomic modalities are one-to-one and must retain every source feature.
  stopifnot(
    identical(
      integration_radix_sort(sn_representatives$feature),
      integration_radix_sort(sn_info$features)
    ),
    identical(
      integration_radix_sort(geomx_representatives$feature),
      integration_radix_sort(geomx_info$features)
    )
  )

  symbols <- list(
    snRNAseq = sn_representatives$symbol,
    GeoMx = geomx_representatives$symbol,
    bulk = bulk_representatives$symbol
  )
  representatives <- list(
    snRNAseq = list(rows = sn_representatives, top = sn_top, feature_col = "gene"),
    GeoMx = list(rows = geomx_representatives, top = geomx_top, feature_col = "symbol"),
    bulk = list(rows = bulk_representatives, top = bulk_top, feature_col = "feature")
  )

  raw <- setNames(vector("list", length(modalities)), modalities)
  standardized <- setNames(vector("list", length(modalities)), modalities)
  standardization <- setNames(vector("list", length(modalities)), modalities)
  for (modality in modalities) {
    spec <- representatives[[modality]]
    raw[[modality]] <- list(
      logFC = integration_effect_matrix(
        spec$top, contrasts, spec$feature_col, spec$rows, "logFC"
      ),
      t = integration_effect_matrix(
        spec$top, contrasts, spec$feature_col, spec$rows, "t"
      )
    )
    standardized[[modality]] <- list()
    standardization[[modality]] <- list()
    for (statistic in c("logFC", "t")) {
      rz <- robust_z(raw[[modality]][[statistic]])
      standardized[[modality]][[statistic]] <- rz$z
      standardization[[modality]][[statistic]] <- rz[c("center", "scale")]
      stopifnot(
        max(abs(
          sweep(
            sweep(rz$z, 1L, rz$scale, "*"),
            1L, rz$center, "+"
          ) - raw[[modality]][[statistic]]
        )) < 1e-9
      )
    }
  }

  union_symbols <- integration_radix_sort(unique(unlist(symbols, use.names = FALSE)))
  membership <- as.integer(
    (union_symbols %in% symbols$snRNAseq) +
      (union_symbols %in% symbols$GeoMx) +
      (union_symbols %in% symbols$bulk)
  )
  names(membership) <- union_symbols
  pairwise <- list(
    snRNAseq_GeoMx = integration_radix_sort(intersect(symbols$snRNAseq, symbols$GeoMx)),
    snRNAseq_bulk = integration_radix_sort(intersect(symbols$snRNAseq, symbols$bulk)),
    GeoMx_bulk = integration_radix_sort(intersect(symbols$GeoMx, symbols$bulk))
  )
  index <- list(
    complete_case = union_symbols[membership == 3L],
    at_least_two = union_symbols[membership >= 2L],
    union = union_symbols,
    pairwise = pairwise,
    membership = membership
  )

  phospho_parent_alt <- list(
    symbols = phospho_representatives$symbol,
    logFC = integration_effect_matrix(
      phospho_top, contrasts, "feature", phospho_representatives, "logFC"
    ),
    t = integration_effect_matrix(
      phospho_top, contrasts, "feature", phospho_representatives, "t"
    ),
    note = paste0(
      "phosphosite -> parent-gene collapse of the SAME TiO2 assay as bulk; ",
      "not a 4th modality"
    )
  )

  counts <- c(
    snRNAseq_features = as.integer(nrow(sn_info$first)),
    GeoMx_features = as.integer(nrow(geomx_info$first)),
    bulk_features = as.integer(nrow(bulk_info$first)),
    phospho_sites = as.integer(nrow(phospho_info$first)),
    snRNAseq_symbol_na = as.integer(sum(is.na(sn_symbol_all))),
    snRNAseq_symbol_blank = as.integer(sum(!is.na(sn_symbol_all) & !nzchar(trimws(sn_symbol_all)))),
    snRNAseq_symbol_duplicates = as.integer(sum(duplicated(sn_symbol_all[sn_valid]))),
    GeoMx_symbol_na = as.integer(sum(is.na(geomx_symbol_all))),
    GeoMx_symbol_blank = as.integer(sum(!is.na(geomx_symbol_all) & !nzchar(trimws(geomx_symbol_all)))),
    GeoMx_symbol_duplicates = as.integer(sum(duplicated(geomx_symbol_all[geomx_valid]))),
    bulk_gene_first_na = as.integer(sum(is.na(bulk_symbol_all))),
    bulk_gene_first_blank = as.integer(sum(!is.na(bulk_symbol_all) & !nzchar(trimws(bulk_symbol_all)))),
    bulk_gene_first_duplicates = as.integer(sum(duplicated(bulk_symbol_all[bulk_valid]))),
    snRNAseq_symbols = as.integer(length(symbols$snRNAseq)),
    GeoMx_symbols = as.integer(length(symbols$GeoMx)),
    bulk_symbols = as.integer(length(symbols$bulk)),
    pairwise_snRNAseq_GeoMx = as.integer(length(pairwise$snRNAseq_GeoMx)),
    pairwise_snRNAseq_bulk = as.integer(length(pairwise$snRNAseq_bulk)),
    pairwise_GeoMx_bulk = as.integer(length(pairwise$GeoMx_bulk)),
    complete_case = as.integer(length(index$complete_case)),
    at_least_two = as.integer(length(index$at_least_two)),
    union = as.integer(length(index$union)),
    phospho_parent_symbols = as.integer(length(phospho_parent_alt$symbols))
  )
  expected_counts <- c(
    snRNAseq_features = 14512L,
    GeoMx_features = 19959L,
    bulk_features = 3379L,
    phospho_sites = 17707L,
    snRNAseq_symbol_na = 0L,
    snRNAseq_symbol_blank = 0L,
    snRNAseq_symbol_duplicates = 0L,
    GeoMx_symbol_na = 0L,
    GeoMx_symbol_blank = 0L,
    GeoMx_symbol_duplicates = 0L,
    bulk_gene_first_na = 15L,
    bulk_gene_first_blank = 0L,
    bulk_gene_first_duplicates = 58L,
    snRNAseq_symbols = 14512L,
    GeoMx_symbols = 19959L,
    bulk_symbols = 3306L,
    pairwise_snRNAseq_GeoMx = 12324L,
    pairwise_snRNAseq_bulk = 3132L,
    pairwise_GeoMx_bulk = 3189L,
    complete_case = 3109L,
    at_least_two = 12427L,
    union = 22241L,
    phospho_parent_symbols = 3019L
  )
  stopifnot(identical(counts, expected_counts))

  result <- list(
    modalities = modalities,
    contrasts = contrasts,
    symbols = symbols,
    raw = raw,
    standardized = standardized,
    standardization = standardization,
    index = index,
    phospho_parent_alt = phospho_parent_alt,
    provenance = list(
      counts = counts,
      representative_rule = paste0(
        "select one feature per symbol once by maximum AveExpr; ties use the ",
        "radix-min source feature id"
      ),
      standardization = paste0(
        "within each modality/statistic/contrast over the full modality symbol set: ",
        "z = (x - median(x)) / mad(x), using stats::mad default constant 1.4826; ",
        "stored centers/scales invert to raw values"
      ),
      modality_note = paste0(
        "three modalities only: snRNAseq, GeoMx, bulk protein-group; the phosphosite ",
        "parent-gene layer is a within-assay alternate of the same TiO2 assay"
      ),
      sources = c(
        snRNAseq = "pb_de_microglia$top joined to symbol_map by ENSMUSG id",
        GeoMx = "geomx_de$primary$top keyed by symbol",
        bulk = "proteome_de_24m$top keyed by gene_first (protein-group primary view)",
        phospho_parent_alt = "phospho_de_24m$top single parent-gene rows only"
      )
    )
  )

  expected_top_names <- c(
    "modalities", "contrasts", "symbols", "raw", "standardized",
    "standardization", "index", "phospho_parent_alt", "provenance"
  )
  stopifnot(
    identical(names(result), expected_top_names),
    identical(result$modalities, c("snRNAseq", "GeoMx", "bulk")),
    identical(names(result$symbols), result$modalities),
    identical(names(result$raw), result$modalities),
    identical(names(result$standardized), result$modalities),
    identical(names(result$standardization), result$modalities),
    !any(grepl("phospho", result$modalities, ignore.case = TRUE)),
    !any(grepl("phospho", names(result$raw), ignore.case = TRUE)),
    !any(grepl("phospho", names(result$standardized), ignore.case = TRUE)),
    identical(result$index$membership[result$index$complete_case],
              setNames(rep.int(3L, length(result$index$complete_case)),
                       result$index$complete_case)),
    identical(result$phospho_parent_alt$note,
              "phosphosite -> parent-gene collapse of the SAME TiO2 assay as bulk; not a 4th modality")
  )

  for (modality in result$modalities) {
    stopifnot(
      identical(result$symbols[[modality]], integration_radix_sort(result$symbols[[modality]])),
      anyDuplicated(result$symbols[[modality]]) == 0L
    )
    for (statistic in c("logFC", "t")) {
      raw_matrix <- result$raw[[modality]][[statistic]]
      z_matrix <- result$standardized[[modality]][[statistic]]
      params <- result$standardization[[modality]][[statistic]]
      stopifnot(
        identical(dimnames(raw_matrix), list(result$contrasts, result$symbols[[modality]])),
        identical(dimnames(z_matrix), dimnames(raw_matrix)),
        all(is.finite(raw_matrix)), all(is.finite(z_matrix)),
        identical(names(params), c("center", "scale")),
        identical(names(params$center), result$contrasts),
        identical(names(params$scale), result$contrasts),
        all(params$scale > 0),
        max(abs(
          sweep(sweep(z_matrix, 1L, params$scale, "*"), 1L, params$center, "+") -
            raw_matrix
        )) < 1e-9
      )
    }
  }
  stopifnot(
    identical(result$index$union, integration_radix_sort(result$index$union)),
    identical(result$index$complete_case, integration_radix_sort(result$index$complete_case)),
    identical(result$index$at_least_two, integration_radix_sort(result$index$at_least_two)),
    all(vapply(result$index$pairwise, function(x) {
      identical(x, integration_radix_sort(x)) && anyDuplicated(x) == 0L
    }, logical(1))),
    identical(result$phospho_parent_alt$symbols,
              integration_radix_sort(result$phospho_parent_alt$symbols)),
    anyDuplicated(result$phospho_parent_alt$symbols) == 0L,
    identical(dimnames(result$phospho_parent_alt$logFC),
              list(result$contrasts, result$phospho_parent_alt$symbols)),
    identical(dimnames(result$phospho_parent_alt$t),
              list(result$contrasts, result$phospho_parent_alt$symbols)),
    !integration_contains_parent(result),
    as.numeric(object.size(result)) < 25 * 1024^2
  )
  result
}


integration_rank_signal <- function(
    singular_values, n_rows, n_cols, rank_cap = 2L,
    relative_cutoff = 0.10) {
  stopifnot(
    is.numeric(singular_values), length(singular_values) > 0L,
    all(is.finite(singular_values)), all(singular_values >= 0),
    is.numeric(relative_cutoff), length(relative_cutoff) == 1L,
    relative_cutoff > 0, relative_cutoff < 1,
    length(rank_cap) == 1L, rank_cap >= 0L
  )
  leading <- singular_values[[1L]]
  numerical_tolerance <- max(n_rows, n_cols) * .Machine$double.eps *
    max(1, leading)
  threshold <- max(relative_cutoff * leading, numerical_tolerance)
  rank_uncapped <- if (leading <= numerical_tolerance) {
    0L
  } else {
    as.integer(sum(singular_values >= threshold))
  }
  rank <- as.integer(min(rank_cap, rank_uncapped))
  total_energy <- sum(singular_values^2)
  energy_fraction <- if (total_energy > 0) {
    singular_values^2 / total_energy
  } else {
    rep.int(0, length(singular_values))
  }
  names(singular_values) <- paste0("sv", seq_along(singular_values))
  names(energy_fraction) <- names(singular_values)

  list(
    singular_values = singular_values,
    relative_singular_values = singular_values / max(leading, numerical_tolerance),
    energy_fraction = energy_fraction,
    cumulative_energy = cumsum(energy_fraction),
    threshold = threshold,
    relative_cutoff = relative_cutoff,
    numerical_tolerance = numerical_tolerance,
    rank_uncapped = rank_uncapped,
    rank = rank,
    hard_cap = as.integer(rank_cap),
    cap_applied = rank_uncapped > rank,
    rule = paste0(
      "retain singular values >= 0.10 times the leading singular value, ",
      "then apply the requested hard rank cap"
    )
  )
}

integration_canonicalize_axes <- function(x) {
  stopifnot(is.matrix(x), is.numeric(x), all(is.finite(x)))
  if (ncol(x) == 0L) return(x)
  for (j in seq_len(ncol(x))) {
    pivot <- which.max(abs(x[, j]))
    if (x[pivot, j] < 0) x[, j] <- -x[, j]
  }
  x
}

integration_random_joint_threshold <- function(
    n_objects, ranks, joint_cap, seed, n_random = 256L,
    probability = 0.95) {
  stopifnot(
    n_objects > 0L, length(ranks) >= 2L, all(ranks > 0L),
    joint_cap > 0L, joint_cap <= min(ranks),
    length(seed) == 1L, is.finite(seed),
    n_random >= 32L, probability > 0, probability < 1
  )
  set.seed(as.integer(seed))
  maximum_stacked <- numeric(n_random)
  maximum_minimum_alignment <- numeric(n_random)

  for (iteration in seq_len(n_random)) {
    random_bases <- lapply(ranks, function(rank) {
      candidate <- matrix(stats::rnorm(n_objects * rank), n_objects, rank)
      qr.Q(qr(candidate), complete = FALSE)[, seq_len(rank), drop = FALSE]
    })
    stacked <- do.call(rbind, lapply(random_bases, t))
    stacked_svd <- svd(stacked, nu = 0L, nv = joint_cap)
    candidate_axes <- stacked_svd$v[, seq_len(joint_cap), drop = FALSE]
    alignment <- vapply(random_bases, function(basis) {
      colSums(crossprod(basis, candidate_axes)^2)
    }, numeric(joint_cap))
    if (joint_cap == 1L) {
      alignment <- matrix(alignment, nrow = 1L)
    }
    maximum_stacked[[iteration]] <- stacked_svd$d[[1L]]
    maximum_minimum_alignment[[iteration]] <- max(apply(alignment, 1L, min))
  }

  empirical_quantile <- function(x) {
    ordered <- sort(x, method = "radix")
    index <- min(length(ordered), ceiling((length(ordered) + 1L) * probability))
    ordered[[index]]
  }

  list(
    stacked_singular_value = empirical_quantile(maximum_stacked),
    minimum_alignment = empirical_quantile(maximum_minimum_alignment),
    seed = as.integer(seed),
    draws = as.integer(n_random),
    probability = probability,
    rule = paste0(
      "fixed-seed independent Gaussian-QR gene subspaces; exact empirical ",
      "order-statistic threshold"
    )
  )
}

integration_joint_basis <- function(
    signal_bases, signal_diagnostics, seed = 8202L, n_random = 256L,
    fallback_angle_degrees = 45) {
  stopifnot(
    is.list(signal_bases), length(signal_bases) >= 2L,
    identical(names(signal_bases), names(signal_diagnostics)),
    all(vapply(signal_bases, is.matrix, logical(1)))
  )
  ranks <- vapply(signal_bases, ncol, integer(1))
  n_objects <- unique(vapply(signal_bases, nrow, integer(1)))
  stopifnot(length(n_objects) == 1L)
  joint_cap <- min(ranks)
  empty_basis <- matrix(
    numeric(0), nrow = n_objects, ncol = 0L,
    dimnames = list(rownames(signal_bases[[1L]]), character())
  )
  if (joint_cap == 0L) {
    return(list(
      basis = empty_basis,
      rank = 0L,
      diagnostics = list(
        stacked_singular_values = numeric(), candidates = data.frame(),
        thresholds = list(
          seed = as.integer(seed), random_draws = as.integer(n_random),
          fallback_fired = FALSE,
          reason = "at least one block selected signal rank zero"
        )
      )
    ))
  }

  stacked <- do.call(rbind, lapply(signal_bases, t))
  stacked_svd <- svd(stacked, nu = 0L, nv = joint_cap)
  candidate_axes <- integration_canonicalize_axes(
    stacked_svd$v[, seq_len(joint_cap), drop = FALSE]
  )
  alignment <- vapply(signal_bases, function(basis) {
    colSums(crossprod(basis, candidate_axes)^2)
  }, numeric(joint_cap))
  if (joint_cap == 1L) alignment <- matrix(alignment, nrow = 1L)
  colnames(alignment) <- names(signal_bases)
  rownames(alignment) <- paste0("joint_candidate_", seq_len(joint_cap))

  random_threshold <- integration_random_joint_threshold(
    n_objects, ranks, joint_cap, seed, n_random
  )
  wedin <- lapply(seq_along(signal_diagnostics), function(i) {
    diagnostic <- signal_diagnostics[[i]]
    rank <- ranks[[i]]
    values <- unname(diagnostic$singular_values)
    signal_value <- values[[rank]]
    noise_value <- if (rank < length(values)) values[[rank + 1L]] else 0
    gap <- signal_value - noise_value
    raw_bound <- noise_value / max(gap, diagnostic$numerical_tolerance)
    bound <- min(1, raw_bound)
    data.frame(
      modality = names(signal_diagnostics)[[i]],
      signal_singular_value = signal_value,
      next_singular_value = noise_value,
      spectral_gap = gap,
      perturbation_bound = bound,
      angle_degrees = asin(bound) * 180 / pi,
      degenerate = !is.finite(raw_bound) || raw_bound >= 1,
      stringsAsFactors = FALSE
    )
  })
  wedin <- do.call(rbind, wedin)
  fallback_fired <- any(wedin$degenerate)
  fallback_alignment <- cos(fallback_angle_degrees * pi / 180)^2

  if (fallback_fired) {
    block_alignment_threshold <- rep.int(fallback_alignment, length(ranks))
    names(block_alignment_threshold) <- names(ranks)
    wedin_stacked_threshold <- NA_real_
    operative_stacked_threshold <- max(
      random_threshold$stacked_singular_value,
      sqrt(length(ranks) * fallback_alignment)
    )
    operative_rule <- paste0(
      "Wedin diagnostic degenerated; fixed ", fallback_angle_degrees,
      "-degree all-block principal-angle cutoff is operative"
    )
  } else {
    block_alignment_threshold <- 1 - wedin$perturbation_bound^2
    names(block_alignment_threshold) <- wedin$modality
    wedin_stacked_threshold <- sqrt(
      max(0, length(ranks) - sum(wedin$perturbation_bound^2))
    )
    operative_stacked_threshold <- max(
      random_threshold$stacked_singular_value,
      wedin_stacked_threshold
    )
    operative_rule <- paste0(
      "fixed-seed random-direction and Wedin perturbation thresholds are ",
      "jointly operative"
    )
  }

  minimum_alignment <- apply(alignment, 1L, min)
  clipped_alignment <- alignment
  clipped_alignment[] <- pmin(1, pmax(0, alignment))
  maximum_angle <- apply(
    acos(sqrt(clipped_alignment)) * 180 / pi,
    1L, max
  )
  comparison_tolerance <- 1e-8
  stacked_pass <- stacked_svd$d[seq_len(joint_cap)] >=
    operative_stacked_threshold - comparison_tolerance
  random_pass <- minimum_alignment >=
    random_threshold$minimum_alignment - comparison_tolerance
  block_pass <- vapply(seq_len(joint_cap), function(j) {
    all(alignment[j, ] >= block_alignment_threshold - comparison_tolerance)
  }, logical(1))
  pass <- stacked_pass & random_pass & block_pass
  selected <- cumprod(as.integer(pass)) == 1L
  joint_rank <- as.integer(sum(selected))
  joint_basis <- if (joint_rank > 0L) {
    candidate_axes[, seq_len(joint_rank), drop = FALSE]
  } else {
    empty_basis
  }
  colnames(joint_basis) <- if (joint_rank > 0L) {
    paste0("joint_", seq_len(joint_rank))
  } else {
    character()
  }

  candidates <- data.frame(
    axis = paste0("joint_candidate_", seq_len(joint_cap)),
    stacked_singular_value = stacked_svd$d[seq_len(joint_cap)],
    stacked_alignment_sum = stacked_svd$d[seq_len(joint_cap)]^2,
    minimum_block_alignment = minimum_alignment,
    maximum_principal_angle_degrees = maximum_angle,
    stacked_threshold_pass = stacked_pass,
    random_alignment_pass = random_pass,
    all_block_angle_pass = block_pass,
    selected = selected,
    stringsAsFactors = FALSE
  )

  list(
    basis = joint_basis,
    rank = joint_rank,
    diagnostics = list(
      stacked_singular_values = stacked_svd$d,
      block_alignment = alignment,
      candidates = candidates,
      thresholds = list(
        random = random_threshold,
        wedin = wedin,
        fallback_fired = fallback_fired,
        fallback_angle_degrees = fallback_angle_degrees,
        block_alignment = block_alignment_threshold,
        wedin_stacked_singular_value = wedin_stacked_threshold,
        operative_stacked_singular_value = operative_stacked_threshold,
        comparison_tolerance = comparison_tolerance,
        operative_rule = operative_rule
      )
    )
  )
}

integration_ajive_decompose <- function(
    blocks, seed = 8202L, n_random = 256L, rank_cap = 2L,
    keep_components = FALSE) {
  stopifnot(
    is.list(blocks), length(blocks) >= 2L, !is.null(names(blocks)),
    all(vapply(blocks, is.matrix, logical(1))),
    all(vapply(blocks, is.numeric, logical(1))),
    all(vapply(blocks, function(x) all(is.finite(x)), logical(1)))
  )
  dimensions <- lapply(blocks, dim)
  stopifnot(
    length(unique(vapply(dimensions, `[[`, integer(1), 1L))) == 1L,
    length(unique(vapply(dimensions, `[[`, integer(1), 2L))) == 1L,
    nrow(blocks[[1L]]) == 5L,
    all(vapply(blocks, function(x) {
      identical(rownames(x), rownames(blocks[[1L]])) &&
        identical(colnames(x), colnames(blocks[[1L]]))
    }, logical(1)))
  )

  block_svds <- lapply(blocks, function(x) {
    svd(x, nu = 0L, nv = min(dim(x)))
  })
  signal_diagnostics <- Map(function(decomposition, x) {
    integration_rank_signal(
      decomposition$d, nrow(x), ncol(x), rank_cap = rank_cap
    )
  }, block_svds, blocks)
  names(signal_diagnostics) <- names(blocks)
  signal_bases <- Map(function(decomposition, diagnostic, x) {
    rank <- diagnostic$rank
    basis <- if (rank > 0L) {
      decomposition$v[, seq_len(rank), drop = FALSE]
    } else {
      matrix(numeric(0), nrow = ncol(x), ncol = 0L)
    }
    basis <- integration_canonicalize_axes(basis)
    rownames(basis) <- colnames(x)
    colnames(basis) <- if (rank > 0L) {
      paste0("signal_", seq_len(rank))
    } else {
      character()
    }
    basis
  }, block_svds, signal_diagnostics, blocks)
  names(signal_bases) <- names(blocks)

  joint <- integration_joint_basis(
    signal_bases, signal_diagnostics, seed = seed, n_random = n_random
  )
  joint_basis <- joint$basis
  joint_rank <- joint$rank
  stopifnot(joint_rank <= min(vapply(signal_diagnostics, `[[`, integer(1), "rank")))

  block_results <- Map(function(x, signal_diagnostic) {
    joint_loadings <- x %*% joint_basis
    joint_component <- if (joint_rank > 0L) {
      joint_loadings %*% t(joint_basis)
    } else {
      matrix(0, nrow(x), ncol(x), dimnames = dimnames(x))
    }
    dimnames(joint_component) <- dimnames(x)
    orthogonal_residual <- x - joint_component
    residual_svd <- svd(
      orthogonal_residual, nu = 0L, nv = min(dim(orthogonal_residual))
    )
    individual_budget <- max(
      0L,
      min(
        signal_diagnostic$rank - joint_rank,
        nrow(x) - joint_rank - 1L
      )
    )
    individual_diagnostic <- integration_rank_signal(
      residual_svd$d, nrow(x), ncol(x), rank_cap = individual_budget
    )
    individual_rank <- individual_diagnostic$rank
    individual_basis <- if (individual_rank > 0L) {
      candidate <- residual_svd$v[, seq_len(individual_rank), drop = FALSE]
      if (joint_rank > 0L) {
        candidate <- candidate - joint_basis %*% crossprod(joint_basis, candidate)
      }
      candidate_qr <- qr(candidate)
      stopifnot(candidate_qr$rank == individual_rank)
      qr.Q(candidate_qr, complete = FALSE)[, seq_len(individual_rank), drop = FALSE]
    } else {
      matrix(numeric(0), nrow = ncol(x), ncol = 0L)
    }
    individual_basis <- integration_canonicalize_axes(individual_basis)
    rownames(individual_basis) <- colnames(x)
    colnames(individual_basis) <- if (individual_rank > 0L) {
      paste0("individual_", seq_len(individual_rank))
    } else {
      character()
    }
    individual_loadings <- orthogonal_residual %*% individual_basis
    individual_component <- if (individual_rank > 0L) {
      individual_loadings %*% t(individual_basis)
    } else {
      matrix(0, nrow(x), ncol(x), dimnames = dimnames(x))
    }
    dimnames(individual_component) <- dimnames(x)
    residual_component <- x - joint_component - individual_component
    dimnames(residual_component) <- dimnames(x)

    total_sum_squares <- sum(x^2)
    component_sum_squares <- c(
      joint = sum(joint_component^2),
      individual = sum(individual_component^2),
      residual = sum(residual_component^2)
    )
    fractions <- component_sum_squares / total_sum_squares
    reconstruction_error <- max(abs(
      x - (joint_component + individual_component + residual_component)
    ))
    basis_cross <- if (joint_rank > 0L && individual_rank > 0L) {
      max(abs(crossprod(joint_basis, individual_basis)))
    } else {
      0
    }
    frobenius_cross <- c(
      joint_individual = sum(joint_component * individual_component),
      joint_residual = sum(joint_component * residual_component),
      individual_residual = sum(individual_component * residual_component)
    )
    relative_frobenius_cross <- frobenius_cross / max(1, total_sum_squares)
    pythagorean_error <- abs(total_sum_squares - sum(component_sum_squares)) /
      max(1, total_sum_squares)

    stopifnot(
      reconstruction_error < 1e-8,
      basis_cross < 1e-8,
      max(abs(relative_frobenius_cross)) < 1e-8,
      pythagorean_error < 1e-8,
      all(fractions >= -1e-12), all(fractions <= 1 + 1e-12),
      abs(sum(fractions) - 1) < 1e-8,
      joint_rank + individual_rank <= nrow(x) - 1L
    )

    out <- list(
      individual_rank = individual_rank,
      residual_rank_budget = as.integer(nrow(x) - joint_rank - individual_rank),
      individual_diagnostic = individual_diagnostic,
      joint_loadings = joint_loadings,
      variance = data.frame(
        component = names(component_sum_squares),
        sum_squares = unname(component_sum_squares),
        fraction = unname(fractions),
        stringsAsFactors = FALSE
      ),
      oracles = list(
        reconstruction_error = reconstruction_error,
        basis_cross = basis_cross,
        frobenius_cross = frobenius_cross,
        relative_frobenius_cross = relative_frobenius_cross,
        pythagorean_error = pythagorean_error
      )
    )
    if (keep_components) {
      out$individual_basis <- individual_basis
      out$components <- list(
        joint = joint_component,
        individual = individual_component,
        residual = residual_component
      )
    }
    out
  }, blocks, signal_diagnostics)
  names(block_results) <- names(blocks)

  variance <- do.call(rbind, lapply(names(block_results), function(modality) {
    data.frame(
      modality = modality,
      block_results[[modality]]$variance,
      row.names = NULL,
      stringsAsFactors = FALSE
    )
  }))
  ranks <- data.frame(
    modality = names(blocks),
    signal_rank_uncapped = vapply(
      signal_diagnostics, `[[`, integer(1), "rank_uncapped"
    ),
    signal_rank = vapply(signal_diagnostics, `[[`, integer(1), "rank"),
    joint_rank = rep.int(joint_rank, length(blocks)),
    individual_rank = vapply(block_results, `[[`, integer(1), "individual_rank"),
    residual_rank_budget = vapply(
      block_results, `[[`, integer(1), "residual_rank_budget"
    ),
    stringsAsFactors = FALSE
  )
  stopifnot(
    all(ranks$signal_rank <= rank_cap),
    all(ranks$joint_rank <= ranks$signal_rank),
    all(ranks$joint_rank + ranks$individual_rank <= nrow(blocks[[1L]]) - 1L),
    all(ranks$residual_rank_budget >= 1L)
  )

  component_order <- c("joint", "individual", "residual")
  across_block <- data.frame(
    component = component_order,
    mean_fraction = vapply(component_order, function(component) {
      mean(variance$fraction[variance$component == component])
    }, numeric(1)),
    energy_weighted_fraction = vapply(component_order, function(component) {
      sum(variance$sum_squares[variance$component == component]) /
        sum(variance$sum_squares)
    }, numeric(1)),
    stringsAsFactors = FALSE
  )
  stopifnot(abs(sum(across_block$energy_weighted_fraction) - 1) < 1e-8)

  joint_loadings <- lapply(block_results, `[[`, "joint_loadings")
  for (modality in names(joint_loadings)) {
    dimnames(joint_loadings[[modality]]) <- list(
      rownames(blocks[[modality]]),
      if (joint_rank > 0L) paste0("joint_", seq_len(joint_rank)) else character()
    )
  }
  diagnostics <- list(
    signal = signal_diagnostics,
    joint = joint$diagnostics,
    individual = lapply(block_results, `[[`, "individual_diagnostic"),
    oracles = lapply(block_results, `[[`, "oracles"),
    rank_rule = signal_diagnostics[[1L]]$rule
  )

  result <- list(
    variance = variance,
    across_block = across_block,
    ranks = ranks,
    joint_gene_scores = joint_basis,
    joint_contrast_loadings = joint_loadings,
    diagnostics = diagnostics
  )
  if (keep_components) {
    result$components <- lapply(block_results, `[[`, "components")
    result$individual_bases <- lapply(block_results, `[[`, "individual_basis")
  }
  result
}

integration_ajive_fixture <- function() {
  set.seed(8203L)
  modalities <- c("snRNAseq", "GeoMx", "bulk")
  contrasts <- integration_contrast_names()
  n_objects <- 240L
  genes <- paste0("fixture_gene_", seq_len(n_objects))
  gene_basis <- qr.Q(
    qr(matrix(stats::rnorm(n_objects * 5L), n_objects, 5L)),
    complete = FALSE
  )
  dimnames(gene_basis) <- list(genes, paste0("gene_axis_", seq_len(5L)))
  loading_seeds <- list(
    cbind(c(4, 2, -1, 3, 1), c(1, -3, 4, 2, -2)),
    cbind(c(-2, 5, 3, 1, -1), c(4, 1, -2, 3, 2)),
    cbind(c(5, -1, 2, 4, 1), c(-2, 3, 1, -1, 5))
  )
  loading_bases <- lapply(loading_seeds, function(x) {
    qr.Q(qr(x), complete = FALSE)[, 1:2, drop = FALSE]
  })

  rank_one_blocks <- setNames(lapply(seq_along(modalities), function(i) {
    signal <- 8 * tcrossprod(loading_bases[[i]][, 1L], gene_basis[, 1L]) +
      5 * tcrossprod(loading_bases[[i]][, 2L], gene_basis[, i + 1L])
    noise <- matrix(stats::rnorm(5L * n_objects), 5L, n_objects) * 1e-10
    out <- signal + noise
    dimnames(out) <- list(contrasts, genes)
    out
  }), modalities)
  rank_one <- integration_ajive_decompose(
    rank_one_blocks, seed = 8204L, n_random = 128L,
    keep_components = TRUE
  )
  rank_one_repeat <- integration_ajive_decompose(
    rank_one_blocks, seed = 8204L, n_random = 128L,
    keep_components = TRUE
  )
  rank_one_angle <- acos(pmin(1, svd(
    crossprod(gene_basis[, 1L, drop = FALSE], rank_one$joint_gene_scores),
    nu = 0L, nv = 0L
  )$d))
  rank_one_error <- max(vapply(seq_along(modalities), function(i) {
    planted <- 8 * tcrossprod(
      loading_bases[[i]][, 1L], gene_basis[, 1L]
    )
    recovered <- rank_one$components[[i]]$joint
    sqrt(sum((recovered - planted)^2) / sum(planted^2))
  }, numeric(1)))
  stopifnot(
    identical(rank_one$ranks$signal_rank, rep.int(2L, 3L)),
    identical(unique(rank_one$ranks$joint_rank), 1L),
    identical(rank_one$ranks$individual_rank, rep.int(1L, 3L)),
    max(rank_one_angle) < 1e-6,
    rank_one_error < 1e-6,
    identical(rank_one$ranks, rank_one_repeat$ranks),
    isTRUE(all.equal(
      rank_one$joint_gene_scores, rank_one_repeat$joint_gene_scores,
      tolerance = 0
    ))
  )

  rank_two_blocks <- setNames(lapply(seq_along(modalities), function(i) {
    signal <- 9 * tcrossprod(loading_bases[[i]][, 1L], gene_basis[, 1L]) +
      6 * tcrossprod(loading_bases[[i]][, 2L], gene_basis[, 2L])
    noise <- matrix(stats::rnorm(5L * n_objects), 5L, n_objects) * 1e-10
    out <- signal + noise
    dimnames(out) <- list(contrasts, genes)
    out
  }), modalities)
  rank_two <- integration_ajive_decompose(
    rank_two_blocks, seed = 8205L, n_random = 128L,
    keep_components = TRUE
  )
  rank_two_angles <- acos(pmin(1, svd(
    crossprod(gene_basis[, 1:2, drop = FALSE], rank_two$joint_gene_scores),
    nu = 0L, nv = 0L
  )$d))
  rank_two_error <- max(vapply(seq_along(modalities), function(i) {
    planted <- 9 * tcrossprod(
      loading_bases[[i]][, 1L], gene_basis[, 1L]
    ) + 6 * tcrossprod(
      loading_bases[[i]][, 2L], gene_basis[, 2L]
    )
    recovered <- rank_two$components[[i]]$joint
    sqrt(sum((recovered - planted)^2) / sum(planted^2))
  }, numeric(1)))
  stopifnot(
    identical(rank_two$ranks$signal_rank, rep.int(2L, 3L)),
    identical(unique(rank_two$ranks$joint_rank), 2L),
    identical(rank_two$ranks$individual_rank, rep.int(0L, 3L)),
    max(rank_two_angles) < 1e-6,
    rank_two_error < 1e-6
  )

  list(
    rank_one = list(
      planted_rank = 1L, recovered_rank = unique(rank_one$ranks$joint_rank),
      maximum_principal_angle = max(rank_one_angle),
      maximum_relative_joint_error = rank_one_error,
      deterministic_repeat = TRUE
    ),
    rank_two = list(
      planted_rank = 2L, recovered_rank = unique(rank_two$ranks$joint_rank),
      maximum_principal_angle = max(rank_two_angles),
      maximum_relative_joint_error = rank_two_error
    )
  )
}

build_integration_decomposition <- function(integration_substrate) {
  seed <- 8202L
  set.seed(seed)
  stopifnot(
    is.list(integration_substrate),
    identical(integration_substrate$modalities, c("snRNAseq", "GeoMx", "bulk")),
    identical(integration_substrate$contrasts, integration_contrast_names()),
    is.list(integration_substrate$standardized),
    identical(names(integration_substrate$standardized),
              integration_substrate$modalities)
  )
  genes <- integration_substrate$index$complete_case
  stopifnot(
    identical(genes, integration_radix_sort(genes)),
    anyDuplicated(genes) == 0L,
    length(genes) == 3109L
  )

  fixture <- integration_ajive_fixture()
  stopifnot(
    fixture$rank_one$recovered_rank == fixture$rank_one$planted_rank,
    fixture$rank_two$recovered_rank == fixture$rank_two$planted_rank
  )

  blocks_for <- function(statistic) {
    out <- setNames(lapply(integration_substrate$modalities, function(modality) {
      source <- integration_substrate$standardized[[modality]][[statistic]]
      stopifnot(
        is.matrix(source), identical(rownames(source), integration_contrast_names()),
        all(genes %in% colnames(source))
      )
      block <- source[, genes, drop = FALSE]
      stopifnot(
        identical(dim(block), c(5L, length(genes))),
        identical(rownames(block), integration_contrast_names()),
        identical(colnames(block), genes),
        all(is.finite(block))
      )
      block
    }), integration_substrate$modalities)
    stopifnot(all(vapply(out, function(x) {
      identical(colnames(x), genes)
    }, logical(1))))
    out
  }

  primary_internal <- integration_ajive_decompose(
    blocks_for("logFC"), seed = seed, n_random = 256L
  )
  sensitivity_internal <- integration_ajive_decompose(
    blocks_for("t"), seed = seed, n_random = 256L
  )

  primary <- c(
    list(statistic = "standardized_logFC"),
    primary_internal
  )
  sensitivity <- list(
    statistic = "standardized_moderated_t",
    variance = sensitivity_internal$variance,
    across_block = sensitivity_internal$across_block,
    ranks = sensitivity_internal$ranks,
    diagnostics = sensitivity_internal$diagnostics
  )
  result <- list(
    primary = primary,
    sensitivity = sensitivity,
    provenance = list(
      universe = "complete_case",
      n_genes = as.integer(length(genes)),
      modalities = integration_substrate$modalities,
      contrasts = integration_substrate$contrasts,
      orientation = "5 contrasts (rows) x shared gene symbols (columns)",
      primary = "per-modality robust-z standardized logFC",
      sensitivity = "per-modality robust-z standardized moderated t",
      method = paste0(
        "pure-R AJIVE-style shared gene-subspace projection using base svd/qr; ",
        "Frobenius-energy joint/individual/residual partition"
      ),
      rank_selection = paste0(
        "deterministic relative singular-value rule with r_k <= 2; joint rank ",
        "uses fixed-seed random-direction plus Wedin diagnostic, with a ",
        "45-degree deterministic principal-angle fallback"
      ),
      rank_budget = "r_J + r_I,k <= 4, leaving at least one of five row dimensions",
      missingness = "complete-case genes only; no imputation",
      interpretation = paste0(
        "descriptive point-estimate decomposition only; genes are not exchangeable ",
        "sampling units and no calibrated cross-modality p-values are produced"
      ),
      pairing = "no per-animal or aliquot pairing is attempted",
      r_jive_reference = paste0(
        "one-time implementation cross-check only, not a dependency; accept <=0.10 ",
        "absolute component-fraction difference and <=0.05 structured-fraction difference"
      ),
      seed = seed
    )
  )

  stopifnot(
    identical(names(result), c("primary", "sensitivity", "provenance")),
    identical(rownames(result$primary$joint_gene_scores), genes),
    ncol(result$primary$joint_gene_scores) ==
      unique(result$primary$ranks$joint_rank),
    identical(names(result$primary$joint_contrast_loadings),
              integration_substrate$modalities),
    all(vapply(result$primary$joint_contrast_loadings, function(x) {
      identical(rownames(x), integration_contrast_names())
    }, logical(1))),
    !integration_contains_parent(result),
    as.numeric(object.size(result)) < 25 * 1024^2
  )
  result
}
