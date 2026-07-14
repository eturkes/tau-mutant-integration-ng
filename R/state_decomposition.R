# P6 state composition-versus-regulation substrate. The heavy annotated Seurat is read once;
# this module emits raw Homeostatic/DAM pseudobulks + exact unit/score/map audits only. No
# Seurat/S4 object, fitted model, or cell-level frame may cross the target boundary.

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

# Recursive isolation guard: target payloads remain base vectors/matrices/data frames/lists.
# Any Seurat/S4/environment/external-pointer reachability means the heavy parent leaked through.
state_substrate_contains_parent <- function(x) {
  if (inherits(x, c("Seurat", "Assay", "Assay5", "DimReduc", "Neighbor")) ||
      isS4(x) || is.environment(x) || typeof(x) == "externalptr") return(TRUE)
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
