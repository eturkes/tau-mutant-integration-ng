# P8 integration is restricted to the P7.3-compliant symbol x five-contrast
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
  columns <- lapply(contrasts, function(contrast) {
    tab <- top[[contrast]]
    idx <- match(representatives$feature, as.character(tab[[feature_col]]))
    stopifnot(!anyNA(idx))
    values <- tab[[statistic]][idx]
    stopifnot(is.numeric(values), all(is.finite(values)))
    values
  })
  out <- do.call(cbind, columns)
  dimnames(out) <- list(representatives$symbol, contrasts)

  # Exact, tolerance-zero reconstruction from every representative source row.
  for (j in seq_along(contrasts)) {
    tab <- top[[contrasts[[j]]]]
    idx <- match(representatives$feature, as.character(tab[[feature_col]]))
    stopifnot(identical(as.vector(out[, j]), as.vector(tab[[statistic]][idx])))
  }
  out
}

robust_z <- function(x) {
  stopifnot(
    is.matrix(x), is.numeric(x), nrow(x) > 0L, ncol(x) > 0L,
    !is.null(rownames(x)), !is.null(colnames(x)), all(is.finite(x))
  )
  center <- apply(x, 2L, stats::median)
  scale <- vapply(
    seq_len(ncol(x)),
    function(j) stats::mad(x[, j], center = center[[j]]),
    numeric(1)
  )
  names(center) <- colnames(x)
  names(scale) <- colnames(x)
  stopifnot(all(is.finite(center)), all(is.finite(scale)), all(scale > 0))

  z <- sweep(sweep(x, 2L, center, "-"), 2L, scale, "/")
  dimnames(z) <- dimnames(x)
  reconstructed <- sweep(sweep(z, 2L, scale, "*"), 2L, center, "+")
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
            sweep(rz$z, 2L, rz$scale, "*"),
            2L, rz$center, "+"
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
        identical(dimnames(raw_matrix), list(result$symbols[[modality]], result$contrasts)),
        identical(dimnames(z_matrix), dimnames(raw_matrix)),
        all(is.finite(raw_matrix)), all(is.finite(z_matrix)),
        identical(names(params), c("center", "scale")),
        identical(names(params$center), result$contrasts),
        identical(names(params$scale), result$contrasts),
        all(params$scale > 0),
        max(abs(
          sweep(sweep(z_matrix, 2L, params$scale, "*"), 2L, params$center, "+") -
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
              list(result$phospho_parent_alt$symbols, result$contrasts)),
    identical(dimnames(result$phospho_parent_alt$t),
              list(result$phospho_parent_alt$symbols, result$contrasts)),
    !integration_contains_parent(result),
    as.numeric(object.size(result)) < 25 * 1024^2
  )
  result
}
