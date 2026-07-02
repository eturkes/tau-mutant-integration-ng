# P3 mechanism helpers: symbol-ranked RNA matrices, decoupleR activity wrappers, and
# fingerprinted OmniPath prior loaders. No package is loaded at source time; every external API
# is required inside the caller so tests/gates fail at the exact contract boundary.

mechanism_tz_preflight <- function(default = "UTC") {
  if (!nzchar(Sys.getenv("TZ"))) Sys.setenv(TZ = default)
  Sys.getenv("TZ")
}

.ensure_dir <- function(path, label = "directory") {
  if (!dir.exists(path)) {
    ok <- dir.create(path, recursive = TRUE, showWarnings = FALSE)
    if (!ok && !dir.exists(path)) stop("failed to create ", label, ": ", path, call. = FALSE)
  }
  if (!dir.exists(path)) stop(label, " does not exist: ", path, call. = FALSE)
  normalizePath(path, mustWork = TRUE)
}

set_mechanism_prior_cache <- function(path = file.path("storage", "cache", "omnipath")) {
  mechanism_tz_preflight()
  if (!requireNamespace("OmnipathR", quietly = TRUE)) stop("OmnipathR is not installed", call. = FALSE)
  cache_path <- .ensure_dir(path, "OmniPath cache")
  OmnipathR::omnipath_set_cachedir(cache_path)
  list(cache_dir = cache_path, tz = Sys.getenv("TZ"))
}

add_symbol_to_top <- function(top_df, symbol_map, gene_col = "gene") {
  stopifnot(is.data.frame(top_df), is.data.frame(symbol_map))
  if (!gene_col %in% names(top_df)) stop("missing gene column: ", gene_col, call. = FALSE)
  stopifnot(all(c("ensembl", "symbol") %in% names(symbol_map)))
  out <- top_df
  map <- symbol_map[!is.na(symbol_map$ensembl) & !is.na(symbol_map$symbol) &
                      symbol_map$ensembl != "" & symbol_map$symbol != "", c("ensembl", "symbol")]
  map <- map[!duplicated(map$ensembl), , drop = FALSE]
  out$symbol <- map$symbol[match(as.character(out[[gene_col]]), map$ensembl)]
  keep <- !is.na(out$symbol) & out$symbol != ""
  attr(out, "symbol_mapping") <- list(
    gene_col = gene_col,
    n_input = nrow(top_df),
    n_mapped = sum(keep),
    n_dropped = sum(!keep),
    n_map_rows = nrow(map)
  )
  out[keep, , drop = FALSE]
}

mechanism_contrasts <- function() {
  pairwise <- names(contrast_definitions)
  c("tau_alone", setdiff(pairwise, "tau_alone"), "interaction")
}

extract_rank_matrix <- function(top_list, symbol_map, stat_col = "t") {
  stopifnot(is.list(top_list), length(top_list) >= 1L, !is.null(names(top_list)), nzchar(stat_col))
  present <- names(top_list)
  if (any(!nzchar(present)) || anyDuplicated(present)) stop("top_list must have unique nonblank names", call. = FALSE)
  if (exists("contrast_definitions", inherits = TRUE)) {
    present <- intersect(mechanism_contrasts(), present)
  }
  if (!length(present)) stop("no named contrasts present", call. = FALSE)
  collapsed <- lapply(present, function(cn) {
    top <- add_symbol_to_top(top_list[[cn]], symbol_map)
    if (!stat_col %in% names(top)) stop("missing stat column: ", stat_col, call. = FALSE)
    stat <- as.numeric(top[[stat_col]])
    ok <- is.finite(stat)
    top <- top[ok, , drop = FALSE]
    stat <- stat[ok]
    if (!nrow(top)) stop("no finite mapped statistics for contrast: ", cn, call. = FALSE)
    idx <- split(seq_along(stat), top$symbol)
    val <- vapply(idx, function(i) {
      best <- i[which.max(abs(stat[i]))]
      stat[best]
    }, numeric(1))
    val[sort(names(val), method = "radix")]
  })
  names(collapsed) <- present
  symbols <- sort(unique(unlist(lapply(collapsed, names), use.names = FALSE)), method = "radix")
  mat <- matrix(NA_real_, nrow = length(symbols), ncol = length(present),
                dimnames = list(symbols, present))
  for (cn in present) mat[names(collapsed[[cn]]), cn] <- collapsed[[cn]]
  attr(mat, "symbol_mapping") <- lapply(top_list[present], function(x)
    attr(add_symbol_to_top(x, symbol_map), "symbol_mapping"))
  mat
}

run_decoupler_matrix <- function(mat, network, minsize = 5L) {
  stopifnot(is.matrix(mat), !is.null(rownames(mat)), !is.null(colnames(mat)),
            is.data.frame(network), all(c("source", "target", "mor") %in% names(network)),
            is.numeric(network$mor), all(network$mor %in% c(-1, 1)), minsize >= 1L)
  if (!requireNamespace("decoupleR", quietly = TRUE)) stop("decoupleR is not installed", call. = FALSE)
  res <- decoupleR::run_ulm(mat, network, .source = "source", .target = "target",
                            .mor = "mor", minsize = minsize, na.rm = TRUE)
  res <- as.data.frame(res, stringsAsFactors = FALSE)
  need <- c("statistic", "source", "condition", "score", "p_value")
  stopifnot(all(need %in% names(res)), all(is.finite(res$score)),
            all(is.finite(res$p_value) | is.na(res$p_value)))
  out <- tibble::as_tibble(res[need])
  attr(out, "has_consensus") <- FALSE
  attr(out, "method") <- "ulm"
  out
}

.normalise_prior_df <- function(x) {
  stopifnot(is.data.frame(x))
  out <- as.data.frame(x, stringsAsFactors = FALSE)
  out[] <- lapply(out, function(col) {
    if (is.factor(col)) as.character(col)
    else if (inherits(col, "integer64")) as.character(col)
    else col
  })
  out <- out[sort(names(out), method = "radix")]
  if (nrow(out)) {
    ord <- do.call(order, c(out, list(method = "radix", na.last = TRUE)))
    out <- out[ord, , drop = FALSE]
  }
  rownames(out) <- NULL
  out
}

.mechanism_pkg_versions <- function(pkgs = c("decoupleR", "OmnipathR", "digest")) {
  stats::setNames(lapply(pkgs, function(pkg) {
    if (requireNamespace(pkg, quietly = TRUE)) as.character(utils::packageVersion(pkg)) else NA_character_
  }), pkgs)
}

.sort_named_list <- function(x) {
  if (!is.list(x)) return(x)
  if (!is.null(names(x))) x <- x[sort(names(x), method = "radix")]
  lapply(x, .sort_named_list)
}

prior_fingerprint <- function(x, query) {
  if (!requireNamespace("digest", quietly = TRUE)) stop("digest is not installed", call. = FALSE)
  payload <- list(data = .normalise_prior_df(x), query = .sort_named_list(query),
                  package_versions = .sort_named_list(.mechanism_pkg_versions()))
  list(
    hash = digest::digest(payload, algo = "sha256"),
    n_rows = nrow(payload$data),
    n_cols = ncol(payload$data),
    query = query,
    package_versions = payload$package_versions
  )
}

.bool_omnipath <- function(x) {
  if (is.logical(x)) return(x)
  tolower(as.character(x)) %in% c("true", "t", "1", "yes")
}

.read_omnipath_tsv <- function(endpoint, query, cache_dir = file.path("storage", "cache", "omnipath")) {
  if (!requireNamespace("digest", quietly = TRUE)) stop("digest is not installed", call. = FALSE)
  cache_dir <- .ensure_dir(cache_dir, "OmniPath REST cache")
  q <- query[sort(names(query), method = "radix")]
  query_text <- paste(vapply(names(q), function(nm)
    paste0(utils::URLencode(nm, reserved = TRUE), "=",
           utils::URLencode(as.character(q[[nm]]), reserved = TRUE)), character(1)),
    collapse = "&")
  url <- paste0(endpoint, "?", query_text)
  key <- digest::digest(list(endpoint = endpoint, query = q), algo = "sha256")
  cached <- file.path(cache_dir, paste0("rest_", key, ".tsv"))
  if (!file.exists(cached)) {
    tmp <- tempfile(fileext = ".tsv")
    utils::download.file(url, tmp, quiet = TRUE, mode = "wb")
    ok <- file.rename(tmp, cached)
    if (!ok) {
      file.copy(tmp, cached, overwrite = TRUE)
      unlink(tmp)
    }
  }
  read.delim(cached, check.names = FALSE, stringsAsFactors = FALSE)
}

.dedupe_signed_prior <- function(out) {
  stopifnot(is.data.frame(out), all(c("source", "target", "mor") %in% names(out)))
  key <- paste(out$source, out$target, sep = "\r")
  conflict <- tapply(out$mor, key, function(v) length(unique(v)) > 1L)
  conflict_keys <- names(conflict)[conflict]
  conflict_row <- key %in% conflict_keys
  clean <- out[!conflict_row, , drop = FALSE]
  duplicate_rows <- sum(duplicated(clean[c("source", "target", "mor")]))
  clean <- clean[!duplicated(clean[c("source", "target", "mor")]), , drop = FALSE]
  rownames(clean) <- NULL
  attr(clean, "dedupe_counts") <- list(
    n_conflicting_pairs = length(conflict_keys),
    n_conflicting_rows = sum(conflict_row),
    n_duplicate_rows = duplicate_rows
  )
  clean
}

standardise_collectri_table <- function(raw, split_complexes = FALSE) {
  need <- c("source", "source_genesymbol", "target_genesymbol", "is_stimulation", "is_inhibition")
  stopifnot(is.data.frame(raw), all(need %in% names(raw)))
  src <- as.character(raw$source_genesymbol)
  is_complex <- grepl("^COMPLEX", as.character(raw$source)) | grepl("_", src, fixed = TRUE)
  if (!split_complexes) {
    src[is_complex & grepl("JUN|FOS", src, ignore.case = TRUE)] <- "AP1"
    src[is_complex & grepl("REL|NFKB", src, ignore.case = TRUE)] <- "NFKB"
  }
  sign_cols <- if (all(c("consensus_stimulation", "consensus_inhibition") %in% names(raw))) {
    c("consensus_stimulation", "consensus_inhibition")
  } else {
    c("is_stimulation", "is_inhibition")
  }
  stim <- .bool_omnipath(raw[[sign_cols[1]]])
  inhib <- .bool_omnipath(raw[[sign_cols[2]]])
  ambiguous <- (stim & inhib) | (!stim & !inhib)
  mor <- ifelse(stim & !inhib, 1, ifelse(inhib & !stim, -1, NA_real_))
  out <- data.frame(source = src, target = as.character(raw$target_genesymbol), mor = mor,
                    stringsAsFactors = FALSE)
  keep <- !is.na(out$source) & !is.na(out$target) & out$source != "" & out$target != "" &
    !is.na(out$mor) & out$mor %in% c(-1, 1)
  out <- out[keep, , drop = FALSE]
  out <- .dedupe_signed_prior(out)
  attr(out, "prior_filter_counts") <- c(list(
    sign_source = paste(sign_cols, collapse = "/"),
    n_raw = nrow(raw),
    n_ambiguous_sign = sum(ambiguous),
    n_dropped_basic = sum(!keep)
  ), attr(out, "dedupe_counts"))
  stopifnot(nrow(out) > 0L, all(out$mor %in% c(-1, 1)))
  out
}

load_collectri_mouse <- function(cache_dir = file.path("storage", "cache", "omnipath"),
                                 try_package = FALSE) {
  cache <- set_mechanism_prior_cache(cache_dir)
  query <- list(endpoint = "https://omnipathdb.org/interactions", datasets = "collectri",
                organisms = 10090L, genesymbols = 1L, format = "tsv")
  endpoint <- query$endpoint
  query_args <- query[setdiff(names(query), "endpoint")]
  pkg_error <- NULL
  raw <- NULL
  if (try_package) {
    raw <- tryCatch(decoupleR::get_collectri(organism = "mouse", split_complexes = FALSE),
                    error = function(e) { pkg_error <<- conditionMessage(e); NULL })
  }
  if (is.null(raw)) {
    raw <- .read_omnipath_tsv(endpoint, query_args, cache$cache_dir)
    net <- standardise_collectri_table(raw, split_complexes = FALSE)
    retrieval <- if (try_package) "omnipath_rest_after_omnipathr_postprocess_failure" else "omnipath_rest"
  } else {
    stopifnot(all(c("source", "target", "mor") %in% names(raw)))
    net <- as.data.frame(raw[c("source", "target", "mor")], stringsAsFactors = FALSE)
    net$mor <- as.numeric(net$mor)
    retrieval <- "decoupleR_get_collectri"
  }
  fp <- prior_fingerprint(net, c(query, list(retrieval = retrieval)))
  attr(net, "provenance") <- c(fp, list(
    cache_dir = cache$cache_dir,
    retrieval = retrieval,
    package_error = pkg_error,
    n_sources = length(unique(net$source)),
    n_targets = length(unique(net$target)),
    source_examples = head(sort(unique(net$source), method = "radix"), 8L),
    filter_counts = attr(net, "prior_filter_counts")
  ))
  net
}

standardise_ksn_table <- function(raw) {
  need <- c("enzyme_genesymbol", "substrate_genesymbol", "residue_type", "residue_offset", "modification")
  stopifnot(is.data.frame(raw), all(need %in% names(raw)))
  src <- trimws(as.character(raw$enzyme_genesymbol))
  sub <- trimws(as.character(raw$substrate_genesymbol))
  aa <- trimws(as.character(raw$residue_type))
  offset <- trimws(as.character(raw$residue_offset))
  mod <- tolower(trimws(as.character(raw$modification)))
  missing_component <- is.na(src) | src == "" | is.na(sub) | sub == "" | is.na(aa) | aa == "" |
    is.na(offset) | offset == "" | tolower(offset) %in% c("na", "nan")
  multi_gene <- grepl("[;,]", src) | grepl("[;,]", sub)
  signed <- mod %in% c("phosphorylation", "dephosphorylation")
  mor <- ifelse(mod == "phosphorylation", 1, ifelse(mod == "dephosphorylation", -1, NA_real_))
  site <- paste0(sub, "_", aa, offset)
  out <- data.frame(source = src, target = site, mor = mor,
                    stringsAsFactors = FALSE)
  keep <- !missing_component & !multi_gene & signed & !is.na(out$mor) & out$mor %in% c(-1, 1)
  out <- out[keep, , drop = FALSE]
  out <- .dedupe_signed_prior(out)
  attr(out, "prior_filter_counts") <- c(list(
    n_raw = nrow(raw),
    n_missing_component = sum(missing_component),
    n_multi_gene = sum(multi_gene & !missing_component),
    n_unsupported_modification = sum(!signed)
  ), attr(out, "dedupe_counts"))
  stopifnot(nrow(out) > 0L, all(out$mor %in% c(-1, 1)))
  out
}

load_omnipath_ksn_mouse <- function(cache_dir = file.path("storage", "cache", "omnipath"),
                                    try_package = FALSE) {
  cache <- set_mechanism_prior_cache(cache_dir)
  query <- list(endpoint = "https://omnipathdb.org/enz_sub", organisms = 10090L,
                genesymbols = 1L, format = "tsv")
  endpoint <- query$endpoint
  query_args <- query[setdiff(names(query), "endpoint")]
  pkg_error <- NULL
  raw <- NULL
  if (try_package) {
    raw <- tryCatch(OmnipathR::enzyme_substrate(organism = 10090L, genesymbols = TRUE),
                    error = function(e) { pkg_error <<- conditionMessage(e); NULL })
  }
  if (is.null(raw)) {
    raw <- .read_omnipath_tsv(endpoint, query_args, cache$cache_dir)
    net <- standardise_ksn_table(raw)
    retrieval <- if (try_package) "omnipath_rest_after_omnipathr_postprocess_failure" else "omnipath_rest"
  } else {
    net <- standardise_ksn_table(raw)
    retrieval <- "OmnipathR_enzyme_substrate"
  }
  fp <- prior_fingerprint(net, c(query, list(retrieval = retrieval)))
  attr(net, "provenance") <- c(fp, list(
    cache_dir = cache$cache_dir,
    retrieval = retrieval,
    package_error = pkg_error,
    n_sources = length(unique(net$source)),
    n_targets = length(unique(net$target)),
    source_examples = head(sort(unique(net$source), method = "radix"), 8L),
    filter_counts = attr(net, "prior_filter_counts")
  ))
  net
}

mechanism_prior_expectations <- function() {
  list(
    collectri = list(hash = "027ee57a61246bff4127d9d36807469713731de552398bb81989a06fd1bc44e6",
                     n_rows = 37096L, n_sources = 1093L, n_targets = 6010L),
    ksn = list(hash = "997b690d5efdfd8bb4424c12a29a80f5a980d8b3404025210e188281d554172d",
               n_rows = 29378L, n_sources = 1397L, n_targets = 13048L),
    ksn_coverage = list(n_matched_sites = 2250L, kinases_passing_minsize = 212L,
                        gsk3b_matched_sites = 245L)
  )
}

assert_mechanism_prior_expectations <- function(collectri = NULL, ksn = NULL, coverage = NULL,
                                                expected = mechanism_prior_expectations()) {
  check_prior <- function(x, exp, label) {
    prov <- attr(x, "provenance")
    stopifnot(is.list(prov), !is.null(prov$hash))
    vals <- list(hash = prov$hash, n_rows = nrow(x),
                 n_sources = prov$n_sources, n_targets = prov$n_targets)
    for (nm in intersect(names(exp), names(vals))) {
      if (!is.na(exp[[nm]]) && !identical(vals[[nm]], exp[[nm]])) {
        stop(label, " prior drift: ", nm, " observed=", vals[[nm]],
             " expected=", exp[[nm]], call. = FALSE)
      }
    }
    TRUE
  }
  if (!is.null(collectri)) check_prior(collectri, expected$collectri, "CollecTRI")
  if (!is.null(ksn)) check_prior(ksn, expected$ksn, "KSN")
  if (!is.null(coverage)) {
    vals <- list(n_matched_sites = coverage$n_matched_sites,
                 kinases_passing_minsize = coverage$kinases_passing_minsize,
                 gsk3b_matched_sites = coverage$gsk3b$matched_sites)
    for (nm in intersect(names(expected$ksn_coverage), names(vals))) {
      if (!is.na(expected$ksn_coverage[[nm]]) &&
          !identical(vals[[nm]], expected$ksn_coverage[[nm]])) {
        stop("KSN coverage drift: ", nm, " observed=", vals[[nm]],
             " expected=", expected$ksn_coverage[[nm]], call. = FALSE)
      }
    }
  }
  TRUE
}

phospho_site_ids <- function(phospho_tbl) {
  stopifnot(is.data.frame(phospho_tbl))
  need <- c("PG.Genes", "PTM.SiteAA", "PTM.SiteLocation")
  stopifnot(all(need %in% names(phospho_tbl)))
  gene <- trimws(as.character(phospho_tbl[["PG.Genes"]]))
  aa <- trimws(as.character(phospho_tbl[["PTM.SiteAA"]]))
  loc <- trimws(as.character(phospho_tbl[["PTM.SiteLocation"]]))
  missing_gene <- is.na(gene) | gene == ""
  multi_gene <- grepl("[;,]", gene)
  missing_site <- is.na(aa) | aa == "" | is.na(loc) | loc == "" | tolower(loc) %in% c("na", "nan")
  keep <- !missing_gene & !multi_gene & !missing_site
  ids <- paste0(gene[keep], "_", aa[keep], loc[keep])
  ids <- ids[!is.na(ids) & ids != ""]
  out <- unique(ids)
  attr(out, "phospho_site_counts") <- list(
    n_rows = nrow(phospho_tbl),
    n_missing_gene = sum(missing_gene),
    n_multi_gene = sum(multi_gene & !missing_gene),
    n_missing_site = sum(missing_site),
    n_kept_rows = length(ids),
    n_unique_sites = length(unique(ids)),
    n_duplicate_rows = length(ids) - length(unique(ids))
  )
  out
}

ksn_coverage_probe <- function(ksn, site_ids, minsize = 5L, kinase = "Gsk3b") {
  stopifnot(is.data.frame(ksn), all(c("source", "target", "mor") %in% names(ksn)),
            is.character(site_ids), minsize >= 1L)
  sites <- unique(site_ids[!is.na(site_ids) & site_ids != ""])
  matched <- ksn[ksn$target %in% sites, , drop = FALSE]
  per_kinase <- sort(table(matched$source), decreasing = TRUE)
  passing <- names(per_kinase)[per_kinase >= minsize]
  source <- unique(as.character(ksn$source))
  case <- data.frame(
    pattern = c("exact_Gsk3b", "upper_GSK3B", "contains_lowercase", "contains_underscore"),
    count = c(sum(source == "Gsk3b"), sum(source == "GSK3B"),
              sum(grepl("[a-z]", source)), sum(grepl("_", source, fixed = TRUE))),
    stringsAsFactors = FALSE
  )
  kinase_rows <- matched[matched$source == kinase, , drop = FALSE]
  list(
    n_phospho_sites = length(sites),
    n_ksn_edges = nrow(ksn),
    n_ksn_sources = length(source),
    n_matched_edges = nrow(matched),
    n_matched_sites = length(unique(matched$target)),
    minsize = minsize,
    kinases_passing_minsize = length(passing),
    passing_sources = passing,
    gsk3b = list(source_present = kinase %in% source,
                 matched_sites = length(unique(kinase_rows$target)),
                 passes_minsize = length(unique(kinase_rows$target)) >= minsize),
    source_case = case,
    top_matched = data.frame(source = names(head(per_kinase, 10L)),
                             n_sites = as.integer(head(per_kinase, 10L)),
                             row.names = NULL, stringsAsFactors = FALSE)
  )
}

# ---- P3-S2: RNA pathway / TF / NF-kB mechanism targets -------------------------------

.empty_df <- function(cols) {
  out <- as.data.frame(stats::setNames(rep(list(logical()), length(cols)), cols),
                       stringsAsFactors = FALSE)
  out[0, , drop = FALSE]
}

collect_rna_rank_matrices <- function(pb_de_microglia, pb_de_substate, symbol_map,
                                      stat_col = "t") {
  stopifnot(is.list(pb_de_microglia), is.list(pb_de_substate), is.data.frame(symbol_map),
            is.list(pb_de_microglia$top), is.list(pb_de_substate$per_substate))
  ranks <- list(
    whole_microglia = list(
      population = "whole_microglia",
      population_type = "whole",
      status = "fit",
      n_cells = pb_de_microglia$n_cells %||% NA_integer_,
      matrix = extract_rank_matrix(pb_de_microglia$top, symbol_map, stat_col = stat_col)
    )
  )
  skipped <- list()
  for (nm in names(pb_de_substate$per_substate)) {
    one <- pb_de_substate$per_substate[[nm]]
    stopifnot(is.list(one), !is.null(one$status))
    if (identical(one$status, "fit")) {
      ranks[[nm]] <- list(
        population = nm,
        population_type = "substate",
        status = "fit",
        substate = nm,
        n_cells = one$n_cells %||% NA_integer_,
        matrix = extract_rank_matrix(one$top, symbol_map, stat_col = stat_col)
      )
    } else {
      skipped[[length(skipped) + 1L]] <- data.frame(
        population = nm,
        population_type = "substate",
        status = as.character(one$status),
        n_cells = one$n_cells %||% NA_integer_,
        reason = as.character(one$reason %||% NA_character_),
        stringsAsFactors = FALSE
      )
    }
  }
  skipped_df <- if (length(skipped)) do.call(rbind, skipped)
                else .empty_df(c("population", "population_type", "status", "n_cells", "reason"))
  rownames(skipped_df) <- NULL
  attr(ranks, "skipped") <- skipped_df
  ranks
}

filter_gene_set_list <- function(x, min_size = 5L, universe = NULL) {
  stopifnot(is.list(x), min_size >= 1L)
  out <- lapply(x, function(v) {
    v <- sort(unique(as.character(v[!is.na(v) & v != ""])), method = "radix")
    if (!is.null(universe)) v <- intersect(v, universe)
    v
  })
  out <- out[lengths(out) >= min_size]
  out[sort(names(out), method = "radix")]
}

nfkb_family_sources <- function(collectri,
                                family = c("NFKB", "Nfkb1", "Nfkb2", "Rel", "Rela", "Relb")) {
  stopifnot(is.data.frame(collectri), "source" %in% names(collectri))
  src <- sort(unique(as.character(collectri$source)), method = "radix")
  norm <- toupper(gsub("[^[:alnum:]]", "", src))
  fam <- toupper(gsub("[^[:alnum:]]", "", family))
  src[norm %in% fam]
}

build_project_gene_sets <- function(collectri,
                                    markers = canonical_microglia_markers,
                                    custom_min_size = 5L) {
  stopifnot(is.data.frame(collectri), all(c("source", "target") %in% names(collectri)))
  nfkb_src <- nfkb_family_sources(collectri)
  if (!length(nfkb_src)) stop("no NF-kB family sources found in CollecTRI prior", call. = FALSE)
  raw <- list(
    DAM = markers$DAM,
    Homeostatic = markers$Homeostatic,
    MHC_APC = markers$MHC_APC,
    IFN = markers$IFN,
    NFkB_CollecTRI_Targets = collectri$target[collectri$source %in% nfkb_src],
    NFkB_Activated_Targets = collectri$target[collectri$source %in% nfkb_src & collectri$mor > 0],
    NFkB_Repressed_Targets = collectri$target[collectri$source %in% nfkb_src & collectri$mor < 0]
  )
  sets <- filter_gene_set_list(raw, min_size = custom_min_size)
  stopifnot(all(c("DAM", "Homeostatic", "MHC_APC", "IFN", "NFkB_CollecTRI_Targets",
                  "NFkB_Activated_Targets", "NFkB_Repressed_Targets") %in% names(sets)))
  attr(sets, "nfkb_sources") <- nfkb_src
  sets
}

.msigdbr_go_sets <- function(subcollection, min_size = 15L) {
  if (!requireNamespace("msigdbr", quietly = TRUE)) stop("msigdbr is not installed", call. = FALSE)
  tbl <- msigdbr::msigdbr(db_species = "MM", species = "Mus musculus",
                          collection = "M5", subcollection = subcollection)
  stopifnot(is.data.frame(tbl), all(c("gs_name", "gene_symbol") %in% names(tbl)))
  raw <- split(tbl$gene_symbol, tbl$gs_name)
  filter_gene_set_list(raw, min_size = min_size)
}

build_mechanism_gene_sets <- function(collectri, go_min_size = 15L, custom_min_size = 5L) {
  project <- build_project_gene_sets(collectri, custom_min_size = custom_min_size)
  go_defs <- c(GO_BP = "GO:BP", GO_CC = "GO:CC", GO_MF = "GO:MF")
  go <- lapply(go_defs, .msigdbr_go_sets, min_size = go_min_size)
  sets <- c(go, list(project = project))
  sizes <- do.call(rbind, lapply(names(sets), function(collection) {
    data.frame(collection = collection, set = names(sets[[collection]]),
               size = lengths(sets[[collection]]), row.names = NULL,
               stringsAsFactors = FALSE)
  }))
  rownames(sizes) <- NULL
  out <- list(
    sets = sets,
    sizes = sizes,
    nfkb_sources = attr(project, "nfkb_sources"),
    thresholds = list(go_min_size = go_min_size, custom_min_size = custom_min_size),
    provenance = list(msigdbr_db_species = "MM", species = "Mus musculus",
                      collections = unname(go_defs), n_sets = nrow(sizes))
  )
  fp <- mechanism_gene_set_fingerprint(out, query = out$provenance)
  out$provenance$gene_set_hash <- fp$hash
  out$provenance$n_gene_set_rows <- fp$n_rows
  out$provenance$package_versions <- fp$package_versions
  assert_mechanism_gene_set_expectations(out)
  out
}

mechanism_gene_set_fingerprint <- function(mechanism_gene_sets, query) {
  if (!requireNamespace("digest", quietly = TRUE)) stop("digest is not installed", call. = FALSE)
  stopifnot(is.list(mechanism_gene_sets), is.list(mechanism_gene_sets$sets))
  rows <- do.call(rbind, lapply(names(mechanism_gene_sets$sets), function(collection) {
    do.call(rbind, lapply(names(mechanism_gene_sets$sets[[collection]]), function(set) {
      data.frame(collection = collection, set = set,
                 gene = sort(unique(mechanism_gene_sets$sets[[collection]][[set]]), method = "radix"),
                 stringsAsFactors = FALSE)
    }))
  }))
  rows <- rows[order(rows$collection, rows$set, rows$gene, method = "radix"), , drop = FALSE]
  rownames(rows) <- NULL
  pkg <- .mechanism_pkg_versions(c("msigdbr", "fgsea", "digest"))
  payload <- list(data = rows, query = .sort_named_list(query), package_versions = .sort_named_list(pkg))
  list(hash = digest::digest(payload, algo = "sha256"),
       n_rows = nrow(rows),
       n_sets = nrow(mechanism_gene_sets$sizes),
       n_cols = ncol(rows),
       package_versions = pkg)
}

mechanism_gene_set_expectations <- function() {
  list(hash = "a9c0842dd34a70b4f88b502bf741526d864d8766fe531140345b9f5089c99a2f",
       n_rows = 840988L,
       n_sets = 6142L)
}

assert_mechanism_gene_set_expectations <- function(mechanism_gene_sets,
                                                   expected = mechanism_gene_set_expectations()) {
  prov <- mechanism_gene_sets$provenance
  vals <- list(hash = prov$gene_set_hash,
               n_rows = prov$n_gene_set_rows,
               n_sets = prov$n_sets)
  for (nm in intersect(names(expected), names(vals))) {
    if (!is.na(expected[[nm]]) && !identical(vals[[nm]], expected[[nm]])) {
      stop("MSigDB/mechanism gene-set drift: ", nm, " observed=", vals[[nm]],
           " expected=", expected[[nm]], call. = FALSE)
    }
  }
  TRUE
}

run_mechanism_tf <- function(pb_de_microglia, pb_de_substate, symbol_map, collectri,
                             minsize = 5L) {
  ranks <- collect_rna_rank_matrices(pb_de_microglia, pb_de_substate, symbol_map)
  activity <- lapply(ranks, function(pop) {
    res <- run_decoupler_matrix(pop$matrix, collectri, minsize = minsize)
    data.frame(
      population = pop$population,
      population_type = pop$population_type,
      statistic = res$statistic,
      source = res$source,
      contrast = res$condition,
      score = as.numeric(res$score),
      p_value = as.numeric(res$p_value),
      direction = as.numeric(res$score),
      method = attr(res, "method") %||% "ulm",
      has_consensus = isTRUE(attr(res, "has_consensus")),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, activity)
  out$fdr <- ave(out$p_value, out$population, out$contrast,
                 FUN = function(p) stats::p.adjust(p, method = "BH"))
  ord <- order(out$population, out$contrast, out$fdr, -abs(out$score), out$source,
               method = "radix", na.last = TRUE)
  out <- out[ord, , drop = FALSE]
  rownames(out) <- NULL
  list(
    activity = out,
    skipped = attr(ranks, "skipped"),
    provenance = list(minsize = minsize, direction = "decoupleR ULM score",
                      fdr_scope = "BH within population x contrast")
  )
}

.capture_mechanism_warnings <- function(expr) {
  warnings <- character()
  value <- withCallingHandlers(
    expr,
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = unique(warnings))
}

.fgsea_one <- function(ranks, pathways, collection, min_size, max_size = Inf) {
  stopifnot(is.numeric(ranks), !is.null(names(ranks)), is.list(pathways),
            min_size >= 1L, is.numeric(max_size), length(max_size) == 1L)
  stats <- ranks[is.finite(ranks) & !is.na(names(ranks)) & names(ranks) != ""]
  stats <- stats[!duplicated(names(stats))]
  stats <- sort(stats, decreasing = TRUE)
  max_eff <- min(if (is.finite(max_size)) max_size else Inf, max(1L, length(stats) - 1L))
  pathways <- lapply(pathways, function(x) intersect(x, names(stats)))
  pathways <- pathways[lengths(pathways) >= min_size & lengths(pathways) <= max_eff]
  if (!length(pathways)) {
    return(list(result = .empty_df(c("pathway", "pval", "padj", "log2err", "ES", "NES", "size")),
                warnings = character()))
  }
  bp <- if (requireNamespace("BiocParallel", quietly = TRUE)) {
    BiocParallel::SerialParam(progressbar = FALSE)
  } else {
    NULL
  }
  cap <- .capture_mechanism_warnings(
    fgsea::fgseaMultilevel(pathways = pathways, stats = stats, minSize = min_size,
                           maxSize = max_eff, eps = 1e-10,
                           scoreType = "std", nproc = 0L, BPPARAM = bp)
  )
  known <- "For some pathways, in reality P-values are less than 1e-10."
  unexpected <- cap$warnings[!startsWith(cap$warnings, known)]
  if (length(unexpected)) {
    stop("unexpected fgsea warning: ", paste(unexpected, collapse = " | "), call. = FALSE)
  }
  res <- as.data.frame(cap$value, stringsAsFactors = FALSE)
  keep <- c("pathway", "pval", "padj", "log2err", "ES", "NES", "size")
  stopifnot(all(keep %in% names(res)))
  res <- res[keep]
  res$pathway <- as.character(res$pathway)
  res$collection <- collection
  res$p_floor_warning <- length(cap$warnings) > 0L
  list(result = res, warnings = cap$warnings)
}

run_mechanism_pathway <- function(pb_de_microglia, pb_de_substate, symbol_map,
                                  mechanism_gene_sets, go_max_size = 500L,
                                  project_max_size = Inf) {
  stopifnot(is.list(mechanism_gene_sets), is.list(mechanism_gene_sets$sets))
  ranks <- collect_rna_rank_matrices(pb_de_microglia, pb_de_substate, symbol_map)
  all_rows <- list()
  warning_rows <- list()
  k <- 0L
  w <- 0L
  for (pop_name in names(ranks)) {
    pop <- ranks[[pop_name]]
    for (collection in names(mechanism_gene_sets$sets)) {
      min_size <- if (identical(collection, "project")) {
        mechanism_gene_sets$thresholds$custom_min_size
      } else {
        mechanism_gene_sets$thresholds$go_min_size
      }
      max_size <- if (identical(collection, "project")) project_max_size else go_max_size
      for (contrast in intersect(mechanism_contrasts(), colnames(pop$matrix))) {
        fg <- .fgsea_one(pop$matrix[, contrast], mechanism_gene_sets$sets[[collection]],
                         collection = collection, min_size = min_size, max_size = max_size)
        if (nrow(fg$result)) {
          k <- k + 1L
          all_rows[[k]] <- transform(
            fg$result,
            population = pop$population,
            population_type = pop$population_type,
            contrast = contrast,
            direction = NES,
            p_value = pval,
            fdr = padj
          )
        }
        if (length(fg$warnings)) {
          w <- w + 1L
          warning_rows[[w]] <- data.frame(
            population = pop$population,
            collection = collection,
            contrast = contrast,
            warning = paste(fg$warnings, collapse = " | "),
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  out <- if (length(all_rows)) do.call(rbind, all_rows)
         else .empty_df(c("pathway", "pval", "padj", "log2err", "ES", "NES", "size", "p_floor_warning",
                          "collection", "population", "population_type", "contrast",
                          "direction", "p_value", "fdr"))
  ord <- order(out$population, out$collection, out$contrast, out$fdr, -abs(out$NES), out$pathway,
               method = "radix", na.last = TRUE)
  out <- out[ord, , drop = FALSE]
  rownames(out) <- NULL
  warn <- if (length(warning_rows)) do.call(rbind, warning_rows)
          else .empty_df(c("population", "collection", "contrast", "warning"))
  list(
    pathway = out,
    skipped = attr(ranks, "skipped"),
    warnings = warn,
    provenance = list(direction = "fgsea NES",
                      fdr_scope = "fgsea padj within population x collection x contrast",
                      go_max_size = go_max_size,
                      project_max_size = project_max_size)
  )
}

.family_best_bh <- function(score, p_value, label) {
  ok <- is.finite(score) & is.finite(p_value) & p_value >= 0 & p_value <= 1 & !is.na(label) & label != ""
  if (!any(ok)) {
    return(list(score = NA_real_, p_value = NA_real_, n = 0L,
                detail = NA_character_, raw_p_value = NA_real_))
  }
  adj <- stats::p.adjust(p_value[ok], method = "BH")
  best <- which.min(adj)
  list(score = score[ok][best],
       p_value = adj[best],
       n = sum(ok),
       detail = as.character(label[ok][best]),
       raw_p_value = p_value[ok][best])
}

nfkb_tf_family_table <- function(tf_activity, nfkb_sources,
                                 contrasts = c("interaction", "tau_in_nlgf")) {
  stopifnot(is.data.frame(tf_activity), length(nfkb_sources) >= 1L,
            all(c("population", "population_type", "source", "contrast", "score", "p_value") %in%
                  names(tf_activity)))
  sub <- tf_activity[tf_activity$source %in% nfkb_sources & tf_activity$contrast %in% contrasts, , drop = FALSE]
  if (!nrow(sub)) {
    return(.empty_df(c("population", "population_type", "contrast", "test", "score", "p_value",
                       "direction", "n_sources", "detail")))
  }
  key <- interaction(sub$population, sub$contrast, drop = TRUE, lex.order = TRUE)
  rows <- lapply(split(sub, key), function(d) {
    agg <- .family_best_bh(d$score, d$p_value, d$source)
    data.frame(
      population = d$population[1],
      population_type = d$population_type[1],
      contrast = d$contrast[1],
      test = "tf_family",
      score = agg$score,
      p_value = agg$p_value,
      raw_p_value = agg$raw_p_value,
      direction = agg$score,
      n_sources = agg$n,
      detail = agg$detail,
      family_members = paste(sort(unique(d$source), method = "radix"), collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

nfkb_target_gsea_table <- function(mechanism_pathway,
                                   components = c(NFkB_Activated_Targets = 1,
                                                  NFkB_Repressed_Targets = -1),
                                   contrasts = c("interaction", "tau_in_nlgf")) {
  pathway <- if (is.list(mechanism_pathway)) mechanism_pathway$pathway else mechanism_pathway
  stopifnot(is.data.frame(pathway),
            all(c("population", "population_type", "collection", "contrast", "pathway",
                  "NES", "pval", "padj", "size", "p_floor_warning") %in% names(pathway)))
  sub <- pathway[pathway$collection == "project" & pathway$pathway %in% names(components) &
                   pathway$contrast %in% contrasts, , drop = FALSE]
  if (!nrow(sub)) {
    return(.empty_df(c("population", "population_type", "contrast", "test", "score", "p_value",
                       "raw_p_value", "direction", "n_sources", "detail", "family_members",
                       "pathway_fdr", "size", "p_floor_warning")))
  }
  sub$signed_score <- sub$NES * unname(components[sub$pathway])
  key <- interaction(sub$population, sub$contrast, drop = TRUE, lex.order = TRUE)
  rows <- lapply(split(sub, key), function(d) {
    agg <- .family_best_bh(d$signed_score, d$pval, d$pathway)
    data.frame(
      population = d$population[1],
      population_type = d$population_type[1],
      contrast = d$contrast[1],
      test = "target_gsea",
      score = agg$score,
      p_value = agg$p_value,
      raw_p_value = agg$raw_p_value,
      direction = agg$score,
      n_sources = NA_integer_,
      detail = agg$detail,
      family_members = paste(sort(unique(d$pathway), method = "radix"), collapse = ";"),
      pathway_fdr = min(stats::p.adjust(d$pval, method = "BH"), na.rm = TRUE),
      size = d$size[match(agg$detail, d$pathway)],
      p_floor_warning = any(d$p_floor_warning),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

call_nfkb_attenuation <- function(rows, alpha = 0.10) {
  stopifnot(is.data.frame(rows), all(c("population", "contrast", "test", "score", "p_value") %in% names(rows)))
  primary_idx <- rows$population == "whole_microglia" & rows$contrast == "interaction" &
    rows$test %in% c("tf_family", "target_gsea")
  primary <- rows[primary_idx, , drop = FALSE]
  stopifnot(nrow(primary) == 2L, all(is.finite(primary$p_value)))
  primary$primary_family_fdr <- stats::p.adjust(primary$p_value, method = "BH")
  rows$primary_family_fdr <- NA_real_
  rows$primary_family_fdr[primary_idx] <- primary$primary_family_fdr
  primary_negative <- all(primary$score < 0)
  primary_positive <- all(primary$score > 0)
  primary_discordant <- any(primary$score < 0) && any(primary$score > 0)
  supported <- primary_negative && any(primary$primary_family_fdr < alpha)
  list(
    rows = rows,
    verdict = list(
      alpha = alpha,
      status = if (supported) "supported" else if (primary_discordant) "discordant"
               else if (primary_negative) "directional_only" else "not_supported",
      supported = supported,
      primary_negative = primary_negative,
      primary_positive = primary_positive,
      primary_discordant = primary_discordant,
      n_primary_supported = sum(primary$score < 0 & primary$primary_family_fdr < alpha),
      rule = "Supported attenuation requires concordant negative whole_microglia interaction primary rows; tau_in_nlgf and substates are supportive."
    )
  )
}

build_nfkb_attenuation <- function(mechanism_tf, mechanism_pathway, mechanism_gene_sets,
                                   alpha = 0.10) {
  stopifnot(is.list(mechanism_tf), is.data.frame(mechanism_tf$activity),
            is.list(mechanism_pathway), is.data.frame(mechanism_pathway$pathway),
            is.list(mechanism_gene_sets), length(mechanism_gene_sets$nfkb_sources) >= 1L)
  tf <- nfkb_tf_family_table(mechanism_tf$activity, mechanism_gene_sets$nfkb_sources)
  gs <- nfkb_target_gsea_table(mechanism_pathway)
  common <- union(names(tf), names(gs))
  fill <- function(d) {
    miss <- setdiff(common, names(d))
    for (m in miss) d[[m]] <- NA
    d[common]
  }
  rows <- rbind(fill(tf), fill(gs))
  if (!"p_floor_warning" %in% names(rows)) rows$p_floor_warning <- FALSE
  rows$p_floor_warning[is.na(rows$p_floor_warning)] <- FALSE
  rows$primary_test <- rows$population == "whole_microglia" & rows$contrast == "interaction" &
    rows$test %in% c("tf_family", "target_gsea")
  rows$supportive_only <- !rows$primary_test
  ord <- order(rows$primary_test, rows$population, rows$contrast, rows$test,
               decreasing = c(TRUE, FALSE, FALSE, FALSE), method = "radix")
  rows <- rows[ord, , drop = FALSE]
  rownames(rows) <- NULL
  called <- call_nfkb_attenuation(rows, alpha = alpha)
  list(
    table = called$rows,
    verdict = called$verdict,
    skipped = mechanism_tf$skipped,
    provenance = list(nfkb_sources = mechanism_gene_sets$nfkb_sources,
                      alpha = alpha,
                      primary_family = c("tf_family", "target_gsea"))
  )
}

# ---- P3-S3: minimal 24M phosphosite DE + kinase activity ------------------------------

phosphosite_probability_column <- function(phospho_tbl) {
  candidates <- c("Phosphosite probability", "Phopshosite probability")
  hit <- intersect(candidates, names(phospho_tbl))
  if (!length(hit)) {
    stop("missing phosphosite probability column; expected one of: ",
         paste(candidates, collapse = ", "), call. = FALSE)
  }
  hit[1]
}

match_24m_intensity_columns <- function(phospho_tbl, sample_key, n_keep = 16L) {
  stopifnot(is.data.frame(phospho_tbl), is.data.frame(sample_key),
            all(c("genotype", "col_stub", "label") %in% names(sample_key)))
  if (nrow(sample_key) != n_keep) {
    stop("sample_key must contain exactly ", n_keep, " 24M rows", call. = FALSE)
  }
  if (anyDuplicated(sample_key$col_stub)) stop("duplicate sample-key stubs", call. = FALSE)
  hits <- match_intensity_columns(names(phospho_tbl), sample_key)
  hits <- hits[!is.na(hits$key_idx), , drop = FALSE]
  if (nrow(hits) != n_keep || anyDuplicated(hits$key_idx) ||
      !setequal(hits$key_idx, seq_len(n_keep))) {
    stop("expected exactly ", n_keep, "/", n_keep, " matched 24M intensity columns", call. = FALSE)
  }
  if (anyDuplicated(hits$stub)) stop("duplicate matched phospho intensity stubs", call. = FALSE)
  hits <- hits[order(hits$key_idx), , drop = FALSE]
  geno <- factor(as.character(hits$genotype), levels = genotype_levels)
  if (anyNA(geno) || !all(as.integer(table(geno)) == n_keep / 4L)) {
    stop("24M phospho columns are not balanced 4/genotype", call. = FALSE)
  }
  meta <- data.frame(
    sample_id = hits$stub,
    column = hits$column,
    stub = hits$stub,
    label = hits$label,
    genotype = geno,
    run_index = seq_len(nrow(hits)),
    stringsAsFactors = FALSE
  )
  rownames(meta) <- meta$sample_id
  list(columns = hits$column, meta = meta,
       matched = hits, n_expected = n_keep, n_matched = nrow(hits))
}

phospho_feature_frame <- function(phospho_tbl) {
  stopifnot(is.data.frame(phospho_tbl))
  need <- c("PG.Genes", "PTM.SiteAA", "PTM.SiteLocation", "PTM.CollapseKey")
  stopifnot(all(need %in% names(phospho_tbl)))
  prob_col <- phosphosite_probability_column(phospho_tbl)
  row_index <- seq_len(nrow(phospho_tbl))
  collapse_key <- trimws(as.character(phospho_tbl[["PTM.CollapseKey"]]))
  collapse_blank <- is.na(collapse_key) | collapse_key == ""
  collapse_trace <- ifelse(collapse_blank, "<blank>", collapse_key)
  feature <- paste0("row", row_index, "|", collapse_trace)
  stopifnot(anyDuplicated(feature) == 0L)

  gene <- trimws(as.character(phospho_tbl[["PG.Genes"]]))
  aa <- trimws(as.character(phospho_tbl[["PTM.SiteAA"]]))
  loc <- trimws(as.character(phospho_tbl[["PTM.SiteLocation"]]))
  missing_gene <- is.na(gene) | gene == ""
  multi_gene <- grepl("[;,]", gene)
  missing_site <- is.na(aa) | aa == "" | is.na(loc) | loc == "" | tolower(loc) %in% c("na", "nan")
  site_id <- ifelse(!missing_gene & !multi_gene & !missing_site, paste0(gene, "_", aa, loc), NA_character_)
  probability <- suppressWarnings(as.numeric(phospho_tbl[[prob_col]]))
  out <- data.frame(
    feature = feature,
    original_row = row_index,
    collapse_key = collapse_key,
    collapse_key_blank = collapse_blank,
    gene = gene,
    site_aa = aa,
    site_location = loc,
    site_id = site_id,
    phosphosite_probability = probability,
    stringsAsFactors = FALSE
  )
  attr(out, "phospho_feature_counts") <- list(
    probability_col = prob_col,
    n_rows = nrow(out),
    n_blank_collapse_key = sum(collapse_blank),
    n_duplicate_collapse_key = sum(duplicated(collapse_key[!collapse_blank])),
    n_missing_gene = sum(missing_gene),
    n_multi_gene = sum(multi_gene & !missing_gene),
    n_missing_site = sum(missing_site),
    n_site_rows = sum(!is.na(site_id)),
    n_unique_sites = length(unique(site_id[!is.na(site_id)]))
  )
  out
}

positive_log2_matrix <- function(mat) {
  stopifnot(is.matrix(mat), is.numeric(mat))
  nonpositive <- !is.na(mat) & mat <= 0
  out <- mat
  out[nonpositive] <- NA_real_
  out <- log2(out)
  attr(out, "log2_counts") <- list(
    n_values = length(mat),
    n_missing_input = sum(is.na(mat)),
    n_nonpositive_to_na = sum(nonpositive),
    n_missing_output = sum(is.na(out))
  )
  out
}

prepare_phospho_24m_matrix <- function(phospho_tbl, sample_key,
                                       min_present = 2L, min_groups = 4L) {
  match <- match_24m_intensity_columns(phospho_tbl, sample_key)
  feat <- phospho_feature_frame(phospho_tbl)
  feature_counts <- attr(feat, "phospho_feature_counts")
  raw <- as.matrix(phospho_tbl[, match$columns, drop = FALSE])
  storage.mode(raw) <- "double"
  rownames(raw) <- feat$feature
  colnames(raw) <- rownames(match$meta)
  log_mat <- positive_log2_matrix(raw)
  norm <- median_normalise(log_mat)
  filt <- prevalence_filter(norm, match$meta$genotype,
                            min_present = min_present, min_groups = min_groups)
  feat <- feat[match(rownames(filt), feat$feature), , drop = FALSE]
  stopifnot(identical(feat$feature, rownames(filt)),
            identical(colnames(filt), rownames(match$meta)))
  list(
    matrix = filt,
    meta = match$meta,
    features = feat,
    matched = match$matched,
    counts = c(feature_counts, attr(log_mat, "log2_counts"),
               list(n_features_raw = nrow(raw),
                    n_features_filtered = nrow(filt),
                    n_filtered_site_rows = sum(!is.na(feat$site_id) & feat$site_id != ""),
                    n_filtered_unique_sites = length(unique(feat$site_id[!is.na(feat$site_id) &
                                                                          feat$site_id != ""])),
                    min_present = min_present,
                    min_groups = min_groups))
  )
}

annotate_phospho_top_tables <- function(top, features) {
  stopifnot(is.list(top), is.data.frame(features), "feature" %in% names(features))
  lapply(top, function(tbl) {
    stopifnot(is.data.frame(tbl), "feature" %in% names(tbl))
    idx <- match(tbl$feature, features$feature)
    if (anyNA(idx)) stop("top table contains feature absent from annotation", call. = FALSE)
    cbind(tbl, features[idx, setdiff(names(features), "feature"), drop = FALSE])
  })
}

run_index_factorial_design <- function(meta, run_col = "run_index") {
  stopifnot(is.data.frame(meta), run_col %in% names(meta))
  fd <- factorial_design(meta, add_batch = FALSE)
  run <- as.numeric(meta[[run_col]])
  if (!all(is.finite(run)) || length(unique(run)) < 2L) {
    return(list(status = "skipped", reason = "run_index is not finite and variable"))
  }
  run_z <- as.numeric(scale(run))
  design <- cbind(fd$design, run_index = run_z)
  if (qr(design)$rank != ncol(design)) {
    return(list(status = "skipped", reason = "additive run-index design is rank deficient"))
  }
  contrasts <- rbind(fd$contrasts, run_index = rep(0, ncol(fd$contrasts)))
  stopifnot(identical(rownames(contrasts), colnames(design)),
            "interaction" %in% colnames(contrasts))
  list(status = "fit", design = design, contrasts = contrasts,
       run_col = run_col, scaled_center = mean(run), scaled_scale = stats::sd(run))
}

run_phospho_de_24m <- function(phospho_tbl, sample_key,
                               min_present = 2L, min_groups = 4L) {
  prep <- prepare_phospho_24m_matrix(phospho_tbl, sample_key,
                                     min_present = min_present, min_groups = min_groups)
  fd <- factorial_design(prep$meta, add_batch = FALSE)
  fit <- fit_limma_log(prep$matrix, fd$design, fd$contrasts)
  top <- annotate_phospho_top_tables(fit$top, prep$features)

  run_index <- run_index_factorial_design(prep$meta)
  if (identical(run_index$status, "fit")) {
    fit_ri <- fit_limma_log(prep$matrix, run_index$design, run_index$contrasts)
    run_index$top <- annotate_phospho_top_tables(fit_ri$top, prep$features)
  }
  list(
    n_samples = ncol(prep$matrix),
    n_features = nrow(prep$matrix),
    meta = prep$meta,
    features = prep$features,
    top = top,
    design = list(add_batch = FALSE, design_cols = colnames(fd$design),
                  contrast_names = colnames(fd$contrasts),
                  residual_df = nrow(fd$design) - qr(fd$design)$rank),
    run_index = run_index,
    filters = prep$counts,
    provenance = list(
      transform = "log2 positive intensities; nonpositive values set to NA",
      normalisation = "sample-wise median shift on log2 scale",
      prevalence = list(min_present = min_present, min_groups = min_groups),
      feature_id = "row<original_row>|<PTM.CollapseKey>"
    )
  )
}

phospho_site_stat_matrix <- function(top_list, stat_col = "t") {
  stopifnot(is.list(top_list), length(top_list) >= 1L)
  present <- intersect(mechanism_contrasts(), names(top_list))
  if (!length(present)) stop("no canonical contrasts present in phospho top tables", call. = FALSE)
  collapsed <- list()
  count_rows <- list()
  for (cn in present) {
    top <- top_list[[cn]]
    need <- c("feature", "site_id", "phosphosite_probability", "original_row", stat_col)
    stopifnot(is.data.frame(top), all(need %in% names(top)))
    stat <- as.numeric(top[[stat_col]])
    ok <- !is.na(top$site_id) & top$site_id != "" & is.finite(stat)
    d <- top[ok, , drop = FALSE]
    d$.stat <- stat[ok]
    if (!nrow(d)) stop("no finite single-site phospho statistics for contrast: ", cn, call. = FALSE)
    site_counts <- table(d$site_id)
    by_site <- split(d, d$site_id)
    vals <- vapply(by_site, function(x) {
      prob <- as.numeric(x$phosphosite_probability)
      prob[!is.finite(prob)] <- -Inf
      ord <- order(-prob, -abs(x$.stat), as.integer(x$original_row), method = "radix")
      x$.stat[ord[1]]
    }, numeric(1))
    vals <- vals[sort(names(vals), method = "radix")]
    collapsed[[cn]] <- vals
    count_rows[[cn]] <- data.frame(
      contrast = cn,
      n_rows_input = nrow(top),
      n_rows_usable = nrow(d),
      n_sites = length(site_counts),
      n_duplicate_sites = sum(site_counts > 1L),
      max_rows_per_site = max(as.integer(site_counts)),
      median_rows_per_site = stats::median(as.integer(site_counts)),
      stringsAsFactors = FALSE
    )
  }
  sites <- sort(unique(unlist(lapply(collapsed, names), use.names = FALSE)), method = "radix")
  mat <- matrix(NA_real_, nrow = length(sites), ncol = length(present),
                dimnames = list(sites, present))
  for (cn in present) mat[names(collapsed[[cn]]), cn] <- collapsed[[cn]]
  stopifnot(anyDuplicated(rownames(mat)) == 0L)
  attr(mat, "collapse_counts") <- do.call(rbind, count_rows)
  mat
}

assert_ksn_activity_coverage <- function(coverage, min_kinases = 50L,
                                         min_matched_sites = 100L,
                                         kinase = "Gsk3b",
                                         require_kinase_minsize = TRUE) {
  stopifnot(is.list(coverage), is.list(coverage$gsk3b))
  if (coverage$n_matched_sites < min_matched_sites) {
    stop("KSN coverage too low: matched_sites=", coverage$n_matched_sites,
         " < ", min_matched_sites, call. = FALSE)
  }
  if (coverage$kinases_passing_minsize < min_kinases) {
    stop("KSN coverage too low: kinases_passing_minsize=",
         coverage$kinases_passing_minsize, " < ", min_kinases, call. = FALSE)
  }
  if (require_kinase_minsize && !isTRUE(coverage$gsk3b$passes_minsize)) {
    stop(kinase, " does not pass KSN minsize; matched_sites=",
         coverage$gsk3b$matched_sites, call. = FALSE)
  }
  TRUE
}

activity_from_site_matrix <- function(site_mat, ksn, minsize = 5L, fit = "primary") {
  stopifnot(is.matrix(site_mat), anyDuplicated(rownames(site_mat)) == 0L)
  net <- ksn[ksn$target %in% rownames(site_mat), , drop = FALSE]
  res <- run_decoupler_matrix(site_mat, net, minsize = minsize)
  out <- data.frame(
    fit = fit,
    statistic = res$statistic,
    source = res$source,
    contrast = res$condition,
    score = as.numeric(res$score),
    p_value = as.numeric(res$p_value),
    direction = as.numeric(res$score),
    method = attr(res, "method") %||% "ulm",
    stringsAsFactors = FALSE
  )
  out$fdr <- ave(out$p_value, out$fit, out$contrast,
                 FUN = function(p) stats::p.adjust(p, method = "BH"))
  ord <- order(out$fit, out$contrast, out$fdr, -abs(out$score), out$source,
               method = "radix", na.last = TRUE)
  out[ord, , drop = FALSE]
}

run_kinase_activity <- function(phospho_de_24m, minsize = 5L,
                                kinase = "Gsk3b",
                                coverage_min_kinases = 50L,
                                coverage_min_sites = 100L) {
  stopifnot(is.list(phospho_de_24m), is.list(phospho_de_24m$top))
  ksn <- load_omnipath_ksn_mouse()
  assert_mechanism_prior_expectations(ksn = ksn)
  primary_mat <- phospho_site_stat_matrix(phospho_de_24m$top)
  coverage <- ksn_coverage_probe(ksn, rownames(primary_mat), minsize = minsize, kinase = kinase)
  assert_ksn_activity_coverage(coverage, min_kinases = coverage_min_kinases,
                               min_matched_sites = coverage_min_sites,
                               kinase = kinase, require_kinase_minsize = TRUE)
  primary <- activity_from_site_matrix(primary_mat, ksn, minsize = minsize, fit = "primary")
  rows <- list(primary)
  sensitivity <- list(status = phospho_de_24m$run_index$status %||% "absent",
                      reason = phospho_de_24m$run_index$reason %||% NA_character_)
  if (identical(phospho_de_24m$run_index$status, "fit") &&
      is.list(phospho_de_24m$run_index$top)) {
    run_mat <- phospho_site_stat_matrix(phospho_de_24m$run_index$top)
    rows[[length(rows) + 1L]] <- activity_from_site_matrix(run_mat, ksn, minsize = minsize,
                                                           fit = "run_index")
    sensitivity$collapse_counts <- attr(run_mat, "collapse_counts")
  }
  activity <- do.call(rbind, rows)
  rownames(activity) <- NULL
  list(
    activity = activity,
    coverage = coverage,
    sensitivity = sensitivity,
    provenance = list(
      minsize = minsize,
      fdr_scope = "BH within fit x contrast",
      kinase_prior = attr(ksn, "provenance"),
      primary_collapse_counts = attr(primary_mat, "collapse_counts"),
      bulk_context = "24M bulk hippocampus phosphoproteomics; not microglia-sorted"
    )
  )
}

build_kinase_mechanism_summary <- function(kinase_activity, kinase = "Gsk3b",
                                           alpha = 0.10) {
  activity <- if (is.list(kinase_activity)) kinase_activity$activity else kinase_activity
  stopifnot(is.data.frame(activity),
            all(c("fit", "source", "contrast", "score", "p_value", "fdr") %in% names(activity)))
  primary <- activity[activity$fit == "primary", , drop = FALSE]
  primary$significant <- is.finite(primary$fdr) & primary$fdr < alpha
  primary$include_reason <- ifelse(primary$source == kinase & primary$significant,
                                   "significant+gsk3b_carry",
                                   ifelse(primary$source == kinase, "gsk3b_carry", "significant"))
  keep <- primary$significant | primary$source == kinase
  out <- primary[keep, , drop = FALSE]
  missing_contrasts <- setdiff(mechanism_contrasts(), out$contrast[out$source == kinase])
  if (length(missing_contrasts)) {
    out <- rbind(out, data.frame(
      fit = "primary", statistic = NA_character_, source = kinase,
      contrast = missing_contrasts, score = NA_real_, p_value = NA_real_,
      direction = NA_real_, method = NA_character_, fdr = NA_real_,
      significant = FALSE, include_reason = "gsk3b_absent",
      stringsAsFactors = FALSE
    ))
  }
  sens <- activity[activity$fit == "run_index", c("source", "contrast", "score", "fdr"), drop = FALSE]
  names(sens) <- c("source", "contrast", "run_index_score", "run_index_fdr")
  out <- merge(out, sens, by = c("source", "contrast"), all.x = TRUE, sort = FALSE)
  out$run_index_sign_concordant <- is.finite(out$score) & is.finite(out$run_index_score) &
    sign(out$score) == sign(out$run_index_score)
  out$run_index_supports <- is.finite(out$run_index_fdr) & out$run_index_fdr < alpha &
    out$run_index_sign_concordant
  out$run_order_confounded <- out$source == kinase & out$significant & !out$run_index_supports
  out$alpha <- alpha
  ord_contrast <- match(out$contrast, mechanism_contrasts())
  ord <- order(ord_contrast, out$source != kinase, out$fdr, -abs(out$score),
               method = "radix", na.last = TRUE)
  out <- out[ord, , drop = FALSE]
  rownames(out) <- NULL
  list(
    table = out,
    coverage = if (is.list(kinase_activity)) kinase_activity$coverage else NULL,
    provenance = list(alpha = alpha,
                      carry_kinase = kinase,
                      rule = "primary significant kinases plus explicit Gsk3b row per contrast")
  )
}
