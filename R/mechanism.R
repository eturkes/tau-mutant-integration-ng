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

extract_rank_matrix <- function(top_list, symbol_map, stat_col = "t") {
  stopifnot(is.list(top_list), length(top_list) >= 1L, !is.null(names(top_list)), nzchar(stat_col))
  present <- names(top_list)
  if (any(!nzchar(present)) || anyDuplicated(present)) stop("top_list must have unique nonblank names", call. = FALSE)
  if (exists("contrast_definitions", inherits = TRUE)) {
    present <- intersect(c("tau_alone", "nlgf_in_maptki", "nlgf_in_p301s", "tau_in_nlgf", "interaction"), present)
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
