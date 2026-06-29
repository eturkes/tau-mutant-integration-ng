# fgsea / pathway-analysis helpers: cached MSigDB loaders, single-contrast
# fgsea wrapper, per-state extension, gene-level cross-modality joins, and
# the focused custom-pathway fgsea used to score snRNAseq-derived gene
# sets against bulk modalities.

get_gobp <- function(cache_path = "storage/cache/msigdb_gobp_mouse.rds") {
  cache_or_run(cache_path, {
    df <- msigdbr::msigdbr(species = "Mus musculus", collection = "C5",
                           subcollection = "GO:BP")
    split(df$gene_symbol, df$gs_name)
  })
}

get_gomf <- function(cache_path = "storage/cache/msigdb_gomf_mouse.rds") {
  cache_or_run(cache_path, {
    df <- msigdbr::msigdbr(species = "Mus musculus", collection = "C5",
                           subcollection = "GO:MF")
    split(df$gene_symbol, df$gs_name)
  })
}

get_gocc <- function(cache_path = "storage/cache/msigdb_gocc_mouse.rds") {
  cache_or_run(cache_path, {
    df <- msigdbr::msigdbr(species = "Mus musculus", collection = "C5",
                           subcollection = "GO:CC")
    split(df$gene_symbol, df$gs_name)
  })
}

# Curated microglia-state gene-set collection (7 mouse-symbol sets) built
# offline by `scripts/build_custom_microglia_states.R` from the original
# Keren-Shaul / Marsh / Friedman / Olah / Sala Frigerio supplements. The
# build script is the single source of truth for set contents and
# provenance (see `storage/cache/custom_microglia_states_provenance.txt`).
# This getter only LOADS the cache so the Rmd pipeline stays a pure
# consumer; if the cache is absent the user is told exactly which script
# to run rather than the helper silently rebuilding.
get_custom_microglia_states <- function(
    cache_path = "storage/cache/custom_microglia_states.rds") {
  if (!file.exists(cache_path)) {
    stop(sprintf(
      "%s not found. Build it once with:\n  Rscript scripts/build_custom_microglia_states.R\n(see storage/cache/custom_microglia_states_provenance.txt for source-paper details).",
      cache_path), call. = FALSE)
  }
  readRDS(cache_path)
}

# Curated tau / amyloid-specific microglia gene-set collection (4 mouse-
# symbol sets: AD1, AD2, LDAM, WAM) built offline by
# `scripts/build_custom_microglia_ad.R` from the Safaiyan 2021 /
# Marschallinger 2020 / Gerrits 2022 supplements. The build script is
# the single source of truth for set contents and provenance (see
# `storage/cache/custom_microglia_ad_provenance.txt`). Same pure-reader
# contract as `get_custom_microglia_states`: the Rmd pipeline is a
# consumer; absent cache => instruct, never silently rebuild.
get_custom_microglia_ad <- function(
    cache_path = "storage/cache/custom_microglia_ad.rds") {
  if (!file.exists(cache_path)) {
    stop(sprintf(
      "%s not found. Build it once with:\n  Rscript scripts/build_custom_microglia_ad.R\n(see storage/cache/custom_microglia_ad_provenance.txt for source-paper details).",
      cache_path), call. = FALSE)
  }
  readRDS(cache_path)
}

# Curated microglia transcriptional-state "module sources" gene-set
# collection (12 mouse-symbol sets representing the published MG0-MG12
# states from Sun, Victor, Mathys et al. 2023 Cell, with MG9 skipped per
# the published 12-state taxonomy) built offline by
# `scripts/build_custom_module_sources.R`. The build script is the single
# source of truth for set contents and provenance (see
# `storage/cache/custom_module_sources_provenance.txt`, which also
# documents the B4 decision-gate deviation explaining why this collection
# is sourced from Sun/Victor 2023 instead of the plan's original Mathys
# 2019 + Olah 2020 defaults). Same pure-reader contract as
# `get_custom_microglia_states` / `get_custom_microglia_ad`: the Rmd
# pipeline is a consumer; absent cache => instruct, never silently
# rebuild.
get_custom_module_sources <- function(
    cache_path = "storage/cache/custom_module_sources.rds") {
  if (!file.exists(cache_path)) {
    stop(sprintf(
      "%s not found. Build it once with:\n  Rscript scripts/build_custom_module_sources.R\n(see storage/cache/custom_module_sources_provenance.txt for source-paper details).",
      cache_path), call. = FALSE)
  }
  readRDS(cache_path)
}

run_fgsea_for_contrast <- function(top_tbl, gene_sets, stat_col = "t",
                                   gene_col = "symbol", min_size = 10,
                                   max_size = 500) {
  ranks <- top_tbl[[stat_col]]
  names(ranks) <- top_tbl[[gene_col]]
  ranks <- ranks[!is.na(names(ranks)) & !is.na(ranks)]
  ranks <- ranks[!duplicated(names(ranks))]
  ranks <- sort(ranks)
  fgsea::fgsea(pathways = gene_sets, stats = ranks,
               minSize = min_size, maxSize = max_size, eps = 0)
}

# Per-modality fgsea driver: deduplicate the gene column, fall back to logFC
# if no `t` column exists, and run `run_fgsea_for_contrast` for each of the
# five canonical contrasts. Used by both `rmd/07_integration.Rmd` and the
# stand-alone build scripts so collection-specific runs share one shim.
prep_t_dedup <- function(top_tbl, gene_col = "symbol") {
  tbl <- top_tbl |> dplyr::filter(!is.na(.data[[gene_col]]),
                                  .data[[gene_col]] != "")
  tbl[!duplicated(tbl[[gene_col]]), ]
}

run_fgsea_per_dataset <- function(de_obj, pathways, gene_col = "symbol",
                                  min_size = 15, max_size = 500,
                                  contrasts = c("nlgf_in_maptki",
                                                "nlgf_in_p301s",
                                                "interaction", "tau_alone",
                                                "tau_in_nlgf")) {
  out <- lapply(contrasts, function(cn) {
    tbl <- prep_t_dedup(de_obj$fit$top[[cn]], gene_col)
    if (!"t" %in% names(tbl)) tbl$t <- tbl$logFC
    run_fgsea_for_contrast(tbl, pathways, stat_col = "t",
                           gene_col = gene_col,
                           min_size = min_size, max_size = max_size)
  })
  names(out) <- contrasts
  out
}

# Join one fgsea per-modality-per-contrast result list into a tidy
# pathway-keyed table: one row per (pathway, contrast) with per-modality
# `nes_<m>` and `padj_<m>` columns. Mirrors the join used in
# `rmd/07_integration.Rmd` for GO BP so script-built collection caches
# yield the same TSV schema.
join_fgsea_results <- function(fgsea_results,
                               contrasts = c("nlgf_in_maptki",
                                             "nlgf_in_p301s",
                                             "interaction", "tau_alone",
                                             "tau_in_nlgf")) {
  modalities <- names(fgsea_results)
  purrr::map_dfr(contrasts, function(cn) {
    extract <- function(name) {
      tbl <- fgsea_results[[name]][[cn]] |> as.data.frame()
      tbl[, c("pathway", "NES", "padj")] |>
        setNames(c("pathway", paste0("nes_", name), paste0("padj_", name)))
    }
    Reduce(function(a, b) dplyr::full_join(a, b, by = "pathway"),
           lapply(modalities, extract)) |>
      dplyr::mutate(contrast = cn)
  })
}

# Run fgsea on every contrast of a per-substate NEBULA fit list (output of
# `fit_nebula_per_state`). Returns a nested list `state -> contrast -> fgseaResult`.
# Uses the same NEBULA z-statistic (`t`) ranking as the whole-microglia fgsea.
run_fgsea_per_state <- function(per_state_fits, pathways,
                                contrasts = c("nlgf_in_maptki", "nlgf_in_p301s",
                                              "interaction", "tau_alone",
                                              "tau_in_nlgf"),
                                gene_col = "symbol",
                                min_size = 10, max_size = 500) {
  out <- list()
  for (st in names(per_state_fits)) {
    fit <- per_state_fits[[st]]
    out[[st]] <- lapply(contrasts, function(cn) {
      run_fgsea_for_contrast(fit$top[[cn]], pathways, stat_col = "t",
                             gene_col = gene_col,
                             min_size = min_size, max_size = max_size)
    })
    names(out[[st]]) <- contrasts
  }
  out
}
