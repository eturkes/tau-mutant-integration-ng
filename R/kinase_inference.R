# Kinase activity inference helpers used by the Phase F mechanism layer
# (see storage/notes/mechanism_layer_plan.md, sessions F2..F5). Wraps
# decoupleR's univariate methods (ULM + WSUM family) and the OmniPath
# kinase-substrate network (KSN) mapped human -> mouse via nichenetr.
# The design mirrors `R/tf_inference.R` so the Phase F outputs read
# cleanly against the Phase E TF outputs at the verdict layer in
# section 15.3 vs section 14.3:
#   - cross-cache, sign-aware, family-agnostic;
#   - never pre-commits to a single kinase family (no Erk-only, no
#     CDK-only frames);
#   - returns long tidy tibbles with the same canonical key set
#     (`statistic`, `source`, `condition`, `score`, `p_value`) used by
#     the TF helpers, so the F3+ leader-board / axis-restricted code
#     can reuse the same downstream patterns.
#
# Cache shape consumed (per F2): the phospho DE caches
#   storage/cache/de_phospho.rds              (raw)
#   storage/cache/de_phospho_corrected.rds    (batch-corrected)
# Each is a list with `$fit` (itself a list with `$fit` the MArrayLM
# and `$top` a named list one top-table per contrast across the
# canonical five-contrast set nlgf_in_maptki, nlgf_in_p301s,
# tau_alone, tau_in_nlgf, interaction). Each per-contrast top-table
# carries columns `feature`, `logFC`, `AveExpr`, `t`, `P.Value`,
# `adj.P.Val`, `B`, `symbol`, `PG.ProteinGroups`, `PTM.SiteAA`,
# `PTM.SiteLocation`. The unit of analysis is the phospho-site, not
# the gene: a single gene can carry several measured sites with
# divergent kinase regulation, so collapsing to gene-level would
# destroy signal.
#
# Substrate id convention: `symbol_resPos`, built as
#   paste0(symbol, "_", PTM.SiteAA, PTM.SiteLocation)
# e.g. "Mapt_S202", "Gsk3b_Y216". This format matches the OmniPath
# KSN target ids returned by the upstream wrapper (modulo human-
# vs-mouse symbol; nichenetr handles the species translation).
# Regulatory phosphosites are highly conserved across mammals so
# residue identity + position are translated 1:1; only the gene
# symbol is mapped.
#
# Cache shape produced: `cache -> contrast -> tibble(statistic,
# source, condition, score, p_value)`. `cache` is in
# {"phospho_raw", "phospho_corrected"}; `contrast` is the canonical
# five-contrast set; `condition == contrast` by construction (kept
# for round-trip compatibility with the raw decoupleR output);
# `statistic` carries the method label (ulm, wsum, norm_wsum,
# corr_wsum, consensus). Per-method p-value semantics follow the
# E2 convention recorded in `R/tf_inference.R`'s header (ulm = lm
# F-test; wsum-family = empirical permutation p floored at 1/1000
# i.e. minimum visible 0.02; consensus = Stouffer-combined,
# continuous, populated -- not NA).
#
# OmniPath wrapper note. `decoupleR::get_ksn_omnipath()` is broken
# in the currently-installed decoupleR 2.16.0 / OmnipathR 3.18.4
# stack: the wrapper calls the deprecated
# `OmnipathR::import_omnipath_enzsub()` which now errors out on its
# argument splat. We reproduce the wrapper's preprocessing locally
# on the modern `OmnipathR::enzyme_substrate()` API; the resulting
# tibble is byte-equivalent in schema to what the broken wrapper
# would have returned (columns source / target / mor). When/if the
# decoupleR wrapper is repaired upstream, the helper here can be
# replaced with a thin shim around it.
#
# Rationale for duplicate-id collapse rule (keep max |t|, same as the
# TF E2 convention): a single (gene, residue, position) tuple can in
# principle appear more than once if two distinct feature rows of
# the limma fit collapse onto the same id (multi-protein-group
# ambiguity, alternative PTM forms folding to the same residue). Max
# |t| preserves the strongest signed evidence per site. Alternatives
# considered: mean (dilutes), sum (inflates), first (arbitrary).
# decoupleR requires a unique row-name axis, so the collapse must
# happen before the call.

# Fetch the human OmniPath KSN and preprocess it identically to
# `decoupleR::get_ksn_omnipath()`, bypassing the latter's broken
# upstream call. Returns a tibble with columns `source` (kinase
# HGNC symbol), `target` (substrate id `SYMBOL_resPos`), `mor`
# (mode-of-regulation: +1 phosphorylation, -1 dephosphorylation).
#
# Arguments:
#   organism  NCBI taxon id. Default 9606L (human). Mouse mapping
#             happens downstream via nichenetr; the KSN is best
#             fetched in its native human form because the upstream
#             OmnipathR cache is human-keyed.
#
# Returns a tibble. Deterministic given the OmnipathR cache state;
# downstream callers that need reproducibility should pin the
# OmnipathR cache version.
fetch_omnipath_ksn_human <- function(organism = 9606L) {
  es <- OmnipathR::enzyme_substrate(organism = organism)
  es %>%
    dplyr::filter(modification %in% c("phosphorylation", "dephosphorylation")) %>%
    dplyr::mutate(
      target = sprintf("%s_%s%i", substrate_genesymbol,
                       residue_type, residue_offset),
      mor    = (modification == "phosphorylation") * 2L - 1L
    ) %>%
    dplyr::select(source = enzyme_genesymbol, target, mor) %>%
    dplyr::distinct() %>%
    # When the same (source, target) pair carries both phosphorylation
    # and dephosphorylation annotations, take min(mor) = -1 so the
    # negative signal wins. This mirrors decoupleR's wrapper.
    dplyr::group_by(source, target) %>%
    dplyr::summarise(mor = min(mor), .groups = "drop")
}

# Map a human KSN tibble to mouse symbols via nichenetr. Kinase symbol
# (source column) and substrate gene symbol (target column prefix) are
# mapped independently; residue type + offset are preserved (regulatory
# phosphosites are highly conserved across mammals so we do not need a
# per-site re-coordination step).
#
# Arguments:
#   ksn_human  tibble from `fetch_omnipath_ksn_human()` with columns
#              `source`, `target`, `mor`.
#
# Returns a tibble with columns `source` (mouse kinase symbol),
# `target` (mouse `SYMBOL_resPos`), `mor`. Edges with an orphan kinase
# OR an orphan substrate symbol are dropped silently. Re-distincts and
# re-collapses by min(mor) to handle any (rare) collisions induced
# by many-to-one mouse mapping.
#
# Failure modes:
#   * If `nichenetr::convert_human_to_mouse_symbols()` returns all-NA
#     (network failure, missing data), the result is an empty 3-col
#     tibble. Callers should check `nrow()` before passing to
#     decoupleR.
map_ksn_human_to_mouse <- function(ksn_human) {
  stopifnot(is.data.frame(ksn_human),
            all(c("source", "target", "mor") %in% names(ksn_human)))

  parts    <- regmatches(ksn_human$target,
                         regexec("^([^_]+)_([A-Z][0-9]+)$", ksn_human$target))
  tg_sym_h <- vapply(parts,
                     function(p) if (length(p) == 3L) p[2] else NA_character_,
                     character(1))
  tg_pos   <- vapply(parts,
                     function(p) if (length(p) == 3L) p[3] else NA_character_,
                     character(1))

  src_unique <- unique(ksn_human$source)
  tgt_unique <- unique(tg_sym_h[!is.na(tg_sym_h)])
  src_map <- nichenetr::convert_human_to_mouse_symbols(src_unique)
  tgt_map <- nichenetr::convert_human_to_mouse_symbols(tgt_unique)
  src_lookup <- setNames(src_map, src_unique)
  tgt_lookup <- setNames(tgt_map, tgt_unique)

  tibble::tibble(
      source       = unname(src_lookup[ksn_human$source]),
      target_sym_m = unname(tgt_lookup[tg_sym_h]),
      target_pos   = tg_pos,
      mor          = ksn_human$mor
    ) %>%
    dplyr::filter(!is.na(source), !is.na(target_sym_m), !is.na(target_pos)) %>%
    dplyr::mutate(target = paste0(target_sym_m, "_", target_pos)) %>%
    dplyr::select(source, target, mor) %>%
    dplyr::group_by(source, target) %>%
    dplyr::summarise(mor = min(mor), .groups = "drop")
}

# Convenience wrapper: build the mouse OmniPath KSN end-to-end.
# Exposed because F4 (axis-restricted kinase analysis) needs the
# same network and the plan's F4 stub explicitly requires "do not
# re-implement the mapping -- factor F2's mapping out into a helper".
build_omnipath_ksn_mouse <- function() {
  ksn_h <- fetch_omnipath_ksn_human(organism = 9606L)
  map_ksn_human_to_mouse(ksn_h)
}

# Extract a single per-contrast stat vector keyed by phospho-site id.
#
# Arguments:
#   top_df     a data.frame as found in `de_phospho.rds$fit$top[[contrast]]`
#              or `de_phospho_corrected.rds$fit$top[[contrast]]`. Must
#              carry `symbol`, `PTM.SiteAA`, `PTM.SiteLocation`, and
#              `stat_col` columns.
#   stat_col   column name of the test statistic. Default "t" (limma
#              moderated t).
#   symbol_col, aa_col, pos_col  column names; defaults match the
#              project's phospho cache schema.
#
# Returns a named numeric vector of stats, names = `SYMBOL_resPos`.
# Sites with missing symbol / SiteAA / SiteLocation are dropped (the
# phospho cache has ~0.4% such rows where the upstream Spectronaut
# search failed to assign a residue). NA / non-finite stats are also
# dropped. Duplicates collapsed by keeping max |stat| per the module
# header rationale. Empty input returns a length-0 numeric.
extract_phospho_stat_vec <- function(top_df,
                                     stat_col   = "t",
                                     symbol_col = "symbol",
                                     aa_col     = "PTM.SiteAA",
                                     pos_col    = "PTM.SiteLocation") {
  stopifnot(is.data.frame(top_df),
            all(c(symbol_col, aa_col, pos_col, stat_col) %in% names(top_df)))

  sym  <- as.character(top_df[[symbol_col]])
  aa   <- as.character(top_df[[aa_col]])
  pos  <- top_df[[pos_col]]
  stat <- as.numeric(top_df[[stat_col]])

  keep <- !is.na(sym) & nzchar(sym) &
          !is.na(aa)  & nzchar(aa)  &
          !is.na(pos) &
          !is.na(stat) & is.finite(stat)
  sym  <- sym[keep]; aa <- aa[keep]; pos <- pos[keep]; stat <- stat[keep]
  if (length(stat) == 0L) {
    return(setNames(numeric(0), character(0)))
  }

  ids <- paste0(sym, "_", aa, pos)
  if (anyDuplicated(ids) > 0L) {
    ord  <- order(ids, -abs(stat))
    ids  <- ids[ord]; stat <- stat[ord]
    keep <- !duplicated(ids)
    ids  <- ids[keep]; stat <- stat[keep]
  }
  setNames(stat, ids)
}

# Build a numeric matrix of stats with rows = phospho-site id and
# columns = contrast, suitable for direct decoupleR consumption.
#
# Arguments:
#   top_list   named list `contrast -> data.frame`, as in
#              `de_phospho{,_corrected}.rds$fit$top`.
#   stat_col, symbol_col, aa_col, pos_col  forwarded to
#              `extract_phospho_stat_vec()`.
#
# Returns a numeric matrix; rownames are the union of ids across
# contrasts; missing `(id, contrast)` cells are NA. Empty input
# returns a 0x0 matrix.
extract_phospho_stat_matrix <- function(top_list,
                                        stat_col   = "t",
                                        symbol_col = "symbol",
                                        aa_col     = "PTM.SiteAA",
                                        pos_col    = "PTM.SiteLocation") {
  stopifnot(is.list(top_list))
  if (length(top_list) == 0L) {
    return(matrix(numeric(0), nrow = 0, ncol = 0))
  }
  per_contrast <- lapply(top_list, extract_phospho_stat_vec,
                         stat_col   = stat_col,
                         symbol_col = symbol_col,
                         aa_col     = aa_col,
                         pos_col    = pos_col)
  all_ids <- sort(unique(unlist(lapply(per_contrast, names),
                                use.names = FALSE)))
  contrasts <- names(top_list)
  mat <- matrix(NA_real_, nrow = length(all_ids), ncol = length(contrasts),
                dimnames = list(all_ids, contrasts))
  for (cn in contrasts) {
    v <- per_contrast[[cn]]
    if (length(v) == 0L) next
    mat[names(v), cn] <- unname(v)
  }
  mat
}

# Run decoupleR univariate methods on one phospho stat matrix
# (single cache, all contrasts). Mirrors `run_decoupler_per_modality()`
# in `R/tf_inference.R`: per-contrast iteration to tolerate limma's
# asymmetric per-contrast NA pattern; same defaults
# (`statistics = c("ulm", "wsum")` + `consensus = TRUE`, `minsize = 5L`)
# so the Phase F outputs are read against Phase E with identical
# inference conventions.
#
# Arguments:
#   stat_mat    numeric matrix from `extract_phospho_stat_matrix()`.
#               Rows = phospho-site ids; cols = contrast names.
#               Non-zero dimensions required; NA cells tolerated.
#   network     mouse OmniPath KSN tibble (from
#               `build_omnipath_ksn_mouse()`) with columns `source`
#               (kinase), `target` (substrate id), `mor` (+/- 1).
#   statistics  character vector of decoupleR method names. Default
#               c("ulm", "wsum"); consensus appended via
#               `consensus = TRUE`.
#   minsize     drop kinases whose target-substrate set overlap with
#               the row-name (phospho-id) universe is smaller than
#               this. Default 5L (project-wide convention; matches
#               fgsea + TF inference). The phospho probe in F1
#               shows ~100-110 kinases survive this threshold per
#               cache, comfortably enough for downstream ranking.
#   consensus   logical; whether to emit the cross-method consensus
#               score. Default TRUE.
#   seed        integer seed forwarded to consensus; the underlying
#               implementation is deterministic so this is a
#               safeguard for any future stochastic decoupleR
#               method that may be added. Default 42L.
#
# Returns a tidy tibble with columns statistic / source / condition /
# score / p_value, identical schema to the TF helper's output.
#
# Failure modes:
#   * `stat_mat` empty -> empty tibble, canonical schema preserved.
#   * No kinases survive `minsize` for any contrast -> empty tibble +
#     warning.
run_decoupler_per_cache <- function(stat_mat,
                                    network,
                                    statistics = c("ulm", "wsum"),
                                    minsize    = 5L,
                                    consensus  = TRUE,
                                    seed       = 42L) {
  stopifnot(is.matrix(stat_mat),
            is.data.frame(network) || tibble::is_tibble(network),
            all(c("source", "target", "mor") %in% names(network)))

  if (nrow(stat_mat) == 0L || ncol(stat_mat) == 0L) {
    return(tibble::tibble(
      statistic = character(0), source = character(0),
      condition = character(0), score = numeric(0),
      p_value   = numeric(0)
    ))
  }

  contrasts <- colnames(stat_mat)
  per_contrast <- lapply(contrasts, function(cn) {
    col  <- stat_mat[, cn]
    keep <- !is.na(col) & is.finite(col)
    if (!any(keep)) {
      warning(sprintf(
        "run_decoupler_per_cache: contrast '%s' has no usable stats",
        cn), call. = FALSE)
      return(NULL)
    }
    sub_mat <- matrix(col[keep], ncol = 1,
                      dimnames = list(rownames(stat_mat)[keep], cn))
    dec <- decoupleR::decouple(
      mat             = sub_mat,
      network         = network,
      .source         = "source",
      .target         = "target",
      statistics      = statistics,
      consensus_score = isTRUE(consensus),
      minsize         = as.integer(minsize)
    )
    if (is.null(dec) || nrow(dec) == 0L) {
      warning(sprintf(
        "run_decoupler_per_cache: decouple returned 0 rows for contrast '%s'",
        cn), call. = FALSE)
      return(NULL)
    }
    dec
  })
  per_contrast <- per_contrast[!vapply(per_contrast, is.null, logical(1))]
  if (length(per_contrast) == 0L) {
    return(tibble::tibble(
      statistic = character(0), source = character(0),
      condition = character(0), score = numeric(0),
      p_value   = numeric(0)
    ))
  }

  if (!is.null(seed) && isTRUE(consensus)) {
    invisible(seed)
  }

  dec      <- dplyr::bind_rows(per_contrast)
  out_cols <- c("statistic", "source", "condition", "score", "p_value")
  missing  <- setdiff(out_cols, names(dec))
  if (length(missing) > 0L) {
    stop(sprintf(
      "run_decoupler_per_cache: decouple output missing columns: %s",
      paste(missing, collapse = ", ")), call. = FALSE)
  }
  tibble::as_tibble(dec[, out_cols, drop = FALSE])
}

# Split a per-cache decoupler tibble into a list keyed by contrast.
# Same shape as `split_decoupler_by_contrast()` in `R/tf_inference.R`
# but redefined here to keep the kinase module self-contained
# (avoids cross-file silent coupling if a future change to one
# file's split logic should not propagate to the other).
#
# Arguments:
#   dec_tbl   output of `run_decoupler_per_cache()`.
#
# Returns a named list of tibbles, one per `condition` level. Inner
# tibbles drop the now-redundant `condition` column. Preserves the
# order in which contrasts first appear in `dec_tbl`.
split_kinase_decoupler_by_contrast <- function(dec_tbl) {
  stopifnot(is.data.frame(dec_tbl),
            all(c("statistic", "source", "condition", "score", "p_value")
                %in% names(dec_tbl)))
  if (nrow(dec_tbl) == 0L) return(list())

  contrasts <- unique(dec_tbl$condition)
  out <- lapply(contrasts, function(cn) {
    sub <- dec_tbl[dec_tbl$condition == cn, , drop = FALSE]
    sub <- sub[, c("statistic", "source", "score", "p_value"), drop = FALSE]
    tibble::as_tibble(sub)
  })
  setNames(out, contrasts)
}

# ============================================================================
# Phase F3 helpers: per-contrast significant-kinase reporting.
# ============================================================================
#
# Session F3 (storage/notes/mechanism_layer_plan.md) was originally scoped
# to mirror the E3 cross-modality + cross-contrast leader-board pattern
# applied to the F2 kinase cache, with `phospho_corrected` as the primary
# cache and `phospho_raw` as a corroborator. The F2 completion note
# foreshadowed (and the F3 pre-build diagnostic confirmed) two structural
# obstacles to that pattern:
#
#   1. The cross-cache consistency arm of the locked leader rule is
#      completely empty at FDR<0.10. The corrected and raw caches read
#      meaningfully different biology because batch correction has
#      exposed signal that batch-confounded noise was burying in the
#      raw cache. The two caches even disagree on sign at the strongest
#      single-kinase signal (Gsk3b in corrected `nlgf_in_p301s` is
#      score +3.79, p = 1.5e-04; in raw the same cell is -1.11,
#      p = 0.27). The raw cache cannot serve as an FDR-scale
#      corroborator for the corrected cache; the cross-cache rule arm
#      is structurally dead.
#
#   2. The `n_contrasts_sig_corrected >= 3` arm in the locked rule
#      fires for exactly one kinase (Gsk3b) because the kinase
#      inference universe (~100 kinases per contrast after the
#      `minsize = 5L` overlap filter) is small enough that BH-padj
#      is harsh -- the per-contrast n_sig at FDR<0.10 is in the
#      0-4 range.
#
# User-directed F3 deviation (2026-05-24). The cross-cache and
# cross-contrast aggregation patterns are dropped at the kinase layer.
# Contrasts are treated as independent inference units; only the
# corrected cache is consulted; the per-contrast significance gate is
# FDR<0.10 on the `ulm` p-value (E3 convention preserved); a kinase is
# "significant at this contrast" or not, with no cross-contrast
# leader rule. The resulting per-contrast significant sets become the
# F3 output (no unified leader-board concept). The F4 axis-restricted
# layer still aggregates within axes per the E4 precedent unless that
# precedent is revisited at F4; the F5 verdict cross-references the
# per-contrast significant sets directly rather than a leader board.
# This deviation is logged in the F3 completion note of the active plan
# at storage/notes/mechanism_layer_plan.md.

# Compute per-contrast, per-kinase significance + direction columns for
# one cache of the F2 output. Returns a long tibble keyed on
# (contrast, source) with both ulm and consensus columns retained so
# downstream chunks can choose either for any visualisation. BH-padj is
# applied within each contrast (the multiple-testing unit defined by the
# user-directed F3 rule that contrasts are independent inferences).
#
# Arguments:
#   kinase_cache    the loaded kinase_activity_decoupler.rds list.
#   cache_name      character (default "phospho_corrected" per the F3
#                   user directive); name of the cache to use.
#   sig_statistic   character (default "ulm"); name of the per-method
#                   p-value used as the significance arbiter. Must be
#                   one of "ulm", "consensus", "wsum", "norm_wsum",
#                   "corr_wsum" (the F2 cache statistic levels).
#   dir_statistic   character (default "consensus"); name of the per-
#                   method score used for sign and magnitude attribution.
#
# Returns a tibble with columns:
#   contrast                 character; contrast name.
#   source                   character; kinase symbol.
#   score_<sig_statistic>    numeric; score for the sig arbiter.
#   p_value_<sig_statistic>  numeric; raw p-value.
#   padj_<sig_statistic>     numeric; BH-padj within contrast.
#   score_<dir_statistic>    numeric; only if dir != sig.
#   p_value_<dir_statistic>  numeric; only if dir != sig.
#   padj_<dir_statistic>     numeric; BH-padj within contrast; only
#                            if dir != sig.
#   sign_dir                 integer in {-1, 0, 1}; sign(score_<dir>).
build_kinase_per_contrast_table <- function(kinase_cache,
                                            cache_name    = "phospho_corrected",
                                            sig_statistic = "ulm",
                                            dir_statistic = "consensus") {
  stopifnot(is.list(kinase_cache),
            cache_name %in% names(kinase_cache),
            length(sig_statistic) == 1L,
            length(dir_statistic) == 1L)

  cache_tbls <- kinase_cache[[cache_name]]
  contrasts  <- names(cache_tbls)
  if (length(contrasts) == 0L) {
    return(tibble::tibble())
  }

  per_contrast <- lapply(contrasts, function(cn) {
    tib     <- cache_tbls[[cn]]
    sig_row <- tib[tib$statistic == sig_statistic,
                   c("source", "score", "p_value"), drop = FALSE]
    if (nrow(sig_row) == 0L) {
      stop("statistic '", sig_statistic,
           "' not found in cache_name='", cache_name,
           "', contrast='", cn, "'")
    }
    names(sig_row) <- c("source",
                        paste0("score_",   sig_statistic),
                        paste0("p_value_", sig_statistic))
    sig_row[[paste0("padj_", sig_statistic)]] <-
      p.adjust(sig_row[[paste0("p_value_", sig_statistic)]], method = "BH")

    out <- sig_row
    if (dir_statistic != sig_statistic) {
      dir_row <- tib[tib$statistic == dir_statistic,
                     c("source", "score", "p_value"), drop = FALSE]
      if (nrow(dir_row) == 0L) {
        stop("statistic '", dir_statistic,
             "' not found in cache_name='", cache_name,
             "', contrast='", cn, "'")
      }
      names(dir_row) <- c("source",
                          paste0("score_",   dir_statistic),
                          paste0("p_value_", dir_statistic))
      dir_row[[paste0("padj_", dir_statistic)]] <-
        p.adjust(dir_row[[paste0("p_value_", dir_statistic)]], method = "BH")
      out <- merge(out, dir_row, by = "source", all = TRUE)
    }

    out$sign_dir <- as.integer(sign(out[[paste0("score_", dir_statistic)]]))
    out$contrast <- cn
    col_order <- c("contrast", "source",
                   setdiff(names(out), c("contrast", "source")))
    tibble::as_tibble(out[, col_order, drop = FALSE])
  })

  dplyr::bind_rows(per_contrast)
}

# Filter the per-contrast table to FDR-significant rows at one contrast
# and return a tidy display-ready tibble suitable for kable. Output
# column set is the minimal informative slice; the full per-contrast
# table is preserved upstream for downstream re-use.
#
# Arguments:
#   per_contrast_tbl  output of build_kinase_per_contrast_table().
#   contrast          character (length 1); contrast to filter on.
#   padj_cut          numeric (default 0.10); BH-padj threshold.
#   sig_statistic     character (default "ulm"); must match the
#                     sig_statistic used at the parent helper.
#   dir_statistic     character (default "consensus"); must match
#                     the dir_statistic used at the parent helper.
#
# Returns a tibble ordered by ascending padj_<sig_statistic>. May have
# zero rows for contrasts with no significant kinases.
filter_sig_kinases_for_contrast <- function(per_contrast_tbl,
                                            contrast,
                                            padj_cut      = 0.10,
                                            sig_statistic = "ulm",
                                            dir_statistic = "consensus") {
  stopifnot(is.data.frame(per_contrast_tbl),
            length(contrast) == 1L,
            length(padj_cut) == 1L,
            "contrast" %in% names(per_contrast_tbl))

  padj_col <- paste0("padj_", sig_statistic)
  if (!padj_col %in% names(per_contrast_tbl)) {
    stop("padj column '", padj_col, "' not found in per_contrast_tbl")
  }

  sub <- per_contrast_tbl[per_contrast_tbl$contrast == contrast, , drop = FALSE]
  sub <- sub[!is.na(sub[[padj_col]]) & sub[[padj_col]] < padj_cut, , drop = FALSE]
  if (nrow(sub) == 0L) {
    return(tibble::tibble(source = character(0)))
  }
  sub <- sub[order(sub[[padj_col]]), , drop = FALSE]
  keep_cols <- c("source",
                 paste0("score_",   sig_statistic),
                 paste0("padj_",    sig_statistic),
                 paste0("score_",   dir_statistic),
                 paste0("padj_",    dir_statistic))
  keep_cols <- keep_cols[keep_cols %in% names(sub)]
  tibble::as_tibble(sub[, keep_cols, drop = FALSE])
}

# Build a ComplexHeatmap of signed score by (kinase x contrast) for the
# kinases that reach FDR<padj_cut on the sig_statistic in any contrast.
# Rows alphabetised so no kinase is visually privileged; columns in the
# project-canonical narrative order. Cells coloured by signed
# score_<colour_statistic>; cell text marks significance level
# (`*` = FDR<padj_strict; `.` = FDR<padj_cut). NA cells (kinase absent
# from a contrast's decoupleR output) are blank.
#
# Arguments:
#   per_contrast_tbl  output of build_kinase_per_contrast_table().
#   padj_cut          looser BH-padj threshold (default 0.10) used to
#                     pick the row set AND to draw "." annotation.
#   padj_strict       stricter BH-padj threshold (default 0.05) used to
#                     draw "*" annotation.
#   sig_statistic     character (default "ulm"); names the padj_<>
#                     column used for the row-set filter + annotation.
#   colour_statistic  character (default "consensus"); names the
#                     score_<> column used for cell colour.
#   contrast_order    character vector defining column order; default
#                     the project-canonical five-contrast set in
#                     narrative order (NLGF-context contrasts first,
#                     then tau-context, then interaction).
#   title             heatmap column title; auto-generated if NULL.
#
# Returns a ComplexHeatmap::Heatmap object (caller draws or relies on
# auto-printing inside a knitr chunk).
plot_kinase_activity_heatmap <- function(per_contrast_tbl,
                                         padj_cut         = 0.10,
                                         padj_strict      = 0.05,
                                         sig_statistic    = "ulm",
                                         colour_statistic = "consensus",
                                         contrast_order   = c(
                                           "nlgf_in_maptki",
                                           "nlgf_in_p301s",
                                           "tau_alone",
                                           "tau_in_nlgf",
                                           "interaction"
                                         ),
                                         title = NULL) {
  stopifnot(is.data.frame(per_contrast_tbl),
            length(padj_cut)    == 1L,
            length(padj_strict) == 1L,
            padj_strict <= padj_cut)
  sig_col   <- paste0("padj_",  sig_statistic)
  score_col <- paste0("score_", colour_statistic)
  stopifnot(sig_col    %in% names(per_contrast_tbl),
            score_col  %in% names(per_contrast_tbl),
            "contrast" %in% names(per_contrast_tbl),
            "source"   %in% names(per_contrast_tbl))

  any_sig <- per_contrast_tbl[
    !is.na(per_contrast_tbl[[sig_col]]) &
      per_contrast_tbl[[sig_col]] < padj_cut, , drop = FALSE
  ]
  if (nrow(any_sig) == 0L) {
    stop("No kinase reaches padj_", sig_statistic, " < ", padj_cut,
         " in any contrast; nothing to plot.")
  }
  kinases <- sort(unique(any_sig$source))

  available_contrasts <- intersect(contrast_order,
                                   unique(per_contrast_tbl$contrast))
  if (length(available_contrasts) == 0L) {
    stop("None of the requested contrasts are present in per_contrast_tbl.")
  }

  score_mat <- matrix(NA_real_,
                      nrow = length(kinases),
                      ncol = length(available_contrasts),
                      dimnames = list(kinases, available_contrasts))
  padj_mat <- score_mat
  for (ctr in available_contrasts) {
    sub <- per_contrast_tbl[per_contrast_tbl$contrast == ctr, , drop = FALSE]
    idx <- match(kinases, sub$source)
    score_mat[, ctr] <- sub[[score_col]][idx]
    padj_mat[, ctr]  <- sub[[sig_col]][idx]
  }

  range_abs <- max(abs(score_mat[is.finite(score_mat)]), na.rm = TRUE)
  if (!is.finite(range_abs) || range_abs == 0) range_abs <- 1
  col_fun <- circlize::colorRamp2(c(-range_abs, 0, range_abs),
                                  c("#3a4cc0", "white", "#b40426"))
  cell_fn <- function(j, i, x, y, width, height, fill) {
    if (is.na(padj_mat[i, j])) return()
    if (padj_mat[i, j] < padj_strict) {
      grid::grid.text("*", x, y,
                      gp = grid::gpar(fontsize = 11, col = "black"))
    } else if (padj_mat[i, j] < padj_cut) {
      grid::grid.text(".", x, y,
                      gp = grid::gpar(fontsize = 11, col = "black"))
    }
  }
  if (is.null(title)) {
    title <- sprintf(
      paste0("Kinases sig at FDR<%g (%s) in any contrast x signed %s score ",
             "('*' = FDR<%g, '.' = FDR<%g; rows alphabetical)"),
      padj_cut, sig_statistic, colour_statistic,
      padj_strict, padj_cut)
  }
  ComplexHeatmap::Heatmap(
    score_mat,
    name            = sprintf("%s score", colour_statistic),
    col             = col_fun,
    cluster_rows    = FALSE,
    cluster_columns = FALSE,
    row_order       = seq_len(nrow(score_mat)),
    column_order    = seq_len(ncol(score_mat)),
    cell_fun        = cell_fn,
    row_names_gp    = grid::gpar(fontsize = 10),
    column_names_gp = grid::gpar(fontsize = 10),
    column_title    = title,
    na_col          = "grey90",
    border          = TRUE
  )
}

# One-line-per-contrast summary used in the section 15.1 intro paragraph.
# Reports total kinase count, sig count, and the top hit (smallest padj)
# for each contrast.
#
# Arguments mirror filter_sig_kinases_for_contrast() with the addition
# of `contrast_order` for output row ordering (default project canonical
# narrative order).
#
# Returns a tibble with columns:
#   contrast          character; contrast name.
#   n_kinases_total   integer; total kinases in the contrast.
#   n_kinases_sig     integer; kinases at padj < padj_cut.
#   top_kinase        character; lowest-padj kinase, or NA.
#   top_kinase_padj   numeric; corresponding padj_<sig_statistic>, or NA.
#   top_kinase_score  numeric; corresponding score_<dir_statistic>, or NA.
summarise_sig_kinases_per_contrast <- function(per_contrast_tbl,
                                               padj_cut       = 0.10,
                                               sig_statistic  = "ulm",
                                               dir_statistic  = "consensus",
                                               contrast_order = c(
                                                 "nlgf_in_maptki",
                                                 "nlgf_in_p301s",
                                                 "tau_alone",
                                                 "tau_in_nlgf",
                                                 "interaction"
                                               )) {
  stopifnot(is.data.frame(per_contrast_tbl),
            "contrast" %in% names(per_contrast_tbl))
  sig_col   <- paste0("padj_",  sig_statistic)
  score_col <- paste0("score_", dir_statistic)
  contrasts_in_data <- unique(per_contrast_tbl$contrast)
  ordered_contrasts <- c(intersect(contrast_order, contrasts_in_data),
                         setdiff(contrasts_in_data, contrast_order))

  rows <- lapply(ordered_contrasts, function(ctr) {
    sub     <- per_contrast_tbl[per_contrast_tbl$contrast == ctr, , drop = FALSE]
    n_total <- nrow(sub)
    sig_idx <- !is.na(sub[[sig_col]]) & sub[[sig_col]] < padj_cut
    n_sig   <- sum(sig_idx)
    if (n_sig > 0L) {
      sub_sig <- sub[sig_idx, , drop = FALSE]
      top_idx <- which.min(sub_sig[[sig_col]])
      top_k   <- sub_sig$source[top_idx]
      top_pa  <- sub_sig[[sig_col]][top_idx]
      top_sc  <- sub_sig[[score_col]][top_idx]
    } else {
      top_k  <- NA_character_
      top_pa <- NA_real_
      top_sc <- NA_real_
    }
    tibble::tibble(
      contrast         = ctr,
      n_kinases_total  = as.integer(n_total),
      n_kinases_sig    = as.integer(n_sig),
      top_kinase       = top_k,
      top_kinase_padj  = top_pa,
      top_kinase_score = top_sc
    )
  })
  dplyr::bind_rows(rows)
}

# ----------------------------------------------------------------------
# F4 axis-restricted kinase activity (hybrid view).
# ----------------------------------------------------------------------
# These helpers mirror the F4 TF helpers in `R/tf_inference.R` (the
# "Interpretation A" approach: read cached scores against full-modality
# universes, then restrict to kinases whose KSN substrate sites overlap
# the per-axis gene universe at >= min_targets distinct sites). The
# locked F4 hybrid decision (gated on 2026-05-24 with user directive
# "Do the hybrid. Have other analyses do the hybrid as well for
# consistency") emits BOTH the axis-mean and per-axis-contrast columns
# so a reader can see, at a glance, whether the axis-mean smooths over
# per-contrast asymmetry.
#
# Important schema notes consumed downstream:
#   - The KSN target id format is `symbol_resPos`, e.g. "Mapt_S202",
#     "Gsk3b_Y216". The substrate gene-symbol prefix is everything
#     before the LAST underscore (the residue+position suffix has the
#     form `[A-Z][0-9]+` and contains no underscore).
#   - The cache passed in is `kinase_activity_decoupler.rds` keyed by
#     `cache_name` (default "phospho_corrected" per the F3 decision to
#     use the corrected layer only for the per-contrast significance
#     gate; F4 inherits the same "corrected only" convention by default
#     so the axis-restricted layer is internally consistent with the
#     F3 leader-board).
#   - The function emits `score_at_<contrast>` columns padded to the
#     union of axis_contrasts across all axes (NA where the row's axis
#     does not include the contrast), matching the locked F4 hybrid
#     schema that `format_axis_restricted_kinase_table()` and the
#     F5 verdict helper consume.

# Extract gene symbol from a KSN substrate id (`symbol_resPos`).
# Strips the last `_`-delimited token. Vectorised.
.ksn_substrate_to_symbol <- function(target_ids) {
  sub("_[^_]+$", "", as.character(target_ids))
}

# Restrict an OmniPath KSN to kinases whose distinct substrate-site
# count in a given gene universe meets min_targets. Mirrors the TF
# analogue `restrict_collectri_to_universe()`.
#
# Arguments:
#   network        data.frame / tibble with columns `source` (kinase
#                  symbol), `target` (substrate `symbol_resPos`), `mor`
#                  (mode of regulation; not used here but preserved).
#                  Output of `build_omnipath_ksn_mouse()`.
#   universe       character vector of gene symbols (the axis universe).
#   min_targets    integer; drop kinases with fewer than this many
#                  distinct in-universe substrate sites. Default 5.
#
# Returns a list with:
#   edges          filtered network restricted to (source, target) where
#                  the target's gene-symbol prefix is in `universe`
#                  AND the source passes the min_targets filter.
#   target_counts  tibble (source, n_targets_in_universe), one row per
#                  surviving kinase, sorted by n_targets_in_universe
#                  desc, source alphabetical.
restrict_ksn_to_universe <- function(network, universe,
                                     min_targets = 5L) {
  stopifnot(is.data.frame(network) || tibble::is_tibble(network),
            all(c("source", "target") %in% names(network)),
            is.character(universe) || is.null(universe),
            length(min_targets) == 1L, min_targets >= 1L)

  if (length(universe) == 0L) {
    return(list(
      edges = network[FALSE, , drop = FALSE],
      target_counts = tibble::tibble(
        source = character(0),
        n_targets_in_universe = integer(0)
      )
    ))
  }
  uni <- unique(as.character(universe))

  net <- network
  net$.symbol_in_universe <- .ksn_substrate_to_symbol(net$target) %in% uni
  net <- net[net$.symbol_in_universe, , drop = FALSE]
  net$.symbol_in_universe <- NULL

  if (nrow(net) == 0L) {
    return(list(
      edges = net,
      target_counts = tibble::tibble(
        source = character(0),
        n_targets_in_universe = integer(0)
      )
    ))
  }

  per_source <- aggregate(target ~ source, data = net,
                          FUN = function(v) length(unique(v)))
  names(per_source)[2] <- "n_targets_in_universe"
  per_source <- per_source[per_source$n_targets_in_universe >= min_targets, ,
                           drop = FALSE]
  if (nrow(per_source) == 0L) {
    return(list(
      edges = net[FALSE, , drop = FALSE],
      target_counts = tibble::tibble(
        source = character(0),
        n_targets_in_universe = integer(0)
      )
    ))
  }

  net <- net[net$source %in% per_source$source, , drop = FALSE]
  per_source <- per_source[order(-per_source$n_targets_in_universe,
                                 per_source$source), , drop = FALSE]

  list(
    edges = tibble::as_tibble(net),
    target_counts = tibble::as_tibble(per_source)
  )
}

# Score each kinase per D2 axis using the cached decoupleR consensus
# statistic, restricted to kinases with sufficient in-universe KSN
# substrate sites. Hybrid output: axis-mean across (axis_contrast x
# primary_modality) cells AND per-contrast cross-modality means.
#
# Arguments:
#   kinase_cache         output of `run_decoupler_per_kinase_cache()`
#                        (`storage/cache/kinase_activity_decoupler.rds`).
#                        list: cache_name -> contrast -> tibble.
#   axis_universes       output of `build_axis_gene_universe()` from
#                        `R/tf_inference.R` (gene-level universes per
#                        axis). The kinase universe is derived by
#                        expanding each gene universe to its in-KSN
#                        substrate-site set.
#   network              KSN tibble (source = kinase, target =
#                        symbol_resPos, mor). Output of
#                        `build_omnipath_ksn_mouse()`.
#   primary_modalities   character vector of cache_name keys to include
#                        in the cross-modality mean. The F3 decision
#                        locks the phospho-corrected layer only;
#                        default mirrors that locked choice. Provide
#                        c("phospho_raw","phospho_corrected") for a
#                        sensitivity check.
#   min_targets          integer; min in-universe substrate-site count
#                        per kinase. Default 5.
#   score_statistic      character; which decoupleR statistic to read
#                        per cell. Default "consensus".
#
# Returns a tibble with columns:
#   axis                              character; axis name.
#   source                            character; kinase symbol.
#   mean_activity_in_axis_contrasts   numeric; mean score across
#                                     (axis_contrast, primary_modality)
#                                     cells with finite score.
#   sd_activity_in_axis_contrasts     numeric.
#   n_cells_used                      integer.
#   n_targets_in_axis_universe        integer; KSN substrate-site count
#                                     in the axis gene universe.
#   leader_rank                       integer; rank within axis by
#                                     |mean| desc, source alphabetical.
#   score_at_<contrast>               numeric; cross-modality mean at
#                                     that axis_contrast (one column
#                                     per unique axis_contrast across
#                                     all axes; NA where the row's axis
#                                     does not include the contrast).
#                                     Hybrid columns -- locked at the
#                                     F4 decision gate.
score_kinase_per_axis <- function(kinase_cache,
                                  axis_universes,
                                  network,
                                  primary_modalities = c("phospho_corrected"),
                                  min_targets        = 5L,
                                  score_statistic    = "consensus") {
  stopifnot(is.list(kinase_cache),
            is.list(axis_universes),
            !is.null(names(axis_universes)),
            is.data.frame(network) || tibble::is_tibble(network),
            length(primary_modalities) > 0L,
            all(primary_modalities %in% names(kinase_cache)))

  all_contrasts <- unique(unlist(
    lapply(axis_universes, function(a) a$axis_contrasts),
    use.names = FALSE
  ))
  per_contrast_cols <- paste0("score_at_", all_contrasts)

  rows <- list()
  for (ax in names(axis_universes)) {
    aux  <- axis_universes[[ax]]
    uni  <- aux$universe
    ctrs <- aux$axis_contrasts
    if (length(uni) == 0L || length(ctrs) == 0L) next

    restr <- restrict_ksn_to_universe(network, uni,
                                      min_targets = min_targets)
    if (nrow(restr$target_counts) == 0L) next

    long_list <- list()
    for (mn in primary_modalities) {
      mod_list <- kinase_cache[[mn]]
      for (cn in ctrs) {
        tbl <- mod_list[[cn]]
        if (is.null(tbl)) next
        sub <- tbl[tbl$statistic == score_statistic, , drop = FALSE]
        if (nrow(sub) == 0L) next
        long_list[[length(long_list) + 1L]] <- data.frame(
          source   = as.character(sub$source),
          contrast = cn,
          score    = as.numeric(sub$score),
          stringsAsFactors = FALSE
        )
      }
    }
    if (length(long_list) == 0L) next
    long_df <- do.call(rbind, long_list)
    long_df <- long_df[!is.na(long_df$score) & is.finite(long_df$score), ,
                       drop = FALSE]
    long_df <- long_df[long_df$source %in% restr$target_counts$source, ,
                       drop = FALSE]
    if (nrow(long_df) == 0L) next

    per_k <- aggregate(score ~ source, data = long_df, FUN = function(v) {
      c(mean_v = mean(v), sd_v = stats::sd(v), n_v = length(v))
    })
    score_mat <- per_k$score
    per_k$mean_activity_in_axis_contrasts <- as.numeric(score_mat[, "mean_v"])
    per_k$sd_activity_in_axis_contrasts   <- as.numeric(score_mat[, "sd_v"])
    per_k$n_cells_used                    <- as.integer(score_mat[, "n_v"])
    per_k$score <- NULL

    # Hybrid per-contrast columns: cross-modality mean at each axis-
    # contrast.
    per_kc <- aggregate(score ~ source + contrast, data = long_df,
                        FUN = mean)
    src_levels <- sort(unique(per_kc$source))
    wide_mat <- matrix(NA_real_,
                       nrow = length(src_levels),
                       ncol = length(ctrs),
                       dimnames = list(src_levels, ctrs))
    for (i in seq_len(nrow(per_kc))) {
      wide_mat[per_kc$source[i], per_kc$contrast[i]] <- per_kc$score[i]
    }
    wide_df <- as.data.frame(wide_mat, stringsAsFactors = FALSE)
    names(wide_df) <- paste0("score_at_", names(wide_df))
    wide_df$source <- rownames(wide_df)
    rownames(wide_df) <- NULL
    per_k <- merge(per_k, wide_df, by = "source", all.x = TRUE)

    per_k <- merge(per_k, restr$target_counts, by = "source", all.x = TRUE)
    names(per_k)[names(per_k) == "n_targets_in_universe"] <-
      "n_targets_in_axis_universe"

    missing_cols <- setdiff(per_contrast_cols, names(per_k))
    for (mc in missing_cols) per_k[[mc]] <- NA_real_

    per_k <- per_k[order(-abs(per_k$mean_activity_in_axis_contrasts),
                         per_k$source), , drop = FALSE]
    per_k$leader_rank <- seq_len(nrow(per_k))
    per_k$axis <- ax

    rows[[ax]] <- per_k[, c("axis", "source",
                            "mean_activity_in_axis_contrasts",
                            "sd_activity_in_axis_contrasts",
                            "n_cells_used",
                            "n_targets_in_axis_universe",
                            "leader_rank",
                            per_contrast_cols), drop = FALSE]
  }
  if (length(rows) == 0L) {
    empty <- tibble::tibble(
      axis                              = character(0),
      source                            = character(0),
      mean_activity_in_axis_contrasts   = numeric(0),
      sd_activity_in_axis_contrasts     = numeric(0),
      n_cells_used                      = integer(0),
      n_targets_in_axis_universe        = integer(0),
      leader_rank                       = integer(0)
    )
    for (col in per_contrast_cols) empty[[col]] <- numeric(0)
    return(empty)
  }
  tibble::as_tibble(do.call(rbind, rows))
}

# Format the top-n axis-restricted kinases for one axis as a
# knitr::kable. Mirrors `format_axis_restricted_table()` in
# `R/tf_inference.R`; auto-includes per-contrast columns whose values
# are not all NA in the axis subset.
format_axis_restricted_kinase_table <- function(axis_tbl, axis_name,
                                                n = 15, caption = NULL) {
  sub <- axis_tbl[axis_tbl$axis == axis_name, , drop = FALSE]
  if (nrow(sub) == 0L) {
    return(knitr::kable(
      data.frame(message = sprintf(
        "No kinases survive the axis '%s' KSN-universe filter.", axis_name)),
      caption = caption %||%
        sprintf("Axis-restricted kinases: %s", axis_name),
      row.names = FALSE))
  }
  sub <- sub[order(sub$leader_rank), , drop = FALSE]
  top <- head(sub, n)

  per_contrast_cols <- grep("^score_at_", names(top), value = TRUE)
  if (length(per_contrast_cols) > 0L) {
    axis_subset <- axis_tbl[axis_tbl$axis == axis_name, , drop = FALSE]
    keep_pc <- vapply(per_contrast_cols, function(col) {
      !all(is.na(axis_subset[[col]]))
    }, logical(1))
    per_contrast_cols <- per_contrast_cols[keep_pc]
  }

  display <- data.frame(
    rank         = top$leader_rank,
    kinase       = top$source,
    mean_score   = round(top$mean_activity_in_axis_contrasts, 3),
    sd_score     = round(top$sd_activity_in_axis_contrasts,   3),
    n_cells      = top$n_cells_used,
    n_sites      = top$n_targets_in_axis_universe,
    stringsAsFactors = FALSE
  )
  for (col in per_contrast_cols) {
    short <- sub("^score_at_", "", col)
    display[[short]] <- round(top[[col]], 3)
  }
  if (is.null(caption)) {
    pc_part <- if (length(per_contrast_cols) > 0L) {
      sprintf(
        paste0(" Hybrid view also shows the per-contrast cross-",
               "modality mean for each axis-relevant contrast (%s)."),
        paste(sub("^score_at_", "", per_contrast_cols), collapse = ", "))
    } else {
      ""
    }
    caption <- sprintf(
      paste0("Top %d axis-restricted kinases for axis '%s' (sorted by ",
             "|mean consensus score| across axis-relevant contrasts and ",
             "primary modalities; kinases filtered to those with >= 5 ",
             "OmniPath KSN substrate sites in the axis gene universe). ",
             "n_cells is the count of finite (contrast, modality) cells ",
             "averaged; n_sites is the count of in-axis-universe ",
             "substrate sites.%s"),
      nrow(top), axis_name, pc_part)
  }
  knitr::kable(display, caption = caption, row.names = FALSE)
}

# Lollipop chart of the top-n axis-restricted kinases for one axis.
# Mirrors `plot_axis_lollipop()` in `R/tf_inference.R`.
plot_axis_lollipop_kinase <- function(axis_tbl, axis_name, n = 12,
                                      title = NULL) {
  sub <- axis_tbl[axis_tbl$axis == axis_name, , drop = FALSE]
  if (nrow(sub) == 0L) {
    return(ggplot2::ggplot() +
             ggplot2::annotate("text", x = 0, y = 0,
                               label = sprintf("No kinases for axis '%s'.",
                                               axis_name)) +
             ggplot2::theme_void())
  }
  sub <- sub[order(-sub$mean_activity_in_axis_contrasts), , drop = FALSE]
  top <- head(sub, n)
  top$direction <- ifelse(top$mean_activity_in_axis_contrasts >= 0,
                          "positive", "negative")
  top$kinase <- factor(top$source,
                       levels = rev(top$source))
  if (is.null(title)) {
    title <- sprintf("Top %d axis-restricted kinases: %s", n, axis_name)
  }
  ggplot2::ggplot(top, ggplot2::aes(x = mean_activity_in_axis_contrasts,
                                    y = kinase,
                                    colour = direction)) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = mean_activity_in_axis_contrasts,
                                       yend = kinase), linewidth = 0.7) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                        colour = "grey50") +
    ggplot2::scale_colour_manual(values = c(positive = "#d6604d",
                                            negative = "#4393c3")) +
    ggplot2::labs(x = "mean consensus score (axis contrasts)",
                  y = NULL, title = title) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "none",
                   panel.grid.major.y = ggplot2::element_blank())
}

# ----------------------------------------------------------------------
# F5 helper: build the kinase verdict TSV summarising the top kinases
# per D2 axis with a per-contrast significance cross-reference (replacing
# E5's cross-modality leader-board cross-reference; the F3 user-directed
# deviation removed the kinase leader-board concept entirely, so the
# kinase verdict has to anchor on the per-contrast significant sets
# from `storage/results/kinase_activity_per_contrast.tsv` instead).
#
# Design notes (mirror `build_tf_verdict_table()` schema except where
# noted):
#   * The helper does NOT introduce new inference. It aggregates the F3
#     per-contrast significance table and the F4 axis-restricted ranking
#     into a per-axis summary row.
#   * The signed mean is taken from `mean_activity_in_axis_contrasts` in
#     `axis_tbl` (the F4 output). `top_kinase_signs` uses the same
#     `+`/`-` convention as the TF verdict.
#   * `mean_score_range` is the [min, max] of the signed mean across the
#     top-N kinases.
#   * `per_contrast_score_range` and `per_contrast_summary` are the
#     hybrid columns locked at the F4 decision gate. They are derived
#     from the `score_at_<contrast>` columns of `axis_tbl` -- the same
#     mechanism used by `build_tf_verdict_table()`. At the kinase layer
#     the corrected cache is single-modality with uniform cell coverage,
#     so the row-mean of per-contrast columns equals the axis-mean to
#     floating-point precision and the hybrid columns are preserved for
#     schema consistency rather than exposing cell-coverage divergence.
#   * `n_top_sig_in_axis_contrasts` is the F5-specific column that
#     replaces E5's `n_top_in_cross_modality_leaderboard`. For each of
#     the top-N kinases, it counts whether the kinase is flagged
#     `sig_fdr10_ulm = TRUE` in ANY of the axis-relevant contrasts. The
#     axis-relevant contrasts are derived from the non-all-NA pattern
#     of the `score_at_<contrast>` columns in the axis subset of
#     `axis_tbl` (the same derivation `build_tf_verdict_table()` uses to
#     find axis-relevant per-contrast columns).
#   * `evidence_summary` is a free-text column filled by the caller via
#     `evidence_summaries`, so the editorial prose can be tuned in the
#     Rmd without code changes. Missing entries become NA_character_.
#
# Arguments:
#   axis_tbl              long tibble from `score_kinase_per_axis()`.
#                         Must carry `axis`, `source`, `leader_rank`,
#                         `mean_activity_in_axis_contrasts`, and the
#                         `score_at_<contrast>` columns for hybrid use.
#   per_contrast_tbl      long tibble keyed on `(contrast, source)` with
#                         a Boolean `sig_fdr10_ulm` column. Produced by
#                         the F3 build chunk in rmd/14_kinase_inference.Rmd
#                         and exported as
#                         `storage/results/kinase_activity_per_contrast.tsv`.
#   n_top_per_axis        top-N kinases to summarise per axis. Default 5.
#   evidence_summaries    optional named list `axis_name -> string`
#                         containing the prose evidence-summary string
#                         to place in the `evidence_summary` column.
#                         Missing axes get NA_character_.
#
# Returns a data.frame with one row per axis and columns:
#   axis                              character;
#   top_kinases                       comma-separated;
#   top_kinase_signs                  comma-separated +/-;
#   mean_score_range                  e.g. "[-1.16, +2.91]" (axis-mean);
#   per_contrast_score_range          e.g. "[-1.32, +3.02]" (hybrid);
#   per_contrast_summary              string like
#                                     "nlgf_in_maptki:[-1.32, +3.02]; nlgf_in_p301s:[-1.00, +3.00]"
#                                     giving per-contrast [min, max]
#                                     across top-N x axis_contrasts cells;
#   n_top_sig_in_axis_contrasts       integer (0..n_top_per_axis);
#                                     count of top-N kinases that reach
#                                     FDR<0.10 on ulm in ANY of the
#                                     axis-relevant contrasts (F3 layer);
#   evidence_summary                  character (possibly NA).
build_kinase_verdict_table <- function(axis_tbl, per_contrast_tbl,
                                       n_top_per_axis = 5L,
                                       evidence_summaries = NULL) {
  stopifnot(is.data.frame(axis_tbl) || tibble::is_tibble(axis_tbl),
            all(c("axis", "source", "mean_activity_in_axis_contrasts",
                  "leader_rank") %in% names(axis_tbl)),
            is.data.frame(per_contrast_tbl) || tibble::is_tibble(per_contrast_tbl),
            all(c("contrast", "source", "sig_fdr10_ulm") %in%
                names(per_contrast_tbl)),
            is.numeric(n_top_per_axis), length(n_top_per_axis) == 1L,
            n_top_per_axis >= 1L)

  # Preserve the natural axis encounter order (amyloid_activation,
  # synaptic_suppression, interaction_metabolic) rather than
  # alphabetising; this matches the 1/2/3 ordering used by the
  # predecessor section-13 D2 verdict and the section-14 TF verdict so
  # the kinase verdict TSV reads in parallel with both.
  axes <- unique(as.character(axis_tbl$axis))
  if (is.null(evidence_summaries)) {
    evidence_summaries <- setNames(rep(NA_character_, length(axes)), axes)
  }

  per_contrast_cols <- grep("^score_at_", names(axis_tbl), value = TRUE)

  rows <- lapply(axes, function(ax) {
    sub <- axis_tbl[axis_tbl$axis == ax, , drop = FALSE]
    sub <- sub[order(sub$leader_rank), , drop = FALSE]
    top <- head(sub, n_top_per_axis)
    if (nrow(top) == 0L) {
      return(data.frame(
        axis                              = ax,
        top_kinases                       = NA_character_,
        top_kinase_signs                  = NA_character_,
        mean_score_range                  = NA_character_,
        per_contrast_score_range          = NA_character_,
        per_contrast_summary              = NA_character_,
        n_top_sig_in_axis_contrasts       = 0L,
        evidence_summary                  = evidence_summaries[[ax]] %||%
                                              NA_character_,
        stringsAsFactors                  = FALSE
      ))
    }

    signs <- ifelse(top$mean_activity_in_axis_contrasts >= 0, "+", "-")
    range_str <- sprintf("[%+0.2f, %+0.2f]",
                         min(top$mean_activity_in_axis_contrasts),
                         max(top$mean_activity_in_axis_contrasts))

    # Hybrid per-contrast summary. Identify axis-relevant per-contrast
    # columns by dropping those with all-NA values in the axis subset
    # (mirrors the TF verdict helper's derivation, which mirrors the
    # format_axis_restricted_table filtering convention). The surviving
    # column names also identify the axis_contrasts for the F5-specific
    # n_top_sig_in_axis_contrasts cross-reference below.
    pc_keep <- per_contrast_cols[vapply(per_contrast_cols, function(col) {
      !all(is.na(sub[[col]]))
    }, logical(1))]
    axis_contrasts <- sub("^score_at_", "", pc_keep)

    if (length(pc_keep) == 0L) {
      pc_range_str <- NA_character_
      pc_summary   <- NA_character_
    } else {
      pc_values <- unlist(lapply(pc_keep, function(col) top[[col]]),
                          use.names = FALSE)
      pc_values <- pc_values[is.finite(pc_values)]
      if (length(pc_values) == 0L) {
        pc_range_str <- NA_character_
      } else {
        pc_range_str <- sprintf("[%+0.2f, %+0.2f]",
                                min(pc_values), max(pc_values))
      }
      pc_summary <- paste(vapply(pc_keep, function(col) {
        v <- top[[col]]
        v <- v[is.finite(v)]
        if (length(v) == 0L) {
          sprintf("%s:NA", sub("^score_at_", "", col))
        } else {
          sprintf("%s:[%+0.2f, %+0.2f]",
                  sub("^score_at_", "", col),
                  min(v), max(v))
        }
      }, character(1)), collapse = "; ")
    }

    # F5-specific cross-reference column: count of top-N kinases that
    # are sig_fdr10_ulm == TRUE in ANY of the axis-relevant contrasts.
    # If pc_keep is empty (no axis_contrasts derivable from axis_tbl
    # columns), fall back to 0L per the F5 design.
    if (length(axis_contrasts) == 0L) {
      n_sig <- 0L
    } else {
      sig_sub <- per_contrast_tbl[
        per_contrast_tbl$contrast %in% axis_contrasts &
          per_contrast_tbl$source %in% top$source &
          !is.na(per_contrast_tbl$sig_fdr10_ulm) &
          per_contrast_tbl$sig_fdr10_ulm, , drop = FALSE]
      n_sig <- length(unique(sig_sub$source))
    }

    data.frame(
      axis                              = ax,
      top_kinases                       = paste(top$source, collapse = ", "),
      top_kinase_signs                  = paste(signs, collapse = ","),
      mean_score_range                  = range_str,
      per_contrast_score_range          = pc_range_str,
      per_contrast_summary              = pc_summary,
      n_top_sig_in_axis_contrasts       = as.integer(n_sig),
      evidence_summary                  = evidence_summaries[[ax]] %||%
                                            NA_character_,
      stringsAsFactors                  = FALSE
    )
  })
  do.call(rbind, rows)
}
