# Agnostic cross-modality pathway-ranking helpers used by the Phase C
# agnostic pathway survey (`rmd/12_pathway_survey.Rmd`). The ranking is
# deliberately metric-driven and sign-symmetric: pathways are ordered by
# how many modalities call them significant with a consistent direction
# of effect, ties broken by mean |NES| across modalities. No pathway
# family receives privileged framing -- this module exists to test
# whether the previously-promoted OXPHOS signal is genuinely dominant or
# one among many. See `storage/notes/pathway_overhaul_plan.md` Phase C
# and the anti-anchoring guardrails at the bottom of that file.
#
# Cache shape consumed throughout: `modality -> contrast -> fgseaResult`
# (a data.table with columns pathway / pval / padj / log2err / ES / NES /
# size / leadingEdge). This is the layout produced by
# `run_fgsea_per_dataset()` in `R/fgsea.R` and by every Phase B build
# script (B1 GO MF + GO CC, B2 custom microglia states, B3 custom
# microglia AD, B4 custom module sources).

# Compute a tidy cross-modality ranking for one contrast from one fgsea
# cache.
#
# Arguments:
#   fgsea_list  named list `modality -> contrast -> fgseaResult`.
#   contrast    contrast name (one of names(fgsea_list[[1]])); the five
#               canonical project contrasts are nlgf_in_maptki,
#               nlgf_in_p301s, interaction, tau_alone, tau_in_nlgf.
#   modalities  character vector of modality names to include; defaults
#               to every name of `fgsea_list`. Used to restrict the
#               ranking to a subset of modalities (e.g. RNA-only).
#   padj_cut   significance threshold for the per-modality call
#              (default 0.05).
#
# Returns a tibble with one row per pathway and the following columns:
#   pathway                            character; pathway name as
#                                      stored in the fgsea cache.
#   n_modalities_sig                   integer; count of modalities
#                                      with padj < padj_cut (any sign).
#   n_modalities_sig_consistent_sign   integer; among the significant
#                                      modalities, the count belonging
#                                      to the larger sign group (0 if
#                                      no significant modalities).
#                                      With (NES = +2.1, +1.8, -1.5)
#                                      and all three significant this
#                                      is 2.
#   mean_abs_nes                       numeric; mean |NES| across
#                                      modalities with finite NES,
#                                      irrespective of significance.
#                                      NA only if every modality is
#                                      missing or NA.
#   sign_consensus                     "+" / "-" / "mixed" / NA. "+"
#                                      if every significant modality
#                                      has positive NES; "-" if every
#                                      significant modality has
#                                      negative NES; "mixed" if signs
#                                      differ; NA if no significant
#                                      modalities. This is the strict
#                                      consensus -- a single
#                                      disagreeing significant call
#                                      flips it to "mixed", even
#                                      though n_modalities_sig_consistent_sign
#                                      still records the larger group.
#   nes_<modality>                     numeric; per-modality NES (NA
#                                      if the pathway was missing from
#                                      that modality's fgsea, e.g.
#                                      below min_size after
#                                      intersection with the modality
#                                      gene universe).
#   padj_<modality>                    numeric; per-modality padj (NA
#                                      if missing).
#   composite_rank                     integer; 1-indexed rank from
#                                      sorting by
#                                      (n_modalities_sig_consistent_sign desc,
#                                       mean_abs_nes desc). Ties broken
#                                      stably by tibble row order.
#
# Pathways present in any modality appear; pathways absent from all
# modalities cannot appear by construction (their cache rows do not
# exist). NAs in nes_<m> / padj_<m> survive faithfully so downstream
# heatmaps can distinguish "missing" from "not significant".
rank_pathways_cross_modality <- function(fgsea_list, contrast,
                                         modalities = NULL,
                                         padj_cut = 0.05) {
  if (is.null(modalities)) modalities <- names(fgsea_list)
  stopifnot(length(modalities) > 0L,
            all(modalities %in% names(fgsea_list)),
            length(contrast) == 1L)

  per_modality <- lapply(modalities, function(mn) {
    tbl <- fgsea_list[[mn]][[contrast]]
    if (is.null(tbl) || nrow(tbl) == 0L) {
      data.frame(pathway = character(0), NES = numeric(0),
                 padj = numeric(0), stringsAsFactors = FALSE)
    } else {
      as.data.frame(tbl)[, c("pathway", "NES", "padj")]
    }
  })
  names(per_modality) <- modalities

  all_pathways <- sort(unique(unlist(lapply(per_modality, `[[`, "pathway"),
                                     use.names = FALSE)))
  if (length(all_pathways) == 0L) {
    empty <- tibble::tibble(
      pathway = character(0),
      n_modalities_sig = integer(0),
      n_modalities_sig_consistent_sign = integer(0),
      mean_abs_nes = numeric(0),
      sign_consensus = character(0)
    )
    for (mn in modalities) {
      empty[[paste0("nes_", mn)]]  <- numeric(0)
      empty[[paste0("padj_", mn)]] <- numeric(0)
    }
    empty$composite_rank <- integer(0)
    return(empty)
  }

  nes_mat  <- matrix(NA_real_, nrow = length(all_pathways),
                     ncol = length(modalities),
                     dimnames = list(all_pathways, modalities))
  padj_mat <- matrix(NA_real_, nrow = length(all_pathways),
                     ncol = length(modalities),
                     dimnames = list(all_pathways, modalities))
  for (mn in modalities) {
    tbl <- per_modality[[mn]]
    if (nrow(tbl) == 0L) next
    idx <- match(tbl$pathway, all_pathways)
    nes_mat[idx, mn]  <- tbl$NES
    padj_mat[idx, mn] <- tbl$padj
  }

  sig_mat <- !is.na(padj_mat) & padj_mat < padj_cut
  n_sig <- rowSums(sig_mat)

  # Per-pathway sign accounting over the SIGNIFICANT subset only. The
  # consistent-sign count is the larger of (positive-sig, negative-sig);
  # the consensus label is strict (any disagreement -> "mixed") so the
  # two columns answer two distinct questions: "how many agree on the
  # dominant direction" vs "do all significant modalities agree on
  # direction".
  pos_sig <- vapply(seq_along(all_pathways), function(i) {
    if (n_sig[i] == 0L) return(0L)
    sum(nes_mat[i, sig_mat[i, ]] > 0, na.rm = TRUE)
  }, integer(1))
  neg_sig <- vapply(seq_along(all_pathways), function(i) {
    if (n_sig[i] == 0L) return(0L)
    sum(nes_mat[i, sig_mat[i, ]] < 0, na.rm = TRUE)
  }, integer(1))

  n_consistent <- pmax(pos_sig, neg_sig)
  sign_consensus <- ifelse(
    n_sig == 0L, NA_character_,
    ifelse(pos_sig > 0L & neg_sig == 0L, "+",
           ifelse(neg_sig > 0L & pos_sig == 0L, "-", "mixed"))
  )

  mean_abs_nes <- apply(abs(nes_mat), 1, function(r) {
    r <- r[is.finite(r)]
    if (length(r) == 0L) NA_real_ else mean(r)
  })

  out <- tibble::tibble(
    pathway                          = all_pathways,
    n_modalities_sig                 = as.integer(n_sig),
    n_modalities_sig_consistent_sign = as.integer(n_consistent),
    mean_abs_nes                     = mean_abs_nes,
    sign_consensus                   = sign_consensus
  )
  for (mn in modalities) {
    out[[paste0("nes_", mn)]]  <- nes_mat[, mn]
    out[[paste0("padj_", mn)]] <- padj_mat[, mn]
  }

  out |>
    dplyr::arrange(dplyr::desc(n_modalities_sig_consistent_sign),
                   dplyr::desc(mean_abs_nes)) |>
    dplyr::mutate(composite_rank = seq_len(dplyr::n()))
}

# Format the top-n rows of a cross-modality ranking as a knitr::kable
# for inline display in the Rmd. The full ranking (every pathway, every
# contrast) is exported to TSV separately; this function only handles
# the visible table.
#
# Arguments:
#   ranking_tbl       output of `rank_pathways_cross_modality()`,
#                     optionally augmented with extra annotation
#                     columns the caller wants to display
#                     (e.g. source_collection for pooled rankings).
#   contrast          contrast string, used in the default caption.
#   n                 top-n rows to display (default 20). If
#                     `nrow(ranking_tbl) <= n` every row is shown.
#   strip_prefix      optional regex passed to `sub()` to strip a
#                     leading collection prefix from pathway names;
#                     e.g. "^GOBP_" for the GO BP collection. Leave
#                     NULL for custom collections whose pathway names
#                     are already human-readable.
#   caption           overrides the default kable caption.
#   include_padj      if TRUE, include per-modality padj columns.
#                     Defaults to FALSE to keep the table narrow; full
#                     padj data is in the exported TSV.
#   extra_cols        optional character vector of column names from
#                     `ranking_tbl` to insert immediately after
#                     `composite_rank` and before `pathway`. Used for
#                     pooled rankings that carry a source_collection
#                     annotation (Phase C4) or a module_id annotation
#                     (Phase C6). Columns must already exist on
#                     `ranking_tbl`.
#
# Returns a knitr::kable object suitable for `print()` inside an
# `results = 'asis'` chunk.
format_ranking_table <- function(ranking_tbl, contrast, n = 20,
                                 strip_prefix = NULL,
                                 caption = NULL,
                                 include_padj = FALSE,
                                 extra_cols = NULL) {
  top <- head(ranking_tbl, n)
  if (!is.null(strip_prefix)) {
    top$pathway <- sub(strip_prefix, "", top$pathway)
  }
  top$pathway <- gsub("_", " ", tolower(top$pathway))

  if (!is.null(extra_cols)) {
    missing <- setdiff(extra_cols, names(top))
    if (length(missing) > 0L) {
      stop(sprintf("format_ranking_table: extra_cols not in ranking_tbl: %s",
                   paste(missing, collapse = ", ")), call. = FALSE)
    }
  }

  nes_cols  <- grep("^nes_",  names(top), value = TRUE)
  padj_cols <- grep("^padj_", names(top), value = TRUE)
  display_cols <- c("composite_rank")
  if (!is.null(extra_cols)) display_cols <- c(display_cols, extra_cols)
  display_cols <- c(display_cols, "pathway",
                    "n_modalities_sig",
                    "n_modalities_sig_consistent_sign",
                    "sign_consensus",
                    "mean_abs_nes",
                    nes_cols)
  if (isTRUE(include_padj)) display_cols <- c(display_cols, padj_cols)

  if (is.null(caption)) {
    caption <- sprintf(
      paste0("Top %d cross-modality pathways for contrast '%s' ",
             "(sorted by sign-consistent significant-modality count, ",
             "then mean |NES|; padj cutoff 0.05)."),
      min(n, nrow(top)), contrast)
  }

  knitr::kable(top[, display_cols], digits = 3, caption = caption)
}

# Build a long tidy ranking from the per-substate fgsea cache.
#
# The per-substate cache produced in Phase B5
# (`storage/cache/fgsea_per_state_results.rds`) has a different shape from
# the whole-microglia caches used by `rank_pathways_cross_modality()`:
# `collection -> substate -> contrast -> fgseaResult`. Each substate's
# pathway universe differs because the snRNAseq nebula run was restricted
# to that substate's microglia, so per-substate top-n lists are the
# natural ranking unit. The substate_breadth column then quantifies how
# concentrated each pathway's signal is across substates -- a pathway
# that lands in the per-substate top-n in only one substate is a
# substate-specific signal, one that lands in three or four is broadly
# active across the microglia compartment.
#
# Arguments:
#   fgsea_substate_cache  named list `collection -> substate -> contrast
#                         -> fgseaResult` as written by the B5 rebuild of
#                         `rmd/02e_snrnaseq_substate_pathway.Rmd`.
#   n_top                 per-substate top-n cutoff used to derive the
#                         `in_top_n` boolean and the `substate_breadth`
#                         count (default 10). The value flows through to
#                         the output column name only via documentation;
#                         the column itself is always called `in_top_n`
#                         and callers can rename it after the fact.
#
# Returns a tibble with one row per (collection, contrast, substate,
# pathway) carrying:
#   collection            collection key (gobp / gomf / gocc /
#                         custom_microglia_states / custom_microglia_ad /
#                         custom_module_sources).
#   contrast              one of nlgf_in_maptki, nlgf_in_p301s,
#                         interaction, tau_alone, tau_in_nlgf.
#   substate              one of homeostatic, DAM, IFN, proliferative.
#   pathway               pathway name as stored in the fgsea cache.
#   NES                   normalised enrichment score for this
#                         (collection, contrast, substate, pathway).
#   padj                  fgsea-adjusted p-value (BH).
#   abs_nes               |NES|; the sort key for the per-substate rank.
#   substate_rank         1-indexed rank within
#                         (collection, contrast, substate) by abs_nes
#                         descending; ties broken stably by input order.
#   in_top_n              TRUE iff substate_rank <= n_top.
#   substate_breadth      integer count of substates where the same
#                         pathway lands in_top_n at the SAME
#                         (collection, contrast). Range 0..4 since the
#                         project has four substates. Pathways with
#                         breadth = 0 cannot appear in any per-substate
#                         top-n; pathways with breadth = 4 are
#                         compartment-wide.
rank_pathways_per_substate <- function(fgsea_substate_cache, n_top = 10L) {
  stopifnot(is.list(fgsea_substate_cache),
            length(fgsea_substate_cache) > 0L,
            n_top >= 1L)
  collections <- names(fgsea_substate_cache)
  rows <- list()
  for (cn in collections) {
    substates <- names(fgsea_substate_cache[[cn]])
    if (length(substates) == 0L) next
    # Contrasts assumed homogeneous across substates within a collection
    # (the B5 builder guarantees this); pick from the first substate.
    contrasts <- names(fgsea_substate_cache[[cn]][[substates[1L]]])
    for (ctr in contrasts) {
      for (ss in substates) {
        tbl <- fgsea_substate_cache[[cn]][[ss]][[ctr]]
        if (is.null(tbl) || nrow(tbl) == 0L) next
        df <- as.data.frame(tbl)[, c("pathway", "NES", "padj"), drop = FALSE]
        df$collection <- cn
        df$contrast   <- ctr
        df$substate   <- ss
        df$abs_nes    <- abs(df$NES)
        df <- df[order(-df$abs_nes), , drop = FALSE]
        df$substate_rank <- seq_len(nrow(df))
        df$in_top_n      <- df$substate_rank <= n_top
        rows[[length(rows) + 1L]] <- df
      }
    }
  }
  if (length(rows) == 0L) {
    return(tibble::tibble(
      collection = character(0), contrast = character(0),
      substate = character(0), pathway = character(0),
      NES = numeric(0), padj = numeric(0), abs_nes = numeric(0),
      substate_rank = integer(0), in_top_n = logical(0),
      substate_breadth = integer(0)
    ))
  }
  long <- do.call(rbind, rows)
  long <- long[, c("collection", "contrast", "substate", "pathway",
                   "NES", "padj", "abs_nes", "substate_rank", "in_top_n"),
               drop = FALSE]
  long <- tibble::as_tibble(long)
  long |>
    dplyr::group_by(.data$collection, .data$contrast, .data$pathway) |>
    dplyr::mutate(substate_breadth = as.integer(sum(.data$in_top_n))) |>
    dplyr::ungroup()
}

# Format the per-substate top-n table for one (collection, contrast,
# substate) cell. Pathways with substate_breadth >= breadth_mark are
# prefixed with a "*" marker so cross-substate consistency is visible
# inline. The display drops abs_nes (redundant given NES) and the
# in_top_n boolean (every row shown is in the top-n by construction).
#
# Arguments:
#   long_tbl         output of `rank_pathways_per_substate()`.
#   collection_key   collection name to filter on.
#   contrast_key     contrast name to filter on.
#   substate_key     substate name to filter on.
#   n_top            top-n cutoff for the table (default 10). Capped at
#                    the number of rows actually present after filtering.
#   breadth_mark     `substate_breadth >= breadth_mark` triggers the
#                    "*" prefix in the pathway display column (default 3).
#   strip_prefix     optional regex passed to `sub()` to strip a leading
#                    collection prefix from pathway names (e.g. "^GOBP_").
#   caption          overrides the default kable caption.
#
# Base R filter/order to sidestep the dplyr NSE name-clash that would
# arise if the function arguments shared names with column names.
format_per_substate_table <- function(long_tbl, collection_key, contrast_key,
                                       substate_key, n_top = 10L,
                                       breadth_mark = 3L,
                                       strip_prefix = NULL, caption = NULL) {
  hits <- long_tbl$collection == collection_key &
          long_tbl$contrast   == contrast_key   &
          long_tbl$substate   == substate_key   &
          long_tbl$in_top_n
  hits[is.na(hits)] <- FALSE
  top <- long_tbl[hits, , drop = FALSE]
  if (nrow(top) == 0L) {
    if (is.null(caption)) {
      caption <- sprintf("No pathways for %s / %s / %s.",
                         collection_key, contrast_key, substate_key)
    }
    return(knitr::kable(
      data.frame(message = "No pathways in this cell."),
      caption = caption, row.names = FALSE))
  }
  top <- top[order(top$substate_rank), , drop = FALSE]
  if (nrow(top) > n_top) top <- top[seq_len(n_top), , drop = FALSE]

  pathway_display <- top$pathway
  if (!is.null(strip_prefix)) {
    pathway_display <- sub(strip_prefix, "", pathway_display)
  }
  pathway_display <- gsub("_", " ", tolower(pathway_display))
  pathway_display <- ifelse(top$substate_breadth >= breadth_mark,
                            paste0("* ", pathway_display),
                            pathway_display)

  out <- data.frame(
    rank             = top$substate_rank,
    pathway          = pathway_display,
    NES              = round(top$NES, 3),
    padj             = signif(top$padj, 3),
    substate_breadth = top$substate_breadth,
    stringsAsFactors = FALSE
  )
  if (is.null(caption)) {
    caption <- sprintf(
      paste0("Top %d pathways by |NES| for %s substate at contrast '%s' ",
             "(* prefix marks pathways with substate_breadth >= %d, i.e. ",
             "appearing in top-%d of at least %d substates)."),
      nrow(out), substate_key, contrast_key, breadth_mark, n_top,
      breadth_mark)
  }
  knitr::kable(out, caption = caption, row.names = FALSE)
}

# Build the unified unbiased leader board from the per-collection
# cross-modality rankings produced in Phase C1-C6. This is the D1 step
# of `storage/notes/pathway_overhaul_plan.md`: collapse the five
# per-collection rankings (gobp, gomf, gocc, custom-pooled, hdwgcna)
# into one row per (collection, pathway), keeping only pathways that
# qualify as "leaders" at at least one contrast under the supplied
# `leader_rule`. The default rule fires when EITHER >=2 modalities
# agree on direction at any single contrast (cross-modality consistent-
# sign breadth) OR >=3 modalities reach significance regardless of sign
# (cross-modality breadth, sign-agnostic; lets mixed-sign hits like
# OXPHOS-at-interaction surface as informative rather than be filtered
# out).
#
# Aggregation columns (one row per `collection`, `pathway`):
#   n_contrasts_leader              integer 1..5; how many contrasts
#                                   pass the leader rule for this
#                                   (collection, pathway).
#   n_contrasts_consistent_sign_ge2 integer 0..5; how many contrasts
#                                   have n_modalities_sig_consistent_sign
#                                   >= 2 (the cross-contrast consistent-
#                                   sign breadth signal).
#   n_contrasts_sig_ge3             integer 0..5; how many contrasts
#                                   have n_modalities_sig >= 3
#                                   (cross-modality breadth at the
#                                   contrast, sign-agnostic).
#   max_consistent_sign             integer; max
#                                   n_modalities_sig_consistent_sign
#                                   across the five contrasts.
#   max_n_modalities_sig            integer; max n_modalities_sig
#                                   across the five contrasts.
#   max_abs_nes                     numeric; max mean_abs_nes across
#                                   the five contrasts (peak effect-
#                                   size magnitude at any contrast).
#   dominant_sign                   "+" / "-" / "mixed". Based on the
#                                   sign_consensus values at LEADER
#                                   contrasts only. "+" if every
#                                   leader contrast has +, "-" if
#                                   every leader contrast has -,
#                                   "mixed" otherwise (including any
#                                   single "mixed" leader). NA only
#                                   if no leader contrasts have a
#                                   computable sign_consensus (which
#                                   cannot happen given the rule).
#   max_substate_breadth            integer 0..4; max substate_breadth
#                                   from `per_state_long` across all
#                                   contrasts for this (collection,
#                                   pathway). 0 when the pathway is
#                                   not in the per-substate cache
#                                   (e.g. hdwgcna modules, which were
#                                   not included in the B5 per-state
#                                   fgsea rebuild).
#   contrasts_summary               character; pipe-delimited string
#                                   listing each LEADER contrast as
#                                   `contrast:n_modalities_sig/sign_consensus`,
#                                   e.g. "interaction:3/mixed |
#                                   tau_in_nlgf:2/-". NA if no leader
#                                   contrasts.
#   leader_score                    numeric composite for primary
#                                   sort:
#                                   `5 * n_contrasts_consistent_sign_ge2
#                                    + n_contrasts_sig_ge3
#                                    + max_consistent_sign / 5`.
#                                   Range ~0 to 31. Designed so
#                                   cross-contrast consistent-sign
#                                   support (worth 5 each, up to 25)
#                                   dominates cross-modality breadth
#                                   at any single contrast (worth 1
#                                   each, up to 5); the third term
#                                   (0 to 1) is a tie-breaker for
#                                   single-contrast peak intensity.
#
# Arguments:
#   rankings_named_list  named list of cross-modality ranking tibbles
#                        (output of `rank_pathways_cross_modality()` for
#                        each contrast, rbound by contrast and exported
#                        to per-collection TSV). The expected names
#                        match the Phase C survey: gobp, gomf, gocc,
#                        custom, hdwgcna. Each tibble must carry
#                        `pathway`, `contrast`,
#                        `n_modalities_sig_consistent_sign`,
#                        `n_modalities_sig`, `sign_consensus`,
#                        `mean_abs_nes`. If a tibble has a
#                        `source_collection` column (the custom-pool
#                        TSV from C4) its values are used as the
#                        collection label for each row, so the unified
#                        output splits the pool back into its
#                        components. Otherwise the list-element name
#                        is used as the collection label for every row.
#                        Callers should pre-remap source_collection
#                        values to match `per_state_long$collection`
#                        if they want the substate-breadth join to fire
#                        for custom-pool entries.
#   per_state_long       optional output of
#                        `rank_pathways_per_substate()` (or its TSV
#                        re-read). Must carry `collection`, `pathway`,
#                        and `substate_breadth`. If NULL, every
#                        `max_substate_breadth` falls back to 0.
#   leader_rule          optional vectorised closure
#                        `function(tbl) -> logical(nrow(tbl))` deciding
#                        whether each (collection, pathway, contrast)
#                        row qualifies as a leader. Default: row's
#                        `n_modalities_sig_consistent_sign >= 2` OR
#                        `n_modalities_sig >= 3`. The closure receives
#                        the full stacked tibble (all collections, all
#                        contrasts) so vectorised operations are the
#                        most efficient pattern; row-by-row evaluation
#                        is supported by accident of R's recycling
#                        rules but is not recommended.
#
# Returns a tibble of leader rows sorted by `leader_score` desc, ties
# broken by `n_contrasts_consistent_sign_ge2` desc, `max_consistent_sign`
# desc, `max_abs_nes` desc, then collection / pathway alphabetical.
build_leader_board <- function(rankings_named_list,
                                per_state_long = NULL,
                                leader_rule    = NULL) {
  stopifnot(is.list(rankings_named_list),
            length(rankings_named_list) > 0L,
            !is.null(names(rankings_named_list)),
            all(nzchar(names(rankings_named_list))))

  if (is.null(leader_rule)) {
    leader_rule <- function(row) {
      row$n_modalities_sig_consistent_sign >= 2L |
      row$n_modalities_sig                  >= 3L
    }
  }

  required_cols <- c("pathway", "contrast", "n_modalities_sig",
                     "n_modalities_sig_consistent_sign",
                     "sign_consensus", "mean_abs_nes")

  # Stack the per-collection inputs into one long tibble with a
  # `collection` label per row. If an input tibble carries a
  # source_collection column (the custom-pool case), use those values
  # as the collection label so the unified output splits the pool back
  # into its three sub-collections. Otherwise the input list name
  # becomes the collection label.
  stacked_list <- lapply(seq_along(rankings_named_list), function(i) {
    cn  <- names(rankings_named_list)[i]
    tbl <- rankings_named_list[[i]]
    missing_cols <- setdiff(required_cols, names(tbl))
    if (length(missing_cols) > 0L) {
      stop(sprintf(
        "build_leader_board: input '%s' missing required columns: %s",
        cn, paste(missing_cols, collapse = ", ")), call. = FALSE)
    }
    if (nrow(tbl) == 0L) return(NULL)
    if ("source_collection" %in% names(tbl)) {
      tbl$collection <- as.character(tbl$source_collection)
    } else {
      tbl$collection <- cn
    }
    tbl
  })
  stacked <- dplyr::bind_rows(stacked_list)
  if (nrow(stacked) == 0L) {
    return(tibble::tibble(
      collection                      = character(0),
      pathway                         = character(0),
      n_contrasts_leader              = integer(0),
      n_contrasts_consistent_sign_ge2 = integer(0),
      n_contrasts_sig_ge3             = integer(0),
      max_consistent_sign             = integer(0),
      max_n_modalities_sig            = integer(0),
      max_abs_nes                     = numeric(0),
      dominant_sign                   = character(0),
      max_substate_breadth            = integer(0),
      contrasts_summary               = character(0),
      leader_score                    = numeric(0)
    ))
  }

  is_leader <- as.logical(leader_rule(stacked))
  if (length(is_leader) != nrow(stacked)) {
    stop("build_leader_board: leader_rule must return a logical vector ",
         "of length nrow(stacked).", call. = FALSE)
  }
  is_leader[is.na(is_leader)] <- FALSE
  stacked$is_leader <- is_leader

  # Per-group reductions for the categorical aggregates. Defined
  # outside dplyr::summarise() to keep that call readable.
  reduce_dominant_sign <- function(sign_consensus, is_leader) {
    sc <- sign_consensus[is_leader]
    sc <- sc[!is.na(sc)]
    if (length(sc) == 0L) return(NA_character_)
    if (all(sc == "+")) return("+")
    if (all(sc == "-")) return("-")
    "mixed"
  }
  reduce_contrasts_summary <- function(contrast, n_modalities_sig,
                                       sign_consensus, is_leader) {
    ix <- which(is_leader)
    if (length(ix) == 0L) return(NA_character_)
    paste(sprintf("%s:%d/%s",
                  contrast[ix],
                  n_modalities_sig[ix],
                  ifelse(is.na(sign_consensus[ix]), "NA",
                         sign_consensus[ix])),
          collapse = " | ")
  }

  per_group <- stacked |>
    dplyr::group_by(collection, pathway) |>
    dplyr::summarise(
      n_contrasts_leader              = as.integer(sum(is_leader, na.rm = TRUE)),
      n_contrasts_consistent_sign_ge2 = as.integer(sum(
        n_modalities_sig_consistent_sign >= 2L, na.rm = TRUE)),
      n_contrasts_sig_ge3             = as.integer(sum(
        n_modalities_sig >= 3L, na.rm = TRUE)),
      max_consistent_sign             = as.integer(suppressWarnings(
        max(n_modalities_sig_consistent_sign, na.rm = TRUE))),
      max_n_modalities_sig            = as.integer(suppressWarnings(
        max(n_modalities_sig, na.rm = TRUE))),
      max_abs_nes                     = suppressWarnings(
        max(mean_abs_nes, na.rm = TRUE)),
      dominant_sign                   = reduce_dominant_sign(
        sign_consensus, is_leader),
      contrasts_summary               = reduce_contrasts_summary(
        contrast, n_modalities_sig, sign_consensus, is_leader),
      .groups = "drop"
    ) |>
    dplyr::filter(n_contrasts_leader >= 1L)

  # Sanitise -Inf produced by max() over an all-NA group (cannot occur
  # given the >=1 leader filter, but be defensive).
  is_neg_inf <- function(x) is.numeric(x) & is.finite(x) == FALSE & x < 0
  per_group$max_consistent_sign[is_neg_inf(per_group$max_consistent_sign)]   <- 0L
  per_group$max_n_modalities_sig[is_neg_inf(per_group$max_n_modalities_sig)] <- 0L
  per_group$max_abs_nes[is_neg_inf(per_group$max_abs_nes)]                   <- NA_real_

  # Join per-state breadth. Pathways absent from per_state_long (e.g.
  # hdwgcna modules) get max_substate_breadth = 0.
  if (!is.null(per_state_long) && nrow(per_state_long) > 0L) {
    ps_required <- c("collection", "pathway", "substate_breadth")
    ps_missing  <- setdiff(ps_required, names(per_state_long))
    if (length(ps_missing) > 0L) {
      stop(sprintf(
        "build_leader_board: per_state_long missing columns: %s",
        paste(ps_missing, collapse = ", ")), call. = FALSE)
    }
    ps_max <- per_state_long |>
      dplyr::group_by(collection, pathway) |>
      dplyr::summarise(
        max_substate_breadth = as.integer(max(substate_breadth, na.rm = TRUE)),
        .groups = "drop"
      )
    per_group <- per_group |>
      dplyr::left_join(ps_max, by = c("collection", "pathway"))
  } else {
    per_group$max_substate_breadth <- NA_integer_
  }
  per_group$max_substate_breadth[is.na(per_group$max_substate_breadth)] <- 0L

  per_group$leader_score <-
    5 * per_group$n_contrasts_consistent_sign_ge2 +
    per_group$n_contrasts_sig_ge3 +
    per_group$max_consistent_sign / 5

  per_group |>
    dplyr::arrange(dplyr::desc(leader_score),
                   dplyr::desc(n_contrasts_consistent_sign_ge2),
                   dplyr::desc(max_consistent_sign),
                   dplyr::desc(max_abs_nes),
                   collection, pathway) |>
    dplyr::select(collection, pathway,
                  n_contrasts_leader,
                  n_contrasts_consistent_sign_ge2,
                  n_contrasts_sig_ge3,
                  max_consistent_sign,
                  max_n_modalities_sig,
                  max_abs_nes,
                  dominant_sign,
                  max_substate_breadth,
                  contrasts_summary,
                  leader_score)
}
