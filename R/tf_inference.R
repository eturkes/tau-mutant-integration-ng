# Transcription factor activity inference helpers used by the Phase E
# mechanism layer (see storage/notes/mechanism_layer_plan.md, sessions
# E2..E5). Wraps decoupleR's univariate methods (ULM + WSUM family) and
# the CollecTRI mouse prior produced by `decoupleR::get_collectri()` /
# OmnipathR. The design mirrors `R/pathway_survey.R`:
#   - cross-modality, sign-aware, family-agnostic;
#   - never pre-commits to a single TF family;
#   - returns long tidy tibbles with a single canonical key set so
#     downstream aggregation is uniform across modalities.
#
# Cache shape consumed: any project DE cache exposing per-contrast
# `$top` tables with at least `symbol` and `t` columns. The three
# supported shapes are:
#   1. nebula whole-microglia:    `cache$top[[contrast]]`
#   2. nebula per-substate:       `cache[[substate]]$top[[contrast]]`
#   3. limma (geomx, proteomics): `cache$fit$top[[contrast]]`
# The five canonical contrasts are nlgf_in_maptki, nlgf_in_p301s,
# interaction, tau_alone, tau_in_nlgf. Phospho is deliberately deferred
# to Phase F because phospho-site stats need protein-level collapse +
# signed re-mapping before they can serve as TF-target evidence.
#
# Cache shape produced: `modality -> contrast -> tibble(statistic,
# source, condition, score, p_value)`. `condition == contrast` by
# construction; the column is retained for round-trip compatibility
# with the raw decoupleR output. `statistic` carries the method label
# (`ulm`, `wsum`, `norm_wsum`, `corr_wsum`, `consensus`). Per-method
# p-values: `ulm` is the linear-model F-test p; `wsum`, `norm_wsum`,
# `corr_wsum` are empirical permutation p-values floored at 1/1000
# (so the minimum visible p_value is 0.02 in practice); `consensus`
# is the cross-method Stouffer-style combined p-value (continuous,
# can dip below the wsum floor). The E2 stub in this file previously
# claimed consensus p_value was NA -- that was a mis-read of an
# earlier decoupleR version; the current 2.16.0 output populates it.
# The corrected behaviour is verified empirically in the E3 smoke test.
#
# Rationale for the duplicate-symbol collapse rule: a single gene
# symbol can map to multiple feature rows (e.g. several Ensembl IDs in
# the snRNAseq nebula output, several ProteinGroups in the proteomics
# limma fit). Choosing the row with the largest |stat| preserves the
# strongest signed evidence per gene. Alternatives considered: mean
# (dilutes signal), sum (inflates), first (arbitrary). decoupleR
# requires a unique gene-symbol axis, so the collapse must happen
# before the call.

# Extract a single per-contrast stat vector keyed by gene symbol.
#
# Arguments:
#   top_df    a data.frame as found in any of the three supported
#             cache shapes' `$top[[contrast]]` slot. Must contain at
#             least `id_col` (default "symbol") and `stat_col`
#             (default "t").
#   stat_col  column name of the test statistic to use as the rank
#             value. Default "t": the post-hoc Wald-style t = logFC/se
#             for nebula and the limma moderated t for geomx /
#             proteomics. Both are signed and on a comparable scale,
#             so a single helper handles all three cache shapes.
#   id_col    column name carrying the gene symbol. Default "symbol".
#             CollecTRI mouse is symbol-indexed, so this is the
#             canonical key.
#
# Returns a named numeric vector of stats, names = gene symbol.
# Duplicates collapsed by keeping max |stat|; NA stats / NA / empty
# symbols dropped. Empty input returns a length-0 numeric vector.
extract_de_stat_vec <- function(top_df,
                                stat_col = "t",
                                id_col   = "symbol") {
  stopifnot(is.data.frame(top_df),
            id_col   %in% names(top_df),
            stat_col %in% names(top_df))

  sym  <- as.character(top_df[[id_col]])
  stat <- as.numeric(top_df[[stat_col]])

  keep <- !is.na(sym) & nzchar(sym) & !is.na(stat) & is.finite(stat)
  sym  <- sym[keep]
  stat <- stat[keep]
  if (length(stat) == 0L) {
    return(setNames(numeric(0), character(0)))
  }

  if (anyDuplicated(sym) > 0L) {
    # Collapse duplicate symbols by keeping the row with the largest
    # |stat|. ord places the largest |stat| within each symbol group
    # first; !duplicated then picks that first occurrence.
    ord <- order(sym, -abs(stat))
    sym  <- sym[ord]
    stat <- stat[ord]
    keep <- !duplicated(sym)
    sym  <- sym[keep]
    stat <- stat[keep]
  }
  setNames(stat, sym)
}

# Build a numeric matrix of stats with rows = gene symbol and columns =
# contrast, suitable for direct decoupleR consumption.
#
# Arguments:
#   top_list  named list `contrast -> data.frame`, as in any of the
#             three supported cache shapes.
#   stat_col, id_col  forwarded to `extract_de_stat_vec()`.
#
# Returns a numeric matrix; rownames are the union of symbols across
# contrasts; missing (symbol, contrast) cells are NA. Empty input
# (no contrasts) returns a 0x0 matrix.
extract_de_stat_matrix <- function(top_list,
                                   stat_col = "t",
                                   id_col   = "symbol") {
  stopifnot(is.list(top_list))
  if (length(top_list) == 0L) {
    return(matrix(numeric(0), nrow = 0, ncol = 0))
  }
  per_contrast <- lapply(top_list, extract_de_stat_vec,
                         stat_col = stat_col, id_col = id_col)
  all_syms <- sort(unique(unlist(lapply(per_contrast, names),
                                 use.names = FALSE)))
  contrasts <- names(top_list)
  mat <- matrix(NA_real_, nrow = length(all_syms), ncol = length(contrasts),
                dimnames = list(all_syms, contrasts))
  for (cn in contrasts) {
    v <- per_contrast[[cn]]
    if (length(v) == 0L) next
    mat[names(v), cn] <- unname(v)
  }
  mat
}

# Run decoupleR univariate methods on one (modality) stat matrix.
#
# Default methods are `ulm` (univariate linear model, the fast workhorse)
# and `wsum` (weighted sum; emits three sub-statistics raw / norm /
# corr). Both are signed and consume the CollecTRI `mor` mode-of-
# regulation column natively. With `consensus = TRUE`, decoupleR
# appends a `consensus` row per (TF, contrast) computed by rank-
# normalisation across methods -- this is the headline activity score
# downstream sessions should rank on.
#
# Arguments:
#   stat_mat    numeric matrix from `extract_de_stat_matrix()`. Must
#               have non-zero dimensions; NA cells are tolerated (the
#               decoupleR backend handles them per contrast).
#   network     CollecTRI prior tibble with columns `source` (TF),
#               `target` (gene symbol), `mor` (+1 / -1). Usually
#               `decoupleR::get_collectri(organism = "mouse")`.
#   statistics  character vector of decoupleR method names. Defaults
#               to `c("ulm", "wsum")`; consensus is added separately
#               via `consensus = TRUE`.
#   minsize     drop TFs whose target-set overlap with the row-name
#               universe is smaller than this. Default 5 (same as
#               the project's fgsea convention). With minsize = 5 on
#               the ~11k-symbol nebula universe roughly half of
#               CollecTRI's 1,114 TFs survive (the rest have too few
#               targets in the modality).
#   consensus   logical; whether to append the cross-method consensus
#               score. Default TRUE -- this is the recommended primary
#               score per the decoupleR vignette.
#   seed        integer seed forwarded to `run_consensus()` (the
#               method uses a permutation-style rank normalisation
#               that benefits from a fixed seed for reproducibility).
#               Default 42L.
#
# Returns a tidy tibble:
#   statistic   character; method label.
#   source      character; TF (CollecTRI source).
#   condition   character; contrast name (= matrix column name).
#   score       numeric; signed TF activity score.
#   p_value     numeric; method-specific p-value (see header for the
#               distribution of behaviours across methods, including
#               the corrected note on consensus -- p_value is populated,
#               not NA).
#
# Failure modes:
#   * `stat_mat` 0 rows or 0 cols -> returns an empty tibble with the
#     canonical schema.
#   * No TFs survive `minsize` (very sparse universe) -> empty tibble
#     plus a warning.
run_decoupler_per_modality <- function(stat_mat,
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

  # Per-contrast iteration. decoupleR::decouple() refuses NA/Inf cells
  # in `mat`, but stat matrices built from limma topTables can have
  # asymmetric per-contrast NAs (a small minority of features fail
  # to estimate for some contrasts but estimate cleanly for others
  # because of small-group covariance structure). Iterating one
  # contrast at a time lets each call see only the non-NA subset for
  # that contrast, preserving every contrast's full measurable gene
  # universe rather than dropping any gene with NA in any contrast.
  # The output shape is identical to the single-call form because
  # decouple()'s `condition` column already carries the contrast name.
  contrasts <- colnames(stat_mat)
  per_contrast <- lapply(contrasts, function(cn) {
    col <- stat_mat[, cn]
    keep <- !is.na(col) & is.finite(col)
    if (!any(keep)) {
      warning(sprintf(
        "run_decoupler_per_modality: contrast '%s' has no usable stats",
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
        "run_decoupler_per_modality: decouple returned 0 rows for contrast '%s'",
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
    # decouple()'s internal consensus_score=TRUE call seeds itself;
    # this branch exists for callers that want to override the seed
    # by passing it through. Currently a no-op safeguard: the default
    # decoupleR implementation is deterministic given a fixed prior
    # and a fixed input matrix, so a seed is not actually consumed.
    invisible(seed)
  }

  dec <- dplyr::bind_rows(per_contrast)
  out_cols <- c("statistic", "source", "condition", "score", "p_value")
  missing  <- setdiff(out_cols, names(dec))
  if (length(missing) > 0L) {
    stop(sprintf(
      "run_decoupler_per_modality: decouple output missing columns: %s",
      paste(missing, collapse = ", ")), call. = FALSE)
  }
  tibble::as_tibble(dec[, out_cols, drop = FALSE])
}

# Split a per-modality decoupler tibble into a list keyed by contrast.
# This is the final step before persisting to the cache shape the
# E3+ sessions read from (`modality -> contrast -> tibble`).
#
# Arguments:
#   dec_tbl   output of `run_decoupler_per_modality()`.
#
# Returns a named list of tibbles, one per `condition` level. The
# inner tibble drops the now-redundant `condition` column but keeps
# `statistic / source / score / p_value`. The list preserves the
# order in which contrasts first appear in `dec_tbl`.
split_decoupler_by_contrast <- function(dec_tbl) {
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

# --------------------------------------------------------------------
# Phase E3 helpers: cross-modality TF activity ranking + unified leader
# board. Mirrors the design of `R/pathway_survey.R::rank_pathways_cross_modality`
# and `::build_leader_board` so the TF inference layer can be read
# alongside the pathway-survey layer with parallel mental models.
#
# Modality split convention (locked at E3):
#   primary_modalities    = c("snrnaseq", "geomx", "proteomics")
#                           These three drive the cross-modality count
#                           `n_modalities_sig` (max = 3 with these
#                           defaults). snrnaseq here means whole-
#                           microglia (not a substate).
#   substate_modalities   = c("snrnaseq_homeostatic", "snrnaseq_DAM",
#                             "snrnaseq_IFN", "snrnaseq_proliferative")
#                           These four drive a parallel breadth axis
#                           `substate_breadth` (max = 4), expressing
#                           how compartment-wide a TF signal is. They
#                           are deliberately excluded from the primary
#                           count to avoid collinearity with the whole-
#                           microglia modality.
# This matches the pathway-survey treatment where the GO/custom
# rankings count one snRNAseq modality and report substate breadth as
# a separate axis on `per_state_long`.
#
# Statistic conventions:
#   sig_statistic   = "ulm" -- the linear-model F-test p-value, BH-
#                     adjusted within each (modality, contrast) to give
#                     per-modality padj. Empirical smoke testing during
#                     E3 development showed that using "consensus" for
#                     significance produced a 1-TF leader board across
#                     the project's three primary modalities at FDR<0.10
#                     (consensus is the Stouffer combination of four
#                     methods including the empirical permutation
#                     wsum family, whose p-values floor at 0.02 and
#                     pull the combined p toward conservatism on small
#                     universes). "ulm" is the canonical frequentist
#                     univariate test in the decoupleR vignette and
#                     gives a meaningfully populated leader board
#                     (~10 leader TFs across the two NLGF contrasts at
#                     FDR<0.10) with face-valid microglial biology
#                     (Nr1h3, Spi1/PU.1, Irf1 as headline drivers).
#   score_statistic = "consensus" -- the cross-method consensus is the
#                     recommended primary score in the decoupleR
#                     vignette; it is signed and tracks ulm at 87%
#                     correlation while smoothing single-method
#                     idiosyncrasies. Used for sign attribution
#                     (sign_consensus) and magnitude (mean_abs_score).
#   padj_cut        = 0.10 -- locked by the mechanism_layer_plan E3 spec.
#                     Slightly more permissive than the pathway-survey's
#                     0.05 because TF activity inference has sparser
#                     per-TF evidence (CollecTRI target-set sizes range
#                     from 5 to several hundred, vs gene-set sizes of
#                     >=15 in fgsea), so the per-modality test power is
#                     lower and the threshold compensates.
# --------------------------------------------------------------------

# Internal: extract a per-modality wide table for one contrast from the
# tf_activity_decoupler cache. Each (modality, contrast) cell contains
# the BH-padj for sig_statistic and the score for score_statistic. The
# return value is the building block for `rank_tfs_cross_modality`.
#
# Arguments:
#   tf_cache         the loaded tf_activity_decoupler.rds list
#                    (modality -> contrast -> tibble).
#   contrast         contrast name.
#   modalities       which modalities to include as columns.
#   sig_statistic    decoupler method whose p_value drives the padj
#                    calculation (default "consensus").
#   score_statistic  decoupler method whose score is reported (default
#                    "consensus"). May equal sig_statistic.
#
# Returns a list with three matrices, all rows = TF (sorted alphabetically),
# columns = modalities (in the order given):
#   score_mat  signed activity score (NA where TF is absent from a
#              modality's cache for the contrast, e.g. dropped by
#              minsize=5).
#   padj_mat   BH-adjusted p-value (NA if absent or non-finite p).
#   pval_mat   raw p_value (kept for diagnostics; not used downstream).
.extract_tf_per_modality <- function(tf_cache, contrast, modalities,
                                     sig_statistic   = "ulm",
                                     score_statistic = "consensus") {
  per_mod <- lapply(modalities, function(mn) {
    if (is.null(tf_cache[[mn]]) || is.null(tf_cache[[mn]][[contrast]])) {
      return(data.frame(source = character(0), padj = numeric(0),
                        score = numeric(0), stringsAsFactors = FALSE))
    }
    tbl <- tf_cache[[mn]][[contrast]]
    sig_rows <- tbl[tbl$statistic == sig_statistic, , drop = FALSE]
    if (nrow(sig_rows) == 0L) {
      return(data.frame(source = character(0), padj = numeric(0),
                        score = numeric(0), stringsAsFactors = FALSE))
    }
    sig_rows$padj <- stats::p.adjust(sig_rows$p_value, method = "BH")
    if (score_statistic == sig_statistic) {
      score_rows <- sig_rows
    } else {
      score_rows <- tbl[tbl$statistic == score_statistic, , drop = FALSE]
    }
    # Merge sig + score on `source`; consensus uses the same source set
    # as ulm/wsum so this is a clean inner join in practice.
    merged <- merge(
      sig_rows[, c("source", "p_value", "padj")],
      score_rows[, c("source", "score")],
      by = "source", all = FALSE
    )
    data.frame(
      source = merged$source,
      padj   = merged$padj,
      score  = merged$score,
      pval   = merged$p_value,
      stringsAsFactors = FALSE
    )
  })
  names(per_mod) <- modalities

  all_tfs <- sort(unique(unlist(lapply(per_mod, `[[`, "source"),
                                use.names = FALSE)))
  empty_mat <- function() matrix(NA_real_, nrow = length(all_tfs),
                                 ncol = length(modalities),
                                 dimnames = list(all_tfs, modalities))
  score_mat <- empty_mat()
  padj_mat  <- empty_mat()
  pval_mat  <- empty_mat()
  for (mn in modalities) {
    tbl <- per_mod[[mn]]
    if (nrow(tbl) == 0L) next
    idx <- match(tbl$source, all_tfs)
    score_mat[idx, mn] <- tbl$score
    padj_mat[idx, mn]  <- tbl$padj
    pval_mat[idx, mn]  <- tbl$pval
  }
  list(score_mat = score_mat, padj_mat = padj_mat, pval_mat = pval_mat)
}

# Compute a tidy cross-modality TF ranking for one contrast.
#
# This is the TF analogue of `rank_pathways_cross_modality()`. Each
# row is one TF; columns mirror the pathway-survey schema so the
# downstream `build_tf_leader_board()` can apply the same leader rule.
#
# Arguments:
#   tf_cache         the loaded tf_activity_decoupler.rds list.
#   contrast         contrast name; must appear in every modality's
#                    inner list.
#   primary_modalities  modalities to count in `n_modalities_sig`.
#                    Default snrnaseq + geomx + proteomics (whole-
#                    microglia, not substates). Substates are surfaced
#                    via `compute_tf_substate_breadth()` instead.
#   padj_cut         significance threshold for the per-modality call
#                    (default 0.10; see the module header for rationale).
#   sig_statistic    statistic whose BH-padj drives sig (default
#                    "ulm" -- the linear-model F-test p-value; see
#                    module header for the rationale behind preferring
#                    ulm over consensus for the sig call).
#   score_statistic  statistic whose signed score drives sign + magnitude
#                    (default "consensus" -- the recommended decoupleR
#                    primary metric for ranking).
#
# Returns a tibble with one row per TF and the following columns:
#   source                              character; TF name.
#   n_modalities_sig                    integer; count of primary
#                                       modalities with padj < padj_cut.
#   n_modalities_sig_consistent_sign    integer; among the significant,
#                                       count of the larger sign group.
#   mean_abs_score                      numeric; mean |score| across
#                                       primary modalities with finite
#                                       score.
#   sign_consensus                      "+" / "-" / "mixed" / NA, defined
#                                       analogously to the pathway-survey
#                                       column.
#   score_<modality>                    numeric; per-modality consensus
#                                       score (NA if TF absent).
#   padj_<modality>                     numeric; per-modality padj
#                                       (NA if TF absent).
#   composite_rank                      integer; 1-indexed rank by
#                                       (n_modalities_sig_consistent_sign
#                                        desc, mean_abs_score desc).
rank_tfs_cross_modality <- function(tf_cache, contrast,
                                    primary_modalities = c("snrnaseq",
                                                           "geomx",
                                                           "proteomics"),
                                    padj_cut = 0.10,
                                    sig_statistic   = "ulm",
                                    score_statistic = "consensus") {
  stopifnot(is.list(tf_cache),
            length(contrast) == 1L,
            length(primary_modalities) > 0L,
            all(primary_modalities %in% names(tf_cache)))

  mats <- .extract_tf_per_modality(tf_cache, contrast, primary_modalities,
                                   sig_statistic = sig_statistic,
                                   score_statistic = score_statistic)
  score_mat <- mats$score_mat
  padj_mat  <- mats$padj_mat
  all_tfs   <- rownames(score_mat)

  if (length(all_tfs) == 0L) {
    empty <- tibble::tibble(
      source = character(0),
      n_modalities_sig = integer(0),
      n_modalities_sig_consistent_sign = integer(0),
      mean_abs_score = numeric(0),
      sign_consensus = character(0)
    )
    for (mn in primary_modalities) {
      empty[[paste0("score_", mn)]] <- numeric(0)
      empty[[paste0("padj_",  mn)]] <- numeric(0)
    }
    empty$composite_rank <- integer(0)
    return(empty)
  }

  sig_mat <- !is.na(padj_mat) & padj_mat < padj_cut
  n_sig <- rowSums(sig_mat)

  # Per-TF sign accounting over the SIGNIFICANT subset only. Mirrors
  # the pathway-survey logic exactly: the consistent-sign count is the
  # larger of (positive-sig, negative-sig); the consensus label is
  # strict (any disagreement -> "mixed").
  pos_sig <- vapply(seq_along(all_tfs), function(i) {
    if (n_sig[i] == 0L) return(0L)
    sum(score_mat[i, sig_mat[i, ]] > 0, na.rm = TRUE)
  }, integer(1))
  neg_sig <- vapply(seq_along(all_tfs), function(i) {
    if (n_sig[i] == 0L) return(0L)
    sum(score_mat[i, sig_mat[i, ]] < 0, na.rm = TRUE)
  }, integer(1))

  n_consistent <- pmax(pos_sig, neg_sig)
  sign_consensus <- ifelse(
    n_sig == 0L, NA_character_,
    ifelse(pos_sig > 0L & neg_sig == 0L, "+",
           ifelse(neg_sig > 0L & pos_sig == 0L, "-", "mixed"))
  )

  mean_abs_score <- apply(abs(score_mat), 1, function(r) {
    r <- r[is.finite(r)]
    if (length(r) == 0L) NA_real_ else mean(r)
  })

  out <- tibble::tibble(
    source                            = all_tfs,
    n_modalities_sig                  = as.integer(n_sig),
    n_modalities_sig_consistent_sign  = as.integer(n_consistent),
    mean_abs_score                    = mean_abs_score,
    sign_consensus                    = sign_consensus
  )
  for (mn in primary_modalities) {
    out[[paste0("score_", mn)]] <- score_mat[, mn]
    out[[paste0("padj_",  mn)]] <- padj_mat[, mn]
  }

  out |>
    dplyr::arrange(dplyr::desc(n_modalities_sig_consistent_sign),
                   dplyr::desc(mean_abs_score)) |>
    dplyr::mutate(composite_rank = seq_len(dplyr::n()))
}

# Compute per-(TF, contrast) substate breadth: how many of the four
# microglia substate modalities call the TF significant at padj_cut on
# sig_statistic. Mirror of the pathway-survey's `per_state_long` but
# keyed on TF.
#
# Arguments:
#   tf_cache             the loaded tf_activity_decoupler.rds list.
#   contrasts            character vector of contrast names. NULL -> use
#                        every contrast that appears under every
#                        substate_modality (intersection).
#   substate_modalities  character vector of substate modality names.
#                        Default the four microglia substates.
#   padj_cut             significance threshold (default 0.10).
#   sig_statistic        statistic whose BH-padj drives sig (default
#                        "ulm"; matches the cross-modality ranking
#                        default so substate breadth uses the same
#                        per-modality sig metric).
#
# Returns a long tibble with columns:
#   source             character; TF name.
#   contrast           character; contrast name.
#   substate_breadth   integer 0..length(substate_modalities); count of
#                      substates calling the TF sig at padj_cut.
#   substates_summary  character; pipe-delimited list of substates that
#                      called the TF sig (e.g. "DAM | IFN").
#
# TFs absent from every substate cache at a contrast yield no row at
# that contrast; downstream joins fill missing rows with breadth = 0.
compute_tf_substate_breadth <- function(tf_cache,
                                        contrasts = NULL,
                                        substate_modalities = c(
                                          "snrnaseq_homeostatic",
                                          "snrnaseq_DAM",
                                          "snrnaseq_IFN",
                                          "snrnaseq_proliferative"
                                        ),
                                        padj_cut      = 0.10,
                                        sig_statistic = "ulm") {
  stopifnot(is.list(tf_cache),
            all(substate_modalities %in% names(tf_cache)))

  if (is.null(contrasts)) {
    inner_sets <- lapply(substate_modalities,
                         function(mn) names(tf_cache[[mn]]))
    contrasts <- Reduce(intersect, inner_sets)
  }

  per_contrast <- lapply(contrasts, function(cn) {
    mats <- .extract_tf_per_modality(tf_cache, cn, substate_modalities,
                                     sig_statistic   = sig_statistic,
                                     score_statistic = sig_statistic)
    padj_mat <- mats$padj_mat
    if (nrow(padj_mat) == 0L) {
      return(tibble::tibble(
        source           = character(0),
        contrast         = character(0),
        substate_breadth = integer(0),
        substates_summary = character(0)
      ))
    }
    sig_mat <- !is.na(padj_mat) & padj_mat < padj_cut
    breadth <- as.integer(rowSums(sig_mat))
    # Strip the "snrnaseq_" prefix from substate column names for the
    # readable summary string. Underscore left in to preserve readability
    # of compound substate names if any are added later.
    short_names <- sub("^snrnaseq_", "", colnames(sig_mat))
    summary_str <- vapply(seq_len(nrow(sig_mat)), function(i) {
      cols <- short_names[sig_mat[i, ]]
      if (length(cols) == 0L) NA_character_ else paste(cols, collapse = " | ")
    }, character(1))
    tibble::tibble(
      source             = rownames(padj_mat),
      contrast           = cn,
      substate_breadth   = breadth,
      substates_summary  = summary_str
    )
  })
  dplyr::bind_rows(per_contrast)
}

# Build the unified TF leader board across contrasts. TF-keyed mirror
# of `R/pathway_survey.R::build_leader_board()`. Returns ONE row per
# TF, with cross-contrast aggregate columns. The plan E3 spec lists
# the columns required: n_modalities_sig (-> max_n_modalities_sig),
# dominant_sign, max_abs_score, agreement_pattern (-> contrasts_summary),
# n_contrasts_consistent_sign_ge2, leader_score.
#
# Arguments:
#   ranking_long           long tibble obtained by rbinding the output
#                          of `rank_tfs_cross_modality()` across
#                          contrasts, with an added `contrast` column
#                          and the per-modality score / padj columns
#                          preserved (so downstream sessions can
#                          recompute axis-restricted statistics from
#                          this tibble alone if needed).
#                          Required columns: source, contrast,
#                          n_modalities_sig, n_modalities_sig_consistent_sign,
#                          sign_consensus, mean_abs_score.
#   substate_breadth_long  optional output of `compute_tf_substate_breadth()`.
#                          Used to populate `max_substate_breadth`.
#                          If NULL, `max_substate_breadth` is 0 for
#                          every TF.
#   leader_rule            optional vectorised closure
#                          `function(tbl) -> logical(nrow(tbl))`.
#                          Default identical to the pathway-survey rule:
#                          `n_modalities_sig_consistent_sign >= 2 OR
#                          n_modalities_sig >= 3`.
#
# Returns a tibble of leader rows, sorted by leader_score desc with the
# same tie-break chain as `build_leader_board`. Columns mirror the
# pathway leader board:
#   source                          TF.
#   n_contrasts_leader              count of contrasts where rule fires.
#   n_contrasts_consistent_sign_ge2 count with consistent-sign >= 2.
#   n_contrasts_sig_ge3             count with n_modalities_sig >= 3.
#   max_consistent_sign             max per-contrast consistent-sign count.
#   max_n_modalities_sig            max per-contrast n_modalities_sig.
#   max_abs_score                   max per-contrast mean_abs_score.
#   dominant_sign                   "+" / "-" / "mixed" / NA across leader
#                                   contrasts.
#   max_substate_breadth            max per-contrast substate breadth.
#   contrasts_summary               "<contrast>:<n_modalities_sig>/<sign_consensus>"
#                                   pipe-delimited across leader contrasts.
#   leader_score                    composite, same formula as pathway:
#                                   `5 * n_contrasts_consistent_sign_ge2 +
#                                    n_contrasts_sig_ge3 +
#                                    max_consistent_sign / 5`.
build_tf_leader_board <- function(ranking_long,
                                  substate_breadth_long = NULL,
                                  leader_rule = NULL) {
  stopifnot(is.data.frame(ranking_long))

  required <- c("source", "contrast", "n_modalities_sig",
                "n_modalities_sig_consistent_sign", "sign_consensus",
                "mean_abs_score")
  missing  <- setdiff(required, names(ranking_long))
  if (length(missing) > 0L) {
    stop(sprintf(
      "build_tf_leader_board: ranking_long missing required columns: %s",
      paste(missing, collapse = ", ")), call. = FALSE)
  }

  if (is.null(leader_rule)) {
    leader_rule <- function(row) {
      row$n_modalities_sig_consistent_sign >= 2L |
      row$n_modalities_sig                  >= 3L
    }
  }

  if (nrow(ranking_long) == 0L) {
    return(tibble::tibble(
      source                          = character(0),
      n_contrasts_leader              = integer(0),
      n_contrasts_consistent_sign_ge2 = integer(0),
      n_contrasts_sig_ge3             = integer(0),
      max_consistent_sign             = integer(0),
      max_n_modalities_sig            = integer(0),
      max_abs_score                   = numeric(0),
      dominant_sign                   = character(0),
      max_substate_breadth            = integer(0),
      contrasts_summary               = character(0),
      leader_score                    = numeric(0)
    ))
  }

  is_leader <- as.logical(leader_rule(ranking_long))
  if (length(is_leader) != nrow(ranking_long)) {
    stop("build_tf_leader_board: leader_rule must return a logical ",
         "vector of length nrow(ranking_long).", call. = FALSE)
  }
  is_leader[is.na(is_leader)] <- FALSE
  ranking_long$is_leader <- is_leader

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

  per_tf <- ranking_long |>
    dplyr::group_by(source) |>
    dplyr::summarise(
      n_contrasts_leader              = as.integer(sum(is_leader,
                                                       na.rm = TRUE)),
      n_contrasts_consistent_sign_ge2 = as.integer(sum(
        n_modalities_sig_consistent_sign >= 2L, na.rm = TRUE)),
      n_contrasts_sig_ge3             = as.integer(sum(
        n_modalities_sig >= 3L, na.rm = TRUE)),
      max_consistent_sign             = as.integer(suppressWarnings(
        max(n_modalities_sig_consistent_sign, na.rm = TRUE))),
      max_n_modalities_sig            = as.integer(suppressWarnings(
        max(n_modalities_sig, na.rm = TRUE))),
      max_abs_score                   = suppressWarnings(
        max(mean_abs_score, na.rm = TRUE)),
      dominant_sign                   = reduce_dominant_sign(
        sign_consensus, is_leader),
      contrasts_summary               = reduce_contrasts_summary(
        contrast, n_modalities_sig, sign_consensus, is_leader),
      .groups = "drop"
    ) |>
    dplyr::filter(n_contrasts_leader >= 1L)

  is_neg_inf <- function(x) is.numeric(x) & is.finite(x) == FALSE & x < 0
  per_tf$max_consistent_sign[is_neg_inf(per_tf$max_consistent_sign)]   <- 0L
  per_tf$max_n_modalities_sig[is_neg_inf(per_tf$max_n_modalities_sig)] <- 0L
  per_tf$max_abs_score[is_neg_inf(per_tf$max_abs_score)]               <- NA_real_

  if (!is.null(substate_breadth_long) && nrow(substate_breadth_long) > 0L) {
    sb_required <- c("source", "contrast", "substate_breadth")
    sb_missing  <- setdiff(sb_required, names(substate_breadth_long))
    if (length(sb_missing) > 0L) {
      stop(sprintf(
        "build_tf_leader_board: substate_breadth_long missing columns: %s",
        paste(sb_missing, collapse = ", ")), call. = FALSE)
    }
    sb_max <- substate_breadth_long |>
      dplyr::group_by(source) |>
      dplyr::summarise(
        max_substate_breadth = as.integer(max(substate_breadth, na.rm = TRUE)),
        .groups = "drop"
      )
    per_tf <- per_tf |>
      dplyr::left_join(sb_max, by = "source")
  } else {
    per_tf$max_substate_breadth <- NA_integer_
  }
  per_tf$max_substate_breadth[is.na(per_tf$max_substate_breadth)] <- 0L

  per_tf$leader_score <-
    5 * per_tf$n_contrasts_consistent_sign_ge2 +
    per_tf$n_contrasts_sig_ge3 +
    per_tf$max_consistent_sign / 5

  per_tf |>
    dplyr::arrange(dplyr::desc(leader_score),
                   dplyr::desc(n_contrasts_consistent_sign_ge2),
                   dplyr::desc(max_consistent_sign),
                   dplyr::desc(max_abs_score),
                   source) |>
    dplyr::select(source,
                  n_contrasts_leader,
                  n_contrasts_consistent_sign_ge2,
                  n_contrasts_sig_ge3,
                  max_consistent_sign,
                  max_n_modalities_sig,
                  max_abs_score,
                  dominant_sign,
                  max_substate_breadth,
                  contrasts_summary,
                  leader_score)
}

# Format the top-n rows of a TF cross-modality ranking as a knitr::kable
# for inline display in the Rmd. Mirror of `format_ranking_table` but
# adapted to TF columns (source instead of pathway; no GO prefix strip).
#
# Arguments:
#   ranking_tbl   output of `rank_tfs_cross_modality()`.
#   contrast      contrast string, used in the default caption.
#   n             top-n rows to display (default 20). If
#                 `nrow(ranking_tbl) <= n` every row is shown.
#   include_padj  include per-modality padj columns (default FALSE to
#                 keep the table narrow).
#   caption       overrides the default caption.
#
# Returns a knitr::kable suitable for `print()` inside an
# `results = 'asis'` chunk.
format_tf_ranking_table <- function(ranking_tbl, contrast, n = 20,
                                    include_padj = FALSE,
                                    caption = NULL) {
  top <- head(ranking_tbl, n)
  score_cols <- grep("^score_", names(top), value = TRUE)
  padj_cols  <- grep("^padj_",  names(top), value = TRUE)
  display_cols <- c("composite_rank", "source",
                    "n_modalities_sig",
                    "n_modalities_sig_consistent_sign",
                    "sign_consensus",
                    "mean_abs_score",
                    score_cols)
  if (isTRUE(include_padj)) display_cols <- c(display_cols, padj_cols)

  if (is.null(caption)) {
    caption <- sprintf(
      paste0("Top %d cross-modality TFs for contrast '%s' (sorted by ",
             "sign-consistent significant-modality count, then mean ",
             "|consensus score|; padj cutoff %.2f on consensus p-value, ",
             "BH-adjusted within modality x contrast)."),
      min(n, nrow(top)), contrast, 0.10)
  }

  knitr::kable(top[, display_cols], digits = 3, caption = caption)
}

# --------------------------------------------------------------------
# Phase E4 helpers: axis-restricted TF activity scoring. The cross-
# modality leader board built in E3 (`build_tf_leader_board()`) ranks
# TFs by their *global* multi-modality activity profile and surfaces
# almost exclusively amyloid-activation drivers at the two NLGF
# contrasts. The D2 verdict in `12_pathway_survey.Rmd` named three
# axes:
#   1. amyloid-driven activation   (+ at NLGF contrasts)
#   2. NLGF-driven synaptic suppression  (- at NLGF contrasts)
#   3. mixed-sign metabolic/translational  (mixed at interaction)
# The E4 axis-restricted scoring asks a different question per axis:
# of the TFs whose CollecTRI target set overlaps the axis gene
# universe at >= min_targets distinct genes, which carry the strongest
# mean activity across that axis's relevant contrasts and across the
# three primary modalities?
#
# Design choice (interpretation A of plan E4): the per-TF activity
# score is read from the existing tf_activity_decoupler.rds cache
# (computed against the full per-modality universe in E2). E4
# restricts the TF SET by universe-overlap, then averages the existing
# (full-universe) consensus scores across (axis-contrasts x primary-
# modalities). This is the decoupleR vignette convention for "targeted"
# inference and avoids the statistical-power cost of re-running
# decoupleR per axis against a 2-4k gene subset. Interpretation B
# (re-run decoupleR with the network's targets clipped to the axis
# universe) is the more aggressive variant; it would change the scores
# but at the cost of much smaller TF target sets per modality and
# weaker power. Interpretation A is the locked choice; the module
# header documents the rationale so a future session can revisit if a
# specific axis demands the rerun.
#
# Axis classification rule (locked at E4):
#   amyloid_activation     <- contrasts_summary matches "nlgf_in_[a-z0-9]+:[0-9]+/[+]"
#   synaptic_suppression   <- contrasts_summary matches "nlgf_in_[a-z0-9]+:[0-9]+/-"
#   interaction_metabolic  <- contrasts_summary matches "interaction:[0-9]+/mixed"
# Multi-label (a leader row may match more than one axis; e.g. DAM_up
# is both axis 1 and axis 3). Orphans (no axis) are dropped from the
# universes. The rule is strictly directional (sign-based) and may
# differ from the predecessor pathway-survey prose in borderline cases
# (e.g. AD2 grouped with axis 1 by curatorial provenance but with
# `nlgf_in_maptki:2/-` lands in axis 2 by the directional rule). The
# directional rule is preferred for the mechanism layer because the
# axis universe should reflect the *biological direction* the axis
# represents (+ amyloid up; - synapse down) -- mixing direction within
# an axis dilutes the TF inference signal.

# Build a gene universe per D2 axis from the unified leader-board TSV
# plus the gene-set caches the survey collections pulled from.
#
# Arguments:
#   leaderboard      a data.frame with columns `collection`, `pathway`,
#                    `contrasts_summary`. Typically
#                    `storage/results/pathway_survey_unified_leaderboard.tsv`
#                    re-read with `readr::read_tsv()`.
#   gene_sets        a named list keyed by leaderboard$collection values.
#                    Supported types per element:
#                      * named list `pathway -> character` (msigdb
#                        collections, custom_microglia_states,
#                        custom_microglia_ad).
#                      * data.frame with `gene_name` (char) and
#                        `module` (char or factor) columns -- this is
#                        the hdWGCNA `modules` slot of
#                        `hdwgcna_microglia.rds`; rows are filtered to
#                        the matching module level.
#                    Collections not present in `gene_sets` raise a
#                    stop().
#   axis_rules       optional named list `axis_name -> regex`. Default
#                    matches the strict directional rule (see module
#                    header). Provide a custom list to experiment with
#                    different axis-membership rules (e.g. content-
#                    aware grouping).
#
# Returns a named list `axis -> list(universe, leader_rows,
# axis_contrasts)`:
#   universe         character vector of gene symbols (union across
#                    every leader row passing the axis rule for the
#                    axis). Sorted alphabetical for reproducibility.
#   leader_rows      data.frame of (collection, pathway, contrasts_summary)
#                    triples that contributed to the universe.
#   axis_contrasts   character vector of contrast names extracted from
#                    `contrasts_summary` strings matching the axis rule.
#                    Sorted; deduplicated. Used downstream by
#                    `score_tf_per_axis()` to define which contrasts'
#                    activity scores enter the per-TF mean.
build_axis_gene_universe <- function(leaderboard, gene_sets,
                                     axis_rules = NULL) {
  stopifnot(is.data.frame(leaderboard),
            all(c("collection", "pathway", "contrasts_summary") %in%
                  names(leaderboard)),
            is.list(gene_sets),
            !is.null(names(gene_sets)))

  if (is.null(axis_rules)) {
    axis_rules <- list(
      amyloid_activation    = "nlgf_in_[a-z0-9]+:[0-9]+/[+]",
      synaptic_suppression  = "nlgf_in_[a-z0-9]+:[0-9]+/-",
      interaction_metabolic = "interaction:[0-9]+/mixed"
    )
  }

  # Look up gene symbols for one (collection, pathway) pair. Handles
  # the two supported gene-set encodings (named list, data.frame).
  get_genes <- function(coll, pwy) {
    if (!coll %in% names(gene_sets)) {
      stop(sprintf("build_axis_gene_universe: gene_sets missing collection '%s'",
                   coll), call. = FALSE)
    }
    src <- gene_sets[[coll]]
    if (is.data.frame(src)) {
      # hdWGCNA-style: module column carries the pathway name; filter
      # rows where as.character(module) == pwy and return gene_name.
      stopifnot(all(c("gene_name", "module") %in% names(src)))
      hits <- src$gene_name[as.character(src$module) == pwy]
      as.character(hits)
    } else if (is.list(src)) {
      g <- src[[pwy]]
      if (is.null(g)) {
        return(character(0))
      }
      as.character(g)
    } else {
      stop(sprintf("build_axis_gene_universe: gene_sets[['%s']] must be a list or data.frame",
                   coll), call. = FALSE)
    }
  }

  # Extract the contrast names from a contrasts_summary string that
  # match a given axis regex. The contrasts_summary string is a
  # pipe-delimited list of `<contrast>:<n_modalities_sig>/<sign>` tokens
  # (built by `build_leader_board()`); the contrast is the prefix
  # before the first colon.
  contrasts_matching <- function(summary_str, rgx) {
    if (is.na(summary_str) || !nzchar(summary_str)) return(character(0))
    tokens <- strsplit(summary_str, "\\s*\\|\\s*")[[1]]
    hit_tokens <- grep(rgx, tokens, value = TRUE)
    if (length(hit_tokens) == 0L) return(character(0))
    sub("^([^:]+):.*$", "\\1", hit_tokens)
  }

  out <- lapply(names(axis_rules), function(ax) {
    rgx <- axis_rules[[ax]]
    flag <- grepl(rgx, leaderboard$contrasts_summary)
    rows <- leaderboard[flag, , drop = FALSE]
    if (nrow(rows) == 0L) {
      return(list(universe = character(0),
                  leader_rows = rows[, c("collection", "pathway",
                                         "contrasts_summary"), drop = FALSE],
                  axis_contrasts = character(0)))
    }
    syms <- character(0)
    for (i in seq_len(nrow(rows))) {
      syms <- union(syms, get_genes(rows$collection[i], rows$pathway[i]))
    }
    ctrs <- unique(unlist(lapply(rows$contrasts_summary,
                                 contrasts_matching, rgx = rgx),
                          use.names = FALSE))
    list(
      universe       = sort(syms),
      leader_rows    = rows[, c("collection", "pathway",
                                "contrasts_summary"), drop = FALSE],
      axis_contrasts = sort(unique(ctrs))
    )
  })
  setNames(out, names(axis_rules))
}

# Restrict a CollecTRI TF-target network to TFs whose targets overlap a
# given gene universe at >= min_targets distinct symbols. Returns the
# filtered edges plus a per-source target-count tibble, both keyed
# consistently with the input network's `source` column.
#
# Arguments:
#   network        a data.frame / tibble with columns `source` (TF),
#                  `target` (gene symbol), `mor` (mode of regulation,
#                  +1 / -1). Typically the output of
#                  `decoupleR::get_collectri(organism = "mouse")`.
#   universe       character vector of gene symbols defining the axis
#                  universe. Duplicates are tolerated; the function
#                  uniques internally.
#   min_targets    integer; drop TFs whose distinct in-universe target
#                  count is below this. Default 5.
#
# Returns a list with two named elements:
#   edges          filtered network tibble, restricted to (source,
#                  target) edges whose target is in `universe` AND
#                  whose source passes the min_targets filter. All
#                  original columns of `network` preserved.
#   target_counts  tibble (source, n_targets_in_universe), one row per
#                  surviving TF, sorted by n_targets_in_universe
#                  descending then source alphabetical. Useful as
#                  metadata for `score_tf_per_axis()` and for the
#                  axis-restricted display tables.
restrict_collectri_to_universe <- function(network, universe,
                                            min_targets = 5L) {
  stopifnot(is.data.frame(network) || tibble::is_tibble(network),
            all(c("source", "target") %in% names(network)),
            is.character(universe))

  if (length(universe) == 0L) {
    return(list(
      edges = network[0L, , drop = FALSE],
      target_counts = tibble::tibble(
        source                  = character(0),
        n_targets_in_universe   = integer(0)
      )
    ))
  }

  uni <- unique(universe)
  in_uni <- network$target %in% uni
  edges_in <- network[in_uni, , drop = FALSE]
  if (nrow(edges_in) == 0L) {
    return(list(
      edges = edges_in,
      target_counts = tibble::tibble(
        source                  = character(0),
        n_targets_in_universe   = integer(0)
      )
    ))
  }

  per_tf <- as.data.frame(table(edges_in$source), stringsAsFactors = FALSE)
  names(per_tf) <- c("source", "n_targets_in_universe")
  per_tf$n_targets_in_universe <- as.integer(per_tf$n_targets_in_universe)
  keep <- per_tf$source[per_tf$n_targets_in_universe >= as.integer(min_targets)]

  edges_keep <- edges_in[edges_in$source %in% keep, , drop = FALSE]
  per_tf_keep <- per_tf[per_tf$source %in% keep, , drop = FALSE]
  per_tf_keep <- per_tf_keep[order(-per_tf_keep$n_targets_in_universe,
                                   per_tf_keep$source), , drop = FALSE]
  list(
    edges         = tibble::as_tibble(edges_keep),
    target_counts = tibble::as_tibble(per_tf_keep)
  )
}

# Score and rank TFs per axis using the existing tf_activity_decoupler
# cache. The score for a (TF, axis) cell is the MEAN of the consensus
# (or chosen `score_statistic`) score across every (axis_contrast,
# primary_modality) cell in the cache where the TF has a non-NA score.
# TFs are filtered to those passing the per-axis universe min_targets
# overlap (see `restrict_collectri_to_universe()`).
#
# The ranking within an axis is by |mean_activity| descending so that
# strong-signal TFs (regardless of sign) surface first. The sign of
# `mean_activity` carries the biological direction (axis 1 leaders
# should be positive; axis 2 leaders negative; axis 3 can go either
# way -- the verdict subsection E5 interprets per axis). |mean| is
# the correct sort key because axes can have either sign, and the
# E4 output is meant to be the input to a sign-aware verdict, not the
# verdict itself.
#
# Arguments:
#   tf_cache             the loaded tf_activity_decoupler.rds list
#                        (modality -> contrast -> tibble).
#   axis_universes       output of `build_axis_gene_universe()`. Each
#                        element must provide `universe` (character)
#                        and `axis_contrasts` (character).
#   network              CollecTRI tibble used to compute per-axis
#                        target overlap. Typically the OmnipathR-
#                        cached `decoupleR::get_collectri(organism = "mouse")`.
#   primary_modalities   modalities whose scores enter the per-(TF, axis)
#                        mean. Defaults to the three primary modalities
#                        used in E3. Substate modalities are excluded
#                        for the same reason E3 excluded them: per-
#                        substate universes are smaller and the
#                        scores are heavily correlated with the whole-
#                        microglia modality.
#   min_targets          forwarded to `restrict_collectri_to_universe()`.
#                        Default 5.
#   score_statistic      decoupleR statistic to read for the per-TF
#                        per-(modality, contrast) score. Default
#                        "consensus" -- E3's primary ranking statistic.
#
# Returns a long tibble with one row per (axis, source/TF):
#   axis                          axis name from names(axis_universes).
#   source                        TF (CollecTRI source).
#   mean_activity_in_axis_contrasts  numeric; mean score across
#                                    (axis_contrast, primary_modality)
#                                    cells with finite score.
#   sd_activity_in_axis_contrasts  numeric; sd of the same cells; useful
#                                  for confidence inspection.
#   n_cells_used                   integer; count of finite (contrast,
#                                  modality) cells that entered the mean.
#                                  Range 1 .. (axis_contrasts x
#                                  primary_modalities).
#   n_targets_in_axis_universe     integer; count of distinct in-universe
#                                  targets in CollecTRI for this TF.
#   leader_rank                    integer 1-indexed; rank within the
#                                  axis by |mean_activity| desc, ties
#                                  broken alphabetically by source.
#   score_at_<contrast>            numeric; mean across `primary_modalities`
#                                  of `score_statistic` cells at this
#                                  contrast (one column per unique
#                                  axis_contrast across `axis_universes`;
#                                  NA where the row's axis does not
#                                  include the contrast). Hybrid columns
#                                  added 2026-05-24 per the F4 decision
#                                  gate (locked to "hybrid": both axis-
#                                  mean and per-contrast). The axis-mean
#                                  semantic above is preserved byte-for-
#                                  byte so previously-saved leader ranks
#                                  remain reproducible. Note: rowMeans of
#                                  the `score_at_<contrast>` columns
#                                  generally does NOT equal
#                                  `mean_activity_in_axis_contrasts` when
#                                  cell coverage differs across
#                                  (modality, contrast) pairs -- the axis
#                                  mean is over ALL finite cells, whereas
#                                  the rowMean of per-contrast columns is
#                                  a mean of per-contrast cross-modality
#                                  means and therefore weights each
#                                  contrast equally. Equality holds when
#                                  every TF has the same modality
#                                  coverage at every contrast; under
#                                  partial coverage the divergence is the
#                                  honest weighting artefact, not a bug.
score_tf_per_axis <- function(tf_cache,
                              axis_universes,
                              network,
                              primary_modalities = c("snrnaseq",
                                                     "geomx",
                                                     "proteomics"),
                              min_targets        = 5L,
                              score_statistic    = "consensus") {
  stopifnot(is.list(tf_cache),
            is.list(axis_universes),
            !is.null(names(axis_universes)),
            is.data.frame(network) || tibble::is_tibble(network),
            length(primary_modalities) > 0L,
            all(primary_modalities %in% names(tf_cache)))

  # Pre-compute the union of axis_contrasts across all axes so the per-
  # contrast score columns are uniform across axis rows (NAs where the
  # row's axis doesn't include the contrast). Column ordering follows
  # first appearance in axis_universes.
  all_contrasts <- unique(unlist(
    lapply(axis_universes, function(a) a$axis_contrasts),
    use.names = FALSE
  ))
  per_contrast_cols <- paste0("score_at_", all_contrasts)

  rows <- list()
  for (ax in names(axis_universes)) {
    aux <- axis_universes[[ax]]
    uni <- aux$universe
    ctrs <- aux$axis_contrasts
    if (length(uni) == 0L || length(ctrs) == 0L) {
      next
    }
    restr <- restrict_collectri_to_universe(network, uni,
                                            min_targets = min_targets)
    if (nrow(restr$target_counts) == 0L) next

    # Gather scores per (TF, modality, contrast) from the cache.
    # Build a long table: source, contrast, score. Modality is implicit
    # through the per-modality loop (not retained as a column) because
    # the downstream hybrid aggregation averages across modalities at
    # each contrast anyway and the axis-mean averages across all cells.
    long_list <- list()
    for (mn in primary_modalities) {
      mod_list <- tf_cache[[mn]]
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
    # Restrict to TFs surviving the universe-overlap filter.
    long_df <- long_df[long_df$source %in% restr$target_counts$source, ,
                       drop = FALSE]
    if (nrow(long_df) == 0L) next

    per_tf <- aggregate(score ~ source, data = long_df, FUN = function(v) {
      c(mean_v = mean(v), sd_v = stats::sd(v), n_v = length(v))
    })
    # aggregate() with a multi-output function returns a matrix in the
    # data column; flatten it.
    score_mat <- per_tf$score
    per_tf$mean_activity_in_axis_contrasts <- as.numeric(score_mat[, "mean_v"])
    per_tf$sd_activity_in_axis_contrasts   <- as.numeric(score_mat[, "sd_v"])
    per_tf$n_cells_used                    <- as.integer(score_mat[, "n_v"])
    per_tf$score <- NULL

    # Hybrid per-contrast columns. Aggregate by (source, contrast) to
    # get the cross-modality mean at each axis_contrast, then pivot to
    # source x contrast wide using base R to avoid an extra dependency.
    per_tc <- aggregate(score ~ source + contrast, data = long_df,
                        FUN = mean)
    src_levels <- sort(unique(per_tc$source))
    wide_mat <- matrix(NA_real_,
                       nrow = length(src_levels),
                       ncol = length(ctrs),
                       dimnames = list(src_levels, ctrs))
    for (i in seq_len(nrow(per_tc))) {
      wide_mat[per_tc$source[i], per_tc$contrast[i]] <- per_tc$score[i]
    }
    wide_df <- as.data.frame(wide_mat, stringsAsFactors = FALSE)
    names(wide_df) <- paste0("score_at_", names(wide_df))
    wide_df$source <- rownames(wide_df)
    rownames(wide_df) <- NULL
    per_tf <- merge(per_tf, wide_df, by = "source", all.x = TRUE)

    # Attach target count and rename to the axis-aware column name used
    # downstream in this function and in the exported TSV. The
    # restrict_collectri_to_universe() helper is axis-agnostic so it
    # emits the neutral `n_targets_in_universe` column; we add the
    # `_axis_` qualifier here, where the universe is by construction
    # the axis universe.
    per_tf <- merge(per_tf, restr$target_counts, by = "source", all.x = TRUE)
    names(per_tf)[names(per_tf) == "n_targets_in_universe"] <-
      "n_targets_in_axis_universe"

    # Add NA columns for any contrast not in this axis, so every row in
    # the final long tibble has the same per-contrast column set.
    missing_cols <- setdiff(per_contrast_cols, names(per_tf))
    for (mc in missing_cols) per_tf[[mc]] <- NA_real_

    # Rank within axis by |mean| desc, source alphabetical.
    per_tf <- per_tf[order(-abs(per_tf$mean_activity_in_axis_contrasts),
                           per_tf$source), , drop = FALSE]
    per_tf$leader_rank <- seq_len(nrow(per_tf))
    per_tf$axis <- ax

    rows[[ax]] <- per_tf[, c("axis", "source",
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

# Format the top-n axis-restricted TFs for one axis as a knitr::kable.
# Used inside the Rmd subsection 14.2 to show one table per axis. Since
# the hybrid retrofit (2026-05-24) the helper auto-includes the per-
# contrast score columns relevant to the axis (any `score_at_<contrast>`
# column whose values are not all NA for the axis subset).
#
# Arguments:
#   axis_tbl     output of `score_tf_per_axis()`.
#   axis_name    one of `unique(axis_tbl$axis)`.
#   n            top-n rows to display (default 15).
#   caption      overrides the default caption.
#
# Returns a knitr::kable suitable for `print()` inside an
# `results = 'asis'` chunk.
format_axis_restricted_table <- function(axis_tbl, axis_name, n = 15,
                                          caption = NULL) {
  sub <- axis_tbl[axis_tbl$axis == axis_name, , drop = FALSE]
  if (nrow(sub) == 0L) {
    return(knitr::kable(
      data.frame(message = sprintf(
        "No TFs survive the axis '%s' universe-overlap filter.", axis_name)),
      caption = caption %||% sprintf("Axis-restricted TFs: %s", axis_name),
      row.names = FALSE))
  }
  sub <- sub[order(sub$leader_rank), , drop = FALSE]
  top <- head(sub, n)

  # Identify the per-contrast columns that actually carry data for this
  # axis (the score_tf_per_axis output pads non-axis contrasts with NA;
  # rendering an all-NA column would clutter the table). We test on the
  # subset of rows being displayed (top n) rather than the full axis
  # subset because the subset is what gets rendered; if a contrast is
  # missing across top n but present for some of the lower-ranked rows
  # that's still informative — relax to checking the full axis subset.
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
    TF           = top$source,
    mean_score   = round(top$mean_activity_in_axis_contrasts, 3),
    sd_score     = round(top$sd_activity_in_axis_contrasts,   3),
    n_cells      = top$n_cells_used,
    n_targets    = top$n_targets_in_axis_universe,
    stringsAsFactors = FALSE
  )
  for (col in per_contrast_cols) {
    # Strip the "score_at_" prefix so column headers stay compact.
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
      paste0("Top %d axis-restricted TFs for axis '%s' (sorted by ",
             "|mean consensus score| across axis-relevant contrasts ",
             "and primary modalities; TFs filtered to those with >= 5 ",
             "CollecTRI targets in the axis gene universe). n_cells is ",
             "the count of finite (contrast, modality) cells averaged.",
             "%s"),
      nrow(top), axis_name, pc_part)
  }
  knitr::kable(display, caption = caption, row.names = FALSE)
}

# Plot a lollipop chart of the top-n axis-restricted TFs for one axis.
# Used inside the Rmd subsection 14.2 to show one chart per axis with
# sign-aware colouring.
#
# Arguments:
#   axis_tbl     output of `score_tf_per_axis()`.
#   axis_name    one of `unique(axis_tbl$axis)`.
#   n            top-n rows to display (default 12).
#   title        plot title; default derives from axis_name.
#
# Returns a ggplot object. Rows are sorted by mean_activity desc within
# axis (signed sort, not absolute), so positive-driver and negative-
# driver TFs are visually separated.
plot_axis_lollipop <- function(axis_tbl, axis_name, n = 12, title = NULL) {
  sub <- axis_tbl[axis_tbl$axis == axis_name, , drop = FALSE]
  if (nrow(sub) == 0L) {
    return(ggplot2::ggplot() +
             ggplot2::annotate("text", x = 0, y = 0,
                               label = sprintf("No TFs for axis '%s'.",
                                               axis_name)) +
             ggplot2::theme_void())
  }
  sub <- sub[order(sub$leader_rank), , drop = FALSE]
  top <- head(sub, n)
  # Re-sort by signed mean for the bar order so positive- and
  # negative-driver TFs are visually separated.
  top <- top[order(top$mean_activity_in_axis_contrasts), , drop = FALSE]
  top$source <- factor(top$source, levels = top$source)
  top$sign <- ifelse(top$mean_activity_in_axis_contrasts >= 0,
                     "positive", "negative")

  if (is.null(title)) {
    title <- sprintf("Axis-restricted TFs: %s (top %d by |mean score|)",
                     axis_name, nrow(top))
  }
  ggplot2::ggplot(top,
                  ggplot2::aes(x = mean_activity_in_axis_contrasts,
                               y = source,
                               colour = sign)) +
    ggplot2::geom_segment(ggplot2::aes(x = 0,
                                       xend = mean_activity_in_axis_contrasts,
                                       y    = source,
                                       yend = source),
                          linewidth = 0.6) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                        colour = "grey60") +
    ggplot2::scale_colour_manual(values = c(positive = "#b40426",
                                            negative = "#3a4cc0"),
                                 name = "mean score sign") +
    ggplot2::labs(title = title,
                  x     = "mean consensus score across axis cells",
                  y     = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
}

# --------------------------------------------------------------------
# Phase E5 helper: build the verdict TSV summarising the top TFs per
# axis with cross-modality leader-board cross-reference, signed mean-
# score range, and a human-curated evidence-summary string. Output is
# concise (one row per axis) so an external reviewer or future session
# can read the bottom-line TF reading per axis without re-running the
# E2/E3/E4 pipeline.
#
# Design notes:
#   * The helper does NOT introduce new inference -- it aggregates the
#     E3 cross-modality leader board and the E4 axis-restricted ranking
#     into a per-axis summary row.
#   * The signed mean is taken from `mean_activity_in_axis_contrasts`
#     in `axis_tbl` and rendered in the `top_TF_signs` column with the
#     same convention as the gene-set verdict in section 13 (`+` =
#     positive activator, `-` = negative repressor at the axis's
#     contrast set).
#   * `mean_score_range` is the [min, max] of the signed mean across
#     the top-N TFs, rendered with explicit signs so the reader sees
#     at a glance whether the top block is unidirectional or mixed.
#   * `n_top_in_cross_modality_leaderboard` counts how many of the
#     axis top-N TFs also reached the E3 cross-modality leader rule
#     (>=2 modalities sig consistent OR all 3 sig at any contrast).
#     For axis 3 this column is almost always zero by construction --
#     the cross-modality leader board concentrates on the two NLGF
#     contrasts and contains no TFs whose leader status comes from
#     the interaction contrast in the current build.
#   * `evidence_summary` is a free-text column filled by the caller
#     via `evidence_summaries`, so the editorial prose can be tuned
#     in the Rmd without code changes. Missing entries become
#     NA_character_.
#
# Arguments:
#   axis_tbl              long tibble from `score_tf_per_axis()` (the
#                         input to `format_axis_restricted_table()`).
#   tf_leaderboard        cross-modality leader-board tibble from
#                         `build_tf_leader_board()`.
#   n_top_per_axis        top-N TFs to summarise per axis. Default 5.
#   evidence_summaries    optional named list `axis_name -> string`
#                         containing the prose evidence-summary string
#                         to place in the `evidence_summary` column.
#                         Missing axes get NA_character_.
#
# Returns a data.frame with one row per axis and columns:
#   axis                                    character;
#   top_TFs                                 comma-separated;
#   top_TF_signs                            comma-separated +/-;
#   mean_score_range                        e.g. "[+2.47, +4.32]" (axis-
#                                           mean across top-N);
#   per_contrast_score_range                e.g. "[+2.31, +5.40]" (the
#                                           min/max of per-contrast
#                                           cross-modality means across
#                                           top-N x axis_contrasts cells;
#                                           added by the 2026-05-24
#                                           hybrid retrofit so a reader
#                                           sees at a glance whether the
#                                           axis-mean smooths over
#                                           heterogeneous per-contrast
#                                           signal);
#   per_contrast_summary                    string like
#                                           "nlgf_in_maptki:[+1.85, +5.40]; nlgf_in_p301s:[+2.21, +3.23]"
#                                           giving per-contrast [min, max]
#                                           of top-N scores per axis_contrast.
#   n_top_in_cross_modality_leaderboard     integer (0..n_top_per_axis);
#   evidence_summary                        character (possibly NA).
build_tf_verdict_table <- function(axis_tbl, tf_leaderboard,
                                   n_top_per_axis = 5L,
                                   evidence_summaries = NULL) {
  stopifnot(is.data.frame(axis_tbl) || tibble::is_tibble(axis_tbl),
            all(c("axis", "source", "mean_activity_in_axis_contrasts",
                  "leader_rank") %in% names(axis_tbl)),
            is.data.frame(tf_leaderboard) || tibble::is_tibble(tf_leaderboard),
            "source" %in% names(tf_leaderboard),
            is.numeric(n_top_per_axis), length(n_top_per_axis) == 1L,
            n_top_per_axis >= 1L)

  # Preserve the natural axis encounter order (amyloid_activation,
  # synaptic_suppression, interaction_metabolic) rather than
  # alphabetising; this matches the 1/2/3 ordering used by the
  # predecessor section-13 D2 verdict so the verdict TSV reads in
  # parallel with the pathway-survey prose.
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
        axis                                   = ax,
        top_TFs                                = NA_character_,
        top_TF_signs                           = NA_character_,
        mean_score_range                       = NA_character_,
        per_contrast_score_range               = NA_character_,
        per_contrast_summary                   = NA_character_,
        n_top_in_cross_modality_leaderboard    = 0L,
        evidence_summary                       = evidence_summaries[[ax]] %||%
                                                  NA_character_,
        stringsAsFactors                       = FALSE
      ))
    }

    signs <- ifelse(top$mean_activity_in_axis_contrasts >= 0, "+", "-")
    range_str <- sprintf("[%+0.2f, %+0.2f]",
                         min(top$mean_activity_in_axis_contrasts),
                         max(top$mean_activity_in_axis_contrasts))
    n_in_lb <- sum(top$source %in% tf_leaderboard$source)

    # Hybrid per-contrast summary. Identify the axis-relevant per-
    # contrast columns by dropping those with all-NA values in the
    # axis subset (consistent with the format_axis_restricted_table
    # filtering convention).
    pc_keep <- per_contrast_cols[vapply(per_contrast_cols, function(col) {
      !all(is.na(sub[[col]]))
    }, logical(1))]
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

    data.frame(
      axis                                   = ax,
      top_TFs                                = paste(top$source, collapse = ", "),
      top_TF_signs                           = paste(signs, collapse = ","),
      mean_score_range                       = range_str,
      per_contrast_score_range               = pc_range_str,
      per_contrast_summary                   = pc_summary,
      n_top_in_cross_modality_leaderboard    = as.integer(n_in_lb),
      evidence_summary                       = evidence_summaries[[ax]] %||%
                                                NA_character_,
      stringsAsFactors                       = FALSE
    )
  })
  do.call(rbind, rows)
}
