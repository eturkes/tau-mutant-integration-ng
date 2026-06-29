#!/usr/bin/env Rscript
# build_biological_model_adjudication.R
#
# H3 deliverable: aggregate the H2 claims ledger by entity and emit the
# adjudication outputs that feed the H4 section-17 chapter and the H5
# upfront mirror.
#
# Inputs:
#   storage/results/biological_model_claims_ledger.tsv (87 rows x 14 cols
#     after Phase L; built by scripts/build_biological_model_ledger.R).
#     Row count: 75 Phase H rows (H2-001..H2-075) + 6 Phase I rows
#     (I-001..I-006 added 2026-05-25 per the per-state NF-kB attenuation
#     test in section 18 of analysis.html) + 2 Phase J rows (J-001..J-002
#     added 2026-06-04 per the causal signalling-network reconstruction in
#     section 20 of analysis.html) + 2 Phase K rows (K-001..K-002 added
#     2026-06-05 per the SCENIC data-driven regulon recovery in section 21)
#     + 2 Phase L rows (L-001..L-002 added 2026-06-05 per the GeoMx
#     spatial-deconvolution tissue-composition layer in section 22).
#
# Outputs:
#   storage/results/biological_model_adjudication.tsv (11 entities x
#     ~18 arithmetic columns; per-entity support / contradict counts
#     overall + per-axis wide breakdown + per-confidence-grade counts).
#   storage/results/biological_model_contest_verdicts.tsv (3 per-axis
#     contests; favoured model per contest with the tie-break rule that
#     was invoked, if any).
#
# Adjudication arithmetic (locked at H1 2026-05-25; see
# storage/notes/biological_model_plan.md Open-questions row 4):
#   - Per-entity: support_count, contradict_count, net_support,
#     support_to_contradict_ratio, n_strong_supports,
#     n_moderate_supports, n_suggestive_supports.
#   - Per-entity x axis: same arithmetic restricted to claims whose
#     `axis` matches; emitted as wide columns
#     support_at_<axis> / contradict_at_<axis> for the four axes
#     amyloid_activation, synaptic_suppression, interaction_metabolic,
#     cross_axis.
#   - Per-contest favoured-model verdict (3 per-axis contests only):
#     model with higher net_support wins; ties broken by
#     n_strong_supports; further ties by n_moderate_supports.
#     Hyp-0 + four cross-axis themes are adjudicated as stand-alone counts
#     in the per-entity table; they do not appear in the verdict table.
#
# Schema invariants enforced on every load (so the adjudication cannot
# be silently corrupted if the H2 ledger is regenerated):
#   - All entity IDs referenced in the ledger's supports_models /
#     contradicts_models columns are members of the locked 11-entity set.
#   - confidence_grade in {Strong, Moderate, Suggestive}.
#   - axis in {amyloid_activation, synaptic_suppression,
#     interaction_metabolic, cross_axis}.

suppressPackageStartupMessages(library(utils))

# ---- locked vocabulary ----------------------------------------------------

# 11 entities, same order as the H2 ledger builder.
ENTITY_IDS <- c(
  "Hyp-1A", "Hyp-1B", "Hyp-2A", "Hyp-2B", "Hyp-3A", "Hyp-3B", "Hyp-0",
  "T-Inflammation", "T-Compartment-suppression",
  "T-Tau-attenuates", "T-Synergy"
)

# Entity catalogue: structural type, contest membership, human-readable
# name. Column-bound to the per-entity adjudication table so that
# downstream H4 rendering can group / filter without re-deriving the
# typology.
ENTITY_CATALOGUE <- data.frame(
  entity_id = c(
    "Hyp-1A", "Hyp-1B",
    "Hyp-2A", "Hyp-2B",
    "Hyp-3A", "Hyp-3B",
    "Hyp-0",
    "T-Inflammation",
    "T-Compartment-suppression",
    "T-Tau-attenuates",
    "T-Synergy"
  ),
  entity_type = c(
    "contest_model", "contest_model",
    "contest_model", "contest_model",
    "contest_model", "contest_model",
    "integrator",
    "theme", "theme", "theme", "theme"
  ),
  contest_id = c(
    "amyloid_activation",   "amyloid_activation",
    "synaptic_suppression", "synaptic_suppression",
    "interaction_metabolic","interaction_metabolic",
    NA_character_,
    NA_character_, NA_character_, NA_character_, NA_character_
  ),
  entity_name = c(
    "Tau-independent amyloid program",
    "Tau attenuates amyloid-driven NF-kB",
    "Classical complement-pruning hypothesis",
    "TREM2 / APP-fragment-mediated clearance",
    "No interaction-specific mechanism",
    "Distinct synergy mechanism (Gsk3b / Myc / axon-guidance)",
    "Cdk5 cross-axis integrator",
    "Additive DAM-amplification across axes",
    "Synaptic compartment suppressed via classical complement",
    "Tau modifies amyloid rather than amplifying it",
    "Qualitatively new mechanism at the tau x amyloid interface"
  ),
  stringsAsFactors = FALSE
)
# Sanity check: ENTITY_CATALOGUE rows match ENTITY_IDS order.
stopifnot(identical(ENTITY_CATALOGUE$entity_id, ENTITY_IDS))

AXES <- c(
  "amyloid_activation",
  "synaptic_suppression",
  "interaction_metabolic",
  "cross_axis"
)

GRADES <- c("Strong", "Moderate", "Suggestive")

# ---- load ledger ----------------------------------------------------------

in_path  <- file.path("storage", "results", "biological_model_claims_ledger.tsv")
if (!file.exists(in_path)) {
  stop(sprintf(
    "Ledger TSV not found at %s; run scripts/build_biological_model_ledger.R first.",
    in_path
  ))
}

# Read tab-separated, no quoting (H2 builder writes with quote = FALSE),
# treat empty string cells as NA (matches the empty-supports semantics
# used by the H2 builder for rows with no supports / contradicts).
ledger <- read.delim(in_path, header = TRUE, sep = "\t",
                     stringsAsFactors = FALSE, na.strings = "",
                     quote = "", comment.char = "",
                     check.names = FALSE)

# ---- input validation (defensive: catches drift if H2 is regenerated) ----

stopifnot(nrow(ledger) > 0)
stopifnot(all(c("claim_id", "axis", "confidence_grade",
                "supports_models", "contradicts_models") %in%
              colnames(ledger)))
stopifnot(!anyDuplicated(ledger$claim_id))
stopifnot(all(ledger$axis              %in% AXES))
stopifnot(all(ledger$confidence_grade  %in% GRADES))

parse_models <- function(s) {
  # Returns a character vector of entity IDs cited in a semicolon-joined
  # cell. NA / empty cell -> character(0).
  if (length(s) != 1L) {
    stop("parse_models expects a single string.")
  }
  if (is.na(s) || s == "") return(character(0))
  ids <- strsplit(s, ";", fixed = TRUE)[[1]]
  ids <- trimws(ids)
  ids[nzchar(ids)]
}

# Validate that every ID in supports_models / contradicts_models is a
# member of ENTITY_IDS. (The H2 builder runs the same check at write
# time; we re-run it on load so a stale ledger cannot leak through.)
walk_validate <- function(col_name) {
  bad <- character(0)
  for (i in seq_len(nrow(ledger))) {
    ids <- parse_models(ledger[[col_name]][i])
    if (!all(ids %in% ENTITY_IDS)) {
      bad <- c(bad, sprintf("%s row %s (%s): %s",
                            col_name, ledger$claim_id[i],
                            ledger$axis[i],
                            paste(setdiff(ids, ENTITY_IDS), collapse = ",")))
    }
  }
  if (length(bad) > 0L) {
    stop("Out-of-vocab entity IDs in ledger:\n  ",
         paste(bad, collapse = "\n  "))
  }
}
walk_validate("supports_models")
walk_validate("contradicts_models")

# ---- core arithmetic helpers ----------------------------------------------

# Per-claim membership matrices: rows = claims, cols = entities, values =
# logical "this claim cites this entity in supports_models /
# contradicts_models". Built once; sliced repeatedly below.
mat_logical <- function(col_name) {
  m <- matrix(FALSE, nrow = nrow(ledger), ncol = length(ENTITY_IDS),
              dimnames = list(ledger$claim_id, ENTITY_IDS))
  for (i in seq_len(nrow(ledger))) {
    ids <- parse_models(ledger[[col_name]][i])
    if (length(ids) > 0L) m[i, ids] <- TRUE
  }
  m
}
supports_mat    <- mat_logical("supports_models")
contradicts_mat <- mat_logical("contradicts_models")

# Per-entity overall counts.
support_count    <- colSums(supports_mat)
contradict_count <- colSums(contradicts_mat)
net_support      <- support_count - contradict_count

# Support-to-contradict ratio with explicit handling of the zero-
# denominator edge case. The realised ledger has no entity with
# (support = 0 AND contradict = 0); the safety branch below returns NA
# for that pathological case rather than NaN to keep the TSV
# self-explanatory.
ratio <- function(s, c) {
  if (c == 0L && s == 0L) return(NA_real_)
  if (c == 0L)            return(Inf)
  s / c
}
support_to_contradict_ratio <- mapply(ratio, support_count, contradict_count,
                                      USE.NAMES = FALSE)

# Per-(entity, grade) support counts.
grade_counts <- function(grade) {
  idx <- which(ledger$confidence_grade == grade)
  if (length(idx) == 0L) return(setNames(integer(length(ENTITY_IDS)), ENTITY_IDS))
  colSums(supports_mat[idx, , drop = FALSE])
}
n_strong_supports     <- grade_counts("Strong")
n_moderate_supports   <- grade_counts("Moderate")
n_suggestive_supports <- grade_counts("Suggestive")

# Per-(entity, axis) support / contradict counts. Wide-emitted below.
axis_counts <- function(mat, axis_name) {
  idx <- which(ledger$axis == axis_name)
  if (length(idx) == 0L) return(setNames(integer(length(ENTITY_IDS)), ENTITY_IDS))
  colSums(mat[idx, , drop = FALSE])
}
support_at <- lapply(AXES, function(a) axis_counts(supports_mat, a))
names(support_at) <- AXES
contradict_at <- lapply(AXES, function(a) axis_counts(contradicts_mat, a))
names(contradict_at) <- AXES

# ---- per-entity adjudication table ----------------------------------------

adjudication <- data.frame(
  entity_id                   = ENTITY_CATALOGUE$entity_id,
  entity_type                 = ENTITY_CATALOGUE$entity_type,
  contest_id                  = ENTITY_CATALOGUE$contest_id,
  entity_name                 = ENTITY_CATALOGUE$entity_name,
  support_count               = unname(support_count),
  contradict_count            = unname(contradict_count),
  net_support                 = unname(net_support),
  support_to_contradict_ratio = support_to_contradict_ratio,
  n_strong_supports           = unname(n_strong_supports),
  n_moderate_supports         = unname(n_moderate_supports),
  n_suggestive_supports       = unname(n_suggestive_supports),
  stringsAsFactors            = FALSE
)

# Append wide per-axis columns in a canonical order: for each axis, the
# support column then the contradict column. Keeps related figures
# adjacent for human readers and for any DT::datatable column-grouping
# the H4 chapter might do.
for (a in AXES) {
  adjudication[[paste0("support_at_",    a)]] <- unname(support_at[[a]])
  adjudication[[paste0("contradict_at_", a)]] <- unname(contradict_at[[a]])
}

# Sort: contest models grouped by contest (Hyp-1A/Hyp-1B, Hyp-2A/Hyp-2B, Hyp-3A/Hyp-3B)
# first; then Hyp-0; then themes ordered by net_support descending. This
# is the canonical reading order for the H4 17.3 sub-table.
contest_block <- adjudication[adjudication$entity_type == "contest_model", ]
contest_block <- contest_block[order(
  match(contest_block$contest_id,
        c("amyloid_activation", "synaptic_suppression", "interaction_metabolic")),
  contest_block$entity_id
), ]
m0_block      <- adjudication[adjudication$entity_type == "integrator", ]
theme_block   <- adjudication[adjudication$entity_type == "theme", ]
theme_block   <- theme_block[order(-theme_block$net_support,
                                   -theme_block$n_strong_supports), ]
adjudication  <- rbind(contest_block, m0_block, theme_block)

out_path_a <- file.path("storage", "results", "biological_model_adjudication.tsv")
write.table(adjudication, out_path_a, sep = "\t", quote = FALSE,
            row.names = FALSE, na = "", fileEncoding = "UTF-8")

# ---- per-axis-contest favoured-model verdict table ------------------------

# Helper: per-entity grade counts as a named integer vector for one
# entity ID. Used by the tie-break logic below.
strong_for    <- function(eid) n_strong_supports[[eid]]
moderate_for  <- function(eid) n_moderate_supports[[eid]]

verdict_row <- function(contest_id, model_a, model_b) {
  s_a <- support_count[[model_a]];    s_b <- support_count[[model_b]]
  c_a <- contradict_count[[model_a]]; c_b <- contradict_count[[model_b]]
  n_a <- s_a - c_a;                   n_b <- s_b - c_b

  # Tie-break order: net_support -> n_strong_supports -> n_moderate_supports.
  favoured <- NA_character_
  tie_rule <- NA_character_
  if (n_a > n_b) {
    favoured <- model_a; tie_rule <- "no_tie_break_needed"
  } else if (n_b > n_a) {
    favoured <- model_b; tie_rule <- "no_tie_break_needed"
  } else {
    # net_support tie: try n_strong_supports
    sa <- strong_for(model_a); sb <- strong_for(model_b)
    if (sa > sb) {
      favoured <- model_a; tie_rule <- "n_strong_supports"
    } else if (sb > sa) {
      favoured <- model_b; tie_rule <- "n_strong_supports"
    } else {
      ma <- moderate_for(model_a); mb <- moderate_for(model_b)
      if (ma > mb) {
        favoured <- model_a; tie_rule <- "n_moderate_supports"
      } else if (mb > ma) {
        favoured <- model_b; tie_rule <- "n_moderate_supports"
      } else {
        favoured <- "TIE"; tie_rule <- "could_not_break_tie"
      }
    }
  }

  data.frame(
    contest_id                  = contest_id,
    model_a_id                  = model_a,
    model_a_name                = ENTITY_CATALOGUE$entity_name[ENTITY_CATALOGUE$entity_id == model_a],
    model_a_support_count       = s_a,
    model_a_contradict_count    = c_a,
    model_a_net_support         = n_a,
    model_a_n_strong_supports   = strong_for(model_a),
    model_b_id                  = model_b,
    model_b_name                = ENTITY_CATALOGUE$entity_name[ENTITY_CATALOGUE$entity_id == model_b],
    model_b_support_count       = s_b,
    model_b_contradict_count    = c_b,
    model_b_net_support         = n_b,
    model_b_n_strong_supports   = strong_for(model_b),
    favoured_model              = favoured,
    favoured_by_margin          = if (favoured == "TIE") NA_integer_
                                  else abs(n_a - n_b),
    tie_break_rule_invoked      = tie_rule,
    stringsAsFactors            = FALSE
  )
}

verdicts <- rbind(
  verdict_row("amyloid_activation",    "Hyp-1A", "Hyp-1B"),
  verdict_row("synaptic_suppression",  "Hyp-2A", "Hyp-2B"),
  verdict_row("interaction_metabolic", "Hyp-3A", "Hyp-3B")
)

out_path_v <- file.path("storage", "results", "biological_model_contest_verdicts.tsv")
write.table(verdicts, out_path_v, sep = "\t", quote = FALSE,
            row.names = FALSE, na = "", fileEncoding = "UTF-8")

# ---- summary to stdout (for log capture + cross-check vs H2 preview) ------

cat(sprintf("Wrote %s (%d entities x %d cols).\n",
            out_path_a, nrow(adjudication), ncol(adjudication)))
cat(sprintf("Wrote %s (%d contests x %d cols).\n",
            out_path_v, nrow(verdicts), ncol(verdicts)))

cat("\nPer-entity headline arithmetic:\n")
preview <- adjudication[, c("entity_id", "entity_type", "contest_id",
                            "support_count", "contradict_count", "net_support",
                            "n_strong_supports", "n_moderate_supports",
                            "n_suggestive_supports")]
print(preview, row.names = FALSE)

cat("\nPer-axis contest verdicts:\n")
v_preview <- verdicts[, c("contest_id",
                          "model_a_id", "model_a_net_support",
                          "model_b_id", "model_b_net_support",
                          "favoured_model", "favoured_by_margin",
                          "tie_break_rule_invoked")]
print(v_preview, row.names = FALSE)

cat("\nPer-axis x per-entity wide breakdown (support_at_<axis>):\n")
wide_preview <- adjudication[, c("entity_id",
                                 paste0("support_at_", AXES))]
print(wide_preview, row.names = FALSE)
