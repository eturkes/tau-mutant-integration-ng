# Handoff: land parked codex P1-S5 review fixes (TRANSIENT -- delete after landing)

WHY: /codex-review of the S5 close (67b7dbc) raised 14 findings, ALL ACCEPTED (none rejected). The fix set was
applied + live-cache-verified EXCEPT 2 small doc edits, but the run overflowed one window before gating. Main
was reverted to 67b7dbc (gate-green); the fixes are PARKED verbatim on branch `wip-codex-p1s5-review`. This
session = land them small. Authoritative codex text (may not persist): scratchpad
`codex-review-tau-mutant-integration-ng-b2e9.review`.

## Task (one bounded EXECUTE)
Restore the parked edits -> apply 2 pending edits -> gate green -> ONE commit
`microglia (p1 s5 review): ... (codex)` -> clean up.

## 1. Restore the parked fixes (5 files; roadmap.md handled in 2c, NOT here)
    git checkout wip-codex-p1s5-review -- R/microglia.R _microglia.qmd tests/test_microglia.R \
      .agent/history.md .agent/memory.md
(_targets.R is NOT parked -- its fix is 2a, applied fresh on HEAD.)

WHAT the parked edits do (already applied; read the diff if a detail is unclear):
- _microglia.qmd (prose accuracy): interaction "smaller effects" -> "sub-threshold-per-contrast" (real |logFC|,
  median ~1.1; min adj.P ~0.17); "three independent" -> "three complementary"; concordance prose from `flagged`
  + DAM-concordant clause; MDE = nominal per-test 80% power (corrected-discovery floor higher); fit/skip
  substates via knitr::combine_words; within-DAM amyloid DE surfaced (dam_amyloid_de); pruning caveat = "did
  not manufacture the signal by discarding DAM-high nuclei"; sccomp captions/caveat gated !is.null(comp$sccomp).
- R/microglia.R: extractor finite/consistency assertions (umap + z finite; prov substate_table == recomputed
  counts; prune$n_retained == ncol) + ~0.5MB comment.
- tests/test_microglia.R: 4 negative tests (Inf umap, NA z, tampered substate_table, wrong n_retained).
- .agent/history.md + .agent/memory.md: "small-effect" -> "stageR-confirmed" wording; concordance run-to-run note.

## 2. Apply the 2 pending edits
### 2a. _targets.R:82 -- last stale ~2MB copy (R/microglia.R already corrected)
old: reads a ~2MB target, not the 612MB Seurat.
new: reads a ~0.5MB target, not the 612MB Seurat.

### 2b. .agent/memory.md test inventory -- the edit that FAILED before (line-wrap gotcha)
"Current set:" lists only 4 tests; the project has 6 (test_composition.R + test_microglia.R both exist). Prior
failure: "test_io" ends one line, "(loader..." starts the next (newline, not space) -> anchor on a SINGLE-LINE
substring. After step 1 restored memory.md the region is ~L331-334. Mechanical edit (anchor all on one line):
old: theme/scale/concordance). They are data-free
new: theme/scale/concordance), test_composition (propeller/sccomp arm + concordance), test_microglia (reprocess + microglia_report_data extractor/guards). They are data-free

### 2c. .agent/roadmap.md -- wording fix + Active-plan reset (at commit time)
- L42 wording: `123 stageR small-effect` -> `123 stageR-confirmed (real |logFC|, sub-threshold per-contrast FDR)`
- Reset "## Active plan" to: (none) -- P1 CLOSED; next = PLAN P2 (restore the pre-park wording, drop the handoff
  pointer). KEEP the dated ledger entry "P1-S5 codex review PARKED" as the permanent record; optionally append
  "landed <date>" to it.

## 3. GATE (quality contract -- green before the commit)
    bash scripts/check.sh
RISK (tests/test_microglia.R Inf-embedding negative test): it builds a umap reduction via
SeuratObject::CreateDimReducObject(embeddings=emb_inf, ...) then expect_error(microglia_report_data(obj_inf),
"finite"). Confirm the Inf raises in the EXTRACTOR's is.finite assertion, NOT earlier at CreateDimReducObject.
If construction rejects Inf (wrong error site/text), set the Inf AFTER construction instead
(obj@reductions$umap@cell.embeddings[1,1] <- Inf) so the extractor is what fails; fix the test, re-gate.

## 4. Commit + clean up
- gate green -> stage all (restored 5 + _targets.R + roadmap.md + memory.md) -> ONE commit:
  microglia (p1 s5 review): de-overclaim interaction effect-size + concordance/pruning/MDE accuracy + extractor guards + negative tests (codex)
- git branch -D wip-codex-p1s5-review
- delete this file (.agent/p1s5_review_handoff.md)

## Verification facts (live-cache-checked -- REUSE, do not re-verify; for gate-debug only)
- F1 interaction: the 123 stageR-confirmed interaction genes ALL carry |log2FC|>lfc (median ~1.14); they fail
  only the stricter STANDALONE per-contrast FDR (min adj.P ~0.17) -> "sub-threshold-per-contrast", never "small
  effect" (came back STRONGER than codex stated).
- Concordance flag fires on propeller-vs-sccomp DISAGREEMENT; the DAM interaction is CONCORDANT (not flagged);
  only 3 sparse-IFN cells flag. sccomp is run-to-run variable -> all concordance prose from `flagged`, sccomp
  captions/caveat gated !is.null(comp$sccomp).
- Pruning: dropping low-DAM nuclei asymmetrically by genotype can only nudge a genotype's DAM% UP; NLGF_MAPTKI
  loses most (~16.5%). Defensible = "did not manufacture the amyloid->DAM signal by discarding DAM-high nuclei"
  (dropped clusters all low-DAM), NOT "cannot be a pruning artefact".
- MDE = NOMINAL per-test 80% power (not multiple-testing-corrected discovery power); interaction df =
  median(fit$df.total)=24 (eBayes-MODERATED, not raw resid df ~9).
- Within-DAM amyloid DE = 160 genes past the large-effect bar (NLGF_P301S vs P301S) -> backs dam_amyloid_de.
