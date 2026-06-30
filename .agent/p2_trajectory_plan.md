# P2 Interaction trajectory — plan

## Scope
Test the PIVOTAL project claim: mutant tau MODULATES amyloid-driven microglial activation as a
PROGRESSION-RATE synergy along the homeostatic->DAM trajectory — the interaction
`(NLGF_P301S - P301S) - (NLGF_MAPTKI - MAPTKI)` that P1 found STATIC-null (propeller composition +
pseudobulk DE both under-powered, sub-threshold) but NOT absent. Build an activation pseudotime, test the
2x2 interaction on PROGRESSION, and DECOMPOSE composition (more DAM cells) from progression (cells further
along the axis). Input = cached `microglia_annotated` (612MB: harmony reduction + per-cell UCell scores +
substate labels + genotype/batch/genotype_batch). Precondition MET (P1 closed) -> NO external/data gate.
Outputs = `microglia_trajectory` + `trajectory_progression` targets + `_trajectory.qmd` report section.

## Stack (DEFAULT; fully on-lock, reuses the factorial machinery)
v1 recon (Arc M — the EXECUTED analysis that found "the one positive orthogonal interaction") + 2026 SOTA
sweep CONVERGE on the same primary. No 2026 method does replication-aware + factorial-interaction +
along-trajectory in one call; the field composes these, and the composition below is the lockable fit.

- TRAJECTORY = `slingshot` (BioCsoft) on harmony[1:15] (15 dims: 30 let cell-cycle pull proliferative in as a
  spurious intermediate — v1 lesson), CLEAN lineage homeostatic->DAM ONLY (omit IFN + proliferative
  confounds), rooted `start.clus="Homeostatic"`. Run on the BATCH-CORRECTED harmony space (not raw PCA).
  ROBUSTNESS anchor = UCell DAM-minus-Homeostatic SCORE-AXIS (zero-dep, assumption-light pseudotime for a
  genuinely linear gradient); record slingshot-vs-score-axis Spearman rho. ROOTING validated POST-HOC
  (DAM_UCell monotone in pt + canonical markers), framed as an ACTIVATION ORDERING, NOT developmental
  time/potency — v1's 3 potency proxies (RNA entropy, n_genes, CytoTRACE2) all REJECTED
  homeostatic-as-most-potent; direction rests on Slingshot monotonicity + marker biology, never a potency claim.

- INTERACTION (headline) = the TRAJECTORY ANALOGUE of the pseudobulk DE: collapse per-cell pseudotime to the
  16 genotype_batch replicates -> REUSE `factorial_design` + `make_contrast_matrix` -> read the `interaction`
  contrast. Respects 16-unit replication EXACTLY (a test on 16 numbers, 12 residual df, interaction on 1 df) —
  the single-cell DE-calibration lesson (Squair/Crowell/muscat: cell-level tests pseudoreplicate, inflate FDR
  10-100x; aggregate-to-replicate restores nominal FDR). `condiments::progressionTest` / `tradeSeq::conditionTest`
  treat cells as INDEPENDENT and expose NO factorial-interaction contrast -> NOT the inferential test (descriptive
  marginal only; dropped from the core to stay lean).

- DECOMPOSITION = Kitagawa/Oaxaca shift-share (v1-validated, exact additive counterfactuals): observed mean_pt =
  sum_s pi_s * mu_s -> composition_cf (vary state fractions pi, hold within-state means mu at global) +
  progression_cf (vary mu, hold pi) -> the interaction contrast on EACH counterfactual channel isolates
  composition-driven vs progression-driven synergy. + within-state means (within_homeostatic = the CLEANEST
  "not composition": even homeostatic-labelled cells advancing cannot be an abundance artefact). v1 result to
  re-derive (R4.6, NOT v1's locked margins): interaction loaded PROGRESSION (~94%, sig), NULL on composition.

- SENSITIVITY = `glmmTMB` beta GLMM `pt01 ~ tau*amyloid + (1|genotype_batch)` on PER-CELL pseudotime (CRAN; TMB =
  C++ template, NOT Stan -> ON-lock) — a replication-aware per-cell confirmation that models the full bounded
  (possibly bimodal) distribution the summary collapses; Wald/LRT on the `tau:amyloid` term. Mirrors the
  composition arm's primary+cross-check pattern, but BOTH arms on-lock here. Singular-fit risk if the interaction
  variance ~ 0 -> fallback rank-normal LMM; degrade gracefully (the limma-summary primary stands alone).

LOCK POSTURE: the ENTIRE primary stack is pure-R from the P3M 2026-06-22 (Bioc 3.23 / CRAN) snapshot — NO Stan,
NO Python, NO GitHub. A leanness + reproducibility WIN over v1 (which carried Python + GitHub off-lock deps).

## Dropped from v1 (bloat; re-derive VALUE not machinery)
- Python triangulation (scanpy PAGA + DPT, CellRank2) — a reproducibility liability off the R snapshot, and it
  buys nothing for snapshot snRNA (no spliced/unspliced for velocity). Replaced by the pure-R score-axis anchor.
- CytoTRACE2 root-validation (GitHub dep) — proxies REJECTED homeostatic-potency anyway; we frame an
  activation ordering and validate rooting by marker monotonicity, no potency claim.
- Arc O gene-level tradeSeq fitGAM interaction — leaned on FRAGILE unexported internals
  (`tradeSeq:::predictGAM` / `:::.getPredictRangeDf`), margin-neutral (110 genes, ~1/4 ambient/ribosomal
  contaminated, Gsk3b/Myc ABSENT). Gene-level interaction is P3 mechanism's job (or stays with the static
  pseudobulk DE). OUT of P2.

## Steps (each one-window-closeable; gate-independent — microglia_annotated already cached)

S1 — Trajectory + pseudotime target.
  NEW `R/trajectory.R` pure helpers: `build_activation_trajectory` (slingshot harmony[1:15], homeostatic->DAM,
  rooted; seed + provenance) | `score_axis_pseudotime` (UCell DAM-minus-Homeostatic) | `trajectory_concordance`
  (Spearman slingshot-vs-score-axis) | provenance (dims/root/pkg versions/rho, mirror reprocess_provenance).
  Target `microglia_trajectory` = COMPACT per-cell frame {cell, genotype_batch, genotype, substate, pseudotime,
  score_axis_pt, DAM_UCell, Homeostatic_UCell} + lineage/provenance lists (~small; NEVER the 612MB object —
  cheap-render invariant). Unit-test pure helpers on a synthetic 2-cluster embedding fixture (lineage present,
  pseudotime monotone, aggregation shapes). LIVE smoke on real microglia_annotated before commit.
  ACCEPTANCE: single clean homeostatic->DAM lineage; pseudotime monotone in DAM_UCell; slingshot-vs-score-axis
  rho recorded (concordance OR honest flag); gate green.

S2 — Progression interaction + decomposition target.
  `pseudotime_per_replicate` (16-unit summaries: mean_pt, median_pt, frac_past[DAM-onset = global DAM-cell
  median pt], q90 leading-edge, within_<state> means) | `fit_trajectory_contrasts` (REUSE factorial_design +
  make_contrast_matrix -> 5 contrasts per measure; pt-scale measures UNtransformed, bounded [0,1] measures
  logit/asin -> SEPARATE fits; NB fit_limma_log's intensity-log is inappropriate for pseudotime -> add a
  transform-arg/`fit_limma_summary` variant OR call limma::lmFit directly on the measures x 16 matrix; reuse
  ONLY the design + contrast builders) | `decompose_progression_vs_composition` (Kitagawa channels) |
  glmmTMB per-cell sensitivity. Target `trajectory_progression` = contrast tables + decomposition channels +
  glmmTMB concordance + provenance. PRE-REGISTER (state thresholds BEFORE applying): primary progression
  endpoint = progression_cf (Kitagawa) + within_homeostatic (composition-robust); frac_past = the interpretable
  bridge metric; mean_pt reported but FLAGGED composition-conflated. FDR primary 0.05, REPORT the 0.05-0.10 band
  + effect + CI (P1 lesson: effect/CI/MDE over a bare binary; re-baseline R4.6 != v1's 0.077/0.013).
  ACCEPTANCE: interaction contrast on every measure + decomposition isolates progression-vs-composition loading +
  glmmTMB `tau:amyloid` directionally concords; honest report if the rebuild diverges from v1; gate green.

S3 — Report + integration.
  `_trajectory.qmd` (trajectory UMAP x pseudotime; pseudotime density/ridge by genotype; interaction FOREST
  [contrasts x measures]; composition-vs-progression decomposition bar; headline progression_cf/frac_past
  interaction; activation-ordering + rooting + transcriptionally-close-substates caveats) + compact
  `trajectory_report_data` extractor (slim frame + the contrast/decomposition tables; cheap-render) + wire into
  index.qmd AFTER _microglia.qmd + update _microglia.qmd's "P2 pointer" -> built. Prose INLINE-COMPUTED from the
  loaded targets (never hardcoded — tracks the cached build). ACCEPTANCE: report renders 0-warning under warn=2;
  headline + decomposition + caveats present; gate green. -> then CLOSE-OUT.

## Reproducibility / risks
- ADD rproject.toml: `slingshot` (BioCsoft; auto-pulls SingleCellExperiment / princurve / TrajectoryUtils —
  pure-R Bioc), `glmmTMB` (CRAN). VERIFY present in the 2026-06-22 snapshot on first `rv add`;
  glmmTMB<->Matrix<->TMB ABI must co-pin (source-compile all 3 from the snapshot -> no binary mismatch warning).
- slingshot determinism: principal-curve is deterministic given embedding + labels (+ seed); the harmony
  embedding is the stable P1 cached input. Record seed + provenance (S1 mirror of reprocess_provenance).
- gate cheap-render: BOTH trajectory targets COMPACT; no qmd tar_loads the 612MB object.
- glmmTMB can degrade (singular fit) like the sccomp arm -> the locked limma-summary + decomposition is the
  standalone primary; the sensitivity is supportive, not load-bearing.

## DECIDED (gate 2026-06-30) — LEAN ON-LOCK
slingshot + UCell score-axis concordance anchor; limma per-replicate-summary interaction + Kitagawa
decomposition + glmmTMB per-cell sensitivity. ENTIRELY pure-R from the pinned snapshot — no Stan/Python/GitHub,
no destiny DPT, no off-lock Lamian arm. The primary is already replication-correct (16-unit summary test); the
score-axis + glmmTMB are the two independent cross-checks. (Rejected: +destiny DPT triangulation, +off-lock
Lamian XCD — both weighed against the lean-rebuild ethos.)
