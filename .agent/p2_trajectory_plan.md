# P2 Interaction trajectory — plan

## Scope
Test the PIVOTAL project claim: mutant tau MODULATES amyloid-driven microglial activation, synergistically
ADVANCING cells along the homeostatic->DAM trajectory — the interaction
`(NLGF_P301S - P301S) - (NLGF_MAPTKI - MAPTKI)`. P1 ESTABLISHED a static COMPOSITIONAL synergy on DAM
(propeller DAM-fraction interaction SIG, fdr_global~0.034; sccomp BORDERLINE, c_fdr~0.044), while the
pseudobulk DE interaction stayed under-powered/sub-threshold (123 stageR-confirmed genes, NO large-effect
call). Composition (MORE DAM cells under the interaction) is therefore ALREADY known -> P2's DISTINCT job =
test whether cells also advance FURTHER along the axis (PROGRESSION) BEYOND that composition shift, and
DECOMPOSE the two channels. MEASURED quantity = position/extent of advance along an inferred CROSS-SECTIONAL
ordering; "rate/acceleration" is the biological INTERPRETATION under the age-matched-snapshot + common-baseline
assumption (stated, NOT measured — no longitudinal time / RNA velocity here). Input = cached
`microglia_annotated` (612MB: harmony reduction + per-cell UCell scores + substate labels +
genotype/batch/genotype_batch). Precondition MET (P1 closed) -> NO external/data gate. Outputs =
`microglia_trajectory` + `trajectory_progression` targets + `_trajectory.qmd` report section.

## Stack (DEFAULT; fully on-lock, reuses the factorial machinery)
v1 recon (Arc M — the EXECUTED analysis that found "the one positive orthogonal interaction") + 2026 SOTA
sweep CONVERGE on the same primary. No 2026 method does replication-aware + factorial-interaction +
along-trajectory in one call; the field composes these, and the composition below is the lockable fit.
(Swept + rejected for the INFERENTIAL test: `condiments::progressionTest` / `tradeSeq::conditionTest` = cells
independent, no factorial interface; Lamian = off-lock; destiny/DPT diffusion + PILOT/scFates = off-lock + add
nothing for a single near-linear axis. Kept descriptively-marginal at most.)

- TRAJECTORY = `slingshot` (BioCsoft) on harmony[1:15] (15 dims: 30 let cell-cycle pull proliferative in as a
  spurious intermediate — v1 lesson), CLEAN lineage homeostatic->DAM ONLY (omit IFN + proliferative
  confounds), rooted `start.clus="Homeostatic"`. Run on the BATCH-CORRECTED harmony space (not raw PCA).
  FORK CHOICES (dims=15, IFN/prolif pruning, DAM-onset threshold) are HIGH-LEVERAGE -> PRE-DECLARE as PRIMARY +
  a FIXED sensitivity table (dims {10,15,20}; all-cells-retained re-fit), NO retuning after contrast inspection.
  CONDITIONING on the H->D lineage is a potential SELECTION effect (if the omitted IFN/prolif fraction differs
  by genotype the summaries shift) -> store n_on_lineage + omitted fraction PER UNIT + an all-retained
  sensitivity. CONCORDANCE check = UCell DAM-minus-Homeostatic SCORE-AXIS (zero-dep, assumption-light) -> record
  Spearman rho; it SHARES the marker system that defines the substate labels, so it catches gross trajectory
  failure, NOT shared-marker bias -> a concordance check, NOT independent robustness. ROOTING validated POST-HOC
  (DAM_UCell monotone in pt + canonical markers), framed as an ACTIVATION ORDERING, NOT developmental
  time/potency — v1's 3 potency proxies (RNA entropy, n_genes, CytoTRACE2) all REJECTED
  homeostatic-as-most-potent; direction rests on Slingshot monotonicity + marker biology, never a potency claim.

- INTERACTION (headline) = the TRAJECTORY ANALOGUE of the pseudobulk DE: collapse per-cell pseudotime to the
  16 genotype_batch replicate SUMMARIES -> REUSE `factorial_design` ALONE (treatment-coded design + its 5
  canonical contrasts incl. `interaction`; do NOT also call `make_contrast_matrix` — they are TWO ALTERNATIVE
  parameterisations of the same column space, pick ONE). Design = 16 x 7 (intercept + tau + nlgf + tau_nlgf +
  3 batch) -> 9 ORDINARY residual df, interaction on 1 df (matches de_pb's documented ~9 df; never call a
  moderated df "residual df"). Aggregating to 16 units REMOVES the dominant cell-level pseudoreplication
  (Squair/Crowell/muscat: cell-level tests inflate FDR 10-100x) but unequal cells/unit -> HETEROSKEDASTIC unit
  means -> weight each unit by inverse summary-variance (sd_pt^2 / n_cells, the trajectory analogue of de_pb's
  voomWQW quality weights); OLS + a replicate-PERMUTATION null as sensitivity. NOT "exact"/"nominal-FDR"
  unconditionally — the unit-level aggregation fixes the unit, the weighting + permutation guard the rest.

- DECOMPOSITION = Kitagawa/Oaxaca shift-share, observed mean_pt_r = sum_s pi_sr * mu_sr (unit r, state s).
  EXACT 3-channel identity (a 2-channel split is NOT exact — it drops a cross term): mean_pt_r =
  composition_cf_r + progression_cf_r + cross_r - const, where composition_cf_r = sum_s pi_sr * mu_bar_s (vary
  pi, hold mu at pooled global), progression_cf_r = sum_s pi_bar_s * mu_sr (vary mu, hold pi), cross_r =
  sum_s (pi_sr - pi_bar_s)(mu_sr - mu_bar_s), const = sum_s pi_bar_s * mu_bar_s. The 2x2 interaction contrast
  L() is LINEAR -> L(mean_pt) = L(composition_cf) + L(progression_cf) + L(cross) EXACTLY (const drops; L
  annihilates any unit-constant) -> report ALL THREE channels, ASSERT reconstruction == L(mean_pt), define
  loading % against the full 3-channel sum. Holds on the RAW ADDITIVE pt scale (pt_raw) ONLY — logit/asin break
  additivity -> Kitagawa + pt-scale fits on pt_raw. + within-state means (within_homeostatic = the CLEANEST
  "not composition": homeostatic-labelled cells advancing cannot be an abundance artefact). Pooled
  pi_bar/mu_bar are CELL-weighted -> inference is conditional-on-anchors; report a replicate-balanced-anchor
  sensitivity. v1 POINT estimate loaded mostly PROGRESSION (~94%); its progression interaction was SUGGESTIVE
  (FDR~0.077, in the 0.05-0.10 band), NOT significant at 0.05 -> re-derive on R4.6, RECONCILE the Kitagawa
  composition channel against P1's SIG DAM-fraction interaction, no pre-claim either way.

- SENSITIVITY = `glmmTMB` beta GLMM on PER-CELL pseudotime (CRAN; TMB = C++ template, NOT Stan -> ON-lock) — a
  replication-aware per-cell confirmation that models the full bounded (possibly bimodal) distribution the
  summary collapses. Beta needs OPEN (0,1) -> squeeze pt01 = (pt*(n-1)+0.5)/n (Smithson-Verkuilen), or `ordbeta`
  for exact 0/1 mass at root/tip. Model `pt01 ~ tau*amyloid + batch + (1|genotype_batch)` (FIXED batch,
  de_pb-consistent); Wald/LRT on `tau:amyloid`. Only 16 clusters -> asymptotics WEAK -> SUPPORTIVE, not
  load-bearing (parametric/cluster bootstrap if leaned on). Singular fit (interaction var ~0) -> fallback
  rank-normal LMM; degrade gracefully (the limma-summary primary stands alone).

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
  Target `microglia_trajectory` = COMPACT per-cell frame {cell, genotype_batch, genotype, substate,
  on_lineage[H/D-membership flag], pt_raw, pt01[Smithson-Verkuilen squeeze], score_axis_pt, DAM_UCell,
  Homeostatic_UCell} + lineage/provenance lists (~small; NEVER the 612MB object — cheap-render invariant).
  Unit-test pure helpers on a synthetic 2-cluster embedding fixture (lineage present, pseudotime monotone,
  aggregation shapes). LIVE smoke on real microglia_annotated before commit.
  ACCEPTANCE: single clean homeostatic->DAM lineage; pseudotime monotone in DAM_UCell; slingshot-vs-score-axis
  rho recorded (concordance OR honest flag); dims {10,15,20} + all-retained sensitivity recorded (PRIMARY = 15);
  per-unit on-lineage/omitted fraction recorded; gate green.

S2 — Progression interaction + decomposition target.
  `pseudotime_per_replicate` (16-unit summaries on pt_raw: mean_pt, median_pt, q90 leading-edge, within_<state>
  means, + per-unit n_cells & sd for inverse-variance weights; frac_past = fraction of cells past DAM-onset
  [pooled DAM-cell median pt_raw], a bounded [0,1] measure) | `fit_trajectory_contrasts` (reuse
  `factorial_design`'s design + 5 contrasts; FIT each measure by inverse-variance-weighted `limma::lmFit` on
  the measures x 16 matrix WITHOUT eBayes [no strength-borrowing across HETEROGENEOUS measures; per-measure
  residual variance, 9 df] -> `contrasts.fit` -> ordinary t + CI; do NOT route through `fit_limma_log` [that is
  limma-TREND + eBayes for MANY log-intensity features, incoherent for a handful of pseudotime endpoints];
  pt-scale measures on pt_raw UNtransformed, bounded [0,1] measures (frac_past, within_state) logit/asin in
  SEPARATE fits) | `decompose_progression_vs_composition` (3 Kitagawa channels comp/prog/cross + reconstruction
  assert) | glmmTMB per-cell sensitivity. WITHIN-STATE FLOOR: mirror P1 `run_pb_de_substate` min-cells/unit ->
  SKIP within_<state> for any state below floor in any unit (no NA/noisy means); ALWAYS store the state x unit
  count + on-lineage table. Target `trajectory_progression` = contrast tables + 3 decomposition channels +
  glmmTMB effect/CI/p + provenance.
  PRE-REGISTER (state BEFORE fitting): PRIMARY progression endpoint = progression_cf (Kitagawa) +
  within_homeostatic (composition-robust); frac_past = the interpretable bridge; mean_pt reported but FLAGGED
  composition-conflated. MULTIPLICITY family = BH across the 2 PRIMARY interaction tests; all other
  measures/contrasts EXPLORATORY with a separately-labelled FDR. Primary FDR 0.05, REPORT the 0.05-0.10
  SUGGESTIVE band + effect + CI (v1's progression interaction sat there ~0.077; R4.6 re-baseline may move it).
  ACCEPTANCE (outcome-INDEPENDENT): interaction contrast computed on every measure; 3-channel decomposition
  reconstructs L(mean_pt) (assert passes); glmmTMB FIT COMPLETES with `tau:amyloid` effect/CI/p RECORDED
  (concordance AND discordance both FLAGGED, neither required for success); divergence from v1 reported
  honestly; gate green.

S3 — Report + integration.
  `_trajectory.qmd` (trajectory UMAP x pseudotime; pseudotime density/ridge by genotype; interaction FOREST
  [contrasts x measures]; composition-vs-progression-vs-cross decomposition bar; headline progression_cf/frac_past
  interaction; activation-ordering + rooting + position-not-rate + lineage-conditioning(omitted-fraction) +
  transcriptionally-close-substates caveats) + compact `trajectory_report_data` extractor (slim frame + the
  contrast/decomposition tables; cheap-render) + wire into index.qmd AFTER _microglia.qmd + update
  _microglia.qmd's "P2 pointer" -> built. Prose INLINE-COMPUTED from the loaded targets (never hardcoded —
  tracks the cached build). ACCEPTANCE: report renders 0-warning under warn=2; headline + 3-channel
  decomposition + caveats present; gate green. -> then CLOSE-OUT.

## Reproducibility / risks
- rproject.toml: `slingshot` (BioCsoft; auto-pulls SingleCellExperiment / princurve / TrajectoryUtils — pure-R
  Bioc) + `glmmTMB` (CRAN) CONFIRMED present in the 2026-06-22 snapshot (slingshot 2.20.0, glmmTMB 1.1.14;
  codex-verified) -> co-pin glmmTMB<->Matrix<->TMB ABI (source-compile from snapshot -> no binary-mismatch warning).
- slingshot determinism: principal-curve is deterministic given embedding + labels (+ seed); the harmony
  embedding is the stable P1 cached input. Record seed + provenance (S1 mirror of reprocess_provenance).
- gate cheap-render: BOTH trajectory targets COMPACT; no qmd tar_loads the 612MB object.
- glmmTMB can degrade (singular fit) like the sccomp arm -> the locked limma-summary + decomposition is the
  standalone primary; the sensitivity is supportive, not load-bearing.

## DECIDED (gate 2026-06-30) — LEAN ON-LOCK
slingshot + UCell score-axis concordance anchor; weighted limma per-replicate-summary interaction (no eBayes) +
3-channel Kitagawa decomposition + glmmTMB per-cell sensitivity. ENTIRELY pure-R from the pinned snapshot — no
Stan/Python/GitHub, no destiny DPT, no off-lock Lamian arm. The primary is replication-correct at the 16-unit
level (weighted summary test); the two cross-checks are the score-axis (shared-marker CONCORDANCE, not
statistically independent) + glmmTMB (per-cell replication-aware, supportive at 16 clusters). (Rejected: +destiny
DPT triangulation, +off-lock Lamian XCD — both weighed against the lean-rebuild ethos.)
