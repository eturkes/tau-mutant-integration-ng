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
`microglia_trajectory` + `trajectory_progression` + `trajectory_glmm_sensitivity` targets + `_trajectory.qmd`
report section.

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

LOCK POSTURE: the ENTIRE stack is on-lock from the P3M 2026-06-22 (Bioc 3.23 / CRAN) snapshot — the S2 primary
is pure-R (no new dep), the S3 glmmTMB arm source-compiles TMB C++ — NO Stan,
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
Each step spec is SELF-CONTAINED (function contracts inline). Resuming mid-plan: read ONLY your step + the
files it names (R/design.R factorial_design, R/trajectory.R tail, tests/helpers.R) — the Scope/Stack/DECIDED
rationale above is read-once, skip it. Implement the contracts as written; they encode decisions already made.

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

S2 split into S2 (pure-R load-bearing primary, NO new dep) + S3 (glmmTMB supportive arm, isolates the source-
compile + ABI verify + fresh-build cost) — TWO INDEPENDENT targets off `microglia_trajectory`, each closeable in
one window. S3's per-cell glmmTMB does NOT depend on S2's per-replicate summary -> clean decouple, no cross-step
edits. The weighted-limma summary + Kitagawa decomposition is the standalone primary; glmmTMB is supportive.

S2 — Progression interaction + decomposition target (pure-R, on-lock, NO new dependency).
  Append to `R/trajectory.R` (after build_activation_trajectory). All non-base calls namespace-qualified
  (stats::, limma::). df.residual = 16 - 7 = 9. Contracts:
  - `derive_batch(genotype_batch, genotype)`: batch = sub(paste0("^", genotype, "_"), "", genotype_batch); assert
    nzchar(batch) AND identical(paste(genotype, batch, sep = "_"), genotype_batch) (round-trip, FAIL-LOUD --
    genotypes contain "_" so a naive strsplit is wrong). S1's cell_frame omits batch; deriving it here avoids
    re-touching the built microglia_trajectory target. Reused by pseudotime_per_replicate + glmmtmb_pt_sensitivity.
  - `pseudotime_per_replicate(cell_frame, lineage_states, dam_state = "DAM", min_within = 10L)`: filter to
    is.finite(pt_raw) (on-lineage); assert substate %in% lineage_states + validate_trajectory_units(unit,geno);
    batch = derive_batch(genotype_batch, genotype) per unit. dam_onset = stats::median(pt_raw[substate==dam_state])
    PRE-DECLARED. per_unit cols {genotype_batch, genotype, batch, n_cells, sd_pt, mean_pt, median_pt, q90,
    frac_past}; frac_past = mean(pt_raw > dam_onset) = the ONLY
    genuinely [0,1] measure. within_<state> = per-unit MEAN pt_raw within each state, on the pt-SCALE and
    UNtransformed (a within-state mean position is NOT in [0,1]; keeping it on pt_raw puts within_homeostatic in
    the SAME additive units as the Kitagawa progression_cf -> the two PRIMARY endpoints share a scale). State x
    unit count matrix `cnt` (assert all >= 1); pi_mat = cnt col-normalised; mu_mat = per-state per-unit mean;
    sd_mat = per-state per-unit SD of pt_raw (parallel to mu_mat -> within-state inverse-variance weights);
    pi_bar = rowSums(cnt)/sum(cnt) (CELL-weighted pooled composition); mu_bar = per-state pooled mean;
    within_skip[state] = any(cnt[state, ] < min_within). Returns list(per_unit, states, units, counts = cnt, pi,
    mu, sd = sd_mat, pi_bar, mu_bar, dam_onset, within_skip, min_within).
  - `ordinary_t_table(fit, contrast_name, conf_level = 0.95)`: ordinary t from a contrasts.fit'd limma fit
    WITHOUT eBayes (topTable REQUIRES eBayes -> compute by hand): coef = fit$coefficients[,c]; se =
    fit$sigma * fit$stdev.unscaled[,c]; df = fit$df.residual; t = coef/se; p = 2*stats::pt(-abs(t), df); CI =
    coef -/+ stats::qt(1-(1-conf)/2, df)*se. Returns df {measure, contrast, coef, se, t, df, p_value, ci_l, ci_r}.
  - `fit_trajectory_contrasts(measure_mat, design, contrasts, weights = NULL, conf_level = 0.95)`: assert
    colnames(measure_mat)==rownames(design), rownames(contrasts)==colnames(design), full-rank, all finite;
    weights (if given) a matrix matching dim(measure_mat), finite & > 0. limma::lmFit(measure_mat, design,
    weights) -> limma::contrasts.fit -> per-contrast ordinary_t_table. Returns list(fit, top). The `interaction`
    contrast == the `tau_nlgf` coefficient (assert in test vs manual OLS). Do NOT route through `fit_limma_log`
    (that = limma-TREND + eBayes for MANY log-intensity features, incoherent for a handful of heterogeneous
    pseudotime endpoints).
  - `kitagawa_channels(pi_mat, mu_mat, pi_bar, mu_bar, tol = 1e-8)`: const = sum(pi_bar*mu_bar); mean_pt =
    colSums(pi_mat*mu_mat); comp_cf = colSums(pi_mat*mu_bar); prog_cf = colSums(pi_bar*mu_mat); cross =
    colSums((pi_mat-pi_bar)*(mu_mat-mu_bar)); assert max|mean_pt-(comp_cf+prog_cf+cross-const)| < tol. Returns
    df per unit {genotype_batch, mean_pt, comp_cf, prog_cf, cross, const}.
  - `decompose_progression_vs_composition(per_rep, design, contrasts, weights = NULL, conf_level = 0.95)`:
    build the 4-channel matrix (rows mean_pt, comp_cf, prog_cf, cross; cols = units) via kitagawa_channels ->
    fit_trajectory_contrasts with ONE shared per-unit weight vector replicated across all 4 channel-rows
    (matrix(weights, nrow, ncol, byrow = TRUE)) -> exact reconstruction preserved AND progression_cf weighted.
    L_int = fit$fit$coefficients[, "interaction"] (a length-4 NAMED vector over the channel-rows); loadings =
    L_int[c("comp_cf","prog_cf","cross")] / L_int["mean_pt"] (NA loadings if abs(L_int["mean_pt"]) < tol). The
    interaction contrast puts ZERO weight on the intercept, so it annihilates any unit-constant (treatment coding
    absorbs the constant entirely into the intercept); ONE shared per-unit weight vector makes the fitted contrast
    identical across channel-rows -> assert recon_resid = |L_int["mean_pt"] - (L_int["comp_cf"]+L_int["prog_cf"]+
    L_int["cross"])| < tol. Cell-weighted anchors = primary; rowMeans(pi)/rowMeans(mu) replicate-balanced anchors =
    sensitivity. Returns list(channels, fit, L_int, loadings, interaction = fit$top$interaction, recon_resid_max,
    balanced).
  - `freedman_lane_interaction(y, design, int_col = "tau_nlgf", weights = NULL, n_perm = 2000L, seed = 42L)`:
    WLS = OLS on weight-scaled data: r = sqrt(weights, or 1 if NULL); yw = r*y, Xw = r*design; t_obs =
    interaction-coef/SE from lm.fit(Xw, yw) (resid df = n - ncol(design)). REDUCED Xw0 = Xw minus int_col ->
    weighted fitted fw0 + residuals ew0 (homoscedastic on the WEIGHTED scale -> exchangeable; permuting RAW
    residuals would NOT be). Each perm: yw* = fw0 + ew0[sample(n)] -> t* = interaction-t from lm.fit(Xw, yw*).
    RNG-PURE (save + restore .Random.seed & RNGkind via on.exit; set.seed(seed) inside). perm_p =
    (1 + sum(|t*| >= |t_obs|)) / (n_perm + 1). Returns list(t_obs, n_perm, perm_p).
  - `run_trajectory_progression(microglia_trajectory, min_within = 10L, n_perm = 2000L, seed = 42L)`: per_rep =
    pseudotime_per_replicate(cell_frame, lineage_states from provenance$lineage_substates) -> meta =
    per_unit[, c("genotype","batch")] (rownames = genotype_batch; assert 4x4 balanced) -> factorial_design(meta)
    (batch now present -> 9 resid df). w_overall = n_cells/sd_pt^2. DIRECT measures {mean_pt, median_pt, q90,
    within_<used>} (used = states NOT in within_skip) -> M = t(per_unit[, direct_rows]); per-endpoint weight
    matrix W (rows = measures, cols = units): mean_pt = n_cells/sd_pt^2; within_<state> =
    cnt[state,]/per_rep$sd[state,]^2 (state-specific inverse-variance); median_pt & q90 = n_cells/sd_pt^2
    (overall-precision proxy -- quantile sampling variance unmodelled, EXPLORATORY) -> weighted fit + an OLS
    (weights = NULL) sensitivity. BOUNDED frac_past -> rbind(frac_past_logit = log((x+0.5)/(n_cells-x+0.5)),
    x = round(frac_past*n_cells); frac_past_asin = asin(sqrt(frac_past))); weights = n_cells -> bounded fit. decompose_
    (weights = w_overall). Freedman-Lane on {progression_cf, within_homeostatic (iff used), frac_past_logit,
    mean_pt}. PRE-REGISTERED primary family = BH across {progression_cf, within_homeostatic} (drop
    within_homeostatic if skipped); exploratory = separate BH. provenance (versions, v1_progression_loading ~
    0.94, v1_progression_fdr ~ 0.077). Postcondition stopifnot: interaction on EVERY measure; recon_resid_max <
    tol; primary BH present. Returns list(per_unit, counts, dam_onset, within_skip, design,
    contrasts = {weighted, ols, bounded}, decomposition, permutation, primary_family, exploratory_family,
    provenance). (NO glmmTMB here -> S3.)
  Target `trajectory_progression` = run_trajectory_progression(microglia_trajectory) (reads the COMPACT S1
  target; pure-R). WITHIN-STATE FLOOR mirrors P1 run_pb_de_substate.
  FIXTURE (ADD to tests/helpers.R -- the prior WIP was reverted; deterministic, NO RNG): `make_trajectory_cell_frame(
  per_state = 6L, adv = c(MAPTKI=0, P301S=0.2, NLGF_MAPTKI=1.0, NLGF_P301S=1.6), dam_extra = 0L)` -> 4x4 genotype x
  batch grid; genotype_batch = paste(genotype, sprintf("batch%02d", 1:4), sep = "_") (so derive_batch round-trips);
  each unit = per_state Homeostatic + (per_state + dam_extra*(genotype=="NLGF_P301S")) DAM cells -- the extra DAM
  mass lands ONLY in the double-mutant interaction cell, a genuine tau:amyloid COMPOSITION interaction (NOT an
  amyloid main effect). pt_raw = 1+adv[g]+ramp (Homeostatic) / 4+adv[g]+ramp (DAM), ramp = (seq_len(n)/n)*0.3;
  pt01 = squeeze_unit_interval(pt_raw); on_lineage = TRUE. DEFAULT adv encodes a PURE within-state interaction =
  (1.6-0.2)-(1.0-0) = 0.4 with dam_extra=0 -> CONSTANT composition -> ~100% progression loading; FLAT adv (all 0)
  + dam_extra>0 -> PURE composition interaction (prog loading ~0, comp ~1). Cols {cell, genotype_batch, genotype,
  substate, on_lineage, pt_raw, pt01} (NO batch col -> exercises derive_batch, matching the real cell_frame).
  Tests on it: per-rep shapes + floor skip (min_within=20, per_state=8 -> all within_skip
  TRUE); ordinary_t/fit_trajectory_contrasts vs manual OLS (t matches; interaction == tau_nlgf coef); kitagawa
  identity + pure-composition / pure-progression; decomposition reconstruction (recon_resid_max < 1e-8) +
  loadings on the pure-progression fixture (DEFAULT adv, per_state=8: prog loading ~ 1, comp ~ 0, interaction
  coef ~ 0.4) AND the pure-composition fixture (FLAT adv, dam_extra>0: comp loading ~ 1, prog ~ 0); derive_batch
  round-trip + fail-loud on a corrupted genotype_batch;
  Freedman-Lane signal (perm_p < 0.05), null (design-orthogonal residual -> perm_p > 0.9), determinism, RNG-
  purity (runif sequence unchanged across a call).
  PRE-REGISTER (BEFORE fitting): PRIMARY = progression_cf (Kitagawa) + within_homeostatic (composition-robust);
  frac_past = the interpretable bridge; mean_pt FLAGGED composition-conflated. Family = BH across the 2 PRIMARY;
  rest EXPLORATORY (separate FDR). Primary FDR 0.05, REPORT the 0.05-0.10 SUGGESTIVE band + effect + CI (v1
  ~0.077; R4.6 re-baseline may move it).
  ACCEPTANCE (outcome-INDEPENDENT): interaction computed on every measure; 3-channel decomposition reconstructs
  L(mean_pt) (assert passes); primary BH family + Freedman-Lane perm_p RECORDED; v1 divergence reported
  honestly; gate green.

S3 — glmmTMB per-cell sensitivity arm (supportive; isolates the new-dependency + fresh-build cost).
  Add `glmmTMB` to rproject.toml (CRAN; TMB = C++ template, NOT Stan -> ON-lock; P3M co-pins
  glmmTMB<->Matrix<->TMB so the source build loads ABI-clean) -> `rv sync` -> verify it library()-loads under
  options(warn=2) with NO binary-mismatch warning. Append `glmmtmb_pt_sensitivity(cell_frame)` to
  `R/trajectory.R`: on-lineage cells (finite pt01); tau/amyloid as integer 0/1 from genotype; batch =
  derive_batch(genotype_batch, genotype), unit = genotype_batch, both as factors. capture() = withCallingHandlers
  MUFFLING + RECORDING both warnings AND messages (the sccomp lesson:
  optimisers report convergence health via message(), not warning() -> the gate's warn=2 + tar_meta would else
  red). Fit pt01 ~ tau*amyloid + batch + (1|unit), glmmTMB::beta_family(); cond = summary(fit)$coefficients$cond; term = intersect(c("tau:amyloid",
  "amyloid:tau"), rownames(cond)) with assert length(term)==1L (interaction row present + unambiguous); est/se/z/p
  from cond[term, ], CI = est -/+ stats::qnorm(0.975)*se, re_sd from VarCorr, singular = re_sd < 1e-4. DEGRADE to
  the rank-normal LMM if the fit is NULL (try-error) OR !fit$sdr$pdHess OR fit$fit$convergence != 0 OR any
  non-finite est/se OR singular. RANK-NORMAL LMM (on-lock, SAME package): rn = stats::qnorm((rank(pt01)-0.5)/
  length(pt01)); glmmTMB::glmmTMB(rn ~ tau*amyloid + batch + (1|unit), family = stats::gaussian()); extract the
  same interaction row; method = "lmm_ranknorm" (else "glmmTMB_beta"). The limma-summary primary stands alone --
  the arm degrades gracefully, never blocks. Returns list(method, term, estimate, se,
  z, p_value, ci_l, ci_r, re_sd, singular, n_cells, warnings, messages). Target `trajectory_glmm_sensitivity` =
  glmmtmb_pt_sensitivity(microglia_trajectory$cell_frame) -- reads the COMPACT per-cell frame, INDEPENDENT of
  trajectory_progression. Smoke-test on the synthetic fixture; then a FRESH `tar_make` of
  trajectory_glmm_sensitivity on the real microglia_trajectory (the gate force-invalidates ONLY `report`, so
  build this target explicitly) to confirm the glmmTMB arm is gate-clean (warnings captured; no tar_meta/render-
  log red). Test: returns a recorded tau:amyloid effect on the synthetic frame; the degrade path is exercised.
  ACCEPTANCE (outcome-INDEPENDENT): glmmTMB (or the degraded LMM) COMPLETES with `tau:amyloid` effect/CI/p
  RECORDED + FLAGGED supportive (concordance AND discordance both fine, neither required); no gate red from the
  captured warnings; gate green.

S4 — Report + integration.
  `_trajectory.qmd` (trajectory UMAP x pseudotime; pseudotime density/ridge by genotype; interaction FOREST
  [contrasts x measures]; composition-vs-progression-vs-cross decomposition bar; headline progression_cf/frac_past
  interaction; glmmTMB tau:amyloid effect/CI annotated SUPPORTIVE; activation-ordering + rooting +
  position-not-rate + lineage-conditioning(omitted-fraction) + transcriptionally-close-substates caveats) +
  compact `trajectory_report_data` extractor (slim frame + the contrast/decomposition tables + the glmmTMB
  effect row; cheap-render) loading BOTH `trajectory_progression` and `trajectory_glmm_sensitivity` + wire into
  index.qmd AFTER _microglia.qmd + update _microglia.qmd's "P2 pointer" -> built. Prose INLINE-COMPUTED from the
  loaded targets (never hardcoded — tracks the cached build). ACCEPTANCE: report renders 0-warning under
  warn=2; headline + 3-channel decomposition + glmmTMB supportive + caveats present; gate green. -> then
  CLOSE-OUT.

## Reproducibility / risks
- rproject.toml: `slingshot` (BioCsoft; auto-pulls SingleCellExperiment / princurve / TrajectoryUtils — pure-R
  Bioc) + `glmmTMB` (CRAN) CONFIRMED present in the 2026-06-22 snapshot (slingshot 2.20.0, glmmTMB 1.1.14;
  codex-verified) -> co-pin glmmTMB<->Matrix<->TMB ABI (source-compile from snapshot -> no binary-mismatch warning).
- slingshot determinism: principal-curve is deterministic given embedding + labels (+ seed); the harmony
  embedding is the stable P1 cached input. Record seed + provenance (S1 mirror of reprocess_provenance).
- gate cheap-render: ALL THREE trajectory targets COMPACT (trajectory_glmm_sensitivity reads the compact
  per-cell frame, not the 612MB object); no qmd tar_loads the 612MB object.
- glmmTMB can degrade (singular fit) like the sccomp arm -> the locked limma-summary + decomposition is the
  standalone primary; the sensitivity is supportive, not load-bearing.

## DECIDED (gate 2026-06-30) — LEAN ON-LOCK
slingshot + UCell score-axis concordance anchor; weighted limma per-replicate-summary interaction (no eBayes) +
3-channel Kitagawa decomposition + glmmTMB per-cell sensitivity. ENTIRELY on-lock from the pinned snapshot (the
S2 primary is pure-R / no new dep; S3 adds source-compiled glmmTMB/TMB) — no Stan/Python/GitHub, no destiny DPT,
no off-lock Lamian arm. The primary is replication-correct at the 16-unit
level (weighted summary test); the two cross-checks are the score-axis (shared-marker CONCORDANCE, not
statistically independent) + glmmTMB (per-cell replication-aware, supportive at 16 clusters). (Rejected: +destiny
DPT triangulation, +off-lock Lamian XCD — both weighed against the lean-rebuild ethos.)
