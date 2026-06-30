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

S2 (the pure-R primary) is SPLIT AGAIN -> S2a + S2b: even the no-new-dep primary (8 fns + fixture + ~6 test
groups + live smoke + gate) overflowed one window TWICE. Split at the VERIFICATION SEAM -- S2a = the estimation
core (per-replicate summary + contrast fit + Kitagawa decomposition; 6 fns, fully verifiable on the deterministic
fixture, NO target) -> S2b = inference + orchestration (Freedman-Lane + run_trajectory_progression + the
`trajectory_progression` target + live smoke; the live-integration half). S3 (glmmTMB supportive arm) stays a
SEPARATE INDEPENDENT target off `microglia_trajectory` (source-compile + ABI verify + fresh build). The
weighted-limma summary + Kitagawa decomposition (S2a math, S2b inference) is the standalone primary; glmmTMB =
supportive. Each of S2a / S2b / S3 closeable in one window.

S2a — Per-replicate summary + contrast fit + Kitagawa decomposition (estimation core; pure-R, NO new dep, NO target yet).
  RESUMING? Read ONLY this S2a block + R/design.R `factorial_design` + the tail of R/trajectory.R (after
  build_activation_trajectory) + tests/helpers.R + tests/test_trajectory.R. SKIP Scope/Stack/other steps.
  Append the 6 estimation fns to `R/trajectory.R` (after build_activation_trajectory). All non-base calls
  namespace-qualified (stats::, limma::). df.residual = 16 - 7 = 9. Contracts:
  - `derive_batch(genotype_batch, genotype)` (VECTORISED; base sub() does NOT vectorise over a length>1 pattern ->
    warns + uses pattern[1] only -> mis-extracts under warn=2, so use literal-prefix string ops, NOT regex):
    prefix = paste0(as.character(genotype), "_"); assert all(startsWith(genotype_batch, prefix)); batch =
    substring(genotype_batch, nchar(prefix) + 1L); assert all(nzchar(batch)) AND identical(paste(genotype, batch,
    sep = "_"), genotype_batch) (round-trip, FAIL-LOUD -- genotypes contain "_" so a naive strsplit is wrong; a
    literal prefix also dodges regex-metachar pitfalls). S1's cell_frame omits batch; deriving it here avoids
    re-touching the built microglia_trajectory target. Reused by pseudotime_per_replicate + glmmtmb_pt_sensitivity.
  - `pseudotime_per_replicate(cell_frame, lineage_states, dam_state = "DAM", min_within = 10L)`: filter to
    is.finite(pt_raw) (on-lineage); assert substate %in% lineage_states + validate_trajectory_units(unit,geno);
    batch = derive_batch(genotype_batch, genotype) per unit. dam_onset = stats::median(pt_raw[substate==dam_state])
    PRE-DECLARED. per_unit cols {genotype_batch, genotype, batch, n_cells, sd_pt, mean_pt, median_pt, q90,
    frac_past}; frac_past = mean(pt_raw > dam_onset) = the ONLY
    genuinely [0,1] measure. within_<lc> (lc = tolower(state) -> within_homeostatic / within_dam; ONE sanitizer, S2b
    reuses it -- NEVER the verbatim-cased state) = per-unit MEAN pt_raw within each state, on the pt-SCALE and
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
    weights (if given) a matrix, identical(dimnames(weights), dimnames(measure_mat)) (limma consumes weights BY
    POSITION -> assert dimnames match, NOT just dims, to catch row/unit drift), finite & > 0. limma::lmFit(measure_mat, design,
    weights) -> cfit = limma::contrasts.fit(...) -> per-contrast ordinary_t_table. Returns list(fit = cfit, top)
    where top = setNames(lapply(colnames(contrasts), function(cn) ordinary_t_table(cfit, cn, conf_level)),
    colnames(contrasts)) -- a NAMED list keyed by contrast (top$interaction = the per-measure interaction table) so
    fit$top$interaction AND fit$fit$coefficients[,"interaction"] both resolve downstream. The `interaction`
    contrast == the `tau_nlgf` coefficient (assert in test vs manual OLS). Do NOT route through `fit_limma_log`
    (that = limma-TREND + eBayes for MANY log-intensity features, incoherent for a handful of heterogeneous
    pseudotime endpoints).
  - `kitagawa_channels(pi_mat, mu_mat, pi_bar, mu_bar, tol = 1e-8)`: assert ALIGNMENT first (every broadcast below is
    POSITION-based -> a silent wrong split on state-order drift): identical(dimnames(pi_mat), dimnames(mu_mat)),
    identical(rownames(pi_mat), names(pi_bar)), identical(names(pi_bar), names(mu_bar)). const = sum(pi_bar*mu_bar);
    mean_pt = colSums(pi_mat*mu_mat); comp_cf = colSums(pi_mat*mu_bar); prog_cf = colSums(pi_bar*mu_mat); cross =
    colSums((pi_mat-pi_bar)*(mu_mat-mu_bar)) (pi_bar/mu_bar recycle DOWN columns -> need the rownames-aligned state
    order asserted above); assert max|mean_pt-(comp_cf+prog_cf+cross-const)| < tol. Returns
    df per unit {genotype_batch, mean_pt, comp_cf, prog_cf, cross, const}.
  - `decompose_progression_vs_composition(per_rep, design, contrasts, weights = NULL, conf_level = 0.95)`:
    build the 4-channel matrix M4 (rows mean_pt, comp_cf, prog_cf, cross; cols = units, dimnames pinned) via
    kitagawa_channels -> fit_trajectory_contrasts with ONE shared per-unit weight vector (weights = a per-unit vector
    NAMED by units) replicated across all 4 channel-rows: W4 = matrix(weights[colnames(M4)], nrow(M4), ncol(M4),
    byrow = TRUE, dimnames = dimnames(M4)) (matrix() DROPS dimnames + index by colnames -> else
    fit_trajectory_contrasts' identical(dimnames(weights), dimnames(M4)) assert fails) -> exact reconstruction preserved AND progression_cf weighted.
    L_int = fit$fit$coefficients[, "interaction"] (a length-4 NAMED vector over the channel-rows); loadings =
    L_int[c("comp_cf","prog_cf","cross")] / L_int["mean_pt"] (NA loadings if abs(L_int["mean_pt"]) < tol). The
    interaction contrast puts ZERO weight on the intercept, so it annihilates any unit-constant (treatment coding
    absorbs the constant entirely into the intercept); ONE shared per-unit weight vector makes the fitted contrast
    identical across channel-rows -> assert recon_resid = |L_int["mean_pt"] - (L_int["comp_cf"]+L_int["prog_cf"]+
    L_int["cross"])| < tol. Cell-weighted anchors = primary; rowMeans(pi)/rowMeans(mu) replicate-balanced anchors =
    sensitivity. Returns list(channels, fit, L_int, loadings, interaction = fit$top$interaction, recon_resid_max,
    balanced).
  FIXTURE (ADD to tests/helpers.R -- the prior WIP was reverted; deterministic, NO RNG): `make_trajectory_cell_frame(
  per_state = 6L, adv = c(MAPTKI=0, P301S=0.2, NLGF_MAPTKI=1.0, NLGF_P301S=1.6), dam_extra = 0L)` -> 4x4 genotype x
  batch grid; genotype_batch = paste(genotype, sprintf("batch%02d", 1:4), sep = "_") (so derive_batch round-trips);
  each unit = per_state Homeostatic + (per_state + dam_extra*(genotype=="NLGF_P301S")) DAM cells -- the extra DAM
  mass lands ONLY in the double-mutant interaction cell, a genuine tau:amyloid COMPOSITION interaction (NOT an
  amyloid main effect). pt_raw = 1+adv[g]+ramp (Homeostatic) / 4+adv[g]+ramp (DAM), ramp = ((seq_len(n)-0.5)/n)*0.3
  (MIDPOINT rule -> block mean EXACTLY 0.15 for ANY block size n, so an unequal DAM count does NOT shift the
  within-state mean -> the pure-composition fixture is EXACTLY pure, prog/cross interaction 0 not merely approximate);
  pt01 = squeeze_unit_interval(pt_raw); on_lineage = TRUE. DEFAULT adv encodes a PURE within-state interaction =
  (1.6-0.2)-(1.0-0) = 0.4 with dam_extra=0 -> CONSTANT composition -> 100% progression loading; FLAT adv (all 0)
  + dam_extra>0 -> PURE composition interaction (prog loading 0, comp 1, EXACT under the midpoint ramp). Optional
  `jitter = 0`: when >0 ADD a deterministic NON-additive perturbation ((gi*bi) %% 5)*jitter to pt_raw (gi/bi =
  genotype/batch indices, NO RNG) -> breaks the saturated design's zero residual (sigma > 0) for S2b's STRUCTURAL
  orchestrator test; default 0 keeps the component tests EXACT-pure. Cols {cell,
  genotype_batch, genotype, substate, on_lineage, pt_raw, pt01} (NO batch col -> exercises derive_batch, matching
  the real cell_frame).
  S2a TESTS (tests/test_trajectory.R; all on the deterministic fixture / make_meta16, warn=2):
  - derive_batch round-trip + fail-loud on a corrupted genotype_batch + a VECTOR call (>1 genotype, exercises the
    vectorised path).
  - per-rep shapes + floor skip (min_within=20, per_state=8 -> all within_skip TRUE).
  - ordinary_t / fit_trajectory_contrasts vs manual OLS: interaction coef == tau_nlgf coef AND t == hand-computed
    coef/se. BANKED -- the manual-OLS response on make_meta16 MUST carry a NON-ADDITIVE term (e.g.
    base + ((gi*bi) %% 5) * 0.07, gi/bi = genotype/batch indices) so the SATURATED genotype+batch design leaves
    residual > 0; a batch-ALIGNED wiggle is fully absorbed -> sigma = 0 -> t = Inf/NaN -> the test is vacuous.
  - kitagawa identity (reconstruction < 1e-8) + pure-composition / pure-progression channels.
  - decomposition reconstruction (recon_resid_max < 1e-8) + loadings: pure-progression fixture (DEFAULT adv,
    per_state=8 -> prog loading ~ 1, comp ~ 0, interaction coef ~ 0.4) AND pure-composition fixture (FLAT adv,
    dam_extra>0 -> comp loading 1, prog 0, EXACT to 1e-8 via the midpoint ramp).
  ACCEPTANCE (S2a): the 6 fns implemented + namespace-qualified; fixture added; the 5 test groups pass at warn=2;
  Kitagawa + decomposition reconstruction asserts pass (< 1e-8). NO target wired (S2b wires it). gate green.

S2b — Progression interaction inference + orchestrator + target (pure-R; the live-integration half).
  RESUMING? Read ONLY this S2b block + the S2a fns now in R/trajectory.R + R/design.R `factorial_design` +
  R/de_pb.R `assert_complete_crossing` + _targets.R (the microglia_trajectory target) + tests/test_trajectory.R.
  SKIP Scope/Stack above.
  Append the 2 inference fns to `R/trajectory.R` (after the S2a fns). namespace-qualified (stats::, limma::, utils::).
  - `freedman_lane_interaction(y, design, int_col = "tau_nlgf", weights = NULL, n_perm = 2000L, seed = 42L)`:
    WLS = OLS on weight-scaled data: r = sqrt(weights, or 1 if NULL); yw = r*y, Xw = r*design. interaction-t helper
    int_t(yv) (FULL model; Xw fixed across perms -> precompute XtXinv = chol2inv(chol(crossprod(Xw))) ONCE -- the
    pivot-FREE (X'X)^-1 (sidesteps qr()'s LINPACK column-pivoting, so XtXinv's index order matches lm.fit's coef
    order for the full-rank design); assert all(is.finite(XtXinv))): f = lm.fit(Xw, yv); p = ncol(Xw); df = nrow(Xw) - p; sigma2 =
    sum(f$residuals^2)/df; j = match(int_col, colnames(design)); t = f$coefficients[j] / sqrt(sigma2 * XtXinv[j, j]).
    t_obs = int_t(yw). REDUCED Xw0 = Xw minus int_col -> weighted fitted fw0 + residuals ew0 (APPROXIMATELY
    exchangeable on the WEIGHTED scale -- weights are ESTIMATED from the same unit summaries, so exchangeability is
    conditional/approximate NOT exact -> this permutation null is a SENSITIVITY, not a nominal-exact test; permuting
    RAW unweighted residuals would be worse still). Each perm: yw* = fw0 + ew0[sample(n)] ->
    t* = int_t(yw*). RNG-PURE: on.exit save+restore .Random.seed & RNGkind; set.seed(seed, kind = "Mersenne-Twister",
    normal.kind = "Inversion", sample.kind = "Rejection") inside (pin all 3 kinds, matching S1). perm_p =
    (1 + sum(|t*| >= |t_obs|)) / (n_perm + 1). Returns list(t_obs, n_perm, perm_p).
  - `run_trajectory_progression(microglia_trajectory, min_within = 10L, n_perm = 2000L, seed = 42L)`: per_rep =
    pseudotime_per_replicate(cell_frame, lineage_states from provenance$lineage_substates) -> meta =
    per_unit[, c("genotype_batch","genotype","batch")] with rownames(meta) = per_unit$genotype_batch ->
    assert_complete_crossing(meta, "genotype_batch") (needs genotype_batch + genotype + batch cols = the 4x4
    balance check, fail-loud) -> factorial_design(meta) (reads genotype/batch + the rownames -> design rownames =
    genotype_batch matching M's colnames; batch present -> 9 resid df). w_overall = n_cells/sd_pt^2 (NAMED by
    genotype_batch). DIRECT measures {mean_pt, median_pt, q90,
    within_<used>} (used = states NOT in within_skip) -> M = t(per_unit[, direct_rows]) (dimnames = direct_rows x
    units); per-endpoint weight matrix W built with IDENTICAL dimnames(W) == dimnames(M) (measures x units, same
    order -- limma weights apply BY POSITION; assert before fitting): mean_pt = n_cells/sd_pt^2; within_<lc> =
    cnt[state,]/per_rep$sd[state,]^2 (state-specific inverse-variance, row keyed within_<tolower(state)> to match
    pseudotime_per_replicate); median_pt & q90 = n_cells/sd_pt^2
    (overall-precision proxy -- quantile sampling variance unmodelled, EXPLORATORY) -> weighted fit + an OLS
    (weights = NULL) sensitivity. BOUNDED frac_past -> rbind(frac_past_logit = log((x+0.5)/(n_cells-x+0.5)),
    x = round(frac_past*n_cells); frac_past_asin = asin(sqrt(frac_past))); weights = n_cells (a precision PROXY:
    EXACT-up-to-constant for the asin VST, only APPROXIMATE for the logit whose delta-method inverse-variance is
    n*p*(1-p) -- both EXPLORATORY bridges, not primary) -> bounded fit. decompose_progression_vs_composition
    (weights = w_overall). Freedman-Lane (EACH call weighted to MATCH its limma fit: progression_cf & mean_pt ->
    w_overall, within_<lc> -> cnt[state,]/sd[state,]^2, frac_past_logit -> n_cells) on {progression_cf,
    within_homeostatic (iff used), frac_past_logit, mean_pt}. PRE-REGISTERED primary family = BH across {progression_cf, within_homeostatic} (drop
    within_homeostatic if skipped); exploratory = separate BH. provenance (versions, v1_progression_loading ~
    0.94, v1_progression_fdr ~ 0.077). Postcondition stopifnot: interaction on EVERY measure; recon_resid_max <
    tol; primary BH present. Returns list(per_unit, counts, dam_onset, within_skip, design,
    contrasts = {weighted, ols, bounded}, decomposition, permutation, primary_family, exploratory_family,
    provenance). (NO glmmTMB here -> S3.)
  Target `trajectory_progression` = run_trajectory_progression(microglia_trajectory) -> ADD to _targets.R AFTER
  the microglia_trajectory target: tar_target(trajectory_progression,
  run_trajectory_progression(microglia_trajectory), format = "qs"). Reads the COMPACT S1 target; pure-R.
  WITHIN-STATE FLOOR mirrors P1 run_pb_de_substate.
  PRE-REGISTER (encoded in run_trajectory_progression, BEFORE fitting): PRIMARY = progression_cf (Kitagawa) +
  within_homeostatic (composition-robust); frac_past = the interpretable bridge; mean_pt FLAGGED
  composition-conflated. Family = BH across the 2 PRIMARY; rest EXPLORATORY (separate FDR). Primary FDR 0.05,
  REPORT the 0.05-0.10 SUGGESTIVE band + effect + CI (v1 ~0.077; R4.6 re-baseline may move it).
  S2b TESTS (tests/test_trajectory.R, warn=2):
  - Freedman-Lane on a SEPARATE deterministic design (make_meta16 -> factorial_design): signal (perm_p < 0.05),
    null (perm_p > 0.9), determinism (same seed -> identical perm_p), RNG-purity (a runif draw advances
    IDENTICALLY across a freedman_lane_interaction call -> the on.exit restore). BANKED FL design -- null:
    e = stats::lm.fit(design, v)$residuals (v deterministic) is design-ORTHOGONAL -> tau_nlgf coef ~ 0 ->
    t_obs ~ 0 -> every |t*| >= |t_obs| -> perm_p ~ 1; signal: y = 2*design[,"tau_nlgf"] + 0.1*e -> perm_p < 0.05.
  - BANKED: do NOT assert INFERENTIAL quantities (t / p / perm_p) from run_trajectory_progression on the EXACT-pure
    fixture -- its zero-residual within-genotype design gives sigma = 0 -> t = Inf/NaN on every measure (a limma
    zero-variance warning would also red warn=2). But DO add a STRUCTURAL orchestrator test on a NON-additive
    fixture variant (make_trajectory_cell_frame(jitter > 0) -> sigma > 0): run the FULL orchestrator + assert
    STRUCTURE only -- every measure carries an interaction row; list fields present; recon_resid_max < tol;
    primary_family BH keys present; the within_<lc> columns wire through. This catches return-shape / field-name /
    wiring bugs CHEAPLY in unit tests, BEFORE the context-heavy live smoke. Outcome assertions still come from (a)
    the S2a component tests, (b) the LIVE smoke below, (c) the gate's tar_make building trajectory_progression
    FRESH on real noisy data.
  LIVE SMOKE (before the gate): Rscript tar_read(microglia_trajectory) -> run_trajectory_progression -> inspect
  primary_family (BH finite), decomposition$recon_resid_max (< 1e-8), permutation perm_p (finite), interaction
  present on EVERY measure.
  ACCEPTANCE (S2b, outcome-INDEPENDENT): interaction computed on every measure; the target builds on real data via
  the gate's tar_make; 3-channel decomposition reconstructs L(mean_pt) on real data (recon_resid_max < 1e-8);
  primary BH family + Freedman-Lane perm_p RECORDED; v1 divergence (loading ~0.94, fdr ~0.077) reported honestly
  in provenance; gate green.

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
  same interaction row; method = "lmm_ranknorm" (else "glmmTMB_beta"). HEALTH-CHECK the fallback LMM with the SAME
  battery (NULL / !pdHess / convergence != 0 / non-finite est|se / singular); if BOTH the beta GLMM AND the
  rank-normal LMM fail, method = "failed" with estimate/se/z/p_value/ci_l/ci_r = NA_real_ + the captured
  warnings/messages (RECORD a failed-supportive result, NEVER error). The limma-summary primary stands alone --
  the arm degrades gracefully, never blocks. Returns list(method, term, estimate, se,
  z, p_value, ci_l, ci_r, re_sd, singular, n_cells, warnings, messages). Target `trajectory_glmm_sensitivity` =
  glmmtmb_pt_sensitivity(microglia_trajectory$cell_frame) -- reads the COMPACT per-cell frame, INDEPENDENT of
  trajectory_progression. Smoke-test on the synthetic fixture; then a FRESH `tar_make` of
  trajectory_glmm_sensitivity on the real microglia_trajectory (the gate force-invalidates ONLY `report`, so
  build this target explicitly) to confirm the glmmTMB arm is gate-clean (warnings captured; no tar_meta/render-
  log red). Test: returns a recorded tau:amyloid effect on the synthetic frame; the degrade path is exercised.
  ACCEPTANCE (outcome-INDEPENDENT): glmmTMB (or the degraded LMM, or a RECORDED method="failed" if both degrade)
  COMPLETES with `tau:amyloid` effect/CI/p RECORDED + FLAGGED supportive (concordance AND discordance both fine,
  neither required); no gate red from the captured warnings; gate green.

S4 — Report + integration. SPLIT (2026-07-01) into S4a (data layer) + S4b (report layer): the combined S4
(extractor + target + qmd + wiring + pointers + test + render-debug + gate + docs + commit + close-out) overflowed
one window. WORK PARKED on branch `wip-p2s4-report` (pre-revert): R/trajectory.R (+`trajectory_report_data` extractor
appended after glmmtmb_pt_sensitivity, BUILT-VALIDATED — the target built 674KB, all fields present), _targets.R
(+`trajectory_report` target after trajectory_glmm_sensitivity), _trajectory.qmd (FULL chapter draft, UNRENDERED).
Each sub-step RESTORES its parked files (`git checkout wip-p2s4-report -- <file>`) → near-mechanical, NO re-derive,
NO re-read of R/ model files. DELETE the branch after S4b lands. Headline numbers live in memory.md P2-S2b/S3 and are
INLINE-COMPUTED in the qmd (never hardcoded). Caveat: main HEAD 863cb75 has no trajectory_report code/object (pruned)
— restore from the branch.

S4a — Compact extractor + target + test (data layer; pure-R, NO render).
  RESUMING? Read ONLY this S4a block + R/microglia.R::microglia_report_data (the MODEL extractor) + the tail of
  R/trajectory.R (where the extractor appends) + tests/test_trajectory.R + tests/helpers.R. SKIP Scope/Stack/other steps.
  RESTORE the parked extractor + target: `git checkout wip-p2s4-report -- R/trajectory.R _targets.R`; `git diff` MUST
  show ONLY the appended fn + the one new target (else inspect).
  `trajectory_report_data(microglia_trajectory, trajectory_progression, trajectory_glmm_sensitivity)` = compact
  extractor reading the 3 COMPACT trajectory targets (NEVER the 612MB Seurat — cheap-render invariant). Returns
  list{cell_frame(genotype/substate/on_lineage/pt_raw/score_axis_pt), interaction(rbind primary+exploratory families,
  +family col), weighted_top(5-contrast mean_pt rows), decomposition(loadings/L_int/interaction), per_unit,
  lineage_per_unit, sensitivity, glmm(the S3 row subset), provenance(dims/root/concordance/loadings/v1/n_perm/seed/
  versions)}. Fail-loud stopifnot on every input's structure + finite/consistency guards (MIRROR microglia_report_data).
  Target `trajectory_report` (format="qs") after trajectory_glmm_sensitivity in _targets.R.
  TEST (add to tests/test_trajectory.R): REUSE the microglia_trajectory stub the S2b structural-orchestrator test
  builds (make_trajectory_cell_frame(jitter>0)); extend it with the per_unit / sensitivity / provenance fields the
  extractor reads → run run_trajectory_progression + glmmtmb_pt_sensitivity on it → trajectory_report_data(...) →
  assert the returned list carries every documented field (+ "progression_cf" in interaction$measure, finite
  coef/fdr, comp_cf/prog_cf/cross in decomposition$loadings) AND a fail-loud guard fires on a malformed input
  (drop a required name → expect_error). The extractor is BUILT-VALIDATED → the only NEW work is this test.
  ACCEPTANCE (S4a): extractor + target restored; test passes at warn=2; the gate builds `trajectory_report` FRESH on
  real data (the force-render still renders index.qmd WITHOUT the trajectory chapter — not wired till S4b → the
  unrendered qmd cannot red the gate); map.md += trajectory_report target + trajectory_report_data fn; gate green;
  commit `trajectory (p2 s4a): compact report-data extractor + target + test`.

S4b — Chapter + wiring (report layer; the render-debug half).
  RESUMING? Read ONLY this S4b block + (after restoring) _trajectory.qmd + _microglia.qmd (MODEL chapter + the 2
  "P2 pointer" sentences) + index.qmd (include list + Overview). SKIP Scope/Stack/other steps AND the R/ model files
  (the qmd is parked-complete — do NOT re-derive it).
  RESTORE the parked chapter: `git checkout wip-p2s4-report -- _trajectory.qmd`. COMPLETE draft: title
  `# The tau-amyloid synergy adds DAM cells rather than advancing them {#sec-trajectory}`; warn=2 setup w/ inline
  helpers (irow/mp_ctr/fmt_p + adaptive glmm_sentence + per-genotype conditioning audit); sections (amyloid axis
  shift / composition-not-progression / per-cell glmmTMB / reconciliation+robustness+concordance / caveats+
  provenance); 3 fig chunks (pseudotime-shift, decomposition, concordance), 2 table chunks (interaction-table,
  conditioning), glmmTMB inline, provenance cat. Prose INLINE-COMPUTED from trajectory_report.
  WIRE: index.qmd → add `{{< include _trajectory.qmd >}}` AFTER `{{< include _microglia.qmd >}}`; extend the Overview
  paragraph to name the trajectory chapter (activation-axis advance; synergy = composition not progression).
  POINTER: _microglia.qmd → rewrite the 2 forward-pointer sentences (locate by "rate question handed to the
  trajectory phase" + "the subject of the next phase") to reflect P2 is BUILT + found the synergy COMPOSITIONAL
  (more DAM cells, NO supported further-advance); cross-ref @sec-trajectory. British English + hyphens.
  DEP: the qmd uses patchwork (`(p_a | p_b)`) — CONFIRMED already in rproject.toml + installed + used by
  _microglia.qmd/_qc.qmd (NO new dep). All other infrastructure (theme_tau, scale_fill_genotype, genotype_levels)
  is confirmed present → render-debug should surface only minor typos (a field name, an aes/scale) if anything.
  RENDER-DEBUG: scripts/check.sh now force-renders index.qmd WITH _trajectory.qmd under warn=2 → fix any chunk
  warning/error. Iterate to 0-warning.
  ACCEPTANCE (S4b): report renders 0-warning under warn=2; headline + 3-channel decomposition + glmmTMB supportive +
  the 5 caveats present; gate green. memory.md += a P2-S4 section (cheap-render invariant for _trajectory.qmd + any
  render gotcha); map.md += _trajectory.qmd include. Commit `trajectory (p2 s4b): trajectory chapter + index wiring +
  microglia pointer`. DELETE branch `wip-p2s4-report`. → CLOSE-OUT (next session = session-prompt CLOSE-OUT mode:
  adversarial plan review, fold P2 digest → history.md, archive the plan, reset Active plan).

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
