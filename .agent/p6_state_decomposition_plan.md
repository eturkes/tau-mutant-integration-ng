# P6 - microglial state composition versus regulation

Status = ACTIVE. Next = S3 exact UCell channel decomposition + verdict.

## Question + choice

Is the `NLGF_P301S` phenotype explained by **more canonical DAM cells**, or does
mutant tau also alter transcription **within** DAM and/or Homeostatic microglia?

Chosen default = decompose the microglial 2x2 response into replicate-level state
occupancy and state-conditional regulation. The phase combines:

1. a beta-binomial Homeostatic/DAM occupancy model;
2. raw-count pseudobulk DE within each state plus a direct paired state-response
   contrast; and
3. an exact reference standardization of the five predeclared UCell programmes
   into composition, within-state, and cross channels.

The phase-defining comparison is the canonical interaction
`(NLGF_P301S - P301S) - (NLGF_MAPTKI - MAPTKI)`. Both within-tau amyloid
contrasts, `tau_in_nlgf`, and the remaining canonical contrast stay in the audit.
This asks whether tau changes the amyloid response without assuming the archived
or current answer.

This replaces the unexecuted cross-cell-type P6 at the user's direction. It is a
tighter test of the current story's unresolved boundary: P2 found an interaction
carried by DAM composition, but did not test whether annotated DAM cells also
change their molecular programme.

Alternatives rejected for this phase:

- State-wise pseudobulk alone: valid within-state DE, but it cannot quantify how
  much the aggregate shift is carried by occupancy versus conditional scores.
- Per-cell mixed models, CoCoA-diff, or deep counterfactual models: cell-level
  inference risks pseudoreplication; matching/causal language is unjustified for
  a post-genotype state label in this cross-sectional experiment; opacity buys
  little with 16 complete replicate units.
- Causal mediation: microglial state is affected by genotype and defined from the
  same transcriptome. P6 estimates statistical standardization, not a natural
  indirect effect or same-cell counterfactual.
- IFN/proliferative primary branches: IFN is sparse (minimum five cells/unit) and
  no cell is Proliferative-dominant. They remain measured programmes, never
  fabricated states.
- Broad GO/network discovery: repeats retired report bloat and invites
  outcome-selected storytelling. Five fixed marker programmes carry the
  programme-level verdict; gene rows remain context.
- `dreamlet`/new dependencies: exactly two complete paired states fit the existing
  `edgeR`/`limma`/`glmmTMB` stack. A new modelling layer would add environment and
  interpretation cost without solving a missing design capability.

## Planning-time feasibility snapshot

Observed from the cached live `microglia_annotated` target on 2026-07-14; every
quantity is rechecked fail-loud in S1 rather than trusted as a promise:

- 23,160 retained microglia: DAM = 11,189; Homeostatic = 11,174; IFN = 797;
  Proliferative-dominant = 0.
- Homeostatic + DAM = 22,363/23,160 (96.56%); coverage across the 16 units ranges
  93.93%-97.91%.
- Both primary states occur in all 16 `genotype_batch` units. DAM cells/unit =
  31-2,151; Homeostatic = 52-1,237. The sparsest full unit has 233 retained cells
  (176 DAM, 52 Homeostatic, 5 IFN).
- All five raw UCell score columns are finite and non-constant across the primary
  cells.
- A live two-state raw-count aggregation/filter smoke retained 13,599
  Homeostatic and 9,148 DAM genes, completed warning-free in 7.6 seconds at
  approximately 2.73 GB peak RSS; the in-memory pseudobulk payload was 14.6 MB.
- The locked beta-binomial formula converged warning-free with convergence code
  zero and positive-definite Hessian on the live 16-row count table.
- Required packages are already locked and loadable; P6 changes neither
  `rproject.toml` nor `rv.lock`.

## Research decision

### A. State occupancy - “more DAM”

- Biological unit = `genotype_batch` (16; four/genotype), not cell.
- Primary state universe is fixed before inference to Homeostatic + DAM. Define
  `DAM_fraction = DAM / (DAM + Homeostatic)`; report excluded-state coverage by
  unit beside it.
- Fit `glmmTMB(cbind(DAM, Homeostatic) ~ 0 + genotype + batch,
  family = betabinomial(link = "logit"))`.
- Extract all five canonical log-odds contrasts and batch-standardized
  probability-scale genotype means/contrasts. Probability estimates average
  predictions equally over the four observed batch levels; analytic delta-method
  covariance propagates through every contrast.
- Primary inference = probability-scale effect/95% CI plus zero-null and
  TREAT-style minimum-effect Wald tests. Meaningful margin = 0.10 absolute
  fraction; the threshold statistic uses `(|effect| - 0.10) / SE` with a
  conservative two-sided tail. Zero-null and minimum-effect BH families each span
  the five canonical contrasts. A “meaningful occupancy shift” requires the
  minimum-effect family, not a large-looking point estimate. Log-odds estimates
  remain a complete model-scale audit.
- Small-sample interaction sensitivity = empirical-logit OLS on
  `log((DAM + 0.5) / (Homeostatic + 0.5))` with the same factorial design plus a
  batch-stratified Freedman-Lane residual permutation generalized from the P2
  machinery. It is sensitivity evidence, not a second chance to choose a
  preferred p-value.
- Diagnostics are build-fatal: warning-clean fit, convergence code zero,
  `pdHess = TRUE`, finite dispersion/VCOV/predictions, full-rank model, and no
  extreme/nonidentified coefficients.
- Wording = retained-nuclei state occupancy. Cell fractions are not tissue
  abundance estimates.

### B. State-conditional raw-count response - “different DAM”

- Aggregate raw RNA counts by `subpopulation x genotype_batch` from the cached
  annotated object. SCT, Harmony, UMAP, and UCell scores never enter count DE.
- Hard inclusion = both Homeostatic and DAM have >=20 cells, positive libraries,
  all 16 aligned units, canonical genotype/batch crossing, and a full-rank design.
- Fit each state separately with the existing project pattern:
  `factorial_design()` (`~ tau + nlgf + tau_nlgf + batch`),
  `edgeR::filterByExpr(min.count = 5)`, `edgeR::normLibSizes`,
  `limma::voomWithQualityWeights`, `lmFit`, robust `eBayes`, and all five
  canonical contrasts.
- Primary gene evidence = raw log2FC/95% CI plus `limma::treat(lfc = 0.5)`;
  BH is explicit within state x contrast. A hit count is descriptive only.
- Formal state specificity uses the intersection of state-filtered genes and
  paired unit-level
  `delta = DAM logCPM - Homeostatic logCPM`. Fit the same factorial design to
  `delta`; canonical contrasts therefore equal
  `(DAM response) - (Homeostatic response)`. Primary precision weight is
  `1 / (1 / w_DAM + 1 / w_Homeostatic)` from the two voom fits; an unweighted
  paired fit is the fixed sensitivity. This direct contrast - never
  “significant here, not there” - carries state-selective claims.
- A synthetic fixture must reconstruct every paired contrast from the two
  separate state coefficients to numerical tolerance before live inference.
- Programme inference uses `limma::mroast`/rotation tests on the weighted
  expression objects for exactly the committed `canonical_microglia_markers`:
  Homeostatic, DAM, IFN, Proliferative, MHC_APC. Feature membership follows the
  committed feature-to-symbol map without collapsing count rows. Use 9,999
  rotations and an RNG-pure fixed seed; FDR spans five programmes within endpoint
  x contrast. Broad GO is outside scope.
- Gene rows contextualize programme findings. No lone gene or count of DE genes
  can trigger the phase verdict.
- The UCell and rotation readouts are two summaries of the same RNA assay, not
  independent replication. Concordance is a representation/scale guard only.

### C. Exact UCell standardization - channel attribution

For every primary-state cell, use the five existing raw UCell scores. Divide each
programme by its pooled Homeostatic+DAM cell-level SD; compare programmes only on
their own standardized scale. `total` is therefore the aggregate score inside the
fixed two-state universe, not all retained microglia. Within unit `u`, state `s`,
programme `p`, define
the cell fraction `pi[u,s]` and state-mean score `mu[u,s,p]`. Equal-unit anchors
prevent large libraries from defining the reference:

```text
pi_bar[s]    = mean_u pi[u,s]
mu_bar[s,p]  = mean_u mu[u,s,p]

total[u,p]       = sum_s pi[u,s] * mu[u,s,p]
anchor[p]        = sum_s pi_bar[s] * mu_bar[s,p]
composition[u,p]  = sum_s (pi[u,s] - pi_bar[s]) * mu_bar[s,p]
within_state[u,p] = sum_s pi_bar[s] * (mu[u,s,p] - mu_bar[s,p])
cross[u,p]        = sum_s (pi[u,s] - pi_bar[s]) *
                          (mu[u,s,p] - mu_bar[s,p])
```

The identity
`total - anchor = composition + within_state + cross` must hold per unit/program
and for every fitted contrast within tolerance `1e-10`. The within-state channel
is conditional on an annotated state; it can include unmeasured within-state
mixture changes and is never described as a causal within-cell response.

- Fit unweighted replicate-unit OLS with `factorial_design()` to total, the three
  channels, within-DAM means, within-Homeostatic means, and paired
  `DAM - Homeostatic` means. Ordinary feature-wise t inference uses the nine
  residual degrees of freedom; no cell pseudoreplication and no empirical-Bayes
  borrowing across only five programmes.
- Report all five contrasts. BH families are separate per endpoint x contrast
  across the five fixed programmes.
- Meaningful score margin = 0.25 pooled cell SD. The same ordinary t distribution
  supplies TREAT-style minimum-effect tests and two one-sided equivalence tests
  with 90% CIs. Minimum-effect and TOST BH families are separate per endpoint x
  contrast across five programmes. Supported-beyond-margin, equivalent-within-
  margin, and unresolved are mutually exclusive evidence states; equivalence
  language requires affirmative TOST evidence.
- A sample-size/precision weighted unit fit is sensitivity only. The equal-unit
  estimate remains primary.
- Report signed channel effects. Channel-share ratios are allowed only when the
  total effect clears its meaningful margin and channels are directionally
  concordant; opposing/near-zero totals remain signed components, not unstable
  percentages.

### D. Bridge + integrated verdict

- Reaggregate pooled Homeostatic+DAM counts and compare all five gene-level
  contrasts with current whole-microglia `pb_de_microglia`: Spearman correlation,
  effect direction, and fixed marker-programme concordance. This quantifies the
  consequence of excluding the 3.44% IFN remainder; it has no biological pass
  threshold and is always reported.
- “More” = a supported meaningful DAM occupancy effect and/or score-composition
  channel in the defining contrast.
- “Different” = a supported meaningful within-DAM UCell programme effect with a
  same-direction raw-count rotation result. A direct DAM-minus-Homeostatic
  response is additionally required for the word **DAM-selective**; otherwise the
  shift is shared/within-state.
- Primary synthesis uses the interaction. `tau_in_nlgf` contextualizes the
  NLGF-bearing genotype difference, while both amyloid arms expose where the
  interaction originates.
- Fixed outcome classes: composition-dominant; composition + state-conditional;
  state-conditional without composition; unresolved. Mixed programme-specific
  outcomes remain mixed rather than forced into one global label.
- If no conditional programme passes, write “no supported within-DAM programme
  shift.” “Equivalent across the tested programmes” is allowed only if all five
  equivalence tests pass. “Identical” is never earned.

Method basis:

- Replicate-level pseudobulk controls single-cell pseudoreplication and false
  discoveries: Squair et al. 2021,
  https://www.nature.com/articles/s41467-021-25960-2.
- Multi-sample, multi-condition state-wise differential-state analysis:
  Crowell et al. 2020,
  https://www.nature.com/articles/s41467-020-19894-4.
- Rank-based per-cell programme scores used by the live object: Andreatta &
  Carmona 2021, https://doi.org/10.1016/j.csbj.2021.06.043.
- `glmmTMB` beta-binomial parameterization:
  https://glmmtmb.github.io/glmmTMB/reference/nbinom2.html.
- `limma` rotation gene-set tests and weighted linear modelling:
  https://bioconductor.org/packages/release/bioc/manuals/limma/man/limma.pdf.

## Locked implementation + report contracts

- New module = `R/state_decomposition.R`; reuse design/contrast, pseudobulk,
  permutation, and P2 exact-decomposition helpers where their contracts match.
- Intended compact target chain:
  `microglia_state_substrate -> microglia_state_response ->
  microglia_state_decomposition -> state_decomposition_figures`.
- `microglia_state_substrate` may retain only compact raw two-state pseudobulks,
  unit metadata/counts, and the necessary score frame; it must never retain or
  serialize the 612 MB Seurat parent.
- `_state_decomposition.qmd` loads only `state_decomposition_figures`; include it
  after `_modality.qmd`. Existing Figure 1-9 sources/order stay unchanged; P6 is
  Figure 10.
- Figure 10 = one conventional publication plate with fixed content: replicate
  DAM occupancy/effect; state-specific amyloid responses plus direct state
  differences; within-DAM fixed-programme effects; and signed
  composition/within-state/cross attribution. All five programmes remain present;
  exact geometry follows result density, never outcome-based feature selection.
- Preserve current visible grammar: numbered heading + figure + folded compact
  code; no visible title, TOC, caption, body prose, table, or global code menu;
  nonblank alt text retained.
- Final figure bundle <=5 MB. QMD/report dependencies cannot reach cell-level
  frames, pseudobulk matrices, or fitted model objects.
- Runtime assertions + deterministic synthetic algebra smokes are primary. Add a
  committed test only when it directly protects Figure 10 or an inferential
  verdict; no broad robustness suite.
- Every execution step ends with a targeted fresh build/smoke and
  `scripts/check.sh`, then one scoped commit.

## Heavy/gated preconditions

- S1 loads only cached `microglia_annotated` (~612 MB), not the 8.95 GB all-cell
  source. Record elapsed time/peak RSS where feasible and prove the parent is not
  retained by the substrate.
- Primary-state coverage must be >=95% overall and >=90% in every unit. Both
  states must have exactly 16 units, >=20 cells/unit, finite scores, positive
  libraries, and canonical full-rank design. Any failure blocks inference rather
  than relaxing a threshold.
- State labels remain the existing cluster-level primary annotation. Per-cell
  programme argmax is diagnostic only and cannot relabel cells for a better
  result.
- `microglia_state_substrate` serialized size ceiling = 25 MB;
  `microglia_state_response` = 60 MB; final qmd bundle = 5 MB.
- Beta-binomial, voom/limma, rotation, and OLS fits must be warning-clean with
  finite coefficients, SEs, CIs, p/q values, and declared multiplicity-family
  identifiers. Missing score programmes/units are explicit failures, not dropped
  rows; raw-count marker coverage is always recorded and a zero-gene set is
  marked untestable rather than fabricated.
- Synthetic gates precede live claims: exact channel identity, contrast
  reconstruction, correct interaction signs, delta-method finite-difference
  agreement, and TOST boundary behavior.
- Dependency/environment files stay byte-stable unless an existing locked
  package unexpectedly fails to expose a required API; that is a stop-and-review
  event, not permission for automatic installation.

## Steps

- [x] **S1 - compact substrate + hard gates.** Implement two-state extraction and
  raw-count aggregation in `R/state_decomposition.R`; add
  `microglia_state_substrate`; force a fresh build. Validate state coverage,
  16-unit alignment, >=20 cells/state/unit, libraries, score variance, gene
  mapping, design rank, timing/RSS, serialization size, and parent-object
  isolation. Acceptance = every heavy data gate passes, substrate <=25 MB, no
  environment or report change.
  DONE 2026-07-14: one sparse membership multiply emits aligned 33,683 x 16 raw
  count matrices for each state plus unit/state UCell means, pooled score scales,
  and exact feature/marker maps. Live gates: 22,363/23,160 primary cells (96.56%),
  unit coverage >=93.93%, state-unit cells 31-2,151, positive libraries, design
  rank 7/residual df 9, all five marker sets fully mapped. Fresh target = 2.3 s
  producer / 5.50 s end-to-end, 2,592,592 KiB peak RSS, 1,431,505 serialized
  bytes / 19,763,288 in-memory bytes, warning/error clean, parent isolation true.
  Both matrices equal the established two-pass pseudobulk implementation element
  for element; all score means/SDs independently reconstruct from the parent.
  `scripts/check.sh` green; dependency lock and report QMD/content unchanged.
- [x] **S2 - occupancy + state response.** Implement beta-binomial occupancy,
  standardized probability contrasts/diagnostics/permutation sensitivity,
  separate Homeostatic/DAM pseudobulk fits, paired state-difference fit, fixed
  programme rotations, and pooled-state bridge. Add
  `microglia_state_response`. Acceptance = five canonical contrasts across every
  declared endpoint; synthetic paired-contrast oracle exact; finite effect/CI/
  family tables; model diagnostics green; target <=60 MB; bridge complete without
  a forced agreement threshold.
  DONE 2026-07-14: `microglia_state_response` carries all five contrasts for
  beta-binomial/model-scale + batch-standardized probability occupancy, empirical-
  logit/permutation sensitivity, Homeostatic/DAM voom+treat, harmonic-weight paired
  DAM-minus-Homeostatic response + unweighted sensitivity, 9,999-rotation fixed
  programmes, and the pooled-state bridge. Live filtered genes = 13,599 Homeostatic /
  9,148 DAM / 9,123 paired / 14,474 pooled. Interaction occupancy = +0.174 absolute
  DAM fraction (95% CI 0.095-0.253): zero-null FDR 1.81e-5, but 0.10-margin FDR
  0.081 = unresolved at 5%; batch-stratified permutation p = 0.021. Interaction
  rotation supports the Homeostatic programme within Homeostatic cells (FDR 0.0055),
  with no DAM-state or direct state-difference programme at FDR <=0.05; S3 remains
  outcome-open. Synthetic paired reconstruction residual = 8.88e-16; probability
  gradient residual = 3.60e-11. Pooled-versus-whole effect rho = 0.982-0.994 across
  contrasts without a pass gate. Fresh target = 10.2 s / 10,407,494 serialized bytes /
  27,221,560 in-memory bytes, warning/error clean, fitted-parent isolation true; full
  report gate green and dependency/report surfaces unchanged.
- [ ] **S3 - exact channel decomposition + verdict.** Implement equal-unit UCell
  standardization, composition/within-state/cross channels, OLS/TOST/sensitivity,
  and the predeclared evidence classifier; add `microglia_state_decomposition`.
  Acceptance = unit- and contrast-level reconstruction residual <=1e-10; ordinary
  df = 9; all five programmes/endpoints/contrasts present; BH/TOST families
  explicit; every verdict traceable to fixed evidence fields and outcome-open
  wording.
- [ ] **S4 - compact Figure 10 payload.** Build `state_decomposition_figures` from
  compact estimates only and implement the publication plate grammar. Acceptance
  = fixed panels/programmes represented including null/unresolved rows, no
  outcome-selected genes, deterministic ordering, bundle <=5 MB, no heavy object
  reachable from the qmd payload, focused plot smoke warning-clean.
- [ ] **S5 - report integration + QA.** Add `_state_decomposition.qmd` after
  `_modality.qmd`; wire report source dependencies; update map/memory only for
  durable shipped facts. Run `scripts/check.sh`, rendered-DOM assertions, and
  Chromium PDF -> PNG inspection. Acceptance = 10 numbered figures and 10
  nonblank alts; Figures 1-9 source/order unchanged; no visible prose/captions/
  tables; no clipping, external refs, warning/error markers, or stale target
  reads; full gate green.

After S5 = CLOSE-OUT: adversarially review plan, code, models, figure, and every
claim boundary; fold durable results into history/memory/map; archive this plan;
reset roadmap Active plan; run the final gate; commit the scoped close-out.

## Definition of done

- The report distinguishes retained-nuclei DAM occupancy from conditional
  molecular response with replicate-unit inference, or labels the evidence
  unresolved.
- Every “more,” “different,” “DAM-selective,” and equivalence statement is backed
  by its declared direct contrast, meaningful margin, uncertainty, and
  multiplicity family.
- The exact aggregate-score decomposition reconstructs every canonical contrast;
  no causal mediation, cell-level p-value, significance-only comparison, or
  gene-hit tally carries the conclusion.
- IFN sparsity and absent Proliferative-dominant cells remain visible coverage
  facts while all five continuous programmes remain tested.
- Figure 10 is compact, conventional, legible, and additive; the stable nine
  figures and dependency lock remain unchanged.
- Documentation matches live wiring/results, and one scoped commit closes every
  execution step plus close-out.
