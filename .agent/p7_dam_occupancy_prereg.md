# P7.4 DAM-occupancy robustness preregistration

Status: **FROZEN at P7.4, before any P7.5 labeling variant is generated or evaluated.**
This protocol applies to existing snRNAseq data only. P7.4 validates the reusable estimator
harness against the established current labeling and uses fabricated counts only for reduced-design
plumbing checks.

## 1. Estimand and endpoint

The primary estimand is the probability-standardized, equal-batch-marginalized DAM-fraction
interaction

`(NLGF_P301S - P301S) - (NLGF_MAPTKI - MAPTKI)`.

It is exactly the `interaction` row of `probability_contrasts` from the beta-binomial occupancy fit.
The live current-label baseline from `microglia_state_response$occupancy` is:

- estimate: **+0.1741141 fraction** (approximately +0.174);
- 95% CI: **[0.0954125, 0.2528158]**;
- zero-null FDR: **1.812987e-05** (approximately 1.81e-5);
- 0.10-margin minimum-effect FDR: **0.08116644**, unresolved at 5%.

## 2. Families and thresholds

- **Zero-null family:** `fdr_zero`, Benjamini-Hochberg adjusted over all five canonical contrasts.
  The interaction is "zero-null supported" if and only if `fdr_zero <= 0.05`.
- **0.10-margin minimum-effect family:** `fdr_minimum`, using margin 0.10 and
  Benjamini-Hochberg adjustment over all five canonical contrasts. The interaction is
  ">0.10 supported" if and only if `fdr_minimum <= 0.05`.

## 3. Estimators

Every labeling variant will be evaluated with all three estimators.

1. **E1 beta-binomial (PRIMARY).** Fit
   `cbind(n_DAM, n_Homeostatic) ~ 0 + genotype + batch` with a beta-binomial logit link,
   then equal-batch-standardize genotype probabilities and evaluate all five contrasts. This is the
   +0.1741141 primary estimand in section 1.
2. **E2 empirical-logit OLS plus permutation.** For every unit, calculate
   `log((n_DAM + 0.5) / (n_Homeostatic + 0.5))`, fit the shared factorial design/contrasts,
   and run a Freedman-Lane interaction test with residual permutations restricted within batch
   strata (`seed = 614`, `n_perm = 9999`). The current-label interaction baseline is coefficient
   **0.9340592** (OLS 95% CI **[0.1658742, 1.702244]**, five-contrast BH FDR **0.02806732**)
   with Freedman-Lane permutation p = **0.021**.
3. **E3 simple-proportion OLS.** Fit raw per-unit `DAM_fraction` through the shared factorial
   design/contrasts, without beta-binomial fitting or probability standardization. The current-label
   interaction baseline is estimate **+0.1733368 fraction**, with five-contrast BH FDR
   **0.01854462**.

## 4. Perturbation axes

Every perturbation is **one-at-a-time from the reference labeling**:
`primary_resolution = 0.4`, `id_floor = 0.15`, `mglike_floor = 0.30`, `tol = 0.10`, and
`amb_floor = 0.10`. Joint grid searches are not introduced.

### A. Clustering resolution

`primary_resolution` in `{0.2, 0.3, 0.4, 0.5, 0.6, 0.8}`, with 0.4 the reference.

### B. Pruning

- vary `id_floor` in `{0.10, 0.15, 0.20}` while holding `mglike_floor = 0.30`;
- vary `mglike_floor` in `{0.20, 0.30, 0.40}` while holding `id_floor = 0.15`;
- add one no-prune variant.

Reference settings are not counted as non-reference variants.

### C. Annotation

- vary `tol` in `{0.05, 0.10, 0.15}` while holding `amb_floor = 0.10`;
- vary `amb_floor` in `{0.05, 0.10, 0.15}` while holding `tol = 0.10`;
- UCell re-scoring is declared an **optional secondary** sensitivity.

Reference settings are not counted as non-reference variants.

### D. Leave one unit out

Sixteen variants, each dropping one `genotype_batch`, leaving 15 units.

### E. Leave one batch out

Four variants, each dropping one batch, leaving 12 units across three batches.

The frozen primary grid therefore contains approximately **34 non-reference labelings**, each
assessed with all three estimators.

## 5. Reduced-design handling

- **Leave one unit out (15 units):** retain all four batch levels; require a full-rank factorial
  design; standardize E1 equally over all present batches.
- **Leave one batch out (12 units):** `droplevels` the absent batch before `factorial_design`;
  construct the standardization grid as genotype times the three present batches.
- The original build-time cell-count and coverage gates in `build_microglia_state_substrate` are
  relaxed for variants only to structural estimability.
- Any non-estimable result -- including rank deficiency, an empty state in a unit,
  non-convergent beta-binomial fitting, or `pdHess != TRUE` -- is retained and recorded as
  `estimator_failed` with its reason. It is never silently dropped.

## 6. Decision and tipping-point rules

These rules are outcome-independent and all outputs are reported unconditionally.

- **Primary direction:** the E1 interaction sign. A **TIPPING** perturbation is any perturbation
  where the E1 interaction is `<= 0`.
- **Zero-null stability:** stable when E1 `fdr_zero <= 0.05`; a change away from that state is a
  **TIPPING** result.
- **0.10-margin behavior:** record every variant that either resolves
  `fdr_minimum <= 0.05` or has `abs(estimate) <= 0.10`.
- **Pre-committed robustness verdict:**
  **ROBUST-POSITIVE** if and only if E1 has `interaction > 0` and `fdr_zero <= 0.05` across
  **all** variants. Otherwise the verdict is **FRAGILE**, with the exact tipping perturbation set
  and the empirical `[min, max]` range of E1 interaction estimates reported.
- Cross-estimator concordance -- whether E2 and E3 signs track E1 -- is reported, but the verdict
  is defined on E1.
- The full variant-by-estimator table, empirical range, and tipping set are reported regardless of
  direction or significance.

## 7. Harness contract

`R/analysis/occupancy_harness.R` exposes:

- `membership_to_unit_coverage()`;
- `occupancy_family()`;
- `occupancy_from_membership()`.

On the current labeling,
`occupancy_family(microglia_state_substrate$unit_coverage, microglia_state_substrate$unit_meta)`
must reproduce `microglia_state_response$occupancy` for E1 `probability_contrasts`,
`probability_means`, `probability_vcov`, and `log_odds`, plus E2 `empirical_logit` and
`permutation`, with `all.equal(..., tolerance = 1e-8)`. The P7.4 validation target
`occupancy_harness_check` enforces this gate and records the E3 current-label baseline.
No labeling variant is generated or evaluated in P7.4.

## 8. Freeze declaration

This specification is **FROZEN at the P7.4 commit**. No variant occupancy was computed or viewed
before the freeze. The harness was validated only against the current labeling; its reduced-design
smoke used clearly fabricated positive counts unrelated to project data. P7.5 will execute this
frozen sweep on the existing snRNAseq data, and the pre-committed verdict and full results will be
reported regardless of direction.
