# Spatial decon follow-up - plan

## Scope
Follow up the closed P4 negative: SpatialDecon abundance was not earned because
no compact reference profile existed. This phase asks one narrow question:

Can the rebuilt project earn a GeoMx tissue-abundance readout, and if so does it
change the final synthesis?

Default route = gated SpatialDecon, broad-first and two-stage-aware:
- build a compact snRNAseq-derived profile from the full reference, not from the
  report bundles;
- run SpatialDecon on GeoMx Q3-normalised expression with Q3-scaled negative
  background;
- keep nuclei-rescaled absolute counts disabled while 42/91 ROIs carry `-1`;
- test log-beta abundance across the same 5 canonical contrasts;
- integrate the result into clearance/cross-modality/synthesis only if the
  profile and fit pass predeclared quality gates.

Out by default: new CCC model, v1 ledger revival, absolute cell counting,
claiming plaque-niche localisation from geometric whole-tissue ROIs, or forcing
microglial subpopulation abundance if the reference profile is collinear.

## Research Digest
Local current state:
- P4 already added the right seams: `geomx_q3_scaled_background`,
  `profile_collinearity`, `geomx_decon_preflight`, and
  `fit_geomx_abundance_de` in `R/crossmodality.R`, with synthetic tests.
- At plan open, live `geomx_de$decon_preflight` = `defer`: `SpatialDecon` was available from
  pinned repos as 1.22.0, Q3/background were usable, nuclei had 42 sentinels, and
  the compact reference profile was still absent.
- `SpatialDecon` is not installed in the current project library. S1 must add it
  to `rproject.toml` / `rv.lock`, smoke-load it warning-clean, and capture any
  package messages like the optimiser-heavy targets do.
- Current full `snrnaseq.rds` RNA rownames are ENSMUSG; symbols live in
  `@misc$geneids` / `symbol_map`. v1's direct-MGI-rowname shortcut is invalid
  for this fresh rebuild.
- Current microglia annotation has coherent Homeostatic, DAM, and IFN clusters;
  Proliferative has no coherent cluster. Do not fabricate a proliferative
  profile from noisy per-cell secondary labels.

v1 archive mining:
- v1 Arc L used SpatialDecon in a two-stage shape: broad profile for total
  microglia, then subpopulation fractions anchored to total microglia. Useful.
- Durable gotchas: Q3-scale negative-probe background, disable nuclei count
  conversion, quantify profile collinearity, residualise spatial autocorrelation
  by genotype/slide, and interpret subpopulation abundance as model-estimated tissue
  composition, not proof.
- Rejected v1 carry-over: ledger rows, full prose chapter shape, hardcoded
  common-gene counts, direct symbol rownames, and treating collinear subpopulation
  output as if it were equally stable as broad-cell output.

Current docs / method sweep:
- SpatialDecon is still the GeoMx-native route: its README states minimal inputs
  as normalised expression, same-scale background, and cell-profile matrix, and
  the package includes custom-profile and GeoMx-background helpers.
  Sources: https://github.com/Nanostring-Biostats/SpatialDecon,
  https://rdrr.io/bioc/SpatialDecon/man/spatialdecon.html,
  https://rdrr.io/bioc/SpatialDecon/man/derive_GeoMx_background.html
- Danaher/Kim 2022 frames SpatialDecon as log-normal regression with background
  modelling for spatial expression deconvolution, designed for region-level
  abundance estimates rather than single-cell placement.
  Source: https://www.nature.com/articles/s41467-022-28020-5
- SpatialDecon's cell-count conversion accepts nuclei counts, but the local
  sentinel pattern makes this inappropriate here. Keep beta/log-abundance scale.
  Source: https://rdrr.io/bioc/SpatialDecon/man/convertCellScoresToCounts.html
- RCTD/spacexr and cell2location are credible spatial decon alternatives, but
  they target spot/pixel count data and heavier platform-effect/Bayesian
  workflows. They are alternatives, not the default for GeoMx WTA.
  Sources: https://github.com/dmcable/spacexr,
  https://bioc.r-universe.dev/spacexr,
  https://www.nature.com/articles/s41587-021-01139-4

## Default Design
Broad-first SpatialDecon with subpopulation attempt guarded separately.

1. Reference profile is the key deliverable.
   Build a compact `geomx_reference_profile` target from the full snRNAseq RDS,
   `microglia_annotated`, `symbol_map`, and `geomx`. It loads the heavy full
   object once, maps ENSMUSG -> symbols, caps cells per class before any dense
   profile build, then stores only profile matrices and QC.

2. Profile levels are explicit.
   Broad profile = non-microglia broad annotations plus pooled retained
   microglia. Subpopulation profile = non-microglia broad annotations plus coherent
   microglia subpopulations that pass min-cell gates. Proliferative is recorded as
   absent unless the current annotation actually contains a coherent cluster.

3. Earning rules are predeclared.
   Broad abundance can be earned if common-gene coverage, per-class cell counts,
   profile collinearity, SpatialDecon fit, finite beta, and warning capture pass.
   Subpopulation abundance is optional and earns only if its own profile/stability
   gates pass; a collinear subpopulation profile becomes a reported negative, not a
   tuning prompt.

4. Deconvolution uses local GeoMx scale correctly.
   Use RNA `data` / Q3-normalised linear expression for `norm`. Build background
   as `NegGeoMean / q_norm_qFactors` broadcast over genes. Do not call
   `convertCellScoresToCounts()` while nuclei sentinels remain.

5. Abundance inference mirrors P4.
   Fit `log(beta + offset)` by slide fixed effect plus
   `duplicateCorrelation(block = bio_unit)` through existing
   `fit_geomx_abundance_de`, with unblocked sensitivity. Add a compact
   residualised spatial-autocorrelation check only as descriptive QC.

6. Report/synthesis become honest to the new state.
   If SpatialDecon remains blocked/not-earned, update the report from
   the old missing-profile negative to the actual reason. If earned, update
   `clearance_axis`, `crossmodality_report`, `synthesis_report`, and report
   prose so the final answer reflects the abundance layer without reviving full
   CCC.

## Alternatives
Alternative A - broad-only SpatialDecon.
Pros: strongest estimability, smallest profile surface, likely cleanest tissue
composition readout. Cons: cannot test DAM tissue subpopulation abundance; only total
microglia and other broad classes.

Alternative B - full two-stage SpatialDecon as primary.
Pros: closest to v1 Arc L and directly asks whether DAM abundance in tissue
corroborates snRNAseq composition. Cons: microglia subpopulations are transcriptionally
close; current Proliferative is absent; high risk of collinear, unstable output.

Alternative C - RCTD/spacexr or cell2location side branch.
Pros: modern spatial decon families with platform-effect or Bayesian modelling.
Cons: less GeoMx-native, heavier dependency/object surface, likely a new phase
rather than a follow-up; not needed until SpatialDecon fails in an informative
way.

Default choice: broad-first SpatialDecon with a gated subpopulation attempt. It gives a
real chance to earn spatial abundance while making "subpopulation not estimable" a
valid outcome.

## Steps
Each step is one closing unit. Run `scripts/check.sh` unless the step is
explicitly docs-only and a lighter check is justified.

### S0 - Route gate [DONE 2026-07-02]
User chose the recommended default: broad-first SpatialDecon with a gated
subpopulation attempt. Proceed to S1 without changing the default design.

Acceptance:
- User chooses default / alternative A / alternative B / alternative C / another
  route.
- If the user chooses away from default, revise this plan and roadmap before
  implementation.

### S1 - Dependency + compact reference profile [DONE 2026-07-02]
Added `SpatialDecon` to the project-local R lock and built the compact reference
profile target.

Result:
- `SpatialDecon` 1.22.0 installed project-locally and smoke-loaded under
  `warn=2` with no warnings/messages.
- `geomx_reference_profile` loads the full snRNAseq RDS once, overlays retained
  microglia subpopulations by barcode, caps at 500 cells/class, and stores only
  broad/subpopulation profiles plus QC/provenance.
- Live target warning-clean: 1.88 MB serialized; broad profile earned
  (15,919 genes x 6 profiles, max |cor| 0.674); subpopulation profile earned
  (16,079 genes x 8 profiles, max |cor| 0.902); Microglia_Proliferative is
  recorded absent.

Contracts:
- Add package via `rproject.toml` (`repository = "BioCsoft"`) + `rv sync`; no
  global library reliance.
- Smoke-load `SpatialDecon` under `warn=2`; capture or fix any startup messages.
- `geomx_reference_profile` loads full `snrnaseq_file` once, maps RNA ENSMUSG
  rows to unique symbols, caps cells per class with a fixed seed, and stores only
  compact profile/QC matrices.
- Reference labels combine full-object `broad_annotations` with
  `microglia_annotated` retained-cell subpopulations by cell barcode.
- QC records common genes, cells per class, absent/under-min classes, max profile
  correlation, condition number, memory estimate, package version, and seed.
- Heavy target uses transient memory/garbage-collection settings and drops the
  full Seurat object before returning.

Acceptance:
- Synthetic tests cover symbol mapping, cap determinism, missing-class handling,
  Proliferative-absent handling, and profile-collinearity gating.
- Fresh `geomx_reference_profile` build warning-clean with `tar_meta` clean.
- Output is compact enough for downstream targets; the heavy full Seurat object
  is not retained.
- `scripts/check.sh` green.

### S2 - SpatialDecon fit + two-stage assembly [DONE 2026-07-02]
Added `geomx_decon` target. It returns the reportable blocked state instead of
silently skipping after S1 earned the profile.

Result:
- GeoMx RNA `data` and Q3-scaled background were aligned and SpatialDecon ran
  warning-clean for both broad and subpopulation profiles.
- Broad arm: blocked because 4/91 AOIs have beta_total=0 after the fit
  (the same unresolved AOIs recur in the subpopulation arm).
- Subpopulation arm: independently attempted and blocked for the same 4 unresolved
  AOIs; two-stage subpopulation assembly is therefore blocked.
- Nuclei absolute rescaling remains disabled (42 sentinels). No abundance claim
  is earned at S2; S3 should pass through the blocked state and surface residual
  diagnostics rather than fitting log-beta DE.

Contracts:
- Use GeoMx RNA `data` layer as linear Q3-normalised `norm`; verify same AOI
  order as metadata and background.
- Background = `geomx_q3_scaled_background(meta)` broadcast to gene x AOI.
- Run broad profile first. Run subpopulation profile only as its own gated arm.
- Capture SpatialDecon warnings/messages; no leaked warning can reach tar_meta or
  the render log.
- Store beta, proportions derived from beta, fit residual/QC, profile gate
  status, unresolved AOI counts, and finite/bounds postconditions.
- Nuclei absolute rescaling remains disabled and recorded.

Acceptance:
- Tests cover SpatialDecon-result normalisation from synthetic beta, two-stage
  assembly with unresolved AOIs, finite/bounds guards, blocked-profile output,
  and warning/message capture.
- Fresh target build warning-clean. If blocked, reasons are precise and flow to
  report prose.
- Broad fit is required before claiming any abundance result; subpopulation failure
  does not invalidate a clean broad result.
- `scripts/check.sh` green.

### S3 - Abundance DE + spatial residual audit [DONE 2026-07-02]
Add `geomx_abundance_de` and a small spatial audit table from `geomx_decon`.

Result:
- Added `geomx_abundance_de`. It fits log-beta abundance only for earned
  deconvolution arms and otherwise returns canonical blocked outputs with empty
  5-contrast top tables.
- Live target warning-clean/tar_meta clean, 5.93 KB serialized. Broad,
  subpopulation, and microglia-subpopulation abundance DE are all blocked because the
  SpatialDecon fit still has 4 unresolved AOIs with beta_total=0.
- Residual audit is available despite blocked abundance: broad and subpopulation
  arms each cover 91 AOIs across 4 slides, using nearest-neighbour summaries of
  genotype-residualised per-AOI RMS residuals. Median RMS residual is ~0.821 for
  broad and subpopulation; descriptive only, no new claim axis.
- `scripts/check.sh` green across 49 current targets/branches.

Contracts:
- Fit log-beta abundance by the existing GeoMx design:
  `~0 + genotype + slide`, `duplicateCorrelation(block = bio_unit)`, 5 canonical
  contrasts, unblocked sensitivity.
- Broad abundance is primary. Subpopulation contrasts are shown only when subpopulation
  decon status is earned.
- Spatial audit = per-slide Moran's I or nearest-neighbour residual summary,
  genotype-residualised where estimable; descriptive only, not a new claim axis.
- Report all broad cell types and all attempted subpopulations, not only microglia.

Acceptance:
- Tests cover abundance-DE structure, canonical contrast coverage, blocked
  decon passthrough, broad/subpopulation gating, and residualised spatial-audit shape.
- Fresh target build warning-clean with `tar_meta` clean.
- Output states whether amyloid/interactions affect total microglia or DAM-like
  abundance, or whether the profile made that unearned.
- `scripts/check.sh` green.

### S4 - Report + synthesis integration [DONE 2026-07-02]
Wire a compact report bundle and update the closed synthesis surfaces to the new
abundance state.

Result:
- Added `spatial_decon_report`, a compact handoff over `geomx_decon`,
  `geomx_abundance_de`, and `geomx_reference_profile`. Live status is
  `blocked`, action `attempted`, with the real reason: 4 unresolved AOIs have
  near-zero total beta. It stores reference QC, arm summaries, unresolved AOIs,
  residual-audit summary, nuclei policy, and provenance, but no beta matrices.
- Rewired `clearance_axis` to accept the compact handoff, so
  `crossmodality_report` and `synthesis_report` are target-derived from the
  attempted SpatialDecon fit instead of the historical P4 preflight.
- Updated `_crossmodality.qmd`, `_synthesis.qmd`, and `index.qmd`: report prose
  now says SpatialDecon abundance is blocked after fitting, residual audit is
  descriptive fit QC, nuclei-rescaled absolute counts are disabled, and full CCC
  is not called.
- Targeted tests, live target rebuild, forced 142-chunk report render, and
  tar_meta warning/error check were green.

Contracts:
- Add `spatial_decon_report` target; qmds read only compact report bundles.
- Add `_spatial_decon.qmd` after `_crossmodality.qmd`, or fold a short section
  into `_crossmodality.qmd` if the result is blocked/not-earned.
- Revise `clearance_axis_data` to accept earned decon targets instead of failing
  when preflight becomes earned.
- Revise `crossmodality_report_data` and `synthesis_report_data` so
  SpatialDecon earned/not-earned status is target-derived, not hardcoded.
- Keep full CCC absent unless a later explicit phase opens it.

Acceptance:
- Report prose is inline-computed from compact targets; no v1 numbers.
- Synthesis evidence table updates the SpatialDecon row to earned / blocked /
  not-earned with the real reason.
- Forced report render warning-clean.
- `scripts/check.sh` green.

### S5 - Follow-up QA pass [DONE 2026-07-02]
Adversarially review the spatial-decon follow-up against the closed P1-P5 story.

Result:
- Accepted/fixed claim-scope wording issues: GeoMx captions now say
  bio-unit-blocked primary model/calls, so they are not confused with the blocked
  SpatialDecon abundance state.
- Accepted/fixed stale-negative wording: the P4 decon preflight reason now
  points to `geomx_reference_profile` / `geomx_decon` follow-up targets, and
  history/memory/spine summaries say SpatialDecon abundance is blocked after
  attempted fit rather than skipped or missing-profile.
- No new figure section was added in the follow-up, so existing hyphenated
  `fig-*` label QA from Figure expansion still applies.
- Full `scripts/check.sh` green; plan is ready for CLOSE-OUT mode.

Contracts:
- Check claim scope: model-estimated tissue abundance, not cell counts; GeoMx
  whole-tissue ROIs, not plaque niches; broad/subpopulation distinction clear.
- Check stale negatives: no remaining missing-profile negative after S1 built
  the reference profile.
- Check figure/caption labels if a new section adds figures.
- Update `.agent/memory.md`, `.agent/map.md`, `.agent/history.md`, and the
  cohesive story only for durable outcomes.

Acceptance:
- Accepted findings fixed before close-out. [DONE]
- Full `scripts/check.sh` green. [DONE]
- Plan ready for CLOSE-OUT mode. [DONE]

## Close-Out
After S5, run the standard close-out:
- adversarially review plan body + shipped code/prose;
- fold durable decisions into `.agent/history.md`;
- archive this plan to `.agent/completed/`;
- reset roadmap Active plan;
- commit `<scope> (spatial decon close): ...`.
