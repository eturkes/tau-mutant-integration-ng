# Figure Expansion Plan

Opened: 2026-07-02. Task: add many more figures to the final standalone report.

## Goal

Increase visual evidence density without changing biological claims or breaking the cheap-render contract.

Current source count = ~16 rendered figure chunks:
`_synthesis` 1, `_qc` 3, `_microglia` 4, `_trajectory` 3, `_mechanism` 2,
`_crossmodality` 3. Target = +20 to +26 additional figures, mostly inline
multi-panel plots inside existing result chapters, so final report lands near 36-42
figure chunks while the synthesis chapter stays answer-first.

## Research Digest

- Repo wiring: final report = `index.qmd` + included `_*.qmd`; force-rendered by
  `scripts/check.sh`; every section uses `options(warn=2)`. Existing compact report
  targets are small enough to support more plotting: `microglia_report` ~1.1MB,
  `trajectory_report` ~0.7MB, `mechanism_report` ~0.13MB, `crossmodality_report`
  ~0.16MB, `synthesis_report` ~45KB.
- Data constraint: qmds load compact report/figure targets only, beyond the existing
  QC exception. Any extra data extraction happens in targets, then qmds load compact
  per-chapter figure data.
- v1 mining: reuse visual idioms, not v1 claims. Good forms: DAM composition bars,
  interaction/null volcanoes, method/cross-modality concordance heatmaps, raw-vs-
  corrected phospho scatter, kinase/NF-kB signed heatmaps, clearance-pair lollipop/
  grid, convergence matrix. Exclude or rewrite v1-only overclaims: Gsk3b support,
  NF-kB attenuation, progression/rate synergy, full CCC, SpatialDecon abundance.
- Web/API check: Quarto supports figure cross-references with `fig` labels and warns
  against underscores in labels; Quarto HTML lightbox can improve large-figure reading;
  `targets` dynamic branching exists, but compact ordinary targets fit this report
  better; current `ggplot2` + `patchwork` is still the right local stack.

Sources consulted:
- Quarto figures/crossrefs/lightbox:
  `https://quarto.org/docs/authoring/figures.html`,
  `https://quarto.org/docs/authoring/cross-references.html`,
  `https://quarto.org/docs/output-formats/html-lightbox-figures.html`
- targets branching:
  `https://books.ropensci.org/targets/dynamic.html`
- ggplot2/patchwork current docs:
  `https://ggplot2.tidyverse.org/news/index.html`,
  `https://patchwork.data-imaginist.com/articles/guides/layout.html`

## Selected Plan: Inline Chapter Expansion

S0 decision: user selected inline expansion on 2026-07-02. Add figures directly
inside the existing report chapters, backed by compact per-chapter figure targets:

- `R/figures.R`: pure figure-data builders + small plotting helpers only where
  repeated grammar warrants it.
- `_targets.R`: compact targets `microglia_figures`, `trajectory_figures`,
  `mechanism_figures`, and `crossmodality_figures`; `synthesis_report` can feed the
  synthesis evidence map directly unless a tiny `synthesis_figures` target proves
  cleaner. Targets read compact report bundles and selected moderate summary tables
  (`pb_de_*`, `mechanism_*`, `crossmodality_*`) as needed. Each target stays qmd-safe;
  expected size <5MB unless justified.
- `_synthesis.qmd`, `_microglia.qmd`, `_trajectory.qmd`, `_mechanism.qmd`,
  `_crossmodality.qmd`: add figure chunks close to their local interpretation.
  Keep prose minimal; captions carry the added interpretation.
- `index.qmd`: no new atlas include. Optionally enable Quarto `lightbox: auto` if
  the rendered offline HTML remains clean and usable.
- `tests/test_figures.R`: schema/finite guards for every inline figure data frame
  that feeds a geom; reject missing anchors and figure-count drift.

Expected figures (+25 baseline, expandable):

Synthesis:
1. Evidence map: claim x source-status tile from `synthesis_report`.

Microglia:
2. Genotype-faceted UMAP by substate.
3. UMAP score triptych: homeostatic / DAM / MHC-APC z.
4. 16-unit substate composition bars by genotype_batch.
5. Substate score distributions by genotype and substate.
6. Composition method concordance/significance grid.
7. Whole-microglia volcano small multiples across all 5 contrasts.
8. Substate DE fit/skip and min-cell audit figure.
9. DAM/Homeostatic within-substate DE paired volcano or effect grid.

Trajectory:
10. Pseudotime density by genotype x substate.
11. Unit-level mean pseudotime vs DAM fraction scatter, batch-labelled.
12. Kitagawa channel forest across all 5 contrasts.
13. Trajectory sensitivity/omitted-fraction audit.

Mechanism:
14. Project pathway heatmap by population x contrast.
15. GO top-pathway dot plot by population/contrast.
16. TF activity lollipop matrix, Myc and NF-kB family highlighted.
17. NF-kB two-primary-row discordance tile, explicitly "not supported".
18. Kinase activity/run-index heatmap, Gsk3b present-but-not-supported highlighted.

Cross-modality:
19. GeoMx volcano small multiples for focus contrasts.
20. GeoMx blocked vs sensitivity support/loss figure.
21. Bulk run-index loss/flip heatmap.
22. Raw vs parent-corrected phospho effect scatter by contrast.
23. Anchor effect/coverage heatmap: clearance, synaptic, complement, tau/Gsk3b.
24. Clearance pair-support grid across focus contrasts.
25. Symbol x modality evidence matrix for top axis symbols.
26. Pathway axis-summary heatmap, broad modality counts + direction call.

## Steps

### S0 - Decision Gate

Acceptance:
- DONE 2026-07-02: user chose inline expansion.
- Plan revised from one atlas chapter to per-chapter inline figure targets; figure
  budget remains +20 to +26.

### S1 - Inline Figure Data Contract

Work:
- Add `R/figures.R` with compact builders:
  `microglia_figure_data`, `trajectory_figure_data`, `mechanism_figure_data`,
  `crossmodality_figure_data`; add `synthesis_figure_data` only if the direct
  `synthesis_report` load becomes awkward.
- Wire compact figure targets after current report bundles and selected non-heavy
  summary targets. Qmds load these compact targets, never raw modality data or the
  612MB Seurat object.
- Add `tests/test_figures.R` with finite/geoms guards, expected slot names, and the
  planned figure manifest by chapter.

Acceptance:
- DONE 2026-07-02: `R/figures.R` + targets `microglia_figures`,
  `trajectory_figures`, `mechanism_figures`, `crossmodality_figures`.
- Fresh target build warning-clean/tar_meta clean. Object sizes all <5MB:
  microglia 2.381MB, trajectory 0.033MB, mechanism 0.107MB, cross-modality
  0.514MB.
- `tests/test_figures.R` proves manifest slot presence and finite geom-fed numeric
  columns on deterministic fixtures. Heavy volcano/scatter shapes are pre-binned.

### S2 - Synthesis + Microglia Inline Figures

Work:
- Implement figures 1-9 in `_synthesis.qmd` and `_microglia.qmd`.
- Captions emphasise composition, under-powered interaction DE, and P1/P5 claim
  wording without changing the synthesis answer.

Acceptance:
- DONE 2026-07-02: `_synthesis.qmd` adds the evidence-map figure from
  `synthesis_report`; `_microglia.qmd` adds the 8 planned microglia figures from
  `microglia_figures`.
- Forced report render warning-clean under `options(warn=2)`; report chunks now
  include the 9 S2 `fig-*` labels.
- Captions preserve closed wording: robust amyloid-to-DAM activation, DAM
  composition interaction, under-powered interaction DE, and composition not
  progression beyond composition.

### S3 - Trajectory + Mechanism Inline Figures

Work:
- Implement figures 10-18 in `_trajectory.qmd` and `_mechanism.qmd`.
- Captions state the trajectory result as composition, not supported progression
  beyond composition.
- Make unsupported mechanism evidence visible: Myc supported; NF-kB attenuation
  discordant/not supported; Gsk3b covered but not recovered; bulk kinase caveats
  explicit.

Acceptance:
- DONE 2026-07-02: `_trajectory.qmd` adds the 4 planned trajectory figures from
  `trajectory_figures` plus guarded `trajectory_report` interaction rows:
  pseudotime density, unit mean-pseudotime vs DAM fraction, channel/decomposition
  forest, and robustness/omission audit.
- DONE 2026-07-02: `_mechanism.qmd` adds the 5 planned mechanism figures from
  `mechanism_figures`: all-population project pathway heatmap, compact GO dot
  plot, Myc/NF-kB-family TF lollipop, NF-kB discordance tile, and kinase/run-index
  heatmap with Gsk3b carried as unsupported unless significant.
- Full `scripts/check.sh` green: tests warn=2, forced 124-chunk report render,
  tar_meta clean, render-log clean; report chunks include the 9 S3 `fig-*` labels.
- Captions preserve closed wording: trajectory interaction = composition/not
  supported progression beyond composition; Myc supported; NF-kB attenuation
  discordant/not supported; Gsk3b covered but not recovered; bulk kinase caveats
  explicit.

### S4 - Cross-Modality Inline Figures

Work:
- Implement figures 19-26 in `_crossmodality.qmd`.
- Pull only compact/prepared data into qmd; any raw/corrected phospho scatter prep
  happens in `crossmodality_figure_data`, not in the qmd.

Acceptance:
- DONE 2026-07-02: `_crossmodality.qmd` adds the 8 planned cross-modality
  figures from `crossmodality_figures`: GeoMx volcanoes, GeoMx
  sensitivity/loss, bulk run-index, raw-vs-corrected phospho, bulk anchor
  heatmap, clearance-pair grid, symbol-modality matrix, and pathway-axis
  heatmap.
- DONE 2026-07-02: captions/visual encodings carry GeoMx blocked-AOI vs
  sensitivity, SpatialDecon defer/no full CCC, bulk-not-microglia-sorted, and
  run-index sensitivity caveats.
- DONE 2026-07-02: symbol/pathway counts use broad `n_modalities_sig`
  (`modality_class` semantics), not layer-level `modality_group`.
- DONE 2026-07-02: full `scripts/check.sh` green: tests warn=2, forced
  140-chunk report render, tar_meta clean, render-log clean.

### S5 - UX, Visual QA, Close-Out

Work:
- Decide and verify `lightbox: auto` vs no lightbox for embedded offline HTML.
- Add cross-reference labels using hyphenated `fig-*` ids, no underscores.
- Run full `scripts/check.sh`.
- Inspect rendered HTML with a small script: figure count, no missing captions, no
  warning/error classes, expected section/include present.
- Update `.agent/map.md`, `.agent/memory.md` only for durable wiring/gotchas; close
  roadmap; archive plan.

Acceptance:
- `scripts/check.sh` green.
- Final report has +20 or more new figure chunks and no stale "figure atlas open"
  or "atlas chapter" wording.
- Close-out review accepts/fixes correctness, claim-honesty, and render-risk issues.

## Alternatives Not Selected

### Alternative A - Compact Figure Atlas

Add one target-backed `_figures.qmd` chapter after synthesis and before QC.

Pros: easiest to keep P1-P5 prose untouched; concentrated visual appendix.
Cons: more section jumping; less local context for each figure.

Choose later only if inline chapters become unwieldy during S2-S4.

### Alternative B - Prebuilt Publication Gallery

Generate PNG/SVG files under a tracked `storage/figures/` or similar path, then
embed them with a Quarto lightbox gallery.

Pros: best for external sharing and visual polish; can use Python/matplotlib for
bespoke layout; HTML render becomes mostly image inclusion.
Cons: introduces regenerated assets and gitignore re-inclusion; harder to keep
figures automatically synced to target data; visual QA burden rises.

Choose if the main deliverable is a slide/paper figure bank, not an analytic HTML.

### Alternative C - Schematic-First Expansion

Build 3-5 large synthesis schematics/convergence figures instead of 20+ analytic
figures.

Pros: strongest narrative polish; compact report length.
Cons: does not satisfy "many more" literally; higher risk of over-smoothing
unsupported/null evidence.

Choose if the goal shifts from evidence density to presentation.

## Figure Claim Rules

- Figures can make existing evidence easier to see; they do not introduce new
  biological inference without a target/test step.
- Unsupported and unearned findings are plotted as first-class outcomes.
- Captions are factual, target-derived where possible, and avoid v1-overclaim words:
  progression/rate/acceleration support, NF-kB attenuation support, Gsk3b recovery,
  full CCC, SpatialDecon abundance.
- Human-facing prose follows repo register: British English; hyphens over em/en
  dashes; concise captions.
