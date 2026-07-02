# Figure Expansion Plan

Opened: 2026-07-02. Task: add many more figures to the final standalone report.

## Goal

Increase visual evidence density without changing biological claims or breaking the cheap-render contract.

Current source count = ~16 rendered figure chunks:
`_synthesis` 1, `_qc` 3, `_microglia` 4, `_trajectory` 3, `_mechanism` 2,
`_crossmodality` 3. Target = +20 to +26 additional figures, mostly atlas-style
multi-panel plots, so final report lands near 36-42 figure chunks while the synthesis
chapter stays answer-first.

## Research Digest

- Repo wiring: final report = `index.qmd` + included `_*.qmd`; force-rendered by
  `scripts/check.sh`; every section uses `options(warn=2)`. Existing compact report
  targets are small enough to support more plotting: `microglia_report` ~1.1MB,
  `trajectory_report` ~0.7MB, `mechanism_report` ~0.13MB, `crossmodality_report`
  ~0.16MB, `synthesis_report` ~45KB.
- Data constraint: qmds must not load heavy Seurat/raw modality targets beyond the
  existing QC exception. Any extra data extraction happens in targets, then qmd loads
  compact atlas data only.
- v1 mining: reuse visual idioms, not v1 claims. Good forms: DAM composition bars,
  interaction/null volcanoes, method/cross-modality concordance heatmaps, raw-vs-
  corrected phospho scatter, kinase/NF-kB signed heatmaps, clearance-pair lollipop/
  grid, convergence matrix. Exclude or rewrite v1-only overclaims: Gsk3b support,
  NF-kB attenuation, progression/rate synergy, full CCC, SpatialDecon abundance.
- Web/API check: Quarto supports figure cross-references with `fig` labels and warns
  against underscores in labels; Quarto HTML lightbox can improve large-figure reading;
  `targets` dynamic branching exists, but a single compact atlas target is simpler for
  this report; current `ggplot2` + `patchwork` is still the right local stack.

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

## Default Plan: Compact Figure Atlas

Add one new target-backed atlas chapter plus small helper infrastructure:

- `R/figures.R`: pure figure-atlas data builders + small plotting helpers only where
  repeated grammar warrants it.
- `_targets.R`: `figure_atlas` target, fed by compact/current analysis targets plus
  selected moderate tables (`pb_de_*`, `mechanism_*`, `crossmodality_*`) as needed.
  Output stays compact and qmd-safe; expected size target <5MB unless justified.
- `_figures.qmd`: new included chapter, probably after `_synthesis.qmd` and before
  `_qc.qmd`, titled "Evidence figure atlas". Minimal prose; captions carry the
  interpretation. Existing P1-P4 chapters remain the audit trail.
- `index.qmd`: include `_figures.qmd`; optionally enable Quarto `lightbox: auto` if
  the rendered offline HTML remains clean and usable.
- `tests/test_figures.R`: schema/finite guards for every atlas data frame that feeds
  a geom; reject missing anchors and figure-count drift.

Expected figures (+22 baseline, expandable):

1. Synthesis evidence map: claim x source-status tile from `synthesis_report`.
2. Microglia genotype-faceted UMAP by substate.
3. UMAP score triptych: homeostatic / DAM / MHC-APC z.
4. 16-unit substate composition bars by genotype_batch.
5. Substate score distributions by genotype and substate.
6. Composition method concordance/significance grid.
7. Whole-microglia volcano small multiples across all 5 contrasts.
8. Substate DE fit/skip and min-cell audit figure.
9. DAM/Homeostatic within-substate DE paired volcano or effect grid.
10. Pseudotime density by genotype x substate.
11. Unit-level mean pseudotime vs DAM fraction scatter, batch-labelled.
12. Kitagawa channel forest across all 5 contrasts.
13. Trajectory sensitivity/omitted-fraction audit.
14. Project pathway heatmap by population x contrast.
15. GO top-pathway dot plot by population/contrast.
16. TF activity lollipop matrix, Myc and NF-kB family highlighted.
17. NF-kB two-primary-row discordance tile, explicitly "not supported".
18. Kinase activity/run-index heatmap, Gsk3b present-but-not-supported highlighted.
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
- User chooses default atlas, inline expansion, publication-gallery, or a hybrid.
- Plan revised if choice changes scope or figure budget.

### S1 - Atlas Data Contract

Work:
- Add `R/figures.R::figure_atlas_data(...)` with compact, named slots:
  `synthesis`, `microglia`, `trajectory`, `mechanism`, `crossmodality`.
- Wire target `figure_atlas` after current report bundles and selected non-heavy
  summary targets. Avoid raw modality / 612MB Seurat reads in qmd.
- Add `tests/test_figures.R` with finite/geoms guards and explicit expected slots.

Acceptance:
- Fresh `tar_make(figure_atlas)` warning-clean.
- `figure_atlas` object size recorded and justified if >5MB.
- Test proves no required figure slot is empty; geom-fed numeric columns finite or
  intentionally handled with `na.rm` plus explicit prose.

### S2 - Microglia + Trajectory Atlas

Work:
- Implement figures 2-13 in `_figures.qmd` or a `_figures_microglia.qmd` include.
- Captions emphasise composition, under-powered interaction DE, and no supported
  progression beyond composition.

Acceptance:
- All figures render under `options(warn=2)`.
- No new claim says "rate", "acceleration", "absence", or "progression synergy"
  beyond the closed P2 wording.
- At least 10 new figures visible in rendered report.

### S3 - Mechanism Atlas

Work:
- Implement figures 14-18.
- Make unsupported mechanism evidence visible: Myc supported; NF-kB attenuation
  discordant/not supported; Gsk3b covered but not recovered; bulk kinase caveats
  explicit.

Acceptance:
- Figure captions derive status from `mechanism_report` / `figure_atlas`, not
  hardcoded current margins.
- NF-kB and Gsk3b plots cannot be visually mistaken for positive support.
- Render stays warning-clean.

### S4 - Cross-Modality Atlas

Work:
- Implement figures 19-26.
- Pull only compact/prepared data into qmd; any raw/corrected phospho scatter prep
  happens in `figure_atlas_data`, not in the qmd.

Acceptance:
- GeoMx repeated-AOI/blocking, SpatialDecon defer, no full CCC, bulk-not-
  microglia-sorted, and run-index sensitivity are visually or caption-explicit.
- Cross-modality counts use `modality_class`, not layer-level `modality_group`.
- Render stays warning-clean.

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
  wording.
- Close-out review accepts/fixes correctness, claim-honesty, and render-risk issues.

## Alternatives

### Alternative A - Inline Chapter Expansion

Add figures directly inside `_microglia`, `_trajectory`, `_mechanism`, and
`_crossmodality`.

Pros: each figure appears next to its local interpretation; fewer section jumps.
Cons: chapters become long and harder to scan; more chances to disturb polished
P1-P5 prose; repeated setup code unless helper discipline is strict.

Choose if the report should read like a full results manuscript rather than a lean
answer plus atlas.

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
