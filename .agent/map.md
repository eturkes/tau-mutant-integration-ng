# Map - live codebase wiring

Current rendered surface (2026-07-14): lean 9-figure report. P6 has compact,
report-disconnected Homeostatic/DAM substrate + response targets; channel/figure/report wiring remains open.
No committed test suite, Python/uv surface, composition/sccomp/CmdStan arm,
retired P1 per-subpopulation DE target, prose-inventory
utility, mechanism/cross-modality/qc/story chapters, or retired agent configs. Historical
science lives in git + `roadmap.md` ledger; this file maps live code plus the rendered surface.

## Bootstrap

Fresh clone:
1. `scripts/install-sysdeps.sh` - apt libs/toolchain for current R stack.
2. `scripts/install-rv.sh` - project R package manager.
3. `scripts/install-quarto.sh` - pinned local Quarto CLI.
4. `rv sync` - `rproject.toml` -> `rv.lock` -> `rv/library`.
5. `scripts/check.sh` - fast force-render of the final HTML target only.

R activation:
- `.Rprofile` sources `rv/scripts/rvr.R` + `rv/scripts/activate.R`.
- Non-interactive R fails loud if project `rv/library` is absent from `.libPaths()`.

## Targets

`_targets.R` sources `R/*.R`, sets pinned `QUARTO_PATH`, and stores heavy/intermediate
objects as `format="qs"`. Expected live target count after P6-S2: 31.

Raw file targets:
- `snrnaseq_file`
- `geomx_file`
- `proteomics_file`
- `phospho_file`
- `sample_key_file`

Loaded modalities:
- `microglia_seurat_raw <- load_snrnaseq(snrnaseq_file)`
- `symbol_map <- build_symbol_map(microglia_seurat_raw)`
- `geomx <- load_geomx(geomx_file)`
- `proteomics <- read_spectronaut_tsv(proteomics_file)`
- `phospho <- read_spectronaut_tsv(phospho_file)`
- `sample_key <- proteomics_sample_meta(sample_key_file)`

P1 microglia:
- `microglia_processed <- reprocess_microglia(microglia_seurat_raw)`
- `microglia_annotated <- annotate_microglia(microglia_processed, symbol_map)`
- `pb_de_microglia <- run_pb_de_microglia(microglia_annotated)`
- `microglia_report <- microglia_report_data(microglia_annotated, symbol_map)`
- `microglia_figures <- microglia_figure_data(microglia_report)`

P6 state decomposition (active; report-disconnected through S2):
- `microglia_state_substrate <- build_microglia_state_substrate(microglia_annotated, symbol_map)`
  emits only two aligned raw-count pseudobulks, unit/state counts + libraries,
  unit/state raw-UCell means + pooled SDs, and exact feature/marker maps. Runtime
  gates fix the Homeostatic/DAM universe, 16 complete units, >=20 cells/state/unit,
  >=95% overall / >=90% per-unit coverage, full-rank design, finite variable scores,
  positive libraries, <=25 MB in-memory + qs-serialized payloads, and no Seurat/S4
  reachability.
- `microglia_state_response <- run_microglia_state_response(microglia_state_substrate, pb_de_microglia)`
  emits beta-binomial + standardized occupancy, empirical-logit/permutation sensitivity,
  state-wise voom/treat gene tables, harmonic-weight direct state differences + unweighted
  sensitivity, five fixed-programme rotations, and the pooled-two-state/whole-MG bridge.
  Runtime algebra/fit/family/size gates keep fitted objects out; live payload = 27.22 MB
  in memory / 10.41 MB serialized. No report target depends on P6 through S2.

P2 trajectory:
- `microglia_trajectory <- build_activation_trajectory(microglia_annotated)`
- `trajectory_progression <- run_trajectory_progression(microglia_trajectory)`
- `trajectory_glmm_sensitivity <- glmmtmb_pt_sensitivity(microglia_trajectory$cell_frame)`
- `trajectory_report <- trajectory_report_data(microglia_trajectory, trajectory_progression, trajectory_glmm_sensitivity)`
- `trajectory_figures <- trajectory_figure_data(trajectory_report)`

Modality context:
- `geomx_de <- run_geomx_de(geomx)` (primary DE + compact sample-heatmap descriptor)
- `proteome_de_24m <- run_proteome_de_24m(proteomics, sample_key)`
- `phospho_de_24m <- run_phospho_de_24m(phospho, sample_key)`
- `modality_scatter_figures <- modality_logfc_scatter_data(pb_de_microglia, symbol_map, geomx_de, proteome_de_24m, phospho_de_24m)`
  carries the first-five-DAM-gene-clustered AOI design/DAM-gene track atlas through
  `descriptive$GeoMx$sample_heatmap`, plus the proteome PCA and phosphoproteome heatmap payloads.
  The amyloid-response scatter uses one shared off-diagonal feature cutoff: `|x-y| >= 3.5`.
  Figure 8 labels all points past the cutoff and draws all facets on one shared
  square coordinate range with a collected line legend for the dotted cutoff bands;
  the functional-category panel scores all same-cutoff selected features and displays categorized rows
  with aggregate `|P301S - MAPTKI| >= 0.5`; role assignment is a priority-ordered GO
  term-family pass that separates complement/MHC, phagocytosis, and chemotaxis before broader
  cell-cell-adhesion/extracellular-matrix/motility and immune residual buckets, and each visible category label
  lists every retained scored feature.
  The bulk context plate combines the proteome sample PCA with the phosphoproteome native heatmap.
  The heatmap selects 20 rows, excludes parent genes `Plcb1` and `Arhgef7`,
  keeps the same effect direction as the top-ranked candidate, collapses exact duplicate log2
  median-normalized profiles to the first ranked representative without label suffixes.

Report:
- `report_sources <- c("_quarto.yml", "index.qmd", "_microglia.qmd", "_trajectory.qmd", "_modality.qmd", R/*.R)`
  so helper-only plot/source edits invalidate `report`.
- `report_extra_files <- c("theme.scss", assets/fonts/*.woff2)`
- `report <- render_report(...)`

## Modules

`R/constants.R`
- Genotype levels, five canonical contrasts, marker sets, data paths.

`R/io.R`
- Load snRNAseq microglia subset, build symbol map, load GeoMx, parse Spectronaut
  exports, parse 24M sample key, match intensity columns.

`R/design.R`
- Shared 2x2 factorial design and cell-means contrast matrix for canonical contrasts.

`R/de_pb.R`
- Pseudobulk counts, limma voom/log helpers, prevalence/median normalization,
  whole-microglia DE, crossing guard.

`R/microglia.R`
- Reprocess, annotate, compact report bundle. `microglia_report_data()` is the only
  report path that reads the heavy annotated Seurat object for microglia figures; it
  emits per-cell UMAP/scores, replicate-unit subpopulation composition, marker dot-plot data,
  prune/provenance.

`R/trajectory.R`
- Activation pseudotime, 16-unit progression/decomposition inference, supportive
  glmmTMB sensitivity, compact trajectory report bundle.

`R/state_decomposition.R`
- P6 compact Homeostatic/DAM substrate extraction + occupancy/state-response
  inference: exact feature/marker mapping, coverage/design/library/size gates,
  beta-binomial standardization, state/delta voom+treat, rotations, bridge, and
  deterministic algebra/model gates. Heavy parents/fits stay out of target payloads.

`R/modality_de.R`
- Lean primary DE for GeoMx, 24M proteome, and 24M phosphosite data. GeoMx also emits
  compact sample-heatmap fields for first-five-DAM-gene-clustered and mean-DAM-rotated
  AOI layout plus compact design/DAM-gene tracks.
  Auxiliary SpatialDecon beta/abundance, run-index, and sensitivity arms stay deleted.

`R/figures.R`
- Compact figure-data builders for rendered slots only:
  `microglia_figure_data()`, `trajectory_figure_data()`,
  `modality_logfc_scatter_data()`.

`R/plot.R`
- Shared report theme, scales, modality and rendered descriptive plot helpers including
  `geomx_sample_heatmap_plot()` and `bulk_modality_context_plot()`.

`R/report.R`
- Quarto render wrapper, embedded-lightbox repair, and report-dir pruning so the final HTML is the only
  user-facing output.

## Report

`index.qmd` includes three qmd fragments. The rendered HTML exposes simple numbered
`Figure 1` ... `Figure 9` headings, but no title/TOC/captions:
- `_microglia.qmd`: subpopulation marker dot plot, vertically stacked subpopulation/DAM UMAPs,
  genotype-faceted subpopulation UMAP, replicate-unit subpopulation composition.
- `_trajectory.qmd`: pseudotime density by genotype/subpopulation.
- `_modality.qmd`: GeoMx AOI metadata-track diagnostic,
  vertically stacked proteome PCA / phosphoproteome heatmap descriptive figure,
  four-method amyloid response scatter, functional-category score panel.

Rendered output = 9 numbered figures plus compact per-figure folded code controls/content in
`report/tau-mutant-integration.html`; `render_report()` removes stale sibling outputs from `report/`.
The browser/tab title is `Tau Mutant Integration`.
Chunk setup uses `options(warn=2)`; data builders pre-filter/guard finite values so report
warnings are treated as real failures.

## Tracked vs Ignored

Tracked live source/config:
- `AGENTS.md`, `.codex/prompts/session.md`, `.agents/skills/session-prompt/SKILL.md`
- `.agent/{memory,map,roadmap,history,archive_digest}.md`
- `_targets.R`, `R/*.R`, `*.qmd`, `_quarto.yml`, `theme.scss`, `assets/`, `scripts/*.sh`
- `rproject.toml`, `rv.lock`, `.Rprofile`, `rv/scripts/*.R`

Ignored/generated/heavy:
- `.git/`, `rv/library/`, `_targets/`, `report/`, `_freeze/`, `.quarto/`
- `storage/data/`, `storage/cache/`, `storage/logs/`, `.agent/completed/`
- project-local env caches and rendered/static artefacts.
