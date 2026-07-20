# Map - live codebase wiring

Current rendered surface (2026-07-16): lean 10-figure report. Closed P6 has compact
Homeostatic/DAM substrate, response, channel, paired gene-atlas, and Figure 10 payload
targets; report integration is complete.
No committed test suite, Python/uv surface, composition/sccomp/CmdStan arm,
retired P1 per-subpopulation DE target, prose-inventory
utility, mechanism/cross-modality/qc/story chapters, or retired agent configs. Historical
science lives in git + `roadmap.md` ledger; this file maps live code plus the rendered surface.

## Bootstrap

Fresh clone:
1. `scripts/bootstrap/sysdeps.sh` - apt libs/toolchain for current R stack.
2. `scripts/bootstrap/rv.sh` - project R package manager.
3. `scripts/bootstrap/quarto.sh` - pinned local Quarto CLI.
4. `rv sync` - `rproject.toml` -> `rv.lock` -> `rv/library`.
5. `scripts/check.sh` - fast force-render of the final HTML plus the DAM-occupancy harness, integration-substrate, integration-decomposition, and integration-concordance guardrails.

R activation:
- `.Rprofile` sources `rv/scripts/rvr.R` + `rv/scripts/activate.R`.
- Non-interactive R fails loud if project `rv/library` is absent from `.libPaths()`.

## Targets

`_targets.R` recursively sources `R/`, sets pinned `QUARTO_PATH`, and stores
heavy/intermediate objects as `format="qs"`. `_targets.yaml` routes the generated
store to `storage/targets/`. Expected live target count: 39.

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

P6 state decomposition (closed; report-integrated):
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
  in memory / 10.41 MB serialized.
- `microglia_state_gene_atlas <- run_microglia_state_gene_atlas(microglia_state_substrate)`
  jointly fits the 32 Homeostatic/DAM pseudobulks with `edgeR::voomLmFit`: unit block
  correlation, sample-quality weights, zero-aware residual df, robust limma eBayes/treat,
  four state/background amyloid effects, two paired DAM-minus-Homeostatic response
  differences, two state interactions, their direct difference, and joint 4-/2-df
  moderated-F families. Live = 14,438 genes, 16 paired units, correlation 0.183,
  1,120 joint-response / zero joint-interaction genes at FDR <=0.05; direct state
  differences yield 122/70 nonzero and 11/7 minimum-effect hits in MAPTKI/P301S;
  48/52 declared markers count-filter passing; 22.80 MB in memory / 7.20 MB serialized.
- `microglia_state_decomposition <- run_microglia_state_decomposition(microglia_state_substrate, microglia_state_response)`
  standardizes five raw-UCell programmes by pooled cell SD and emits exact equal-unit
  total/composition/within-state/cross channels, both state means + paired differences,
  ordinary OLS zero/minimum/TOST families, fixed cell-count-weighted sensitivity, compact
  S2 evidence, and the predeclared interaction classifier. Unit/contrast algebra, 9-df,
  family completeness, TOST boundary, size, and parent-isolation gates are runtime-fatal;
  live payload = 0.20 MB in memory / 0.054 MB serialized.
- `occupancy_harness_check <- check_occupancy_harness(microglia_state_substrate, microglia_state_response)`
  is the P7.4 non-report validation leaf: exact current-label E1/E2 reproduction, E3 baseline,
  complete aggregation round-trip including retained-denominator filler/coverage, strict E2/E3
  output schemas and finiteness, fabricated estimable 15-/12-unit smoke, and a separate fabricated
  all-estimator-failure fixture. `scripts/check.sh` invalidates and rebuilds it on every report gate.
- `occupancy_robustness <- run_occupancy_robustness(microglia_processed, microglia_annotated,
  microglia_state_substrate, symbol_map)` is the P7.5 non-report leaf: frozen 35-variant/105-estimator
  sweep, exact reference anchors, wide interaction table, E1 range/tipping/margin summaries,
  failure inventory, cross-estimator sign concordance, and outcome-independent verdict. Compact
  base-R payload only; explicitly built, absent from report ancestry and `scripts/check.sh`.
- `state_decomposition_figures <- state_decomposition_figure_data(microglia_state_response, microglia_state_gene_atlas)`
  emits the Figure 10 contract: accepted 16-unit occupancy; 208 fixed response rows =
  52 unique declared genes x four state/background amyloid effects, including 95% CI;
  and 28,876 transcriptome-wide paired state-difference rows = 14,438 genes x two tau
  backgrounds. All declared genes remain represented independently of outcomes; B2m's
  duplicate membership is collapsed. Payload = 4.45 MB in memory / 1.19 MB serialized,
  deterministic + parent-isolated. `state_decomposition_figure_plot()` draws occupancy,
  two direct DAM-minus-Homeostatic response maps, and 48 shared-scale per-gene factorial
  small multiples; six direct minimum-effect marker/background rows carry labels + CI,
  and four count-filter failures are listed. The report names only this compact leaf.

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
  `descriptive$GeoMx$sample_heatmap`, plus the TiO2 phospho protein-group PCA and phosphosite heatmap
  payloads (historical `Proteome`/`Phospho` keys; same TiO2 assay at two reporting levels).
  The amyloid-response scatter uses one shared off-diagonal feature cutoff: `|x-y| >= 3.5`.
  Figure 8 labels all points past the cutoff and draws all facets on one shared
  square coordinate range with a collected line legend for the dotted cutoff bands;
  the functional-category panel scores all same-cutoff selected features and displays categorized rows
  with aggregate `|P301S - MAPTKI| >= 0.5`; role assignment is a priority-ordered GO
  term-family pass that separates complement/MHC, phagocytosis, and chemotaxis before broader
  cell-cell-adhesion/extracellular-matrix/motility and immune residual buckets, and each visible category label
  lists every retained scored feature.
  The bulk context plate combines the TiO2 phospho protein-group-sum sample PCA with the same assay's phosphosite native heatmap.
  The heatmap selects 20 rows, excludes parent genes `Plcb1` and `Arhgef7`,
  keeps the same effect direction as the top-ranked candidate, collapses exact duplicate log2
  median-normalized profiles to the first ranked representative without label suffixes.
- `integration_substrate <- build_integration_substrate(pb_de_microglia, symbol_map, geomx_de, proteome_de_24m, phospho_de_24m)`
  is the P8.1 compact, parent-isolated NON-report leaf over exactly three modalities
  (`snRNAseq`, `GeoMx`, bulk protein-group). It stores raw and invertible robust-z `logFC`/moderated-`t`
  matrices on the five canonical contrasts; indexes 3,109 complete-case, 12,427 >=2-modality, and
  22,241 union symbols (pairwise 12,324/3,132/3,189); and keeps the 3,019 phosphosite parent-gene
  collapse as a same-TiO2-assay alternate, never a fourth modality. `scripts/check.sh` invalidates
  and rebuilds this self-validating leaf on every gate.
- `integration_decomposition <- build_integration_decomposition(integration_substrate)` is the P8.2
  compact, parent-isolated NON-report leaf. It decomposes the aligned 3,109 complete-case genes into
  joint/individual/residual Frobenius-energy shares on standardized logFC, with standardized
  moderated-t sensitivity, deterministic capped ranks and angle diagnostics, joint gene scores, and
  per-modality five-contrast loadings. In-builder reconstruction, orthogonality, rank-budget, planted
  rank-1/rank-2, and size oracles are runtime-fatal; `scripts/check.sh` rebuilds it on every gate.
- `integration_concordance <- build_integration_concordance(integration_substrate)` is the P8.3
  compact, parent-isolated NON-report leaf. It stores the common-3,109 raw-logFC Spearman 3 x 5
  matrix/tidy correlations, Pearson/moderated-t and per-pair coverage sensitivities, exact directional
  counts, and deterministic corrected-`phyper(q - 1, ...)` RRHO maxima. Bootstrap calibration is
  deferred; all outputs are descriptive-only. `scripts/check.sh` rebuilds it on every gate.

Report:
- `report_sources <- c("_quarto.yml", "index.qmd", sections/*.qmd, R/**/*.R)`
  so helper-only plot/source edits invalidate `report`.
- `report_extra_files <- c("assets/theme.scss", assets/fonts/*.woff2)`
- `report <- render_report(...)`

## Modules

`R/core/constants.R`
- Genotype levels, five canonical contrasts, marker sets, data paths.

`R/core/io.R`
- Load snRNAseq microglia subset, build symbol map, load GeoMx, parse Spectronaut
  exports, parse 24M sample key, match intensity columns.

`R/core/design.R`
- Shared 2x2 factorial design and cell-means contrast matrix for canonical contrasts.

`R/core/utils.R`
- Small shared operators and report-selection helpers.

`R/core/spine.R`
- Pinned-stack provenance target.

`R/analysis/de_pb.R`
- Pseudobulk counts, limma voom/log helpers, prevalence/median normalization,
  whole-microglia DE, crossing guard.

`R/analysis/microglia.R`
- Reprocess, annotate, compact report bundle. `microglia_report_data()` is the only
  report path that reads the heavy annotated Seurat object for microglia figures; it
  emits per-cell UMAP/scores, replicate-unit subpopulation composition, marker dot-plot data,
  prune/provenance.

`R/analysis/trajectory.R`
- Activation pseudotime, 16-unit progression/decomposition inference, supportive
  glmmTMB sensitivity, compact trajectory report bundle.

`R/analysis/state_decomposition.R`
- P6 compact Homeostatic/DAM substrate, occupancy/state response, paired multivariate gene
  atlas, and exact UCell-channel inference: feature/marker mapping,
  coverage/design/library/size gates, beta-binomial standardization, state/delta voom+treat,
  paired blocked voomLmFit + joint moderated F, rotations, bridge, ordinary
  OLS/TOST/weighted sensitivity, fixed classifier, and deterministic algebra/model gates.
  Heavy parents/fits stay out of target payloads.

`R/analysis/occupancy_harness.R`
- `membership_to_unit_coverage()` -> `occupancy_family()` -> `occupancy_from_membership()`;
  `check_occupancy_harness()` validates generalized E1 beta-binomial/E2 empirical-logit/E3
  simple-proportion inference, exact current-label reproduction, and fabricated reduced-design smoke.

`R/analysis/occupancy_sweep.R`
- `run_occupancy_robustness()` executes the frozen P7.5 reference + 34 one-at-a-time resolution,
  pruning, annotation, LOU, and LOBO variants from cached objects; records all three estimator
  interactions/failures, exact reference determinism, tipping/margin sets, sign concordance, and the
  all-variant E1 verdict. Returned target is compact and parent-isolated; Seurat objects/fits stay out.

`R/analysis/modality_de.R`
- Lean primary DE for GeoMx and the 24M TiO2 phospho assay: the historical `proteome_*` tokens
  are protein-group sums of TiO2 phospho-PTM rows (`PG.ProteinGroups`), `phospho_*` is the
  phosphosite view of the same acquisition. GeoMx also emits compact sample-heatmap fields for
  first-five-DAM-gene-clustered and mean-DAM-rotated AOI layout plus compact design/DAM-gene tracks.
  Both bulk DE producers carry a compact `$run_order_sensitivity` (rank-5 additive mean-centered
  `run_index`, 11 residual df) as an integrity record only; the genotype-blocked acquisition means
  run order is aliased with genotype, so figures keep the primary no-batch `$top`. Auxiliary
  SpatialDecon beta/abundance and broad-family sensitivity arms stay deleted.

`R/analysis/integration.R`
- `build_integration_substrate()` harmonizes snRNAseq, GeoMx, and bulk protein-group effects by
  gene symbol across the five canonical contrasts. `robust_z()` stores invertible per-modality,
  per-statistic, per-contrast median/MAD scaling; deterministic max-`AveExpr` representatives,
  exact source reconstruction, overlap counts, compactness, and parent isolation are runtime-fatal
  guards. The phosphosite parent-gene payload is explicitly a within-assay alternate.
- `build_integration_decomposition()` uses `integration_rank_signal()`,
  `integration_random_joint_threshold()`, `integration_joint_basis()`, and
  `integration_ajive_decompose()` for the pure-R AJIVE-style split; `integration_ajive_fixture()`
  self-tests planted joint ranks 1 and 2 with near-zero-noise reconstruction.
- `build_integration_concordance()` uses `integration_modality_pairs()`,
  `integration_hypergeom_upper()`, `integration_rrho_thresholds()`,
  `integration_rrho_maximum()`, and `integration_directional_overlap()` for common-universe and
  pairwise raw-effect concordance plus compact four-quadrant directional-overlap summaries.

`R/report/figures.R`
- Compact figure-data builders for rendered slots only:
  `microglia_figure_data()`, `trajectory_figure_data()`,
  `modality_logfc_scatter_data()`, `state_decomposition_figure_data()`.

`R/report/plot.R`
- Shared report theme, scales, modality and rendered descriptive plot helpers including
  `geomx_sample_heatmap_plot()`, `bulk_modality_context_plot()`, and the report-integrated
  `state_decomposition_figure_plot()`.

`R/report/render.R`
- Quarto render wrapper, embedded-lightbox repair, and report-dir pruning so the final HTML is the only
  user-facing output.

## Report

`index.qmd` includes four qmd fragments. The rendered HTML exposes simple numbered
`Figure 1` ... `Figure 10` headings, but no title/TOC/captions:
- `sections/microglia.qmd`: subpopulation marker dot plot, vertically stacked subpopulation/DAM UMAPs,
  genotype-faceted subpopulation UMAP, replicate-unit subpopulation composition.
- `sections/trajectory.qmd`: pseudotime density by genotype/subpopulation.
- `sections/modality.qmd`: GeoMx AOI metadata-track diagnostic,
  vertically stacked TiO2 phospho protein-group PCA / phosphosite heatmap descriptive figure,
  four-panel amyloid response scatter, functional-category score panel.
- `sections/state-decomposition.qmd`: compact two-tier plate with retained-state occupancy,
  transcriptome-wide two-state interaction geometry, and ungrouped line-profile fields
  spanning all declared marker genes.

Rendered output = 10 numbered figures plus compact per-figure folded code controls/content in
`report/tau-mutant-integration.html`; `render_report()` removes stale sibling outputs from `report/`.
The browser/tab title is `Tau Mutant Integration`.
Chunk setup uses `options(warn=2)`; data builders pre-filter/guard finite values so report
warnings are treated as real failures.

## Tracked vs Ignored

Tracked live source/config:
- `README.md`, `CLAUDE.md`
- `.claude/settings.json`, `.claude/commands/session-prompt.md`
- `.serena/.gitignore`, `.serena/project.yml`
- `.agent/{memory.md,map.md,roadmap.md,history.md,archive_digest.md,context.sh}` +
  `.agent/completed/`
- `_targets.R`, `_targets.yaml`, `_quarto.yml`, `index.qmd`, `sections/`, `R/`,
  `assets/`, `scripts/`
- `rproject.toml`, `rv.lock`, `.Rprofile`, `rv/scripts/*.R`

Ignored/generated/heavy:
- `.git/`, `rv/library/`, `report/`, `_freeze/`, `.quarto/`
- `storage/{data,cache,logs,qa,targets}/` and legacy/default `_targets/`
- project-local env caches, Claude local/runtime state, and rendered/static artefacts.
