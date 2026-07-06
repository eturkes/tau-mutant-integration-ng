# Map - codebase wiring (grows with the build)

SEVEN-figure report: snRNAseq microglia (P1) + activation trajectory (P2) + a four-method
amyloid-response logFC scatter + off-diagonal functional-group score panel (2026-07-06 adds:
re-wire the GeoMx / proteome / phospho PRIMARY DE via the lean `R/modality_de.R`). Shows STRUCTURE
(what calls what, what lives where);
`memory.md` holds the WHY (facts/gotchas/decisions), `roadmap.md` the trajectory. The
mechanism / cross-modality / qc / story chapters + their targets + `R/mechanism.R` +
`R/crossmodality.R` + their tests stay DELETED (roadmap Ledger 2026-07-06). Residual dead code
(figures.R story/mechanism/crossmodality builders + `.fig_*` helpers + qc_figure_data, plot.R
concordance_plot, constants rbc) stays UNWIRED -> not in this map. The io.R geomx/proteome/phospho
loaders + data_paths geomx/proteomics/phospho/sample_key are RE-WIRED (modality DE). Keep current:
load order, data -> module -> output flow, cache producer -> consumer pairs.

## P0 spine

### Bootstrap (fresh clone -> green pipeline)
`scripts/install-sysdeps.sh`   # apt: build-essential gfortran libglpk40
  -> `scripts/install-rv.sh`       # pinned rv -> ~/.local/bin (MUST be on PATH)
  -> `scripts/install-quarto.sh`   # pinned quarto -> tools/quarto/<ver>/ + bin/ wrapper
  -> `rv sync`                     # rproject.toml -> rv.lock -> rv/library (R pkgs)
  -> `uv sync`                     # pyproject.toml -> uv.lock -> .venv (Python)
  -> `scripts/check.sh`            # canonical green-check: sync + tests + FORCE-render report + zero-fault enforce
       (or `Rscript -e 'targets::tar_make()'` to just build the DAG)

### R session activation (every R/Rscript launched in project root)
`.Rprofile`
  -> `rv/scripts/rvr.R`       # rv helper fns
  -> `rv/scripts/activate.R`  # Sys.which("rv") -> shells `rv info` -> options(repos) + .libPaths(rv/library)
  -> guard (in .Rprofile): non-interactive stop() unless rv/library in .libPaths()
     (rv off PATH / `rv info` fail / R-version mismatch -> fail loud, no silent global-lib fallback)

### Pipeline (targets DAG)
`_targets.R`
  - QUARTO_PATH = tools/quarto/bin/quarto + file.exists() preflight stop()  # pinned CLI; no PATH fallback
  - Stan preflight (OPTIONAL): prepend tools/rlib-stan + set CMDSTAN iff tools/cmdstan exists (scripts/install-cmdstan.sh)
       -> sccomp cross-check arm; absent -> propeller-only (fresh clone still green)
  - tar_option_set(packages="quarto", memory="transient", garbage_collection=TRUE, trust_timestamps=TRUE)
  - tar_source("R")                                     # loads every R/*.R pure fn
  R/ pure fns:
   constants.R (genotype_levels, contrast_definitions, microglia_identity_markers,
      canonical_microglia_markers, microglia_substate_levels, contam_signatures; data_paths -- all 5
      entries wired (snrnaseq + geomx/proteomics/phospho/sample_key [modality DE]); rbc_marker_symbols dead) |
      utils.R (`%||%`, write_tsv_safe) | io.R (load_snrnaseq + build_symbol_map + load_geomx /
      read_spectronaut_tsv / proteomics_sample_meta / match_intensity_columns / normalise_ptm_stub
      ALL wired [snRNAseq + modality DE]) | spine.R (spine_versions)
   + design.R: factorial_design (treatment ~tau+nlgf+tau_nlgf[+batch]) + make_contrast_matrix
      (cell-means ~0+genotype) -> the 5 canonical contrasts; two equivalent parameterisations |
      de_pb.R: pseudobulk_counts/build_pseudobulk (replicate=genotype_batch; `cells=` -> per-substate
      subset pre-aggregation), fit_limma_voom (voomWQW default + confint) / fit_limma_log (log-intensity),
      median_normalise, prevalence_filter, de_pseudobulk, stage_wise_test (stageR screen df1=rank=3 +
      Holm confirm), interaction_power, assert_complete_crossing.
   + plot.R: theme_tau (ggplot base theme; base_family="" -> device font, warning-free; installs
      saturated-but-controlled ggplot discrete defaults) + manual scale_colour/fill_genotype (+ scale_color_
      alias; limits/breaks=genotype_levels, drop=FALSE) + manual microglia-substate / tau-background / binary /
      direction scales + richer continuous scales (`scale_fill_rwb`, `scale_colour_rwb`; signed panels pass
      midpoint=0; count panels use a neutral sequential gradient). concordance_plot retained UNWIRED (P4-only);
      modality_interaction_scatter WIRED -> the four-method scatter (per-modality amyloid logFC panel: dashed y=x
      identity + zero crosshairs + OLS trend + within-method Q99 |x-y| off-diagonal repel labels + coord_equal 1:1);
      functional_group_score_plot WIRED -> Figure 7 category-score dumbbell facets (per modality/free-y, primary
      GO-BP keyword-union role or explicit fallback category per within-method Q99 Figure 6 off-diagonal gene/protein;
      phosphoproteomics uses Figure 6 parent-protein mean points, segment colour = P301S-MAPTKI, fill = dedicated
      high-separation MAPTKI/P301S endpoint pair, bottom guides stacked, size = scored items; no enrichment/FDR display).
      Report visual identity = theme.scss.
   + (P1-S1) microglia.R: reprocess_microglia (SCT-v2/glmGamPoi -> Harmony[batch] -> Louvain multi-res ->
      UMAP; seeds+threads -> @misc$reprocess_provenance; strips stale reduction-coord/cluster meta shadows) +
      marker_mean_by_cluster (post-Harmony substate-separation check; map symbols->ensembl first) +
      reprocess_thread_env (thread provenance snapshot).
   + (P1-S2) microglia.R: annotate_microglia (UCell-score identity+substate+aux+contam on SCT data ->
      drop clear contaminant clusters -> calibrated cluster-argmax substate labels; guards required QC meta;
      @misc$microglia_prune[qc_rationale+DAM/lineage medians, separation margins] + $substate_provenance)
      ORCHESTRATES pure helpers: marker_sets_to_ensembl (symbol sets->present ensembl, error empty/near-empty
      min_n, n_used) | zscale_signatures (per-signature z; non-finite->0) | assign_substate (argmax +
      ambiguous/unassigned buckets, eps-floor) | cluster_mean_z (per-cluster mean z) | flag_contaminant_clusters
      (id_med<0.15 OR mglike_frac<0.30, finite-guard).
   + (P1-S3) composition.R: test_composition (orchestrator) -> composition_results. composition_counts (per-sample
      [genotype_batch] x substate count table; drops globally-empty levels; covariate-constancy fail-loud) |
      run_propeller (CELL-MEANS ~0+genotype+batch via make_contrast_matrix -- speckle PropRatio needs per-genotype
      mean coefs; asserts balanced crossed design; getTransformedProps -> propeller.ttest per contrast; logit PRIMARY + asin) |
      run_sccomp (OFF-lock sccomp ~0+genotype+(1|batch), cores-only [source-verified]; withCallingHandlers muffles+
      records BOTH R warnings->attr AND cmdstanr message() sampler notes ->messages attr [incl. "Warning:"-prefixed
      divergence lines -- else they red the gate's ^Warning: scan on a FRESH build]; final-fit diagnostic_summary
      (quiet=TRUE) [divergent/treedepth/E-BFMI] ->diagnostics attr, hardened to NULL on API drift; c_rhat structurally NA) |
      sccomp_backend_ready (PROJECT-LOCAL CmdStan gate -> no global leak) | composition_concordance (sign/sig cross-
      method flag; every present method must cover all keys -> fail-loud). orchestrator: backend present + sccomp error
      -> LOUD (allow_sccomp_failure to skip). propeller LOCKED primary; sccomp OPTIONAL off-lock cross-check
      (scripts/install-cmdstan.sh -> tools/rlib-stan + tools/cmdstan).
   + (P1-S4) de_pb.R: run_pb_de_microglia (whole-MG) / run_pb_de_substate (per-substate, min-cell floor fit-or-skip)
      -> de_pseudobulk (build_pseudobulk -> factorial_design -> fit_limma_voom[voomWQW+robust eBayes] -> stage_wise_test
      -> interaction_power). stage_wise_test = stageR screen (omnibus F, df1=rank=3) + modified-Holm per-contrast confirm.
      dam_direction = amyloid->DAM concordance vs v1.
   + (P1-S5) microglia.R: microglia_report_data(microglia_annotated, symbol_map) -- extracts a COMPACT report
      frame from the 612MB microglia_annotated (per-cell {umap_1/2, genotype, substate, *_UCell_z} cell_frame +
      n_cells + prune/provenance summaries; asserts finite z + non-NA factors) PLUS the substate-marker dot-plot
      slot via substate_marker_panel (per-substate mean_expr+pct_expr + per-gene z off SCT `data`; maps marker
      symbols -> SCT ensembl via symbols_to_ensembl; SOLE heavy-object access) -> the `microglia_report` target
      keeps _microglia.qmd (+ every gate force-render) reading ~0.5MB, NOT the 612MB Seurat.
   + (P2-S1) trajectory.R: build_activation_trajectory (orchestrator -> microglia_trajectory) -- slingshot on
      harmony[1:15] of microglia_annotated, FORCED single Homeostatic->DAM lineage (2 substate super-clusters), IFN/
      Proliferative omitted (on_lineage flag + NA pt); compact per-cell frame + per-unit omitted fraction + dims
      {10,15,20}+all-retained sensitivity + score-axis concordance + provenance. ORCHESTRATES pure helpers:
      run_slingshot_lineage (slingshot fit -> DAM-terminal lineage pt; branch-safe) | score_axis_pseudotime (raw
      DAM_UCell-Homeostatic_UCell) | squeeze_unit_interval (Smithson-Verkuilen -> open (0,1) for the S3 beta GLMM) |
      trajectory_concordance (Spearman pt vs score-axis) | trajectory_provenance (pkg versions/seed/RNG/threads).
      helpers.R += make_trajectory_embedding.
   + (P2-S2a/b) trajectory.R: run_trajectory_progression (orchestrator -> trajectory_progression) -- collapses
      microglia_trajectory on-lineage pt to 16 genotype_batch summaries (pseudotime_per_replicate) -> factorial_design
      (9 resid df) -> weighted/OLS/bounded limma interaction fits (fit_trajectory_contrasts) + exact 3-channel Kitagawa
      shift-share (decompose_progression_vs_composition -> kitagawa_channels) + Freedman-Lane perm null
      (freedman_lane_interaction; WLS-as-OLS, pivot-free chol2inv, RNG-pure, SENSITIVITY not nominal-exact);
      PRE-REGISTERED primary BH {progression_cf, within_homeostatic} vs separate exploratory BH; mean_pt FLAGGED
      composition-conflated. S2a estimation primitives (pure, crossing-agnostic): derive_batch | pseudotime_per_replicate
      | ordinary_t_table | fit_trajectory_contrasts | kitagawa_channels | within_state_col. Pure-R; reads the COMPACT
      S1 target (never the 612MB Seurat).
   + (P2-S3) trajectory.R: glmmtmb_pt_sensitivity (-> trajectory_glmm_sensitivity) -- per-cell beta GLMM
      pt01 ~ tau*amyloid + batch + (1|unit) on the COMPACT S1 cell_frame; SUPPORTIVE (composition-conflated, like
      mean_pt -> corroborates the position shift, NOT progression). Degrade cascade: nonestimable / singular RE /
      !pdHess / non-convergence / degenerate se|z|p -> rank-normal LMM (same gate) -> RECORDED method="failed"
      (NA + fail_reason). A FIT failure never throws; MALFORMED INPUT fails LOUD. Records n_units + fail_reason.
      Helpers: .capture_quietly (muffle+record warnings AND messages) + .fit_health_ok (PURE unit-tested battery
      gate) + .fit_pt_interaction (fit + Wald-row-by-POSITION + positional-integrity guard).
   + (P2-S4a) trajectory.R: trajectory_report_data (-> trajectory_report target) -- bundles the 3 COMPACT
      trajectory targets (microglia_trajectory / trajectory_progression / trajectory_glmm_sensitivity, NEVER the
      612MB Seurat) into ONE ~0.34MB object: slim per-cell plotting frame + interaction table (primary+exploratory
      families; coef/CI/perm_p/FDR -- the comp_cf/progression_cf/cross decomposition channels are ROWS here) +
      per_unit + lineage_per_unit + sensitivity + glmm 13-name subset + provenance (incl. the 3 decomposition
      loadings). NO separate `decomposition` field (the qmd reads loadings from provenance + per-channel coefs from
      interaction rows). Two guard layers: up-front stopifnot validates the 3 INPUTS' schema; render-cleanliness
      POSTCONDITIONS validate the ASSEMBLED bundle (col-EXISTENCE before is.finite; measure uniqueness; finite
      coef/ci/p/fdr/perm_p; weighted mean_pt present+finite on all 5 contrasts; glmm/provenance scalars). Pure: no RNG/IO.
   + (modality DE) modality_de.R: run_geomx_de / run_proteome_de_24m / run_phospho_de_24m -> the
      geomx_de / proteome_de_24m / phospho_de_24m targets = PRIMARY per-contrast topTables (logFC keyed by the 5
      canonical contrasts) for the 3 non-snRNAseq modalities, restored LEAN from the deleted P4 crossmodality/
      mechanism (auxiliary sensitivity / run-index / decon-preflight arms NOT restored). GeoMx = voom+TMM + slide
      fixed effect + bio-unit duplicateCorrelation (geomx_count_matrix / geomx_meta / geomx_slide_design /
      .fit_geomx_voom / .geomx_top_tables -> $primary$top); proteome = protein-group-summed log2 median-normalised
      limma-trend (protein_group_features / aggregate_proteome_raw / prepare_proteome_24m_matrix /
      .limma_log_de_from_matrix -> $top); phospho = log2 median-normalised phosphosite limma-trend
      (phospho_feature_frame / positive_log2_matrix / prepare_phospho_24m_matrix -> $top). match_24m_bulk_columns
      (16/16 balanced 4/genotype, shared by both bulk arms) + reuses fit_limma_log / median_normalise /
      prevalence_filter (de_pb.R) + factorial_design(add_batch=FALSE) + io.R loaders.
   + figures.R (LIVE surface after teardown): figure_manifest (11 hyphenated fig-* ids: microglia 7 + trajectory 4;
      report renders 7 figures -- 5 manifest + fig-modality-amyloid-effect + fig-modality-functional-scores, NOT
      in the vestigial manifest) + compact
      inline builders microglia_figure_data / trajectory_figure_data / modality_logfc_scatter_data (qmd-ready slots,
      finite geom guards, pre-binned/top-row reductions; modality_logfc_scatter_data = per-modality
      {feature, label, gene_symbols, y=nlgf_in_maptki, x=nlgf_in_p301s, interaction=x-y} logFC-pair frames,
      key-aligned, plus `groups$summary` = primary functional-category aggregate scores over empirical Figure 6
      off-diagonal labels per method; phospho scatter collapses finite sites to parent-protein mean points) + visual_reduction_slot_map
      + visual_slot_coverage (gate-wired vestigial prose-slot check, memory.md relic note) + generic `.fig_*` geom
      helpers. Dead story/mechanism/crossmodality/qc builders remain in-file (roadmap Ledger) -- UNWIRED, not mapped.
   + report.R: render_report (-> report target) calls quarto::quarto_render(quiet=FALSE), then
      repair_embedded_lightbox (rewrites Quarto embedded-lightbox anchors from absent local
      `index_files/figure-html/*.png` hrefs to the already embedded data-URI img src; fails loud if a local
      lightbox href has no embedded image shape).
  targets (31):
  - `spine` <- spine_versions()  [R/spine.R]                       # R + core-pkg version provenance df
  - raw inputs (format="file"; mtime-tracked via trust_timestamps):
       snrnaseq_file <- data_paths$snrnaseq  |  geomx_file / proteomics_file / phospho_file / sample_key_file (modality DE)
  - modalities (format="qs"):
       microglia_seurat_raw <- load_snrnaseq(snrnaseq_file)            # RNA-only microglia 33683 x 26104
       symbol_map           <- build_symbol_map(microglia_seurat_raw)  # {ensembl,symbol} 33683 x 2
       geomx      <- load_geomx(geomx_file)                            # GeoMx WTA Seurat 19963 x 91 AOIs
       proteomics <- read_spectronaut_tsv(proteomics_file)            # 24M proteome PTM export tibble
       phospho    <- read_spectronaut_tsv(phospho_file)               # 24M phospho PTM export tibble
       sample_key <- proteomics_sample_meta(sample_key_file)          # 16-run 24M key {genotype, col_stub}
  - P1 microglia core (format="qs"; consumes the snRNAseq modality):
       microglia_processed  <- reprocess_microglia(microglia_seurat_raw)  # SCT+pca+harmony+12 clusters@0.4+umap (687MB)
       microglia_annotated  <- annotate_microglia(microglia_processed, symbol_map)  # UCell substates + prune {6,7,8,11}; 23160 cells, 612MB
       composition_results  <- test_composition(microglia_annotated)  # propeller(logit+asin) [+gated sccomp] x 5 contrasts; counts+stats+concordance (small qs)
       pb_de_microglia      <- run_pb_de_microglia(microglia_annotated, symbol_map)  # whole-MG pseudobulk DE x 5 contrasts (voomWQW+stageR) + DAM concordance (4.7MB)
       pb_de_substate       <- run_pb_de_substate(microglia_annotated)  # per-substate DE: Homeo+DAM fit, IFN/Prolif skip (min-cell floor); cell_counts (7MB)
       microglia_report     <- microglia_report_data(microglia_annotated, symbol_map)  # compact report frame (umap+substate+z) + substate_markers + prune/provenance; ~0.5MB (keeps gate render cheap)
       microglia_figures    <- microglia_figure_data(microglia_report, composition_results, pb_de_microglia, pb_de_substate, symbol_map)  # microglia figure slots, pre-binned DE volcanoes, composition aliases; ~429KB qs live
  - P2 interaction trajectory (format="qs"; consumes microglia_annotated):
       microglia_trajectory <- build_activation_trajectory(microglia_annotated)  # slingshot H->D pseudotime; per-cell frame + per-unit omitted-frac + sensitivity + concordance; serialized ~0.8MB / in-mem ~3.3MB
       trajectory_progression <- run_trajectory_progression(microglia_trajectory)  # 16-unit pseudotime summaries -> weighted/ols/bounded interaction fits + 3-channel Kitagawa decompose + Freedman-Lane null; primary BH {progression_cf, within_homeostatic}; reads COMPACT S1 target (no 612MB load)
       trajectory_glmm_sensitivity <- glmmtmb_pt_sensitivity(microglia_trajectory$cell_frame)  # per-cell beta GLMM tau:amyloid (degrade -> rank-normal LMM -> method="failed"); supportive, INDEPENDENT of trajectory_progression; ~0.3KB
       trajectory_report    <- trajectory_report_data(microglia_trajectory, trajectory_progression, trajectory_glmm_sensitivity)  # bundles the 3 compact targets -> one ~0.34MB render object; two-layer guards (input schema + assembled-bundle postconditions); keeps gate force-render cheap
       trajectory_figures   <- trajectory_figure_data(trajectory_report, composition_results)  # pseudotime density, DAM-fraction join, Kitagawa forest, decomposition/concordance/audit slots; ~18KB live
  - four-method amyloid-response scatter (format="qs"; consumes the 3 modalities + P1 pb_de_microglia/symbol_map):
       geomx_de        <- run_geomx_de(geomx)                          # GeoMx primary voom DE x 5 contrasts ($primary$top; 19959 genes kept)
       proteome_de_24m <- run_proteome_de_24m(proteomics, sample_key)  # 24M proteome limma-trend x 5 contrasts ($top; 3379 groups)
       phospho_de_24m  <- run_phospho_de_24m(phospho, sample_key)      # 24M phosphosite limma-trend x 5 contrasts ($top; 17707 rows)
       modality_scatter_figures <- modality_logfc_scatter_data(pb_de_microglia, symbol_map, geomx_de, proteome_de_24m, phospho_de_24m)  # 4 compact per-modality logFC-pair frames (y=nlgf_in_maptki, x=nlgf_in_p301s) + off-diagonal functional-group scores; compact
  - report_sources <- c("_quarto.yml", "index.qmd", "_microglia.qmd", "_trajectory.qmd", "_modality.qmd")  # file target; explicit qmd invalidation
    report_extra_files <- c("theme.scss", "assets/code-tools-fix.html", assets/fonts/*.woff2)  # file target; explicit theme/font/after-body invalidation
    `report` <- render_report(report_sources, report_extra_files, microglia_report, microglia_figures, trajectory_figures, modality_scatter_figures)  # ONE offline HTML; quarto_render quiet=FALSE -> Quarto/Pandoc warnings reach the gate log; post-render repairs embedded-lightbox hrefs to data URIs
       reads `_quarto.yml` (type default; render index.qmd; output report/; lang en-GB; freeze false)
            -> `index.qmd` (format html, embed-resources, lightbox=auto, theme=theme.scss; no prose body;
                no author metadata; execute.echo=true + code-fold + code-tools -> chunk code shown as
                collapsed <details> folds; include-after-body=assets/code-tools-fix.html re-binds
                Show/Hide All Code (stock selector misses the code-copy scaffold); immediately includes the 3 chapters)
                                                          --{{< include >}}--> `_microglia.qmd`
               (microglia chapter {#sec-microglia}: setup `options(warn=2)` -> chunk warnings fail the render;
                tar_load microglia_report + microglia_figures [compact] -> "Substate landscape" (substate-marker
                per-gene-z dot plot `fig-microglia-substate-markers`, substate + DAM-score UMAPs `fig-microglia-umap`,
                genotype-faceted substate UMAP `fig-microglia-umap-substate`) + "Amyloid expands the DAM compartment"
                (replicate-unit stacked substate bars `fig-microglia-unit-composition`). 4 of the 7 figures.)
                                                          --{{< include >}}--> `_trajectory.qmd`
               (trajectory chapter {#sec-trajectory}: setup `options(warn=2)`; tar_load trajectory_figures ->
                genotype x substate pseudotime density `fig-trajectory-pt-density` (geom_area,
                facet_grid(substate~genotype)). The 5th figure. @sec cross-refs resolve across the full doc.)
                                                          --{{< include >}}--> `_modality.qmd`
               (modality chapter {#sec-modality}: setup `options(warn=2)`; tar_load modality_scatter_figures ->
                four-panel amyloid-response scatter `fig-modality-amyloid-effect` (modality_interaction_scatter x4
                via patchwork::wrap_plots; per method y=logFC nlgf_in_maptki, x=logFC nlgf_in_p301s, dashed y=x
                identity + OLS + within-method Q99 |x-y| off-diagonal labels; phospho = parent-protein mean points)
                + off-diagonal functional-category score facets
                `fig-modality-functional-scores` (functional_group_score_plot; within-method Q99 Figure 6 off-diagonal
                genes/proteins, phosphoproteomics uses displayed parent-protein mean points, rows carry primary
                role/fallback categories + leading score labels; connected points show aggregate NLGF_MAPTKI and NLGF_P301S scores,
                segment colour = P301S-MAPTKI). The 6th and 7th
                figures.)
       `theme.scss` = deep-blue/teal/slate chrome + IBM Plex (9 woff2 in assets/fonts/, base64-inlined offline)
       + figure-output overflow override (prevents print/PDF scrollbar chrome over figures).
       Figure labels: every captioned figure chunk uses a hyphenated `fig-*` id + `fig-cap` + `fig-alt` (7 total).

### Report prose inventory (vestigial after teardown)
`scripts/prose_inventory.py` (stdlib Python; non-DAG utility) + `.agent/prose_replacement_manifest.tsv`:
  now gate-wired ONLY through test_figures.R's `visual_slot_coverage()` assertion (the figure/schematic
  prose-replacement set was emptied during caption-only curation -> passes trivially). Retained because the
  test still sources the manifest. Full prose-reduction history -> roadmap Ledger.

### Tests
`tests/test_*.R` each: source the R/ files it exercises + `tests/helpers.R` (expect_error, make_meta16,
make_fake_seurat = synthetic Seurat fixtures, make_trajectory_embedding = synthetic slingshot embedding),
run stopifnot checks (fail-loud, no testthat dep), print `ok - <name>`. Run from project root:
`Rscript tests/test_<x>.R`.
  - test_design.R : 5-contrast exact weights + factorial==cell-means equivalence (property)
  - test_composition.R : composition_counts shapes/empty-drop/constancy-guard + propeller direction (logit+asin) + balance-guard + concordance (incl. completeness fail-loud) + sccomp-gate logical
  - test_microglia.R : reprocess/annotate pure-helper + synthetic-Seurat fixtures (S1/S2) + microglia_report_data extract/guards (S5)
  - test_de_pb.R  : pseudobulk -> 16 cols, median/prevalence, fit_limma_voom/log smokes (S3) + cells= subset,
                    de_pseudobulk/stageR matrix/interaction MDE, run_pb_de_substate fit-or-skip, dam_direction (S4)
  - test_io.R     : io contract tests (pure helpers + loader fail-loud asserts on tempfiles)
  - test_plot.R   : device-free -- theme_tau/scale_*_genotype/substate/background/binary/direction/rwb +
                    concordance_plot + modality_interaction_scatter (7-layer y=x panel, coord_equal 1:1, finite filter) +
                    functional_group_score_plot (warning-free aggregate-score dumbbell facets) class/wiring
  - test_figures.R : figure_manifest (11) + microglia_figure_data / trajectory_figure_data /
                    modality_logfc_scatter_data (y=nlgf_in_maptki / x=nlgf_in_p301s axis mapping + key-align +
                    per-modality labels + gene-symbol tokens + off-diagonal functional-group scoring + finite-drop /
                    missing-contrast fail-loud) builder contracts +
                    visual_slot_coverage (manifest-driven slot coverage) + synthetic finite guards for qmd geom inputs
  - test_modality_de.R : restored DE pure helpers -- positive_log2_matrix (nonpositive->NA before log2),
                    protein_group_features + aggregate_proteome_raw (group sum, present->NA), geomx_slide_design
                    (full-rank cell-means + 5 canonical contrasts, <2-slide fail), match_24m_bulk_columns (16/16 balanced fail-loud)
  - test_report.R  : repair_embedded_lightbox rewrite/no-op/fail-loud cases for embedded Quarto lightbox anchors
  - test_trajectory.R : (P2-S1) score-axis/squeeze/concordance + run_slingshot_lineage (single H->D + branched DAM-terminal) + provenance; (P2-S2a) derive_batch + per-replicate summary + contrast fit + Kitagawa exact-pure reconstruction; (P2-S2b) freedman_lane_interaction (null/signal/determinism/RNG-purity/weighted) + run_trajectory_progression structure on the jitter>0 non-additive fixture [sources R/de_pb.R for assert_complete_crossing]; (P2-S3) .fit_health_ok degrade-gate branches + glmmtmb_pt_sensitivity (beta success / singular->failed / nonestimable->failed / unknown-genotype fail-loud); (P2-S4a) trajectory_report_data field/measure/contrast presence + finite inference + malformed-input expect_error cases

### Quality gate -- `scripts/check.sh` (fail-loud, `set -euo pipefail`; `CHECK_SKIP_SYNC=1` skips env sync)
  1. `rv sync` + `uv sync`  2. loop `tests/test_*.R` (each `options(warn=2)` -> stray warnings = errors)
  3. FORCE-render report (`tar_invalidate(any_of("report")); tar_make()`, tee'd to a log, `if !`-wrapped) ->
     two render-time catches in the SOURCES: every included `_*.qmd` setup `options(warn=2)`
     (R chunk warning -> render error) + `render_report()` uses `quarto_render(quiet=FALSE)` (Quarto/Pandoc `[WARNING]` -> log)
  4. enforce zero-fault: (a) `tar_meta(error,warnings)` all-NA scoped to manifest names + dynamic branches;
     (b) anchored render-log grep (`^[WARNING]`/`^Warning:`/...), exit 0=fault / 1=clean / >=2=infra.
Any error/warning/log-hit -> non-zero exit. Detail (force-render rationale, warn=2, scoping, anchored grep,
negative tests) -> memory.md Quality gate.

### Config: tracked vs regenerated
tracked : rproject.toml rv.lock | pyproject.toml uv.lock .python-version | _targets.R R/*.R tests/*.R |
          _quarto.yml index.qmd _microglia.qmd _trajectory.qmd _modality.qmd theme.scss assets/code-tools-fix.html assets/fonts/*.woff2 |
          .Rprofile rv/scripts/*.R rv/.gitignore | scripts/install-*.sh scripts/prose_inventory.py |
          AGENTS.md .agents/skills/** .codex/prompts/*.md
regen   : rv/library _targets/ report/ _freeze/ .quarto/ .venv tools/  (gitignored + read-economy skip);
          sccomp_draws_files/ (sccomp per-chain CSV draws at build CWD; gitignored)
