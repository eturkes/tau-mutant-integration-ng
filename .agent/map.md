# Map - codebase wiring (grows with the build)

Fresh rebuild; modules land per phase. Shows STRUCTURE (what calls what, what lives
where); `memory.md` holds the WHY (facts/gotchas/decisions). Keep current: load order,
the data -> module -> output flow, and any cache producer -> consumer pairs.

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
  R/ pure fns (S2): constants.R (genotype_levels/colours, contrast_definitions, marker lists,
      rbc_marker_symbols, data_paths) | utils.R (`%||%`, write_tsv_safe) | io.R (loaders) | spine.R
   + (S3) design.R: factorial_design (treatment ~tau+nlgf+tau_nlgf[+batch]) + make_contrast_matrix
      (cell-means ~0+genotype) -> the 5 canonical contrasts; two equivalent parameterisations |
      de_pb.R: pseudobulk_counts/build_pseudobulk (replicate=genotype_batch; `cells=` -> per-substate
      subset pre-aggregation), fit_limma_voom (voomWQW default + confint) / fit_limma_log (log-intensity),
      median_normalise, prevalence_filter. S3 = machinery only; P1-S4 wires the DE targets.
   + (S4) plot.R: theme_tau (ggplot base theme; base_family="" -> device font, warning-free) +
      scale_colour/fill_genotype (+ scale_color_ alias; limits/breaks=genotype_levels, drop=FALSE) +
      concordance_plot (two-effect scatter, P4 cross-modality). Report visual identity = theme.scss.
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
      (id_med<0.15 OR mglike_frac<0.30, finite-guard). constants.R now carries microglia_identity_markers (pan QC) +
      canonical_microglia_markers (Homeostatic/DAM/IFN/Proliferative+MHC_APC) + microglia_substate_levels +
      contam_signatures (Oligo/Neuron/Astro).
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
      -> interaction_power). stage_wise_test = stageR screen (omnibus F, df1=rank=3) + modified-Holm per-contrast confirm. dam_direction =
      amyloid->DAM concordance vs v1. stageR added (rproject.toml, BioCsoft).
   + (P1-S5) microglia.R: microglia_report_data (extracts a COMPACT report frame from the 612MB
      microglia_annotated -- per-cell {umap_1/2, genotype, substate, *_UCell_z} cell_frame + n_cells +
      prune/provenance summaries; asserts finite z + non-NA factors) -> the `microglia_report` target keeps
      _microglia.qmd (+ every gate force-render) reading ~0.5MB, NOT the 612MB Seurat.
   + (P2-S1) trajectory.R: build_activation_trajectory (orchestrator -> microglia_trajectory) -- slingshot on
      harmony[1:15] of microglia_annotated, FORCED single Homeostatic->DAM lineage (2 substate super-clusters), IFN/
      Proliferative omitted (on_lineage flag + NA pt); compact per-cell frame + per-unit omitted fraction + dims
      {10,15,20}+all-retained sensitivity + score-axis concordance + provenance. ORCHESTRATES pure helpers:
      run_slingshot_lineage (slingshot fit -> DAM-terminal lineage pt; branch-safe) | score_axis_pseudotime (raw
      DAM_UCell-Homeostatic_UCell) | squeeze_unit_interval (Smithson-Verkuilen -> open (0,1) for the S3 beta GLMM) |
      trajectory_concordance (Spearman pt vs score-axis) | trajectory_provenance (pkg versions/seed/RNG/threads).
      rproject.toml += slingshot (BioCsoft; pulls princurve/TrajectoryUtils). helpers.R += make_trajectory_embedding.
  targets:
  - `spine` <- spine_versions()  [R/spine.R]            # R + core-pkg version provenance df
  - input files (format="file"): snrnaseq_file/geomx_file/proteomics_file/phospho_file/sample_key_file
       = data_paths$* (storage/data/*); change-tracked by mtime (trust_timestamps)
  - modalities (format="qs"; P1-P5 read via tar_load/tar_read):
       microglia_seurat_raw <- load_snrnaseq(snrnaseq_file)            # RNA-only microglia 33683 x 26104
       symbol_map           <- build_symbol_map(microglia_seurat_raw)  # {ensembl,symbol} 33683 x 2
       geomx                <- load_geomx(geomx_file)                  # Seurat 19963 x 91 AOIs
       proteomics           <- read_spectronaut_tsv(proteomics_file)   # tibble 45972 x 30
       phospho              <- read_spectronaut_tsv(phospho_file)      # tibble 64328 x 81
       sample_key           <- proteomics_sample_meta(sample_key_file) # tibble 16 x 4 (24M timepoint)
  - P1 microglia core (format="qs"; consumes the snRNAseq modality):
       microglia_processed  <- reprocess_microglia(microglia_seurat_raw)  # SCT+pca+harmony+12 clusters@0.4+umap (687MB)
       microglia_annotated  <- annotate_microglia(microglia_processed, symbol_map)  # UCell substates + prune {6,7,8,11}; 23160 cells, 612MB
       composition_results  <- test_composition(microglia_annotated)  # propeller(logit+asin) [+gated sccomp] x 5 contrasts; counts+stats+concordance (small qs)
       pb_de_microglia      <- run_pb_de_microglia(microglia_annotated, symbol_map)  # whole-MG pseudobulk DE x 5 contrasts (voomWQW+stageR) + DAM concordance (4.7MB)
       pb_de_substate       <- run_pb_de_substate(microglia_annotated)  # per-substate DE: Homeo+DAM fit, IFN/Prolif skip (min-cell floor); cell_counts (7MB)
       microglia_report     <- microglia_report_data(microglia_annotated)  # compact report frame (umap+substate+z) + prune/provenance; ~0.5MB (keeps gate render cheap)
  - P2 interaction trajectory (format="qs"; consumes microglia_annotated):
       microglia_trajectory <- build_activation_trajectory(microglia_annotated)  # slingshot H->D pseudotime; per-cell frame + per-unit omitted-frac + sensitivity + concordance; serialized ~0.8MB / in-mem ~3.3MB
  - `report` <- tar_quarto(path=".", quiet=FALSE, extra_files=c("theme.scss", assets/fonts/*.woff2))  # ONE offline HTML; quiet=FALSE -> Quarto/Pandoc warnings reach the gate log
       reads `_quarto.yml` (type default; render index.qmd; output _report/; lang en-GB; freeze false)
            -> `index.qmd` (format html, embed-resources, theme=theme.scss) --{{< include >}}--> `_qc.qmd`
               (QC-sanity chapter: setup `options(warn=2)` -> chunk warnings fail the render; tar_load 4
                modalities + sample_key -> dims, 16x16 design bijection, bounds)
                                                          --{{< include >}}--> `_microglia.qmd`
               (P1 microglia chapter: setup `options(warn=2)`; tar_load microglia_report + composition_results +
                pb_de_microglia + pb_de_substate + symbol_map -> substate UMAP, composition forest/table,
                amyloid->DAM volcano + DE counts, under-powered interaction + P2 pointer, Thrupp + dropout caveats)
       `theme.scss` = crimson colours (#B0344D) + IBM Plex (9 woff2 in assets/fonts/, base64-inlined offline)

### Tests (S3; gate-wired at S5)
`tests/test_*.R` each: source the R/ files it exercises + `tests/helpers.R` (expect_error,
make_meta16, make_fake_seurat = synthetic Seurat fixtures, make_trajectory_embedding = synthetic
slingshot embedding), run stopifnot checks (fail-loud, no testthat dep), print `ok - <name>`. Run
from project root: `Rscript tests/test_<x>.R`.
  - test_design.R : 5-contrast exact weights + factorial==cell-means equivalence (property)
  - test_composition.R : composition_counts shapes/empty-drop/constancy-guard + propeller direction (logit+asin) + balance-guard + concordance (incl. completeness fail-loud) + sccomp-gate logical
  - test_microglia.R : reprocess/annotate pure-helper + synthetic-Seurat fixtures (S1/S2) + microglia_report_data extract/guards (S5)
  - test_de_pb.R  : pseudobulk -> 16 cols, median/prevalence, fit_limma_voom/log smokes (S3) + cells= subset,
                    de_pseudobulk/stageR matrix/interaction MDE, run_pb_de_substate fit-or-skip, dam_direction (S4)
  - test_io.R     : io contract tests (pure helpers + loader fail-loud asserts on tempfiles)
  - test_plot.R   : device-free -- theme_tau/scale_*_genotype/concordance_plot class + wiring checks
  - test_trajectory.R : score-axis/squeeze/concordance + run_slingshot_lineage (single H->D + branched DAM-terminal) + provenance (P2-S1)

### Quality gate (S5; review-hardened)
`scripts/check.sh` (fail-loud, `set -euo pipefail`; `CHECK_SKIP_SYNC=1` skips env sync):
  1. `rv sync` + `uv sync`  2. loop `tests/test_*.R` (each `options(warn=2)` -> stray warnings = errors)
  3. FORCE-render report (`tar_invalidate(any_of("report")); tar_make()`, tee'd to a log, `if !`-wrapped) ->
     two render-time catches in the SOURCES: every section qmd (`_qc.qmd`, `_microglia.qmd`) setup `options(warn=2)`
     (R chunk warning -> render error) + `_targets.R` `tar_quarto(quiet=FALSE)` (Quarto/Pandoc `[WARNING]` -> log)
  4. enforce zero-fault: (a) `tar_meta(error,warnings)` all-NA scoped to manifest names + dynamic branches;
     (b) anchored render-log grep (`^[WARNING]`/`^Warning:`/...), exit 0=fault / 1=clean / >=2=infra.
Any error/warning/log-hit -> non-zero exit. Detail (force-render rationale, warn=2, scoping, anchored grep,
negative tests) -> memory.md Quality gate.

### Config: tracked vs regenerated
tracked : rproject.toml rv.lock | pyproject.toml uv.lock .python-version | _targets.R R/*.R tests/*.R |
          _quarto.yml index.qmd _qc.qmd _microglia.qmd theme.scss assets/fonts/*.woff2 | .Rprofile rv/scripts/*.R
          rv/.gitignore | scripts/install-*.sh
regen   : rv/library _targets/ _report/ _freeze/ .quarto/ .venv tools/  (gitignored + deny-Read);
          sccomp_draws_files/ (sccomp per-chain CSV draws at build CWD; gitignored)
