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
  - tar_option_set(packages="quarto", memory="transient", garbage_collection=TRUE, trust_timestamps=TRUE)
  - tar_source("R")                                     # loads every R/*.R pure fn
  R/ pure fns (S2): constants.R (genotype_levels/colours, contrast_definitions, marker lists,
      rbc_marker_symbols, data_paths) | utils.R (`%||%`, write_tsv_safe) | io.R (loaders) | spine.R
   + (S3) design.R: factorial_design (treatment ~tau+nlgf+tau_nlgf[+batch]) + make_contrast_matrix
      (cell-means ~0+genotype) -> the 5 canonical contrasts; two equivalent parameterisations |
      de_pb.R: pseudobulk_counts/build_pseudobulk (replicate=genotype_batch), fit_limma_voom
      (counts) / fit_limma_log (log-intensity), median_normalise, prevalence_filter.
      S3 = machinery only -> NO new targets; P1+ wires the DE targets (consumes design + de_pb).
   + (S4) plot.R: theme_tau (ggplot base theme; base_family="" -> device font, warning-free) +
      scale_colour/fill_genotype (+ scale_color_ alias; limits/breaks=genotype_levels, drop=FALSE) +
      concordance_plot (two-effect scatter, P4 cross-modality). Report visual identity = theme.scss.
   + (P1-S1) microglia.R: reprocess_microglia (SCT-v2/glmGamPoi -> Harmony[batch] -> Louvain multi-res ->
      UMAP; seeds+threads -> @misc$reprocess_provenance; strips stale reduction-coord/cluster meta shadows) +
      marker_mean_by_cluster (post-Harmony substate-separation check; map symbols->ensembl first) +
      reprocess_thread_env (thread provenance snapshot).
   + (P1-S2) microglia.R: annotate_microglia (UCell-score identity+substate+aux+contam on SCT data ->
      drop clear contaminant clusters -> calibrated cluster-argmax substate labels; @misc$microglia_prune +
      $substate_provenance) ORCHESTRATES pure helpers: marker_sets_to_ensembl (symbol sets->present ensembl,
      error-on-empty, n_used) | zscale_signatures (per-signature z; zero-var->0) | assign_substate (argmax +
      ambiguous/unassigned buckets) | cluster_mean_z (per-cluster mean z) | flag_contaminant_clusters
      (id_med<0.15 OR mglike_frac<0.30). constants.R now carries microglia_identity_markers (pan QC) +
      canonical_microglia_markers (Homeostatic/DAM/IFN/Proliferative+MHC_APC) + microglia_substate_levels +
      contam_signatures (Oligo/Neuron/Astro).
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
  - `report` <- tar_quarto(path=".", quiet=FALSE, extra_files=c("theme.scss", assets/fonts/*.woff2))  # ONE offline HTML; quiet=FALSE -> Quarto/Pandoc warnings reach the gate log
       reads `_quarto.yml` (type default; render index.qmd; output _report/; lang en-GB; freeze false)
            -> `index.qmd` (format html, embed-resources, theme=theme.scss) --{{< include >}}--> `_qc.qmd`
               (QC-sanity chapter: setup `options(warn=2)` -> chunk warnings fail the render; tar_load 4
                modalities + sample_key -> dims, 16x16 design bijection, bounds)
       `theme.scss` = crimson colours (#B0344D) + IBM Plex (9 woff2 in assets/fonts/, base64-inlined offline)

### Tests (S3; gate-wired at S5)
`tests/test_*.R` each: source the R/ files it exercises + `tests/helpers.R` (expect_error,
make_meta16, make_fake_seurat = synthetic Seurat fixtures), run stopifnot checks (fail-loud,
no testthat dep), print `ok - <name>`. Run from project root: `Rscript tests/test_<x>.R`.
  - test_design.R : 5-contrast exact weights + factorial==cell-means equivalence (property)
  - test_de_pb.R  : pseudobulk -> 16 cols, median/prevalence, fit_limma_voom/log smokes
  - test_io.R     : io contract tests (pure helpers + loader fail-loud asserts on tempfiles)
  - test_plot.R   : device-free -- theme_tau/scale_*_genotype/concordance_plot class + wiring checks

### Quality gate (S5; review-hardened)
`scripts/check.sh` (fail-loud, `set -euo pipefail`; `CHECK_SKIP_SYNC=1` skips env sync):
  1. `rv sync` + `uv sync`  2. loop `tests/test_*.R` (each `options(warn=2)` -> stray warnings = errors)
  3. FORCE-render report (`tar_invalidate(any_of("report")); tar_make()`, tee'd to a log, `if !`-wrapped) ->
     two render-time catches in the SOURCES: `_qc.qmd` `options(warn=2)` (R chunk warning -> render error) +
     `_targets.R` `tar_quarto(quiet=FALSE)` (Quarto/Pandoc `[WARNING]` -> log)
  4. enforce zero-fault: (a) `tar_meta(error,warnings)` all-NA scoped to manifest names + dynamic branches;
     (b) anchored render-log grep (`^[WARNING]`/`^Warning:`/...), exit 0=fault / 1=clean / >=2=infra.
Any error/warning/log-hit -> non-zero exit. Detail (force-render rationale, warn=2, scoping, anchored grep,
negative tests) -> memory.md Quality gate.

### Config: tracked vs regenerated
tracked : rproject.toml rv.lock | pyproject.toml uv.lock .python-version | _targets.R R/*.R tests/*.R |
          _quarto.yml index.qmd _qc.qmd theme.scss assets/fonts/*.woff2 | .Rprofile rv/scripts/*.R
          rv/.gitignore | scripts/install-*.sh
regen   : rv/library _targets/ _report/ _freeze/ .quarto/ .venv tools/  (gitignored + deny-Read)
