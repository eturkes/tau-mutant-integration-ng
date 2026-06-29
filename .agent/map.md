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
  -> `Rscript -e 'targets::tar_make()'`   # build the DAG

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
  - `report` <- tar_quarto(path=".", extra_files=c("theme.scss", assets/fonts/*.woff2))  # ONE offline HTML
       reads `_quarto.yml` (type default; render index.qmd; output _report/; lang en-GB; freeze false)
            -> `index.qmd` (format html, embed-resources, theme=theme.scss) --{{< include >}}--> `_qc.qmd`
               (QC-sanity chapter: tar_load 4 modalities + sample_key -> dims, 16x16 design bijection, bounds)
       `theme.scss` = crimson colours (#B0344D) + IBM Plex (9 woff2 in assets/fonts/, base64-inlined offline)

### Tests (S3; gate-wired at S5)
`tests/test_*.R` each: source the R/ files it exercises + `tests/helpers.R` (expect_error,
make_meta16, make_fake_seurat = synthetic Seurat fixtures), run stopifnot checks (fail-loud,
no testthat dep), print `ok - <name>`. Run from project root: `Rscript tests/test_<x>.R`.
  - test_design.R : 5-contrast exact weights + factorial==cell-means equivalence (property)
  - test_de_pb.R  : pseudobulk -> 16 cols, median/prevalence, fit_limma_voom/log smokes
  - test_io.R     : io contract tests (pure helpers + loader fail-loud asserts on tempfiles)
  - test_plot.R   : device-free -- theme_tau/scale_*_genotype/concordance_plot class + wiring checks
S5 `scripts/check.sh` loops `tests/test_*.R` (non-zero exit on any failure).

### Config: tracked vs regenerated
tracked : rproject.toml rv.lock | pyproject.toml uv.lock .python-version | _targets.R R/*.R tests/*.R |
          _quarto.yml index.qmd _qc.qmd theme.scss assets/fonts/*.woff2 | .Rprofile rv/scripts/*.R
          rv/.gitignore | scripts/install-*.sh
regen   : rv/library _targets/ _report/ _freeze/ .quarto/ .venv tools/  (gitignored + deny-Read)
