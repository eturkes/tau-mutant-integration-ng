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
  - `book` <- tar_quarto(path=".")                      # renders the Quarto book
       reads `_quarto.yml` (type book; render `*.qmd` + `!rv/`; output _book/; freeze false)
            -> `index.qmd` (+ future analysis chapters reading the modalities)

### Config: tracked vs regenerated
tracked : rproject.toml rv.lock | pyproject.toml uv.lock .python-version |
          _targets.R R/*.R _quarto.yml index.qmd | .Rprofile rv/scripts/*.R rv/.gitignore |
          scripts/install-*.sh
regen   : rv/library _targets/ _book/ _freeze/ .quarto/ .venv tools/  (gitignored + deny-Read)
