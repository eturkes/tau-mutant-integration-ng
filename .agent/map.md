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
  - tar_option_set(packages = "quarto")
  - tar_source("R")                                     # loads every R/*.R pure fn
  - target `spine` <- spine_versions()   [`R/spine.R`]  # R + core-pkg version provenance df
  - target `book`  <- tar_quarto(path=".")              # renders the Quarto book
       reads `_quarto.yml` (type book; render `*.qmd` + `!rv/`; output _book/; freeze false)
            -> `index.qmd` (+ future analysis chapters)

### Config: tracked vs regenerated
tracked : rproject.toml rv.lock | pyproject.toml uv.lock .python-version |
          _targets.R R/*.R _quarto.yml index.qmd | .Rprofile rv/scripts/*.R rv/.gitignore |
          scripts/install-*.sh
regen   : rv/library _targets/ _book/ _freeze/ .quarto/ .venv tools/  (gitignored + deny-Read)
