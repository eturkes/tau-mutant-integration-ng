# targets pipeline entrypoint. Pure functions in R/ load via tar_source() (the DAG orders
# execution -- no manual loader). Data + heavy producers store format="qs" (qs2 backend).
library(targets)
library(tarchetypes)

# Point the `quarto` R package (used by tar_quarto / quarto_inspect) at the pinned,
# project-local Quarto CLI. Resolve the bin DIR (real) but leave the `quarto` symlink
# unresolved -> stable handle across version bumps. Fail loud if absent: a missing path
# lets the quarto pkg silently fall back to a PATH binary (unpinned) -> repro hole.
local({
  quarto_bin <- file.path(normalizePath("tools/quarto/bin", mustWork = FALSE), "quarto")
  if (!file.exists(quarto_bin)) {
    stop("pinned Quarto missing at ", quarto_bin, " -- run scripts/install-quarto.sh", call. = FALSE)
  }
  Sys.setenv(QUARTO_PATH = quarto_bin)
})

# memory="transient" + gc: release the ~8G snRNAseq load + its 340MB subset between targets.
# trust_timestamps: detect raw-input change by mtime/size, not by re-hashing the 8G file.
tar_option_set(
  packages = "quarto",
  memory = "transient",
  garbage_collection = TRUE,
  trust_timestamps = TRUE
)

tar_source("R")

list(
  # reproducibility-spine self-check: pinned-stack provenance via a tar_source()'d function
  tar_target(spine, spine_versions()),

  # --- raw input files (registered for DAG change-tracking; paths = data_paths, R/constants.R) ---
  tar_target(snrnaseq_file,   data_paths$snrnaseq,   format = "file"),
  tar_target(geomx_file,      data_paths$geomx,      format = "file"),
  tar_target(proteomics_file, data_paths$proteomics, format = "file"),
  tar_target(phospho_file,    data_paths$phospho,    format = "file"),
  tar_target(sample_key_file, data_paths$sample_key, format = "file"),

  # --- analysis-ready modalities (P1-P5 read these via tar_load; qs2 serialization) ---
  tar_target(microglia_seurat_raw, load_snrnaseq(snrnaseq_file),          format = "qs"),
  tar_target(symbol_map,           build_symbol_map(microglia_seurat_raw), format = "qs"),
  tar_target(geomx,                load_geomx(geomx_file),                 format = "qs"),

  # --- P1 snRNAseq microglia core ---
  # S1: reprocess (SCT-v2 + glmGamPoi) -> Harmony(batch) -> cluster (Louvain, multi-res) -> UMAP.
  # Heavy build; works on the 340MB RNA-counts subset (not the 8G load). qs2 serializes the Assay5
  # SCT object fine. memory="transient"+gc (global default) release it between targets.
  tar_target(microglia_processed,  reprocess_microglia(microglia_seurat_raw), format = "qs"),

  # S2: UCell substate scoring (identity + Homeostatic/DAM/IFN/Proliferative + MHC_APC aux + contam)
  # -> drop clear contaminant clusters -> calibrated argmax substate labels on the clean population.
  tar_target(microglia_annotated,  annotate_microglia(microglia_processed, symbol_map), format = "qs"),

  tar_target(proteomics,           read_spectronaut_tsv(proteomics_file),  format = "qs"),
  tar_target(phospho,              read_spectronaut_tsv(phospho_file),     format = "qs"),
  tar_target(sample_key,           proteomics_sample_meta(sample_key_file), format = "qs"),

  # Standalone HTML report render (path = project root with _quarto.yml; renders index.qmd, which
  # pulls in _qc.qmd via {{< include >}}). extra_files: quarto inspection tracks the .qmd target
  # deps but NOT the theme or its inlined fonts -> list theme.scss + the IBM Plex woff2 so editing
  # either reinvalidates the report (list.files keeps the 9 faces in sync with assets/fonts/).
  tar_quarto(
    report, path = ".",
    quiet = FALSE,   # surface Quarto/Pandoc warnings to the render log so the gate can scan them; default TRUE hides them
    extra_files = c("theme.scss",
                    list.files("assets/fonts", pattern = "\\.woff2$", full.names = TRUE))
  )
)
