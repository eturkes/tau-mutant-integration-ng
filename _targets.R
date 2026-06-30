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

# Optional Stan backend for the sccomp composition cross-check (P1-S3). Prepend the off-lock,
# project-local cmdstanr library + point CMDSTAN at the compiled tree, but ONLY when both exist
# (provisioned by scripts/install-cmdstan.sh). Absent -> not added -> R/composition.R degrades to
# propeller-only (the reproducible primary), so a fresh clone still builds + the gate stays green.
local({
  stan_lib <- "tools/rlib-stan"
  cmdstan  <- Sys.glob(file.path("tools", "cmdstan", "cmdstan-*"))
  if (dir.exists(stan_lib) && length(cmdstan) >= 1L) {
    .libPaths(c(normalizePath(stan_lib), .libPaths()))
    Sys.setenv(CMDSTAN = normalizePath(cmdstan[1]))
  }
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

  # S3: substate composition across the 5 canonical contrasts. propeller (logit + asin) = locked
  # primary; sccomp (Bayesian, random batch) = optional cross-check, run iff the Stan backend above
  # is present, else recorded as skipped. Small result (count tables + per-contrast stats + concordance).
  tar_target(composition_results,  test_composition(microglia_annotated), format = "qs"),

  # S4: pseudobulk limma-voom DE across the 5 canonical contrasts. voomWithQualityWeights +
  # robust eBayes; stageR family-screen (omnibus F -> per-contrast Holm confirm). Whole-microglia
  # = the headline amyloid->DAM activation programme (+ DAM-direction concordance vs v1);
  # per-substate = fit-or-skip by a min-cell floor (Homeostatic/DAM fit, IFN/Proliferative thin).
  tar_target(pb_de_microglia, run_pb_de_microglia(microglia_annotated, symbol_map), format = "qs"),
  tar_target(pb_de_substate,  run_pb_de_substate(microglia_annotated),              format = "qs"),

  # S5: compact report-data extraction. Pulls the per-cell plotting frame (UMAP + substate +
  # activation z-scores) + the small prune/provenance summaries out of the heavy annotated object
  # so _microglia.qmd (and every force-rendered gate run) reads a ~0.5MB target, not the 612MB Seurat.
  tar_target(microglia_report, microglia_report_data(microglia_annotated), format = "qs"),

  # --- P2 interaction trajectory ---
  # S1: activation pseudotime. slingshot on the harmony embedding (dims 1:15), forced single
  # Homeostatic->DAM lineage (2 substate super-clusters), IFN/Proliferative omitted (off-lineage
  # flag, not deleted). Compact per-cell frame (pt_raw/pt01/score-axis + on-lineage flag) +
  # per-unit omitted fraction + a fixed sensitivity (dims 10 & 20 vs primary 15, + all-retained)
  # + score-axis concordance. Reads the 612MB annotated object once; stores a compact target
  # (serialized ~0.8MB, in-memory ~3.3MB; cheap-render).
  tar_target(microglia_trajectory, build_activation_trajectory(microglia_annotated), format = "qs"),

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
