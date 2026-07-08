# targets pipeline entrypoint. Pure functions in R/ load via tar_source() (the DAG orders
# execution -- no manual loader). Data + heavy producers store format="qs" (qs2 backend).
library(targets)

# Point the `quarto` R package (used by render_report -> quarto_render) at the pinned,
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

  # --- analysis-ready modalities (P1-P2 read these via tar_load; qs2 serialization) ---
  tar_target(microglia_seurat_raw, load_snrnaseq(snrnaseq_file),          format = "qs"),
  tar_target(symbol_map,           build_symbol_map(microglia_seurat_raw), format = "qs"),

  # GeoMx spatial + 24M bulk proteome / phosphosite loads feed the cross-tau amyloid-response
  # scatter (below). Small tables (parsed Spectronaut exports + a Seurat object), not the 8G snRNAseq.
  tar_target(geomx,      load_geomx(geomx_file),               format = "qs"),
  tar_target(proteomics, read_spectronaut_tsv(proteomics_file), format = "qs"),
  tar_target(phospho,    read_spectronaut_tsv(phospho_file),    format = "qs"),
  tar_target(sample_key, proteomics_sample_meta(sample_key_file), format = "qs"),

  # --- P1 snRNAseq microglia core ---
  # S1: reprocess (SCT-v2 + glmGamPoi) -> Harmony(batch) -> cluster (Louvain, multi-res) -> UMAP.
  # Heavy build; works on the 340MB RNA-counts subset (not the 8G load). qs2 serializes the Assay5
  # SCT object fine. memory="transient"+gc (global default) release it between targets.
  tar_target(microglia_processed,  reprocess_microglia(microglia_seurat_raw), format = "qs"),

  # S2: UCell substate scoring (identity + Homeostatic/DAM/IFN/Proliferative + MHC_APC aux + contam)
  # -> drop clear contaminant clusters -> calibrated argmax substate labels on the clean population.
  tar_target(microglia_annotated,  annotate_microglia(microglia_processed, symbol_map), format = "qs"),

  # Whole-microglia pseudobulk limma-voom DE across the 5 canonical contrasts. The report
  # uses these topTables for the snRNAseq panel in the four-method amyloid-response scatter.
  tar_target(pb_de_microglia, run_pb_de_microglia(microglia_annotated), format = "qs"),

  # S5: compact report-data extraction. Pulls the per-cell plotting frame (UMAP + substate +
  # activation z-scores) + the small prune/provenance summaries out of the heavy annotated object
  # so _microglia.qmd (and every force-rendered gate run) reads a ~0.5MB target, not the 612MB Seurat.
  tar_target(microglia_report, microglia_report_data(microglia_annotated, symbol_map), format = "qs"),
  tar_target(microglia_figures,
             microglia_figure_data(microglia_report),
             format = "qs"),

  # --- P2 interaction trajectory ---
  # S1: activation pseudotime. slingshot on the harmony embedding (dims 1:15), forced single
  # Homeostatic->DAM lineage (2 substate super-clusters), IFN/Proliferative omitted (off-lineage
  # flag, not deleted). Compact per-cell frame (pt_raw/pt01/score-axis + on-lineage flag) +
  # per-unit omitted fraction + a fixed sensitivity (dims 10 & 20 vs primary 15, + all-retained)
  # + score-axis concordance. Reads the 612MB annotated object once; stores a compact target
  # (serialized ~0.8MB, in-memory ~3.3MB; cheap-render).
  tar_target(microglia_trajectory, build_activation_trajectory(microglia_annotated), format = "qs"),

  # S2b: progression-interaction inference. Collapse on-lineage pseudotime to the 16 genotype_batch
  # summaries -> factorial design (9 resid df) -> weighted / OLS / bounded interaction fits + the
  # exact 3-channel Kitagawa decomposition + a Freedman-Lane permutation null; PRE-REGISTERED primary
  # BH family {progression_cf, within_homeostatic}. Pure-R off the COMPACT S1 target (no 612MB load).
  tar_target(trajectory_progression, run_trajectory_progression(microglia_trajectory), format = "qs"),

  # S3: per-cell glmmTMB beta-GLMM sensitivity for the tau:amyloid interaction on bounded pseudotime
  # (pt01 ~ tau*amyloid + batch + (1|unit)). SUPPORTIVE (16 clusters -> weak asymptotics); degrades
  # singular RE -> rank-normal LMM -> a RECORDED method="failed", never blocking the limma-summary
  # primary. Reads the COMPACT S1 per-cell frame (not the 612MB Seurat), INDEPENDENT of trajectory_progression.
  tar_target(trajectory_glmm_sensitivity, glmmtmb_pt_sensitivity(microglia_trajectory$cell_frame), format = "qs"),

  # S4: compact report-data extraction. Bundles the per-cell pseudotime frame + the interaction /
  # decomposition tables + the glmmTMB supportive row from the three (already compact) trajectory
  # targets into one small object, so _trajectory.qmd (and every force-rendered gate run) tar_loads
  # a single compact target -- no 612MB Seurat, no heavy re-read (all three inputs are compact).
  tar_target(trajectory_report, trajectory_report_data(microglia_trajectory, trajectory_progression,
                                                       trajectory_glmm_sensitivity), format = "qs"),
  tar_target(trajectory_figures,
             trajectory_figure_data(trajectory_report),
             format = "qs"),

  # --- cross-tau amyloid-response (four-method logFC scatter) ---
  # Primary per-contrast DE for the three non-snRNAseq modalities (R/modality_de.R): GeoMx
  # voom+TMM with a slide fixed effect + bio-unit duplicateCorrelation; 24M bulk proteome and
  # phosphosite limma-trend on log2 median-normalised, prevalence-filtered intensities (no batch).
  # Each returns per-contrast topTables with logFC keyed by the 5 canonical contrasts.
  tar_target(geomx_de,        run_geomx_de(geomx),                       format = "qs"),
  tar_target(proteome_de_24m, run_proteome_de_24m(proteomics, sample_key), format = "qs"),
  tar_target(phospho_de_24m,  run_phospho_de_24m(phospho, sample_key),   format = "qs"),

  # Compact per-modality amyloid-response logFC pairs: y = nlgf_in_maptki (amyloid on the tau-KO
  # background), x = nlgf_in_p301s (amyloid on the mutant-tau background). The interaction contrast is
  # x - y (nlgf_in_p301s - nlgf_in_maptki), so off-diagonal distance ranks its magnitude. Also carries
  # Functional-category scores over within-method Q99 off-diagonal genes/proteins from the
  # amyloid-effect scatter (phosphoproteomics scatter points are parent-protein means; unmapped
  # symbols kept in fallback categories).
  # Reads only the compact DE topTables (no heavy object).
  tar_target(modality_scatter_figures,
             modality_logfc_scatter_data(pb_de_microglia, symbol_map, geomx_de,
                                         proteome_de_24m, phospho_de_24m),
             format = "qs"),

  # Standalone HTML report render. Source-file targets make report invalidation explicit so
  # caption-only post-render repair can run inside the same `report` target. The render still
  # depends on all compact qmd inputs declared below, and quiet=FALSE keeps Quarto/Pandoc
  # warnings in the gate log. repair_embedded_lightbox() rewrites Quarto's local lightbox
  # PNG hrefs to the already embedded image data URIs, preserving the single offline HTML.
  tar_target(
    report_sources,
    c("_quarto.yml", "index.qmd", "_microglia.qmd", "_trajectory.qmd", "_modality.qmd",
      sort(list.files("R", pattern = "\\.R$", full.names = TRUE))),
    format = "file"
  ),
  tar_target(
    report_extra_files,
    c("theme.scss", "assets/code-tools-fix.html",
      list.files("assets/fonts", pattern = "\\.woff2$", full.names = TRUE)),
    format = "file"
  ),
  tar_target(
    report,
    render_report(
      report_sources = report_sources,
      report_extra_files = report_extra_files,
      microglia_report = microglia_report,
      microglia_figures = microglia_figures,
      trajectory_figures = trajectory_figures,
      modality_scatter_figures = modality_scatter_figures
    ),
    format = "file"
  )
)
