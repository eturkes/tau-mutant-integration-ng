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
  tar_target(microglia_figures,
             microglia_figure_data(microglia_report, composition_results,
                                   pb_de_microglia, pb_de_substate, symbol_map),
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
             trajectory_figure_data(trajectory_report, composition_results),
             format = "qs"),

  # --- P3 mechanism ---
  # S2 RNA mechanism targets. Reuse the cached/fingerprinted direct-mouse CollecTRI prior from S1,
  # then derive focused mouse gene sets, decoupleR ULM TF activity, fgsea preranked pathway scores,
  # and the sign-aware predeclared NF-kB attenuation gate. Inputs are the existing replicate-corrected
  # pseudobulk DE targets (whole microglia + fit substates); skipped substates travel as metadata.
  tar_target(mechanism_collectri, {
    net <- load_collectri_mouse()
    assert_mechanism_prior_expectations(collectri = net)
    net
  }, format = "qs"),
  tar_target(mechanism_gene_sets, build_mechanism_gene_sets(mechanism_collectri), format = "qs"),
  tar_target(mechanism_tf, run_mechanism_tf(pb_de_microglia, pb_de_substate, symbol_map,
                                           mechanism_collectri), format = "qs"),
  tar_target(mechanism_pathway, run_mechanism_pathway(pb_de_microglia, pb_de_substate,
                                                     symbol_map, mechanism_gene_sets), format = "qs"),
  tar_target(nfkb_attenuation, build_nfkb_attenuation(mechanism_tf, mechanism_pathway,
                                                     mechanism_gene_sets), format = "qs"),

  tar_target(proteomics,           read_spectronaut_tsv(proteomics_file),  format = "qs"),
  tar_target(phospho,              read_spectronaut_tsv(phospho_file),     format = "qs"),
  tar_target(sample_key,           proteomics_sample_meta(sample_key_file), format = "qs"),
  tar_target(qc_figures,
             qc_figure_data(microglia_seurat_raw, geomx, proteomics, phospho, sample_key),
             format = "qs"),

  # S3 kinase layer. Minimal 24M bulk-phosphosite DE (16 sample-key-matched runs, no batch
  # term) feeds decoupleR ULM activity over the fingerprinted direct-mouse OmniPath KSN. The
  # summary carries significant kinases plus Gsk3b for every contrast, with run-index sensitivity
  # columns so report prose can downgrade run-order-dependent support.
  tar_target(phospho_de_24m, run_phospho_de_24m(phospho, sample_key), format = "qs"),
  tar_target(kinase_activity, run_kinase_activity(phospho_de_24m), format = "qs"),
  tar_target(kinase_mechanism_summary, build_kinase_mechanism_summary(kinase_activity), format = "qs"),

  # S4 report bundle. Selects compact TF/pathway/NF-kB/kinase highlights plus P1/P2
  # anchors for synthesis. No heavy Seurat object is read or stored.
  tar_target(mechanism_report, mechanism_report_data(mechanism_tf, mechanism_pathway,
                                                     nfkb_attenuation, kinase_mechanism_summary,
                                                     composition_results, trajectory_report),
             format = "qs"),
  tar_target(mechanism_figures, mechanism_figure_data(mechanism_report), format = "qs"),

  # --- P4 cross-modality ---
  # S1 GeoMx spatial DE: RNA count layer (explicit, despite the object's SCT default) ->
  # edgeR TMM + limma-voom, slide fixed effect, duplicateCorrelation block = genotype:bio_rep,
  # plus unblocked and bio-unit-collapsed sensitivities. Deconvolution is NOT run here; the
  # preflight records Q3/background/nuclei/reference/package feasibility for the S3 gate.
  tar_target(geomx_de, run_geomx_de(geomx), format = "qs"),

  # Spatial decon follow-up S1: compact full-reference profile. Loads the full snRNAseq RDS
  # once, overlays retained microglia substates from microglia_annotated, caps cells per
  # class before profile averaging, and stores only broad/substate profile matrices + QC.
  tar_target(geomx_reference_profile,
             geomx_reference_profile_data(snrnaseq_file, microglia_annotated, symbol_map, geomx),
             format = "qs"),

  # Spatial decon follow-up S2: run SpatialDecon on GeoMx RNA/data Q3-normalised
  # expression with Q3-scaled negative-probe background. Broad profile is primary;
  # substate fit is gated separately and assembled back onto broad Microglia beta.
  # Nuclei-based absolute count conversion stays disabled while nuclei sentinels remain.
  tar_target(geomx_decon,
             run_geomx_decon(geomx, geomx_reference_profile),
             format = "qs"),

  # Spatial decon follow-up S3: if SpatialDecon earns finite beta, fit log-beta abundance
  # DE on the same GeoMx slide + bio-unit design; otherwise pass through the blocked state
  # while retaining residual spatial diagnostics from the attempted fit.
  tar_target(geomx_abundance_de,
             run_geomx_abundance_de(geomx_decon, geomx),
             format = "qs"),
  tar_target(spatial_decon_report,
             spatial_decon_report_data(geomx_decon, geomx_abundance_de,
                                       geomx_reference_profile),
             format = "qs"),

  # S2 bulk proteome + corrected phospho: 24M 16-run sample-key-matched bulk hippocampus
  # layers. Proteome aggregates raw positive rows to PG.ProteinGroups before log2/median
  # normalisation/limma-trend. Corrected phospho subtracts matched parent-protein log2
  # intensity from the existing phosphosite layer, then refits limma-trend. Summary is compact
  # and report-oriented; P3's raw phospho target is reused, not duplicated.
  tar_target(proteome_de_24m, run_proteome_de_24m(proteomics, sample_key), format = "qs"),
  tar_target(phospho_corrected_24m,
             run_phospho_corrected_24m(phospho, sample_key, proteome_de_24m), format = "qs"),
  tar_target(bulk_omics_summary,
             bulk_omics_summary_data(proteome_de_24m, phospho_de_24m, phospho_corrected_24m),
             format = "qs"),

  # Spatial-composition gate + clearance-axis CCC-lite. The follow-up report bundle carries
  # the actual SpatialDecon attempted/blocked state into cross-modality + synthesis surfaces.
  tar_target(clearance_axis,
             clearance_axis_data(pb_de_microglia, pb_de_substate, symbol_map, geomx_de,
                                 bulk_omics_summary, mechanism_gene_sets,
                                 spatial_decon_report),
             format = "qs"),

  # S4 integrated gene/pathway divergence view. Evidence rows collapse duplicated RNA symbols,
  # protein groups, and phosphosites to one symbol x contrast x modality_group row while retaining
  # feature/site counts. Pathway/divergence leaves nulls and mixed signs explicit for the S5 report.
  tar_target(crossmodality_table,
             crossmodality_table_data(pb_de_microglia, pb_de_substate, symbol_map, geomx_de,
                                      proteome_de_24m, phospho_de_24m,
                                      phospho_corrected_24m, mechanism_tf,
                                      kinase_mechanism_summary),
             format = "qs"),
  tar_target(crossmodality_pathway,
             crossmodality_pathway_data(crossmodality_table, mechanism_gene_sets,
                                        mechanism_pathway),
             format = "qs"),
  tar_target(crossmodality_divergence,
             crossmodality_divergence_data(crossmodality_table, crossmodality_pathway,
                                           clearance_axis),
             format = "qs"),

  # S5 compact report bundle. Selects GeoMx, bulk, clearance, pathway, and divergence
  # slices for the chapter so _crossmodality.qmd loads one small object, never the full
  # GeoMx/proteome/phospho targets or the 10MB harmonised evidence table.
  tar_target(crossmodality_report,
             crossmodality_report_data(geomx_de, bulk_omics_summary, clearance_axis,
                                       crossmodality_divergence, crossmodality_pathway),
             format = "qs"),
  tar_target(crossmodality_figures,
             crossmodality_figure_data(crossmodality_report, geomx_de, bulk_omics_summary,
                                       phospho_de_24m, phospho_corrected_24m,
                                       crossmodality_table),
             format = "qs"),

  # Story plates: compact synthesis figures over already-built compact result bundles.
  # No new inference; this target just pre-assembles publication-grade front-of-report
  # plotting data so the caption-only report can tell the coherent story first.
  tar_target(story_figures,
             story_figure_data(qc_figures, composition_results, pb_de_microglia,
                               trajectory_report, mechanism_report,
                               crossmodality_report, crossmodality_figures),
             format = "qs"),

  # Standalone HTML report render. Source-file targets make report invalidation explicit so
  # caption-only post-render repair can run inside the same `report` target. The render still
  # depends on all compact qmd inputs declared below, and quiet=FALSE keeps Quarto/Pandoc
  # warnings in the gate log. repair_embedded_lightbox() rewrites Quarto's local lightbox
  # PNG hrefs to the already embedded image data URIs, preserving the single offline HTML.
  tar_target(
    report_sources,
    c("_quarto.yml", "index.qmd", "_qc.qmd", "_microglia.qmd",
      "_trajectory.qmd", "_mechanism.qmd", "_crossmodality.qmd"),
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
      qc_figures = qc_figures,
      microglia_report = microglia_report,
      composition_results = composition_results,
      pb_de_microglia = pb_de_microglia,
      pb_de_substate = pb_de_substate,
      symbol_map = symbol_map,
      microglia_figures = microglia_figures,
      trajectory_report = trajectory_report,
      trajectory_figures = trajectory_figures,
      mechanism_report = mechanism_report,
      mechanism_figures = mechanism_figures,
      crossmodality_report = crossmodality_report,
      crossmodality_figures = crossmodality_figures,
      story_figures = story_figures
    ),
    format = "file"
  )
)
