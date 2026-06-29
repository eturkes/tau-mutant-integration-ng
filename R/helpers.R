# Shared utilities for tau-mutant microglia integration analysis.
# All paths are resolved relative to the project root; the Rmd sets
# knitr::opts_knit$set(root.dir = project_root) before sourcing.
#
# This file is a thin loader: project-wide library() calls plus
# source()s of every domain file in dependency order. Add new R/*.R
# files to the source order below as the project grows.

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(readr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(edgeR)
  library(limma)
  library(fgsea)
  library(msigdbr)
  library(ComplexHeatmap)
  library(circlize)
  library(BiocParallel)
  library(nebula)
  library(nichenetr)  # convert_human_to_mouse_symbols uses lazy-data
                      # geneinfo_human that only resolves when the
                      # package is attached (not via ::). Needed by
                      # R/kinase_inference.R::map_ksn_human_to_mouse.
})

# Source order is dependency-driven, not alphabetical:
#   constants.R       must come first (every other file references at
#                     least one symbol defined here).
#   utils.R           defines cache_or_run, %||%, write_tsv_safe, isTRUE_vec.
#   io.R              defines symbols_to_ensembl (used by microglia.R).
#   design.R          defines factorial_design (used by hdwgcna.R,
#                     de_pb.R-callers).
#   de_pb.R           pseudobulk + bulk DE methods.
#   de_sc.R           single-cell DE (NEBULA, glmmTMB).
#   fgsea.R           fgsea / pathway analysis helpers + collection
#                     getters.
#   pathway_survey.R  agnostic cross-modality pathway-ranking helpers
#                     consumed by rmd/12_pathway_survey.Rmd (depends on
#                     the fgsea cache shape produced by fgsea.R).
#   tf_inference.R    transcription-factor activity inference helpers
#                     used by the Phase E mechanism layer
#                     (storage/notes/mechanism_layer_plan.md, sessions
#                     E2..E5). Wraps decoupleR + CollecTRI mouse prior;
#                     consumes the same DE cache shapes as fgsea.R, so
#                     it is sourced after both fgsea.R and
#                     pathway_survey.R.
#   kinase_inference.R kinase activity inference helpers used by the
#                     Phase F mechanism layer (mechanism_layer_plan.md,
#                     sessions F2..F5). Wraps decoupleR + OmniPath KSN
#                     (mouse-mapped via nichenetr); consumes the
#                     phospho DE caches (de_phospho.rds and
#                     de_phospho_corrected.rds). Sourced after
#                     tf_inference.R to keep the mechanism-layer
#                     files grouped, and because the kinase
#                     module-header references the TF module-header
#                     conventions for reading the Phase F outputs
#                     against the Phase E outputs side-by-side.
#                     Self-contained: does not call any function from
#                     tf_inference.R; the cross-reference is
#                     documentation-only.
#   ccc_inference.R   cell-cell communication inference helpers used
#                     by the Phase G mechanism layer
#                     (mechanism_layer_plan.md, sessions G2..G5).
#                     Wraps the Python LIANA+ framework (Dimitrov 2024)
#                     via reticulate against a project-local virtual
#                     env at `.venv/bin/python`. Consumes the project-
#                     wide Seurat object (seurat_full_processed.rds);
#                     produces a per-genotype + per-contrast cache
#                     (storage/cache/liana_output.rds) intended to
#                     rbind cleanly with the CellChat and MultiNicheNet
#                     per-contrast tibbles at G3. Sourced after
#                     kinase_inference.R to keep the three mechanism-
#                     layer modules grouped in plan-phase order
#                     (E -> F -> G). Self-contained: does not call any
#                     function from tf_inference.R or
#                     kinase_inference.R; the cross-reference is
#                     documentation-only.
#   causal_network.R  causal signalling-network reconstruction helpers
#                     (narrative arc J; causal_network_plan.md, J2..J5).
#                     Assembles the OmniPath signed-directed PKN and the
#                     CARNIVAL inputs that wire the Phase F kinase
#                     UPSTREAM layer to the Phase E TF DOWNSTREAM layer.
#                     Reuses tf_inference.R (.extract_tf_per_modality,
#                     build_axis_gene_universe), so it is sourced after
#                     the three E/F/G mechanism modules to sit as their
#                     capstone. Reads the kinase cache directly (does
#                     not call kinase_inference.R).
#   scenic.R          data-driven SCENIC regulon contrast modelling
#                     (narrative arc K; scenic_regulons_plan.md, K4).
#                     Reuses tf_inference.R (extract_de_stat_matrix,
#                     run_decoupler_per_modality, split_decoupler_by_contrast)
#                     for the controlled CollecTRI->SCENIC network swap and
#                     design.R/de_pb.R for the AUCell pseudobulk factorial,
#                     so it is sourced after causal_network.R.
#   spatial_decon.R   GeoMx tissue-composition deconvolution (narrative
#                     arc L; spatial_deconvolution_plan.md, L2). Wraps
#                     SpatialDecon (create_profile_matrix + spatialdecon)
#                     to turn the snRNAseq reference + the GeoMx-WTA object
#                     into per-ROI cell-type / microglial-substate
#                     abundances, then reuses design.R (factorial_design)
#                     + de_pb.R (fit_limma_log) for the locked 5-contrast
#                     factorial on log-abundance. Self-contained otherwise;
#                     sourced after scenic.R as the next display-only arc.
#   trajectory.R      microglial activation-trajectory / dynamics layer
#                     (narrative arc M; trajectory_plan.md, M2). Wraps
#                     Slingshot (pseudotime on the harmony embedding, the 4
#                     `state`s as clusters, root = homeostatic) + a root-free
#                     potency cross-check, summarises pseudotime per
#                     genotype_batch replicate, and reuses design.R
#                     (factorial_design) + de_pb.R (fit_limma_log) for the
#                     locked 5-contrast factorial on that summary plus a
#                     progression-vs-composition decomposition. The project's
#                     first DYNAMICS readout (all prior arcs are static).
#                     Self-contained; sourced after spatial_decon.R as the
#                     next display-only arc.
#   celltype_specificity.R  cross-cell-type specificity layer (narrative arc
#                     N; celltype_specificity_plan.md, N2). Re-runs the EXACT
#                     microglia DE treatment (subset pseudobulk limma-voom +
#                     single-cell NEBULA) symmetrically over 6 units (5 broad
#                     non-microglial types + microglia-as-whole) under NATIVE
#                     and power-MATCHED (balanced-downsample) regimes to ask
#                     whether the tau x amyloid interaction is microglia-
#                     specific. Reuses de_pb.R (build_pseudobulk,
#                     fit_limma_voom), de_sc.R (fit_nebula_microglia), and
#                     design.R (make_contrast_matrix); the project's first
#                     CROSS-CELL-TYPE readout. Self-contained; sourced after
#                     trajectory.R as the next display-only arc.
#   microglia.R       substate annotation + RBC contaminant flag.
#   hdwgcna.R         hdWGCNA pipeline + module DE.
#   plot.R            ggplot helpers (uses %||% from utils.R).
#   report.R          shared bslib HTML report theme (report_theme()); also
#                     wired self-contained into each report's YAML, so it must
#                     stand alone at parse time (no cross-file deps).
source("R/constants.R")
source("R/utils.R")
source("R/io.R")
source("R/design.R")
source("R/de_pb.R")
source("R/de_sc.R")
source("R/fgsea.R")
source("R/pathway_survey.R")
source("R/tf_inference.R")
source("R/kinase_inference.R")
source("R/ccc_inference.R")
source("R/causal_network.R")
source("R/scenic.R")
source("R/spatial_decon.R")
source("R/trajectory.R")
source("R/celltype_specificity.R")
source("R/microglia.R")
source("R/hdwgcna.R")
source("R/plot.R")
source("R/report.R")
