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
  R/ pure fns (S2): constants.R (genotype_levels, contrast_definitions, marker lists,
      rbc_marker_symbols, data_paths) | utils.R (`%||%`, write_tsv_safe) | io.R (loaders) | spine.R
   + (S3) design.R: factorial_design (treatment ~tau+nlgf+tau_nlgf[+batch]) + make_contrast_matrix
      (cell-means ~0+genotype) -> the 5 canonical contrasts; two equivalent parameterisations |
      de_pb.R: pseudobulk_counts/build_pseudobulk (replicate=genotype_batch; `cells=` -> per-substate
      subset pre-aggregation), fit_limma_voom (voomWQW default + confint) / fit_limma_log (log-intensity),
      median_normalise, prevalence_filter. S3 = machinery only; P1-S4 wires the DE targets.
   + (S4) plot.R: theme_tau (ggplot base theme; base_family="" -> device font, warning-free) +
      default ggplot scale_colour/fill_genotype (+ scale_color_ alias; limits/breaks=genotype_levels, drop=FALSE) +
      RWB continuous heatmap scales (`scale_fill_rwb`, `scale_colour_rwb`; signed panels pass midpoint=0) +
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
   + (P2-S2a/b) trajectory.R: run_trajectory_progression (orchestrator -> trajectory_progression) -- collapses
      microglia_trajectory on-lineage pt to 16 genotype_batch summaries (pseudotime_per_replicate) -> factorial_design
      (9 resid df) -> weighted/OLS/bounded limma interaction fits (fit_trajectory_contrasts) + exact 3-channel Kitagawa
      shift-share (decompose_progression_vs_composition -> kitagawa_channels) + Freedman-Lane perm null
      (freedman_lane_interaction; WLS-as-OLS, pivot-free chol2inv, RNG-pure, SENSITIVITY not nominal-exact);
      PRE-REGISTERED primary BH {progression_cf, within_homeostatic} vs separate exploratory BH; mean_pt FLAGGED
      composition-conflated. S2a estimation primitives (pure, crossing-agnostic): derive_batch | pseudotime_per_replicate
      | ordinary_t_table | fit_trajectory_contrasts | kitagawa_channels | within_state_col (col-name sanitizer S2a
      writes / S2b reads). Pure-R, NO new dep; reads the COMPACT S1 target (never the 612MB Seurat).
   + (P2-S3) trajectory.R: glmmtmb_pt_sensitivity (-> trajectory_glmm_sensitivity) -- per-cell beta GLMM
      pt01 ~ tau*amyloid + batch + (1|unit) on the COMPACT S1 cell_frame; SUPPORTIVE (composition-conflated, like
      mean_pt -> corroborates the position shift, NOT progression). Degrade cascade: nonestimable / singular RE /
      !pdHess / non-convergence / degenerate se|z|p -> rank-normal LMM (same gate) -> RECORDED method="failed"
      (NA + fail_reason). A FIT failure never throws; MALFORMED INPUT (missing cols / boundary pt01 / unknown
      genotype / broken genotype_batch) fails LOUD. Records n_units (asymptotics basis, NOT asserted) + fail_reason.
      Helpers: .capture_quietly (muffle+record warnings AND messages, sccomp lesson) + .fit_health_ok (PURE
      unit-tested battery gate: pdHess & converged & finite est & se>0 & finite z & valid p & non-singular) +
      .fit_pt_interaction (fit + Wald-row-by-POSITION + positional-integrity guard [z==est/se, p==2*pnorm(-|z|)]).
      rproject.toml += glmmTMB (CRAN -> P3M trixie BINARY, ABI-clean; namespace-qualified, NOT in tar_option_set packages).
   + (P2-S4a) trajectory.R: trajectory_report_data (-> trajectory_report target) -- bundles the 3 COMPACT
      trajectory targets (microglia_trajectory / trajectory_progression / trajectory_glmm_sensitivity, NEVER the
      612MB Seurat) into ONE ~0.34MB object: slim per-cell plotting frame + interaction table (primary+exploratory
      families; coef/CI/perm_p/FDR -- the comp_cf/progression_cf/cross decomposition channels are ROWS here) +
      per_unit + lineage_per_unit + sensitivity + glmm 13-name subset + provenance (incl. the 3 decomposition
      loadings) -> keeps _trajectory.qmd (+ every gate force-render) compact. NO separate `decomposition` field
      (S4b, codex 955): the qmd reads loadings from provenance + per-channel coefs from interaction rows -> a
      decomposition bundle would only duplicate those two live sources.
      Two guard layers: up-front stopifnot validates the 3 INPUTS' schema; render-cleanliness POSTCONDITIONS
      validate the ASSEMBLED bundle -- interaction col-EXISTENCE asserted BEFORE is.finite (a dropped rbind-sourced
      col fails loud, not vacuously via all(is.finite(NULL))==TRUE) + measure uniqueness + finite coef/ci/p/fdr/perm_p;
      weighted mean_pt coef/CI/p_value present+finite on all 5 canonical contrasts; per_unit/sensitivity non-empty +
      finite; glmm method-enum + finite estimate/CI/p/re_sd; provenance fin1/int1/str1 scalars + logical concordant;
      labelled genotype/substate. Pure: no RNG/IO. (_trajectory.qmd render layer = the report include chain below, S4b DONE.)
   + (P3-S1) mechanism.R: pure mechanism-contract helpers (NO source-time package loads). add_symbol_to_top /
      extract_rank_matrix convert pseudobulk topTables to symbol x contrast statistic matrices with duplicate symbols
      collapsed by max |stat|; run_decoupler_matrix wraps decoupleR::run_ulm into canonical
      {statistic,source,condition,score,p_value}; set_mechanism_prior_cache pins OmniPath cache under
      storage/cache/omnipath and sets TZ=UTC before OmnipathR load (container /etc/localtime warning guard);
      prior_fingerprint hashes sorted prior table + query args + package versions; mechanism_prior_expectations +
      assert_mechanism_prior_expectations pin S1 live prior/coverage drift gates; load_collectri_mouse /
      load_omnipath_ksn_mouse default to official OmniPath REST endpoints (cached TSV) because OmnipathR 4.0.0
      postprocessing currently fails on live schema missing ncbi_tax_id, while try_package=TRUE remains an explicit
      drift probe; standardise_collectri_table / standardise_ksn_table validate direct mouse source-target-mor priors
      (consensus signs, unresolved ambiguous rows dropped, conflicting duplicate signed pairs dropped);
      phospho_site_ids builds `SYMBOL_AApos` from {PG.Genes, PTM.SiteAA, PTM.SiteLocation} with drop counts;
      ksn_coverage_probe reports matched sites, minsize-pass kinases, Gsk3b coverage, and source-case diagnostics.
   + (P3-S2) mechanism.R: RNA mechanism orchestration. collect_rna_rank_matrices converts `pb_de_microglia` +
      fit `pb_de_substate` topTables to symbol x contrast t-stat matrices (whole_microglia + Homeostatic + DAM;
      skipped IFN/Proliferative metadata carried). build_mechanism_gene_sets pulls native mouse MSigDB GO BP/CC/MF
      via msigdbr plus project sets (DAM/Homeostatic/MHC_APC/IFN/NF-kB union + activated/repressed CollecTRI
      targets) and pins a sorted gene-set payload hash. run_mechanism_tf -> decoupleR ULM activity, FDR BH within
      population x contrast. run_mechanism_pathway -> fgseaMultilevel NES, FDR/padj within population x collection
      x contrast; GO maxSize=500, project sets uncapped; known sub-1e-10 fgsea p-floor warning captured, any other
      fgsea warning fails. build_nfkb_attenuation -> NF-kB TF family best source after BH + sign-aware target GSEA
      (activated NES, repressed -NES); support gate only whole_microglia interaction, requires concordant negative
      primary rows after BH across {tf_family,target_gsea}; discordant signs -> `discordant`; tau_in_nlgf/substates
      supportive-only.
   + (P3-S3) mechanism.R: phospho/kinase orchestration. match_24m_intensity_columns asserts exactly 16/16
      sample-key-matched 24M phospho columns, 4/genotype, no duplicate stubs. phospho_feature_frame builds traceable
      feature ids `row<idx>|<PTM.CollapseKey>` + biological site ids `SYMBOL_AApos`, recording blank/multi-gene/
      missing-site counts; positive_log2_matrix converts nonpositive intensities to NA + counts them.
      run_phospho_de_24m -> log2 + median-normalise + prevalence filter (2 samples in all 4 genotypes) ->
      factorial_design(add_batch=FALSE) -> limma-trend 5 contrasts + additive run-index sensitivity fit.
      phospho_site_stat_matrix collapses duplicate biological sites per contrast by highest phosphosite probability,
      then max |t|, then original row. run_kinase_activity loads drift-gated KSN, gates coverage, runs decoupleR ULM
      for primary + run-index fits; build_kinase_mechanism_summary keeps significant kinases plus Gsk3b on every
      contrast with run-index support/confounding columns.
   + (P3-S4) mechanism.R: mechanism_report_data (-> mechanism_report target) -- bundles compact S2/S3 mechanism
      results + P1/P2 anchors into one ~26KB report object: project gene-set rows + top GO survey rows + TF
      highlights (top rows plus key Myc/NF-kB/v1 candidates) + NF-kB gate table/verdict + kinase summary/coverage
      + DAM composition interaction + trajectory composition/progression anchors. Guard layer checks required cols
      and build-fatal chapter anchors (Myc whole interaction, 2 NF-kB primary rows, Gsk3b all contrasts, DAM
      composition interaction, trajectory anchors). No heavy Seurat object.
   + (P4-S1) crossmodality.R: GeoMx DE + decon preflight. geomx_count_matrix explicitly reads RNA/counts
      (GeoMx object default is SCT), drops empty genes, records non-integer residues, and only rounds integer-ish
      matrices. geomx_meta aligns AOI metadata (genotype, slide_rep, bio_unit=genotype:bio_rep, roi/SampleID,
      ROI XY, Q3, negative background, nuclei). fit_geomx_de -> edgeR TMM + limma-voom with `~0+genotype+slide`,
      duplicateCorrelation(block=bio_unit), robust eBayes, canonical 5 contrasts; stores unblocked AOI and
      bio-unit-collapsed sensitivities separately. geomx_decon_preflight records SpatialDecon pinned-repo
      availability + Q3/background/nuclei/reference/profile/memory feasibility; no SpatialDecon install/run in P4-S1.
   + (Spatial decon follow-up S1) crossmodality.R: geomx_reference_profile_data -> geomx_reference_profile.
      Heavy-once full-reference loader: read full snrnaseq_file RNA/counts, map ENSMUSG->symbol via symbol_map,
      overlay retained microglia substates by barcode from microglia_annotated, cap cells/class under fixed seed,
      full-library-normalise selected cells, sparse-average expression over GeoMx-overlap genes, and return only
      broad/substate profile matrices + QC.
      Helpers: spatialdecon_package_info (warn/message capture), reference_label_sets, reference_gene_map,
      reference_select_cells, reference_profile_from_counts, profile_condition_number.
   + (Spatial decon follow-up S2) crossmodality.R: run_geomx_decon -> geomx_decon. Reads GeoMx RNA/data
      (Q3-normalised linear expression) via geomx_norm_matrix, broadcasts geomx_background_matrix
      (NegGeoMean / q_norm_qFactors) over genes, and calls SpatialDecon::spatialdecon under .capture_spatialdecon
      for independent broad + substate arms. normalise_spatialdecon_result stores beta, beta-derived proportions,
      unresolved AOI counts, and residual QC; unresolved beta totals block the arm but retain diagnostics.
      assemble_microglia_substate_abundance anchors substate fractions to broad Microglia beta only if both arms fit.
      Live S2 target is warning-clean but blocked: both arms have 4 unresolved AOIs, so no abundance claim yet.
   + (Spatial decon follow-up S3) crossmodality.R: run_geomx_abundance_de -> geomx_abundance_de. Fits
      fit_geomx_abundance_de (log beta + slide fixed + duplicateCorrelation bio-unit block + unblocked sensitivity)
      only for decon arms with status fit; blocked arms return canonical empty 5-contrast top tables plus reasons.
      geomx_spatial_residual_audit joins SpatialDecon residual QC to GeoMx coordinates and emits per-slide
      nearest-neighbour summaries of genotype-residualised AOI RMS residuals (descriptive QC only). Live S3 target
      is warning-clean but abundance-blocked by the same 4 unresolved AOIs; broad/substate residual audit stored.
   + (Spatial decon follow-up S4) crossmodality.R: spatial_decon_report_data -> spatial_decon_report. Compact
      report handoff over geomx_decon + geomx_abundance_de + geomx_reference_profile: reference QC, decon/abundance
      arm summaries, unresolved AOIs, residual-audit summary, nuclei policy, provenance; NO beta matrices. Live
      status is blocked/action attempted because 4 AOIs have beta_total=0; residual QC remains descriptive.
   + (P4-S2) crossmodality.R: bulk proteome + corrected phospho. match_24m_bulk_columns asserts exact 16-run
      sample-key matching for proteome/phospho-style exports. protein_group_features + aggregate_proteome_raw +
      prepare_proteome_24m_matrix sum raw positive intensities to `PG.ProteinGroups` (NO zero-imputation),
      log2/median-normalise/prevalence-filter, and carry gene-symbol/raw-row provenance. run_proteome_de_24m ->
      limma-trend + additive run-index sensitivity. prepare_phospho_corrected_24m_matrix reuses P3 phosphosite
      prep, asserts identical sample order, subtracts matched filtered parent-protein log2 intensity by exact
      `PG.ProteinGroups`, drops no-parent rows with counts, and re-filters. run_phospho_corrected_24m -> limma-
      trend + run-index. bulk_omics_summary_data compacts feature/significance/run-index/anchor coverage for S4/S5.
   + (P4-S3) crossmodality.R: spatial-composition gate + clearance-axis CCC-lite. geomx_q3_scaled_background
      divides negative-probe background by Q3 factor; profile_collinearity records max abs profile correlation;
      geomx_decon now carries blocked/fit beta + residual diagnostics; fit_geomx_abundance_de is the S3
      log-beta abundance DE path (slide fixed + duplicateCorrelation bio-unit block, unblocked sensitivity).
      clearance_axis_data -> clearance_axis now accepts optional spatial_decon_report; without it an earned preflight
      still fails loud, with it the real follow-up status flows downstream. It harmonises measured Apoe/Trem2/App/Cd74/Pros1/Mertk,
      complement, and synaptic anchors across microglia RNA, GeoMx, and bulk layers with a conservative
      pair-support verdict.
   + (P4-S4) crossmodality.R: integration helpers. crossmodality_table_data -> crossmodality_table harmonises
      snRNAseq whole/substate RNA, GeoMx, proteome, raw/corrected phospho, TF activity, and kinase summary rows to
      one symbol x contrast x modality_group x feature_type table, collapsing duplicate RNA/protein/phosphosite
      symbols by best FDR then abs(statistic) while retaining feature/site counts + missingness reasons. Carries
      both modality_group (layer-level evidence) and modality_class (broad count semantics). crossmodality_pathway_data
      selects project sets + top RNA-mechanism GO sets, scores each selected set from ranked modality statistics,
      and summarises n_modalities_present/sig + n_evidence_groups_present/sig + sign consistency. crossmodality_divergence_data
      focuses {interaction,nlgf_in_maptki,nlgf_in_p301s,tau_in_nlgf}, keeps mixed signs explicit, and outputs compact
      symbol/pathway/clearance highlights for S5.
   + (P4-S5) crossmodality.R: crossmodality_report_data -> crossmodality_report compact chapter bundle. Selects
      GeoMx DE counts/top rows, bulk feature/significance/run-index/anchor slices, clearance/decon verdicts, divergence
      symbol/pathway highlights, and axis-level pathway summaries. Guard layer validates every field _crossmodality.qmd
      reads; render target stays small and does not load GeoMx/proteome/phospho or the 10MB evidence table.
   + (Figure expansion S1) figures.R: figure_manifest (25 hyphenated fig-* ids) + compact inline
      figure-data builders: microglia_figure_data / trajectory_figure_data / mechanism_figure_data /
      crossmodality_figure_data. Builders emit qmd-ready slots, finite geom guards, and pre-binned/top-row
      reductions for heavy shapes (whole/substate volcanoes, GeoMx volcanoes, raw-vs-corrected phospho) so
      later qmd chunks tar_load compact figure targets rather than raw/heavy analysis tables.
   + (Figure-caption-only S4) report.R: render_report (-> report target) calls quarto::quarto_render
      with quiet=FALSE, then repair_embedded_lightbox. Repair rewrites Quarto embedded-lightbox anchors from
      absent local `index_files/figure-html/*.png` hrefs to the already embedded data-URI img src values; fails
      loud if a local lightbox href has no embedded image shape.
   + (Prose-to-figures S2) figures.R: visual_reduction_slot_map + visual_slot_coverage + qc_figure_data.
      Adds compact visual-grammar contracts without qmd rewrites: QC slots from already materialised modality
      targets and alias board slots inside the existing chapter figure targets. Coverage test =
      every S1 manifest `figure`/`schematic` slot has a compact source.
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
       qc_figures           <- qc_figure_data(microglia_seurat_raw, geomx, proteomics, phospho, sample_key)  # Prose-to-figures S2 compact QC visual slots; qmd-safe 4.37KB live
  - P1 microglia core (format="qs"; consumes the snRNAseq modality):
       microglia_processed  <- reprocess_microglia(microglia_seurat_raw)  # SCT+pca+harmony+12 clusters@0.4+umap (687MB)
       microglia_annotated  <- annotate_microglia(microglia_processed, symbol_map)  # UCell substates + prune {6,7,8,11}; 23160 cells, 612MB
       composition_results  <- test_composition(microglia_annotated)  # propeller(logit+asin) [+gated sccomp] x 5 contrasts; counts+stats+concordance (small qs)
       pb_de_microglia      <- run_pb_de_microglia(microglia_annotated, symbol_map)  # whole-MG pseudobulk DE x 5 contrasts (voomWQW+stageR) + DAM concordance (4.7MB)
       pb_de_substate       <- run_pb_de_substate(microglia_annotated)  # per-substate DE: Homeo+DAM fit, IFN/Prolif skip (min-cell floor); cell_counts (7MB)
       microglia_report     <- microglia_report_data(microglia_annotated)  # compact report frame (umap+substate+z) + prune/provenance; ~0.5MB (keeps gate render cheap)
       microglia_figures    <- microglia_figure_data(microglia_report, composition_results, pb_de_microglia, pb_de_substate, symbol_map)  # Figure expansion S1 + Prose-to-figures S2: qmd-ready microglia figure slots, pre-binned DE volcanoes, summary/composition aliases; ~430KB qs live
  - P2 interaction trajectory (format="qs"; consumes microglia_annotated):
       microglia_trajectory <- build_activation_trajectory(microglia_annotated)  # slingshot H->D pseudotime; per-cell frame + per-unit omitted-frac + sensitivity + concordance; serialized ~0.8MB / in-mem ~3.3MB
       trajectory_progression <- run_trajectory_progression(microglia_trajectory)  # S2b: 16-unit pseudotime summaries -> weighted/ols/bounded interaction fits + 3-channel Kitagawa decompose + Freedman-Lane null; primary BH {progression_cf, within_homeostatic}; reads COMPACT S1 target (no 612MB load)
       trajectory_glmm_sensitivity <- glmmtmb_pt_sensitivity(microglia_trajectory$cell_frame)  # S3: per-cell beta GLMM tau:amyloid (degrade -> rank-normal LMM -> method="failed"); supportive, INDEPENDENT of trajectory_progression; ~0.3KB
       trajectory_report    <- trajectory_report_data(microglia_trajectory, trajectory_progression, trajectory_glmm_sensitivity)  # S4a: bundles the 3 compact targets -> one ~0.34MB render object (slim cell_frame + interaction table [decomposition channels as rows] + weighted_top + per_unit + lineage_per_unit + sensitivity + glmm row + provenance [incl. decomposition loadings]); two-layer guards (input schema + assembled-bundle postconditions: col-existence-before-finiteness, weighted mean_pt p on all 5 contrasts, glmm/provenance scalars); keeps gate force-render cheap
       trajectory_figures   <- trajectory_figure_data(trajectory_report, composition_results)  # Figure expansion S1 + Prose-to-figures S2: pseudotime density, DAM-fraction join, Kitagawa forest, decomposition/concordance/logic-board slots; ~18KB live
  - P3 mechanism (format="qs"; consumes compact P1 pseudobulk DE + S1 prior helpers + minimal phospho/sample-key layer):
       mechanism_collectri <- load_collectri_mouse() + assert_mechanism_prior_expectations(collectri)  # cached direct-mouse CollecTRI prior; drift-gated
       mechanism_gene_sets <- build_mechanism_gene_sets(mechanism_collectri)  # GO_BP/GO_CC/GO_MF native mouse MSigDB + project sets incl. signed NF-kB target arms; hash-pinned
       mechanism_tf        <- run_mechanism_tf(pb_de_microglia, pb_de_substate, symbol_map, mechanism_collectri)  # decoupleR ULM TF activity x populations x 5 contrasts; skipped substates carried
       mechanism_pathway   <- run_mechanism_pathway(pb_de_microglia, pb_de_substate, symbol_map, mechanism_gene_sets)  # fgsea NES x populations x {GO/project} x contrasts; GO maxSize=500
       nfkb_attenuation    <- build_nfkb_attenuation(mechanism_tf, mechanism_pathway, mechanism_gene_sets)  # sign-aware primary whole-MG interaction gate + supportive rows
       phospho_de_24m      <- run_phospho_de_24m(phospho, sample_key)  # minimal 24M bulk phosphosite limma-trend DE; no batch; run-index sensitivity stored
       kinase_activity     <- run_kinase_activity(phospho_de_24m)  # decoupleR ULM over direct-mouse KSN; KSN coverage gate + primary/run-index activity tables
       kinase_mechanism_summary <- build_kinase_mechanism_summary(kinase_activity)  # significant kinases + explicit Gsk3b rows with run-index support flags
       mechanism_report    <- mechanism_report_data(mechanism_tf, mechanism_pathway, nfkb_attenuation, kinase_mechanism_summary, composition_results, trajectory_report)  # S4: one compact report object (~26KB) for _mechanism.qmd; no heavy Seurat
       mechanism_figures   <- mechanism_figure_data(mechanism_report)  # Figure expansion S1 + Prose-to-figures S2: pathway/GO/TF/NF-kB/kinase figure slots + status/project/TF aliases; ~25KB live
  - P4 cross-modality:
       geomx_de <- run_geomx_de(geomx)  # S1: GeoMx RNA/counts DE; slide fixed + duplicateCorrelation bio-unit block primary; unblocked/collapsed sensitivities; decon preflight status/reasons, no SpatialDecon run
       geomx_reference_profile <- geomx_reference_profile_data(snrnaseq_file, microglia_annotated, symbol_map, geomx)  # spatial-decon follow-up S1: full-snRNAseq compact broad/substate reference profile + QC; serialized 1.88MB live
       geomx_decon <- run_geomx_decon(geomx, geomx_reference_profile)  # spatial-decon follow-up S2: SpatialDecon broad/substate beta + proportions + residual QC; live blocked by 4 unresolved AOIs, no abundance claim
       geomx_abundance_de <- run_geomx_abundance_de(geomx_decon, geomx)  # spatial-decon follow-up S3: abundance-DE pass-through + residual audit; live blocked, 5.93KB compact target
       spatial_decon_report <- spatial_decon_report_data(geomx_decon, geomx_abundance_de, geomx_reference_profile)  # S4: tiny compact handoff for report; status blocked/action attempted; no beta matrices
       proteome_de_24m <- run_proteome_de_24m(proteomics, sample_key)  # S2: protein-group bulk proteome limma-trend + run-index; raw positive rows summed before log2
       phospho_corrected_24m <- run_phospho_corrected_24m(phospho, sample_key, proteome_de_24m)  # S2: phosphosite minus matched parent protein, re-filter/refit; raw phospho target reused from P3
       bulk_omics_summary <- bulk_omics_summary_data(proteome_de_24m, phospho_de_24m, phospho_corrected_24m)  # S2 compact feature/FDR/run-index/anchor summary (~23KB)
       clearance_axis <- clearance_axis_data(pb_de_microglia, pb_de_substate, symbol_map, geomx_de, bulk_omics_summary, mechanism_gene_sets, spatial_decon_report)  # compact measured-axis table; carries target-derived SpatialDecon blocked state + CCC-lite pair support
       crossmodality_table <- crossmodality_table_data(pb_de_microglia, pb_de_substate, symbol_map, geomx_de, proteome_de_24m, phospho_de_24m, phospho_corrected_24m, mechanism_tf, kinase_mechanism_summary)  # S4 harmonised symbol evidence table (~10MB; 337k rows live), broad modality_class + layer-level modality_group
       crossmodality_pathway <- crossmodality_pathway_data(crossmodality_table, mechanism_gene_sets, mechanism_pathway)  # S4 selected gene-set x modality-class scoring (~108KB live)
       crossmodality_divergence <- crossmodality_divergence_data(crossmodality_table, crossmodality_pathway, clearance_axis)  # S4 compact divergence summary (~1.9MB live), mixed signs + highlights for S5
       crossmodality_report <- crossmodality_report_data(geomx_de, bulk_omics_summary, clearance_axis, crossmodality_divergence, crossmodality_pathway)  # S5 compact report object (~23KB qs live); _crossmodality.qmd reads this plus crossmodality_figures
       crossmodality_figures <- crossmodality_figure_data(crossmodality_report, geomx_de, bulk_omics_summary, phospho_de_24m, phospho_corrected_24m)  # Figure expansion S1 + Prose-to-figures S2: GeoMx/phospho heavy tables reduced to binned/top-row plot data + status/count aliases; ~60KB qs live
  - report_sources <- c("_quarto.yml", "index.qmd", "_qc.qmd", "_microglia.qmd", "_trajectory.qmd", "_mechanism.qmd", "_crossmodality.qmd")  # file target; explicit qmd invalidation
    report_extra_files <- c("theme.scss", assets/fonts/*.woff2)  # file target; explicit theme/font invalidation
    `report` <- render_report(report_sources, report_extra_files, qc_figures, microglia_report, composition_results, pb_de_microglia, pb_de_substate, symbol_map, microglia_figures, trajectory_report, trajectory_figures, mechanism_report, mechanism_figures, crossmodality_report, crossmodality_figures)  # ONE offline HTML; quarto_render quiet=FALSE -> Quarto/Pandoc warnings reach the gate log; post-render repairs embedded-lightbox hrefs to data URIs
       reads `_quarto.yml` (type default; render index.qmd; output _report/; lang en-GB; freeze false)
            -> `index.qmd` (format html, embed-resources, lightbox=auto, theme=theme.scss; no prose body;
                no author metadata; execute.echo=false keeps visible path code-free;
                immediately includes report chapters)
                                                          --{{< include >}}--> `_qc.qmd`
               (caption-only QC chapter: setup `options(warn=2)` -> chunk warnings fail the render;
                tar_load `qc_figures` only -> modality/GeoMx/sample-key grid, 16-cell genotype-batch
                heatmap, depth/fraction histograms, hidden compact-contract checks, metric bounds status)
                                                          --{{< include >}}--> `_microglia.qmd`
               (caption-only microglia chapter: setup `options(warn=2)`; tar_load microglia_report +
                composition_results + pb_de_microglia + pb_de_substate + symbol_map + microglia_figures ->
                summary board, substate UMAPs/score maps, composition shift + unit composition, composition
                forest/concordance, amyloid + all-contrast volcanoes, substate fit/DE audit, pruning audit.
                Hidden finite-check replaces visible sccomp stdout. Claims remain source-derived:
                amyloid->DAM, DAM composition interaction, interaction DE under-powered not absent.)
                                                          --{{< include >}}--> `_trajectory.qmd`
               (caption-only trajectory chapter, {#sec-trajectory}: setup `options(warn=2)`; tar_load
                trajectory_report + trajectory_figures [compact targets] -> logic board, pseudotime shift/density,
                3-channel decomposition, unit DAM-fraction/mean-pt scatter, Kitagawa forest, score-axis
                concordance, robustness/omission audit, method-status board. Headline = synergy adds DAM cells;
                progression and within-homeostatic advance are explicitly not supported; per-cell arm is
                position-only support.)
                                                          --{{< include >}}--> `_mechanism.qmd`
               (caption-only mechanism chapter, {#sec-mechanism}: setup `options(warn=2)`; tar_load
                mechanism_report + mechanism_figures [compact targets] -> mechanism status board, project/GO
                pathway panels, Myc/NF-kB-family TF panels, NF-kB discordance tile, kinase/run-index heatmap.
                Live read = Myc supported, NF-kB discordant/not supported, Gsk3b not recovered; kinase caveat =
                24M bulk hippocampus, not microglia-sorted, genotype-blocked run order.)
                                                          --{{< include >}}--> `_crossmodality.qmd`
               (caption-only cross-modality chapter, {#sec-crossmodality}: setup `options(warn=2)`;
                tar_load crossmodality_report + crossmodality_figures [compact targets] -> status board,
                GeoMx counts/volcano/sensitivity, 24M bulk proteome/phospho counts + run-index + correction,
                anchor heatmap, clearance grid, symbol matrix, pathway heatmap. Modality wording keeps bulk
                hippocampus != microglia-sorted, GeoMx AOIs repeated, SpatialDecon attempted but blocked by
                unresolved AOIs, and CCC-lite != full CCC.)
       `theme.scss` = crimson colours (#B0344D) + IBM Plex (9 woff2 in assets/fonts/, base64-inlined offline)
       Figure labels: every captioned figure chunk uses a hyphenated `fig-*` id. Last rendered HTML QA:
       Figure-caption-only S4, 2026-07-03 (strict rendered main path pass; 48 figures / 48 captions /
       48 nonblank alts; 48 data-URI lightbox hrefs; 0 local figure refs, duplicate IDs, visible body
       paragraphs/tables/stdout/text outputs, code UI blocks, or warning/error markers; full gate green
       after forced 109-chunk render; tar_meta clean across 52 current targets/branches).

### Report prose inventory (Prose-to-figures S1)
`scripts/prose_inventory.py` (stdlib Python; no env deps, non-DAG utility):
  - parses `index.qmd` + `_*.qmd`; skips YAML, executable code bodies, hidden
    setup/helper code, and ordinary source comments.
  - keeps Quarto caption metadata (`fig-cap` / `tbl-cap` / `fig-alt`) as
    report-facing text.
  - emits chapter baseline summary + block-level TSV manifest with qmd, block id,
    line, kind, section, word count, disposition, target slot, label, and text.
  - S1 command:
    `python3 scripts/prose_inventory.py --manifest .agent/prose_replacement_manifest.tsv --summary-only`
  - Figure-caption-only S1 strict gate:
    `uv run python scripts/prose_inventory.py --strict --html _report/index.html --summary-only`
    -> expected red until conversion is complete; source allows only headings/captions,
       rendered HTML fails visible body paragraphs, tables, text-only `.cell-output-display`,
       and stdout provenance.
    -> baseline 5,111 words / 119 blocks; 33 headings kept as navigation;
    86/86 prose/caption blocks assigned non-keep dispositions; selected target
    >=55% reduction (<=2,300 counted words), stretch <=1,800.
  - S5 final: same command -> 1,164 words / 117 blocks (77% reduction, stretch
    met); manifest remains the compact source of qmd block/disposition/slot state.

### Tests (S3; gate-wired at S5)
`tests/test_*.R` each: source the R/ files it exercises + `tests/helpers.R` (expect_error,
make_meta16, make_fake_seurat = synthetic Seurat fixtures, make_trajectory_embedding = synthetic
slingshot embedding), run stopifnot checks (fail-loud, no testthat dep), print `ok - <name>`. Run
from project root: `Rscript tests/test_<x>.R`.
  - test_design.R : 5-contrast exact weights + factorial==cell-means equivalence (property)
  - test_composition.R : composition_counts shapes/empty-drop/constancy-guard + propeller direction (logit+asin) + balance-guard + concordance (incl. completeness fail-loud) + sccomp-gate logical
  - test_crossmodality.R : (P4-S1) GeoMx RNA/count extraction + meta alignment + slide rank guard +
                    duplicateCorrelation primary + unblocked/collapsed sensitivity status +
                    malformed metadata + decon preflight defer/block/earned reasons; (spatial-decon S2/S3)
                    background broadcast + SpatialDecon-result normalisation + unresolved-AOI blocked diagnostics +
                    warning/message capture + two-stage assembly + abundance-DE earned/blocked orchestration +
                    canonical empty top tables + residualised nearest-neighbour spatial audit; (P4-S2) 16-run bulk matching +
                    protein aggregation/no-imputation + parent-protein correction/sample-order guard +
                    missing-parent counts + run-index summary + duplicate/multi-gene provenance; (P4-S3) Q3
                    background scaling + profile-collinearity + abundance-DE design shape + clearance-axis
                    earned/not-earned classification + fail-loud decon-earned guard + spatial_decon_report compact
                    blocked-state handoff; (P4-S4) duplicate RNA/protein/
                    phosphosite collapse + missingness + canonical-contrast guard + pathway modality-score invariants
                    + divergence mixed-sign / clearance-highlight checks; (P4-S5) compact crossmodality_report_data
                    bundle structure + malformed-input schema guards
  - test_microglia.R : reprocess/annotate pure-helper + synthetic-Seurat fixtures (S1/S2) + microglia_report_data extract/guards (S5)
  - test_de_pb.R  : pseudobulk -> 16 cols, median/prevalence, fit_limma_voom/log smokes (S3) + cells= subset,
                    de_pseudobulk/stageR matrix/interaction MDE, run_pb_de_substate fit-or-skip, dam_direction (S4)
  - test_io.R     : io contract tests (pure helpers + loader fail-loud asserts on tempfiles)
  - test_plot.R   : device-free -- theme_tau/scale_*_genotype/scale_*_rwb/concordance_plot class + wiring checks
  - test_figures.R : (Figure expansion S1) 26-figure manifest + compact figure builder contracts for microglia /
                    trajectory / mechanism / cross-modality; (Prose-to-figures S2) QC/report-visual builders,
                    per-chapter board aliases, S1 manifest figure/schematic slot coverage; synthetic finite guards
                    for qmd geom inputs
  - test_report.R  : (Figure-caption-only S4) repair_embedded_lightbox rewrite/no-op/fail-loud cases for
                    embedded Quarto lightbox anchors
  - test_trajectory.R : (P2-S1) score-axis/squeeze/concordance + run_slingshot_lineage (single H->D + branched DAM-terminal) + provenance; (P2-S2a) derive_batch + per-replicate summary + contrast fit + Kitagawa exact-pure reconstruction; (P2-S2b) freedman_lane_interaction (null/signal/determinism/RNG-purity/weighted) + run_trajectory_progression structure on the jitter>0 non-additive fixture [sources R/de_pb.R for assert_complete_crossing]; (P2-S3) .fit_health_ok degrade-gate branches + glmmtmb_pt_sensitivity: beta success (glmmTMB_beta, n_units=16) / singular->failed / non-estimable->failed (fail_reason + captured messages) / unknown-genotype fail-loud; (P2-S4a) trajectory_report_data field/measure/contrast presence + finite interaction inference on the jitter>0 fixture (per-unit DAM composition VARIED so comp_cf/cross fdr stay non-degenerate) + a positive finite-bundle assertion + 13 malformed-input expect_error cases (up-front schema + assembled-bundle postconditions: col-drop/measure/weighted-p/per_unit/sensitivity/glmm-enum/provenance-scalar/NaN-endpoint; fixed=TRUE patterns match the TRUNCATED stopifnot deparse prefix)
  - test_mechanism.R  : (P3-S1) symbol mapping/drop counts + duplicate max-|stat| collapse + decoupleR ULM schema +
                    cache/TZ preflight + prior fingerprint determinism + synthetic CollecTRI/KSN standardisers +
                    phosphosite ID/drop counts + KSN coverage/Gsk3b assertions + prior drift assertion helper;
                    (P3-S2) RNA rank-matrix population/skipped propagation + gene-set filtering/sign splitting +
                    NF-kB source extraction + TF/pathway FDR/direction conventions + p-floor propagation +
                    tau_in_nlgf-supportive-only / discordant / supported attenuation rules; (P3-S3) 16-column
                    phospho match + log2 nonpositive guard + site-id/drop counts + add_batch=FALSE design +
                    run-index sensitivity design + duplicate biological-site collapse + KSN coverage gate +
                    explicit Gsk3b summary carry-through; (P3-S4) mechanism_report_data compact bundle shape +
                    top-row capping + no heavy cell_frame + fail-loud dropped Myc/Gsk3b/trajectory anchors

### Quality gate (S5; review-hardened)
`scripts/check.sh` (fail-loud, `set -euo pipefail`; `CHECK_SKIP_SYNC=1` skips env sync):
  1. `rv sync` + `uv sync`  2. loop `tests/test_*.R` (each `options(warn=2)` -> stray warnings = errors)
  3. FORCE-render report (`tar_invalidate(any_of("report")); tar_make()`, tee'd to a log, `if !`-wrapped) ->
     two render-time catches in the SOURCES: every included `_*.qmd` setup `options(warn=2)`
     (R chunk warning -> render error) + `render_report()` uses `quarto_render(quiet=FALSE)` (Quarto/Pandoc `[WARNING]` -> log)
  4. enforce zero-fault: (a) `tar_meta(error,warnings)` all-NA scoped to manifest names + dynamic branches;
     (b) anchored render-log grep (`^[WARNING]`/`^Warning:`/...), exit 0=fault / 1=clean / >=2=infra.
Any error/warning/log-hit -> non-zero exit. Detail (force-render rationale, warn=2, scoping, anchored grep,
negative tests) -> memory.md Quality gate.

### Config: tracked vs regenerated
tracked : rproject.toml rv.lock | pyproject.toml uv.lock .python-version | _targets.R R/*.R tests/*.R |
          _quarto.yml index.qmd _qc.qmd _microglia.qmd _trajectory.qmd _mechanism.qmd _crossmodality.qmd theme.scss assets/fonts/*.woff2 | .Rprofile rv/scripts/*.R
          rv/.gitignore | scripts/install-*.sh scripts/prose_inventory.py | AGENTS.md .agents/skills/** .codex/prompts/*.md
regen   : rv/library _targets/ _report/ _freeze/ .quarto/ .venv tools/  (gitignored + read-economy skip);
          sccomp_draws_files/ (sccomp per-chain CSV draws at build CWD; gitignored)
