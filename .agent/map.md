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
  R/ pure fns (S2): constants.R (genotype_levels/colours, contrast_definitions, marker lists,
      rbc_marker_symbols, data_paths) | utils.R (`%||%`, write_tsv_safe) | io.R (loaders) | spine.R
   + (S3) design.R: factorial_design (treatment ~tau+nlgf+tau_nlgf[+batch]) + make_contrast_matrix
      (cell-means ~0+genotype) -> the 5 canonical contrasts; two equivalent parameterisations |
      de_pb.R: pseudobulk_counts/build_pseudobulk (replicate=genotype_batch; `cells=` -> per-substate
      subset pre-aggregation), fit_limma_voom (voomWQW default + confint) / fit_limma_log (log-intensity),
      median_normalise, prevalence_filter. S3 = machinery only; P1-S4 wires the DE targets.
   + (S4) plot.R: theme_tau (ggplot base theme; base_family="" -> device font, warning-free) +
      scale_colour/fill_genotype (+ scale_color_ alias; limits/breaks=genotype_levels, drop=FALSE) +
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
      availability + Q3/background/nuclei/reference/profile/memory feasibility; no SpatialDecon install/run in S1.
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
      fit_geomx_abundance_de is the un-wired future SpatialDecon beta/log-abundance DE path (slide fixed +
      duplicateCorrelation bio-unit block, unblocked sensitivity). clearance_axis_data -> clearance_axis target:
      fails loud if geomx_de$decon_preflight becomes `earned` before decon targets exist; otherwise records the
      intentional decon skip and harmonises measured Apoe/Trem2/App/Cd74/Pros1/Mertk, complement, and synaptic
      anchors across microglia RNA, GeoMx, and bulk layers with a conservative pair-support verdict.
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
   + (P5-S1) synthesis.R: synthesis_report_data -> synthesis_report compact closing bundle. Reads only
      microglia_report, trajectory_report, mechanism_report, and crossmodality_report; builds headline strings,
      a 10-row descriptive evidence table, status counts, unsupported/open rows, and tiny source highlights. Guard
      layer checks the synthesis anchors (amyloid->DAM, DAM composition, trajectory comp/prog rows, Myc, NF-kB,
      Gsk3b, GeoMx/bulk caveats, SpatialDecon, clearance pairs incl. empty earned set) and rejects ledger-like
      score columns. `.synthesis_clean_text` normalises raw caveat strings before table output, so the compact
      target itself stays free of stale phase-step wording.
   + (Figure expansion S1) figures.R: figure_manifest (26 hyphenated fig-* ids) + compact inline
      figure-data builders: microglia_figure_data / trajectory_figure_data / mechanism_figure_data /
      crossmodality_figure_data. Builders emit qmd-ready slots, finite geom guards, and pre-binned/top-row
      reductions for heavy shapes (whole/substate volcanoes, GeoMx volcanoes, raw-vs-corrected phospho) so
      later qmd chunks tar_load compact figure targets rather than raw/heavy analysis tables.
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
       composition_results  <- test_composition(microglia_annotated)  # propeller(logit+asin) [+gated sccomp] x 5 contrasts; counts+stats+concordance (small qs)
       pb_de_microglia      <- run_pb_de_microglia(microglia_annotated, symbol_map)  # whole-MG pseudobulk DE x 5 contrasts (voomWQW+stageR) + DAM concordance (4.7MB)
       pb_de_substate       <- run_pb_de_substate(microglia_annotated)  # per-substate DE: Homeo+DAM fit, IFN/Prolif skip (min-cell floor); cell_counts (7MB)
       microglia_report     <- microglia_report_data(microglia_annotated)  # compact report frame (umap+substate+z) + prune/provenance; ~0.5MB (keeps gate render cheap)
       microglia_figures    <- microglia_figure_data(microglia_report, composition_results, pb_de_microglia, pb_de_substate, symbol_map)  # Figure expansion S1: qmd-ready microglia figure slots, pre-binned DE volcanoes; 2.381MB in memory / 420KB qs live
  - P2 interaction trajectory (format="qs"; consumes microglia_annotated):
       microglia_trajectory <- build_activation_trajectory(microglia_annotated)  # slingshot H->D pseudotime; per-cell frame + per-unit omitted-frac + sensitivity + concordance; serialized ~0.8MB / in-mem ~3.3MB
       trajectory_progression <- run_trajectory_progression(microglia_trajectory)  # S2b: 16-unit pseudotime summaries -> weighted/ols/bounded interaction fits + 3-channel Kitagawa decompose + Freedman-Lane null; primary BH {progression_cf, within_homeostatic}; reads COMPACT S1 target (no 612MB load)
       trajectory_glmm_sensitivity <- glmmtmb_pt_sensitivity(microglia_trajectory$cell_frame)  # S3: per-cell beta GLMM tau:amyloid (degrade -> rank-normal LMM -> method="failed"); supportive, INDEPENDENT of trajectory_progression; ~0.3KB
       trajectory_report    <- trajectory_report_data(microglia_trajectory, trajectory_progression, trajectory_glmm_sensitivity)  # S4a: bundles the 3 compact targets -> one ~0.34MB render object (slim cell_frame + interaction table [decomposition channels as rows] + weighted_top + per_unit + lineage_per_unit + sensitivity + glmm row + provenance [incl. decomposition loadings]); two-layer guards (input schema + assembled-bundle postconditions: col-existence-before-finiteness, weighted mean_pt p on all 5 contrasts, glmm/provenance scalars); keeps gate force-render cheap
       trajectory_figures   <- trajectory_figure_data(trajectory_report, composition_results)  # Figure expansion S1: pseudotime density, DAM-fraction join, Kitagawa forest, audit slots; 0.033MB live
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
       mechanism_figures   <- mechanism_figure_data(mechanism_report)  # Figure expansion S1: pathway/GO/TF/NF-kB/kinase figure slots; 0.107MB live
  - P4 cross-modality:
       geomx_de <- run_geomx_de(geomx)  # S1: GeoMx RNA/counts DE; slide fixed + duplicateCorrelation bio-unit block primary; unblocked/collapsed sensitivities; decon preflight status/reasons, no SpatialDecon run
       proteome_de_24m <- run_proteome_de_24m(proteomics, sample_key)  # S2: protein-group bulk proteome limma-trend + run-index; raw positive rows summed before log2
       phospho_corrected_24m <- run_phospho_corrected_24m(phospho, sample_key, proteome_de_24m)  # S2: phosphosite minus matched parent protein, re-filter/refit; raw phospho target reused from P3
       bulk_omics_summary <- bulk_omics_summary_data(proteome_de_24m, phospho_de_24m, phospho_corrected_24m)  # S2 compact feature/FDR/run-index/anchor summary (~23KB)
       clearance_axis <- clearance_axis_data(pb_de_microglia, pb_de_substate, symbol_map, geomx_de, bulk_omics_summary, mechanism_gene_sets)  # S3 compact measured-axis table (~37KB); decon skipped/defer recorded; Apoe_Trem2 currently the only earned CCC-lite pair
       crossmodality_table <- crossmodality_table_data(pb_de_microglia, pb_de_substate, symbol_map, geomx_de, proteome_de_24m, phospho_de_24m, phospho_corrected_24m, mechanism_tf, kinase_mechanism_summary)  # S4 harmonised symbol evidence table (~10MB; 337k rows live), broad modality_class + layer-level modality_group
       crossmodality_pathway <- crossmodality_pathway_data(crossmodality_table, mechanism_gene_sets, mechanism_pathway)  # S4 selected gene-set x modality-class scoring (~108KB live)
       crossmodality_divergence <- crossmodality_divergence_data(crossmodality_table, crossmodality_pathway, clearance_axis)  # S4 compact divergence summary (~1.9MB live), mixed signs + highlights for S5
       crossmodality_report <- crossmodality_report_data(geomx_de, bulk_omics_summary, clearance_axis, crossmodality_divergence, crossmodality_pathway)  # S5 compact report object (~23KB qs live); _crossmodality.qmd reads this plus crossmodality_figures
       crossmodality_figures <- crossmodality_figure_data(crossmodality_report, geomx_de, bulk_omics_summary, phospho_de_24m, phospho_corrected_24m)  # Figure expansion S1: GeoMx/phospho heavy tables reduced to binned/top-row plot data; 0.514MB live
  - P5 synthesis:
       synthesis_report <- synthesis_report_data(microglia_report, trajectory_report, mechanism_report, crossmodality_report)  # S1 compact read-only synthesis (~4.8KB live); no crossmodality_divergence/heavy targets; descriptive status rows only
  - `report` <- tar_quarto(path=".", quiet=FALSE, extra_files=c("theme.scss", assets/fonts/*.woff2))  # ONE offline HTML; quiet=FALSE -> Quarto/Pandoc warnings reach the gate log
       reads `_quarto.yml` (type default; render index.qmd; output _report/; lang en-GB; freeze false)
            -> `index.qmd` (format html, embed-resources, theme=theme.scss) --{{< include >}}--> `_synthesis.qmd`
               (P5 synthesis chapter, {#sec-synthesis}: setup `options(warn=2)`; tar_load synthesis_report [ONE
                compact target] -> answer-first paragraph + status-count figure + S2 claim-source evidence map +
                compact evidence table + unsupported/unearned paragraph; no heavy target reads, no ledger scoring)
                                                          --{{< include >}}--> `_qc.qmd`
               (QC-sanity chapter: setup `options(warn=2)` -> chunk warnings fail the render; tar_load 4
                modalities + sample_key -> dims, 16x16 design bijection, bounds)
                                                          --{{< include >}}--> `_microglia.qmd`
               (P1 microglia chapter: setup `options(warn=2)`; tar_load microglia_report + composition_results +
                pb_de_microglia + pb_de_substate + symbol_map + microglia_figures -> original substate UMAP,
                composition forest/table, amyloid->DAM volcano + S2 inline figures {genotype-faceted UMAP, score
                triptych/distributions, 16-unit composition, concordance grid, all-contrast volcanoes, substate
                fit audit, within-substate DE counts}, under-powered interaction + @sec-trajectory pointer, Thrupp
                + dropout caveats)
                                                          --{{< include >}}--> `_trajectory.qmd`
               (P2 trajectory chapter, {#sec-trajectory}: setup `options(warn=2)`; tar_load trajectory_report +
                trajectory_figures [compact targets] -> pseudotime-shift + composition-not-progression 3-channel
                decomposition + per-cell glmmTMB supportive + score-axis concordance + S3 inline figures
                {pseudotime density, unit DAM-fraction/mean-pt scatter, channel/decomposition forest,
                robustness/omission audit} + 5 caveats/provenance; headline = synergy adds DAM cells, no
                supported further-advance; inference numbers inline-computed from trajectory_report, never hardcoded
                [fixed design constants -- resid df, sensitivity dims -- stated as text])
                                                          --{{< include >}}--> `_mechanism.qmd`
               (P3 mechanism chapter, {#sec-mechanism}: setup `options(warn=2)`; tar_load mechanism_report +
                mechanism_figures [compact targets] -> pathway survey + TF activity + NF-kB attenuation gate +
                Gsk3b/kinase support + S3 inline figures {all-population project pathway heatmap, GO dot plot,
                Myc/NF-kB-family TF lollipop, NF-kB discordance tile, kinase/run-index heatmap} + synthesis/caveats.
                Live read = Myc supported, NF-kB discordant/not supported, Gsk3b not recovered; kinase caveat =
                24M bulk hippocampus, not microglia-sorted, genotype-blocked run order.)
                                                          --{{< include >}}--> `_crossmodality.qmd`
               (P4 cross-modality chapter, {#sec-crossmodality}: setup `options(warn=2)`; tar_load
                crossmodality_report + crossmodality_figures [compact targets] -> GeoMx spatial DE, 24M bulk
                proteome/phospho + run-index caveats, decon skip + clearance-axis CCC-lite, integrated pathway/
                symbol divergence, and S4 inline figures {GeoMx volcanoes, sensitivity/loss, bulk run-index,
                raw-vs-corrected phospho, anchor heatmap, clearance grid, symbol matrix, pathway heatmap}.
                Modality wording keeps bulk hippocampus != microglia-sorted, GeoMx AOIs repeated,
                SpatialDecon skipped/defer, and CCC-lite != full CCC.)
       `theme.scss` = crimson colours (#B0344D) + IBM Plex (9 woff2 in assets/fonts/, base64-inlined offline)

### Tests (S3; gate-wired at S5)
`tests/test_*.R` each: source the R/ files it exercises + `tests/helpers.R` (expect_error,
make_meta16, make_fake_seurat = synthetic Seurat fixtures, make_trajectory_embedding = synthetic
slingshot embedding), run stopifnot checks (fail-loud, no testthat dep), print `ok - <name>`. Run
from project root: `Rscript tests/test_<x>.R`.
  - test_design.R : 5-contrast exact weights + factorial==cell-means equivalence (property)
  - test_composition.R : composition_counts shapes/empty-drop/constancy-guard + propeller direction (logit+asin) + balance-guard + concordance (incl. completeness fail-loud) + sccomp-gate logical
  - test_crossmodality.R : (P4-S1) GeoMx RNA/count extraction + meta alignment + slide rank guard +
                    duplicateCorrelation primary + unblocked/collapsed sensitivity status +
                    malformed metadata + decon preflight defer/block/earned reasons; (P4-S2) 16-run bulk matching +
                    protein aggregation/no-imputation + parent-protein correction/sample-order guard +
                    missing-parent counts + run-index summary + duplicate/multi-gene provenance; (P4-S3) Q3
                    background scaling + profile-collinearity + abundance-DE design shape + clearance-axis
                    earned/not-earned classification + fail-loud decon-earned guard; (P4-S4) duplicate RNA/protein/
                    phosphosite collapse + missingness + canonical-contrast guard + pathway modality-score invariants
                    + divergence mixed-sign / clearance-highlight checks; (P4-S5) compact crossmodality_report_data
                    bundle structure + malformed-input schema guards
  - test_microglia.R : reprocess/annotate pure-helper + synthetic-Seurat fixtures (S1/S2) + microglia_report_data extract/guards (S5)
  - test_de_pb.R  : pseudobulk -> 16 cols, median/prevalence, fit_limma_voom/log smokes (S3) + cells= subset,
                    de_pseudobulk/stageR matrix/interaction MDE, run_pb_de_substate fit-or-skip, dam_direction (S4)
  - test_io.R     : io contract tests (pure helpers + loader fail-loud asserts on tempfiles)
  - test_plot.R   : device-free -- theme_tau/scale_*_genotype/concordance_plot class + wiring checks
  - test_figures.R : (Figure expansion S1) 26-figure manifest + compact figure builder contracts for microglia /
                    trajectory / mechanism / cross-modality; synthetic finite guards for qmd geom inputs
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
  - test_synthesis.R : (P5-S1) synthesis_report_data compact bundle shape + status enum validation +
                    missing-anchor failures + empty earned-pair handling + no-ledger-column invariant

### Quality gate (S5; review-hardened)
`scripts/check.sh` (fail-loud, `set -euo pipefail`; `CHECK_SKIP_SYNC=1` skips env sync):
  1. `rv sync` + `uv sync`  2. loop `tests/test_*.R` (each `options(warn=2)` -> stray warnings = errors)
  3. FORCE-render report (`tar_invalidate(any_of("report")); tar_make()`, tee'd to a log, `if !`-wrapped) ->
     two render-time catches in the SOURCES: every included `_*.qmd` setup `options(warn=2)`
     (R chunk warning -> render error) + `_targets.R` `tar_quarto(quiet=FALSE)` (Quarto/Pandoc `[WARNING]` -> log)
  4. enforce zero-fault: (a) `tar_meta(error,warnings)` all-NA scoped to manifest names + dynamic branches;
     (b) anchored render-log grep (`^[WARNING]`/`^Warning:`/...), exit 0=fault / 1=clean / >=2=infra.
Any error/warning/log-hit -> non-zero exit. Detail (force-render rationale, warn=2, scoping, anchored grep,
negative tests) -> memory.md Quality gate.

### Config: tracked vs regenerated
tracked : rproject.toml rv.lock | pyproject.toml uv.lock .python-version | _targets.R R/*.R tests/*.R |
          _quarto.yml index.qmd _synthesis.qmd _qc.qmd _microglia.qmd _trajectory.qmd _mechanism.qmd _crossmodality.qmd theme.scss assets/fonts/*.woff2 | .Rprofile rv/scripts/*.R
          rv/.gitignore | scripts/install-*.sh | AGENTS.md .agents/skills/** .codex/prompts/*.md
regen   : rv/library _targets/ _report/ _freeze/ .quarto/ .venv tools/  (gitignored + read-economy skip);
          sccomp_draws_files/ (sccomp per-chain CSV draws at build CWD; gitignored)
