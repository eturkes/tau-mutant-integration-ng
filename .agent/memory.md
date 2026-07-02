# Memory - standing contract (read every session)

Durable facts / decisions / gotchas surviving all plans. Codex-only repo:
`AGENTS.md` is the canonical instruction file; prompt templates live under
`.codex/prompts/`; per-run Codex review output lives under ignored `.codex/runs/`.
Companions: `roadmap.md`
(direction), `map.md` (wiring), `history.md` (new decision digests),
`archive_digest.md` (v1 reference). This is a FRESH streamlined rebuild; v1 lives
on branch `archive`.

## Project
Integrate snRNAseq + GeoMx spatial + bulk proteomics + bulk phosphoproteomics across
4 mouse AD genotypes to read how amyloid (NLGF) reshapes microglia under different
tau backgrounds. Genotypes (2x2): MAPTKI (tau-KO baseline), P301S (mutant tau),
NLGF_MAPTKI (amyloid alone), NLGF_P301S (amyloid + tau). Design = tau (MAPTKI vs
P301S) x amyloid (-/+ NLGF) + batch. Divergence = interaction
(NLGF_P301S - P301S) - (NLGF_MAPTKI - MAPTKI).

## Raw data (storage/data = symlink -> host Documents tree, shared/external read-only copy;
gitignored via /storage/* + Codex read-economy skip; Rscript can read it when a task requires live data.
Shapes VERIFIED live in S2 via the R/io.R loaders; data is immutable so numbers are stable.)
- snrnaseq.rds 8.3G: Seurat; full = RNA 33683 genes (Assay5, ENSMUSG rownames) + SCT 28299
  (active assay). GOTCHA: dim(obj)=SCT(28299) but @misc$geneids (33683 symbols) aligns to RNA,
  NOT SCT -> build_symbol_map uses RNA rownames (asserts unique non-NA symbols; positional alignment
  rests on the v1 @misc$geneids contract, length-checked only). broad_annotations=="Microglia" -> 26104 cells.
  meta PRECOMPUTED (P1/QC consume, don't recompute): genotype (canonical), batch (batch01-04),
  genotype_batch (16-lvl factor, 4x4 fully crossed all-nonzero), sex, nCount/nFeature_RNA,
  percent_mt/ribo/malat1/contam, doublets. load_snrnaseq drops SCT+reductions -> 340MB RNA-counts
  microglia subset (the qs target microglia_seurat_raw).
- geomx.rds 22M: Seurat 19963 genes x 91 AOIs; genotype already canonical (n=20/20/28/23). Spatial
  design cols bio_rep/tech_rep/slide_rep/roi/segment/SampleID; NO batch/genotype_batch (differs from snRNAseq).
- proteomics_nonfiltered_nonnormalised.tsv 15M: Spectronaut PTM export, 45972 x 30. Cols 1-14 = PTM
  annotation (PG.*, Gene-pSite, PTM.SiteAA/Location, "Phopshosite probability"[sic typo]); 15-30 = 16
  intensity `Naoto-Hippo_TiO2_DIA_NN.raw.PTM.Quantity` (24M set only). peptide->protein-group sum = P4.
- phosphoproteomics_*.tsv 35.5M: Spectronaut PTM, 64328 x 81. 1-14 annotation; 15-81 = 67 intensity
  `*.PTM.Quantity` (Naoto 01-26 + Set6 01-41). The phospho target stores ALL 67 (nothing dropped at
  load); only the first 16 (Naoto 01-16) are the 24M timepoint -> the 24M column-subset is a P4 step.
- proteomics_sample_key.csv: 67 rows {`File name`, `Sample/Condtion`[sic]}; proteomics_sample_meta
  n_keep=16 (24M, 4 geno x 4 reps; asserts balanced reps + exact labels + unique join keys).
  normalise_ptm_stub (shared by the key producer + match_intensity_columns) strips `.PTM.Quantity`
  then trailing `.raw` -> both files' intensity cols collapse to one run stub, .raw discrepancy handled.
  match_intensity_columns is helper-only in P0 (NA = non-sample col); P4 wires it + asserts 16/16.
  read_spectronaut_tsv stop()s on any parse problem + strips readr spec/problems attrs (stale ->
  bad_weak_ptr after qs restore) -> stores a plain tibble.
Missing vs v1 (do NOT re-acquire unless a phase explicitly needs them): cisTarget
mm10 (SCENIC), SEA-AD h5ads (human validation) - both are v1 bloat, out of scope.

## Carried scientific decisions (v1-validated; treat as defaults, revisit per phase)
- Pseudobulk replicate = genotype_batch (16 ids, 4/genotype); batch in every design.
- Proteomics: peptide -> protein-group by sum; treat as total proteome.
- Phospho: first 16 samples = 24M timepoint; bulk hippocampus, NOT microglia-sorted
  (restate "not microglia-sorted" in any kinase prose).
- Microglia: re-cluster on the subset, drop only clear outliers (no over-pruning).
- Substates (v1 02a): homeostatic / DAM / IFN / proliferative.

## snRNAseq microglia reprocess (P1-S1, built) -- `R/microglia.R` -> target `microglia_processed`
- `reprocess_microglia`: SCT-v2(glmGamPoi, regress percent_mt+percent_contam) -> RunPCA(npcs=30) ->
  RunHarmony(batch ONLY) -> FindNeighbors+RunUMAP(harmony dims 1:20) -> FindClusters(Louvain algo 1, res
  {0.2,0.4,0.6}). Output (687MB qs, DefaultAssay SCT): reductions pca/harmony/umap, fresh SCT_snn_res.* +
  `microglia_clusters` (=res 0.4 primary, 12 clusters, Idents). Re-run STABLE: observed cluster
  ARI=1.0 under the recorded thread config -- NOT bitwise-guaranteed (clusters come from the seeded
  PCA->Harmony->SNN->Louvain chain; UMAP is viz-only, irrelevant; Annoy/BLAS threads unpinned per
  up-to-tolerance). seed 42 + RNGkind + thread snapshot in
  @misc$reprocess_provenance (read provenance$primary_cluster_col / $resolutions, NOT a grep -- robust to
  stale cols). Post-Harmony marker separation CONFIRMED (homeostatic/DAM/IFN/prolif distinct argmax clusters
  -> batch-only Harmony did NOT wash out substate biology). marker_mean_by_cluster(symbols mapped ->ensembl
  first; SCT/RNA rownames are ensembl) = the reusable separation helper (S2 reuses).
- `microglia_seurat_raw` carries STALE upstream meta from snrnaseq.rds processing: reduction COORDS as columns
  (pca1/2, umap1/2) + old clusters (SCT_snn_res.0.01, seurat_clusters) -> reprocess STRIPS the reduction-coord
  + non-computed-resolution shadows. PRESERVED for S2: `allen_labels` (v1 fine annotation -> reconcile), cell-
  cycle (S.Score/G2M.Score/Phase/cc_diff -> QC confound check), nCount/nFeature_SCT (recomputed by SCT anyway).
- Seurat-ecosystem GATE gotcha (forward to S2/S3): some fns emit a once-per-session WARNING (RunUMAP "default
  method changed") that lands in tar_meta(warnings) -> FAILS the zero-fault gate. Silence via the fn's OWN
  option (here Seurat.warn.umap.uwot=FALSE), never blanket suppressWarnings -- keeps every real warning's gate
  signal. SCTransform also trips future's 500MiB globals cap on 26k cells -> set options(future.globals.maxSize).
  Check each new pkg (UCell/sccomp) for the same patterns. harmony 2.0 dropped the v1 `assay.use` arg (assay
  implicit in reduction.use) -> v1 recipe args can be stale under the P3M-2026 pkg versions; verify signatures.
- P1-S1 codex-review hardening (durable): reprocess_microglia ends with BUILD-TIME postconditions (3 reductions
  + microglia_clusters factor + Idents match + no reduction shadows + only fresh *_snn_res.* + provenance) -> a
  silent recipe regression fails tar_make (the warn=2 unit tests skip the heavy body); the marker-SEPARATION
  argmax assertion is deferred to S2 QC (needs its curated ensembl sets). load_snrnaseq also clears
  @graphs/@neighbors (v1 source carries stale SCT_nn/SCT_snn). Matrix is used directly but INTENTIONALLY
  undeclared in rproject.toml (recommended pkg, ABI-coupled to R, always present via SeuratObject Imports ->
  an explicit CRAN decl only risks a mismatched reinstall).
- 5 canonical contrasts everywhere: tau_alone, nlgf_in_maptki, nlgf_in_p301s,
  tau_in_nlgf, interaction (factorial 2x2 + batch).
- State thresholds before applying them; present the axes with no pre-privileged winner.

## snRNAseq microglia annotate (P1-S2, built) -- `R/microglia.R::annotate_microglia` -> `microglia_annotated`
- UCell scoring: AddModuleScore_UCell(assay="SCT", slot="data", ncores=1) = DETERMINISTIC (rank-based, no RNG);
  `slot=`layer is GATE-SAFE on Assay5 (no deprecation warn -- verified 0 warnings). All markers map 100% to SCT
  ensembl rownames. SCT assay = 21333 genes (SCTransform sparsity-filtered from 33683 RNA); subset(cells=)
  preserves genes (only cells drop). Build 55s/612MB. Pure helpers unit-tested; heavy body smoke-tested live.
- Marker constants (constants.R) restructured: microglia_identity_markers (pan, state-INDEP QC: Csf1r/C1q/Ctss/
  Fcrls/Hexb/Tyrobp) -- NEVER use homeostatic markers for "is a microglia" (DAM downregulates P2ry12/Tmem119) +
  canonical_microglia_markers (Homeostatic / DAM[s1+s2 MERGED, one Apoe-Trem2 programme] / IFN / Proliferative +
  MHC_APC aux) + microglia_substate_levels (argmax set) + contam_signatures (Oligo/Neuron/Astro).
- PRUNE rule uses RAW identity-vs-best-contam, NOT z. z-based argmax FAILS: ambient oligo/neuron/astro is pervasive
  background -> z-centering destroys the absolute "is this a microglia at all". Drop if id_med<0.15 OR
  mglike_frac<0.30 (frac cells with raw identity > best raw contam) -- a 2-D OR rule (a cluster dropped on one axis
  can sit high on the other) catching BOTH low-identity AND neuron-doublet contaminants. id_floor is LOAD-BEARING
  (a kept cluster sat only 0.0076 over it) -> re-derive per dataset. Exact margins + cluster ids + DAM/lineage
  medians stored @misc$microglia_prune (separation, qc_rationale); dropped 2944/26104 (11.3%) -> 23160 retained.
- SUBSTATE: cluster-PRIMARY is authoritative; per-cell secondary is NOISY (z-argmax over-calls sparse IFN/Prolif)
  -> report cluster-level proportions, not per-cell. z-scale 4 substate sigs on RETAINED cells -> cluster-mean-z
  argmax; unassigned if best z<=0; ambiguous if top-two both>amb_floor(0.10) AND within tol(0.10) -- amb_floor
  guards a noise runner-up (cluster 2: weak-homeo 0.12 vs prolif-noise 0.07 -> Homeostatic, not ambiguous).
  Result Homeo{0,2,5}=11174, DAM{1,3,4}=11189, IFN{9,10}=797, Proliferative 0 clusters (no prolif-dominant
  cluster -- biologically fine; per-cell secondary over-calls it 4038, sparse-marker noise that cluster-mean
  averages out -> "0 Proliferative" = no coherent population, NOT 0 cells), 0 ambiguous/unassigned.
- HEADLINE confirmed descriptively (amyloid -> homeostatic->DAM): MAPTKI 3521H/415D -> NLGF_MAPTKI 2641H/4484D ->
  NLGF_P301S 1612H/5885D (S3 composition tests it; S4 DE the programme).
- CAVEAT for S3 (genotype-associated QC dropout): cluster 6 (worst QC, dropped) is 86% NLGF_MAPTKI -> prune
  removes more NLGF_MAPTKI low-quality cells (retained genotype frac 0.176/0.170/0.317/0.337 vs original
  0.167/0.168/0.337/0.329). Dropped clusters = low-identity {6,7,11} or neuron-doublets {8}; NONE DAM-high
  (DAM med 0.099-0.129 vs kept DAM clusters 0.195-0.225, stored) -> amyloid->DAM cannot be a prune artifact;
  REPORT the asymmetry, never hide it. Full stats in @misc$microglia_prune (qc_rationale, separation,
  dropped_by_genotype) +
  $substate_provenance (cluster_mean_z, substate_table, n_used, thresholds).

## snRNAseq microglia composition (P1-S3, built) -- `R/composition.R::test_composition` -> `composition_results`
- METHOD-STACK REVERSAL vs plan (reproducibility-driven): plan had sccomp PRIMARY + propeller sensitivity; BUILT
  the opposite. propeller (speckle) = LOCKED primary (logit) + sensitivity (asin), reproducible from the P3M
  snapshot. sccomp = OPTIONAL OFF-lock cross-check. WHY: CmdStan is a compiled C++ tree off the P3M snapshot ->
  rv cannot lock it -> a Bayesian arm cannot be THE reproducible call. The locked guarantee is propeller.
- propeller NEEDS a CELL-MEANS design (THE S3 build bug, fixed): speckle::propeller.ttest derives PropRatio by
  raising per-GENOTYPE mean-proportion coefficients to the contrast powers (apply(coef^cont, 2, prod)). A
  treatment/factorial design (Intercept+effects, factorial_design()) makes PropRatio a meaningless ratio of effect
  coefs AND HARD-CRASHES ("dim(X) must have a positive length") when a contrast loads ONE coefficient (tau_alone ->
  apply collapses to a vector). FIX: model.matrix(~ 0 + genotype + batch), genotype cols renamed to bare levels,
  contrasts via make_contrast_matrix() (cell-means, proven == factorial in test_design.R). t-stat is batch-ADJUSTED
  (full design); PropRatio is the MARGINAL genotype mean ratio (speckle refits raw props on ONLY the contrast
  genotype cols, batch dropped) -> for `interaction` a raw-scale ratio-of-ratios whose SIGN may differ from the
  logit t -> read direction from t. The t-stat alone would survive a factorial design; only PropRatio forces cell-means.
  run_propeller asserts the balanced fully-crossed design (table(geno,batch)==1) since PropRatio's marginal reading needs it.
- DISCORDANCE RULE (pre-declared): propeller-LOGIT stands as THE call; where asin or sccomp differ in effect SIGN
  or significance for a (contrast, substate), FLAG + report, never average. composition_concordance() keys
  (contrast,substate); dir=sign(t|c_effect), sig=FDR<alpha across available methods; sccomp absent -> propeller-
  only comparison (no flag from a missing arm). Each PRESENT method MUST cover every base (contrast,substate) key ->
  a missing/partial/empty arm fails LOUD (no false-green); run_sccomp asserts nonempty output.
- Batch random-vs-fixed asymmetry (intentional): propeller/limma FIXED batch (de_pb-consistent); sccomp RANDOM
  batch intercept (~ 0 + genotype + (1|batch)) -> priors regularise the 4-level batch. Unit = genotype_batch (16).
  sccomp cell-means -> colon-free contrasts (no backtick hazard).
- sccomp arm internals (P1-S3 build detail -- closed; full saga in git + archived plan): `cores`-ONLY fit_model call
  (the real parallelism knob; `parallel_chains` would collide via `...`); HMC health read from the final cmdstanr
  fit's $diagnostic_summary via attr(fit,"fit") (sccomp_test's c_rhat/c_ess are all-NA for contrast rows), hardened to
  NULL on API drift, RECORDED not gated. DURABLE (also in the Quality-gate section, forwarded to the P2-S3 glmmTMB
  arm): cmdstanr reports divergences via message() carrying a literal "Warning:" prefix, NOT warning() -> a
  warning-only handler misses it AND it reaches the tee'd log -> the gate's `^Warning:` scan reddens; the
  cached-target blind spot (check.sh rebuilds ONLY `report`) masked it until a forced fresh build. FIX: run_sccomp's
  withCallingHandlers now muffles+records BOTH warnings AND messages -> capture both for any heavy/optimiser arm.
- OFF-lock backend (scripts/install-cmdstan.sh, idempotent): cmdstanr from the Stan r-universe -> tools/rlib-stan
  (SEPARATE lib so `rv sync` never prunes it) + CmdStan compiled -> tools/cmdstan. _targets.R prepends the lib +
  sets CMDSTAN iff BOTH exist. sccomp_backend_ready() requires the PROJECT-LOCAL tools/rlib-stan + tools/cmdstan/
  cmdstan-* dirs (not just any CmdStan on .libPaths/CMDSTAN) -> a GLOBAL CmdStan can't activate the arm on a fresh
  clone -> propeller-only, fresh clone stays green. Orchestrator: backend present + sccomp STRUCTURAL error -> LOUD
  (allow_sccomp_failure=TRUE downgrades to skip); cmdstanr/CmdStan version + path RECORDED in provenance$sccomp_backend.
  sccomp dumps per-chain CSV draws to ./sccomp_draws_files/ at the build CWD -> gitignored (regenerable).
- STATUS: propeller primary + pure helpers UNIT-VALIDATED (tests/test_composition.R green at warn=2: count
  shapes/empty-drop/order/correctness, covariate-constancy + balanced-crossed + concordance-completeness fail-loud,
  propeller direction DAM-up + Homeostatic-down + prop_ratio>1 for logit AND asin, concordance flagging, sccomp gate
  logical). CODEX-REVIEWED + hardened (8 findings: concordance false-green, sccomp-error-now-loud, PropRatio
  interaction-scale doc, project-local backend gate, c_R_k_hat surfaced, balance guard, cores-comment fix, backend
  provenance). sccomp cores/chains API SOURCE-VERIFIED (cores-only call correct).
- LIVE-RUN VERIFIED 2026-06-30 (full scripts/check.sh green, incl. a FORCED fresh composition_results rebuild;
  cores fix holds, no parallel_chains collision; headline + run numbers -> history.md + provenance, drift-prone).
  sccomp final fit ~2-3% divergent (run-to-run), healthy E-BFMI/treedepth -> recorded NOT gated, treat SUPPORTIVE
  not definitive: the few-level (1|batch) funnel is the PLAUSIBLE divergence source but NOT localized by the summary
  diagnostics (do not claim it shown; adapt_delta the lever if a later phase hardens the Bayesian arm). Point
  estimates corroborate propeller; amyloid->DAM STRONG + concordant both methods, interaction DAM-positive sig
  under propeller / near-boundary under sccomp -> propeller-primary stands, static synergy -> P2. Concordance flags
  = sparse-IFN sign/sig noise (n=797) + the DAM-interaction cell only on sccomp draws >0.05 -> report INLINE-computes
  the live set (propeller-primary authoritative, never averaged).

## snRNAseq microglia pseudobulk DE (P1-S4, built) -- `R/de_pb.R` -> `pb_de_microglia` + `pb_de_substate`
- DE on RAW RNA counts (33683 ENSMUSG genes), NOT SCT -- microglia_annotated keeps BOTH assays; aggregate by
  genotype_batch (16 replicate units). SCT clusters/scores; pseudobulk does inference (counts). Whole-MG +
  per-substate, 5 canonical contrasts each. orchestrators: run_pb_de_microglia / run_pb_de_substate -> de_pseudobulk
  (generic: build_pseudobulk -> factorial_design -> fit_limma_voom -> stage_wise_test -> interaction_power).
- fit_limma_voom UPGRADED: voomWithQualityWeights (quality_weights=TRUE default) -- sample weights down-weight
  empirically-noisy units (correlates with, != low cell count); topTable confint=TRUE -> CI.L/CI.R (effect-size
  reporting); robust eBayes kept. quality_weights=FALSE -> plain voom (proteome/phospho may opt out later).
- pseudobulk_counts/build_pseudobulk gained `cells=` -> restrict to a cell subset (one substate) BEFORE aggregation,
  AVOIDS Seurat::subset (its benign messages would dirty heavy build logs). Caller pre-guards estimability.
- stageR (1.34.0, BioCsoft, on-snapshot -> LOCKED reproducible layer; pulled no heavy deps -- SummarizedExperiment
  already present). stage_wise_test: pScreen = fit$F.p.value (moderated omnibus F = "any genotype effect"; the 5
  contrasts span the rank-3 genotype subspace; limma's classifyTestsF keys df1 off the contrast-covariance RANK
  (df1=3, NOT ncol=5 -- VERIFIED: F.p.value == pf(F, 3, df.total); the moderated-F VALUE is mildly basis-dependent
  under rank deficiency, but its null df -- the screen's validity -- is right), pConfirmation = per-contrast p,
  method="holm" = stageR's MODIFIED post-screen Holm (NOT plain p.adjust Holm: folds in the OFDR screen scaling;
  FWER-valid under ARBITRARY dependence -> correct despite the rank-deficient contrast family), alpha=OFDR.
  getAdjustedPValues -> matrix [genes x {padjScreen, 5 contrasts}]; confirmation NA where the gene fails the screen
  -> a contrast is sig at OFDR iff its column <= alpha. GATE NOTE: getAdjustedPValues (NOT
  stageWiseAdjustment) emits ONE message() restating the fixed-OFDR caveat -- informational/deterministic, text has
  NO ^Warning/^WARN anchor (so it would not red the gate even on a fresh build reaching the log) -> suppressMessages
  it for clean logs. These targets stay cached under check.sh (only `report` is force-invalidated) so they don't
  re-run during the gate anyway.
- interaction_power: median per-gene posterior SE (sqrt(s2.post)*stdev.unscaled[,"interaction"]) -> NOMINAL minimum
  detectable log2FC at 80% power for the MEDIAN gene, via qt(median df.total). Per-test nominal-t power -- NOT
  stageR/OFDR/BH discovery power, NOT gene-specific; context for the threshold count, never a bare "0 genes".
- MIN-CELL FLOOR (per-substate fit-or-skip): fit iff EVERY genotype_batch unit has >= min_cells (default 10) of the
  substate (0-cell units fail too) -> full estimable factorial design; else SKIP -> descriptive-only. cell_counts
  (substate x 16-unit) table ALWAYS stored (report dropout asymmetry). Real argmax substates only. Both run_pb_de_*
  assert a COMPLETE genotype x batch crossing up front (assert_complete_crossing: n_units == prod covariate levels)
  -> an absent unit fails LOUD, never a silent <16-unit sub-design. LIVE: Homeostatic
  (min 52) + DAM (min 31) FIT; IFN (min 5, 15/16 units pass) + Proliferative (0) SKIP.
- LIVE RESULTS (2026-06-30, R 4.6 re-baseline NOT v1's margins; full counts -> history.md + the cached targets,
  drift-prone): whole-MG amyloid->DAM STRONG (DAM markers amyloid-UP, v1-concordant at the DE level). INTERACTION
  0 large-effect genes BUT 123 stageR-confirmed (|logFC|>0.5 median ~1.1, failing only the stricter standalone
  per-contrast FDR) + MDE@80%~0.92 log2FC -> UNDER-POWERED not absent (report the 123 + MDE/CI, synergy -> P2;
  no-LARGE-DE != no-DE, absence-of-evidence != evidence-of-absence). Per-substate: Homeostatic + DAM FIT (IFN/Prolif
  skip); fewer cells -> larger MDE (honest), amyloid still shifts WITHIN DAM.
- tests/test_de_pb.R extended (warn=2 clean): cells= subset + validation (all-present/no-dup, partial-bad fail-loud),
  SCREEN df1=rank=3 calibration (F.p.value==pf(F,3,df) != pf(F,5,df)), de_pseudobulk structure + CI cols + stageR
  matrix colnames + finite interaction MDE, run_pb_de_substate fit/skip statuses + 4x16 cell_counts + skip-reason +
  tighter-floor skips DAM, COMPLETE-CROSSING guard (absent unit -> fail-loud both orchestrators), dam_direction shape
  + zero-marker NA. Synthetic Seurat (RNG-free), voomWQW warning-free.

## snRNAseq microglia report (P1-S5, built) -- `_microglia.qmd` + `R/microglia.R::microglia_report_data` -> `microglia_report`
- CHEAP-RENDER INVARIANT (the S5 design constraint): the gate FORCE-renders the report every run reading CACHED
  targets -> a section qmd must NEVER tar_load the 612MB microglia_annotated. microglia_report_data extracts a
  COMPACT frame (per-cell {umap_1/2, genotype[factor], substate[factor], 3 *_UCell_z} cell_frame + n_cells + the
  small $prune/$provenance lists) -> `microglia_report` ~0.5MB. _microglia.qmd tar_loads microglia_report +
  composition_results + pb_de_microglia + pb_de_substate + symbol_map (all compact) -> render ~12s. The extractor
  ASSERTS finite z + non-NA genotype/substate (a downstream NA -> a ggplot "removed rows" warning -> warn=2 gate
  fail) -> the plotting data is render-clean by construction.
- _microglia.qmd setup `options(warn=2)` + tar_source(); included in index.qmd AFTER _qc.qmd (section flow in map.md).
  RENDER-WARNING avoidance (warn=2 INSIDE the qmd): ggrepel on the small labelled DAM subset uses max.overlaps=Inf
  (else "unlabeled data points"); volcano y = -log10(pmax(P.Value, 1e-300)) (caps zero-p -> Inf -> drop-warning);
  geom_text_repel seed=42L (deterministic layout).
- PROSE = INLINE-COMPUTED from the loaded targets, NEVER hardcoded -> tracks the actual cached build (sccomp
  FDR/divergences run-to-run variable -> a hardcoded value drifts). Three adversarial-review fixes (numbers must be
  RIGHT, not just present), traced live 2026-06-30:
  (1) interaction df: describe the RAW design df qualitatively ("few replicate degrees of freedom"), show ONLY the
      computed pbm$interaction$df (=24 = eBayes-MODERATED df.total) as the "effective total" -- df.total != residual
      df (~9); never label 24 "residual df".
  (2) sccomp interaction: a setup `int_sentence` reports the ACTUAL sccomp FDR + 95% CI bounds (or "not re-tested ...
      propeller-only" when sccomp NULL); never assert "CI spanning zero" (run-to-run).
  (3) dropout caveat: back "no dropped cluster DAM-high" against the DAM-ENRICHED ceiling kept_dam_max (~0.225), NOT
      the kept FLOOR -- dropped clusters' DAM medians (<=0.13) OVERLAP kept HOMEOSTATIC clusters (~0.10) so
      drop_dam_max<kept_dam_min is FALSE; report the dropout GENOTYPE (drop_geno=NLGF_MAPTKI via which.max rowSums)
      vs CLUSTER (drop_clu=6 via colSums) each correctly labelled (old prose mislabelled the genotype "cluster").

## snRNAseq microglia activation trajectory (P2-S1, built) -- `R/trajectory.R` -> target `microglia_trajectory`
- build_activation_trajectory: slingshot 2.20.0 (returns a PseudotimeOrdering; version-stable accessors
  slingPseudotime / slingLineages / slingCurves) on harmony[1:15] of microglia_annotated (NO recompute -- reads the
  cached harmony embedding + per-cell UCell scores + substate labels). PRIMARY = a FORCED single Homeostatic->DAM
  lineage: clusterLabels = the 2 substate super-clusters (MST on 2 nodes = one edge = one lineage, clean BY
  CONSTRUCTION; the principal curve still fits the full 15-D cell cloud) -> NO arbitrary root-cluster pick.
  start.clus="Homeostatic" orients pt (homeostatic low). IFN (797 cells, 3.4%) + Proliferative (0) OMITTED from the
  lineage (on_lineage flag + NA pt), never deleted -> per-unit omitted fraction reported.
- pseudotime = an ACTIVATION ORDERING (position/extent of advance), NOT developmental time/potency -- direction rests
  on slingshot rooting + DAM-marker monotonicity (validated post-hoc), never a potency claim. FORK CHOICES (dims=15,
  H->D-only, 2-cluster labels) PRE-DECLARED, not retuned after contrast inspection.
- LIVE (2026-06-30, R4.6): single clean lineage; mean pt Homeo 22.3 < DAM 36.2 (direction correct);
  Spearman(pt,DAM_UCell)=+0.56, (pt,Homeo_UCell)=-0.40 (rooting validated, recorded). Score-axis concordance rho=0.62
  (n=22363) -- MODERATE-LARGE positive, clears the 0.5 gross-failure floor; EXPECTED below 1 (slingshot = transcriptome
  geometry vs score-axis = DAM-minus-Homeo marker contrast: related, NOT identical -> the score-axis is a CONCORDANCE
  anchor sharing the marker system, NOT statistically-independent robustness). `concordant` flag RECORDED not gated.
  SENSITIVITY highly robust: dims-10 rho 0.990 / dims-20 0.997 / all-retained 0.990 vs primary. all-retained (IFN
  included, dims 15) = a single PATH Homeostatic->IFN->DAM (IFN intermediate, not a branch).
- SELECTION-EFFECT audit (the lineage-conditioning risk): omitted (IFN) fraction BALANCED across genotypes
  (low single-digit %, exact in per_unit table) -> conditioning on H->D barely skews interaction; report, never hide. (16 units)
  stores n_cells / n_on_lineage / omitted_frac.
- score_axis_pt = RAW DAM_UCell - Homeostatic_UCell (assumption-light: both UCell [0,1] rank scores, no population
  z-centring). pt01 = Smithson-Verkuilen squeeze (min-max scale -> (y*(n-1)+0.5)/n) into the OPEN (0,1) for the S3
  beta GLMM; off-lineage NA preserved. Target = COMPACT list {cell_frame 23160x10, per_unit, lineage,
  sensitivity df, provenance(versions/seed/RNG/threads + dims/rho/dam_pt_rho/omitted)}; serialized ~0.8MB / in-memory
  ~3.3MB -- NOT the 612MB Seurat (build
  reads it once, ~218s for 4 slingshot fits; one-time, the gate force-renders cached targets). 0 build warnings.
- slingshot GOTCHA: its covariance-scaled MST distance (.dist_clusters_scaled -> solve(s1+s2)) ERRORS "system is
  computationally singular" on a near-degenerate per-cluster covariance -> synthetic TEST fixtures need FULL-RANK
  well-conditioned blobs (distinct-frequency sinusoids, comparable amplitude), NOT tiny collinear ripples (real
  harmony clusters are full-rank, no issue). run_slingshot_lineage extracts the DAM-TERMINAL lineage (>2 clusters may
  branch -> pick the terminal-in-DAM lineage, longest if tied; off-lineage cells keep slingshot's NA).
- rproject.toml: slingshot (BioCsoft; pulled CRAN princurve [needs Rcpp] + Bioc TrajectoryUtils, SingleCellExperiment already
  present). glmmTMB (per-cell GLMM sensitivity) deferred to P2-S3 where it is actually loaded (ABI-warning handling
  belongs with the live load test). Pure helpers UNIT-tested (tests/test_trajectory.R, warn=2 clean) on
  make_trajectory_embedding (helpers.R: 2-cluster + with_ifn branch) + validate_trajectory_units metadata guards;
  orchestrator covered by the gate tar_make (full-Seurat unit test DEFERRED). run_slingshot_lineage is RNG-PURE
  (pins seed + all 3 kinds, restores caller .Random.seed on exit); rooting now GATED both ways (dam_pt_rho>0,
  homeo_pt_rho<0). validate_trajectory_units = fail-loud per-unit audit (no-NA/blank unit, geno in levels, one geno/unit).

## snRNAseq microglia trajectory estimation core (P2-S2a, built) -- `R/trajectory.R` (6 pure fns, NO target)
- Collapse on-lineage per-cell pt -> 16 genotype_batch summaries -> factorial_design + 5 contrasts, ordinary t at
  9 resid df (NO eBayes). Fns: derive_batch, pseudotime_per_replicate, ordinary_t_table, fit_trajectory_contrasts,
  kitagawa_channels, decompose_progression_vs_composition + within_state_col(state)=paste0("within_",tolower(state))
  = the ONE col sanitizer S2a writes + S2b reads back (never hand-case). EXACT 3-channel Kitagawa shift-share
  (comp/prog/cross) holds on RAW pt_raw ONLY (logit/asin break additivity); the interaction contrast is
  intercept-free so it annihilates the unit-constant -> L(mean_pt)=L(comp)+L(prog)+L(cross), reconstruction <1e-8
  VERIFIED. Coefficient-reconstruction exactness needs ONE shared per-unit weight vector replicated across the 4
  channel-rows (same WLS operator each row); differing per-row weights generally break it. Parametric weighted
  contrast SE: limma::contrasts.fit derives a multi-coef contrast SE from the UNWEIGHTED corr ((X'X)^-1) -> exact
  only for a single-coef contrast OR a balanced design (interaction=tau_nlgf is single-coef -> always exact; the
  real 16-unit factorial is balanced -> exact, but a dropped unit unbalances it). fit_trajectory_contrasts now
  OVERRIDES stdev.unscaled with the exact per-feature normal-equation value (chol2inv(chol(X'W_iX))) when weighted
  -> ordinary t exact for EVERY contrast independent of balance. The 16-unit/9-df complete crossing is enforced by
  S2b's run_trajectory_progression (assert_complete_crossing); the S2a primitives stay crossing-agnostic (full-rank
  guard only).
- GATE gotcha (limma 1-row fit): limma::lmFit DROPS the single feature rowname (coefficients rowname NULL) ->
  fit_trajectory_contrasts RESTORES rownames(coefficients/stdev.unscaled) <- rownames(measure_mat) so
  ordinary_t_table keys the measure labels (forward to ANY 1-feature limma fit).
- GATE gotcha (tapply 1-D array): tapply over a SINGLE factor returns a 1-D ARRAY (carries `dim`) -> array*matrix
  is "non-conformable arrays" (both have dims) whereas a plain vector recycles -> coerce single-factor tapply
  results to a plain named vector (as.numeric+setNames) before broadcasting against a matrix (the mu_bar pooled-mean
  bug). rowSums/colSums already return plain vectors (no issue).
- FIXTURE make_trajectory_cell_frame (helpers.R; S2b/S3 reuse): 4x4 geno x batch, midpoint ramp ((i-0.5)/n*0.3)
  -> each block mean EXACTLY 0.15 for ANY n -> the pure-composition fixture is EXACTLY pure (unequal DAM count does
  not shift a within-state mean). DEFAULT adv -> pure within-state interaction 0.4 (prog loading 1); FLAT adv +
  dam_extra>0 -> pure composition (comp loading 1); jitter>0 = NON-additive ((gi*bi)%%5)*jitter -> sigma>0 for S2b's
  structural orchestrator test. NO batch col -> exercises derive_batch.

## snRNAseq microglia trajectory progression inference (P2-S2b, built) -- `R/trajectory.R` +2 fns -> target `trajectory_progression`
- run_trajectory_progression(microglia_trajectory): COMPACT S1 target -> pseudotime_per_replicate -> factorial_design (9 resid
  df) -> 3 limma fits {weighted direct (mean_pt/median_pt/q90/within_<used>, SHARED n/sd^2 weights = mean-precision heuristic, NOT quantile-exact for median_pt/q90 + chosen for EXACT Kitagawa reconstruction so progression_cf = shared-weight DECOMPOSITION test, not channel-specific IV inference), ols sensitivity, bounded
  (frac_past logit+asin VST, w=n_cells, EXPLORATORY)} + decompose_progression_vs_composition + freedman_lane_interaction on
  {progression_cf, within_homeostatic, frac_past_logit, mean_pt}. PRE-REGISTERED primary BH family = {progression_cf,
  within_homeostatic}; rest separate exploratory BH; mean_pt FLAGGED composition-conflated. Pure-R, NO new dep.
- freedman_lane_interaction = WLS-as-OLS permutation null (permute REDUCED-model weighted residuals; pivot-free XtXinv once;
  RNG-pure -- pins seed+3 kinds, restores caller stream on.exit, RNGkind-THEN-.Random.seed). SENSITIVITY not nominal-exact
  (weights ESTIMATED from same summaries -> approximate exchangeability); each FL call weighted to MATCH its limma fit.
- HEADLINE (R4.6 re-baseline; EXACT numbers DRIFT-PRONE -> S4 inline-computes, NEVER hardcode; the QUALITATIVE finding is durable):
  the 2x2 interaction RAISES mean pseudotime position (mean_pt coef>0, ~p 0.04) BUT the 3-channel Kitagawa shows this is
  COMPOSITION (comp_cf loading ~+1.25, the SIG channel ~fdr 0.025) NOT progression (prog_cf loading NEGATIVE ~-0.55, NS, perm_p
  ~0.18; cross ~+0.29; loadings sum 1.0, recon 6.7e-15). PRE-REGISTERED primary {progression_cf, within_homeostatic} BOTH NS ->
  no statistically SUPPORTED progression-beyond-composition signal (negative/NS estimate, NOT proven absence -- absence of
  evidence). DIVERGES from v1 (~0.94 progression loading): R4.6 attributes the interaction advance MAINLY to MORE DAM CELLS
  (composition -- CONFIRMS P1's sig DAM-fraction interaction), with no DETECTED cells-advancing-FURTHER signal. mean_pt
  alone (+2, p~0.04) MISLEADS as "advance" -> the decomposition is LOAD-BEARING (read comp/prog channels, never bare mean_pt). S4
  frames composition-confirmed + progression-null honestly.
- tests: source R/de_pb.R (assert_complete_crossing dep); FL null (design-orthogonal resid -> t_obs~0 -> perm_p~1) / signal
  (2*tau_nlgf col -> perm_p<0.05) / determinism / RNG-purity / weighted; structural orchestrator on the jitter>0 fixture
  (sigma>0 dodges the limma zero-variance warn under warn=2 -- the EXACT-pure fixture saturates -> sigma=0 -> Inf t).
- REVIEW-HARDENED (codex, 8 findings accepted): LOCAL 16-unit/9-df full-4x4 crossing guard in the orchestrator --
  assert_complete_crossing (shared P1 fn) only checks n_units==prod(OBSERVED levels), so a dropped batch (4x3=12u/6df)
  PASSES it; guard LOCALLY (leave the shared fn untouched -> no heavy-pb-target invalidation). min_within floor 1->2
  (within-state sd needs >=2 cells/unit) + explicit within-sd finite+positive guard; provenance += planned_primary /
  primary_within_skipped (primary = analyzable subset if root state too thin, not silently 1-row). FL lm-oracle TEST
  GOTCHA: a design-orthogonal y_det -> near-collinear WEIGHTED fit -> |t|~5e15 -> all.equal(tolerance=) REQUIRED (abs()<eps
  fails though relative match is 4e-16 when |t| is huge). Weight HONESTY: w_overall=n/sd^2 = shared mean-precision (chosen
  for exact Kitagawa reconstruction), NOT channel-specific IV nor quantile-precision -> progression_cf = shared-weight
  DECOMPOSITION test; "no progression synergy" softened to "no statistically SUPPORTED signal" (absence of evidence).

## snRNAseq microglia trajectory glmmTMB sensitivity (P2-S3, built) -- `R/trajectory.R::glmmtmb_pt_sensitivity` -> target `trajectory_glmm_sensitivity`
- SUPPORTIVE per-cell beta GLMM (replication-aware confirmation modelling the full bounded distribution the 16-unit
  summary collapses): pt01 ~ tau*amyloid + batch + (1|unit), glmmTMB::beta_family(); FIXED batch (de_pb-consistent),
  (1|unit) random intercept. On-lineage = finite pt01; tau/amyloid integer 0/1 from genotype (matching factorial_design);
  batch=derive_batch, unit=genotype_batch factors. Reads the COMPACT microglia_trajectory$cell_frame (NOT 612MB); target
  INDEPENDENT of trajectory_progression. Returns list(method, term, estimate, se, z, p_value, ci_l, ci_r, re_sd,
  singular, n_cells, n_units, fail_reason, warnings, messages); n_units = genotype_batch clusters present (asymptotics
  basis, RECORDED not asserted -> a dropped unit is honestly reported, not silently reframed as 16).
- DEGRADE cascade (graceful, NEVER blocks the limma-summary primary): health gate = .fit_health_ok() PURE helper
  (pdHess & convergence==0 & finite est & finite se>0 & finite z & valid p in [0,1] & non-singular re_sd>=1e-4) -> on any
  fail (incl. nonestimable interaction = est NA from a rank-deficient column drop, or NULL-on-error) -> rank-normal LMM
  rn=qnorm((rank(pt01)-0.5)/n) ~ same formula, gaussian() (SAME package on-lock, SAME gate) -> if BOTH fail, method="failed"
  + NA effect + fail_reason (per-arm error/nonestimable/singular/nonconverge) + captured warnings/messages. A FIT failure
  NEVER throws; MALFORMED INPUT (missing cols / boundary pt01 / unknown genotype not in genotype_levels / broken
  genotype_batch) fails LOUD via stopifnot (surfaces an upstream break, not masks it as failed). method in {glmmTMB_beta,
  lmm_ranknorm, failed}.
- Wald row read by POSITION (cond cols fixed order Estimate/Std.Error/z/Pr(>|z|)) -> a column-NAME drift can't silently
  mis-extract; a positional-integrity guard (z==est/se, p==2*pnorm(-|z|), tol 1e-5) catches the converse (a column-ORDER
  change) by degrading; term=intersect(c("tau:amyloid","amyloid:tau"), rownames(cond)) asserted length-1.
  .capture_quietly = withCallingHandlers muffling+recording BOTH warnings AND messages (the sccomp lesson: TMB optimisers
  report health via message() too -> a fresh build would red warn=2/tar_meta/^Warning: scan) -> 0 leaked to the gate.
- glmmTMB on CRAN -> P3M serves the trixie BINARY (the plan's "source-compile TMB" framing was off; CRAN=binary, ABI-built
  together at the snapshot -> loads 0-warning under warn=2; no co-pin handwork needed). glmmTMB 1.1.14 / TMB 1.9.21 /
  Matrix 1.7.5. Namespace-qualified (glmmTMB::), NOT in tar_option_set packages (like slingshot); rproject.toml + rv.lock co-pin.
- LIVE (2026-07-01, R4.6; EXACT numbers DRIFT-PRONE -> S4 inline-computes, the QUALITATIVE read is durable):
  method=glmmTMB_beta on 22363 on-lineage cells, tau:amyloid est +0.123 (logit), se 0.050, z 2.44, p~0.015, CI
  [0.024,0.222] excludes 0, re_sd 0.047 (non-singular). Build 1.5s, tar_meta NA/NA. COMPOSITION-CONFLATED (the per-cell
  MEAN-position analogue of mean_pt; corroborates the position-shift interaction + P1's propeller DAM-fraction composition
  signal) -> NOT progression-specific; the S2b Kitagawa decomposition stays LOAD-BEARING for composition-vs-progression.
  Supportive only; concordance AND discordance both fine.
- tests (tests/test_trajectory.R, warn=2, deterministic nlminb / NO RNG; invisible(loadNamespace("glmmTMB")) at head = load
  WITHOUT attaching, so an accidental unqualified prod call still fails): .fit_health_ok branches unit-tested directly (each
  FALSE arm, no optimiser coaxing); jitter=0.3 fixture -> non-singular RE -> beta fit (assert method=="glmmTMB_beta",
  n_units=16, finite effect+CI, fail_reason NA); DEFAULT fixture (identical within-genotype units) -> singular both arms ->
  failed (fail_reason "singular", 0 msgs); no-amyloid subset -> tau:amyloid rank-deficient, glmmTMB DROPS it (captured
  MESSAGE not exception) -> failed (fail_reason "nonestimable", msg "rank-deficient", n_units=8); unknown genotype -> fail loud.

## snRNAseq microglia trajectory report (P2-S4, built) -- `_trajectory.qmd` + `R/trajectory.R::trajectory_report_data` -> `trajectory_report`
- CHEAP-RENDER INVARIANT (mirrors P1-S5): _trajectory.qmd tar_loads ONE compact target (`trajectory_report` ~0.34MB),
  NEVER the 612MB Seurat -> gate force-render stays cheap. trajectory_report_data bundles the 3 COMPACT trajectory
  targets (microglia_trajectory / trajectory_progression / trajectory_glmm_sensitivity) into: slim per-cell cell_frame
  (genotype/substate/on_lineage/pt_raw/score_axis_pt) + interaction table (primary+exploratory BH families;
  coef/CI/perm_p/FDR) + weighted_top (5 contrasts) + per_unit + lineage_per_unit + sensitivity + glmm 13-name subset +
  provenance. Drift-prone inference numbers (coefs/p/loadings/rho, R version, unit count) INLINE-COMPUTED from
  trajectory_report (R4.6 re-baseline), NEVER hardcoded; only fixed design constants (9 resid df, sensitivity dims
  10/20) stated as text (they change only with the design/config, not the data).
- EXTRACTOR GUARD-BAR (S4a, mirrors microglia_report_data; the qmd render-cleanliness contract): up-front input-schema
  stopifnot on EVERY nested field the body reads + assembled-bundle postconditions with COL-EXISTENCE-BEFORE-FINITENESS
  (a dropped col fails HERE, never vacuously -- all(is.finite(NULL))==TRUE lets a missing col slip). Guards: finite
  coef/CI/p/fdr across BOTH interaction families; perm_p finite on the 2 INLINED rows {mean_pt, progression_cf}; each of
  the 5 canonical weighted_top contrasts has a finite mean_pt coef/CI (feeds p_ctr geom_pointrange); glmm 13-name set +
  provenance scalars finite/int-valued/string by role. The finite-p/fdr assertion is an INTENTIONAL build-fatal
  data-quality gate (validated non-degenerate on CURRENT data, NOT a proven universal); a future degenerate EXPLORATORY
  endpoint would red it -> then add graceful "-" formatting + exempt that row.
- FIXTURE COMPOSITION-DEGENERACY (baked into the S4a guard test, durable): CONSTANT per-unit composition (12/12) zeroes
  the comp_cf/cross residual variance -> se=0 -> NaN p -> NaN fdr -> the finite-fdr postcondition reds. The test VARIES
  per-unit composition (drop k%%3 DAM cells/unit, each >=10 retained) so comp_cf/cross are non-degenerate.
- decomposition NOT re-bundled (S4b, codex 955): the qmd draws its loadings figure+prose from provenance
  (composition/progression/cross_loading, guarded finite) + per-channel coefs from the comp_cf/progression_cf/cross ROWS
  of the interaction table -> a `decomposition` field only DUPLICATED those two live sources (dead figure-shaped output).
  Kept: `glmm` (full 13-name mirror of the S3 target) + `per_unit` (full 16-unit summary) as faithful RESULT records
  despite a few unread fields -- distinct from decomposition (a redundant duplicate of data used elsewhere).
- RENDER GOTCHA (ggplot2 4.0.3, generalises to any qmd): scale_*_gradient(`trans=`) DEPRECATED at 3.5.0 -> use
  `transform=` (else a lifecycle deprecation WARNING -> warn=2 render error -> red gate). The concordance geom_bin2d fill
  uses transform="log10". Every other trajectory geom already renders clean under _microglia.qmd's 4.0.3 pass.
- WIRING: index.qmd includes _trajectory.qmd AFTER _microglia.qmd; the Overview names it. `@sec-trajectory` cross-refs
  resolve across the full rendered doc regardless of include order. _microglia.qmd's 2 forward-pointers rewritten (P2
  BUILT; synergy = DAM composition, NO supported further-advance; cross-ref @sec-trajectory). tar_quarto auto-detects the
  trajectory_report tar_load edge THROUGH the include; NO new extra_files (theme/fonts shared), NO new dep (patchwork
  already used by _microglia.qmd/_qc.qmd).
- CHAPTER (title "The tau-amyloid synergy adds DAM cells rather than advancing them"): amyloid-axis-shift ->
  composition-not-progression 3-channel decomposition -> per-cell glmmTMB supportive (composition-conflated) ->
  reconciliation/robustness/concordance -> 5 caveats+provenance (activation-ordering not time; position not rate;
  lineage-conditioning balanced; cell-weighted anchors; transcriptionally-close substates). Full render = 52 chunks
  0-warning, report ~3.6MB, 23 targets.

## Environment (project-local; NO Docker, NO system-wide installs)
- Run as eturkes:eturkes (single-user Distrobox) -> files land user-owned, NO chown
  needed (v1's `chown rstudio:rstudio` was a rocker artefact, obsolete).
- R is 4.6.0 here (v1 pinned 4.5.2 / Bioc 3.22) -> numbers WILL differ; we
  re-baseline, never reproduce v1's locked margins (18/12/55).
- Stack (P0 built): **rv** (R pkgs) + **uv** (Python) + project-local **Quarto** + **targets** DAG,
  P3M-pinned (snapshot 2026-06-22, CRAN+Bioc same date). No bitwise guarantee (no Docker):
  targets=pipeline, rv.lock/uv.lock=versions, P3M=pinning. Fresh-clone bootstrap ORDER:
  `scripts/install-sysdeps.sh` -> `install-rv.sh` + `install-quarto.sh` -> `rv sync` -> `uv sync` -> `tar_make()`.
- **rv MUST be on PATH** (~/.local/bin, like uv): `.Rprofile`->`rv/scripts/activate.R` finds rv via
  `Sys.which("rv")` + shells `rv info` to set `.libPaths(rv/library)`; a tools/-only rv breaks
  activation. Pinned (version+sha256) in `install-rv.sh`. `.Rprofile` runs activate.R then a
  fail-loud guard: non-interactive `stop()` unless `rv/library` is in `.libPaths()` (catches
  rv-off-PATH / `rv info` fail / R-version-mismatch safe-mode -> NO silent global-lib fallback;
  re-add the guard if rv regenerates `.Rprofile`). NO repos override (base-R `install.packages`
  would write off-lock; use `rv add` OR edit rproject.toml + `rv sync`; Bioc pkgs need the
  repository qualifier e.g. `{ name = "glmGamPoi", repository = "BioCsoft" }`). rv installs to
  `rv/library/<R-ver>/<arch>/<distro>` (NOT `rv/library/<pkg>`) -> check a pkg via Rscript
  `requireNamespace`, not `test -d rv/library/<pkg>`.
- Repos (`rproject.toml` -> `rv sync` -> `rv.lock`): CRAN = plain `p3m.dev/cran/<date>` (rv inserts
  `__linux__/trixie` -> binary; do NOT hardcode the binary path); Bioc 3.23 =
  `p3m.dev/bioconductor/<date>/packages/3.23/{bioc,data/annotation,data/experiment,workflows}` +
  `force_source` (source-only on Debian). `tar_source` lives in `targets`; the `quarto` R pkg finds the
  pinned CLI via `QUARTO_PATH` (`_targets.R`; a `file.exists` preflight `stop()`s if missing -> no
  silent PATH-quarto fallback); `_quarto.yml`: render `*.qmd`+`!rv/` (rv#332), `freeze:false`
  (targets owns caching -> no stale `_freeze` divergence).
- Sysdeps (`scripts/install-sysdeps.sh`; `rv sysdeps` returns [] on trixie -> useless): build-essential
  + gfortran (Bioc source compiles) + libglpk40 (libglpk.so.40 for the igraph binary). Re-derive any new
  missing lib: `ldd`-scan `rv/library/**/*.so` for "not found".
- Python: uv `.venv` (gitignored), `pyproject.toml`+`uv.lock`+`.python-version` 3.13 (empty deps in P0); SOTA per phase.
- Heavy installs/compute: expect long runs; smoke-test helpers via `Rscript -e '...'`
  against the live data BEFORE any full run.

## Quality gate (P0-S5) -- `scripts/check.sh` (self-documenting header), fail-loud zero-fault, run from anywhere
tar_make's exit is NOT enough: it returns 0 on CAPTURED target warnings + is blind to the report's SEPARATE
knitr/Quarto render. Layered: env (`rv sync`+`uv sync`; skip `CHECK_SKIP_SYNC=1`) | tests warn=2 | FORCE-render
report tee'd to a log | `tar_meta` error+warnings all-NA (scoped to current targets+branches) | anchored render-log
grep. CHEAP (~12s: reads cached ~0.3GB targets, does NOT re-run the heavy load_snrnaseq build). Durable lessons:
- EVERY section qmd setup carries `options(warn=2)` -> a chunk warning halts the render (else it renders SILENTLY
  into the HTML, never reaching the log/tar_meta). `tar_quarto(report, quiet=FALSE)` lets Quarto/Pandoc `[WARNING]`
  lines reach the log (default quiet=TRUE hid them).
- CACHED-TARGET BLIND SPOT (load-bearing): check.sh force-invalidates ONLY `report`; every other target stays
  CACHED during the gate -> a heavy/off-lock arm's warnings never re-surface unless forced. Any "RECORDED not
  gated" claim MUST be re-verified on a FRESH build of that target (tar_invalidate it + report), not a cached one
  (the sccomp message() hole hid exactly here -- see P1-S3; same risk for the P2-S3 glmmTMB arm).
- render-log scan = `command grep -nE` (GNU, not the rg-fff shadow), ANCHORED forms only (`^[WARNING]` `^Warning:`
  `^Warning in ` `^Warning messages?:` `^WARN`) so benign "warn" substrings can't false-red; exit 0=fault/1=clean/
  >=2=infra, sits in an `if`-cond (exempt from set -e + ERR trap). NEGATIVE-TESTED both ways (chunk warning -> exit
  1; benign lines ignored). Residual: out-of-band edits to `_report/` output aren't scanned (moot -- re-rendered each run).
- Committed tests `tests/test_*.R` = plain `stopifnot` fail-loud scripts (zero new deps), source the R/ files they
  exercise + `tests/helpers.R` (deterministic synthetic fixtures, NO RNG/clock; expect_error[+pattern]), print
  `ok - <x>`, non-zero exit on fail. Set: test_design, test_de_pb, test_io, test_plot, test_composition,
  test_microglia, test_trajectory. Data-free; per-module live-data smoke-test still happens once before commit.
- Reproducible: fresh clone -> bootstrap order (map.md) -> `scripts/check.sh` green end-to-end.

## Operational
- Prose register: British English; human-facing report/figure text uses hyphens over
  em/en-dashes (commas or parentheses for asides, colons for restatements). R `#`
  comments + code stay exempt (LLM-facing).
- Read economy set lives in `AGENTS.md`: avoid generated/heavy trees (`storage/**`,
  `_targets/`, `_report/`, `.quarto/`, `rv/library/`, `tools/`, `.venv/`, fonts,
  completed plans) unless the task requires them. For generated outputs, prefer compact
  targets, target metadata, or runtime indirection over large direct reads. storage/data
  is a symlink now -> `git check-ignore` on a path UNDER it fatals (`beyond a symbolic
  link`); probe `storage/data` itself.
- New `R/*.R` = pure functions loaded by targets via `tar_source()`; the DAG orders
  execution (no manual dependency loader). Heavy producers = `tar_target`s storing
  `format="qs"`/`"file"`; Python steps = `uv run <script>` as `tar_file` targets.
- Test/DE gotchas (S3; harness itself described under Quality gate): edgeR 4.x ->
  `edgeR::normLibSizes` (calcNormFactors deprecated, emits a message); `limma::makeContrasts`
  needs ALL design columns to be syntactically valid R names (cell-means form: rename genotype
  columns to bare levels, build batch from a named factor; an inline `factor()` yields invalid
  names). Sort character keys with `method = "radix"` (locale-independent -> reproducible
  column/level order across machines). factorial_design fails loud if add_batch=TRUE but the
  batch column is absent -> a batch-less modality (GeoMx) MUST pass add_batch=FALSE. de_pb fitters
  fail loud on misaligned inputs (limma fits/contrasts by POSITION, only warns on name mismatch):
  both assert identical(colnames(data), rownames(design)) + identical(rownames(contrasts),
  colnames(design)) + full-rank design (qr rank == ncol); factorial_design guarantees full-rank
  output; pseudobulk_counts asserts the meta.data<->counts cell alignment Seurat maintains.
- targets serialization (verified S2): `format="qs"` works via the **qs2** backend (the
  `qs` pkg is NOT installed; qs2 serializes Seurat Assay5 objects fine). `format="file_fast"`
  is deprecated -> use `format="file"` + `tar_option_set(trust_timestamps=TRUE)` (set in
  `_targets.R`) so the 8G snrnaseq input is change-detected by mtime/size, not re-hashed each
  run. Heavy seurat target: `memory="transient"` + `garbage_collection=TRUE` release the load.
- Codex-only project config: keep `CLAUDE.md`, `.claude/`, and `.serena/` untracked/ignored.
  Session/review reusable prompts belong in `.codex/prompts/`; run outputs belong in
  `.codex/runs/` (ignored).

## Reports (Quarto; built P0-S4)
- Report = ONE self-contained OFFLINE HTML: a standalone `format: html` doc (`index.qmd`,
  `embed-resources: true` inlines CSS/JS/fonts) + modular `{{< include _section.qmd >}}` (leading
  `_` -> not rendered standalone). NOT a Quarto book (multi-file + `Could not fetch resource
  ./<sibling>.html` nav warnings under embed-resources -> caught by the gate's render-log scan once quiet=FALSE). tar_quarto
  still detects `tar_load`s inside EACH included `_*.qmd` (verified: edges through both `_qc.qmd` and `_microglia.qmd`
  -- the latter pulls microglia_report/composition_results/pb_de_microglia/pb_de_substate/symbol_map; the gate's
  force-render + tar_meta scope exercise them); list the theme/css in `tar_quarto(extra_files=)` (inspection misses them).
- Theme + fonts SHIP via theme.scss (Quarto 1.9.38, verified live on the production render 2026-06-29):
  scss:defaults COLOUR vars ($primary/$link-color #B0344D, $code-color #3F5A6B) + the IBM Plex stack
  ($font-family-sans-serif/$headings-font-family/$font-family-monospace) + 9 `@font-face` (scss:rules)
  with a relative `url("assets/fonts/<n>.woff2") format("woff2")`. Quarto base64-INLINES each woff2 into
  the embedded CSS under embed-resources -> ONE offline file (render PROVED: 9 faces inlined, magic d09GMg,
  0 external). The 9 woff2 are COMMITTED (assets/fonts/; avoid direct reads via `AGENTS.md` read economy);
  list them in `tar_quarto(extra_files=)` -- inspection misses them, `list.files("assets/fonts",
  pattern="woff2", full.names=TRUE)` keeps the list in sync. ggplot panels keep `theme_tau(base_family="")`
  (device font) so figures stay decoupled from this chrome + warning-free.
- Theme-CSS DETECTION gotcha: with an `@font-face url()` in the theme, Quarto URL-encodes the WHOLE compiled
  theme CSS into a `data:text/css,...` URI -> to verify embedding, match ENCODED tokens (`#B0344D`->`%23B0344D`,
  woff2 magic -> `d09GMg`); RAW `.count` reads ~0 theme-side and URLdecoding the ~1MB blob is very slow. Figures
  embed as `data:image/png` base64, so their colours are not raw either.
- Quarto caches the Sass compile in `.quarto/` -> a theme edit is invisible until cleared; `.quarto`
  is generated/heavy -> clear via runtime indirection (R `unlink(".quarto", recursive=TRUE)` in the render
  script), not a direct large-tree read. Inspect output HTML via a tiny python/R script building the path
  ("_"+"report") + `.count` substrings (a `#hex` regex proved flaky).

## Codex workflow
- Fresh session: start from `.codex/prompts/session.md` or manually follow its load order.
- Review: `scripts/codex-review.sh [focus]` runs read-only generic `codex exec` with
  `.codex/prompts/reviewer.md` because Codex CLI 0.142.5 rejects custom stdin prompts
  with `review --uncommitted`; accept/reject findings explicitly, fix accepted ones, then
  commit one scoped unit.
- Headroom: `.agent/context.sh` scans newest Codex JSONL session for this cwd and parses
  latest `token_count` (`last_token_usage.input_tokens`, `model_context_window`); override
  the window with `CODEX_CONTEXT_WINDOW`.
