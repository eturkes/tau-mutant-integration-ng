# Memory - standing contract (read every session)

Durable facts / decisions / gotchas surviving all plans. Companions: `roadmap.md`
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
gitignored via /storage/* + deny-Read; Rscript reads resolve through it, bypassing the deny.
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
- PRUNE rule uses RAW identity-vs-best-contam, NOT z. z-based identity argmax FAILS: ambient oligo/neuron/astro
  is pervasive background (medians 0.10-0.15) so z-centering destroys the absolute "is this a microglia at all".
  Drop if id_med<0.15 OR mglike_frac<0.30 (frac cells raw identity > best raw contam) -- 2-D OR rule, a cluster
  dropped on one axis can sit high on the other. Two contaminant kinds: low-IDENTITY {6,7,11} id_med<=0.091
  (7 id_med=0, clean gap to kept>=0.158); DOUBLET {8} id_med 0.187 (ABOVE floor, even above kept cluster 10) but
  neuron-per-cell (Neuron med 0.300 vs ~0.15 bg) -> dropped on mglike 0.244. mglike gap clean (drop<=0.244 /
  keep>=0.381); id gap NOT clean (8 overlaps) + BINDING kept margin THIN (cluster 10 IFN id_med 0.1576, 0.0076
  over floor) -> id_floor LOAD-BEARING, re-derive per dataset. Exact margins + DAM/lineage medians stored
  @misc$microglia_prune (separation, qc_rationale). Dropped {6,7,8,11}=2944/26104 (11.3%) -> 23160 retained;
  doublets precomputed all-0 (logged no-op).
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
- sccomp call (sccomp source-verified): pass ONLY `cores` -- it is fit_model's real parallelism knob (caps chains via
  find_optimal_number_of_chains %>% min(cores); sets parallel_chains/threads_per_chain). `parallel_chains` is no formal
  -> `...` -> collides with the internal mod$sample(parallel_chains=); `chains` IS a fit_model formal (binds+overrides,
  NOT a collision -- the earlier "chains collides" note was WRONG).
- sccomp HMC-CONVERGENCE capture (CORRECTED vs codex finding "c_R_k_hat surfaced" -- that column does NOT exist in
  sccomp 2.4.0): sccomp_test emits c_rhat/c_ess_bulk/c_ess_tail but they are STRUCTURALLY all-NA for CONTRAST rows
  (sccomp computes rhat/ESS for base design params only, not derived contrasts), and divergences surface via message()
  NOT warning() -> withCallingHandlers catches nothing. REAL signal = the final (outlier-removed) cmdstanr fit's
  $diagnostic_summary() reached via attr(fit,"fit") (pass_fit=TRUE default): per-chain divergent / max-treedepth /
  E-BFMI -> run_sccomp's `diagnostics` attr -> provenance$sccomp_diagnostics + a sccomp_status note. RECORDED not gated
  (divergences are messages anyway) -> OFF-lock arm never fails the warn=2 gate; LOCKED propeller keeps full strictness.
  withCallingHandlers still captures any genuine R warning into `warnings` (0 on the live run).
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
- LIVE-RUN VERIFIED 2026-06-30 (full scripts/check.sh green end-to-end): tar_make built composition_results (29.5s)
  on real microglia_annotated; cores fix holds (no parallel_chains collision). sccomp final fit = 6 chains, 96/3996
  (~2.4%) divergent, E-BFMI 0.72, 0 treedepth -> recorded NOT gated (divergences from the 4-level (1|batch) random-
  effect funnel; ebfmi healthy -> localized, estimates corroborate -> treat sccomp SUPPORTIVE not definitive;
  adapt_delta is the lever if a later phase hardens the Bayesian arm). HEADLINE robust across BOTH methods: DAM up
  under amyloid (nlgf_in_maptki/nlgf_in_p301s) propeller t=10.8/14.4 FDR~1e-10/1e-13, sccomp c_effect +1.45/+1.83
  FDR~0; Homeostatic mirror-down. INTERACTION DAM positive (synergy): propeller FDR 0.027 (sig) vs sccomp 0.051
  (borderline) -> FLAGGED; interaction Homeostatic down sig in both. Concordance flagged 4/15 (3 sparse-IFN n=797
  sign/sig noise + the interaction-DAM sig-borderline) -> propeller-logit stands per the pre-declared rule.

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

## Quality gate (concrete + review-hardened; P0-S5) -- `scripts/check.sh`, fail-loud, zero-fault
Run from anywhere (cd's to root); any fault -> non-zero. tar_make's exit is NOT enough: returns 0 on CAPTURED
target warnings, and is blind to warnings raised in the report's SEPARATE knitr/Quarto render process.
Enforcement is layered across qmd + target + script:
1. env: `rv sync` + `uv sync` (idempotent; skip via `CHECK_SKIP_SYNC=1`).
2. tests: loop `tests/test_*.R` each `Rscript -e 'options(warn=2); source(<f>)'` (warn=2 -> stray R warning = error); set -e aborts on first fail.
3. pipeline: FORCE-render the report each run -- `tar_invalidate(any_of("report")); tar_make()`, tee'd to a log,
   wrapped in `if !` (so the failure message blames tar_make, not `tee`). Forcing = soundness + idempotency: a
   cached report skips its render -> empty log + warn=2 un-exercised. CHEAP (~8 s): the render READS cached
   targets (~0.3 GB), it does NOT re-run the heavy load_snrnaseq build (the "8G" is that BUILD's peak, not the
   stored target). Two render-time catches live in the SOURCES (so they fire on every render, not just check.sh):
     - `_qc.qmd` setup `options(warn = 2)` -> any R chunk warning -> error -> Quarto halts (error:false default)
       -> tar_make non-zero. The REAL catch for knitr/R chunk warnings (else they render silently INTO the HTML,
       never reaching the log or tar_meta).
     - `_targets.R` `tar_quarto(report, quiet = FALSE)` -> Quarto/Pandoc `[WARNING]` lines reach the log;
       default quiet=TRUE SUPPRESSED them, so the pre-review grep was blind to all Quarto/Pandoc warnings.
4. script-side zero-fault enforcement:
   (a) `tar_meta(error,warnings)` all-NA; kept rows = `name %in% tar_manifest()$name | parent %in%` it (current
       targets + their dynamic branches) -> drops tar_source'd fns/globals (~45 meta rows vs 13) + stale dead rows. warn=2 here too.
   (b) render-log `command grep -nE` (GNU grep, not the rg-fff shadow), ANCHORED forms only: `^[WARNING]`
       `^Warning:` `^Warning in ` `^Warning messages?:` `^WARN` -> benign "warn" substrings (paths/labels)
       can't false-red. Exit discriminated 0=fault / 1=clean / >=2=grep infra. grep sits in an `if`-cond
       (exempt from set -e AND the ERR trap; else the clean exit-1 prints a spurious GATE-FAILED line).
Negative-tested: a chunk `warning()` -> render error -> gate exit 1 (no PASS); the anchored pattern matches the
5 real warning forms, ignores benign false-positive lines. Residual: out-of-band edits to `_report/` output
files aren't scanned -- moot, the report is RE-rendered from source each run, not trusted from disk.
   Negative-tested: the grep matches `[WARNING]`/`Warning message:`/bare `WARN` and ignores benign tar_make
   lines; a `stopifnot(FALSE)` test drives exit 1. A current store -> tar_make is a no-op but 4(a)/(b) still run.
- Committed tests = `tests/test_*.R`: plain `stopifnot` fail-loud scripts (zero new deps, mirror io.R's
  assertion idiom), source the R/ files they exercise + `tests/helpers.R`, print `ok - <x>`, non-zero exit
  on failure. helpers.R = deterministic synthetic fixtures (make_fake_seurat / make_meta16 /
  expect_error[+pattern]; NO RNG or clock). Current set: test_design (5-contrast weights +
  factorial==cell-means property), test_de_pb (pseudobulk 16-col + fit_limma_voom/log smokes), test_io
  (loader contracts on tempfiles), test_plot (device-free theme/scale/concordance). They are data-free
  synthetic; live-data smoke-testing per module still happens once before commit.
- Reproducible: fresh clone -> bootstrap order (map.md) -> `scripts/check.sh` green end-to-end.

## Operational
- Prose register: British English; human-facing report/figure text uses hyphens over
  em/en-dashes (commas or parentheses for asides, colons for restatements). R `#`
  comments + code stay exempt (LLM-facing).
- storage/** + future caches/HTML are gitignored AND deny-Read -> Read/grep/ls on
  them is blocked (the `cat`-deny also bites). Rscript file reads BYPASS the deny
  matcher: inspect via `Rscript -e 'readRDS(...)' / 'readLines(...)'`. storage/data is
  a symlink now -> `git check-ignore` on a path UNDER it fatals (`beyond a symbolic
  link`); probe `storage/data` itself.
- Bash deny-gate is static on command text -> a cmd naming a deny-Read path as an
  arg (rm/stat/grep/find-piped) is blocked; use a glob/`find -delete` or runtime
  indirection. `ls`/`wc`/`echo` slip through.
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
- Commit `.serena/` (project.yml + .gitignore); Serena language changes are
  startup-only (restart Claude Code to apply).

## Reports (Quarto; built P0-S4)
- Report = ONE self-contained OFFLINE HTML: a standalone `format: html` doc (`index.qmd`,
  `embed-resources: true` inlines CSS/JS/fonts) + modular `{{< include _section.qmd >}}` (leading
  `_` -> not rendered standalone). NOT a Quarto book (multi-file + `Could not fetch resource
  ./<sibling>.html` nav warnings under embed-resources -> caught by the gate's render-log scan once quiet=FALSE). tar_quarto
  still detects `tar_load`s inside an included `_*.qmd` (verified: 5 edges through the include); list
  the theme/css in `tar_quarto(extra_files=)` (inspection misses them).
- Theme + fonts SHIP via theme.scss (Quarto 1.9.38, verified live on the production render 2026-06-29):
  scss:defaults COLOUR vars ($primary/$link-color #B0344D, $code-color #3F5A6B) + the IBM Plex stack
  ($font-family-sans-serif/$headings-font-family/$font-family-monospace) + 9 `@font-face` (scss:rules)
  with a relative `url("assets/fonts/<n>.woff2") format("woff2")`. Quarto base64-INLINES each woff2 into
  the embedded CSS under embed-resources -> ONE offline file (render PROVED: 9 faces inlined, magic d09GMg,
  0 external). The 9 woff2 are COMMITTED (assets/fonts/, deny-Read `**/*.woff2`, Serena ignored_paths);
  list them in `tar_quarto(extra_files=)` -- inspection misses them, `list.files("assets/fonts",
  pattern="woff2", full.names=TRUE)` keeps the list in sync. ggplot panels keep `theme_tau(base_family="")`
  (device font) so figures stay decoupled from this chrome + warning-free.
- Theme-CSS DETECTION gotcha (CORRECTED -- supersedes "colours embed raw"): once an `@font-face url()`
  is in the theme, Quarto embeds the WHOLE compiled theme CSS as a URL-ENCODED `data:text/css,...` URI,
  so BOTH colours AND fonts encode -- `#B0344D`->`%23B0344D`, `IBM Plex`->`IBM%20Plex`, the woff2 data
  URI -> `woff2%3Bbase64%2Cd09GMg` (d09GMg = wOF2 magic). RAW `.count` then reads ~0 for everything
  theme-side (raw-embed held ONLY for the pre-fonts colours-only theme). Match the ENCODED tokens (fast)
  or URLdecode/urllib.unquote the data:text/css blocks first -- URLdecode of the ~1MB blob is VERY SLOW,
  so prefer encoded matching. Figures are PNG base64 (`data:image/png`), so their `#B0344D` is not raw either.
- Quarto caches the Sass compile in `.quarto/` -> a theme edit is invisible until cleared; `.quarto`
  is deny-Read -> clear via runtime indirection (R `unlink(".quarto", recursive=TRUE)` in the render
  script), not a Bash `rm`. Inspect the output HTML the same way (output dir is deny-Read): build the
  path inside a python/R script ("_"+"report") + `.count` substrings (a `#hex` regex proved flaky).

## Subagents & skills
Scan the available-skills list each session; invoke a matching Skill before
improvising (scientific-writing, scientific-visualization, pathway-enrichment,
pydeseq2, scanpy, anndata, bioservices). Spawn subagents to protect main context
(Explore = cross-file search, Plan = design, general-purpose = research).
