# P1 plan - snRNAseq microglia core

Active-plan working doc. Posture BUILD: streamlined SOTA-informed rebuild, ONE cohesive
story, KISS. Companions: roadmap.md (direction), memory.md (contract), map.md (wiring),
archive_digest.md (v1 ref). Closed -> fold to history.md, archive here. Gate decisions
locked 2026-06-29 (see "Resolved decisions"); fold the durable ones to memory.md at close.

## Scope
From `microglia_seurat_raw` (RNA counts, 33683 x 26104, already broad-annotated Microglia)
produce, as targets DAG nodes feeding ONE incremental report section:
1. reprocessed (SCT) + integrated (Harmony) + clustered microglia object;
2. substate labels (homeostatic / DAM / IFN / proliferative; aux MHC/APC score);
3. compositional test of substate shifts across the 5 contrasts (sccomp + propeller);
4. pseudobulk DE across the 5 contrasts (whole-MG + per-substate) -> the robust
   amyloid->DAM activation programme (headline #1) + an honest STATIC-NULL interaction
   that forward-points to P2 (trajectory);
5. a `_microglia.qmd` report section + CLOSE-OUT.

Out of P1 (deferred / dropped): trajectory + interaction-as-rate -> P2; TF/kinase
mechanism -> P3; cross-modality -> P4. SINGLE-CELL DE dropped entirely (gate). v1 bloat
dropped: SCENIC, spatial-decon, gene-level dynamics, celltype-specificity, CARNIVAL,
hdWGCNA, the ledger/contest machinery.

## Resolved decisions (gate 2026-06-29)
1. Single-cell DE = DROPPED. Pseudobulk limma-voom is the SOLE DE inference (Squair 2021 /
   Murphy-Skene 2022 / Crowell 2020: cell-level DE = pseudoreplication -> FDR inflation).
   v1's NEBULA static interaction was NULL anyway -> nothing lost. SUPERSEDES de_pb.R's
   "de_sc.R" forward-ref comment (fix that comment in S4). No nebula dep; no S5-single-cell step.
2. Composition = sccomp PRIMARY + propeller-logit cross-check (concordance = robustness at n~=4).
   Accept the cmdstanr/CmdStan install weight.
3. Normalisation = SCTransform v2 (v1 continuity). Restores v1's reprocess recipe (SCT +
   glmGamPoi, regress percent_mt+percent_contam); substate scoring on the SCT assay. v1 Harmonised
   over c("batch","sex") (archive rmd/02a:21); REVISED to batch-only (sex perfectly aliased w/ batch -
   batch01/03 male, 02/04 female, archive_digest - so batch-only is equivalent-or-finer, sex is absorbed).
Folded into default (SOTA-clear, no gate): UCell over AddModuleScore (per-signature calibrated before
argmax - raw UCell scores aren't cross-signature comparable); Harmony over BATCH ONLY (never integrate
over genotype/amyloid -> DAM is biology); 4 substates by CALIBRATED argmax + ambiguous/unassigned bucket
+ aux MHC/APC score; voomWithQualityWeights + eBayes(robust=TRUE); per-contrast BH base FDR + stageR
screen-confirm ACROSS the 5-contrast family (a contrast-FAMILY tool; degenerate on one contrast alone).

## SOTA research carried (2026; mandated by AGENTS.md) -> memory.md at close
- Pseudobulk = sole condition/population DE inference; cell-level DE FDR-inflated (decision 1).
- UCell (rank-based, robust to dropout/depth/batch) > AddModuleScore (control-bin, unstable on sparse nuclei).
  Scores NOT cross-signature calibrated (depend on signature size/coherence) -> z-scale per signature before
  argmax; assign at CLUSTER level (primary) + per-cell argmax (secondary); ambiguous (top-two within tol /
  all sub-null) -> unassigned bucket, never force-assign.
- Harmony over batch only; never over genotype/amyloid; verify DAM/IFN/prolif still separate post-Harmony
  (else lower theta / go uncorrected). scVI = overkill at 26k/4-batch.
- Composition: sccomp (Beta-binom, sum-constrained, factorial+interaction+RE, outlier-robust; Mangiola 2023)
  + propeller-logit (Phipson 2022, best n=3-5) cross-check; aggregate to 16 sample units, never nuclei.
- Interaction in pseudobulk: diff-of-diffs ~2x SE, ~9 resid df -> eBayes cross-gene shrinkage STABILISES the
  variance (does NOT recover replication/power -> the interaction stays under-powered; report effect-size/CI,
  not just p), robust=TRUE, NOT treat(); voomWithQualityWeights down-weights EMPIRICALLY-noisy units
  (correlates with, not equal to, low cell count); batch FIXED (4 crossed levels too few for
  duplicateCorrelation); per-contrast BH = base FDR, stageR adds a family-level screen-confirm across the 5
  contrasts (degenerate on one contrast).
- CARRY caveat (Thrupp 2020, CellRep): snRNA depletes ~18% of DAM activation genes (Apoe/Spp1/Cd74/B2m/Cst3)
  -> attenuated DAM expected; SCORE don't threshold; DE on RAW counts regardless. CellBender NOT applicable
  (needs raw CellRanger matrix w/ empty droplets; we have the 26k filtered subset) -> regress
  percent_mt+percent_contam in SCT + flag contaminant clusters instead.
- QC confounds (precomputed in meta, memory: doublets, percent_ribo, percent_malat1): audit by
  cluster/genotype/substate (v1 trajectory flagged ribo/ambient confounding); prune/down-weight w/ logged
  rationale. Locked SCT recipe regresses percent_mt+percent_contam; adding ribo/malat1 to vars.to.regress is
  a FALLBACK iff a cluster is confounded (don't silently mutate the locked recipe).
- Canonical mouse markers (expand constants.R per these): homeostatic P2ry12/P2ry13/Cx3cr1/Tmem119/Hexb/
  Sall1/Selplg/Siglech/Olfml3/Gpr34; DAM-s1(Trem2-indep) Tyrobp/Apoe/B2m/Ctsb/Ctsd/Fth1/Lyz2; DAM-s2(Trem2-
  dep) Trem2/Cst7/Lpl/Cd9/Itgax/Clec7a/Spp1/Gpnmb/Igf1/Axl/Cd63; IFN Ifit1/2/3/Irf7/Oasl2/Isg15/Mx1/Ifitm3/
  Usp18/Bst2/Rsad2/Stat1; proliferative Mki67/Top2a/Birc5/Mcm5/Stmn1/Cenpa; MHC/APC(aux) Cd74/H2-Aa/H2-Ab1/
  H2-Eb1/H2-K1/Tap1/B2m. DAM~MGnD~ARM~WAM = ONE Apoe-Trem2 convergent programme (don't treat as orthogonal);
  ARM=DAM+MHC (Sala Frigerio 2019, App-NL-G-F -> most NLGF-relevant).

## Science framing (the spine for P1)
- HEADLINE (robust): amyloid (NLGF) drives homeostatic->DAM activation; microglia-led. P1 evidence =
  pseudobulk nlgf_in_maptki / nlgf_in_p301s DE (DAM_up genes top-ranked, direction matches v1) +
  composition (DAM proportion up with amyloid) + UCell DAM score shift. Antigen-presentation thread =
  aux MHC/APC score (part of the activation axis).
- INTERACTION (honest, outcome-OPEN): run the `interaction` contrast. v1 PRIOR = ~0 genes (matched-power
  null) but acceptance is PROCESS not outcome -> report the result WITH a power/effect-size statement (CI or
  MDE on the diff-of-diffs); a "static-null" call must be BACKED by that, never asserted from "0 genes passed
  FDR" (absence of evidence != evidence of absence). If genes DO surface, report them, don't dismiss as noise.
  Forward-point: the tau x amyloid synergy is a progression-RATE effect -> tested in P2 (trajectory). No
  pre-privileged axis; thresholds stated before applying.

## Method stack (locked)
- Reprocess (SCT-v2; v1 recipe, Harmony REVISED batch-only): SCTransform(method="glmGamPoi", vst.flavor="v2",
  vars.to.regress=c("percent_mt","percent_contam")) -> RunPCA(npcs=30) -> RunHarmony(group.by="batch")
  -> FindNeighbors+RunUMAP(reduction="harmony", dims=1:20) -> FindClusters(Louvain, multi-res, pick).
  assay="SCT". Seed+RNGkind+thread-count fixed & recorded -> reproducible UP TO TOLERANCE, NOT bitwise (memory
  contract; multithreaded Harmony/UMAP/RcppParallel). NO genotype/amyloid in integration.
- Substates: UCell::AddModuleScore_UCell(assay="SCT") on {Homeostatic, DAM, IFN, Proliferative} -> z-scale
  per signature -> CLUSTER-level assignment primary, per-cell argmax secondary (raw UCell not cross-signature
  calibrated); ambiguous (top-two within tol / all sub-null) -> unassigned bucket; cross-tab vs de-novo
  clusters + v1 labels (reconcile); aux MHC/APC + *_contam + rbc + doublet/ribo scores -> flag/prune
  contaminant clusters (logged rationale, not silent).
- Composition: propeller (speckle) PRIMARY (logit) + asin SENSITIVITY -- LOCKED/reproducible from the P3M
  snapshot; cell-means ~0+genotype+batch (speckle PropRatio needs per-genotype mean coefs), batch FIXED (de_pb-
  consistent). sccomp(~0+genotype+(1|batch), RANDOM batch) = OPTIONAL OFF-lock cross-check gated on the CmdStan
  backend (unlockable compiled C++ -> can't be the reproducible primary; REVERSED from the original sccomp-primary
  plan). 16 sample units (genotype_batch); 5 contrasts incl. interaction. Discordance rule (pre-declared):
  propeller-logit call STANDS; asin/sccomp sign-or-significance differences flagged+reported, never averaged. Batch
  random(sccomp)-vs-fixed(propeller) asymmetry intentional (priors regularise few-level batch); state it.
- Pseudobulk DE: build_pseudobulk(replicate=genotype_batch) on RNA counts -> extend fit_limma_voom
  (voomWithQualityWeights + eBayes(robust=TRUE)) across 5 contrasts; whole-MG always; per-substate ONLY where
  every genotype_batch unit clears a PRE-DECLARED min-cell floor (e.g. >=10/unit) for that substate, else skip
  -> descriptive-only + log the cell-count table. per-contrast BH = base FDR + stageR family-screen across the
  5 contrasts. filterByExpr(group=genotype) per level.

## Steps (each closeable in one window; acceptance per step; [HEAVY]=compute, [DEP]=new pkgs+gate)
Each step: write pure fn + synthetic unit test (gate-independent) -> smoke-test on live data ->
wire target + full run -> verify quality gate (scripts/check.sh) before AND after.

- **S1 reprocess + cluster** [HEAVY][DEP: harmony, glmGamPoi]. **[DONE 2026-06-29]**
  R/microglia.R::reprocess_microglia() (SCT recipe above) + marker_mean_by_cluster (separation check); target
  `microglia_processed` (qs: SCT+pca+harmony+clusters+umap; 12 clusters @res0.4, 687MB). harmony+glmGamPoi
  added. ACCEPT met: reductions {pca,harmony,umap}+cluster factor present; post-Harmony marker separation
  confirmed (distinct argmax homeostatic/DAM/IFN/prolif); re-run STABLE (observed ARI=1.0, recorded threads,
  not bitwise-guaranteed; seed 42 +
  RNGkind + thread snapshot in @misc$reprocess_provenance); gate green. Pkg-drift fixes vs v1: harmony 2.0 drops
  assay.use; future.globals.maxSize raised for SCT; Seurat.warn.umap.uwot=FALSE (else UMAP notice -> tar_meta
  warning -> gate fail). Stale upstream meta shadows (pca1/umap1, SCT_snn_res.0.01) stripped.

- **S2 substate annotation + QC prune** [DEP: UCell]. **[DONE 2026-06-30]**
  constants.R restructured (microglia_identity_markers pan-QC + canonical_microglia_markers
  Homeostatic/DAM/IFN/Proliferative+MHC_APC + microglia_substate_levels + contam_signatures);
  microglia.R::annotate_microglia (UCell score -> prune -> calibrated argmax) + pure helpers
  marker_sets_to_ensembl/zscale_signatures/assign_substate/cluster_mean_z/flag_contaminant_clusters; target
  `microglia_annotated` (612MB). UCell added (BioCsoft). ACCEPT met: every retained cell labelled-or-bucketed
  (postcondition !anyNA); per-signature z-calibration; per-substate self-enrichment asserted (build-time);
  contaminant clusters {6,7,8,11}=2944 cells dropped w/ @misc$microglia_prune rationale (id_med<0.15 OR
  mglike_frac<0.30, thresholds in natural gaps; doublets all-0 no-op); substate x genotype table in
  $substate_provenance; Thrupp noted (constants+code); gate green. Key finding: amyloid->homeostatic->DAM
  confirmed descriptively; genotype-associated QC dropout (cluster 6 = 86% NLGF_MAPTKI) carried to S3 caveat.
  Gotcha: z-based prune FAILS (ambient contam pervasive -> use RAW identity-vs-contam); per-cell substate noisy
  (cluster-level primary authoritative).

- **S3 composition** [DEP: speckle + sccomp(off-lock cmdstanr/CmdStan)]. **[DONE 2026-06-30]**
  R/composition.R::test_composition (propeller primary + gated sccomp) across 5 contrasts; target
  `composition_results`; tests/test_composition.R (synthetic count table, KNOWN amyloid->DAM direction). speckle +
  sccomp via rv; cmdstanr/CmdStan OFF-lock via scripts/install-cmdstan.sh (project-local tools/, gitignored).
  LIVE-VERIFIED: tar_make built composition_results (29.5s) on real microglia_annotated + full scripts/check.sh
  green; sccomp ran (final fit ~2.4% divergent, E-BFMI 0.72, recorded). Fixed latent diagnostics bug (codex
  c_R_k_hat column absent in sccomp 2.4.0 + per-contrast c_rhat structurally NA) -> capture fit$diagnostic_summary.
  ACCEPT (revised): propeller(logit) PRIMARY + asin sensitivity, sccomp OPTIONAL gated cross-check (reproducibility
  reversal: CmdStan unlockable); per-contrast proportion estimates incl. interaction; amyloid->DAM quantified+
  tested; discordance per the pre-declared rule (propeller-logit stands); batch random-vs-fixed asymmetry justified;
  replicate=genotype_batch; gate green. ALL MET (live-verified 2026-06-30).

- **S4 pseudobulk DE** [DEP: stageR(Bioc)].
  Extend R/de_pb.R (voomWithQualityWeights + eBayes(robust=TRUE); stageR helper; FIX the stale "de_sc.R"
  forward-ref comment -> single-cell DE dropped); targets `pb_de_microglia` (whole-MG) + `pb_de_substate`
  (per-substate) x 5 contrasts. rv add stageR. Extend tests/test_de_pb.R.
  ACCEPT: contrasts computed + topTables stored per contrast x level; pre-stated thresholds applied
  (FDR<0.05 base, |logFC|>0.5 for "sig"); amyloid->DAM direction-concordance with v1 REPORTED (whatever it
  shows, not required to "surface"); interaction reported WITH a power/effect-size statement (CI/MDE on the
  diff-of-diffs) - null claimed only if backed by that; per-substate fit-or-skip by the min-cell rule +
  cell-count table stored; gate green.

- **S5 report section + CLOSE-OUT** [render].
  `_microglia.qmd` (UCell substate UMAP, proportions, composition results, amyloid->DAM DE programme,
  static-null interaction + P2 forward-pointer, Thrupp caveat); include in index.qmd; add to
  tar_quarto(extra_files) if needed; figures via theme_tau + scale_*_genotype (British English, hyphens).
  New tests wired into gate. Update memory (SOTA decisions, gotchas) + map (P1 wiring) + history (P1 digest);
  archive this plan; reset roadmap Active plan.
  ACCEPT: report renders 0-warning with microglia section; full gate green end-to-end; docs updated.

## New dependencies (rv add per step)
harmony + glmGamPoi (S1) | UCell (S2) | sccomp + speckle (S3; sccomp pulls cmdstanr -> CmdStan compile) |
stageR (S4). Re-derive any new sysdep via ldd-scan of new .so (memory).

## Risks / watch
- Harmony over-correction washing out DAM -> verify marker separation; reduce theta or go uncorrected.
- sccomp/CmdStan install weight (Stan compile) -> if it fights the project-local pin, fall back to
  propeller-only + record the blocker (decision was sccomp+propeller; degrade only if truly stuck).
- SCT DE-conservative on nuclei (SOTA caveat the user accepted for v1 continuity) -> DE is on RAW counts
  (pseudobulk) so this hits clustering/scoring not DE; watch that SCT clustering still resolves substates.
- Per-substate pseudobulk thin for rare states (IFN/proliferative) -> PRE-DECLARED min-cell floor per unit
  gates fit-or-skip (not just filterByExpr); report n; substate DE is secondary to whole-MG.
- 8G load peak is the load_snrnaseq BUILD (cached); reprocess works on the 340MB subset -> moderate.
- Keep gate green each step; near 80% context -> drive to clean state + close out.
