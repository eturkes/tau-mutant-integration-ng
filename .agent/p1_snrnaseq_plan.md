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
   glmGamPoi, regress percent_mt+percent_contam); substate scoring on the SCT assay.
Folded into default (SOTA-clear, no gate): UCell over AddModuleScore; Harmony over BATCH ONLY
(sex perfectly aliased w/ batch -> batch absorbs it; never integrate over genotype/amyloid ->
DAM is biology); 4 argmax substates + aux MHC/APC score; voomWithQualityWeights + eBayes(robust=TRUE)
+ stageR for the interaction.

## SOTA research carried (2026; mandated by AGENTS.md) -> memory.md at close
- Pseudobulk = sole condition/population DE inference; cell-level DE FDR-inflated (decision 1).
- UCell (rank-based, robust to dropout/depth/batch) > AddModuleScore (control-bin, unstable on sparse nuclei).
- Harmony over batch only; never over genotype/amyloid; verify DAM/IFN/prolif still separate post-Harmony
  (else lower theta / go uncorrected). scVI = overkill at 26k/4-batch.
- Composition: sccomp (Beta-binom, sum-constrained, factorial+interaction+RE, outlier-robust; Mangiola 2023)
  + propeller-logit (Phipson 2022, best n=3-5) cross-check; aggregate to 16 sample units, never nuclei.
- Interaction in pseudobulk: diff-of-diffs ~2x SE, ~9 resid df -> eBayes cross-gene shrinkage is the rescue,
  robust=TRUE, NOT treat(); voomWithQualityWeights down-weights low-cell units; batch FIXED (4 crossed
  levels too few for duplicateCorrelation); stageR screen-confirm FDR (+ per-contrast BH).
- CARRY caveat (Thrupp 2020, CellRep): snRNA depletes ~18% of DAM activation genes (Apoe/Spp1/Cd74/B2m/Cst3)
  -> attenuated DAM expected; SCORE don't threshold; DE on RAW counts regardless. CellBender NOT applicable
  (needs raw CellRanger matrix w/ empty droplets; we have the 26k filtered subset) -> regress
  percent_mt+percent_contam in SCT + flag contaminant clusters instead.
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
- INTERACTION (static-null, honest): run the `interaction` contrast; EXPECT ~0 genes (v1 matched-power
  null). State it plainly + forward-point: the tau x amyloid synergy is a progression-RATE effect ->
  tested in P2 (trajectory). Present axes with no pre-privileged winner; thresholds stated before applying.

## Method stack (locked)
- Reprocess (SCT, v1 recipe): SCTransform(method="glmGamPoi", vst.flavor="v2",
  vars.to.regress=c("percent_mt","percent_contam")) -> RunPCA(npcs=30) -> RunHarmony(group.by="batch")
  -> FindNeighbors+RunUMAP(reduction="harmony", dims=1:20) -> FindClusters(Louvain, multi-res, pick).
  assay="SCT". Seed-fixed (deterministic). NO genotype/amyloid in integration.
- Substates: UCell::AddModuleScore_UCell(assay="SCT") on {Homeostatic, DAM, IFN, Proliferative} ->
  per-cell argmax = substate; cross-tab vs de-novo clusters (reconcile); aux MHC/APC + *_contam + rbc
  scores -> flag/prune contaminant clusters (logged rationale, not silent).
- Composition: sccomp(~tau*amyloid + (1|batch)) PRIMARY + propeller-logit cross-check; 16 sample units
  (genotype_batch); report 5 contrasts incl. interaction; concordance = the call.
- Pseudobulk DE: build_pseudobulk(replicate=genotype_batch) on RNA counts -> extend fit_limma_voom
  (voomWithQualityWeights + eBayes(robust=TRUE)) across 5 contrasts; whole-MG + per-substate (where cells
  allow); stageR screen-confirm FDR (+ per-contrast BH). filterByExpr per level.

## Steps (each closeable in one window; acceptance per step; [HEAVY]=compute, [DEP]=new pkgs+gate)
Each step: write pure fn + synthetic unit test (gate-independent) -> smoke-test on live data ->
wire target + full run -> verify quality gate (scripts/check.sh) before AND after.

- **S1 reprocess + cluster** [HEAVY][DEP: harmony, glmGamPoi].
  R/microglia.R::reprocess_microglia() (SCT recipe above); target `microglia_processed` (qs:
  SCT+pca+harmony+clusters+umap). rv add harmony glmGamPoi. Smoke-test reprocess on the live 26k object
  before the full target run.
  ACCEPT: reductions {pca,harmony,umap} + cluster factor present; DAM/IFN/homeostatic marker separation
  confirmed post-Harmony (else lower correction); deterministic re-run (fixed seed); gate green.

- **S2 substate annotation + QC prune** [DEP: UCell].
  Expand `canonical_microglia_markers` (constants.R) per SOTA marker lists above (+ MHC_APC aux);
  microglia.R::score_substates (UCell, SCT assay) + label argmax + prune_contaminant_clusters; target
  `microglia_annotated`. rv add UCell. Unit-test argmax + prune logic on synthetic fixture.
  ACCEPT: every retained cell labelled; per-substate UCell enrichment asserted (DAM cells high DAM score
  etc.); contaminant clusters dropped w/ logged counts+rationale; substate x genotype proportion table;
  Thrupp attenuation noted; gate green.

- **S3 composition** [DEP: sccomp(+cmdstanr/CmdStan, HEAVY install) + speckle].
  R/composition.R::test_composition (sccomp + propeller) across 5 contrasts; target `composition_results`.
  rv add sccomp speckle. cmdstanr downloads+compiles CmdStan (heavy; needs C++ toolchain - have
  build-essential; may need to `cmdstanr::install_cmdstan()` to a project-local path). Unit-test the
  contrast wiring on a synthetic count table.
  ACCEPT: per-contrast proportion estimates incl. interaction; amyloid->DAM shift quantified+tested;
  sccomp<->propeller concordance reported; replicate=genotype_batch; gate green.

- **S4 pseudobulk DE** [DEP: stageR(Bioc)].
  Extend R/de_pb.R (voomWithQualityWeights + eBayes(robust=TRUE); stageR helper; FIX the stale "de_sc.R"
  forward-ref comment -> single-cell DE dropped); targets `pb_de_microglia` (whole-MG) + `pb_de_substate`
  (per-substate) x 5 contrasts. rv add stageR. Extend tests/test_de_pb.R.
  ACCEPT: amyloid->DAM programme surfaces (DAM_up genes top in nlgf contrasts, direction matches v1
  qualitatively); interaction reported static-null honestly; topTables stored per contrast x level;
  thresholds stated (FDR<0.05 base, |logFC|>0.5 for "sig" calls); gate green.

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
- Per-substate pseudobulk thin for rare states (IFN/proliferative) -> filterByExpr + report n; substate DE
  is secondary to whole-MG.
- 8G load peak is the load_snrnaseq BUILD (cached); reprocess works on the 340MB subset -> moderate.
- Keep gate green each step; near 80% context -> drive to clean state + close out.
</content>
