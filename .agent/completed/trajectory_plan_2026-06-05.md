# Microglial activation-trajectory plan (arc M): dynamics layer

Adds the project's first **dynamics / progression** readout. Every prior arc reads
a **static** quantity: expression DE (02b), pathway/module enrichment (D), TF /
kinase / CCC / causal-topology / SCENIC activity (E/F/G/J/K), or — most recently —
cell **composition** (L, SpatialDecon). None reads the **ordering** of microglia
along an activation continuum. Here we infer a microglial activation trajectory
(homeostatic → DAM / IFN / proliferative) from the snRNAseq itself, assign each
cell a pseudotime, and test whether **amyloid (NLGF) advances microglia along the
homeostatic→DAM activation axis**, and whether the **tau background modulates that
advance** — the interaction `(NLGF_P301S−P301S)−(NLGF_MAPTKI−MAPTKI)` re-expressed
as a *progression* effect rather than an expression, activity, or composition one.

Why this is an independent readout, not a restatement of 02b/L: static substate DE
asks *which genes change within a fixed cluster set*; composition (L) asks *how many
cells fall in each discrete cluster*. A pseudotime asks *how far along the
activation program each cell sits on a continuum* — including within-cluster
advancement that neither DE nor discrete composition can see. The honest risk
(guardrail #1) is that a per-replicate mean-pseudotime shift could be **entirely a
composition shift** (more DAM cells ⇒ higher mean pseudotime), in which case it
reduces to L; the design therefore *decomposes* progression from reshuffling and
reports both. A trajectory that turns out additive (null interaction) — as the L
composition and human interaction layers already were — is a finding, not a failure.

Readout = per-cell pseudotime on a validated homeostatic→DAM lineage + a factorial
test (incl. interaction) on a per-**replicate** (genotype_batch, the locked unit)
pseudotime summary + an explicit progression-vs-composition decomposition + a
cross-method pseudotime concordance check. Does NOT touch the locked mouse 2×2
analysis except via an explicit, gated ledger-feed decision (M5).

## STATUS
- [x] M1 (GATE): DONE — gate answer = **Maximal triangulation** (Slingshot + PAGA/DPT + CellRank2 + tradeSeq, CytoTRACE2 root-check; user). Toolchain installed (slingshot/tradeSeq/princurve/TrajectoryUtils Bioc 3.22; cellrank 2.3.0 in .venv; scanpy 1.12.1 already present; CytoTRACE2-R from GitHub). Anti-anchoring design locked. End-to-end smoke OK (per-replicate power + slingshot homeostatic→DAM ordering). See M1 COMPLETION NOTE.
- [x] M2: DONE — `R/trajectory.R` built (10 pure fns: trajectory+pseudotime, root-potency compute+validate, per-replicate summary, progression-vs-composition decomposition, factorial fit, cross-method concordance, 3 plot fns), wired into `R/helpers.R` after `spatial_decon.R`. All 10 smoke-tested on the live cache (exit 0); 2 robustness bugs found+fixed. See M2 COMPLETION NOTE.
- [x] M3: DONE — `scripts/build_trajectory.R` (+ `scripts/build_trajectory_python.py` scanpy/CellRank lane) → `storage/cache/trajectory.rds` (15 slots) + 4 results TSVs. Validated (pt monotone, 16 reps, 5 contrasts/measure). **Headline: NLGF advances microglia along homeostatic→DAM pseudotime; tau supra-additively amplifies it (interaction logFC +2.46, FDR 0.077) — carried by PROGRESSION (+2.30 sig) not COMPOSITION (+0.13 null), so non-redundant with arc-L.** See M3 COMPLETION NOTE.
- [x] M4: DONE — `rmd/22_trajectory.Rmd` (§23, display-only) built + wired as `child-trajectory` before `child-session`. 5 subsections (23.1 topology+root, 23.2 density, 23.3 factorial forest+heatmap, 23.4 progression-vs-composition decomposition+concordance, 23.5 verdict). Reads `trajectory.rds` only; reuses the 3 R/trajectory.R plot fns + inline ggplot. Smoke-tested (12 builds, exit 0) then full knit: **0 error / 0 warning**, renders as §23 (session→§24, offset confirmed). Writes `results/trajectory_verdict.tsv`. See M4 COMPLETION NOTE.
- [x] M5 (GATE): DONE — gate answer (user, 2026-06-06) = **option B: 2 margin-neutral Suggestive rows, layer="dynamics"**. M-001 (amyloid_activation) supports Hyp-1A;Hyp-1B;T-Inflammation (mirrors K-001/L-001); M-002 (interaction_metabolic) records the POSITIVE+significant trajectory interaction thematically — supports **T-Synergy only**, no Hyp-3A/Hyp-3B support/contradict — so contests stay **18/12/55**. Built phase_m + extended layer-validation set ("dynamics"); rebuilt ledger (89 rows) + adjudication; **verified contests held 18/12/55**; synced rmd/16 §17 prose + fixed stale analysis.Rmd headline (was frozen at Phase J/83 rows; now 89, themes 30/30); re-knit **0 error / 0 warning**. See M5 COMPLETION NOTE.
- [ ] M6: close out — refresh `map.md` + `completed/DIGEST.md`, archive plan with date suffix, final clean knit + commit.

## Locked facts (verified at plan time 2026-06-05; do not re-derive)
- **Substrate object:** `storage/cache/microglia_seurat_processed.rds` — Seurat,
  **26,104 cells**, default assay **SCT**; assays RNA + SCT, layers counts/data/
  scale.data only. Reductions present: **pca, harmony, umap** (harmony = the
  batch-corrected embedding; use it for trajectory, NOT raw pca). `broad_annotations`
  = all "Microglia" (the `allen_labels` neuronal entries are noisy Allen-transfer
  artefacts on microglia, ignore). Meta has `genotype`, `batch`, `genotype_batch`
  (16 ids = the locked replication unit; ~1,631 microglia/replicate), `seurat_clusters`
  (13), `SCT_snn_res.{0.01,0.2,0.4,0.6}`, cell-cycle (`S.Score`/`G2M.Score`/`Phase`).
- **Substate labels are NOT cached.** The 4-level `state` factor
  {homeostatic, DAM, IFN, proliferative} is produced at runtime by
  `label_microglia_states(seurat_obj, symbol_map)` (R/microglia.R) via AddModuleScore
  argmax of 4 canonical signature scores. The build script MUST re-apply it (needs
  `symbol_map` ← `storage/cache/snrnaseq_symbol_map.rds`; `canonical_ids` ← constants.R).
  Deterministic; identical to the per-state NEBULA caches' assignment. (Verify the
  state×genotype_batch table for per-replicate power at M1 — IFN/proliferative may be
  thin in some replicates; the homeostatic→DAM lineage is the powered headline.)
- **RNA velocity is OUT.** No spliced/unspliced layers exist (RNA assay = counts/
  data/scale.data only) and snRNA spliced capture is unreliable; pseudotime only.
  Re-quantification from BAMs is not in scope.
- **R trajectory packages MISSING** (install at M1 per chosen tool): `slingshot`,
  `princurve`, `tradeSeq`, `CytoTRACE2`, `monocle3`, `destiny`. `SingleCellExperiment`
  is OK. Bioc = 3.22.
- **Python lane ready:** project `.venv` has **scanpy 1.12.1** + anndata (PAGA via
  `sc.tl.paga`, DPT via `sc.tl.dpt` — zero install). `cellrank` is NOT installed
  (needs `pip install cellrank` only if the CellRank2 option is chosen). Reuse the
  SCENIC-style microglia→h5ad export pattern if a Python pseudotime is built.
- **Replication / design (locked, project-standard):** unit = genotype_batch (16);
  headline interaction model on the per-replicate pseudotime summary via
  `de_pb.R::fit_limma_log` `~ 0 + genotype + batch` → `design.R::make_contrast_matrix`
  → the 5 canonical contrasts (`nlgf_in_maptki`, `nlgf_in_p301s`, `interaction`,
  `tau_alone`, `tau_in_nlgf`); FDR threshold 0.10. Section: rmd/22 → **§23**
  (verify the file≠§ offset at M4).

## M1 — design lock + env setup + smoke test (DECISION GATE at start)
**Gate:** choose the primary trajectory tool + headline readout. Present these
before AskUserQuestion (default first):
  1. **(DEFAULT) Slingshot pseudotime-shift.** R `slingshot` on the harmony
     embedding with the 4 `state`s as clusters, root = homeostatic; lineage topology
     (MST) + per-cell per-lineage pseudotime. Cross-check pseudotime with scanpy
     PAGA+DPT (free). Validate root with CytoTRACE2 (root-free potency). Headline =
     per-replicate **mean pseudotime on the homeostatic→DAM lineage** → limma 2×2+batch
     interaction. Rationale: R-native to the knit, cluster-based curves suit 4 defined
     substates, light installs, multi-method honesty (Slingshot + PAGA + CytoTRACE2).
  2. **CellRank2 fate-probability (SOTA alt).** scanpy PAGA+DPT + `cellrank` v2
     velocity-free (CytoTRACE/pseudotime kernel) → DAM-fate **absorption
     probabilities**; headline = per-replicate DAM-fate-probability interaction.
     More 2024-SOTA + probabilistic (macrostates, fate maps), velocity-free; but
     heavier Python dep (`pip install cellrank`), less R-native, readout shifts from
     pseudotime to fate-probability.
  3. **Maximal triangulation.** Slingshot + PAGA/DPT + CellRank2 + tradeSeq
     differential-dynamics (4-genotype condition smoothers + interaction contrast).
     Fullest house-style triangulation; most installs/compute; risk of over-build.
**Anti-anchoring locks (fix BEFORE inspecting which contrast moves; mirror K1/L1):**
root = homeostatic (validated, not asserted); embedding = harmony; lineage of
interest = homeostatic→DAM (IFN/proliferative lineages reported but secondary);
significance = FDR<0.10; replication unit = genotype_batch; contrast set = the 5
canonical. Velocity excluded.
**Setup + smoke:** install chosen toolchain; `chown rstudio:rstudio`; apply
`label_microglia_states`, print state×genotype_batch counts (power check), build a
minimal trajectory + pseudotime end-to-end on the live cache (exit 0), confirm a
sane homeostatic→DAM ordering. Record the gate answer + smoke result in the M1 note.

### M1 COMPLETION NOTE (2026-06-05)
**Gate answer:** Maximal triangulation — Slingshot (R, primary per-lineage
pseudotime) + scanpy PAGA/DPT (independent topology + pseudotime) + CellRank2
(velocity-free fate probabilities) + tradeSeq (differential dynamics), with
CytoTRACE2 root-free potency for root validation. The widest method set; the
smoke-test (below) shows WHY it is warranted here (close substates ⇒ method-
sensitive topology, so cross-method concordance is itself a result).

**Toolchain installed (all this session):** Bioc 3.22 — `slingshot`, `tradeSeq`,
`princurve`, `TrajectoryUtils` (exit 0, requireNamespace OK). `.venv` —
`cellrank` 2.3.0 (`import cellrank` OK); `scanpy` 1.12.1 already present (PAGA via
`sc.tl.paga`, DPT via `sc.tl.dpt`). `CytoTRACE2`-R installed from GitHub
(`digitalcytometry/cytotrace2`, subdir `cytotrace2_r`; exit 0, namespace loads;
pulled `HiClimR`). Documented fallback chain if it ever breaks: cellrank's
`CytoTRACEKernel` (free with cellrank) → per-cell RNA Shannon entropy /
n_expressed_genes as the potency proxy for root validation. Velocity
stays OUT (no spliced layers; snRNA capture unreliable).

**Anti-anchoring design LOCKED (before inspecting which contrast moves):** root =
homeostatic (to be *validated* by potency at M2, not asserted); embedding = harmony
(26104×30; **test 10–15-dim truncation at M2** — see implication 2); lineage of
interest = homeostatic→DAM (IFN/proliferative secondary); FDR<0.10; replication
unit = genotype_batch (16); contrast set = the 5 canonical.

**Smoke-test results (live cache).** (a) `label_microglia_states(s, symbol_map)`
re-derives the 4-level `state` deterministically (not cached): totals homeostatic
13550 / DAM 8774 / IFN 1760 / proliferative 2020. Per-replicate power: the headline
homeostatic→DAM lineage is well-powered in ALL 16 ids (min homeostatic 117, min DAM
91); IFN/proliferative are thin only in NLGF_P301S_batch01 (15 / 19) → reported but
secondary (guardrail #7). genotype×batch = 4×4 = 16 confirmed; harmony = 26104×30.
(b) Slingshot end-to-end on a 6000-cell subsample (38 s): topology Lineage1
homeostatic→proliferative→DAM, Lineage2 homeostatic→proliferative→IFN; root ordering
CORRECT (mean pseudotime on L1: homeostatic 37.0 < proliferative 38.2 < IFN 37.8 <
DAM 43.9 — homeostatic lowest, DAM highest).

**Design implications carried to M2/M3 (deviations/refinements):**
1. **Substates are transcriptionally close** (pseudotime spread only ~37→44 on the
   DAM lineage; matches the L-arc profile9 κ=35 collinearity). Consequence: the
   per-replicate mean-pseudotime shift will be substantially composition-driven ⇒
   the **progression-vs-composition decomposition (guardrail #1) and cross-method
   concordance (guardrail #6) are load-bearing, not optional**. Expect a modest /
   possibly null interaction (like L composition + human) — that is a valid finding.
2. **Cell-cycle confound.** Proliferative appears as a spurious *intermediate* on
   the full-30-dim harmony curves (proliferative = cycling; `S.Score`/`G2M.Score`/
   `Phase` exist). M2/M3: test fewer harmony dims (10–15), keep the headline on
   homeostatic→DAM, and treat proliferative as a cycling **side**-state (consider
   excluding it from the primary lineage and/or noting the cycle axis). Do NOT let a
   cycle-driven axis masquerade as activation progression.
3. Root-validation tool fallback chain documented above (implication: CytoTRACE2-R
   is preferred but not a hard dependency).

## M2 — `R/trajectory.R` helpers
Pure fns (no globals; smoke-testable in isolation). Proposed surface (finalise to
the chosen tool): `build_microglia_trajectory()` (harmony + states → lineages +
per-cell, per-lineage pseudotime), `validate_root_potency()` (CytoTRACE2 or entropy
→ confirm homeostatic = least-differentiated), `pseudotime_per_replicate()`
(genotype_batch × lineage → mean pseudotime + fraction-past-DAM-threshold +
within-state mean), `fit_trajectory_contrasts()` (reuse `fit_limma_log` +
`make_contrast_matrix` on the per-replicate summary), `decompose_progression_vs_composition()`
(partition the per-replicate shift into between-state composition vs within-state
advancement; guardrail #1), `pseudotime_concordance()` (Slingshot vs PAGA/DPT
Spearman), plot fns (trajectory on UMAP, pseudotime density per genotype,
interaction forest). Source in `R/helpers.R` AFTER `spatial_decon.R` (read its top
for dependency order). Smoke-test each fn against live caches via `Rscript -e`
before any knit.

### M2 COMPLETION NOTE (2026-06-05)
**Built `R/trajectory.R`** (411 lines, pure fns, no globals at source time;
sourced in `R/helpers.R` after `spatial_decon.R`, before `microglia.R`; doc block
added to helpers.R's source-order comment). Final surface (10 fns):
`build_microglia_trajectory(embedding, clusters, start_clus, end_clus,
omit_clusters, dam_state, seed)` runs Slingshot on a low-dim embedding with the 4
states as clusters; returns per-cell x lineage pseudotime/weights (NA off-lineage
or omitted), lineages, mst, the `sds` PseudotimeOrdering (for embedCurves/plots),
`terminal_states`, and the auto-identified `dam_lineage`. `cell_potency(counts,
method=entropy|n_genes)` (root-free proxy; higher=more potent) split OUT from
`validate_root_potency(potency, states, root)` (per-state mean-potency ranks +
`root_is_most_potent`/`delta_to_next`) so CytoTRACE2 (computed in M3) and the
entropy/n_genes proxies all flow through the SAME validator. `pseudotime_per_
replicate(pt, states, replicate, weights, dam_state, dam_threshold, state_levels,
replicate_meta)` returns per-id mean_pt + frac_past (default threshold = global
median pt of DAM cells) + within_<state> means. `decompose_progression_vs_
composition()` returns per-id observed / composition_cf / progression_cf (all
pseudotime scale, fit together at M3). `fit_trajectory_contrasts(summary_mat, meta,
transform, padj_cut)` reuses `factorial_design`+`fit_limma_log` (mirrors
`fit_aucell_contrasts`), carries `se`(=logFC/t) for the §23 forest CIs.
`pseudotime_concordance(pt_df)` pairwise Spearman. Plots: `plot_trajectory_umap`,
`plot_pseudotime_density` (project palette), `plot_interaction_forest` (5-contrast
logFC +/- 1.96 se).

**Smoke (8,000-cell stratified subsample, 15 harmony dims, omit proliferative;
/tmp/smoke_traj.R, exit 0):** all 10 fns validated. Two robustness bugs found+
fixed: (1) `decompose` produced NaN when a `state_level` was entirely absent on
the lineage (pi_bar=0, mu_bar=NA, 0*NA poisons the sum) -> now intersects
state_levels with on-lineage observed states; (2) `fit_trajectory_contrasts` now
drops non-finite measure rows (limma-trend eBayes rejects an NA/Inf covariate) with
a naming warning. Also swapped `geom_errorbarh` (soft-deprecated, ggplot2 4.0.2)
for `geom_errorbar(orientation="y")` so M4 knits warning-free (verified via
ggplot_build). `trajectory.R` sources standalone with 0 warnings (the 1 helpers.R
source-warning is the pre-existing nichenetr `e1071::element` namespace notice). No
function-name collisions with later-sourced files.

**Findings carried to M3 (NOT final results -- subsample diagnostics):**
1. **Topology is method-sensitive (reinforces M1).** 15 dims + omit proliferative
   gave a single lineage homeostatic->**IFN**->DAM (IFN now the spurious
   intermediate, as proliferative was at 30 dims). Root ordering still correct
   (mean pt DAM 53.7 > IFN 50.0 > homeostatic 46.2; DAM terminal). M3 policy must
   explore n_dims AND which side-states to omit, and likely report a CLEAN
   homeostatic->DAM variant (omit BOTH proliferative+IFN, or `end_clus`) alongside
   the all-states topology (guardrail #7). The fns are deliberately un-opinionated:
   policy (n_dims, omit_clusters, end_clus) lives in the M3 build script, fixed
   before inspecting which contrast moves (guardrail #8).
2. **Crude potency disagrees with the asserted root (guardrail #2 is live).** Both
   entropy and n_genes rank homeostatic only 3rd (IFN>DAM>homeostatic>prolif) --
   the expected confound (activated microglia induce many response genes ->
   inflated crude potency). CytoTRACE2 (M3) is the load-bearing arbiter; if it ALSO
   disagrees, report honestly: microglial "homeostatic" is a mature RESTING state,
   not a stem state, so the axis is ACTIVATION ordering, not developmental potency
   (ties to guardrail #3, pseudotime != time). Do not force the root to validate.
3. **Decomposition discriminates channels (the machinery works).** On the
   subsample the interaction was carried by `progression_cf` (logFC 2.99, sig) not
   `composition_cf` (0.31, ns); 2nd-order interaction residual small (max |.| 0.76
   vs a ~3.3 mean shift, so the additive partition holds). This is ONLY a mechanics
   demonstration on a confounded subsample/topology -- the real read is M3's.
4. **/tmp/smoke_traj.R is disposable** (outside the repo; not committed).

## M3 — `scripts/build_trajectory.R` → cache + TSVs (out-of-knit)
Idempotent. Reads microglia_seurat_processed + snrnaseq_symbol_map; applies
`label_microglia_states`; builds trajectory (+ optional Python PAGA/DPT/CellRank2
sub-step exporting microglia→h5ad, reuse SCENIC export pattern); computes per-cell
pseudotime, root validation, per-replicate summaries, factorial contrasts, the
progression-vs-composition decomposition, cross-method concordance. Writes
`storage/cache/trajectory.rds` (multi-slot list: lineages, per-cell pseudotime,
potency, per-replicate summary, contrast fits, decomposition, concordance, params)
+ results TSVs (`trajectory_pseudotime_by_genotype.tsv`, `trajectory_contrasts.tsv`,
`trajectory_progression_decomposition.tsv`, `trajectory_method_concordance.tsv`).
`chown`. Validate: pseudotime finite/monotone vs state means, all 16 replicates
summarised, contrast table has 5 rows × {logFC,t,P,adj.P}.

### M3 COMPLETION NOTE (2026-06-06)
**Built** `scripts/build_trajectory.R` (thin out-of-knit orchestrator, house
build pattern) + `scripts/build_trajectory_python.py` (scanpy PAGA/DPT + CellRank2
lane, SCENIC-style R→Python CSV bridge). Outputs (all chowned rstudio:rstudio):
`storage/cache/trajectory.rds` (15 slots: traj_clean/traj_all Slingshot runs,
per_cell [26104×14], potency, root_validation, per_replicate, decomposition,
summary_matrix, contrasts, concordance, by_genotype, differential_dynamics,
python_provenance, replicate_meta, params) + 4 TSVs `trajectory_contrasts.tsv`
(45=9 measures×5 contrasts), `_pseudotime_by_genotype.tsv` (36), `_progression_
decomposition.tsv` (16 reps), `_method_concordance.tsv` (6). Runtime ~12 min
(CytoTRACE2-dominated). Locked params: 15 harmony dims, n_neighbors 30, n_dcs 10,
seed 1, clean lineage omits IFN+proliferative, dam_threshold 33.6, design
`~ 0 + genotype + batch`, FDR<0.10.

**HEADLINE RESULT (the real read; M2 subsample was only mechanics).** Amyloid
(NLGF) strongly advances microglia along the homeostatic→DAM activation pseudotime:
`nlgf_in_maptki` mean_pt logFC +7.87 (p 2.7e-7), `nlgf_in_p301s` +10.32 (p 1.3e-8);
tau alone does NOT (`tau_alone` −1.20, ns). The **tau (P301S) background supra-
additively amplifies the amyloid advance**: interaction logFC +2.46, FDR 0.077
(≈ 10.32−7.87). The progression-vs-composition decomposition (guardrail #1) is the
load-bearing result: the interaction is carried almost entirely by **progression**
(`progression_cf` +2.30, FDR 0.077 — genuine within-state advancement, ~94%) and
**not composition** (`composition_cf` +0.13, FDR 0.81 — cluster reshuffling, null).
So the dynamics layer reveals a tau×amyloid effect that is invisible to discrete
composition — and the null composition-interaction is exactly consistent with arc-L
(SpatialDecon) and the human/mouse interaction layers, confirming this readout is
NON-REDUNDANT, not a restatement. Corroborating channels: `frac_past` interaction
+0.19 (FDR 0.013, strongest; absolutes 4%→40%[MAPTKI bg]→59%[P301S bg] of cells
past the DAM-threshold) and `within_homeostatic` +2.68 (sig) > `within_DAM` +1.71
(ns) — the amplification pushes even homeostatic-labelled cells rightward. The
amyloid MAIN effects, by contrast, split across BOTH channels (composition +1.76/
+1.90 sig AND progression +6.69/+8.99 sig), i.e. amyloid both makes more DAM cells
AND advances them; tau only adds the progression piece. mean_pt absolutes: MAPTKI
25.3 / P301S 24.1 / NLGF_MAPTKI 33.1 / NLGF_P301S 34.4 (highest). `observed`≡
`mean_pt` to the digit (law-of-total-expectation consistency check passed).

**Two honest findings reported, not forced (guardrails #2/#3, #6).** (1) Root
potency: all THREE proxies reject homeostatic as most-potent — entropy rank 3/4,
n_genes 3/4, CytoTRACE2 4/4 (the load-bearing arbiter), all `most_potent=FALSE`.
As M2 anticipated, microglial "homeostatic" is a mature RESTING state, so the axis
is an ACTIVATION ordering, not developmental potency/time. The script reports the
disagreement; the trajectory direction rests on the validated Slingshot
homeostatic→DAM monotonicity (DAM mean pt 35.1 > homeostatic 27.7) + canonical
marker biology, NOT on a potency claim. (2) Cross-method concordance: Slingshot-
clean vs Slingshot-all ρ=0.98 (nested), Slingshot vs **CellRank ρ=0.57** (an
independent fate-probability method corroborates the ordering), but Slingshot vs
**DPT ρ=−0.09** (and DPT-vs-CellRank −0.22). DPT is the outlier: diffusion
pseudotime has weak dynamic range on this graded activation continuum and is pulled
by the peripheral IFN/proliferative geometry (mean DPT ranks them latest). So 2 of
3 independent methods corroborate; the 3rd is a known-weak estimator here — exactly
the kind of disagreement maximal triangulation (M1 gate) is meant to surface.

**Deviations from the M3 spec / design (all deliberate, justified here).**
(a) **CellRank route changed.** The first build hit a hard wall: the planned
automatic GPCCA macrostate discovery (`compute_schur(20)` + `compute_macrostates
(n_states=[3,8])`) densified to an O(n³) real-Schur on the 26k transition matrix
(7–8 BLAS threads, 32 GB, >45 min, killed). Replaced with the **principled cheap
route**: the biology fixes the terminal set, so set terminal states MANUALLY to
{DAM, IFN} and solve the SPARSE absorption-probability system (55 s, all 26104
cells finite). This is more honest, not a shortcut — it uses the known end-states
rather than re-deriving them. Auto discovery retained behind `--cellrank-auto` with
a SPARSE Krylov Schur + bounded single n_states so it can never regress to the
dense monster. Added a DPT-only checkpoint CSV write (a later CellRank stall can no
longer cost the computed DPT). (b) **tradeSeq GATED OFF** by default
(`differential_dynamics=NULL`); the per-cell/per-replicate progression readout is
the headline and tradeSeq's gene-level differential dynamics is a separate, slower
question — enable with `--tradeseq` at M4 only if §23 wants the gene panel.
(c) **igraph** installed in `.venv` for scanpy PAGA; PAGA wrapped non-blocking
(a failure never costs DPT/CellRank). PAGA succeeded (state order returned).
(d) **Two-scale limma:** pseudotime-scale measures (mean_pt, observed, composition_
cf, progression_cf, within_*) and bounded [0,1] measures (frac_past, cellrank_dam_
fate) fit in SEPARATE eBayes calls — mixing scales would distort limma-trend's
shared variance prior. (e) **DPT computed all-states** (it must, to feed the
CellRank PseudotimeKernel); quarantined to the concordance sidebar — every
downstream result (decomposition, per-replicate, contrasts) uses Slingshot-clean.
**For M4:** plot Slingshot-clean pt as the headline; present DPT with the honest
weak-concordance caveat; the forest plot's lead panel should be the interaction
decomposed by channel (progression vs composition), since that IS the result.

## M4 — `rmd/22_trajectory.Rmd` (§23, display-only)
Subsections: 23.1 trajectory topology + root validation (UMAP + CytoTRACE2);
23.2 pseudotime distribution per genotype (density/ridge by 2×2); 23.3 factorial
interaction test (forest of the 5 contrasts on the homeostatic→DAM lineage);
23.4 progression-vs-composition decomposition + cross-method concordance; 23.5
sign-aware verdict prose. Reads `trajectory.rds` only (recomputes nothing); plot
fns live in the chapter's own R module / `R/trajectory.R` so they render inside the
shared knit session (mirror §15/§20). Wire `child-trajectory` before
`child-session` in analysis.Rmd. Verify file≠§ offset. `chown` root-owned outputs;
`grep -c 'class="error"' analysis.html` and `class="warning"` both 0.

### M4 COMPLETION NOTE (2026-06-06)
**Built** `rmd/22_trajectory.Rmd` (display-only, ~250 lines, mirrors the rmd/21
spatial-decon verdict-chapter template) and wired `child-trajectory` between
`child-spatial-deconvolution` and `child-session` in `analysis.Rmd`. The chapter
`readRDS`-loads `trajectory.rds` in its first chunk and **recomputes nothing**;
all prose scalars are read from the cache via two unique-row lookup closures
(`ct_get(measure, contrast)`, `bg_get(measure, genotype)`) so the text can never
drift from the numbers. Five subsections per the M4 spec: **23.1** topology + root
validation (two UMAPs via `plot_trajectory_umap` — state + clean pseudotime — plus
a 3-proxy potency rank table); **23.2** per-cell pseudotime density by genotype
(`plot_pseudotime_density`); **23.3** the 5-contrast forest on mean_pt
(`plot_interaction_forest`) + a 9-measure × 5-contrast logFC heatmap; **23.4** the
load-bearing progression-vs-composition decomposition (a faceted
composition|progression forest of the amyloid + interaction contrasts) + the
cross-method Spearman concordance heatmap; **23.5** sign-aware verdict prose +
`results/trajectory_verdict.tsv` (3-axis, same schema as `spatial_decon_verdict.tsv`).

**Render mechanics / house-rule compliance.** The 3 plot fns already live in
`R/trajectory.R` (sourced via helpers.R in the shared knit session, so they render
in-session like §15/§20 — no scientific-visualization). Inline plots use bare
`ggplot()` (ggplot2 is attached globally) and the project `diverging` palette;
forests use `geom_errorbar(orientation="y")` (NOT the soft-deprecated
`geom_errorbarh`, which the M2 note flagged for ggplot2 4.0.2). Section offset
**verified**: pre-insert, session was §23; post-insert the trajectory chapter is
**§23** and session shifted to §24 (file 22 → §23, the documented file≠§ offset).
Smoke-tested every chunk + every plot build against the live cache
(`/tmp/smoke_rmd22.R`, disposable, 12 OK builds, exit 0) BEFORE the knit; the only
warning is the pre-existing nichenetr `e1071::element` namespace notice at
source-time. Full knit (5 min): `grep -c 'class="error"'` = **0**,
`class="warning"` = **0**. `chown rstudio:rstudio` on `rmd/22_trajectory.Rmd` and
on the root-written knit outputs; also restored `analysis.Rmd` ownership (the Edit
tool rewrote it as root). `analysis.html`, `storage/results/*` and `omnipathr-log/`
are all gitignored (confirmed via `git check-ignore`), so no stray artefacts.

**Biology surfaced (display only — the result is M3's; this chapter presents it
honestly).** The headline reads cleanly in-document: amyloid advances microglia
along the homeostatic→DAM activation pseudotime in both tau backgrounds
(`nlgf_in_maptki` +7.87, `nlgf_in_p301s` +10.32, both FDR<1e-5), tau alone does
not (−1.20 ns), and the P301S tau background supra-additively amplifies the amyloid
advance (interaction +2.46, FDR 0.077). The decomposition figure is the chapter's
load-bearing visual: the interaction loads on **progression** (+2.30, FDR 0.077,
94% of total) and **not composition** (+0.13, FDR 0.81), with the amyloid main
effects splitting across BOTH channels — making the dynamics interaction
non-redundant with the arc-L composition null. Root validation is reported
honestly (homeostatic ranks last/near-last across all 3 potency proxies → the axis
is an ACTIVATION ordering, not developmental potency) and cross-method concordance
is shown in full with DPT as the flagged weak outlier (2/3 methods corroborate).

**No deviations from the M4 spec.** The chapter is purely display/read-only as
specified; the verdict-TSV write is the established house pattern (§20/§21 each
write a `*_verdict.tsv`). Next: **M5 is a GATE** (ledger-feed decision) — ended
the session cleanly here per the execution model; M5 needs user confirmation.

### M5 COMPLETION NOTE (2026-06-06)
**Gate answer (user): option B — 2 margin-neutral Suggestive rows, `layer="dynamics"`.**
Fed the trajectory evidence into the biological-model ledger as a new evidence
class without disturbing any contest verdict.

**What was built.** `scripts/build_biological_model_ledger.R` gained a `phase_m`
block (rows **M-001**, **M-002**) and `"dynamics"` was added to the `layer`
validation `stopifnot` set (the only builder change needed — the adjudication
script ignores `layer`, per the L precedent). M-001 (axis amyloid_activation):
the amyloid MAIN effect advances homeostatic→DAM pseudotime in BOTH tau
backgrounds (mean-pt logFC +7.87 / +10.32, FDR<1e-5; frac-past-DAM 4%/4%→40%/59%);
supports Hyp-1A;Hyp-1B;T-Inflammation, contradicts none — margin-neutral by the
K-001/L-001 both-amyloid-models construction. M-002 (axis interaction_metabolic):
the project's FIRST positive+significant interaction in a non-expression layer
(interaction mean-pt +2.46 FDR 0.077; frac-past +0.19 FDR 0.013), carried by
within-state PROGRESSION (+2.30, ~94%) not composition (+0.13 null, the arc-L
channel); recorded THEMATICALLY against **T-Synergy only**, no Hyp-3A/Hyp-3B
support/contradict. Both Suggestive (single inferred ordering; activation-not-
potency root; DPT cross-method outlier — guardrail #1).

**Rebuild + verification.** Ran `build_biological_model_ledger.R` (→ 89-row
`biological_model_claims_ledger.tsv`, layer dist now includes `dynamics=2`) then
`build_biological_model_adjudication.R` (→ adjudication + contest_verdicts TSVs).
**Contests held exactly 18 / 12 / 55** (Hyp-1B 33−15=18, Hyp-2B 8−(−4)=12,
Hyp-3B 28−(−27)=55) — confirmed in both the TSV and the rendered §17 prose. Theme
net-support: T-Inflammation 29→30, T-Synergy 29→30 (each +1 from one Phase M row);
all other entities unchanged.

**Prose sync (rmd/16 §17).** Updated: intro row-count 87→89 + Phase M clause;
datatable comment + `lengthMenu` 87→89; axis-1 margin arithmetic (32−14)→(33−15)
and the K-001/L-001 both-supporting-rows list extended with M-001 (Hyp-1A
supporting rows 23→24); cross-entity ranking (Hyp-1B 33, T-Synergy 30,
T-Inflammation 30, Hyp-1A net 15); per-axis-anchoring percentages (27/30 for both
top themes, still 90%); §17.5 verdict theme ranking (30/30); and a new Phase M
paragraph after the Phase L one (full M-001/M-002 description + the option-B gate
rationale + the Suggestive caveats).

**Deviation (deliberate, +scope):** also corrected the `analysis.Rmd` "Headline
integrated biological model" section, which had drifted — it was frozen at Phase J
("83 atomic evidence claims", themes T-Synergy 29 / T-Inflammation 26); Phase K
and Phase L were never synced there. Phase M widened that gap into an internal
inconsistency with §17 inside the same rendered HTML, so I brought the two factual
drift points current (89 rows + full phase parenthetical through Phase M; theme
ranking 30/30). The biology narrative and the 18/12/55 margins in that headline
were already correct and untouched. Flagging because the M5 plan scoped prose sync
to "rmd/16 §17" only; this is the accuracy-over-tight-scope call (CLAUDE.md).

**Knit gate: 0 error / 0 warning.** Full re-knit; analysis.html shows §23
trajectory / §24 session, ledger 89 rows × 14 cols, `dynamics=2`, margins 18/12/55.
Chowned the three Edit-touched files (analysis.Rmd, rmd/16, ledger script) +
storage/cache/trajectory/* + the two notes back to rstudio:rstudio. Next: **M6**
close-out (map.md, DIGEST.md, archive plan, final knit, commit).

## M5 — ledger-feed decision (DECISION GATE)
Same shape as K6/L5. Default (present + ≥1 alt): **add 2 Suggestive rows** in a new
`layer="dynamics"` evidence class — one for the amyloid-driven homeostatic→DAM
**progression** (supports BOTH amyloid models Hyp-1A;Hyp-1B + T-Inflammation, like
K-001/L-001 ⇒ margin-neutral), one for the **interaction-on-trajectory** result
(if positive: support Hyp-3B;T-Synergy + contradict Hyp-3A; if null: empty/empty
honest non-corroboration, like L-002). Alternatives: 1 row only; or no feed
(standalone chapter). Anti-double-count: scope rows to *trajectory progression /
dynamics*, explicitly distinct from the substate expression DE (already counted)
and the discrete composition (L). Suggestive (not Strong) BECAUSE pseudotime is a
single inferred ordering partly confounded with composition (guardrail #1). The
11-entity set + grade definitions stay untouched; `layer="dynamics"` added to the
builder's layer validation only (the adjudication script ignores `layer`, per L).
If feeding: rebuild ledger + adjudication, **verify contests stay 18/12/55** (unless
a deliberate, user-approved move), sync rmd/16 §17 prose, re-knit clean.

## M6 — close out
Refresh `map.md` (new rmd/22 row, R/trajectory.R entry, trajectory.rds cache row,
ledger phase-M note if M5 fed rows) + `completed/DIGEST.md` (arc-M digest in the
established schema) + archive this file to `completed/trajectory_plan_2026-06-05.md`.
Final clean knit; commit.

## Execution model (per-step loop; mirrors the project standard)
0. Check for a fitting Skill/subagent first (scanpy/scvi-tools/anndata for the
   Python lane; statistical-analysis for the interaction model; scientific-
   visualization is NOT used for in-knit plots — they must render in the shared R
   session, per §15/§20).
1. Smoke-test new helpers / cache-readers against live caches via `Rscript -e`
   BEFORE knitting (knit ≈ 3 min).
2. Mark the step DONE in STATUS + write a multi-paragraph COMPLETION NOTE (what was
   built, paths, biology, explicit deviations).
3. `chown rstudio:rstudio` new files + the knit's root-owned outputs.
4. Re-knit; `grep -c 'class="error"' analysis.html` = 0 (same for warning).
5. Commit locally (imperative subject <70 chars + HEREDOC body + Co-Authored-By
   trailer). Gates need user confirmation via AskUserQuestion.
Design-heavy steps (M1 lock, M2 metric design, M3 cache shape) get a fresh session;
chain only trivial pattern-extensions. End cleanly at a gate, a design fork, or low
context.

## Anti-anchoring guardrails
1. **Progression ≠ composition.** A per-replicate mean-pseudotime shift may be pure
   reshuffling into the (already-known, L-arc) larger DAM cluster. ALWAYS report the
   decomposition (between-state composition vs within-state advancement) and the
   within-DAM interaction; if the shift is fully composition-driven, say so (it
   reduces to L, not a new finding).
2. **Root validated, not asserted.** Confirm homeostatic = least-differentiated with
   a root-free potency method (CytoTRACE2 / RNA entropy); report disagreement.
3. **Pseudotime ≠ time.** It is an activation-state ordering, not chronological or
   developmental time; frame as "activation progression" throughout.
4. **Null interaction is a finding.** L composition + human interaction were both
   null/under-identified; an additive (non-synergistic) trajectory drive is a valid,
   reportable outcome — do not fish for significance or retune the lineage.
5. **Batch.** Use the harmony embedding + batch in the limma design; check the
   trajectory's leading axis is activation, not a batch axis.
6. **Multi-method honesty.** Report cross-method pseudotime Spearman; discordance is
   a finding, not smoothed (mirror K's head-to-head honesty).
7. **Report all states / all lineages.** Absence of an IFN or proliferative effect
   is evidence (anti-anchoring); the homeostatic→DAM lineage is the powered headline
   but is not the only one reported.
8. **Locks fixed before inspection** (M1): root, embedding, lineage-of-interest,
   threshold, unit, contrast set — all before seeing which contrast moves.
