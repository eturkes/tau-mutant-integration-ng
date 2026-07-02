# Roadmap - forward index (read first each session)

Posture · cohesive story · active plan · phased backlog. Companions in `.agent/`:
`memory.md` (standing contract), `map.md` (wiring, grows with the build),
`history.md` (new-project decision digests), `archive_digest.md` (v1 reference -
mine for promising threads; full v1 code on branch `archive`).

## Posture: BUILD (fresh)
v1 (23 rmd / 11 arcs D-O / 92-row ledger / 3 prose-heavy reports) is archived on
branch `archive`. `main` rebuilds from scratch: leaner pipeline, ONE cohesive
story, minimal prose, project-local env, no Docker. Re-derive value, drop bloat.

## Cohesive story (the spine - proposed, revisable)
Data = 4 mouse AD genotypes in a 2x2 (tau: MAPTKI ~ tau-KO vs P301S mutant; amyloid:
-/+ NLGF) x 4 modalities (snRNAseq, GeoMx spatial, bulk proteome, bulk phospho).
Question = how amyloid reshapes microglia under each tau background. Divergence =
interaction (NLGF_P301S - P301S) - (NLGF_MAPTKI - MAPTKI).
Durable findings mined from v1 (the headline to rebuild around):
1. Amyloid drives a microglial homeostatic->DAM activation programme (DAM / MG-M2 /
   antigen-presentation); microglia-led, corroborated across modalities.
2. Mutant tau MODULATES the amyloid response (the interaction) - clearest as a
   compositional DAM-cell synergy: amyloid adds more DAM cells on the mutant-tau
   background than additivity predicts. P2 does NOT detect a supported
   progression-beyond-composition signal along the activation trajectory.
3. Mechanism converges on Gsk3b / Myc signalling; tau attenuates amyloid-driven NF-kB.
4. Secondary axis: amyloid-linked synaptic suppression + TREM2 / APP clearance.
Streamlined OUT (v1 bloat): the 11-arc ledger + contest machinery, the margin-neutral
corroboration arcs (SCENIC, spatial-decon, celltype-specificity, gene-level dynamics),
the human-validation layer, the capstone convergence matrix, the heavy prose.

## Active plan: P3 Mechanism
Plan: `.agent/p3_mechanism_plan.md` (opened 2026-07-02 after user confirmed default P3).
Next `$session-prompt` mode = EXECUTE: implement S3.

Steps:
- [x] S1 dependencies + API contracts.
- [x] S2 RNA pathway + TF + NF-kB targets.
- [ ] S3 minimal phosphosite DE + kinase activity.
- [ ] S4 mechanism report + integration.

## Backlog - phased build (each phase = closeable increments; mine archive_digest per phase)
- P0 Foundations [DONE 2026-06-29]: project-local env (rv for R + uv .venv for Python), shared
  helpers (io / design+contrasts / plot theme), data load + QC sanity, 2x2
  factorial + 5 contrasts, concrete quality gate.
- P1 snRNAseq microglia core [DONE 2026-06-30]: reprocess (SCT) + substates (homeostatic / DAM /
  IFN / proliferative, UCell), composition (propeller primary + sccomp cross-check), pseudobulk DE across
  contrasts -> the robust amyloid->DAM activation programme (3-way confirmed: composition + DE + UCell score).
  Interaction = no large-effect DE, under-powered NOT absent (reported w/ MDE/CI + 123 stageR-confirmed, real |logFC| but sub-threshold per-contrast FDR) ->
  synergy handed to P2 and resolved as DAM-cell composition rather than progression/acceleration. Digest -> history.md.
- P2 Interaction trajectory [DONE 2026-07-02 -> `.agent/completed/p2_trajectory_plan_2026-07-02.md`]:
  activation pseudotime (homeostatic->DAM); amyloid advances the activation axis, but the tau x amyloid
  interaction decomposes to more DAM cells, not supported progression beyond composition.
- P3 Mechanism: focused pathway/module survey; TF (decoupleR / CollecTRI) + kinase
  (decoupleR / OmniPath) -> Gsk3b / Myc; NF-kB attenuation check.
- P4 Cross-modality: GeoMx spatial DE (+ light deconvolution if it earns it),
  proteome + phospho DE, CCC for the synaptic/clearance axis, integrated divergence view.
- P5 Synthesis: ONE lean report - cohesive narrative + compact evidence table (no
  ledger / contest machinery).

## Ledger (trajectory)
- 2026-06-29 archived v1 -> branch `archive`; opened fresh orphan `main`; reset
  `.agent` docs + initial agent config; reframed history as `archive_digest.md`; drafted
  this streamlined phase plan.
- 2026-06-29 S4 DONE (report engine closed): pivoted Quarto book -> ONE standalone offline HTML
  (`index.qmd` format:html embed-resources + `{{< include _qc.qmd >}}` + `theme.scss`); theme = crimson
  colours (#B0344D) + IBM Plex via 9 `@font-face` (relative `url("assets/fonts/<n>.woff2")` -> Quarto
  base64-inlines each woff2 under embed-resources). Render PROVED: 9 faces inlined OFFLINE (d09GMg magic),
  0 external loads, 0 error/0 warning (`tar_meta` all-NA), QC bounds pass (16x16 genotype-batch bijection).
  9 woff2 COMMITTED (assets/fonts/, avoid direct reads per AGENTS read economy), listed in
  `tar_quarto(extra_files=)`. Added `R/plot.R` + device-free `tests/test_plot.R`; de-staled map.md
  (book -> standalone report wiring) + S1 plan refs. Reasoning reversal (reviewed): the earlier
  "IBM Plex can't work via theme.scss" was a MEASUREMENT artifact -- once `@font-face url()` is present
  the whole theme CSS URL-encodes, so raw greps for `IBM Plex`/colours read ~0; match encoded tokens.
- 2026-06-29 S5 DONE -> P0 CLOSED: concrete quality gate `scripts/check.sh` (fail-loud, zero-fault) =
  `rv sync` + `uv sync` | loop `tests/test_*.R` at warn=2 | `tar_make()` tee'd to a log | enforce
  `tar_meta(error,warnings)` all-NA SCOPED to `tar_manifest()$name` (drops functions/globals + stale dead
  rows) + render-log grep for Quarto/pandoc/knitr `warning`/`warn` (knitr's separate process bypasses
  tar_meta). Negative-tested (grep pattern + `stopifnot(FALSE)` -> exit 1); full gate green end-to-end.
  Locked memory.md gate (provisional -> concrete); finalised map.md wiring; folded P0 digest -> history.md;
  archived plan -> `.agent/completed/`. P0 foundation complete; next = P1.
- 2026-06-29 P0-S5 gate review-hardened (post-close): the pre-review render-log grep was BLIND to
  render warnings -- `tar_quarto` defaults quiet=TRUE (suppresses Quarto/Pandoc warnings from the log) and
  knitr renders R chunk warnings INTO the HTML (never to log or tar_meta). Fixes: `_qc.qmd` setup
  `options(warn=2)` (chunk warning -> render error -> tar_make fails); `_targets.R` `tar_quarto(quiet=FALSE)`
  (warnings reach the log); check.sh FORCE-renders the report each run (closes the cached-clean blind spot;
  cheap -- reads cached ~0.3G, not the 8G load_snrnaseq build), anchored render-log grep with exit-code
  discrimination (no false-reds/false-greens), tar_meta scoped to manifest + dynamic branches, `if !`-wrapped
  tar_make. Negative-tested (chunk warning -> red; clean -> green); report still renders 0-warning under the
  stricter checks. Docs corrected (memory/map/history).
- 2026-06-29 P1 OPENED -> `.agent/p1_snrnaseq_plan.md`. Research: mined archive_digest + v1 archive branch
  (Explore) + 2026 SOTA web sweep (4-agent) -> reconciled v1 vs SOTA. Decision gate (3 forks): single-cell
  DE DROPPED (pseudobulk sole inference, Squair/Murphy-Skene; v1 NEBULA interaction was null); composition =
  sccomp + propeller; normalisation = SCT-v2 (v1 continuity, user choice over SOTA logNorm lean). Locked
  defaults: UCell scoring, Harmony batch-only (sex aliased), 4 argmax substates + aux MHC/APC, voomQW +
  robust eBayes + stageR. Key reframing from mining: the tau x amyloid interaction is STATIC-NULL (v1 matched-
  power null even single-cell) -> the synergy is a position/extent effect (rate only under the age-matched snapshot assumption), deferred to P2 (trajectory); P1
  nails the robust amyloid->DAM headline + substates + composition. Carry Thrupp 2020 (snRNA under-detects
  ~18% DAM activation genes). 5 steps S1-S5; next = S1 reprocess+cluster.
- 2026-06-29 P1-S1 DONE -> `microglia_processed` target + `R/microglia.R` (reprocess_microglia +
  marker_mean_by_cluster) + `tests/test_microglia.R`. SCT-v2/glmGamPoi -> Harmony(batch-only) -> Louvain
  multi-res {0.2,0.4,0.6}, primary=0.4 (12 clusters) -> UMAP, on the live 26k subset (138s, 687MB qs). Three
  pkg-drift fixes vs v1 recipe: harmony 2.0 dropped `assay.use` (assay implicit in reduction.use); SCTransform
  trips future's 500MiB globals cap (raise it); RunUMAP's "default changed" NOTICE is a WARNING that would fail
  the gate via tar_meta -> silenced by Seurat's own `Seurat.warn.umap.uwot=FALSE` (other warnings still
  surface). Also strips stale upstream meta shadows (pca1/umap1 coords-as-cols, SCT_snn_res.0.01). Acceptance:
  reductions {pca,harmony,umap}+cluster factor present; post-Harmony marker separation confirmed (distinct
  argmax per substate); re-run STABLE (observed ARI=1.0 under recorded threads; clusters derive from the
  seeded PCA->Harmony->SNN->Louvain chain, UMAP is viz-only; NOT bitwise-guaranteed per up-to-tolerance); gate green.
- 2026-06-30 P1-S2..S4 DONE (annotate / composition / pseudobulk-DE; outcomes folded to history.md + memory.md P1
  sections, full detail in the archived plan) -> P1-S5 DONE -> P1 CLOSED. S5 = `_microglia.qmd` microglia chapter
  (substate UMAP, composition forest/table, amyloid->DAM volcano + DE counts, under-powered-interaction + P2
  pointer, Thrupp + dropout caveats) reading a COMPACT `microglia_report` target (microglia_report_data extracts
  umap+substate+z+prune/provenance from the 612MB annotated -> the gate's force-render stays cheap ~12s). Prose
  INLINE-COMPUTED from targets (never hardcoded). Adversarial self-review fixed 3 prose-accuracy issues (interaction
  df.total-vs-residual-df conflation; sccomp CI reported factually, not asserted "spanning zero"; dropout caveat
  backed against the DAM-ENRICHED ceiling -- dropped clusters' DAM medians overlap kept HOMEOSTATIC, so the floor
  comparison was false). Full gate green end-to-end (6 tests warn=2 + force-render 0-warning + tar_meta clean across
  19 targets). Folded P1 digest -> history.md; updated memory (P1-S5 + cheap-render + warn=2-per-section) + map
  (microglia_report target + _microglia.qmd include + microglia_report_data fn); archived plan ->
  `.agent/completed/p1_snrnaseq_plan_2026-06-30.md`; reset Active plan to (none). Next = PLAN P2.
- 2026-06-30 P1-S5 review PARKED (context-overflow remediation): Review of the S5 close raised 14
  findings, ALL accepted (interaction "small-effect" over-claim -> "sub-threshold-per-contrast"; concordance /
  pruning / MDE prose accuracy; extractor finite+consistency guards + 4 negative tests; stale ~2MB figure). The
  fix set was applied + live-cache-verified except 2 small doc edits, but the run overflowed one window before
  gating -> reverted main to 67b7dbc, parked the fixes on branch `wip-codex-p1s5-review`, wrote
  `.agent/p1s5_review_handoff.md` so the next session lands it small (restore + 2 pending edits + gate + one
  `microglia (p1 s5 review): ... (codex)` commit). No main code changed; the S5 close stands gate-green.
  LANDED 2026-06-30 (next commit): restored the parked fixes + applied the 2 pending doc edits, gate green,
  branch + handoff deleted -- P1-S5 review CLOSED.
- 2026-06-30 P2 OPENED -> `.agent/p2_trajectory_plan.md`. User confirmed P2 (vs P3/P4 reorder). Research =
  2 parallel agents: Explore mined v1's EXECUTED trajectory analysis (Arc M -- "the one positive orthogonal
  interaction" = synergistic homeostatic->DAM acceleration; slingshot on harmony[1:15] + limma on 16-replicate
  pseudotime summaries via the SAME factorial machinery + Kitagawa composition/progression decomposition ->
  interaction loaded progression ~94% sig, null composition; Arc-O gene tradeSeq margin-neutral/contaminated/
  fragile-internals); general-purpose swept 2026 SOTA (condiments/tradeSeq pseudoreplicate + no factorial
  interaction interface -> descriptive only; per-replicate-summary + factorial = the replication-correct
  inferential route; glmmTMB beta-GLMM = on-lock per-cell sensitivity, TMB not Stan). CONVERGED. Decision gate ->
  LEAN ON-LOCK (rejected +destiny DPT, +off-lock Lamian). Full stack on-lock from the pinned snapshot (S2 pure-R
  primary; S3 source-compiled glmmTMB/TMB). Plan = 3 steps (later split to 4; see below). next = EXECUTE S1.
- 2026-06-30 P2-S1 DONE -> `R/trajectory.R` (build_activation_trajectory + pure helpers) + target
  `microglia_trajectory` + `tests/test_trajectory.R` + rproject.toml slingshot. slingshot on harmony[1:15], FORCED
  single Homeostatic->DAM lineage (2 substate super-clusters -> clean by construction), IFN/Prolif omitted
  (on_lineage flag + NA pt). LIVE: direction correct (mean pt Homeo 22.3 < DAM 36.2; Spearman pt-vs-DAM +0.56);
  score-axis concordance rho=0.62 (moderate-large, clears 0.5 gross-failure floor; slingshot=transcriptome vs
  score-axis=marker contrast -> related-not-identical); sensitivity robust (dims-10/20 + all-retained rho 0.99 vs
  primary); IFN omitted-fraction balanced across genotypes (low single-digit % -> conditioning barely skews interaction).
  Compact target ~3.4MB, 0 build warnings, gate green. Outcomes -> memory.md.
- 2026-06-30 P2-S2 SPLIT (no code shipped; a prior session ran out of window writing the original combined
  S2). Old S2 bundled the pure-R primary (weighted-limma 16-unit interaction + 3-channel Kitagawa) with the
  glmmTMB supportive arm (source-compiled new dep + ABI verify + fresh real-data build) + all tests -> two
  windows of work. Reverted that session, re-split into S2 (pure-R, NO new dep -> trajectory_progression) + S3
  (glmmTMB -> trajectory_glmm_sensitivity, INDEPENDENT target off microglia_trajectory$cell_frame), report -> S4.
  Plan steps now carry inline function contracts (the decided design) so the next session implements, not
  re-derives. Next open = S2.
- 2026-06-30 P2-S2 SPLIT AGAIN -> S2a + S2b (no code shipped; this session re-authored the full S2 [8 fns +
  fixture + tests] from the inline contracts but ran out of window before testing/gating -> reverted to 60774f9).
  LESSON: even the pure-R S2 alone (6 estimation fns + 2 inference fns + fixture + ~6 test groups + live smoke +
  gate) exceeds one window. Re-split at the VERIFICATION SEAM: S2a = per-replicate summary + contrast fit +
  Kitagawa decomposition (derive_batch / pseudotime_per_replicate / ordinary_t_table / fit_trajectory_contrasts /
  kitagawa_channels / decompose_progression_vs_composition + fixture + 5 unit-test groups, fully verifiable on the
  deterministic fixture, NO target); S2b = Freedman-Lane + run_trajectory_progression orchestrator +
  trajectory_progression target + FL tests + live smoke + gate (the live-integration half). Banked the expensive
  verification insights into each step spec (OLS-vs-manual needs a NON-additive term or sigma=0; NO
  full-orchestrator fixture test -- the exact-pure fixture has sigma=0 -> t=Inf). Next open = S2a.
- 2026-06-30 P2-S2a/S2b plan reviewed (PRE-implementation, no code): 11 findings, ALL accepted -> contracts
  HARDENED (within_<lc> col naming pinned S2a<->S2b; run_trajectory_progression meta now feeds
  assert_complete_crossing + factorial_design correctly; decompose weight-matrix dimnames set; fit_trajectory_contrasts
  top = named-list-by-contrast; pivot-free FL chol2inv(chol(crossprod)); per-endpoint FL weights specified; structural
  orchestrator test on a non-additive fixture jitter added). Core math re-verified CORRECT (Kitagawa identity,
  const-annihilation under weighted treatment coding, ordinary-t, FL null/signal, midpoint ramp). Next open = S2a.
- 2026-06-30 P2-S2a DONE -> `R/trajectory.R` +6 estimation fns (`derive_batch`, `pseudotime_per_replicate`,
  `ordinary_t_table`, `fit_trajectory_contrasts`, `kitagawa_channels`, `decompose_progression_vs_composition`) +
  `within_state_col` col-name sanitizer; `tests/helpers.R` `make_trajectory_cell_frame` (midpoint-ramp exact-pure
  fixture, default=pure-progression / flat-adv+dam_extra=pure-composition); `tests/test_trajectory.R` +5 groups.
  Kitagawa identity + intercept-free interaction-contrast reconstruction EXACT (<1e-8, ONE shared per-unit weight
  vector across the 4 channel-rows). 2 gate gotchas fixed -> memory.md (limma drops the 1-row coef/stdev.unscaled
  rowname; tapply over a single factor returns a non-conformable 1-D array). NO target/wiring (pure fns). Gate green.
  Next open = S2b.
- 2026-06-30 P2-S2b DONE -> `R/trajectory.R` +2 inference fns (`freedman_lane_interaction`,
  `run_trajectory_progression`) -> target `trajectory_progression` (`_targets.R`, reads COMPACT
  microglia_trajectory, NO 612MB load) + `tests/test_trajectory.R` +2 groups (source R/de_pb.R for
  assert_complete_crossing). Weighted/OLS/bounded limma interaction fits + 3-channel Kitagawa decompose +
  Freedman-Lane perm null; PRE-REGISTERED primary BH {progression_cf, within_homeostatic} (detail -> memory.md).
  HEADLINE (R4.6 re-baseline; numbers DRIFT-PRONE): the interaction RAISES mean_pt (p~0.04) BUT Kitagawa
  attributes it to COMPOSITION (comp_cf SIG fdr~0.025) NOT progression (prog_cf NEGATIVE/NS, perm_p~0.18; recon
  6.7e-15); primary family BOTH NS -> no statistically supported progression-beyond-composition signal (neg/NS
  estimate, NOT proven absence). DIVERGES v1 (~0.94 prog loading): R4.6 says the advance is MAINLY MORE
  DAM CELLS (confirms P1's sig DAM-fraction interaction), no detected further-advance -> defensible
  negative on the distinct progression question (acceptance
  outcome-independent). FL = SENSITIVITY not nominal-exact. Gate green (21 targets, built fresh on real data,
  tar_meta + render clean). CODEX-REVIEWED (8 findings all accepted, gate stays green): local 16u/9df crossing
  guard (shared assert_complete_crossing passes a complete sub-rectangle), min_within>=2 + within-sd guard,
  planned_primary/primary_within_skipped provenance, FL lm-oracle equality test (relative tol), + claim-honesty
  softening (shared-weight decomposition test; "no SUPPORTED progression signal" not "no synergy"). Next open = S3
  (glmmTMB per-cell sensitivity).
- 2026-07-01 P2-S3 DONE -> `R/trajectory.R` +`glmmtmb_pt_sensitivity` (+ helpers `.capture_quietly` /
  `.fit_pt_interaction`) -> target `trajectory_glmm_sensitivity` (`_targets.R`, reads COMPACT
  microglia_trajectory$cell_frame, INDEPENDENT of trajectory_progression) + `tests/test_trajectory.R` +3 groups +
  rproject.toml/rv.lock += glmmTMB. SUPPORTIVE per-cell beta GLMM pt01 ~ tau*amyloid + batch + (1|unit); degrade
  cascade singular/!pdHess/non-converge/non-finite/error -> rank-normal LMM -> RECORDED method="failed" (NA, never
  errors); Wald row by POSITION; .capture_quietly muffles+records warnings AND messages (sccomp lesson) -> 0 gate
  leak (detail -> memory.md). glmmTMB on CRAN = P3M trixie BINARY (plan's "source-compile" framing off), ABI-clean
  0-warning load (1.1.14/TMB 1.9.21/Matrix 1.7.5). LIVE (R4.6, DRIFT-PRONE): method=glmmTMB_beta on 22363 cells,
  tau:amyloid +0.123 logit (p~0.015, CI excludes 0), non-singular; 1.5s, tar_meta NA/NA. COMPOSITION-CONFLATED (the
  per-cell mean_pt analogue -> corroborates the position shift, NOT progression; S2b Kitagawa stays load-bearing).
  Full gate green (39 tests, 22 targets, render clean). Next open = S4 (report + integration) -> then CLOSE-OUT.
- 2026-07-01 P2-S4 SPLIT -> S4a + S4b (no code shipped on main; the done work PARKED on branch wip-p2s4-report,
  main CODE reverted to 863cb75 [restructure docs committed on top -> HEAD 1f8f3f2], orphan trajectory_report object pruned, gate green on clean main). The combined S4
  (trajectory_report_data extractor + trajectory_report target + _trajectory.qmd chapter + index wiring +
  _microglia.qmd pointer + test + render-debug + gate + docs + commit + close-out) overflowed one window. PARKED
  (built-validated pre-revert): the extractor + target (R/trajectory.R + _targets.R; the target built 674KB, all
  fields) + the FULL _trajectory.qmd draft (title "synergy adds DAM cells rather than advancing them", inline-
  computed prose, 5 sections / 3 figs / 2 tables / glmmTMB inline / provenance; UNRENDERED). Re-split at the
  data/report seam: S4a = compact extractor + trajectory_report target + a structure/guard test (pure-R, no
  render); S4b = restore the qmd + wire index/_microglia pointer + render-debug to 0-warning + docs. [SUPERSEDED by
  the 2nd-split entry below -- S4a restores its 3 files from wip-p2s4a-hardened; ONLY _trajectory.qmd comes from
  wip-p2s4-report.] DELETE both branches after S4b lands. Next open = S4a.
- 2026-07-01 P2-S4a guard-hardening PARKED (context-overflow remediation, 2nd split): S4a's NEW work -- harden the
  restored extractor's stopifnot to guard EVERY nested field the body reads (mirror microglia_report_data) + the
  structure/guard test -- required heavy re-derivation (enumerate the exact nested-field guards against the real
  cached targets; discover + fix a fixture composition-degeneracy [CONSTANT 12/12 per-unit composition -> zero
  comp_cf/cross residual variance -> NaN fdr -> VARY per-unit DAM drop k%%3, each >=10 retained]; validate a fresh
  674KB trajectory_report build) and overflowed one window before the full gate. Remediation (mirrors the P1-S5
  review + P2-S4 parks): reverted main to the gated d15e6c7 (clean; orphan trajectory_report object pruned; gate
  RE-VERIFIED green), parked the VALIDATED work (hardened extractor+guards + trajectory_report target + guard test;
  trajectory test passes warn=2, target built fresh on real data) on branch wip-p2s4a-hardened (bae0cf2), and
  REWROTE S4a into a MECHANICAL restore (git checkout the 3 files + gate + map.md + commit -- no re-derivation, no
  large-file reads, no fixture debugging). S4b's qmd stays on wip-p2s4-report (restore ONLY _trajectory.qmd there;
  its pre-hardening extractor is superseded). Next open = S4a (mechanical).
- 2026-07-01 P2-S4a review (parked code + handoff): review found guard-vs-use gaps -- `interaction$perm_p`
  (qmd table + inlined prose) + each canonical `weighted_top` mean_pt row (feeds the p_ctr geom_pointrange -> a
  dropped/NaN row warns under warn=2, red gate) were unguarded, a finite-fdr postcondition overclaim, and a stale
  restore-from-wip-p2s4-report sentence. Hardened the parked extractor (guard perm_p existence + inlined-row
  finiteness + each weighted_top mean_pt finite coef/CI; document finite-p/fdr as an INTENTIONAL build-fatal
  data-quality gate, not a proven universal) + 2 negative tests; re-validated (test warn=2 + fresh trajectory_report
  build 25ms, tar_meta NA); AMENDED the branch bae0cf2 -> 7139b09 (now +234: R/trajectory.R +157 / _targets.R +7 /
  test +70) + synced the S4a spec (hash/counts/explicit fresh-build cmd) + marked the stale sentence superseded.
  S4a stays MECHANICAL.
- 2026-07-01 P2-S4a DONE (mechanical restore, as planned -- no re-derivation): reconciled first (main d15e6c7..HEAD
  did NOT touch the 3 files -- only .agent docs), then `git checkout wip-p2s4a-hardened -- R/trajectory.R _targets.R
  tests/test_trajectory.R` -> diff HEAD = EXACTLY +234 (R/trajectory.R +157 / _targets.R +7 / test +70), the parked
  hardened extractor + `trajectory_report` target + guard test. Force-built the new leaf FRESH (tar_invalidate +
  tar_make trajectory_report, 29ms/338kB, tar_meta error+warnings NA -> guards pass on the REAL cached targets, no
  false-red) then full gate GREEN (23 targets now, up from 22; tests warn=2 incl. the +70-line guard test; force-
  render clean; render-log clean). map.md += trajectory_report_data fn + trajectory_report target + the S4a test
  line. Both wip branches LEFT for S4b to delete. Next open = S4b (restore _trajectory.qmd + wire index + rewrite
  the microglia pointer + render-debug to 0-warning).
- 2026-07-01 P2-S4b DONE -> P2 STEPS COMPLETE. Restored ONLY `_trajectory.qmd` from wip-p2s4-report (its pre-
  hardening R/trajectory.R + _targets.R superseded by S4a's), wired `{{< include _trajectory.qmd >}}` after
  _microglia.qmd + extended index Overview to name the trajectory chapter (activation-axis advance; synergy =
  composition not progression), rewrote the 2 _microglia.qmd forward-pointers (P2 built; synergy is compositional
  -- more DAM cells, no supported further-advance; cross-ref @sec-trajectory). Chapter (title "the tau-amyloid
  synergy adds DAM cells rather than advancing them", {#sec-trajectory}) reads the COMPACT trajectory_report,
  prose ALL inline-computed. RESOLVED the S4a review low:955 (trd$decomposition dead output): DROPPED it -- the qmd
  draws loadings from prov$*_loading (guarded finite) + per-channel coefs from the interaction comp_cf/
  progression_cf/cross rows, so a decomposition field only duplicated two live sources (figure-shaped dead output);
  updated R/trajectory.R extractor (guard-names + postconditions) + test to match. RENDER gotcha -> memory.md:
  ggplot2 4.0.3 deprecates scale_*_gradient(trans=) at 3.5.0 -> transform= (lifecycle warning reds warn=2);
  caught proactively, rendered 0-warning first try. Force-built trajectory_report fresh (44ms/338kB, tar_meta NA)
  then full gate GREEN (52 render chunks 0-warning, report 3.61MB/14.1s, tar_meta clean 23 targets, render-log
  clean; test_trajectory.R all groups warn=2). Docs: memory.md += P2-S4 built section; map.md += _trajectory.qmd
  include chain + @sec-trajectory pointer. Both wip branches (wip-p2s4a-hardened, wip-p2s4-report) DELETED. Next
  mode = CLOSE-OUT (adversarial plan review; fold P2 digest -> history.md; archive plan; reset Active plan; revise
  Cohesive-story finding #2 wording per the R4.6 re-baseline).
- 2026-07-01 P2-S4b review (of d5c4ccb): 8 findings (0 high / 3 med / 5 low), ALL accepted, gate stays green.
  HONESTY (med): intro's bare "the answer is no" overclaimed proven absence -> "does not detectably / no
  statistically supported evidence of further advance"; dropped "faster advance" rate wording (chapter disavows
  rate) across _trajectory.qmd + _microglia.qmd; index Overview += "mainly". MED2: the provenance cat printed
  "random-effect SD NA, singular = NA" on a FAILED glmm (only the prose sentence branched on method=="failed") ->
  added a glmm_prov branch that prints fail_reason instead. MED3: report-layer extractor now backs the qmd's
  "reconstruction residual" print + "three loadings sum to one" claim (recon_resid_max < 1e-6 + loadings-sum within
  1e-5; run_trajectory_progression still gates the exact < 1e-8, this catches a corrupted provenance copy). LOW:
  inlined the two data-derived hardcodes (16 units -> glmm$n_units, R 4.6 -> prov$r_version), left 9 resid df +
  sensitivity dims 10/20 as fixed DESIGN constants + scoped the memory/map "inline-computed" claim to match; added
  lineage_per_unit value/bounds postconditions (finite positive n_cells, on-lineage in [0,n_cells], known genotype)
  + folded v1_progression_loading/fdr into the finite provenance guard; dropped the un-surfaced/un-guarded
  "replicate-balanced re-anchoring is carried alongside" caveat clause (computed upstream, never surfaced here);
  test now asserts !("decomposition" %in% names(trd)) to lock the S4a-955 drop. Review CHECKED CLEAN: no qmd reads
  trd$decomposition, comp_cf/progression_cf/cross labels align across extractor/test/qmd, include placement +
  #sec-trajectory unique, no residual ggplot trans=. Rebuilt trajectory_report fresh on real data (new postconds
  pass, no false-red) + re-render 0-warning + tests warn=2. Commit `trajectory (p2 s4b review): ... (codex)`.
- 2026-07-02 infra Codex-only adaptation: retired tracked `CLAUDE.md`, `.claude/`, and `.serena/` project config;
  made `AGENTS.md` the sole canonical agent instruction surface; added `.codex/prompts/session.md`; rewired
  `.agent/context.sh` to Codex JSONL sessions; updated live memory/map/history wording from Claude deny-Read to
  Codex read-economy.
- 2026-07-02 infra Codex skill wrapper: added repo-scoped `.agents/skills/session-prompt` per Codex docs
  (`.agents/skills` repo discovery; explicit `$...` invocation). `$session-prompt` reads
  `.codex/prompts/session.md`; fresh `codex exec` smoke checks confirmed it appears in initial skill context.
- 2026-07-02 infra review wrapper retired: removed redundant review skills/prompts/script; session prompt now
  uses direct self-checks plus `.agent/context.sh` headroom tracking.
- 2026-07-02 P2 CLOSED: close-out review of plan body + shipped trajectory code/prose found no remaining blocker
  after the S4b review fixes; the only forward-state correction was the spine wording. Folded P2 digest ->
  history.md; archived plan -> `.agent/completed/p2_trajectory_plan_2026-07-02.md`; reset Active plan to none.
  Cohesive-story finding #2 now says mutant tau's amyloid-response modulation is DAM-cell COMPOSITION, not
  supported progression/acceleration. Full `scripts/check.sh` rerun green (tests warn=2, force-render 52 chunks,
  tar_meta clean across 23 targets, render-log clean). Next = PLAN P3 by default, with user confirmation.
- 2026-07-02 P3 OPENED -> `.agent/p3_mechanism_plan.md`. User confirmed default P3 (vs P4 reorder). Research =
  repo/map/history + v1 archive mining + 2026 web/API sweep. Default plan = lean mechanism rebuild: RNA pathway/TF
  + targeted NF-kB from existing pseudobulk targets, minimal 24M phosphosite DE solely for OmniPath kinase activity,
  compact mechanism report. Explicitly OUT vs v1: CCC/topology/SCENIC/human/ledger bloat; P4 keeps full
  cross-modality interpretation. Key S1 precondition = add/smoke-test `decoupleR`, `OmnipathR`, `fgsea`, `msigdbr`
  under the repo lock; prefer direct mouse `OmnipathR::enzyme_substrate(organism=10090)` over v1's off-lock
  `nichenetr` mapping, but fail loud on KSN coverage before analysis. CODEX-REVIEWED before commit (6 findings,
  all accepted): project-local prior cache + prior hashes/counts; explicit KSN gene-symbol validation; live S1
  KSN coverage against current phospho IDs; duplicate phosphosite collapse before decoupleR; run-order sensitivity
  for genotype-blocked 24M phospho samples; NF-kB attenuation primary = negative `interaction`, with
  `tau_in_nlgf` supportive only. Next = EXECUTE S1.
- 2026-07-02 P3-S1 DONE -> rproject.toml/rv.lock += decoupleR, OmnipathR, fgsea, msigdbr, digest;
  `R/mechanism.R` + `tests/test_mechanism.R` implement symbol-ranked matrices, decoupleR ULM wrapper, prior
  cache/fingerprint + expectation helpers, CollecTRI/KSN normalisers, phosphosite IDs, and KSN coverage. LIVE API smoke:
  `OmnipathR` loads warning-clean only after `TZ=UTC`; its CollecTRI + enzyme_substrate postprocessors fail on the
  current server schema (`ncbi_tax_id` missing), so loaders default to official OmniPath REST endpoints cached under
  `storage/cache/omnipath` (still direct mouse OmniPath, no nichenetr fallback). Hardened signs/components:
  CollecTRI uses consensus sign columns, drops 2,449 ambiguous rows; KSN drops 542 conflicting sign pairs.
  Observed REST priors pinned in code: CollecTRI 37,096 edges / 1,093 TFs / 6,010 targets hash
  027ee57a61246bff4127d9d36807469713731de552398bb81989a06fd1bc44e6; KSN 29,378 edges / 1,397 enzymes /
  13,048 sites hash 997b690d5efdfd8bb4424c12a29a80f5a980d8b3404025210e188281d554172d. Real phospho-site probe:
  64,328 rows -> 44,896 unique single-gene site IDs; KSN overlap 2,250 sites, 212 kinases pass minsize>=5, Gsk3b
  present + 245 matched sites. Next = S2.
- 2026-07-02 P3-S2 DONE -> RNA mechanism targets: `mechanism_collectri` (prior drift gate), `mechanism_gene_sets`
  (native mouse MSigDB GO BP/CC/MF + project sets DAM/Homeostatic/MHC_APC/IFN/NF-kB union + activated/repressed
  targets), `mechanism_tf` (decoupleR ULM CollecTRI on whole microglia + fit substates), `mechanism_pathway`
  (fgsea preranked GSEA), and `nfkb_attenuation` (predeclared primary-family gate). Live build clean in tar_meta;
  GO fgsea capped at maxSize=500 per fgsea manual runtime guidance, project sets uncapped; known sub-1e-10 fgsea
  p-value floor notices are recorded as provenance and any other fgsea warning fails loud; MSigDB/project gene-set
  payload hash pinned. CODEX-REVIEWED (7 findings all accepted): sign-aware NF-kB target-GSEA (activated NES,
  repressed -NES), conservative family-best BH for NF-kB TF/components, concordant-negative support gate, p-floor
  propagation, central contrast order. Live outcome: Myc is the top DAM interaction TF (negative, significant);
  NF-kB attenuation is DISCORDANT, not supported (target-GSEA negative, TF-family positive). Full `scripts/check.sh`
  green. Next = S3.

## Context ledger (per work-unit session)
Retro-recorded from session transcripts (this metric was meant to be logged per unit at the time, but
wasn't). Value = FINAL window occupancy = exactly what `.agent/context.sh` reports at session close (last
assistant turn: input + cache_creation + cache_read, over the window; format `N% used/window`); every row
is /200K (all recorded sessions ran on the 200K cap). One row per work SESSION (step + its review in one row;
commits chronological). Compacted sessions read BELOW their in-session peak (final != high-water; e.g. the
P1-S4 row is 45% final vs an 88% peak) -> read each as end-state, not total workload. Omitted: 9 no-commit
sessions with no unit recorded here (research / orient / failed-prompt / remote-push help, the temporary
milestone report-refresh, and this bookkeeping session).
```text
ctx (final)      date   work unit  (session step [+ review]; commits chronological)
---------------  -----  ---------------------------------------------------------------------
51% 103K/200K    06-29  git init + upstream-config merge; seed v1 -> archive (82130a5, archive branch)
95% 191K/200K    06-29  scaffold rebuild: fresh main + reset .agent/config (586f691)
30% 60K/200K     06-29  memory: storage/data symlink provenance (222a3ab)
46% 92K/200K     06-29  P0 plan + review hardening (f765085, f612ede)
20% 39K/200K     06-29  infra: review rollout-path fix (d0e870b)
59% 117K/200K    06-29  P0-S1 spine: env + targets/quarto scaffold +review (37216a8, 39bb692)
62% 123K/200K    06-29  P0-S2 data: io loaders + 4 modality targets +review (a8312b9, 66aa492)
88% 176K/200K    06-29  P0-S3 design: 5-contrast + pseudobulk/bulk DE +review (e6b183f, 466c65f)
68% 136K/200K    06-29  P0-S4 handoff: standalone-HTML report decision +review (b4f1496, a2c4db1)
42% 85K/200K     06-29  P0-S4 report checkpoint: offline-HTML pivot +review (2cbf238, 4951c63)
73% 147K/200K    06-29  P0-S4 report: IBM Plex theme + plot tests, S4 close +review (9365c82, d29e38b)
67% 134K/200K    06-29  P0-S5 gate scripts/check.sh + P0 close +review (998950d, 961dd22, b53f22f)
53% 106K/200K    06-29  P1 plan +review (0555e55, 0316bcb)
72% 145K/200K    06-29  P1-S1 microglia reprocess+cluster +review (6c69986, 1706fa9)
84% 168K/200K    06-30  P1-S2 substate annotate + prune +review (612dbd0, 5992e39)
88% 176K/200K    06-30  P1-S3 composition propeller + sccomp +review (7713248, cdff4bc)
76% 151K/200K    06-30  infra: gitignore .tokensave/ index (7a2518b)
79% 158K/200K    06-30  P1-S3 composition close: live sccomp arm +review (e27c41e, 7ca25f9)
54% 109K/200K    06-30  infra: upstream agent config sync (a8593ec)
45% 90K/200K     06-30  P1-S4 pseudobulk DE + stageR +review (64bf2c2, 94d4505)
55% 110K/200K    06-30  P1-S5 microglia report + P1 close + handoff +review (67b7dbc, 7a5e4d3, cf974bd)
67% 134K/200K    06-30  P2 plan +review (a6d58e2, 3fca9c2)
59% 118K/200K    06-30  P2-S1 slingshot pseudotime +review (d79219d, 3a1ff11)
83% 167K/200K    06-30  P2 plan: S2 -> S2/S3 split + memory trim +reviews (1d8d6c3, 7f96505, 20be54f, 3507b84, 60774f9)
86% 173K/200K    06-30  P2 plan: S2 -> S2a/S2b split +review (a62fc40, c54d854)
67% 134K/200K    06-30  P2-S2a estimation core +review (d92ace9, a99312b)
87% 174K/200K    06-30  P2-S2b interaction inference + FL null +review (3a8e2ef, 2c9e8b1)
81% 161K/200K    07-01  P2-S3 glmmTMB per-cell sensitivity +review (14fe8aa, 863cb75)
72% 145K/200K    07-01  P2-S4 split (S4a/S4b) +review (1f8f3f2, d15e6c7)
88% 177K/200K    07-01  P2-S4a handoff: park hardened extractor +review (0818d5d, 9111d69)
21% 42K/200K     07-01  P2-S4a extractor+target restore +review (b82b6f2, 15be78b)
71% 143K/200K    07-01  P2-S4b chapter + wiring (P2 steps complete) +review (d5c4ccb, cbbd129)
```
