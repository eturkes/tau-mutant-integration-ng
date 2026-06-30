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
   SYNERGISTIC acceleration along the activation trajectory (v1's one positive,
   significant orthogonal interaction).
3. Mechanism converges on Gsk3b / Myc signalling; tau attenuates amyloid-driven NF-kB.
4. Secondary axis: amyloid-linked synaptic suppression + TREM2 / APP clearance.
Streamlined OUT (v1 bloat): the 11-arc ledger + contest machinery, the margin-neutral
corroboration arcs (SCENIC, spatial-decon, celltype-specificity, gene-level dynamics),
the human-validation layer, the capstone convergence matrix, the heavy prose.

## Active plan: P2 interaction trajectory -> `.agent/p2_trajectory_plan.md` (OPENED 2026-06-30)
Tests the PIVOTAL claim P1 deferred: tau x amyloid synergy as PROGRESSION (extent-of-advance) along the
homeostatic->DAM activation trajectory. P1 ESTABLISHED a static COMPOSITIONAL synergy (propeller DAM-fraction
interaction sig; DE under-powered) -> P2 isolates the PROGRESSION channel BEYOND that composition shift. Stack
DECIDED = LEAN ON-LOCK (user gate): slingshot (harmony 15-dim) + UCell score-axis concordance anchor; weighted
per-replicate (16-unit) pseudotime-summary interaction through EXISTING factorial_design (no eBayes; 9 resid
df); 3-channel Kitagawa composition/progression/cross decomposition; glmmTMB per-cell sensitivity (supportive).
S2 primary = pure-R/no new dep; full P2 stack on-lock from the pinned snapshot (S3 = source-compiled glmmTMB/TMB;
no Stan/Python/GitHub). Converged: v1 Arc M (the executed analysis, found the one +ve orthogonal interaction)
+ 2026 SOTA sweep. Dropped v1 bloat: Python triangulation, CytoTRACE2, fragile Arc-O gene-dynamics.
5 steps (S1 trajectory+pseudotime -> S2a estimation-core + S2b interaction-inference [pure-R primary] -> S3
glmmTMB per-cell sensitivity [new-dep arm] -> S4 report). S1 + S2a DONE; next open = S2b. (S2 split out S3, then S2
itself split into S2a/S2b, on 2026-06-30 to fit one window each — see ledger.)

## Backlog - phased build (each phase = closeable increments; mine archive_digest per phase)
- P0 Foundations [DONE 2026-06-29]: project-local env (rv for R + uv .venv for Python), shared
  helpers (io / design+contrasts / plot theme), data load + QC sanity, 2x2
  factorial + 5 contrasts, concrete quality gate.
- P1 snRNAseq microglia core [DONE 2026-06-30]: reprocess (SCT) + substates (homeostatic / DAM /
  IFN / proliferative, UCell), composition (propeller primary + sccomp cross-check), pseudobulk DE across
  contrasts -> the robust amyloid->DAM activation programme (3-way confirmed: composition + DE + UCell score).
  Interaction = no large-effect DE, under-powered NOT absent (reported w/ MDE/CI + 123 stageR-confirmed, real |logFC| but sub-threshold per-contrast FDR) ->
  synergy = trajectory position/extent effect to P2 (rate interpretation only under the age-matched snapshot
  assumption). Digest -> history.md.
- P2 Interaction trajectory [OPENED 2026-06-30 -> p2_trajectory_plan.md]: activation pseudotime
  (homeostatic->DAM); test amyloid advance + tau x amyloid progression synergy, decomposed composition vs
  progression.
- P3 Mechanism: focused pathway/module survey; TF (decoupleR / CollecTRI) + kinase
  (decoupleR / OmniPath) -> Gsk3b / Myc; NF-kB attenuation check.
- P4 Cross-modality: GeoMx spatial DE (+ light deconvolution if it earns it),
  proteome + phospho DE, CCC for the synaptic/clearance axis, integrated divergence view.
- P5 Synthesis: ONE lean report - cohesive narrative + compact evidence table (no
  ledger / contest machinery).

## Ledger (trajectory)
- 2026-06-29 archived v1 -> branch `archive`; opened fresh orphan `main`; reset
  `.agent` docs + Claude config; reframed history as `archive_digest.md`; drafted
  this streamlined phase plan.
- 2026-06-29 S4 DONE (report engine closed): pivoted Quarto book -> ONE standalone offline HTML
  (`index.qmd` format:html embed-resources + `{{< include _qc.qmd >}}` + `theme.scss`); theme = crimson
  colours (#B0344D) + IBM Plex via 9 `@font-face` (relative `url("assets/fonts/<n>.woff2")` -> Quarto
  base64-inlines each woff2 under embed-resources). Render PROVED: 9 faces inlined OFFLINE (d09GMg magic),
  0 external loads, 0 error/0 warning (`tar_meta` all-NA), QC bounds pass (16x16 genotype-batch bijection).
  9 woff2 COMMITTED (assets/fonts/, deny-Read `**/*.woff2` + Serena ignored_paths), listed in
  `tar_quarto(extra_files=)`. Added `R/plot.R` + device-free `tests/test_plot.R`; de-staled map.md
  (book -> standalone report wiring) + S1 plan refs. Reasoning reversal (codex-reviewed): the earlier
  "IBM Plex can't work via theme.scss" was a MEASUREMENT artifact -- once `@font-face url()` is present
  the whole theme CSS URL-encodes, so raw greps for `IBM Plex`/colours read ~0; match encoded tokens.
- 2026-06-29 S5 DONE -> P0 CLOSED: concrete quality gate `scripts/check.sh` (fail-loud, zero-fault) =
  `rv sync` + `uv sync` | loop `tests/test_*.R` at warn=2 | `tar_make()` tee'd to a log | enforce
  `tar_meta(error,warnings)` all-NA SCOPED to `tar_manifest()$name` (drops functions/globals + stale dead
  rows) + render-log grep for Quarto/pandoc/knitr `warning`/`warn` (knitr's separate process bypasses
  tar_meta). Negative-tested (grep pattern + `stopifnot(FALSE)` -> exit 1); full gate green end-to-end.
  Locked memory.md gate (provisional -> concrete); finalised map.md wiring; folded P0 digest -> history.md;
  archived plan -> `.agent/completed/`. P0 foundation complete; next = P1.
- 2026-06-29 P0-S5 gate review-hardened (codex-review, post-close): the pre-review render-log grep was BLIND to
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
- 2026-06-30 P1-S5 codex review PARKED (context-overflow remediation): /codex-review of the S5 close raised 14
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
- 2026-06-30 P2-S2a/S2b plan codex-reviewed (PRE-implementation, no code): 11 findings, ALL accepted -> contracts
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
