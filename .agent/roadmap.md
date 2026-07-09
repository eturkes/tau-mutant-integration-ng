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
3. [TORN DOWN 2026-07-06] Mechanism (Myc-linked DAM interaction; NF-kB attenuation discordant/not supported;
   Gsk3b not recovered in 24M bulk phospho) -- chapter + targets + `R/mechanism.R` deleted from report + pipeline;
   science in git history + Ledger.
4. [TORN DOWN 2026-07-06] Cross-modality (amyloid-response spine + synaptic/clearance axis; focused Apoe-Trem2;
   SpatialDecon abundance blocked; full CCC not called) -- chapter + targets + `R/crossmodality.R` deleted; science
   in git history + Ledger.
REPORT SCOPE (current 2026-07-09): the rendered report = NINE visible figures with simple numbered
`Figure 1` ... `Figure 9` headings plus compact per-figure folded code only (no visible document title,
TOC, captions, body prose, tables, or global code menu):
microglia (P1) + trajectory (P2) +
one GeoMx modality-native figure (the former Figure 10 sample heatmap), one retained bulk modality-native
figure combining the proteome PCA and phosphoproteome heatmap, and two modality-context figures (four-method amyloid-response logFC
scatter; functional-group aggregate scores over the scatter's off-diagonal genes/proteins). The pipeline loads
snRNAseq plus lean GeoMx/proteome/phospho primary-DE targets solely for the modality figures; GeoMx carries only
the retained sample-heatmap descriptor as native report payload. The dedicated mechanism/cross-modality/qc/story
chapters, targets, R modules, tests, Python/uv surface, composition/sccomp target, per-subpopulation pseudobulk,
stageR layer, prose-inventory utility, and retired GeoMx exploratory/native panels remain deleted or retired from
the live DAG (Ledger 2026-07-08). 29 targets remain. P3/P4/P5 + the figure-expansion passes below stay as
historical DONE records (this roadmap holds the trajectory); their report chapters no longer exist.
Streamlined OUT (v1 bloat): the 11-arc ledger + contest machinery, the margin-neutral
corroboration arcs (SCENIC, spatial-decon, celltype-specificity, gene-level dynamics),
the human-validation layer, the capstone convergence matrix, the heavy prose.

## Active plan: none

No active plan. Next `$session-prompt` = PLAN: confirm the next backlog phase with
user before writing the plan.

## Backlog - phased build (each phase = closeable increments; mine archive_digest per phase)
- P0 Foundations [DONE 2026-06-29; live env leaned 2026-07-07]: project-local rv R env,
  project-local Quarto, shared helpers (io / design+contrasts / plot theme), data load + QC sanity,
  2x2 factorial + 5 contrasts, report-render quality gate. Historical Python/uv/test-loop gate was removed
  during the lean iteration cut.
- P1 snRNAseq microglia core [DONE 2026-06-30; live DAG leaned 2026-07-07]: reprocess (SCT) +
  subpopulations (homeostatic / DAM / IFN / proliferative, UCell), whole-microglia pseudobulk DE, and direct
  replicate-unit composition bars in the compact report bundle. Historical composition/sccomp/stageR/per-subpopulation
  DE machinery is retired from the live path. Synergy is carried into P2 as DAM-cell composition rather than
  progression/acceleration. Digest -> history.md.
- P2 Interaction trajectory [DONE 2026-07-02 -> `.agent/completed/p2_trajectory_plan_2026-07-02.md`]:
  activation pseudotime (homeostatic->DAM); amyloid advances the activation axis, but the tau x amyloid
  interaction decomposes to more DAM cells, not supported progression beyond composition.
- P3 Mechanism [DONE 2026-07-02 -> `.agent/completed/p3_mechanism_plan_2026-07-02.md`]:
  focused RNA pathway / TF + targeted NF-kB + minimal 24M bulk-phosphosite kinase. Myc supported;
  NF-kB attenuation discordant / not supported; Gsk3b not recovered.
- P4 Cross-modality [DONE 2026-07-02 -> `.agent/completed/p4_cross_modality_plan_2026-07-02.md`]:
  GeoMx spatial DE, 24M proteome + raw/corrected phospho, targeted clearance-axis CCC-lite, integrated
  divergence view. Strengthens amyloid-response / synaptic-clearance axes; interaction stays mostly
  microglia-composition-specific; SpatialDecon/full CCC not earned; bulk run-index sensitivity load-bearing.
- P5 Synthesis [DONE 2026-07-02 -> `.agent/completed/p5_synthesis_plan_2026-07-02.md`]:
  ONE lean report - cohesive narrative + compact evidence table (no ledger /
  contest machinery).
- Figure expansion [DONE 2026-07-02 -> `.agent/completed/figure_expansion_plan_2026-07-02.md`]:
  post-report visual-density pass. Inline chapter expansion backed by compact
  per-chapter figure targets; +26 planned figures landed without changing claims.
- Spatial decon follow-up [DONE 2026-07-02 -> `.agent/completed/spatial_decon_followup_plan_2026-07-02.md`]:
  gated GeoMx tissue-abundance follow-up to the P4 "SpatialDecon not earned"
  row. S0 chose the default: install/run SpatialDecon only after a compact
  full-reference profile earns it, keep nuclei absolute counts disabled, report
  broad abundance primary and microglia subpopulations only if stable. Final status:
  profile earned, SpatialDecon ran, abundance DE blocked by unresolved AOIs;
  residual audit is report QC; synthesis says blocked, not skipped.
- Prose-to-figures reduction [DONE 2026-07-03 -> `.agent/completed/prose_to_figures_plan_2026-07-03.md`]:
  user feedback: final output remains too prose-heavy even after Figure
  expansion. Plan opens a visual-first reduction pass: prose inventory,
  replacement manifest, compact visual data slots/schematics, chapter conversion,
  and before/after prose-density QA. S0 chose aggressive inline visual conversion;
  S1 measured the baseline and replacement manifest; S2 built visual grammar /
  compact data slots; S3 converted overview+synthesis; S4 converted result
  chapters; S5 closed with final report prose 5,111 -> 1,164 words (77%
  reduction), rendered HTML QA green, and full gate green.
- Figure-caption-only report [DONE 2026-07-03 -> `.agent/completed/figure_caption_only_plan_2026-07-03.md`]:
  user requested the extreme endpoint: no prose, just figures and captions.
  Plan target = rendered main path with headings + figures + captions only;
  paragraphs/tables/provenance text removed or figure-encoded; `fig-alt`
  retained for accessibility. S1-S5 closed: rendered main path has 48 figures /
  48 captions / 48 alt attributes, no visible body prose/tables/provenance/code
  UI, strict gate green, claim-parity review no accepted blockers, and full gate
  green.
- Box-figure curation [DONE 2026-07-03]:
  user feedback: too many figures were just boxes. Removed pure status/logic/
  checklist/matrix figures from QC, result, mechanism, and cross-modality
  chapters; retained data-rich journal-relevant plots. Rendered main path now has
  31 figures / 31 captions, no removed-label hits, strict caption-only HTML QA
  green, and full gate green.
- Figure story layout [DONE 2026-07-03]:
  user feedback: add/shape figures and make the layout tell more story while
  staying conventional for the field. Added a data-backed 2x2 design/sample
  support figure, reshaped the DAM composition panel around a direct tau x
  amyloid DAM-response plot, tagged composite panels, added a trajectory
  DAM-fraction trend, and replaced NF-kB status-box logic with a primary-score
  lollipop. Rendered main path now has 33 figures / 33 captions, strict
  caption-only HTML QA green, and full gate green.
- Visual maturity pass [DONE 2026-07-03]:
  user feedback: current colour scheme and figure types felt juvenile. Replaced
  crimson/bright plot defaults with deep-ink chrome, muted
  graphite/teal/ochre/wine categorical accents, neutral count-density fills,
  and teal/paper/wine signed fills. Reworked oversized design circles into a
  compact 2x2 tile matrix, replaced the Myc/NF-kB lollipop with a TF focus
  heatmap, made the NF-kB primary gate signed bars, removed batch text labels
  from the trajectory scatter, and muted binary/direction count panels. Rendered
  main path remains 33 figures / 33 captions / 33 nonblank alts, duplicate IDs
  0; full gate green.
- Color saturation pass [DONE 2026-07-03]:
  user feedback: the matured palette was too dull. Raised chroma without
  returning to toy-like defaults: deep-blue/teal/slate HTML chrome,
  steel-blue/teal/amber/cranberry genotype accents, blue/cranberry/amber/violet
  subpopulation accents, stronger blue/paper/cranberry signed fills, and a richer
  blue-teal sequential count gradient. Figure grammar from the visual maturity
  pass is unchanged. Rendered main path remains 33 figures / 33 captions / 33
  nonblank alts, duplicate IDs 0; full gate green.
- Four-modality integration figures [DONE 2026-07-03]:
  user requested more figures integrating the four modalities. Added a compact
  evidence-table reduction to `crossmodality_figures` and three visible
  cross-modality panels: FDR support by assay family, pathway-axis support
  across assay families, and selected axis-symbol modality tiles. Rendered main
  path now has 36 figures / 36 captions / 36 nonblank alts, duplicate IDs 0;
  full gate green.
- Figure elegance pass [DONE 2026-07-03]:
  user feedback: too many figures still looked blocky (bar plots, squares).
  Replaced visible QMD `geom_col`/`geom_tile`/`geom_rect`/`geom_bin2d` grammar with
  node/stem design, distribution traces, circular density dots, bubble/dot
  matrices, contours, and point-stems; updated captions/alts accordingly.
  Added theme print-overflow override after Chromium PDF QA exposed figure
  scrollbar chrome. Rendered main path remains 36 figures / 36 captions / 36
  role-img elements; full gate green; Chromium PDF contact sheet clean.
- Story plate synthesis [DONE 2026-07-03]:
  user asked to comb through the results and tease out a coherent scientific
  story, then create publication-grade supporting figures. Added `story_figures`
  and `_story.qmd` before the result chapters: a core-evidence plate (replicate
  DAM response, signed whole-MG DE counts, interaction decomposition) and a
  mechanism/integration plate (pathway-axis modality support, Myc/NF-kB/Gsk3b
  triage, measured clearance-pair support). Strict caption-only gate counts the
  new chapter; Chromium PDF story-page QA clean; full gate green.
- Field-convention figure pass [DONE 2026-07-03]:
  user feedback: some figures did not look conventional for the field. Web/literature
  check supported point volcanoes/scatters, stacked composition bars, and tile
  heatmaps over the prior circular-density/bubble grammar. Revised visible
  design/composition/matrix/volcano/scatter panels accordingly, while retaining
  compact target-backed figure data.
- Cross-modality narrative figure pass [DONE 2026-07-03]:
  user feedback: integrated cross-modality figures are useful but dashboard-like.
  Replaced generic assay/contrast dashboards with journal-style evidence plates
  plus a closing model while preserving the earned claim boundaries. Final report
  has 37 figures / 37 captions / 37 nonblank alts; strict caption-only HTML QA,
  Chromium PDF content-page QA, and full `scripts/check.sh` green. Residual:
  Chromium emits a trailing blank PDF page, but content pages have no clipped
  legends or unreadable cross-modality labels.
- Plot variety pass [DONE 2026-07-03]:
  user feedback: too many repeated stem-style charts. Replaced repeated visible stem
  grammars with a mixed set of diverging count bars, point-only effect plates,
  dot-matrix support counts, signed heatmaps, bubble audit, and signed loading
  bars. Remaining visible `geom_segment` use is the closing-model arrows. Full
  `scripts/check.sh` green; Chromium PDF QA checked changed pages.
- Figure elaboration pass [DONE 2026-07-03]:
  user feedback: many figures remained overly simplistic. Enriched the retained
  conventional report figures without changing claims: QC design now includes
  replicate-unit support; story mechanism/cross-modality plate shows a signed
  evidence ladder plus all clearance-pair rows; subpopulation and trajectory audits
  expose replicate-unit distributions; trajectory decomposition uses a signed
  reconstruction path; NF-kB gate shows primary/supportive signed rows; GeoMx
  sensitivity shows supported-row totals plus lost/gained/flip changes. Targeted
  render green; Chromium PDF spot checks clean.
- Conventional figure cleanup [DONE 2026-07-03]:
  user feedback: some figures still looked strange relative to field norms. Kept
  claims/data fixed and replaced the remaining custom grammars with standard
  biology-paper forms: score distributions -> violin+box plots, pruning audit ->
  faceted bars, trajectory reconstruction path -> signed contribution bars,
  NF-kB stem audit -> score heatmap, GeoMx sensitivity bubble/segment hybrid ->
  count heatmaps, and story mechanism stems -> point evidence. Targeted render
  green; Chromium PDF spot checks clean.
- GeoMx exploratory figures [DONE 2026-07-08 -> `.agent/completed/geomx_exploratory_figures_plan_2026-07-08.md`]:
  user requested web-searched exploratory GeoMx figure inventory + roadmap sessions. Added
  nine compact GeoMx exploratory figures in separate sessions: AOI QC atlas, normalization/RLE,
  ordination, gene detection, sample heatmap, spatial program overlays, contrast diagnostics,
  ROI/segment replicate audit, and decon feasibility/status. All figures are descriptive,
  exclude no AOIs, change no DE model, and keep SpatialDecon abundance blocked/not claimed.

## Ledger (trajectory)
- 2026-07-09 Figure 7 vertical stack DONE (ad hoc user review): changed the bulk context plate
  from horizontal PCA/heatmap panels to a centered vertical stack. The figure device is now
  15.6x19.6 with row spacers, preserving the old child plot allocations: PCA 6.084x9.8,
  phosphosite heatmap 9.516x9.8.
- 2026-07-09 Figure 2 vertical stack DONE (ad hoc user review): changed the two microglia UMAP
  subplots from a horizontal patchwork to a vertical A/B stack and resized the figure device from
  12.4x4.4 to 6.2x8.8 so each child plot keeps its prior subplot footprint.
- 2026-07-09 Figure 7/8 bulk context merge DONE (ad hoc user task): combined the former
  Figure 7 proteome PCA and former Figure 8 phosphosite heatmap into one Figure 7 bulk context
  plate, removed both proteome and phosphoproteome volcano plots from the live report payload,
  and renumbered the four-method scatter / functional-group panels to Figures 8/9.
- 2026-07-09 Figure 10 right-edge clipping fix DONE (ad hoc user review): padded the
  functional-group x-scale so the maximum-score endpoint no longer sits on the panel boundary.
- 2026-07-09 Figure 10 point-size expansion DONE (ad hoc user review): kept the
  existing minimum functional-group point size and raised only the maximum size so
  larger scored-item buckets read more prominently.
- 2026-07-09 Figure 10 accepted label polish DONE (ad hoc user review): visible buckets now use
  `Cell-Cell Adhesion` for the leukocyte adhesion GO family and `Extracellular Matrix` for the
  matrix/adhesion GO family.
- 2026-07-09 Figure 10 label revision DONE (ad hoc user review): replaced disliked compact labels
  with literal shorter buckets (`Leukocyte Adhesion`, `Matrix`) and split `Chemotaxis / Phagocytosis`
  into `Phagocytosis` for `Camk1d` and `Chemotaxis` for `Mpp1`.
- 2026-07-09 Figure 10 label compaction DONE (ad hoc user review): shortened long visible
  GO term-family labels (`Complement / MHC`, `Leukocyte Trafficking`, `ECM / Adhesion`) and
  capped each retained feature list at two lines by increasing items per line when needed.
- 2026-07-09 Figure 10 role-label hardening DONE (ad hoc user review): split
  phagocytosis/chemotaxis out of the complement/antigen bucket and replaced the vague cell-adhesion
  row with clearer GO term-family buckets (`Leukocyte Adhesion / Migration`, `ECM / Substrate
  Adhesion`, `Cell Motility / Cytoskeleton`). Current visible examples: `Camk1d` =
  chemotaxis/phagocytosis, `Cr1l` = complement/antigen, `Icam2` = leukocyte adhesion/migration.
- 2026-07-09 Figure 10 category cleanup DONE (ad hoc user task): replaced the broad-first
  `Microglial Activation` role bucket with priority-ordered GO-BP role buckets. Narrow process
  buckets now claim selected off-diagonal features before a broad immune/inflammatory residual
  bucket, so Figure 10 no longer presents assay-specific opposite-sign feature sets under the
  same microglial-activation label.
- 2026-07-09 fast-iteration cut: `scripts/check.sh` now invalidates/builds only `report`
  (no env sync, no all-target metadata/log scan), and `render_report()` resets/prunes `report/`
  so the sole user-facing output is `report/tau-mutant-integration.html`.
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
  defaults: UCell scoring, Harmony batch-only (sex aliased), 4 argmax subpopulations + aux MHC/APC, voomQW +
  robust eBayes + stageR. Key reframing from mining: the tau x amyloid interaction is STATIC-NULL (v1 matched-
  power null even single-cell) -> the synergy is a position/extent effect (rate only under the age-matched snapshot assumption), deferred to P2 (trajectory); P1
  nails the robust amyloid->DAM headline + subpopulations + composition. Carry Thrupp 2020 (snRNA under-detects
  ~18% DAM activation genes). 5 steps S1-S5; next = S1 reprocess+cluster.
- 2026-06-29 P1-S1 DONE -> `microglia_processed` target + `R/microglia.R` (reprocess_microglia +
  marker_mean_by_cluster) + `tests/test_microglia.R`. SCT-v2/glmGamPoi -> Harmony(batch-only) -> Louvain
  multi-res {0.2,0.4,0.6}, primary=0.4 (12 clusters) -> UMAP, on the live 26k subset (138s, 687MB qs). Three
  pkg-drift fixes vs v1 recipe: harmony 2.0 dropped `assay.use` (assay implicit in reduction.use); SCTransform
  trips future's 500MiB globals cap (raise it); RunUMAP's "default changed" NOTICE is a WARNING that would fail
  the gate via tar_meta -> silenced by Seurat's own `Seurat.warn.umap.uwot=FALSE` (other warnings still
  surface). Also strips stale upstream meta shadows (pca1/umap1 coords-as-cols, SCT_snn_res.0.01). Acceptance:
  reductions {pca,harmony,umap}+cluster factor present; post-Harmony marker separation confirmed (distinct
  argmax per subpopulation); re-run STABLE (observed ARI=1.0 under recorded threads; clusters derive from the
  seeded PCA->Harmony->SNN->Louvain chain, UMAP is viz-only; NOT bitwise-guaranteed per up-to-tolerance); gate green.
- 2026-06-30 P1-S2..S4 DONE (annotate / composition / pseudobulk-DE; outcomes folded to history.md + memory.md P1
  sections, full detail in the archived plan) -> P1-S5 DONE -> P1 CLOSED. S5 = `_microglia.qmd` microglia chapter
  (subpopulation UMAP, composition forest/table, amyloid->DAM volcano + DE counts, under-powered-interaction + P2
  pointer, Thrupp + dropout caveats) reading a COMPACT `microglia_report` target (microglia_report_data extracts
  umap+subpopulation+z+prune/provenance from the 612MB annotated -> the gate's force-render stays cheap ~12s). Prose
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
  single Homeostatic->DAM lineage (2 subpopulation super-clusters -> clean by construction), IFN/Prolif omitted
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
  targets), `mechanism_tf` (decoupleR ULM CollecTRI on whole microglia + fit subpopulations), `mechanism_pathway`
  (fgsea preranked GSEA), and `nfkb_attenuation` (predeclared primary-family gate). Live build clean in tar_meta;
  GO fgsea capped at maxSize=500 per fgsea manual runtime guidance, project sets uncapped; known sub-1e-10 fgsea
  p-value floor notices are recorded as provenance and any other fgsea warning fails loud; MSigDB/project gene-set
  payload hash pinned. CODEX-REVIEWED (7 findings all accepted): sign-aware NF-kB target-GSEA (activated NES,
  repressed -NES), conservative family-best BH for NF-kB TF/components, concordant-negative support gate, p-floor
  propagation, central contrast order. Live outcome: Myc is the top DAM interaction TF (negative, significant);
  NF-kB attenuation is DISCORDANT, not supported (target-GSEA negative, TF-family positive). Full `scripts/check.sh`
  green. Next = S3.
- 2026-07-02 P3-S3 DONE -> minimal 24M bulk-phosphosite kinase layer: `phospho_de_24m` (16/16 sample-key matched
  runs, 4/genotype; log2 positive intensities with 1,213 nonpositive -> NA; median-normalised; prevalence >=2 samples
  in all 4 genotypes; limma-trend no batch + additive run-index sensitivity), `kinase_activity` (direct-mouse KSN
  drift-gated + coverage-gated; decoupleR ULM primary + run-index), `kinase_mechanism_summary` (significant kinases +
  explicit Gsk3b rows every contrast). Live target clean: 64,328 raw rows -> 17,707 filtered DE rows -> 12,938
  filtered unique single-gene site IDs; KSN coverage 1,164 matched sites / 123 kinases minsize>=5 / Gsk3b passes with
  169 sites. Live qualitative read: Gsk3b is covered but NOT significant for interaction or tau_in_nlgf; several
  tau_in_nlgf primary kinases (CAMK-family etc.) pass FDR<0.10 but weaken under run-index adjustment -> S4 must carry
  the bulk/not-microglia + genotype-blocked run-order caveats. Full `scripts/check.sh` green (31 current
  targets/branches; force-render report clean). Next = S4 report + integration.
- 2026-07-02 P3-S4 DONE -> mechanism report + integration: `mechanism_report_data` bundles compact pathway/TF/NF-kB/
  kinase highlights plus P1/P2 anchors into `mechanism_report` (~26KB; no heavy Seurat), `_mechanism.qmd` included after
  trajectory. Live report read: Myc is the strongest rebuilt interaction TF signal; NF-kB attenuation gate is discordant
  and unsupported; Gsk3b is covered but not significant in interaction/tau_in_nlgf; tau_in_nlgf kinase hits weaken under
  run-index sensitivity, so the bulk hippocampus/not-microglia + genotype-blocked run-order caveats are explicit. Full
  `scripts/check.sh` green (64 render chunks, report ~3.82MB, tar_meta/render-log clean). Self-review accepted/fixed
  2 lows: Myc-specific wording (not Myc-family) + zero-row TF-top selection edge case. Next = CLOSE-OUT.
- 2026-07-02 P3 CLOSED: close-out review of plan body + shipped mechanism code/prose found no blocker. One low
  report-prose robustness issue was accepted/fixed: the NF-kB paragraph now branches on the actual gate status instead
  of hardcoding the current discordant sign pattern. Folded P3 digest -> history.md; archived plan ->
  `.agent/completed/p3_mechanism_plan_2026-07-02.md`; reset Active plan to none. Cohesive-story finding #3 now states
  the rebuilt mechanism asymmetry (Myc supported; NF-kB/Gsk3b not recovered), not the v1 mechanism headline.
- 2026-07-02 P4 OPENED -> `.agent/p4_cross_modality_plan.md`. User confirmed default next phase. Research =
  repo/map/history + v1 archive mining + current-method sweep. Default plan = lean on-lock cross-modality rebuild:
  replicate-aware GeoMx DE with deconvolution gated, 24M proteome + raw/corrected phospho, targeted clearance-axis
  CCC-lite, and an integrated gene/pathway divergence view. Explicit alternatives: proteomics-heavy
  QFeatures/msqrob2PTM, full off-lock CCC stack, or no deconvolution. Next = EXECUTE S1.
- 2026-07-02 P4-S1 DONE -> `R/crossmodality.R` + `tests/test_crossmodality.R` + target `geomx_de`. GeoMx DE =
  explicit RNA count layer (object default is SCT), edgeR TMM + limma-voom, slide fixed effect,
  duplicateCorrelation block `genotype:bio_rep`, plus unblocked and bio-unit-collapsed sensitivities. Live build
  warning-clean: 91 AOIs / 15 bio-units / 19,959 genes kept; duplicateCorrelation consensus small positive; both
  sensitivities fit. Counts layer has 351 non-integer residues (max ~0.5) -> recorded, not rounded. Decon preflight
  = `defer`: SpatialDecon is pinned-repo available, Q3/background usable, but 42 nuclei sentinels disable absolute
  nuclei rescaling and no compact reference profile is built in S1. Full `scripts/check.sh` green (33 current
  targets/branches). Next = S2.
- 2026-07-02 P4-S2 DONE -> `proteome_de_24m`, `phospho_corrected_24m`, `bulk_omics_summary`.
  Proteome = exact 16-run sample-key match, raw positive row-sum to `PG.ProteinGroups`, log2, median normalise,
  prevalence filter, limma-trend + additive run-index sensitivity. Corrected phospho = P3 raw phosphosite 24M layer
  minus matched parent-protein log2 intensity, then re-filter/refit; P3 raw `phospho_de_24m` reused. Live build
  warning-clean/tar_meta clean: proteome 3,379 protein groups; corrected phospho 15,477 sites; 15,647 phosphosite rows
  had matched filtered parents and 2,059 lacked a filtered parent. Run-index sensitivity is harsh for most primary bulk
  hits, so downstream support must downgrade run-order-dependent signals. Next = S3.
- 2026-07-02 P4-S3 DONE -> `clearance_axis` + gated spatial-composition helpers. SpatialDecon remains intentionally
  skipped: live GeoMx preflight status is still `defer` (42 nuclei sentinels disable absolute rescaling; no compact
  profile was built), and `clearance_axis_data()` now fails loud if that preflight ever becomes `earned` before
  `geomx_decon` / `geomx_abundance_de` targets exist. Added Q3-scaled background, profile-collinearity, and
  log-beta abundance-DE helpers + tests. Live clearance-axis target warning-clean/tar_meta clean: 1,400 measured
  anchor rows, 166 synaptic GO-set rows, all 15 dictionary anchors measured somewhere; CCC-lite verdict is earned
  only for `Apoe_Trem2` in `nlgf_in_p301s` via coherent supported GeoMx + snRNAseq microglia evidence. No full CCC
  method is called. Next = S4.
- 2026-07-02 P4-S4 DONE -> `crossmodality_table`, `crossmodality_pathway`, `crossmodality_divergence`.
  Harmonised snRNAseq/GeoMx/proteome/raw+corrected phospho/TF/kinase rows to one symbol-evidence table, preserving
  feature/site collapse provenance + missingness. Pathway target scores selected project + RNA-mechanism GO sets from
  ranked modality statistics; divergence target focuses the four story contrasts and keeps mixed signs explicit.
  Count-honesty hardening: `modality_class` drives `n_modalities_present/sig`; `modality_group` is layer-level evidence
  only. Live fresh targets warning-clean/tar_meta clean (~10MB table, ~108KB pathway, ~1.9MB divergence). Full
  `scripts/check.sh` green. Next = S5 report chapter.
- 2026-07-02 P4-S5 DONE -> P4 STEPS COMPLETE. Added compact `crossmodality_report` (~23KB qs live) via
  `crossmodality_report_data`, `_crossmodality.qmd` include after mechanism, index overview wording, and S5 tests.
  Chapter reads only the compact target and covers GeoMx spatial DE, 24M bulk proteome/raw+corrected phospho with
  run-index caveats, then-skipped decon + clearance-axis CCC-lite, integrated pathway/symbol divergence, and P4 synthesis.
  Live read: GeoMx/bulk strongest in amyloid contrasts; interaction is smaller outside microglia composition/trajectory;
  SpatialDecon was skipped/defer at P4 close (superseded by follow-up S4's blocked-fit state); CCC-lite earns only
  Apoe_Trem2 in `nlgf_in_p301s`; no full CCC. Full
  `scripts/check.sh` green (tests warn=2, forced 82-chunk report render, tar_meta clean across 41 targets/branches,
  render-log clean). Next mode = CLOSE-OUT.
- 2026-07-02 P4 CLOSED: close-out review of plan body + shipped cross-modality code/prose found one low report-prose
  robustness issue; accepted/fixed by deriving earned clearance-pair prose from `crossmodality_report` instead of
  hardcoding the current Apoe-Trem2 row. Folded P4 digest -> history.md; archived plan ->
  `.agent/completed/p4_cross_modality_plan_2026-07-02.md`; reset Active plan to none. Cohesive-story finding #4 now
  states cross-modality corroboration + measured Apoe-Trem2 support, with SpatialDecon/full CCC not earned. Final
  close-out gate green. Then = PLAN P5 Synthesis, now closed.
- 2026-07-02 P5 OPENED -> `.agent/p5_synthesis_plan.md`. User confirmed default final phase. Research =
  current report/DAG wiring + history/map + v1 capstone archive + current Quarto/targets/MOFA2 docs. Default plan =
  compact read-only synthesis target, upfront `_synthesis.qmd`, and final lean report pass. Explicitly OUT:
  v1 convergence/ledger/contest machinery, MOFA/meta-analysis, human/SCENIC/topology/full CCC/decon side arcs, and
  new biological inference. Next = EXECUTE S1.
- 2026-07-02 P5-S1 DONE -> `R/synthesis.R` + `tests/test_synthesis.R` + target `synthesis_report`.
  `synthesis_report_data` reads ONLY compact report bundles (`microglia_report`, `trajectory_report`,
  `mechanism_report`, `crossmodality_report`; no `crossmodality_divergence`) and returns a ~4.8KB qs object:
  headline, 10-row descriptive evidence table, status counts, unsupported/open rows, and tiny source highlights.
  Guarded anchors: amyloid->DAM support, DAM composition interaction, trajectory comp/prog rows, Myc, NF-kB gate,
  Gsk3b rows, GeoMx/bulk run-index caveats, SpatialDecon status, and earned-pair empty/not-empty handling. No
  support/contradict/net-score/ledger columns. Synthetic tests + fresh target build green; manifest/raw deps show only
  the four compact report targets + the synthesis function. Next = S2 synthesis chapter + report wiring.
- 2026-07-02 P5-S2 DONE -> `_synthesis.qmd` + `index.qmd` report wiring. The synthesis chapter is included
  immediately after Overview and tar_loads ONLY `synthesis_report`; it renders an answer-first paragraph, status-count
  bar plot, compact evidence table, and unsupported/unearned paragraph from the compact target. Overview now states
  the final answer up front and treats P1-P4 chapters as audit trail (no "final synthesis still open" wording).
  Docs: map.md include chain + memory.md render contract. Next = S3 lean report pass.
- 2026-07-02 P5-S3 DONE -> P5 STEPS COMPLETE. Lean report pass removed stale P3/P4 forward pointers, cleaned
  required report-source stale search (`P5` / `final synthesis` / `still open` / `before the final synthesis`),
  tightened trajectory/mechanism/synthesis wording around progression-vs-rate, and moved raw caveat cleanup into
  `synthesis_report_data` so the compact evidence table itself no longer carries `value(s)` or
  `deconvolution deferred to S3`. Adversarial review accepted/fixed the raw-caveat leak + acceleration-wording
  caveat. Full gate GREEN after final fixes (tests warn=2, forced 88-chunk report render, tar_meta clean across
  42 current targets/branches, render-log clean). Next mode = CLOSE-OUT.
- 2026-07-02 P5 CLOSED: close-out review of plan body + shipped synthesis/report code/prose found no remaining
  blocker. Folded P5 digest -> history.md; archived plan ->
  `.agent/completed/p5_synthesis_plan_2026-07-02.md`; reset Active plan to none. Cohesive-story spine unchanged by
  P5; the upfront synthesis chapter now states it as the final compact answer. Final close-out gate green.
- 2026-07-02 Figure expansion OPENED -> `.agent/figure_expansion_plan.md`. Research =
  current qmd/target inventory, compact report-target capacity, v1 figure-form mining, and current Quarto/targets/
  ggplot2/patchwork docs. Default = add a compact `figure_atlas` target + `_figures.qmd` atlas chapter with ~20-26
  additional figures; keep existing P1-P5 claims closed and synthesis answer-first. Alternatives recorded for inline
  chapter expansion, prebuilt publication gallery, or schematic-first expansion. Next = S0 user-choice gate before
  implementation.
- 2026-07-02 Figure expansion S0 DONE: user chose inline chapter expansion. Revised
  plan away from a standalone atlas chapter toward compact per-chapter figure targets
  plus inline qmd additions; figure budget stays +20-26. Next = S1 data contract.
- 2026-07-02 Figure expansion S1 DONE -> `R/figures.R` + four compact targets:
  `microglia_figures`, `trajectory_figures`, `mechanism_figures`,
  `crossmodality_figures`; `tests/test_figures.R` pins the 26-figure manifest and
  finite/slot guards. Fresh target build warning-clean/tar_meta clean; object sizes
  all <5MB (microglia 2.381MB, trajectory 0.033MB, mechanism 0.107MB,
  crossmodality 0.514MB). Heavy plot shapes are pre-binned/top-row compressed
  (volcanoes, GeoMx, raw-vs-corrected phospho) so qmds can stay compact. Next = S2
  synthesis + microglia inline figures.
- 2026-07-02 Figure expansion S2 DONE -> `_synthesis.qmd` + `_microglia.qmd`
  inline visual-density pass. Added 9 labelled `fig-*` chunks: synthesis
  claim-source evidence map plus 8 microglia figures (genotype-faceted UMAP,
  score triptych/distributions, 16-unit composition, composition concordance,
  all-contrast whole-MG volcanoes, subpopulation fit audit, within-subpopulation DE counts).
  Forced report render warning-clean (106 chunks). Claims unchanged: robust
  amyloid-to-DAM activation, DAM composition interaction, under-powered
  interaction DE, and no supported progression beyond composition. Next = S3
  trajectory + mechanism inline figures.
- 2026-07-02 Figure expansion S3 DONE -> `_trajectory.qmd` + `_mechanism.qmd`
  inline visual-density pass. Added 9 labelled `fig-*` chunks: trajectory
  pseudotime density, unit DAM-fraction/mean-pt scatter, channel/decomposition
  forest, robustness/omission audit; mechanism project-pathway heatmap, GO dot
  plot, Myc/NF-kB-family TF lollipop, NF-kB discordance tile, and kinase/run-index
  heatmap. Full `scripts/check.sh` green (tests warn=2, forced 124-chunk report
  render, tar_meta/render-log clean). Claims unchanged: trajectory interaction =
  DAM composition, no supported progression beyond composition; Myc supported;
  NF-kB attenuation discordant/not supported; Gsk3b covered but not recovered.
  Next = S4 cross-modality inline figures.
- 2026-07-02 Figure expansion S4 DONE -> `_crossmodality.qmd` inline
  visual-density pass. Added 8 labelled `fig-*` chunks: GeoMx volcanoes,
  GeoMx sensitivity/loss, bulk run-index, raw-vs-corrected phospho, bulk anchor
  heatmap, clearance-pair grid, symbol-modality matrix, and pathway-axis
  heatmap. Full `scripts/check.sh` green: tests warn=2, forced 140-chunk
  report render, tar_meta clean across 46 current targets/branches, render-log
  clean; report 8.24MB.
  Claim guards unchanged: GeoMx AOIs are blocked/repeated observations;
  SpatialDecon remains deferred; no full CCC is called; bulk layers are 24M
  hippocampus, not microglia-sorted; run-index sensitivity is load-bearing;
  symbol/pathway modality counts use broad modality classes. Next = S5 UX,
  visual QA, and close-out.
- 2026-07-02 Figure expansion S5 DONE -> CLOSED: enabled Quarto `lightbox: auto`
  in the embedded offline HTML and normalised every captioned figure chunk to a
  hyphenated `fig-*` label. Full `scripts/check.sh` green after forced 140-chunk
  render. Rendered HTML QA green: 42 figure blocks, 42 captions, 42 source
  `fig-*` labels, lightbox assets embedded, expected sections present, 0 external
  resource refs, 0 warning/error markers. Folded digest -> history.md; archived
  plan -> `.agent/completed/figure_expansion_plan_2026-07-02.md`; reset Active
  plan to none.
- 2026-07-02 Spatial decon follow-up OPENED -> `.agent/spatial_decon_followup_plan.md`.
  User requested roadmap units. Research = current P4 decon seams + v1 Arc L
  mining + current SpatialDecon/RCTD docs. Default plan = broad-first SpatialDecon
  with a compact full-snRNAseq reference profile, Q3-scaled background, no nuclei
  absolute counts while 42/91 sentinels remain, abundance DE/spatial residual
  audit, and report/synthesis rewiring only if earned. Active S0 route gate awaits
  user choice; no implementation yet.
- 2026-07-02 Spatial decon follow-up S0/S1 DONE: user chose the recommended
  broad-first SpatialDecon route with gated subpopulation attempt. Added
  SpatialDecon 1.22.0 to the repo lock; built `geomx_reference_profile` from the
  full snRNAseq RDS + retained microglia subpopulations + GeoMx gene overlap. Live
  compact target warning-clean/tar_meta clean: broad 15,919 genes x 6 profiles
  earned, subpopulation 16,079 genes x 8 profiles earned, Proliferative absent
  recorded, serialized 1.88 MB. Next = S2 decon fit.
- 2026-07-02 Spatial decon follow-up S2 DONE -> `geomx_decon`. Added
  SpatialDecon fit adapter around GeoMx RNA/data + Q3-scaled background,
  warning/message capture, beta/proportion/residual-QC normalisation, independent
  broad/subpopulation arms, and two-stage microglia-subpopulation assembly. Live target
  warning-clean/tar_meta clean but NOT earned: broad and subpopulation arms both
  blocked by the same 4/91 unresolved AOIs with beta_total=0
  (E03/E04/E05/G12 on DSP-1001660019825-A); 87 AOIs resolved, residual QC stored,
  nuclei absolute rescaling disabled (42 sentinels), no abundance claim. Full
  `scripts/check.sh` green across 48 current targets/branches. Next = S3
  blocked-abundance passthrough + residual audit.
- 2026-07-02 Spatial decon follow-up S3 DONE -> `geomx_abundance_de`. Added the
  abundance-DE orchestrator and nearest-neighbour residual audit. Live target
  warning-clean/tar_meta clean, 5.93 KB: broad/subpopulation/microglia-subpopulation DE all
  blocked by the same 4 unresolved SpatialDecon AOIs; canonical 5-contrast empty
  top tables preserve downstream shape. Residual audit earned descriptively for
  broad and subpopulation arms (91 AOIs x 4 slides; genotype-residualised RMS
  nearest-neighbour summaries), no abundance claim. Full `scripts/check.sh` green
  across 49 current targets/branches. Next = S4 report/synthesis integration.
- 2026-07-02 Spatial decon follow-up S4 DONE -> report/synthesis integration.
  Added compact `spatial_decon_report` (~1.5 KB live) over `geomx_decon`,
  `geomx_abundance_de`, and `geomx_reference_profile`; no beta matrices enter
  report bundles. `clearance_axis` now consumes the report handoff, so
  `crossmodality_report` and `synthesis_report` derive SpatialDecon status from
  the attempted fit: blocked/action attempted because 4 AOIs have beta_total=0;
  residual audit remains descriptive fit QC; nuclei-rescaled absolute counts and
  full CCC remain absent. `_crossmodality.qmd`, `_synthesis.qmd`, and `index.qmd`
  now say blocked after fitting, not skipped/reference-absent. Targeted tests,
  live target rebuild, forced 142-chunk report render, and tar_meta check were
  warning-clean. Next = S5 follow-up QA pass.
- 2026-07-02 Spatial decon follow-up S5 DONE -> steps complete. Adversarial QA
  accepted/fixed stale/ambiguous wording: GeoMx figure captions now distinguish
  the bio-unit-blocked primary DE model from the blocked SpatialDecon abundance
  state; the historical decon-preflight reason now points to the follow-up
  targets; memory/history/spine say SpatialDecon abundance is attempted but
  blocked, not merely skipped or missing-profile. No new figure section was
  added, so fig-label QA stayed unchanged. Full `scripts/check.sh` green. Next
  mode = CLOSE-OUT.
- 2026-07-02 Spatial decon follow-up CLOSED: close-out review of plan body +
  shipped code/prose found no accepted blocker. Folded digest -> history.md;
  archived plan -> `.agent/completed/spatial_decon_followup_plan_2026-07-02.md`;
  reset Active plan to none. Cohesive-story spine already carries the final
  blocked-fit SpatialDecon state. Final close-out gate green. Next = PLAN only
  after user confirms a new roadmap unit.
- 2026-07-03 Prose-to-figures reduction OPENED (later archived ->
  `.agent/completed/prose_to_figures_plan_2026-07-03.md`). User feedback: output report remains far
  too prose-heavy. Research = current qmd/report inventory, completed Figure
  expansion infrastructure, v1 archive report-shape mining, and current Quarto
  figure/crossref/diagram/tabset docs. Default plan = aggressive inline visual
  conversion: measure prose-only baseline, classify every human-facing block,
  replace paragraphs with figures/schematics/status matrices/caption microcopy,
  then verify before/after prose reduction plus rendered HTML QA. Alternatives
  recorded for visual-abstract-first, publication PNG gallery, or collapsible
  audit slimming. Next = S0 user route gate before report-source edits.
- 2026-07-03 Prose-to-figures reduction S0 DONE: user selected aggressive
  inline visual conversion. Provisional reduction floor = >=40% prose-only word
  count reduction, with S1 allowed to raise it after baseline inventory. No
  report sources edited before route selection. Next = S1 prose inventory and
  replacement manifest.
- 2026-07-03 Prose-to-figures reduction S1 DONE -> `scripts/prose_inventory.py`
  + `.agent/prose_replacement_manifest.tsv`. Baseline = 5,111 prose-only words
  across 119 human-facing blocks in `index.qmd` + `_*.qmd`; section headings
  are 33 kept navigation blocks, leaving 86 prose/caption blocks. All 86 have
  non-keep dispositions (`caption` 42, `figure` 32, `collapsed_audit` 7,
  `schematic` 5) and target slots. Selected reduction floor raised
  to >=55% (<=2,300 counted words; stretch <=1,800). Checked with
  `python3 -m py_compile scripts/prose_inventory.py` + inventory command; no
  report-source edits, so full render gate deferred to S2+ implementation.
  Next = S2 visual grammar and compact data contract.
- 2026-07-03 Prose-to-figures reduction S2 DONE -> compact visual grammar/data
  contract. Added `qc_figures` (QC modality/design/metric slots) and
  `report_visuals` (report spine schematic, synthesis visual abstract/source
  matrix, unsupported status grid, caveat glyphs, chapter evidence boards);
  enriched existing chapter figure targets with board/alias slots
  (microglia summary/composition, trajectory logic/decomposition/concordance,
  mechanism status, cross-modality status/count aliases). `visual_reduction_slot_map`
  + `visual_slot_coverage` cover every S1 `figure`/`schematic` disposition.
  Live build warning-clean and compact: `qc_figures` 4.37 KB, `report_visuals`
  4.28 KB; enriched chapter targets remain small. Focused figure test green.
  Next = S3 synthesis and overview conversion.
- 2026-07-03 Prose-to-figures reduction S3 DONE -> `index.qmd` + `_synthesis.qmd`
  visual-first conversion. Overview now draws `fig-report-spine-schematic` from
  compact `report_visuals` instead of prose chapter previews. Synthesis now loads
  `synthesis_report` + `report_visuals`, keeps one source-derived answer
  sentence, replaces the compact evidence table/status prose with
  `fig-synthesis-visual-abstract`, `fig-synthesis-evidence-map`, and
  `fig-synthesis-status`. Local counted prose: Overview+synthesis 305 -> 46
  words (85% reduction, clears >=55% floor). Forced render + full
  `scripts/check.sh` green. Next = S4 result chapter conversion.
- 2026-07-03 Prose-to-figures reduction S4 DONE -> result chapter conversion.
  `_microglia.qmd`, `_trajectory.qmd`, `_mechanism.qmd`, and
  `_crossmodality.qmd` now use compact visual boards, existing figure targets,
  and one-line captions instead of table-heavy explanatory prose. Counted words:
  microglia 1,438 -> 251, trajectory 1,317 -> 186, mechanism 755 -> 148,
  cross-modality 958 -> 195; whole report now 1,164 words (77% below S1).
  Preserved visible null/blocked states: progression not supported, NF-kB
  discordant/not supported, Gsk3b not recovered, SpatialDecon abundance blocked,
  full CCC absent, bulk run-index caveat. Manifest writer now emits `.` for
  empty cells to avoid trailing whitespace. Full `scripts/check.sh` green. Next =
  S5 visual QA and close-out.
- 2026-07-03 Prose-to-figures reduction S5 DONE -> CLOSED. Final source
  inventory = 1,164 counted words / 117 blocks vs 5,111 / 119 baseline (77%
  reduction; clears >=55% floor and <=1,800 stretch). Rendered HTML QA:
  49 figures / 49 captions, no >32-word captions, no duplicate IDs, no broken
  internal anchors, no external href/src refs, lightbox present, no visible
  warning/error markers, no underscored rendered fig-* IDs. Claim-parity review
  accepted no blockers: nulls/blocked states remain target-derived in boards,
  captions, or visible panels. Full `scripts/check.sh` green. Folded digest ->
  history.md; archived plan -> `.agent/completed/prose_to_figures_plan_2026-07-03.md`;
  reset Active plan to none. Next = PLAN only after user confirms a new roadmap
  unit.
- 2026-07-03 Report top sections REMOVED by direct user request. Deleted
  `_synthesis.qmd`, `R/synthesis.R`, `tests/test_synthesis.R`, `synthesis_report`,
  plus the Overview body, `report_visuals` target/helper, and report-spine manifest
  slot. `index.qmd` is now YAML + includes only; rendered order is QC -> result
  chapters. Local mechanism/cross-modality "Synthesis" tails were renamed to status
  headings.
- 2026-07-03 Figure-caption-only report OPENED -> `.agent/figure_caption_only_plan.md`.
  User request = extreme prose elimination: no prose, just figures/captions.
  Research = current qmd inventory, completed prose-to-figures plan, top-section
  removal commits, Quarto figure/code docs, WCAG non-text-content guidance.
  Current source inventory = 1,119 words / 109 blocks: paragraph 33/594,
  caption 45/388, headings 31/137. Selected route = strict visible-path
  conversion: headings + figures + captions only; no visible markdown prose,
  `knitr::kable()` tables, or `cat()` provenance; `fig-alt` remains source-level
  accessibility text. Next = S1 strict inventory gate.
- 2026-07-03 Figure-caption-only S2 DONE: `_qc.qmd` now loads only compact
  `qc_figures` and renders headings + 5 figures/captions: modality/GeoMx/sample-key
  grid, 16-cell genotype-by-batch heatmap, depth and fraction histograms, and
  structural/bounds status tiles. Removed QC markdown prose, `knitr::kable()`
  tables, and visible `cat()` provenance. QC source strict gate passes; rendered
  strict HTML has no QC table/stdout blockers, with remaining blockers reserved
  for result chapters + YAML author metadata. Full `scripts/check.sh` green.
  Next = S3 result chapter body-prose removal.
- 2026-07-03 Figure-caption-only S3 DONE: result chapters now have zero
  paragraph/list/table source blockers: headings + figure chunks/captions only.
  Removed chapter opener/status-tail prose from `_microglia.qmd`, `_trajectory.qmd`,
  `_mechanism.qmd`, `_crossmodality.qmd`; hid microglia sccomp diagnostics as a
  render-fatal check; converted trajectory provenance/per-cell status from stdout
  into a figure later removed by curation. Captions pass S3 rule (53 captions, max 12 words,
  median 8). Strict rendered HTML is down to the known YAML author paragraph only;
  no result-body paragraphs, tables, text-only outputs, or stdout provenance remain.
  Full `scripts/check.sh` green. Next = S4 caption/alt/HTML QA + author metadata.
- 2026-07-03 Figure-caption-only S4 DONE: removed YAML author metadata and visible
  code UI, added/verified `fig-alt` for every captioned figure (48 `fig-cap` /
  48 `fig-alt`; caption max 12 words), and replaced the report target with
  `render_report()` so Quarto renders quietly-false through the DAG then repairs
  embedded lightbox anchors from missing local `index_files/figure-html/*.png`
  hrefs to data URIs. Rendered strict HTML QA green: 48 figures/captions/img alts,
  48 data-URI lightbox links, 0 local figure refs, 0 duplicate IDs, 0 visible
  paragraphs/tables/stdout/text-only outputs, 0 code UI, 0 warning/error markers.
  Full `scripts/check.sh` green across 52 current targets/branches. Next = S5
  claim-parity close-out.
- 2026-07-03 Figure-caption-only S5 CLOSED: adversarial claim-parity review
  found no accepted blockers (removed prose was redundant with headings, encoded
  in plotted boards/status panels/captions, or intentionally dropped audit text).
  Strict gate green: source has 0 paragraphs/lists/tables and rendered HTML has
  0 body prose/tables/stdout/text-only outputs. DOM QA: 48 figures / 48 captions /
  48 nonblank alts, 48 data-URI lightbox hrefs, 0 local figure hrefs, 0 duplicate
  IDs, 0 external refs, 0 code UI. Full `scripts/check.sh` green. Archived plan
  to `.agent/completed/figure_caption_only_plan_2026-07-03.md`; reset Active plan
  to none. Next = PLAN only after user confirms a new roadmap direction.
- 2026-07-03 Box-figure curation DONE (ad hoc user task): removed visible pure
  box/status figures from the caption-only report even where that omits prepared
  audit information. Dropped QC modality/design/bounds boxes, microglia summary +
  composition-concordance tiles, trajectory logic/method boards, mechanism status/
  project-set/NF-kB-discordance boxes, and cross-modality status/run-index/anchor/
  clearance/integrated-divergence matrices. Retained 31 data-rich figures:
  distributions, UMAPs, compositions, forests, volcanoes, densities, scatter,
  dot/lollipop plots, and quantitative heatmaps. `figure_manifest()` now pins 18
  curated expansion ids; `visual_slot_coverage()` accepts the empty figure/
  schematic prose-replacement set after manifest regeneration. Strict rendered
  HTML QA green (31 captions, 0 body prose/tables/stdout/text-only outputs,
  0 removed-label hits); full `scripts/check.sh` green across 52 targets/branches.
- 2026-07-03 Figure story layout DONE (ad hoc user task): kept caption-only
  surface and conventional figure classes, but reshaped the story. Added
  `qc_figures$study_design` + `fig-qc-study-design` (2x2 genotype grid +
  modality support bars); expanded `fig-microglia-composition-shift` to a
  3-panel composition / tau x amyloid DAM-response / DAM-score figure; added
  panel tags to composite microglia/trajectory figures; added an lm trend to the
  trajectory DAM-fraction scatter; added `fig-mechanism-nfkb-primary` as a
  primary-score lollipop instead of a status box. Strict rendered QA: 33 figures /
  33 captions / 33 nonblank alts, 0 duplicate IDs, 0 external refs, 0 visible body
  prose/tables/stdout/text outputs; full `scripts/check.sh` green across
  52 targets/branches.
- 2026-07-03 Visual maturity pass DONE (ad hoc user task): addressed user
  feedback that the report felt juvenile. Shared plot layer now pins a restrained
  graphite/teal/ochre/wine palette, binary/direction helper scales, and neutral
  sequential count fills; theme.scss uses deep-ink/teal/slate chrome. Render
  layer changes: QC design circles -> tile matrix; mechanism TF lollipop ->
  `fig-mechanism-tf-focus` heatmap; NF-kB primary lollipop -> signed bars;
  trajectory DAM-fraction batch text -> shape coding; direct-labelled count bars
  muted. Visual PDF QA inspected opening + mechanism pages; focused tests green;
  DOM QA = 33 figures / 33 captions / 33 nonblank alts / 0 duplicate IDs; full
  `scripts/check.sh` green across 52 targets/branches.
- 2026-07-03 Color saturation pass DONE (ad hoc user task): addressed user
  feedback that the visual maturity pass went too dull. Palette raised to
  steel-blue / teal / amber / cranberry / violet accents while preserving the
  report's journal-style plot grammar. Visual PDF QA inspected opening +
  mechanism pages; DOM QA = 33 figures / 33 captions / 33 nonblank alts /
  0 duplicate IDs; full `scripts/check.sh` green across 52 targets/branches.
- 2026-07-03 Cross-modality narrative figure pass CLOSED: replaced the generic
  cross-modality dashboard panels with target-backed evidence plates for
  amyloid-response, synaptic/clearance context, and interaction boundaries; kept
  GeoMx/phospho boundary panels; added the final closing model; and fixed print
  legend overrun in the trajectory DAM-fraction scatter plus the cross-modality
  evidence plates. Strict caption-only HTML QA = 37 figures / 37 captions /
  37 nonblank alts, zero body prose/tables/stdout blockers, no stale dashboard
  ids, no duplicate ids, and no local/external refs. Chromium PDF content-page
  QA clean for trajectory, cross-modality plates, GeoMx, phospho, and closing
  model; Chromium still emits one trailing blank page. Full `scripts/check.sh`
  green across 53 current targets/branches. Active plan reset to none.
- 2026-07-03 infra agent-surface revert: restored the Claude/Serena project
  surface (`CLAUDE.md`, `.claude/commands/`, `.claude/settings.json`,
  `.serena/`) and slash-command session/review entry points; `.agent/context.sh`
  again reads Claude Code transcripts. Removed the Codex-only prompt/skill
  surface while preserving the current roadmap/report state.
- 2026-07-04 infra repo hygiene: renamed the Quarto `output-dir` `_report/` -> `report/` (the ONE
  free-choice name; the other leading-`_` names -- `_quarto.yml`, `_targets.R`/`_targets/`, `_*.qmd` --
  are tool contracts, KEPT). Synced `_quarto.yml`, `R/report.R` html_file, `.gitignore` `/report/`,
  map.md + memory.md; `report` is `format="file"` so the target value
  just re-points -- no DAG change. Redirected OmnipathR API logs -> `storage/cache/omnipathr-log/` via
  `.Rprofile` `options(omnipathr.logdir=)` set before load (was littering repo-root `./omnipathr-log/`).
  Purged stray gitignored working-tree litter (`sccomp_draws_files/` ~174M, root `omnipathr-log/`, old
  `_report/`; ~271M freed) + deleted stale `.agent/crossmodality_narrative_manifest.md` (completed-S1 audit
  manifest, no live refs). No committed junk existed -- all tracked files legit, clutter was gitignored
  caches. New memory.md "Repo layout" section locks the rule. Gate green before + after.
- 2026-07-06 Report figure-cut / chapter teardown DONE (user task, confirmed full-teardown scope): permanent
  pivot to a FIVE-figure microglia+trajectory report. Cut the rendered document to exactly old figures
  2,3,4,7,16 (renumber 1-5): fig-microglia-subpopulation-markers, fig-microglia-umap, fig-microglia-umap-subpopulation,
  fig-microglia-unit-composition, fig-trajectory-pt-density. DELETED the qc/mechanism/cross-modality chapters +
  every one of their targets + R modules + tests: `_qc.qmd`/`_mechanism.qmd`/`_crossmodality.qmd`, `R/mechanism.R`,
  `R/crossmodality.R`, `tests/test_mechanism.R`, `tests/test_crossmodality.R`; dropped all geomx/proteome/phospho +
  mechanism + cross-modality + spatial-decon + story + qc targets from `_targets.R` (and the story/qc/mechanism/
  crossmodality builders from `R/figures.R`). `_targets.R` now = 19 targets (microglia P1 + trajectory P2 pipelines
  intact, cached, gate-cheap). `figure_manifest()` trimmed 24->11 rows (microglia 7 + trajectory 4);
  `tests/test_figures.R` updated in lockstep (11-row assertion). `index.qmd` = YAML + `{{< include _microglia.qmd >}}`
  + `{{< include _trajectory.qmd >}}`; `render_report()`/`report_sources` trimmed. Full `scripts/check.sh` green:
  render processed exactly 5 figure chunks, microglia_figures/trajectory_figures rebuilt cheaply, tar_meta clean
  across 19 targets, render-log clean. Findings 3 (mechanism) + 4 (cross-modality) science preserved in git history
  + this Ledger, dropped from report + pipeline. Residual dead code, deferred as a separate hygiene unit: (a) io.R geomx/proteome/phospho loaders +
  constants.R data_paths + rbc_marker_symbols -- dead but KEPT (test_io-covered) for possible modality re-add
  (memory.md "Raw data" note); (b) figures.R story/mechanism/crossmodality builders + their `.fig_*` helpers +
  qc_figure_data -- pure dead weight referencing permanently-deleted targets, NOT test-covered, removable; plus
  plot.R concordance_plot (P4-only). The two kept builders (microglia/trajectory) still emit unrendered manifest
  slots (11 built, 5 rendered).
- 2026-07-06 infra Codex project restore: reversed the `f946306` Claude-workflow surface while preserving
  newer report/pipeline state. Tracked contract is again Codex-only: `AGENTS.md` canonical, `$session-prompt`
  via `.agents/skills/session-prompt` + `.codex/prompts/session.md`, `.agent/context.sh` reads Codex JSONL
  sessions, and `CLAUDE.md` / `.claude/` / `.serena/` are retired untracked/ignored state. Synced memory/map/
  history wording from deny-Read/Serena back to Codex read-economy.
- 2026-07-06 Off-diagonal pathway figure DONE (ad hoc user task): added `fig-modality-offdiag-pathways`
  immediately after Figure 6. `modality_logfc_scatter_data` now carries gene-symbol tokens, signed
  interaction `x-y`, and `pathways$summary`: mouse MSigDB GO Biological Process overlap over the top 250
  unique off-diagonal genes/proteins per method by |x-y| (descriptive, not a restored mechanism chapter).
  `_modality.qmd` renders a modality x GO-BP bubble matrix (size = overlap count, fill = mean x-y, row
  subtitles = leading hit genes; FDR support now marked by an offset asterisk, not a hard-to-see ring). Target
  remains compact (~1.7MB); report now has 7 captions/alts.
- 2026-07-06 Figure 7 functional-score rework DONE (ad hoc user task): replaced the GO-BP enrichment-style
  bubble matrix with `fig-modality-functional-scores`. `modality_logfc_scatter_data` now returns `groups$summary`:
  broad functional-role gene sets assembled from GO-BP keyword unions over the same genes/proteins labelled in
  Figure 6 (top 12 display labels per method by |x-y| after duplicate-label collapse), scored by aggregate
  NLGF_MAPTKI amyloid logFC, aggregate NLGF_P301S amyloid logFC, and `delta=x-y`. The plot is a per-modality
  dumbbell score facet (connected MAPTKI/P301S points; segment colour = P301S-MAPTKI; point size = Figure 6
  labels) and displays no enrichment/FDR result.
- 2026-07-06 Figure 7 phosphosite parent-gene scoring DONE (ad hoc user follow-up; superseded by parent-protein
  mean scatter collapse below): phosphosite labels still enter
  through the Figure 6 label rule, but group scoring substitutes the best-fit parent gene (`gene` when usable,
  otherwise site_id prefix) and collapses duplicate labelled phosphosites to the highest-|x-y| site per parent.
  Legend/caption now say scored genes/proteins/items rather than raw Figure 6 labels.
- 2026-07-06 Figure 6/7 empirical off-diagonal categorization DONE (ad hoc user follow-up; threshold superseded by
  next ledger item): replaced fixed top-12 Figure 6 labels with empirical threshold labels after finite filtering
  and duplicate display-label collapse. Reworked Figure 7 from shared functional-group rows to
  modality-specific/free-y primary categories over those empirical outliers: GO-BP keyword-union roles first,
  explicit fallback categories for unmapped/predicted/olfactory labels, phosphosites scored through parent genes.
- 2026-07-06 Figure 6/7 standardized off-diagonal cutoff DONE (ad hoc user correction; live threshold/counts
  superseded by parent-protein mean scatter collapse below): replaced the per-modality
  empirical threshold with one pooled Q99.8 cutoff over finite `|x-y|` distances from all four Figure 6 methods
  (`|x-y| >= 4.480587` live). Live labels now reflect shared absolute distance: snRNAseq 0, GeoMx 0, Proteome 18,
  Phospho 87; Figure 7 scores Proteome 18 + Phospho 59 parent-gene/protein items after phosphosite collapse.
- 2026-07-06 Figure 6 phosphoproteomics readability DONE (ad hoc user correction; threshold/counts superseded by
  Q99 label relaxation below): phospho DE remains site-level,
  but `modality_logfc_scatter_data` collapses finite phosphosite logFC pairs to best-fit parent-protein means before
  pooled cutoff, labels, and Figure 7 scoring. Compact target at the time: phospho 17707 finite sites -> 3092 parent
  proteins; pooled Q99.8 `|x-y| >= 3.307601`; Figure 6 labels snRNAseq 5 / GeoMx 0 / Proteome 57 / Phospho 20;
  Figure 7 rows snRNAseq 4 / GeoMx 0 / Proteome 8 / Phospho 5.
- 2026-07-06 Figure 6 Q99 label relaxation DONE (ad hoc user correction): standardized pooled cutoff relaxed from
  Q99.8 to Q99 to capture more transcriptomic off-diagonal points while retaining the parent-protein phospho collapse.
  Live compact target: pooled Q99 `|x-y| >= 2.009746`; Figure 6 labels snRNAseq 65 / GeoMx 0 / Proteome 247 /
  Phospho 97; Figure 7 rows snRNAseq 7 / GeoMx 0 / Proteome 9 / Phospho 9.
- 2026-07-06 Figure 7 endpoint colour separation DONE (ad hoc user correction): kept global tau-background
  palette unchanged, but Figure 7 aggregate-score endpoint fills now use a dedicated high-separation blue/burnt-orange
  pair instead of the similar blue/green pair; stacked bottom guides to prevent the size legend from clipping.
- 2026-07-06 Figure 6/7 within-method Q99 DONE (ad hoc user correction; supersedes pooled-Q99 label selection):
  empirical off-diagonal labels and Figure 7 scored categories now use each method's own Q99 `|x-y|` cutoff, while
  preserving phospho parent-protein collapse. Live compact target: cutoffs snRNAseq 1.688 / GeoMx 0.774 /
  Proteome 3.899 / Phospho 2.776; Figure 6 labels snRNAseq 145 / GeoMx 199 / Proteome 34 / Phospho 31; Figure 7
  rows snRNAseq 10 / GeoMx 10 / Proteome 8 / Phospho 7. High-label scatter panels use smaller labels/points.
- 2026-07-07 Non-snRNAseq modality landscapes DONE (ad hoc user correction): replaced the 3 descriptive
  GeoMx/proteome/phosphoproteome prose paragraphs with 3 standalone target-derived figures:
  `fig-modality-geomx-landscape`, `fig-modality-proteome-landscape`, `fig-modality-phospho-landscape`. Each combines
  the modality amyloid-effect scatter (standalone label cap 30 for readability) with same-modality Q99 off-diagonal
  functional-category scores. Rendered report now has 10 unique captioned figure containers; full `scripts/check.sh`
  green; Chromium PDF QA pages 5-6 clean.
- 2026-07-07 Non-snRNAseq modality-native figures DONE (ad hoc user correction; supersedes the generic landscapes
  above): kept the same 3 chunk ids but replaced their contents with assay-typical descriptive figures. GeoMx now
  shows a slide-faceted spatial AOI score plate; proteome now
  shows 16-run PCA + a protein-level NLGF_P301S-vs-P301S volcano; phospho now shows a phosphosite-level
  NLGF_P301S-vs-P301S volcano + top-site z-score heatmap. `modality_scatter_figures` now carries `$descriptive`;
  `geomx_de` carries compact `$spatial`; `phospho_de_24m` keeps its filtered matrix for the heatmap. Targeted tests
  green; report render green; Chromium PDF QA pages 5-6 clean.
- 2026-07-07 GeoMx informative plate DONE (ad hoc user correction): refined only
  `fig-modality-geomx-landscape`. The GeoMx panel is now a spatial + quantitative plate: slide-faceted AOI map
  (score fill, genotype ring, AOI-area size, high-magnitude ROI labels), AOI-score distribution by genotype, and
  top score-gene driver logFCs for the two amyloid-background contrasts. `geomx_de$spatial$aoi` now carries
  `segment` + `aoi_area`; plot/tests require `$spatial$genes`. Targeted tests green; report render green; Chromium
  PDF page 5 clean.
- 2026-07-07 Lean iteration infrastructure cut DONE (user task): optimized for faster report iteration over
  robustness. Removed committed tests, Python/uv config, prose inventory, CmdStan/sccomp composition arm,
  stageR layer, per-subpopulation pseudobulk, and dead figure/plot/constants helpers that did not feed the current
  10-figure report. `microglia_report_data()` now emits replicate-unit composition directly; `microglia_figures`
  and `trajectory_figures` carry only rendered slots; `_targets.R` is the 29-target report DAG. `scripts/check.sh`
  is now a lean report gate: optional `rv sync`, forced `report` rebuild, target metadata scan, render-log scan.
  `AGENTS.md`, session prompt/skill, memory, and map now state the lean posture.
- 2026-07-08 GeoMx exploratory figures PLAN OPENED (direct user task): web search
  found common GeoMx exploratory figures in Bruker GeoScript/User Manual, standR workflow/manual, and
  SpatialDecon docs: NGS/segment QC, normalization/RLE, PCA/MDS/loadings, gene detection,
  heatmaps, spatial overlays, volcano/MA, ROI/segment structure, and decon/status diagnostics.
  Opened the plan now archived at `.agent/completed/geomx_exploratory_figures_plan_2026-07-08.md`;
  roadmap Active plan implemented one figure per session (S1-S9). Next = EXECUTE S1.
- 2026-07-08 GeoMx exploratory figures S1 DONE -> `fig-geomx-qc-atlas`. Added
  `geomx_qc_descriptor()` into `geomx_de` and `geomx_qc_atlas_plot()` into `_modality.qmd`.
  The atlas shows library size, detected genes, nuclei, AOI area, negative background, and Q3 factor
  by slide/genotype/segment plus flag counts. Live QC: 91 AOIs, 53 with >=1 descriptive flag;
  detected genes are constant across AOIs and correctly produce 0 low-gene flags; nuclei sentinels =
  42. No exclusion or claim change. Targeted `geomx_de` / `modality_scatter_figures` builds and report
  render green; next = EXECUTE S2.
- 2026-07-08 GeoMx exploratory figures S2 DONE -> `fig-geomx-normalization-rle`. Added
  `geomx_normalization_descriptor()` into `geomx_de`, passed it through
  `modality_scatter_figures$descriptive$GeoMx$normalization`, and rendered
  `geomx_normalization_rle_plot()` in `_modality.qmd`. The figure shows raw vs TMM logCPM
  per-AOI quantiles, TMM RLE quantiles by slide/genotype/segment, Q3 factor vs negative-control
  background, and the saved voom mean-variance trend from the primary slide-adjusted,
  bio-unit-blocked fit. Live data: 91 AOIs, 19,959/19,963 genes kept by
  `filterByExpr(min.count=5)`, Q3/background Spearman rho = 0.994. No AOI exclusion or claim
  change. Hardened `report_sources` so R helper edits invalidate `report`. Focused target/render
  build green; Chromium PDF page-6 QA clean; next = EXECUTE S3.
- 2026-07-08 GeoMx exploratory figures S3 DONE -> `fig-geomx-ordination`. Added
  `geomx_ordination_descriptor()` into `geomx_de`, passed it through
  `modality_scatter_figures$descriptive$GeoMx$ordination`, and rendered
  `geomx_ordination_plot()` in `_modality.qmd`. The figure shows slide-faceted PCA and
  classical MDS over TMM-logCPM top-variable filter-passing genes, a PCA scree curve, and
  PC1/PC2 signed loading genes. Live data: 91 AOIs, 19,959/19,963 genes kept by
  `filterByExpr(min.count=5)`, 2,000 variable genes used; PC1 = 19.36%, PC2 = 6.94%.
  No AOI exclusion or claim change. Focused target build green and warn=2 plot smoke-render
  clean; next = EXECUTE S4.
- 2026-07-08 GeoMx exploratory figures S4 DONE -> `fig-geomx-gene-detection`. Added
  `geomx_gene_detection_descriptor()` into `geomx_de`, passed it through
  `modality_scatter_figures$descriptive$GeoMx$gene_detection`, and rendered
  `geomx_gene_detection_plot()` in `_modality.qmd`. The figure shows mean TMM logCPM vs
  AOI detection fraction, the existing low-coverage `filterByExpr(min.count=5)` decision,
  microglia identity/homeostatic/DAM marker measurability, and highest-detected WTA genes.
  Live data: 91 AOIs, 19,959/19,963 genes filter-passing, 4 low-coverage genes, marker
  genes present/pass = Microglia 8/8, Homeostatic 10/10, DAM 18/18. No AOI exclusion,
  DE change, or report-claim change. Focused target/render build green; full
  `TZ=UTC CHECK_SKIP_SYNC=1 scripts/check.sh` green; Chromium PDF page-8 QA clean.
  Next = EXECUTE S5.
- 2026-07-08 GeoMx exploratory figures S5 DONE -> `fig-geomx-sample-heatmap`. Added
  `geomx_sample_heatmap_descriptor()` into `geomx_de`, passed it through
  `modality_scatter_figures$descriptive$GeoMx$sample_heatmap`, and rendered
  `geomx_sample_heatmap_plot()` in `_modality.qmd`. The figure shows clustered AOI
  tracks for genotype/slide/segment/bio-unit/ROI, the signed GeoMx amyloid-response
  score, and a capped top-variable-gene row-z heatmap. Live data: 91 AOIs,
  19,959/19,963 genes filter-passing, 40 top-variable genes displayed, z clipped
  at +/-2.5. No AOI exclusion, DE change, or report-claim change. Focused target
  build green; rendered report has 15 figures / 15 nonblank alts; full
  `CHECK_SKIP_SYNC=1 scripts/check.sh` green; Chromium PDF page-9 QA clean.
  Next = EXECUTE S6.
- 2026-07-08 GeoMx exploratory figures S6 DONE -> `fig-geomx-spatial-program-overlays`.
  Added `geomx_spatial_program_descriptor()` into `geomx_de`, passed it through
  `modality_scatter_figures$descriptive$GeoMx$spatial_programs`, and rendered
  `geomx_spatial_program_overlay_plot()` in `_modality.qmd`. The figure shows
  coordinate-only AOI maps for Homeostatic/DAM/IFN/MHC_APC signatures plus Apoe/Trem2
  single genes, with genotype median/IQR summaries. Live data: 91 AOIs x 6 programs
  (546 rows), 19,959/19,963 genes filter-passing, scored features = Homeostatic 10/10,
  DAM 18/18, IFN 12/12, MHC/APC 7/7, Apoe 1/1, Trem2 1/1. No AOI exclusion, DE change,
  or report-claim change. Focused target/render build green; rendered report has
  16 figures / 16 nonblank alts; full `scripts/check.sh` green; Chromium PDF page-10
  QA clean. Next = EXECUTE S7.
- 2026-07-08 GeoMx exploratory figures S7 DONE -> `fig-geomx-contrast-diagnostics`.
  Added `geomx_contrast_diagnostic_descriptor()` into `geomx_de`, passed it through
  `modality_scatter_figures$descriptive$GeoMx$contrast_diagnostics`, and rendered
  `geomx_contrast_diagnostic_plot()` in `_modality.qmd`. The figure shows five
  canonical contrast volcano facets, five matching MA facets, signed FDR-support
  counts, and top interaction genes with confidence intervals. Live data: 99,795
  contrast-feature rows (19,959 genes x 5 contrasts), 20 deterministic labels,
  support counts at FDR <= 0.10 = tau_alone 4 up/0 down, nlgf_in_maptki 5,258
  up/1,559 down, nlgf_in_p301s 1,386 up/541 down, tau_in_nlgf 126 up/270 down,
  interaction 45 up/117 down. No new inference model, AOI exclusion, DE change,
  or report-claim change. Focused target build green; full `scripts/check.sh`
  green; rendered HTML has 17 figures / 17 captions; Chromium PDF page-11 QA
  clean. Next = EXECUTE S8.
- 2026-07-08 GeoMx exploratory figures S8 DONE -> `fig-geomx-roi-segment-replicates`.
  Added `geomx_roi_replicate_descriptor()` into `geomx_de`, passed it through
  `modality_scatter_figures$descriptive$GeoMx$roi_replicates`, and rendered
  `geomx_roi_replicate_plot()` in `_modality.qmd`. The figure shows genotype-by-
  bio-replicate AOI support, AOI counts per duplicateCorrelation block, and
  AOI-pair expression-correlation distributions by same-bio-unit / same-genotype
  different-unit / different-genotype relationship. Live data: 91 AOIs, 15/16
  expected bio-units present (MAPTKI:1 absent), one segment level (`Segment 1`),
  all 15 present bio-units are repeated-observation blocks, block sizes 2-7 AOIs,
  duplicateCorrelation consensus = 0.0085, 19,959/19,963 genes filter-passing,
  2,000 top-variable genes used for AOI-pair correlations, pair counts = 248 /
  763 / 3,084. No paired segment-difference panel is drawn because all AOIs use
  one segment level. No AOI exclusion, DE change, or report-claim change. Focused
  target build and report render green; full `CHECK_SKIP_SYNC=1 scripts/check.sh`
  green; rendered HTML has 18 figure containers / 18 captions / 18 nonblank image
  alts. Next = EXECUTE S9.
- 2026-07-08 GeoMx exploratory figures S9 DONE -> `fig-geomx-decon-feasibility`.
  Added `geomx_decon_feasibility_descriptor()` into `geomx_de`, passed it through
  `modality_scatter_figures$descriptive$GeoMx$decon_feasibility`, and rendered
  `geomx_decon_feasibility_plot()` in `_modality.qmd`. The figure shows candidate
  marker-set coverage, AOI precondition/blocker coordinates, genotype-level
  blocker counts, and marker-coherence proxy residuals. Live data: 91 AOIs,
  19,959/19,963 genes filter-passing, 8/8 marker components covered with >=2
  filter-passing marker genes, no live `SpatialDecon` dependency or
  reference-profile/beta/abundance-DE target, and AOI input-status bins = 37 no
  local blocker / 3 low-input tail / 9 background-Q3 tail / 42 absolute-count
  blocked by nuclei sentinel. The figure is a blocked diagnostic, not a
  deconvolved cell-abundance claim. Focused target build green; isolated plot
  render green under `warn=2`; full `CHECK_SKIP_SYNC=1 scripts/check.sh` green;
  rendered HTML has 19 captions / 19 nonblank image alts; Chromium PDF page-13
  QA clean. Next = CLOSE-OUT active plan.
- 2026-07-08 GeoMx native figure curation DONE (ad hoc user task): retained only former
  Figure 10 (`fig-geomx-sample-heatmap`) as the GeoMx panel under "Non-snRNAseq
  modality-native figures"; removed visible GeoMx QC/normalization/ordination/gene-detection/
  spatial-program/contrast/ROI/decon diagnostics plus `fig-modality-geomx-landscape`. Kept the
  proteome and phosphoproteome native figures plus the two modality-context figures. `geomx_de`
  and `modality_scatter_figures` now carry only `sample_heatmap` as the GeoMx descriptive
  payload; rendered surface expected = 10 figures, target count unchanged at 29.
- 2026-07-09 Report HTML polish DONE (ad hoc user task): added simple visible `Figure 1` ...
  `Figure 10` headings across the 10 qmd figure chunks, kept `toc: false`, and shrank the native
  Quarto `details.code-fold > summary` control plus expanded folded-code text in `theme.scss`.
  Report surface remains title-free, caption-free, body-prose-free, and global-code-menu-free.
- 2026-07-09 Figure 9/10 off-diagonal tightening DONE (ad hoc user task): tightened the four-method
  amyloid-response scatter label rule from within-method `max(q0.99, top-24)` to `max(q0.99, top-8)`.
  The functional-group score panel still scores the exact Figure 9 selected features, now with 10 displayed
  categorized rows and no rows with `|P301S - MAPTKI| < 0.5` under the live sensitivity check.
- 2026-07-09 Figure 9/10 modality-aware threshold correction DONE (ad hoc user task): replaced the
  uniform top-8 off-diagonal label budget with modality-specific budgets: snRNAseq 16, GeoMx 4,
  proteome 24, phosphoproteome 24. Figure 10 now filters displayed aggregate rows to
  `|P301S - MAPTKI| >= 0.5`, so near-overlap aggregate categories stay out while the stronger
  proteomic/phosphoproteomic and snRNAseq outliers return.
- 2026-07-09 Figure 9/10 shared-cutoff correction DONE (ad hoc user task): replaced modality-specific
  Figure 9 feature thresholds with one round shared cutoff, `|x-y| >= 2.0`, chosen just above the
  live GeoMx maximum so GeoMx contributes no off-diagonal selected features. Figure 9 caps displayed
  text labels at 24 per panel for readability only; Figure 10 scores all shared-cutoff selected
  features and keeps the global aggregate display floor `|P301S - MAPTKI| >= 0.5`.
- 2026-07-09 Figure 9 shared-scale cutoff bands DONE (ad hoc user task): Figure 9 now draws all four
  amyloid-response scatter facets on one shared square coordinate range, so the two dotted `|x-y|=2.0`
  cutoff bands sit the same visual distance from the center identity diagonal in every modality.
- 2026-07-09 Figure 9 threshold legend DONE (ad hoc user task): added a collected line legend for
  the dotted off-diagonal cutoff bands, labelled `threshold: |x-y| >= 2 log2FC`.
- 2026-07-09 Figure 9 threshold legend sizing DONE (ad hoc user task): enlarged the collected
  threshold legend text/key to match the report axis-label scale.
- 2026-07-09 Figure 9 threshold text-only DONE (ad hoc user task): removed the dotted swatch from
  the threshold note and set its text between the axis-label scale and the larger 1.35x draft.
- 2026-07-09 Figure 10 title polish DONE (ad hoc user task): centered the functional-group
  plot guide stack, shifted the inner plot title right to align with it, enlarged the title,
  and added vertical title padding; the special centered Figure 10 heading also gets slightly
  larger padded spacing.

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
84% 168K/200K    06-30  P1-S2 subpopulation annotate + prune +review (612dbd0, 5992e39)
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
