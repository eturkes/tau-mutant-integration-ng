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

## Active plan: P1 snRNAseq microglia core -> `.agent/p1_snrnaseq_plan.md` (opened 2026-06-29)
5 steps: S1 reprocess+cluster (SCT+Harmony-batch-only+Louvain) | S2 substates (UCell calibrated-argmax: homeostatic/
DAM/IFN/proliferative + aux MHC/APC) + contaminant prune | S3 composition (sccomp + propeller) | S4
pseudobulk DE (whole-MG + per-substate, voomQW + robust eBayes + stageR) | S5 report `_microglia.qmd` + close.
Gate-locked (2026-06-29, SOTA-researched): single-cell DE DROPPED (pseudobulk sole inference; supersedes
de_pb.R de_sc forward-ref); composition = sccomp+propeller; normalisation = SCT-v2 (v1 continuity). Science:
headline = amyloid->DAM (robust, microglia-led); the tau x amyloid interaction = v1-prior null -> P1 reports it
WITH a power/effect-size stmt (outcome-open), synergy = a rate effect handed to P2. Carry Thrupp 2020 caveat
(snRNA under-detects ~18% DAM genes; score, not threshold).
S1 DONE 2026-06-29: `microglia_processed` (SCT-v2/glmGamPoi -> Harmony[batch] -> Louvain 12 clusters @res0.4 ->
UMAP); marker separation confirmed post-Harmony (homeostatic/DAM/IFN/prolif distinct argmax); re-run STABLE
(observed ARI=1.0, recorded threads); gate green. **NEXT:** S2 substate annotation + QC prune (rv add UCell).

## Backlog - phased build (each phase = closeable increments; mine archive_digest per phase)
- P0 Foundations [DONE 2026-06-29]: project-local env (rv for R + uv .venv for Python), shared
  helpers (io / design+contrasts / plot theme), data load + QC sanity, 2x2
  factorial + 5 contrasts, concrete quality gate.
- P1 snRNAseq microglia core [ACTIVE]: reprocess (SCT) + substates (homeostatic / DAM /
  IFN / proliferative, UCell), composition (sccomp), pseudobulk DE across contrasts -> the
  amyloid->DAM activation programme. Static interaction = v1-prior null (reported w/ power stmt) -> rate effect to P2.
- P2 Interaction trajectory: activation pseudotime (homeostatic->DAM); test amyloid
  advance + tau x amyloid progression synergy, decomposed composition vs progression.
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
  power null even single-cell) -> the synergy is a progression-RATE effect, deferred to P2 (trajectory); P1
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
