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

## Active plan: `p0_foundations_plan.md` - P0 Foundations (open; **S5 = last step**). Stack: targets + rv
+ uv + project-local Quarto, P3M-pinned. Reports = ONE offline self-contained HTML (standalone `index.qmd`
+ `{{< include _qc.qmd >}}` + `theme.scss` + LOCAL IBM Plex; 9 woff2 COMMITTED, base64-inlined). S1-S4 done.
**NEXT (S5 = close P0):** runnable quality gate `scripts/check.sh` = `rv sync` + `uv sync` + `tar_make()`,
then ENFORCE zero-fault rather than trust `tar_make`'s exit -- assert `tar_meta(fields=c("error",
"warnings"))` all-NA + loop `tests/test_*.R` + grep the Quarto render log for `Warning`/`WARN`; non-zero
exit on any hit (use tar_make's DEFAULT reporter, not the deprecated `reporter="summary"`). Then lock
memory.md's provisional quality gate to this concrete gate, finalise map.md wiring, and close P0 (archive
plan -> `.agent/completed/`, reset Active plan to none). Detail: the plan's S5 step.

## Backlog - phased build (each phase = closeable increments; mine archive_digest per phase)
- P0 Foundations: project-local env (rv for R + uv .venv for Python), shared
  helpers (io / design+contrasts / plot theme), data load + QC sanity, 2x2
  factorial + 5 contrasts, concrete quality gate.
- P1 snRNAseq microglia core: QC, microglia subset + substates (homeostatic / DAM /
  IFN / proliferative), DE (pseudobulk + single-cell) across contrasts -> the
  amyloid->DAM activation programme + the gene/pathway interaction.
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
