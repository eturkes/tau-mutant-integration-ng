# Memory - standing contract

Claude Code repo. Canonical instructions = `CLAUDE.md`; session entry =
`/session-prompt` from `.claude/commands/session-prompt.md`. Use `roadmap.md`
for trajectory/history, `map.md` for live wiring, `history.md` for decision
digests, `archive_digest.md`/branch `archive` for v1 mining.

## Current Scope

Goal: integrate snRNAseq + GeoMx spatial + the 24M TiO2 phospho assay at
protein-group-sum and phosphosite reporting levels across 4 AD mouse genotypes.
Design = 2x2 tau x amyloid:
- MAPTKI = wild-type humanized tau (human MAPT knock-in; NOT tau-KO)
- P301S = mutant humanized tau (base-edited MAPT^P301S;Int10+3)
- NLGF_MAPTKI = amyloid (App^NL-G-F) + WT humanized tau
- NLGF_P301S = amyloid (App^NL-G-F) + mutant humanized tau

Canonical interaction = `(NLGF_P301S - P301S) - (NLGF_MAPTKI - MAPTKI)`.

Construct provenance (P7.1; G1 literature-verified 2026-07-18) — the `MAPTKI`/`P301S` factor
TOKENS are stable identifiers kept unchanged across code; only the described MEANING was
corrected from the earlier wrong "tau-KO baseline":
- MAPTKI = Saito humanized wild-type MAPT knock-in: endogenous murine Mapt replaced with WT
  human genomic MAPT, so mice express all six 3R/4R human tau isoforms under native regulation
  and lack MOUSE tau but are NOT tau-null (a true Mapt-KO tau-negative genotype is not in this study).
- P301S = Watamura base-edited MAPT^P301S;Int10+3G>A allele ("+3" = the Int10+3 splice mutation,
  NOT an age label; 24M/30M are separate age tokens); also a humanized-MAPT knock-in.
- Four groups form a clean 2x2 with NO tau-null arm; the tau factor compares WT humanized tau
  vs mutant (P301S;Int10+3) humanized tau, optionally stratified by App^NL-G-F amyloid.
- Sources: Saito 2014 (App^NL-G-F KI); Saito 2019 (WT humanized MAPT-KI, App/MAPT-KI double KI);
  Watamura 2025 (seven base-edited human-MAPT alleles; MAPT-KI is not Mapt-KO); Benzow 2024
  (separate H1/H2 MAPT-GR series, distinct from the Saito MAPT-KI); Morito 2025 (P301S;Int10+3
  with App^NL-G-F).
- STANDING external-record request (user-supplied; not in-repo): breeding/allele records and
  per-animal provenance (MAPTKI/P301S animal IDs, promoter/locus, genotyping confirmation) to
  confirm the literature-inferred construct at the animal level.

Live report scope (2026-07-16): 10 visible figures, 4 included qmd fragments, 34 report-ancestry
targets (37 live: + `occupancy_harness_check`, `occupancy_robustness`, and `integration_substrate`
non-report leaves, outside report ancestry).
Rendered HTML artifact = `report/tau-mutant-integration.html`; the report directory is pruned after each
render so that HTML is the only user-facing output. Browser/tab title =
`Tau Mutant Integration`. Visible surface = simple numbered figure headings (`Figure 1` ...
`Figure 10`) + figures + per-figure folded code only; no visible document title, TOC,
captions, body prose, tables, or global code-tools menu. Folded-code summaries and expanded
code blocks are intentionally compact via `assets/theme.scss`.
Retired infrastructure remains absent: committed tests, Python/uv files,
composition/sccomp/CmdStan target, P1 per-subpopulation DE target, prose inventory,
stageR layer, mechanism/crossmodality/qc/story chapters and modules. Retained
non-snRNAseq modality-native set = GeoMx sample heatmap (former Figure 10) + one bulk context plate
combining the TiO2 phospho protein-group-sum PCA and phosphosite heatmap; the historical Proteome/Phospho
volcano plots are removed from the live report. The other GeoMx exploratory/native panels are historical only. Historical
claims remain in git + `roadmap.md`; do not treat them as live pipeline contracts.

P6 state decomposition closed 2026-07-14 and is report-integrated. Compact
`microglia_state_substrate` contains aligned Homeostatic/DAM raw-count pseudobulks,
16-unit state/coverage/library metadata, five-programme raw-UCell unit/state means +
pooled SDs, and exact feature/marker maps. `microglia_state_response` contains only
occupancy/state/delta/rotation/bridge inference tables + diagnostics.
`microglia_state_decomposition` contains unit-level standardized score channels,
ordinary/weighted inference, compact S2 evidence, and the fixed classifier; it is now fully
report-disconnected, and no score aggregate reaches Figure 10.
`microglia_state_gene_atlas` jointly fits 32 state pseudobulks with paired-unit
`edgeR::voomLmFit`, sample-quality weights, robust limma eBayes/treat, nine explicit gene
contrasts (four responses, two paired state differences, three interactions), and joint
four-response/two-interaction moderated-F families. `state_decomposition_figures` is the
1.19 MB Figure 10 leaf: accepted raw-unit DAM occupancy, 52 unique declared genes x four
response effects with 95% CI, and 14,438 genes x two paired state-difference backgrounds.
Figure 10 renders direct DAM-minus-Homeostatic response maps plus 48 shared-scale gene-level
factorial small multiples; six direct minimum-effect marker/background rows carry labels +
CI, B2m is display-collapsed, and four below-filter genes are listed. The two-tier plate is
18.8 x 12.2 inches; marker selection remains fixed and programme-free.
`sections/state-decomposition.qmd` loads only this leaf after the stable Figures 1-9;
no Seurat/S4 parent or fitted model crosses the report boundary.

## Data

`storage/data` is a symlink to external host data; read-only for this repo, gitignored.
Use loaders instead of ad hoc parsing.

Raw data facts:
- `snrnaseq.rds`: Seurat, RNA 33683 ENSMUSG genes, full object ~8.3G. Loader keeps
  broad_annotations == Microglia and drops SCT/reductions before reprocessing.
- snRNAseq metadata: genotype, batch, genotype_batch (16 fully crossed units), sex,
  QC metrics, doublets.
- `geomx.rds`: GeoMx WTA Seurat, 19963 genes x 91 AOIs; genotype canonical; spatial
  design cols include bio_rep/tech_rep/slide_rep/roi/segment/SampleID.
- GeoMx source gap: 91 of ideal 112 AOIs → 21 missing across 7
  (`genotype`, `bio_rep`) cells; only `NLGF_MAPTKI` is complete (28/28), while
  `MAPTKI` `bio_rep = 1` is wholly absent. `bio_rep == slide_rep`; one segment.
  Current loaders/descriptors exclude no AOIs; the gap predates this repo's live GeoMx
  model/report path. Full per-cell accounting → `.agent/p7_g3_provenance_dossier.md`.
- The `proteomics_*.tsv` input is itself a TiO2 phospho-PTM Spectronaut export: it
  carries phospho-PTM annotation columns plus 16 `.raw.PTM.Quantity` columns and is
  summed by `PG.ProteinGroups` for the historical `proteome_*` layer. Both that
  protein-group-sum layer and the phosphosite `phosphoproteomics_*.tsv` layer come
  from the same TiO2 `Naoto-Hippo_TiO2_DIA` acquisition; the report uses the 24M
  16-run subset from `proteomics_sample_key.csv`. NO global-proteome file exists in
  `storage/data/`; `proteome_*` code tokens are historical stable names, not a claim
  of an independent global-proteome assay.
- Bulk provenance contract: report bulk = 16/67 manifest runs, exactly the balanced
  Naoto 24M 2x2; Naoto 30M + all Set6 arms are unused. Full manifest/disposition →
  `.agent/p7_g3_provenance_dossier.md`.
- Cross-modality gate: no in-repo animal/aliquot crosswalk. Paired/animal-level claims
  (`DIABLO`, `MOFA`, mediation, within-animal correlation) are PROHIBITED until the user
  supplies a validated cross-assay crosswalk; gene-symbol harmonization/meta-comparison
  remains allowed. Current report is compliant.

## Standing external-record requests (user-only; do not block)

Items the report cannot self-supply; itemized in the G-dossiers, banked here so they survive
plan archival:
- G1: breeding/allele confirmation records for the humanized-tau lines.
- G2: vendor TiO2 enrichment template + phospho acquisition method.
- G3: animal/aliquot cross-assay crosswalk; per-animal survival/attrition, litter, cage, sex,
  batch; GeoMx raw files (DCC/PKC), neg-probe/scan/mask metadata, and the missing-AOI disposition.

## Scientific Spine

Durable headline:
- Amyloid drives microglia homeostatic -> DAM activation.
- Mutant tau modulates amyloid response mostly through DAM-cell composition, not a
  supported progression-beyond-composition shift along the activation trajectory. The
  DAM-composition interaction (+0.174) is directionally robust (positive in all 35 P7.5
  variants) but its zero-null significance is resolution-fragile (P7.5 verdict = FRAGILE):
  read this as a robust-sign, significance-caveated claim.
- Non-snRNAseq modalities provide context figures, not resurrected mechanism or
  cross-modality chapters.

Subpopulations: Homeostatic, DAM, IFN, Proliferative. Current coherent clusters contain
Homeostatic/DAM/IFN; no Proliferative-dominant cluster in the built annotated object.

## Analysis Contracts

Replicate unit = `genotype_batch` (16 ids, 4/genotype). Batch belongs in snRNAseq
designs when estimable. GeoMx uses slide/bio-unit structure; 24M bulk omits batch.

Five canonical contrasts everywhere:
- `tau_alone`
- `nlgf_in_maptki`
- `nlgf_in_p301s`
- `tau_in_nlgf`
- `interaction`

Cross-modality integration (P8):
- `integration_substrate` is the symbol x five-contrast effect-size substrate over exactly three
  modalities: snRNAseq, GeoMx, and bulk protein-group. It stores raw `logFC`/moderated-`t` plus
  invertible per-modality/per-contrast median-MAD robust-z matrices.
- Symbol coverage is complete-case 3,109 / >=2-modality 12,427 / union 22,241; pairwise overlaps
  are snRNAseq-GeoMx 12,324, snRNAseq-bulk 3,132, and GeoMx-bulk 3,189.
- The 3,019-symbol phosphosite parent-gene collapse is a within-assay alternate of the SAME TiO2
  assay, NOT a modality; bulk remains one modality and no per-animal pairing is attempted.
- `integration_substrate` is a compact, parent-isolated, self-validating non-report leaf rebuilt by
  `scripts/check.sh`; exact source reconstruction, invertibility, deterministic representatives,
  fixed counts, and the <25 MiB size ceiling are runtime-fatal.

Microglia reprocess:
- SCT-v2/glmGamPoi on RNA counts, regress percent_mt + percent_contam.
- Harmony integrates batch only; genotype/amyloid are biology, not correction terms.
- Louvain clusters at multiple resolutions; primary stable column = `microglia_clusters`.
- Re-run is seed-deterministic up to tolerance, not bitwise.

Microglia annotate:
- UCell scoring on SCT data; cluster-level subpopulation assignment is authoritative.
- Contaminant pruning uses raw identity-vs-contaminant evidence, not z-scaled subpopulation
  scores.
- Prune asymmetry exists: dropped low-quality clusters are enriched for NLGF_MAPTKI but
  are not DAM-high. Report caveat factually.

Microglia composition:
- No live composition inference target. Replicate-unit stacked bars are derived inside
  `microglia_report_data()` from annotated Seurat metadata and emitted as
  `unit_composition`.

Microglia state decomposition (P6 closed):
- Primary universe = existing cluster labels Homeostatic + DAM; S1 coverage =
  22,363/23,160 cells (96.56%), >=93.93% in every unit; both states have all 16
  units and >=31 cells/unit.
- `microglia_state_substrate` uses raw RNA counts only, keeps both 33,683 x 16
  state matrices column-aligned to `unit_meta`, and records all five fixed marker
  programmes without collapsing feature rows.
- `microglia_state_response` uses 16 units/9 residual df throughout. Live genes =
  13,599 Homeostatic / 9,148 DAM / 9,123 paired; direct response uses harmonic voom
  precision and carries an unweighted sensitivity. Interaction DAM occupancy = +0.174
  fraction (95% CI 0.095-0.253): zero-null supported, 0.10 minimum-effect FDR = 0.081
  (unresolved at 5%); empirical-logit batch-stratified permutation p = 0.021.
  Fixed raw-count rotations support an interaction Homeostatic programme within
  Homeostatic cells, not a DAM or direct state-difference programme at FDR <=0.05.
  S3 exact UCell channels reconstruct at <=1.24e-15 and use ordinary 9-df OLS with
  0.25 pooled-SD minimum/equivalence margins; cell-count/harmonic-count WLS is sensitivity.
  Integrated interaction verdict = unresolved. All five composition-score channels are
  equivalent within the margin, but occupancy minimum-effect FDR = 0.081; no within-DAM
  programme has meaningful UCell + same-direction significant rotation support, and all-five
  within-DAM equivalence is unearned. `DAM-selective` requires concordant meaningful direct
  DAM-minus-Homeostatic evidence in both score and raw-count programme representations.
- Pooled Homeostatic+DAM versus whole-MG effect correlations = 0.982-0.994 across
  contrasts; descriptive exclusion bridge only, with no agreement threshold.
- Paired multivariate gene atlas = 14,438 filter-passing genes, within-unit state
  correlation 0.183, sample weights 0.526-1.431, and residual df 2-18 after zero-aware
  adjustment. Joint four-response moderated F supports 1,120 genes at BH 5%; joint
  Homeostatic/DAM interaction F supports zero. The combined-model interaction effects agree
  with established separate-state fits (Pearson 0.992/0.993; median absolute delta
  0.037/0.044 log2FC for Homeostatic/DAM).

DAM-occupancy robustness (P7.5 executed):
- `R/analysis/occupancy_harness.R` maps membership to the occupancy family: E1 beta-binomial
  probability standardization (primary), E2 empirical-logit OLS/permutation, E3 raw-proportion OLS.
- Frozen preregistration = `.agent/p7_dam_occupancy_prereg.md`; execution dossier =
  `.agent/p7_dam_occupancy_robustness_results.md`.
- `occupancy_robustness <- run_occupancy_robustness(microglia_processed, microglia_annotated,
  microglia_state_substrate, symbol_map)` is the compact, parent-isolated non-report leaf: exact
  reference substrate/membership/re-annotation anchors; 35-row reference + 34 one-at-a-time variant
  table; 105 E1/E2/E3 attempts; E1 range/tipping/margin summaries; failure inventory; sign concordance;
  frozen verdict. Live payload is compact, orders under the 256 KiB audit budget;
  `audit$serialized_bytes` (~5,769) self-measures before that field is populated, so it is a
  pre-write estimate a few bytes off the authoritative on-disk `tar_meta(bytes)`, not an exact size.
- Pre-committed verdict = **FRAGILE**: all E1 estimates remain positive, but resolutions 0.5 and 0.6
  cross to `fdr_zero > 0.05`; E1 empirical range = [0.1090994, 0.2306348]. No estimator failed;
  E2/E3 signs track E1 in 35/35 variants. E1 `fdr_minimum` spans [7.369534e-05, 1.0], resolves in
  four variants, and no variant has `abs(estimate) <= 0.10`.
- Frozen-rule details: every variant, including an E1 failure, must pass for ROBUST-POSITIVE;
  direction tipping is `estimate <= 0`; margin reporting includes both `fdr_minimum <= 0.05` and
  `abs(estimate) <= 0.10`. Optional UCell re-scoring was not run (declared optional secondary).
- `occupancy_harness_check` reproduces the established +0.174 current-label family exactly and
  exercises reduced designs only with fabricated counts. Its Layer-A round-trip includes non-primary
  filler rows so `n_retained` and coverage reproduce exactly; E2/E3 success requires exact schemas
  and finite inferential output, and reduced-design smoke requires E2/E3 `ok` plus a separate fabricated
  non-estimable fixture with all three failures recorded. `scripts/check.sh` invalidates and rebuilds
  this guardrail together with `integration_substrate`; `occupancy_robustness` stays an explicitly
  built non-report leaf.

Microglia DE:
- Live target = `pb_de_microglia` only.
- Raw RNA pseudobulk by genotype_batch, limma voom with quality weights, topTables for
  the four-panel amyloid-response scatter.
- No stageR, per-subpopulation DE, MDE, or DAM-direction helper in the live DAG.

Trajectory:
- Slingshot on Harmony embedding with forced Homeostatic -> DAM lineage.
- IFN/Proliferative are omitted from the lineage frame via flags, not deleted from the
  biological record.
- Primary inference uses 16-unit pseudotime summaries plus Kitagawa-style decomposition;
  glmmTMB per-cell model is supportive and degrades to recorded failure instead of
  blocking the report.

Modality context:
- `R/analysis/modality_de.R` restores only primary DE needed for report figures:
  GeoMx voom/TMM with slide fixed effect + duplicateCorrelation plus the retained sample
  heatmap descriptor; the historical `proteome` layer is the TiO2 phospho-PTM report
  summed to `PG.ProteinGroups`, while `phospho` is the phosphosite view of the same
  `Naoto-Hippo_TiO2_DIA` acquisition. Both use limma-trend on log2 median-normalized
  24M intensities; no independent global-proteome assay exists.
- The 24M run order is genotype-blocked (01-04 MAPTKI, 05-08 P301S, 09-12
  NLGF_MAPTKI, 13-16 NLGF_P301S), so between-genotype effects cannot be separated
  from acquisition order/batch. Primary bulk figures continue to use the no-batch
  `$top` fits. Each bulk target also stores compact `$run_order_sensitivity` from a
  rank-5 design with mean-centered continuous `run_index` (11 residual df); it captures
  only within-acquisition linear drift, does not fix the alias, and keeps the interaction
  aggregate shift approximately zero by design. The bulk layers are context-only.
- `geomx_de$sample_heatmap` is descriptive only: the retained GeoMx modality-native figure is a compact AOI track
  atlas with AOI columns average-linkage clustered by the displayed first five DAM genes from the prior full-row order,
  then dendrogram-rotated by mean displayed DAM z-score
  (`B2m`, `Apoe`, `Ctsd`, `Tyrobp`, `Trem2`),
  with genotype/tau/amyloid context above and top legends for genotype plus shared tau/amyloid no-versus-yes colors.
  Spatial/QC, bio/slide replicate, tech-replicate, ROI, signature, and non-DAM marker tracks are omitted; ROI exactly encodes
  genotype block + tech_rep. It excludes no AOIs and changes no DE model.
- Four-panel amyloid-response scatter uses one shared off-diagonal feature cutoff:
  `|x-y| >= 3.5`. Figure 8 labels all points past the cutoff and draws all four facets on one shared square coordinate
  range so the dotted cutoff bands have the same visual distance from the identity line; its line
  legend labels the dotted cutoff bands as `threshold: |x-y| >= 3.5 log2FC`. Figure 9
  scores all shared-cutoff selected features and displays priority-ordered GO term-family
  role/fallback buckets with aggregate `|P301S - MAPTKI| >= 0.5`; complement/MHC, phagocytosis,
  and chemotaxis are separate buckets, lipid/endolysosome/synaptic buckets precede broader
  cell-cell-adhesion, extracellular-matrix, motility/cytoskeleton, and broad immune/inflammatory residual buckets,
  predicted/unannotated + other-annotated no-role buckets
  are excluded, and each visible category label lists every retained scored feature in that category.
  Current visible Figure 9 facets = snRNAseq + bulk phospho (protein-group) + bulk phospho (site);
  GeoMx has no retained categorized group rows under the live shared-cutoff/filter rule.
- Bulk context plate combines the TiO2 phospho protein-group-sum sample PCA (historical
  `Proteome` payload key) with the same assay's phosphosite heatmap. The heatmap keeps
  20 top mutant-tau amyloid phosphosite rows after excluding
  parent genes `Plcb1` and `Arhgef7`, keeping the same effect direction as the top-ranked candidate,
  and silently collapsing exact duplicate log2 median-normalized profiles to the first ranked
  representative.
- Retired GeoMx QC/normalization/ordination/gene-detection/spatial-program/contrast/ROI/decon
  figures are ledger history, not live report/path contracts.
- Auxiliary SpatialDecon beta/abundance and broad mechanism/cross-modality target families
  stay deleted. Bulk run-order sensitivity is now a compact field inside the two primary DE
  targets only; it is not a target family and is not wired into any figure.

Report:
- `sections/{microglia,trajectory,modality,state-decomposition}.qmd`;
  rendered artifact =
  `report/tau-mutant-integration.html`; visible HTML = simple numbered figure headings +
  figures plus per-chunk folded code. Alt text stays in image attributes; captions/TOC/visible
  title stay absent.
- qmd chunks set `options(warn=2)`. Any warning during render is a failure.
- Compact report targets must keep qmd renders off the 612MB annotated Seurat object.
- Legend hygiene: `geom_text` / `ggrepel` layers that inherit mapped aesthetics must use
  `show.legend = FALSE` unless text glyphs are intentionally part of the key; otherwise ggplot
  can add the default `a` text glyph beside point keys.

## Environment And Gate

Live stack = rv-managed R + project-local Quarto + targets. No committed Python/uv
surface unless a future report-producing step earns it.

Iteration mode = rapid + failure-tolerant. User owns figure inspection. Claude Code skips
Chromium/PDF/PNG visual QA and optional test/check gates; run only the minimum command
needed to produce a requested artifact or diagnose a concrete failure.

Fresh bootstrap:
1. `scripts/bootstrap/sysdeps.sh`
2. `scripts/bootstrap/rv.sh`
3. `scripts/bootstrap/quarto.sh`
4. `rv sync`
5. `scripts/check.sh`

`scripts/check.sh`:
- invalidates and rebuilds `report` plus the non-report `occupancy_harness_check` and
  `integration_substrate` validation leaves
- does not sync the environment or scan all target metadata/logs
- relies on qmd `options(warn=2)`, target/harness failures, and `render_report()`'s self-contained/pruned HTML assertion
- functional green does not currently imply byte stability: repeated renders vary only in Figure 8
  `plot-modality-amyloid-effect`; `modality_interaction_scatter()` fixes ggrepel seed 42 but also sets
  wall-clock `max.time=3`, so the time-bounded label layout changes PNG pixels across runs. P7.5
  observed 10 figures throughout but distinct hashes/sizes on consecutive clean renders; fixing this
  requires a separately authorized `R/report/plot.R` edit (P7.5 explicitly leaves report code untouched).

Run `rv sync` manually after dependency changes; use `scripts/check.sh` for fast local report iteration.

## Maintenance

Keep `.gitignore` aligned with generated artefacts. Current deleted infrastructure should
stay absent unless it directly accelerates or protects the final report path.

Layout contract: `.agent/` = tracked project state; `.claude/` = shared Claude Code
settings + commands (local/runtime files remain ignored); `.serena/` = tracked project
LSP config. Report fragments live in `sections/`; R code is grouped under
`R/{core,analysis,report}/`; generated pipeline/QA state lives in `storage/{targets,qa}/`.

When adding a dependency, prefer current project-local R/Quarto path first; update
`rproject.toml` + `rv.lock` together. Avoid global library/tool leakage.

Project state sources:
- durable project fact -> `.agent/memory.md`
- live wiring -> `.agent/map.md`
- trajectory/ledger -> `.agent/roadmap.md`
- workflow entry -> `.claude/commands/session-prompt.md`
- LSP languages/read exclusions -> `.serena/project.yml`
