# P8 — Cross-modality symbol/effect-size integration (plan)

Posture: BUILD. Opened after P7 REVIEWED (2026-07-19). Direction = "cross-modality paired
integration" (user-selected at the direction gate); constrained to the VALID substitute the user
then confirmed: **symbol/effect-size joint integration WITHOUT a crosswalk**.

Hard fact settled at planning: per-animal paired methods (DIABLO / MOFA-on-matched-samples /
mediation / within-animal correlation) are INFEASIBLE here — the only identifier shared across
snRNAseq (16 units), GeoMx (91 AOIs), bulk TiO2 (16 runs) is `genotype` (4 levels); no
animal/aliquot crosswalk exists. They stay PROHIBITED by the P7.3 gate and are NOT attempted.
Integration operates ONLY on the shared gene-symbol × 5-contrast EFFECT-SIZE space (P7.3-compliant:
symbol-level harmonization is the allowed side of the gate).

## Answer being built
Quantify how the amyloid×tau response is SHARED vs MODALITY-SPECIFIC across the three assays, at
gene and pathway resolution, honestly (signal is modest by construction). Extends spine finding #1
(amyloid→DAM) cross-modally; independently corroborates finding #2 (the tau×amyloid interaction is
mostly microglia-composition-specific → expected near-zero cross-modality concordance).

## Precondition — MET (real-input evidence, this planning session)
Confirmed by reading the cached DE targets (`Rscript` + `targets::tar_read_raw`, store
`storage/targets`), no fabricated pairing:
- Four producers emit 5 per-contrast topTables (`tau_alone, nlgf_in_maptki, nlgf_in_p301s,
  tau_in_nlgf, interaction`; `R/core/design.R`) over the full retained feature set, each with
  `logFC` + moderated `t` + `P.Value`/`adj.P.Val` and a symbol column:
  - snRNAseq `pb_de_microglia$top[[c]]` — 14,512 genes; symbol = `gene` (ENSEMBL) → map via
    `symbol_map` (ENSMUSG→symbol, 33,683 rows, 0% NA on this set).
  - GeoMx `geomx_de$primary$top[[c]]` — 19,959 genes; symbol = `symbol` (direct).
  - bulk-PG `proteome_de_24m$top[[c]]` — 3,379 protein groups; symbol = `gene_first`
    (multi-symbol list in `gene_symbols`).
  - bulk-site `phospho_de_24m$top[[c]]` — 17,707 sites; parent symbol = `gene`; `site_id`=`Gene_Ssite`.
- Shared symbols: snR∩GeoMx = 12,324; **complete-case (all 3 modalities) = 3,109**;
  **≥2-modality = 12,427**; total distinct = 22,241.
- Scale: per-contrast logFC sd — snRNAseq 0.29–0.48, GeoMx 0.16–0.40, **bulk 0.73–1.12** (~2–4×).
  This scale gap matters ONLY for the variance-based DECOMPOSITION (else bulk dominates joint
  variance) → per-modality-per-contrast standardization there. It does NOT matter for CONCORDANCE:
  Spearman is rank-based and Pearson is scale-invariant, so the concordance layer uses RAW logFC.
- Concordance preview (Spearman on 3,109): amyloid contrasts strongest (`nlgf_in_p301s`
  snR↔GeoMx 0.19; `nlgf_in_maptki` 0.13), `interaction` near-zero everywhere (0.07/−0.06/−0.05) →
  structured but modest signal, matching the spine.
- Joint-SVD preview (standardized concat, 3,109×15): real structure (PC1 = 22% var, interpretable
  amyloid-weighted loadings) BUT concatenated PCA MIXES shared + modality-specific variation (PC1
  signs conflict across assays: GeoMx interaction −0.27 vs bulk +0.36) → a principled
  joint-vs-individual method is required, not naive concatenation.

Precondition MET. No gated external input (the whole point of the symbol-level route). All units OPEN.

## Integrity constraint (load-bearing)
- Bulk `proteome_de_24m` (protein-group sum) and `phospho_de_24m` (phosphosite) are ONE TiO2 assay
  at two aggregation levels (identical `Naoto-Hippo_TiO2_DIA` 16-run meta; G2 dossier). Integration
  treats **bulk = one modality** (protein-group gene-level = primary bulk view). Phosphosite
  parent-gene collapse is an OPTIONAL within-assay alternate, NEVER a fourth "modality" → no
  double-counting.
- Three integration modalities: `snRNAseq`, `GeoMx`, `bulk` (protein-group).

## Honesty posture (load-bearing — set by the method research)
Effect-size integration here operates on per-contrast POINT ESTIMATES over UNMATCHED modalities; it
cannot recover replicate-level uncertainty, and genes are NOT exchangeable sampling units (coexpression,
shared samples, gene-specific SE/coverage). Consequences the whole milestone honours:
- **All three layers are DESCRIPTIVE** — the decomposition, the concordance network, and the pathway
  consensus. A cross-gene label permutation for a cross-modality concordance/enrichment p is
  ANTI-CONSERVATIVE and INVALID (Goeman & Bühlmann 2007; the CAMERA inter-gene-correlation point);
  from cached topTables alone NO design-valid calibrated cross-modality p exists. Report effect sizes,
  ρ, and variance-explained with explicit n + coverage; never a gene-permutation p.
- The ONLY design-valid calibration is **per-modality resampling of BIOLOGICAL UNITS + DE refit +
  bootstrap ρ** (each modality's own units — snRNAseq pseudobulk units, GeoMx bio-units, bulk runs —
  resampled, its DE refit, ρ recomputed → a bootstrap CI). Feasible from the cached per-unit matrices
  (NOT the 8GB Seurat), but heavier than topTables-only; it is an OPTIONAL P8.3 sensitivity
  (feasibility-gated), with the descriptive ρ as the honest floor.
- Primary statistic = **logFC** (the effect) for BOTH concordance (Spearman primary, Pearson
  sensitivity) and the decomposition (standardized). **moderated-t is NOT an effect** — it conflates
  logFC with assay/gene precision, so it enters only as a SECONDARY "evidence-statistic concordance"
  view, never the primary. Both stored (raw + standardized) in P8.1.
- Signal is modest by construction; a near-null joint component or weak ρ is a VALID reportable
  outcome. No statement asserts strong cross-modality agreement the data does not show.

## Scope
IN: harmonized effect-size substrate → joint-vs-individual decomposition → cross-modality concordance
network → pathway/functional consensus → report figure section.
OUT: per-animal pairing / DIABLO / MOFA / mediation / within-animal correlation; any new external
data; resurrecting the deleted mechanism/cross-modality chapters or their target families;
SpatialDecon; run-order refits (bulk stays context-only per memory).

## Method decisions (locked; research- + real-input-grounded)

### Primary joint decomposition = reimplemented pure-R AJIVE (joint vs individual vs residual)
The scientific question — "how much of the amyloid×tau response is SHARED across assays vs
modality-specific" — requires an explicit joint/individual/residual variance split. No maintained R
package delivers that cleanly for outcome-free, genes-as-shared-dimension blocks:
- MOFA2 (BioC 1.22) → needs Python `mofapy2` via reticulate/basilisk (reintroduces the Python
  surface the lean pass removed) AND its own authors flag D=5 « the D<15 weak-view warning, factors
  are gene embeddings ignoring replicate SE → exploratory-only in its own framing. REJECTED.
- RGCCA/SGCCA (CRAN 3.0.3) → genuinely unsupervised, native blockwise NA, but pulls `caret` (heavy
  transitive) AND maximizes cross-block CORRELATION, not joint/individual separation; sparsity is
  pointless at p=5. Not primary.
- mixOmics 6.36 `block.pls/spls` → REQUIRE an outcome (`Y`/`indY`); DIABLO strictly supervised.
  REJECTED (not outcome-free).
- MCIA `omicade4` 1.52 (pinned Bioc 3.23 snapshot, light: ade4/made4) → maximizes a weighted sum of
  squared block↔consensus COVARIANCES (co-inertia = a ±-sign-ambiguous consensus axis), not
  joint/individual; quiet since 2020. Not the right decomposition.
- Concatenated PCA → SVD preview proved it conflates joint + individual. REJECTED as primary.
- `r_jive` (CRAN r.jive 2.4) IS current (all 13 CRAN checks OK 2026-07-19, no compilation, light
  deps gplots/abind, GPL-3) and does the actual JIVE joint/individual/residual split (perm rank
  selection, row-center + block-Frobenius scaling, `orthIndiv=TRUE`) → adopted as the **packaged
  CROSS-CHECK reference** for P8.2 (validate the reimplementation's variance split against it), NOT a
  committed dep. Caveat: its wrapper passes `ncomp=min(dim)` (full rank) and relies on `SVDmiss` fill
  — moot on our complete-case (no NA) but a reason not to trust its rank behaviour blindly. AJIVE
  (Feng/Jiang/Hannig/Marron, JMVA 2018, arXiv:1704.02060) handles CORRELATED individual subspaces
  better than JIVE — relevant here — but its only package (`py_jive`/`jive`) is DEPRECATED (successor
  `mvdr` unmaintained). → Per CLAUDE.md (reject library-availability-driven choices; reimplement
  compact cores), P8.2 REIMPLEMENTS the AJIVE decomposition in pure R (base `svd`/`qr`), cross-checked
  against r.jive.

Orientation + algorithm (both Explore finders, reconciled):
- Input = **standardized logFC** (per-modality-per-contrast robust-z, so neither bulk's ~2–4× scale
  nor any single contrast dominates the variance split). moderated-t is a documented sensitivity.
- Blocks are oriented **features × shared-objects = 5 contrasts (rows) × genes (columns)**; genes are
  the common objects. "Joint" = the shared **gene-score subspace** (common column space) — a
  gene-direction all three modalities express; each modality's expression of a joint component is its
  own 5-contrast loading vector (this reveals WHICH contrasts, expected amyloid, drive the shared
  axis). "Joint" does NOT mean equal contrast loadings across modalities.
- Per-block SVD → initial signal rank `r_k`. **p=5 hard constraint**: the random-direction/Wedin
  perturbation bound samples an r-dim space orthogonal to an r-dim space (needs both dims ≥ 2r), so
  at 5 features it is defensible only for `r_init ≤ 2`; enforce `r_J + r_{I,k} ≤ 5` and leave ≥1
  residual dimension (no saturation). Deterministic scree/variance threshold; report singular values
  + selection diagnostics.
- Joint rank = # shared directions whose principal angle between block signal row-spaces beats a
  permutation/random-direction bound (fixed RNG seed → deterministic — this is a NUMERICAL rank
  threshold on the row-spaces, NOT an inferential gene-label null). Reconstruct
  `X_k = J_k + A_k + E_k`; report per-block variance-explained for J/A/E, joint gene scores, and each
  block's joint-component 5-contrast loadings. Interaction expected to load individual/residual. This
  is a DESCRIPTIVE variance partition, not a test of "sharedness".
- Fallback if the reimplemented joint-rank proves fragile: fixed-rank joint projection with reported
  robustness, and/or the r.jive cross-check as the fallback engine.

### Cross-modality concordance = pairwise logFC-Spearman network + directional overlap (DESCRIPTIVE)
- Statistic: **Spearman ρ of RAW logFC** (primary; Pearson as a scale-invariant sensitivity), per
  modality-pair × per contrast = the 3×5 grid. Plus a SECONDARY moderated-t Spearman "evidence-
  statistic concordance" view (does the evidence pattern agree, distinct from the effect).
- Universe: the **common 3,109-gene complete-case set is PRIMARY** for edge comparability (all 3
  edges on one gene set → directly comparable ρ); the per-pair ≥2-modality overlaps are a coverage
  SENSITIVITY. NEVER zero-impute a missing gene.
- Inference: report ρ **DESCRIPTIVELY** (magnitude, sign, n) — NO cross-gene label permutation p
  (anti-conservative; genes not exchangeable). OPTIONAL design-valid calibration = the per-modality
  unit-resample + DE-refit bootstrap CI on ρ (honesty posture above; feasibility-gated).
- Directional overlap: reimplement the rank–rank hypergeometric directly (RRHO2's stratified
  four-quadrant geometry is the design reference, but its GitHub v1.0 is dormant since 2020 AND both
  RRHO/RRHO2 sources carry a `phyper` N+1 off-by-one → audit/reimplement, no dep). DESCRIPTIVE map:
  the raw pixel p is uncalibrated (up to N² nested dependent tests, max selected) and design-valid
  resampling is infeasible cross-modality, so the overlap illustrates direction only.
- Expected: amyloid contrasts concordant, interaction near-zero.

### Pathway/functional consensus = `msigdbr` GO-BP + ≥2-modality directional consensus (DESCRIPTIVE)
Default (lean, no new dep): `msigdbr` (already locked) mouse GO-BP (+ optional project
DAM/Homeostatic/IFN/MHC_APC sets). Per modality × contrast, a compact set-level DIRECTIONAL score
(mean standardized effect of set genes, or a preranked statistic). Consensus = sets scoring same-sign
above threshold in ≥2 modalities for a contrast, with an explicit minimum-modality-coverage count
(guards single-modality pass-through). Reported DESCRIPTIVELY: a competitive gene-set (gene-permutation)
null is anti-conservative (CAMERA's whole point), and a valid rotation/CAMERA null needs per-sample
residuals + design (NOT topTables) → no per-modality enrichment p is claimed by default; an optional
CAMERA-style inter-gene-correlation-corrected enrichment is available ONLY if refitting from the
cached per-sample matrices (P8.3 feasibility call).
Surveyed packaged alternatives (documented; P8.3 makes the final reimplement-vs-adopt call — the
compact reimplement is the default):
- **ActivePathways** 2.0.6 (CRAN, no-compile, GPL-3; quiet ~1yr) — closest structural fit: per-gene
  multi-modal P-value matrix → Brown covariance-aware merge (of the MODALITIES, not the genes) →
  hypergeometric enrichment, directional **DPM** (2024) encoding relative up/down and cancelling
  discordant modalities (NA → P 1 / direction 0). Packaged fallback; note its Brown step corrects
  cross-MODALITY dependence, not the cross-GENE non-exchangeability, so its enrichment p inherits the
  same competitive-null caveat.
- **mitch** 1.24 (BioC, no-compile, maintained Jan 2026) — MANOVA over a signed rank matrix, models
  contrast-column covariance, per-axis signed `s`; best for multi-CONTRAST within ONE gene universe,
  less natural for cross-MODALITY consensus (would conflate modality × contrast axes).
- **multiGSEA** 1.22 (BioC, delegates fgsea) — per-layer fgsea + Stouffer/Fisher/Edgington
  p-combination, but DROPS NES sign (needs manual sign gating). Heterogeneous universes.
- (`fgsea` re-add — cut in the lean pass — remains an option for a battle-tested preranked core.)

### Guardrails threaded through every unit
- Symbol-level only; no pairing; P7.3-compliant. Bulk = one modality (no double-count).
- logFC (effect) primary everywhere; moderated-t only as a secondary evidence-concordance view.
  Standardization is a DECOMPOSITION-variance device only; concordance uses raw logFC (scale-invariant).
- Deterministic outputs (fixed seed for the AJIVE numerical rank threshold + any optional bootstrap);
  compact targets read the CACHED DE topTables (no 8GB Seurat, no 612MB annotated); the optional ρ
  bootstrap may reach the cached per-unit matrices only. Oracle-tested (exact reconstruction from raw
  DE targets; decomposition identity + row-space orthogonality; synthetic-rank recovery fixtures).
- **All three layers DESCRIPTIVE**: no design-valid calibrated cross-modality p from topTables
  (gene-exchangeability fails); report variance-explained / ρ / directional scores with explicit n +
  coverage; interaction stays "mostly modality-specific"; a near-null shared component is valid.

## Units (5; sequential; all gate-independent — the original P8.3 "concordance + pathway" was SIZE-CHECK-split at the concordance/pathway seam, 2026-07-20; live DONE/OPEN status in `.agent/roadmap.md`)
New module `R/analysis/integration.R`; new fragment `sections/integration.qmd`; new-target prefix
`integration_` (distinct from the deleted `crossmodality_*` chapter). Report is figure-only; the new
section appends AFTER `sections/modality.qmd`.

- **P8.1 — Harmonized effect-size substrate.** `integration_substrate` (compact `qs` target) +
  `R/analysis/integration.R` core + a focused project-local oracle test (re-adding a committed test
  for the new analysis module, matching the P6 units; the lean pass dropped committed tests but the
  oracle-heavy math warrants one). Harmonize the 3 modalities to gene symbols (snRNAseq
  ENSEMBL→symbol via `symbol_map`, dup-symbol collapse by max `AveExpr`; GeoMx direct; bulk
  `gene_first`, dup collapse), emit per-modality **[5 contrasts × genes]** logFC (primary) AND
  moderated-t (secondary) matrices — RAW plus per-modality-per-contrast standardized (robust z; the
  standardized copy is consumed by the decomposition, concordance uses raw) — the complete-case
  (3,109) and ≥2-modality (12,427) index sets with per-pair overlaps, phosphosite parent alternate
  recorded, and full provenance/coverage.
  Acceptance: exact reconstruction of each modality's per-contrast logFC/t from the raw DE targets
  (tol 0); shared-set counts reproduce 3,109 / 12,427; standardization invertible; bulk single-modality;
  target compact + parent-isolated; no report change; gate green.

- **P8.2 — Joint-vs-individual decomposition.** `integration_decomposition` (compact target) +
  AJIVE reimplementation in `R/analysis/integration.R` + oracle tests + r.jive cross-check. Reads the
  `integration_substrate` complete-case **standardized logFC** (5×genes) blocks (moderated-t = a
  sensitivity). Emit joint/individual/residual variance per block, joint gene scores, per-block joint
  contrast-loadings, chosen ranks (`r_init ≤ 2`, `r_J+r_{I,k} ≤ 5`, ≥1 residual) + selection diagnostics.
  Acceptance: `X_k = J_k + A_k + E_k` to <1e-8; joint orthogonal to individual row-spaces; deterministic
  ranks under fixed seed; synthetic fixture with planted joint rank recovers it; variance split agrees
  with an r.jive run on the same blocks (within a documented method tolerance); compact + parent-isolated;
  no report change; gate green. MAIN SIZE-CHECK before dispatch (decomposition math + fixtures + cross-check).

- **P8.3 — Concordance network.** `integration_concordance` (compact target) + tests. Pairwise
  **logFC** Spearman on the common 3,109-universe (primary; Pearson + per-pair ≥2-modality overlap +
  moderated-t = sensitivities), 3 modality-pairs × 5 contrasts, reported DESCRIPTIVELY — NO
  gene-permutation p; OPTIONAL per-modality unit-resample + DE-refit bootstrap CI on ρ
  (feasibility-gated — DEFER with recorded rationale if it cannot reach cached per-unit matrices
  cheaply) — plus the reimplemented RRHO-style directional-overlap (phyper off-by-one audited out,
  descriptive). Reads ONLY `integration_substrate` (raw logFC/t + index sets); the optional bootstrap
  may reach cached per-unit matrices only. Acceptance: ρ + overlap counts reproduce on the real
  substrate (tol 0); any bootstrap/seeded output deterministic under a fixed seed; overlap map audited
  free of the phyper off-by-one; compact + parent-isolated; no report change; gate green.

- **P8.4 — Pathway consensus.** `integration_pathway` (compact target) + tests. `msigdbr` mouse GO-BP
  (+ optional project DAM/Homeostatic/IFN/MHC_APC sets), per modality × contrast a compact set-level
  DIRECTIONAL score (mean standardized effect over set genes, coverage-gated) + ≥2-modality
  coverage-gated DESCRIPTIVE consensus (same-sign above threshold in ≥2 modalities per contrast, with
  an explicit minimum-modality-coverage count guarding single-modality pass-through). The
  reimplement-vs-ActivePathways decision is made here — DEFAULT = lean reimplement over the
  already-locked `msigdbr` (no new dep, per CLAUDE.md); ActivePathways / fgsea re-add only if the
  reimplement proves inadequate, with recorded rationale. Competitive-null gene-permutation p is
  anti-conservative → NO calibrated enrichment p; report scores + coverage descriptively. Reads ONLY
  `integration_substrate` (standardized effect for scoring) + `msigdbr`. Acceptance: set scores +
  consensus membership reproduce on the real substrate (tol 0); consensus deterministic under a fixed
  seed; coverage gating correct; compact + parent-isolated; no report change; gate green.

- **P8.5 — Report integration.** `integration_figures` (compact leaf) + `sections/integration.qmd`
  (append after modality in `index.qmd`) + wire into `report`/`report_sources`/`render_report()`
  args + memory/map update. New numbered figures (Figure 11+): (a) joint/individual/residual variance
  + joint contrast-loading heatmap + top joint genes (framed descriptive/exploratory); (b) per-contrast
  cross-modality concordance heatmap (logFC Spearman, 3 pairs × 5 contrasts, descriptive ρ) +
  directional-overlap; (c) pathway consensus (concordant GO-BP for the amyloid contrasts). Acceptance:
  `scripts/check.sh` green (warn=2, self-contained pruned HTML); new figures render with fold +
  `fig-alt`; no disturbance to Figures 1–10; DAG count updated in memory/map. Extends the curated
  report from 10 to ~12–13 figures (intentional scope growth per user direction).

Sizing: each unit is a focused module (~200–500 lines + oracle test) over CACHED compact DE
topTables (no heavy Seurat load; the optional ρ bootstrap adds cached per-unit matrices) → fits one
compaction-free 272K Agent window with the ~72K reserve (P6 comparable units peaked ~119–193K impl).
`main=`/`impl=` recorded per unit at close.

## Risks / watch-items
- Reimplemented AJIVE joint-rank at p=5 is the trickiest piece → P8.2 pins it with a synthetic
  planted-rank fixture + a fixed seed + `r_init ≤ 2` / ≥1-residual constraints, AND cross-checks the
  variance split against r.jive; fallback = fixed-rank joint projection with reported sensitivity (or
  r.jive as the engine) if the angle threshold is fragile.
- Only 5 contrast columns per block → low per-block rank (≤4, effectively ≤2 joint); keep ranks small,
  report selection diagnostics; sparsity/SGCCA unnecessary at p=5.
- INFERENCE HONESTY (load-bearing): the integration is DESCRIPTIVE — no calibrated cross-modality p is
  claimed from cached topTables (a cross-gene permutation is anti-conservative; genes not
  exchangeable). The only design-valid calibration is the optional per-modality unit-resample + DE-refit
  bootstrap on ρ. logFC (effect) is primary; moderated-t conflates effect + precision → secondary
  evidence-concordance view only. Every cross-modality statement carries its variance/ρ + n/coverage.
- Re-adding a committed test file reverses the lean-pass "no committed tests" posture for the new
  module only; justified by the oracle-heavy math. Keep tests project-local + runnable; the gate stays
  the report-render path unless a unit wires tests in.

## Closeout expectations
Per unit: MAIN restates scope, SIZE-CHECKs, dispatches one bypass Agent (accepted scope + locations
+ acceptance + gates), inspects diff, reruns decisive gates on real inputs, records `main=`/`impl=`,
sets unit DONE, commits `<scope> (M8.<u>): …`. Milestone → IMPLEMENTED when all units DONE → then
MILESTONE-REVIEW. Memory/roadmap/map updated at each implemented-unit close.
