# GeoMx spatial deconvolution plan (arc L): tissue composition layer

Adds the project's first **cell-composition** readout. Every prior arc (D pathway,
E/F/G TF/kinase/CCC, H ledger, I NF-κB, human, J causal topology, K SCENIC) reads
**expression or activity**; none reads **how much** of each cell type / microglial
substate is present in tissue. The GeoMx is the only location-aware modality yet is
currently reduced to bulk-ish genotype DE (rmd/04). Here we deconvolve each of the
91 GeoMx ROIs into cell-type + microglial-substate **abundances** using the
snRNAseq as a reference (SpatialDecon — the GeoMx-WTA-native method), then test
whether amyloid reshapes microglial composition in tissue and whether **tau
background modulates that reshaping** (the interaction, as a composition effect),
with ROI XY coordinates enabling a spatial-gradient check.

Why this is an independent readout, not a restatement of the snRNAseq: the snRNAseq
substate proportions are sorting/loading-biased (DAM already shows a large
amyloid-driven jump in the reference counts — 859/847 → 3520/3548 under NLGF — but
that is a sorted-nuclei number, not tissue composition). Deconvolving the unsorted
GeoMx tissue gives a composition estimate that the snRNAseq cannot. Framed as
**orthogonal corroboration / honest reporting of what resolves**, NOT proof; it does
NOT touch the locked mouse 2×2 analysis except via an explicit, gated ledger-feed
decision (L5).

Readout = per-ROI cell-type/substate abundance + factorial contrasts (incl.
interaction) on microglial abundance + spatial gradient. A substate that does not
deconvolve cleanly is a finding (reference-collinearity caveat), not hidden.

## STATUS
- [x] L1 (GATE): DONE — granularity = **two-stage** (user choice); other 3 design choices locked to project defaults; SpatialDecon pipeline smoke-tested end-to-end (exit 0). See L1 COMPLETION NOTE.
- [x] L2: DONE — `R/spatial_decon.R` (6 pure fns) built + wired in helpers.R after scenic.R; smoke-tested end-to-end on live caches (exit OK, 113s). 3 correctness fixes vs L1 smoke (Q3-scaled background, dropped nuclei cell_counts, two-stage NA handling). See L2 COMPLETION NOTE.
- [x] L3: DONE — `scripts/build_spatial_deconvolution.R` → `spatial_decon.rds` (13 slots) + 3 results TSVs; built out-of-knit (exit 0, ~2.5 min, RSS ~20 G). Validation passed (18,512 common genes; stage1 6×91 + stage2 4×91 finite; consistency 1.0000). See L3 COMPLETION NOTE.
- [x] L4: DONE — `rmd/21_spatial_deconvolution.Rmd` (§22, 5 subsections, 7 in-knit ggplot figs + 4 tables, display-only) wired as `child-spatial-deconvolution` before `child-session`; writes `results/spatial_decon_verdict.tsv`. Clean knit (0 error / 0 warning); §22.1–22.5 + session→§23 verified. See L4 COMPLETION NOTE.
- [x] L5 (GATE): DONE — gate answer = **Feed 2 rows** (user). Added Phase-L block L-001 (amyloid DAM tissue-compartment expansion) + L-002 (null composition interaction = honest non-corroboration) to `build_biological_model_ledger.R` in a new `layer="composition"` evidence class; both margin-neutral (contests **verified 18/12/55**, ledger 85→87). Synced rmd/16 §17 prose to the rebuilt tallies; clean re-knit (0 error / 0 warning); arc closed (map + DIGEST + archive + commit). See L5 COMPLETION NOTE.

## Locked facts (verified at plan time 2026-06-05; do not re-derive)
- **Tool:** `SpatialDecon` **1.20.1** installed (Bioc 3.22). The NanoString GeoMx-WTA
  deconvolution method (Danaher 2022): log-normal regression of normalised
  expression on a cell-profile matrix + per-ROI negative-probe background.
  `spacexr`/`Giotto`/`MuSiC`/`DWLS`/`BisqueRNA`/`SCDC`/`granulator` are NOT installed
  (SpatialDecon is purpose-built for GeoMx; alternatives would need install + are
  Visium/bulk-oriented).
- **Reference object:** `storage/cache/seurat_full_processed.rds` — Seurat, **286,285
  cells × 33,683 genes**, assay RNA (layers counts + data; default RNA). `cell_type`
  is the **9-level** working identity: Astrocyte 24894, Neuronal 162954,
  Oligodendrocyte 47988, OPC 8153, Vascular 16192, **Microglia_DAM 8774,
  Microglia_homeostatic 13550, Microglia_IFN 1760, Microglia_proliferative 2020**
  (built by `scripts/build_seurat_full.R`: broad_annotations minus "Unknown", with
  Microglia replaced by the 4 §02a substate labels). `cell_type × genotype` shows
  amyloid-driven DAM expansion in the reference (sorting-biased — the motivation for
  tissue deconvolution, not a result to reproduce).
- **Reference gene IDs are MGI symbols** (`Xkr4`, `Gm1992`, …; 33,683, unique) — NOT
  Ensembl. (L1 smoke-test correction: differs from the Ensembl-keyed microglia object;
  all 33,683 are in `symbol_map$symbol`.) → **NO mapping needed**; intersect rownames
  directly with the GeoMx symbols = **18,512 common genes** (`create_profile_matrix`
  filtering keeps ~16,746 for the profile). `snrnaseq_symbol_map` is NOT used by this arc.
- **Target object:** `storage/data/geomx.rds` — Seurat, **19,963 genes × 91 ROIs**,
  assay RNA (counts + data). `data` = vendor **Q3-normalised** = `counts /
  q_norm_qFactors`. **Single AOI/segment per ROI** (`segment`/`aoi`/`tags` all
  1-valued) — geometric whole-tissue ROIs, NO morphological/plaque/Iba1 targeting.
  Per-ROI **negative-probe background** = `NegGeoMean_Mm_R_NGS_WTA_v1.0` (median 6.30,
  range 1.00–53.9) + `NegGeoSD_…` — exactly SpatialDecon's background input. XY =
  `ROI Coordinate X` / `ROI Coordinate Y` (+ per-slide `Scan Offset X/Y`); `nuclei`
  per ROI (median 38; **-1 = missing sentinel**, handle). 4 slides (`slide name`).
  genotype: MAPTKI 20, P301S 20, NLGF_MAPTKI 28, NLGF_P301S 23.
- **Existing GeoMx design to mirror:** rmd/04 fits `~ 0 + genotype + slide` (raw
  counts + voom + TMM). Reuse the SAME `genotype + slide` design for the abundance
  contrasts so the interaction stays cross-modality-comparable.
- **5 canonical contrasts everywhere** via `design.R::make_contrast_matrix`:
  `nlgf_in_maptki, nlgf_in_p301s, interaction, tau_alone, tau_in_nlgf`. Interaction =
  `(NLGF_P301S − P301S) − (NLGF_MAPTKI − MAPTKI)`, here on **log cell-type
  abundance**. 3 axes (no pre-privileged winner): amyloid_activation &
  synaptic_suppression use {nlgf_in_maptki, nlgf_in_p301s}; interaction_metabolic
  uses {interaction}.
- **Host:** 8 cores, 62 GB RAM, ~1 TB free. seurat_full load is ~1.9 G on disk /
  larger in RAM — L3 is the heavy step (load + profile build), run out-of-knit.

## Arc pattern (mirror J/K exactly)
Heavy compute OUT of the knit in a build script → light cache → display-only Rmd
that recomputes nothing. New `R/spatial_decon.R` sourced via `R/helpers.R` (after
`scenic.R`; reuses `design.R` factorial + `io.R` symbol mapping + `de_pb.R`
`fit_limma_log` for log-abundance limma-trend). Cache `storage/cache/spatial_decon.rds`.
Results `storage/results/spatial_decon_*.tsv`. Chapter `rmd/21_spatial_deconvolution.Rmd`
renders **§22** (file≠§ offset: rmd/20→§21, so rmd/21→§22; verify post-knit).

## Execution model (per step)
1. Check for a fitting Skill/subagent first (Explore for cross-file search;
   general-purpose for research; the `exploratory-data-analysis` /
   `scientific-visualization` skills may help L4 figures). Subagents use the largest
   model.
2. Smoke-test new helper / cache-reading code against live caches via `Rscript -e`
   BEFORE knitting (knits ~3 min; shape checks are cheap).
3. Mark the step DONE in STATUS + append a multi-paragraph completion note (what was
   built, file paths, biology, explicit deviations) under that step.
4. chown rstudio:rstudio all new files + root-owned knit outputs (analysis.html,
   storage/results/*.tsv, new storage/cache/*).
5. Re-knit (only from L4 on; L2–L3 build the cache outside the knit):
   `Rscript -e 'rmarkdown::render("analysis.Rmd", quiet=TRUE)'`; verify
   `grep -c 'class="error"' analysis.html` == 0 and same for `class="warning"`.
6. Commit locally: imperative subject <70 chars + HEREDOC body + Co-Authored-By
   trailer (this session's model).
Give L1, L3, L4, L5 fresh sessions (design-heavy or heavy compute); L2 may chain
onto L1 if context is low and the design is locked.

## Anti-anchoring guardrails (spatial-specific; enforce every step)
- **Lock the reference granularity + design in L1 BEFORE inspecting which abundances
  move.** Choosing the granularity after seeing whether DAM "comes out" is anchoring.
  State the policy, then honour it.
- **Report ALL 9 (or chosen-level) cell-type abundances even-handedly** + the
  deconvolution QC honestly (profile-matrix collinearity / condition number,
  per-ROI residuals, fit R²). A substate that does NOT resolve is reported as a
  reference-collinearity finding, never dropped to tidy the story.
- **No axis pre-privileged:** report amyloid + synaptic + interaction effects on
  abundance even-handedly; do not foreground the interaction.
- **State thresholds before applying:** FDR 0.10 (project standard), log-abundance
  offset, Moran's I permutation p, any deconvolution inclusion/QC cutoff,
  per-cell-type capping for the profile (if used). Record SpatialDecon defaults used.
- **Caveats stated at EVERY interpretation:** (a) the snRNAseq reference proportions
  are sorting/loading-biased and the deconvolution is NOT validated against
  ground-truth histology → abundances are model estimates; (b) the 4 microglial
  substates are transcriptionally similar (IFN/proliferative are small + close to
  DAM/homeostatic) → per-substate abundance may be collinear/unstable; (c) GeoMx
  ROIs are whole-tissue geometric samples with NO plaque/morphology targeting →
  shifts are regional-bulk composition, not plaque-niche; (d) single-nucleus
  reference vs GeoMx-WTA platform/segment-length differences.

---

## L1 (GATE): design lock + end-to-end smoke-test

**Goal.** Lock the four design choices below, then smoke-test that the SpatialDecon
pipeline runs end-to-end on a small subset (profile-matrix build → background →
`spatialdecon` → a `beta` abundance matrix), so L2/L3 build on a verified path.

**Decision gate at start — the ONE user-gated fork is reference granularity.** The
other three choices have a clear best option (locked as default below, alternative
noted but not gated, for cross-modality consistency).

1. **Reference granularity (USER GATE).**
   - **Default — 9-level direct:** build the profile at full granularity (5 broad +
     4 microglial substates); deconvolve all 9 at once. Pro: reads DAM/substate
     abundance per ROI directly (the headline). Con: the 4 substates are
     transcriptionally similar → collinear columns in the profile → per-substate
     abundance may be unstable. Mitigation: report profile collinearity + residuals;
     if substates do not resolve, the instability IS the finding (guardrail).
   - **Alt A — 6-level collapsed:** 5 broad + 1 pooled Microglia. Robust total-microglia
     abundance; forgoes per-substate (loses the DAM headline). Simplest.
   - **Alt B — two-stage:** stage-1 6-level (stable total microglia) → stage-2
     sub-deconvolve the microglial fraction into the 4 substates. Best-of-both;
     most complex; substate estimates conditional on stage-1.
   Present default + both alts, then `AskUserQuestion`.
2. **Abundance-testing model (LOCK default; alt noted).** Default: `log(beta + offset)`
   → `de_pb.R::fit_limma_log` (limma-trend) with `~ 0 + genotype + slide` →
   `make_contrast_matrix` → 5 contrasts. Mirrors every other modality's interaction
   contrast exactly (consistency > novelty). Alt (noted, not gated): CLR on
   `prop_of_all` + linear, or a compositional model (scCODA-style) — defer unless L4
   shows the limma reading is compositionally misleading.
3. **Spatial-autocorrelation method (LOCK default; alt noted).** Default: Moran's I
   per slide on per-ROI abundance using XY coords (k-nearest or inverse-distance
   weights), permutation p. 4 slides handled separately (offsets differ). Alt:
   variogram / Geary's C — noted, not gated.
4. **SpatialDecon I/O (LOCK; method-canonical).** Profile from snRNAseq **RNA counts**
   averaged per `cell_type` via `SpatialDecon::create_profile_matrix` (Ensembl→MGI
   mapped, intersected with GeoMx symbols). GeoMx input = **Q3 `data` layer** +
   background matrix from `NegGeoMean` (per-ROI geomean broadcast over genes, the
   `derive_GeoMx_background` pattern). Optional per-type cell cap for profile speed
   (state it if used). `nuclei` (drop -1 sentinels) feeds absolute-abundance scaling.

**Smoke-test (before declaring L1 done).** On a subsampled reference (e.g. ≤2k
cells/type) + the full GeoMx: build the profile, derive background, run
`spatialdecon`, confirm a `beta` matrix (9 or chosen-level × 91 ROIs) with finite
values + a sane fit. Record the locked decisions + smoke-test result in the L1 note.

### L1 COMPLETION NOTE (2026-06-05)
**Gate answer = two-stage** (user). The other three design choices are locked to the
project-consistent defaults above (limma-trend factorial on log-abundance; Moran's I
per slide; SpatialDecon Q3-data + neg-probe background). Validated end-to-end on a
1200-cell/type subsample of `seurat_full` + the full 91-ROI GeoMx (exit 0, ~2 min,
`/tmp/smoke_spatial_decon.R`, throwaway):

- **Gene handling (CORRECTED):** `seurat_full` rownames are MGI symbols, not Ensembl
  → drop the planned `symbol_map` mapping; `intersect(rownames(ref), rownames(geomx))`
  = 18,512; `create_profile_matrix(geneList=common, minCellNum=15, minGenes=100,
  normalize=TRUE, outDir=NULL)` keeps 16,746 profile genes. (My first smoke run
  wrongly mapped symbol→Ensembl, gutting the matrix to 248 genes / ~1 nonzero per
  cell — the bug that surfaced the keying fact.)
- **SpatialDecon I/O (validated):** GeoMx input = Q3 `data` layer (linear);
  background = `NegGeoMean_Mm_R_NGS_WTA_v1.0` per ROI broadcast across genes
  (`sweep(matrix(1,G,N), 2, negmean, "*")`); `spatialdecon(norm, bg, X, align_genes=
  TRUE)` → `$beta` (celltypes × 91 ROI), all finite. (`derive_GeoMx_background` exists
  but expects neg-probe ROWS in `norm`; the WTA object stores the geomean in metadata,
  so the manual broadcast is the right path.)
- **Two-stage mechanic (validated, locked):** Stage-1 = 6-level profile (5 broad +
  pooled Microglia) → `beta6`; total-microglia = `beta6["Microglia",]` (robust). Stage-2
  = 9-level profile → `beta9`; within-microglia substate fractions `f_s =
  beta9[substate,]/Σ_substates`, anchored to absolute abundance `f_s ×
  beta6["Microglia",]`. Stage-1 total ≡ Σ stage-2 substates (**Spearman 1.000**), so
  the anchoring is near-exact. Both layers reported; stage-2 caveated (conditional +
  substate collinearity).
- **Biology preview (capped, UNTESTED — do not over-read; L3 does the real stats):**
  stage-1 total-microglia mean abundance MAPTKI 0.5 / P301S 1.7 / NLGF_MAPTKI 6.9 /
  NLGF_P301S 6.0 (sharp amyloid-driven tissue expansion). Stage-2 within-microglia DAM
  fraction 0.050 / 0.015 / 0.214 / **0.324** — DAM rises under amyloid AND is higher on
  the P301S tau background (NLGF_P301S > NLGF_MAPTKI): a candidate tau×amyloid
  interaction the two-stage split isolates from the total-microglia signal (which leans
  the other way). Confirms the pipeline yields sensible, testable output.
- **L3 to-do flagged by smoke:** (a) `create_profile_matrix` densifies internally
  (2.7 GiB at cap 1200) → for the full build CAP cells per type (~3000/type; the
  profile is a per-type mean, so capping is lossless for the mean and avoids a ~38 GB
  densification at 286k cells). (b) `res$prop_of_all` returned NaN in my quick
  per-genotype summary (orientation/shape differs from assumed) → in L2 derive
  proportions from `beta` directly, or inspect `prop_of_all` structure before using it.

## L2: `R/spatial_decon.R` helper

**Goal.** A small set of pure functions, sourced in `R/helpers.R` after `scenic.R`.
Smoke-test each against live caches via `Rscript -e` before any knit.
- `build_reference_profile(sc, cell_type_col, geomx_symbols, cap=3000L)` →
  gene × cell-type profile (wraps `create_profile_matrix`, `geneList=intersect(rownames,
  geomx_symbols)`, `minCellNum=15, minGenes=100, normalize=TRUE, outDir=NULL`; caps
  cells/type for the densification; reports coverage + collinearity/condition number).
  NO symbol mapping — `seurat_full` rownames are already MGI symbols (L1 fact).
- `derive_geomx_background(geomx)` → gene × ROI matrix = per-ROI `NegGeoMean` broadcast
  across genes (manual sweep — `SpatialDecon::derive_GeoMx_background` wants neg-probe
  rows the WTA object lacks).
- `run_spatialdecon(geomx, profile, background)` → list(beta, prop_of_all, yhat,
  resid, sigmas, fit_qc). Wraps `SpatialDecon::spatialdecon(align_genes=TRUE)` on the
  Q3 `data` layer. (Inspect `prop_of_all` shape; derive proportions from `beta` if it
  is awkward — L1 to-do.)
- `combine_two_stage(beta6, beta9)` → stage-1 total-microglia + 5 broad (from beta6);
  stage-2 within-microglia substate fractions anchored to the stage-1 total (from
  beta9). Returns both layers (validated: Σ-9L ≡ 6L microglia, Spearman 1.000).
- `fit_abundance_contrasts(abund, meta, offset)` → reuse `design.R` factorial +
  `de_pb.R::fit_limma_log`; per-contrast top tables (cell_type × stat) for the 5
  contrasts. Run on BOTH the stage-1 6-level abundances and the stage-2 substate layer.
- `spatial_autocorrelation(abund, coords, slide)` → per-slide Moran's I + perm p per
  cell type.
Keep functions side-effect-free (no file writes); L3 orchestrates + writes.

### L2 COMPLETION NOTE (2026-06-05)
Built `R/spatial_decon.R` (6 pure, side-effect-free functions + private
`.morans_I`) and wired `source("R/spatial_decon.R")` into `R/helpers.R` after
`scenic.R`, before `microglia.R`, with a matching dependency-order doc comment.
Functions: `build_reference_profile`, `derive_geomx_background`,
`run_spatialdecon`, `combine_two_stage`, `fit_abundance_contrasts`,
`spatial_autocorrelation`. Each reuses the locked machinery (`factorial_design`
from design.R + `fit_limma_log` from de_pb.R for the 5-contrast log-abundance
factorial; `%||%` from utils.R). Smoke-tested the WHOLE chain end-to-end on a
cap=600/type `seurat_full` subsample + the full 91-ROI GeoMx (profile6+profile9
→ background → spatialdecon×2 → combine_two_stage → fit_abundance_contrasts on
both layers → spatial_autocorrelation raw+residualised): exit "ALL FUNCTIONS OK",
113s, every shape/finiteness assertion passed. The two-stage NA fix was
additionally verified with a synthetic unit check.

**Three correctness improvements over the L1 smoke (which was a runs-check, not a
correctness lock).** (1) **Q3-scaled background.** Verified on the live object that
`data == counts / q_norm_qFactors` exactly (data/counts ratio constant across all
19,963 genes per ROI, sd ~6e-17; qFactor 0.10–7.49, median 1.04). SpatialDecon's
`norm` and `bg` must share a scale, so `derive_geomx_background` divides
`NegGeoMean` by `q_norm_qFactors` per ROI to land in the same Q3 space as the
`data` layer — the L1 smoke broadcast the RAW NegGeoMean (a scale mismatch up to
~7×). (2) **`nuclei`/`cell_counts` dropped.** 42/91 ROIs carry the `-1` missing
sentinel, so passing `cell_counts` to `spatialdecon` would NA half the ROIs; beta
stays on its native `normalize=TRUE` abundance scale (absolute-count rescaling is
abandoned, not just deferred). (3) **`combine_two_stage` unresolved-ROI handling.**
An ROI whose 9-level substate betas are all zero gives a 0/0 within-microglia
fraction; those columns are set NA (a zero cannot be split — honest) and counted in
`n_unresolved`, and `consistency_spearman` is computed over complete columns
(returns exactly 1.000, confirming Σ-stage2 ≡ stage-1 microglia total on resolved
ROIs).

**Two design choices made at L2.** (a) **Moran's I is hand-rolled** (no new dep):
`spdep` is absent and `ape::Moran.I` returns a normal-approximation p, not the
permutation p the plan locked, so `.morans_I` + a 999-permutation one-sided
(positive-autocorrelation) test is used. (b) **Optional `covar` residualisation in
`spatial_autocorrelation`.** The slide×genotype crosstab shows slides MIX genotypes
(e.g. "290724 1st slide" carries all four), so a raw per-slide Moran's I is
confounded by within-slide genotype composition. The smoke made this concrete: raw
Moran's I is strongly positive and significant (Neuronal/Astrocyte/Vascular
I≈0.3–0.46, p=0.005) but collapses to near-null once genotype is residualised out
(Astrocyte I=0.12 p=0.045, the rest n.s.) — most of the apparent spatial structure
IS genotype composition. L4 must report the residualised result as the honest
spatial-gradient read and caveat the raw. `offset` for `log(beta+offset)` defaults
to 1 (a one-cell pseudocount; `normalize=TRUE` betas are on a cell-abundance scale);
L3 records the realised value.

**Biology/QC surfaced (preview only — L3 runs the real stats).** The profile QC
makes the guardrail's substate-collinearity concern visible and quantified, not
hidden: profile6 is well-conditioned (κ=4.5, max off-diagonal column cor 0.84) but
profile9 is collinear (κ=22.9, **max column cor 0.978**), and SpatialDecon pushes
Microglia_IFN to ~0 across ROIs (constant → logFC=0,t=0 in the interaction fit) —
the expected "small, transcriptionally-close substate does not resolve cleanly"
finding, reported not dropped. Direction sanity: the stage-2 DAM interaction logFC
is positive (+0.13 at cap=600), consistent with the L1 preview (DAM higher on the
P301S tau background under amyloid). **Memory:** `build_reference_profile`
pre-subsets the count matrix to the 18,512 common genes and caps cells/type BEFORE
`create_profile_matrix` (which densifies internally), halving the densification vs
the L1 smoke that passed all 33,683 genes (2.7 GiB at cap 1200). One benign stderr
message ("more than one class 'dist' … spam/BiocGenerics") appears when
`stats::dist` runs under loaded `spam`; it is cosmetic and out-of-knit (L3 only), so
it cannot affect analysis.html. No knit at L2 (cache is built out-of-knit at L3).

## L3: `scripts/build_spatial_deconvolution.R` → cache + TSVs

**Goal.** Idempotent build script (params embedded), run out-of-knit. Read
`seurat_full_processed.rds` + `geomx.rds`; build BOTH profiles (6-level pooled +
9-level substates, cap ~3000 cells/type); derive background; run `spatialdecon`
twice → `beta6`, `beta9`; `combine_two_stage` → stage-1 (total microglia + 5 broad)
+ stage-2 (anchored substate) abundance layers; fit the 5 factorial contrasts on
each layer; spatial autocorrelation → write `storage/cache/spatial_decon.rds` =
list(profile6, profile9, beta6, beta9, stage1, stage2, fit_qc, abundance_contrasts,
spatial_autocorr, params, meta) + results TSVs: `spatial_decon_abundance_by_genotype.tsv`,
`spatial_decon_contrasts.tsv`, `spatial_decon_spatial_autocorr.tsv`. Validate object
shape + log coverage/QC (profile collinearity, common-gene count 18,512, per-ROI fit).
chown outputs. (Heavy: 1.9 G reference load; cap keeps `create_profile_matrix`
densification tractable — smoke hit 2.7 GiB at cap 1200.)

### L3 COMPLETION NOTE (2026-06-05)
Built `scripts/build_spatial_deconvolution.R` (idempotent, `--overwrite` guard;
mirrors the `build_scenic_contrasts.R` I/O-wrapper pattern — all compute
out-of-knit, rmd/21 will be display-only). It sources `R/helpers.R`, loads
`seurat_full_processed.rds` + `geomx.rds`, derives a 6-level `cell_type6` column
(collapse `Microglia_*`→`Microglia`), builds **both** reference profiles
(`build_reference_profile` cap=3000/type, seed=1), derives the Q3-scaled
background, runs `spatialdecon` ×2, two-stage-assembles, fits the locked
5-contrast 2×2 factorial on each layer, runs per-slide Moran's I (raw +
genotype-residualised), and writes `storage/cache/spatial_decon.rds` (13 slots:
profile6/9, beta6/9, stage1/2, two_stage, fit_qc, abundance_contrasts,
abundance_by_genotype, spatial_autocorr, meta, params) + 3 TSVs
(`spatial_decon_{abundance_by_genotype,contrasts,spatial_autocorr}.tsv`). Ran
clean (exit 0, ~2.5 min, peak RSS ~20 G on the 62 G host); all outputs
chowned rstudio:rstudio. The only stderr noise is the expected out-of-knit trio
(nichenetr e1071/ggplot2 import, create_profile_matrix dense-coercion warnings,
the cosmetic spam/BiocGenerics `dist` message) — none touch analysis.html.

**Deconvolution QC (guardrail: report collinearity even-handedly).** Common
genes ref∩GeoMx = **18,512** (matches the L1 lock exactly). After
`create_profile_matrix` filtering: profile6 = 16,990 genes × 6 types, **κ=3.9,
max|cor|=0.853** (well-conditioned); profile9 = 17,097 genes × 9 types, **κ=35.1,
max|cor|=0.989** — the substate columns are collinear, the predicted
"transcriptionally-close substates → unstable betas" finding made quantitative,
reported not hidden. Per-ROI fit `resid_rmse=0.835` for both deconvolutions
(identical to 3 dp because the broad-type fit dominates the global residual; the
microglia 4-way split barely moves it). The two-stage consistency check is exact
on resolved ROIs (Σ-stage2 ≡ stage1 microglia total, **Spearman 1.000**).

**New fact surfaced at full scale: `n_unresolved = 5`.** 5/91 ROIs have all-zero
microglial substate betas in the 9-level fit (no resolvable microglia to split),
so `combine_two_stage` sets their stage-2 substate columns NA (a zero cannot be
split — honest). Two downstream consequences, both handled in-script: (a) the
stage-2 contrast fit drops those ROIs per limma's row-wise NA handling (it still
returned a clean 4-substate × 5-contrast table); (b) the spatial-autocorrelation
input is `rbind(stage1, stage2)` **only when n_unresolved==0**, else stage1 only
— with 5 unresolved here, the per-slide Moran's I is reported for the **6 stage-1
broad types only** (an NA substate column would NaN the statistic). L4 must state
that the spatial gradient is a broad-type read; substate spatial structure is not
estimable on this ROI set.

**Biology (honest, sign-aware; the real stats, not the smoke preview).** Amyloid
(NLGF) MAIN effects dominate and are strong; the interaction is NOT significant.
Within-microglia the headline is a **DAM↑ / proliferative↓ substate swap under
amyloid**: DAM fraction 0.18 (MAPTKI) / 0.18 (P301S) → **0.79 (NLGF_MAPTKI) /
0.93 (NLGF_P301S)**, DAM abundance logFC +0.40 (`nlgf_in_maptki`) / +0.53
(`nlgf_in_p301s`), both **FDR<0.001**; proliferative is the mirror image
(0.82/0.82→0.21/0.07; logFC −0.45/−0.59, FDR<0.001). Pooled-microglia tissue
**fraction** rises (0.37/0.45→0.72/0.75) while its absolute beta stays ~flat
(broad types dilute), so the composition shift is relative — the two-stage split
is what isolates the DAM redistribution from the flat total. The broad types
(Astrocyte/Vascular/Neuronal/Oligodendrocyte) fall under amyloid in both
backgrounds (relative dilution). **Interaction: 0/10 cell types significant**
(best adj.P=0.25, Astrocyte); the DAM interaction point estimate is positive
(+0.13, i.e. DAM higher on the P301S tau background under amyloid — same
direction as the L1 preview and the raw means 0.93>0.79) but t=0.98, **p=0.33** —
the composition layer does **not** corroborate the tau×amyloid interaction at
FDR 0.10. IFN and homeostatic flatline at interaction (logFC=t=0; SpatialDecon
pushes the small, collinear IFN substate to ~0) — the κ=35 collinearity made
concrete in the contrast table. **Spatial:** raw per-slide Moran's I is broadly
positive (14/24 sig) but collapses to near-null once genotype is residualised out
(**1/24 sig**: Astrocyte on "290724 1st slide", I=0.10 p=0.043) — the apparent
spatial structure IS within-slide genotype composition (L2 finding reproduced at
full scale). The honest spatial read is therefore **no robust within-genotype
gradient**.

**Caveat to carry (L4/L5):** "proliferative" being the dominant *baseline*
microglial substate in the deconvolution is plausibly a reference-profile
artifact (its profile may absorb generic/cycling-agnostic microglia), not a
literal claim that most baseline microglia proliferate — flag in §22.1/§22.3
alongside the standard sorting-bias / no-histology-ground-truth / whole-tissue-ROI
caveats. **This pre-stages the L5 gate:** the clean, significant, new-evidence-class
result is **amyloid-driven DAM tissue expansion** (candidate Suggestive
corroboration of the amyloid-activation axis / Hyp-1), explicitly **NOT** the
interaction (composition interaction is null) — any ledger feed must be scoped to
the amyloid axis and to *tissue composition* to avoid double-counting the snRNAseq
DE already in the ledger.

**Cache-shape deviation from the plan spec (additive, benign):** the cache carries
two slots beyond the listed set — `two_stage` (full `combine_two_stage` output:
microglia_total, substate_fractions, n_unresolved, consistency) and
`abundance_by_genotype` (the per-genotype mean abundance+prop summary, 40 rows) —
both convenience inputs for rmd/21; nothing is removed. `meta` is a clean 91×6
ROI frame (roi, genotype, slide, x, y, nuclei) so rmd/21 + the design avoid the
`"slide name"` space. `params` records every locked threshold (cap 3000,
minCellNum 15, minGenes 100, offset 1, padj 0.10, n_perm 999, seed 1, design
"genotype + slide", Q3-scaled background, cell_counts_used=FALSE).

## L4: `rmd/21_spatial_deconvolution.Rmd` (→ §22), display-only

**Goal.** New chapter reading ONLY `spatial_decon.rds` (recompute nothing). Wire as
`child-spatial-deconvolution` before `child-session` in analysis.Rmd. Subsections:
- 22.1 deconvolution QC (profile collinearity, per-ROI fit/residuals, coverage).
- 22.2 composition by genotype (stacked/abundance plots, all cell types even-handed).
- 22.3 factorial contrasts on microglial abundance incl. the interaction (DAM the
  focus, but all substates reported); sign-aware.
- 22.4 spatial gradient (Moran's I per slide; XY abundance maps).
- 22.5 verdict (multi-axis, no pre-privileged winner; writes
  `results/spatial_decon_verdict.tsv`).
Plot helpers live in the chapter's own R module or `plot.R` so they render inside
the shared knit session (mirror §20's in-module plotting; do NOT use the Python
visualization skill for in-knit figures). Knit-verify 0 error / 0 warning.

### L4 COMPLETION NOTE (2026-06-05)
Built `rmd/21_spatial_deconvolution.Rmd` as a display-only chapter that reads
ONLY `storage/cache/spatial_decon.rds` (recomputes nothing — mirrors the §20
arc-K pattern) and wired it as `child-spatial-deconvolution` between
`child-scenic-regulons` and `child-session` in `analysis.Rmd`. It renders **§22**
(confirmed: rmd/20→§21, rmd/21→§22; session info shifted to §23) with five
subsections matching the plan spec exactly: 22.1 deconvolution QC, 22.2 tissue
composition by genotype, 22.3 factorial contrasts on abundance (incl. the
interaction), 22.4 spatial gradient, 22.5 verdict. All figures are in-knit
ggplot2 (7 figs) built inside the shared knit session — no Python visualization
skill, no new `plot.R` helpers needed (the chapter's plots are one-off, so inline
chunks match §20). The verdict chunk writes `storage/results/spatial_decon_verdict.tsv`
(3 axis rows). Full knit clean: `grep -c 'class="error"'` = 0 and
`class="warning"` = 0; 4m53s.

**Smoke-tested before the full knit** (execution-model step 2): rendered rmd/21
alone through a throwaway temp parent (`/tmp`, since cleaned) with
`knit_root_dir` = project root and `results_dir` pointed at `/tmp` to keep
`storage/results` pristine — exit 0, 0 error / 0 warning blocks, 7 `<img>`,
verdict TSV numbers byte-correct. This caught the inline-R + figure-build paths
exactly as the real knit would, at ~10s vs the 5-min knit.

**Chapter design choices.** (a) **Prose syncs to the cache via inline R**: the
load chunk computes every headline scalar (κ, DAM/proliferative logFC+FDR, pooled-
microglia fractions, interaction nullity, Moran's raw-vs-residualised counts) and
the narrative references them with `r sprintf(...)`, so the text cannot drift from
the cache (the §20 discipline). (b) **Even-handed ordering** per the guardrail:
broad classes and substates are `sort()`ed (alphabetical, neutral — DAM is not
foregrounded by position), all 6 broad + 4 substates reported in every table/fig.
(c) **Figures**: profile9 between-type collinearity heatmap (visualises the 0.97–
0.99 substate block / κ=35), stacked broad + stacked substate composition bars
(each sums to 1), a stage-faceted contrast heatmap (logFC fill, * = FDR<0.10), an
interaction forest with approximate Wald 95% intervals (all cross zero — the null
made visual), a per-slide XY abundance map (viridis), and a raw→residualised
Moran's I paired-collapse plot. (d) **Verdict TSV** keeps the project's 3-axis
frame for cross-modality comparability: amyloid_activation = corroborated (new
evidence class, tissue composition); interaction_metabolic = direction-concordant
but null; synaptic_suppression = no abundance endpoint (honest non-finding,
mirroring §14/§20).

**Biology rendered (the L3 stats, now in-document, honest + sign-aware).**
Headline = amyloid-driven **DAM↑ / proliferative↓ substate swap** in tissue: DAM
within-microglia fraction 18%/18% → 79%/93% (NLGF_MAPTKI/NLGF_P301S), DAM
abundance logFC +0.40 (FDR 3.5e-5) / +0.53 (FDR 2.4e-7); proliferative the mirror
(−0.45/−0.59). Pooled-microglia tissue fraction rises 37%/45% → 72%/75% while
absolute β is ~flat (broad-type dilution). **Interaction null** (0/10 sig; DAM
+0.13, t=0.98, p=0.33 — same sign as the snRNAseq but not significant; best
Astrocyte FDR 0.25). **Spatial**: raw Moran's I 14/21 sig collapses to 1/21 after
genotype-residualisation (only Astrocyte on "290724 1st slide", I=0.10, p=0.043) —
no robust within-genotype gradient. Caveats carried at every interpretation:
sorting-biased reference / no histology ground truth; substate collinearity
(κ=35) → IFN+homeostatic collinearity-suppressed to 0, only DAM+proliferative
trustworthy; "proliferative" baseline dominance plausibly a profile artefact;
whole-tissue geometric ROIs (no plaque targeting) → regional-bulk not plaque-niche;
sn-vs-WTA platform difference.

**No deviations from the plan spec.** chowned rmd/21, analysis.html,
spatial_decon_verdict.tsv, analysis_files/ to rstudio:rstudio. **This pre-stages
the L5 gate (next, GATE):** the clean new-evidence-class result for a ledger feed
is the *amyloid* DAM tissue expansion (Suggestive, scoped to tissue composition to
avoid double-counting the snRNAseq DE), explicitly NOT the interaction (composition
interaction is null).

## L5 (GATE): ledger feed + arc close

### L5 GATE ANSWER (2026-06-05) = Feed 2 rows
User chose **Feed 2 rows** over the 1-row default and the no-feed alternative. The
fresh L5 session executes this exact path:

1. **Read first** (do not guess the schema): `scripts/build_biological_model_ledger.R`
   — the hardcoded `row()` calls grouped by phase (H/I/J/K). Use the **Phase K rows
   K-001 (Spi1) / K-002 (Rel)** as the template (the SCENIC §21→§17 corroboration
   feed): Suggestive, existence-not-activity, margin-neutral. Mirror their `row()`
   argument shape exactly.
2. **Add a Phase L block (2 rows):**
   - **L-001 — amyloid DAM tissue-compartment expansion.** Evidence class = *tissue
     composition* (GeoMx SpatialDecon two-stage), axis = **amyloid_activation /
     Hyp-1**, scope = existence of the composition shift (DAM↑/proliferative↓ in
     both backgrounds; logFC +0.40/+0.53, FDR <1e-4), **Suggestive**. Explicitly
     distinct from the snRNAseq expression DE already counted — it corroborates the
     amyloid axis from the *unsorted tissue*, a new evidence class. Cite §22.
   - **L-002 — null composition interaction (honest non-corroboration).** The
     tau×amyloid interaction does NOT surface as a composition effect (0/10 cell
     types sig; DAM interaction +0.13, direction-concordant with the snRNAseq but
     p=0.33). Record as a new-evidence-class **negative** result on the
     interaction_metabolic / divergence axis (existence-of-non-finding), not a
     corroboration. This is the honest-discordance counterpart to L-001.
3. **Margin-neutrality is REQUIRED for both** (mirror K-001/K-002): after
   `build_biological_model_ledger.R` + `build_biological_model_adjudication.R`
   rebuild, verify the locked contests stay **18/12/55** (read
   `results/biological_model_contest_verdicts.tsv` / the adjudication TSV). If a row
   would perturb a contest, engineer it margin-neutral (scope/strength wording) as
   K-001/K-002 did. Do NOT relitigate the locked margins.
4. **Re-knit** `analysis.Rmd` (the §17 ledger + §16 model chapters read the rebuilt
   TSVs) → verify 0 error / 0 warning. chown outputs.
5. **Arc close:** refresh `map.md` (add the rmd/21→§22 pipeline row + the
   `spatial_decon.rds`→21 / `spatial_decon_verdict.tsv` wiring + the new Phase-L
   ledger rows in the biological-model note); add a DIGEST.md entry for arc L;
   archive this plan to `completed/spatial_deconvolution_plan_2026-06-05.md`. Commit.

(Caveat to encode in the row prose, per the L3/L4 guardrails: composition is a
model estimate — sorting-biased reference, no histology ground truth, κ=35 substate
collinearity, whole-tissue non-plaque ROIs. That model-dependence is WHY L-001 is
Suggestive, not Strong.)

**Decision gate at end (K6-style).** Composition/abundance is a NEW evidence class
(all 85 ledger rows are expression/activity). If the deconvolution shows amyloid-driven
microglial-abundance expansion with tau-background modulation in tissue, that is a
candidate Suggestive corroboration of the amyloid-activation axis (Hyp-1) and possibly
the interaction. Options to present (default + alts), scoped to avoid double-counting
the snRNAseq DE already in the ledger:
- **Default — add margin-aware Suggestive row(s)** scoped to *tissue composition*
  (abundance), explicitly distinct from the expression footprints already counted;
  engineer margin-neutrality if the row would otherwise perturb the locked contests
  (mirror K-001/K-002).
- **Alt — standalone chapter, no ledger feed** (if the composition signal is weak,
  unstable, or compositionally ambiguous).
Then close the arc: refresh `map.md` + `completed/DIGEST.md`, archive this plan to
`completed/spatial_deconvolution_plan_<YYYY-MM-DD>.md`.

### L5 COMPLETION NOTE (2026-06-05)
Executed the gate answer (**Feed 2 rows**) exactly. Read
`scripts/build_biological_model_ledger.R` and mirrored the Phase-K K-001/K-002
`row()` shape (14 args) for a new `phase_l <- rbind(row("L-001", ...),
row("L-002", ...))` block inserted before `# ---- ASSEMBLE + WRITE`, appended
`phase_l` to the `ledger <- rbind(...)` call, and extended the `layer` validation
vector with **`"composition"`** (a genuinely new evidence class — all 85 prior
rows are expression/activity). Verified the adjudication script consumes only
`claim_id`/`axis`/`confidence_grade`/`supports_models`/`contradicts_models` and
**ignores `layer`**, so only the builder's validation needed extending.

- **L-001** (axis `amyloid_activation`, layer `composition`, Suggestive): amyloid
  DAM tissue-compartment expansion (within-microglia DAM fraction 0.18/0.18 →
  0.79/0.93; DAM abundance logFC +0.40 `nlgf_in_maptki` FDR 3.5e-5 / +0.53
  `nlgf_in_p301s` FDR 2.4e-7; proliferative the mirror; pooled-microglia tissue
  fraction 0.37/0.45 → 0.72/0.75). `supports="Hyp-1A;Hyp-1B;T-Inflammation"`,
  `contradicts=""`. Margin-neutral by the K-001 mechanic: supporting BOTH amyloid
  models lifts each net by one without moving the amyloid margin.
- **L-002** (axis `interaction_metabolic`, layer `composition`, Suggestive): the
  null composition interaction (0/10 cell types FDR<0.10, best Astrocyte adj.P
  0.25; Microglia_DAM interaction logFC +0.13, t=0.98, p=0.33; IFN + homeostatic
  collinearity-suppressed to 0). `supports=""`, `contradicts=""` — an honest
  non-corroboration / existence-of-non-finding that feeds NO model (precedent: the
  empty/empty context rows). Margin-neutral by touching no contest.

Rebuilt both TSVs and **verified the locked contests EXACTLY** —
`amyloid_activation: Hyp-1B by margin 18 (14 vs 32); synaptic_suppression: Hyp-2B
by margin 12 (-4 vs 8); interaction_metabolic: Hyp-3B by margin 55 (-27 vs 28)` —
row count 87, `composition` layer = 2, L-001/L-002 supports/contradicts correct.
Then synced the `rmd/16_biological_model.Rmd` §17 hardcoded prose to the rebuilt
tallies (the §17 dynamic tables read the TSVs and would otherwise contradict the
prose): ledger 85→87 + Phase-L clauses; lengthMenu 85→87; amyloid margin
`(31-13)`→`(32-14)` (still 18); Hyp-1A 22→23 supporting rows (K-001 + L-001 each
support both amyloid models); T-Inflammation 25/28→26/29 (axis-1 share 89%→90%);
cross-entity ranking Hyp-1B 31→32, T-Inflammation 28→29 reordered ABOVE Hyp-3B
(net-29 tie broken by Strong: T-Synergy 9 > T-Inflammation 6, both above Hyp-3B
28); per-axis-anchoring 89-94%→90-94%; added a Phase-L paragraph to 17.5 mirroring
the Phase-K one. Verified every hardcoded number against a direct ledger tally
(Hyp-1B 32/11-Strong, T-Inflammation 29/6, T-Synergy 29/9, Hyp-3B 28/8, Hyp-1A 23
sup/9 con/14 net/6 Strong, Hyp-2B 8, Hyp-0 7, Hyp-2A −4, Hyp-3A −27) and confirmed
the §17 table sort (`desc(net_support), desc(n_strong_supports)`) matches the prose
ordering. Re-knit `analysis.Rmd` clean (0 `class="error"` / 0 `class="warning"`,
4m53s); `analysis.html` self-contained 35 MB, chowned rstudio:rstudio.

Arc closed: refreshed `map.md` (pipeline row rmd/21→§22; verified-anchor `21→§22`;
`spatial_decon.rds` cache row + `seurat_full_processed` consumer; spatial-lane
scripts note; `R/spatial_decon.R` in the helper source order; Phase-L extension of
the biological-model TSV note); added the `spatial_deconvolution_plan_2026-06-05`
DIGEST entry + arc-L clause in the narrative-arc intro; archived this plan.
**No deviations from the L5 gate answer.** summary.Rmd is OUT of scope (a separate
standalone, not rebuilt by the analysis.Rmd knit; its pre-existing "margin 53"
interaction citation is stale vs the live 55 but was not introduced by arc L).
