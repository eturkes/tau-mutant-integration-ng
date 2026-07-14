# P6 - cross-cell-type response specificity

Status = ACTIVE. Next = S1 environment + compact all-cell pseudobulk substrate.

## Question + choice

Does amyloid remain microglia-led across the major brain cell classes, and is the
tau-dependent amyloid divergence (`interaction`) actually microglia-specific?

User delegated phase selection. Chosen default = re-baseline the six source-level
`broad_annotations` classes with replicate-correct pseudobulk, formal cross-cell-type
response contrasts, and one compact report plate. This directly stress-tests the
microglia-centred spine inside the same randomized 2x2 experiment.

Alternatives considered:

- Human validation: high translational value, but SEA-AD source files are absent;
  human Thal/Braak axes are observational + correlated; the live mouse mechanism
  targets that backed v1 signatures are retired. Defer until explicitly desired.
- Gene-level trajectory dynamics: potentially names interaction genes, but v1's
  110-gene result was ~1/4 ambient/ribosomal-confounded and depended on per-cell
  tradeSeq inference; current P2 supports composition, not further progression.
  Defer unless a trajectory-specific question justifies that risk.
- Broad-cell composition: useful validation, but nuclei capture/loading makes raw
  class proportions non-tissue abundance. Keep cell counts as precision/QC only.

## Research decision

Primary method = all-cell-type pseudobulk with `dreamlet` 1.10.0:

- biological unit = existing `genotype_batch` (16; 4/genotype);
- raw RNA counts aggregated by `broad_annotations x genotype_batch`;
- six fixed classes = Astrocyte, Microglia, Neuronal, Oligodendrocyte, OPC,
  Vascular; source labels stay fixed;
- per-class model = `~ 0 + genotype + batch`; five canonical contrasts;
- `dreamlet` precision weights retain all cells while modelling unequal observed
  cell counts/depth; source counts never become abundance evidence;
- primary calls = raw frequentist effects/CIs + study-wide BH per contrast;
  `getTreat(lfc = 0.5)` tests scientifically non-trivial effects;
- formal specificity = stacked pseudobulk model
  `~ 0 + stackedAssay:genotype + stackedAssay:batch + (1|genotype_batch)`;
  direct `(Microglia contrast) - (other-class contrast)` effects plus joint
  heterogeneity tests;
- `mashr`/`compositePosteriorTest()` = supportive shared-vs-specific prioritizer
  only. `dreamlet` documents cross-class overshrinkage, so posterior means never
  replace raw effects or generate the phase verdict;
- correlation-aware `zenith_gsa()` = contextual pathway readout; gene/pathway hit
  counts never define specificity by themselves.
- fixed project programmes = `canonical_microglia_markers` exactly as committed
  at plan open: Homeostatic, DAM, IFN, Proliferative, MHC_APC. Full mouse GO-BP
  stays exploratory; Figure 10 selects set names from, per story contrast, the two
  lowest-FDR positive + two lowest-FDR negative set/class rows, deduplicates exact
  memberships, then unions/caps at 12 and displays every selected set across all
  six classes with unsupported rows retaining their FDR state.

Why this supersedes v1 Arc N: v1 used native + one seed of K=289 cell-count
downsampling and NEBULA. That equalized neither transcript depth nor inferential
precision, discarded most neuronal data, and made significance tallies carry the
claim. P6 uses every biological unit, directly tests response differences, and
separates effect magnitude from discovery power. v1 outcomes are hypotheses only;
P6 thresholds + readouts are fixed before live P6 results.

Method/tooling basis:

- Squair et al. 2021: replicate-level pseudobulk controls false discoveries from
  cell-level pseudoreplication (Nature Communications 12:5692;
  https://doi.org/10.1038/s41467-021-25960-2).
- Hoffman et al. 2023: `dreamlet` precision-weighted pseudobulk supports complex
  designs, unequal cells, and cross-cell-type inference
  (https://doi.org/10.1093/bioinformatics/btad116).
- `dreamlet` manual: `stackAssays()` supports direct differential-response
  contrasts; `run_mash()` explicitly recommends MASH for specificity
  prioritization and warns that shared-effect shrinkage can over-shrink class
  effects (Bioconductor 3.23 package manual).
- Urbut et al. 2019: multivariate adaptive shrinkage learns shared/specific effect
  patterns across conditions (Nature Genetics 51:187-195;
  https://doi.org/10.1038/s41588-018-0268-8).
- Locked-repo probe 2026-07-14: the pinned 2026-06-22 repositories expose
  `dreamlet` 1.10.0, `variancePartition` 1.42.0, `zenith` 1.14.0,
  `mashr` 0.2.79, and `muscat` 1.26.0. `muscat` is unnecessary for the live path.

## Locked contracts

### Data

- `storage/data/snrnaseq.rds` = 8.95 GB compressed; read once inside an aggregation
  target. Never serialize a second full all-cell Seurat target.
- Extract RNA counts + required metadata, build compact pseudobulk, release parent.
- S1 validates exact six-class membership, all 16 units/class, canonical genotype
  crossing, one genotype/batch per unit, cell-count and UMI-depth distributions,
  positional feature/symbol mapping, and no missing design fields.
- Source broad labels enable a fair class comparison. Current pruned/reprocessed
  microglia remain the live P1 analysis; P6 carries an explicit bridge comparison.
- Sex remains omitted because it is aliased with batch in this experiment.

### Inference + claims

- Raw counts only; SCT/Harmony/UMAP never enter DE.
- Separate-class fit estimates all five canonical contrasts. Story focus =
  `nlgf_in_maptki`, `nlgf_in_p301s`, `interaction`; remaining contrasts stay in
  the compact audit payload.
- Specificity means a formally different response, not “significant here and not
  there.” Direct class-difference effects/CIs carry that claim.
- Meaningful gene-scale effect margin = 0.5 log2FC, fixed before results.
- Non-significance = unresolved. Similar/shared wording requires affirmative
  evidence: multiplicity-controlled two-one-sided equivalence at +/-0.5 log2FC
  where estimable and/or supportive MASH sharing probability, always beside raw
  90% equivalence CIs and 95% effect CIs.
- Study-wide BH is separate per canonical contrast across evaluated gene x class
  rows. Direct Microglia-vs-class difference tests use one declared family per
  contrast. Pathway correction spans set x class within contrast.
- Per-class discovery counts remain descriptive. Cell/UMI counts remain precision
  diagnostics. Neither becomes a biological ranking.
- Microglia bridge = compare P6 source-label effects with live pruned
  `pb_de_microglia`; report all five Spearman correlations + anchor directions.
  A synthetic contrast oracle is the fail-loud algebra gate; live divergence is
  interpreted as source-label/pruning sensitivity, not silently forced to agree.
- `mashr` is supportive because six-class shrinkage can pull effects toward a
  common value. Raw `dreamlet` coefficients, formal stacked contrasts, and CIs are
  load-bearing.
- Outcome-open wording: P6 may validate microglia enrichment, show broader/shared
  response, reveal a distinct non-microglial interaction, or remain underpowered.

### Report + scope

- New module = `R/celltype_specificity.R`.
- New targets, intended compact chain:
  `celltype_pseudobulk -> celltype_response -> celltype_specificity -> celltype_figures`.
- `_celltype.qmd` loads only `celltype_figures`; append after `_modality.qmd` so
  existing Figure 1-9 numbering stays stable and P6 becomes Figure 10.
- Figure 10 = one publication-style multi-panel plate: sampling/precision context,
  amyloid-response effects by class, formal interaction-difference effects, and
  shared/specific programme context. Exact grammar follows live result density.
- Preserve visible surface: numbered heading + figure + folded compact code; no
  title/TOC/caption/body prose/table/global code menu; nonblank alt text retained.
- Report render stays off the 8.95 GB source, pseudobulk matrices, and model fits.
  `celltype_figures` target size ceiling = 5 MB.
- Current lean posture holds: runtime guards + targeted smoke commands; add a
  committed test only if it directly protects Figure 10 or the inferential call.
- Every execution step runs its targeted fresh build/smoke and `scripts/check.sh`;
  before S5 the latter protects the existing nine-figure path while the new leaf
  targets remain outside `report`.

## Heavy/gated preconditions

- Observed host headroom at planning = 30 GiB RAM / 30 GiB swap, ~16 GiB RAM
  available. S1 loads/aggregates the 8.95 GB RDS alone, records elapsed time + peak
  RSS where feasible, and proves the full object is released before serialization.
- Add direct locked dependencies only after a clean resolver/load smoke. Update
  `rproject.toml` + `rv.lock` together; `rv sync` is mandatory for S1.
- `celltype_pseudobulk` fresh build must complete before downstream method work.
  Acceptance = compact serialized object <=250 MB, six assays/classes, 16 aligned
  columns each, finite positive libraries, and no dropped genotype/batch cell.
- `dreamlet` processing must retain all 16 units and the full model in every class;
  formula drops/errors are build-fatal.
- Planning-time fixed-design oracle for the locked stacked formula = 96 rows x 42
  columns, rank 42; S3 rechecks rank on live metadata before adding the unit random
  intercept.
- Stacked contrast algebra is proved on a synthetic fixture before live
  interpretation: each direct coefficient equals separately reconstructed
  `(class A contrast) - (class B contrast)` within tolerance.

## Steps

- [ ] **S1 - environment + compact substrate.** Add/sync `dreamlet` + direct
  namespaces actually called; warning-clean load/version smoke. Implement
  all-cell-to-pseudobulk adapter; add `celltype_pseudobulk`; force fresh build;
  validate data/design/memory/size contracts. Acceptance = all heavy preconditions
  above pass; no report change yet.
- [ ] **S2 - primary class-wise response.** Implement warning-fatal `dreamlet`
  processing + five-contrast fit + raw/topTreat extraction + study-wide families;
  add `celltype_response`; run synthetic coefficient checks and live microglia
  bridge. Acceptance = six classes x five contrasts, finite coefficient/SE/CI
  tables, all units/formulas retained, no target warnings, explicit bridge audit.
- [ ] **S3 - formal specificity + sharing.** Stack assays; fit the locked
  unit-correlated class-specific genotype/batch model; derive five
  Microglia-vs-other direct differences for each amyloid contrast + interaction
  (15 total) and one joint heterogeneity test per contrast.
  Run interaction MASH separately; emit raw, equivalence, and posterior classes
  without collapsing unresolved rows. Add `celltype_specificity`. Acceptance =
  synthetic oracle exact, live contrasts finite/identified, multiplicity families
  explicit, raw-vs-shrunk provenance complete.
- [ ] **S4 - pathway context + Figure 10 data.** Run correlation-aware fixed project
  programmes + mouse GO-BP context with deterministic display selection; build
  `celltype_figures` from compact tables only. Acceptance = predeclared sets always
  represented (including null/unmeasured states), GO selection reproducible,
  figure bundle <=5 MB, no heavy upstream object reachable from qmd.
- [ ] **S5 - report integration + QA.** Add `_celltype.qmd` after existing report
  includes; render one Figure 10 plate; update report source/argument wiring,
  map/memory as warranted. Run `scripts/check.sh`, rendered-DOM assertions, and
  Chromium PDF/PNG visual inspection. Acceptance = 10 numbered figures, 10
  nonblank alts, existing Figure 1-9 unchanged, no visible prose/captions/tables,
  no clipping/external refs/warning/error markers, gate green.

After S5 = CLOSE-OUT: adversarially review plan + code + figure claim boundaries;
fold durable decisions into history/memory/map; archive this plan; reset roadmap
Active plan; run final gate; commit scoped close-out.

## Definition of done

- The report answers both scope questions with replicate-unit inference and direct
  class-response comparisons, or explicitly reports unresolved evidence.
- Every specificity/shared statement is backed by the declared effect/CI family;
  no significance-only contrast, discovery-count ranking, or posterior-only call.
- v1 Arc N is independently re-baselined; P6 reports its own result without trying
  to reproduce the archived K=289/NEBULA verdict.
- Figure 10 is compact, conventional, legible, and additive: current nine figures
  remain byte/source-stable except unavoidable global numbering-independent render
  metadata.
- Full source/env/report documentation matches live wiring and one scoped commit
  closes each executed step.
