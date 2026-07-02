# P4 Cross-modality - plan

## Scope
Build the cross-modality evidence layer around the rebuilt P1-P3 story:
amyloid drives microglial DAM activation, mutant tau modulates the amyloid
response mainly by adding DAM cells, and the rebuilt mechanism layer supports a
Myc-linked RNA signal while downgrading NF-kB attenuation and Gsk3b.

P4 asks what the non-snRNAseq modalities add:
- GeoMx spatial WTA: spatial/tissue transcriptomic DE and, only if it earns it,
  light deconvolution for tissue cell-state abundance.
- Bulk proteome + phosphoproteome: total protein, raw phosphosite, and
  protein-corrected phosphosite effects over the same 24M 16-run sample set.
- Synaptic / clearance axis: focused CCC-lite around APP/TREM2/APOE/CD74/MERTK
  biology if measured support exists; no full off-lock CCC stack by default.
- Integration: one compact gene/pathway/axis view showing concordance,
  divergence, nulls, and caveats across modalities, without reviving the v1
  claims ledger / contest machinery.

Inputs already built: `geomx`, `proteomics`, `phospho`, `sample_key`,
`pb_de_microglia`, `pb_de_substate`, `composition_results`,
`trajectory_report`, `mechanism_report`, `mechanism_gene_sets`,
`phospho_de_24m`, `kinase_activity`, `kinase_mechanism_summary`, `symbol_map`.

Outputs: P4 modality targets, one compact `crossmodality_report`, and
`_crossmodality.qmd`. P5 remains the final lean synthesis chapter.

## Research Digest
Local / v1:
- Current GeoMx target = 19,963 genes x 91 AOIs; genotype counts 20/20/28/23,
  4 slides, one segment, 15 genotype:bio_rep biological units, and repeated
  AOIs per unit. Negative-probe and Q3 metadata are present:
  `NegGeoMean_Mm_R_NGS_WTA_v1.0` and `q_norm_qFactors`. `nuclei` has many `-1`
  sentinels, so any deconvolution must stay on beta/log-abundance scale rather
  than nuclei-rescaled absolute counts.
- v1 GeoMx DE used raw counts + limma-voom + slide fixed effect. v1 spatial
  deconvolution used SpatialDecon and found useful composition evidence but
  carried load-bearing gotchas: Q3-scale the negative-probe background, drop
  nuclei rescaling, expect collinear microglia substate profiles, and report
  residualised spatial autocorrelation rather than raw slide-confounded Moran's I.
- v1 bulk proteomics summed peptide/PTM rows to parent protein groups, median
  normalised log2 intensities, and used limma-trend. v1 phospho added a
  phospho-minus-parent-protein correction; P3 has already built the raw 24M
  phosphosite DE and KSN kinase activity leaf.
- v1 integration found three durable axes: amyloid-driven activation
  (microglia/DAM/MHC), NLGF-linked synaptic suppression, and mixed-sign
  interaction/metabolic signals. P4 should rebuild the useful cross-modality
  view, not the 11-arc ledger / contest layer.
- v1 CCC triangulated CellChat + MultiNicheNet + LIANA and surfaced
  `Apoe_Trem2` and `App_Cd74` for the synaptic/clearance axis, but those tools
  were heavy/off-lock and used the full snRNAseq reference. In this rebuild,
  full CCC must be an explicit alternative, not a hidden default.

Current-method sweep (2026-07-02):
- `standR` is the current Bioconductor GeoMx workflow; its published workflow
  recommends limma-voom with `duplicateCorrelation` to account for repeated
  measurements / individuals. Source:
  https://academic.oup.com/nar/article/52/1/e2/7416375
- SpatialDecon remains the GeoMx-native deconvolution option and is available
  on the pinned Bioconductor 3.23 repos (`SpatialDecon` 1.22.0 locally
  available). It explicitly targets mixed-cell abundance estimation from GeoMx
  / spatial expression and supports custom profile matrices plus background
  estimation. Sources: https://www.nature.com/articles/s41467-022-28020-5,
  https://github.com/Nanostring-Biostats/SpatialDecon
- Proteomics SOTA alternatives are on-lock: `msqrob2`, `DEqMS`, `MSstats`,
  `MSstatsPTM`, `QFeatures`, and `PhosR` are available from the pinned Bioc
  repo. `msqrob2PTM` specifically distinguishes differential PTM abundance
  from differential PTM usage / parent-protein-corrected signal. Source:
  https://pubmed.ncbi.nlm.nih.gov/38154689/
- Default still stays limma-trend for the 24M bulk layer: the inputs are already
  Spectronaut-summary intensities, n=16, the existing helpers are warning-clean,
  and P4's question is cross-modality consistency rather than a new proteomics
  methods benchmark. `msqrob2PTM` becomes the reasoned alternative if the default
  phospho/protein-correction diagnostics are unstable.
- CCC SOTA has moved toward multi-condition frameworks (LIANA+, MultiNicheNet,
  CellChat v2), but local pinned-R availability check found `liana`, `CellChat`,
  and `multinichenetr` absent from CRAN/Bioc repos while LIANA+ is Python /
  scverse. Sources: https://www.nature.com/articles/s41556-024-01469-w,
  https://github.com/saeyslab/multinichenetr,
  https://liana-py.readthedocs.io/

Local availability probe:
`SpatialDecon`, `standR`, `msqrob2`, `DEqMS`, `PhosR`, `QFeatures`,
`MSstats`, and `MSstatsPTM` are on the pinned repos; `liana`, `CellChat`,
and `multinichenetr` are not.

## Default Design
Lean, on-lock, replicate-aware:

1. GeoMx DE is primary spatial transcriptomics.
   Use the raw `counts` layer with edgeR TMM + limma-voom. Model genotype
   effects with slide fixed effect and `duplicateCorrelation` block =
   `interaction(genotype, bio_rep)` to account for repeated AOIs from the same
   biological unit. Sensitivities: AOI-level slide-fixed fit without blocking
   and biological-unit-collapsed fit if the collapsed design is estimable. The
   report surfaces block-vs-unblocked concordance and downgrades claims that
   exist only in the unblocked AOI fit.

2. SpatialDecon is gated, not assumed.
   S1 only checks data/reference feasibility and pinned-repo availability.
   Install/load `SpatialDecon` in S3 only if the gate says a real
   deconvolution attempt is warranted. Deconvolution "earns it" if all are
   true: Q3-normalised GeoMx data and Q3-scaled negative background are usable,
   a reference profile can be built without a large full-object target, profile
   collinearity is below a predeclared threshold, at least broad-cell or
   microglia-vs-nonmicroglia abundance is estimable, and fresh target build
   stays warning-clean. If not, P4 records the failed precondition and reports
   GeoMx DE only.

3. Bulk proteome and phospho stay 24M, sample-key matched, no imputation.
   `proteome_de_24m`: exactly 16 matched runs, protein-group aggregation by sum,
   positive log2, median normalisation, prevalence filter, `fit_limma_log()`.
   `phospho_de_24m` from P3 remains the raw phosphosite layer. Add
   `phospho_corrected_24m`: phosphosite log2 intensity minus matched parent
   protein log2 intensity, refit limma-trend, carry matched-parent counts and
   missing-parent drop counts. Carry run-index sensitivity from P3 and add it
   for proteome/corrected phospho where full-rank.

4. CCC-lite is targeted to the synaptic/clearance axis.
   Default P4 does not install off-lock CellChat/MultiNicheNet/LIANA. Instead,
   build a focused `clearance_axis` table from measured evidence: Apoe, Trem2,
   App, Cd74, Pros1, Mertk, complement and synaptic marker sets across
   microglia RNA, GeoMx, proteome, and phospho/protein-corrected layers.
   Report it as ligand-receptor / clearance-axis support only if both sides of
   a pair are measured with coherent direction in at least two modalities or one
   modality plus a strong microglia anchor. Otherwise report "not earned" rather
   than inflating it to CCC.

5. Integration is descriptive and auditable.
   Build a harmonised symbol x contrast table with effect sizes, p/FDR, modality,
   feature provenance, and broad evidence class. Phosphosites aggregate to gene
   by max |t| / max |logFC| with site counts retained. Pathway integration reuses
   `mechanism_gene_sets` for GO/project sets and ranks by modality count,
   consistent sign, and named-axis membership. Nulls and mixed signs are
   first-class outcomes.

## Alternatives
Alternative A - Proteomics-heavy default (`QFeatures` + `msqrob2` /
`msqrob2PTM`).
Pros: most modern MS-specific modelling, explicit differential PTM usage.
Cons: substantially more dependency and object-surface complexity for n=16
Spectronaut-summary data; higher chance of spending P4 on workflow translation
rather than the cross-modality story. Good fallback if limma/protein-correction
diagnostics are unstable.

Alternative B - Full CCC stack now (LIANA+ / MultiNicheNet / CellChat).
Pros: closest to current CCC SOTA and v1's LR layer; better for APP/TREM2
claims if the user wants a dedicated communication arc. Cons: off-lock or Python
stack, likely requires a new full snRNAseq broad-cell reference target, and can
recreate v1 bloat. Default keeps this out unless targeted CCC-lite fails and the
synaptic/clearance axis becomes the load-bearing P4 question.

Alternative C - No deconvolution.
Pros: keeps P4 entirely on expression/protein/phosphosite statistics and avoids
heavy reference-building. Cons: loses the only tissue-composition readout outside
snRNAseq. Default is a gate: try only if the preconditions are clean.

Default choice: lean on-lock P4. It covers the backlog, respects the rebuild's
anti-bloat posture, and makes every off-lock or heavy branch opt-in.

## Steps
Each step is one closing unit. Resuming mid-plan: read this plan's Scope + your
step + `.agent/memory.md` relevant sections + files named in the step. Run
`scripts/check.sh` unless explicitly docs-only.

### S1 - GeoMx DE + deconvolution preflight [OPEN]
Add `R/crossmodality.R` and `tests/test_crossmodality.R` with GeoMx helpers.
Add `geomx_de` target. Do not install `SpatialDecon` in this step; record
pinned-repo availability and data/reference feasibility in the preflight
provenance, then let S3 decide whether to add the package and run it.

Contracts:
- `geomx_count_matrix(geomx)` reads RNA assay counts explicitly (the object
  defaults to SCT), coerces integer-ish counts,
  drops empty genes, and records non-integer residues.
- `geomx_meta(geomx)` returns rows aligned to AOIs with `genotype`, `slide`,
  `bio_unit = paste(genotype, bio_rep, sep=":")`, `roi`, `SampleID`, coordinates,
  Q3 factor, negative background, and nuclei sentinel counts.
- `fit_geomx_de()` primary = voom + TMM + duplicateCorrelation(block=bio_unit)
  + slide fixed effect + robust eBayes + 5 canonical contrasts. Sensitivities:
  unblocked slide-fixed fit and collapsed biological-unit fit if full-rank.
- Every top table carries `symbol`, `contrast`, `logFC`, `P.Value`,
  `adj.P.Val`, `t`, `CI.L`, `CI.R`, and provenance for primary/sensitivity.
- Decon preflight returns `status in {earned,defer,blocked}` plus precise reasons:
  background scale, nuclei sentinel, reference availability, profile collinearity
  if tested, and estimated memory footprint.

Acceptance:
- Synthetic tests cover count/meta alignment, slide design rank guard,
  duplicateCorrelation branch, sensitivity status, and malformed GeoMx metadata.
- Live `geomx_de` fresh build warning-clean with `tar_meta` clean.
- If decon status is not `earned`, the reason is specific enough for report text.
- `scripts/check.sh` green.

### S2 - Bulk proteome + corrected phospho [PENDING]
Add targets: `proteome_de_24m`, `phospho_corrected_24m`, and
`bulk_omics_summary`.

Contracts:
- Match exactly 16/16 sample-key runs for proteome and phospho; assert balanced
  4/genotype and identical sample order for phospho-protein correction.
- Proteome feature = `PG.ProteinGroups`; aggregate rows by sum on raw positive
  intensities, then log2, median-normalise, prevalence-filter, limma-trend.
  Carry first/unique gene symbol(s) and row-count-per-protein provenance.
- Corrected phospho subtracts matched parent protein log2 intensity per sample;
  drops sites without a matched parent with counts; re-applies prevalence filter.
- Add additive run-index sensitivity where rank permits; downgrade support if a
  key contrast flips or loses support under run-index adjustment.
- Reuse P3 raw `phospho_de_24m`; do not duplicate that target.

Acceptance:
- Unit tests cover protein aggregation, sample-order match, no-imputation policy,
  parent-protein correction, missing-parent counts, run-index sensitivity, and
  duplicate/multi-gene provenance.
- Fresh targets warning-clean; `tar_meta` clean.
- Live summary reports feature counts, missingness, significant counts, and
  Gsk3b/tau/synaptic anchors without hardcoded margins.
- `scripts/check.sh` green.

### S3 - Spatial composition + clearance-axis CCC-lite [PENDING]
If S1 decon preflight status is `earned`, add `geomx_decon` and
`geomx_abundance_de`. Always add `clearance_axis` as a focused measured-axis
table; it may return status `not_earned`.

Contracts:
- SpatialDecon reference path avoids storing a full 8G-derived object. If a full
  snRNAseq reference is needed, build only the compact profile matrix and record
  cell caps, common genes, cell-type labels, and memory/RSS.
- Q3-normalised GeoMx expression and Q3-scaled negative background must be on the
  same scale. Nuclei-rescaled absolute abundance stays disabled while `-1`
  sentinels remain.
- Abundance DE uses log(beta + offset) with slide fixed effect and bio-unit
  blocking/sensitivity where estimable; profile collinearity and unresolved
  abundance are reported.
- Clearance-axis table includes measured rows for Apoe/Trem2/App/Cd74/Pros1/Mertk,
  complement, synaptic GO/project sets, and modality support counts. It never
  calls CCC if evidence is single-sided or single-modality only.

Acceptance:
- Tests cover Q3 background scaling, profile-collinearity gate, abundance design
  shape, and clearance-axis support/not-earned classification.
- Fresh targets warning-clean; if decon is skipped, skip is intentional and
  reportable.
- `scripts/check.sh` green.

### S4 - Integrated gene/pathway divergence view [PENDING]
Add targets: `crossmodality_table`, `crossmodality_pathway`, and
`crossmodality_divergence`.

Contracts:
- Harmonise symbols across microglia RNA, GeoMx, proteome, phospho raw,
  phospho-corrected, pathway/TF/kinase summaries. Preserve feature provenance:
  protein groups, phosphosite counts, corrected/raw status, and source target.
- For each gene x contrast x modality, carry effect, p/FDR, statistic, sign,
  significance flag, and missingness reason.
- Pathway summary reuses `mechanism_gene_sets`; modalities contribute ranked
  statistics appropriate to their layer. Score rows by n_modalities_present,
  n_modalities_sig, n_consistent_sign, and axis membership, not by a sealed
  contest margin.
- Divergence view focuses on `interaction`, `nlgf_in_maptki`,
  `nlgf_in_p301s`, and `tau_in_nlgf`; report mixed-sign rows explicitly.

Acceptance:
- Tests cover duplicate symbol/site aggregation, missingness, sign-consistency
  scoring, pathway score invariants, and no silent drop of canonical contrasts.
- Fresh targets warning-clean and compact enough for report rendering.
- `scripts/check.sh` green.

### S5 - Cross-modality report chapter [PENDING]
Add `crossmodality_report` target and `_crossmodality.qmd`; include the chapter
after `_mechanism.qmd` in `index.qmd`.

Chapter shape:
1. GeoMx spatial DE: amyloid/DAM concordance and interaction caveats.
2. Bulk proteome + phospho: total-protein, raw phosphosite, corrected
   phosphosite, run-order caveats, kinase/protein correction reconciliation.
3. Spatial composition / clearance-axis status: decon result or explicit
   precondition failure; APP/TREM2/APOE/CD74 evidence if earned.
4. Integrated divergence: compact heatmap/table of concordant, mixed, and null
   cross-modality rows.
5. P4 synthesis: what cross-modality evidence strengthens, weakens, or leaves
   open for P5.

Contracts:
- `_crossmodality.qmd` loads only `crossmodality_report`; no heavy Seurat,
  GeoMx, proteomics, or phospho target in the render.
- All prose margins inline-computed. No v1 locked numbers or hardcoded p-values.
- Wording keeps modality scope honest: bulk hippocampus is not microglia-sorted;
  GeoMx AOIs are repeated spatial observations; deconvolution is abundance
  estimation, not direct cell counting.

Acceptance:
- `crossmodality_report_data()` guards every field the qmd reads.
- Full report renders with 0 warnings under `options(warn=2)`.
- `scripts/check.sh` green.
