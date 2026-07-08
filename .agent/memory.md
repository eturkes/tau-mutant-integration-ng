# Memory - standing contract

Codex-only repo. Canonical instructions = `AGENTS.md`; session entry = `$session-prompt`
backed by `.codex/prompts/session.md`; repo skill = `.agents/skills/session-prompt/`.
Use `roadmap.md` for trajectory/history, `map.md` for live wiring, `history.md` for
decision digests, `archive_digest.md`/branch `archive` for v1 mining.

## Current Scope

Goal: integrate snRNAseq + GeoMx spatial + 24M proteome + 24M phosphoproteome across
4 AD mouse genotypes. Design = 2x2 tau x amyloid:
- MAPTKI = tau-KO baseline
- P301S = mutant tau
- NLGF_MAPTKI = amyloid alone
- NLGF_P301S = amyloid + mutant tau

Canonical interaction = `(NLGF_P301S - P301S) - (NLGF_MAPTKI - MAPTKI)`.

Live report scope (2026-07-08): 17 figures, 3 qmd chapters, expected 29 targets. Infrastructure
that does not directly feed the final analysis document is removed: committed tests, Python/uv
files, composition/sccomp/CmdStan target, per-substate pseudobulk, prose inventory, stageR
layer, mechanism/crossmodality/qc/story chapters and modules. The GeoMx QC atlas,
normalization/RLE, ordination, gene-detection, sample-heatmap, spatial-program, and contrast-diagnostic figures are compact
figures riding on `geomx_de`, not a restored QC chapter. Historical claims remain in git + `roadmap.md`;
do not treat them as live pipeline contracts.

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
- Proteome/phospho TSVs: Spectronaut PTM exports. Current report uses the 24M 16-run
  subset from `proteomics_sample_key.csv`.

## Scientific Spine

Durable headline:
- Amyloid drives microglia homeostatic -> DAM activation.
- Mutant tau modulates amyloid response mostly through DAM-cell composition, not a
  supported progression-beyond-composition shift along the activation trajectory.
- Non-snRNAseq modalities provide context figures, not resurrected mechanism or
  cross-modality chapters.

Substates: Homeostatic, DAM, IFN, Proliferative. Current coherent clusters contain
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

Microglia reprocess:
- SCT-v2/glmGamPoi on RNA counts, regress percent_mt + percent_contam.
- Harmony integrates batch only; genotype/amyloid are biology, not correction terms.
- Louvain clusters at multiple resolutions; primary stable column = `microglia_clusters`.
- Re-run is seed-deterministic up to tolerance, not bitwise.

Microglia annotate:
- UCell scoring on SCT data; cluster-level substate assignment is authoritative.
- Contaminant pruning uses raw identity-vs-contaminant evidence, not z-scaled substate
  scores.
- Prune asymmetry exists: dropped low-quality clusters are enriched for NLGF_MAPTKI but
  are not DAM-high. Report caveat factually.

Microglia composition:
- No live composition inference target. Replicate-unit stacked bars are derived inside
  `microglia_report_data()` from annotated Seurat metadata and emitted as
  `unit_composition`.

Microglia DE:
- Live target = `pb_de_microglia` only.
- Raw RNA pseudobulk by genotype_batch, limma voom with quality weights, topTables for
  the four-method amyloid-response scatter.
- No stageR, per-substate DE, MDE, or DAM-direction helper in the live DAG.

Trajectory:
- Slingshot on Harmony embedding with forced Homeostatic -> DAM lineage.
- IFN/Proliferative are omitted from the lineage frame via flags, not deleted from the
  biological record.
- Primary inference uses 16-unit pseudotime summaries plus Kitagawa-style decomposition;
  glmmTMB per-cell model is supportive and degrades to recorded failure instead of
  blocking the report.

Modality context:
- `R/modality_de.R` restores only primary DE needed for report figures:
  GeoMx voom/TMM with slide fixed effect + duplicateCorrelation; proteome/phospho
  limma-trend on log2 median-normalized 24M intensities.
- `geomx_de$qc` is descriptive only: AOI library/detected-gene/nuclei/area/background/Q3
  flags are visual diagnostics and exclude no AOIs.
- `geomx_de$normalization` is descriptive only: raw/TMM logCPM quantiles, TMM RLE
  quantiles, Q3 factor vs negative-control background, and saved voom trend. It excludes no
  AOIs and does not alter the primary GeoMx model.
- `geomx_de$ordination` is descriptive only: deterministic TMM-logCPM PCA/classical MDS on
  top-variable filter-passing genes plus PC1/PC2 top loadings. It excludes no AOIs and does
  not alter the primary GeoMx model.
- `geomx_de$gene_detection` is descriptive only: mean TMM logCPM vs raw-count detection
  fraction, existing `filterByExpr(min.count=5)` status, microglia identity/homeostatic/DAM
  marker measurability, and highest-detected genes. It excludes no AOIs and adds no new
  gene filter beyond the primary GeoMx model.
- `geomx_de$sample_heatmap` is descriptive only: top-variable filter-passing TMM-logCPM
  genes, row-z heatmap clipped for display, deterministic average-linkage AOI/gene
  clustering, genotype/slide/segment/bio-unit/ROI tracks, and the existing signed GeoMx
  amyloid-response score. It excludes no AOIs and changes no DE model.
- `geomx_de$spatial_programs` is descriptive only: coordinate-only AOI maps for
  Homeostatic/DAM/IFN/MHC_APC signatures plus Apoe/Trem2 single genes, scored as
  TMM-logCPM row-z summaries over filter-passing genes. Tissue images are absent from
  the live report path, so the figure must say coordinate-only. It excludes no AOIs and
  changes no DE model.
- `geomx_de$contrast_diagnostics` is descriptive only: volcano/MA rows, deterministic top
  labels, signed support counts, and an interaction top-gene panel from the existing
  primary GeoMx top tables. It adds no inference model, excludes no AOIs, and changes no
  DE model.
- Auxiliary deconvolution, run-index sensitivity, and broad mechanism/cross-modality
  target families stay deleted.

Report:
- `_microglia.qmd`, `_trajectory.qmd`, `_modality.qmd`; no body prose beyond headings,
  captions, alt text, and code display controls.
- qmd chunks set `options(warn=2)`. Any warning during render is a failure.
- Compact report targets must keep qmd renders off the 612MB annotated Seurat object.

## Environment And Gate

Live stack = rv-managed R + project-local Quarto + targets. No committed Python/uv
surface unless a future report-producing step earns it.

Fresh bootstrap:
1. `scripts/install-sysdeps.sh`
2. `scripts/install-rv.sh`
3. `scripts/install-quarto.sh`
4. `rv sync`
5. `scripts/check.sh`

`scripts/check.sh`:
- runs `rv sync` unless `CHECK_SKIP_SYNC=1`
- invalidates and rebuilds `report`
- scans `targets::tar_meta()` warnings/errors/messages
- scans render log for warning/error anchors

Use `CHECK_SKIP_SYNC=1 scripts/check.sh` after a known-good sync for faster local iteration.

## Maintenance

Keep `.gitignore` aligned with generated artefacts. Current deleted infrastructure should
stay absent unless it directly accelerates or protects the final report path.

When adding a dependency, prefer current project-local R/Quarto path first; update
`rproject.toml` + `rv.lock` together. Avoid global library/tool leakage.

Instruction/prompts are active source:
- durable project fact -> `.agent/memory.md`
- live wiring -> `.agent/map.md`
- trajectory/ledger -> `.agent/roadmap.md`
- workflow prompt change -> `.codex/prompts/session.md` and matching skill docs.
