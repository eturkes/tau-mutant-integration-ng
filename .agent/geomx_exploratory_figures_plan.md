# GeoMx exploratory figures plan

User task 2026-07-08: web-search exploratory figures for GeoMx, write them down, and
update roadmap so each figure is implemented in its own session.

## Source Digest

Primary/current sources searched 2026-07-08:
- Bruker GeoScript Hub:
  https://brukerspatialbiology.com/products/geomx-digital-spatial-profiler/geoscript-hub/
  - GeoMx-validated script classes include SpatialDecon, Dimension Reduction, Volcano Plot,
    Cell-Type Contouring, General Heatmap, SpatialOmicsOverlay.
  - Dimension Reduction outputs PCA/tSNE/UMAP segment scatter plus PCA variance/loadings.
  - General Heatmap supports labelled publication heatmaps; SpatialOmicsOverlay overlays
    gene counts, annotations, pathway scores, and cell-type proportions on GeoMx images.
- Bruker GeoMx DSP Data Analysis User Manual:
  https://www.brukerspatialbiology.com/wp-content/uploads/2022/06/MAN-10154-01-GeoMx-DSP-Data-Analysis-User-Manual.pdf
  - NGS QC uses line plot + study summary over read-processing outcomes; segment QC checks
    nuclei/area and warning-tagged segments; PCA is a core visualization for segment
    similarity/clustering and target loadings.
- Davis Lab standR workflow:
  https://davislaboratory.github.io/GeoMXAnalysisWorkflow/articles/GeoMXAnalysisWorkflow.html
  - Recommended GeoMx WTA workflow uses SpatialExperiment, ROI QC, gene QC, RLE, PCA,
    pairwise PCA, batch/biology evaluation stats, limma-voom, MA plots, and downstream
    visualization.
- standR manual:
  https://bioconductor.posit.co/packages/3.23/bioc/manuals/standR/man/standR.pdf
  - Relevant figure functions: `drawPCA`, `plotMDS`, `plotPairPCA`, `plotPCAbiplot`,
    `plotGeneQC`, `plotRLExpr`, `plotROIQC`, `plotSampleInfo`, `plotScreePCA`.
- SpatialDecon docs:
  https://github.com/Nanostring-Biostats/SpatialDecon
  and
  https://bioconductor.statistik.tu-dortmund.de/packages/3.18/bioc/vignettes/SpatialDecon/inst/doc/SpatialDecon_vignette_NSCLC.html
  - SpatialDecon estimates mixed cell abundance for GeoMx; exploratory displays include
    abundance/proportion bars, `florets` spatial cockscomb plots, collapsed abundance
    heatmaps, reverse-decon residual heatmaps, and fit-dependency scatter.

Local seam:
- Live report has eight GeoMx figures after S7: `fig-geomx-qc-atlas` = AOI QC atlas;
  `fig-geomx-normalization-rle` = raw/TMM logCPM, RLE, Q3/background, voom trend;
  `fig-geomx-ordination` = slide-faceted PCA/MDS + scree + PC1/PC2 loadings;
  `fig-geomx-gene-detection` = WTA gene detectability, existing filterByExpr decision,
  microglia marker measurability, highest-detected genes;
  `fig-geomx-sample-heatmap` = clustered top-variable-gene row-z heatmap with
  genotype/slide/segment/bio-unit/ROI tracks and signed-response score;
  `fig-geomx-spatial-program-overlays` = coordinate-only biology program maps;
  `fig-geomx-contrast-diagnostics` = volcano/MA/support diagnostics over the five
  canonical GeoMx contrasts;
  `fig-modality-geomx-landscape` = slide-faceted AOI signed-score map + genotype score
  distribution + top score-gene drivers.
- `geomx_de$spatial` already carries compact AOI fields:
  `slide`, `roi`, `segment`, `sample_id`, `genotype`, `aoi_area`, `x_coord`, `y_coord`,
  `q3_factor`, `neg_background`, `nuclei`, `signed_response_score`.
- `geomx_de$qc` carries compact S1 QC data: long AOI metrics, per-slide/segment flag counts,
  thresholds, and provenance. `geomx_de$normalization` carries compact S2 raw/TMM
  distribution quantiles, TMM RLE quantiles, Q3/background AOI data, and saved voom trend.
  `geomx_de$ordination` carries compact S3 PCA/MDS/scree/loading data.
  `geomx_de$gene_detection` carries compact S4 gene-detection/marker-measurability data.
  `geomx_de$sample_heatmap` carries compact S5 clustered top-variable-gene heatmap data.
  `geomx_de$spatial_programs` carries compact S6 coordinate-only biology program data.
  `geomx_de$contrast_diagnostics` carries compact S7 volcano/MA/support data. Report
  chunks read these descriptors through `modality_scatter_figures`.
- Current lean DAG removed deconvolution targets. Any decon figure must either re-earn a
  compact decon/status target in its own session or render a truthful blocked/feasibility
  diagnostic, not an abundance claim.

## Default Plan: One Figure Per Session

Acceptance shared by every session:
- Adds exactly one new visible GeoMx exploratory figure to `_modality.qmd` or a GeoMx-only
  include, with one `fig-*` id, caption, and nonblank `fig-alt`.
- Reads compact GeoMx data targets in report chunks; heavy `geomx` object access stays in
  target builders.
- Preserves current report claims unless the figure reveals a concrete issue; any
  changed claim is target-derived and recorded in roadmap.
- Runs a focused target/render check; run full `scripts/check.sh` when code touches shared
  helpers, target wiring, or report-wide behavior.

S1 [DONE 2026-07-08] - `fig-geomx-qc-atlas`
- Purpose: AOI/segment QC first-look figure.
- Panels: library size / detected genes / nuclei / area / negative background / q3 factor,
  stratified by slide, genotype, segment, and warning-like sentinel states.
- Source basis: Bruker NGS line/summary + segment QC; standR `plotROIQC`.
- Local data: counts from `geomx_count_matrix()`, metadata from `geomx_meta()`.
- Acceptance: flags low-depth / nuclei-sentinel / small-area AOIs without excluding samples
  or changing DE.
- Status: landed as compact `geomx_de$qc` + visible `_modality.qmd` figure. Live flags =
  53/91 AOIs with >=1 descriptive flag: 42 nuclei sentinels, 5 low-library, 5 small-area,
  5 high-background, 7 high-Q3, 0 low-gene flags because detected genes are constant. No
  AOI exclusion, DE change, or report-claim change.

S2 [DONE 2026-07-08] - `fig-geomx-normalization-rle`
- Purpose: normalization/background sanity figure.
- Panels: raw logCPM distribution vs TMM/logCPM distribution, RLE by slide/genotype, q3
  factor vs negative background, optional voom mean-variance trend if cheap.
- Source basis: standR `plotRLExpr`; Bruker normalization/negative-control QC.
- Acceptance: makes slide-driven technical spread visible and states model still uses
  limma-voom with slide fixed effect + duplicateCorrelation.
- Status: landed as compact `geomx_de$normalization` + visible `_modality.qmd` figure. Live
  data = 91 AOIs, 19,959/19,963 genes kept by `filterByExpr(min.count=5)`, raw/TMM per-AOI
  logCPM quantiles, TMM RLE quantiles, Q3/background AOI scatter, saved voom trend from the
  primary fit; Q3 factor vs negative-control background Spearman rho = 0.994. No AOI
  exclusion, DE change, or report-claim change.

S3 [DONE 2026-07-08] - `fig-geomx-ordination`
- Purpose: sample similarity / biology-vs-slide ordination figure.
- Panels: PCA or MDS AOI scatter colored by genotype and shaped/faceted by slide/segment,
  variance explained, top loading genes for PC1/PC2.
- Source basis: Bruker Dimension Reduction; standR `drawPCA`, `plotMDS`,
  `plotPCAbiplot`.
- Acceptance: precomputes deterministic PCA/MDS inside target data; report draw is pure.
- Status: landed as compact `geomx_de$ordination` + visible `_modality.qmd` figure. Live
  data = 91 AOIs, 19,959/19,963 genes kept by `filterByExpr(min.count=5)`, 2,000
  top-variable genes used for TMM-logCPM PCA/MDS; PC1 = 19.36%, PC2 = 6.94% variance
  explained. No AOI exclusion, DE change, or report-claim change.

S4 [DONE 2026-07-08] - `fig-geomx-gene-detection`
- Purpose: gene-level QC/detectability figure.
- Panels: mean expression vs detection fraction, low-coverage filter boundary, labelled
  high-detection microglia/DAM/homeostatic marker genes, and top detected-driver summary.
- Source basis: standR `plotGeneQC`; GeoMx WTA negative-probe handling.
- Acceptance: documents which marker/pathway genes are measurable before interpretation.
- Status: landed as compact `geomx_de$gene_detection` + visible `_modality.qmd` figure.
  Live data = 91 AOIs, 19,959/19,963 genes pass `filterByExpr(min.count=5)`, 4
  low-coverage genes, marker genes present/pass = Microglia 8/8, Homeostatic 10/10,
  DAM 18/18. The figure documents measurability only; no AOI exclusion, DE change, or
  report-claim change.

S5 [DONE 2026-07-08] - `fig-geomx-sample-heatmap`
- Purpose: sample/gene structure heatmap.
- Panels: AOI correlation or top-variable-gene z-score heatmap with annotation tracks for
  genotype, slide, segment, bio-unit/ROI, and signed-response score.
- Source basis: Bruker General Heatmap; standard GeoMx exploratory clustering.
- Acceptance: keeps feature rows capped and precomputed so embedded report size stays lean.
- Status: landed as compact `geomx_de$sample_heatmap` + visible `_modality.qmd` figure.
  Live data = 91 AOIs, 19,959/19,963 genes pass `filterByExpr(min.count=5)`, 40
  top-variable genes displayed as row-z scores clipped at +/-2.5, with deterministic
  average-linkage AOI/gene clustering. The figure documents sample/gene structure only;
  no AOI exclusion, DE change, or report-claim change.

S6 [DONE 2026-07-08] - `fig-geomx-spatial-program-overlays`
- Purpose: spatial small-multiple overlay beyond the existing single signed-score plate.
- Panels: AOI coordinate maps for a compact set of biology-first scores/genes, e.g. DAM
  score, homeostatic score, IFN score, Apoe/Trem2 or top amyloid-response drivers.
- Source basis: Bruker SpatialOmicsOverlay concept; local coordinate-only fallback because
  OME-TIFF images are not in the live repo.
- Acceptance: explicitly labels coordinate-only status if tissue images are unavailable.
- Status: landed as compact `geomx_de$spatial_programs` + visible `_modality.qmd`
  figure. Live data = 91 AOIs x 6 programs (546 AOI-program rows), 19,959/19,963
  genes pass `filterByExpr(min.count=5)`, scored features = Homeostatic 10/10,
  DAM 18/18, IFN 12/12, MHC/APC 7/7, Apoe 1/1, Trem2 1/1. Coordinate-only status is
  explicit because tissue images are not in the live report path. No AOI exclusion, DE
  change, or report-claim change.

S7 [DONE 2026-07-08] - `fig-geomx-contrast-diagnostics`
- Purpose: GeoMx-only DE diagnostic figure.
- Panels: five canonical contrast volcano/MA small multiples, top labels, signed support
  counts, and interaction emphasis.
- Source basis: Bruker Volcano Plot; standR limma-voom + MA visualization.
- Acceptance: uses existing `geomx_de$primary$top`; no new inference model.
- Status: landed as compact `geomx_de$contrast_diagnostics` + visible
  `_modality.qmd` figure. Live data = 99,795 contrast-feature rows
  (19,959 genes x 5 canonical contrasts), 20 deterministic labels, and
  FDR <= 0.10 support counts: tau_alone 4 up/0 down; nlgf_in_maptki
  5,258 up/1,559 down; nlgf_in_p301s 1,386 up/541 down; tau_in_nlgf
  126 up/270 down; interaction 45 up/117 down. The figure reuses the
  existing primary GeoMx top tables; no new inference model, AOI exclusion,
  DE change, or report-claim change.

S8 - `fig-geomx-roi-segment-replicates`
- Purpose: ROI/segment/repeated-observation structure figure.
- Panels: paired ROI/segment score differences when repeated segments exist; otherwise
  bio-unit/slide support map, AOI-per-unit counts, and duplicateCorrelation block audit.
- Source basis: GeoMx segment/ROI model; Bruker line/summary visualizations.
- Acceptance: clarifies pseudo-replication risk and why bio-unit blocking is load-bearing.

S9 - `fig-geomx-decon-feasibility`
- Purpose: deconvolution/status exploratory figure without overstating abundance.
- Panels: reference-overlap/coverage, blocked-AOI beta-total or nuclei-sentinel map, and
  residual/fit QC; if decon is re-earned in-session, add abundance proportions/florets.
- Source basis: Bruker SpatialDecon/SpatialOmicsOverlay; SpatialDecon `TIL_barplot`,
  `florets`, collapsed abundance heatmap, reverse-decon residual diagnostics.
- Acceptance: if live decon remains absent or blocked, the figure is a blocked diagnostic
  and says so in caption; no cell-abundance claim.

## Alternative

Bundle these into 2-3 dashboard sessions (`QC+normalization`, `ordination+heatmaps`,
`DE+spatial+decon`). Rejected for this roadmap because the user asked for one session per
figure and prior figure work in this repo has overflowed when sessions bundle too many
plot contracts.

## Active-State Notes

Default execution order is QC -> normalization -> ordination -> gene detectability ->
heatmap -> spatial overlays -> DE diagnostics -> replicate structure -> decon/status.
Later sessions may drop or merge a figure only after an implemented predecessor proves a
figure redundant or unsupported by available local data; record that decision in roadmap.
