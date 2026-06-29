# P0 Foundations - plan

Lay the reusable, reproducible foundation every downstream phase (P1-P5) builds on.
Fresh `main` rebuild; v1 foundation lives on branch `archive` (re-derive, do NOT
reproduce its locked numbers - R 4.6 re-baselines). Companions: `roadmap.md`
(direction), `memory.md` (contract), `archive_digest.md` + `git show archive:<path>`
(v1 reference). v1 foundation files mapped: `R/{constants,utils,io,design,de_pb,
de_sc,plot,report,microglia,helpers}.R`, `rmd/{01_data,02a_snrnaseq_qc}.Rmd`.

## Scope (P0 = foundations only)
IN: project-local reproducible env (R + Python); shared helpers (constants, utils,
io, design+contrasts, plot+report theme, loader); data load -> 4 analysis-ready
modalities; QC **sanity** pass (loads + sane metrics, light figures); 2x2 factorial +
5-contrast machinery; concrete runnable quality gate.
OUT (-> P1+): full snRNAseq microglia reprocess (SCT/Harmony/cluster-prune), microglia
substate assignment, single-cell DE (NEBULA/glmmTMB), pseudobulk DE *results*. P0
builds the contrast + pseudobulk *machinery* and smoke-tests it; it runs no analysis.

## Data (live, verified through `storage/data` symlink -> host Documents)
snrnaseq.rds 8.3G (Seurat; `broad_annotations=="Microglia"` subset) | geomx.rds 22M
(WTA spatial, ~91 ROIs) | proteomics_*.tsv 15M (peptide) | phosphoproteomics_*.tsv
**35.5M** (phosphosite) | proteomics_sample_key.csv 3.3K (67 rows; `File name` +
`Sample/Condtion` [sic]; TiO2/PTM phospho run map; **24M = rows 1-16** = 4 genotypes x
4 reps). Design: tau(MAPTKI|P301S) x amyloid(-/+NLGF) + batch; replicate =
genotype_batch (16 ids); 5 contrasts {tau_alone, nlgf_in_maptki, nlgf_in_p301s,
tau_in_nlgf, interaction}; `interaction = (NLGF_P301S-P301S)-(NLGF_MAPTKI-MAPTKI)`.

## STACK (gate resolved 2026-06-29)
**Quarto book** (reports) + **targets** DAG (`tarchetypes::tar_quarto`) + **rv** (R pkgs;
A2-ai Rust, declarative `rproject.toml` + `rv.lock`) + **uv** (Python), versions pinned
via **P3M dated snapshot**. User picked rv over renv (agent-native, mirrors uv; accepts
pre-1.0 + manual Bioc, cheap fallback). rv specifics: hand-wire Bioc **3.23** (R 4.6)
repos (no auto R<->Bioc coupling); exclude `rv/library` in `_quarto.yml` (bug #332).
SOTA rationale (cross-checked): Rmd = Posit maintenance-only; targets content-hash skip
pays off on the 8G object + multi-hour fits; rv = one declarative model across R+Python.

Cross-cutting (both stacks): pin versions via **P3M dated snapshot**. CRAN = binary
(trixie, ~7-8x faster); Bioc 3.23 = **likely source-only** via P3M for rv (the ~7-8x is
CRAN-only; prove the Bioc repo path at S1). Source fallback preserves **R-package
versions only** (system toolchain/libs separate -> needs Debian sysdeps: compilers, BLAS,
-dev headers; some packages may fail to build). Repro layers, NO Docker, NO bitwise
guarantee: targets=pipeline correctness; renv/rv-lock + uv.lock=versions; P3M=pinning.
Flagged for later: **BPCells** (Seurat-v5 on-disk, relieves 8G RAM ceiling) at S2/P1;
**qs2** `format="qs"` for targets serialization.

## Steps (each closeable in one ~200K window; acceptance gates the close)

### S1 - Reproducibility spine + env scaffold  [GATE-DEPENDENT; long installs]  -- DONE 2026-06-29
- `rv init` -> `rproject.toml` + `rv/library` + `rv/.gitignore` (NOT the lock; `rv sync`
  produces `rv.lock`). Repos in order: CRAN P3M dated snapshot (binary, trixie) + Bioc
  **3.23** (R 4.6) component repos (BioCsoft/ann/exp/workflows; likely **source-only** via
  P3M for rv -> prove at S1). Pin core stack: Seurat/SeuratObject, limma, edgeR,
  ggplot2/ggrepel/patchwork, BiocParallel, readr/dplyr/tibble/tidyr/stringr/purrr,
  targets/tarchetypes/qs2.
- uv: `uv init` + `.venv` + pyproject (minimal; no Python compute in P0, establish venv).
- Quarto: install pinned quarto-cli **project-local** (versioned dir + wrapper on a fixed
  path + checksum; NOT a system/PATH binary -> honour the no-system-install contract);
  `_quarto.yml` book skeleton (**exclude `rv/library`**).
- Orchestration: `_targets.R` + `R/` pure functions via `tar_source()` (DAG orders
  execution -> no manual dependency loader; supersedes v1's `helpers.R`).
- **Accept:** fresh `rv sync` + `uv sync` rebuild env; `library()` smoke-load of core
  stack passes; the **project-local** quarto wrapper reports the pinned version (not a
  PATH binary); `tar_manifest()` clean; lockfiles committed; any Debian sysdeps needed
  for source builds captured in an install log.
  **.gitignore** updated for Quarto (`_freeze/`, `.quarto/`, `/_book/`) + `_targets/` + `rv/`.

### S2 - constants + utils + io; data load  [HEAVY: 8G object, peak RAM ~9-10G]  -- DONE 2026-06-29
<!-- Deviations: load_proteomics+load_phospho merged -> read_spectronaut_tsv (identical bodies, KISS);
     dead args dropped (col_pattern, allow_raw); genotype_batch already PRECOMPUTED on the object
     (assert n_distinct==16 baked into load_snrnaseq, not derived). @misc$geneids aligns to RNA(33683)
     not the active SCT assay(28299). qs2 backend for format="qs"; trust_timestamps for the 8G input. -->

- `R/constants.R`: genotype_levels, genotype_colours (hex), contrast_definitions,
  canonical_microglia_markers, rbc markers, data paths.
- `R/utils.R`: `%||%`, `write_tsv_safe` (targets' store supersedes v1's `cache_or_run`).
- `R/io.R`: load_snrnaseq (readRDS -> subset Microglia -> drop SCT/reductions -> gc),
  build symbol_map (from `@misc$geneids`), load_geomx, load_proteomics, load_phospho,
  `proteomics_sample_meta(n_keep=16)`, match_intensity_columns, symbols_to_ensembl.
- Data-load target/module: materialise microglia_seurat_raw, symbol_map, geomx,
  proteomics, phospho, sample_key (store as `format="qs"`/`"file"` targets).
- **Accept:** each io fn returns documented-shape object vs LIVE data (Rscript
  smoke-test); microglia subset built; assert `dplyr::n_distinct(genotype_batch)==16`
  (one id per genotype x batch cell, no missing, all cells nonzero) -- not bare `table()` counts.

### S3 - design + contrasts + pseudobulk machinery  [GATE-INDEPENDENT; cheap]
- `R/design.R`: `factorial_design()` (~tau+nlgf+tau_nlgf[+batch]) + `make_contrast_matrix()`
  (cell-means form). BOTH must reproduce the same 5 named contrasts (key by name).
- `R/de_pb.R`: pseudobulk_counts, build_pseudobulk (replicate=genotype_batch),
  fit_limma_voom (min_count=5), fit_limma_log (proteomics/phospho), median_normalise,
  prevalence_filter (>=3 present in >=2 groups). Defer edgeR/deseq2/dream until a phase
  needs them (KISS).
- **Accept:** UNIT TEST - contrast matrix reproduces exact weights for all 5 (esp.
  `interaction = tau_nlgf`; `nlgf_in_p301s = nlgf+tau_nlgf`; etc.) on a synthetic 16-row
  meta; pseudobulk on S2 live object -> 16 genotype_batch columns.

### S4 - plot + report theme; QC-sanity chapter  [GATE-DEPENDENT engine; light render]
- `R/plot.R`: base ggplot theme fn + `concordance_plot()`.
- Report theme: Quarto `_brand.yml` / theme scss + IBM Plex (carry v1 bslib palette:
  primary/link = #B0344D). British English; human-facing prose uses hyphens not em/en-dashes.
- QC-sanity chapter (.qmd): load data; print dims + genotype x batch crosstab;
  QC metric distributions (nFeature/nCount/percent_mt/percent_contam from raw object's
  existing cols); 2-3 sanity figures. NO microglia reprocess (P1).
- **Accept:** chapter renders clean (0 error / 0 warning; enforced per the S5 gate),
  self-contained HTML; QC metrics within **predeclared numeric bounds asserted in-chapter**
  (declare bounds in code, no pre-privileged threshold, no eyeballed "sane").

### S5 - concrete quality gate + lock + close  [GATE-INDEPENDENT; light]
- Runnable gate (concrete): `scripts/check.sh` wraps `rv sync` + `uv sync` + `tar_make()`,
  then **enforces** zero-fault rather than trusting `tar_make()`'s exit (it returns 0 even
  with captured warnings): assert `tar_meta(fields=c("error","warnings"))` all NA across
  targets; set `options(warn=2)` in steps that tolerate it; grep the Quarto render log for
  `Warning`/`WARN`. Any error/warning/log-hit -> non-zero exit.
- Lock memory.md quality gate (provisional -> concrete, referencing the runnable gate);
  update map.md wiring (loader order, data->module flow, cache producer->consumer).
- **Accept:** gate runs green end-to-end on a clean state.

## Risks / caveats
- P3M trixie/Bioc binaries may be absent -> long source compiles at S1 (R-package
  versions stay pinned; system toolchain/libs are a separate apt concern -> may need sysdeps).
- 8G snRNAseq load: peak RAM ~9-10G; engineer gc; consider BPCells/qs at S2 if tight.
- targets x 8G object: use `format="qs"`/`"file"` + `memory="transient"`; heap may not
  return to OS (peak-RAM real, engineer it).
- rv: hand-wire Bioc 3.23 repos; `_quarto.yml` excludes `rv/library` (#332); prefer
  `rv remove` (exists in v0.22), hand-edit TOML only for what rv can't express; pre-1.0
  -> if it bites, fall back to renv.

## Close-out (when all S done)
Adversarial review (correctness/scope/token-efficiency/obsolescence); fold decisions
into history.md; archive plan -> `.agent/completed/p0_foundations_plan_YYYY-MM-DD.md`;
reset roadmap Active plan to (none); commit `<scope> (p0 close): ...`.
