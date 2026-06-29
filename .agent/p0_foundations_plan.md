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
     not the active SCT assay(28299). qs2 backend for format="qs"; trust_timestamps for the 8G input.
     Review (codex) hardening: refactor-then-assert + 4x4 bijection check (load_snrnaseq); fail-loud
     geomx genotype; read_spectronaut_tsv stop()s on parse problems + strips readr attrs (qs bad_weak_ptr);
     shared normalise_ptm_stub; balanced/unique-key asserts (proteomics_sample_meta); unique-symbol assert
     (build_symbol_map). Deferred: col-map DAG validation target -> P4; committed io contract tests -> S3. -->

- `R/constants.R`: genotype_levels, genotype_colours (hex), contrast_definitions,
  canonical_microglia_markers, rbc markers, data paths.
- `R/utils.R`: `%||%`, `write_tsv_safe` (targets' store supersedes v1's `cache_or_run`).
- `R/io.R`: load_snrnaseq (readRDS -> subset Microglia -> drop SCT/reductions -> gc),
  build symbol_map (from `@misc$geneids`), load_geomx, read_spectronaut_tsv (merged load_proteomics+load_phospho),
  `proteomics_sample_meta(n_keep=16)`, match_intensity_columns, symbols_to_ensembl.
- Data-load target/module: materialise microglia_seurat_raw, symbol_map, geomx,
  proteomics, phospho, sample_key (store as `format="qs"`/`"file"` targets).
- **Accept:** each io fn returns documented-shape object vs LIVE data (Rscript
  smoke-test); microglia subset built; assert `dplyr::n_distinct(genotype_batch)==16`
  (one id per genotype x batch cell, no missing, all cells nonzero) -- not bare `table()` counts.

### S3 - design + contrasts + pseudobulk machinery  [GATE-INDEPENDENT; cheap]  -- DONE 2026-06-29
<!-- Deviations: both design fns mirror v1 (algebra unchanged); canonical contrast column order
     (tau_alone, nlgf_in_maptki, nlgf_in_p301s, tau_in_nlgf, interaction). de_pb: dropped the dead
     `group` arg (fit_limma_voom filters via design; fit_limma_log never used it); edgeR
     calcNormFactors -> normLibSizes (4.x); pseudobulk_counts via Matrix::rowSums; build_pseudobulk
     asserts covariate-constancy within sample + 1:1 meta<->count alignment. edgeR-QLF/DESeq2/dream
     deferred (KISS, P1+). Fail-loud guards (out-of-level genotype, absent level column, batch
     requested but batch col absent -> batch-less modalities pass add_batch=FALSE); deterministic
     radix sort for pseudobulk/prevalence group order. Tests = stopifnot harness (no testthat) in tests/ + tests/helpers.R; the io contract tests
     deferred from S2 review landed here (test_io.R). Acceptance met: test_design (exact weights +
     factorial==cell-means equivalence (estimator-map equality -> proven for ANY response) +
     hand-computed contrast values); live microglia_seurat_raw (26104 cells) -> build_pseudobulk =
     16 genotype_batch columns. Review (codex) hardening: fail-loud input contracts on the de_pb
     fitters (name alignment + full rank) + pseudobulk_counts/build_pseudobulk (cell alignment,
     no-NA covariates) + prevalence_filter (no-NA group) + factorial_design (full rank, no-NA
     batch); expect_error gains a message-pattern check; pseudobulk sum verified for all 16 groups. -->
- `R/design.R`: `factorial_design()` (~tau+nlgf+tau_nlgf[+batch]) + `make_contrast_matrix()`
  (cell-means form). BOTH must reproduce the same 5 named contrasts (key by name).
- `R/de_pb.R`: pseudobulk_counts, build_pseudobulk (replicate=genotype_batch),
  fit_limma_voom (min_count=5), fit_limma_log (proteomics/phospho), median_normalise,
  prevalence_filter (>=3 present in >=2 groups). Defer edgeR/deseq2/dream until a phase
  needs them (KISS).
- **Accept:** UNIT TEST - contrast matrix reproduces exact weights for all 5 (esp.
  `interaction = tau_nlgf`; `nlgf_in_p301s = nlgf+tau_nlgf`; etc.) on a synthetic 16-row
  meta; pseudobulk on S2 live object -> 16 genotype_batch columns.

### S4 - plot + report theme; QC-sanity chapter  [GATE-DEPENDENT engine; light render]  -- IN PROGRESS (resume here)
<!-- S4 STATUS 2026-06-29 (mid-step; UNCOMMITTED WIP in tree; compaction handoff):
  DECISION (user-picked): report = ONE self-contained OFFLINE HTML file, NOT a Quarto book.
    Verified Quarto fact: `embed-resources` STRIPS the brand/bootstrap theme in a `book` project
    (+ a "Could not fetch resource ./<chapter>.html" nav WARNING); it applies cleanly without a book.
    -> drop the book; render ONE standalone `format: html` doc (embed-resources works there); modular
    sources via `{{< include _section.qmd >}}` (leading `_` => not rendered standalone); theme via
    `theme.scss` (Bootstrap Sass vars) NOT `_brand.yml`; IBM Plex woff2 BUNDLED LOCALLY (offline;
    embed inlines them) NOT brand `source: google` (= a view-time Google-CDN dep).
  DONE + working: `R/plot.R` (smoke-tested 0-warn, namespace-qualified => sources w/o ggplot2 attached:
    theme_tau(base_family="" -> no missing-font warning) + scale_colour/fill_genotype() [+ scale_color_
    alias] + concordance_plot() [v1 port, ends theme_tau(), geom_smooth formula spelt out to mute its
    msg]). `qc.qmd` RENDERS 0-error, all 17 chunks run, every in-chapter sanity bound passes (bounds are
    definitional rails: counts>=1 & <1e6, features<=nCount<=n_genes, %s in 0-100, +structural 4x4-nonzero
    /canonical-levels/16 ids). Targets store already built => reads microglia_seurat_raw (33683 x 26104)
    from qs, NO 8G reload.
  WIP TREE (uncommitted): R/plot.R (keep) | qc.qmd (keep content -> RENAME `_qc.qmd`, becomes the include)
    | _brand.yml (SUPERSEDED -> delete) | _quarto.yml (edited for book -> REWORK). _book/ artefacts gitignored.
  REMAINING (fresh session, in order):
    1. Bundle IBM Plex woff2 -> `assets/fonts/` (latin: Sans 400/500/600/700+400i, Serif 600/700, Mono
       400/600) from a PINNED fontsource ver via jsDelivr; record the version; git-add (small binaries).
    2. `theme.scss`: scss:defaults {$primary,$link-color=#B0344D; $code-color=#3F5A6B; $font-family-base/
       $headings-font-family/$font-family-monospace = IBM Plex Sans/Serif/Mono} + scss:rules {@font-face
       url(assets/fonts/*.woff2)}. Delete `_brand.yml`.
    3. RENAME qc.qmd -> `_qc.qmd` (its setup already does library + tar_source + tar_load + theme_set).
       `index.qmd` = the ONE report: front-matter format.html {theme: theme.scss, embed-resources: true,
       code-fold: true, code-tools: true} + overview prose + `{{< include _qc.qmd >}}`.
    4. `_quarto.yml`: drop `type: book` + `book:` block -> `project:{type: default, render: [index.qmd]}`;
       move format.html opts (or keep in index front-matter). output-dir: rename `_book`->`_report` AND
       update .gitignore (`/_book/`->`/_report/`) + .claude/settings.json deny (`Read(./_book/**)`->
       `Read(./_report/**)`) -- deny globs "never trust by eye": Read-test 1 must-block + 1 must-read after.
       (Cheaper alt: keep the `_book` dir name; decide.)
    5. `_targets.R`: `tar_quarto(book, path=".")` -> `tar_quarto(report, path=".")`. VERIFY tar_quarto
       detects the `tar_load` INSIDE included `_qc.qmd` (tar_manifest/tar_network); else add deps explicitly.
    6. Render (pinned quarto): assert offline (0 external http loads), #B0344D + IBM Plex + @font-face
       data-URIs INLINED, 0 WARNING in log, bounds stopifnot pass. Inspect output via runtime-indirection
       python (output dir is deny-Read: path must be ABSENT from the command text).
    7. `tests/test_plot.R` (device-free: source constants+utils+plot; class-check theme_tau/scale_*_genotype/
       concordance_plot, NO draw) so S5 check.sh loops it. Run all tests/test_*.R green.
    8. Update map.md (plot.R + report wiring: index.qmd -include-> _qc.qmd; theme.scss + assets/fonts;
       report target). Fold the report-architecture DECISION into memory.md (durable, outlives P0).
       Mark S4 done (plan + roadmap). Commit `report (p0 s4): ...`.
-->
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
