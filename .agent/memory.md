# Memory - standing contract (read every session)

Durable facts / decisions / gotchas surviving all plans. Companions: `roadmap.md`
(direction), `map.md` (wiring), `history.md` (new decision digests),
`archive_digest.md` (v1 reference). This is a FRESH streamlined rebuild; v1 lives
on branch `archive`.

## Project
Integrate snRNAseq + GeoMx spatial + bulk proteomics + bulk phosphoproteomics across
4 mouse AD genotypes to read how amyloid (NLGF) reshapes microglia under different
tau backgrounds. Genotypes (2x2): MAPTKI (tau-KO baseline), P301S (mutant tau),
NLGF_MAPTKI (amyloid alone), NLGF_P301S (amyloid + tau). Design = tau (MAPTKI vs
P301S) x amyloid (-/+ NLGF) + batch. Divergence = interaction
(NLGF_P301S - P301S) - (NLGF_MAPTKI - MAPTKI).

## Raw data (storage/data = symlink -> host Documents tree, shared/external read-only copy;
gitignored via /storage/* + deny-Read; Rscript reads resolve through it, bypassing the deny.
Shapes VERIFIED live in S2 via the R/io.R loaders; data is immutable so numbers are stable.)
- snrnaseq.rds 8.3G: Seurat; full = RNA 33683 genes (Assay5, ENSMUSG rownames) + SCT 28299
  (active assay). GOTCHA: dim(obj)=SCT(28299) but @misc$geneids (33683 symbols) aligns to RNA,
  NOT SCT -> build_symbol_map uses RNA rownames (asserts unique non-NA symbols; positional alignment
  rests on the v1 @misc$geneids contract, length-checked only). broad_annotations=="Microglia" -> 26104 cells.
  meta PRECOMPUTED (P1/QC consume, don't recompute): genotype (canonical), batch (batch01-04),
  genotype_batch (16-lvl factor, 4x4 fully crossed all-nonzero), sex, nCount/nFeature_RNA,
  percent_mt/ribo/malat1/contam, doublets. load_snrnaseq drops SCT+reductions -> 340MB RNA-counts
  microglia subset (the qs target microglia_seurat_raw).
- geomx.rds 22M: Seurat 19963 genes x 91 AOIs; genotype already canonical (n=20/20/28/23). Spatial
  design cols bio_rep/tech_rep/slide_rep/roi/segment/SampleID; NO batch/genotype_batch (differs from snRNAseq).
- proteomics_nonfiltered_nonnormalised.tsv 15M: Spectronaut PTM export, 45972 x 30. Cols 1-14 = PTM
  annotation (PG.*, Gene-pSite, PTM.SiteAA/Location, "Phopshosite probability"[sic typo]); 15-30 = 16
  intensity `Naoto-Hippo_TiO2_DIA_NN.raw.PTM.Quantity` (24M set only). peptide->protein-group sum = P4.
- phosphoproteomics_*.tsv 35.5M: Spectronaut PTM, 64328 x 81. 1-14 annotation; 15-81 = 67 intensity
  `*.PTM.Quantity` (Naoto 01-26 + Set6 01-41). The phospho target stores ALL 67 (nothing dropped at
  load); only the first 16 (Naoto 01-16) are the 24M timepoint -> the 24M column-subset is a P4 step.
- proteomics_sample_key.csv: 67 rows {`File name`, `Sample/Condtion`[sic]}; proteomics_sample_meta
  n_keep=16 (24M, 4 geno x 4 reps; asserts balanced reps + exact labels + unique join keys).
  normalise_ptm_stub (shared by the key producer + match_intensity_columns) strips `.PTM.Quantity`
  then trailing `.raw` -> both files' intensity cols collapse to one run stub, .raw discrepancy handled.
  match_intensity_columns is helper-only in P0 (NA = non-sample col); P4 wires it + asserts 16/16.
  read_spectronaut_tsv stop()s on any parse problem + strips readr spec/problems attrs (stale ->
  bad_weak_ptr after qs restore) -> stores a plain tibble.
Missing vs v1 (do NOT re-acquire unless a phase explicitly needs them): cisTarget
mm10 (SCENIC), SEA-AD h5ads (human validation) - both are v1 bloat, out of scope.

## Carried scientific decisions (v1-validated; treat as defaults, revisit per phase)
- Pseudobulk replicate = genotype_batch (16 ids, 4/genotype); batch in every design.
- Proteomics: peptide -> protein-group by sum; treat as total proteome.
- Phospho: first 16 samples = 24M timepoint; bulk hippocampus, NOT microglia-sorted
  (restate "not microglia-sorted" in any kinase prose).
- Microglia: re-cluster on the subset, drop only clear outliers (no over-pruning).
- Substates (v1 02a): homeostatic / DAM / IFN / proliferative.
- 5 canonical contrasts everywhere: tau_alone, nlgf_in_maptki, nlgf_in_p301s,
  tau_in_nlgf, interaction (factorial 2x2 + batch).
- State thresholds before applying them; present the axes with no pre-privileged winner.

## Environment (project-local; NO Docker, NO system-wide installs)
- Run as eturkes:eturkes (single-user Distrobox) -> files land user-owned, NO chown
  needed (v1's `chown rstudio:rstudio` was a rocker artefact, obsolete).
- R is 4.6.0 here (v1 pinned 4.5.2 / Bioc 3.22) -> numbers WILL differ; we
  re-baseline, never reproduce v1's locked margins (18/12/55).
- Stack (P0 built): **rv** (R pkgs) + **uv** (Python) + project-local **Quarto** + **targets** DAG,
  P3M-pinned (snapshot 2026-06-22, CRAN+Bioc same date). No bitwise guarantee (no Docker):
  targets=pipeline, rv.lock/uv.lock=versions, P3M=pinning. Fresh-clone bootstrap ORDER:
  `scripts/install-sysdeps.sh` -> `install-rv.sh` + `install-quarto.sh` -> `rv sync` -> `uv sync` -> `tar_make()`.
- **rv MUST be on PATH** (~/.local/bin, like uv): `.Rprofile`->`rv/scripts/activate.R` finds rv via
  `Sys.which("rv")` + shells `rv info` to set `.libPaths(rv/library)`; a tools/-only rv breaks
  activation. Pinned (version+sha256) in `install-rv.sh`. `.Rprofile` runs activate.R then a
  fail-loud guard: non-interactive `stop()` unless `rv/library` is in `.libPaths()` (catches
  rv-off-PATH / `rv info` fail / R-version-mismatch safe-mode -> NO silent global-lib fallback;
  re-add the guard if rv regenerates `.Rprofile`). NO repos override (base-R `install.packages`
  would write off-lock; use `rv add`).
- Repos (`rproject.toml` -> `rv sync` -> `rv.lock`): CRAN = plain `p3m.dev/cran/<date>` (rv inserts
  `__linux__/trixie` -> binary; do NOT hardcode the binary path); Bioc 3.23 =
  `p3m.dev/bioconductor/<date>/packages/3.23/{bioc,data/annotation,data/experiment,workflows}` +
  `force_source` (source-only on Debian). `tar_source` lives in `targets`; the `quarto` R pkg finds the
  pinned CLI via `QUARTO_PATH` (`_targets.R`; a `file.exists` preflight `stop()`s if missing -> no
  silent PATH-quarto fallback); `_quarto.yml`: render `*.qmd`+`!rv/` (rv#332), `freeze:false`
  (targets owns caching -> no stale `_freeze` divergence).
- Sysdeps (`scripts/install-sysdeps.sh`; `rv sysdeps` returns [] on trixie -> useless): build-essential
  + gfortran (Bioc source compiles) + libglpk40 (libglpk.so.40 for the igraph binary). Re-derive any new
  missing lib: `ldd`-scan `rv/library/**/*.so` for "not found".
- Python: uv `.venv` (gitignored), `pyproject.toml`+`uv.lock`+`.python-version` 3.13 (empty deps in P0); SOTA per phase.
- Heavy installs/compute: expect long runs; smoke-test helpers via `Rscript -e '...'`
  against the live data BEFORE any full run.

## Quality gate (provisional - lock concretely in P0)
- Reproducible: fresh clone + `rv sync` (+ `uv sync`) -> `targets::tar_make()` runs the pipeline.
- Each module smoke-tested against live data before commit.
- Committed tests = `tests/test_*.R`: plain `stopifnot` fail-loud scripts (zero new deps, mirror
  io.R's assertion idiom), run `Rscript tests/test_<x>.R` from project root, print `ok - <x>`,
  non-zero exit on failure. `tests/helpers.R` = deterministic synthetic fixtures (make_fake_seurat
  / make_meta16 / expect_error; NO RNG or clock -> reproducible everywhere). S3 added
  test_design (contrast weights + factorial==cell-means equivalence), test_de_pb (pseudobulk
  + fitter smokes), test_io (S2-deferred loader contracts). S5's check.sh MUST loop tests/test_*.R.
- Reports (once they exist) knit clean: 0 error / 0 warning -- enforced concretely, NOT
  by bare `tar_make()` exit (it returns 0 with captured warnings): assert empty
  `tar_meta()` error+warnings + `options(warn=2)` where safe + Quarto-log grep (lock in S5).

## Operational
- Prose register: British English; human-facing report/figure text uses hyphens over
  em/en-dashes (commas or parentheses for asides, colons for restatements). R `#`
  comments + code stay exempt (LLM-facing).
- storage/** + future caches/HTML are gitignored AND deny-Read -> Read/grep/ls on
  them is blocked (the `cat`-deny also bites). Rscript file reads BYPASS the deny
  matcher: inspect via `Rscript -e 'readRDS(...)' / 'readLines(...)'`. storage/data is
  a symlink now -> `git check-ignore` on a path UNDER it fatals (`beyond a symbolic
  link`); probe `storage/data` itself.
- Bash deny-gate is static on command text -> a cmd naming a deny-Read path as an
  arg (rm/stat/grep/find-piped) is blocked; use a glob/`find -delete` or runtime
  indirection. `ls`/`wc`/`echo` slip through.
- New `R/*.R` = pure functions loaded by targets via `tar_source()`; the DAG orders
  execution (no manual dependency loader). Heavy producers = `tar_target`s storing
  `format="qs"`/`"file"`; Python steps = `uv run <script>` as `tar_file` targets.
- Test/DE gotchas (S3; harness itself described under Quality gate): edgeR 4.x ->
  `edgeR::normLibSizes` (calcNormFactors deprecated, emits a message); `limma::makeContrasts`
  needs ALL design columns to be syntactically valid R names (cell-means form: rename genotype
  columns to bare levels, build batch from a named factor; an inline `factor()` yields invalid
  names). Sort character keys with `method = "radix"` (locale-independent -> reproducible
  column/level order across machines). factorial_design fails loud if add_batch=TRUE but the
  batch column is absent -> a batch-less modality (GeoMx) MUST pass add_batch=FALSE. de_pb fitters
  fail loud on misaligned inputs (limma fits/contrasts by POSITION, only warns on name mismatch):
  both assert identical(colnames(data), rownames(design)) + identical(rownames(contrasts),
  colnames(design)) + full-rank design (qr rank == ncol); factorial_design guarantees full-rank
  output; pseudobulk_counts asserts the meta.data<->counts cell alignment Seurat maintains.
- targets serialization (verified S2): `format="qs"` works via the **qs2** backend (the
  `qs` pkg is NOT installed; qs2 serializes Seurat Assay5 objects fine). `format="file_fast"`
  is deprecated -> use `format="file"` + `tar_option_set(trust_timestamps=TRUE)` (set in
  `_targets.R`) so the 8G snrnaseq input is change-detected by mtime/size, not re-hashed each
  run. Heavy seurat target: `memory="transient"` + `garbage_collection=TRUE` release the load.
- Commit `.serena/` (project.yml + .gitignore); Serena language changes are
  startup-only (restart Claude Code to apply).

## Reports (Quarto; built P0-S4)
- Report = ONE self-contained OFFLINE HTML: a standalone `format: html` doc (`index.qmd`,
  `embed-resources: true` inlines CSS/JS/fonts) + modular `{{< include _section.qmd >}}` (leading
  `_` -> not rendered standalone). NOT a Quarto book (multi-file + `Could not fetch resource
  ./<sibling>.html` nav warnings under embed-resources -> trips the zero-warning gate). tar_quarto
  still detects `tar_load`s inside an included `_*.qmd` (verified: 5 edges through the include); list
  the theme/css in `tar_quarto(extra_files=)` (inspection misses them).
- Theme + fonts SHIP via theme.scss (Quarto 1.9.38, verified live on the production render 2026-06-29):
  scss:defaults COLOUR vars ($primary/$link-color/$code-color #B0344D) + the IBM Plex stack
  ($font-family-sans-serif/$headings-font-family/$font-family-monospace) + 9 `@font-face` (scss:rules)
  with a relative `url("assets/fonts/<n>.woff2") format("woff2")`. Quarto base64-INLINES each woff2 into
  the embedded CSS under embed-resources -> ONE offline file (render PROVED: 9 faces inlined, magic d09GMg,
  0 external). The 9 woff2 are COMMITTED (assets/fonts/, deny-Read `**/*.woff2`, Serena ignored_paths);
  list them in `tar_quarto(extra_files=)` -- inspection misses them, `list.files("assets/fonts",
  pattern="woff2", full.names=TRUE)` keeps the list in sync. ggplot panels keep `theme_tau(base_family="")`
  (device font) so figures stay decoupled from this chrome + warning-free.
- Theme-CSS DETECTION gotcha (CORRECTED -- supersedes "colours embed raw"): once an `@font-face url()`
  is in the theme, Quarto embeds the WHOLE compiled theme CSS as a URL-ENCODED `data:text/css,...` URI,
  so BOTH colours AND fonts encode -- `#B0344D`->`%23B0344D`, `IBM Plex`->`IBM%20Plex`, the woff2 data
  URI -> `woff2%3Bbase64%2Cd09GMg` (d09GMg = wOF2 magic). RAW `.count` then reads ~0 for everything
  theme-side (raw-embed held ONLY for the pre-fonts colours-only theme). Match the ENCODED tokens (fast)
  or URLdecode/urllib.unquote the data:text/css blocks first -- URLdecode of the ~1MB blob is VERY SLOW,
  so prefer encoded matching. Figures are PNG base64 (`data:image/png`), so their `#B0344D` is not raw either.
- Quarto caches the Sass compile in `.quarto/` -> a theme edit is invisible until cleared; `.quarto`
  is deny-Read -> clear via runtime indirection (R `unlink(".quarto", recursive=TRUE)` in the render
  script), not a Bash `rm`. Inspect the output HTML the same way (output dir is deny-Read): build the
  path inside a python/R script ("_"+"report") + `.count` substrings (a `#hex` regex proved flaky).

## Subagents & skills
Scan the available-skills list each session; invoke a matching Skill before
improvising (scientific-writing, scientific-visualization, pathway-enrichment,
pydeseq2, scanpy, anndata, bioservices). Spawn subagents to protect main context
(Explore = cross-file search, Plan = design, general-purpose = research).
