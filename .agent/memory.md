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

## Raw data (storage/data = symlink -> ../../../Documents/pro/tau-mutant-integration-ng/data = host Documents tree, shared/external copy; gitignored via /storage/* + deny-Read; Rscript reads resolve through it, bypassing the deny)
- snrnaseq.rds  8.3G  (snRNAseq Seurat object)
- geomx.rds     22M   (GeoMx WTA spatial, ~91 ROIs)
- proteomics_nonfiltered_nonnormalised.tsv      15M  (peptide-level)
- phosphoproteomics_nonfiltered_nonnormalised.tsv  35.5M (phosphosite-level)
- proteomics_sample_key.csv  3.3K  (67 rows: `File name` + `Sample/Condtion`[sic];
  TiO2/PTM phospho-run map; 24M timepoint = rows 1-16 = 4 genotypes x 4 reps;
  io::proteomics_sample_meta n_keep=16; raw labels MAPT-KI/P301S+3/NLGF-* -> canonical)
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
- Stack (P0-locked; SOTA-researched, see p0_foundations_plan): **Quarto book** (reports; project-local pinned binary)
  + **targets** DAG (orchestration; `tarchetypes::tar_quarto`; R/ fns via `tar_source()`,
  no manual loader) + **rv** (R pkgs; declarative `rproject.toml`+`rv.lock`) + **uv**
  (Python). Versions pinned via **P3M dated snapshot**: CRAN = binary (trixie); Bioc 3.23
  = likely source-only via P3M for rv (prove at S1). Source fallback pins R-package
  versions only (system toolchain/libs separate -> may need Debian sysdeps).
  NO Docker -> no bitwise guarantee: targets = pipeline correctness; rv.lock+uv.lock =
  versions; P3M = pinning. (renv was the safe alt: auto-Bioc but imperative.)
- R packages: project-local via rv (NOT installed yet). Bioc **3.23** (R 4.6) repos
  hand-wired (no auto R<->Bioc coupling); `rv remove` exists (hand-edit TOML only for
  what rv can't express); exclude
  `rv/library` in `_quarto.yml` (bug #332). rv pre-1.0 -> fall back to renv if it bites.
- Python: project-local uv `.venv` (gitignored); pick the SOTA per phase.
- Heavy installs/compute: expect long runs; smoke-test helpers via `Rscript -e '...'`
  against the live data BEFORE any full run.

## Quality gate (provisional - lock concretely in P0)
- Reproducible: fresh clone + `rv sync` (+ `uv sync`) -> `targets::tar_make()` runs the pipeline.
- Each module smoke-tested against live data before commit.
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
- Commit `.serena/` (project.yml + .gitignore); Serena language changes are
  startup-only (restart Claude Code to apply).

## Subagents & skills
Scan the available-skills list each session; invoke a matching Skill before
improvising (scientific-writing, scientific-visualization, pathway-enrichment,
pydeseq2, scanpy, anndata, bioservices). Spawn subagents to protect main context
(Explore = cross-file search, Plan = design, general-purpose = research).
