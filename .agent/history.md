# Decision digests (new project)

Compact per-phase decision summaries for the fresh `main` rebuild - read these, not
the full archived plans. v1 digests live in `archive_digest.md` (full v1 plans on
branch `archive`).

## P0 Foundations -- closed 2026-06-29 (-> `.agent/completed/p0_foundations_plan_2026-06-29.md`)

Built the reproducible spine for P1-P5: project-local env, shared pure-fn modules, 4 analysis-ready
modalities, a QC-sanity report, 2x2 factorial + 5-contrast machinery, a concrete quality gate. NO analysis
(P1+). Live operational facts -> memory.md; wiring -> map.md. Below = decisions + rejected alternatives.

- STACK: targets DAG + rv (R) + uv (Python) + project-local Quarto, P3M-pinned (snapshot 2026-06-22). rv
  over renv (agent-native, declarative, mirrors uv; accepts pre-1.0 + hand-wired Bioc 3.23). targets over
  plain Rmd (content-hash skip pays on the 8G object + multi-hour fits). NO Docker -> repro =
  targets(pipeline) + locks(versions) + P3M(pinning), explicitly NO bitwise guarantee. R 4.6 / Bioc 3.23 ->
  RE-BASELINE (numbers differ from v1's R 4.5.2; never reproduce v1's locked margins).
- REPORT: pivoted Quarto BOOK -> ONE standalone offline HTML (`index.qmd` format:html embed-resources +
  `{{< include _qc.qmd >}}`). Book rejected: multi-file nav emits `Could not fetch resource ./<sibling>.html`
  under embed-resources -> trips the zero-warning gate. theme.scss = crimson #B0344D + IBM Plex (9 woff2
  base64-inlined offline; woff2 COMMITTED, deny-Read).
- MODULES: pure fns via `tar_source` (the DAG orders execution -> no manual loader; supersedes v1 helpers.R):
  constants/utils/io/design/de_pb/plot/spine. 6 modalities materialized as qs2 targets. 5 canonical contrasts
  via TWO equivalent parameterisations (factorial `~tau+nlgf+tau_nlgf[+batch]` AND cell-means `~0+genotype`),
  proven equal by a property test. de_pb = limma-voom/log only; edgeR-QLF/DESeq2/dream deferred (KISS). S3 =
  machinery only (no DE results).
- QUALITY GATE: `scripts/check.sh` -- concrete, fail-loud, zero-fault (4 phases; detail -> memory.md).
  Enforces MORE than tar_make's exit (it returns 0 on CAPTURED warnings): tar_meta all-NA (scoped to the
  manifest) + render-log warning grep. Tests = stopifnot harness + synthetic fixtures (no testthat dep).

Verification (honest): each module smoke-tested vs LIVE data at its step (S2 shapes incl. the 8G load, S3
16-col pseudobulk, S4 clean render); `scripts/check.sh` green end-to-end on the populated store (S5). The
one-shot FRESH-CLONE rebuild is the design contract, verified incrementally across S1-S5, NOT re-run at S5
(heavy 8G reload). Gate boundary: the render-log grep only scans renders that occur DURING a gate run.

Deferred -> P1+: microglia reprocess (SCT/Harmony/cluster-prune), substate assignment, single-cell DE
(NEBULA/glmmTMB), pseudobulk DE RESULTS, edgeR-QLF/DESeq2/dream, BPCells (8G RAM relief), col-map DAG
validation (P4). Out of scope (v1 bloat): cisTarget mm10/SCENIC, SEA-AD human validation.
