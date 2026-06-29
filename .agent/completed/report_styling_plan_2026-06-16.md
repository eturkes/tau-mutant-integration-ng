# Report styling + prose pass (UI/UX policy alignment)

Intent: bring the 3 human-facing HTML reports into line with the CLAUDE.md UI/UX
bullet — "unique fonts, cohesive colors and themes ... human-facing [prose] ...
human-like and clear of LLM-isms, smells, and cliches" + the concrete prose rules
the user added 2026-06-16: "prefer to use hyphens over other kinds of dashes,
enumerate flexibly, and vary comparative constructions". Audit (2026-06-16) found all three
reports render with stock rmarkdown `html_document` defaults (Bootstrap 3
"default" theme, system sans-serif, pygments highlight); zero intentional document
theming; the only styling present is per-plot ggplot `theme_bw()`. Scope = the 3
top-level reports only (`analysis.Rmd`, `summary.Rmd`, `synthesis.Rmd`); child
`rmd/*.Rmd` inherit the parent theme and need no per-file change.

## STATUS
- [DONE] S0  aesthetic direction — RESOLVED 2026-06-16: custom bslib bs_theme(v5)
- [DONE] S1  shared theme helper — R/report.R::report_theme(), sourced last in helpers.R
- [DONE] S2  wired theme into all 3 report YAML headers (self-contained !expr)
- [DONE] S3  human-facing prose dash/LLM-ism sweep (em/en → hyphen/comma/colon) — 2026-06-16
- [DONE] S4  cold-knit all 3 (incl publication_mode); 0 error/0 warning each   — 2026-06-16

## S0 — aesthetic direction (DECISION GATE at start)
Default (recommended): a single `bslib::bs_theme(version = 5, ...)` object —
Bootstrap 5 base, a distinctive Google-font pairing (humanist serif headings +
clean sans body + a mono for code), and a cohesive palette anchored to the
project's own identity (derive an accent from the `genotype_colours` constant so
report chrome and figures share one visual language). Centralise it so the three
reports stay identical and drift-proof.
Alternatives: (B) a stock bootswatch preset (`theme: litera/flatly/cosmo`) — one
line, clean, but shared with countless sites and weak on the "unique" requirement,
limited font control; (C) hand-written CSS via `includes: in_header` — maximal
control, maximal maintenance + drift risk. Resolve via AskUserQuestion before S1.

## S1 — shared theme helper
Add `report_theme()` (new `R/report.R`, sourced last in `R/helpers.R`) returning
the agreed bslib object; or, if the gate picks B/C, the matching YAML snippet / CSS
file. Load fonts via `bslib::font_google(..., local = TRUE)` so the html stays
self-contained and renders offline. Smoke-test: `Rscript -e 'bslib::bs_theme(...)'`
builds without error before any knit.

## S2 — wire into reports
Replace each report's `output: html_document:` theme slot with `theme: !expr
report_theme()` (bslib) — note html_document under a bslib theme auto-upgrades to
Bootstrap 5; verify `toc_float`, `code_folding: hide`, and `df_print: paged` still
behave. Keep every existing key (number_sections, fig sizes, params). No change to
child rmds; `publication_mode` and the locked section-17 content stay untouched.

## S3 — prose sweep
Scan only HUMAN-FACING narrative prose (report body text, section intros, the
synthesis conclusion, plus rendered strings: kable captions, verdict text,
`cat()`/`sprintf()` output) for LLM-isms / smells / cliches ("delve", "it is worth
noting", "in conclusion", hedge-stacking, empty transitions, "rich tapestry", etc.)
and tighten to the dry, precise register CLAUDE.md models. Apply the concrete
2026-06-16 rules: (a) prefer hyphens over other kinds of dashes — recast em-dash
asides as commas or parentheses (the user's own register), convert en-dashes
(ranges, compounds, citation page-spans) to hyphens; leave maths notation (×, →,
⊣, superscripts, numeric minus signs) untouched; (b) enumerate flexibly; (c) vary
comparative constructions. Leave R `#` comments and code (LLM-facing) and the
locked scientific numbers / claims / verdicts untouched, byte-for-byte. British
English throughout. Surfaces by dash density: `synthesis.Rmd` (flagship, ~98
sites), `rmd/16,11,17,02c,12,15` + `analysis.Rmd` (human-facing strings only);
`summary.Rmd` / `summary_chunks/` carry none.

## S4 — verify + close
Knit all three cold (`rmarkdown::render` each); `grep -c 'class="error"'` and
`'class="warning"'` both 0 per report; eyeball the rendered theme (fonts load,
palette cohesive, TOC float intact, code folding works, figures legible on the new
background). chown rstudio:rstudio the html + any new R file. One scoped commit.

## Execution model
S0 gate first (user). S1->S2 are a tight pair (theme plumbing) for one session; S3
(prose) is independent and can be its own session; S4 closes. Per step: mark DONE
in STATUS + a completion note, chown, knit-verify, commit.

## Anti-anchoring guardrails
- The section-17 claims ledger, contest margins (18/12/55), and all scientific
  numbers are LOCKED — styling / prose passes must not alter a single figure or
  verdict.
- `publication_mode` masking must still work after the theme swap; test a
  `params = list(publication_mode = TRUE)` knit of summary at least once.
- Theme is chrome only: do NOT restyle the ggplot figures' data ink to match; at
  most align figure accent colours with the palette if trivial.
- "Unique" is not "ornate". Favour legibility and a restrained scientific
  aesthetic over decoration; the audience is the user and collaborators.
- Keep reports self-contained (fonts cached locally, not live CDN pulls) so the
  gitignored html still renders offline.

## Completion notes
S1+S2 done 2026-06-16 (one session, theme plumbing pair). Built
`R/report.R::report_theme()` — a single `bslib::bs_theme(version = 5)` returning
the IBM Plex superfamily (Serif headings / Sans body / Mono code, each
`font_google(local = TRUE)` so fonts embed and the html renders offline) with
`primary`/`link-color` = NLGF_P301S crimson `#B0344D` (the divergence-endpoint
genotype, the palette anchor) and a steel `code-color` `#3F5A6B` so inline code
reads distinct from links. Sourced last in `R/helpers.R` (after plot.R).

Wiring (S2): each of the 3 reports gained one YAML line under `html_document:`,
`theme: !expr 'local({ source("R/report.R", local = TRUE); report_theme() })'`.
The self-contained `local({source;…})` form is deliberate — a throwaway probe
confirmed rmarkdown evaluates output-format `!expr` at front-matter parse time in
a fresh `render` session that has sourced nothing, so a bare `report_theme()`
(relying on a pre-sourced global) would NOT resolve under the documented re-knit
command. Sourcing report.R inline inside the expr makes it robust while keeping
ONE shared object (drift-proof). Every pre-existing YAML key was preserved.

Deviation from plan default: the accent is the inlined literal `#B0344D` (not a
live read of `genotype_colours[["NLGF_P301S"]]`) because report.R must stand alone
at parse time, before constants.R loads; cross-referenced in a comment, palette is
locked so stable. KISS over conditional-sourcing.

Verification: all 3 reports knit clean (analysis 9 min; summary default AND
`publication_mode = TRUE`; synthesis) — 0 errors, 0 warnings each; BS5 upgrade,
all 3 IBM Plex families embedded, crimson + code-folding + toc_float + paged
tables all intact in analysis.html. publication_mode masking guardrail passes
under the new theme. chown'd R/report.R + the html outputs (and corrected lone
root-owned R/helpers.R drift).

Remaining: S3 (user-facing prose LLM-ism/cliche sweep) is independent — own
session. S4 (final cold knit + close) runs after S3.

S3+S4 done 2026-06-16 (one session, prose sweep + close), triggered by the
CLAUDE.md UI/UX bullet gaining concrete rules: "prefer to use hyphens over other
kinds of dashes, enumerate flexibly, and vary comparative constructions" + the
"human-facing" terminology and the added "smells". The high-signal, mechanically
verifiable part is the dash rule, which is what S3 executed.

synthesis.Rmd (flagship): 59 em-dashes recast (commas for short asides;
parentheses where the aside itself carried commas; colons for a restatement /
list) and 39 en-dashes (numeric ranges, compounds, citation page-spans) → hyphens
— matching the user's own register in CLAUDE.md (commas/parentheses for asides,
hyphens for compounds, no em/en-dashes). analysis.Rmd (1 site) + 11 child
rmd/*.Rmd human-facing prose AND rendered strings (kable captions, verdict text,
sprintf("%s - %s") separators, a ggplot axis label) swept; R `#` comments left
untouched (LLM-facing, exempt) so only rmd/11's 5 comment em-dashes remain
repo-wide. Both pre-built figure scripts (build_capstone_figure.py,
build_synthesis_model_figure.py) had their baked-in PNG text swept (titles,
panel headers, band captions → colons; the Apoe-Trem2/App-Cd74 node label →
hyphens) and BOTH PNGs regenerated + visually confirmed clean.

Content-safety proof: every edited file's alphanumeric stream is byte-identical
to HEAD (`diff <(tr -cd '[:alnum:]')`), so not one number, claim, citation or
gene name changed — only punctuation. The child-rmd sweep was delegated to an
Opus subagent under that same alnum-invariant guard and re-verified here per file.
British English and the locked §17 numbers untouched throughout.

LLM-ism/cliche + comparative component: the report prose was already in the dry
register CLAUDE.md models; a scan surfaced only 2 low-severity borderline items in
rmd/16 (a "closes that gap" setup phrase; one quantified hedge-stack), both judged
load-bearing and left. Comparative-construction variation applied minimally — the
existing "rather than" usages are precise scientific contrasts, not filler, so
forcing synonyms would risk meaning drift for no gain.

S4: cold-knit analysis (~9 min) + synthesis + summary (default AND
publication_mode) — 0 error / 0 warning each (confirmed via Rscript readLines, the
deny-Read-immune count method; the 3 HTMLs are gitignored so none are committed).
publication_mode masking still behaves under the theme. chowned the edited sources
and regenerated PNGs back to rstudio:rstudio. Plan S0-S4 all DONE → archived to
completed/; durable summary folded into DIGEST.md.
