# Shared HTML report theme. ONE bslib object reused by all three top-level
# reports (analysis.Rmd, summary.Rmd, synthesis.Rmd) so their chrome stays
# identical and drift-proof. Each report's YAML wires it self-contained via
#   theme: !expr 'local({ source("R/report.R", local = TRUE); report_theme() })'
# which resolves at front-matter parse time even in a fresh `rmarkdown::render`
# session (the documented re-knit form sources nothing first). report_theme is
# always called by name, so this file's statement order is unconstrained.
#
# Design (CLAUDE.md UI/UX policy: unique fonts + cohesive theme + project palette):
#   * Bootstrap 5 base (html_document auto-upgrades under a bslib theme).
#   * IBM Plex superfamily — Serif headings / Sans body / Mono code — one coherent
#     type family: distinctive without being ornate. Fonts embed locally
#     (local = TRUE) so the self-contained html renders offline.
#   * Accent anchored to the project's own identity: NLGF_P301S crimson, the
#     genotype the divergence contrast culminates in. Drives links + active TOC.
#     Inline code takes a neutral steel (a darkened MAPTKI blue from the same
#     ramp) so code reads distinct from links.
#
# The accent is the literal value of genotype_colours[["NLGF_P301S"]] (R/constants.R).
# It is inlined rather than read from that constant because this runs at YAML-parse
# time, before helpers.R/constants.R are sourced; the palette constants are locked,
# so the cross-reference is stable. Keep the two in sync if the palette ever moves.

report_theme <- function() {
  accent <- "#B0344D"  # == genotype_colours[["NLGF_P301S"]] — crimson palette anchor
  bslib::bs_theme(
    version      = 5,
    heading_font = bslib::font_google("IBM Plex Serif", wght = c(500, 600, 700),
                                      ital = 0, local = TRUE),
    base_font    = bslib::font_google("IBM Plex Sans",  wght = c(400, 500, 600, 700),
                                      ital = c(0, 1), local = TRUE),
    code_font    = bslib::font_google("IBM Plex Mono",  wght = c(400, 600),
                                      ital = 0, local = TRUE),
    primary      = accent,
    "link-color" = accent,
    "code-color" = "#3F5A6B"   # steel (darkened MAPTKI blue) — code != link
  )
}
