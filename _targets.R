# targets pipeline entrypoint. Pure functions in R/ load via tar_source() (the DAG orders
# execution -- no manual loader). Heavy producers (P1+) store format="qs" via qs2.
library(targets)
library(tarchetypes)

# Point the `quarto` R package (used by tar_quarto / quarto_inspect) at the pinned,
# project-local Quarto CLI; resolved from the project root on every tar_manifest/tar_make.
Sys.setenv(QUARTO_PATH = normalizePath("tools/quarto/bin/quarto", mustWork = FALSE))

tar_option_set(packages = "quarto")

tar_source("R")

list(
  # reproducibility-spine self-check: pinned-stack provenance via a tar_source()'d function
  tar_target(spine, spine_versions()),
  # Quarto book render (path = project root holding _quarto.yml); chapters grow over P1-P5
  tar_quarto(book, path = ".")
)
