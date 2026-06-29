#!/usr/bin/env bash
# Render analysis.Rmd to HTML. Designed to be launched once and run to completion
# in the background. All output goes to knit.log next to the Rmd.
set -euo pipefail
cd /home/rstudio/tau-mutant-integration-ng
exec Rscript -e 'rmarkdown::render("analysis.Rmd", output_file = "analysis.html")'
