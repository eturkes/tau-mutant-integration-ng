# Reproducibility-spine self-check: a pure function recording the pinned stack, loaded by
# targets via tar_source(). Proves the recursive R/ loader wiring; P1+ adds the analysis modules
# (constants, utils, io, design, de_pb, plot ...) alongside it.

spine_versions <- function() {
  pkgs <- c("Seurat", "targets", "limma", "edgeR", "qs2")
  data.frame(
    component = c("R", pkgs),
    version = c(
      as.character(getRversion()),
      vapply(pkgs, function(p) as.character(utils::packageVersion(p)), character(1))
    ),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}
