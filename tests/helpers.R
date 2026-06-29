# Shared test fixtures + assertions. Source AFTER R/constants.R (uses genotype_levels).
# Fixtures are fully deterministic (no RNG, no clock) so tests are reproducible everywhere.
# Convention: each tests/test_*.R sources the R/ files it exercises + this helper, runs
# stopifnot-based checks (fail loud, non-zero exit on first failure), prints "ok - <name>".

# Assert that evaluating `expr` raises an error (for fail-loud contract checks). With `pattern`
# (a literal substring), also require the error message to contain it -> the error fired for the
# INTENDED reason, not an incidental earlier failure.
expect_error <- function(expr, pattern = NULL) {
  e <- try(expr, silent = TRUE)
  stopifnot(inherits(e, "try-error"))
  if (!is.null(pattern)) {
    stopifnot(grepl(pattern, conditionMessage(attr(e, "condition")), fixed = TRUE))
  }
}

# Canonical 16-row pseudobulk-style metadata: 4 genotypes x 4 batches, fully crossed.
make_meta16 <- function() {
  data.frame(
    genotype  = rep(genotype_levels, each = 4L),
    batch     = rep(paste0("batch0", 1:4), times = 4L),
    row.names = paste0("s", 1:16),
    stringsAsFactors = FALSE
  )
}

# Deterministic counts matrix (genes x cells), all entries >= 1 (no all-zero features/cells).
# Returned sparse (dgCMatrix) so CreateSeuratObject ingests it without a coercion warning.
.fake_counts <- function(n_genes, cell_names) {
  ens <- paste0("ENSMUSG", sprintf("%08d", seq_len(n_genes)))
  m <- 1 + outer(seq_len(n_genes), seq_along(cell_names),
                 function(i, j) (i * 3L + j * 7L) %% 11L)
  m <- matrix(as.double(m), nrow = n_genes, dimnames = list(ens, cell_names))
  Matrix::Matrix(m, sparse = TRUE)
}

# Minimal synthetic Seurat object for io / pseudobulk contract tests.
#   with_broad : add broad_annotations == "Microglia" on the design cells.
#   n_other    : append this many non-microglia ("Neuron") cells (subset out by load_snrnaseq).
#   with_sct   : add an SCT assay (load_snrnaseq drops it).
# Microglia cells form a fully-crossed genotypes x batches grid, cells_per cells each.
make_fake_seurat <- function(genotypes = genotype_levels,
                             batches = paste0("batch0", 1:4),
                             cells_per = 4L, n_genes = 30L,
                             with_broad = FALSE, with_sct = FALSE, n_other = 0L) {
  combos <- expand.grid(genotype = genotypes, batch = batches, stringsAsFactors = FALSE)
  meta <- combos[rep(seq_len(nrow(combos)), each = cells_per), , drop = FALSE]
  meta$genotype_batch <- paste(meta$genotype, meta$batch, sep = "_")
  if (with_broad) meta$broad_annotations <- "Microglia"
  if (n_other > 0L && with_broad) {
    other <- data.frame(
      genotype = rep(genotypes, length.out = n_other),
      batch    = rep(batches, length.out = n_other),
      genotype_batch = NA_character_,
      broad_annotations = "Neuron",
      stringsAsFactors = FALSE
    )
    meta <- rbind(meta, other)
  }
  rownames(meta) <- paste0("cell", seq_len(nrow(meta)))
  counts <- .fake_counts(n_genes, rownames(meta))
  obj <- SeuratObject::CreateSeuratObject(counts = counts, meta.data = meta)
  obj@misc$geneids <- paste0("Sym", seq_len(n_genes))
  if (with_sct) obj[["SCT"]] <- SeuratObject::CreateAssay5Object(counts = counts)
  obj
}
