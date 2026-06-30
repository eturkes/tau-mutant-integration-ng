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

# Deterministic reduced-dim embedding + cluster labels for slingshot trajectory tests (no
# RNG). Each cluster is a well-conditioned, full-rank blob: D1 = an along-axis gradient, the
# other dims = distinct-frequency sinusoids (mutually independent, comparable amplitude -> no
# singular per-cluster covariance, which slingshot's scaled MST distance would choke on).
# Homeostatic (D1 ~ 0) and DAM (D1 ~ 2) are well separated along D1 -> a clean Homeostatic->
# DAM lineage with monotone pseudotime. with_ifn appends an "IFN" blob displaced HIGH on D2
# near the homeostatic end -> the MST branches, so the DAM-terminal lineage must be selected
# and IFN cells get NA pt.
make_trajectory_embedding <- function(n_per = 80L, n_dims = 4L, with_ifn = FALSE) {
  stopifnot(n_dims >= 2L)
  blob <- function(n, d1_base, d2_base = 0) {
    i <- seq_len(n)
    do.call(cbind, lapply(seq_len(n_dims), function(j) {
      if (j == 1L) d1_base + (i / n) * 0.6            # gradient along the activation axis
      else if (j == 2L) d2_base + 0.4 * sin(i * 1.1)  # branch dimension
      else 0.4 * sin(i * (0.7 * j + 0.3))             # distinct freqs -> full-rank covariance
    }))
  }
  m <- rbind(blob(n_per, 0), blob(n_per, 2.0))
  lab <- c(rep("Homeostatic", n_per), rep("DAM", n_per))
  if (with_ifn) {
    ki <- n_per %/% 2L
    m <- rbind(m, blob(ki, 0.3, d2_base = 2.0))       # high on D2 near the homeostatic end
    lab <- c(lab, rep("IFN", ki))
  }
  rownames(m) <- paste0("c", seq_len(nrow(m))); colnames(m) <- paste0("D", seq_len(n_dims))
  list(embedding = m, labels = lab)
}

# Deterministic per-cell trajectory frame (NO RNG) for the S2a estimation-core tests -- mirrors
# the real microglia_trajectory$cell_frame columns (no batch col -> exercises derive_batch). 4x4
# genotype x batch grid; each unit = per_state Homeostatic + (per_state + dam_extra iff the
# double-mutant) DAM cells. pt_raw = 1 + adv[g] + ramp (Homeostatic) / 4 + adv[g] + ramp (DAM),
# ramp = ((seq_len(n) - 0.5)/n)*0.3 -> the MIDPOINT rule makes EACH block mean EXACTLY 0.15 for
# ANY block size, so an unequal DAM count does NOT shift a within-state mean (the pure-composition
# fixture is EXACTLY pure). DEFAULT adv encodes a pure WITHIN-STATE interaction
# (1.6-0.2)-(1.0-0) = 0.4 at constant composition -> 100% progression loading; FLAT adv + dam_extra>0
# -> a pure COMPOSITION interaction (extra DAM mass only in NLGF_P301S) -> 100% composition loading,
# exact. jitter>0 adds a deterministic NON-additive ((gi*bi) %% 5)*jitter per unit (gi/bi = genotype
# /batch indices) -> breaks the saturated design's zero residual (sigma>0) for S2b's structural test.
make_trajectory_cell_frame <- function(per_state = 6L,
                                       adv = c(MAPTKI = 0, P301S = 0.2,
                                               NLGF_MAPTKI = 1.0, NLGF_P301S = 1.6),
                                       dam_extra = 0L, jitter = 0) {
  stopifnot(setequal(names(adv), genotype_levels))
  batches <- sprintf("batch%02d", 1:4)
  ramp <- function(n) ((seq_len(n) - 0.5) / n) * 0.3       # block mean EXACTLY 0.15 for any n
  rows <- list(); k <- 0L
  for (gi in seq_along(genotype_levels)) {
    g <- genotype_levels[gi]
    for (bi in seq_along(batches)) {
      n_dam <- per_state + dam_extra * (g == "NLGF_P301S")
      jit   <- ((gi * bi) %% 5L) * jitter                  # NON-additive in (gi,bi) -> sigma>0 when jitter>0
      sub   <- c(rep("Homeostatic", per_state), rep("DAM", n_dam))
      pt    <- c(1 + adv[[g]] + ramp(per_state), 4 + adv[[g]] + ramp(n_dam)) + jit
      rows[[length(rows) + 1L]] <- data.frame(
        cell = paste0("c", k + seq_along(sub)),
        genotype_batch = paste(g, batches[bi], sep = "_"),
        genotype = g, substate = sub, on_lineage = TRUE, pt_raw = pt,
        row.names = NULL, stringsAsFactors = FALSE)
      k <- k + length(sub)
    }
  }
  cf <- do.call(rbind, rows)
  cf$genotype <- factor(cf$genotype, levels = genotype_levels)
  cf$substate <- factor(cf$substate, levels = sort(unique(cf$substate)))
  cf$pt01 <- squeeze_unit_interval(cf$pt_raw)
  cf[c("cell", "genotype_batch", "genotype", "substate", "on_lineage", "pt_raw", "pt01")]
}
