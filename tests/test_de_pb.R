# S3 acceptance: pseudobulk + bulk DE machinery. Pseudobulk of a synthetic 4x4 Seurat -> 16
# genotype_batch columns (mirrors the live-object check run separately on microglia_seurat_raw);
# build_pseudobulk alignment + covariate-constancy guard; median_normalise / prevalence_filter
# logic; and end-to-end smokes that fit_limma_voom / fit_limma_log run and return 5 named
# per-contrast topTables (catches namespace / signature regressions; runs no real analysis).

source("R/constants.R")
source("R/design.R")
source("R/de_pb.R")
source("tests/helpers.R")

canonical <- c("tau_alone", "nlgf_in_maptki", "nlgf_in_p301s", "tau_in_nlgf", "interaction")
meta16 <- make_meta16()

# --- pseudobulk_counts: 16 columns + sum correctness ------------------------------------
obj <- make_fake_seurat(cells_per = 4L, n_genes = 30L)      # 64 cells, 16 genotype_batch levels
pc  <- pseudobulk_counts(obj, "genotype_batch")
stopifnot(ncol(pc) == 16L, nrow(pc) == 30L,
          identical(sort(colnames(pc), method = "radix"), colnames(pc)))   # deterministic (radix) column order
counts_m <- SeuratObject::GetAssayData(obj, assay = "RNA", layer = "counts")
for (gb in colnames(pc)) {                                   # verify the sum for EVERY group, not only the first
  cells <- rownames(obj@meta.data)[obj@meta.data$genotype_batch == gb]
  stopifnot(isTRUE(all.equal(unname(pc[, gb]), unname(Matrix::rowSums(counts_m[, cells, drop = FALSE])))))
}

# --- build_pseudobulk: alignment + covariate-constancy fail-loud ------------------------
pb <- build_pseudobulk(obj, "genotype_batch", c("genotype", "batch"))
stopifnot(ncol(pb$counts) == 16L, nrow(pb$meta) == 16L,
          identical(colnames(pb$counts), rownames(pb$meta)),
          all(c("genotype", "batch") %in% colnames(pb$meta)))
obj_bad <- obj
i_bad   <- which(obj_bad@meta.data$genotype_batch == colnames(pc)[1])[1]
obj_bad@meta.data$genotype[i_bad] <- "P301S"                 # genotype now varies within one pseudobulk sample
expect_error(build_pseudobulk(obj_bad, "genotype_batch", c("genotype", "batch")), "n_unique")

# --- median_normalise: exact shift ------------------------------------------------------
m  <- matrix(c(1, 2, 3, 11, 12, 13), nrow = 3)               # column medians 2 and 12, global 7
mn <- median_normalise(m)
stopifnot(isTRUE(all.equal(mn, matrix(c(6, 7, 8, 6, 7, 8), nrow = 3))))

# --- prevalence_filter: presence logic + rowname guard ----------------------------------
mat <- rbind(g1 = c(1, 2, 3, 4, 5, 6, 7, 8),                 # present in all -> both groups >=3 -> keep
             g2 = c(1, NA, NA, NA, NA, NA, NA, NA))          # A present 1, B present 0 -> drop
grp <- rep(c("A", "B"), each = 4L)
f   <- prevalence_filter(mat, grp, min_present = 3L, min_groups = 2L)
stopifnot(identical(rownames(f), "g1"), ncol(f) == 8L)
mat2 <- rbind(keep = c(1, 2, 3, 4, 5, 6, 7, 8), c(1, 2, 3, 4, 5, 6, 7, 8))
rownames(mat2) <- c("keep", "")                              # empty-name row dropped despite full prevalence
stopifnot(identical(rownames(prevalence_filter(mat2, grp)), "keep"))
expect_error(prevalence_filter(mat, c(NA, "A", "A", "A", "B", "B", "B", "B")), "anyNA")   # NA group -> fail loud

# --- fit_limma_voom smoke: 5 named topTables --------------------------------------------
fd     <- factorial_design(meta16)
counts <- 30 + outer(seq_len(200), seq_len(16), function(i, j) (i * 5L + j * 13L) %% 19L)
rownames(counts) <- paste0("ENS", seq_len(200)); colnames(counts) <- rownames(meta16)
vfit <- fit_limma_voom(counts, fd$design, fd$contrasts, min_count = 5)   # quality_weights = TRUE (default)
stopifnot(identical(names(vfit$top), canonical), vfit$kept >= 1L)
for (tt in vfit$top) {
  stopifnot(is.data.frame(tt), nrow(tt) == vfit$kept,                    # confint=TRUE -> CI.L/CI.R for effect-size reporting
            all(c("gene", "logFC", "CI.L", "CI.R", "P.Value", "adj.P.Val") %in% names(tt)))
}
# fail loud when counts columns are not aligned to the design rows (limma fits by position, not name)
expect_error(fit_limma_voom(counts[, 16:1], fd$design, fd$contrasts, min_count = 5), "colnames(counts)")

# --- fit_limma_log smoke (proteomics-style, no batch): 5 named topTables ----------------
meta_p <- data.frame(genotype = rep(genotype_levels, each = 4L), row.names = paste0("p", 1:16))
fdp    <- factorial_design(meta_p, add_batch = FALSE)
logmat <- 10 + outer(seq_len(100), seq_len(16), function(i, j) ((i * 2L + j * 3L) %% 7L)) / 3
rownames(logmat) <- paste0("PROT", seq_len(100)); colnames(logmat) <- rownames(meta_p)
lfit <- fit_limma_log(median_normalise(logmat), fdp$design, fdp$contrasts)
stopifnot(identical(names(lfit$top), canonical))
for (tt in lfit$top) {
  stopifnot(is.data.frame(tt), nrow(tt) == 100L,
            all(c("feature", "logFC", "P.Value", "adj.P.Val") %in% names(tt)))
}

# ============================== P1-S4: DE orchestration ==================================
# Synthetic microglia object: 16 fully-crossed units x 12 cells; a deterministic substate split
# (6 Homeostatic / 5 DAM / 1 IFN per unit) -> exercises the per-substate min-cell fit-or-skip.
s4obj <- make_fake_seurat(cells_per = 12L, n_genes = 200L)
.uidx <- ave(seq_len(ncol(s4obj)), s4obj$genotype_batch, FUN = seq_along)
s4obj$microglia_substate <- factor(
  ifelse(.uidx <= 6L, "Homeostatic", ifelse(.uidx <= 11L, "DAM", "IFN")),
  levels = microglia_substate_levels)

# --- cells= subsetting: pseudobulk only the DAM cells, still 16 units --------------------
dam_cells <- colnames(s4obj)[s4obj$microglia_substate == "DAM"]
pcsub <- pseudobulk_counts(s4obj, "genotype_batch", cells = dam_cells)
stopifnot(ncol(pcsub) == 16L, nrow(pcsub) == 200L)
expect_error(pseudobulk_counts(s4obj, "genotype_batch", cells = "no_such_cell"), "any(sel)")

# --- de_pseudobulk: 16-sample fit, 5 CI topTables, stageR matrix, finite interaction MDE -
dpb <- de_pseudobulk(s4obj)
stopifnot(dpb$n_samples == 16L, length(dpb$lib_size) == 16L, dpb$kept >= 1L,
          identical(names(dpb$top), canonical),
          all(c("gene", "CI.L", "CI.R") %in% names(dpb$top$interaction)))
sp <- dpb$stageR$stage_padj                                   # genes x {padjScreen, 5 contrasts}
stopifnot(is.matrix(sp), identical(colnames(sp), c("padjScreen", canonical)),
          nrow(sp) == dpb$kept, dpb$stageR$n_screened >= 0L,
          dpb$thresholds$fdr == 0.05, dpb$thresholds$lfc == 0.5)
stopifnot(is.finite(dpb$interaction$mde), dpb$interaction$mde > 0,
          dpb$interaction$contrast == "interaction", is.finite(dpb$interaction$df))

# --- run_pb_de_substate: floor=3 -> Homeostatic+DAM FIT, IFN(1/unit)+Proliferative(0) SKIP
rsub <- run_pb_de_substate(s4obj, min_cells = 3L)
status <- vapply(rsub$per_substate, function(x) x$status, character(1))
stopifnot(identical(unname(status[c("Homeostatic", "DAM")]), c("fit", "fit")),
          identical(unname(status[c("IFN", "Proliferative")]), c("skipped", "skipped")),
          identical(dim(rsub$cell_counts), c(4L, 16L)), rsub$min_cells == 3L,
          identical(names(rsub$per_substate$DAM$top), canonical),          # fit branch carries DE result
          grepl("floor", rsub$per_substate$IFN$reason, fixed = TRUE))      # skip branch carries a reason
# tighter floor skips DAM too (some unit has 5 < 6)
stopifnot(run_pb_de_substate(s4obj, min_cells = 6L)$per_substate$DAM$status == "skipped")

# --- dam_direction: per-amyloid-contrast direction summary over the marker genes -----------
s4map <- data.frame(ensembl = rownames(s4obj),
                    symbol = paste0("Sym", seq_len(nrow(s4obj))), stringsAsFactors = FALSE)
dconc <- dam_direction(dpb$top, s4map, markers = c("Sym10", "Sym20", "Sym30"))
stopifnot(identical(names(dconc), c("nlgf_in_maptki", "nlgf_in_p301s")),
          dconc$nlgf_in_maptki$n_markers_in_fit >= 1L,
          is.finite(dconc$nlgf_in_maptki$frac_up), is.finite(dconc$nlgf_in_maptki$mean_logFC))

cat("ok - test_de_pb\n")
