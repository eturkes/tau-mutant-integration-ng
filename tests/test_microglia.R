# S1 acceptance (gate-independent slice): reprocess_microglia input-guard fail-loud paths +
# marker_mean_by_cluster correctness on a synthetic fixture. The heavy SCT/Harmony/cluster body
# of reprocess_microglia is validated separately by a live smoke-test on microglia_seurat_raw
# (needs glmGamPoi-scale data); here we only exercise the cheap guards that short-circuit BEFORE
# any compute, plus the fully-deterministic separation helper.

source("R/constants.R")
source("R/microglia.R")
source("tests/helpers.R")

# --- reprocess_microglia: input guards short-circuit before SCTransform --------------------
obj <- make_fake_seurat(cells_per = 4L, n_genes = 30L)   # RNA assay, genotype/batch meta; NO percent_* cols

expect_error(reprocess_microglia(list()), "Seurat")                                  # not a Seurat object
expect_error(reprocess_microglia(obj, batch_col = "nope"), "batch_col")              # missing batch column
expect_error(reprocess_microglia(obj, regress_vars = "percent_mt"), "regress_vars")  # regress var absent from meta

obj_qc <- obj
obj_qc$percent_mt <- 1; obj_qc$percent_contam <- 1                                   # satisfy the regress-var guard
expect_error(reprocess_microglia(obj_qc, primary_resolution = 0.9), "primary_resolution")  # primary not in resolutions
expect_error(reprocess_microglia(obj_qc, dims = seq_len(40L)), "npcs")               # dims exceed npcs

# --- marker_mean_by_cluster: exact means, deterministic order, absent-marker handling ------
obj_m <- make_fake_seurat(cells_per = 4L, n_genes = 30L)
obj_m$microglia_clusters <- factor(rep(c("c1", "c2"), length.out = ncol(obj_m)))
rn   <- rownames(obj_m[["RNA"]])
sets <- list(A = rn[1:3], B = rn[4:6])
mm   <- marker_mean_by_cluster(obj_m, sets, cluster_col = "microglia_clusters",
                               assay = "RNA", layer = "counts")
stopifnot(is.matrix(mm), nrow(mm) == 2L, ncol(mm) == 2L,
          identical(rownames(mm), c("c1", "c2")), identical(colnames(mm), c("A", "B")))

# exact value: per-cell mean over set A's genes, averaged within cluster c1
counts    <- as.matrix(SeuratObject::GetAssayData(obj_m, assay = "RNA", layer = "counts"))
set_mean_A <- colMeans(counts[rn[1:3], , drop = FALSE])
c1_cells   <- colnames(obj_m)[obj_m$microglia_clusters == "c1"]
stopifnot(isTRUE(all.equal(unname(mm["c1", "A"]), mean(set_mean_A[c1_cells]))))

# absent markers dropped (set A keeps only the present id); a set with NO present marker errors
sets2 <- list(A = c(rn[1], "NOT_A_GENE"), B = rn[4:6])
mm2   <- marker_mean_by_cluster(obj_m, sets2, cluster_col = "microglia_clusters",
                                assay = "RNA", layer = "counts")
stopifnot(identical(attr(mm2, "markers_used")$A, rn[1]))
expect_error(marker_mean_by_cluster(obj_m, list(A = c("X", "Y"), B = rn[4:6]),
                                    cluster_col = "microglia_clusters",
                                    assay = "RNA", layer = "counts"))

cat("ok - test_microglia\n")
