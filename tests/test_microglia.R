# S1 acceptance (gate-independent slice): reprocess_microglia input-guard fail-loud paths +
# marker_mean_by_cluster correctness on a synthetic fixture. The heavy SCT/Harmony/cluster body
# of reprocess_microglia is validated separately by a live smoke-test on microglia_seurat_raw
# (needs glmGamPoi-scale data); here we only exercise the cheap guards that short-circuit BEFORE
# any compute, plus the fully-deterministic separation helper.

source("R/constants.R")
source("R/io.R")          # symbols_to_ensembl (marker_sets_to_ensembl depends on it)
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
expect_error(reprocess_microglia(obj_qc, dims = integer(0)), "dims")                 # empty dims (no heavy run)
expect_error(reprocess_microglia(obj_qc, npcs = c(30L, 40L)), "npcs")               # non-scalar npcs
expect_error(reprocess_microglia(obj_qc, resolutions = c(0.2, NA)), "resolutions")   # non-finite resolution
obj_na <- obj_qc; obj_na@meta.data$batch[1] <- NA
expect_error(reprocess_microglia(obj_na), "batch_col")                               # NA in grouping metadata

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
                                    assay = "RNA", layer = "counts"), "no present marker")

# --- P1-S2 pure helpers: substate argmax + contaminant prune (gate-independent) -------------

# marker_sets_to_ensembl: map present, drop misses (unmapped OR absent from assay), record n_used
sm2 <- data.frame(symbol = c("A","B","C","D"), ensembl = c("E1","E2","E3","E4"), stringsAsFactors = FALSE)
es  <- marker_sets_to_ensembl(list(S1 = c("A","B","Z"), S2 = c("C","D")), sm2,
                              present_ids = c("E1","E2","E3","E4","E9"))
stopifnot(identical(unname(es$S1), c("E1","E2")), identical(unname(es$S2), c("E3","E4")),
          identical(attr(es, "n_used"), c(S1 = "2/3", S2 = "2/2")))
expect_error(marker_sets_to_ensembl(list(S1 = c("Z","Y")), sm2, c("E1")), "no present ensembl")
expect_error(marker_sets_to_ensembl(list(S1 = "A"), sm2, present_ids = c("E2","E3")), "no present ensembl")
# min_n floor: a set reduced to a SINGLE present id is rejected (a 1-gene signature fabricates a score)
expect_error(marker_sets_to_ensembl(list(S1 = c("A","Z","Y")), sm2, present_ids = c("E1","E9")), "near-empty")
stopifnot(length(marker_sets_to_ensembl(list(S1 = c("A","Z")), sm2,           # min_n=1 permits a singleton
                                        present_ids = c("E1","E9"), min_n = 1L)$S1) == 1L)

# zscale_signatures: per-column z (mean 0, sd 1); a zero-variance column -> 0
zz <- zscale_signatures(matrix(c(1,2,3, 5,5,5), nrow = 3, dimnames = list(NULL, c("v","const"))))
stopifnot(abs(mean(zz[,"v"])) < 1e-9, abs(sd(zz[,"v"]) - 1) < 1e-9, all(zz[,"const"] == 0),
          identical(colnames(zz), c("v","const")))
# non-finite -> 0: a stray NA cell is zeroed (no-enrichment) while the column's non-NA cells keep
# valid z (scale() na.rm's the NA from the mean/sd); a clean column is unaffected. No NA ever leaks.
znf <- zscale_signatures(matrix(c(1,2,NA, 1,2,3), nrow = 3, dimnames = list(NULL, c("a","b"))))
stopifnot(all(is.finite(znf)), znf[3,"a"] == 0, abs(znf[1,"a"] + 0.7071068) < 1e-6,
          abs(sd(znf[,"b"]) - 1) < 1e-9)

# assign_substate: clear argmax / unassigned (all<=0) / ambiguous (two real & close) / noise runner-up kept
zm <- rbind(
  clearDAM = c(Homeostatic = -0.5, DAM =  1.2, IFN = -0.2, Proliferative = -0.3),  # -> DAM
  allneg   = c(Homeostatic = -0.5, DAM = -0.2, IFN = -0.9, Proliferative = -0.4),  # -> unassigned
  tie      = c(Homeostatic =  0.55, DAM = 0.50, IFN = -0.2, Proliferative = -0.3), # gap .05<tol, both>floor -> ambiguous
  noise2nd = c(Homeostatic =  0.40, DAM = 0.05, IFN = -0.2, Proliferative = -0.3)  # runner-up <floor -> Homeostatic
)
lab <- assign_substate(zm, tol = 0.10, amb_floor = 0.10)
stopifnot(identical(unname(lab), c("DAM","unassigned","ambiguous","Homeostatic")),
          identical(names(lab), rownames(zm)))
# eps: a numerically-zero best (floating-point residual of a degenerate centred mean) -> unassigned,
# not force-labelled to an arbitrary argmax
fp <- assign_substate(rbind(noise = c(Homeostatic = 1e-12, DAM = 5e-13, IFN = -1e-12, Proliferative = 0)))
stopifnot(identical(unname(fp), "unassigned"))
# ambiguity boundary is STRICT (< tol): an exact tol gap is NOT ambiguous (uses exactly-representable
# halves to dodge the 0.3-0.2 != 0.1 float trap); a gap just under tol IS ambiguous
bnd <- assign_substate(rbind(
  atol = c(Homeostatic = 1.0, DAM = 0.5, IFN = -0.2, Proliferative = -0.3),   # gap == tol(0.5) -> Homeostatic
  utol = c(Homeostatic = 1.0, DAM = 0.6, IFN = -0.2, Proliferative = -0.3)),  # gap 0.4 < tol     -> ambiguous
  tol = 0.5, amb_floor = 0.10)
stopifnot(identical(unname(bnd), c("Homeostatic", "ambiguous")))
# finite-input contract guard (z-scale coerces non-finite upstream; the helper still fails loud on NA)
expect_error(assign_substate(rbind(x = c(Homeostatic = NA, DAM = 1, IFN = 0, Proliferative = 0))), "finite")

# cluster_mean_z: exact per-cluster column means; rows in (dropped-empty) level order
cmz <- cluster_mean_z(matrix(c(0,2, 4,6, 10,10), ncol = 2, byrow = TRUE,
                             dimnames = list(NULL, c("p","q"))), c("a","a","b"))
stopifnot(identical(rownames(cmz), c("a","b")), cmz["a","p"] == 2, cmz["b","q"] == 10)

# flag_contaminant_clusters: drop low-identity OR low-mglike; keep clean
st <- data.frame(id_med = c(0.30, 0.05, 0.40), mglike_frac = c(0.60, 0.50, 0.20),
                 row.names = c("keep","lowid","lowmgl"))
stopifnot(setequal(flag_contaminant_clusters(st, id_floor = 0.15, mglike_floor = 0.30),
                   c("lowid","lowmgl")),
          length(flag_contaminant_clusters(st[1, , drop = FALSE])) == 0L)
# finite-stat guard: an NA stat fails loud rather than returning an NA-named "ghost" drop
expect_error(flag_contaminant_clusters(
  data.frame(id_med = c(0.3, NA), mglike_frac = c(0.6, 0.5), row.names = c("a","b"))), "finite")

# annotate_microglia: input guards short-circuit before any UCell compute
sm3   <- data.frame(symbol = "A", ensembl = "E1", stringsAsFactors = FALSE)
obj_a <- make_fake_seurat(cells_per = 4L, n_genes = 30L, with_sct = TRUE)
obj_a$microglia_clusters <- factor(rep(c("c1","c2"), length.out = ncol(obj_a)))
expect_error(annotate_microglia(list(), sm3), "Seurat")
expect_error(annotate_microglia(obj_a, sm3, cluster_col = "nope"), "cluster_col")
expect_error(annotate_microglia(obj_a, "notadf"), "symbol_map")
expect_error(annotate_microglia(obj_a, sm3,
             marker_sets = list(Homeostatic = "A", DAM = "A", IFN = "A", Proliferative = "A")), "MHC_APC")
expect_error(annotate_microglia(obj_a, sm3), "percent_contam")   # required QC meta absent -> fail at the gate

# --- microglia_report_data: compact extraction from the annotated object -------------------
# Build a minimal "annotated" fixture (umap reduction + primary substate + z-scores + the two
# @misc summaries) and confirm the extractor returns a cell-aligned plotting frame + passes the
# small summaries through, with the heavy-object guards firing loud.
obj_r <- make_fake_seurat(cells_per = 4L, n_genes = 30L)         # 16 cells, genotype/batch meta
nc    <- ncol(obj_r)
emb   <- matrix(as.double(seq_len(nc * 2L)), nrow = nc,
                dimnames = list(colnames(obj_r), c("UMAP_1", "UMAP_2")))
obj_r[["umap"]] <- SeuratObject::CreateDimReducObject(embeddings = emb, key = "UMAP_", assay = "RNA")
obj_r$microglia_substate <- factor(rep(c("Homeostatic", "DAM"), length.out = nc),
                                   levels = c(microglia_substate_levels, "ambiguous", "unassigned"))
for (z in c("Homeostatic_UCell_z", "DAM_UCell_z", "MHC_APC_UCell_z"))
  obj_r@meta.data[[z]] <- as.double(seq_len(nc))
obj_r@misc$microglia_prune <- list(n_retained = nc, n_dropped = 0L)
obj_r@misc$substate_provenance <- list(
  substate_table = table(genotype = obj_r$genotype, substate = obj_r$microglia_substate))

rd <- microglia_report_data(obj_r)
stopifnot(
  is.data.frame(rd$cell_frame), nrow(rd$cell_frame) == nc,
  all(c("umap_1", "umap_2", "genotype", "substate",
        "Homeostatic_UCell_z", "DAM_UCell_z", "MHC_APC_UCell_z") %in% names(rd$cell_frame)),
  identical(levels(rd$cell_frame$genotype), genotype_levels),    # canonical genotype factor preserved
  is.factor(rd$cell_frame$substate),
  isTRUE(all.equal(rd$cell_frame$umap_1, as.double(seq_len(nc)))),   # umap col 1 carried through, order-stable
  rd$n_cells == nc, identical(rd$prune$n_retained, nc),
  !is.null(rd$provenance$substate_table)
)
# guards: missing umap / missing z column / missing substate column fail loud (no silent partial frame)
obj_nou <- obj_r; obj_nou[["umap"]] <- NULL
expect_error(microglia_report_data(obj_nou), "umap")
expect_error(microglia_report_data(obj_r, z_cols = "NOPE_UCell_z"), "z_cols")
expect_error(microglia_report_data(obj_r, substate_col = "nope"), "substate_col")

cat("ok - test_microglia\n")
