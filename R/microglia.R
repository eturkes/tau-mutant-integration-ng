# snRNAseq microglia reprocessing + clustering (P1-S1). Pure functions consumed by the P1
# microglia targets; the DAG orders execution. Recipe carried from v1 (SCT-v2 + glmGamPoi
# backend + Harmony), REVISED to integrate over BATCH ONLY -- sex is perfectly aliased with
# batch (batch01/03 male, batch02/04 female), so batch-only absorbs sex (equivalent-or-finer).
# Genotype/amyloid are NEVER integrated over: amyloid-driven DAM activation is the biology we
# measure, not a batch effect to remove. All non-base calls are namespace-qualified (targets
# attaches only `quarto`).

# Snapshot the thread-control environment for the reproducibility provenance record. Threads
# are NOT pinned: we concede non-bitwise reproducibility (multithreaded Harmony/UMAP/BLAS), and
# the inference-relevant outputs (PCA/Harmony/Louvain clusters) are seed-deterministic given
# the SNN graph. Recorded so a re-run's threading is auditable.
reprocess_thread_env <- function() {
  thread_vars <- c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                   "VECLIB_MAXIMUM_THREADS", "RCPP_PARALLEL_NUM_THREADS")
  c(as.list(Sys.getenv(thread_vars)), list(detected_cores = parallel::detectCores()))
}

# Reprocess the RNA-counts microglia subset into an SCT-normalised, Harmony-integrated,
# multi-resolution-clustered Seurat object:
#   SCTransform(glmGamPoi, vst.flavor="v2", regress percent_mt + percent_contam)
#     -> RunPCA(npcs) -> RunHarmony(batch_col) -> FindNeighbors + RunUMAP(harmony, dims)
#     -> FindClusters(Louvain, algorithm 1) at each `resolutions`.
# The primary cluster column (resolution `primary_resolution`) is copied to a stable
# `microglia_clusters` column (decoupled from the SCT_snn_res.* naming) and set as Idents; the
# other resolutions stay in meta for S2 reconciliation. Seeds, RNGkind and the thread snapshot
# go in @misc$reprocess_provenance -- re-runs reproduce UP TO TOLERANCE, NOT bitwise (memory
# contract). Returns the processed Seurat object (DefaultAssay = SCT).
reprocess_microglia <- function(seurat_obj,
                                batch_col = "batch",
                                regress_vars = c("percent_mt", "percent_contam"),
                                npcs = 30L, dims = seq_len(20L),
                                resolutions = c(0.2, 0.4, 0.6),
                                primary_resolution = 0.4,
                                seed = 42L) {
  stopifnot(
    inherits(seurat_obj, "Seurat"),
    "RNA" %in% SeuratObject::Assays(seurat_obj),
    is.numeric(npcs), length(npcs) == 1L, npcs >= 2L, npcs == round(npcs),
    is.numeric(dims), length(dims) >= 1L, all(dims >= 1L), all(dims == round(dims)), max(dims) <= npcs,
    is.numeric(resolutions), length(resolutions) >= 1L, all(is.finite(resolutions)), all(resolutions > 0),
    length(primary_resolution) == 1L, primary_resolution %in% resolutions,
    batch_col %in% colnames(seurat_obj@meta.data),
    all(regress_vars %in% colnames(seurat_obj@meta.data)),
    !anyNA(seurat_obj@meta.data[[batch_col]]),                                  # complete grouping/regression
    all(vapply(regress_vars, function(v) !anyNA(seurat_obj@meta.data[[v]]), logical(1)))
  )
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")   # pin the stream; Seurat's internal set.seed inherits the current kind
  set.seed(seed)
  # SCTransform dispatches per-gene work through future.apply; its closures (~0.6G on 26k cells)
  # exceed future's 500MiB default globals cap -> raise it (our controlled data, bounded headroom).
  # Seurat.warn.umap.uwot=FALSE silences RunUMAP's once-per-session "default method changed" NOTICE
  # (uwot is pinned below) -- this is Seurat's own gating option, so EVERY other warning still
  # surfaces and the zero-fault gate keeps its signal.
  old_opt <- options(future.globals.maxSize = 8 * 1024^3, Seurat.warn.umap.uwot = FALSE)
  on.exit(options(old_opt), add = TRUE)

  obj <- seurat_obj
  SeuratObject::DefaultAssay(obj) <- "RNA"
  obj <- Seurat::SCTransform(
    obj, assay = "RNA", vst.flavor = "v2", method = "glmGamPoi",
    vars.to.regress = regress_vars, seed.use = seed, verbose = FALSE
  )
  obj <- Seurat::RunPCA(obj, assay = "SCT", npcs = npcs, seed.use = seed, verbose = FALSE)
  obj <- harmony::RunHarmony(   # harmony 2.0 API: assay is implicit in reduction.use (PCA on SCT); no assay.use arg
    obj, group.by.vars = batch_col,
    reduction.use = "pca", reduction.save = "harmony", verbose = FALSE
  )
  obj <- Seurat::FindNeighbors(obj, reduction = "harmony", dims = dims, verbose = FALSE)
  obj <- Seurat::RunUMAP(obj, reduction = "harmony", dims = dims, seed.use = seed,
                         umap.method = "uwot", verbose = FALSE)
  for (res in resolutions) {
    obj <- Seurat::FindClusters(obj, resolution = res, algorithm = 1L,   # 1 = Louvain (no leidenalg python dep)
                                random.seed = seed, verbose = FALSE)
  }

  # Clear stale meta "shadows" of the very outputs this function regenerates: old reduction
  # COORDINATES carried in as columns (pca1/umap1...; reprocess writes reductions, not these meta
  # cols, so the leftovers would masquerade as current) + any pre-existing cluster column at a
  # resolution we did not compute (e.g. an upstream SCT_snn_res.0.01). Fresh SCT_snn_res.<our
  # resolutions> + seurat_clusters stay; QC / annotation / cell-cycle / design columns are untouched.
  fresh_res <- paste0("SCT_snn_res.", resolutions)
  stale <- setdiff(grep("^(pca|umap|tsne|harmony)[0-9]+$|^[A-Za-z0-9]+_snn_res\\.[0-9.]+$",
                        colnames(obj@meta.data), value = TRUE), fresh_res)   # anchored: exact shadow shapes only
  if (length(stale)) obj@meta.data[stale] <- NULL

  primary_col <- paste0("SCT_snn_res.", primary_resolution)
  stopifnot(primary_col %in% colnames(obj@meta.data))
  obj$microglia_clusters <- factor(obj@meta.data[[primary_col]])
  SeuratObject::Idents(obj) <- "microglia_clusters"

  obj@misc$reprocess_provenance <- list(
    seed = seed, rngkind = RNGkind(), threads = reprocess_thread_env(),
    npcs = npcs, dims = dims, resolutions = resolutions,
    primary_resolution = primary_resolution, primary_cluster_col = primary_col,
    batch_col = batch_col, regress_vars = regress_vars,
    r_version = as.character(getRversion())
  )
  # Build-time postconditions -- the gate's warn=2 unit tests skip this heavy body, so the S1
  # acceptance INVARIANTS are pinned HERE (a regression -> tar_make stops, gate red). Structural
  # only (not the empirical "12 clusters" count, which legitimately shifts with npcs/dims tuning):
  stopifnot(
    all(c("pca", "harmony", "umap") %in% SeuratObject::Reductions(obj)),
    is.factor(obj$microglia_clusters), nlevels(obj$microglia_clusters) >= 2L,
    identical(as.character(SeuratObject::Idents(obj)), as.character(obj$microglia_clusters)),
    !any(grepl("^(pca|umap|tsne|harmony)[0-9]+$", colnames(obj@meta.data))),       # no reduction-coord shadows
    setequal(grep("_snn_res\\.", colnames(obj@meta.data), value = TRUE), fresh_res),  # only fresh cluster cols
    length(obj@misc$reprocess_provenance) >= 1L
  )
  obj
}

# Mean per-cell expression of each marker SET across the levels of a cluster column -- the
# post-Harmony substate-separation check (homeostatic/DAM/IFN must stay distinguishable; if
# they collapse, Harmony over-corrected -> lower theta or go uncorrected). Markers are matched
# to assay rownames AS GIVEN (map symbols -> ensembl upstream via symbols_to_ensembl; the
# SCT/RNA assays carry ensembl rownames). Returns a clusters x marker-set matrix: the per-cell
# mean across each set's genes, averaged within each cluster's cells, on `assay`/`layer`.
# Absent markers are dropped; a set left with no present marker errors (a silent empty set
# would fabricate NaN separation). attr "markers_used" records the retained ids. Reused by S2.
marker_mean_by_cluster <- function(seurat_obj, marker_sets,
                                   cluster_col = "microglia_clusters",
                                   assay = "SCT", layer = "data") {
  stopifnot(is.list(marker_sets), length(marker_sets) >= 1L,
            !is.null(names(marker_sets)),
            cluster_col %in% colnames(seurat_obj@meta.data))
  expr     <- SeuratObject::GetAssayData(seurat_obj, assay = assay, layer = layer)
  clusters <- as.character(seurat_obj@meta.data[[cluster_col]])
  stopifnot(length(clusters) == ncol(expr), !anyNA(clusters))
  u       <- unique(clusters)                           # deterministic row order, human-legible:
  nums    <- suppressWarnings(as.numeric(u))             # numeric-like labels -> numeric sort (0,1,..,10,11),
  lvls    <- if (!anyNA(nums)) u[order(nums)] else sort(u, method = "radix")   # else locale-independent radix
  present <- lapply(marker_sets, function(m) intersect(m, rownames(expr)))
  empty   <- names(present)[vapply(present, length, integer(1)) == 0L]
  if (length(empty))                                    # a silent empty set would fabricate NaN separation
    stop("marker_mean_by_cluster: no present marker in set(s): ", paste(empty, collapse = ", "),
         " -- map symbols -> ensembl upstream (assay rownames are ensembl)")
  out <- vapply(present, function(genes) {
    set_mean <- Matrix::colMeans(expr[genes, , drop = FALSE])   # per-cell mean across the set's genes
    vapply(lvls, function(g) mean(set_mean[clusters == g]), numeric(1))
  }, numeric(length(lvls)))
  if (!is.matrix(out)) out <- matrix(out, nrow = length(lvls))   # one cluster / one set -> re-wrap to matrix
  dimnames(out) <- list(lvls, names(marker_sets))
  stopifnot(all(is.finite(out)))                        # no NA/NaN leaked into the separation matrix
  attr(out, "markers_used") <- present
  out
}
