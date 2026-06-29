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

# --- P1-S2: UCell substate annotation + contaminant-cluster prune --------------------------

# Map named marker-symbol sets -> present ensembl ids via symbol_map (assay rownames are ensembl).
# Drops symbols with no ensembl hit AND ids absent from `present_ids`; a set left empty errors (a
# silent empty signature would fabricate UCell scores). Returns the named ensembl-id list, with
# attr "n_used" = "<kept>/<requested>" per set for the provenance record. Pure.
marker_sets_to_ensembl <- function(marker_sets, symbol_map, present_ids) {
  stopifnot(is.list(marker_sets), length(marker_sets) >= 1L, !is.null(names(marker_sets)),
            is.data.frame(symbol_map), all(c("symbol", "ensembl") %in% colnames(symbol_map)))
  ens <- lapply(marker_sets, function(s) intersect(symbols_to_ensembl(s, symbol_map), present_ids))
  empty <- names(ens)[vapply(ens, length, integer(1)) == 0L]
  if (length(empty))
    stop("marker_sets_to_ensembl: no present ensembl id for set(s): ", paste(empty, collapse = ", "),
         " -- map symbols -> ensembl upstream (assay rownames are ensembl)")
  n_used <- vapply(names(marker_sets),
                   function(nm) sprintf("%d/%d", length(ens[[nm]]), length(marker_sets[[nm]])), character(1))
  attr(ens, "n_used") <- n_used
  ens
}

# Calibrate a units x signatures score matrix: z-scale each signature (column) so sets of
# different size/coherence become comparable before argmax (raw UCell is NOT cross-signature
# comparable). A zero-variance column -> z undefined -> set to 0 (no enrichment). Pure.
zscale_signatures <- function(score_mat) {
  stopifnot(is.matrix(score_mat), ncol(score_mat) >= 1L)
  z <- scale(score_mat)
  z[is.nan(z)] <- 0
  attr(z, "scaled:center") <- NULL
  attr(z, "scaled:scale")  <- NULL
  dimnames(z) <- dimnames(score_mat)
  z
}

# Argmax substate assignment with explicit ambiguous/unassigned buckets (never force-assign). For
# each row (a cell or a cluster) of a z-scaled signature matrix:
#   unassigned : best signature z <= 0 (no positive enrichment for ANY substate)
#   ambiguous  : the top-two are BOTH > amb_floor (genuinely enriched) AND within tol (co-dominant)
#   else       : the argmax signature name.
# amb_floor guards against a noise-level runner-up (z just above 0, e.g. a sparse signature with a
# single stray gene) faking ambiguity -- ambiguity requires TWO real competitors. Returns a
# row-named character vector. Pure.
assign_substate <- function(z_mat, tol = 0.10, amb_floor = 0.10) {
  stopifnot(is.matrix(z_mat), ncol(z_mat) >= 2L, !is.null(colnames(z_mat)))
  apply(z_mat, 1L, function(r) {
    o <- order(r, decreasing = TRUE)
    if (r[o[1]] <= 0) return("unassigned")
    if (r[o[2]] > amb_floor && (r[o[1]] - r[o[2]]) < tol) return("ambiguous")
    colnames(z_mat)[o[1]]
  })
}

# Aggregate a per-cell z matrix to per-cluster MEAN z -- cluster-level (primary) assignment
# averages cell z within each cluster. Returns clusters x signatures, rows in the cluster factor's
# (dropped-empty) level order. Pure.
cluster_mean_z <- function(z_mat, clusters) {
  stopifnot(is.matrix(z_mat), length(clusters) == nrow(z_mat))
  cl  <- droplevels(factor(clusters))
  out <- apply(z_mat, 2L, function(x) tapply(x, cl, mean))
  if (!is.matrix(out))
    out <- matrix(out, nrow = nlevels(cl), dimnames = list(levels(cl), colnames(z_mat)))
  out
}

# Decide which clusters are non-microglial contamination from per-cluster RAW UCell stats. A
# cluster is dropped when microglia identity is essentially ABSENT (id_med < id_floor) OR a
# contaminant signature beats identity in most of its cells (mglike_frac < mglike_floor;
# doublet/ambient dominated). Both defaults sit in observed natural gaps for this data (real
# microglia clusters had id_med >= 0.158 & mglike_frac >= 0.38; contaminants <= 0.091 & <= 0.24)
# -> conservative, drops only clear outliers (no over-pruning). Pure; returns dropped cluster names.
flag_contaminant_clusters <- function(stats, id_floor = 0.15, mglike_floor = 0.30) {
  stopifnot(is.data.frame(stats), all(c("id_med", "mglike_frac") %in% colnames(stats)),
            !is.null(rownames(stats)))
  rownames(stats)[stats$id_med < id_floor | stats$mglike_frac < mglike_floor]
}

# Annotate the reprocessed microglia: UCell-score identity + substate + aux + contaminant
# signatures, drop clear contaminant clusters, assign substates on the clean population.
#   1. UCell (rank-based; robust to dropout/depth/batch) scores every signature on the SCT `layer`
#      (ncores=1 -> deterministic, no parallel warning). Raw scores -> meta `<sig>_UCell`.
#   2. Prune: per cluster, raw microglia-identity median + fraction of cells whose identity beats
#      its best contaminant signature -> flag_contaminant_clusters -> DROP (subset out). Doublets
#      are precomputed (0 here -> logged no-op). Per-cluster QC rationale (id/mglike/contam/ribo/
#      malat1/dropped), dropped ids+counts, and the dropped genotype x cluster table go in
#      @misc$microglia_prune -- nothing hidden (the dropout is mildly genotype-associated; reported).
#   3. On RETAINED cells: z-scale the 4 substate signatures (calibrate per signature on the clean
#      population) -> cluster-mean-z argmax = PRIMARY label (broadcast to its cells); per-cell-z
#      argmax = SECONDARY (diagnostic; noisier for sparse states). MHC_APC z = a continuous aux
#      axis (ARM = DAM + MHC). Ambiguous/unassigned buckets, never force-assigned.
# Thrupp 2020 caveat carried: snRNA depletes ~18% of DAM genes -> DAM is SCORED not thresholded,
# and the broad DAM signature absorbs the dropout. Returns the pruned, annotated obj (DefaultAssay
# SCT). Heavy body (UCell on the full subset) is smoke-tested live; pure helpers are unit-tested.
annotate_microglia <- function(seurat_obj, symbol_map,
                               marker_sets      = canonical_microglia_markers,
                               identity_markers = microglia_identity_markers,
                               contam_sets      = contam_signatures,
                               substate_levels  = microglia_substate_levels,
                               cluster_col = "microglia_clusters",
                               assay = "SCT", layer = "data",
                               id_floor = 0.15, mglike_floor = 0.30,
                               tol = 0.10, amb_floor = 0.10) {
  stopifnot(
    inherits(seurat_obj, "Seurat"), assay %in% SeuratObject::Assays(seurat_obj),
    cluster_col %in% colnames(seurat_obj@meta.data),
    is.list(marker_sets), all(substate_levels %in% names(marker_sets)), "MHC_APC" %in% names(marker_sets),
    is.list(contam_sets), length(contam_sets) >= 1L, !is.null(names(contam_sets)),
    is.character(identity_markers), length(identity_markers) >= 1L,
    is.data.frame(symbol_map), all(c("symbol", "ensembl") %in% colnames(symbol_map))
  )
  present  <- rownames(SeuratObject::GetAssayData(seurat_obj, assay = assay, layer = layer))
  all_sets <- c(list(Microglia_identity = identity_markers), marker_sets, contam_sets)
  ens_sets <- marker_sets_to_ensembl(all_sets, symbol_map, present)

  # UCell raw scores (ties.method default "average"; missing_genes="skip" -- sets are pre-filtered
  # to present ids so nothing is missing). slot=`layer` is gate-safe on Assay5 (no deprecation warn).
  obj <- UCell::AddModuleScore_UCell(seurat_obj, features = ens_sets, assay = assay, slot = layer,
                                     ncores = 1L, name = "_UCell", missing_genes = "skip")
  stopifnot(all(paste0(names(all_sets), "_UCell") %in% colnames(obj@meta.data)))

  # --- prune: per-cluster raw identity vs best contaminant ---
  md  <- obj@meta.data
  cl0 <- droplevels(factor(md[[cluster_col]]))
  id_raw     <- md[["Microglia_identity_UCell"]]
  contam_raw <- do.call(pmax, md[paste0(names(contam_sets), "_UCell")])
  stats <- data.frame(
    n           = as.integer(table(cl0)),
    id_med      = tapply(id_raw, cl0, median),
    mglike_frac = tapply(id_raw > contam_raw, cl0, mean),
    row.names   = levels(cl0)
  )
  drop_clusters <- flag_contaminant_clusters(stats, id_floor, mglike_floor)
  keep <- !(as.character(md[[cluster_col]]) %in% drop_clusters)
  stopifnot(any(keep), length(drop_clusters) < nlevels(cl0))   # never prune everything

  qc_rationale <- data.frame(
    n           = stats$n,
    id_med      = round(stats$id_med, 3),
    mglike_frac = round(stats$mglike_frac, 3),
    pct_contam  = round(tapply(md$percent_contam, cl0, mean), 2),
    pct_ribo    = round(tapply(md$percent_ribo,   cl0, mean), 2),
    pct_malat1  = round(tapply(md$percent_malat1, cl0, mean), 2),
    dropped     = levels(cl0) %in% drop_clusters,
    row.names   = levels(cl0)
  )
  if (is.numeric(md$doublets)) qc_rationale$doublet_mean <- round(tapply(md$doublets, cl0, mean), 4)
  prune_log <- list(
    qc_rationale = qc_rationale, dropped = drop_clusters,
    n_dropped = sum(!keep), n_retained = sum(keep),
    thresholds = c(id_floor = id_floor, mglike_floor = mglike_floor),
    dropped_by_genotype = table(genotype = md$genotype[!keep], cluster = droplevels(cl0[!keep]))
  )
  obj <- subset(obj, cells = colnames(obj)[keep])
  obj@meta.data[[cluster_col]] <- droplevels(factor(obj@meta.data[[cluster_col]]))

  # --- substate assignment on the clean population (z-scale calibrated on retained cells) ---
  z_sub <- zscale_signatures(as.matrix(obj@meta.data[paste0(substate_levels, "_UCell")]))
  colnames(z_sub) <- substate_levels
  for (s in substate_levels) obj@meta.data[[paste0(s, "_UCell_z")]] <- z_sub[, s]
  obj@meta.data[["MHC_APC_UCell_z"]] <-
    as.numeric(zscale_signatures(as.matrix(obj@meta.data["MHC_APC_UCell"])))

  cl  <- droplevels(factor(obj@meta.data[[cluster_col]]))
  cmz <- cluster_mean_z(z_sub, cl)
  cluster_label <- assign_substate(cmz, tol = tol, amb_floor = amb_floor)          # PRIMARY (cluster)
  lvls <- c(substate_levels, "ambiguous", "unassigned")
  obj@meta.data[["microglia_substate"]] <- factor(cluster_label[as.character(cl)], levels = lvls)
  obj@meta.data[["microglia_substate_percell"]] <-                                  # SECONDARY (cell)
    factor(assign_substate(z_sub, tol = tol, amb_floor = amb_floor), levels = lvls)

  # --- postconditions: labelled-or-bucketed + self-consistent enrichment ---
  sub <- obj@meta.data[["microglia_substate"]]
  stopifnot(
    !anyNA(sub),                                                       # every retained cell bucketed
    nlevels(cl) == nrow(cmz),                                          # no dropped-cluster ghosts
    all(c("pca", "harmony", "umap") %in% SeuratObject::Reductions(obj))# reductions survive subset
  )
  for (s in intersect(levels(droplevels(sub)), substate_levels)) {     # S-labelled cells score higher on S
    sc <- obj@meta.data[[paste0(s, "_UCell")]]
    stopifnot(mean(sc[sub == s]) > mean(sc[sub != s]))
  }

  obj@misc$microglia_prune <- prune_log
  obj@misc$substate_provenance <- list(
    cluster_mean_z = cmz, cluster_label = cluster_label,
    substate_table = table(genotype = obj$genotype, substate = sub),
    n_used = attr(ens_sets, "n_used"),
    thresholds = c(id_floor = id_floor, mglike_floor = mglike_floor, tol = tol, amb_floor = amb_floor),
    assay = assay, layer = layer
  )
  obj
}
