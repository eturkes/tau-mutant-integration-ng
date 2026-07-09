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
# post-Harmony subpopulation-separation check (homeostatic/DAM/IFN must stay distinguishable; if
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

# --- P1-S2: UCell subpopulation annotation + contaminant-cluster prune --------------------------

# Map named marker-symbol sets -> present ensembl ids via symbol_map (assay rownames are ensembl).
# Drops symbols with no ensembl hit AND ids absent from `present_ids`. A set reduced to 0 errors
# ("no present ensembl"); a set reduced below `min_n` errors too -- a near-empty signature (e.g. the
# Thrupp 2020 ~18% snRNA DAM dropout taken to an extreme, or a wrong-organism marker list) would
# fabricate an unreliable UCell score that argmax then trusts. Returns the named ensembl-id list,
# with attr "n_used" = "<kept>/<requested>" per set for the provenance record. Pure.
marker_sets_to_ensembl <- function(marker_sets, symbol_map, present_ids, min_n = 2L) {
  stopifnot(is.list(marker_sets), length(marker_sets) >= 1L, !is.null(names(marker_sets)),
            is.data.frame(symbol_map), all(c("symbol", "ensembl") %in% colnames(symbol_map)),
            is.numeric(min_n), length(min_n) == 1L, min_n >= 1L)
  ens <- lapply(marker_sets, function(s) intersect(symbols_to_ensembl(s, symbol_map), present_ids))
  n   <- vapply(ens, length, integer(1))
  if (any(n == 0L))
    stop("marker_sets_to_ensembl: no present ensembl id for set(s): ",
         paste(names(ens)[n == 0L], collapse = ", "),
         " -- map symbols -> ensembl upstream (assay rownames are ensembl)")
  if (any(n < min_n))
    stop("marker_sets_to_ensembl: under ", min_n, " present ensembl id(s) (near-empty signature) for set(s): ",
         paste(sprintf("%s (%d)", names(ens)[n < min_n], n[n < min_n]), collapse = ", "),
         " -- widen the set or map more symbols (a single-gene signature is unreliable)")
  n_used <- vapply(names(marker_sets),
                   function(nm) sprintf("%d/%d", length(ens[[nm]]), length(marker_sets[[nm]])), character(1))
  attr(ens, "n_used") <- n_used
  ens
}

# Calibrate a units x signatures score matrix: z-scale each signature (column) so sets of
# different size/coherence become comparable before argmax (raw UCell is NOT cross-signature
# comparable). Any non-finite result -> 0 (no enrichment): a zero-variance column gives NaN, and
# a stray NA/Inf in the input would otherwise propagate to argmax as a silent mis-label. Pure.
zscale_signatures <- function(score_mat) {
  stopifnot(is.matrix(score_mat), ncol(score_mat) >= 1L)
  z <- scale(score_mat)
  z[!is.finite(z)] <- 0
  attr(z, "scaled:center") <- NULL
  attr(z, "scaled:scale")  <- NULL
  dimnames(z) <- dimnames(score_mat)
  z
}

# Argmax subpopulation assignment with explicit ambiguous/unassigned buckets (never force-assign). For
# each row (a cell or a cluster) of a z-scaled signature matrix:
#   unassigned : best signature z <= eps (no positive enrichment for ANY subpopulation; eps also absorbs
#                floating-point zeros, e.g. a degenerate single-cluster object whose centred z ~ 1e-17)
#   ambiguous  : the top-two are BOTH > amb_floor (genuinely enriched) AND within tol (co-dominant)
#   else       : the argmax signature name.
# amb_floor guards against a noise-level runner-up (z just above 0, e.g. a sparse signature with a
# single stray gene) faking ambiguity -- ambiguity requires TWO real competitors. A weak-but-clear
# argmax (best z in (eps, amb_floor]) IS assigned: argmax must pick a winner where one positively
# enriched signature dominates. Requires finite input (z-scale upstream coerces non-finite -> 0).
# Returns a row-named character vector. Pure.
assign_subpopulation <- function(z_mat, tol = 0.10, amb_floor = 0.10, eps = 1e-8) {
  stopifnot(is.matrix(z_mat), ncol(z_mat) >= 2L, !is.null(colnames(z_mat)), all(is.finite(z_mat)))
  apply(z_mat, 1L, function(r) {
    o <- order(r, decreasing = TRUE)
    if (r[o[1]] <= eps) return("unassigned")
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
            !is.null(rownames(stats)),
            all(is.finite(stats$id_med)), all(is.finite(stats$mglike_frac)))   # NA -> drop NA-named ghosts
  rownames(stats)[stats$id_med < id_floor | stats$mglike_frac < mglike_floor]
}

# Annotate the reprocessed microglia: UCell-score identity + subpopulation + aux + contaminant
# signatures, drop clear contaminant clusters, assign subpopulations on the clean population.
#   1. UCell (rank-based; robust to dropout/depth/batch) scores every signature on the SCT `layer`
#      (ncores=1 -> deterministic, no parallel warning). Raw scores -> meta `<sig>_UCell`.
#   2. Prune: per cluster, raw microglia-identity median + fraction of cells whose identity beats
#      its best contaminant signature -> flag_contaminant_clusters -> DROP (subset out). Doublets
#      are precomputed (0 here -> logged no-op). Per-cluster QC rationale (id/mglike/DAM + per-lineage
#      contaminant medians/contam/ribo/malat1/dropped), kept-vs-dropped separation margins, dropped
#      ids+counts, and the dropped genotype x cluster table go in @misc$microglia_prune -- nothing
#      hidden (the dropout is mildly genotype-associated; reported, and the dropped clusters are NOT
#      DAM-high so the amyloid->DAM headline cannot be a pruning artifact).
#   3. On RETAINED cells: z-scale the 4 subpopulation signatures (calibrate per signature on the clean
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
                               subpopulation_levels  = microglia_subpopulation_levels,
                               cluster_col = "microglia_clusters",
                               assay = "SCT", layer = "data",
                               id_floor = 0.15, mglike_floor = 0.30,
                               tol = 0.10, amb_floor = 0.10) {
  stopifnot(
    inherits(seurat_obj, "Seurat"), assay %in% SeuratObject::Assays(seurat_obj),
    cluster_col %in% colnames(seurat_obj@meta.data),
    is.list(marker_sets), all(subpopulation_levels %in% names(marker_sets)), "MHC_APC" %in% names(marker_sets),
    is.list(contam_sets), length(contam_sets) >= 1L, !is.null(names(contam_sets)),
    is.character(identity_markers), length(identity_markers) >= 1L,
    is.data.frame(symbol_map), all(c("symbol", "ensembl") %in% colnames(symbol_map)),
    all(c("percent_contam", "percent_ribo", "percent_malat1", "genotype") %in%   # prune-log/provenance inputs:
          colnames(seurat_obj@meta.data))                                        # fail at the gate, not mid-compute
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
    id_med      = round(stats$id_med, 4),                       # 4dp: the binding kept margin is ~0.008
    mglike_frac = round(stats$mglike_frac, 4),
    DAM_med     = round(tapply(md$DAM_UCell, cl0, median), 3),  # dropped clusters are NOT DAM-high (audit)
    pct_contam  = round(tapply(md$percent_contam, cl0, mean), 2),
    pct_ribo    = round(tapply(md$percent_ribo,   cl0, mean), 2),
    pct_malat1  = round(tapply(md$percent_malat1, cl0, mean), 2),
    dropped     = levels(cl0) %in% drop_clusters,
    row.names   = levels(cl0)
  )
  # Per-cluster contaminant-lineage medians -> evidence the lineage call for each dropped cluster
  # (a neuron-doublet cluster shows an elevated Neuron median over the pervasive ambient background).
  qc_rationale <- cbind(qc_rationale,
    as.data.frame(sapply(names(contam_sets),
      function(nm) round(tapply(md[[paste0(nm, "_UCell")]], cl0, median), 3))))
  if (is.numeric(md$doublets)) qc_rationale$doublet_mean <- round(tapply(md$doublets, cl0, mean), 4)
  # Exact separation between kept and dropped clusters at the two thresholds -- makes the "thresholds
  # sit in a gap" claim reproducible from the target. The prune is an OR rule (drop if id_med OR
  # mglike below floor), so a cluster dropped on one axis can sit high on the other (a neuron-doublet
  # cluster keeps a decent id_med but loses on mglike). binding_kept_margin = the smallest distance any
  # KEPT cluster has to the boundary; a small value (this data: ~0.008, the low-identity IFN cluster)
  # flags that the exact id_floor is load-bearing and should be re-derived per dataset, not assumed.
  kept_lvl   <- rownames(stats)[!(rownames(stats) %in% drop_clusters)]
  separation <- c(
    min_kept_id_med = min(stats[kept_lvl, "id_med"]),
    min_kept_mglike = min(stats[kept_lvl, "mglike_frac"]),
    max_drop_id_med = if (length(drop_clusters)) max(stats[drop_clusters, "id_med"]) else NA_real_,
    max_drop_mglike = if (length(drop_clusters)) max(stats[drop_clusters, "mglike_frac"]) else NA_real_,
    binding_kept_margin = min(min(stats[kept_lvl, "id_med"]) - id_floor,
                              min(stats[kept_lvl, "mglike_frac"]) - mglike_floor)
  )
  prune_log <- list(
    qc_rationale = qc_rationale, dropped = drop_clusters,
    n_dropped = sum(!keep), n_retained = sum(keep),
    thresholds = c(id_floor = id_floor, mglike_floor = mglike_floor),
    separation = separation,
    dropped_by_genotype = table(genotype = md$genotype[!keep], cluster = droplevels(cl0[!keep]))
  )
  obj <- subset(obj, cells = colnames(obj)[keep])
  obj@meta.data[[cluster_col]] <- droplevels(factor(obj@meta.data[[cluster_col]]))

  # --- subpopulation assignment on the clean population (z-scale calibrated on retained cells) ---
  z_sub <- zscale_signatures(as.matrix(obj@meta.data[paste0(subpopulation_levels, "_UCell")]))
  colnames(z_sub) <- subpopulation_levels
  for (s in subpopulation_levels) obj@meta.data[[paste0(s, "_UCell_z")]] <- z_sub[, s]
  obj@meta.data[["MHC_APC_UCell_z"]] <-
    as.numeric(zscale_signatures(as.matrix(obj@meta.data["MHC_APC_UCell"])))

  cl  <- droplevels(factor(obj@meta.data[[cluster_col]]))
  cmz <- cluster_mean_z(z_sub, cl)
  cluster_label <- assign_subpopulation(cmz, tol = tol, amb_floor = amb_floor)          # PRIMARY (cluster)
  lvls <- c(subpopulation_levels, "ambiguous", "unassigned")
  obj@meta.data[["microglia_subpopulation"]] <- factor(cluster_label[as.character(cl)], levels = lvls)
  obj@meta.data[["microglia_subpopulation_percell"]] <-                                  # SECONDARY (cell)
    factor(assign_subpopulation(z_sub, tol = tol, amb_floor = amb_floor), levels = lvls)

  # --- postconditions: every retained cell labelled-or-bucketed; reductions survive the subset ---
  sub <- obj@meta.data[["microglia_subpopulation"]]
  stopifnot(
    !anyNA(sub),                                                       # every retained cell bucketed
    all(c("pca", "harmony", "umap") %in% SeuratObject::Reductions(obj))# reductions survive subset
  )
  # Self-CONSISTENCY guard, not independent validation: the labels derive from these same z-scaled
  # scores, so this only catches a sign/indexing inversion (s-labelled cells must out-score non-s
  # cells on the raw s signature). It cannot detect wrong marker definitions or overfit thresholds.
  for (s in intersect(levels(droplevels(sub)), subpopulation_levels)) {
    sc <- obj@meta.data[[paste0(s, "_UCell")]]
    if (any(sub != s)) stopifnot(mean(sc[sub == s]) > mean(sc[sub != s]))   # skip when s is the only state
  }

  obj@misc$microglia_prune <- prune_log
  obj@misc$subpopulation_provenance <- list(
    cluster_mean_z = cmz, cluster_label = cluster_label,
    subpopulation_table = table(genotype = obj$genotype, subpopulation = sub),
    n_used = attr(ens_sets, "n_used"),
    thresholds = c(id_floor = id_floor, mglike_floor = mglike_floor, tol = tol, amb_floor = amb_floor),
    assay = assay, layer = layer
  )
  obj
}

# --- P1-S5: compact report-data extraction (keeps the gate render cheap) --------------------

# Per-subpopulation marker-expression panel -- the "genes that DEFINE each subpopulation" dot-plot data.
# For every marker SET (signature) and each of its genes present in `assay`, compute the mean
# expression and the fraction of cells expressing (layer value > 0) WITHIN each subpopulation group.
# Symbols map -> ensembl via symbol_map (the SCT/RNA rownames are ensembl); genes absent from the
# assay drop, a signature left under `min_present` errors (a near-empty set would misrepresent the
# state). Subpopulation groups default to the set names (Homeostatic/DAM/IFN) and must each carry >= 1
# cell in `subpopulation_col`. Returns a long data.frame {signature, gene(symbol), ensembl, subpopulation,
# n_cells, mean_expr, pct_expr}; signature/gene/subpopulation are factors ordered as given so the qmd
# renders a clean signature-blocked, marker-ordered dot grid. All numeric finite, pct in [0,1].
# Pure: no RNG, no I/O. Feeds microglia_report$subpopulation_markers -> the Subpopulation-landscape dot plot.
subpopulation_marker_panel <- function(seurat_obj, symbol_map, marker_sets,
                                  subpopulations = names(marker_sets),
                                  subpopulation_col = "microglia_subpopulation",
                                  assay = "SCT", layer = "data", min_present = 2L) {
  stopifnot(
    inherits(seurat_obj, "Seurat"),
    is.data.frame(symbol_map), all(c("symbol", "ensembl") %in% colnames(symbol_map)),
    is.list(marker_sets), length(marker_sets) >= 1L, !is.null(names(marker_sets)),
    is.character(subpopulations), length(subpopulations) >= 1L, !anyDuplicated(subpopulations),
    subpopulation_col %in% colnames(seurat_obj@meta.data)
  )
  expr        <- SeuratObject::GetAssayData(seurat_obj, assay = assay, layer = layer)
  present_ids <- rownames(expr)
  sub         <- as.character(seurat_obj@meta.data[[subpopulation_col]])
  stopifnot(length(sub) == ncol(expr))
  miss <- setdiff(subpopulations, unique(sub))
  if (length(miss))                                    # a requested state with no cells would give NaN means
    stop("subpopulation_marker_panel: subpopulation(s) absent from ", subpopulation_col, ": ",
         paste(miss, collapse = ", "))
  pieces <- lapply(names(marker_sets), function(sig) {
    hit <- symbols_to_ensembl(marker_sets[[sig]], symbol_map)   # named: names=symbol, values=ensembl (input order)
    hit <- hit[hit %in% present_ids]                            # keep only genes present in the assay
    if (length(hit) < min_present)                              # a single-gene panel row misrepresents the signature
      stop("subpopulation_marker_panel: under ", min_present, " present gene(s) for signature '", sig,
           "' -- widen the set or map more symbols")
    set_mat <- expr[unname(hit), , drop = FALSE]                # genes x cells for this signature
    do.call(rbind, lapply(subpopulations, function(st) {
      m <- set_mat[, sub == st, drop = FALSE]
      data.frame(
        signature = sig,
        gene      = names(hit),
        ensembl   = unname(hit),
        subpopulation  = st,
        n_cells   = ncol(m),
        mean_expr = as.numeric(Matrix::rowMeans(m)),
        pct_expr  = as.numeric(Matrix::rowMeans(m > 0)),
        stringsAsFactors = FALSE
      )
    }))
  })
  out <- do.call(rbind, pieces)
  out$signature <- factor(out$signature, levels = names(marker_sets))
  out$subpopulation  <- factor(out$subpopulation, levels = subpopulations)
  out$gene      <- factor(out$gene, levels = unique(out$gene))   # signature-then-marker order preserved
  rownames(out) <- NULL
  stopifnot(
    all(is.finite(out$mean_expr)), all(is.finite(out$pct_expr)),
    all(out$pct_expr >= 0 & out$pct_expr <= 1), all(out$n_cells > 0L)
  )
  out
}

# Extract ONLY what _microglia.qmd plots from the ~612MB annotated Seurat, so a
# report render reads one compact target instead of the full object. The bundle
# carries the per-cell UMAP/score frame, the unit-composition bars, and the
# subpopulation marker panel. Pure: no RNG, no I/O.
microglia_report_data <- function(seurat_obj, symbol_map,
                                  subpopulation_col = "microglia_subpopulation",
                                  z_cols = c("Homeostatic_UCell_z", "DAM_UCell_z", "MHC_APC_UCell_z"),
                                  marker_sets = canonical_microglia_markers[c("Homeostatic", "DAM", "IFN")],
                                  marker_subpopulations = c("Homeostatic", "DAM", "IFN"),
                                  marker_layer = "data") {
  stopifnot(
    inherits(seurat_obj, "Seurat"),
    is.data.frame(symbol_map),
    "umap" %in% SeuratObject::Reductions(seurat_obj),
    subpopulation_col %in% colnames(seurat_obj@meta.data),
    "genotype" %in% colnames(seurat_obj@meta.data),
    "genotype_batch" %in% colnames(seurat_obj@meta.data),
    "batch" %in% colnames(seurat_obj@meta.data),
    all(z_cols %in% colnames(seurat_obj@meta.data)),
    !is.null(seurat_obj@misc$microglia_prune),
    !is.null(seurat_obj@misc$subpopulation_provenance)
  )
  md  <- seurat_obj@meta.data
  emb <- SeuratObject::Embeddings(seurat_obj, "umap")
  stopifnot(identical(rownames(emb), rownames(md)), ncol(emb) >= 2L)   # umap cell-aligned to meta
  sub <- md[[subpopulation_col]]
  observed_subpopulations <- unique(as.character(sub))
  sub_levels <- intersect(c(microglia_subpopulation_levels, "ambiguous", "unassigned"), observed_subpopulations)
  if (!length(sub_levels)) sub_levels <- sort(observed_subpopulations, method = "radix")
  cell_frame <- data.frame(
    umap_1   = as.numeric(emb[, 1]),                      # as.numeric strips cell names -> no row.names inference
    umap_2   = as.numeric(emb[, 2]),
    genotype = factor(as.character(md$genotype), levels = genotype_levels),
    subpopulation = factor(as.character(sub), levels = sub_levels),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  cell_frame[z_cols] <- md[, z_cols, drop = FALSE]         # append the activation z-scores by name (cell-aligned)
  rownames(cell_frame) <- NULL
  unit_levels <- sort(unique(as.character(md$genotype_batch)), method = "radix")
  unit_tab <- as.data.frame(
    table(
      genotype_batch = factor(as.character(md$genotype_batch), levels = unit_levels),
      subpopulation = factor(as.character(sub), levels = sub_levels)
    ),
    stringsAsFactors = FALSE
  )
  names(unit_tab)[3L] <- "n_cells"
  unit_tab$n_cells <- as.numeric(unit_tab$n_cells)
  unit_meta <- unique(data.frame(
    genotype_batch = as.character(md$genotype_batch),
    genotype = factor(as.character(md$genotype), levels = genotype_levels),
    batch = as.character(md$batch),
    stringsAsFactors = FALSE
  ))
  stopifnot(nrow(unit_meta) == length(unit_levels), !anyNA(unit_meta$genotype),
            !anyDuplicated(unit_meta$genotype_batch))
  unit_comp <- merge(unit_tab, unit_meta, by = "genotype_batch", all.x = TRUE, sort = FALSE)
  totals <- stats::aggregate(n_cells ~ genotype_batch, data = unit_comp, FUN = sum)
  names(totals)[2L] <- "unit_total"
  unit_comp <- merge(unit_comp, totals, by = "genotype_batch", all.x = TRUE, sort = FALSE)
  unit_comp$proportion <- ifelse(unit_comp$unit_total > 0,
                                 unit_comp$n_cells / unit_comp$unit_total, NA_real_)
  unit_comp$genotype <- factor(as.character(unit_comp$genotype), levels = genotype_levels)
  unit_comp <- unit_comp[order(match(unit_comp$genotype, genotype_levels),
                               unit_comp$genotype_batch, unit_comp$subpopulation,
                               method = "radix"), , drop = FALSE]
  rownames(unit_comp) <- NULL
  prov  <- seurat_obj@misc$subpopulation_provenance
  prune <- seurat_obj@misc$microglia_prune
  # cell_frame is the single source of truth -> assert the passed-through summaries AGREE with it
  # (catches drift between the S2/S3 provenance table and this S5 per-cell frame).
  sub_counts <- table(factor(as.character(sub), levels = colnames(prov$subpopulation_table)))
  stopifnot(
    !anyNA(cell_frame$genotype), !anyNA(cell_frame$subpopulation),   # every cell placed (annotate guarantees it)
    all(is.finite(cell_frame$umap_1)), all(is.finite(cell_frame$umap_2)),   # finite coords -> no ggplot drop
    all(vapply(z_cols, function(z) all(is.finite(cell_frame[[z]])), logical(1))),
    all(is.finite(unit_comp$n_cells)), all(is.finite(unit_comp$unit_total)),
    all(is.finite(unit_comp$proportion)), !anyNA(unit_comp$genotype),
    !anyNA(unit_comp$batch),
    identical(as.integer(colSums(prov$subpopulation_table)), as.integer(sub_counts)),  # provenance == per-cell counts
    isTRUE(prune$n_retained == ncol(seurat_obj))               # retained count == frame rows
  )
  # Genes-that-define-each-subpopulation dot-plot data (compact: ~40 genes x 3 subpopulations). Computed
  # here (the sole point with the heavy annotated object + its SCT layer) so the qmd reads it off
  # the ~0.5MB target like cell_frame -- no extra heavy load at render.
  subpopulation_markers <- subpopulation_marker_panel(
    seurat_obj, symbol_map, marker_sets = marker_sets, subpopulations = marker_subpopulations,
    subpopulation_col = subpopulation_col, assay = "SCT", layer = marker_layer)
  list(cell_frame = cell_frame, n_cells = ncol(seurat_obj), prune = prune, provenance = prov,
       subpopulation_markers = subpopulation_markers, unit_composition = unit_comp)
}
