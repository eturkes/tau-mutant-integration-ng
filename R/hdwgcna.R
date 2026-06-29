# hdWGCNA pipeline for unbiased microglia co-expression modules, plus
# module-level DE (whole + per-substate) and hypergeometric enrichment of
# module gene sets against named references.
#
# All cells are treated as a single "microglia" cell type (the wgcna_group
# meta column). Metacells are pooled within state x genotype strata so the
# resulting metacells are biologically homogeneous and condition-aware,
# while module construction itself runs across the union of metacells.
# Eigengenes are batch-corrected (Harmony) at the ME stage. Downstream
# DE on modules aggregates MEs to the same genotype_batch pseudosamples
# used everywhere else in the project.

run_hdwgcna_pipeline <- function(seurat_obj,
                                 wgcna_name          = "microglia",
                                 group_by            = c("state", "genotype"),
                                 metacell_k          = 25,
                                 metacell_min_cells  = 50,
                                 metacell_max_shared = 10,
                                 metacell_reduction  = "harmony",
                                 gene_fraction       = 0.05,
                                 batch_var           = "batch",
                                 network_type        = "signed",
                                 tom_type            = "signed",
                                 min_module_size     = 50,
                                 merge_cut           = 0.2,
                                 deep_split          = 4,
                                 soft_power          = NULL,
                                 soft_power_grid     = c(1:10, seq(12, 30, by = 2)),
                                 sft_r2_target       = 0.8,
                                 seed                = 1L,
                                 verbose             = TRUE) {
  # End-to-end hdWGCNA: SetupForWGCNA -> MetacellsByGroups ->
  # NormalizeMetacells -> SetDatExpr -> TestSoftPowers (auto soft_power
  # if NULL) -> ConstructNetwork -> ModuleEigengenes -> ModuleConnectivity
  # -> ResetModuleNames. Returns the augmented Seurat plus a side-channel
  # summary with modules, hub genes, cell-level MEs and the power table.
  set.seed(seed)
  suppressPackageStartupMessages({
    library(Seurat)
    library(hdWGCNA)
    library(WGCNA)
  })
  try(WGCNA::allowWGCNAThreads(
    nThreads = max(1L, parallel::detectCores() - 2L)
  ), silent = TRUE)

  DefaultAssay(seurat_obj) <- "RNA"
  rna_layers <- SeuratObject::Layers(seurat_obj[["RNA"]])
  if (!("data" %in% rna_layers)) {
    if (verbose) message("[hdwgcna] RNA data layer missing; running NormalizeData ...")
    seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  }

  # Cell-level single-group label used by ModuleConnectivity below. Not
  # propagated to the metacell object by MetacellsByGroups; SetDatExpr
  # therefore uses one of the group_by columns instead.
  seurat_obj$wgcna_group <- "microglia"

  if (verbose) message("[hdwgcna] SetupForWGCNA (gene_select=fraction, ",
                       "fraction=", gene_fraction, ") ...")
  seurat_obj <- SetupForWGCNA(
    seurat_obj,
    gene_select = "fraction",
    fraction    = gene_fraction,
    wgcna_name  = wgcna_name
  )

  if (verbose) message("[hdwgcna] MetacellsByGroups k=", metacell_k,
                       " group_by=", paste(group_by, collapse = "/"),
                       " reduction=", metacell_reduction, " ...")
  seurat_obj <- MetacellsByGroups(
    seurat_obj  = seurat_obj,
    group.by    = group_by,
    ident.group = group_by[1L],
    k           = metacell_k,
    max_shared  = metacell_max_shared,
    min_cells   = metacell_min_cells,
    reduction   = metacell_reduction,
    assay       = "RNA",
    slot        = "counts"
  )
  seurat_obj <- NormalizeMetacells(seurat_obj)

  # Select all metacells. The metacell object only carries the columns
  # passed to MetacellsByGroups (group_by), so we use group_by[1L] with
  # its full set of levels as the selection key.
  selection_col <- group_by[1L]
  metacell_obj  <- hdWGCNA::GetMetacellObject(seurat_obj, wgcna_name)
  selection_levels <- unique(metacell_obj[[selection_col]][[1L]])
  if (verbose) message("[hdwgcna] SetDatExpr (group.by=", selection_col,
                       ", group_name=", paste(selection_levels, collapse = "/"),
                       ") ...")
  seurat_obj <- SetDatExpr(
    seurat_obj,
    group_name    = selection_levels,
    group.by      = selection_col,
    use_metacells = TRUE,
    assay         = "RNA",
    slot          = "data"
  )

  if (verbose) message("[hdwgcna] TestSoftPowers ...")
  seurat_obj <- TestSoftPowers(
    seurat_obj,
    powers      = soft_power_grid,
    networkType = network_type
  )
  pt <- GetPowerTable(seurat_obj)
  chosen_power <- soft_power
  if (is.null(chosen_power)) {
    pass <- pt[!is.na(pt$SFT.R.sq) & pt$SFT.R.sq >= sft_r2_target, , drop = FALSE]
    chosen_power <- if (nrow(pass) > 0) min(pass$Power) else
      pt$Power[which.max(pt$SFT.R.sq)]
  }
  if (verbose) message("[hdwgcna] Chosen soft_power = ", chosen_power)

  if (verbose) message("[hdwgcna] ConstructNetwork ...")
  tom_outdir <- file.path(tempdir(), paste0("hdwgcna_", wgcna_name, "_tom"))
  dir.create(tom_outdir, showWarnings = FALSE, recursive = TRUE)
  seurat_obj <- ConstructNetwork(
    seurat_obj,
    soft_power     = chosen_power,
    tom_outdir     = tom_outdir,
    tom_name       = wgcna_name,
    minModuleSize  = min_module_size,
    mergeCutHeight = merge_cut,
    deepSplit      = deep_split,
    networkType    = network_type,
    TOMType        = tom_type,
    overwrite_tom  = TRUE,
    randomSeed     = seed
  )

  # ModuleEigengenes with group.by.vars (Harmony batch correction at the
  # ME level) requires scale.data on the WGCNA features. Scale only those
  # genes for speed.
  wgcna_features <- hdWGCNA::GetWGCNAGenes(seurat_obj, wgcna_name)
  if (verbose) message("[hdwgcna] ScaleData on ", length(wgcna_features),
                       " WGCNA features ...")
  seurat_obj <- ScaleData(seurat_obj, features = wgcna_features, verbose = FALSE)

  if (verbose) message("[hdwgcna] ModuleEigengenes (group.by.vars=",
                       batch_var, ") ...")
  seurat_obj <- ModuleEigengenes(
    seurat_obj,
    group.by.vars = batch_var,
    assay         = "RNA"
  )

  if (verbose) message("[hdwgcna] ModuleConnectivity ...")
  seurat_obj <- ModuleConnectivity(
    seurat_obj,
    group.by   = "wgcna_group",
    group_name = "microglia",
    harmonized = TRUE,
    assay      = "RNA",
    slot       = "data"
  )

  seurat_obj <- ResetModuleNames(seurat_obj, new_name = "MG-M")

  modules    <- GetModules(seurat_obj, wgcna_name = wgcna_name)
  hub_genes  <- GetHubGenes(seurat_obj, n_hubs = 25, wgcna_name = wgcna_name)
  MEs_cells  <- GetMEs(seurat_obj, harmonized = TRUE, wgcna_name = wgcna_name)

  module_meta <- modules |>
    dplyr::group_by(module) |>
    dplyr::summarise(
      n_genes = dplyr::n(),
      colour  = unique(color),
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(n_genes))

  list(
    seurat_obj  = seurat_obj,
    soft_power  = chosen_power,
    power_table = pt,
    modules     = modules,
    hub_genes   = hub_genes,
    MEs_cells   = MEs_cells,
    module_meta = module_meta,
    wgcna_name  = wgcna_name
  )
}

# Aggregate cell-level module eigengenes to genotype_batch pseudosamples
# and fit limma with the project-wide 2x2 factorial design + batch fixed
# effect. Returns the canonical 5-contrast tibble schema, one row per
# (module, contrast) pair.
fit_module_de <- function(MEs_cells,
                          meta,
                          id_col       = "genotype_batch",
                          genotype_col = "genotype",
                          batch_col    = "batch",
                          drop_modules = c("grey", "MG-M0")) {
  stopifnot(nrow(MEs_cells) == nrow(meta))
  if (!is.null(rownames(MEs_cells)) && !is.null(rownames(meta))) {
    stopifnot(all(rownames(MEs_cells) == rownames(meta)))
  }
  drop_modules <- intersect(drop_modules, colnames(MEs_cells))
  ME <- as.matrix(MEs_cells[, setdiff(colnames(MEs_cells), drop_modules), drop = FALSE])

  id_vec <- meta[[id_col]]
  ME_pb <- do.call(rbind, lapply(split(seq_len(nrow(ME)), id_vec), function(idx) {
    colMeans(ME[idx, , drop = FALSE])
  }))
  storage.mode(ME_pb) <- "double"

  sample_meta <- unique(meta[, c(id_col, genotype_col, batch_col), drop = FALSE])
  rownames(sample_meta) <- sample_meta[[id_col]]
  sample_meta <- sample_meta[rownames(ME_pb), , drop = FALSE]

  fd <- factorial_design(
    sample_meta,
    genotype_col = genotype_col,
    batch_col    = batch_col
  )
  des       <- fd$design
  contrasts <- fd$contrasts

  fit  <- limma::lmFit(t(ME_pb), des)
  fit2 <- limma::contrasts.fit(fit, contrasts)
  fit2 <- limma::eBayes(fit2)

  per_contrast <- lapply(colnames(contrasts), function(ct) {
    tt <- limma::topTable(fit2, coef = ct, number = Inf, sort.by = "none")
    tibble::tibble(
      module    = rownames(tt),
      logFC     = tt$logFC,
      t         = tt$t,
      P.Value   = tt$P.Value,
      adj.P.Val = tt$adj.P.Val,
      B         = tt$B
    ) |> dplyr::arrange(P.Value)
  })
  names(per_contrast) <- colnames(contrasts)

  list(
    fit         = fit2,
    top         = per_contrast,
    ME_pb       = ME_pb,
    sample_meta = sample_meta,
    design      = des,
    contrasts   = contrasts,
    n_modules   = ncol(ME_pb),
    n_samples   = nrow(ME_pb)
  )
}

# Hypergeometric enrichment of each module's gene set against a named list
# of reference gene sets (gene symbols). Universe defaults to the union of
# module genes (i.e. the genes that survived the gene-selection filter).
module_enrichment <- function(modules, gene_sets,
                              gene_col   = "gene_name",
                              module_col = "module",
                              universe   = NULL,
                              drop_modules = "grey") {
  if (is.null(universe)) universe <- unique(modules[[gene_col]])
  modules <- modules[!(modules[[module_col]] %in% drop_modules), , drop = FALSE]
  # Cast to character so split() drops unused factor levels (otherwise the
  # grey level survives as a zero-row group and emits spurious rows).
  per_module <- split(modules[[gene_col]],
                      as.character(modules[[module_col]]))

  out <- vector("list", length = length(per_module) * length(gene_sets))
  k <- 1L
  for (mod_name in names(per_module)) {
    mod_genes <- intersect(per_module[[mod_name]], universe)
    n_mod <- length(mod_genes)
    for (gs_name in names(gene_sets)) {
      gs_genes <- intersect(gene_sets[[gs_name]], universe)
      n_gs <- length(gs_genes)
      ov <- intersect(mod_genes, gs_genes)
      n_ov <- length(ov)
      n_universe <- length(universe)
      # phyper: P(X >= n_ov) under hypergeometric null
      p <- if (n_mod == 0 || n_gs == 0) 1 else
        phyper(n_ov - 1L, n_gs, n_universe - n_gs, n_mod,
               lower.tail = FALSE)
      odds <- if (n_mod > 0 && n_gs > 0) {
        (n_ov / max(n_mod, 1)) / max(n_gs / n_universe, .Machine$double.eps)
      } else NA_real_
      out[[k]] <- tibble::tibble(
        module     = mod_name,
        gene_set   = gs_name,
        n_module   = n_mod,
        n_gene_set = n_gs,
        n_overlap  = n_ov,
        n_universe = n_universe,
        odds_ratio = odds,
        P.Value    = p,
        overlap    = paste(sort(ov), collapse = ";")
      )
      k <- k + 1L
    }
  }
  res <- dplyr::bind_rows(out)
  res$adj.P.Val <- p.adjust(res$P.Value, method = "BH")
  res |> dplyr::arrange(P.Value)
}

# Per-state module eigengene differential expression. Splits cell-level MEs
# by `state_col`, pseudobulks within each state to `id_col` (= one mean ME
# per genotype_batch within that state), and reuses `fit_module_de()` for
# the 2x2 factorial + batch limma fit. Returns a named list state -> result
# matching `fit_module_de()`'s schema, plus an extra `cells_used` count.
#
# All 16 genotype_batch ids must be populated within each state (true for
# this dataset's microglia; min 15 cells in the rarest combination). The
# function checks for missing ids and skips degenerate states with a
# message rather than erroring.
fit_module_de_per_state <- function(MEs_cells, meta,
                                     state_col    = "state",
                                     id_col       = "genotype_batch",
                                     genotype_col = "genotype",
                                     batch_col    = "batch",
                                     drop_modules = c("grey", "MG-M0"),
                                     min_ids      = 8L,
                                     verbose      = TRUE) {
  stopifnot(nrow(MEs_cells) == nrow(meta))
  stopifnot(state_col %in% names(meta))
  if (!is.null(rownames(MEs_cells)) && !is.null(rownames(meta))) {
    stopifnot(all(rownames(MEs_cells) == rownames(meta)))
  }
  states <- as.character(unique(meta[[state_col]]))
  states <- sort(states[!is.na(states) & nzchar(states)])
  res <- vector("list", length(states))
  names(res) <- states
  for (st in states) {
    sub_idx  <- which(meta[[state_col]] == st)
    sub_meta <- meta[sub_idx, , drop = FALSE]
    sub_ME   <- MEs_cells[sub_idx, , drop = FALSE]
    n_ids    <- length(unique(sub_meta[[id_col]]))
    if (verbose) {
      message(sprintf("  state=%s: %d cells across %d %s ids",
                      st, length(sub_idx), n_ids, id_col))
    }
    if (n_ids < min_ids) {
      message(sprintf("  state=%s: only %d ids (< %d); skipping",
                      st, n_ids, min_ids))
      res[[st]] <- NULL
      next
    }
    res[[st]] <- tryCatch(
      fit_module_de(MEs_cells = sub_ME,
                    meta       = sub_meta,
                    id_col     = id_col,
                    genotype_col = genotype_col,
                    batch_col  = batch_col,
                    drop_modules = drop_modules),
      error = function(e) {
        message(sprintf("  state=%s fit failed: %s",
                        st, conditionMessage(e)))
        NULL
      }
    )
    if (!is.null(res[[st]])) {
      res[[st]]$cells_used <- length(sub_idx)
      res[[st]]$state      <- st
    }
  }
  res
}
