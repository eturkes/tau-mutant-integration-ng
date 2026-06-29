# Microglia-specific annotation helpers: substate scoring via canonical
# markers, and the RBC-contamination flag for hub displays.

# Annotate a (sub-)microglia Seurat object with a discrete `state` factor
# in {homeostatic, DAM, IFN, proliferative} via argmax of four AddModuleScore
# columns. Identical to the snrnaseq-substate chunk's logic so per-state
# NEBULA caches stay in lock-step with the Rmd state assignments.
#
# Args:
#   seurat_obj: Seurat object (uses SCT assay).
#   symbol_map: data frame with `ensembl` and `symbol` columns matching the
#               object's row names.
#   nbin, ctrl: AddModuleScore parameters; defaults match the Rmd.
#
# Returns the Seurat object with new columns:
#   score_homeostatic, score_DAM, score_IFN, score_proliferative, state.
label_microglia_states <- function(seurat_obj, symbol_map,
                                   nbin = 12, ctrl = 50) {
  canonical_ids <- lapply(canonical_microglia_markers, symbols_to_ensembl,
                          symbol_map = symbol_map)
  hom_ids  <- intersect(canonical_ids$Microglia,     rownames(seurat_obj))
  dam_ids  <- intersect(canonical_ids$DAM,           rownames(seurat_obj))
  ifn_ids  <- intersect(canonical_ids$IFN,           rownames(seurat_obj))
  prol_ids <- intersect(canonical_ids$Proliferative, rownames(seurat_obj))
  prev_assay <- SeuratObject::DefaultAssay(seurat_obj)
  SeuratObject::DefaultAssay(seurat_obj) <- "SCT"
  seurat_obj <- Seurat::AddModuleScore(seurat_obj,
    features = list(hom_ids, dam_ids, ifn_ids, prol_ids),
    name = "module_", assay = "SCT", nbin = nbin, ctrl = ctrl
  )
  SeuratObject::DefaultAssay(seurat_obj) <- prev_assay
  score_cols   <- paste0("module_", 1:4)
  state_labels <- c("homeostatic", "DAM", "IFN", "proliferative")
  names(seurat_obj@meta.data)[match(score_cols,
    names(seurat_obj@meta.data))] <- paste0("score_", state_labels)
  ms <- seurat_obj@meta.data[, paste0("score_", state_labels)]
  seurat_obj$state <- factor(state_labels[apply(ms, 1, which.max)],
                             levels = state_labels)
  seurat_obj
}

# Mark module hub rows that match a list of contaminant symbols. Returns
# the hub data frame with an `is_likely_rbc` logical column appended and,
# optionally, a `display_rank` integer that ranks within each module after
# excluding the contaminants. Useful for surfacing publication-ready hub
# lists while keeping the raw module assignment intact.
flag_contaminant_hubs <- function(hub_df, contaminants = rbc_marker_symbols,
                                  symbol_col = "symbol",
                                  module_col = "module",
                                  rank_col   = "kME") {
  stopifnot(all(c(symbol_col, module_col, rank_col) %in% names(hub_df)))
  contaminants <- toupper(contaminants)
  hub_df$is_likely_rbc <- toupper(hub_df[[symbol_col]]) %in% contaminants
  # Per-module display rank: rank by descending `rank_col` over non-flagged
  # rows only. Flagged rows keep `display_rank = NA_integer_`.
  hub_df$display_rank <- NA_integer_
  for (m in unique(hub_df[[module_col]])) {
    idx <- which(hub_df[[module_col]] == m & !hub_df$is_likely_rbc)
    if (length(idx) == 0L) next
    ord <- idx[order(-hub_df[[rank_col]][idx])]
    hub_df$display_rank[ord] <- seq_along(ord)
  }
  hub_df
}
