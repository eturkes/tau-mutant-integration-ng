# --------------------------------------------------------------------
# Cell-type specificity of the tau x amyloid interaction (plan arc N).
# Pure functions that ask the project's first CROSS-CELL-TYPE question:
# is the locked interaction (NLGF_P301S-P301S)-(NLGF_MAPTKI-MAPTKI) -- and
# the amyloid / tau main effects -- microglia-SPECIFIC, or do astrocytes /
# oligodendrocytes / OPCs / neurons / vascular cells show comparable
# factorial responses? Every prior arc (D..M) reads ONLY microglia; this
# layer re-runs the EXACT microglia DE treatment (rmd/02b: pseudobulk
# limma-voom + single-cell NEBULA) symmetrically over 6 units (5 broad
# non-microglial types + microglia-as-whole, the 4 substates pooled).
#
# The I/O driver is scripts/build_celltype_specificity.R (step N3); the
# display chapter is rmd/23 (N4). NO function here writes to disk.
#
# THE CONFOUND these helpers are built around: cell count varies ~25x
# across units (Neuronal ~163k vs the IFN-thin microglia), so pseudobulk
# precision -> power is asymmetric and a specificity claim from NATIVE
# counts alone is confounded. The load-bearing control is the MATCHED
# regime: downsample_balanced() equalises cells per unit x replicate to a
# common K (the min-over-units of the per-replicate minimum) BEFORE refit,
# so the headline specificity verdict rests on power-equalised fits. NATIVE
# fits are reported too but flagged power-confounded (guardrail).
#
# Locked design wired in (N0/N1, verified against the live cache):
#  * Substrate = seurat_full_processed.rds (286,285 cells); metadata carries
#    genotype / batch / sex / genotype_batch (16 ids, 4/genotype) / cell_type
#    (9-level: Microglia_<state> for microglia, broad label otherwise).
#  * 6 units = Astrocyte, Neuronal, Oligodendrocyte, OPC, Vascular, and
#    Microglia (the 4 Microglia_* labels pooled). All 6 carry all 16
#    replicates, so pseudobulk is well-powered for every unit.
#  * TWO estimators (pseudobulk + NEBULA, the N1 full stack) x TWO power
#    regimes (NATIVE + MATCHED), all 6 units, the locked 5 canonical
#    contrasts. Both estimators mirror rmd/02b EXACTLY so the 6 units are
#    method-comparable (subset_pseudobulk_de == 02b's pseudobulk chunk;
#    subset_nebula_de reuses the GENERIC fit_nebula_microglia per subset).
#
# Anti-anchoring guardrails enforced downstream (state at EVERY readout):
# all 6 units reported SYMMETRICALLY with no privileged ordering
# (alphabetical); the specificity verdict MUST cite the MATCHED result (a
# NATIVE-only headline is confounded and disallowed); a non-microglial
# interaction hit -- especially Astrocyte -- is a FINDING, not noise to
# explain away; microglia-restriction VALIDATES the focus, a broader
# response REFRAMES (not refutes) it.
# --------------------------------------------------------------------

# Canonical contrast order used by every table/plot here (interaction last,
# the headline column). Mirrors the project-wide 5-contrast set.
specificity_contrasts <- c("tau_alone", "nlgf_in_maptki", "nlgf_in_p301s",
                           "tau_in_nlgf", "interaction")

# Pool the 9-level `cell_type` into the 6 comparison UNITS: every
# Microglia_<state> label collapses to "Microglia"; the 5 broad
# non-microglial labels pass through unchanged. Returns a character vector
# the same length as `cell_type`.
celltype_unit_labels <- function(cell_type) {
  ct <- as.character(cell_type)
  ifelse(grepl("^Microglia", ct), "Microglia", ct)
}

# Per-unit minimum cells-per-replicate: a (unit) named integer vector whose
# value is the smallest cell count that unit has in any single genotype_batch
# replicate. The downsample target K is min() of this vector (the largest
# common per-replicate depth every unit can supply). Used to FIX K from data
# before any fit (anti-anchoring: K is data-determined, not tuned to results).
min_cells_per_replicate <- function(sc, unit_col = "unit",
                                    id_col = "genotype_batch") {
  tab <- table(as.character(sc@meta.data[[unit_col]]),
               as.character(sc@meta.data[[id_col]]))
  apply(tab, 1, min)
}

# --- Estimator A: pseudobulk limma-voom on one unit (mirrors rmd/02b) -----
#
# Subset `sc` to the cells of `unit`, aggregate raw counts per genotype_batch
# (16 pseudo-samples), and fit limma-voom with the design
#   ~ 0 + genotype + batch
# and the 5 canonical contrasts -- byte-for-byte the rmd/02b snrnaseq-de
# chunk, so the 6 units are directly comparable to the microglia headline.
# Returns list(meta, design, contrasts, fit, top, n_cells, n_genes); `top`
# is the per-contrast tibble list (lifted to the top level so pseudobulk and
# NEBULA share the `$top[[contrast]]` accessor downstream).
subset_pseudobulk_de <- function(sc, unit, unit_col = "unit",
                                 min_count = 10, symbol_map = NULL) {
  cells <- colnames(sc)[as.character(sc@meta.data[[unit_col]]) == unit]
  stopifnot(length(cells) > 0)
  sub <- subset(sc, cells = cells)

  pb <- build_pseudobulk(sub, sample_col = "genotype_batch",
                         covariate_cols = c("genotype", "batch", "sex"))
  meta <- pb$meta
  meta$genotype <- factor(as.character(meta$genotype), levels = genotype_levels)
  meta$batch    <- factor(meta$batch)
  design <- model.matrix(~ 0 + genotype + batch, data = meta)
  colnames(design) <- sub("^genotype", "", colnames(design))
  colnames(design) <- gsub("[^A-Za-z0-9_]", "_", colnames(design))
  cm  <- make_contrast_matrix(design)
  fit <- fit_limma_voom(pb$counts, group = meta$genotype, design = design,
                        contrasts = cm, min_count = min_count)
  # seurat_full_processed is SYMBOL-keyed (rownames are MGI symbols), so when
  # no ensembl->symbol map is supplied the `gene` column IS the symbol. Pass
  # a map only for an ensembl-keyed object.
  for (cn in names(fit$top)) {
    fit$top[[cn]] <- if (!is.null(symbol_map)) {
      dplyr::left_join(fit$top[[cn]], symbol_map, by = c("gene" = "ensembl"))
    } else {
      dplyr::mutate(fit$top[[cn]], symbol = .data$gene)
    }
  }
  list(meta = meta, design = design, contrasts = cm, fit = fit,
       top = fit$top, n_cells = ncol(sub), n_genes = fit$kept)
}

# --- Estimator B: single-cell NEBULA on one unit (mirrors rmd/02b) --------
#
# Subset `sc` to the cells of `unit` and fit the GENERIC NEBULA factorial
# (fit_nebula_microglia is verified unit-agnostic despite the name: it builds
# the 2x2 design from `genotype` internally). Returns the standard
# fit_nebula_microglia shape: list(fit, top, design, n_genes, n_cells) with
# `$top[[contrast]]` matching the pseudobulk wrapper's accessor.
subset_nebula_de <- function(sc, unit, unit_col = "unit",
                             min_cell_frac = 0.01,
                             ncore = max(1L, parallel::detectCores() - 2L),
                             symbol_map = NULL) {
  cells <- colnames(sc)[as.character(sc@meta.data[[unit_col]]) == unit]
  stopifnot(length(cells) > 0)
  sub <- subset(sc, cells = cells)
  fit <- fit_nebula_microglia(sub, id_col = "genotype_batch",
                              genotype_col = "genotype", assay = "RNA",
                              layer = "counts", min_cell_frac = min_cell_frac,
                              ncore = ncore, symbol_map = symbol_map)
  # Symbol-keyed substrate: with no map the `gene` column already holds the
  # symbol (assemble_nebula_top left `symbol` NA). Mirror the pseudobulk path.
  if (is.null(symbol_map)) {
    for (cn in names(fit$top)) fit$top[[cn]]$symbol <- fit$top[[cn]]$gene
  }
  fit
}

# --- Power-equalising downsample (the MATCHED regime) ---------------------
#
# Sample (without replacement) exactly K cells from every unit x genotype_batch
# stratum, so all 6 units carry the same per-replicate depth and pseudobulk /
# NEBULA precision is no longer cell-count-confounded. Strata with <= K cells
# are taken whole (should not occur when K is the data-derived min). A single
# set.seed(seed) before a fixed-order (sorted unit, sorted id) loop makes the
# draw fully reproducible; K and seed are recorded by the caller. Returns the
# downsampled Seurat object.
downsample_balanced <- function(sc, K, seed = 1L, unit_col = "unit",
                                id_col = "genotype_batch", units = NULL) {
  meta  <- sc@meta.data
  u     <- as.character(meta[[unit_col]])
  g     <- as.character(meta[[id_col]])
  cells <- colnames(sc)
  if (is.null(units)) units <- sort(unique(u))
  ids <- sort(unique(g))

  set.seed(seed)
  keep <- vector("list", length(units) * length(ids))
  i <- 0L
  for (unit in units) {
    for (gid in ids) {
      i <- i + 1L
      idx <- which(u == unit & g == gid)
      if (length(idx) == 0L) next
      take <- if (length(idx) <= K) idx else sample(idx, K)
      keep[[i]] <- cells[take]
    }
  }
  keep_cells <- unlist(keep, use.names = FALSE)
  subset(sc, cells = keep_cells)
}

# --- R1 / R2: significant-gene tally --------------------------------------
#
# Per (unit x contrast) significant-gene counts at FDR<0.05 and FDR<0.10,
# up/down split, plus cells/genes tested. `fits_by_unit` is a named list
# unit -> fit object (either wrapper's output; both expose $top / $n_cells).
# `regime` ("native"/"matched") and `estimator` ("pseudobulk"/"nebula") are
# stamped as columns so tallies stack across the full 2x2x6 grid. The
# `interaction` rows are the specificity headline; the verdict reads the
# MATCHED ones.
specificity_tally <- function(fits_by_unit, regime, estimator,
                              contrasts = specificity_contrasts) {
  purrr::imap_dfr(fits_by_unit, function(fo, unit) {
    purrr::map_dfr(contrasts, function(cn) {
      tt <- fo$top[[cn]]
      sig05 <- tt$adj.P.Val < 0.05
      sig10 <- tt$adj.P.Val < 0.10
      data.frame(
        regime = regime, estimator = estimator, unit = unit, contrast = cn,
        n_cells        = fo$n_cells %||% NA_integer_,
        n_genes_tested = nrow(tt),
        n_sig_05  = sum(sig05, na.rm = TRUE),
        n_up_05   = sum(sig05 & tt$logFC > 0, na.rm = TRUE),
        n_down_05 = sum(sig05 & tt$logFC < 0, na.rm = TRUE),
        n_sig_10  = sum(sig10, na.rm = TRUE),
        n_up_10   = sum(sig10 & tt$logFC > 0, na.rm = TRUE),
        n_down_10 = sum(sig10 & tt$logFC < 0, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    })
  })
}

# --- R3a: interaction-logFC concordance vs the microglia reference --------
#
# Genome-wide Spearman of `contrast` logFC between the `reference` unit
# (microglia) and every other unit, on shared genes, per estimator. A LOW rho
# means a unit's interaction response is genuinely distinct from microglia's;
# a HIGH rho means it tracks microglia. Returns long rows tagged
# comparison="vs_microglia".
interaction_concordance <- function(fits_by_unit, estimator, regime,
                                    reference = "Microglia",
                                    contrast = "interaction") {
  if (!reference %in% names(fits_by_unit)) {
    warning(sprintf("reference unit '%s' absent; concordance skipped.", reference))
    return(data.frame())
  }
  ref <- fits_by_unit[[reference]]$top[[contrast]] |>
    dplyr::transmute(gene, logFC_ref = logFC)
  others <- setdiff(names(fits_by_unit), reference)
  purrr::map_dfr(others, function(u) {
    cmp <- fits_by_unit[[u]]$top[[contrast]] |>
      dplyr::transmute(gene, logFC_u = logFC) |>
      dplyr::inner_join(ref, by = "gene") |>
      dplyr::filter(is.finite(logFC_ref), is.finite(logFC_u))
    rho <- if (nrow(cmp) > 2)
      suppressWarnings(cor(cmp$logFC_ref, cmp$logFC_u, method = "spearman"))
      else NA_real_
    data.frame(regime = regime, estimator = estimator,
               comparison = "vs_microglia", unit = u, contrast = contrast,
               n_shared = nrow(cmp), spearman = rho, stringsAsFactors = FALSE)
  })
}

# --- R3c: cross-estimator (pseudobulk vs NEBULA) concordance per unit ------
#
# The rmd/02b "pseudobulk vs NEBULA power gain" panel generalised to all 6
# units: Spearman of `contrast` logFC between the two estimators within each
# unit, on shared genes. Returns long rows tagged comparison="pb_vs_nebula"
# (estimator="cross") so they stack with interaction_concordance().
cross_estimator_concordance <- function(pb_by_unit, neb_by_unit, regime,
                                        contrast = "interaction") {
  units <- intersect(names(pb_by_unit), names(neb_by_unit))
  purrr::map_dfr(units, function(u) {
    pb <- pb_by_unit[[u]]$top[[contrast]]  |> dplyr::transmute(gene, logFC_pb = logFC)
    nb <- neb_by_unit[[u]]$top[[contrast]] |> dplyr::transmute(gene, logFC_nb = logFC)
    cmp <- dplyr::inner_join(pb, nb, by = "gene") |>
      dplyr::filter(is.finite(logFC_pb), is.finite(logFC_nb))
    rho <- if (nrow(cmp) > 2)
      suppressWarnings(cor(cmp$logFC_pb, cmp$logFC_nb, method = "spearman"))
      else NA_real_
    data.frame(regime = regime, estimator = "cross",
               comparison = "pb_vs_nebula", unit = u, contrast = contrast,
               n_shared = nrow(cmp), spearman = rho, stringsAsFactors = FALSE)
  })
}

# --- R4: per-gene specificity class ---------------------------------------
#
# For every gene significant at `contrast` (FDR<fdr) in ANY unit, record the
# set of units in which it is significant and classify it:
#   microglia_unique      -- sig in microglia only
#   non_microglial_unique -- sig in >=1 non-microglial unit but NOT microglia
#   shared                -- sig in microglia AND >=1 other unit
# Returns one row per gene (gene, symbol, units_sig string, n_units_sig,
# in_microglia, class) plus the regime/estimator/contrast stamps. An empty
# data.frame is returned when nothing is significant.
specificity_class_table <- function(fits_by_unit, regime, estimator,
                                    contrast = "interaction", fdr = 0.10,
                                    microglia_unit = "Microglia") {
  units <- names(fits_by_unit)
  sig_by_unit <- lapply(units, function(u) {
    tt <- fits_by_unit[[u]]$top[[contrast]]
    tt$gene[which(tt$adj.P.Val < fdr)]
  })
  names(sig_by_unit) <- units
  all_sig <- unique(unlist(sig_by_unit, use.names = FALSE))
  if (length(all_sig) == 0)
    return(data.frame(regime = character(), estimator = character(),
                      contrast = character(), gene = character(),
                      symbol = character(), units_sig = character(),
                      n_units_sig = integer(), in_microglia = logical(),
                      class = character(), stringsAsFactors = FALSE))

  membership <- vapply(units, function(u) all_sig %in% sig_by_unit[[u]],
                       logical(length(all_sig)))
  rownames(membership) <- all_sig
  in_mg    <- if (microglia_unit %in% units) membership[, microglia_unit] else rep(FALSE, length(all_sig))
  others   <- setdiff(units, microglia_unit)
  in_other <- if (length(others)) rowSums(membership[, others, drop = FALSE]) > 0
              else rep(FALSE, length(all_sig))
  cls <- ifelse(in_mg & !in_other, "microglia_unique",
         ifelse(!in_mg & in_other, "non_microglial_unique", "shared"))

  units_sig <- apply(membership, 1, function(r) paste(units[r], collapse = ","))

  # Resolve a symbol per gene by coalescing the symbol column across units.
  sym_map <- do.call(rbind, lapply(units, function(u) {
    tt <- fits_by_unit[[u]]$top[[contrast]]
    if ("symbol" %in% names(tt)) unique(tt[, c("gene", "symbol")]) else NULL
  }))
  sym <- if (!is.null(sym_map)) sym_map$symbol[match(all_sig, sym_map$gene)] else NA_character_

  data.frame(regime = regime, estimator = estimator, contrast = contrast,
             gene = all_sig, symbol = sym, units_sig = units_sig,
             n_units_sig = as.integer(rowSums(membership)),
             in_microglia = in_mg, class = cls,
             row.names = NULL, stringsAsFactors = FALSE)
}

# --- Sanity gate: microglia subset vs the canonical microglia DE ----------
#
# Assert the pooled-microglia unit from seurat_full reproduces the canonical
# microglia DE (de_snrnaseq / de_snrnaseq_nebula). `reference_top` is the
# per-contrast top list of the canonical fit (de_snrnaseq$fit$top for
# pseudobulk, de_snrnaseq_nebula$top for NEBULA). The two are keyed in
# DIFFERENT spaces -- seurat_full is symbol-keyed, the canonical caches are
# ensembl-keyed with a `symbol` column -- so both sides are bridged to SYMBOL
# space (NA symbols dropped, duplicates collapsed by mean logFC) before the
# join. Returns per-contrast Spearman(logFC) on shared symbols; the build
# asserts rho > `min_rho`. Not expected to be 1.0 -- seurat_full microglia may
# differ slightly from microglia_seurat_processed, and many->one symbol
# collapse adds noise.
microglia_crosscheck <- function(mg_fit, reference_top, estimator,
                                 contrasts = specificity_contrasts,
                                 min_rho = 0.95) {
  to_symbol_space <- function(tbl) {
    sym <- tbl$symbol %||% tbl$gene
    d <- data.frame(symbol = as.character(sym), logFC = tbl$logFC,
                    stringsAsFactors = FALSE)
    d <- d[!is.na(d$symbol) & d$symbol != "" & is.finite(d$logFC), , drop = FALSE]
    dplyr::summarise(dplyr::group_by(d, symbol),
                     logFC = mean(logFC), .groups = "drop")
  }
  purrr::map_dfr(contrasts, function(cn) {
    a <- to_symbol_space(mg_fit$top[[cn]])    |> dplyr::rename(logFC_a = logFC)
    b <- to_symbol_space(reference_top[[cn]]) |> dplyr::rename(logFC_b = logFC)
    cmp <- dplyr::inner_join(a, b, by = "symbol")
    rho <- if (nrow(cmp) > 2)
      suppressWarnings(cor(cmp$logFC_a, cmp$logFC_b, method = "spearman"))
      else NA_real_
    data.frame(estimator = estimator, contrast = cn, n_shared = nrow(cmp),
               spearman = rho, passes = !is.na(rho) & rho > min_rho,
               stringsAsFactors = FALSE)
  })
}

# --- Orchestrator: assemble R1-R4 from the full fit grid ------------------
#
# `fits` is the nested grid produced by scripts/build_celltype_specificity.R:
#   fits$<regime>$<estimator>$<unit> = wrapper output
# with regime in {native, matched} and estimator in {pseudobulk, nebula}.
# Returns a named list of tidy tables:
#   tally                   -- R1+R2 stacked (split into _native/_matched TSVs by the caller)
#   interaction_concordance -- R3a (vs_microglia, per estimator) + R3c (pb_vs_nebula) stacked
#   specificity_class       -- R4 per-gene classes (interaction, both regimes x estimators)
assemble_specificity_tables <- function(fits,
                                        contrasts = specificity_contrasts,
                                        class_contrast = "interaction",
                                        class_fdr = 0.10) {
  regimes    <- names(fits)
  estimators <- c("pseudobulk", "nebula")

  tally <- purrr::map_dfr(regimes, function(rg) {
    purrr::map_dfr(estimators, function(es) {
      if (is.null(fits[[rg]][[es]])) return(NULL)
      specificity_tally(fits[[rg]][[es]], regime = rg, estimator = es,
                        contrasts = contrasts)
    })
  })

  conc_ref <- purrr::map_dfr(regimes, function(rg) {
    purrr::map_dfr(estimators, function(es) {
      if (is.null(fits[[rg]][[es]])) return(NULL)
      interaction_concordance(fits[[rg]][[es]], estimator = es, regime = rg)
    })
  })
  conc_cross <- purrr::map_dfr(regimes, function(rg) {
    if (is.null(fits[[rg]]$pseudobulk) || is.null(fits[[rg]]$nebula)) return(NULL)
    cross_estimator_concordance(fits[[rg]]$pseudobulk, fits[[rg]]$nebula, regime = rg)
  })

  spec_cls <- purrr::map_dfr(regimes, function(rg) {
    purrr::map_dfr(estimators, function(es) {
      if (is.null(fits[[rg]][[es]])) return(NULL)
      specificity_class_table(fits[[rg]][[es]], regime = rg, estimator = es,
                              contrast = class_contrast, fdr = class_fdr)
    })
  })

  list(tally = tally,
       interaction_concordance = dplyr::bind_rows(conc_ref, conc_cross),
       specificity_class = spec_cls)
}

# --- R5: pathway-level interaction enrichment per unit --------------------
#
# The pathway-level analogue of R1/R2: instead of asking "how many GENES move
# at the interaction in each unit", ask "how many PATHWAYS are enriched at the
# interaction in each unit, and are the microglial ones cell-type-restricted?"
# Reuses R/fgsea.R's run_fgsea_for_contrast on each unit's interaction ranking
# (the NEBULA z-stat / limma t, keyed by symbol; logFC fallback if no `t`).
# Returns a long (regime x estimator x unit x pathway) frame -- the build
# filters it for the TSVs and the cache stores the microglia-significant
# slice for rmd/23. Set a seed before calling (fgseaMultilevel samples).
specificity_pathway_fgsea <- function(fits_by_unit, gene_sets, regime,
                                      estimator, contrast = "interaction",
                                      min_size = 15, max_size = 500) {
  purrr::imap_dfr(fits_by_unit, function(fo, unit) {
    tt <- fo$top[[contrast]]
    if (!"t" %in% names(tt)) tt$t <- tt$logFC      # logFC fallback (mirrors run_fgsea_per_dataset)
    fg <- run_fgsea_for_contrast(tt, gene_sets, stat_col = "t",
                                 gene_col = "symbol",
                                 min_size = min_size, max_size = max_size) |>
      as.data.frame()
    if (!nrow(fg)) return(NULL)
    data.frame(regime = regime, estimator = estimator, unit = unit,
               contrast = contrast, pathway = fg$pathway, NES = fg$NES,
               pval = fg$pval, padj = fg$padj, size = fg$size,
               stringsAsFactors = FALSE)
  })
}

# R5 tally: per (regime x estimator x unit) count of pathways enriched at the
# interaction at padj<0.05 and <0.10, split by NES sign (up = NES>0). The
# pathway-level headline; the verdict reads the MATCHED rows. `pathway_long`
# is one or more specificity_pathway_fgsea() frames row-bound.
specificity_pathway_tally <- function(pathway_long) {
  if (!nrow(pathway_long)) return(pathway_long)
  pathway_long |>
    dplyr::group_by(regime, estimator, unit, contrast) |>
    dplyr::summarise(
      n_path_tested = dplyr::n(),
      n_path_05  = sum(padj < 0.05, na.rm = TRUE),
      n_path_up_05   = sum(padj < 0.05 & NES > 0, na.rm = TRUE),
      n_path_down_05 = sum(padj < 0.05 & NES < 0, na.rm = TRUE),
      n_path_10  = sum(padj < 0.10, na.rm = TRUE),
      n_path_up_10   = sum(padj < 0.10 & NES > 0, na.rm = TRUE),
      n_path_down_10 = sum(padj < 0.10 & NES < 0, na.rm = TRUE),
      .groups = "drop")
}

# R5 cross-unit view: for every pathway enriched in MICROGLIA at the
# interaction (padj<fdr), report its NES / padj in ALL units, per
# regime x estimator. Directly answers "is the microglial interaction
# enrichment cell-type-restricted?" -- a pathway with n_units_sig==1 is
# microglia-restricted; high n_units_sig means a shared programme.
# `n_units_sig` counts units (incl microglia) at padj<fdr for that pathway.
microglia_pathway_cross_unit <- function(pathway_long, fdr = 0.10,
                                         microglia_unit = "Microglia") {
  if (!nrow(pathway_long)) return(pathway_long[0, ])
  purrr::map_dfr(split(pathway_long, ~ regime + estimator, drop = TRUE), function(d) {
    mg <- d[d$unit == microglia_unit & is.finite(d$padj) & d$padj < fdr, , drop = FALSE]
    if (!nrow(mg)) return(NULL)
    sel <- d[d$pathway %in% mg$pathway, , drop = FALSE]
    n_sig <- tapply(sel$padj < fdr, sel$pathway, sum, na.rm = TRUE)
    sel$n_units_sig <- as.integer(n_sig[sel$pathway])
    sel$mg_restricted <- sel$n_units_sig == 1L
    sel
  })
}

# --- Plots (display-only; assembled in rmd/23) ----------------------------

# R1/R2 tally heatmap: unit x contrast tiles coloured by significant-gene
# count at the chosen FDR, faceted by estimator, one regime per call. Units
# default to alphabetical order (no privileged ordering -- guardrail).
plot_specificity_tally_heatmap <- function(tally, regime,
                                           fdr = c("05", "10"),
                                           contrasts = specificity_contrasts) {
  fdr  <- match.arg(fdr)
  col  <- paste0("n_sig_", fdr)
  d <- tally |> dplyr::filter(.data$regime == !!regime)
  d$value    <- d[[col]]
  d$unit     <- factor(d$unit, levels = sort(unique(d$unit)))
  d$contrast <- factor(d$contrast, levels = contrasts)
  ggplot(d, aes(unit, contrast, fill = value)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = value), size = 3) +
    scale_fill_gradient(low = "grey92", high = "#b40426", name = "n sig") +
    facet_wrap(~ estimator) +
    labs(title = sprintf("Significant genes (FDR<0.%s) -- %s regime", fdr, regime),
         x = NULL, y = NULL) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# R3b curated interaction-logFC heatmap: rows = a headline gene set (symbols),
# cols = the 6 units, fill = interaction log2FC, for one regime x estimator.
# `genes` is a vector of gene SYMBOLS (the DAM programme + M-002 progression
# drivers + H-phase interaction genes are assembled by the caller). Genes
# absent from a unit's table show as NA tiles.
plot_interaction_logfc_heatmap <- function(fits_by_unit, genes, regime,
                                           estimator, contrast = "interaction") {
  units <- names(fits_by_unit)
  long <- purrr::map_dfr(units, function(u) {
    tt <- fits_by_unit[[u]]$top[[contrast]]
    tt <- tt[!is.na(tt$symbol) & tt$symbol %in% genes, c("symbol", "logFC")]
    if (!nrow(tt)) return(NULL)
    tt <- dplyr::group_by(tt, symbol) |> dplyr::summarise(logFC = logFC[1], .groups = "drop")
    data.frame(unit = u, symbol = tt$symbol, logFC = tt$logFC,
               stringsAsFactors = FALSE)
  })
  long$unit   <- factor(long$unit, levels = sort(unique(long$unit)))
  long$symbol <- factor(long$symbol, levels = rev(intersect(genes, unique(long$symbol))))
  lim <- max(abs(long$logFC), na.rm = TRUE)
  ggplot(long, aes(unit, symbol, fill = logFC)) +
    geom_tile(colour = "white") +
    scale_fill_gradient2(low = "#3b4cc0", mid = "white", high = "#b40426",
                         midpoint = 0, limits = c(-lim, lim),
                         name = "log2FC") +
    labs(title = sprintf("Interaction log2FC over headline genes -- %s / %s",
                         regime, estimator), x = NULL, y = NULL) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# R3c cross-estimator scatter: pseudobulk vs NEBULA `contrast` log2FC, one
# facet per unit (the 02b power-gain panel over all 6 units), one regime.
plot_cross_estimator_panel <- function(pb_by_unit, neb_by_unit, regime,
                                       contrast = "interaction",
                                       units = NULL) {
  if (is.null(units)) units <- sort(intersect(names(pb_by_unit), names(neb_by_unit)))
  long <- purrr::map_dfr(units, function(u) {
    pb <- pb_by_unit[[u]]$top[[contrast]]  |> dplyr::transmute(gene, logFC_pb = logFC)
    nb <- neb_by_unit[[u]]$top[[contrast]] |> dplyr::transmute(gene, logFC_nb = logFC)
    cmp <- dplyr::inner_join(pb, nb, by = "gene") |>
      dplyr::filter(is.finite(logFC_pb), is.finite(logFC_nb))
    if (!nrow(cmp)) return(NULL)
    cmp$unit <- u
    cmp
  })
  long$unit <- factor(long$unit, levels = units)
  ggplot(long, aes(logFC_pb, logFC_nb)) +
    geom_hline(yintercept = 0, colour = "grey80") +
    geom_vline(xintercept = 0, colour = "grey80") +
    geom_abline(slope = 1, intercept = 0, colour = "#b40426", linetype = "dashed") +
    geom_point(alpha = 0.2, size = 0.4) +
    facet_wrap(~ unit, scales = "free") +
    labs(title = sprintf("pseudobulk vs NEBULA %s log2FC -- %s regime",
                         contrast, regime),
         x = "pseudobulk log2FC", y = "NEBULA log2FC") +
    theme_bw()
}

# R4 specificity-class bar: gene counts per class (microglia_unique / shared /
# non_microglial_unique), faceted by estimator, one regime per call.
plot_specificity_class <- function(spec_class, regime) {
  d <- spec_class |> dplyr::filter(.data$regime == !!regime)
  if (!nrow(d)) return(ggplot() + theme_void() +
    labs(title = sprintf("No interaction-significant genes (%s regime)", regime)))
  lv <- c("microglia_unique", "shared", "non_microglial_unique")
  d$class <- factor(d$class, levels = lv)
  counts <- dplyr::count(d, estimator, class, .drop = FALSE)
  ggplot(counts, aes(class, n, fill = class)) +
    geom_col() +
    geom_text(aes(label = n), vjust = -0.3, size = 3) +
    scale_fill_manual(values = c(microglia_unique = "#b40426",
                                 shared = "#9e9ac8",
                                 non_microglial_unique = "#3b4cc0"),
                      guide = "none") +
    facet_wrap(~ estimator) +
    labs(title = sprintf("Interaction-gene specificity class -- %s regime", regime),
         x = NULL, y = "genes") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 20, hjust = 1))
}
