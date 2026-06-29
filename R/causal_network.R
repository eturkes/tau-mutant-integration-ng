# Causal signalling-network reconstruction helpers (narrative arc J;
# plan storage/notes/causal_network_plan.md, sessions J2..J5). The
# mechanism layer (E/F/G) produced endpoint activities -- kinases
# (Phase F), TFs (Phase E), LR pairs (Phase G) -- but no causal wiring
# between them. This module assembles the inputs CARNIVAL needs to
# reconstruct, per contrast, the directed signed signalling subnetwork
# that connects the kinase/phospho UPSTREAM layer to the
# TF/transcription DOWNSTREAM layer over the OmniPath prior-knowledge
# network (PKN), solved by ILP.
#
# Design (J1 gate, locked 2026-06-04):
#   1. Standard per-contrast CARNIVAL: kinases seed `perturbations`
#      (signed inputObj), TF activities are `measurements` (signed
#      measObj). 5 canonical contrasts. (NOT inverse, NOT cosmosR.)
#   2. TF complexes are resolved by reading the split-complexes cache
#      (tf_activity_decoupler_split.rds, get_collectri(split_complexes=
#      TRUE)); sources are then single-protein mouse symbols, so the
#      PKN mapping is the identity and NF-kB resolves to Nfkb1/Rela.
#   3. The PKN is pruned to microglia-expressed symbols (the microglial
#      readout drives the question). Bulk-phospho kinase seeds therefore
#      enter a microglia-filtered network -- a cross-compartment bridge
#      that MUST be caveated wherever kinase-seeded paths are read.
#
# CARNIVAL I/O formats (verified against the package toy data, 2026-06-04):
#   * priorKnowledgeNetwork: data.frame(source<chr>, interaction<int in
#     {-1,+1}>, target<chr>). Column name is `interaction`, NOT `sign`.
#   * perturbations (inputObj): named numeric vector, values are the
#     perturbation SIGN in {-1,+1}; names are node ids (mouse symbols).
#   * measurements (measObj): named numeric vector, signed CONTINUOUS
#     activity; names are node ids.
#   * result: list(weightedSIF = df(Node1, Sign, Node2, Weight),
#     nodesAttributes = df(Node, NodeType[P/M/inferred], ZeroAct, UpAct,
#     DownAct, AvgAct), sifAll, attributesAll). Weight and AvgAct are
#     the consensus across all optimal solutions (0..100 / -100..100),
#     so reading them satisfies the multi-solution-averaging guardrail.
#
# Cache shapes consumed (both: modality -> contrast -> tibble(statistic,
# source, score, p_value); see causal_network_plan.md "Grounded facts"):
#   * tf_activity_decoupler_split.rds  TF measObj source (single-protein
#                                      mouse symbols after split).
#   * kinase_activity_decoupler.rds    kinase inputObj source (mouse
#                                      symbols; phospho_raw /
#                                      phospho_corrected modalities).
# Convention (matches R/tf_inference.R): sign/magnitude from the
# `consensus` statistic; significance from BH-adjusted `ulm` p_value,
# threshold 0.10, computed within modality x contrast. The sign of the
# activity score IS the direction (no separate sign column exists).
#
# Reuses R/tf_inference.R (sourced earlier in R/helpers.R):
#   .extract_tf_per_modality()  per-modality score/padj matrices.
#   build_axis_gene_universe()  the 3 D2 axis gene universes (X.2).
# This module is otherwise self-contained.

# --------------------------------------------------------------------
# PKN assembly
# --------------------------------------------------------------------

# Assemble the OmniPath signed-directed prior-knowledge network in the
# CARNIVAL `(source, interaction, target)` format, in mouse gene-symbol
# space.
#
# An edge is kept iff it is directed AND carries exactly one of
# stimulation / inhibition (so the sign is unambiguous); the sign is
# +1 for stimulation, -1 for inhibition. Self-loops are dropped (ILP
# does not use them). A (source, target) pair that survives with BOTH
# signs is sign-ambiguous in a single-sign PKN and is dropped wholesale
# (count recorded in the `n_conflicts_dropped` attribute).
#
# Arguments:
#   organism      NCBI taxon id passed to OmnipathR. Default 10090
#                 (mouse); the J1 feasibility probe confirmed the mouse
#                 PKN (~18k signed-directed edges) covers the endpoints,
#                 so no human-ortholog fallback is needed.
#   expr_filter   optional character vector of gene symbols. When
#                 supplied, only edges whose BOTH endpoints are in the
#                 set are kept (the microglia-expressed prune, J1 dec.3).
#                 NULL leaves the PKN unfiltered.
#   drop_self_loops  logical; drop source == target edges. Default TRUE.
#
# Returns a tibble(source<chr>, interaction<int>, target<chr>), distinct
# rows, ready to pass straight to `runVanillaCarnival(priorKnowledge
# Network = ...)`. Attributes: `n_conflicts_dropped` (sign-ambiguous
# pairs removed), `n_edges_raw` (rows before filtering).
build_omnipath_pkn <- function(organism        = 10090,
                               expr_filter      = NULL,
                               drop_self_loops  = TRUE) {
  # `import_omnipath_interactions` was renamed to `omnipath_interactions`
  # in OmnipathR 3.x (the former is now deprecated); prefer the new name
  # when present so the build log stays clean, falling back otherwise.
  has_new <- "omnipath_interactions" %in% getNamespaceExports("OmnipathR")
  ia <- if (has_new) OmnipathR::omnipath_interactions(organism = organism)
        else         OmnipathR::import_omnipath_interactions(organism = organism)
  ia <- tibble::as_tibble(ia)
  req <- c("is_directed", "is_stimulation", "is_inhibition",
           "source_genesymbol", "target_genesymbol")
  missing <- setdiff(req, names(ia))
  if (length(missing) > 0L) {
    stop(sprintf("build_omnipath_pkn: OmniPath output missing columns: %s",
                 paste(missing, collapse = ", ")), call. = FALSE)
  }
  n_edges_raw <- nrow(ia)

  # Coerce the directionality flags to logical (OmnipathR may return
  # logical or 0/1 integer); NA flags are treated as FALSE.
  as_lgl0 <- function(x) { v <- as.logical(x); v[is.na(v)] <- FALSE; v }
  is_dir  <- as_lgl0(ia$is_directed)
  is_stim <- as_lgl0(ia$is_stimulation)
  is_inh  <- as_lgl0(ia$is_inhibition)
  keep    <- is_dir & xor(is_stim, is_inh)

  pkn <- tibble::tibble(
    source      = as.character(ia$source_genesymbol)[keep],
    interaction = ifelse(is_stim[keep], 1L, -1L),
    target      = as.character(ia$target_genesymbol)[keep]
  )

  ok <- !is.na(pkn$source) & nzchar(pkn$source) &
        !is.na(pkn$target) & nzchar(pkn$target)
  pkn <- pkn[ok, , drop = FALSE]
  if (isTRUE(drop_self_loops)) {
    pkn <- pkn[pkn$source != pkn$target, , drop = FALSE]
  }
  pkn <- dplyr::distinct(pkn, source, interaction, target)

  # A (source, target) carrying both +1 and -1 is sign-ambiguous; drop
  # both directions rather than arbitrarily choosing a sign.
  conflict <- pkn |>
    dplyr::group_by(source, target) |>
    dplyr::summarise(n_sign = dplyr::n_distinct(interaction),
                     .groups = "drop") |>
    dplyr::filter(n_sign > 1L)
  if (nrow(conflict) > 0L) {
    pkn <- dplyr::anti_join(pkn, conflict, by = c("source", "target"))
  }

  if (!is.null(expr_filter)) {
    uni <- unique(as.character(expr_filter))
    pkn <- pkn[pkn$source %in% uni & pkn$target %in% uni, , drop = FALSE]
  }

  pkn <- tibble::as_tibble(pkn)
  attr(pkn, "n_conflicts_dropped") <- nrow(conflict)
  attr(pkn, "n_edges_raw")         <- n_edges_raw
  pkn
}

# Derive the microglia-expressed gene-symbol universe for the PKN prune
# (J1 dec.3). The snRNAseq objects are keyed by Ensembl gene id, while
# the PKN and the activity caches are in mouse-symbol space, so this
# maps Ensembl -> symbol through the project symbol_map.
#
# Arguments:
#   seurat       a Seurat object (the microglia subset). Genes detected
#                (count > 0) in at least `min_cells` cells of the chosen
#                assay are called expressed.
#   symbol_map   data.frame with `ensembl` and `symbol` columns
#                (snrnaseq_symbol_map.rds).
#   assay        assay to read counts from. Default "RNA".
#   min_cells    minimum number of cells with a non-zero count for a
#                gene to count as expressed. Default 10. (Absolute count
#                rather than a fraction so the threshold is independent
#                of the subset size; a future session may switch to a
#                fraction by passing `min_cells = ceiling(frac * ncol)`.)
#
# Returns a sorted character vector of unique mouse gene symbols.
microglia_expressed_symbols <- function(seurat, symbol_map,
                                         assay     = "RNA",
                                         min_cells = 10L) {
  stopifnot(is.data.frame(symbol_map),
            all(c("ensembl", "symbol") %in% names(symbol_map)))
  counts <- tryCatch(
    SeuratObject::GetAssayData(seurat, assay = assay, layer = "counts"),
    error = function(e)
      SeuratObject::GetAssayData(seurat, assay = assay, slot = "counts")
  )
  n_cells_detected <- Matrix::rowSums(counts > 0)
  expressed_ens <- rownames(counts)[n_cells_detected >= as.integer(min_cells)]
  syms <- symbol_map$symbol[match(expressed_ens, symbol_map$ensembl)]
  syms <- syms[!is.na(syms) & nzchar(syms)]
  # The symbol_map falls back to the Ensembl id as the "symbol" for genes
  # with no assigned symbol; drop those so the returned universe is true
  # gene-symbol space (they would never match the symbol-keyed PKN anyway).
  syms <- syms[!grepl("^ENSMUSG", syms)]
  sort(unique(as.character(syms)))
}

# Per-contrast reachability prune of the PKN (J3 tractability fix,
# user-approved 2026-06-04 -- a deliberate, documented deviation from
# the J1-locked "solve over the full microglia-expressed PKN" design).
# The full microglia-pruned PKN (~10.7k edges / ~3.8k nodes) yields a
# CARNIVAL ILP that neither lpSolve nor cbc could solve to optimality in
# a workable wall-clock (lpSolve >6 min on one contrast; cbc reported no
# feasible solution within a 600 s limit). Restricting the PKN, per
# contrast, to the nodes that lie on a directed seed->target walk of
# length <= `max_path` collapses the ILP to a few-hundred-edge problem
# cbc proves optimal in seconds (nlgf_in_p301s @ L=3: 807 edges/153
# nodes; cbc Optimal in ~3.6 s).
#
# Keep a node iff fwd_dist(seed -> node) + bwd_dist(node -> target) <=
# max_path, computed by multi-source BFS (forward from the seeds over
# source->target, backward from the targets over the reversed edges).
# This is SOLUTION-PRESERVING up to depth max_path: every node on any
# seed->target path of length <= max_path satisfies the inequality, so
# the induced subgraph retains every such path intact (a node at
# forward-position d1 on a length-L path has fwd <= d1 and bwd <= L-d1,
# hence fwd+bwd <= L). Paths LONGER than max_path hops are excluded by
# construction -- `max_path` is a tunable mechanistic horizon, set to 3
# to match the 1-3 hop kinase->TF wiring J1 anticipated; it is recorded
# in the cache manifest and any read of the networks MUST caveat that
# deeper indirect routes are out of scope.
#
# Arguments:
#   pkn       PKN tibble(source, interaction, target) from
#             build_omnipath_pkn().
#   seeds     character vector of perturbation node ids (kinase symbols,
#             names(inputObj)); intersected with PKN nodes internally.
#   targets   character vector of measurement node ids (TF symbols,
#             names(measObj)); intersected with PKN nodes internally.
#   max_path  maximum seed->target walk length (hops) a node may lie on
#             to be kept. Default 3L.
#
# Returns the induced sub-PKN tibble (same columns). When either side
# has no node in the PKN no seed->target path can exist, so a 0-row PKN
# is returned (the solve then records an honest empty network). Carries
# attributes `max_path`, `n_seeds_in`, `n_seeds_kept`, `n_targets_in`,
# `n_targets_kept` (kept = endpoint survives the prune).
restrict_pkn_to_reachable <- function(pkn, seeds, targets, max_path = 3L) {
  stopifnot(is.data.frame(pkn),
            all(c("source", "interaction", "target") %in% names(pkn)))
  max_path <- as.integer(max_path)
  nodes_in <- unique(c(pkn$source, pkn$target))
  seeds    <- intersect(unique(as.character(seeds)),   nodes_in)
  targets  <- intersect(unique(as.character(targets)), nodes_in)

  annotate <- function(p, kept_nodes) {
    attr(p, "max_path")       <- max_path
    attr(p, "n_seeds_in")     <- length(seeds)
    attr(p, "n_targets_in")   <- length(targets)
    attr(p, "n_seeds_kept")   <- length(intersect(seeds,   kept_nodes))
    attr(p, "n_targets_kept") <- length(intersect(targets, kept_nodes))
    p
  }
  if (length(seeds) == 0L || length(targets) == 0L) {
    return(annotate(pkn[0L, , drop = FALSE], character(0)))
  }

  # multi-source BFS hop-distance over a directed edge list (from -> to)
  bfs_dist <- function(from, to, sources, nodes) {
    d <- stats::setNames(rep(Inf, length(nodes)), nodes)
    d[sources] <- 0
    adj <- split(to, from)          # from -> vector of to
    frontier <- sources; k <- 0L
    while (length(frontier) > 0L) {
      k  <- k + 1L
      nb <- unique(unlist(adj[frontier], use.names = FALSE))
      nb <- nb[!is.na(nb)]
      nb <- nb[!is.finite(d[nb])]   # unvisited only
      if (length(nb) == 0L) break
      d[nb] <- k
      frontier <- nb
    }
    d
  }
  fwd <- bfs_dist(pkn$source, pkn$target, seeds,   nodes_in)  # from seeds
  bwd <- bfs_dist(pkn$target, pkn$source, targets, nodes_in)  # to targets

  on_path <- is.finite(fwd[nodes_in]) & is.finite(bwd[nodes_in]) &
             (fwd[nodes_in] + bwd[nodes_in] <= max_path)
  keep    <- nodes_in[on_path]
  pruned  <- pkn[pkn$source %in% keep & pkn$target %in% keep, , drop = FALSE]
  annotate(pruned, keep)
}

# --------------------------------------------------------------------
# CARNIVAL input formatters (from the activity caches)
# --------------------------------------------------------------------

# Build the signed TF measurement vector (measObj) for one contrast from
# the split-complexes TF activity cache. A TF enters the measurements
# iff it is significant (BH-adjusted `ulm` p_value < fdr) in at least
# one of `modalities`; its value is the mean `consensus` score across
# the modalities where it is significant (signed, continuous -- the
# CARNIVAL measObj convention).
#
# Arguments:
#   tf_cache          the loaded TF activity list. MUST be the
#                     split-complexes variant (tf_activity_decoupler_
#                     split.rds) so `source` ids are single-protein
#                     mouse symbols mapping directly onto the PKN.
#   contrast          contrast name (one of the 5 canonical contrasts).
#   modalities        modalities to draw measurements from. Default the
#                     three primary TF modalities (whole-microglia
#                     snRNAseq + GeoMx + proteomics) -- the same set the
#                     E3/E4 leader board ranks on, so the verdict TFs
#                     are recoverable. A strictly-microglial measObj is
#                     obtained by passing modalities = "snrnaseq".
#   fdr               BH significance threshold. Default 0.10 (project
#                     convention for activity inference).
#   sig_statistic     statistic whose BH p drives significance. Default
#                     "ulm" (matches tf_inference.R).
#   score_statistic   statistic whose signed score is the measurement.
#                     Default "consensus".
#
# Returns a named numeric vector (names = TF symbol, values = signed
# mean consensus score), dropping any zero/non-finite value. Empty
# named numeric(0) when no TF is significant.
tf_meas_from_cache <- function(tf_cache, contrast,
                               modalities      = c("snrnaseq", "geomx",
                                                   "proteomics"),
                               fdr             = 0.10,
                               sig_statistic   = "ulm",
                               score_statistic = "consensus") {
  stopifnot(is.list(tf_cache), length(contrast) == 1L,
            all(modalities %in% names(tf_cache)))
  mats <- .extract_tf_per_modality(tf_cache, contrast, modalities,
                                   sig_statistic   = sig_statistic,
                                   score_statistic = score_statistic)
  score_mat <- mats$score_mat
  padj_mat  <- mats$padj_mat
  if (nrow(score_mat) == 0L) return(stats::setNames(numeric(0), character(0)))

  sig_mat <- !is.na(padj_mat) & padj_mat < fdr
  keep    <- rowSums(sig_mat) > 0L
  if (!any(keep)) return(stats::setNames(numeric(0), character(0)))

  idx  <- which(keep)
  vals <- vapply(idx, function(i) {
    s <- score_mat[i, sig_mat[i, ]]
    s <- s[is.finite(s)]
    if (length(s) == 0L) NA_real_ else mean(s)
  }, numeric(1))
  names(vals) <- rownames(score_mat)[idx]
  vals[is.finite(vals) & vals != 0]
}

# Build the signed kinase perturbation vector (inputObj) for one
# contrast from the kinase activity cache. A kinase enters the
# perturbations iff significant (BH-adjusted `ulm` p < fdr) in the
# chosen modality; its value is the SIGN of its `consensus` score
# (CARNIVAL inputObj convention: the perturbation direction in
# {-1,+1}).
#
# Arguments:
#   kin_cache         the loaded kinase activity list.
#   contrast          contrast name.
#   modality          which kinase modality to seed from. Default
#                     "phospho_corrected" (protein-corrected sites; the
#                     standing project choice). "phospho_raw" available
#                     for sensitivity.
#   fdr               BH significance threshold. Default 0.10.
#   sig_statistic     statistic whose BH p drives significance. Default
#                     "ulm".
#   score_statistic   statistic whose sign drives the perturbation
#                     direction. Default "consensus".
#   as_sign           if TRUE (default) values are sign(score) in
#                     {-1,+1}; if FALSE the continuous consensus score
#                     is returned (for sensitivity / weighted runs).
#
# Returns a named numeric vector (names = kinase symbol). Empty named
# numeric(0) when no kinase is significant.
kinase_input_from_cache <- function(kin_cache, contrast,
                                     modality        = "phospho_corrected",
                                     fdr             = 0.10,
                                     sig_statistic   = "ulm",
                                     score_statistic = "consensus",
                                     as_sign         = TRUE) {
  stopifnot(is.list(kin_cache), modality %in% names(kin_cache))
  if (is.null(kin_cache[[modality]][[contrast]])) {
    return(stats::setNames(numeric(0), character(0)))
  }
  tbl <- kin_cache[[modality]][[contrast]]
  sig <- tbl[tbl$statistic == sig_statistic, , drop = FALSE]
  if (nrow(sig) == 0L) return(stats::setNames(numeric(0), character(0)))
  sig$padj <- stats::p.adjust(sig$p_value, method = "BH")
  sco <- tbl[tbl$statistic == score_statistic,
             c("source", "score"), drop = FALSE]
  m <- merge(sig[, c("source", "padj")], sco, by = "source", all = FALSE)
  m <- m[is.finite(m$padj) & m$padj < fdr &
         is.finite(m$score) & m$score != 0, , drop = FALSE]
  if (nrow(m) == 0L) return(stats::setNames(numeric(0), character(0)))
  v <- if (isTRUE(as_sign)) sign(m$score) else m$score
  stats::setNames(as.numeric(v), as.character(m$source))
}

# --------------------------------------------------------------------
# CARNIVAL runner + result tidying
# --------------------------------------------------------------------

.cn_empty_edges <- function() {
  tibble::tibble(source = character(0), interaction = integer(0),
                 target = character(0), weight = numeric(0))
}
.cn_empty_nodes <- function() {
  tibble::tibble(node = character(0), node_type = character(0),
                 activity = numeric(0), up = numeric(0),
                 down = numeric(0), zero = numeric(0))
}

# Run standard per-contrast CARNIVAL and return a tidy result.
#
# Measurements and perturbations are first intersected with the PKN node
# universe (CARNIVAL ignores ids absent from the network); the dropped
# ids are recorded in `meta` because an endpoint stranded by the
# microglia-expressed prune is itself a finding. A node present in both
# inputObj and measObj is kept only as a perturbation (it cannot be both
# the seed and the readout). lpSolve `workdir`/`outputFolder` are
# redirected to a tempdir with `keepLPFiles = FALSE` (J1 housekeeping
# trap: the defaults pollute the project root).
#
# Arguments:
#   pkn               PKN tibble from build_omnipath_pkn().
#   measObj           named signed measurement vector (TF activities).
#   inputObj          named signed perturbation vector (kinases).
#                     Standard CARNIVAL requires a non-empty inputObj;
#                     when it is empty/NULL the function records an
#                     honest empty network (status "no_perturbations_in_pkn")
#                     rather than failing.
#   solver            "lpSolve" (default; the in-R free solver) /
#                     "cbc" / "cplex".
#   beta_weight       network-size penalty (CARNIVAL betaWeight). NULL
#                     keeps the solver default.
#   contrast          contrast label stored in `meta`.
#   fdr               FDR used to build the inputs, stored in `meta`.
#   carnival_options  optional pre-built options list (overrides solver/
#                     beta_weight/tempdir handling).
#
# Returns list(nodes, edges, meta):
#   edges  tibble(source, interaction, target, weight) -- recovered
#          edges (weightedSIF rows with weight > 0).
#   nodes  tibble(node, node_type[input/measured/inferred], activity
#          (AvgAct), up, down, zero) -- all nodesAttributes rows.
#   meta   list(contrast, status, solver, beta_weight, fdr,
#          n_pkn_nodes, n_pkn_edges, n_meas_in, n_meas_dropped,
#          meas_dropped, n_pert_in, n_pert_dropped, pert_dropped,
#          n_solutions, n_edges_zero_dropped).
run_carnival_for_contrast <- function(pkn, measObj, inputObj = NULL,
                                      solver           = "lpSolve",
                                      beta_weight      = NULL,
                                      contrast         = NA_character_,
                                      fdr              = NA_real_,
                                      carnival_options = NULL) {
  stopifnot(is.data.frame(pkn),
            all(c("source", "interaction", "target") %in% names(pkn)))
  nodes_in_pkn <- unique(c(pkn$source, pkn$target))

  meas_use  <- measObj[names(measObj) %in% nodes_in_pkn]
  meas_drop <- setdiff(names(measObj), names(meas_use))
  if (is.null(inputObj)) inputObj <- stats::setNames(numeric(0), character(0))
  pert_use  <- inputObj[names(inputObj) %in% nodes_in_pkn]
  pert_drop <- setdiff(names(inputObj), names(pert_use))
  # A node cannot be both perturbed and measured; keep it as a seed.
  meas_use  <- meas_use[!names(meas_use) %in% names(pert_use)]

  base_meta <- list(
    contrast       = contrast,
    solver         = solver,
    beta_weight    = beta_weight,
    fdr            = fdr,
    n_pkn_nodes    = length(nodes_in_pkn),
    n_pkn_edges    = nrow(pkn),
    n_meas_in      = length(meas_use),
    n_meas_dropped = length(meas_drop),
    meas_dropped   = meas_drop,
    n_pert_in      = length(pert_use),
    n_pert_dropped = length(pert_drop),
    pert_dropped   = pert_drop,
    n_solutions    = 0L,
    n_edges_zero_dropped = 0L
  )

  empty_with <- function(status) {
    base_meta$status <- status
    list(nodes = .cn_empty_nodes(), edges = .cn_empty_edges(),
         meta = base_meta)
  }
  if (length(pert_use) == 0L) return(empty_with("no_perturbations_in_pkn"))
  if (length(meas_use) == 0L) return(empty_with("no_measurements_in_pkn"))

  if (is.null(carnival_options)) {
    wd <- tempfile("carnival_")
    dir.create(wd, showWarnings = FALSE)
    carnival_options <- switch(
      solver,
      lpSolve = CARNIVAL::defaultLpSolveCarnivalOptions(),
      cbc     = CARNIVAL::defaultCbcSolveCarnivalOptions(),
      cplex   = CARNIVAL::defaultCplexCarnivalOptions(),
      stop(sprintf("run_carnival_for_contrast: unsupported solver '%s'",
                   solver), call. = FALSE)
    )
    carnival_options$outputFolder <- wd
    carnival_options$workdir      <- wd
    carnival_options$keepLPFiles  <- FALSE
    if (!is.null(beta_weight)) carnival_options$betaWeight <- beta_weight
  }

  res <- tryCatch(
    suppressMessages(CARNIVAL::runVanillaCarnival(
      perturbations          = pert_use,
      measurements           = meas_use,
      priorKnowledgeNetwork  = as.data.frame(pkn),
      carnivalOptions        = carnival_options
    )),
    error = function(e) {
      warning(sprintf("run_carnival_for_contrast(%s): CARNIVAL error: %s",
                      contrast, conditionMessage(e)), call. = FALSE)
      NULL
    }
  )
  if (is.null(res) || is.null(res$weightedSIF)) {
    return(empty_with("solver_error_or_empty"))
  }

  sif <- as.data.frame(res$weightedSIF, stringsAsFactors = FALSE)
  edges_all <- tibble::tibble(
    source      = as.character(sif$Node1),
    interaction = as.integer(sif$Sign),
    target      = as.character(sif$Node2),
    weight      = as.numeric(sif$Weight)
  )
  n_zero <- sum(edges_all$weight == 0, na.rm = TRUE)
  edges  <- edges_all[is.finite(edges_all$weight) & edges_all$weight > 0,
                      , drop = FALSE]

  att <- as.data.frame(res$nodesAttributes, stringsAsFactors = FALSE)
  type_map <- c(P = "input", S = "input", M = "measured", T = "measured")
  raw_type <- as.character(att$NodeType)
  node_type <- unname(type_map[raw_type])
  node_type[is.na(node_type) | !nzchar(raw_type)] <- "inferred"
  nodes <- tibble::tibble(
    node      = as.character(att$Node),
    node_type = node_type,
    activity  = as.numeric(att$AvgAct),
    up        = as.numeric(att$UpAct),
    down      = as.numeric(att$DownAct),
    zero      = as.numeric(att$ZeroAct)
  )

  base_meta$status               <- "ok"
  base_meta$n_solutions          <- length(res$sifAll)
  base_meta$n_edges_zero_dropped <- as.integer(n_zero)
  base_meta$beta_weight          <- carnival_options$betaWeight %||% beta_weight
  list(nodes = nodes, edges = edges, meta = base_meta)
}

# --------------------------------------------------------------------
# Network summaries + axis restriction
# --------------------------------------------------------------------

# One-row summary of a tidy CARNIVAL network (counts active nodes and
# recovered edges, and the fraction of supplied measurements recovered
# as active). "Active" = AvgAct != 0; recovered edge = weight > 0
# (already enforced in run_carnival_for_contrast).
#
# Returns a one-row tibble: contrast, status, n_nodes_active,
# n_edges, n_input_active, n_measured_recovered, n_inferred_active,
# n_meas_in, frac_meas_recovered, solver, beta_weight.
summarise_network <- function(net) {
  stopifnot(is.list(net), all(c("nodes", "edges", "meta") %in% names(net)))
  nodes <- net$nodes
  meta  <- net$meta
  active <- if (nrow(nodes) == 0L) logical(0) else nodes$activity != 0
  n_meas_in <- meta$n_meas_in %||% NA_integer_
  n_meas_rec <- sum(active & nodes$node_type == "measured")
  tibble::tibble(
    contrast             = meta$contrast %||% NA_character_,
    status               = meta$status %||% NA_character_,
    n_nodes_active       = sum(active),
    n_edges              = nrow(net$edges),
    n_input_active       = sum(active & nodes$node_type == "input"),
    n_measured_recovered = n_meas_rec,
    n_inferred_active    = sum(active & nodes$node_type == "inferred"),
    n_meas_in            = n_meas_in,
    frac_meas_recovered  = if (is.na(n_meas_in) || n_meas_in == 0L) NA_real_
                           else n_meas_rec / n_meas_in,
    solver               = meta$solver %||% NA_character_,
    beta_weight          = meta$beta_weight %||% NA_real_
  )
}

# Restrict a tidy network to a node set (X.2 axis-restricted view).
# Pass an axis gene universe from build_axis_gene_universe() to see how
# much of the recovered network falls within a D2 axis.
#
# Arguments:
#   net        a tidy network list(nodes, edges, meta).
#   node_set   character vector of node symbols (e.g.
#              axis_universes$interaction_metabolic$universe).
#   mode       "induced" (default): keep nodes in node_set and edges
#              with BOTH endpoints in node_set. "incident": keep edges
#              with AT LEAST ONE endpoint in node_set, plus the nodes
#              they touch (the 1-hop neighbourhood -- use this to retain
#              signalling connectors that route a path through a node
#              outside the axis gene set).
#
# Returns a tidy network list with `meta$restricted_to` and
# `meta$restrict_mode` annotated.
restrict_network_to_axis <- function(net, node_set,
                                      mode = c("induced", "incident")) {
  mode <- match.arg(mode)
  stopifnot(is.list(net), all(c("nodes", "edges") %in% names(net)),
            is.character(node_set))
  set <- unique(as.character(node_set))
  edges <- net$edges
  if (nrow(edges) > 0L) {
    if (mode == "induced") {
      edges <- edges[edges$source %in% set & edges$target %in% set,
                     , drop = FALSE]
    } else {
      edges <- edges[edges$source %in% set | edges$target %in% set,
                     , drop = FALSE]
    }
  }
  keep_nodes <- if (mode == "induced") set
                else union(set, unique(c(edges$source, edges$target)))
  nodes <- net$nodes[net$nodes$node %in% keep_nodes, , drop = FALSE]

  meta <- net$meta %||% list()
  meta$restricted_to <- length(set)
  meta$restrict_mode <- mode
  list(nodes = tibble::as_tibble(nodes),
       edges = tibble::as_tibble(edges), meta = meta)
}

# ---------------------------------------------------------------------------
# Display helper (J4, rmd/19 -> section 20). Render one reconstructed
# signed-directed CARNIVAL subnetwork as a node-link diagram. Nodes are
# shaped by CARNIVAL role (kinase seed = triangle, intermediate = circle,
# measured TF = square) and filled by inferred activity sign (up = red,
# down = blue); edges are coloured by interaction sign (activation = green,
# inhibition = red) with an arrowhead for direction. Only ACTIVE nodes
# (activity != 0, i.e. those on a recovered path) are drawn. Returns NULL
# for an empty network so the caller can emit an honest "no network" note
# instead of a blank panel.
#
# The `stress` layout (graphlayouts, ggraph's robust default) tolerates the
# Gsk3b<->Mapk14 feedback cycles these networks contain without the
# DAG-breaking messages a layered layout would emit; arrowheads carry the
# directionality a layered layout would otherwise encode by vertical rank.
# `max.overlaps = Inf` forces every node label to render (ggrepel would
# otherwise silently drop labels on the busier ~30-node contrasts).
plot_causal_network <- function(net, title = NULL, layout = "stress",
                                seed = 1L) {
  stopifnot(is.list(net), all(c("nodes", "edges") %in% names(net)))
  edges <- net$edges
  if (is.null(edges) || nrow(edges) == 0L) {
    return(NULL)
  }
  used  <- union(edges$source, edges$target)
  nodes <- net$nodes[net$nodes$node %in% used, , drop = FALSE]

  node_df <- data.frame(
    name = nodes$node,
    role = factor(
      nodes$node_type,
      levels = c("input", "inferred", "measured"),
      labels = c("kinase (seed)", "intermediate", "TF (measured)")),
    dir = factor(
      ifelse(nodes$activity > 0, "up",
             ifelse(nodes$activity < 0, "down", "none")),
      levels = c("up", "down", "none")),
    stringsAsFactors = FALSE)

  edge_df <- data.frame(
    from = edges$source,
    to   = edges$target,
    sign = factor(ifelse(edges$interaction > 0, "activation", "inhibition"),
                  levels = c("activation", "inhibition")),
    stringsAsFactors = FALSE)

  g <- tidygraph::tbl_graph(nodes = node_df, edges = edge_df, directed = TRUE)

  set.seed(seed)  # `stress` layout is deterministic given a seed
  p <- ggraph::ggraph(g, layout = layout) +
    ggraph::geom_edge_link(
      ggplot2::aes(colour = sign),
      arrow   = grid::arrow(length = grid::unit(2.4, "mm"), type = "closed"),
      end_cap = ggraph::circle(4.5, "mm"),
      start_cap = ggraph::circle(4.5, "mm"),
      width = 0.6, alpha = 0.85) +
    ggraph::geom_node_point(
      ggplot2::aes(shape = role, fill = dir),
      size = 5, colour = "grey25", stroke = 0.5) +
    ggraph::geom_node_text(
      ggplot2::aes(label = name), repel = TRUE, size = 3,
      max.overlaps = Inf, seed = seed) +
    ggplot2::scale_shape_manual(
      values = c("kinase (seed)" = 24, "intermediate" = 21,
                 "TF (measured)" = 22),
      drop = FALSE, name = "node role") +
    ggplot2::scale_fill_manual(
      values = c(up = "#d73027", down = "#4575b4", none = "grey80"),
      drop = FALSE, name = "activity") +
    ggraph::scale_edge_colour_manual(
      values = c(activation = "#1a9850", inhibition = "#d73027"),
      drop = FALSE, name = "edge sign") +
    ggplot2::guides(
      fill = ggplot2::guide_legend(override.aes = list(shape = 21))) +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(legend.position = "right")

  if (!is.null(title)) {
    p <- p + ggplot2::ggtitle(title)
  }
  p
}
