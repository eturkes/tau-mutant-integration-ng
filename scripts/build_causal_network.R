#!/usr/bin/env Rscript
# Build the per-contrast CARNIVAL causal signalling networks for the
# causal-network layer (narrative arc J; plan
# storage/notes/causal_network_plan.md, session J3). The mechanism layer
# (E/F/G) produced endpoint activities -- kinases (Phase F), TFs (Phase
# E) -- but no causal wiring between them. This script reconstructs, per
# contrast, the directed signed signalling subnetwork connecting the
# kinase/phospho UPSTREAM layer to the TF/transcription DOWNSTREAM layer
# over the microglia-pruned OmniPath PKN, solved by ILP (CARNIVAL).
# Heavy compute lives here, OUTSIDE the knit (project convention in
# scripts/build_*.R); rmd/19 (J4) readRDS-loads the resulting cache.
#
# Design (J1 gate, locked 2026-06-04 -- see plan):
#   1. Standard per-contrast CARNIVAL: kinases seed `perturbations`
#      (signed inputObj), TF activities are `measurements` (signed
#      measObj). 5 canonical contrasts. (NOT inverse, NOT cosmosR.)
#   2. TF measObj read from the split-complexes cache
#      (tf_activity_decoupler_split.rds) so sources are single-protein
#      mouse symbols mapping directly onto the PKN (NF-kB -> Nfkb1/Rela).
#   3. PKN pruned to microglia-expressed symbols (snRNAseq microglia
#      subset). Bulk-phospho kinase seeds therefore enter a
#      microglia-filtered network -- a cross-compartment bridge caveated
#      wherever kinase-seeded paths are read (J4 X.3).
#   4. Per-contrast reachability prune at L=max_path (J3 deviation,
#      user-approved 2026-06-04). The full microglia-expressed PKN
#      (~10.7k edges) is ILP-intractable: lpSolve ran >6 min on one
#      contrast; cbc found no feasible solution in 600 s. Restricting
#      the PKN per contrast to nodes on a directed seed->target walk of
#      <= max_path hops (restrict_pkn_to_reachable) makes cbc prove
#      optimality in seconds. Solution-preserving up to depth max_path;
#      longer indirect routes are out of scope (recorded in params).
#
# Thresholds STATED upfront (plan anti-anchoring guardrail), recorded in
# the cache `params`/`pkn_meta` for the manifest:
#   solver       cbc (COIN-OR ILP; PuLP-bundled binary auto-located in
#                     the project .venv, override with --solver-path)
#   betaWeight   0.2  (CARNIVAL network-size sparsity penalty)
#   FDR          0.10 (project activity-inference convention; BH within
#                      modality x contrast on the `ulm` p_value)
#   TF measObj   modalities snrnaseq + geomx + proteomics (the E3/E4
#                leader-board set; sign/magnitude from `consensus`)
#   kinase input phospho_corrected (standing project choice), sign only
#   PKN filter   microglia-expressed (>= min_cells with non-zero count
#                in the PROCESSED microglia subset)
#   max_path     3 (reachability horizon; node kept iff on a <= max_path
#                   hop seed->target walk -- see design pt.4)
#
# Idempotent: skips the cache write if the output already exists unless
# `--overwrite` is passed. A subset run (`--contrasts a,b`) prints
# results + timings but NEVER writes/overwrites the cache (it is a
# probe); the cache is only written when all 5 contrasts ran in one
# invocation.
#
# CLI:
#   --overwrite          rebuild even if the cache exists
#   --contrasts a,b,c    restrict to a subset (probe / re-run; no write)
#   --solver cbc|lpSolve solver (default cbc)
#   --solver-path PATH   cbc binary (default: auto-located in .venv)
#   --beta W             override betaWeight (default 0.2)
#   --max-path N         reachability horizon in hops (default 3)
#   --threads N          cbc threads (default 1; >1 risks a cbc segfault)
#   --timelimit S        cbc CPU-second limit (default 600)
#
# Inputs (rds, under storage/cache/):
#   tf_activity_decoupler_split.rds   TF measObj source (split complexes)
#   kinase_activity_decoupler.rds     kinase inputObj source
#   microglia_seurat_processed.rds    expr filter (PROCESSED subset)
#   snrnaseq_symbol_map.rds           Ensembl -> symbol for the filter
#
# Output (storage/cache/causal_network.rds): a named list ->
#   $<contrast>  list(nodes, edges, meta) per the run_carnival_for_contrast
#                tidy schema (5 entries; honestly-empty nets keep a typed
#                meta$status rather than being dropped)
#   $summary     per-contrast one-row summaries (summarise_network) rbind
#   $pkn_meta    list(organism, n_edges, n_nodes, n_conflicts_dropped,
#                     n_edges_raw, n_expr_symbols, min_cells)
#   $params      list(solver, solver_path, beta_weight, max_path,
#                     threads, timelimit, fdr, tf_modalities,
#                     kinase_modality, contrasts, built_from)
#
# Runtime: with the reachability prune, cbc proves each contrast optimal
# in seconds (the unpruned PKN was ILP-intractable -- design pt.4). The
# script logs per-contrast wall-clock; a solve approaching --timelimit
# (600 s) signals the pruned ILP grew unexpectedly and the network may
# be a non-optimal incumbent, so watch that figure when reading results.

suppressPackageStartupMessages({
  library(dplyr); library(tibble)
})
setwd("/home/rstudio/tau-mutant-integration-ng")  # anchor relative paths

# ---- CLI --------------------------------------------------------------
args        <- commandArgs(trailingOnly = TRUE)
overwrite   <- "--overwrite" %in% args
arg_val <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[[i + 1L]]
}
CONTRASTS_ALL <- c("nlgf_in_maptki", "nlgf_in_p301s", "interaction",
                   "tau_alone", "tau_in_nlgf")
sub_arg      <- arg_val("--contrasts", NULL)
run_contrasts <- if (is.null(sub_arg)) {
  CONTRASTS_ALL
} else {
  trimws(strsplit(sub_arg, ",", fixed = TRUE)[[1]])
}
bad <- setdiff(run_contrasts, CONTRASTS_ALL)
if (length(bad) > 0L) {
  stop(sprintf("unknown contrast(s): %s (valid: %s)",
               paste(bad, collapse = ","),
               paste(CONTRASTS_ALL, collapse = ",")), call. = FALSE)
}
is_full_run <- setequal(run_contrasts, CONTRASTS_ALL)

SOLVER      <- arg_val("--solver", "cbc")
SOLVER_PATH <- arg_val("--solver-path", NULL)
BETA_WEIGHT <- as.numeric(arg_val("--beta", "0.2"))
MAX_PATH    <- as.integer(arg_val("--max-path", "3"))   # reachability horizon
# threads=1 by default: cbc 2.10.3 (the PuLP-bundled 2019 build) segfaults
# mid branch-and-bound in multi-threaded mode (observed at threads=6 on
# the L=3 nlgf_in_p301s ILP). Single-threaded is stable AND deterministic
# -- the pruned problems are tiny (<= ~800 edges), so it stays fast.
THREADS     <- as.integer(arg_val("--threads", "1"))    # cbc only
TIMELIMIT   <- as.numeric(arg_val("--timelimit", "600")) # cbc only (s)
FDR         <- 0.10
TF_MODS     <- c("snrnaseq", "geomx", "proteomics")
KIN_MOD     <- "phospho_corrected"
MIN_CELLS   <- 10L

# cbc ships with PuLP under the project .venv; the x86_64 build is the
# `i64` one (the arm64 sibling errors "Exec format error" on this host).
# Auto-locate it when the user did not pass --solver-path so the default
# `--solver cbc` run is turnkey.
locate_cbc <- function() {
  cand <- Sys.glob(file.path(
    ".venv/lib/python*/site-packages/pulp/solverdir/cbc/linux/i64/cbc"))
  cand <- cand[file.exists(cand)]
  if (length(cand) > 0L) return(cand[[1L]])
  p <- Sys.which("cbc")
  if (nzchar(p)) return(unname(p))
  stop("could not auto-locate a cbc binary; pass --solver-path <cbc>",
       call. = FALSE)
}
if (SOLVER == "cbc" && is.null(SOLVER_PATH)) SOLVER_PATH <- locate_cbc()

source("R/helpers.R")  # pulls in R/causal_network.R + tf_inference.R

cache_dir <- "storage/cache"
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
out_path  <- file.path(cache_dir, "causal_network.rds")

if (is_full_run && file.exists(out_path) && !overwrite) {
  cat(sprintf("[build_causal_network] cache exists, skipping: %s\n", out_path))
  cat("Pass --overwrite to rebuild.\n")
  quit(save = "no", status = 0)
}

cat(sprintf("[build_causal_network] solver=%s betaWeight=%.3g fdr=%.2g max_path=%d%s\n",
            SOLVER, BETA_WEIGHT, FDR, MAX_PATH,
            if (SOLVER == "cbc")
              sprintf(" threads=%d timelimit=%gs", THREADS, TIMELIMIT) else ""))
if (SOLVER == "cbc") cat(sprintf("[build_causal_network] cbc: %s\n", SOLVER_PATH))
cat(sprintf("[build_causal_network] contrasts (%s): %s\n",
            if (is_full_run) "full" else "SUBSET/probe -- no cache write",
            paste(run_contrasts, collapse = ", ")))

# ---- per-run CARNIVAL options (uniform workdir + keepLPFiles=FALSE) ----
# Built explicitly here (not delegated to run_carnival_for_contrast's
# internal default) so the betaWeight/workdir handling is visible and the
# cbc solver-path escalation needs no helper change. The J1 housekeeping
# trap: lpSolve defaults pollute the project root with .lp files, so the
# workdir/outputFolder are redirected to a tempdir and keepLPFiles=FALSE.
make_carnival_options <- function() {
  wd <- tempfile("carnival_j3_")
  dir.create(wd, showWarnings = FALSE)
  opts <- switch(
    SOLVER,
    lpSolve = CARNIVAL::defaultLpSolveCarnivalOptions(),
    cbc     = CARNIVAL::defaultCbcSolveCarnivalOptions(),
    cplex   = CARNIVAL::defaultCplexCarnivalOptions(),
    stop(sprintf("unsupported solver '%s'", SOLVER), call. = FALSE)
  )
  opts$outputFolder <- wd
  opts$workdir      <- wd
  opts$keepLPFiles  <- FALSE
  opts$betaWeight   <- BETA_WEIGHT
  if (SOLVER == "cbc") {
    if (is.null(SOLVER_PATH)) {
      stop("--solver cbc requires --solver-path <cbc binary>", call. = FALSE)
    }
    opts$solverPath <- SOLVER_PATH
    opts$threads    <- THREADS
    opts$timelimit  <- TIMELIMIT  # cbc returns best incumbent if hit
  }
  opts
}

# ---- load inputs ------------------------------------------------------
cat("[build_causal_network] loading activity caches...\n")
tf_split   <- readRDS(file.path(cache_dir, "tf_activity_decoupler_split.rds"))
kin        <- readRDS(file.path(cache_dir, "kinase_activity_decoupler.rds"))
symbol_map <- readRDS(file.path(cache_dir, "snrnaseq_symbol_map.rds"))

cat("[build_causal_network] loading PROCESSED microglia subset (~730M)...\n")
t0 <- proc.time()
seurat <- readRDS(file.path(cache_dir, "microglia_seurat_processed.rds"))
cat(sprintf("  loaded in %.1fs\n", (proc.time() - t0)[3]))

# ---- microglia-expressed PKN (J1 dec.3) -------------------------------
cat(sprintf("[build_causal_network] microglia-expressed filter (min_cells=%d)...\n",
            MIN_CELLS))
expr <- microglia_expressed_symbols(seurat, symbol_map, min_cells = MIN_CELLS)
rm(seurat); invisible(gc())
cat(sprintf("  %d microglia-expressed symbols\n", length(expr)))

cat("[build_causal_network] assembling microglia-pruned OmniPath PKN...\n")
pkn <- build_omnipath_pkn(organism = 10090, expr_filter = expr)
pkn_nodes <- unique(c(pkn$source, pkn$target))
cat(sprintf("  PKN: %d edges, %d nodes (raw %s; %s sign-conflicts dropped)\n",
            nrow(pkn), length(pkn_nodes),
            attr(pkn, "n_edges_raw"), attr(pkn, "n_conflicts_dropped")))

# ---- per-contrast CARNIVAL --------------------------------------------
nets <- list()
for (ct in run_contrasts) {
  meas <- tf_meas_from_cache(tf_split, ct, modalities = TF_MODS, fdr = FDR)
  inp  <- kinase_input_from_cache(kin, ct, modality = KIN_MOD, fdr = FDR)
  cat(sprintf("\n[%s] measObj=%d TF(s) (%d in PKN), inputObj=%d kinase(s) (%d in PKN)\n",
              ct, length(meas), sum(names(meas) %in% pkn_nodes),
              length(inp), sum(names(inp) %in% pkn_nodes)))

  # Reachability prune (design pt.4): collapse the intractable full PKN
  # to the nodes on a <= max_path-hop kinase->TF walk for this contrast.
  pkn_ct <- restrict_pkn_to_reachable(pkn, seeds = names(inp),
                                      targets = names(meas), max_path = MAX_PATH)
  cat(sprintf("[%s] reachability prune L=%d: %d edges, %d nodes (seeds %d/%d, targets %d/%d kept)\n",
              ct, MAX_PATH, nrow(pkn_ct),
              length(unique(c(pkn_ct$source, pkn_ct$target))),
              attr(pkn_ct, "n_seeds_kept"),   attr(pkn_ct, "n_seeds_in"),
              attr(pkn_ct, "n_targets_kept"), attr(pkn_ct, "n_targets_in")))

  t0 <- proc.time()
  net <- run_carnival_for_contrast(
    pkn_ct, measObj = meas, inputObj = inp,
    solver = SOLVER, beta_weight = BETA_WEIGHT,
    contrast = ct, fdr = FDR,
    carnival_options = make_carnival_options()
  )
  dt <- (proc.time() - t0)[3]
  cat(sprintf("[%s] status=%s  %.1fs  edges=%d  active_nodes=%d  meas_recovered=%d/%d\n",
              ct, net$meta$status, dt, nrow(net$edges),
              sum(net$nodes$activity != 0),
              sum(net$nodes$activity != 0 & net$nodes$node_type == "measured"),
              net$meta$n_meas_in))
  net$meta$wall_seconds   <- as.numeric(dt)
  net$meta$max_path       <- MAX_PATH
  net$meta$n_seeds_kept   <- attr(pkn_ct, "n_seeds_kept")
  net$meta$n_targets_kept <- attr(pkn_ct, "n_targets_kept")
  nets[[ct]] <- net
}

# ---- summary + assemble -----------------------------------------------
summary_tbl <- dplyr::bind_rows(lapply(nets, summarise_network))
cat("\n[build_causal_network] per-contrast summary:\n")
print(as.data.frame(summary_tbl))

if (!is_full_run) {
  cat("\n[build_causal_network] SUBSET/probe run -- cache NOT written.\n")
  cat("Re-run without --contrasts to build the full cache.\n")
  quit(save = "no", status = 0)
}

out <- nets
out$summary  <- summary_tbl
out$pkn_meta <- list(
  organism            = 10090L,
  n_edges             = nrow(pkn),
  n_nodes             = length(pkn_nodes),
  n_conflicts_dropped = attr(pkn, "n_conflicts_dropped"),
  n_edges_raw         = attr(pkn, "n_edges_raw"),
  n_expr_symbols      = length(expr),
  min_cells           = MIN_CELLS
)
out$params <- list(
  solver          = SOLVER,
  solver_path     = if (SOLVER == "cbc") SOLVER_PATH else NA_character_,
  beta_weight     = BETA_WEIGHT,
  max_path        = MAX_PATH,
  threads         = if (SOLVER == "cbc") THREADS else NA_integer_,
  timelimit       = if (SOLVER == "cbc") TIMELIMIT else NA_real_,
  fdr             = FDR,
  tf_modalities   = TF_MODS,
  kinase_modality = KIN_MOD,
  contrasts       = CONTRASTS_ALL,
  built_from      = c("tf_activity_decoupler_split.rds",
                      "kinase_activity_decoupler.rds",
                      "microglia_seurat_processed.rds",
                      "snrnaseq_symbol_map.rds")
)

cat(sprintf("\n[build_causal_network] writing cache: %s\n", out_path))
saveRDS(out, out_path, compress = "xz")
Sys.chmod(out_path, mode = "0644")
cat("[build_causal_network] done.\n")
