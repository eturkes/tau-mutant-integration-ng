# P4-S1 GeoMx contracts: count/meta extraction, slide + block design, replicate-aware
# limma-voom fit, sensitivity status, and deconvolution preflight provenance.

source("R/constants.R")
source("R/utils.R")
source("R/io.R")
source("R/design.R")
source("R/de_pb.R")
source("R/mechanism.R")
source("R/crossmodality.R")
source("tests/helpers.R")

canonical <- c("tau_alone", "nlgf_in_maptki", "nlgf_in_p301s", "tau_in_nlgf", "interaction")

make_fake_geomx <- function(n_genes = 80L, n_bio = 3L, aoi_per = 2L) {
  combos <- expand.grid(genotype = genotype_levels,
                        bio_rep = seq_len(n_bio),
                        aoi_rep = seq_len(aoi_per),
                        stringsAsFactors = FALSE)
  combos$slide_rep <- ((match(combos$genotype, genotype_levels) + combos$bio_rep + combos$aoi_rep) %% 4L) + 1L
  combos$roi <- paste0("roi", seq_len(nrow(combos)))
  combos$SampleID <- paste(combos$genotype, combos$bio_rep, combos$aoi_rep, sep = "_")
  combos$`ROI Coordinate X` <- seq_len(nrow(combos)) * 10
  combos$`ROI Coordinate Y` <- seq_len(nrow(combos)) * 5
  combos$q_norm_qFactors <- 0.8 + (seq_len(nrow(combos)) %% 5L) / 10
  combos$NegGeoMean_Mm_R_NGS_WTA_v1.0 <- 4 + (seq_len(nrow(combos)) %% 7L)
  combos$nuclei <- ifelse(seq_len(nrow(combos)) %% 11L == 0L, -1, 30 + seq_len(nrow(combos)))
  rownames(combos) <- paste0("AOI", seq_len(nrow(combos)))

  gidx <- match(combos$genotype, genotype_levels)
  counts <- 20 + outer(seq_len(n_genes), seq_len(nrow(combos)),
                       function(i, j) (i * 7L + j * 11L) %% 37L)
  counts[seq_len(12L), ] <- counts[seq_len(12L), ] +
    rep(c(0, 3, 5, 9)[gidx], each = 12L)
  counts[1, ] <- 0L
  rownames(counts) <- paste0("Gene", seq_len(n_genes))
  colnames(counts) <- rownames(combos)
  obj <- SeuratObject::CreateSeuratObject(counts = Matrix::Matrix(counts, sparse = TRUE),
                                          meta.data = combos,
                                          min.cells = 0, min.features = 0)
  obj[["SCT"]] <- SeuratObject::CreateAssay5Object(counts = Matrix::Matrix(counts + 1000L, sparse = TRUE))
  SeuratObject::DefaultAssay(obj) <- "SCT"
  attr(obj, "raw_counts") <- counts
  obj
}

gx <- make_fake_geomx()
raw_counts <- attr(gx, "raw_counts")

# --- count/meta extraction: explicit RNA counts despite SCT default, zero-gene drop -----
cnt <- geomx_count_matrix(gx)
cp <- attr(cnt, "geomx_count_provenance")
stopifnot(nrow(cnt) == nrow(raw_counts) - 1L, ncol(cnt) == ncol(raw_counts),
          identical(colnames(cnt), colnames(raw_counts)),
          cnt[1, 1] == raw_counts[2, 1],                    # proves RNA counts, not SCT+1000
          cp$n_genes_dropped_empty == 1L,
          cp$n_non_integer == 0L, cp$coerced_integer,
          identical(storage.mode(cnt), "integer"))

meta <- geomx_meta(gx)
mp <- attr(meta, "geomx_meta_provenance")
stopifnot(identical(rownames(meta), colnames(cnt)),
          all(c("genotype", "slide", "bio_unit", "roi", "SampleID", "x", "y",
                "q3_factor", "neg_background", "nuclei") %in% names(meta)),
          nlevels(meta$bio_unit) == 4L * 3L,
          mp$nuclei_sentinel_count == sum(meta$nuclei < 0))

# malformed metadata fails at the named missing column, not later in model fitting
gx_bad <- gx
gx_bad@meta.data$bio_rep <- NULL
expect_error(geomx_meta(gx_bad), "bio_rep")

# --- slide design rank guard -------------------------------------------------------------
fd <- geomx_slide_design(meta)
stopifnot(identical(colnames(fd$contrasts), canonical),
          all(genotype_levels %in% colnames(fd$design)))
meta_rank_bad <- meta
meta_rank_bad$slide <- factor(as.character(meta_rank_bad$genotype))
expect_error(geomx_slide_design(meta_rank_bad), "rank-deficient")

# --- replicate-aware fit: duplicateCorrelation primary + sensitivity branches ------------
de <- fit_geomx_de(cnt, meta, min_count = 1L)
stopifnot(de$primary$status == "fit",
          de$primary$duplicate_correlation$used,
          is.finite(de$primary$duplicate_correlation$consensus_correlation),
          identical(names(de$primary$top), canonical),
          de$sensitivity$unblocked$status == "fit",
          !de$sensitivity$unblocked$duplicate_correlation$used,
          de$sensitivity$collapsed_bio_unit$status == "fit",
          de$sensitivity$collapsed_bio_unit$n_bio_units == 12L)
for (tt in de$primary$top) {
  stopifnot(is.data.frame(tt), nrow(tt) == de$primary$kept,
            all(c("symbol", "contrast", "logFC", "P.Value", "adj.P.Val", "t", "CI.L", "CI.R") %in% names(tt)),
            length(unique(tt$contrast)) == 1L)
}

miss <- meta$genotype != "NLGF_P301S"
skip <- fit_geomx_collapsed_sensitivity(cnt[, miss, drop = FALSE], meta[miss, , drop = FALSE], min_count = 1L)
stopifnot(skip$status == "skipped", grepl("rank-deficient", skip$reason, fixed = TRUE))

one_unit <- meta$bio_rep == "1"
sat <- fit_geomx_collapsed_sensitivity(cnt[, one_unit, drop = FALSE], meta[one_unit, , drop = FALSE], min_count = 1L)
stopifnot(sat$status == "skipped", grepl("no residual degrees", sat$reason, fixed = TRUE))

# --- decon preflight: specific defer/block/earned reasons --------------------------------
sp_ok <- list(package = "SpatialDecon", available = TRUE, version = "0.0.0",
              repos = "synthetic", error = NA_character_, warnings = character(), messages = character())
pf <- geomx_decon_preflight(meta, cnt, spatialdecon = sp_ok)
bg_scaled <- geomx_q3_scaled_background(meta)
stopifnot(pf$status == "defer",
          any(grepl("reference profile", pf$reasons, fixed = TRUE)),
          pf$nuclei$n_sentinel > 0L,
          isTRUE(all.equal(as.numeric(bg_scaled), meta$neg_background / meta$q3_factor,
                           tolerance = 1e-12)),
          identical(attr(bg_scaled, "scale"),
                    "negative-probe background divided by q_norm_qFactors"),
          is.finite(pf$memory$estimated_peak_mb))

meta_bad_q3 <- meta
meta_bad_q3$q3_factor[1] <- 0
pf_bad <- geomx_decon_preflight(meta_bad_q3, cnt, spatialdecon = sp_ok)
stopifnot(pf_bad$status == "blocked",
          any(grepl("Q3 normalisation", pf_bad$reasons, fixed = TRUE)))

profile_bad <- cbind(A = 1:6, B = (1:6) * 2, C = c(6, 1, 5, 2, 4, 3))
pf_col <- geomx_decon_preflight(meta, cnt, profile = profile_bad, spatialdecon = sp_ok)
stopifnot(pf_col$status == "blocked",
          any(grepl("collinearity", pf_col$reasons, fixed = TRUE)))

profile_ok <- cbind(A = c(1, 2, 3, 4, 5, 6),
                    B = c(6, 5, 3, 1, 2, 4),
                    C = c(2, 6, 1, 5, 3, 4))
pc_ok <- profile_collinearity(profile_ok)
pf_ok <- geomx_decon_preflight(meta, cnt, profile = profile_ok, spatialdecon = sp_ok)
stopifnot(pf_ok$status == "earned", pf_ok$reference$profile_ok,
          is.finite(pc_ok$max_abs_correlation),
          pc_ok$n_profiles == 3L)

# --- abundance DE design: log beta + slide fixed effect + bio-unit blocking ------------
abund <- rbind(
  microglia = 0.20 + rep(c(0.01, 0.02, 0.04, 0.06), each = 3L * 2L),
  astro = 0.30 + (seq_len(nrow(meta)) %% 5L) / 100,
  neuron = 0.40 - rep(c(0.00, 0.01, 0.02, 0.03), each = 3L * 2L)
)
colnames(abund) <- rownames(meta)
ab_de <- fit_geomx_abundance_de(abund, meta)
stopifnot(ab_de$status == "fit",
          ab_de$duplicate_correlation$used,
          is.finite(ab_de$duplicate_correlation$consensus_correlation),
          identical(names(ab_de$top), canonical),
          nrow(ab_de$top$interaction) == nrow(abund),
          all(c("feature", "contrast", "logFC", "P.Value", "adj.P.Val", "t", "CI.L", "CI.R") %in%
                names(ab_de$top$interaction)),
          identical(ab_de$sensitivity$unblocked$status, "fit"))
expect_error(fit_geomx_abundance_de(abund[, -1, drop = FALSE], meta), "identical")

# --- P4-S2 bulk proteome + corrected phospho ------------------------------------------
bulk_stubs <- paste0("run", sprintf("%02d", 1:16))
bulk_key <- data.frame(
  file_name = paste0(bulk_stubs, ".PTM.Quantity"),
  label = rep(c("MAPT-KI_24M", "P301S+3_24M", "NLGF-MAPT-KI_24M", "NLGF-P301S+3_24M"), each = 4L),
  genotype = factor(rep(genotype_levels, each = 4L), levels = genotype_levels),
  col_stub = bulk_stubs,
  stringsAsFactors = FALSE
)
bulk_effect <- rep(c(0, 8, 15, 28), each = 4L) + rep(c(0, 2, -1, 1), times = 4L)
make_bulk_values <- function(n, base = 100) {
  m <- outer(seq_len(n), seq_len(16), function(i, j) base + i * 19 + bulk_effect[j] +
               ((i * j) %% 7L))
  storage.mode(m) <- "double"
  m
}

prot_ann <- data.frame(
  PG.ProteinGroups = c("P.GSK3B", "P.GSK3B", "P.MAPT", "P.APP", "P.SYN", "P.OTHER1",
                       "P.OTHER2", "P.OTHER3", "P.OTHER4", "P.OTHER5", ""),
  PG.Genes = c("Gsk3b", "Gsk3b", "Mapt", "App", "Syn1;Syp", "Other1",
               "Other2", "Other3", "Other4", "Other5", "Blank"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
prot_vals <- as.data.frame(make_bulk_values(nrow(prot_ann), base = 200))
names(prot_vals) <- paste0(bulk_stubs, ".raw.PTM.Quantity")
prot_vals[1, 1] <- 100
prot_vals[2, 1] <- 20
prot_vals[1, 2] <- 0
prot_vals[2, 2] <- NA_real_
bulk_proteomics <- cbind(prot_ann, prot_vals)

matched_bulk <- match_24m_bulk_columns(bulk_proteomics, bulk_key, modality = "proteome")
stopifnot(length(matched_bulk$columns) == 16L,
          identical(rownames(matched_bulk$meta), bulk_stubs),
          identical(matched_bulk$columns, names(prot_vals)))
expect_error(match_24m_bulk_columns(bulk_proteomics[, -ncol(bulk_proteomics)], bulk_key,
                                    modality = "proteome"),
             "expected exactly")

pfeat <- protein_group_features(bulk_proteomics)
pcounts <- attr(pfeat, "protein_group_counts")
stopifnot(pfeat["P.GSK3B", "n_raw_rows"] == 2L,
          pfeat["P.SYN", "gene_symbols"] == "Syn1;Syp",
          pcounts$n_missing_protein_group == 1L,
          pcounts$n_duplicate_rows_by_group == 1L)
agg <- aggregate_proteome_raw(bulk_proteomics, matched_bulk$columns, pfeat)
stopifnot(agg["P.GSK3B", matched_bulk$columns[1]] == 120,
          is.na(agg["P.GSK3B", matched_bulk$columns[2]]))

prot_prep <- prepare_proteome_24m_matrix(bulk_proteomics, bulk_key, min_present = 1L, min_groups = 4L)
stopifnot(ncol(prot_prep$matrix) == 16L,
          identical(colnames(prot_prep$matrix), rownames(prot_prep$meta)),
          all(c("P.GSK3B", "P.MAPT", "P.APP") %in% rownames(prot_prep$matrix)))
prot_de <- run_proteome_de_24m(bulk_proteomics, bulk_key, min_present = 1L, min_groups = 4L)
stopifnot(prot_de$n_samples == 16L,
          identical(names(prot_de$top), canonical),
          identical(prot_de$run_index$status, "fit"),
          is.list(prot_de$run_index$top))

phos_ann <- data.frame(
  PG.ProteinGroups = c("P.GSK3B", "P.MAPT", "P.APP", "P.SYN", "P.MISSING",
                       "P.OTHER1", "P.OTHER2", "P.OTHER3"),
  PG.Genes = c("Gsk3b", "Mapt", "App", "Syn1", "OtherMissing",
               "Other1", "Other2", "Other3"),
  PTM.SiteAA = c("S", "T", "S", "Y", "S", "S", "T", "Y"),
  PTM.SiteLocation = c(9, 231, 12, 4, 1, 2, 3, 4),
  PTM.CollapseKey = paste0("ck", 1:8),
  `Phosphosite probability` = c(0.9, 0.95, 0.8, 0.7, 0.6, 0.85, 0.83, 0.81),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
phos_vals <- as.data.frame(make_bulk_values(nrow(phos_ann), base = 350))
names(phos_vals) <- paste0(bulk_stubs, ".PTM.Quantity")
bulk_phospho <- cbind(phos_ann, phos_vals)

raw_phos <- run_phospho_de_24m(bulk_phospho, bulk_key, min_present = 1L, min_groups = 4L)
corr_prep <- prepare_phospho_corrected_24m_matrix(bulk_phospho, bulk_key, prot_de,
                                                  min_present = 1L, min_groups = 4L)
stopifnot(corr_prep$counts$n_parent_matched == 7L,
          corr_prep$counts$n_parent_not_in_filtered_proteome == 1L,
          corr_prep$counts$n_missing_corrected_output > 0L,
          "parent_protein_group" %in% names(corr_prep$features),
          !"row5|ck5" %in% rownames(corr_prep$matrix))
phos_prep <- prepare_phospho_24m_matrix(bulk_phospho, bulk_key, min_present = 1L, min_groups = 4L)
fid <- "row1|ck1"
sid <- rownames(corr_prep$meta)[1]
stopifnot(all.equal(corr_prep$matrix[fid, sid],
                    phos_prep$matrix[fid, sid] - prot_de$matrix["P.GSK3B", sid],
                    tolerance = 1e-12))

bad_prot <- prot_de
bad_prot$meta <- bad_prot$meta[16:1, , drop = FALSE]
expect_error(prepare_phospho_corrected_24m_matrix(bulk_phospho, bulk_key, bad_prot,
                                                  min_present = 1L, min_groups = 4L),
             "sample order")

corr_de <- run_phospho_corrected_24m(bulk_phospho, bulk_key, prot_de,
                                     min_present = 1L, min_groups = 4L)
stopifnot(corr_de$n_samples == 16L,
          corr_de$n_features == nrow(corr_prep$matrix),
          identical(names(corr_de$top), canonical),
          identical(corr_de$run_index$status, "fit"),
          is.list(corr_de$run_index$top))

sig <- bulk_significant_counts(prot_de$top, "proteome")
ri <- bulk_run_index_summary(prot_de$top, prot_de$run_index, "proteome")
stopifnot(nrow(sig) == length(canonical), nrow(ri) == length(canonical),
          all(is.finite(sig$n_features)),
          all(ri$status == "fit"))

summary <- bulk_omics_summary_data(prot_de, raw_phos, corr_de)
stopifnot(nrow(summary$feature_counts) == 3L,
          nrow(summary$significant_counts) == 3L * length(canonical),
          summary$feature_counts$n_missing_output[summary$feature_counts$layer == "phospho_corrected"] ==
            corr_de$filters$n_missing_corrected_output,
          all(c("Gsk3b", "Mapt", "Syn1") %in% summary$anchor_coverage$symbol),
          any(summary$anchors$anchor_symbols == "Gsk3b"),
          any(grepl("synaptic", summary$anchors$anchor_class, fixed = TRUE)))

# --- P4-S3 clearance-axis measured table + conservative CCC-lite verdict ---------------
axis_symbols <- c("Apoe", "Trem2", "App", "Cd74", "Pros1", "Mertk", "Syn1", "C1qb")
axis_map <- data.frame(ensembl = paste0("ens", seq_along(axis_symbols)),
                       symbol = axis_symbols, stringsAsFactors = FALSE)
make_axis_top <- function(id_col = "gene", symbol_col = NULL, geomx = FALSE,
                          fdr_pair = 0.01) {
  lfc <- c(Apoe = 1.2, Trem2 = 0.9, App = -0.5, Cd74 = 0.7,
           Pros1 = 0.2, Mertk = 0.4, Syn1 = -1.0, C1qb = 0.8)
  fdr <- rep(0.5, length(lfc)); names(fdr) <- names(lfc)
  fdr[c("Apoe", "Trem2")] <- fdr_pair
  base <- data.frame(logFC = unname(lfc), t = unname(lfc) * 3,
                     P.Value = pmin(unname(fdr) / 2, 0.99),
                     adj.P.Val = unname(fdr), stringsAsFactors = FALSE)
  if (geomx) {
    base <- cbind(symbol = names(lfc), base)
  } else {
    base <- cbind(gene = axis_map$ensembl[match(names(lfc), axis_map$symbol)], base)
  }
  stats::setNames(lapply(canonical, function(cn) base), canonical)
}
axis_pb <- list(top = make_axis_top())
axis_sub <- list(per_substate = list(
  Homeostatic = list(status = "skipped"),
  DAM = list(status = "skipped"),
  IFN = list(status = "skipped"),
  Proliferative = list(status = "skipped")
))
axis_geomx <- list(
  primary = list(top = make_axis_top(geomx = TRUE)),
  decon_preflight = pf
)
axis_sets <- list(sets = list(
  GO_BP = list(GO_SYNAPTIC_SIGNALING = c("Syn1", "Syp", "Snap25", "Other")),
  project = list(DAM = c("Apoe", "Trem2", "Tyrobp"))
))
ca <- clearance_axis_data(axis_pb, axis_sub, axis_map, axis_geomx, summary, axis_sets)
apo <- ca$pair_support[ca$pair_support$pair == "Apoe_Trem2" &
                         ca$pair_support$contrast == "nlgf_in_maptki", , drop = FALSE]
pros <- ca$pair_support[ca$pair_support$pair == "Pros1_Mertk" &
                          ca$pair_support$contrast == "nlgf_in_maptki", , drop = FALSE]
stopifnot(ca$spatial_decon$status == "defer",
          identical(ca$verdict$ccc_called, FALSE),
          ca$verdict$status == "earned",
          nrow(ca$coverage) == nrow(clearance_axis_dictionary()),
          ca$coverage$measured[match("Trem2", ca$coverage$symbol)],
          nrow(ca$synaptic_gene_sets) == 1L,
          apo$status == "earned",
          grepl("GeoMx_spatial", apo$coherent_supported_modalities, fixed = TRUE),
          grepl("snRNAseq_microglia", apo$coherent_supported_modalities, fixed = TRUE),
          pros$status == "not_earned")
axis_geomx_earned <- axis_geomx
axis_geomx_earned$decon_preflight$status <- "earned"
expect_error(clearance_axis_data(axis_pb, axis_sub, axis_map, axis_geomx_earned,
                                 summary, axis_sets),
             "add geomx_decon")

cat("ok - test_crossmodality\n")
