# S4 acceptance: plot layer. Device-free structural checks on the shared ggplot helpers -- no
# graphics device is opened and nothing is drawn. Verifies theme_tau is a ggplot theme, the
# genotype colour/fill scales pin the restrained project palette + domain, the generic ggplot
# discrete defaults inherit the same palette after theme_tau(), microglia-substate scales keep
# the biology colours stable, binary/direction helpers avoid bright defaults, and concordance_plot
# assembles the expected layered ggplot without rendering.

source("R/constants.R")
source("R/utils.R")
source("R/plot.R")
source("tests/helpers.R")

# --- %||% (R/utils.R) underpins concordance_plot's label defaults -----------------------
stopifnot(identical("a" %||% "b", "a"), identical(NULL %||% "b", "b"))

# --- theme_tau: a ggplot theme object, args plumb through -------------------------------
stopifnot(inherits(theme_tau(), "theme"), inherits(theme_tau(), "gg"),
          inherits(theme_tau(base_size = 14, base_family = "serif"), "theme"),
          theme_tau(base_size = 14)$text$size == 14)   # base_size plumbs through (value, not just class)

# --- palette contract: saturated defaults + genotype domain pinned to the canonical four -
sc <- scale_colour_genotype()
sf <- scale_fill_genotype()
sc_default <- ggplot2::scale_colour_discrete()
sf_default <- ggplot2::scale_fill_discrete()
stopifnot(
  identical(tau_discrete_colours,
            c("#3F5F7F", "#0B7A75", "#C8841C", "#A63A50",
              "#2F7EA8", "#7D5CB8", "#8A6A32", "#6F7782")),
  identical(genotype_colours,
            c(MAPTKI = "#3F5F7F", P301S = "#0B7A75",
              NLGF_MAPTKI = "#C8841C", NLGF_P301S = "#A63A50")),
  identical(getOption("ggplot2.discrete.colour"), tau_discrete_scale_types),
  identical(getOption("ggplot2.discrete.fill"), tau_discrete_scale_types),
  identical(sc$aesthetics, "colour"), identical(sf$aesthetics, "fill"),
  identical(sc$limits, genotype_levels), identical(sc$breaks, genotype_levels),  # all 4 in domain + legend
  identical(sf$limits, genotype_levels), identical(sf$breaks, genotype_levels),
  identical(sc$palette(length(genotype_levels)), genotype_colours),
  identical(sf$palette(length(genotype_levels)), genotype_colours),
  identical(sc_default$fallback_palette(4), tau_discrete_colours[1:4]),
  identical(sf_default$fallback_palette(4), tau_discrete_colours[1:4]),
  identical(scale_color_genotype, scale_colour_genotype)                          # US-spelling alias = same fn
)

# --- substate scales: saturated biology colours, canonical domain pinned ----------------
ss_c <- scale_colour_substate(breaks = c("Homeostatic", "DAM", "IFN"))
ss_f <- scale_fill_substate(breaks = c("DAM", "IFN"))
stopifnot(
  identical(substate_colours,
            c(Homeostatic = "#2F78A0", DAM = "#A63A50",
              IFN = "#C8841C", Proliferative = "#7D5CB8")),
  identical(ss_c$aesthetics, "colour"), identical(ss_f$aesthetics, "fill"),
  identical(ss_c$limits, microglia_substate_levels),
  identical(ss_f$limits, microglia_substate_levels),
  identical(ss_c$breaks, c("Homeostatic", "DAM", "IFN")),
  identical(ss_f$breaks, c("DAM", "IFN")),
  identical(ss_c$palette(length(microglia_substate_levels)), substate_colours),
  identical(ss_f$palette(length(microglia_substate_levels)), substate_colours),
  identical(scale_color_substate, scale_colour_substate)
)

# --- binary/direction helpers: stable states, no ggplot bright defaults -----------------
bg_c <- scale_colour_tau_background()
bin_c <- scale_colour_tau_binary()
bin_f <- scale_fill_tau_binary()
dir_f <- scale_fill_direction()
stopifnot(
  identical(tau_background_colours, c(MAPTKI = "#3F5F7F", P301S = "#0B7A75")),
  identical(tau_binary_colours, c(`FALSE` = "#7C838A", `TRUE` = "#0B6F7E")),
  identical(tau_direction_colours, c(down = "#2F78A0", up = "#A63A50")),
  identical(bg_c$limits, names(tau_background_colours)),
  identical(bin_c$palette(2), tau_binary_colours),
  identical(bin_f$palette(2), tau_binary_colours),
  identical(dir_f$palette(2), tau_direction_colours),
  identical(scale_color_tau_background, scale_colour_tau_background),
  identical(scale_color_tau_binary, scale_colour_tau_binary)
)

# --- heatmap scales: neutral sequential counts + teal/wine signed panels ---------------
rwb_fill_seq <- scale_fill_rwb(name = "x")
rwb_fill_mid <- scale_fill_rwb(name = "x", midpoint = 0)
rwb_col_mid <- scale_colour_rwb(name = "x", midpoint = 0)
stopifnot(
  identical(unname(rwb_colours), c("#2F78A0", "#F8F5ED", "#A63A50")),
  identical(tau_sequential_colours, c(low = "#F1EEE5", high = "#1F6F8B")),
  inherits(rwb_fill_seq, "ScaleContinuous"), inherits(rwb_fill_mid, "ScaleContinuous"),
  inherits(rwb_col_mid, "ScaleContinuous"),
  identical(rwb_fill_seq$aesthetics, "fill"), identical(rwb_fill_mid$aesthetics, "fill"),
  identical(rwb_col_mid$aesthetics, "colour"),
  identical(scale_color_rwb, scale_colour_rwb)
)

# --- concordance_plot: expected layered ggplot; finite-row filter; correlations reported -
df <- data.frame(
  gene = paste0("g", 1:20),
  x    = seq(-2, 2, length.out = 20),
  y    = seq(-2, 2, length.out = 20) * 1.5 + rep(c(-0.1, 0.1), 10),   # deterministic, correlated
  stringsAsFactors = FALSE
)
p <- concordance_plot(df, "x", "y", x_lab = "x effect", y_lab = "y effect", title = "T")
stopifnot(
  inherits(p, "ggplot"), inherits(p, "gg"),
  length(p$layers) == 5L,                              # hline, vline, point, smooth (lm), text_repel
  identical(p$labels$x, "x effect"),                   # explicit labs honoured
  identical(p$labels$title, "T"),
  grepl("Spearman", p$labels$subtitle, fixed = TRUE),  # both correlations computed into the subtitle
  grepl("Pearson",  p$labels$subtitle, fixed = TRUE),
  grepl("n = 20",   p$labels$subtitle, fixed = TRUE)
)
# default labels fall back to the column names via %||%
p_def <- concordance_plot(df, "x", "y")
stopifnot(identical(p_def$labels$x, "x"), identical(p_def$labels$y, "y"))
# non-finite rows are dropped before plotting (NA + Inf -> 20 finite of 22)
df_bad <- rbind(df, data.frame(gene = c("gNA", "gInf"), x = c(NA, Inf), y = c(1, 1),
                               stringsAsFactors = FALSE))
p_bad <- concordance_plot(df_bad, "x", "y")
stopifnot(nrow(p_bad$data) == 20L, grepl("n = 20", p_bad$labels$subtitle, fixed = TRUE))

# --- modality_interaction_scatter: 7-layer amyloid-response panel; y=x line; finite filter;
#     symmetric coord_equal; axis labels pin y=tau-KO effect and x=mutant-tau effect ----------
mdf <- data.frame(
  label = paste0("f", 1:20),
  x     = seq(-2, 2, length.out = 20),
  y     = seq(-2, 2, length.out = 20) + rep(c(-0.05, 0.05), 10),   # near-diagonal, correlated
  stringsAsFactors = FALSE
)
mp <- modality_interaction_scatter(mdf, title = "M")
stopifnot(
  inherits(mp, "ggplot"), inherits(mp, "gg"),
  length(mp$layers) == 7L,                             # hline, vline, abline(y=x), point, smooth, point(top), text_repel
  isTRUE(mp$coordinates$ratio == 1),                   # coord_equal (1:1 aspect) -> diagonal reads at 45 deg
  identical(mp$labels$title, "M"),
  grepl("NLGF_P301S", mp$labels$x, fixed = TRUE),      # x-axis = amyloid effect on the mutant-tau background
  grepl("NLGF_MAPTKI", mp$labels$y, fixed = TRUE),     # y-axis = amyloid effect on the tau-KO background
  grepl("Spearman", mp$labels$subtitle, fixed = TRUE),
  grepl("Pearson",  mp$labels$subtitle, fixed = TRUE),
  grepl("n = 20",   mp$labels$subtitle, fixed = TRUE)
)
# non-finite rows dropped before plotting (NA + Inf -> 20 finite of 22)
mdf_bad <- rbind(mdf, data.frame(label = c("fNA", "fInf"), x = c(NA, Inf), y = c(1, 1),
                                 stringsAsFactors = FALSE))
mp_bad <- modality_interaction_scatter(mdf_bad)
stopifnot(nrow(mp_bad$data) == 20L, grepl("n = 20", mp_bad$labels$subtitle, fixed = TRUE))
cat("ok - modality_interaction_scatter: 7-layer y=x scatter, finite filter, correlations\n")

# --- plate_support_matrix: cross-modality bubble matrix; builds warning-free with & without a
#     not-observed layer. ggplot_build forces scale training so a shape scale left dangling when
#     `missing` is empty would raise "no shared levels" -> under warn=2 that is an error here. ----
mods <- c("snRNAseq", "GeoMx", "phospho")
feats <- c("Cst7", "Apoe")
psm_measured <- data.frame(
  modality_label = factor(c("snRNAseq", "GeoMx", "snRNAseq", "GeoMx"), levels = mods),
  feature_label_plot = factor(c("Apoe", "Apoe", "Cst7", "Cst7"), levels = feats),
  plot_effect = c(1.5, 0.9, -0.8, 2.1),
  plot_status = c("supported/earned", "measured, not supported",
                  "supported/earned", "measured, not supported"),
  stringsAsFactors = FALSE)
psm_missing <- data.frame(
  modality_label = factor(c("phospho", "phospho"), levels = mods),
  feature_label_plot = factor(c("Apoe", "Cst7"), levels = feats),
  plot_effect = NA_real_,
  plot_status = c("not observed", "blocked"),
  stringsAsFactors = FALSE)
p_full <- plate_support_matrix(psm_measured, psm_missing, effect_name = "effect")
p_nomiss <- plate_support_matrix(psm_measured, missing = NULL)
p_emptymiss <- plate_support_matrix(psm_measured, psm_missing[0, , drop = FALSE])
stopifnot(
  inherits(p_full, "ggplot"), inherits(p_nomiss, "ggplot"),
  length(p_full$layers) == 2L,     # measured bubbles + not-observed glyphs
  length(p_nomiss$layers) == 1L,   # measured bubbles only, no dangling shape layer
  length(p_emptymiss$layers) == 1L,
  inherits(ggplot2::ggplot_build(p_full)$plot, "ggplot"),      # trains scales, no warning at warn=2
  inherits(ggplot2::ggplot_build(p_nomiss)$plot, "ggplot"),
  inherits(ggplot2::ggplot_build(p_emptymiss)$plot, "ggplot")
)

# --- plate_pair_matrix: pair x contrast count matrix; builds warning-free incl. an all-not-earned
#     panel (only the FALSE ring level present). --------------------------------------------------
ppm <- data.frame(
  x = factor(c("amyloid on P301S", "amyloid on MAPTKI"),
             levels = c("amyloid on MAPTKI", "amyloid on P301S")),
  y = factor(c("Apoe-Trem2", "App-Cd74"), levels = c("App-Cd74", "Apoe-Trem2")),
  count = c(2L, 1L), earned = c(TRUE, FALSE), stringsAsFactors = FALSE)
p_pair <- plate_pair_matrix(ppm)
ppm0 <- transform(ppm, count = c(0L, 0L), earned = c(FALSE, FALSE))
p_pair0 <- plate_pair_matrix(ppm0)
stopifnot(
  inherits(p_pair, "ggplot"), length(p_pair$layers) == 2L,     # bubbles + count text
  inherits(ggplot2::ggplot_build(p_pair)$plot, "ggplot"),
  inherits(ggplot2::ggplot_build(p_pair0)$plot, "ggplot")
)

cat("ok - test_plot\n")
