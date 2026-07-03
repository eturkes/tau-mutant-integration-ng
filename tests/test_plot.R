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

# --- palette contract: muted defaults + genotype domain pinned to the canonical four ----
sc <- scale_colour_genotype()
sf <- scale_fill_genotype()
sc_default <- ggplot2::scale_colour_discrete()
sf_default <- ggplot2::scale_fill_discrete()
stopifnot(
  identical(tau_discrete_colours,
            c("#56616D", "#3F6F6A", "#9B7A3C", "#7A4052",
              "#647C8A", "#8A6F83", "#7C745F", "#8A8A84")),
  identical(genotype_colours,
            c(MAPTKI = "#56616D", P301S = "#3F6F6A",
              NLGF_MAPTKI = "#9B7A3C", NLGF_P301S = "#7A4052")),
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

# --- substate scales: muted biology colours, canonical domain pinned --------------------
ss_c <- scale_colour_substate(breaks = c("Homeostatic", "DAM", "IFN"))
ss_f <- scale_fill_substate(breaks = c("DAM", "IFN"))
stopifnot(
  identical(substate_colours,
            c(Homeostatic = "#5E7483", DAM = "#7A4052",
              IFN = "#9B7A3C", Proliferative = "#8A6F83")),
  identical(ss_c$aesthetics, "colour"), identical(ss_f$aesthetics, "fill"),
  identical(ss_c$limits, microglia_substate_levels),
  identical(ss_f$limits, microglia_substate_levels),
  identical(ss_c$breaks, c("Homeostatic", "DAM", "IFN")),
  identical(ss_f$breaks, c("DAM", "IFN")),
  identical(ss_c$palette(length(microglia_substate_levels)), substate_colours),
  identical(ss_f$palette(length(microglia_substate_levels)), substate_colours),
  identical(scale_color_substate, scale_colour_substate)
)

# --- binary/direction helpers: stable muted states, no ggplot bright defaults -----------
bg_c <- scale_colour_tau_background()
bin_c <- scale_colour_tau_binary()
bin_f <- scale_fill_tau_binary()
dir_f <- scale_fill_direction()
stopifnot(
  identical(tau_background_colours, c(MAPTKI = "#56616D", P301S = "#3F6F6A")),
  identical(tau_binary_colours, c(`FALSE` = "#8A8A84", `TRUE` = "#315E6F")),
  identical(tau_direction_colours, c(down = "#5E7483", up = "#7A4052")),
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
  identical(unname(rwb_colours), c("#4F6D7A", "#F7F5F0", "#7A4052")),
  identical(tau_sequential_colours, c(low = "#F0EEE8", high = "#3F5F6F")),
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

cat("ok - test_plot\n")
