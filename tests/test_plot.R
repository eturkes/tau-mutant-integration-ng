# S4 acceptance: plot layer. Device-free structural checks on the shared ggplot helpers -- no
# graphics device is opened and nothing is drawn. Verifies theme_tau is a ggplot theme, the
# genotype colour/fill scales use ggplot defaults while pinning the canonical domain, the US-spelling alias is the same
# function, and concordance_plot assembles the expected layered ggplot (dropping non-finite rows,
# reporting correlations) without rendering. Mirrors the stopifnot idiom of the other tests/test_*.R.

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

# --- genotype scales: ggplot defaults + domain pinned to the canonical four -------------
sc <- scale_colour_genotype()
sf <- scale_fill_genotype()
sc_default <- ggplot2::scale_colour_discrete()
sf_default <- ggplot2::scale_fill_discrete()
stopifnot(
  identical(sc$aesthetics, "colour"), identical(sf$aesthetics, "fill"),
  identical(sc$limits, genotype_levels), identical(sc$breaks, genotype_levels),  # all 4 in domain + legend
  identical(sf$limits, genotype_levels), identical(sf$breaks, genotype_levels),
  identical(sc$fallback_palette(length(genotype_levels)), sc_default$fallback_palette(length(genotype_levels))),
  identical(sf$fallback_palette(length(genotype_levels)), sf_default$fallback_palette(length(genotype_levels))),
  identical(scale_color_genotype, scale_colour_genotype)                          # US-spelling alias = same fn
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
