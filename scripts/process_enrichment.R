library("lattice")
library("RColorBrewer")
library("gridExtra")
library("DescTools")
library("dplyr")

source("./scripts/_theme.R")


## Constants ===================================================================

# Vectors to change strings to factors in tables frame. The names attribute
# contains the string and the actual value will be the name of the level.
CNSTRCTS <- c(
  "Wt" = "WT",
  "Palm" = "Delta-Palm",
  "Deltapalm" = "Delta-Palm",
  "Palml200" = "Delta-Palm_L200",
  "Pdzb" = "Aka[PDZB]",
  "pdzb2" = "Gli[PDZB]",
  "Pdzb3" = "Aka[PDZB]Gli[PDZB]",
  "L200" = "Aka[L200]",
  "Isoformb" = "Isoform B",
  "Isoformc" = "Isoform C",
  "Isoformd" = "isoform D",
  "Isoforme" = "Isoform E",
  "Isoformf" = "Isoform F"
)
DRIVERS <- c(
  "Endo" = "Endogenous",
  "69b" = "69B-Gal4"
)
PROTEINS <- c(
  "Nrg" = "Nrg",
  "M6" = "M6",
  "Aka" = "Aka",
  "Gli" = "Gli",
  "Scrib" = "Scribble"
)
EXPRESSION <- c(
  "Hetero" = "Heterozygous",
  "Homo" = "Homozygous"
)


## Read and process intensities ================================================

enrichment_per_cell <- list.files("./data", "_measurements.csv", recursive = TRUE, full.names = TRUE) %>%
  lapply(read.csv) %>%
  Reduce(f = rbind) %>%
  mutate(
    Enrichment = (TCJ - BG) / (BCJ - BG),
    Construct = factor(Cons, levels = names(CNSTRCTS), labels = CNSTRCTS),
    Embryo = Emb,
    Cell = Cell,
    Date = as.Date(as.character(Date), format = "%Y%m%d"),
    Protein = factor(Prot, levels = names(PROTEINS), labels = PROTEINS),
    Driver = factor(Driver, levels = names(DRIVERS), labels = DRIVERS),
    Expression = factor(Exp, levels = names(EXPRESSION), labels = EXPRESSION)
  ) %>%
  select(-Cons, -Emb, -Cell, -Date, -Prot, -Exp) %>%
  filter(Enrichment <= 12)


# Variable names in the constructed environment are prefixed with a dot to
# emphasize that these names are not used to construct the column names of
# the final intensities frame.
enrichment_per_embryo <- enrichment_per_cell %>%
  group_by(Protein, Construct, Embryo, Driver, Expression) %>%
  summarise(
    Mean_Enrichment = mean(Enrichment),
    Mean_Intensity = mean(TCJ),
    N = n()
  )


enrichment_per_embryo_summary <- enrichment_per_embryo %>%
  group_by(Protein, Construct, Driver, Expression) %>%
  summarise(
    Mean = mean(Mean_Enrichment),
    Variance = var(Mean_Enrichment),
    Deviation = sd(Mean_Enrichment),
    n = sum(N),
    N = n()
  ) %>%
  mutate(Cells_per_embryo = n / N)


## Export clean intensities ====================================================

write.csv(enrichment_per_cell, "output/data/enrichment_per_cell.csv")
write.csv(enrichment_per_embryo, "output/data/enrichment_per_embryo.csv")
write.csv(enrichment_per_embryo_summary, "output/data/enrichment_summary.csv")

openxlsx::write.xlsx(
  list(enrichment_per_cell = enrichment_per_cell,
       enrichment_per_embryo = enrichment_per_embryo,
       summary = enrichment_per_embryo_summary),
  "output/data/enrichment.xlsx"
)


## Statistical analysis: Palmitoylation ========================================

# Are groups normally distributed?
enrichment_per_embryo %>%
  filter(
    Protein == "M6",
    Construct %in% c("WT", "Delta-Palm")
  ) %>%
  split(~ Driver + Expression + Construct, drop = TRUE) %>%
  lapply("[[", "Mean_Enrichment") %>%
  lapply(shapiro.test)


# Do groups exhibit equal variances?
enrichment_per_embryo |>
  filter(
    Protein == "M6",
    Construct %in% c("WT", "Delta-Palm"),
    Driver != "69B-Gal4"
  ) %>%
  split(~ Driver + Expression, drop = TRUE) %>%
  lapply(\(x) var.test(Mean_Enrichment ~ Construct, x))


# t-Test for all but the overexpression samples.
enrichment_per_embryo |>
  filter(
    Protein == "M6",
    Construct %in% c("WT", "Delta-Palm"),
    Driver != "69B-Gal4"
  ) %>%
  split(~ Driver + Expression, drop = TRUE) %>%
  lapply(\(x) t.test(Mean_Enrichment ~ Construct, x, var.equal = TRUE))


# Wilcoxon-test for the overexpression groups.
enrichment_per_embryo |>
  filter(
    Protein == "M6",
    Construct %in% c("WT", "Delta-Palm"),
    Driver == "69B-Gal4"
  ) %>%
  split(~ Driver + Expression, drop = TRUE) %>%
  lapply(\(x) wilcox.test(Mean_Enrichment ~ Construct, x))

enrichment_per_embryo |>
  filter(
    Protein == "Aka",
    Construct %in% c("WT", "Delta-Palm"),
    Driver != "69B-Gal4"
  ) %>%
  split(~ Driver + Expression, drop = TRUE) %>%
  lapply(\(x) wilcox.test(Mean_Enrichment ~ Construct, x))


## Statistical analysis: Isoforms ==============================================

# Are groups normally distributed?
enrichment_per_embryo %>%
  filter(
    Protein == "M6",
    Driver == "69B-Gal4",
    Construct != "Delta-Palm"
  ) %>%
  split(~ Construct, drop = TRUE) %>%
  lapply("[[", "Mean_Enrichment") %>%
  lapply(shapiro.test)

# Pairwise Wilcoxon-Test
wilcoxon_isoforms <- enrichment_per_embryo %>%
  filter(
    Protein == "M6",
    Driver == "69B-Gal4",
    Construct != "Delta-Palm"
  ) %>%
  with(pairwise.wilcox.test(Mean_Enrichment, Construct, p.adjust = "none"))

wilcoxon_isoforms[["p.value"]][, 1] |>
  p.adjust(method = "holm") |>
  print()

wilcoxon_nrg <- enrichment_per_embryo %>%
  filter(
    Protein %in% c("Nrg", "M6"),
    Driver == "69B-Gal4",
    Construct != "Delta-Palm"
  ) %>%
  with(pairwise.wilcox.test(Mean_Enrichment, Protein:Construct, p.adj = "none"))

wilcoxon_nrg[["p.value"]][, 1] %>%
  p.adjust(method = "holm") %>%
  round(4)


## Statistical analysis: Scribble ==============================================

enrichment_per_embryo %>%
  filter(Protein =="Scribble" & Mean_Enrichment <= 6) %>%
  split(~ Construct, drop = TRUE) %>%
  lapply("[[", "Mean_Enrichment") %>%
  lapply(shapiro.test)

bartlett.test(
  Mean_Enrichment ~ Construct,
  data = enrichment_per_embryo,
  subset = Protein == "Scribble"
)

enrichment_per_embryo %>%
  filter(Protein == "Scribble") %>%
  with(pairwise.t.test(Mean_Enrichment, Construct, p.adj = "none")) %>%
  `[[`("p.value") %>%
  round(3)


## Plot intensities ============================================================

# Plot for transgene enrichment (WT vs. Palmitoylation)
trellis.device(
  device = "pdf",
  file = "output/plots/enrichment_transgene.pdf",
  width = 24 / 25,
  height = 35 / 25
)

bwplot(
  Mean_Enrichment ~ Construct,
  data = enrichment_per_embryo,
  groups = Construct,
  subset = Protein == "M6"
    & Driver == "69B-Gal4"
    & Construct %in% c("WT", "Delta-Palm"),
  panel = function(x, y, ...) {
    panel.grid(-1, -1)
    panel.abline(h = 1.9, lty = "dashed")

    panel.bwplot(x, y, ..., outlier = FALSE)
    panel.xyplot(as.integer(x) + runif(length(x), .26, .4), y, ...,
                 pch = 16)
    tapply(y, x, FUN = mean) |> panel.xyplot(x = 1:2, cex = .5)
    # panel.density(..., adjust = 1)

    panel.text(c(1, 2), c(1.4, 1.4), tapply(y, x, FUN = length))

    panel.lines(c(1, 2), c(11.3, 11.3), col = "#000000")
    panel.text(1.5, 11.8, expression(italic(p < 0.001)))
  },
  box.ratio = .5,
  xlab = NULL,
  ylab = "Enrichment",
  scales = list(
    x = list(labels = c("WT", expression(Delta * "Palm"))),
    y = list(alternating = 1)
  ),
  col = paste0(pal, "40"),
  fill = paste0(pal, "40"),
  do.out = FALSE,
  ylim = c(.8, 12.5),
  par.settings = theme,
)

dev.off()


# Plot for CRISPR enrichment homozygous (WT vs. Palmitoylation)
trellis.device(
  device = "pdf",
  file = "output/plots/enrichment_homozygous.pdf",
  width = 24 / 25,
  height = 35 / 25
)

bwplot(
  Mean_Enrichment ~ Construct,
  data = enrichment_per_embryo,
  groups = Construct,
  subset = Protein == "M6"
    & Driver == "Endogenous"
    & Expression == "Homozygous"
    & Construct %in% c("WT", "Delta-Palm"),
  panel = function(x, y, ...) {
    panel.grid(-1, -1)
    panel.abline(h = 1.9, lty = "dashed")

    panel.bwplot(x, y, ..., outlier = FALSE)
    panel.xyplot(as.integer(x) + runif(length(x), .26, .4), y, ...,
                 pch = 16)
    tapply(y, x, FUN = mean) |> panel.xyplot(x = 1:2, cex = .5)
    # panel.density(..., adjust = 1)

    panel.text(c(1, 2), c(1.4, 1.4), tapply(y, x, FUN = length))

    panel.lines(c(1, 2), c(11.3, 11.3), col = "#000000")
    panel.text(1.5, 11.8, expression(italic(p %~~% 0.102)))
  },
  box.ratio = .5,
  xlab = NULL,
  ylab = "Enrichment",
  scales = list(
    x = list(labels = c("WT", expression(Delta * "Palm"))),
    y = list(alternating = 1)
  ),
  col = paste0(pal, "40"),
  fill = paste0(pal, "40"),
  do.out = FALSE,
  ylim = c(.8, 12.5),
  par.settings = theme,
)

dev.off()


# Plot for CRISPR enrichment heterozygous (WT vs. Palmitoylation)
trellis.device(
  device = "pdf",
  file = "output/plots/enrichment_heterozygous.pdf",
  width = 24 / 25,
  height = 35 / 25
)

bwplot(
  Mean_Enrichment ~ Construct,
  data = enrichment_per_embryo,
  groups = Construct,
  subset = Protein == "M6"
    & Driver == "Endogenous"
    & Expression == "Heterozygous"
    & Construct %in% c("WT", "Delta-Palm"),
  panel = function(x, y, ...) {
    panel.grid(-1, -1)
    panel.abline(h = 1.9, lty = "dashed")

    panel.bwplot(x, y, ..., outlier = FALSE)
    panel.xyplot(as.integer(x) + runif(length(x), .26, .4), y, ...,
                 pch = 16)
    tapply(y, x, FUN = mean) |> panel.xyplot(x = 1:2, cex = .5)
    # panel.density(..., adjust = 1)

    panel.text(c(1, 2), c(1.4, 1.4), tapply(y, x, FUN = length))

    panel.lines(c(1, 2), c(11.3, 11.3), col = "#000000")
    panel.text(1.5, 11.8, expression(italic(p < 0.001)))
  },
  box.ratio = .5,
  xlab = NULL,
  ylab = "Enrichment",
  scales = list(
    x = list(labels = c("WT", expression(Delta * "Palm"))),
    y = list(alternating = 1)
  ),
  col = paste0(pal, "40"),
  fill = paste0(pal, "40"),
  do.out = FALSE,
  ylim = c(.8, 12.5),
  par.settings = theme,
)

dev.off()


# Plot for CRISPR enrichment homozygous Aka (WT vs. Palmitoylation)
trellis.device(
  device = "pdf",
  file = "output/plots/enrichment_anakonda.pdf",
  width = 24 / 25,
  height = 35 / 25
)

bwplot(
  Mean_Enrichment ~ Construct,
  data = enrichment_per_embryo,
  groups = Construct,
  subset = Protein == "Aka"
    & Driver == "Endogenous"
    & Expression == "Homozygous"
    & Construct %in% c("WT", "Delta-Palm"),
  panel = function(x, y, ...) {
    panel.grid(-1, -1)
    panel.abline(h = 1.9, lty = "dashed")

    panel.bwplot(x, y, ..., outlier = FALSE)
    panel.xyplot(as.integer(x) + runif(length(x), .26, .4), y, ...,
                 pch = 16)
    tapply(y, x, FUN = mean) |> panel.xyplot(x = 1:2, cex = .5)
    # panel.density(..., adjust = 1)

    panel.text(c(1, 2), c(1.4, 1.4), tapply(y, x, FUN = length))

    panel.lines(c(1, 2), c(11.3, 11.3), col = "#000000")
    panel.text(1.5, 11.8, expression(italic(p %~~% 0.55)))
  },
  box.ratio = .5,
  xlab = NULL,
  ylab = "Enrichment",
  scales = list(
    x = list(at = c(1, 2), labels = c("WT", expression(Delta * "Palm"))),
    y = list(alternating = 1)
  ),
  col = paste0(pal, "40"),
  fill = paste0(pal, "40"),
  do.out = FALSE,
  xlim = c(.4, 2.6),
  ylim = c(.8, 12.5),
  par.settings = theme
)

dev.off()


# Plot for CRISPR enrichment homozygous Aka (WT vs. Palmitoylation)
trellis.device(
  device = "pdf",
  file = "output/plots/enrichment_gliotactin.pdf",
  width = 38 / 25,
  height = 53.5 / 25
)

bwplot(
  Mean_Enrichment ~ Construct,
  data = enrichment_per_embryo,
  groups = Construct,
  subset = Protein == "Gli"
    & Driver == "Endogenous"
    & Expression == "Homozygous"
    & Construct %in% c("WT", "Delta-Palm"),
  panel = function(x, y, ...) {
    panel.grid(-1, -1)
    panel.abline(h = 1.9)

    panel.bwplot(x, y, ..., outlier = FALSE)
    panel.xyplot(as.integer(x) + runif(length(x), .26, .4), y, ...,
                 pch = 16)

    panel.text(c(1, 2), c(1.2, 1.2), paste0("n=", tapply(y, x, FUN = length)))

    panel.lines(c(1, 2), c(11.5, 11.5), col = "#000000")
    panel.text(1.5, 11.8, expression(italic(p < 0.001)))
  },
  box.ratio = .5,
  xlab = NULL,
  ylab = "Enrichment [a.u.]",
  scales = list(
    x = list(at = c(1, 2), labels = c("WT", expression(Delta * "Palm"))),
    y = list(alternating = 1)
  ),
  col = paste0(pal, "40"),
  fill = paste0(pal, "40"),
  do.out = FALSE,
  xlim = c(.4, 2.6),
  ylim = c(.8, 12.5),
  par.settings = theme,
)

dev.off()


# Plot for different isoforms.
trellis.device(
  device = "pdf",
  file = "output/plots/enrichment_isoforms_2.pdf",
  width = 68 / 25,
  height = 63 / 25
)

bwplot(
  Mean_Enrichment ~ interaction(Protein, Construct, Driver, lex.order = TRUE),
  data = enrichment_per_embryo,
  subset = (
    Protein %in% c("Nrg", "M6") &
    Construct != "Delta-Palm" &
    (Driver == "69B-Gal4" |
      (Driver == "Endogenous" & Expression == "Homozygous"))
  ),
  panel = function(x, y, ...) {
    panel.grid(-1, -1)
    panel.abline(h = 1.9, lty = "dashed")

    panel.bwplot(x, y, ..., outlier = FALSE)
    tapply(y, x, FUN = mean) |> panel.xyplot(x = 1:8, cex = .5)
    panel.xyplot(as.integer(x) + runif(length(x), .26, .4), y, ...,
                 pch = 16)
    # panel.density(..., adjust = 1)

    panel.text(1:8, rep(1.2, 7), paste0("n=", tapply(y, x, FUN = length)))
    (tapply(y, x, FUN = max) + 1.5) |>
      panel.text(x = 1:8, cex = 1, labels = c(
        "", "***", "***", "*", "***", "***", "***", "***"
    ))
  },
  box.ratio = .5,
  xlab = NULL,
  ylab = "Enrichment",
  scales = list(
    x = list(
      labels = c(
        "69B>Nrg-GFP",
        "GFP::M6",
        "69B>M6-GFP",
        "69B>M6(B)-GFP",
        "69B>M6(C)-GFP",
        "69B>M6(D)-GFP",
        "69B>M6(E)-GFP",
        "69B>M6(F)-GFP"
      ),
      rot = 30
    ),
    y = list(alternating = 1)
  ),
  col = "#00000040",
  # fill = paste0(pal, "40"),
  do.out = FALSE,
  ylim = c(.8, 12.5),
  par.settings = theme,
)

dev.off()


# Plot for CRISPR enrichment homozygous Aka (WT vs. Palmitoylation)
trellis.device(
  device = "pdf",
  file = "output/plots/enrichment_scribble.pdf",
  width = 52 / 25,
  height = 42 / 25
)

bwplot(
  Mean_Enrichment ~ Construct,
  data = filter(enrichment_per_embryo, Protein == "Scribble") %>% droplevels(),
  groups = Construct,
  # subset = Protein == "Scribble",
  panel = function(x, y, ...) {
    panel.grid(-1, -1)
    # panel.abline(h = 1.9, lty = "dashed")

    panel.bwplot(x, y, ...)
    panel.xyplot(as.integer(x) + .33, y, jitter.x = TRUE, amount = .1, ...)
    tapply(y, x, FUN = mean) |> panel.xyplot(x = 1:3, cex = .5)

    panel.text(1:3, rep(.2, 3), tapply(y, x, FUN = length))

    panel.text(2, 6.1, expression(italic(p %~~% 0.01)))
    panel.text(3, 6.1, expression(italic(p < 0.001)))
  },
  box.ratio = .4,
  # main = "Scribble",
  xlab = NULL,
  ylab = "Enrichment",
  col = "#00000040",
  fill = "#ffffff",
  do.out = FALSE,
  # xlim = c(.4, 2.6),
  ylim = c(-.2, 7.2),
  scales = list(
    x = list(
      at = 1:3,
      labels = c(
        "Control",
        expression(italic("aka"^{Delta*"PDZB"})),
        expression(italic("aka"^"L200"))
      ),
      rot = 30
    ),
    y = list(alternating = 1, tck=c(1, 0))
  ),
  par.settings = theme
)

dev.off()


# Plot for CRISPR enrichment homozygous Aka (WT vs. Palmitoylation)
trellis.device(
  device = "pdf",
  file = "output/plots/enrichment_anakonda_pdzb.pdf",
  width = 38.5 / 25,
  height = 60 / 25
)

bwplot(
  Mean_Enrichment ~ Construct | Protein,
  data = filter(
    enrichment_per_embryo,
    Protein %in% c("Aka", "Gli"),
    Construct %in% c("WT", "Delta-Palm", "Aka[PDZB]", "Aka[PDZB]Gli[PDZB]")
  ) %>% droplevels(),
  groups = Construct,
  # subset = Protein == "Scribble",
  panel = function(x, y, ...) {
    panel.grid(-1, -1)
    panel.abline(h = 1.9, lty = "dashed")

    print(x)

    panel.bwplot(x, y, ...)
    panel.xyplot(as.integer(x) + .33, y, jitter.x = TRUE, amount = .1, ...)
    tapply(y, x, FUN = mean) |> panel.xyplot(x = 1:4, cex = .5)

    panel.text(1:4, rep(1.4, 4), tapply(y, x, FUN = length))

    # panel.lines(c(1, 2), c(5.8, 5.8), col = "#000000")
    # panel.text(1.5, 6.1, expression(italic(p %~~% 0.01)))
  },
  box.ratio = .4,
  main = "Anakonda",
  xlab = NULL,
  ylab = "Enrichment",
  col = paste0(pal, "40"),
  fill = paste0(pal, "40"),
  do.out = FALSE,
  # xlim = c(.4, 2.6),
  ylim = c(.8, 11.8),
  scales = list(
    x = list(
      at = 1:4,
      labels = c(
        "Control",
        expression("M6["*Delta*"Palm]"),
        expression("aka["*Delta*"PDZB]"),
        expression("aka["*Delta*"PDZB] Gli["*Delta*"PDZB]")
      ),
      rot = 30
    ),
    y = list(alternating = 1, tck=c(1, 0))
  ),
  par.settings = theme
)

dev.off()