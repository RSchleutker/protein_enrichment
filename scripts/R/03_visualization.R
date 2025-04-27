## Plot intensities ============================================================

# Plot for transgene enrichment (WT vs. Palmitoylation)
plot_enrichment_transgene <- bwplot(
  Mean_Enrichment ~ Construct,
  data = data_enrichment_per_embryo,
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
) |> print()


# Plot for CRISPR enrichment homozygous (WT vs. Palmitoylation)
plot_enrichment_homozygous <- bwplot(
  Mean_Enrichment ~ Construct,
  data = data_enrichment_per_embryo,
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
) |> print()


# Plot for CRISPR enrichment heterozygous (WT vs. Palmitoylation)
plot_enrichment_heterozygous <- bwplot(
  Mean_Enrichment ~ Construct,
  data = data_enrichment_per_embryo,
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
) |> print()


# Plot for CRISPR enrichment homozygous Aka (WT vs. Palmitoylation)
plot_enrichment_anakonda <- bwplot(
  Mean_Enrichment ~ Construct,
  data = data_enrichment_per_embryo,
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
) |> print()


# Plot for CRISPR enrichment homozygous Aka (WT vs. Palmitoylation)
plot_enrichment_gliotactin <- bwplot(
  Mean_Enrichment ~ Construct,
  data = data_enrichment_per_embryo,
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
) |> print()

attr(plot_enrichment_gliotactin, "width") <- 38
attr(plot_enrichment_gliotactin, "height") <- 53.5


# Plot for different isoforms.
plot_enrichment_isoforms_2 <- bwplot(
  Mean_Enrichment ~ interaction(Protein, Construct, Driver, lex.order = TRUE),
  data = data_enrichment_per_embryo,
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
) |> print()

attr(plot_enrichment_isoforms_2, "width") <- 68
attr(plot_enrichment_isoforms_2, "height") <- 63


# Plot for CRISPR enrichment homozygous Aka (WT vs. Palmitoylation)
plot_enrichment_scribble <- bwplot(
  Mean_Enrichment ~ Construct,
  data = filter(data_enrichment_per_embryo, Protein == "Scribble") %>% droplevels(),
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
) |> print()

attr(plot_enrichment_scribble, "width") <- 52
attr(plot_enrichment_scribble, "height") <- 42


# Plot for CRISPR enrichment homozygous Aka (WT vs. Palmitoylation)
plot_enrichment_anakonda_pdzb <- bwplot(
  Mean_Enrichment ~ Construct | Protein,
  data = filter(
    data_enrichment_per_embryo,
    Protein %in% c("Aka", "Gli"),
    Construct %in% c("WT", "Delta-Palm", "Aka[PDZB]", "Aka[PDZB]Gli[PDZB]")
  ) %>% droplevels(),
  groups = Construct,
  # subset = Protein == "Scribble",
  panel = function(x, y, ...) {
    panel.grid(-1, -1)
    panel.abline(h = 1.9, lty = "dashed")
    
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
) |> print()

attr(plot_enrichment_anakonda_pdzb, "width") <- 38.5
attr(plot_enrichment_anakonda_pdzb, "height") <- 60


## Plot expression data ========================================================
plot_expression_relation <- xyplot(
  Mean_Enrichment ~ Mean_Expression,
  data = data_combined,
  groups = Construct == "Delta-Palm",
  construct = data_combined$Construct,
  panel = \(x, y, construct, ...) {
    panel.grid(-1, -1)
    panel.xyplot(x, y, pch = 16, ...)
    panel.lines(x, predict(model, data.frame(Mean_Expression = x)), col = "black")
    ltext(27, 4.9, expression(R^2 %~~% 0.9483), adj = c(1, .5))
    ltext(27, 4.5, expression(Adj.~R^2 %~~% 0.9354), adj = c(1, .5))
    ltext(
      x, y,
      c(expression(Delta*Palm), "B", "C", "D", "E", "F", "WT"),
      pos = c(1, 3, 1, 1, 3, 1, 1)
    )
  },
  # col = c("#000000", "#ff0000"),
  xlab = "Expression [a.u.]",
  ylab = "Enrichment",
  ylim = c(1.7, 5.3),
  par.settings = theme
) |> print()

attr(plot_expression_relation, "width") <- 65
attr(plot_expression_relation, "height") <- 45.8


trellis.device(
  device = "pdf",
  file = "output/plots/expression_levels.pdf",
  width = 65 / 25,
  height = 51 / 25
)

plot_expression_levels <- bwplot(
  Mean ~ Construct,
  data = data_expression,
  panel = \(x, y, ...) {
    panel.grid(-1, -1)
    panel.bwplot(x, y, ..., outlier = FALSE)
    # panel.xyplot(as.integer(x) + runif(length(x), .26, .4), y, ...,
    #             pch = 16)
    tapply(y, x, FUN = mean) |> panel.xyplot(x = 1:7, cex = .5)
    panel.text(1:7, rep(0, 7), paste0("n=", tapply(y, x, FUN = length)))
  },
  box.ratio = .8,
  xlab = NULL,
  ylab = "Expression [a.u.]",
  scales = list(
    x = list(
      labels = c(
        "69B>M6-GFP",
        expression("69B>M6"^{Delta * "Palm"}*"-GFP"),
        "69B>M6(B)-GFP",
        "69B>M6(C)-GFP",
        "69B>M6(D)-GFP",
        "69B>M6(E)-GFP",
        "69B>M6(F)-GFP"
      ),
      rot = 30
    ),
    y = list(
      alternating = 1,
      limits = c(-3, 39)
    )
  ),
  col = "#00000040",
  par.settings = theme
) |> print()

attr(plot_enrichment_gliotactin, "width") <- 65
attr(plot_enrichment_gliotactin, "height") <- 51



## Export plots ================================================================

# Loop through all objects in the global environment and consider all those with
# names starting with 'plot_[...]'.
for (obj in ls(pattern = "plot_")) {
  current_plot <- get(obj)
  
  width <- attr(current_plot, "Width")
  height <- attr(current_plot, "height")
  
  if(is.null(width)) width <- WIDTH
  if(is.null(height)) height <- HEIGHT
  
  for (fmt in c("pdf", "png")) {
    device <- get(fmt)
    filename <- file.path("output", "plots", paste0(obj, ".", fmt))
    
    if (fmt == "pdf")
      device(file = filename, width / 25, height / 25)
    else
      device(file = filename, width, height, units = "mm", res = DPI)
    
    trellis.par.set(theme)
    print(current_plot)
    
    dev.off()
    
  }
}
