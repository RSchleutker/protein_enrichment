library("lattice")
library("RColorBrewer")
library("gridExtra")
library("DescTools")
library("dplyr")

source("src/R/theme.R")


CNSTRCTS <- c(
  "wt" = "WT",
  "palm" = "Delta-Palm",
  "deltapalm" = "Delta-Palm",
  "palml200" = "Delta-Palm_L200",
  "pdzb" = "Delta-PDZB",
  "pdzb2" = "Double-Delta-PDZB",
  "isoformb" = "Isoform B",
  "isoformc" = "Isoform C",
  "isoformd" = "isoform D",
  "isoforme" = "Isoform E",
  "isoformf" = "Isoform F"
)
DRIVERS <- c(
  "endo" = "Endogenous",
  "69b" = "69B-Gal4"
)
PROTEINS <- c(
  "nrg" = "Nrg",
  "m6" = "M6",
  "aka" = "Aka",
  "gli" = "Gli"
)
EXPRESSION <- c(
  "hetero" = "Heterozygous",
  "homo" = "Homozygous"
)


expression <- read.csv("data/expression/expression.csv") %>%
  rename(
    Date = date,
    Protein = prot,
    Driver = driver,
    Construct = cons,
    Expression = exp,
    Embryo = emb
  ) %>%
  mutate(
    Date = as.Date(as.character(Date), format = "%Y%m%d"),
    Protein = factor(Protein, levels = names(PROTEINS), labels = PROTEINS),
    Driver = factor(Driver, levels = names(DRIVERS), labels = DRIVERS),
    Expression = factor(Expression, levels = names(EXPRESSION), labels = EXPRESSION),
    Construct = factor(Construct, levels = names(CNSTRCTS), labels = CNSTRCTS)
  )

expression_summary <- expression %>%
  group_by(Protein, Driver, Construct) %>%
  summarise(
    Mean_Expression = mean(Mean)
  )

enrichment_summary <- enrichment_per_embryo %>%
  filter(Protein == "M6", Driver == "69B-Gal4") %>%
  group_by(Protein, Driver, Construct) %>%
  summarise(Mean_Enrichment = mean(Mean_Enrichment))

combined <- merge(enrichment_summary, expression_summary)

model <- lm(Mean_Enrichment ~ Mean_Expression, data = combined, subset = Construct != "Delta-Palm")
summary(model)


trellis.device(
  device = "pdf",
  file = "output/plots/expression_enrichment_relation.pdf",
  width = 65 / 25,
  height = 45.8 / 25
)

xyplot(
  Mean_Enrichment ~ Mean_Expression,
  data = combined,
  groups = Construct == "Delta-Palm",
  construct = combined$Construct,
  panel = \(x, y, construct, ...) {
    print(construct)
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
)

dev.off()


trellis.device(
  device = "pdf",
  file = "output/plots/expression_levels.pdf",
  width = 65 / 25,
  height = 51 / 25
)

bwplot(
  Mean ~ Construct,
  data = expression,
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
)

dev.off()
