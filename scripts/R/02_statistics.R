## Statistical analysis: Palmitoylation ========================================

# Are groups normally distributed?
data_enrichment_per_embryo %>%
  filter(
    Protein == "M6",
    Construct %in% c("WT", "Delta-Palm")
  ) %>%
  split(~ Driver + Expression + Construct, drop = TRUE) %>%
  lapply("[[", "Mean_Enrichment") %>%
  lapply(shapiro.test)


# Do groups exhibit equal variances?
data_enrichment_per_embryo |>
  filter(
    Protein == "M6",
    Construct %in% c("WT", "Delta-Palm"),
    Driver != "69B-Gal4"
  ) %>%
  split(~ Driver + Expression, drop = TRUE) %>%
  lapply(\(x) var.test(Mean_Enrichment ~ Construct, x))


# t-Test for all but the overexpression samples.
data_enrichment_per_embryo |>
  filter(
    Protein == "M6",
    Construct %in% c("WT", "Delta-Palm"),
    Driver != "69B-Gal4"
  ) %>%
  split(~ Driver + Expression, drop = TRUE) %>%
  lapply(\(x) t.test(Mean_Enrichment ~ Construct, x, var.equal = TRUE))


# Wilcoxon-test for the overexpression groups.
data_enrichment_per_embryo |>
  filter(
    Protein == "M6",
    Construct %in% c("WT", "Delta-Palm"),
    Driver == "69B-Gal4"
  ) %>%
  split(~ Driver + Expression, drop = TRUE) %>%
  lapply(\(x) wilcox.test(Mean_Enrichment ~ Construct, x))

data_enrichment_per_embryo |>
  filter(
    Protein == "Aka",
    Construct %in% c("WT", "Delta-Palm"),
    Driver != "69B-Gal4"
  ) %>%
  split(~ Driver + Expression, drop = TRUE) %>%
  lapply(\(x) wilcox.test(Mean_Enrichment ~ Construct, x))


## Statistical analysis: Isoforms ==============================================

# Are groups normally distributed?
data_enrichment_per_embryo %>%
  filter(
    Protein == "M6",
    Driver == "69B-Gal4",
    Construct != "Delta-Palm"
  ) %>%
  split(~ Construct, drop = TRUE) %>%
  lapply("[[", "Mean_Enrichment") %>%
  lapply(shapiro.test)

# Pairwise Wilcoxon-Test
wilcoxon_isoforms <- data_enrichment_per_embryo %>%
  filter(
    Protein == "M6",
    Driver == "69B-Gal4",
    Construct != "Delta-Palm"
  ) %>%
  with(pairwise.wilcox.test(Mean_Enrichment, Construct, p.adjust = "none"))

wilcoxon_isoforms[["p.value"]][, 1] |>
  p.adjust(method = "holm") |>
  print()

wilcoxon_nrg <- data_enrichment_per_embryo %>%
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

data_enrichment_per_embryo %>%
  filter(Protein =="Scribble" & Mean_Enrichment <= 6) %>%
  split(~ Construct, drop = TRUE) %>%
  lapply("[[", "Mean_Enrichment") %>%
  lapply(shapiro.test)

bartlett.test(
  Mean_Enrichment ~ Construct,
  data = data_enrichment_per_embryo,
  subset = Protein == "Scribble"
)

data_enrichment_per_embryo %>%
  filter(Protein == "Scribble") %>%
  with(pairwise.t.test(Mean_Enrichment, Construct, p.adj = "none")) %>%
  `[[`("p.value") %>%
  round(3)


model <- lm(
  Mean_Enrichment ~ Mean_Expression,
  data = data_combined,
  subset = Construct != "Delta-Palm"
)
summary(model)