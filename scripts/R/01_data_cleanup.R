## Read and process intensities ================================================

data_enrichment_per_cell <- list.files("data", "_measurements.csv", recursive = TRUE, full.names = TRUE) %>%
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
data_enrichment_per_embryo <- data_enrichment_per_cell %>%
  group_by(Protein, Construct, Embryo, Driver, Expression) %>%
  summarise(
    Mean_Enrichment = mean(Enrichment),
    Mean_Intensity = mean(TCJ),
    N = n()
  )


data_enrichment_per_embryo_summary <- data_enrichment_per_embryo %>%
  group_by(Protein, Construct, Driver, Expression) %>%
  summarise(
    Mean = mean(Mean_Enrichment),
    Variance = var(Mean_Enrichment),
    Deviation = sd(Mean_Enrichment),
    n = sum(N),
    N = n()
  ) %>%
  mutate(Cells_per_embryo = n / N)


## Read and process expression data ============================================

data_expression <- read.csv("data/expression/expression.csv") %>%
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

data_expression_summary <- data_expression %>%
  group_by(Protein, Driver, Construct) %>%
  summarise(
    Mean_Expression = mean(Mean)
  )

data_enrichment_summary <- data_enrichment_per_embryo %>%
  filter(Protein == "M6", Driver == "69B-Gal4") %>%
  group_by(Protein, Driver, Construct) %>%
  summarise(Mean_Enrichment = mean(Mean_Enrichment))

data_combined <- merge(data_enrichment_summary, data_expression_summary)


## Export clean intensities ====================================================

write.csv(data_enrichment_per_cell, file.path("output", "tables", "enrichment_per_cell.csv"))
write.csv(data_enrichment_per_embryo, file.path("output", "tables", "enrichment_per_embryo.csv"))
write.csv(data_enrichment_per_embryo_summary, file.path("output", "tables", "enrichment_summary.csv"))

openxlsx::write.xlsx(
  list(enrichment_per_cell = data_enrichment_per_cell,
       enrichment_per_embryo = data_enrichment_per_embryo,
       summary = data_enrichment_per_embryo_summary),
  file.path("output", "tables", "enrichment.xlsx")
)
