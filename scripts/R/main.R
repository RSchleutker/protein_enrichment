library("lattice")
library("RColorBrewer")
library("gridExtra")
library("DescTools")
library("dplyr")


source_folder <- file.path("scripts", "R")

## Constants ===================================================================

WIDTH <- 24
HEIGHT <- 35
DPI <- 300

# Vectors to change strings to factors in tables frame. The names attribute
# contains the string and the actual value will be the name of the level.
CNSTRCTS <- c(
  "Wt" = "WT",
  "wt" = "WT",
  "Palm" = "Delta-Palm",
  "palm" = "Delta-Palm",
  "Deltapalm" = "Delta-Palm",
  "deltapalm" = "Delta-Palm",
  "Palml200" = "Delta-Palm_L200",
  "palml200" = "Delta-Palm_L200",
  "Pdzb" = "Aka[PDZB]",
  "pdzb" = "Delta-PDZB",
  "pdzb2" = "Gli[PDZB]",
  "pdzb2" = "Double-Delta-PDZB",
  "Pdzb3" = "Aka[PDZB]Gli[PDZB]",
  "L200" = "Aka[L200]",
  "Isoformb" = "Isoform B",
  "Isoformc" = "Isoform C",
  "Isoformd" = "isoform D",
  "Isoforme" = "Isoform E",
  "Isoformf" = "Isoform F",
  "isoformb" = "Isoform B",
  "isoformc" = "Isoform C",
  "isoformd" = "isoform D",
  "isoforme" = "Isoform E",
  "isoformf" = "Isoform F"
)

DRIVERS <- c(
  "Endo" = "Endogenous",
  "endo" = "Endogenous",
  "69b" = "69B-Gal4"
)

PROTEINS <- c(
  "Nrg" = "Nrg",
  "M6" = "M6",
  "Aka" = "Aka",
  "Gli" = "Gli",
  "Scrib" = "Scribble",
  "nrg" = "Nrg",
  "m6" = "M6",
  "aka" = "Aka",
  "gli" = "Gli"
)
EXPRESSION <- c(
  "Hetero" = "Heterozygous",
  "Homo" = "Homozygous",
  "hetero" = "Heterozygous",
  "homo" = "Homozygous"
)


## Load helper functions =======================================================

source(file.path(source_folder, "helpers", "_theme.R"))


# Read and process data --------------------------------------------------------
source(file.path(source_folder, "01_data_cleanup.R"))

# Run statistical analysis -----------------------------------------------------
source(file.path(source_folder, "02_statistics.R"))

# Create and export plots ------------------------------------------------------
source(file.path(source_folder, "03_visualization.R"))
