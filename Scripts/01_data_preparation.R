###############################################################################################
###               Multi-trait selection defines a plant-architecture ideotype               ###
###                           for soybean in the Russian Far East                           ###
###                                                                                         ###
###                                       ----------                                        ###
###                               Script 01: Data preparation                               ###
###                                                                                         ###
###                            R script by Yawar Habib, Skoltech                            ###
###           Center for Bio and Medical Technologies, ORCID 0009-0008-1812-3910            ###
###############################################################################################

######### PURPOSE #############################################################################
## Load the analysis dataset, set the factor structure, define the canonical trait names,
## labels and abbreviations used throughout, and verify the design before any model is fitted.

######### DATA ################################################################################
## data/Architecture_63genotypes_2018_2020_clean.csv
## 63 soybean accessions evaluated at the Far Eastern Research Institute of Agriculture,
## Vostochnoye, Khabarovsk Krai, over three consecutive seasons (2018-2020) in a randomised
## complete block design with two to three replications.
##
## Scripts 02-08 are self-contained and each load this same file. Run this script first to
## confirm the dataset is intact.

######### PREPARATION OF THE WORKING INTERFACE IN R ###########################################
### I. Set working directory
# Set the working directory to the root of this repository, i.e. the folder that contains
# the data/, Scripts/ and output/ folders.
# On RStudio: tab 'Session' -> Set Working Directory -> Choose Directory.

### II. R packages needed for the analysis
# dplyr, tidyr

### III. Initialisation of the working space
# To erase all graphs
graphics.off()
# To erase objects from the working space - Clean up of the memory
rm(list = ls())
gc()

library(dplyr)
library(tidyr)

## Analysis dataset
Architecture_data <- read.csv("data/Architecture_63genotypes_2018_2020_clean.csv",
                              header = TRUE) %>%
  mutate(
    Year        = as.factor(Year),
    Replication = as.factor(Replication),
    Genotype    = as.factor(Genotype)
  )

head(Architecture_data)
str(Architecture_data)


######### CANONICAL TRAIT DEFINITIONS #########################################################
architecture_traits <- c("Height", "First_pod_height", "No_branches",
                         "No_ProdNodes_MainStem", "No_ProdNodes_LateralBranches",
                         "Vegetation_period")

trait_labels <- c(
  Height                       = "Plant height",
  First_pod_height             = "First pod height",
  No_branches                  = "Number of branches",
  No_ProdNodes_MainStem        = "Productive nodes on main stem",
  No_ProdNodes_LateralBranches = "Productive nodes on lateral branches",
  Vegetation_period            = "Days to maturity"
)

trait_abbreviations <- c(
  Height                       = "PH",
  First_pod_height             = "FPH",
  No_branches                  = "NB",
  No_ProdNodes_MainStem        = "PNMS",
  No_ProdNodes_LateralBranches = "PNLB",
  Vegetation_period            = "DM"
)


######### DESIGN CHECKS #######################################################################

## 1. Observations, genotypes and years
nrow(Architecture_data)
nlevels(Architecture_data$Genotype)
sort(as.integer(levels(Architecture_data$Year)))

## 2. Replications per genotype-year combination
replication_check_per_genotype_per_year <- Architecture_data %>%
  count(Genotype, Year) %>%
  pivot_wider(names_from = Year, values_from = n, values_fill = 0)
replication_check_per_genotype_per_year

## 3. Every genotype present in all three years with at least two replications
all(replication_check_per_genotype_per_year[, -1] >= 2)

## 4. Rows accounted for by the replication table
nrow(Architecture_data) == sum(replication_check_per_genotype_per_year[, -1])

## 5. No missing values in the six traits
colSums(is.na(Architecture_data[, architecture_traits]))

## 6. Replications per cell and their harmonic mean. The harmonic mean is used
##    as r in the genotype-mean phenotypic variance (script 04) and accounts
##    for the unbalanced replication.
replicates_per_cell <- Architecture_data %>%
  count(Genotype, Year, name = "reps") %>%
  pull(reps)

table(replicates_per_cell)
number_of_replicates <- length(replicates_per_cell) / sum(1 / replicates_per_cell)
number_of_replicates

number_of_years <- nlevels(Architecture_data$Year)
number_of_years


cat("\n=== Dataset verified:",
    nrow(Architecture_data), "observations,",
    nlevels(Architecture_data$Genotype), "genotypes,",
    number_of_years, "years, harmonic mean replications =",
    round(number_of_replicates, 4), "===\n")
