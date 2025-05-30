#### Loading packages ####

library(here) #Operating system agnostic file paths with .Rproj
library(stringr) #String manipulation
library(dplyr) #Data wrangling with tidyverse syntax
library(tidytable) #Brings the data.table package's fast operations to dplyr syntax
library(tidylog) #Optional, prints detailed information about changes in the dataframe during data wrangling with tidyverse functions

#### Setting file paths ####

mdd.synonymy.path <- here("data",
                          "raw",
                          "mdd_taxa_v2.1_raw.csv")

mdd.functions.path <- here("scripts",
                           "functions",
                           "cleaning_mdd.R")

#Optional
mdd.main.path <- here("data",
                      "raw",
                      "mdd_species_v2.1_raw.csv")

#### Loading data and functions####

mdd.synonymy.data <- read.csv(mdd.synonymy.path,
                              na.strings = "")

source(mdd.functions.path)

#Optional
mdd.main.data <- read.csv(mdd.main.path,
                              na.strings = c("NA", ""))

rm(mdd.synonymy.path,
   mdd.main.path,
   mdd.functions.path)
gc()

#### Preparing function arguments ####

select.group_col <- "MDD_order"
select.group_name <- "Carnivora" #Uncomment bellow to check available groups
#available_orders <- unique(mdd.synonymy.data$MDD_order)
#available_families <- distinct(select(mdd.synonymy.data, MDD_order, MDD_family))
#available_genera <- distinct(select(mdd.synonymy.data, MDD_order, MDD_family, MDD_genus))

select.original_rank_col <- "MDD_original_rank"
select.original_rank_name <- c("species", "synonym_species") #Uncomment bellow to check available original ranks
#available_ranks <- unique(mdd.synonymy.data$MDD_original_rank)

select.species_col <- "MDD_species"
select.original_combination_col <- "MDD_original_combination"

select.cols_to_keep <- c("MDD_species_id",
                         "MDD_syn_ID",
                         "MDD_order",
                         "MDD_family",
                         "MDD_genus",
                         "MDD_species",
                         "MDD_validity",
                         "MDD_original_combination",
                         "MDD_author",
                         "MDD_year",
                         "MDD_authority_parentheses") #Uncomment bellow to check available columns
#colnames(mdd.synonymy.data)
select.new_name_separator <- "_"

#### Filtering data ####

mdd.synonymy.data <- mdd.filter_synonyms(data = mdd.synonymy.data,
                                         species_col = select.species_col,
                                         group_col = select.group_col,
                                         group_name = select.group_name,
                                         original_rank_col = select.original_rank_col,
                                         original_rank_name = select.original_rank_name,
                                         cols_to_keep = select.cols_to_keep)
gc()

#### Original rank == c("species", "synonym_species") ####

mdd.synonyms <- mdd.format_synonyms(data = mdd.synonymy.data,
                                    species_col = select.species_col,
                                    original_combination_col = select.original_combination_col,
                                    sep = select.new_name_separator) #Generates a list with 2 items: 1) mdd_syns_final = correctly formatted names; 2) mdd_syns_manual = names in need of manual cleaning
gc()

#### Optional step - Continent ####

mdd.synonyms$final <- mdd.main.data %>%
  select(id,
         continentDistribution) %>%
  rename(MDD_species_id = id,
         MDD_continent = continentDistribution) %>%
  inner_join(mdd.synonyms$final, .)

#### Write CSV ####

write.csv(mdd.synonyms$final,
          file = here("data",
                      "processed",
                      "mdd_taxa_species-names_v2.1.csv"),
          row.names = FALSE)

rm(list = ls())

gc()
