##### Packages #####
library(here)
library(stringr)
library(dplyr)
library(tidytable)
library(fuzzyjoin)
library(tidylog)

##### Loading data #####
nowTaxNames <- read.csv(here("data",
                              "processed",
                              "now_occs_species-names_2025-05-06.csv"))

mddTaxNames <- read.csv(here("data",
                             "processed",
                             "mdd_taxa_species-names_v2.1.csv"))

##### Subsetting data #####

mddTaxNames_species <- mddTaxNames %>%
  filter(MDD_actual_rank == "species") %>%
  distinct()

mddTaxNames_acceptedSpecies <- mddTaxNames_species %>%
  filter(MDD_validity == "species") %>%
  select(MDD_family,
         MDD_species,
         MDD_original_combination,
         MDD_author,
         MDD_year,
         MDD_authority_parentheses) %>%
  distinct()

rm(mddTaxNames)
gc()

##### Loading functions #####

source(here("scripts",
            "harmonization_source-functions.R"))

##### Pre-harmonization #####
#### Setting up default arguments ####

b1_suffix <- "mdd"
b1_col <- "MDD_species"
b2_suffix <- "now"
b2_col <- "now_accepted_name"
select_name_separator <- "_"

#### Merging ####

exact_valid <- harmonize_exact_match(base1_dtf = mddTaxNames_acceptedSpecies,
                                     base2_dtf = nowTaxNames)
exact_valid$exact_summary

exact_valid_eval <- exact_valid$exact_found %>%
  mutate(accepted_match = TRUE) %>%
  mutate(match_notes = paste0(match_notes, " (mdd accepted, now accepted)"))

evaluated_pairs <- exact_valid_eval

finalSynonymy <- add_to_synonymy(data.frame(),
                                 exact_valid_eval)

mddTaxNames_unmatched <- mddTaxNames_species %>%
  anti_join(filter(exact_valid_eval, accepted_match == TRUE))

nowTaxNames_unmatched <- nowTaxNames%>%
  anti_join(filter(exact_valid_eval, accepted_match == TRUE))

rm(exact_valid,
   exact_valid_eval)
gc()

##### Harmonization - MDD synonymy x PBDB accepted names #####
#### Setting up default arguments ####

b1_suffix <- "mdd"
b1_col <- "MDD_original_combination"
b2_suffix <- "now"
b2_col <- "now_accepted_name"

#### Exact ####
exact_mddSyns <- harmonize_exact_match(base1_dtf = mddTaxNames_unmatched,
                                       base2_dtf = nowTaxNames_unmatched)
exact_mddSyns$exact_summary

