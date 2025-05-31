# _Info -------------------------------------------------------------------
##
## Title: harmonization_now-mdd.R
## Purpose:
##
## Author: Thaís G. P. Faria
## Github: https://github.com/Thais-Faria
## Date created: 2025-05-31
## Copyright (c) Thaís G. P. Faria, 2025
## License: GNU General Public License version 3
##
# _Set up -----------------------------------------------------------------
# __Packages --------------------------------------------------------------

library(here)
library(stringr)
library(dplyr)
library(tidytable)
library(fuzzyjoin)
library(tidylog)

# __File paths ------------------------------------------------------------

mdd.path <- here("data",
                   "processed",
                   "mdd_taxa_species-names_v2.1.csv")

now.path <- here("data",
                   "processed",
                   "now_occs_species-names_2025-05-06.csv")

functions.path <- here("scripts",
                       "harmonization_source-functions.R")

# __Loading functions -----------------------------------------------------

source(functions.path)

# __Loading data ----------------------------------------------------------

mddTaxNames <- read.csv(mdd.path)

nowOccNames <- read.csv(now.path)

# __Sub-setting data ------------------------------------------------------

# MDD
mddTaxNames_epithet <- mddTaxNames %>%
  filter(MDD_actual_rank == "epithet") %>%
  distinct()

mddTaxNames_subspecies <- mddTaxNames %>%
  filter(MDD_actual_rank == "subspecies") %>%
  distinct()

mddTaxNames_species <- mddTaxNames %>%
  filter(MDD_actual_rank == "species") %>%
  distinct()

mddTaxNames_acceptedSpecies <- mddTaxNames_species %>%
  filter(MDD_validity == "species") %>%
  select(MDD_species_id,
         MDD_syn_ID,
         MDD_family,
         MDD_species,
         MDD_original_combination,
         MDD_author,
         MDD_year,
         MDD_authority_parentheses,
         MDD_continent) %>%
  distinct()

# NOW
nowOccNames_0to2max <- nowOccNames %>%
  filter(MAX_AGE <= 2)

nowOccNames_2to5max <- nowOccNames %>%
  filter(MAX_AGE > 2 & MAX_AGE <= 5)

nowOccNames_0to1min <- nowOccNames %>%
  filter(MAX_AGE > 5 & MIN_AGE <= 1)

nowOccNames_other <- nowOccNames %>%
  filter(MAX_AGE > 5 & MIN_AGE > 1)

# __Cleanup ---------------------------------------------------------------

rm(mdd.path,
   now.path,
   functions.path,
   nowOccNames)
   gc()

# _Harmonization ----------------------------------------------------------
# __MDD accepted species X NOW accepted species ---------------------------
# ___Set up default function arguments ------------------------------------

b1_suffix <- "mdd"
b1_col <- "MDD_species"
b2_suffix <- "now"
b2_col <- "now_accepted_name"
select_name_separator <- "_"
select_column_order <- c(
  "MDD_syn_ID",
  "SIDNUM",
  "MDD_species",
  "MDD_original_combination",
  "now_accepted_name",
  "string_distance",
  "match_notes",
  "MDD_family",
  "FAMILY",
  "MDD_continent",
  "now_continent",
  "MIN_AGE",
  "MAX_AGE"
)

# ___Max age up to 2 ma ---------------------------------------------------
# ____Exact ---------------------------------------------------------------

exact_valid_0to2max <- harmonize_exact_match(base1_dtf = mddTaxNames_acceptedSpecies,
                                             base2_dtf = nowOccNames_0to2max)
exact_valid_0to2max$exact_summary

exact_valid_0to2max_eval <- exact_valid_0to2max$exact_found %>%
  mutate(accepted_match = TRUE) %>%
  mutate(match_notes = paste0(match_notes, " (mdd accepted, now accepted)"))

evaluated_pairs <- exact_valid_0to2max_eval
finalSynonymy <- add_to_synonymy(data.frame(),
                                 exact_valid_0to2max_eval)

nowOccNames_unmatched_0to2max <- nowOccNames_0to2max %>%
  anti_join(filter(evaluated_pairs, accepted_match == TRUE))

rm(exact_valid_0to2max,
   exact_valid_0to2max_eval)
gc()

# ____Fuzzy ---------------------------------------------------------------

fuzzy_valid_0to2max <- harmonize_fuzzy_match(base1_dtf = mddTaxNames_acceptedSpecies,
                                             base2_dtf = nowOccNames_unmatched_0to2max,
                                             min_dist = 0, max_dist = 4)

fuzzy_valid_0to2max_eval <- fuzzy_valid_0to2max$min0max4_match_found %>%
  mutate(accepted_match = case_when(
    MDD_species == "Mustela_eversmanii" & now_accepted_name == "Mustela_eversmanni" & SIDNUM == "25444" ~ TRUE,
    TRUE ~ FALSE
  )) %>%
  mutate(match_notes = paste0(match_notes, " (mdd accepted, now accepted)"))

evaluated_pairs <- distinct(bind_rows(evaluated_pairs, fuzzy_valid_0to2max_eval))
finalSynonymy <- add_to_synonymy(finalSynonymy,
                                 fuzzy_valid_0to2max_eval)

nowOccNames_unmatched_0to2max <- nowOccNames_0to2max %>%
  anti_join(filter(evaluated_pairs, accepted_match == TRUE))

rm(nowOccNames_0to2max,
   fuzzy_valid_0to2max,
   fuzzy_valid_0to2max_eval)
gc()

# ___Max age up to 5 ma ---------------------------------------------------
# ____Exact ---------------------------------------------------------------

exact_valid_2to5max <- harmonize_exact_match(base1_dtf = mddTaxNames_acceptedSpecies,
                                             base2_dtf = nowOccNames_2to5max)
exact_valid_2to5max$exact_summary

exact_valid_2to5max_eval <- exact_valid_2to5max$exact_found %>%
  mutate(accepted_match = TRUE) %>%
  mutate(match_notes = paste0(match_notes, " (mdd accepted, now accepted)"))

evaluated_pairs <- distinct(bind_rows(evaluated_pairs, exact_valid_2to5max_eval))
finalSynonymy <- add_to_synonymy(finalSynonymy,
                                 exact_valid_2to5max_eval)

nowOccNames_unmatched_2to5max <- nowOccNames_2to5max %>%
  anti_join(filter(evaluated_pairs, accepted_match == TRUE))

rm(exact_valid_2to5max,
   exact_valid_2to5max_eval)
gc()

# ____Fuzzy ---------------------------------------------------------------

fuzzy_valid_2to5max <- harmonize_fuzzy_match(base1_dtf = mddTaxNames_acceptedSpecies,
                                             base2_dtf = nowOccNames_unmatched_2to5max,
                                             min_dist = 0, max_dist = 4)

fuzzy_valid_2to5max_eval <- fuzzy_valid_2to5max$min0max4_match_found %>%
  mutate(accepted_match = case_when(
    MDD_syn_ID == "100003947" & SIDNUM == "29262" ~ TRUE,
    TRUE ~ FALSE
  )) %>%
  mutate(match_notes = paste0(match_notes, " (mdd accepted, now accepted)"))

evaluated_pairs <- distinct(bind_rows(evaluated_pairs, fuzzy_valid_2to5max_eval))
finalSynonymy <- add_to_synonymy(finalSynonymy,
                                 fuzzy_valid_2to5max_eval)

nowOccNames_unmatched_2to5max <- nowOccNames_2to5max %>%
  anti_join(filter(evaluated_pairs, accepted_match == TRUE))

rm(nowOccNames_2to5max,
   fuzzy_valid_2to5max,
   fuzzy_valid_2to5max_eval)
gc()


# ___Min age up to 1 ma ---------------------------------------------------
# ____Exact ---------------------------------------------------------------

exact_valid_0to1min <- harmonize_exact_match(base1_dtf = mddTaxNames_acceptedSpecies,
                                             base2_dtf = nowOccNames_0to1min)
exact_valid_0to1min$exact_summary

exact_valid_0to1min_eval <- exact_valid_0to1min$exact_found %>%
  mutate(accepted_match = TRUE) %>%
  mutate(match_notes = paste0(match_notes, " (mdd accepted, now accepted)"))

evaluated_pairs <- distinct(bind_rows(evaluated_pairs, exact_valid_0to1min_eval))
finalSynonymy <- add_to_synonymy(finalSynonymy,
                                 exact_valid_0to1min_eval)

nowOccNames_unmatched_0to1min <- nowOccNames_0to1min %>%
  anti_join(filter(evaluated_pairs, accepted_match == TRUE))

rm(exact_valid_0to1min,
   exact_valid_0to1min_eval)
gc()

# ____Fuzzy ---------------------------------------------------------------

fuzzy_valid_0to1min <- harmonize_fuzzy_match(base1_dtf = mddTaxNames_acceptedSpecies,
                                             base2_dtf = nowOccNames_unmatched_0to1min,
                                             min_dist = 0, max_dist = 4)

fuzzy_valid_0to1min_eval <- fuzzy_valid_0to1min$min0max4_match_found %>%
  mutate(accepted_match = case_when(
    MDD_syn_ID == "100003754" & SIDNUM == "28792" ~ TRUE,
    TRUE ~ FALSE
  )) %>%
  mutate(match_notes = paste0(match_notes, " (mdd accepted, now accepted)"))

evaluated_pairs <- distinct(bind_rows(evaluated_pairs, fuzzy_valid_0to1min_eval))
finalSynonymy <- add_to_synonymy(finalSynonymy,
                                 fuzzy_valid_0to1min_eval)

nowOccNames_unmatched_0to1min <- nowOccNames_0to1min %>%
  anti_join(filter(evaluated_pairs, accepted_match == TRUE))

rm(nowOccNames_0to1min,
   fuzzy_valid_0to1min,
   fuzzy_valid_0to1min_eval)
gc()


# ___Other ages -----------------------------------------------------------
# ____Exact ---------------------------------------------------------------

exact_valid_other <- harmonize_exact_match(base1_dtf = mddTaxNames_acceptedSpecies,
                                             base2_dtf = nowOccNames_other)
exact_valid_other$exact_summary

exact_valid_other_eval <- exact_valid_other$exact_found %>%
  mutate(accepted_match = case_when(
    MDD_syn_ID == "100006414" & SIDNUM == "24488" ~ FALSE, #Age range incompatible with an extant species
    MDD_syn_ID == "100006096" & SIDNUM == "84395" ~ TRUE,
    TRUE ~ FALSE
  )) %>%
  mutate(match_notes = paste0(match_notes, " (mdd accepted, now accepted)"))

evaluated_pairs <- distinct(bind_rows(evaluated_pairs, exact_valid_other_eval))
finalSynonymy <- add_to_synonymy(finalSynonymy,
                                 exact_valid_other_eval)

nowOccNames_unmatched_other <- nowOccNames_other %>%
  anti_join(filter(evaluated_pairs, accepted_match == TRUE))

rm(exact_valid_other,
   exact_valid_other_eval)
gc()

# ____Fuzzy ---------------------------------------------------------------

fuzzy_valid_other <- harmonize_fuzzy_match(base1_dtf = mddTaxNames_acceptedSpecies,
                                             base2_dtf = nowOccNames_unmatched_other,
                                             min_dist = 0, max_dist = 4)

fuzzy_valid_other_eval <- fuzzy_valid_other$min0max4_match_found %>%
  mutate(accepted_match = FALSE) %>%
  mutate(match_notes = paste0(match_notes, " (mdd accepted, now accepted)"))

evaluated_pairs <- distinct(bind_rows(evaluated_pairs, fuzzy_valid_other_eval))
finalSynonymy <- add_to_synonymy(finalSynonymy,
                                 fuzzy_valid_other_eval)

nowOccNames_unmatched_other <- nowOccNames_other %>%
  anti_join(filter(evaluated_pairs, accepted_match == TRUE))

rm(nowOccNames_other,
   mddTaxNames_acceptedSpecies,
   fuzzy_valid_other,
   fuzzy_valid_other_eval)

# __MDD synonym species X NOW accepted species ----------------------------
# ___Set up default function arguments ------------------------------------

b1_suffix <- "mdd"
b1_col <- "MDD_original_combination"
b2_suffix <- "now"
b2_col <- "now_accepted_name"

# ___Max age up to 2 ma ---------------------------------------------------
# ____Exact ---------------------------------------------------------------

exact_syn_0to2max <- harmonize_exact_match(base1_dtf = mddTaxNames_species,
                                           base2_dtf = nowOccNames_unmatched_0to2max)
exact_syn_0to2max$exact_summary

exact_syn_0to2max_eval <- exact_syn_0to2max$exact_found %>%
  mutate(accepted_match = case_when(
    MDD_syn_ID == "100033238" & SIDNUM == "33077" ~ FALSE, #MDD synonym refers to domestic dog, and NOW occurrence is between 10k and 2ma. It's unclear whether this is a match.
    MDD_syn_ID == "100003738" & SIDNUM == "85148" ~ FALSE, #NOW comments on this occurrence calls it "Canis anthus primaevus" and MDD synonym refers to Cuon alpinus
    TRUE ~ TRUE
  )) %>%
  mutate(match_notes = paste0(match_notes, " (mdd synonym, now accepted)"))

evaluated_pairs <- distinct(bind_rows(evaluated_pairs, exact_syn_0to2max_eval))
finalSynonymy <- add_to_synonymy(finalSynonymy,
                                 exact_syn_0to2max_eval)

nowOccNames_unmatched_0to2max <- nowOccNames_unmatched_0to2max %>%
  anti_join(filter(evaluated_pairs, accepted_match == TRUE))

rm(exact_syn_0to2max,
   exact_syn_0to2max_eval)
gc()

# ____Fuzzy ---------------------------------------------------------------

fuzzy_syn_0to2max <- harmonize_fuzzy_match(base1_dtf = mddTaxNames_species,
                                           base2_dtf = nowOccNames_unmatched_0to2max,
                                           min_dist = 1, max_dist = 4)

fuzzy_syn_0to2max_eval <- fuzzy_syn_0to2max$min1max4_match_found %>%
  mutate(accepted_match = case_when(
    MDD_syn_ID == "100003743" & SIDNUM == "30020" ~ TRUE,
    MDD_syn_ID == "100034881" & SIDNUM == "85096" ~ TRUE,
    MDD_syn_ID == "100003743" & SIDNUM == "29615" ~ TRUE,
    MDD_syn_ID == "100003743" & SIDNUM == "29815" ~ TRUE,
    TRUE ~ FALSE
  )) %>%
  mutate(match_notes = paste0(match_notes, " (mdd synonym, now accepted)"))

evaluated_pairs <- distinct(bind_rows(evaluated_pairs, fuzzy_syn_0to2max_eval))
finalSynonymy <- add_to_synonymy(finalSynonymy,
                                 fuzzy_syn_0to2max_eval)

nowOccNames_unmatched_0to2max <- nowOccNames_unmatched_0to2max %>%
  anti_join(filter(evaluated_pairs, accepted_match == TRUE))

rm(fuzzy_syn_0to2max,
   fuzzy_syn_0to2max_eval)
gc()


# ___Max age up to 5 ma ---------------------------------------------------
# ____Exact ---------------------------------------------------------------

exact_syn_2to5max <- harmonize_exact_match(base1_dtf = mddTaxNames_species,
                                           base2_dtf = nowOccNames_unmatched_2to5max)
exact_syn_2to5max$exact_summary

exact_syn_2to5max_eval <- exact_syn_2to5max$exact_found %>%
  mutate(accepted_match = case_when(
    MDD_syn_ID == "100003634" & SIDNUM == "23952" ~ FALSE, #NOW comments on this occurrence calls it "Canis lupus variabilis" (described by Pei Wenzhong, 1934), not the same as the MDD Canis variabilis Wied-Neuwied, 1841
    MDD_syn_ID == "100003610" & SIDNUM == "82073" ~ FALSE, #NOW locality is Arizona, which likely means it refers to Canis thooides Tedford 2009, not the same as the MDD synonym Canis thooides Hilzheimer, 1906
    TRUE ~ TRUE
  )) %>%
  mutate(match_notes = paste0(match_notes, " (mdd synonym, now accepted)"))

evaluated_pairs <- distinct(bind_rows(evaluated_pairs, exact_syn_2to5max_eval))
finalSynonymy <- add_to_synonymy(finalSynonymy,
                                 exact_syn_2to5max_eval)

nowOccNames_unmatched_2to5max <- nowOccNames_unmatched_2to5max %>%
  anti_join(filter(evaluated_pairs, accepted_match == TRUE))

rm(exact_syn_2to5max,
   exact_syn_2to5max_eval)
gc()

# ____Fuzzy ---------------------------------------------------------------

fuzzy_syn_2to5max <- harmonize_fuzzy_match(base1_dtf = mddTaxNames_species,
                                           base2_dtf = nowOccNames_unmatched_2to5max,
                                           min_dist = 1, max_dist = 4)

fuzzy_syn_2to5max_eval <- fuzzy_syn_2to5max$min1max4_match_found %>%
  mutate(accepted_match = case_when(
    MDD_syn_ID == "100050545" & SIDNUM == "25319" ~ TRUE,
    TRUE ~ FALSE
  )) %>%
  mutate(match_notes = paste0(match_notes, " (mdd synonym, now accepted)"))

evaluated_pairs <- distinct(bind_rows(evaluated_pairs, fuzzy_syn_2to5max_eval))
finalSynonymy <- add_to_synonymy(finalSynonymy,
                                 fuzzy_syn_2to5max_eval)

nowOccNames_unmatched_2to5max <- nowOccNames_unmatched_2to5max %>%
  anti_join(filter(evaluated_pairs, accepted_match == TRUE))

rm(fuzzy_syn_2to5max,
   fuzzy_syn_2to5max_eval)
gc()





# ___Min age up 1ma -------------------------------------------------------
# ____Exact ---------------------------------------------------------------

exact_syn_0to1min <- harmonize_exact_match(base1_dtf = mddTaxNames_species,
                                           base2_dtf = nowOccNames_unmatched_0to1min)
exact_syn_0to1min$exact_summary

exact_syn_0to1min_eval <- exact_syn_0to1min$exact_found %>%
  mutate(accepted_match = TRUE) %>%
  mutate(match_notes = paste0(match_notes, " (mdd synonym, now accepted)"))

evaluated_pairs <- distinct(bind_rows(evaluated_pairs, exact_syn_0to1min_eval))
finalSynonymy <- add_to_synonymy(finalSynonymy,
                                 exact_syn_0to1min_eval)

nowOccNames_unmatched_0to1min <- nowOccNames_unmatched_0to1min %>%
  anti_join(filter(evaluated_pairs, accepted_match == TRUE))

rm(exact_syn_0to1min,
   exact_syn_0to1min_eval)
gc()

# ____Fuzzy ---------------------------------------------------------------

fuzzy_syn_0to1min <- harmonize_fuzzy_match(base1_dtf = mddTaxNames_species,
                                           base2_dtf = nowOccNames_unmatched_0to1min,
                                           min_dist = 1, max_dist = 4)

fuzzy_syn_0to1min_eval <- fuzzy_syn_0to1min$min1max4_match_found %>%
  mutate(accepted_match = FALSE) %>%
  mutate(match_notes = paste0(match_notes, " (mdd synonym, now accepted)"))

evaluated_pairs <- distinct(bind_rows(evaluated_pairs, fuzzy_syn_0to1min_eval))
finalSynonymy <- add_to_synonymy(finalSynonymy,
                                 fuzzy_syn_0to1min_eval)

nowOccNames_unmatched_0to1min <- nowOccNames_unmatched_0to1min %>%
  anti_join(filter(evaluated_pairs, accepted_match == TRUE))

rm(fuzzy_syn_0to1min,
   fuzzy_syn_0to1min_eval)
gc()




# ___Other ages -----------------------------------------------------------
# ____Exact ---------------------------------------------------------------

exact_syn_other <- harmonize_exact_match(base1_dtf = mddTaxNames_species,
                                         base2_dtf = nowOccNames_unmatched_other)
exact_syn_other$exact_summary

exact_syn_other_eval <- exact_syn_other$exact_found %>%
  mutate(accepted_match = case_when(
    MDD_syn_ID == "100034868" & SIDNUM == "24488" ~ FALSE,
    TRUE ~ TRUE
  )) %>%
  mutate(match_notes = paste0(match_notes, " (mdd synonym, now accepted)"))

evaluated_pairs <- distinct(bind_rows(evaluated_pairs, exact_syn_other_eval))
finalSynonymy <- add_to_synonymy(finalSynonymy,
                                 exact_syn_other_eval)

nowOccNames_unmatched_other <- nowOccNames_unmatched_other %>%
  anti_join(filter(evaluated_pairs, accepted_match == TRUE))

rm(exact_syn_other,
   exact_syn_other_eval)
gc()

# ____Fuzzy ---------------------------------------------------------------

fuzzy_syn_other <- harmonize_fuzzy_match(base1_dtf = mddTaxNames_species,
                                         base2_dtf = nowOccNames_unmatched_other,
                                         min_dist = 1, max_dist = 4)

fuzzy_syn_other_eval <- fuzzy_syn_other$min1max4_match_found %>%
  mutate(accepted_match = FALSE) %>%
  mutate(match_notes = paste0(match_notes, " (mdd synonym, now accepted)"))

evaluated_pairs <- distinct(bind_rows(evaluated_pairs, fuzzy_syn_other_eval))
finalSynonymy <- add_to_synonymy(finalSynonymy,
                                 fuzzy_syn_other_eval)

nowOccNames_unmatched_other <- nowOccNames_unmatched_other %>%
  anti_join(filter(evaluated_pairs, accepted_match == TRUE))

rm(fuzzy_syn_other,
   fuzzy_syn_other_eval)
gc()







a <- finalSynonymy %>%
  #filter(accepted_match == TRUE) %>%
  group_by(now_accepted_name, MDD_species) %>%
  summarise()



