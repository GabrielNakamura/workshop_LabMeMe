# _Info -------------------------------------------------------------------
##
## Title: harmonization_pbdb-mdd.R
## Purpose:
##
## Author: Thaís G. P. Faria
## Github: https://github.com/Thais-Faria
## Date created: 2025-06-01
## Copyright (c) Thaís G. P. Faria, 2025
## License: GNU General Public License version 3
##
# _Set up -----------------------------------------------------------------
# __Packages --------------------------------------------------------------

library(here) #Operating system agnostic file paths with .Rproj
library(stringr) #String manipulation
library(dplyr) #Data wrangling with tidyverse syntax
library(tidytable) #Brings the data.table package's fast operations to dplyr syntax
library(fuzzyjoin) #Merging data frames with fuzzy matching
library(tidylog) #Optional, prints detailed information about changes in the dataframe during data wrangling with tidyverse functions

# __File paths ------------------------------------------------------------

pbdbTaxNames.path <- here("data",
                          "processed",
                          "pbdb_taxa_species-names_2025-06-01.csv")

mddTaxNames.path <- here("data",
                         "processed",
                         "mdd_taxa_species-names_v2.1.csv")

evaluated.path <- here("data",
                       "processed",
                       "mdd_pbdb_evaluated-pairs.csv")

functions.path <- here("scripts",
                       "harmonization_source-functions.R")

# __Loading functions -----------------------------------------------------

source(functions.path)

# __Loading data ----------------------------------------------------------

mddTaxNames <- read.csv(mddTaxNames.path)

pbdbTaxNames <- read.csv(pbdbTaxNames.path)

evaluated_pairs <- read.csv(evaluated.path)

# __Sub-setting data ------------------------------------------------------

#PBDB

pbdbTaxNames.allSpecies <- pbdbTaxNames %>%
  select(pbdb_accepted_no,
         pbdb_taxon_no,
         pbdb_family,
         pbdb_accepted_name,
         pbdb_accepted_author_name,
         pbdb_accepted_author_year,
         pbdb_accepted_author_parentheses)

pbdbTaxNames.acceptedSpecies <- pbdbTaxNames %>%
  filter(pbdb_taxon_no == pbdb_accepted_no) %>%
  select(pbdb_accepted_no,
         pbdb_taxon_no,
         pbdb_family,
         pbdb_accepted_name,
         pbdb_accepted_author_name,
         pbdb_accepted_author_year,
         pbdb_accepted_author_parentheses) %>%
  distinct()

#MDD

mddTaxNames.allspecies <- mddTaxNames %>%
  filter(MDD_actual_rank == "species") %>%
  distinct()

mddTaxName.acceptedSpecies <- mddTaxNames.allspecies %>%
  filter(MDD_validity == "species") %>%
  distinct()

# Evaluated pairs

evaluated.keys <- evaluated_pairs %>%
  select(MDD_syn_ID, pbdb_taxon_no, accepted_match) %>%
  distinct()

# __Cleanup ---------------------------------------------------------------

rm(mddTaxNames,
   pbdbTaxNames,
   functions.path,
   pbdbTaxNames.path,
   mddTaxNames.path)
gc()

# _Harmonization ----------------------------------------------------------
# __MDD accepted species X PBDB accepted species --------------------------
# ___Set up default function arguments ------------------------------------

b1_suffix <- "mdd"
b1_col <- "MDD_species"
b2_suffix <- "pbdb"
b2_col <- "pbdb_accepted_name"
select_name_separator <- "_"
select_column_order <- c(
  "MDD_syn_ID",
  "pbdb_taxon_no",
  "MDD_species",
  "MDD_root_name",
  "MDD_original_combination",
  "pbdb_accepted_name",
  "string_distance",
  "match_notes",
  "MDD_family",
  "pbdb_family",
  "MDD_continent"
)
select.base1_id_col <- "MDD_syn_ID"
select.base2_id_col <- "pbdb_taxon_no"

# ___Exact ----------------------------------------------------------------

exact_valid <- harmonize_exact_match(base1_dtf = mddTaxName.acceptedSpecies,
                                     base2_dtf = pbdbTaxNames.acceptedSpecies)
exact_valid$exact_summary

exact_valid_eval <- exact_valid$exact_found %>%
  remove_already_evaluated() %>%
  mutate(accepted_match = case_when(
    MDD_author == "Boddaert" & pbdb_accepted_author_name == "Schreber" ~ FALSE,
    MDD_author == "Boddaert" & pbdb_accepted_author_name == "Gmelin" ~ FALSE,
    MDD_author == "Brongniart" & pbdb_accepted_author_name == "Cuvier" ~ FALSE,
    MDD_author == "A. G. Desmarest" & pbdb_accepted_author_name == "Lesson" ~ FALSE,
    MDD_author == "Griffith" & pbdb_accepted_author_name == "Smith" ~ FALSE,
    MDD_author == "I. Geoffroy Saint-Hilaire" & pbdb_accepted_author_name == "Cuvier" ~ FALSE,
    MDD_author == "Lesson" & pbdb_accepted_author_name == "Wood-Jones" ~ FALSE,
    MDD_author == "Matschie" & pbdb_accepted_author_name == "Hermann" ~ FALSE,
    TRUE ~ TRUE
  )) %>%
  mutate(match_notes = paste0(match_notes, " (mdd accepted, pbdb accepted)"))

evaluated_pairs <- bind_rows(evaluated_pairs, exact_valid_eval)

checkLater_pairs <- exact_valid_eval %>%
  inner_join(mddTaxName.acceptedSpecies) %>%
  filter(accepted_match == FALSE)

#finalSynonymy <- add_to_synonymy(data.frame(),
#                                 exact_valid_eval)

mddTaxNames_unmatched <- mddTaxNames.allspecies %>%
  anti_join(filter(exact_valid_eval, accepted_match == TRUE))

pbdbTaxNames_unmatched <- pbdbTaxNames.acceptedSpecies %>%
  anti_join(filter(exact_valid_eval, accepted_match == TRUE))

rm(exact_valid,
   exact_valid_eval,
   pbdbTaxNames.acceptedSpecies)
gc()


# ___Fuzzy ----------------------------------------------------------------
# ____min1max1 ------------------------------------------------------------

fuzzy_valid_min1max1 <- harmonize_fuzzy_match(base1_dtf = mddTaxName.acceptedSpecies,
                                              base2_dtf = pbdbTaxNames_unmatched,
                                              min_dist = 1, max_dist = 1)

fuzzy_valid_min1max1_eval <- fuzzy_valid_min1max1$min1max1_match_found %>%
  remove_already_evaluated() %>%
  distinct() %>%
  mutate(accepted_match = TRUE) %>%
  mutate(match_notes = paste0(match_notes, " (mdd accepted, pbdb accepted)"))

evaluated_pairs <- bind_rows(evaluated_pairs, fuzzy_valid_min1max1_eval) %>%
  distinct()

#finalSynonymy <- add_to_synonymy(finalSynonymy,
#                                 fuzzy_valid_min1max1_eval)

pbdbTaxNames_unmatched <- pbdbTaxNames_unmatched %>%
  anti_join(filter(fuzzy_valid_min1max1_eval, accepted_match == TRUE))

rm(fuzzy_valid_min1max1,
   fuzzy_valid_min1max1_eval)
gc()

# ____min2max2 ------------------------------------------------------------

fuzzy_valid_min2max2 <- harmonize_fuzzy_match(base1_dtf = mddTaxName.acceptedSpecies,
                                              base2_dtf = pbdbTaxNames_unmatched,
                                              min_dist = 2, max_dist = 2)

fuzzy_valid_min2max2_eval <- fuzzy_valid_min2max2$min2max2_match_found %>%
  remove_already_evaluated() %>%
  distinct() %>%
  mutate(accepted_match = case_when(
    MDD_syn_ID == "100004513" & pbdb_taxon_no == "224096" ~ FALSE,
    MDD_syn_ID == "100003617" & pbdb_taxon_no == "44879" ~ FALSE,
    MDD_syn_ID == "100004159" & pbdb_taxon_no == "224047" ~ FALSE,
    MDD_syn_ID == "100004185" & pbdb_taxon_no == "104159" ~ FALSE,
    MDD_syn_ID == "100005777" & pbdb_taxon_no == "157455" ~ FALSE,
    MDD_syn_ID == "100004629" & pbdb_taxon_no == "46524" ~ FALSE,
    MDD_syn_ID == "100003904" & pbdb_taxon_no == "235498" ~ FALSE,
    TRUE ~ TRUE
  )) %>%
  mutate(match_notes = paste0(match_notes, " (mdd accepted, pbdb accepted)"))

evaluated_pairs <- bind_rows(evaluated_pairs, fuzzy_valid_min2max2_eval) %>%
  distinct()

#finalSynonymy <- add_to_synonymy(finalSynonymy,
#                                 fuzzy_valid_min2max2_eval)

pbdbTaxNames_unmatched <- pbdbTaxNames_unmatched %>%
  anti_join(filter(fuzzy_valid_min2max2_eval, accepted_match == TRUE))

rm(fuzzy_valid_min2max2,
   fuzzy_valid_min2max2_eval)
gc()

# ____min3max3 ------------------------------------------------------------

fuzzy_valid_min3max3 <- harmonize_fuzzy_match(base1_dtf = mddTaxName.acceptedSpecies,
                                              base2_dtf = pbdbTaxNames_unmatched,
                                              min_dist = 3, max_dist = 3)

fuzzy_valid_min3max3_eval <- fuzzy_valid_min3max3$min3max3_match_found %>%
  remove_already_evaluated() %>%
  distinct() %>%
  mutate(accepted_match = case_when(
    MDD_syn_ID == "100006325" & pbdb_taxon_no == "237210" ~ TRUE,
    MDD_syn_ID == "100006203" & pbdb_taxon_no == "105014" ~ TRUE,
    TRUE ~ FALSE
  )) %>%
  mutate(match_notes = paste0(match_notes, " (mdd accepted, pbdb accepted)"))

evaluated_pairs <- bind_rows(evaluated_pairs, fuzzy_valid_min3max3_eval) %>%
  distinct()

#finalSynonymy <- add_to_synonymy(finalSynonymy,
#                                 fuzzy_valid_min3max3_eval)

pbdbTaxNames_unmatched <- pbdbTaxNames_unmatched %>%
  anti_join(filter(fuzzy_valid_min3max3_eval, accepted_match == TRUE))

rm(fuzzy_valid_min3max3,
   fuzzy_valid_min3max3_eval)
gc()

# ____min4max4 ------------------------------------------------------------

fuzzy_valid_min4max4 <- harmonize_fuzzy_match(base1_dtf = mddTaxName.acceptedSpecies,
                                              base2_dtf = pbdbTaxNames_unmatched,
                                              min_dist = 4, max_dist = 4)

fuzzy_valid_min4max4_eval <- fuzzy_valid_min4max4$min4max4_match_found %>%
  remove_already_evaluated() %>%
  distinct() %>%
  mutate(accepted_match = case_when(
    MDD_syn_ID == "100005246" & pbdb_taxon_no == "233625" ~ TRUE,
    MDD_syn_ID == "100003957" & pbdb_taxon_no == "365928" ~ TRUE,
    TRUE ~ FALSE
  )) %>%
  mutate(match_notes = paste0(match_notes, " (mdd accepted, pbdb accepted)"))

evaluated_pairs <- bind_rows(evaluated_pairs, fuzzy_valid_min4max4_eval) %>%
  distinct()

#finalSynonymy <- add_to_synonymy(finalSynonymy,
#                                 fuzzy_valid_min4max4_eval)

pbdbTaxNames_unmatched <- pbdbTaxNames_unmatched %>%
  anti_join(filter(fuzzy_valid_min4max4_eval, accepted_match == TRUE))

rm(fuzzy_valid_min4max4,
   fuzzy_valid_min4max4_eval)
gc()

# ____min0max4 ------------------------------------------------------------

fuzzy_valid_min0max4 <- harmonize_fuzzy_match(base1_dtf = mddTaxName.acceptedSpecies,
                                              base2_dtf = pbdbTaxNames_unmatched,
                                              min_dist = 0, max_dist = 4)

fuzzy_valid_min0max4_eval <- fuzzy_valid_min0max4$min0max4_match_found %>%
  remove_already_evaluated() %>%
  distinct() %>%
  mutate(accepted_match = case_when(
    MDD_syn_ID == "100005246" & pbdb_taxon_no == "233625" ~ TRUE,
    MDD_syn_ID == "100003957" & pbdb_taxon_no == "365928" ~ TRUE,
    TRUE ~ FALSE
  )) %>%
  mutate(match_notes = paste0(match_notes, " (mdd accepted, pbdb accepted)"))

evaluated_pairs <- bind_rows(evaluated_pairs, fuzzy_valid_min0max4_eval) %>%
  distinct()

#finalSynonymy <- add_to_synonymy(finalSynonymy,
#                                 fuzzy_valid_min4max4_eval)

pbdbTaxNames_unmatched <- pbdbTaxNames_unmatched %>%
  anti_join(filter(fuzzy_valid_min0max4_eval, accepted_match == TRUE))

rm(fuzzy_valid_min0max4,
   fuzzy_valid_min0max4_eval)
gc()

# __MDD synonym species X PBDB accepted species ---------------------------
# ___Set up default function arguments ------------------------------------

b1_suffix <- "mdd"
b1_col <- "MDD_original_combination"
b2_suffix <- "pbdb"
b2_col <- "pbdb_accepted_name"

# ___Exact ----------------------------------------------------------------

exact_mddSyns <- harmonize_exact_match(base1_dtf = mddTaxNames.allspecies,
                                       base2_dtf = pbdbTaxNames_unmatched)
exact_mddSyns$exact_summary

exact_mddSyns_eval <- exact_mddSyns$exact_found %>%
  remove_already_evaluated() %>%
  mutate(accepted_match = case_when(
    MDD_original_combination == "Canis_thooides" & pbdb_accepted_name == "Canis_thooides" ~ FALSE,
    MDD_original_combination == "Canis_variabilis" & pbdb_accepted_name == "Canis_variabilis" ~ FALSE,
    MDD_species == "Cerdocyon_thous" & MDD_original_combination == "Lycalopex_vetulus" & pbdb_accepted_name == "Lycalopex_vetulus" ~ FALSE,
    MDD_original_combination == "Hyaena_gigantea" & pbdb_accepted_name == "Hyaena_gigantea" & pbdb_accepted_author_name == "Schlosser" ~ FALSE,
    MDD_species == "Herpestes_sanguineus" & MDD_original_combination == "Herpestes_fuscus" & pbdb_accepted_name == "Herpestes_fuscus" ~ FALSE,
    MDD_species == "Leopardus_geoffroyi" & MDD_original_combination == "Felis_pardoides" & pbdb_accepted_name == "Felis_pardoides" ~ FALSE,
    MDD_species == "Nyctereutes_procyonoides" & MDD_original_combination == "Nyctereutes_sinensis" & pbdb_accepted_name == "Nyctereutes_sinensis" ~ FALSE,
    MDD_species == "Panthera_pardus" & MDD_original_combination == "Felis_prisca" & pbdb_accepted_name == "Felis_prisca" ~ FALSE,
    MDD_syn_ID == "100031829" & pbdb_taxon_no == "232911" ~ FALSE,
    MDD_syn_ID == "100029377" & pbdb_taxon_no == "232911" ~ FALSE,
    MDD_syn_ID == "100028187" & pbdb_taxon_no == "232911" ~ FALSE,
    TRUE ~ TRUE
  )) %>%
  mutate(match_notes = paste0(match_notes, " (mdd synonyms, pbdb accepted)")) %>%
  select(-MDD_order,
         -MDD_genus,
         -MDD_actual_rank,
         -MDD_validity)

evaluated_pairs <- bind_rows(evaluated_pairs, exact_mddSyns_eval) %>%
  distinct()

pbdbTaxNames_unmatched <- pbdbTaxNames_unmatched %>%
  anti_join(filter(evaluated_pairs, accepted_match == TRUE))

#finalSynonymy <- add_to_synonymy(finalSynonymy, exact_mddSyns_eval)

rm(exact_mddSyns_eval)

gc()
# ___Fuzzy ----------------------------------------------------------------
# ____min1max1 ------------------------------------------------------------

fuzzy_min1max1 <- harmonize_fuzzy_match(base1_dtf = mddTaxNames.allspecies,
                                        base2_dtf = pbdbTaxNames_unmatched,
                                        min_dist = 1, max_dist = 1)

fuzzy_min1max1_eval <- fuzzy_min1max1$min1max1_match_found %>%
  remove_already_evaluated() %>%
  mutate(accepted_match = case_when(
    MDD_species == "Vulpes_chama" & pbdb_accepted_name == "Canis_cana" ~ FALSE,
    MDD_syn_ID == "100029968" & pbdb_taxon_no == "104168" ~ FALSE,
    TRUE ~ TRUE
  ),
  match_notes = paste0(match_notes, " (mdd synonyms, pbdb accepted)"))

evaluated_pairs <- bind_rows(evaluated_pairs,
                             fuzzy_min1max1_eval) %>%
  distinct()

pbdbTaxNames_unmatched <- pbdbTaxNames_unmatched %>%
  anti_join(filter(evaluated_pairs, accepted_match == TRUE))

#finalSynonymy <- add_to_synonymy(finalSynonymy, fuzzy_min1max1_eval)

rm(exact_mddSyns,
   fuzzy_min1max1_eval,
   fuzzy_min1max1)
gc()


# ____min2max2 ------------------------------------------------------------

fuzzy_min2max2 <- harmonize_fuzzy_match(base1_dtf = mddTaxNames.allspecies,
                                        base2_dtf = pbdbTaxNames_unmatched,
                                        min_dist = 2, max_dist = 2)

fuzzy_min2max2_eval <- fuzzy_min2max2$min2max2_match_found %>%
  remove_already_evaluated() %>%
  mutate(accepted_match = case_when(
    MDD_original_combination == "Arctocephalus_tasmanicus" & pbdb_accepted_name == "Arctocephalus_tasmanica" ~ TRUE,
    MDD_original_combination == "Profelis_temminckii" & pbdb_accepted_name == "Pardofelis_temminckii" ~ TRUE,
    MDD_original_combination == "Chaeffia_adusta" & pbdb_accepted_name == "Schaeffia_adusta" ~ TRUE,
    MDD_original_combination == "Pagophilus_groenlandicus" & pbdb_accepted_name == "Pagophilus_groenlandica" ~ TRUE,
    MDD_original_combination == "Poecilictis_libyca" & pbdb_accepted_name == "Poecilictis_lybica" ~ TRUE,
    MDD_syn_ID == "100031118" & pbdb_taxon_no == "231845" ~ TRUE,
    MDD_syn_ID == "100053318" & pbdb_taxon_no == "232925" ~ TRUE,
    MDD_syn_ID == "100029959" & pbdb_taxon_no == "471809" ~ TRUE,
    MDD_syn_ID == "100031552" & pbdb_taxon_no == "433170" ~ TRUE,
    MDD_syn_ID == "100033591" & pbdb_taxon_no == "233974" ~ TRUE,
    MDD_syn_ID == "100046105" & pbdb_taxon_no == "224043" ~ TRUE,
    MDD_syn_ID == "100028590" & pbdb_taxon_no == "224043" ~ TRUE,
    MDD_syn_ID == "100032235" & pbdb_taxon_no == "224043" ~ TRUE,
    MDD_syn_ID == "100048291" & pbdb_taxon_no == "231413" ~ TRUE,
    MDD_syn_ID == "100003696" & pbdb_taxon_no == "433177" ~ TRUE,
    MDD_syn_ID == "100058800" & pbdb_taxon_no == "148879" ~ TRUE,
    TRUE ~ FALSE
  ),
  match_notes = paste0(match_notes, " (mdd synonyms, pbdb accepted)"))

evaluated_pairs <- distinct(bind_rows(evaluated_pairs, fuzzy_min2max2_eval))

pbdbTaxNames_unmatched <- pbdbTaxNames_unmatched %>%
  anti_join(filter(evaluated_pairs, accepted_match == TRUE))

#finalSynonymy <- add_to_synonymy(finalSynonymy, fuzzy_min2max2_eval)

rm(fuzzy_min2max2, fuzzy_min2max2_eval)
gc()

# ____min3max3 ------------------------------------------------------------

fuzzy_min3max3 <- harmonize_fuzzy_match(base1_dtf = mddTaxNames.allspecies,
                                        base2_dtf = pbdbTaxNames_unmatched,
                                        min_dist = 3, max_dist = 3)

fuzzy_min3max3_eval <- fuzzy_min3max3$min3max3_match_found %>%
  remove_already_evaluated() %>%
  anti_join(evaluated_pairs, by = c("pbdb_accepted_name" = "pbdb_accepted_name")) %>%
  mutate(accepted_match = case_when(
    MDD_family != pbdb_family ~ FALSE
  ), accepted_match = case_when(
    MDD_original_combination == "Pusa_hispida" & pbdb_accepted_name == "Phoca_hispida" ~ TRUE,
    MDD_syn_ID == "100053429" & pbdb_taxon_no == "224062" ~ TRUE,
    MDD_syn_ID == "100046071" & pbdb_taxon_no == "48139" ~ TRUE,
    MDD_syn_ID == "100048565" & pbdb_taxon_no == "157458" ~ TRUE,
    MDD_syn_ID == "100039307" & pbdb_taxon_no == "98693" ~ TRUE,
    MDD_syn_ID == "100006600" & pbdb_taxon_no == "52609" ~ TRUE,
    MDD_syn_ID == "100034865" & pbdb_taxon_no == "81237" ~ TRUE,
    TRUE ~ FALSE
  ),
  match_notes = paste0(match_notes, " (mdd synonyms, pbdb accepted)"))

#finalSynonymy <- add_to_synonymy(finalSynonymy, fuzzy_min3max3_eval)
evaluated_pairs <- distinct(bind_rows(evaluated_pairs, fuzzy_min3max3_eval))

pbdbTaxNames_unmatched <- pbdbTaxNames_unmatched %>%
  anti_join(filter(evaluated_pairs, accepted_match == TRUE))

rm(fuzzy_min3max3, fuzzy_min3max3_eval)
gc()

# ____min4max4 ------------------------------------------------------------

fuzzy_min4max4 <- harmonize_fuzzy_match(base1_dtf = mddTaxNames.allspecies,
                                        base2_dtf = pbdbTaxNames_unmatched,
                                        min_dist = 4, max_dist = 4)

fuzzy_min4max4_eval <- fuzzy_min4max4$min4max4_match_found %>%
  remove_already_evaluated() %>%
  mutate(accepted_match = case_when(
    MDD_original_combination == "Aeluropus_fovealis" & pbdb_accepted_name == "Ailuropoda_fovealis" ~ TRUE,
    MDD_syn_ID == "100058195" & pbdb_taxon_no == "224062" ~ TRUE,
    MDD_syn_ID == "100046080" & pbdb_taxon_no == "233976" ~ TRUE,
    TRUE ~ FALSE
  ),
  match_notes = paste0(match_notes, " (mdd synonyms, pbdb accepted)"))


#finalSynonymy <- add_to_synonymy(finalSynonymy, fuzzy_min4max4_eval)
evaluated_pairs <- distinct(bind_rows(evaluated_pairs, fuzzy_min4max4_eval))

pbdbTaxNames_unmatched <- pbdbTaxNames_unmatched %>%
  anti_join(filter(evaluated_pairs, accepted_match == TRUE))

rm(fuzzy_min4max4,
   fuzzy_min4max4_eval)
gc()

# ____min0max4 ------------------------------------------------------------

fuzzy_min0max4 <- harmonize_fuzzy_match(base1_dtf = mddTaxNames.allspecies,
                                        base2_dtf = pbdbTaxNames_unmatched,
                                        min_dist = 0, max_dist = 4)

fuzzy_min0max4_eval <- fuzzy_min0max4$min0max4_match_found %>%
  remove_already_evaluated() %>%
  mutate(accepted_match = case_when(
    MDD_original_combination == "Aeluropus_fovealis" & pbdb_accepted_name == "Ailuropoda_fovealis" ~ TRUE,
    MDD_syn_ID == "100058195" & pbdb_taxon_no == "224062" ~ TRUE,
    MDD_syn_ID == "100046080" & pbdb_taxon_no == "233976" ~ TRUE,
    TRUE ~ FALSE
  ),
  match_notes = paste0(match_notes, " (mdd synonyms, pbdb accepted)"))


#finalSynonymy <- add_to_synonymy(finalSynonymy, fuzzy_min4max4_eval)
evaluated_pairs <- distinct(bind_rows(evaluated_pairs, fuzzy_min0max4_eval))

pbdbTaxNames_unmatched <- pbdbTaxNames_unmatched %>%
  anti_join(filter(evaluated_pairs, accepted_match == TRUE))

rm(fuzzy_min0max4,
   fuzzy_min0max4_eval)
gc()










#### Fuzzy part 2 ####

fuzzy_min1max4 <- harmonize_fuzzy_match(base1_dtf = mddTaxNames.allspecies,
                                        base2_dtf = fuzzy_min4max4$pbdb_min4max4_match_failed,
                                        min_dist = 1, max_dist = 4)

fuzzy_min1max4_eval <- fuzzy_min1max4$min1max4_match_found %>%
  mutate(accepted_match = FALSE,
         match_notes = paste0(match_notes, " (mdd synonyms, pbdb accepted)"))


evaluated_pairs <- distinct(bind_rows(evaluated_pairs, fuzzy_min1max4_eval))
finalSynonymy <- add_to_synonymy(finalSynonymy, fuzzy_min1max4_eval)

pbdbTaxNames_unmatched <- fuzzy_min1max4$pbdb_min1max4_match_failed %>%
  inner_join(distinct(select(pbdbTaxNames, pbdb_accepted_name, pbdb_alternative_synonym)))

rm(fuzzy_min4max4,
   fuzzy_min1max4,
   fuzzy_min1max4_eval)
gc()

#### Extra ####

b1_suffix <- "mdd"
b1_col <- "MDD_original_combination"
b2_suffix <- "pbdb"
b2_col <- "pbdb_alternative_synonym"

pbdb_alternativeNames <- create_alternative_with_subgenus(data = pbdbTaxNames_unmatched,
                                                          name_col = "pbdb_accepted_name") %>%
  rbind(., create_alternative_with_subgenus(data = pbdbTaxNames_unmatched,
                                            name_col = "pbdb_alternative_synonym")) %>%
  select(-pbdb_alternative_synonym) %>%
  distinct() %>%
  rename(pbdb_alternative_synonym = del_or_swap_subgen) %>%
  rbind(filter(pbdbTaxNames_unmatched, !str_detect(pbdb_accepted_name, pattern = "\\("))) %>%
  filter(pbdb_accepted_name != pbdb_alternative_synonym)

fuzzy_min0max4 <- harmonize_fuzzy_match(base1_dtf = mddTaxNames.allspecies,
                                        base2_dtf = pbdb_alternativeNames,
                                        min_dist = 0, max_dist = 4)

fuzzy_min0max4_eval <- fuzzy_min0max4$min0max4_match_found %>%
  mutate(accepted_match = case_when(
    MDD_family != pbdb_family ~ FALSE
  ), accepted_match = case_when(
    MDD_original_combination == "Lynx_montanus" & pbdb_alternative_synonym == "Lynx_montana" ~ TRUE,
    MDD_original_combination == "Felis_montana" & pbdb_alternative_synonym == "Lynx_montana" ~ TRUE,
    MDD_original_combination == "Morunga_elephantina" & pbdb_alternative_synonym == "Mirounga_elephantinus" ~ TRUE,
    MDD_original_combination == "Arctocephalus_gracilis" & pbdb_alternative_synonym == "Arctocephalus_gracilis" ~ TRUE,
    MDD_original_combination == "Arctocephalus_delalandii" & pbdb_alternative_synonym == "Arctocephalus_delalandii" ~ TRUE,
    MDD_original_combination == "Lynx_vulgaris" & pbdb_alternative_synonym == "Lynx_vulgaris" ~ TRUE,
    MDD_original_combination == "Ursus_optimus" & pbdb_alternative_synonym == "Ursus_optimus" ~ TRUE,
    TRUE ~ FALSE
  ), match_notes = paste0(match_notes, " (mdd synonyms, pbdb alternative)")) %>%
  select(-MDD_order,
         -MDD_genus,
         -MDD_actual_rank,
         -MDD_validity,
         -pbdb_alternative_synonym) %>%
  distinct()

evaluated_pairs <- distinct(bind_rows(evaluated_pairs, fuzzy_min0max4_eval))
finalSynonymy <- add_to_synonymy(finalSynonymy, fuzzy_min0max4_eval)

rm(pbdb_alternativeNames,
   fuzzy_min0max4,
   fuzzy_min0max4_eval)
gc()

##### Harmonization - MDD accepted names x PBDB synonymy #####
#### Setting up default arguments and data ####

b1_suffix <- "mdd"
b1_col <- "MDD_species"
b2_suffix <- "pbdb"
b2_col <- "pbdb_taxon_name"

pbdbTaxNames_unmatched <- pbdbTaxNames %>%
  filter(!pbdb_accepted_name %in% finalSynonymy$pbdb_accepted_name) %>%
  filter(!pbdb_taxon_name %in% evaluated_pairs$pbdb_accepted_name)

#### Exact + Fuzzy ####
fuzzy_min0max4 <- harmonize_fuzzy_match(base1_dtf = mddTaxName.acceptedSpecies,
                                        base2_dtf = pbdbTaxNames_unmatched,
                                        min_dist = 0, max_dist = 4)

fuzzy_min0max4_eval <- fuzzy_min0max4$min0max4_match_found %>%
  mutate(accepted_match = case_when(
    MDD_family != pbdb_family ~ FALSE
  ), accepted_match = case_when(
    MDD_species == "Genetta_poensis" & pbdb_taxon_name == "Genetta_poensis" ~ TRUE,
    MDD_species == "Viverricula_indica" & pbdb_taxon_name == "Viverricula_indica" ~ TRUE,
    TRUE ~ FALSE
  ), match_notes = paste0(match_notes, " (mdd accepted, pbdb synonyms)")) %>%
  distinct()

evaluated_pairs <- distinct(bind_rows(evaluated_pairs, fuzzy_min0max4_eval))
finalSynonymy <- add_to_synonymy(finalSynonymy, fuzzy_min0max4_eval)

rm(fuzzy_min0max4_eval)
gc()

#### Extra ####

b1_suffix <- "mdd"
b1_col <- "MDD_species"
b2_suffix <- "pbdb"
b2_col <- "pbdb_alternative_synonym"

pbdb_alternativeNames <- create_alternative_with_subgenus(data = fuzzy_min0max4$pbdb_min0max4_match_failed,
                                                          name_col = "pbdb_taxon_name") %>%
  rbind(., create_alternative_with_subgenus(data = fuzzy_min0max4$pbdb_min0max4_match_failed,
                                            name_col = "pbdb_alternative_synonym")) %>%
  select(-pbdb_alternative_synonym) %>%
  distinct() %>%
  rename(pbdb_alternative_synonym = del_or_swap_subgen) %>%
  rbind(filter(pbdbTaxNames_unmatched, !str_detect(pbdb_accepted_name, pattern = "\\("))) %>%
  filter(pbdb_accepted_name != pbdb_alternative_synonym)

fuzzy_min0max4 <- harmonize_fuzzy_match(base1_dtf = mddTaxName.acceptedSpecies,
                                        base2_dtf = pbdb_alternativeNames,
                                        min_dist = 0, max_dist = 4)

fuzzy_min0max4_eval <- fuzzy_min0max4$min0max4_match_found %>%
  mutate(accepted_match = FALSE,
         match_notes = paste0(match_notes, " (mdd synonyms, pbdb alternative")) %>%
  select(-pbdb_alternative_synonym) %>%
  distinct()

evaluated_pairs <- distinct(bind_rows(evaluated_pairs, fuzzy_min0max4_eval))
finalSynonymy <- add_to_synonymy(finalSynonymy, fuzzy_min0max4_eval)

write.csv(evaluated_pairs, here("data", "processed", "mdd_pbdb_evaluated-pairs.csv"))
