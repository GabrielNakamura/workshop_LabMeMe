##### Packages #####
library(here)
library(stringr)
library(dplyr)
library(tidytable)
library(fuzzyjoin)
library(tidylog)

##### Loading data #####
pbdbTaxNames <- read.csv(here("data",
                              "processed",
                              "pbdb_taxa_species-names_2025-05-06.csv"))

mddTaxNames <- read.csv(here("data",
                             "processed",
                             "mdd_taxa_species-names_v2.1.csv"))

##### Subsetting data #####
pbdbTaxNames_acceptedSpecies <- pbdbTaxNames %>%
  select(pbdb_family,
         pbdb_accepted_name,
         pbdb_accepted_author_name,
         pbdb_accepted_author_year,
         pbdb_accepted_author_parentheses) %>%
  distinct()

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
b2_suffix <- "pbdb"
b2_col <- "pbdb_accepted_name"
select_name_separator <- "_"

#### Merging ####

exact_valid <- harmonize_exact_match(base1_dtf = mddTaxNames_acceptedSpecies,
                                     base2_dtf = pbdbTaxNames_acceptedSpecies)
exact_valid$exact_summary

exact_valid_eval <- exact_valid$exact_found %>%
  inner_join(mddTaxNames_acceptedSpecies) %>%
  distinct() %>%
  inner_join(pbdbTaxNames_acceptedSpecies) %>%
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

evaluated_pairs <- exact_valid_eval

checkLater_pairs <- exact_valid_eval %>%
  inner_join(mddTaxNames_acceptedSpecies) %>%
  filter(accepted_match == FALSE)

finalSynonymy <- add_to_synonymy(data.frame(),
                                 exact_valid_eval)

mddTaxNames_unmatched <- mddTaxNames_species %>%
  anti_join(filter(exact_valid_eval, accepted_match == TRUE))

pbdbTaxNames_unmatched <- pbdbTaxNames_acceptedSpecies %>%
  anti_join(filter(exact_valid_eval, accepted_match == TRUE))

rm(exact_valid,
   exact_valid_eval,
   pbdbTaxNames_acceptedSpecies)
gc()

##### Harmonization - MDD synonymy x PBDB accepted names #####
#### Setting up default arguments ####

b1_suffix <- "mdd"
b1_col <- "MDD_original_combination"
b2_suffix <- "pbdb"
b2_col <- "pbdb_accepted_name"

#### Exact ####
exact_mddSyns <- harmonize_exact_match(base1_dtf = mddTaxNames_unmatched,
                                       base2_dtf = pbdbTaxNames_unmatched)
exact_mddSyns$exact_summary

exact_mddSyns_eval <- exact_mddSyns$exact_found %>%
  mutate(accepted_match = case_when(
    MDD_original_combination == "Canis_thooides" & pbdb_accepted_name == "Canis_thooides" ~ FALSE,
    MDD_original_combination == "Canis_variabilis" & pbdb_accepted_name == "Canis_variabilis" ~ FALSE,
    MDD_species == "Cerdocyon_thous" & MDD_original_combination == "Lycalopex_vetulus" & pbdb_accepted_name == "Lycalopex_vetulus" ~ FALSE,
    MDD_original_combination == "Hyaena_gigantea" & pbdb_accepted_name == "Hyaena_gigantea" & pbdb_accepted_author_name == "Schlosser" ~ FALSE,
    MDD_species == "Herpestes_sanguineus" & MDD_original_combination == "Herpestes_fuscus" & pbdb_accepted_name == "Herpestes_fuscus" ~ FALSE,
    MDD_species == "Leopardus_geoffroyi" & MDD_original_combination == "Felis_pardoides" & pbdb_accepted_name == "Felis_pardoides" ~ FALSE,
    MDD_species == "Nyctereutes_procyonoides" & MDD_original_combination == "Nyctereutes_sinensis" & pbdb_accepted_name == "Nyctereutes_sinensis" ~ FALSE,
    MDD_species == "Panthera_pardus" & MDD_original_combination == "Felis_prisca" & pbdb_accepted_name == "Felis_prisca" ~ FALSE,
    TRUE ~ TRUE
  )) %>%
  mutate(match_notes = paste0(match_notes, " (mdd synonyms, pbdb accepted)")) %>%
  select(-MDD_order,
         -MDD_genus,
         -MDD_actual_rank,
         -MDD_validity)

evaluated_pairs <- bind_rows(evaluated_pairs, exact_mddSyns_eval) %>%
  distinct()

finalSynonymy <- add_to_synonymy(finalSynonymy, exact_mddSyns_eval)

rm(exact_mddSyns_eval,
   mddTaxNames_unmatched,
   pbdbTaxNames_unmatched)

gc()

#### Fuzzy part 1 ####

fuzzy_min1max1 <- harmonize_fuzzy_match(base1_dtf = exact_mddSyns$mdd_exact_failed,
                                        base2_dtf = exact_mddSyns$pbdb_exact_failed,
                                        min_dist = 1, max_dist = 1)

fuzzy_min1max1_eval <- fuzzy_min1max1$min1max1_match_found %>%
  mutate(accepted_match = case_when(
    MDD_species == "Vulpes_chama" & pbdb_accepted_name == "Canis_cana" ~ FALSE,
    TRUE ~ TRUE
  ),
  match_notes = paste0(match_notes, " (mdd synonyms, pbdb accepted)")) %>%
  select(-MDD_order,
         -MDD_genus,
         -MDD_actual_rank,
         -MDD_validity)

evaluated_pairs <- bind_rows(evaluated_pairs,
                             fuzzy_min1max1_eval) %>%
  distinct()

finalSynonymy <- add_to_synonymy(finalSynonymy, fuzzy_min1max1_eval)

rm(exact_mddSyns,
   fuzzy_min1max1_eval)
gc()

fuzzy_min2max2 <- harmonize_fuzzy_match(base1_dtf = fuzzy_min1max1$mdd_min1max1_match_failed,
                                        base2_dtf = fuzzy_min1max1$pbdb_min1max1_match_failed,
                                        min_dist = 2, max_dist = 2)

fuzzy_min2max2_eval <- fuzzy_min2max2$min2max2_match_found %>%
  anti_join(evaluated_pairs, by = c("pbdb_accepted_name" = "pbdb_accepted_name")) %>%
  mutate(accepted_match = case_when(
    MDD_original_combination == "Arctocephalus_tasmanicus" & pbdb_accepted_name == "Arctocephalus_tasmanica" ~ TRUE,
    MDD_original_combination == "Profelis_temminckii" & pbdb_accepted_name == "Pardofelis_temminckii" ~ TRUE,
    MDD_original_combination == "Chaeffia_adusta" & pbdb_accepted_name == "Schaeffia_adusta" ~ TRUE,
    MDD_original_combination == "Pagophilus_groenlandicus" & pbdb_accepted_name == "Pagophilus_groenlandica" ~ TRUE,
    MDD_original_combination == "Poecilictis_libyca" & pbdb_accepted_name == "Poecilictis_lybica" ~ TRUE,
    TRUE ~ FALSE
  ),
  match_notes = paste0(match_notes, " (mdd synonyms, pbdb accepted)")) %>%
  select(-MDD_order,
         -MDD_genus,
         -MDD_actual_rank,
         -MDD_validity)

evaluated_pairs <- distinct(bind_rows(evaluated_pairs, fuzzy_min2max2_eval))
finalSynonymy <- add_to_synonymy(finalSynonymy, fuzzy_min2max2_eval)

rm(fuzzy_min1max1, fuzzy_min2max2_eval)
gc()

fuzzy_min3max3 <- harmonize_fuzzy_match(base1_dtf = fuzzy_min2max2$mdd_min2max2_match_failed,
                                        base2_dtf = fuzzy_min2max2$pbdb_min2max2_match_failed,
                                        min_dist = 3, max_dist = 3)

fuzzy_min3max3_eval <- fuzzy_min3max3$min3max3_match_found %>%
  anti_join(evaluated_pairs, by = c("pbdb_accepted_name" = "pbdb_accepted_name")) %>%
  mutate(accepted_match = case_when(
    MDD_family != pbdb_family ~ FALSE
  ), accepted_match = case_when(
    MDD_original_combination == "Pusa_hispida" & pbdb_accepted_name == "Phoca_hispida" ~ TRUE,
    TRUE ~ FALSE
  ),
  match_notes = paste0(match_notes, " (mdd synonyms, pbdb accepted)")) %>%
  select(-MDD_order,
         -MDD_genus,
         -MDD_actual_rank,
         -MDD_validity)

finalSynonymy <- add_to_synonymy(finalSynonymy, fuzzy_min3max3_eval)
evaluated_pairs <- distinct(bind_rows(evaluated_pairs, fuzzy_min3max3_eval))

rm(fuzzy_min2max2, fuzzy_min3max3_eval)
gc()

fuzzy_min4max4 <- harmonize_fuzzy_match(base1_dtf = fuzzy_min3max3$mdd_min3max3_match_failed,
                                        base2_dtf = fuzzy_min3max3$pbdb_min3max3_match_failed,
                                        min_dist = 4, max_dist = 4)

fuzzy_min4max4_eval <- fuzzy_min4max4$min4max4_match_found %>%
  anti_join(evaluated_pairs, by = c("pbdb_accepted_name" = "pbdb_accepted_name")) %>%
  mutate(accepted_match = case_when(
    MDD_family != pbdb_family ~ FALSE
  ), accepted_match = case_when(
    MDD_original_combination == "Aeluropus_fovealis" & pbdb_accepted_name == "Ailuropoda_fovealis" ~ TRUE,
    TRUE ~ FALSE
  ),
  match_notes = paste0(match_notes, " (mdd synonyms, pbdb accepted)")) %>%
  select(-MDD_order,
         -MDD_genus,
         -MDD_actual_rank,
         -MDD_validity)

finalSynonymy <- add_to_synonymy(finalSynonymy, fuzzy_min4max4_eval)
evaluated_pairs <- distinct(bind_rows(evaluated_pairs, fuzzy_min4max4_eval))

rm(fuzzy_min3max3,
   fuzzy_min4max4_eval)
gc()

#### Fuzzy part 2 ####

fuzzy_min1max4 <- harmonize_fuzzy_match(base1_dtf = mddTaxNames_species,
                                        base2_dtf = fuzzy_min4max4$pbdb_min4max4_match_failed,
                                        min_dist = 1, max_dist = 4)

fuzzy_min1max4_eval <- fuzzy_min1max4$min1max4_match_found %>%
  mutate(accepted_match = FALSE,
         match_notes = paste0(match_notes, " (mdd synonyms, pbdb accepted)")) %>%
  select(-MDD_order,
         -MDD_genus,
         -MDD_actual_rank,
         -MDD_validity)

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

fuzzy_min0max4 <- harmonize_fuzzy_match(base1_dtf = mddTaxNames_species,
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
fuzzy_min0max4 <- harmonize_fuzzy_match(base1_dtf = mddTaxNames_acceptedSpecies,
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

fuzzy_min0max4 <- harmonize_fuzzy_match(base1_dtf = mddTaxNames_acceptedSpecies,
                                        base2_dtf = pbdb_alternativeNames,
                                        min_dist = 0, max_dist = 4)

fuzzy_min0max4_eval <- fuzzy_min0max4$min0max4_match_found %>%
  mutate(accepted_match = FALSE,
         match_notes = paste0(match_notes, " (mdd synonyms, pbdb alternative")) %>%
  select(-pbdb_alternative_synonym) %>%
  distinct()

evaluated_pairs <- distinct(bind_rows(evaluated_pairs, fuzzy_min0max4_eval))
finalSynonymy <- add_to_synonymy(finalSynonymy, fuzzy_min0max4_eval)

