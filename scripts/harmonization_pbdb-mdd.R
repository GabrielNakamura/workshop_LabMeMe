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
                           "pbdb_taxa_species-names_2025-05-06.csv")) %>%
  rename(pbdb_accepted_name = accepted_name,
         pbdb_taxon_name = taxon_name,
         pbdb_family = family)

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

#### Merging ####

exact_valid <- harmonize_exact_match(base1_dtf = mddTaxNames_acceptedSpecies,
                                     base2_dtf = pbdbTaxNames_acceptedSpecies)
exact_valid$exact_summary

finalSynonymy_acceptedAll <- exact_valid$exact_found %>%
  inner_join(mddTaxNames_acceptedSpecies) %>%
  distinct() %>%
  inner_join(pbdbTaxNames_acceptedSpecies)

pbdbTaxNames_unmatched <- anti_join(pbdbTaxNames,
                                    exact_valid$exact_found) %>%
  distinct()

##### Harmonization - MDD synonymy x PBDB accepted names #####
#### Setting up default arguments ####

b1_suffix <- "mdd"
b1_col <- "MDD_original_combination"
b2_suffix <- "pbdb"
b2_col <- "pbdb_accepted_name"

#### Exact ####
exact_mddSyns <- harmonize_exact_match(base1_dtf = mddTaxNames_species,
                                       base2_dtf = exact_valid$pbdb_exact_failed)
exact_mddSyns$exact_summary

finalSynonymy <- inner_join(mddTaxNames_species,
                            exact_mddSyns$exact_found) %>%
  distinct() %>%
  inner_join(pbdbTaxNames_acceptedSpecies) %>%
  distinct()


finalSynonymy <- exact$exact_found

#### Fuzzy part 1 ####

fuzzy_min1max1 <- harmonize_fuzzy_match(base1_dtf = exact$mdd_exact_failed,
                                        base2_dtf = exact$pbdb_exact_failed,
                                        min_dist = 1, max_dist = 1, delim = "_")

fuzzy_min1max1_eval <- fuzzy_min1max1$min1max1_match_found %>%
  mutate(accepted_match = TRUE)

finalSynonymy <- add_to_synonymy(finalSynonymy, fuzzy_min1max1_eval)
evaluated_pairs <- distinct(fuzzy_min1max1_eval)

rm(exact, fuzzy_min1max1_eval)
gc()

fuzzy_min2max2 <- harmonize_fuzzy_match(base1_dtf = fuzzy_min1max1$mdd_min1max1_match_failed,
                                        base2_dtf = fuzzy_min1max1$pbdb_min1max1_match_failed,
                                        min_dist = 2, max_dist = 2, delim = "_")

fuzzy_min2max2_eval <- fuzzy_min2max2$min2max2_match_found %>%
  anti_join(evaluated_pairs) %>%
  mutate(accepted_match = case_when(
    MDD_original_combination == "Arctocephalus_tasmanicus" & pbdb_taxon_name == "Arctocephalus_tasmanica" ~ TRUE,
    MDD_original_combination == "Athylax_paludosus" & pbdb_taxon_name == "Athylax_paludinosus" ~ TRUE,
    MDD_original_combination == "Bassaris_astuta" & pbdb_taxon_name == "Bassaris_astutus" ~ TRUE,
    MDD_original_combination == "Mephitis_leuconota" & pbdb_taxon_name == "Mephitis_leuconotus" ~ TRUE,
    MDD_original_combination == "Mephitis_mesoleuca" & pbdb_taxon_name == "Mephitis_mesoleucus" ~ TRUE,
    MDD_original_combination == "Crocotta_kibonotensis" & pbdb_taxon_name == "Crocotta_kibotonensis" ~ TRUE,
    MDD_original_combination == "Hemigalea_zebra" & pbdb_taxon_name == "Hemigale_zebre" ~ TRUE,
    MDD_original_combination == "Hemigalea_derbiana" & pbdb_taxon_name == "Hemigalea_derbianus" ~ TRUE,
    MDD_original_combination == "Mongo_ichneumon" & pbdb_taxon_name == "Mungos_ichneumon" ~ TRUE,
    MDD_original_combination == "Chaeffia_adusta" & pbdb_taxon_name == "Schaeffia_adusta" ~ TRUE,
    MDD_original_combination == "Lutra_congica" & pbdb_taxon_name == "Lutra_congicus" ~ TRUE,
    MDD_original_combination == "Paradoxurus_musschenbroekii" & pbdb_taxon_name == "Paradoxurus_musschenbroecki" ~ TRUE,
    MDD_original_combination == "Poecilictis_libyca" & pbdb_taxon_name == "Poecilictis_lybica" ~ TRUE,
    MDD_original_combination == "Phoca_fetida" & pbdb_taxon_name == "Phoca_foetica" ~ TRUE,
    MDD_original_combination == "Canis_ruppelii" & pbdb_taxon_name == "Canis_rueppellii" ~ TRUE,
   TRUE ~ FALSE
  ))

finalSynonymy <- add_to_synonymy(finalSynonymy, fuzzy_min2max2_eval)
evaluated_pairs <- distinct(bind_rows(evaluated_pairs, fuzzy_min2max2_eval))

rm(fuzzy_min1max1, fuzzy_min2max2_eval)
gc()

fuzzy_min3max3 <- harmonize_fuzzy_match(base1_dtf = fuzzy_min2max2$mdd_min2max2_match_failed,
                                        base2_dtf = fuzzy_min2max2$pbdb_min2max2_match_failed,
                                        min_dist = 3, max_dist = 3, delim = "_")

fuzzy_min3max3_eval <- fuzzy_min3max3$min3max3_match_found %>%
  anti_join(evaluated_pairs) %>%
  mutate(accepted_match = case_when(
    MDD_original_combination == "Profelis_temmincki" & pbdb_taxon_name == "Pardofelis_temminckii" ~ TRUE,
    MDD_original_combination == "Herpestes_aegyptius" & pbdb_taxon_name == "Herpestides_aegypticus" ~ TRUE,
    MDD_original_combination == "Eunothocyon_parvidens" & pbdb_taxon_name == "Nothocyon_parvidens" ~ TRUE,
    MDD_original_combination == "Eunothocyon_urostictus" & pbdb_taxon_name == "Nothocyon_urostictus" ~ TRUE,
    MDD_original_combination == "Calocephalus_vitulinus" & pbdb_taxon_name == "Callocephalus_vitulina" ~ TRUE,
    MDD_original_combination == "Thalassarctos_polaris" & pbdb_taxon_name == "Thalarctos_polaris" ~ TRUE,
    MDD_original_combination == "Zalophus_californicus" & pbdb_taxon_name == "Zalophus_californiana" ~ TRUE,
    TRUE ~ FALSE
  ))

finalSynonymy <- add_to_synonymy(finalSynonymy, fuzzy_min3max3_eval)
evaluated_pairs <- distinct(bind_rows(evaluated_pairs, fuzzy_min3max3_eval))

rm(fuzzy_min2max2, fuzzy_min3max3_eval)
gc()

fuzzy_min4max4 <- harmonize_fuzzy_match(base1_dtf = fuzzy_min3max3$mdd_min3max3_match_failed,
                                        base2_dtf = fuzzy_min3max3$pbdb_min3max3_match_failed,
                                        min_dist = 4, max_dist = 4, delim = "_")

fuzzy_min4max4_eval <- fuzzy_min4max4$min4max4_match_found %>%
  anti_join(evaluated_pairs) %>%
  mutate(accepted_match = case_when(
    MDD_original_combination == "Badiofelis_badia" & pbdb_taxon_name == "Profelis_badia" ~ TRUE,
    MDD_original_combination == "Pelagios_monachus" & pbdb_taxon_name == "Pelagocyon_monachus" ~ TRUE,
    MDD_original_combination == "Mustela_campestris" & pbdb_taxon_name == "Martes_campestris" ~ TRUE,
    MDD_original_combination == "Felis_augusta" & pbdb_taxon_name == "Uncia_augusta" ~ TRUE,
    TRUE ~ FALSE
  )) #Euarctos_floridanus Tremarctos_floridanus

finalSynonymy <- add_to_synonymy(finalSynonymy, fuzzy_min4max4_eval)
evaluated_pairs <- distinct(bind_rows(evaluated_pairs, fuzzy_min4max4_eval))

rm(fuzzy_min3max3, fuzzy_min4max4, fuzzy_min4max4_eval)
gc()



#### Final combination of valid names ####

a <- distinct(select(mddTaxNames_species,
                     MDD_family,
                     MDD_species,
                     MDD_original_combination)) %>%
  inner_join(finalSynonymy) %>%
  inner_join((distinct(select(pbdbTaxNames,
                              pbdb_accepted_name,
                              pbdb_taxon_name)))) %>%
  select(MDD_species, pbdb_accepted_name, everything())

b <- distinct(select(a, MDD_species, pbdb_accepted_name))
b$pbdb_accepted_name <- str_replace_all(b$pbdb_accepted_name, pattern = " ", replacement = "_")
b$equal <- b$MDD_species == b$pbdb_accepted_name


