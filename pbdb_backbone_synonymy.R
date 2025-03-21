# Carregando pacotes

library(groundhog)
groundhog.day <- "2025-03-01"
groundhog.packages <- c("here",
                        "readr",
                        "tidyr",
                        "dplyr",
                        "stringr",
                        "forcats",
                        "tidylog")

groundhog.library(groundhog.packages, 
                  groundhog.day)
rm(groundhog.day,
   groundhog.packages)

# Limpeza inicial

pbdb_250308_carn_raw_dtf <- read_csv(here("01_data", 
                                          "01_raw-data",
                                          "pbdb_data_250318_carnivora_namao.csv"), 
                                     skip = 15)

pbdb_colnames_vct <- c("pbdb_original_id" = "orig_no",
                       "pbdb_taxon_id" = "taxon_no",
                       "pbdb_taxon_status_record_type" = "record_type",
                       "pbdb_taxon_status_flagged" = "flags",
                       "pbdb_taxon_rank" = "taxon_rank",
                       "pbdb_taxon_name" = "taxon_name",
                       "pbdb_taxon_authority" = "taxon_attr",
                       "pbdb_common_name" = "common_name",
                       "pbdb_status_difference" = "difference",
                       "pbdb_accepted_id" = "accepted_no",
                       "pbdb_accepted_rank" = "accepted_rank",
                       "pbdb_accepted_name" = "accepted_name",
                       "pbdb_accepted_parent_id" = "parent_no",
                       "pbdb_accepted_parent_name" = "parent_name",
                       #"pbdb_accepted_immediate_parent_id" = "immpar_no",
                       #"pbdb_accepted_immediate_parent_name" = "immpar_name",
                       "pbdb_reference_id" = "reference_no",
                       "pbdb_accepted_status_extinct" = "is_extant",
                       "pbdb_accepted_count_occurences" = "n_occs",
                       "pbdb_accepted_info_firstapp_max_ma" = "firstapp_max_ma",
                       "pbdb_accepted_info_firstapp_min_ma" = "firstapp_min_ma",
                       "pbdb_accepted_info_lastapp_max_ma" = "lastapp_max_ma",
                       "pbdb_accepted_info_lastapp_min_ma" = "lastapp_min_ma",
                       "pbdb_accepted_info_early_interval" = "early_interval",
                       "pbdb_accepted_info_late_interval" = "late_interval",
                       "pbdb_accepted_info_taxonsize" = "taxon_size",
                       "pbdb_accepted_info_extansize" = "extant_size",
                       "pbdb_accepted_sist_phylum_name" = "phylum",
                       #"pbdb_accepted_sist_phylum_id" = "phylum_no",
                       "pbdb_accepted_sist_class_name" = "class",
                       #"pbdb_accepted_sist_class_id" = "class_no",
                       "pbdb_accepted_sist_order_name" = "order",
                       #"pbdb_accepted_sist_order_id" = "order_no",
                       "pbdb_accepted_sist_family_name" = "family",
                       #"pbdb_accepted_sist_family_id" = "family_no",
                       "pbdb_accepted_sist_genus_name" = "genus",
                       #"pbdb_accepted_sist_genus_id" = "genus_no",
                       #"pbdb_accepted_sist_subgenus_id" = "subgenus_no",
                       "pbdb_accepted_sist_type_name" = "type_taxon",
                       #"pbdb_accepted_sist_type_id" = "type_taxon_no",
                       "pbdb_accepted_count_orders" = "n_orders",
                       "pbdb_accepted_count_families" = "n_families",
                       "pbdb_accepted_count_genera" = "n_genera",
                       "pbdb_accepted_count_species" = "n_species")

info_typo_fix_vct <- c("Cusafont" = "Crusafont",
                       "Geoffroy-Saint-Hillaire" = "Geoffroy-Saint-Hilaire",
                       "Pockock" = "Pocock",
                       "Rossenbüller" = "Rosenmüller",
                       "Solunias" = "Solounias", #https://sjpp.springeropen.com/articles/10.1007/s13358-012-0042-y
                       "Zimmerman " = "Zimmermann ",
                       "Kretsoi" = "Kretzoi") 

pbdb_250308_carn_raw_dtf <- pbdb_250308_carn_raw_dtf %>%
  rename(all_of(pbdb_colnames_vct)) %>%
  mutate(across(where(is.character), 
                ~str_replace_all(., info_typo_fix_vct)))

accepted_taxa_attr <- pbdb_250308_carn_raw_dtf %>%
  select(pbdb_taxon_name,
         pbdb_accepted_name,
         pbdb_taxon_authority,
         pbdb_taxon_id,
         pbdb_accepted_id
  ) %>%
  mutate(pbdb_taxon_authority_copy = pbdb_taxon_authority) %>%
  filter(pbdb_accepted_name == pbdb_taxon_name) %>%
  rename(pbdb_accepted_authority = pbdb_taxon_authority_copy) %>%
  mutate(
    pbdb_accepted_authority_parentheses = case_when(
      str_detect(pbdb_accepted_authority,
                 pattern = "\\("
      ) ~ "yes",
      !str_detect(pbdb_accepted_authority,
                  pattern = "\\("
      ) ~ "no"
    ),
    pbdb_accepted_authority_parentheses = as.factor(pbdb_accepted_authority_parentheses),
    accepted_copy = pbdb_accepted_authority,
    accepted_copy = str_replace_all(accepted_copy,
                                    pattern = "\\(|\\)|,|\\.",
                                    replacement = ""
    ),
    accepted_copy = str_replace_all(accepted_copy,
                                    pattern = "  ",
                                    replacement = " "
    )
  ) %>%
  distinct()

#Filtrando autoridades que tem ano e criando coluna de ano
accepted_taxa_attr_format <- accepted_taxa_attr %>% 
  filter(str_detect(accepted_copy, 
                    pattern = "[:digit:]")) %>% 
  mutate(accepted_copy = stringi::stri_replace_last_fixed(accepted_copy, 
                                                          pattern = " ",
                                                          replacement = "~")) %>%
  separate(accepted_copy, 
           into = c("pbdb_accepted_author_name_format", 
                    "pbdb_accepted_author_year"), 
           sep = "~") %>%
  
  distinct()

#Juntando com autoridades sem ano e formatando nomes dos autores
accepted_taxa_attr_format <- full_join(accepted_taxa_attr, accepted_taxa_attr_format) %>%
  mutate(pbdb_accepted_author_name_format = case_when(
    is.na(pbdb_accepted_author_name_format) ~ pbdb_accepted_authority,
    TRUE ~ pbdb_accepted_author_name_format)) %>%
  mutate(pbdb_accepted_author_name_format = stringi:: stri_trans_general(pbdb_accepted_author_name_format, "Latin-ASCII"),
         pbdb_accepted_author_name_format = toupper(pbdb_accepted_author_name_format),
         pbdb_accepted_author_name_format = str_replace_all(pbdb_accepted_author_name_format,
                                                            pattern = " ET AL",
                                                            replacement = "_ET-AL"),
         pbdb_accepted_author_name_format = str_replace_all(pbdb_accepted_author_name_format,
                                                            pattern = " AND ",
                                                            replacement = "_"),
         pbdb_accepted_author_name_format = str_replace_all(pbdb_accepted_author_name_format,
                                                            pattern = " AND_",
                                                            replacement = "_"),
         pbdb_accepted_author_name_format = str_replace_all(pbdb_accepted_author_name_format,
                                                            pattern = " ",
                                                            replacement = "-")) %>%
  select(-accepted_copy,
         -pbdb_taxon_name,
         -pbdb_taxon_id,
         -pbdb_taxon_authority) %>%
  distinct()

#Juntando com resto da tabela, formantando nomes de taxa e família
pbdb_250308_carn_raw_dtf <- pbdb_250308_carn_raw_dtf %>%
  full_join(., accepted_taxa_attr_format) %>%
  mutate(pbdb_accepted_name = str_replace_all(pbdb_accepted_name, 
                                              pattern = " ", 
                                              replacement = "_"),
         pbdb_taxon_name = str_replace_all(pbdb_taxon_name, 
                                           pattern = " ", 
                                           replacement = "_"),
         pbdb_accepted_sist_family_name = toupper(pbdb_accepted_sist_family_name)) %>%
  distinct()

pbdb_backbone_carn_250308 <- pbdb_250308_carn_raw_dtf %>%
  select(pbdb_taxon_name,
         pbdb_taxon_authority,
         pbdb_accepted_name,
         pbdb_accepted_rank,
         pbdb_accepted_authority,
         pbdb_accepted_author_name_format,
         pbdb_accepted_author_year,
         pbdb_accepted_author_year,
         pbdb_accepted_authority_parentheses)

write.csv(pbdb_backbone_carn_250308,
          file = here("01_data",
                      "02_clean-data",
                      "pbdb_backbone_carn_250308.csv"))
