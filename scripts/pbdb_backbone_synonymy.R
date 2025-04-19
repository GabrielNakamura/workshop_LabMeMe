#### Carregando pacotes com controle de versão ####

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

#### Limpeza inicial ####
# Carregar arquivo #

pbdb_taxa <- read.csv(here("data",
                             "raw",
                             "pbdb_taxa_2025-04-03_raw.csv"), skip = 17)

# Vetores de correção #

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
                       "pbdb_reference_author" = "ref_author",
                       "pbdb_reference_year" = "ref_pubyr",
                       "pbdb_reference_id" = "reference_no",
                       "pbdb_accepted_status_extinct" = "is_extant",
                       "pbdb_accepted_info_count_occurences" = "n_occs",
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
                       "pbdb_accepted_info_count_orders" = "n_orders",
                       "pbdb_accepted_info_count_families" = "n_families",
                       "pbdb_accepted_info_count_genera" = "n_genera",
                       "pbdb_accepted_info_count_species" = "n_species",
                       "pbdb_accepted_info_environment" = "taxon_environment", #apesar do nome "taxon", uma análise dos dados indica que os valores são consistentes com accepted_name, não com os taxon_name
                       "pbdb_accepted_info_environment_basis" = "environment_basis",
                       "pbdb_accepted_info_motility" = "motility",
                       "pbdb_accepted_info_life_habit" = "life_habit",
                       "pbdb_accepted_info_vision" = "vision",
                       "pbdb_accepted_info_diet" = "diet",
                       "pbdb_accepted_info_reproduction" = "reproduction",
                       "pbdb_accepted_info_ontogeny" = "ontogeny",
                       "pbdb_accepted_info_ecospace_comments" = "ecospace_comments",
                       "pbdb_accepted_info_composition" = "composition",
                       "pbdb_accepted_info_architecture" = "architecture",
                       "pbdb_accepted_info_thickness" = "thickness",
                       "pbdb_accepted_info_reinforcement" = "reinforcement")

info_typo_fix_vct <- c("Cusafont" = "Crusafont",
                       "Geoffroy-Saint-Hillaire" = "Geoffroy-Saint-Hilaire",
                       "Pockock" = "Pocock",
                       "Rossenbüller" = "Rosenmüller",
                       "Solunias" = "Solounias", #https://sjpp.springeropen.com/articles/10.1007/s13358-012-0042-y
                       "Zimmerman " = "Zimmermann ",
                       "Kretsoi" = "Kretzoi")


# Renomeando colunas e corrigindo erros de digitação #
pbdb_taxa <- pbdb_taxa %>%
  rename(all_of(pbdb_colnames_vct)) %>%
  mutate(across(where(is.character),
                ~str_replace_all(., info_typo_fix_vct)))

#### Criação e formatação de colunas de autores aceitos ####

# O PBDB não tem uma coluna de autoria dos accepted_name, só dos taxon_name, porque todo accepted_name
# também é um taxon_name com autoria. A ideia aqui é filtrar os taxon_name que também são accepted_name e
# usar a autoria do taxon_name pra criar a coluna de accepted_authority. O primeiro passo é marcar que
# nomes de autor têm parênteses em uma nova coluna para começar a padronização do nome dos autores.

accepted_taxa_attr <- pbdb_taxa %>%
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
      ) ~ "TRUE",
      !str_detect(pbdb_accepted_authority,
                  pattern = "\\("
      ) ~ "FALSE"
    ),
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

# Filtrando autoridades que tem ano e criando coluna de ano
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

# Juntando com autoridades sem ano e formatando nomes dos autores aceitos para retirar pontos,
# vírgulas, espaços, separar autores por "_" e nomes compostos por "-", tirar acentos e deixar
# tudo em maiúscula. Tudo isso é feito para padronizar ao máximo as strings dos nomes de autores
# e facilitar os processos de join para a harmonização taxonômica com outras bases de dados.

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
pbdb_taxa <- pbdb_taxa %>%
  full_join(., accepted_taxa_attr_format) %>%
  mutate(pbdb_accepted_name = str_replace_all(pbdb_accepted_name,
                                              pattern = " ",
                                              replacement = "_"),
         pbdb_taxon_name = str_replace_all(pbdb_taxon_name,
                                           pattern = " ",
                                           replacement = "_"),
         pbdb_accepted_sist_family_name = toupper(pbdb_accepted_sist_family_name)) %>%
  distinct()

#### Criação e formatação de colunas de autores aceitos ####

original_taxa <- pbdb_taxa %>%
  select(pbdb_taxon_id,
        pbdb_taxon_authority) %>%
  mutate(pbdb_taxon_authority_copy = pbdb_taxon_authority,
    pbdb_taxon_authority_parentheses = case_when(
      str_detect(pbdb_taxon_authority,
                 pattern = "\\("
      ) ~ "TRUE",
      !str_detect(pbdb_taxon_authority,
                  pattern = "\\("
      ) ~ "FALSE"
    ),
    pbdb_taxon_authority_copy = str_replace_all(pbdb_taxon_authority_copy,
                                    pattern = "\\(|\\)|,|\\.",
                                    replacement = ""
    ),
    pbdb_taxon_authority_copy = str_replace_all(pbdb_taxon_authority_copy,
                                    pattern = "  ",
                                    replacement = " "
    )
  ) %>%
  distinct()




pbdb_backbone_carn_250308 <- pbdb_taxa %>%
  select(pbdb_taxon_rank,
         pbdb_accepted_rank,
         pbdb_taxon_name,
         pbdb_accepted_name,
         pbdb_taxon_authority,
         pbdb_accepted_authority,
         pbdb_accepted_author_name_format,
         pbdb_accepted_author_year,
         pbdb_accepted_author_year,
         pbdb_accepted_authority_parentheses)

write.csv(pbdb_backbone_carn_250308,
          file = here("01_data",
                      "02_clean-data",
                      "pbdb_backbone_carn_250308.csv"))

