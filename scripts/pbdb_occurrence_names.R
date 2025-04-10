# Anotação: A coluna accepted_attr deveria conter as autoridades taxonômicas dos
# nomes aceitos, mas veio vazia por alguma razão. Caso os dados sejam atualizados
# no futuro, esse script poderá ser expandido para formatar a coluna de autores
# e anos.

# Carregando pacotes e dados

library(groundhog)
groundhog.day <- "2025-04-01"
groundhog.packages <- c("here",
                        "tidyr",
                        "dplyr",
                        "stringr",
                        "forcats",
                        "tidylog")

groundhog.library(groundhog.packages,
                  groundhog.day)
rm(groundhog.day,
   groundhog.packages)

pbdb_250308_carn_raw_dtf <- read.csv(here("data",
                                          "raw",
                                          "pbdb_occs_2025-04-03_raw.csv"),
                                     skip = 19)

# Preparação dos dados para a harmonização

pbdb_colnames_vct <- c("pbdb_identified_id" = "identified_no",
                       "pbdb_identified_rank" = "identified_rank",
                       "pbdb_identified_name" = "identified_name",
                       "pbdb_accepted_authority" = "accepted_attr",
                       "pbdb_status_difference" = "difference",
                       "pbdb_accepted_id" = "accepted_no",
                       "pbdb_accepted_rank" = "accepted_rank",
                       "pbdb_accepted_name" = "accepted_name",
                       "pbdb_accepted_sist_phylum_name" = "phylum",
                       "pbdb_accepted_sist_phylum_id" = "phylum_no",
                       "pbdb_accepted_sist_class_name" = "class",
                       "pbdb_accepted_sist_class_id" = "class_no",
                       "pbdb_accepted_sist_order_name" = "order",
                       "pbdb_accepted_sist_order_id" = "order_no",
                       "pbdb_accepted_sist_family_name" = "family",
                       "pbdb_accepted_sist_family_id" = "family_no",
                       "pbdb_accepted_sist_genus_name" = "genus",
                       "pbdb_accepted_sist_genus_id" = "genus_no",
                       "pbdb_accepted_sist_subgenus_id" = "subgenus_no",
                       "pbdb_accepted_sist_subgenus_name" = "subgenus_name",
                       "pbdb_info_subgenus_resolution" = "subgenus_reso",
                       "pbdb_info_primary_resolution" = "primary_reso",
                       "pbdb_info_species_resolution" = "species_reso")

pbdb_occs_names <- pbdb_250308_carn_raw_dtf %>%
  select(identified_name,
         difference,
         accepted_name,
         identified_rank,
         accepted_rank,
         accepted_attr,
         phylum,
         class,
         order,
         family,
         genus,
         subgenus_name,
         subgenus_reso,
         primary_reso,
         species_reso,
         phylum_no,
         class_no,
         order_no,
         family_no,
         genus_no,
         subgenus_no,
         identified_no,
         accepted_no) %>%
  mutate(across(where(is.character),
                ~na_if(., "")),
         family_no = as.integer(na_if(family_no, "NF"))) %>%
  rename(all_of(pbdb_colnames_vct))

# Salvando o arquivo
write.csv(pbdb_occs_names,
          row.names = FALSE,
          file = here("data",
                      "processed",
                      "pbdb_backbone_carn_250403.csv"))
