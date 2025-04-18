# Anotação: A coluna accepted_attr deveria conter as autoridades taxonômicas dos
# nomes aceitos, mas veio vazia por alguma razão. Caso os dados sejam atualizados
# no futuro, esse script poderá ser expandido para formatar a coluna de autores
# e anos.

# Carregando pacotes e dados
library(here)
library(groundhog)
groundhog.day <- "2025-04-01"
groundhog.packages <- c("here",
                        "dplyr",
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
         family_no = as.integer(na_if(family_no, "NF")))

# Salvando o arquivo
write.csv(pbdb_occs_names,
          row.names = FALSE,
          file = here("data",
                      "processed",
                      "pbdb_occs_species-names_2025-04-03.csv"))
