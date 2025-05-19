# Anotação: A coluna accepted_attr deveria conter as autoridades taxonômicas dos
# nomes aceitos, mas veio vazia por alguma razão. Caso os dados sejam atualizados
# no futuro, esse script poderá ser expandido para formatar a coluna de autores
# e anos.

# Carregando pacotes e dados
library(here)
library(groundhog)
groundhog.day <- "2025-04-01"
groundhog.packages <- c("here",
                        "tidyr",
                        "dplyr",
                        "stringr",
                        "tidylog")

groundhog.library(groundhog.packages,
                  groundhog.day)
rm(groundhog.day,
   groundhog.packages)

pbdb_250506_carn_raw_dtf <- read.csv(here("data",
                                          "raw",
                                          "pbdb_occs_2025-05-06_raw.csv"),
                                     skip = 19)

# Preparação dos dados para a harmonização

pbdb_occs_names <- pbdb_250506_carn_raw_dtf %>%
  select(identified_name, #Selecionando colunas com informação taxonômica
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
  filter(accepted_rank == "species") %>% #Filtrando apenas nomes pertencentes a espécies
  distinct() %>% #Compactand a tabela para conter apenas linhas distintas, não repetidas
  mutate(across(where(is.character), #Substituindo células vazias na tabela toda por NA
                ~na_if(., "")),
         across(.cols = c(identified_name, accepted_name, genus), #Substituindo os espaços por _ nas colunas de nome identificado, aceito e gênero
                ~str_replace_all(., pattern = " ", replacement = "_")),
         family_no = as.integer(na_if(family_no, "NF")), #Substituindo valor "NF" na coluna "family_no" para NA, e transformando a colunas para o tipo "integer", pois ela corresponde ao número de registro das famílias válidas na base
         split_name = accepted_name, #Copiando coluna de nomes aceitos para substituir o gênero do nome aceito pelo gênero da coluna genus
         split_name = str_replace_all(split_name, #Retirando o nome do subgênero para não ficar repetido com a coluna genus
                                      pattern = "\\(.*\\)_",
                                      replacement = "")) %>%
  separate_wider_delim(cols = split_name, #Separando cópia do nome aceito para juntar com a coluna genus
                       names = c("old_genus", "rest"),
                       delim = "_",
                       too_few = "align_start",
                       too_many = "merge") %>%
  mutate(internal_synonym = paste(genus, rest, sep = "_")) %>% #Criando a coluna de sinônimos internos do pbdb
  relocate(internal_synonym, .after = accepted_name) %>% #Reorganizando colunas
  rename_with(~paste0("pbdb_", .x, recycle0 = TRUE)) #Adicionando o sufixo "pbdb_" no nome de todas as colunas

# Salvando o arquivo
write.csv(pbdb_occs_names,
          row.names = FALSE,
          file = here("data",
                      "processed",
                      "pbdb_occs_species-names_2025-05-06.csv"))
