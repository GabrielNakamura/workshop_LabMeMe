#### Carregando pacotes com controle de versão ####

library(groundhog)
groundhog.day <- "2025-04-01"
groundhog.packages <- c("here",
                        "dplyr",
                        "stringr",
                        "tidylog")

groundhog.library(groundhog.packages,
                  groundhog.day)
rm(groundhog.day,
   groundhog.packages)

#### Limpeza inicial ####

# Nesse pedaço, eu estou lendo o arquivo baixado direto do NOW com separação por tabs, por isso
# ele é tsv e não csv, mas tenho quase certeza que um arquivo baixado pelo script do Carlos
# "mod1_read_wrangle_now_raw.Rmd" funcionaria com esse script normalmente, só precisaria mudar
# a função de read.delim pra read.csv.

now_data <- read.delim(here("data",
                            "raw",
                            "now_occs_2025-05-06_raw.csv"),
                       na.strings = "\\N")


# A função abaixo vai servir para substituir todos os valores especificados
#(NA, "", "\\N", "sp.", "Gen.", "gen.", "incertae sedis") por "indet.", sendo
#apenas necessário especificar o nome da coluna no parâmetro col_name. A função
#`str_replace_all()` seria mais simples de usar, mas ela não consegue substituir as
#células vazias da tabela por razões que desconheço, enquanto a `case_when()` encaixa
#bem na pipeline do tidyverse, consegue substituir as células vazias e funciona bem
#dentro da função `across()` (usada mais abaixo) para checar todas as colunas da tabela
#de uma vez. Certamente alguma função do R base deve fazer algo equivalente com
#maior eficiência, então se alguém quiser escrever de outra forma, é bem-vind@.

indet_fix <- function(col_name) {
  case_when(
    is.na(col_name) | col_name == "" | col_name == "\\N" | col_name == "sp." | col_name == "Gen." | col_name == "gen." | col_name == "incertae sedis" ~ "indet.",
    TRUE ~ col_name
  )
}

# O pipeline a seguir seleciona as colunas de interesse taxonômico e realiza as operações
#de formatação e padronização dos dados. A parte mais complexa dela é o uso sequencial da
#função `case_when()` (dentro da função `mutate()`) para criar a coluna de ranking taxonômico
#do nome aceito, a fim de facilitar a filtragem dos nomes do NOW como é feito com o PBDB.
#Especificamente, a função `case_when()` preenche os valores da nova coluna `now_accepted_rank`
#de forma condicional aos valores das outras colunas com informação taxonômica, como uma
#operação if-else, seguindo a sintaxe `condição verdadeira ~ "string a ser adicionada na coluna"`.
#Como as linhas da coluna que não correspondem à nenhuma condição verdadeira são transformadas
#em NA, a sintaxe `TRUE ~ nome da coluna` deve ser usada para preservar eventuais valores
#anteriores.

now_taxa <- now_data %>%
  select(SUBCLASSORSUPERORDER, #Selecionando colunas de interesse taxonômico
         ORDER,
         SUBORDERORSUPERFAMILY,
         FAMILY,
         SUBFAMILY,
         GENUS,
         SPECIES) %>%
  distinct() %>% #Filtrando linhas com valores distintos (considerando todas as colunas)
  mutate(across(everything(), #Aplicando a função indet_fix em todas as colunas da tabela
                ~indet_fix(.)),
         now_accepted_name = paste(GENUS, SPECIES, sep = "_")) %>% #Criando uma coluna com o nome completo das espécies separado por "_"
  mutate(now_accepted_rank = case_when( #Usando várias vezes a função case_when(), ver detalhes acima
    SPECIES != "indet." ~ "species",
    SPECIES == "indet." & GENUS != "indet." ~ "genus",
    SPECIES == "indet." & GENUS == "indet." ~ "Other"
  ),
  now_accepted_rank = case_when(
    now_accepted_rank == "Other" & SUBFAMILY != "indet." ~ "subfamily_tribe",
    TRUE ~ now_accepted_rank
  ),
  now_accepted_rank = case_when(
    now_accepted_rank == "Other" & FAMILY != "indet." ~ "family",
    TRUE ~ now_accepted_rank
  ),
  now_accepted_rank = case_when(
    now_accepted_rank == "Other" & SUBORDERORSUPERFAMILY != "indet." ~ "suborder_superfamily",
    TRUE ~ now_accepted_rank
  ),
  now_accepted_rank = case_when(
    now_accepted_rank == "Other" & ORDER != "indet." ~ "order",
    TRUE ~ now_accepted_rank
  ),
  now_accepted_rank = case_when(
    now_accepted_rank == "Other" & SUBCLASSORSUPERORDER != "indet." ~ "subclass_superorder",
    TRUE ~ now_accepted_rank
  ),
  now_accepted_name = str_replace(now_accepted_name, #Tirando "_indet."s desnecessário na coluna now_accepted_name
                                  pattern = "_indet\\.",
                                  replacement = ""))

# Filtrando só as espécies para o workshop de Carnivora

now_species <- now_taxa %>%
  filter(now_accepted_rank == "species")

# Salvando o arquivo
write.csv(now_species,
          file = here("data",
                      "processed",
                      "now_occs_species-names_2025-05-06.csv"),
          row.names = FALSE)



