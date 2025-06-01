#### Loading packages ####

library(here) #Operating system agnostic file paths with .Rproj
library(stringr) #String manipulation
library(dplyr) #Data wrangling with tidyverse syntax
library(tidytable) #Brings the data.table package's fast operations to dplyr syntax
library(tidylog) #Optional, prints detailed information about changes in the dataframe during data wrangling with tidyverse functions
library(countrycode) #Needed to summarise country names into continents

#### Setting file paths ####

now.occurrences.path <- here("data",
                             "raw",
                             "now_occs_2025-05-06_raw.csv")

now.functions.path <- here("scripts",
                           "functions",
                           "cleaning_now.R")

#### Loading data and functions####

now.occurrences.data <- read.delim(now.occurrences.path,
                                   na.strings = "\\N")

source(now.functions.path)

rm(now.occurrences.path,
   now.functions.path)
gc()

#### Preparing function arguments ####

select.cols_to_keep <- c("SIDNUM",
                         "SUBCLASSORSUPERORDER",
                         "ORDER",
                         "SUBORDERORSUPERFAMILY",
                         "FAMILY",
                         "SUBFAMILY",
                         "GENUS",
                         "SPECIES") #Uncomment bellow to check available columns
#colnames(now.occurrences.data)

select.taxon_id_col <- "SIDNUM"
select.species_col <- "SPECIES"
select.genus_col <- "GENUS"
select.subfamily_col <- "SUBFAMILY"
select.family_col <- "FAMILY"
select.suborder_superfamily_col <- "SUBORDERORSUPERFAMILY"
select.order_col <- "ORDER"
select.subclass_superorder_col <- "SUBCLASSORSUPERORDER"
select.min_age_col <- "MIN_AGE"
select.max_age_col <- "MAX_AGE"
select.country_col <- "COUNTRY"

select.new_name_separator <- "_"

#### Filtering and fixing indet. values ####

now.taxa <- now.occurrences.data %>%
  select(all_of(select.cols_to_keep)) %>%
  distinct() %>%
  mutate(across(is.character,
                ~now.fix_indet(.)))

#### Adding summarised data ####

now.taxa <- now.summarise_ages_by_id(data = now.occurrences.data,
                                     taxon_id = select.taxon_id_col,
                                     min_age_col = select.min_age_col,
                                     max_age_col = select.max_age_col) %>%
  inner_join(now.taxa, .)

#### Creating new columns ####

now.taxa <- now.create_species_name_col(data = now.taxa,
                                        new_col_name = "now_accepted_name",
                                        species_col = select.species_col,
                                        genus_col = select.genus_col,
                                        sep = select.new_name_separator)

now.taxa <- now.create_accepted_rank_col(data = now.taxa,
                                         new_col_name = "now_accepted_rank",
                                         species_col = select.species_col,
                                         genus_col = select.genus_col,
                                         subfamily_col = select.subfamily_col,
                                         family_col = select.family_col,
                                         suborder_superfamily_col = select.suborder_superfamily_col,
                                         order_col = select.order_col,
                                         subclass_superorder_col = select.subclass_superorder_col)








#### Old comments to translate ####
#### Limpeza inicial

# Nesse pedaço, eu estou lendo o arquivo baixado direto do NOW com separação por tabs, por isso
# ele é tsv e não csv, mas tenho quase certeza que um arquivo baixado pelo script do Carlos
# "mod1_read_wrangle_now_raw.Rmd" funcionaria com esse script normalmente, só precisaria mudar
# a função de read.delim pra read.csv.

# A função abaixo vai servir para substituir todos os valores especificados
#(NA, "", "\\N", "sp.", "Gen.", "gen.", "incertae sedis") por "indet.", sendo
#apenas necessário especificar o nome da coluna no parâmetro col_name. A função
#`str_replace_all()` seria mais simples de usar, mas ela não consegue substituir as
#células vazias da tabela por razões que desconheço, enquanto a `case_when()` encaixa
#bem na pipeline do tidyverse, consegue substituir as células vazias e funciona bem
#dentro da função `across()` (usada mais abaixo) para checar todas as colunas da tabela
#de uma vez. Certamente alguma função do R base deve fazer algo equivalente com
#maior eficiência, então se alguém quiser escrever de outra forma, é bem-vind@.

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
         SPECIES,
         SIDNUM) %>%
  distinct() %>% #Filtrando linhas com valores distintos (considerando todas as colunas)
  mutate(across(is.character, #Aplicando a função indet_fix em todas as colunas da tabela
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

# Preparando alguns resumos de dados sobre as ocorrências, já que o NOW não tem a autoria dos nomes
age_by_sidnum <- now_data %>%
  select(MIN_AGE, MAX_AGE, SIDNUM) %>%
  distinct() %>%
  group_by(SIDNUM) %>%
  summarise(MAX_AGE = max(MAX_AGE),
            MIN_AGE = min(MIN_AGE))

continent_by_sidnum <- now_data %>%
  select(COUNTRY,
         SIDNUM) %>%
  distinct() %>%
  mutate(now_continent = countrycode::countrycode(COUNTRY,
                                                  origin = "country.name.en",
                                                  destination = "continent")) %>%
  select(-COUNTRY) %>%
  distinct() %>%
  group_by(SIDNUM) %>%
  mutate(now_continent = paste0(now_continent,
                                collapse = " | ")) %>%
  distinct()


# Filtrando só as espécies para o workshop de Carnivora

now_species <- now_taxa %>%
  filter(now_accepted_rank == "species")

now_species <- now_species %>%
  inner_join(age_by_sidnum) %>%
  inner_join(continent_by_sidnum)

# Salvando o arquivo
write.csv(now_species,
          file = here("data",
                      "processed",
                      "now_occs_species-names_2025-05-06.csv"),
          row.names = FALSE)


