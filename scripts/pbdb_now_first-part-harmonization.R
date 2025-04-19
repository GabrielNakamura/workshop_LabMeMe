#Informações:
# Para este script, você precisa de tabelas de 2 bases de dados com nomes aceitos
# de ocorrências e o rank taxonômico de cada nome. No momento, ele está escrito
# para a harmonizar as bases PBDB e NOW, mas ele pode ser ajustado no futuro para
# atender a quaisquer 2 bases.
#
# A harmonização é focada em gênero, espécie e subspécie, que são os ranks com
# nomes mais mutáveis e mais úteis em estudos de macroecologia e macroevolução.
#
# Os nomes de gênero, epíteto de espécie e epíteto de subspécie são retirados
# diretamente do nome aceito, e não de colunas específicas das bases. Por exemplo,
# o PBDB tem uma coluna específica para o gênero que frequentemente não condiz
# com o gênero do nome aceito, então aqui é usado o gênero do nome aceito.
#
# Como o NOW não tem nome de autoria, a autoria não foi checada aqui, mas isso
# pode ser expandido no futuro

#### Carregando pacotes com controle de versão ####

# O pacote groundhog instala as versões dos outros pacotes disponíveis no CRAN
# no dia especificado.

library(groundhog)
groundhog.day <- "2025-04-01"
groundhog.packages <- c("here",
                        "dplyr",
                        "stringr",
                        "fuzzyjoin",
                        "tidylog")

groundhog.library(groundhog.packages,
                  groundhog.day)
rm(groundhog.day,
   groundhog.packages)

#### Carregar arquivos e selecionar colunas ####

pbdb_taxa <- read.csv(here("data",
                           "processed",
                           "pbdb_occs_species-names_2025-04-03.csv")) %>%
  select(pbdb_accepted_name,
         pbdb_internal_synonym,
         pbdb_accepted_rank) %>%
  distinct()

now_taxa <- read.csv(here("data",
                          "processed",
                          "now_occs_species-names_2025-04-03.csv")) %>%
  select(now_accepted_name,
         now_accepted_rank) %>%
  distinct()

#### Preparar colunas com cada parte do nome aceito

separate_name_parts <- function(x, #your dataframe
                                col_accepted_name, #occurrence accepted name column
                                col_accepted_rank, #accepted name taxonomic rank column
                                ranks_to_keep = c("genus", "species", "subgenus", "subspecies"), #up to genus, as they are written in your dataframe
                                source_suffix #suffix to add to column names so you know from which base the taxonomic names come from
                                ){

    dtf <- x %>%
           filter(eval(parse(text = col_accepted_rank)) %in% ranks_to_keep) %>%
           separate_wider_delim(cols = col_accepted_name,
                                delim = "_",
                                names = c("suffix_genus", "suffix_other"),
                                too_many = "merge",
                                too_few = "align_start",
                                cols_remove = FALSE) %>%
           separate_wider_delim(cols = "suffix_other",
                                delim = ")",
                                names = c("suffix_subgenus", "suffix_other"),
                                too_few = "align_end") %>%
           mutate(suffix_subgenus = str_replace_all(suffix_subgenus,
                                                    pattern = "\\(",
                                                    replacement = ""),
                  suffix_other = str_replace(suffix_other,
                                             pattern = "^_",
                                             replacement = "")) %>%
           separate_wider_delim(cols = "suffix_other",
                                delim = "_",
                                names = c("suffix_speciesEpiteth", "suffix_subspeciesEpiteth"),
                                too_few = "align_start")

  colnames(dtf) <- str_replace_all(colnames(dtf),
                                   pattern = "suffix",
                                   replacement = source_suffix)
  return(dtf)

  }

pbdb_taxa <- separate_name_parts(pbdb_taxa,
                                 col_accepted_name = "pbdb_accepted_name",
                                 col_accepted_rank = "pbdb_accepted_rank",
                                 source_suffix = "pbdb")

now_taxa <- separate_name_parts(now_taxa,
                                col_accepted_name = "now_accepted_name",
                                col_accepted_rank = "now_accepted_rank",
                                ranks_to_keep = c("genus", "subgenus", "species", "subspecies"),
                                source_suffix = "now")

#### União das bases ####

##### União exata #####

harmonize_exact_match <- function(base1_dtf,
                                  base2_dtf,
                                  base1_suffix,
                                  base2_suffix,
                                  base1_col,
                                  base2_col) {

  by_cols <- base2_col
  names(by_cols) <- base1_col

  select_cols <- c(base1_col, base2_col)

  match_exact <- full_join(base1_dtf, base2_dtf,
                           by = by_cols,
                           keep = TRUE) %>%
    mutate(base1_present = !is.na(eval(parse(text = base1_col))),
           base2_present = !is.na(eval(parse(text = base2_col)))) %>%
    distinct()

  exact_match_found <- match_exact %>%
    filter(base1_present == TRUE & base2_present == TRUE) %>%
    mutate(string_distance = 0,
           match_notes = "exact match") %>%
    select(all_of(select_cols),
           string_distance,
           match_notes)

  base1_failed_match <- match_exact %>%
    filter(base1_present == FALSE) %>%
    select(any_of(colnames(base2_dtf)))

  base2_failed_match <- match_exact %>%
    filter(base2_present == FALSE) %>%
    select(any_of(colnames(base1_dtf)))

  result <- list(exact_match_found = exact_match_found,
                 base1_failed_match = base1_failed_match,
                 base2_failed_match = base2_failed_match)

  names(result) <- str_replace_all(names(result),
                                   c("base1" = base1_suffix,
                                     "base2" = base2_suffix))

  return(result)

  #adicionar parte de checar se o sufixo já vem no nome das colunas de nomes ou se precisa colocar
  #adicionar o sumariozinho aqui dentro também
}

exact_match <- harmonize_exact_match(base1_dtf = pbdb_taxa,
                                     base2_dtf = now_taxa,
                                     base1_suffix = "pbdb",
                                     base2_suffix = "now",
                                     base1_col = "pbdb_accepted_name",
                                     base2_col = "now_accepted_name")


##### União com 1 ou 2 letras diferentes #####

match_likely_typos <- stringdist_full_join(select(filter(match_exact_failed,
                                    now_present == FALSE),
                                    pbdb_accepted_name,
                                    pbdb_accepted_rank),
                      select(filter(match_exact_failed,
                                    pbdb_present == FALSE),
                             now_accepted_name,
                             now_accepted_rank),
                      by = accepted_name_cols_vct,
                      distance_col = "string_distance",
                      max_dist = 2) %>%
  mutate(pbdb_present = !is.na(pbdb_accepted_name),
         now_present = !is.na(now_accepted_name)) %>%
  distinct()

match_likely_typos_found <- match_likely_typos %>%
  filter(pbdb_present == TRUE & now_present == TRUE) %>%
  select(-pbdb_present,
         -now_present) %>%
  left_join(select(filter(match_exact_failed,
                          pbdb_present == TRUE),
                   starts_with("pbdb_"))) %>%
  left_join(select(filter(match_exact_failed,
                          now_present == TRUE),
                   starts_with("now_"))) %>%
  mutate(genus_dif = case_when(stringdist::stringdist(pbdb_genus, now_genus) > 0 ~ "different spelling on genus"),
         subgenus_dif = case_when(stringdist::stringdist(pbdb_subgenus, now_subgenus) > 0 ~ "different spelling on subgenus"),
         speciesEpi_dif = case_when(stringdist::stringdist(pbdb_speciesEpiteth, now_speciesEpiteth) > 0 ~ "different spelling on species epiteth"),
         subspeciesEpi_dif = case_when(stringdist::stringdist(pbdb_subspeciesEpiteth, now_subspeciesEpiteth) > 0 ~ "different spelling on subspecies epiteth"),
         match_notes = paste(genus_dif,
                             subgenus_dif,
                             speciesEpi_dif,
                             subspeciesEpi_dif,
                             sep = ","),
         match_notes = str_replace_all(match_notes,
                                       pattern = "NA,|,NA",
                                       replacement = "")) %>%
  select(pbdb_accepted_name,
         pbdb_accepted_rank,
         now_accepted_name,
         now_accepted_rank,
         string_distance,
         match_notes)

match_likely_typos_failed <- match_likely_typos %>%
  filter(pbdb_present == FALSE | now_present == FALSE)

##### União com 3 ou 4 letras diferentes #####

match_unlikely_typos <- stringdist_full_join(select(filter(match_likely_typos_failed,
                                                           now_present == FALSE),
                                                    pbdb_accepted_name,
                                                    pbdb_accepted_rank),
                                             select(filter(match_exact_failed,
                                                           pbdb_present == FALSE),
                                                    now_accepted_name,
                                                    now_accepted_rank),
                                             by = accepted_name_cols_vct,
                                             distance_col = "string_distance",
                                             max_dist = 4) %>%
  mutate(pbdb_present = !is.na(pbdb_accepted_name),
         now_present = !is.na(now_accepted_name)) %>%
  distinct()

match_unlikely_typos_found <- match_unlikely_typos %>%
  filter(pbdb_present == TRUE & now_present == TRUE) %>%
  select(-pbdb_present,
         -now_present) %>%
  left_join(select(filter(match_exact_failed,
                          pbdb_present == TRUE),
                   starts_with("pbdb_"))) %>%
  left_join(select(filter(match_exact_failed,
                          now_present == TRUE),
                   starts_with("now_"))) %>%
  mutate(genus_dif = case_when(stringdist::stringdist(pbdb_genus, now_genus) > 0 ~ "different spelling on genus"),
         subgenus_dif = case_when(stringdist::stringdist(pbdb_subgenus, now_subgenus) > 0 ~ "different spelling on subgenus"),
         speciesEpi_dif = case_when(stringdist::stringdist(pbdb_speciesEpiteth, now_speciesEpiteth) > 0 ~ "different spelling on species epiteth"),
         subspeciesEpi_dif = case_when(stringdist::stringdist(pbdb_subspeciesEpiteth, now_subspeciesEpiteth) > 0 ~ "different spelling on subspecies epiteth"),
         match_notes = paste(genus_dif,
                             subgenus_dif,
                             speciesEpi_dif,
                             subspeciesEpi_dif,
                             sep = ","),
         match_notes = str_replace_all(match_notes,
                                       pattern = "NA,|,NA",
                                       replacement = "")) %>%
  select(pbdb_accepted_name,
         pbdb_accepted_rank,
         now_accepted_name,
         now_accepted_rank,
         string_distance,
         match_notes)

#### Resultado geral ####

match_exact_summary_stats <- function(){

  summ1 <- match_exact$pbdb_accepted_name %>%
    na.omit() %>%
    unique() %>%
    length()

  summ2 <- match_exact$now_accepted_name %>%
    na.omit() %>%
    unique() %>%
    length()

  summ3 <- match_exact %>%
    filter(pbdb_present == TRUE & now_present == TRUE) %>%
    nrow()

  summ_table <- data.frame("pbdb_n" = c(summ1, summ3, NA),
                           "pbdb_%" = c(NA, NA, NA),
                           "now_n" = c(summ2, summ3, NA),
                           "now_%" = c(NA, NA, NA))
  summ_table[3,] <- summ_table[1,] - summ_table[2,]
  summ_table[2] <- summ_table[1]/summ1
  summ_table[4] <- summ_table[3]/summ2

  return(summ_table)

}

fuzzymatch_summary_stats <- function(distance_between_strings) {

  summ1 <- bind_rows(filter(match_likely_typos, pbdb_present == TRUE & now_present == TRUE),
                     filter(match_unlikely_typos, pbdb_present == TRUE & now_present == TRUE)) %>%
    filter(string_distance == distance_between_strings)

  summ2 <- length(summ1$string_distance)

  summ3 <- length(unique(summ1$pbdb_accepted_name))

  summ4 <- length(unique(summ1$now_accepted_name))

  summ5 <- summ1 %>%
    group_by(pbdb_accepted_name) %>%
    summarise(n = n()) %>%
    pull(n) %>%
    mean()

  summ6 <- summ1 %>%
    group_by(now_accepted_name) %>%
    summarise(n = n()) %>%
    pull(n) %>%
    mean()

  return(c("Distância entre as strings" = distance_between_strings,
           "Número de pares sugeridos" = summ2,
           "Número de pbdb_accepted_name's distintos que aparecem nos pares" = summ3,
           "Número de now_accepted_name's distintos que aparecem nos pares" = summ4,
           "Média do número de now_accepted_name's propostos para cada pbdb_accepted_name" = summ5,
           "Média do número de pbdb_accepted_name's propostos para cada now_accepted_name" = summ6))
}

# Algumas estatísticas sobre as tabelas
match_exact_summary_stats <- match_exact_summary_stats()

fuzzymatch_summary_stats(1)
fuzzymatch_summary_stats(2)
fuzzymatch_summary_stats(3)
fuzzymatch_summary_stats(4)

View(rbind(fuzzymatch_summary_stats(1),
      fuzzymatch_summary_stats(2)
))

# As tabelas de nomes
match_exact_found #Correspondências exatas
match_exact_failed #Nomes sem correspondências exatas
match_likely_typos_found #Correspondências com 1 ou 2 letras diferentes (derivada de match_exact_failed)
match_likely_typos_failed #Nomes que não tiveram nenhuma correspondência com um máximo de 2 letras diferentes.
match_unlikely_typos_found #Correspondências com 3 ou 4 letras diferentes para checagem (match_likely_typos_failed)



