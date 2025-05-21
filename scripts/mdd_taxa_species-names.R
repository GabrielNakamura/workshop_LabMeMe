library(here) #Operating system agnostic file paths with .Rproj
library(stringr) #String manipulation
library(dplyr) #Data wrangling with tidyverse syntax
library(tidytable) #Brings the data.table package's fast operations to dplyr syntax
#library(tidylog) #Optional, prints detailed information about changes in the dataframe during data wrangling with tidyverse functions

mdd_data <- read.csv(here("data",
                          "raw",
                          "mdd_taxa_v2.1_raw.csv"),
                     na.strings = "")

#### Arguments ####

select_group_name <- "Carnivora" #Uncomment bellow to check available groups
  #available_orders <- unique(mdd_data$MDD_order)
  #available_families <- distinct(select(mdd_data, MDD_order, MDD_family))
  #available_genera <- distinct(select(mdd_data, MDD_order, MDD_family, MDD_genus))
select_original_rank <- c("species", "synonym_species") #Uncomment bellow to check available original ranks
  #available_ranks <- unique(mdd_data$MDD_original_rank)
select_cols_to_keep <- c("MDD_order",
                         "MDD_family",
                         "MDD_genus",
                         "MDD_species",
                         "MDD_original_combination",
                         "MDD_author",
                         "MDD_year",
                         "MDD_authority_parentheses") #Uncomment bellow to check available columns
  #colnames(MDD_data)
select_species_col <- "MDD_species"
select_group_col <- "MDD_order"
select_original_rank_col <- "MDD_original_rank"
select_original_combination_col <- "MDD_original_combination"

#### Filtering data ####

filter_mdd_synonyms <- function(
    data = mdd_data,
    species_col = select_species_col,
    group_col = select_group_col,
    group_name = select_group_name,
    original_rank_col = select_original_rank_col,
    original_rank_name = select_original_rank,
    columns_to_keep = select_cols_to_keep
){
  colnames <- colnames(data)

  if (!species_col %in% colnames) {stop("\n `species_col` not found. Please provide a valid column name.\n Tip: The standard MDD synonyms table since v2.0 has MDD_species.")}
  if (!group_col %in% colnames & group_col != "none") {stop("\n `group_col` not found. Please provide one valid column name or 'none'.\n Tip: The standard MDD synonyms table since v2.0 usually has MDD_order, MDD_family and MDD_genus.")}
  if (!original_rank_col %in% colnames) {stop("\n `original_rank_col` not found. Please provide one valid column name or 'none'.\n Tip: The standard MDD synonyms table since v2.0 usually has MDD_original_rank.")}
  if (all(is.na(original_rank_name)) | is.null(original_rank_name)) {stop("\n Please provide at least one valid original rank. \n Tip: Check spelling, or the available options on your data with some code like the following `unique(mdd_data$MDD_original_rank)`.")}

  mdd_data <- filter(data, str_detect(string = eval(parse(text = species_col)), pattern = "incertae_sedis", negate = TRUE))

  if (group_col != "none"){
    mdd_data <- filter(mdd_data, eval(parse(text = group_col)) %in% group_name)
  }

  if(nrow(mdd_data) == 0){
    warning("\n `group_name` not found. Returning an empty dataframe. \n Tip: Check if the column provided to `group_col` is the correct one for your desired group.")
    return(mdd_data)
  }

  mdd_data <- filter(mdd_data, eval(parse(text = original_rank_col)) %in% original_rank_name)

  if(nrow(mdd_data) == 0){
    warning("\n None of the strings provided for `original_rank_name` were found. Returning an empty dataframe.")
    return(mdd_data)
  }

  mdd_data <- mdd_data %>%
    select(all_of(columns_to_keep)) %>%
    distinct()

  return(mdd_data)

}

mdd_data <- filter_mdd_synonyms()

#### Original rank == c("species", "synonym_species") ####

format_mdd_synonyms <- function(
    data = mdd_data,
    species_col = select_species_col,
    synonyms_col = select_original_combination_col
){
  colnames <- colnames(data)

  if (!species_col %in% colnames) {stop("\n `species_col` not found. Please provide a valid column name.\n Tip: The standard MDD synonyms table since v2.0 has MDD_species.")}
  if (!synonyms_col %in% colnames) {stop("\n `synonyms_col` not found. Please provide a valid column name.\n Tip: The standard MDD synonyms table since v2.0 has MDD_original_combination")}

  mdd_data <- data %>%
    mutate(!!select_species_col := str_replace_all(.data[[select_species_col]],
                                                   pattern = " ",
                                                   replacement = "_"),
           !!select_original_combination_col := str_replace_all(.data[[select_original_combination_col]],
                                                                c("-" = "",
                                                                  " le | Le " = " le",
                                                                  " von | Von " = " von",
                                                                  " van | Van " = " van",
                                                                  '"' = '')),
           !!select_original_combination_col := stringi::stri_trans_general(.data[[select_original_combination_col]],
                                                                            id = "Latin-ASCII"))

  mdd_syns_manual <- mdd_data %>%
    filter(str_detect(.data[[select_original_combination_col]],
                      pattern = "\\.|\\?| [:alpha:]* [:upper:]|[:digit:]| forma | varietas | Varietas |\\]$|\\[sic\\]|\\[.* .*\\]"))

  mdd_data <- mdd_data %>%
    anti_join(mdd_syns_manual) %>%
    mutate(!!select_original_combination_col := str_replace_all(.data[[select_original_combination_col]],
                                                                c("\\[" = "\\(",
                                                                  "\\]\\.|\\]" = "\\)")))

  mdd_data_epithets <- mdd_data %>%
    filter(str_detect(.data[[select_original_combination_col]],
                      pattern = " |\\.|\\([:upper:]",
                      negate = TRUE)) %>%
    mutate(!!select_original_combination_col := str_to_lower(.data[[select_original_combination_col]]),
           MDD_actual_rank = "epithet")

  mdd_data_subgen <- mdd_data %>%
    filter(str_detect(.data[[select_original_combination_col]],
                      pattern = "\\([:upper:]") &
             str_detect(.data[[select_original_combination_col]],
                        pattern = "\\.|\\?| [:alpha:]* [:upper:]",
                        negate = TRUE)) %>%
    mutate(!!select_original_combination_col := str_replace(str_replace(str_replace(.data[[select_original_combination_col]],
                                                                                    pattern = " ",
                                                                                    replacement = "_"),
                                                                        pattern = " ",
                                                                        replacement = "~"),
                                                            pattern = " ",
                                                            replacement = "-")) %>%
    separate_wider_delim(cols = .data[[select_original_combination_col]],
                         names = c("col1", "col2"),
                         delim = "~",
                         too_few = "align_start",
                         too_many = "merge") %>%
    mutate(col2 = str_to_lower(col2)) %>%
    unite(col = !!select_original_combination_col,
          col1,
          col2,
          sep = "~",
          remove = TRUE,
          na.rm = TRUE)

  mdd_data_species <- mdd_data %>%
    filter(str_detect(.data[[select_original_combination_col]],
                      pattern = " ") &
             str_detect(.data[[select_original_combination_col]],
                        pattern = "\\.|\\([:upper:]|\\?| [:alpha:]* [:upper:]|[:digit:]",
                        negate = TRUE)) %>%
    mutate(!!select_original_combination_col := str_replace(str_replace(.data[[select_original_combination_col]],
                                                                        pattern = " ",
                                                                        replacement = "_"),
                                                            pattern = " ",
                                                            replacement = "-")) %>%
    bind_rows(., filter(mdd_data_subgen,
                        str_detect(.data[[select_original_combination_col]],
                                   pattern = "~",
                                   negate = TRUE))) %>%
    mutate(!!select_original_combination_col := str_to_sentence(str_replace_all(.data[[select_original_combination_col]],
                                                                                pattern = "\\(|\\)",
                                                                                replacement = "")))

  mdd_data_subsp <- mdd_data_species %>%
    filter(str_detect(.data[[select_original_combination_col]],
                      pattern = "-")) %>%
    mutate(!!select_original_combination_col := str_replace_all(.data[[select_original_combination_col]],
                                                                c(" " = "_",
                                                                  "\\(|\\)" = ""))) %>%
    bind_rows(filter(mdd_data_subgen,
                     str_detect(.data[[select_original_combination_col]],
                                pattern = "-"))) %>%
    mutate(!!select_original_combination_col := str_replace_all(.data[[select_original_combination_col]],
                                                                pattern = "-|~",
                                                                replacement = "_"),
           MDD_actual_rank = case_when(
             str_detect(.data[[select_original_combination_col]], pattern = "\\(") ~ "subspecies (with subgenus)",
             TRUE ~ "subspecies"))

  mdd_data_subgen <- mdd_data_subgen %>%
    filter(!str_detect(.data[[select_original_combination_col]],
                       pattern = "-") &
             str_detect(.data[[select_original_combination_col]],
                        pattern = "~")) %>%
    mutate(!!select_original_combination_col := str_replace_all(.data[[select_original_combination_col]],
                                                                pattern = "~",
                                                                replacement = "_"),
           MDD_actual_rank = "species (with subgenus)")

  mdd_data_species <- mdd_data_species %>%
    filter(!str_detect(.data[[select_original_combination_col]],
                       pattern = "-")) %>%
    mutate(MDD_actual_rank = "species")

  mdd_syns_final <- bind_rows(mdd_data_epithets,
                              mdd_data_species,
                              mdd_data_subgen,
                              mdd_data_subsp)

  return(list(mdd_syns_final = mdd_syns_final,
              mdd_syns_manual = mdd_syns_manual))

}

mdd_syns <- format_mdd_synonyms() #Generates a list with 2 items: 1) mdd_syns_final = correctly formatted names; 2) mdd_syns_manual = names in need of manual cleaning

#### Write CSV ####

write.csv(mdd_syns$mdd_syns_final,
          file = here("data",
                      "processed",
                      "mdd_taxa_species-names_v2.1.csv"),
          row.names = FALSE)
