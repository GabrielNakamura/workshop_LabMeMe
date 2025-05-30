mdd.filter_synonyms <- function(
    data = mdd.synonymy.data,
    species_col = select.species_col,
    group_col = select.group_col,
    group_name = select.group_name,
    original_rank_col = select.original_rank_col,
    original_rank_name = select.original_rank_name,
    cols_to_keep = select.cols_to_keep
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
    select(all_of(cols_to_keep)) %>%
    distinct()

  return(mdd_data)

}

mdd.format_synonyms <- function(
    data = mdd.synonymy.data,
    species_col = select.species_col,
    original_combination_col = select.original_combination_col,
    sep = select.new_name_separator
){
  colnames <- colnames(data)

  if (!species_col %in% colnames) {stop("\n `species_col` not found. Please provide a valid column name.\n Tip: The standard MDD synonyms table since v2.0 has MDD_species.")}
  if (!original_combination_col %in% colnames) {stop("\n `original_combination_col` not found. Please provide a valid column name.\n Tip: The standard MDD synonyms table since v2.0 has MDD_original_combination")}

  mdd_data <- data %>%
    mutate(!!species_col := str_replace_all(.data[[species_col]],
                                                   pattern = " ",
                                                   replacement = "_"),
           !!original_combination_col := str_replace_all(.data[[original_combination_col]],
                                                                c("-" = "",
                                                                  " le | Le " = " le",
                                                                  " von | Von " = " von",
                                                                  " van | Van " = " van",
                                                                  '"' = '')),
           !!original_combination_col := stringi::stri_trans_general(.data[[original_combination_col]],
                                                                            id = "Latin-ASCII"))

  mdd_syns_manual <- mdd_data %>%
    filter(str_detect(.data[[original_combination_col]],
                      pattern = "\\.|\\?| [:alpha:]* [:upper:]|[:digit:]| forma | varietas | Varietas |\\]$|\\[sic\\]|\\[.* .*\\]"))

  mdd_data <- mdd_data %>%
    anti_join(mdd_syns_manual) %>%
    mutate(!!original_combination_col := str_replace_all(.data[[original_combination_col]],
                                                                c("\\[" = "\\(",
                                                                  "\\]\\.|\\]" = "\\)")))

  mdd_data_epithets <- mdd_data %>%
    filter(str_detect(.data[[original_combination_col]],
                      pattern = " |\\.|\\([:upper:]",
                      negate = TRUE)) %>%
    mutate(!!original_combination_col := str_to_lower(.data[[original_combination_col]]),
           MDD_actual_rank = "epithet")

  mdd_data_subgen <- mdd_data %>%
    filter(str_detect(.data[[original_combination_col]],
                      pattern = "\\([:upper:]") &
             str_detect(.data[[original_combination_col]],
                        pattern = "\\.|\\?| [:alpha:]* [:upper:]",
                        negate = TRUE)) %>%
    mutate(!!original_combination_col := str_replace(str_replace(str_replace(.data[[original_combination_col]],
                                                                                    pattern = " ",
                                                                                    replacement = "_"),
                                                                        pattern = " ",
                                                                        replacement = "~"),
                                                            pattern = " ",
                                                            replacement = "-")) %>%
    separate_wider_delim(cols = .data[[original_combination_col]],
                         names = c("col1", "col2"),
                         delim = "~",
                         too_few = "align_start",
                         too_many = "merge") %>%
    mutate(col2 = str_to_lower(col2)) %>%
    unite(col = !!original_combination_col,
          col1,
          col2,
          sep = "~",
          remove = TRUE,
          na.rm = TRUE)

  mdd_data_species <- mdd_data %>%
    filter(str_detect(.data[[original_combination_col]],
                      pattern = " ") &
             str_detect(.data[[original_combination_col]],
                        pattern = "\\.|\\([:upper:]|\\?| [:alpha:]* [:upper:]|[:digit:]",
                        negate = TRUE)) %>%
    mutate(!!original_combination_col := str_replace(str_replace(.data[[original_combination_col]],
                                                                        pattern = " ",
                                                                        replacement = "_"),
                                                            pattern = " ",
                                                            replacement = "-")) %>%
    bind_rows(., filter(mdd_data_subgen,
                        str_detect(.data[[original_combination_col]],
                                   pattern = "~",
                                   negate = TRUE))) %>%
    mutate(!!original_combination_col := str_to_sentence(str_replace_all(.data[[original_combination_col]],
                                                                                pattern = "\\(|\\)",
                                                                                replacement = "")))

  mdd_data_subsp <- mdd_data_species %>%
    filter(str_detect(.data[[original_combination_col]],
                      pattern = "-")) %>%
    mutate(!!original_combination_col := str_replace_all(.data[[original_combination_col]],
                                                                c(" " = "_",
                                                                  "\\(|\\)" = ""))) %>%
    bind_rows(filter(mdd_data_subgen,
                     str_detect(.data[[original_combination_col]],
                                pattern = "-"))) %>%
    mutate(!!original_combination_col := str_replace_all(.data[[original_combination_col]],
                                                                pattern = "-|~",
                                                                replacement = "_"),
           MDD_actual_rank = case_when(
             str_detect(.data[[original_combination_col]], pattern = "\\(") ~ "subspecies (with subgenus)",
             TRUE ~ "subspecies"))

  mdd_data_subgen <- mdd_data_subgen %>%
    filter(!str_detect(.data[[original_combination_col]],
                       pattern = "-") &
             str_detect(.data[[original_combination_col]],
                        pattern = "~")) %>%
    mutate(!!original_combination_col := str_replace_all(.data[[original_combination_col]],
                                                                pattern = "~",
                                                                replacement = "_"),
           MDD_actual_rank = "species (with subgenus)")

  mdd_data_species <- mdd_data_species %>%
    filter(!str_detect(.data[[original_combination_col]],
                       pattern = "-")) %>%
    mutate(MDD_actual_rank = "species")

  mdd_syns_final <- bind_rows(mdd_data_epithets,
                              mdd_data_species,
                              mdd_data_subgen,
                              mdd_data_subsp)

  if (sep != "_"){

    mdd_syns_final <- mutate(mdd_syns_final,
                             across(everything(),
                                    ~str_replace_all(., pattern = "_",
                                                     replacement = sep)))

  }

  return(list(final = mdd_syns_final,
              manual = mdd_syns_manual))

}
