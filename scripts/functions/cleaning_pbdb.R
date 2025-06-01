pbdb.simple_replace_strings <- function(data = pbdbTaxa.data,
                                        replace_list) {

  if(!is.data.frame(data)){stop()}
  dt_tbl <- data

  check_replacementList(data = dt_tbl, object = replace_list)

  dt_tbl <- dt_tbl %>%
    mutate(across(all_of(replace_list$cols),
                  ~str_squish(.)),
           across(all_of(replace_list$cols),
                  ~str_replace_all(., replace_list$old_new_pairs)))

  return(dt_tbl)

}

pbdb.format_attribution <- function(data = pbdbTaxa.data,
                                    attr_col = select.taxon_attribution_column,
                                    taxon_name_col = select.taxon_name_column,
                                    accepted_name_col = select.accepted_name_column,
                                    create_accepted_attribution_cols = TRUE,
                                    author_name_letter_case = c("upper_case", "lower_case", "capitalized")){

  if (!is.data.frame(data)) {stop('"data" must be a dataframe')}
  if (!attr_col %in% colnames(data)) {stop("Column ", attr_col, " not found")}
  match.arg(author_name_letter_case)
  check_trueOrFalse(create_accepted_attribution_cols)

  dtf <- data %>%
    mutate(taxon_author_parentheses = case_when(
      str_detect(eval(parse(text = attr_col)),
                 pattern = "\\(") ~ TRUE,
      !str_detect(eval(parse(text = attr_col)),
                  pattern = "\\(") ~ FALSE)) %>%
    mutate(attr_format = str_squish(str_replace_all(eval(parse(text = attr_col)),
                                                    pattern = "\\(|\\)|,|\\.",
                                                    replacement = ""))) %>%
    mutate(attr_format = case_when(
      str_detect(attr_format,
                 pattern = "[:digit:]") ~ stringi::stri_replace_last_fixed(attr_format,
                                                                           pattern = " ",
                                                                           replacement = "~"),
      !str_detect(attr_format,
                  pattern = "[:alpha:]") ~ "No-author-info",
      is.na(attr_format) ~ "No-author-info",
      TRUE ~ attr_format
    )) %>%
    separate_wider_delim(attr_format,
                         names = c("taxon_author_name",
                                   "taxon_author_year"),
                         delim = "~",
                         too_few = "align_start") %>%
    mutate(taxon_author_name = stringi:: stri_trans_general(taxon_author_name,
                                                            "Latin-ASCII"),
           taxon_author_name = toupper(taxon_author_name),
           taxon_author_name = str_replace_all(taxon_author_name,
                                               pattern = " ET AL",
                                               replacement = "_ET-AL"),
           taxon_author_name = str_replace_all(taxon_author_name,
                                               pattern = " AND ",
                                               replacement = "_"),
           taxon_author_name = str_replace_all(taxon_author_name,
                                               pattern = " AND_",
                                               replacement = "_"),
           taxon_author_name = str_replace_all(taxon_author_name,
                                               pattern = " ",
                                               replacement = "-"))

  if(author_name_letter_case == "lower_case"){
    dtf$taxon_author_name <- str_to_lower(dtf$taxon_author_name)
  }

  if(author_name_letter_case == "capitalized"){
    dtf$taxon_author_name <- str_replace_all(dtf$taxon_author_name, pattern = "_", replacement = " ")
    dtf$taxon_author_name <- str_to_title(dtf$taxon_author_name)
    dtf$taxon_author_name <- str_replace_all(dtf$taxon_author_name, pattern = "Et-Al", replacement = "et-al")
    dtf$taxon_author_name <- str_replace_all(dtf$taxon_author_name, pattern = " ", replacement = "_")
  }


  if(create_accepted_attribution_cols == TRUE){
    if (!taxon_name_col %in% colnames(data)) {
      warning("Column ", taxon_name_col, " not found, returning a dataframe without accepted attribution columns.")
      return(dtf)
    }
    if (!accepted_name_col %in% colnames(data)) {
      warning("Column ", accepted_name_col, " not found, returning a dataframe without accepted attribution columns.")
      return(dtf)
    }

    accepted_authors <- dtf %>%
      filter(taxon_name == accepted_name) %>%
      select(accepted_name,
             taxon_author_parentheses,
             taxon_author_name,
             taxon_author_year) %>%
      distinct()

    colnames(accepted_authors) <- str_replace_all(colnames(accepted_authors),
                                                  pattern = "taxon",
                                                  replacement = "accepted")

    dtf <- full_join(dtf, accepted_authors) %>%
      distinct()
    #Ele continua criado duplicates aqui, checar o caso de Arctocephalus forsteri

  }

  return(dtf)

}

pbdb.format_alternative_synonym_from_genus <- function(data = pbdbTaxa.data,
                                                       accepted_name_col = select.accepted_name_column,
                                                       genus_col = select.genus_column,
                                                       delim = set.separator_to_fix$old_new_pairs){

  if (!is.data.frame(data)) {stop('"data" must be a dataframe')}
  if (!accepted_name_col %in% colnames(data)) {stop("Column ", accepted_name_col, " not found")}
  if (!genus_col %in% colnames(data)) {stop("Column ", genus_col, " not found")}
  #Adicionar uma checagem da presença do delim nos nomes, se não ele dá um error críptico no separate_wider_delim

  dtf <- data %>%
    select(!!accepted_name_col,
           !!genus_col) %>%
    distinct() %>%
    mutate(accepted_name_copy = eval(parse(text = accepted_name_col))) %>%
    separate_wider_delim(accepted_name_copy,
                         delim = delim,
                         names = c("accepted_genus",
                                   "rest"),
                         too_many = "merge") %>%
    mutate(alternative_synonym = paste(eval(parse(text = genus_col)),
                                       rest ,
                                       sep = delim),
           genus_difference = eval(parse(text = accepted_name_col)) != alternative_synonym) %>%
    select(-accepted_genus,
           -rest) %>%
    inner_join(data,.)

  return(dtf)
}
