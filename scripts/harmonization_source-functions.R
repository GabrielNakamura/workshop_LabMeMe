separate_name_parts <- function(data,
                                name_col,
                                base_suffix,
                                delim) {

  #Separate genus part from the rest
  dtf <- data %>%
    separate_wider_delim(cols = name_col,
                         delim = delim,
                         names = c("suffix_genus", "suffix_other"),
                         too_many = "merge",
                         too_few = "align_start",
                         cols_remove = FALSE)

  #Separate subgenus part from the rest and clean markers
  dtf <- dtf %>%
    separate_wider_delim(cols = "suffix_other",
                         delim = ")",
                         names = c("suffix_subgenus", "suffix_other"),
                         too_few = "align_end") %>%
    mutate(suffix_subgenus = str_replace_all(suffix_subgenus,
                                             pattern = "\\(",
                                             replacement = ""),
           suffix_other = str_replace(suffix_other,
                                      pattern = paste0("^", delim),
                                      replacement = ""))

  #Separate species and subspecies epithets
  dtf <- dtf %>%
    separate_wider_delim(cols = "suffix_other",
                         delim = delim,
                         names = c("suffix_speciesEpiteth",
                                   "suffix_subspeciesEpiteth"),
                         too_few = "align_start")

  colnames(dtf) <- str_replace_all(colnames(dtf),
                                   pattern = "suffix",
                                   replacement = base_suffix)

  return(dtf)
}

delete_subgenus <- function(data,
                            name_col,
                            delim) {

  dtf <- data %>%
    filter(str_detect(eval(parse(text = name_col)), pattern = "\\(")) %>%
    mutate(col_copy = eval(parse(text = name_col)),
           col_copy = str_replace_all(col_copy,
                                      pattern = paste0(delim,
                                                       "\\(.*\\)"),
                                      replacement = ""))
  colnames(dtf) <- str_replace(colnames(dtf),
                               pattern = "col_copy",
                               replacement = paste0(name_col, "_delSubgen"))
  if (nrow(dtf) == 0) {
    warning("No taxonomic names containing a subgenus portion enclosed in parentheses found, returning an empty dataframe")
  }

  return(dtf)
}

swap_subgenus <- function(data,
                          name_col,
                          delim) {


  dtf <- data %>%
    filter(str_detect(eval(parse(text = name_col)), pattern = "\\("))

  dtf <- dtf %>%
    mutate(col_copy = eval(parse(text = name_col)),
           col_copy = str_replace_all(col_copy,
                                      c(".*\\(" = "",
                                        "\\)" = "")))

  colnames(dtf) <- str_replace(colnames(dtf),
                               pattern = "col_copy",
                               replacement = paste0(name_col, "_swapSubgen"))

  if (nrow(dtf) == 0) {
    warning("No taxonomic names containing a subgenus portion enclosed in parentheses found, returning an empty dataframe")
  }
  return(dtf)
}

harmonize_exact_match <- function(base1_dtf,
                                  base2_dtf,
                                  base1_suffix = b1_suffix,
                                  base2_suffix = b2_suffix,
                                  base1_col = b1_col,
                                  base2_col = b2_col) {

  by_cols <- base2_col
  names(by_cols) <- base1_col

  select_cols <- c(base1_col, base2_col)

  match_exact <- full_join(base1_dtf, base2_dtf,
                           by = by_cols,
                           keep = TRUE) %>%
    mutate(base1_present = !is.na(eval(parse(text = base1_col))),
           base2_present = !is.na(eval(parse(text = base2_col)))) %>%
    distinct()

  exact1_found <- match_exact  %>%
    filter(base1_present == TRUE & base2_present == TRUE) %>%
    mutate(string_distance = 0,
           match_notes = "exact match") %>%
    select(all_of(select_cols),
           string_distance,
           match_notes)

  base1_failed_match <- match_exact %>%
    filter(base2_present == FALSE) %>%
    select(any_of(colnames(base1_dtf)))

  base2_failed_match <- match_exact %>%
    filter(base1_present == FALSE) %>%
    select(any_of(colnames(base2_dtf)))

  summ1 <- match_exact[[base1_col]] %>%
    na.omit() %>%
    unique() %>%
    length()

  summ2 <- match_exact[[base2_col]] %>%
    na.omit() %>%
    unique() %>%
    length()

  summ3 <- match_exact %>%
    filter(base1_present == TRUE & base2_present == TRUE) %>%
    nrow()

  summ_table <- data.frame(base1_count = c(summ1, summ3, NA),
                           base1_proportion = c(NA, NA, NA),
                           base2_count = c(summ2, summ3, NA),
                           base2_proportion = c(NA, NA, NA))
  summ_table[3,] <- summ_table[1,] - summ_table[2,]
  summ_table[2] <- summ_table[1]/summ1
  summ_table[4] <- summ_table[3]/summ2

  colnames(summ_table) <- str_replace_all(colnames(summ_table),
                                          c("base1" = base1_suffix,
                                            "base2" = base2_suffix))

  row.names(summ_table) <- c("unique_names",
                             "matched_names",
                             "unmatched_names")

  result <- list(exact_summary = summ_table,
                 exact_found = exact1_found,
                 base1_exact_failed = base1_failed_match,
                 base2_exact_failed = base2_failed_match)

  names(result) <- str_replace_all(names(result),
                                   c("base1" = base1_suffix,
                                     "base2" = base2_suffix))

  return(result)

}

harmonize_fuzzy_match <- function(base1_dtf,
                                  base2_dtf,
                                  base1_suffix = b1_suffix,
                                  base2_suffix = b2_suffix,
                                  base1_col = b1_col,
                                  base2_col = b2_col,
                                  min_dist,
                                  max_dist,
                                  delim) {

  by_cols <- base2_col
  names(by_cols) <- base1_col

  select_cols <- c(base1_col, base2_col)

  base1_names <- base1_dtf[base1_col]
  base2_names <- base2_dtf[base2_col]

  match_fuzzy <- stringdist_full_join(base1_names, base2_names,
                                      by = by_cols,
                                      distance_col = "string_distance",
                                      max_dist = max_dist) %>%
    mutate(base1_present = !is.na(eval(parse(text = base1_col))),
           base2_present = !is.na(eval(parse(text = base2_col)))) %>%
    distinct()

  fuzzy_match_found <- match_fuzzy %>%
    filter(string_distance >= min_dist & string_distance <= max_dist) %>%
    select(-base1_present,
           -base2_present)

  fuzzy_match_found <- separate_name_parts(fuzzy_match_found,
                                           name_col = base1_col,
                                           base_suffix = "base1",
                                           delim = delim)

  fuzzy_match_found <- separate_name_parts(fuzzy_match_found,
                                           name_col = base2_col,
                                           base_suffix = "base2",
                                           delim = delim)

  fuzzy_match_found <- fuzzy_match_found %>%
    mutate(genus_dif = case_when(stringdist::stringdist(base1_genus, base2_genus) > 0 ~ "different spelling on genus"),
           subgenus_dif = case_when(stringdist::stringdist(base1_subgenus, base2_subgenus) > 0 ~ "different spelling on subgenus"),
           speciesEpi_dif = case_when(stringdist::stringdist(base1_speciesEpiteth, base2_speciesEpiteth) > 0 ~ "different spelling on species epiteth"),
           subspeciesEpi_dif = case_when(stringdist::stringdist(base1_subspeciesEpiteth, base2_subspeciesEpiteth) > 0 ~ "different spelling on subspecies epiteth"))

  fuzzy_match_found <- fuzzy_match_found %>%
    mutate(match_notes = paste(genus_dif,
                               subgenus_dif,
                               speciesEpi_dif,
                               subspeciesEpi_dif,
                               sep = ","),
           match_notes = str_replace_all(match_notes,
                                         pattern = "NA,|,NA",
                                         replacement = "")) %>%
    select(all_of(select_cols),
           string_distance,
           match_notes) %>%
    mutate(match_notes = case_when(string_distance == 0 ~ "exact match",
                                   TRUE ~ match_notes))

  base1_dist_failed_match <- match_fuzzy %>%
    filter(base2_present == FALSE) %>%
    select(any_of(colnames(base1_dtf))) %>%
    left_join(base1_dtf, by = base1_col)

  base2_dist_failed_match <- match_fuzzy %>%
    filter(base1_present == FALSE) %>%
    select(any_of(colnames(base2_dtf))) %>%
    left_join(base2_dtf, by = base2_col)

  result <- list(dist_match_found = fuzzy_match_found,
                 base1_dist_match_failed = base1_dist_failed_match,
                 base2_dist_match_failed = base2_dist_failed_match)

  names(result) <- str_replace_all(names(result),
                                   c("base1" = base1_suffix,
                                     "base2" = base2_suffix))

  names(result) <- str_replace_all(names(result),
                                   pattern = "dist",
                                   replacement = paste0("min",
                                                        min_dist,
                                                        "max",
                                                        max_dist))

  return(result)

}

add_to_synonymy <- function(synonymy_dtf,
                            eval_dtf){

  dtf <- distinct(bind_rows(synonymy_dtf, select(filter(eval_dtf, accepted_match == TRUE), -accepted_match)))

  return(dtf)
}
