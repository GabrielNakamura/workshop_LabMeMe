now.fix_indet <- function(col_name) {
  case_when(
    is.na(col_name) | col_name == "" | col_name == "\\N" | col_name == "sp." | col_name == "Gen." | col_name == "gen." | col_name == "incertae sedis" ~ "indet.",
    TRUE ~ col_name
  )
}

now.create_species_name_col <- function(data = now.taxa,
                                        new_col_name = "now_accepted_name",
                                        species_col = select.species_col,
                                        genus_col = select.genus_col,
                                        sep = select.new_name_separator){

  dtf <- data %>%
    mutate(new_col = paste(eval(parse(text = genus_col)), eval(parse(text = species_col)), sep = sep),
           new_col = str_replace(new_col,
                                 pattern = "_indet\\.",
                                 replacement = "")) %>%
    rename(!!new_col_name := new_col)

  return(dtf)

}


now.create_accepted_rank_col <- function(data = now.taxa,
                                         new_col_name = "now_accepted_rank",
                                         species_col = select.species_col,
                                         genus_col = select.genus_col,
                                         subfamily_col = select.subfamily_col,
                                         family_col = select.family_col,
                                         suborder_superfamily_col = select.suborder_superfamily_col,
                                         order_col = select.order_col,
                                         subclass_superorder_col = select.subclass_superorder_col
){

  dtf <- data %>%
    mutate(new_col = case_when(
      eval(parse(text = species_col)) != "indet." ~ "species",
      eval(parse(text = species_col)) == "indet." & eval(parse(text = genus_col)) != "indet." ~ "genus",
      eval(parse(text = species_col))  == "indet." & eval(parse(text = genus_col)) == "indet." ~ "Other"
    ),
    new_col = case_when(
      new_col == "Other" & eval(parse(text = subfamily_col)) != "indet." ~ "subfamily",
      TRUE ~ new_col
    ),
    new_col = case_when(
      new_col == "subfamily" & str_detect(eval(parse(text = subfamily_col)), pattern = "ini$") ~ "tribe",
      TRUE ~ new_col
    ),
    new_col = case_when(
      new_col == "Other" & eval(parse(text = family_col)) != "indet." ~ "family",
      TRUE ~ new_col
    ),
    new_col = case_when(
      new_col == "Other" & eval(parse(text = suborder_superfamily_col)) != "indet." ~ "suborder",
      TRUE ~ new_col
    ),
    new_col = case_when(
      new_col == "suborder" & str_detect(eval(parse(text = suborder_superfamily_col)), pattern = "oidea$") ~ "superfamily",
      TRUE ~ new_col
    ),
    new_col = case_when(
      new_col == "Other" & eval(parse(text = order_col)) != "indet." ~ "order",
      TRUE ~ new_col
    ),
    new_col = case_when(
      new_col == "Other" & eval(parse(text = subclass_superorder_col)) != "indet." ~ "subclass_superorder",
      TRUE ~ new_col
    )) %>%
    rename(!!new_col_name := new_col)

  return(dtf)

}

now.summarise_ages_by_id <- function(data = now.occurrences.data,
                                     taxon_id = select.taxon_id_col,
                                     min_age_col = select.min_age_col,
                                     max_age_col = select.max_age_col){

  dtf <- data %>%
    select(all_of(c(taxon_id,
                    min_age_col,
                    max_age_col))) %>%
    distinct()%>%
    group_by_at(taxon_id) %>%
    summarise(!!min_age_col := min(eval(parse(text = min_age_col))),
              !!max_age_col := max(eval(parse(text = max_age_col))))

  return(dtf)
}
