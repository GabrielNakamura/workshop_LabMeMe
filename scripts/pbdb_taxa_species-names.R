#### Loading packages ####

library(here) #Operating system agnostic file paths with .Rproj
library(stringr) #String manipulation
library(dplyr) #Data wrangling with tidyverse syntax
library(tidytable) #Brings the data.table package's fast operations to dplyr syntax
#library(tidylog) #Optional, prints detailed information about changes in the dataframe during data wrangling with tidyverse functions

#### Loading files ####

file_path <- here("data",
                  "raw",
                  "pbdb_taxa_2025-05-06_raw.csv")

pbdb_data <- read.csv(file_path,
                      na.strings = "",
                      skip = 17)

#### Loading support functions ####

check_trueOrFalse <- function(arg){
  msg <- paste0("`", deparse(substitute(arg)), "` must be either TRUE or FALSE")
  if(is.list(unclass(arg))){stop(msg)}
  if(length(arg) != 1){stop(msg)}
  if(!is.logical(arg)){stop(msg)}
  if(is.na(arg)){stop(msg)}
}

check_replacementList <- function(data, object){
  if(is.null(object$cols)){stop()}
  if(is.list(unclass(object$cols))){stop()}
  if(!is.character(object$cols)){stop()}
  cols <- na.omit(object$cols)
  if(!all(cols %in% colnames(data))){stop()}

  if(is.null(object$old_new_pairs)){stop()}
  if(is.list(unclass(object$old_new_pairs))){stop()}
  if(!is.character(object$old_new_pairs)){stop()}
  pairs <- na.omit(object$old_new_pairs)
  if(length(pairs) < 1){stop()}
  pair_names <- names(pairs)
  pair_names <- pair_names[pair_names != ""]
  if(length(pair_names) != length(pairs)){stop()}
}

#### Setting up arguments for functions ####

#colnames(pbdb_data) #Uncomment to check available columns
select_cols_to_keep <- c("flags",
                         "taxon_no",
                         "taxon_rank",
                         "taxon_name",
                         "taxon_attr",
                         "difference",
                         "accepted_no",
                         "accepted_rank",
                         "accepted_name",
                         "parent_name",
                         "immpar_name",
                         "phylum",
                         "class",
                         "order",
                         "family",
                         "genus")

#unique(pbdb_data$accepted_rank) #Uncomment to check available ranks
select_rank_to_filter <- c("subspecies",
                           "species")

select_accepted_name_column <- "accepted_name"
select_taxon_name_column <- "taxon_name"
select_taxon_attribution_column <- "taxon_attr"
select_genus_column <- "genus"

#String replacement guides
separator_to_fix <- list(cols = c("taxon_name",
                                  "accepted_name"),
                         old_new_pairs = c(" " = "_"))

typos_to_fix <- list(cols = "taxon_attr",
                     old_new_pairs = c("Bones" = "Bonis",
                                       "Croziet" = "Croizet",
                                       "Cusafont" = "Crusafont",
                                       "Eherenberg" = "Ehrenberg",
                                       "Ferrusquı" = "Ferrusquia",
                                       "Geoffroy-Saint-Hillaire" = "Geoffroy-Saint-Hilaire",
                                       "Hendy" = "Hendey",
                                       "Jiangzou" = "Jiangzuo",
                                       "Kretsoi" = "Kretzoi",
                                       "Linneaus" = "Linnaeus",
                                       "Nillson" = "Nilsson",
                                       "Pertényi" = "Petényi",
                                       "Pockock" = "Pocock",
                                       "Rossenbüller" = "Rosenmüller",
                                       "Shultz" = "Schultz",
                                       "Solunias" = "Solounias", #https://sjpp.springeropen.com/articles/10.1007/s13358-012-0042-y
                                       "Teilard" = "Teilhard",
                                       "Zimmerman" = "Zimmermann"))

##### Cleaning, formatting and wrangling #####
#### Filtering data ####

pbdb_data <- pbdb_data %>%
  filter(accepted_rank %in% select_rank_to_filter) %>%
  select(all_of(select_cols_to_keep)) %>%
  distinct()

gc()

#### Replacing strings ####

simple_replace_strings <- function(data = pbdb_data,
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

pbdb_data <- simple_replace_strings(replace_list = separator_to_fix)
pbdb_data <- simple_replace_strings(replace_list = typos_to_fix)

gc()

#### Format author names ####

format_pbdb_attribution <- function(data = pbdb_data,
                             attr_col = select_taxon_attribution_column,
                             taxon_name_col = select_taxon_name_column,
                             accepted_name_col = select_accepted_name_column,
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

    dtf <- full_join(dtf, accepted_authors)

  }

  return(dtf)

}

pbdb_data <- format_pbdb_attribution(author_name_letter_case = "capitalized")

gc()

#### Create accepted name with genus column ####

format_alternative_synonym_from_genus <- function(data = pbdb_data,
                                                  accepted_name_col = select_accepted_name_column,
                                                  genus_col = select_genus_column,
                                                  delim = "_"){

  if (!is.data.frame(data)) {stop('"data" must be a dataframe')}
  if (!accepted_name_col %in% colnames(data)) {stop("Column ", accepted_name_col, " not found")}
  if (!genus_col %in% colnames(data)) {stop("Column ", genus_col, " not found")}

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

pbdb_data <- format_alternative_synonym_from_genus()

gc()

#### Save data ####

write.csv(x = pbdb_data,
          file = here("data",
                      "processed",
                      "pbdb_taxa_species-names_2025-05-06.csv"),
          row.names = FALSE)

rm(list = ls())

gc()
