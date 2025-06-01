# _Info -------------------------------------------------------------------
##
## Title: pbdb_occs_species-names.R
## Purpose:
##
## Author: Thaís G. P. Faria
## Github: https://github.com/Thais-Faria
## Date created: 2025-06-01
## Copyright (c) Thaís G. P. Faria, 2025
## License: GNU General Public License version 3
##
# _Set up -----------------------------------------------------------------
# __Packages --------------------------------------------------------------

library(here) #Operating system agnostic file paths with .Rproj
library(stringr) #String manipulation
library(dplyr) #Data wrangling with tidyverse syntax
library(tidytable) #Brings the data.table package's fast operations to dplyr syntax
#library(tidylog) #Optional, prints detailed information about changes in the dataframe during data wrangling with tidyverse functions

# __File paths ------------------------------------------------------------

pbdbOccs.path <- here("data",
                      "raw",
                      "pbdb_occs_2025-06-01_raw.csv")

functions.1.path <- here("scripts",
                         "functions",
                         "checking_general.R")

functions.2.path <- here("scripts",
                         "functions",
                         "cleaning_pbdb.R")

# __Loading functions -----------------------------------------------------

source(functions.1.path)
source(functions.2.path)

# __Loading data ----------------------------------------------------------

pbdbOccs.data <- read.csv(pbdbOccs.path,
                          skip = 19,
                          na.strings = "")

# __Cleanup ---------------------------------------------------------------

rm(pbdbOccs.path,
   functions.1.path,
   functions.2.path)
gc()

# _Main -------------------------------------------------------------------
# __Setting up function arguments -----------------------------------------

#colnames(pbdbOccs.data) #Uncomment to check available columns
select.cols_to_keep <- c("identified_no",
                         "accepted_no",
                         "identified_name",
                         "difference",
                         "accepted_name",
                         "identified_rank",
                         "accepted_rank",
                         "accepted_attr",
                         "phylum",
                         "class",
                         "order",
                         "family",
                         "genus",
                         "subgenus_name",
                         "subgenus_reso",
                         "primary_reso",
                         "species_reso",
                         "phylum_no",
                         "class_no",
                         "order_no",
                         "family_no",
                         "genus_no",
                         "subgenus_no")

#unique(pbdbOccs.data$accepted_rank) #Uncomment to check available ranks
select.accepted_rank_name <- c(
  #"subspecies",
  "species"
)

select.accepted_name_column <- "accepted_name"
select.genus_column <- "genus"

#String replacement guides
set.separator_to_fix <- list(cols = c("identified_name",
                                      "accepted_name",
                                      "genus"),
                             old_new_pairs = c(" " = "_"))

# __Filtering data --------------------------------------------------------

pbdbOccs.data <- pbdbOccs.data %>%
  filter(accepted_rank %in% select.accepted_rank_name) %>%
  select(all_of(select.cols_to_keep)) %>%
  distinct()

gc()

# __Cleaning and wrangling data -------------------------------------------
# ___Replacing strings ----------------------------------------------------

pbdbOccs.data <- pbdb.simple_replace_strings(data = pbdbOccs.data,
                                             replace_list = set.separator_to_fix)

pbdbOccs.data <- pbdbOccs.data %>%
  mutate(family_no = as.integer(na_if(family_no, "NF")))

gc()

# ___Formatting author names ----------------------------------------------

# Notes: A coluna accepted_attr deveria conter as autoridades taxonômicas dos
#        nomes aceitos, mas veio vazia por alguma razão. Caso os dados sejam
#        atualizados no futuro, esse script poderá ser expandido para formatar
#        a coluna de autores e anos.

# ___Create alternative synonym from accepted name and genus column -------

pbdbOccs.data <- pbdb.format_internal_synonym_from_genus(data = pbdbOccs.data)

gc()

# ___Adding suffix to column names ----------------------------------------

colnames(pbdbOccs.data) <- paste0("pbdb_", colnames(pbdbOccs.data))

# __Saving data -----------------------------------------------------------

write.csv(x = pbdbOccs.data,
          file = here("data",
                      "processed",
                      "pbdb_occs_species-names_2025-06-01.csv"),
          row.names = FALSE)

# _Cleanup ----------------------------------------------------------------

#rm(list = ls())
#gc()
