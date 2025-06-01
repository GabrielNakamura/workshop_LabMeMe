# _Info -------------------------------------------------------------------
##
## Title: pbdb_taxa_species-names.R
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
library(tidylog) #Optional, prints detailed information about changes in the dataframe during data wrangling with tidyverse functions

# __File paths ------------------------------------------------------------

pbdbTaxa.path <- here("data",
                      "raw",
                      "pbdb_taxa_2025-06-01_raw.csv")

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

pbdbTaxa.data <- read.csv(pbdbTaxa.path,
                          na.strings = "",
                          skip = 17)

# __Cleanup ---------------------------------------------------------------

rm(pbdbTaxa.path,
   functions.1.path,
   functions.2.path)
gc()

# _Main -------------------------------------------------------------------
# __Setting up function arguments -----------------------------------------

#colnames(pbdbTaxa.data) #Uncomment to check available columns
select.cols_to_keep <- c("flags",
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

#unique(pbdbTaxa.data$accepted_rank) #Uncomment to check available ranks
select.accepted_rank_name <- c(
  #"subspecies",
  "species"
)

select.accepted_name_column <- "accepted_name"
select.taxon_name_column <- "taxon_name"
select.taxon_attribution_column <- "taxon_attr"
select.genus_column <- "genus"

#String replacement guides
set.separator_to_fix <- list(cols = c("taxon_name",
                                      "accepted_name"),
                             old_new_pairs = c(" " = "_"))

set.typos_to_fix <- list(cols = "taxon_attr",
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

# __Filtering data --------------------------------------------------------

pbdbTaxa.data <- pbdbTaxa.data %>%
  filter(accepted_rank %in% select.accepted_rank_name) %>%
  select(all_of(select.cols_to_keep)) %>%
  distinct()

gc()

# __Cleaning and wrangling data -------------------------------------------
# ___Replacing strings ----------------------------------------------------

pbdbTaxa.data <- pbdb.simple_replace_strings(replace_list = set.separator_to_fix)
pbdbTaxa.data <- pbdb.simple_replace_strings(replace_list = set.typos_to_fix)

gc()

# ___Formatting author names ----------------------------------------------

pbdbTaxa.data <- pbdb.format_attribution(author_name_letter_case = "capitalized")

gc()


# ___Create alternative synonym from accepted name and genus column -------

pbdbTaxa.data <- pbdb.format_alternative_synonym_from_genus()

gc()

# ___Adding suffix to column names ----------------------------------------

colnames(pbdbTaxa.data) <- paste0("pbdb_", colnames(pbdbTaxa.data))


# __Saving data -----------------------------------------------------------

write.csv(x = pbdbTaxa.data,
          file = here("data",
                      "processed",
                      "pbdb_taxa_species-names_2025-06-01.csv"),
          row.names = FALSE)

# _Cleanup ----------------------------------------------------------------

#rm(list = ls())
#gc()
