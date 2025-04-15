#### Carregando pacotes com controle de versão ####

library(groundhog)
groundhog.day <- "2025-04-01"
groundhog.packages <- c("here",
                        "tidyr",
                        "dplyr",
                        "stringr",
                        "forcats",
                        "tidylog")

groundhog.library(groundhog.packages, 
                  groundhog.day)
rm(groundhog.day,
   groundhog.packages)

#### Limpeza inicial ####
# Carregar arquivo #

now_data <- read.delim(here("data", 
                            "raw",
                            "now_occs_2025-04-03_raw.tsv"),
                       na.strings = "\\N")

now_colnames_vct <- c("now_accepted_sist_subclass_superorder" = "SUBCLASSORSUPERORDER",
                      "now_accepted_sist_order" = "ORDER",
                      "now_accepted_sist_suborder_superfamily" = "SUBORDERORSUPERFAMILY",
                      "now_accepted_sist_family" = "FAMILY",
                      "now_accepted_sist_subfamily_tribe" = "SUBFAMILY",
                      "now_accepted_sist_genus" = "GENUS",
                      "now_accepted_sist_epithet" = "SPECIES")

indet_fix <- function(col_name) {
  case_when(
    col_name == "" | col_name == "\\N" | col_name == "sp." | col_name == "Gen." | col_name == "gen." | col_name == "incertae sedis" ~ "indet.",
    TRUE ~ col_name
  )
}

now_taxa <- now_data %>%
  rename(all_of(now_colnames_vct)) %>%
  select(now_accepted_sist_subclass_superorder,
         now_accepted_sist_order,
         now_accepted_sist_suborder_superfamily,
         now_accepted_sist_family,
         now_accepted_sist_subfamily_tribe,
         now_accepted_sist_genus,
         now_accepted_sist_epithet) %>%
  distinct() %>%
  mutate(across(everything(),
                ~indet_fix(.))) %>%
  mutate(now_accepted_sist_subclass_superorder = indet_fix(now_accepted_sist_subclass_superorder),
         now_accepted_sist_order = indet_fix(now_accepted_sist_order),
         now_accepted_sist_suborder_superfamily = indet_fix(now_accepted_sist_suborder_superfamily),
         now_accepted_sist_family = indet_fix(now_accepted_sist_family),
         now_accepted_sist_subfamily_tribe = indet_fix(now_accepted_sist_subfamily_tribe),
         now_accepted_sist_genus = indet_fix(now_accepted_sist_genus),
         now_accepted_sist_epithet = indet_fix(now_accepted_sist_epithet),
         now_accepted_name = paste(now_accepted_sist_genus, now_accepted_sist_epithet, sep = "_")) %>%
  mutate(across(everything(),
                ~str_replace_all(., pattern = "[:space:]",
                                 replacement = ""))) %>%
  mutate(now_accepted_rank = case_when(
    now_accepted_sist_epithet != "indet." ~ "species",
    now_accepted_sist_epithet == "indet." & now_accepted_sist_genus != "indet." ~ "genus",
    now_accepted_sist_epithet == "indet." & now_accepted_sist_genus == "indet." ~ "Other"
  ),
  now_accepted_rank = case_when(
    now_accepted_rank == "Other" & now_accepted_sist_subfamily_tribe != "indet." ~ "subfamily_tribe",
    TRUE ~ now_accepted_rank
  ),
  now_accepted_rank = case_when(
    now_accepted_rank == "Other" & now_accepted_sist_family != "indet." ~ "family",
    TRUE ~ now_accepted_rank
  ),
  now_accepted_rank = case_when(
    now_accepted_rank == "Other" & now_accepted_sist_suborder_superfamily != "indet." ~ "suborder_superfamily",
    TRUE ~ now_accepted_rank
  ),
  now_accepted_rank = case_when(
    now_accepted_rank == "Other" & now_accepted_sist_order != "indet." ~ "order",
    TRUE ~ now_accepted_rank
  ),
  now_accepted_rank = case_when(
    now_accepted_rank == "Other" & now_accepted_sist_subclass_superorder != "indet." ~ "subclass_superorder",
    TRUE ~ now_accepted_rank
  ),
  now_accepted_name = str_replace(now_accepted_name, 
                                  pattern = "_indet\\.", 
                                  replacement = ""))

# Salvando o backbone
write.csv(now_taxa,
          file = here("data",
                      "processed",
                      "now_occs_names.csv"))



