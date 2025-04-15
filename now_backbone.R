#### Carregando pacotes com controle de versão ####

library(groundhog)
groundhog.day <- "2025-03-01"
groundhog.packages <- c("here",
                        "readr",
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

now_data <- read.csv(here("01_data", 
                          "01_raw-data",
                          "now_data.csv"), 
                     sep = ";",
                     strip.white = TRUE)

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
    now_accepted_sist_epithet != "indet." ~ "Species",
    now_accepted_sist_epithet == "indet." & now_accepted_sist_genus != "indet." ~ "Genus",
    now_accepted_sist_epithet == "indet." & now_accepted_sist_genus == "indet." ~ "Other"
  ),
  now_accepted_rank = case_when(
    now_accepted_rank == "Other" & now_accepted_sist_subfamily_tribe != "indet." ~ "Subfamily_tribe",
    TRUE ~ now_accepted_rank
  ),
  now_accepted_rank = case_when(
    now_accepted_rank == "Other" & now_accepted_sist_family != "indet." ~ "Family",
    TRUE ~ now_accepted_rank
  ),
  now_accepted_rank = case_when(
    now_accepted_rank == "Other" & now_accepted_sist_suborder_superfamily != "indet." ~ "Suborder_Superfamily",
    TRUE ~ now_accepted_rank
  ),
  now_accepted_rank = case_when(
    now_accepted_rank == "Other" & now_accepted_sist_order != "indet." ~ "Order",
    TRUE ~ now_accepted_rank
  ),
  now_accepted_rank = case_when(
    now_accepted_rank == "Other" & now_accepted_sist_subclass_superorder != "indet." ~ "Subclass_Superorder",
    TRUE ~ now_accepted_rank
  ),
  now_accepted_name = str_replace(now_accepted_name, 
                                  pattern = "_indet\\.", 
                                  replacement = ""))

# Salvando o backbone
write.csv(now_taxa,
          file = here("01_data",
                      "02_clean-data",
                      "now_backbone_carn.csv"))



