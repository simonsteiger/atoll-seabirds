
library("brms")
library("broom")
library("dplyr")
library("tidymodels")
library("magrittr")
library("stringr")

pop <- read.csv("data/atoll_seabird_populations_25Feb.csv")
envs <- read.csv("data/seabird_atolls_envs_25Feb.csv")

laridae <- c("Anous", "Gygis", "Onychoprion", "Sternula", "Hydroprogne", "Sterna", "Thalasseus")
phaethontidae <- "Phaethon"
fregatidae <- "Fregata"
sulidae <- "Sula"

pop_recode <- pop %>% 
  dplyr::select(-starts_with("X")) %>% 
  dplyr::mutate(
    across(c(where(is.character), -atoll), as.numeric)
  ) %>% 
  tidyr::pivot_longer(!atoll, names_to = "species", values_to = "presence") %>% 
  dplyr::mutate(
    presence = ifelse(!is.na(presence), TRUE, FALSE),
    family = dplyr::case_when(
      str_detect(species, paste0(laridae, collapse = "|")) ~ "laridae",
      str_detect(species, paste0(phaethontidae, collapse = "|")) ~ "phaethontidae",
      str_detect(species, paste0(fregatidae, collapse = "|")) ~ "fregatidae",
      str_detect(species, paste0(sulidae, collapse = "|")) ~ "sulidae",
      .default = NULL
    )
  ) %>% 
  drop_na(family)

joined <- left_join(pop_recode, envs, by = "atoll")



