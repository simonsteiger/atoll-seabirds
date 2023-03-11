
library("brms")
library("broom")
library("dplyr")
library("tidyr")
library("magrittr")
library("stringr")

pop <- read.csv("data/atoll_seabird_populations_10Mar.csv")
envs <- read.csv("data/seabird_atolls_envs_10Mar.csv")

pop_recode <- pop %>% 
  dplyr::select(-starts_with("X")) %>% 
  dplyr::mutate(
     across(c(where(is.character), -atoll), as.numeric)
   ) %>% 
  tidyr::pivot_longer(!atoll, names_to = "species", values_to = "presence") %>% 
  dplyr::mutate(
    presence = ifelse(!is.na(presence), TRUE, FALSE),
  )

#' @export
joined <- left_join(pop_recode, envs, by = "atoll")
