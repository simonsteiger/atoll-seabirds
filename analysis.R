
library("brms")
library("broom")
library("dplyr")
library("tidyr")
library("magrittr")
library("stringr")

# source("functions.R")

fam_filter <- function(.data, v) {
  .data %>% 
    filter(family == v)
}

wrap_brm <- function(.data, formula) {
  brm(formula = formula,
      family = bernoulli(link = "logit"),
      data = .data,
      warmup = 1000,
      iter = 2000,
      chains = 4,
      thin = 1
  )
}

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
  tidyr::drop_na(family)

joined <- left_join(pop_recode, envs, by = "atoll")

birdfamilies <- c("phaethontidae", "laridae", "fregatidae", "sulidae")

list_df <- set_names(birdfamilies) %>% 
  map(~fam_filter(joined, v = .x))
  
simple_mod <- presence ~ number_islets*species + land_area_sqkm*species + human_population*species

set.seed(1252)
list_fit <- set_names(birdfamilies) %>% 
  map(~wrap_brm(.data = list_df[[.x]], formula = simple_mod))

