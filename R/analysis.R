box::use(
  dp = dplyr,
  tdr = tidyr,
  ts = tidyselect,
)

pop <- read.csv("data/atoll_seabird_populations_11Mar.csv")
envs_jicp <- read.csv("data/envs_jicp.csv")

pop_recode <- pop %>% 
  dp$select(-starts_with("X")) %>% 
  dp$mutate(
     dp$across(c(ts$where(is.character), -atoll), as.numeric)
   ) %>% 
  tdr$pivot_longer(!atoll, names_to = "species", values_to = "presence") %>% 
  dp$mutate(
    presence = ifelse(!is.na(presence), TRUE, FALSE),
  )

#' @export
joined <- dp$left_join(pop_recode, envs_jicp, by = "atoll") %>% 
  tdr$pivot_wider(names_from = "species", values_from = "presence") %>% 
  dp$select(-X)
