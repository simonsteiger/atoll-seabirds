box::use(
  ut = utils,
  dp = dplyr,
  ts = tidyselect,
  tbl = tibble,
  magrittr[`%>%`],
  tdr = tidyr,
)

box::use(
  R/load,
  cp = R/clean_copernicus,
  ji = R/clean_jisao,
)

predictions <- tbl$tibble(read.csv("data/predictpresence.csv"))

observed <- tbl$tibble(read.csv("data/atoll_seabird_populations_29Jul.csv")) %>%
  dp$mutate(dp$across(ts$where(is.numeric), \(x) ifelse(is.na(x), 0, 1)))

obs_long <- tdr$pivot_longer(observed, !atoll, names_to = "species", values_to = "presence")
pred_long <- tdr$pivot_longer(predictions, !atoll, names_to = "species", values_to = "presence")
full_long <- dp$bind_rows(obs_long, pred_long)
full_wide <- tdr$pivot_wider(full_long, names_from = "species", values_from = "presence")

write.csv(full_wide, file = "data/full_presence.csv")

write.csv(cp$cp_data, file = "data/cp_data.csv")

write.csv(out, file = "data/out.csv")

write.csv(ji$ji_data, file = "data/ji_data.csv")

write.csv(ji$out, file = "data/out_ji.csv")

ji_cp <- dp$left_join(ji$out, cp$out, by = "atoll") %>% 
  dp$select(-ts$contains(c("lat", "long")))

envs_jicp <- dp$left_join(load$envs, ji_cp, by = "atoll")

write.csv(envs_jicp, file = "data/envs_jicp.csv")
