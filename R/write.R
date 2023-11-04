box::use(
  ut = utils,
  dp = dplyr,
  ts = tidyselect,
  tbl = tibble,
  magrittr[`%>%`],
)

box::use(
  R/load,
  cp = R/clean_copernicus,
  ji = R/clean_jisao,
)

pcs <- tbl$tibble(read.csv("data/jl_envscores.csv")) %>%
  dp$distinct(atoll, .keep_all = TRUE)

predictions <- tbl$tibble(read.csv("data/newpreds.csv"))

observed <- tbl$tibble(read.csv("data/atoll_seabird_populations_29Jul.csv")) %>%
  dp$mutate(dp$across(ts$where(is.numeric), \(x) ifelse(is.na(x), 0, 1)))

full_presence <- dp$bind_rows(observed, predictions) %>%
  dp$left_join(pcs, by = "atoll")

# TODO pivot longer -> merge by atoll, species -> pivot wider

out <- dp$left_join(cp$out, full_presence)

write.csv(out, file = "data/full_presence.csv")

write.csv(cp$cp_data, file = "data/cp_data.csv")

write.csv(out, file = "data/out.csv")

write.csv(ji$ji_data, file = "data/ji_data.csv")

write.csv(ji$out, file = "data/out_ji.csv")

ji_cp <- dp$left_join(ji$out, cp$out, by = "atoll") %>% 
  dp$select(-ts$contains(c("lat", "long")))

envs_jicp <- dp$left_join(load$envs, ji_cp, by = "atoll")

write.csv(envs_jicp, file = "data/envs_jicp.csv")
