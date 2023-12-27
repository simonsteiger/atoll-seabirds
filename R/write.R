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

pcs <- tbl$tibble(read.csv("data/jl_envscores.csv")) %>%
  dp$distinct(atoll, .keep_all = TRUE)

write.csv(cp$cp_data, file = "data/cp_data.csv")

write.csv(ji$ji_data, file = "data/ji_data.csv")

write.csv(ji$out, file = "data/out_ji.csv")

ji_cp <- dp$left_join(ji$out, cp$out, by = "atoll") %>% 
  dp$select(-ts$contains(c("lat", "long")))

envs_jicp <- dp$left_join(load$envs, ji_cp, by = "atoll")

write.csv(envs_jicp, file = "data/envs_jicp.csv")
