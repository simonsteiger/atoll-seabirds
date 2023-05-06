box::use(
  ut = utils,
  dp = dplyr,
  ts = tidyselect,
)

box::use(
  R/load,
  cp = R/clean_copernicus,
  ji = R/clean_jisao,
)

write.csv(cp$cp_data, file = "data/cp_data.csv")

write.csv(cp$out, file = "data/out.csv")

write.csv(ji$ji_data, file = "data/ji_data.csv")

write.csv(ji$out, file = "data/out_ji.csv")

ji_cp <- dp$left_join(ji$out, cp$out, by = "atoll") %>% 
  dp$select(-ts$contains(c("lat", "long")))

envs_jicp <- dp$left_join(load$envs, ji_cp, by = "atoll")

write.csv(envs_jicp, file = "data/envs_jicp.csv")
