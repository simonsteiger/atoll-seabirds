box::use(
  ut = utils,
  dp = dplyr,
  ts = tidyselect,
  tbl = tibble,
  magrittr[`%>%`],
  tdr = tidyr,
  here,
)

box::use(
  R / wrangle / load,
  cp = R / wrangle / clean_copernicus,
  ji = R / wrangle / clean_jisao,
)

write.csv(cp$cp_data, file = here$here("data/remotesensing/cp_data.csv"))

write.csv(ji$ji_data, file = here$here("data/remotesensing/ji_data.csv"))

write.csv(ji$out, file = here$here("data/remotesensing/out_ji.csv"))

ji_cp <- dp$left_join(ji$out, cp$out, by = "atoll") %>%
  dp$select(-ts$contains(c("lat", "long")))

envs_jicp <- dp$left_join(load$envs, ji_cp, by = "atoll")

write.csv(envs_jicp, file = here$here("data/remotesensing/envs_jicp.csv"))
