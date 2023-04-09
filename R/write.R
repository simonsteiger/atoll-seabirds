box::use(
  ut = utils
)

box::use(
  cp = R/clean_copernicus,
  ji = R/clean_jisao,
)

write.csv(cp$cp_data, file = "data/cp_data.csv")

write.csv(cp$out, file = "data/out.csv")

write.csv(ji$ji_data, file = "data/ji_data.csv")

write.csv(ji$out, file = "data/out_ji.csv")
