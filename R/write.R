box::use(
  ut = utils
)

box::use(
  cop = R/clean_copernicus
)

write.csv(cop$cp_data, file = "data/cp_data.csv")
