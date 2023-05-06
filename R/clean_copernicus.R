box::use(
  lub = lubridate,
  dp = dplyr,
  tdr = tidyr,
  tbl = tibble,
  st = stats,
  pr = purrr,
  ts = tidyselect,
  magrittr[`%>%`],
)

box::use(
  R/load,
  R/global,
  fn = R/functions,
)

params <- tbl$tibble(
  data = load$list_matrix,
  name = global$names,
  lat = c(rep(list(load$latitude), 3), rep(list(load$latitude_fine), 2)),
  long = c(rep(list(load$longitude), 3), rep(list(load$longitude_fine), 2))
  )

#' @export
cp_data <- 
  pr$pmap(
    params,
    fn$convert_to_tibble, 
    .progress = "Converting matrices to tibbles..."
    ) %>% 
  pr$set_names(global$names) %>% 
  pr$reduce(dp$left_join, by = c("latitude", "longitude")) %>%
  dp$mutate(
    dp$across(global$names, log1p)
  )

#' @export
out <- 
  pr$map(
  seq_len(nrow(load$envs_trans_coord)), # nrow load$envs_trans_coord = number of atolls
  ~ fn$filter_copernicus(cp_data, load$envs_trans_coord[.x, ]$lat, load$envs_trans_coord[.x, ]$long),
  .progress = "Summarising variables per atoll..."
  ) %>%
  pr$list_rbind() %>%
  dp$mutate(
    atoll = load$envs_trans_coord$atoll,
    lat = load$envs_trans_coord$lat,
    long = load$envs_trans_coord$long
  )

#' @export
envs <- dp$left_join(load$envs_trans_coord, out, by = "atoll", suffix = c("", ".dupl")) %>% 
  dp$select(-ts$ends_with(".dupl"))

# What was that for?
# res <- data %>% 
#   dplyr::filter(
#     latitude >= target_lat - 1,
#     latitude <= target_lat + 1,
#     longitude >= target_long - 1,
#     longitude <= target_long + 1
#   ) %>% 
#   summarise(mean = mean(nppv, na.rm = TRUE))
