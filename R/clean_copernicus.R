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
  fn = R/functions,
)

# fix negative longitudes
#' @export
envs <- load$envs %>% 
  dp$mutate(
    long = ifelse(long < 0, long + 360, long)
  )

#' @export
cp_data <- pr$set_names(load$names) %>% 
  pr$map(~ fn$convert_to_tibble(load$list_matrix[[.x]], .x)) %>% 
  pr$reduce(dp$left_join, by = c("latitude", "longitude"))
  

#' @export
out <- 
  pr$map(
  seq_len(nrow(envs)), # nrow envs = number of atolls
  ~ fn$filter_copernicus(cp_data, envs[.x, ]$lat, envs[.x, ]$long)
  ) %>%
  pr$list_rbind() %>%
  dp$mutate(
    atoll = envs$atoll,
    lat = envs$lat,
    long = envs$long
  )

#' @export
cop_envs <- dp$left_join(load$envs, out, by = "atoll", suffix = c("", ".dupl")) %>% 
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
