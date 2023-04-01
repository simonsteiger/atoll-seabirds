box::use(
  lub = lubridate,
  dp = dplyr,
  tdr = tidyr,
  tbl = tibble,
  st = stats,
  pr = purrr,
  magrittr[`%>%`],
)

box::use(
  R/load
)

# fix negative longitudes
#' @export
envs <- load$envs %>% 
  dp$mutate(
    long = ifelse(long < 0, long + 360, long)
  )

#' @export
cp_data <- tbl$as_tibble(load$nc_mean) %>%
  dp$bind_cols(longitude = load$longitude) %>% 
  tdr$pivot_longer(-longitude, names_to = "latitude", values_to = "nppv") %>% 
  dp$mutate(
    latitude = as.numeric(latitude)
    )

filter_copernicus <- function(cp_data, lat, long) {
  cp_data %>% 
    dp$filter(
      latitude >= lat - 1,
      latitude <= lat + 1,
      longitude >= long - 1,
      longitude <= long + 1
    ) %>% 
    dp$summarise(
      mean_nppv = mean(nppv, na.rm = TRUE),
      sd_nppv = st$sd(nppv, na.rm = TRUE)
      )
}

#' @export
out <- pr$map(
  seq_len(nrow(envs)),
  ~ filter_copernicus(cp_data, envs[.x, ]$lat, envs[.x, ]$long)
  ) %>%
  pr$list_rbind() %>% 
  dp$mutate(
    atoll = envs$atoll,
    lat = envs$lat,
    long = envs$long,
    .before = mean_nppv
  )

# What was that for?
# res <- data %>% 
#   dplyr::filter(
#     latitude >= target_lat - 1,
#     latitude <= target_lat + 1,
#     longitude >= target_long - 1,
#     longitude <= target_long + 1
#   ) %>% 
#   summarise(mean = mean(nppv, na.rm = TRUE))
