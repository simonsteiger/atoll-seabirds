box::use(
  tbl = tibble,
  dp = dplyr,
  tdr = tidyr,
  st = stats,
  magrittr[`%>%`, set_colnames]
)

box::use(
  R/load,
)

#' @export
convert_to_tibble <- function(data, name) {
  data %>% 
    # maybe use set_colnames first, then as_tibble?
    tbl$as_tibble() %>%
    set_colnames(load$latitude) %>%
    dp$bind_cols(longitude = load$longitude) %>% 
    tdr$pivot_longer(-longitude, names_to = "latitude", values_to = name) %>% 
    dp$mutate(
      latitude = as.numeric(latitude)
    )
}

#' @export
filter_copernicus <- function(cp_data, lat, long) {
  cp_data %>% 
    dp$filter(
      latitude >= lat - 1,
      latitude <= lat + 1,
      longitude >= long - 1,
      longitude <= long + 1
    ) %>% 
    dp$summarise(
      across(c(load$names), ~ mean(.x, na.rm = TRUE), .names = "mean_{.col}")
    ) 
}