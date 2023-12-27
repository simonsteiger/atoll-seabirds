box::use(
  tbl = tibble,
  dp = dplyr,
  tdr = tidyr,
  st = stats,
  magrittr[`%>%`, set_colnames]
)

box::use(
  R / global,
)

#' @export
convert_to_tibble <- function(data, name, lat, long) {
  data %>%
    # maybe use set_colnames first, then as_tibble?
    tbl$as_tibble() %>%
    set_colnames(lat) %>%
    dp$bind_cols(longitude = long) %>%
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
      dp$across(c(global$names), ~ mean(.x, na.rm = TRUE), .names = "{.col}")
    )
}

#' @export
filter_jisao <- function(ji_data, lat, long) {
  ji_data %>%
    dp$mutate(
      dif_lat = abs(latitude - lat),
      dif_long = abs(longitude - long)
    ) %>%
    dp$filter(
      dif_lat == min(dif_lat) & dif_long == min(dif_long)
    ) %>%
    dp$select(-c(dif_lat, dif_long))
}
