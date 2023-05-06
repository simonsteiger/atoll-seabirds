box::use(
  nc = ncdf4,
  ut = utils,
  pr = purrr,
  ab = abind,
  dp = dplyr,
  ts = tidyselect,
  magrittr[`%>%`],
)

box::use(
  R/global,
  R/load,
  fn = R/functions,
)

data_nc <- nc$nc_open("data/preci_anom_cor.nc")

latitude <- as.vector(data_nc$dim[[1]]$vals)

longitude <- as.vector(data_nc$dim[[2]]$vals)

whole_array <- nc$ncvar_get(data_nc, data_nc$var[[1]])

precip_anom <- whole_array[, , 2] # cm/month is second slice

#' @export
ji_data <- fn$convert_to_tibble(precip_anom, "precip_anom", lat = latitude, long = longitude)

#' @export
out <- 
  pr$map(
    seq_len(nrow(load$envs)), # nrow load$envs = number of atolls
    ~ fn$filter_jisao(ji_data, load$envs_trans_coord[.x, ]$lat, load$envs_trans_coord[.x, ]$long)
  ) %>%
  pr$list_rbind() %>%
  dp$mutate(
    atoll = load$envs$atoll,
    lat = load$envs$lat,
    long = load$envs$long
  )

#' @export
envs <- dp$left_join(load$envs, out, by = "atoll", suffix = c("", ".dupl")) %>% 
  dp$select(-ts$ends_with(".dupl"))

