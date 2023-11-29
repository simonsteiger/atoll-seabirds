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
  R / global
)

#' @export
list_nc <- pr$set_names(global$names_load) %>%
  pr$map(~ nc$nc_open(paste0("data/copernicus_", .x, ".nc")))

#' @export
latitude <- as.vector(list_nc$nppv$dim[[3]]$vals)

#' @export
longitude <- as.vector(list_nc$nppv$dim[[4]]$vals)

#' @export
latitude_fine <- as.vector(list_nc$temp_1$dim[[3]]$vals)

#' @export
longitude_fine <- as.vector(list_nc$temp_1$dim[[4]]$vals)

list_var <- pr$set_names(global$names_load) %>%
  pr$map(~ nc$ncvar_get(list_nc[[.x]], list_nc[[.x]]$var[[1]]), .progress = "Getting variables...")

# absolute values
list_abs <- pr$set_names(global$names_load) %>%
  pr$map(~ abs(list_var[[.x]]), .progress = "Calculating absolute values...")

list_matrix <- pr$set_names(global$names_load) %>%
  pr$map(~ apply(list_abs[[.x]], c(1, 2), mean), .progress = "Calculating means...")

list_matrix$temp <- ab$abind(list_matrix$temp_1, list_matrix$temp_2, along = 3) %>%
  apply(., c(1, 2), mean)

list_matrix$velo <- ab$abind(
  list_matrix$velo_1,
  list_matrix$velo_2,
  list_matrix$velo_3,
  list_matrix$velo_4,
  along = 3
) %>%
  apply(., c(1, 2), mean, na.rm = TRUE)

#' @export
list_matrix <- list_matrix[-c(4:9)]

#' @export
envs <- ut$read.csv("data/seabird_atolls_envs.csv")

#' @export
envs_trans_coord <- envs %>%
  dp$mutate(
    long = ifelse(long < 0, long + 360, long)
  )

#' @export
pop <- ut$read.csv("data/atoll_seabird_populations.csv")

# Assert that all counts are integers
if (!all(dp$summarise(pop, dp$across(ts$where(is.numeric), is.integer)))) {
  stop("Population counts must be integers!")
}
