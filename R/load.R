box::use(
  nc = ncdf4,
  ut = utils,
  pr = purrr,
  magrittr[`%>%`],
)

box::use(
  R/global
)

#' @export
list_nc <- pr$set_names(global$names) %>%
  pr$map(~ nc$nc_open(paste0("data/copernicus_", .x, ".nc")))

#' @export
latitude <- as.vector(list_nc$nppv$dim[[3]]$vals)

#' @export
longitude <- as.vector(list_nc$nppv$dim[[4]]$vals)

list_var <- pr$set_names(global$names) %>% 
  pr$map(~ nc$ncvar_get(list_nc[[.x]], list_nc[[.x]]$var[[1]]))

#' @export
list_matrix <- pr$set_names(global$names) %>% 
  pr$map(~ apply(list_var[[.x]], c(1, 2), mean))

# list_matrix[[1]] <- colnames(latitude)
# list_matrix[[2]] <- colnames(latitude)
# list_matrix[[3]] <- colnames(latitude)

#' @export
envs <- ut$read.csv("data/seabird_atolls_envs_10Mar.csv")

#' @export
pop <- ut$read.csv("data/atoll_seabird_populations_10Mar.csv")
