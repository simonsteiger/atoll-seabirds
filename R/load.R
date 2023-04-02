box::use(
  nc = ncdf4,
  ut = utils,
)


#' @export
nc_file <- nc$nc_open("data/copernicus_nppv.nc")

#' @export
latitude <- as.vector(nc_file$dim[[3]]$vals)

#' @export
longitude <- as.vector(nc_file$dim[[4]]$vals)

x <- nc$ncvar_get(nc_file, nc_file$var[[1]])

matrix <- apply(x, c(1, 2), mean)

colnames(matrix) <- latitude

#' @export
nc_mean <- matrix

#' @export
envs <- ut$read.csv("data/seabird_atolls_envs_10Mar.csv")

#' @export
pop <- ut$read.csv("data/atoll_seabird_populations_10Mar.csv")