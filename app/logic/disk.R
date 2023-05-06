box::use(
  ut = utils,  
  dp = dplyr,
  magrittr[`%>%`],
)

box::use(
  
)


#' @export
out_ji <- ut$read.csv("data/out_ji.csv")

#' @export
cp_data <- ut$read.csv("data/cp_data.csv")

#' @export
out <- ut$read.csv("data/out.csv")

#' @export
ji_data <- ut$read.csv("data/ji_data.csv")

ggplot() + 
  geom_raster(aes(longitude, latitude, fill = precip_anom), ji_data) + 
  geom_point(aes(longitude, latitude, color = precip_anom), out_ji) +
  scale_fill_viridis_c() +
  scale_color_viridis_c(option = "magma")
