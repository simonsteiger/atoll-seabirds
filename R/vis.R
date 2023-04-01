box::use(
  sf,
  gg = ggplot2,
  dp = dplyr,
  magrittr[`%>%`],
  lf = leaflet,
  gir = ggiraph,
)

box::use(
  cop = R/clean_copernicus,
  R/load,
)

# sf_data <- sf$st_as_sf(cop$cp_data, coords = c("longitude", "latitude")) %>% 
#   dp$filter(!is.na(nppv))

# gg$ggplot() +
#   gg$geom_sf(data = sf_data) +
#   gg$geom_point(data = cop$envs, gg$aes(x = long, y = lat), color = "red", size = 0.1)

p <- 
  gg$ggplot() +
  gg$geom_raster(data = cop$cp_data, gg$aes(longitude, latitude, fill = log(nppv))) +
  gir$geom_point_interactive(
    data = cop$out,
    gg$aes(
      x = long, 
      y = lat,
      color = mean_nppv,
      tooltip = paste0("name: ", cop$out$atoll, "\nmean nppv: ", round(cop$out$mean_nppv, 2))
      ),
    size = 1,
    alpha = 0.8) +
  gg$scale_fill_viridis_c(option = "mako") +
  gg$scale_color_viridis_c(option = "inferno", begin = 0.4) +
  gg$coord_quickmap() +
  gg$theme_minimal()

gir$girafe(ggobj = p)

lf$leaflet() %>% 
  lf$addProviderTiles(lf$providers$Stamen.TonerLite,
                   options = lf$providerTileOptions(noWrap = TRUE)
  ) %>%
  lf$addMarkers(data = cop$envs %>% dp$select(lat, long))
