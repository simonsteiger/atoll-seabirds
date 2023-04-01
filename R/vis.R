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

css_default_hover <- gir$girafe_css_bicolor(primary = "#1bf702", secondary = "#1bf702")
css_box <- "font-family:Arial;padding:7px;background-color:#ffffff;color:#000000;border-radius:3px;"

gir$set_girafe_defaults(
  opts_hover = gir$opts_hover(css = css_default_hover),
  opts_zoom = gir$opts_zoom(min = 1, max = 4),
  opts_tooltip = gir$opts_tooltip(css = css_box),
  opts_sizing = gir$opts_sizing(rescale = TRUE),
  opts_toolbar = gir$opts_toolbar(saveaspng = FALSE, position = "bottom", delay_mouseout = 5000)
)

p <- 
  gg$ggplot() +
  gg$geom_raster(data = cop$cp_data, gg$aes(longitude, latitude, fill = log(nppv))) +
  gir$geom_point_interactive(
    data = cop$out,
    gg$aes(
      x = long, 
      y = lat,
      color = mean_nppv,
      tooltip = paste0("<b>Atoll</b> ", cop$out$atoll, "\n<b>NPPV</b> ", round(cop$out$mean_nppv, 2)),
      data_id = atoll
      ),
    size = 0.5,
    alpha = 0.8,
    hover_nearest = TRUE
    ) +
  gg$scale_fill_viridis_c(option = "mako") +
  gg$scale_color_viridis_c(option = "inferno", begin = 0.4) +
  gg$coord_quickmap() +
  gg$theme_minimal() +
  gg$theme(legend.position = "bottom") +
  gg$guides(
    colour = gg$guide_colorbar(title.position = "bottom", title.hjust = 0.5, label.position = "top"),
    fill = gg$guide_colorbar(title.position = "bottom", title.hjust = 0.5, label.position = "top")
    )

gir$girafe(ggobj = p)

lf$leaflet() %>% 
  lf$addProviderTiles(lf$providers$Stamen.TonerLite,
                   options = lf$providerTileOptions(noWrap = TRUE)
  ) %>%
  lf$addMarkers(data = cop$envs %>% dp$select(lat, long))


# sf_data <- sf$st_as_sf(cop$cp_data, coords = c("longitude", "latitude")) %>% 
#   dp$filter(!is.na(nppv))

# gg$ggplot() +
#   gg$geom_sf(data = sf_data) +
#   gg$geom_point(data = cop$envs, gg$aes(x = long, y = lat), color = "red", size = 0.1)
