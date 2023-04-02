box::use(
  gg = ggplot2,
  dp = dplyr,
  magrittr[`%>%`],
  gir = ggiraph,
  sh = shiny,
  pal = palettes,
)

box::use(
  app/logic/disk,
  col = app/logic/colors,
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

#' @export
atoll_map <- function(fill) {
  
  pal_name <- switch(fill,
    "nppv" = "mako",
    "chl" = "vangogh3",
    "phyc" = "okeeffe2"
  )
  
  p <- sh$reactive(
    gg$ggplot() +
      gg$geom_raster(data = disk$cp_data, gg$aes(longitude, latitude, fill = .data[[fill]])) +
      gir$geom_point_interactive(
        data = disk$out,
        gg$aes(
          x = long, 
          y = lat,
          color = .data[[fill]],
          tooltip = paste0(
            "<b>Atoll</b> ", disk$out$atoll, "\n<b>", toupper(fill),"</b> ", round(disk$out[[fill]], 2)
          ),
          data_id = atoll
        ),
        size = 0.5,
        alpha = 0.8,
        hover_nearest = TRUE
      ) +
      pal$scale_fill_palette_c(
        col$get_palette(1000, pal_name = pal_name),
        na.value = "grey80"
      ) +
      gg$scale_color_viridis_c(option = "inferno", begin = 0.4) +
      gg$coord_quickmap() +
      gg$theme_minimal() +
      gg$theme(legend.position = "bottom") +
      gg$guides(
        colour = gg$guide_colorbar(
          title = paste0("Atoll vicinity ", fill),
          title.position = "bottom",
          title.hjust = 0.5,
          label.position = "top"
          ),
        fill = gg$guide_colorbar(
          title = paste0("Oceanic ", fill),
          title.position = "bottom", 
          title.hjust = 0.5, 
          label.position = "top"
          )
      )
  )
  
  gir$girafe(ggobj = p())
}
