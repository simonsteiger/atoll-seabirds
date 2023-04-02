box::use(
  pal = palettes,
  vir = viridisLite,
)

pal_list <- pal$pal_palette(
  mako = vir$mako(100),
  # inspired by met_palette VanGogh3, added darker start and changed last color
  vangogh3 = c("#0a0f07", "#192813", "#1E3D14", "#1F5B25", "#3C7C3D", "#669D62", "#9CC184", "#C2D184", "#E6F2D0"),
  # inspired by met_palette OKeeffe2, added several dark tones
  okeeffe2 = c("#0f0301", "#381911", "#61291a", "#92351E", "#B9563F", "#D37750", "#E69C6B", "#ECB27D", "#F2C88F", "#FBE3C2")
)

#' @export
get_palette <- function(n, pal_name) {
  pillars <- pal_list[[pal_name]]
  pal$pal_ramp(pillars, n = n, interpolate = "spline")
}
