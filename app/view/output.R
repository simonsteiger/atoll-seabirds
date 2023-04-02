box::use(
  sh = shiny,
  swd = shinyWidgets,
  bsl = bslib,
  gir = ggiraph,
  htl = htmltools,
  e4r = echarts4r,
  emo,
)

box::use(
  R / global,
  app / logic / disk,
  app / logic / draw,
  inp = app / logic / input_ui,
)

#' @export
ui <- function(id) {
  ns <- sh$NS(id)
  sh$tagList(
    bsl$card(
      class = "m-8 shadow",
      bsl$card_header(class = "bg-info", "Map"),
      bsl$card_body_fill(
        gir$girafeOutput(ns("map"), width = "100%")
      )
    ),
    bsl$card(
      class = "m-8 shadow",
      bsl$card_header(class = "bg-info", "Histogram"),
      bsl$card_body_fill(
        e4r$echarts4rOutput(ns("hist"), width = "100%")
      )
    )
  )
}
