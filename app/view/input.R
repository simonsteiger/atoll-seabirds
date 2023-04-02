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
  
      bsl$card(
        class = "m-8 shadow",
        bsl$card_header(class = "bg-info", "Controls"),
        bsl$card_body_fill(
          inp$pick_var(ns("variable")),
          inp$slider_range(ns("slider"))
        
      
    )
  )
}