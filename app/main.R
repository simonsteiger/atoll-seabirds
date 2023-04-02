box::use(
  sh = shiny,
  bsl = bslib,
)

box::use(
  app/view/map,
  app/logic/theme,
)

#' @export
ui <- function(id) {
  ns <- sh$NS(id)
  bsl$page(
    theme = theme$light,
    map$ui(ns("map"))
  )
}

#' @export
server <- function(id) {
  sh$moduleServer(id, function(input, output, session) {
    map$server("map")
  })
}
