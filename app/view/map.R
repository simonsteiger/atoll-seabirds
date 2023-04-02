box::use(
  sh = shiny,
  swd = shinyWidgets,
  bsl = bslib,
  gir = ggiraph,
  bsi = bsicons,
)

box::use(
  R / global,
  app / logic / disk,
  app / logic / draw,
)

#' @export
ui <- function(id) {
  ns <- sh$NS(id)
  sh$tagList(
    bsl$card(
      bsl$card_header("Controls"),
      bsl$card_body_fill(
        swd$prettyRadioButtons(
          ns("variable"),
          label = "Pick variable",
          choices = global$names,
          status = "primary",
          shape = "round",
          fill = TRUE,
          outline = TRUE,
          plain = FALSE,
          animation = "smooth",
          icon = sh$icon("check"),
          inline = TRUE
        )
      )
    ),
    bsl$card(
      full_screen = TRUE,
      bsl$card_header("Map"),
      bsl$card_body(
        gir$girafeOutput(ns("map"), height = "800px")
      )
    )
  )
}

#' @export
server <- function(id) {
  sh$moduleServer(
    id,
    function(input, output, session) {
      output$map <- gir$renderGirafe(draw$atoll_map(input$variable))
    }
  )
}
