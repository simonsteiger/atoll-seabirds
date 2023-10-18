box::use(
  sh = shiny,
  swd = shinyWidgets,
  bsl = bslib,
  gir = ggiraph,
  htl = htmltools,
  dp = dplyr,
  magrittr[`%>%`],
  emo,
)

box::use(
  R / global,
  app / logic / disk,
  app / logic / draw,
)

#' @export
ui <- function(id) {
  ns <- sh$NS(id)

  bsl$card(
    bsl$card_header(class = "bg-info", "Map"),
    bsl$card_body_fill(
      gir$girafeOutput(ns("map"), width = "650px")
    )
  )
}

#' @export
server <- function(id) {
  sh$moduleServer(
    id,
    function(input, output, session) {
      # sh$observe({
      #   min_val <- round(min(disk$cp_data[[input$variable]], na.rm = TRUE), 2)
      #   max_val <- round(max(disk$cp_data[[input$variable]], na.rm = TRUE), 2)
# 
      #   sh$updateSliderInput(
      #     session = session,
      #     "slider",
      #     min = min_val,
      #     max = max_val,
      #     value = c(min_val, max_val)
      #   )
      # })

      #filtered_atolls <- sh$reactive(
      #  disk$out %>%
      #    dp$filter(.data[[input$variable]] %>% dp$between(min(input$slider), max(input$slider)))
      #)

      output$map <- gir$renderGirafe(draw$atoll_map(disk$out, input$variable))
    }
  )
}
