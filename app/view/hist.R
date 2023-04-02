box::use(
  sh = shiny,
  swd = shinyWidgets,
  bsl = bslib,
  gir = ggiraph,
  htl = htmltools,
  dp = dplyr,
  magrittr[`%>%`],
  emo,
  e4r = echarts4r,
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
    bsl$card_header(class = "bg-info", "Histogram"),
    bsl$card_body_fill(
      e4r$echarts4rOutput(ns("hist"), width = "100%")
    )
  )
}

#' @export
server <- function(id) {
  sh$moduleServer(
    id,
    function(input, output, session) {
      sh$observe({
        min_val <- round(min(disk$cp_data[[input$variable]], na.rm = TRUE), 2)
        max_val <- round(max(disk$cp_data[[input$variable]], na.rm = TRUE), 2)

        sh$updateSliderInput(
          session = session,
          "slider",
          min = min_val,
          max = max_val,
          value = c(min_val, max_val)
        )
      })

      filtered_atolls <- sh$reactive(
        disk$out %>%
          dp$filter(.data[[input$variable]] %>% dp$between(min(input$slider), max(input$slider)))
      )

      output$hist <- e4r$renderEcharts4r(draw$atoll_hist(filtered_atolls(), input$variable))
    }
  )
}
