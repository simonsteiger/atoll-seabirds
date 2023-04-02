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
  inp = app / logic / input_ui,
)

#' @export
ui <- function(id) {
  ns <- sh$NS(id)
  sh$tagList(
    sh$div(
      class = "d-flex flex-column p-4 align-items-center",
      sh$tags$h1(paste0(emo$ji("fish"), " Atoll Explorer ", emo$ji("island")), class = "m-4"),
      bsl$card(
        class = "m-8 shadow",
        bsl$card_header(class = "bg-info", "Controls"),
        bsl$card_body_fill(
          class = "bg-fade",
          inp$pick_var(ns("variable")),
          inp$slider_range(ns("slider"))
        )
      ),
      bsl$card(
        class = "m-8 shadow",
        full_screen = TRUE,
        bsl$card_header(class = "bg-info", "Map"),
        bsl$card_body_fill(
          gir$girafeOutput(ns("map"), width = "100%")
        )
    )
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
      
      output$map <- gir$renderGirafe(draw$atoll_map(filtered_atolls(), input$variable))
    }
  )
}
