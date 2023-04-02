box::use(
  sh = shiny,
  bsl = bslib,
  emo,
)

box::use(
  app/view/map,
  app/logic/theme,
  app / view / input,
  app / view / map,
  app / view / hist,
)

#' @export
ui <- function(id) {
  ns <- sh$NS(id)
  bsl$page(
    theme = theme$light,
    sh$div(
      class = "d-flex flex-column p-4 align-items-center",
      sh$tags$h1(paste0(emo$ji("fish"), " Atoll Explorer ", emo$ji("island")), class = "m-4"),
      bsl$layout_column_wrap(
        width = 1/2, height = 300,
        input$ui(ns("main")),
        hist$ui(ns("main"))
        )
    ),
    sh$div(
      class = "d-flex flex-column pb-4 align-items-center",
      bsl$layout_column_wrap(
        width = 1, height = 600,
        map$ui(ns("main"))
      )  
    )
    
    # sh$div(
    #   class = "d-flex flex-column p-4 align-items-center",
    #   
    # )
  )
}

#' @export
server <- function(id) {
  sh$moduleServer(id, function(input, output, session) {
    map$server("main")
    hist$server("main")
  })
}
