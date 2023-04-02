box::use(
  sh = shiny,
  swd = shinyWidgets,
)

box::use(
  R/global,
)

pick_var <- function(id) {
  swd$prettyRadioButtons(
    id,
    label = "Pick variable",
    choices = global$names,
    status = "primary",
    shape = "round",
    fill = TRUE,
    outline = FALSE,
    plain = FALSE,
    animation = "smooth",
    icon = sh$icon("check"),
    inline = TRUE
  ) 
}

slider_range <- function(id) {
  sh$sliderInput(
    id,
    label = "Show atolls with values between",
    value = c(0.2, 0.7),
    min = 0,
    max = 1,
    step = 0.2
  )
}