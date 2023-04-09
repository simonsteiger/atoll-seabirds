box::use(
  ut = utils,  
)

box::use(
  
)

#' @export
cp_data <- ut$read.csv("data/cp_data.csv")

#' @export
out <- ut$read.csv("data/out_cp.csv")

#' @export
ji_data <- ut$read.csv("data/ji_data.csv")

#' @export
out_ji <- ut$read.csv("data/out_ji.csv")