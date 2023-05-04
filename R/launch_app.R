#' Launch Shiny App
#'
#' @param name The name of the app to run
#' @param ... arguments to pass to shiny::runApp
#'
#' @export
#'
launch_app <- function(name = "app.R") {
  source(system.file(paste0("apps/", name), package = "ChromoCorrect"), chdir = FALSE)$value
}
