#' Launch Shiny App
#'
#' @param name The name of the app to run
#' @param ... arguments to pass to shiny::runApp
#'
#' @export
#'
launch_app <- function(name = "app.R", wd = getwd()) {
  appDir <- system.file(paste0("apps/", name), package = "ChromoCorrect")
  if (appDir == "") stop("The shiny app ", name, " does not exist")
  shiny::runApp(appDir, wd = wd)
}
