#' Get path to parsemsf example
#'
#' From Hadley Wickham's readr package. parsemsf comes bundled with a number of sample
#' files in its \code{inst/extdata} directory. This function make them easy to access
#'
#' @param path Name of file
#' @export
#' @keywords internal
#' @examples
#' parsemsf_example("test_db.msf")
parsemsf_example <- function(path) {
  system.file("extdata", path, package = "parsemsf", mustWork = TRUE)
}
