#' Print an `spline_mixture` Object
#' @description Print call and input for a `spline_mixture()` object
#'
#' @param x the result of a call to `spline_mixture()`
#' @param digits the number of significant digits to print
#' @param ... for additional print arguments
#'
#' @returns invisibly returns the `spling_mixture()` object
#' that is provided as an argument
#'
#' @export
print.spline_mixture <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
 cat("\n")
 cat("Family:", x$family$family,"\n")
 cat("Link function:", x$family$link, "\n")
 cat("\n")

 cat("Call:\n")
 print(x$call$call, quote = F, digits = digits)
 cat("\n")

 invisible(x)
}
