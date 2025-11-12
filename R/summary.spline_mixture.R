#' Summarize an `spline_mixture` Object
#' @description Summarize results from a `spline_mixture()` object
#'
#' @param x the result of a call to `spline_mixture()`
#' @param digits the number of significant digits to print
#' @param ... for additional summary arguments
#'
#' @returns invisibly returns the `spline_mixture()` object
#' that is provided as an argument
#'
#' @export
summary.spline_mixture <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
 cat("Call:\n")
 print(x$call$call, quote = F, digits = digits)
 cat("\n")

 cat("Estimated Error Variance (Posterior Mean):\n")
 print(x$par_sigmasq$rate / (x$par_sigmasq$shape - 1), quote = F, digits = digits)
 cat("\n")

 cat("R-squared:\n")
 print(x$R2, quote = F, digits = digits)
 cat("\n")

 cat("Estimated Mismatch Rate (Posterior Mean):\n")
 print(x$par_alpha$a / (x$par_alpha$a + x$par_alpha$b) , quote = F, digits = digits)
 cat("\n")

 invisible(x)
}
