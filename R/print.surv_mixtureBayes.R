#' Print a brief summary of a surv_mixtureBayes model
#'
#' @param x An object of class \code{surv_mixtureBayes}.
#' @param digits Minimum number of significant digits to show.
#' @param ... Further arguments (unused).
#' @export
print.surv_mixtureBayes <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("Call:\n")
  print(x$call, quote = FALSE)
  cat("\n")
  cat("Coefficients (component 1 posterior means):\n")
  # Posterior mean of component 1 coefficients
  coef_mean <- apply(x$estimates$coefficients, 2, mean)
  print(format(signif(coef_mean, digits)), quote = FALSE, print.gap = 2)
  cat("\n")
  invisible(x)
}
