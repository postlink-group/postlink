#' Print a `surv_mixtureBayes` Object
#' @description
#' Print the model call and posterior mean coefficients for component 1.
#'
#' @param x an object of class \code{surv_mixtureBayes}.
#' @param digits the number of significant digits to print.
#' @param ... for additional print arguments (unused).
#'
#' @returns invisibly returns the \code{surv_mixtureBayes} object that is provided.
#'
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
