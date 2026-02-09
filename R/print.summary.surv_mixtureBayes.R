#' Print a `summary.surv_mixtureBayes` Object
#' @description
#' Print a formatted summary for a \code{surv_mixtureBayes} fit, including
#' posterior mean and uncertainty summaries for both mixture components.
#'
#' @param x a \code{summary.surv_mixtureBayes} object.
#' @param digits the number of significant digits to print.
#' @param ... for additional print arguments (unused).
#'
#' @returns invisibly returns the \code{summary.surv_mixtureBayes} object that is provided.
#'
#' @export

print.summary.surv_mixtureBayes <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("Call:\n")
  print(x$call, quote = FALSE)
  cat("\nFamily:", x$family, "\n\n")

  # Component 1 (usually correct matches in linkage context)
  cat("(For Correct Matches):\n")
  cat("Outcome Model Coefficients:\n")
  stats::printCoefmat(x$coefficients, digits = digits, signif.stars = FALSE)
  cat("\n")
  if(x$family == "gamma") {
    cat("Shape parameter (component 1):\n")
    cat(format(signif(x$shape, digits)), "\n\n")
  } else if(x$family == "weibull") {
    cat("Shape parameter (component 1):\n")
    cat(format(signif(x$shape1, digits)), "\n")
    cat("Scale parameter (component 1):\n")
    cat(format(signif(x$scale1, digits)), "\n\n")
  }

  # Component 2 (mismatches)
  cat("(For Mismatches):\n")
  cat("Outcome Model Coefficients:\n")
  stats::printCoefmat(x$m.coefficients, digits = digits, signif.stars = FALSE)
  cat("\n")
  if(x$family == "gamma") {
    cat("Shape parameter (component 2):\n")
    cat(format(signif(x$m.shape, digits)), "\n")
  } else if(x$family == "weibull") {
    cat("Shape parameter (component 2):\n")
    cat(format(signif(x$shape2, digits)), "\n")
    cat("Scale parameter (component 2):\n")
    cat(format(signif(x$scale2, digits)), "\n")
  }

  invisible(x)
}
