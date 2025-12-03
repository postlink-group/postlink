#' Print a brief summary of a glm_mixtureBayes model
#'
#' @param x An object of class \code{glm_mixtureBayes}.
#' @param digits Minimum number of significant digits to show.
#' @param ... Further arguments (unused).
#' @return The input \code{x}, invisibly.
#' @export
print.glm_mixtureBayes <- function(x, digits = max(3L, getOption("digits") - 3L),...){
  cat("Call:\n")
  print(x$call, quote = F, digits = digits)
  cat("\n")
  
  cat("Coefficients:", sep="\n")
  vals <- apply(x$estimates$coefficients, 2, mean)
  print(format(signif(vals, digits)), print.gap = 2, quote = F)
  cat("\n")
  
  invisible(x)
}
