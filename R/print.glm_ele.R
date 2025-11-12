#' Print a `glm_ele` Object
#' @description Print call and outcome model coefficients from a `glm_ele()` object
#'
#' @param x the result of a call to `glm_ele()`
#' @param digits the number of significant digits to print
#' @param ... for additional print arguments
#'
#' @returns invisibly returns the `glm_ele()` object that is provided as an argument
#' print(fit)
#'
#' @export
print.glm_ele <- function(x, digits = max(3L, getOption("digits") - 3L),...){
 cat("Call:\n")
 print(x$call, quote = F, digits = digits)
 cat("\n")

 cat("Coefficients:", sep="\n")
 print(format(signif(x$coefficients, digits)), print.gap = 2, quote = F)
 cat("\n")

 invisible(x)
}
