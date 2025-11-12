#' Print a `coxph_ele` Object
#' @description Print call and outcome model coefficients from a `coxph_ele()` object
#'
#' @param x the result of a call to `coxph_ele()`
#' @param digits the number of significant digits to print
#' @param ... for additional print arguments
#'
#' @returns invisibly returns the `coxph_ele()` object that is provided as an argument
#'
#' @export
print.coxph_ele <- function(x, digits = max(3L, getOption("digits") - 3L),...){
 cat("Call:\n")
 print(x$call, quote = F, digits = digits)
 cat("\n")

 printCoefmat(cbind(coef = x$coefficients, "exp(coef)" = exp(x$coefficients)),
              quote=F, digits = digits, has.Pvalue = FALSE)
 cat("\n")

 invisible(x)
}
