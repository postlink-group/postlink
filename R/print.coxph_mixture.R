#' Print a `coxph_mixture` Object
#' @description Print call and outcome model coefficients from a `coxph_mixture()` object
#'
#' @param x the result of a call to `coxph_mixture()`
#' @param digits the number of significant digits to print
#' @param ... for additional print arguments
#'
#' @returns invisibly returns the `coxph_mixture()` object that is provided as an argument
#'
#' @export
print.coxph_mixture <- function(x, digits = max(3L, getOption("digits") - 3L),...){
 cat("Call:\n")
 print(x$call, quote = F, digits = digits)
 cat("\n")

 #printCoefmat(cbind(coef = x$coefficients, "exp(coef)" = exp(x$coefficients)),
#               quote=F, digits = digits, has.Pvalue = FALSE)
 l <- length(x$coefficients)
 l2 <- l + 1

 zval <- x$coefficients/x$standard.errors[1:l]
 pval <- 2 * (1 - pnorm(abs(zval)))

 e <- length(x$standard.errors)
 TAB <- cbind(x$coefficients, exp(x$coefficients),
              x$standard.errors[1:l],
              zval, pval)
 colnames(TAB) <- c("coef", "exp(coef)", "se(coef)", "z", "p")
 rownames(TAB) <- rownames(TAB)
 print(TAB, quote = F, digits = digits)
 cat("\n")

 invisible(x)
}
