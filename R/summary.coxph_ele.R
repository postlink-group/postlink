#' Summarize a `coxph_ele` Object
#' @description Summarize results from a `coxph_ele()` object
#'
#' @param object the result of a call to `coxph_ele()`
#' @param ... for additional summary arguments
#'
#' @returns a list of results from the function called depending on the "family" specified.
#' \item{call}{the matched call}
#' \item{coefficients}{a matrix with the outcome model's coefficient estimates, standard errors, t or z values, and p-values}
#'
#' @export
summary.coxph_ele <- function(object,...){
 l <- length(object$coefficients)
 l2 <- l + 1

 zval <- object$coefficients/object$standard.errors[1:l]
 pval <- 2 * (1 - pnorm(abs(zval)))

 e <- length(object$standard.errors)
 TAB <- cbind(object$coefficients, exp(object$coefficients),
              object$standard.errors[1:l],
              zval, pval)
 colnames(TAB) <- c("coef", "exp(coef)", "se(coef)", "z value", "Pr(>|z|)")
 rownames(TAB) <- substring(rownames(TAB), first=2)

 object <- list(call = object$call,
                coefficients = TAB)

 class(object)    <- "summary.coxph_ele"
 object
}
