#' Summarize a `coxph_mixture` Object
#' @description Summarize results from a `coxph_mixture()` object
#'
#' @param object the result of a call to `coxph_mixture()`
#' @param ... for additional summary arguments
#'
#' @returns a list of results from the function called depending on the "family" specified.
#' \item{call}{the matched call}
#' \item{family}{the assumed type of (outcome) regression model}
#' \item{coefficients}{a matrix with the outcome model's coefficient estimates, standard errors, t or z values, and p-values}
#' \item{m.coefficients}{a matrix with the correct match model's coefficient estimates and standard errors}
#' \item{avgcmr}{the average correct match rate among all records}
#' \item{match.prob}{the posterior correct match probabilities for observations given parameter estimates}
#'
#' @export
summary.coxph_mixture <- function(object,...){
 l <- length(object$coefficients)
 l2 <- l + 1

  zval <- object$coefficients/object$standard.errors[1:l]
  pval <- 2 * (1 - pnorm(abs(zval)))

  e <- length(object$standard.errors)
  TAB <- cbind(object$coefficients, exp(object$coefficients),
               object$standard.errors[1:l],
               zval, pval)
  colnames(TAB) <- c("coef", "exp(coef)", "se(coef)", "z value", "Pr(>|z|)")
  rownames(TAB) <- rownames(TAB)
  TAB2 <- cbind(object$m.coefficients, object$standard.errors[l2:e])
  colnames(TAB2) <- c("Estimate","Std. Error")

 object <- list(call = object$call, family = object$family,
                coefficients = TAB, m.coefficients = TAB2,
                avgcmr = mean(object$match.prob), match.prob = object$hs)

 class(object)    <- "summary.coxph_mixture"
 object
}
