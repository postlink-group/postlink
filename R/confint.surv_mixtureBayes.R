#' Credible intervals for surv_mixtureBayes coefficients
#'
#' @param object A \code{surv_mixtureBayes} model object.
#' @param parm Not used (returns intervals for all coefficients of component 1).
#' @param level Probability level for the intervals (default 0.95).
#' @param ... Not used.
#' @return A matrix with two columns (lower and upper bounds) and one row per coefficient (component 1).
#' @export
confint.surv_mixtureBayes <- function(object, parm, level = 0.95, ...) {
 alpha <- 1 - level
 coef_draws <- object$estimates$coefficients  # S x p matrix for comp1
 ci_mat <- t(apply(coef_draws, 2, stats::quantile, probs = c(alpha/2, 1 - alpha/2)))
 colnames(ci_mat) <- c(
  paste0(round(alpha/2*100, 1), "%"),
  paste0(round((1 - alpha/2)*100, 1), "%")
 )
 return(ci_mat)
}

