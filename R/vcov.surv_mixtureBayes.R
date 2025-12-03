#' Covariance matrix of coefficient estimates for surv_mixtureBayes
#'
#' @param object A \code{surv_mixtureBayes} model object.
#' @param ... Not used.
#' @return Posterior covariance matrix of component 1's coefficient vector.
#' @export
vcov.surv_mixtureBayes <- function(object, ...) {
  coef_draws <- object$estimates$coefficients  # S x p (comp1 coefficients)
  return(cov(coef_draws))
}
