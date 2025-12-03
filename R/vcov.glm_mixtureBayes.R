#' Covariance matrix of coefficient estimates for glm_mixtureBayes
#'
#' @param object A \code{glm_mixtureBayes} model object.
#' @param ... Not used.
#' @return Posterior covariance matrix of component 1's coefficient vector.
#' @export
vcov.glm_mixtureBayes <- function(object,...){
  return(cov(object$estimates$coefficients))
}