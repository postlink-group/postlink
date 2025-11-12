#' Variance-Covariance Matrix from a `glm_mixture` Object
#' @description
#' Extract the variance-covariance matrix in the `glm_mixture` object for
#' the model parameter estimates.
#'
#' @param object an object of class `glm_mixture`.
#' @param ... for additional `vcov` arguments.
#'
#' @returns Return the variance-covariance matrix for the model parameter estimates.
#'
#' @export
vcov.glm_mixture <- function(object,...){
 return(object$vcov)
}

