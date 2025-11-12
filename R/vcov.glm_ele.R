#' Variance-Covariance Matrix from a `glm_ele` Object
#' @description
#' Extract the variance-covariance matrix in the `glm_ele` object for
#' the model parameter estimates.
#'
#' @param object an object of class "glm_ele".
#' @param ... for additional `vcov` arguments.
#'
#' @returns Return the variance-covariance matrix for the model parameter estimates.
#'
#' @export
vcov.glm_ele <- function(object,...){
 return(object$vcov)
}
