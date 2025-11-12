#' Variance-Covariance Matrix from a `coxph_mixture` Object
#' @description
#' Extract the variance-covariance matrix in the `coxph_mixture` object for
#' the model parameter estimates.
#'
#' @param object an object of class "coxph_mixture".
#' @param ... for additional `vcov` arguments.
#'
#' @returns Return the variance-covariance matrix for the model parameter estimates.
#'
#' @export
vcov.coxph_mixture <- function(object,...){
 return(object$vcov)
}
