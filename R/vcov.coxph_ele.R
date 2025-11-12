#' Variance-Covariance Matrix from a `coxph_ele` Object
#' @description
#' Extract the variance-covariance matrix in the `coxph_ele` object for
#' the model parameter estimates.
#'
#' @param object an object of class "coxph_ele".
#' @param ... for additional `vcov` arguments.
#'
#' @returns Return the variance-covariance matrix for the model parameter estimates.
#'
#' @export
vcov.coxph_ele <- function(object,...){
 return(object$vcov)
}
