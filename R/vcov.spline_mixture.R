#' Variance-Covariance Matrix from a `spline_mixture` Object
#' @description Extract the variance-covariance matrix
#'
#' @param x the result of a call to `spline_mixture()`
#' @param ... for additional vcov arguments
#'
#' @returns variance-covariance matrix
#'
#' @export
vcov.spline_mixture <- function(x,...){
 return(x$par_coef$Sigma)
}
