#' Variance-Covariance Matrix from a `ctable_mixture` Object
#' @description
#' Extract the variance-covariance matrix for the estimated cell probabilities
#' in the `ctable_mixture` object.
#'
#' @param object an object of class "ctablemixture".
#' @param ... for additional `vcov` arguments.
#'
#' @returns A matrix of the estimated covariances between the cell probabilities.
#' This should have row and column names corresponding to the row and column
#' variables’ levels of the cells in the underlying two-way contingency table.
#'
#' @export
vcov.ctablemixture <- function(object,...){
 return(object$vcov_phat)
}
