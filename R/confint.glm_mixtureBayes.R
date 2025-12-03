#' Credible intervals for glm_mixtureBayes coefficients
#'
#' @param object A \code{glm_mixtureBayes} model object.
#' @param parm Not used (returns intervals for all coefficients of component 1).
#' @param level Probability level for the intervals (default 0.95).
#' @param ... Not used.
#' @return A matrix with two columns (lower and upper bounds) and one row per coefficient (component 1).
#' @export

confint.glm_mixtureBayes <- function(object, level = 0.95,...){
  alpha <- 1-level
  vals <- t(apply(object$estimates$coefficients, 2, 
                function(x) quantile(x, probs = c(alpha/2, 1-alpha/2))))
  return(vals)              
}
