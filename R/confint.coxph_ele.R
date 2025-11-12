#' Confidence intervals from a `coxph_ele` Object.
#' @description
#' Compute confidence intervals for the outcome model coefficients.
#'
#' @param object an object of class "coxph_ele".
#' @param level the required confidence level.
#' @param ... for additional `confint` arguments.
#'
#' @returns a data frame with the upper and lower confidence limits for each
#' outcome model coefficient. These limits are labelled as (1-level)/2 and
#' 1 - (1-level)/2 (e.g., 2.5% and 97.5% for level 0.95).
#'
#' @export
confint.coxph_ele <- function(object, level = 0.95,...){
 coef <- object$coefficients
 ses <- object$standard.errors
 alpha <- 1 - level
 vals <- data.frame(coef - ses * qnorm(1 - alpha/2),
                    coef + ses * qnorm(1 - alpha/2))
 colnames(vals) <- c(paste(alpha/2*100, "%"),
                     paste((1 - alpha/2)*100, "%"))
 return(vals)
}
