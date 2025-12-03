#' Confidence intervals from a `ctable_mixture` Object.
#' @description
#' Compute confidence intervals for the cell probabilities from the fitted model.
#'
#' @param object an object of class "ctablemixture".
#' @param level the required confidence level.
#' @param ... for additional `confint` arguments.
#'
#' @returns a data frame with the upper and lower confidence limits for each
#' cell probability. These limits are labelled as (1-level)/2 and
#' 1 - (1-level)/2 (e.g., 2.5\% and 97.5\% for level 0.95).
#'
#' @export
confint.ctablemixture <- function(object, level = 0.95,...){
 phat <- c(t(object$phat))
 ses <- sqrt(diag(object$vcov_phat))
 alpha <- 1 - level
 vals <- data.frame(expand.grid(colnames(object$phat), rownames(object$phat)),
                    phat - ses * qnorm(1 - alpha/2),
                    phat + ses * qnorm(1 - alpha/2))
 colnames(vals) <- c(names(dimnames(object$phat))[2:1], paste(alpha/2*100, "%"),
                     paste((1 - alpha/2)*100, "%"))
 return(vals)
}
