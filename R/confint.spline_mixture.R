#' Confidence Intervals from a `spline_mixture` Object
#' @description Return confidence intervals for coefficients
#' from a `spline_mixture()` object
#'
#' @param x the result of a call to `spline_mixture()`
#' @param level the desired confidence level (default is 0.95)
#' @param ... for additional confint arguments
#'
#' @returns a matrix with the upper and lower confidence limits for each
#' outcome model coefficient. These limits are labelled as (1-level)/2 and
#' 1 - (1-level)/2 (e.g., 2.5\% and 97.5\% for level 0.95).
#'
#' @export
confint.spline_mixture <- function(x, level = 0.95){
 coef <- c(x$par_coef$mu)
 ses <- sqrt(diag(x$par_coef$Sigma))
 alpha <- 1 - level

 df <- nrow(x$B) - length(coef)
 vals <- cbind(lower = coef - ses * qt(1 - alpha/2, df),
               upper = coef + ses * qt(1 - alpha/2, df))
 colnames(vals) <- c(paste(alpha/2*100, "%"),
                     paste((1 - alpha/2)*100, "%"))
 return(vals)
}
