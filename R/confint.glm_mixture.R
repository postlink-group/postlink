#' Confidence intervals from a `glm_mixture` Object.
#' @description
#' Compute confidence intervals for the outcome model coefficients.
#'
#' @param object an object of class `glm_mixture`.
#' @param level the required confidence level.
#' @param ... for additional `confint` arguments.
#'
#' @returns a data frame with the upper and lower confidence limits for each
#' outcome model coefficient. These limits are labelled as (1-level)/2 and
#' 1 - (1-level)/2 (e.g., 2.5\% and 97.5\% for level 0.95).
#'
#' @export
confint.glm_mixture <- function(object, level = 0.95,...){
 coef <- object$coefficients
 ses <- object$standard.errors[1:length(coef)]
 alpha <- 1 - level
 t <- qt(1 - alpha/2, df = object$wfit$df.residual)
 vals <- data.frame(coef - ses * t,
                    coef + ses * t)
 colnames(vals) <- c(paste(alpha/2*100, "%"),
                     paste((1 - alpha/2)*100, "%"))
 return(vals)
}
