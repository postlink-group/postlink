#' Summarize a surv_mixtureBayes model fit
#'
#' @param object An object of class \code{surv_mixtureBayes}.
#' @param ... Not used.
#' @return An object of class \code{"summary.surv_mixtureBayes"}, which is printed with a custom method.
#' @export
summary.surv_mixtureBayes <- function(object, ...) {

 # Calculate posterior summary for coefficients of each component
 coef_draws1 <- object$estimates$coefficients   # S x p matrix
 coef_draws2 <- object$estimates$m.coefficients

 # Credible intervals for coefficients
 ci1 <- t(apply(coef_draws1, 2, stats::quantile, probs = c(0.025, 0.975)))
 ci2 <- t(apply(coef_draws2, 2, stats::quantile, probs = c(0.025, 0.975)))

 TAB1 <- cbind(
  Estimate    = apply(coef_draws1, 2, mean),
  `Std. Error`= apply(coef_draws1, 2, stats::sd),
  `2.5 %`     = ci1[, 1],
  `97.5 %`    = ci1[, 2]
 )

 TAB2 <- cbind(
  Estimate    = apply(coef_draws2, 2, mean),
  `Std. Error`= apply(coef_draws2, 2, stats::sd),
  `2.5 %`     = ci2[, 1],
  `97.5 %`    = ci2[, 2]
 )

 # Prepare output list
 res <- list(
  call          = object$call,
  family        = object$family,
  coefficients  = TAB1,
  m.coefficients= TAB2
 )

 # Family-specific parameters
 if (object$family == "gamma") {
  res$shape   <- mean(object$estimates$shape)
  res$m.shape <- mean(object$estimates$m.shape)
 }

 if (object$family == "weibull") {
  res$shape  <- c(comp1 = mean(object$estimates$shape),
                  comp2 = mean(object$estimates$m.shape))
  res$scale  <- c(comp1 = mean(object$estimates$scale),
                  comp2 = mean(object$estimates$m.scale))
  res$shape1 <- mean(object$estimates$shape)
  res$shape2 <- mean(object$estimates$m.shape)
  res$scale1 <- mean(object$estimates$scale)
  res$scale2 <- mean(object$estimates$m.scale)
 }

 class(res) <- "summary.surv_mixtureBayes"
 res
}
