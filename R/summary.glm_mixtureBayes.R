#' Summarize a glm_mixtureBayes model fit
#'
#' @param object An object of class \code{glm_mixtureBayes}.
#' @param ... Not used.
#' @return An object of class \code{"summary.glm_mixtureBayes"}, which is printed with a custom method.
#' @export
summary.glm_mixtureBayes <- function(object, ...) {

 ci_beta1 <- t(apply(object$estimates$coefficients, 2,
                     function(x) stats::quantile(x, probs = c(0.025, 0.975))))
 ci_beta2 <- t(apply(object$estimates$m.coefficients, 2,
                     function(x) stats::quantile(x, probs = c(0.025, 0.975))))

 TAB1 <- cbind(
  Estimates   = apply(object$estimates$coefficients, 2, mean),
  `Std. Error`= apply(object$estimates$coefficients, 2, stats::sd),
  `2.5 %`     = ci_beta1[, 1],
  `97.5 %`    = ci_beta1[, 2]
 )

 TAB2 <- cbind(
  Estimates   = apply(object$estimates$m.coefficients, 2, mean),
  `Std. Error`= apply(object$estimates$m.coefficients, 2, stats::sd),
  `2.5 %`     = ci_beta2[, 1],
  `97.5 %`    = ci_beta2[, 2]
 )

 if (object$family %in% c("gaussian", "gamma")) {
  TAB3 <- cbind(
   Estimate    = mean(object$estimates$dispersion),
   `Std. Error`= stats::sd(object$estimates$dispersion)
  )
  TAB4 <- cbind(
   Estimate    = mean(object$estimates$m.dispersion),
   `Std. Error`= stats::sd(object$estimates$m.dispersion)
  )
  rownames(TAB3) <- ""
  rownames(TAB4) <- ""
 }

 out <- list(
  call          = object$call,
  family        = object$family,
  coefficients  = TAB1,
  m.coefficients= TAB2
 )

 if (object$family %in% c("gaussian", "gamma")) {
  out$dispersion    <- TAB3
  out$m.dispersion  <- TAB4
 }

 class(out) <- "summary.glm_mixtureBayes"
 out
}
