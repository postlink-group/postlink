#' Predictions from a `coxph_mixture` Object
#' @description Obtain predictions from a `coxph_mixture()` object using
#' `predict.coxph()`.
#'
#' @param object the result of a call to `coxph_mixture()`
#' @param newdata optional new data to obtain predictions for. The original data
#' is used by default.
#' @param type For the "cox" family, the choices are the linear predictor ("lp"),
#' the risk score exp(lp) ("risk"), the expected number of events given the
#' covariates and follow-up time ("expected"), and the terms of the linear
#' predictor ("terms"). The survival probability for a subject is equal to
#' exp(-expected).
#' @param terms the terms when type = "terms". By default, all terms are included.
#' @param collapse optional vector of subject identifiers. If specified, the
#' output will contain one entry per subject rather than one entry per observation.
#' @param na.action a function for what to do with missing values in `newdata`.
#' The default is to predict "NA".
#' @param reference Reference for centering predictions.
#' Available options are c("strata" - default, "sample", "zero").
#' @param ... for future predict arguments
#'
#' @returns a vector or matrix of predictions based on arguments specified.
#'
#' @export
predict.coxph_mixture <- function(object, newdata,
                                  type=c("lp", "risk", "expected", "terms", "survival"),
                                  na.action=na.pass, terms = NULL, collapse = NULL,
                                  reference=c("strata", "sample", "zero"),...){
  coxph_object <- coxph(fit$wfit$y ~ fit$wfit$x, weights = fit$wfit$weights)

  type <- match.arg(type)
  reference <- match.arg(reference)
  terms <- ifelse(is.null(terms), names(coxph_object$assign), terms)

  if(is.null(collapse)){
    if(missing(newdata)){
      predictions <- predict(object = coxph_object,
                            type = type, se.fit = FALSE, na.action = na.action,
                            terms = terms, reference = reference)
    } else{
      predictions <- predict(object = coxph_object, newdata = newdata,
                            type = type, se.fit = FALSE, na.action = na.action,
                            terms = terms, reference = reference)
    }
  } else{
    if(missing(newdata)){
      predictions <- predict(object = coxph_object,
                            type = type, se.fit = FALSE, na.action = na.action,
                            terms = terms, collapse = collapse,
                            reference = reference)
    } else{
      predictions <- predict(object = coxph_object, newdata = newdata,
                            type = type, se.fit = FALSE, na.action = na.action,
                            terms = terms, collapse = collapse,
                            reference = reference)
    }
  }

  return(predictions)
}
