#' Predictions from a `coxph_ele` Object
#' @description Obtain predictions from a `coxph_ele()` object.
#'
#' @param object the result of a call to `coxph_ele()`
#' @param newdata optional new data to obtain predictions for. The original data
#' is used by default.
#' @param type For the "cox" family, the choices are the linear predictor ("lp"),
#' and the risk score exp(lp) ("risk").
#' @param ... for future predict arguments
#'
#' @returns a vector of predictions based on arguments specified.
#'
#' @export
predict.coxph_ele <- function(object, newdata, type = c("lp", "risk")){

 if (!missing(newdata)){
  coef <- object$coefficients
  mf <- model.frame(fit$call[["formula"]], data = newdata)
  Z <- mf[,-1]
  lp <- c(X %*% coef)
 } else{
  lp <- object$linear.predictors
 }

 colnames(lp) <- rownames(object$coefficients)
 type <- match.arg(type)
 if(type == "lp"){
  return(lp)
 } else if(type == "risk"){
  return(exp(lp))
 } else{
  stop("Error: type should be lp or risk")
 }

}
