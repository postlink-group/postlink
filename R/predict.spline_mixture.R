#' Predictions from a `spline_mixture` Object
#' @description Return predicted values from an `spline_mixture()` object
#'
#' @param x the result of a call to `spline_mixture()`
#' @param ... for additional predict arguments
#'
#' @returns matrix of predictions w/ se and point-wise confidence bands
#'
#' @export
predict.spline_mixture <- function(x, newdata, se.fit = FALSE,
                                   interval = c("none", "confidence"),
                                   level = 0.95, ...){
 vals <- list()
 if(!missing(newdata)){
 gamobj <- gam(formula = x$call$formula, data = newdata, fit = FALSE)
 B <- as.matrix(gamobj$X)
 vals$vals <- B  %*% x$par_coef$mu
 } else{
 B <- x$B
 vals$vals <- B %*% x$par_coef$mu
 }

 if(se.fit){
  var <- B %*% x$par_coef$Sigma %*% t(B)
  vals$se <- sqrt(diag(var))
 }

 interval <- match.arg(interval)
 if(interval == "confidence"){
  alpha <- 1-level
  vals$low <- vals$vals - qnorm(1 - alpha/2)*vals$se
  vals$high <- vals$vals + qnorm(1 - alpha/2)*vals$se
  names(vals)[3:4] <- c(paste(alpha/2*100, "%"),
                        paste((1 - alpha/2)*100, "%"))
 }

 return(vals)
}
