#' @noRd
print.summary.glm_ele <- function(x, digits = max(3L, getOption("digits") - 3L),
                                      signif.stars = getOption("show.signif.stars"),...){
 cat("Call:", sep="\n")
 print(x$call,quote=F)
 cat(" ", sep="\n")
 cat("Family:", x$family, " ", sep="\n")

 cat("Outcome Model Coefficients:", sep="\n")
 for(i in 1:length(x$coefficients)){
 print(names(x$coefficients)[i])
 printCoefmat(x$coefficients[[i]],quote=F, digits = digits,
              signif.stars = signif.stars)
 cat(" ", sep="\n")
 }

 if (x$family == "gamma"){
  cat("Dispersion: ", format(signif(as.numeric(x$dispersion), digits)))
  cat("\n")
  cat("\n")
 }

 if(x$family == "poisson" | x$family == "binomial"){
  cat("Dispersion: ", 1)
  cat("\n")
  cat("\n")
 }

 if (x$family == "gaussian"){
  cat("Dispersion:", sep="\n")
  print(x$dispersion, quote=F, digits = digits)
  cat(" ", sep="\n")
 }

 invisible(x)
}
