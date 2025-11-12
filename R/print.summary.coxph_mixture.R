#' @noRd
#' @export
print.summary.coxph_mixture <- function(x, digits = max(3L, getOption("digits") - 3L),
                                     signif.stars = getOption("show.signif.stars"),...){
 cat("Call:", sep="\n")
 print(x$call,quote=F)
 cat(" ", sep="\n")
 cat("Family:", x$family, " ", sep="\n")
 
 cat("Outcome Model Coefficients:", sep="\n")
 printCoefmat(x$coefficients,quote=F, digits = digits,
              signif.stars = signif.stars)
 cat(" ", sep="\n")
 
 cat("Correct Match Model Coefficients:", sep="\n")
 printCoefmat(x$m.coefficients,quote=F, digits = digits,
              signif.stars = signif.stars)
 cat(" ", sep="\n")
 
 cat("Average Correct Match Rate: ", format(signif(x$avgcmr, digits)))
 cat(" ", sep="\n")
 
 invisible(x)
}