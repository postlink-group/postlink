#' @noRd
print.summary.coxph_ele <- function(x, digits = max(3L, getOption("digits") - 3L),
                                        signif.stars = getOption("show.signif.stars"),...){
 cat("Call:", sep="\n")
 print(x$call,quote=F)
 cat(" ", sep="\n")

 cat("Outcome Model Coefficients:", sep="\n")
 printCoefmat(x$coefficients,quote=F, digits = digits,
              signif.stars = signif.stars)
 cat(" ", sep="\n")

 invisible(x)
}
