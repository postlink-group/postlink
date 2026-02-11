#' Analysis of Contingency Tables with Linkage Error Adjustment
#'
#' \code{plctable} constructs a contingency table from linked data and fits a model
#' adjusting for linkage errors.
#'
#' @param formula a formula object with the left and right hand sides specifying the column
#' and row variable of the flat table, respectively.
#' @param adjustment An object inheriting from class \code{"adjustment"}, or a
#'   \code{list} containing the necessary parameter specifications.
#' @param subset An optional vector specifying a subset of observations to be used.
#' @param na.action A function which indicates what should happen when the data contain NAs.
#' @param exclude Vector of values to be excluded when forming the table (passed to \code{xtabs}).
#' @param control A list of parameters for controlling the linkage error
#' adjustment process.
#' @param ... Additional arguments passed to the internal engine.
#'
#' @return A fitted model object tailored for contingency table analysis.
#'
#' @export
plctable <- function(formula,
                     adjustment,
                     subset,
                     na.action,
                     exclude = c(NA, NaN),
                     control = list(),
                     ...) {

 cl <- match.call()

 if (missing(data)) data <- NULL

 # Data Retrieval
 if (missing(adjustment) || (!inherits(adjustment, "adjustment") && !is.list(adjustment))) {
  stop("'adjustment' must be a valid adjustment object or a list.", call. = FALSE)
 }

 data_linked <- NULL
 if (!is.null(adjustment$data_ref) && exists("data", envir = adjustment$data_ref)) {
  data_linked <- adjustment$data_ref$data
 } else if (is.list(adjustment) && !is.null(adjustment$data) && is.data.frame(adjustment$data)) {
  data_linked <- adjustment$data
 }

 # Create Contingency Table
 xt_call <- call("xtabs", formula = formula, exclude = exclude)

 # Priority: Adjustment Data > Argument Data
 if (!is.null(data_linked)) {
  xt_call$data <- data_linked
 } else if (!missing(data)) {
  xt_call$data <- substitute(data)
 }

 if (!missing(subset)) xt_call$subset <- substitute(subset)
 if (!missing(na.action)) xt_call$na.action <- substitute(na.action)

 ftable <- eval(xt_call, envir = parent.frame())

 # Dispatch
 fit <- fitctable(ftable = ftable,
                  adjustment = adjustment,
                  control = control,
                  ...)

 # Post-Processing & Class Assignment
 fit$call <- cl
 fit$formula <- formula

 # Append package-level class "plctable" to the engine-level class
 class(fit) <- c(class(fit), "plctable")

 return(fit)
}
