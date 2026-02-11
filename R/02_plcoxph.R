#' Fit Cox Proportional Hazards Models with Linkage Error Adjustment
#'
#' \code{plcoxph} fits Cox proportional hazards models to linked data, adjusting for
#' potential mismatch errors. It serves as a wrapper around the internal \code{fitcoxph}
#' engine, ensuring compatibility with the \code{\link[survival]{coxph}} syntax.
#'
#' @param formula A formula object, with the response on the left of a ~ operator,
#'   and the terms on the right. The response must be a survival object as returned
#'   by the \code{\link[survival]{Surv}} function.
#' @param adjustment An object inheriting from class \code{"adjustment"}, or a
#'   \code{list} containing the necessary parameter specifications.
#' @param subset An optional vector specifying a subset of observations.
#' @param na.action A function for handling NAs.
#' @param model Logical; if \code{TRUE}, the model frame is returned.
#' @param x,y Logical; if \code{TRUE}, the model matrix (\code{x}) and response
#'   (\code{y}) are returned. Defaults are \code{FALSE} and \code{FALSE}.
#' @param control A list of parameters for controlling the linkage error
#' adjustment process.
#' @param ... Additional arguments passed to the internal engine.
#'
#' @return A fitted model object containing the \code{call} and requested design components.
#'
#' @seealso \code{\link{plglm}}, \code{\link[survival]{coxph}}
#' @export
plcoxph <- function(formula,
                    adjustment,
                    subset,
                    na.action,
                    model = TRUE,
                    x = FALSE,
                    y = FALSE,
                    control = list(),
                    ...) {

 cl <- match.call()

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

 # Model Frame Construction
 mf <- match.call(expand.dots = FALSE)
 m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0L)
 mf <- mf[c(1L, m)]
 mf$drop.unused.levels <- TRUE
 mf[[1L]] <- quote(stats::model.frame)

 # Override if adjustment object has data
 if (!is.null(data_linked)) {
  mf$data <- data_linked
 }

 mf <- eval(mf, parent.frame())

 # Extract X and Y
 Y_obj <- stats::model.response(mf)
 if (!inherits(Y_obj, "Surv")) {
  stop("Response must be a 'Surv' object.", call. = FALSE)
 }

 mt <- attr(mf, "terms")
 X_mat <- stats::model.matrix(mt, mf)
 if (attr(mt, "intercept") == 1) X_mat <- X_mat[, -1, drop = FALSE]

 # Dispatch
 fit <- fitcoxph(x = X_mat,
                 y = Y_obj,
                 adjustment = adjustment,
                 control = control,
                 ...)

 # Post-Processing
 fit$call <- cl
 if (model) fit$model <- mf
 if (x) fit$x <- X_mat
 if (y) fit$y <- Y_obj

 class(fit) <- c(class(fit), "plcoxph", "coxph")

 return(fit)
}
