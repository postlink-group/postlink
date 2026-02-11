#' Fit Generalized Linear Models with Linkage Error Adjustment
#'
#' \code{plglm} fits generalized linear models (GLMs) to linked data, incorporating
#' adjustments for linkage error as specified in the provided \code{adjustment} object.
#' It mimics the interface of \code{\link[stats]{glm}} to ensure familiarity for users.
#'
#' @param formula A symbolic description of the model to be fitted.
#' @param family A description of the error distribution and link function to be used
#'   in the model. This can be a character string naming a family function, a family
#'   function or the result of a call to a family function.
#' @param adjustment An object inheriting from class \code{"adjustment"} (e.g.,
#'   \code{adjELE}, \code{adjMixture}), or a \code{list} containing the necessary
#'   data and parameter specifications.
#' @param subset An optional vector specifying a subset of observations to be used
#'   in the fitting process.
#' @param na.action A function which indicates what should happen when the data
#'   contain NAs. The default is set by the \code{na.action} setting of \code{options},
#'   and is \code{\link[stats]{na.fail}} if that is unset.
#' @param model Logical; if \code{TRUE} (default), the model frame is returned.
#' @param x,y Logical; if \code{TRUE}, the model matrix (\code{x}) and response
#'   vector (\code{y}) are returned. Default is \code{FALSE} for both.
#' @param control A list of parameters for controlling the linkage error
#' adjustment process.
#' @param ... Additional arguments passed to the underlying fitting function.
#'
#' @return An object representing the fitted model. The exact class depends on the
#'   class of the \code{adjustment} object and the specific internal engine invoked
#'   (e.g., \code{fitglm}). The object includes the component \code{call} and,
#'   optionally, \code{model}, \code{x}, and \code{y}.
#'
#' @details
#' This function attempts to extract the linked data from the \code{adjustment}
#' object. It supports both reference-based storage (via \code{data_ref}) and
#' direct list components (\code{adjustment$data}). If the data is not present
#' (e.g., NULL), the function will attempt to resolve variables from the
#' environment of the \code{formula}.
#'
#' It applies the standard \code{model.frame} processing steps (formula parsing,
#' subsetting, NA handling) and dispatches the resulting design matrix and
#' response vector to the appropriate \code{fitglm} method.
#'
#' @examples
#' \dontrun{
#' # Method 1: Using a constructor (Recommended)
#' adj_obj <- adjELE(linked.data = my_data, m.rate = 0.1)
#' fit1 <- plglm(y ~ x1 + x2, family = binomial, adjustment = adj_obj)
#'
#' # Method 2: Manual list specification
#' # Note: You must ensure the list structure matches what the internal engine expects
#' adj_list <- list(data = my_data, m.rate = 0.1)
#' class(adj_list) <- "adjELE" # Manually assign class for S3 dispatch
#' fit2 <- plglm(y ~ x1 + x2, family = binomial, adjustment = adj_list)
#' }
#'
#' @seealso \code{\link{adjELE}}, \code{\link{adjMixture}}, \code{\link{adjMixBayes}}
#' @export
plglm <- function(formula,
                  family = gaussian,
                  adjustment,
                  subset,
                  na.action,
                  model = TRUE,
                  x = FALSE,
                  y = FALSE,
                  control = list(),
                  ...) {

 # Capture Call
 cl <- match.call()

 # Validate Family (Robust Lookup)
 if (is.character(family)) {
  family <- tryCatch({
   get(family, mode = "function", envir = parent.frame())
  }, error = function(e) NULL)
 }
 if (is.function(family)) {
  family <- family()
 }
 if (is.null(family$family)) {
  stop("'family' not recognized", call. = FALSE)
 }

 # Data Retrieval from Adjustment Object
 if (missing(adjustment) || (!inherits(adjustment, "adjustment") && !is.list(adjustment))) {
  stop("'adjustment' must be a valid adjustment object or a list.", call. = FALSE)
 }

 data_linked <- NULL
 if (!is.null(adjustment$data_ref) && exists("data", envir = adjustment$data_ref)) {
  data_linked <- adjustment$data_ref$data
 } else if (is.list(adjustment) && !is.null(adjustment$data) && is.data.frame(adjustment$data)) {
  data_linked <- adjustment$data
 }

 # Standard Model Frame Construction
 mf <- match.call(expand.dots = FALSE)
 # Added "data" to match() so it is captured from arguments if provided
 m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0L)
 mf <- mf[c(1L, m)]
 mf$drop.unused.levels <- TRUE
 mf[[1L]] <- quote(stats::model.frame)

 # 1. If adjustment has data, force it (overrides 'data' arg).
 # 2. If adjustment has NO data, 'mf$data' remains what was passed in 'data=' arg.
 # 3. If 'data=' arg is also missing, model.frame uses formula environment.
 if (!is.null(data_linked)) {
  mf$data <- data_linked
 }

 # Evaluate model frame
 mf <- eval(mf, parent.frame())

 # Extract Model Matrix (X) and Response (Y)
 mt <- attr(mf, "terms")
 X_mat <- stats::model.matrix(mt, mf)
 Y_vec <- stats::model.response(mf, "any")

 # Dispatch to Internal Engine
 fit <- fitglm(x = X_mat,
               y = Y_vec,
               family = family,
               adjustment = adjustment,
               control = control,
               ...)

 # Post-Processing
 fit$call <- cl
 if (model) fit$model <- mf
 if (x) fit$x <- X_mat
 if (y) fit$y <- Y_vec

 # class hierarchy
 # 1. Engine Class (e.g., "glmMixture") - preserved from fitglm
 # 2. Package Class ("plglm") - added here
 # 3. Standard Classes ("glm", "lm") - added here for compatibility
 class(fit) <- c(class(fit), "plglm", "glm", "lm")

 return(fit)
}
