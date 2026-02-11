#' Constructor for Secondary Analysis of Linked Data Based on Bayesian Mixture Modeling
#'
#' Specifies the linked data and information on the underlying record linkage
#' process for a mismatch error adjustment using a Bayesian framework based on
#' mixture modeling as developed by Gutman et al. (2016). This framework uses
#' a mixture model for pairs of linked records whose two components reflect distributions
#' conditional on match status, i.e., correct match or mismatch.
#' Posterior inference is carried out via data augmentation or multiple
#' imputation.
#'
#' @param linked.data A data.frame containing the linked dataset.
#'
#' @return An object of class \code{c("adjMixBayes", "adjustment")}. To minimize
#' memory overhead, the underlying \code{linked.data} is stored by reference
#' within an environment inside this object.
#'
#' @details
#' Explicit provision of \code{linked.data} is strongly recommended for
#' reproducibility and to ensure the adjustment object fully encapsulates
#' the necessary data for downstream model fitting.
#'
#' #' @examples
#' # Example: Using the included brfss demonstration dataset
#' data(brfss, package = "postlink")
#'
#' adj_object <- adjMixBayes(linked.data = brfss)
#'
#' @note
#' The reference below discusses the implemented framework in more detail.
#'
#' @references
#' Gutman, R., Sammartino, C., Green, T., & Montague, B. (2016). Error adjustments for file
#' linking methods using encrypted unique client identifier (eUCI) with application to recently
#' released prisoners who are HIV+. \emph{Statistics in Medicine}, 35(1), 115–129. \doi{10.1002/sim.6586}
#'
#' @export
adjMixBayes <- function(linked.data = NULL){
 # 1. Validate linked.data
 if (!is.null(linked.data)) {
  if (is.environment(linked.data) || is.list(linked.data)) {
   linked.data <- tryCatch(as.data.frame(linked.data), error = function(e) {
    stop("'linked.data' must be a data.frame or coercible to one.", call. = FALSE)
   })
  } else if (!is.data.frame(linked.data)) {
   stop("'linked.data' must be a data.frame, list, or environment.", call. = FALSE)
  }
 }

 # 2. Construct and Return the S3 Object with Reference Semantics
 data_ref <- new.env(parent = emptyenv())
 data_ref$data <- linked.data

 out <- structure(
  list(
   data_ref = data_ref
  ),
  class = c("adjMixBayes", "adjustment")
 )

 return(out)
}
