#' Secondary Analysis Constructor Based on Bayesian Mixture Modeling
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
#' @param priors A named \code{list} (or \code{NULL}) of prior specifications.
#'   Because the Stan models are pre-compiled, these strings are parsed into numeric
#'   hyperparameters and passed to the model's data block. Any missing entries are
#'   automatically filled with symmetric defaults dynamically during the model fitting phase.
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
#' @examples
#' # Simulate true data
#' set.seed(123)
#' n <- 200
#' x <- rnorm(n)
#' true_y <- 1.5 + 2 * x + rnorm(n)
#'
#' # Induce linkage mismatch errors
#' # Assume we have a match probability for each record
#' match_prob <- rbeta(n, 8, 2)
#' is_mismatch <- rbinom(n, 1, 1 - match_prob)
#'
#' obs_y <- true_y
#' mismatch_idx <- which(is_mismatch == 1)
#' if(length(mismatch_idx) > 1) {
#'   obs_y[mismatch_idx] <- sample(obs_y[mismatch_idx])
#' }
#' linked_data <- data.frame(y = obs_y, x = x, match_prob = match_prob)
#'
#' # Construct the Bayesian mixture adjustment object
#' adj_bayes <- adjMixBayes(
#'   linked.data = linked_data,
#'   priors = list(theta = "beta(2, 2)") # Optional: Override default
#' )
#'
#' class(adj_bayes)
#'
#' @seealso
#' * [plglm()] for generalized linear regression modeling
#' * [plsurvreg()] for parametric survival modeling
#'
#' @references
#' Gutman, R., Sammartino, C., Green, T., & Montague, B. (2016). Error adjustments for file
#' linking methods using encrypted unique client identifier (eUCI) with application to recently
#' released prisoners who are HIV+. \emph{Statistics in Medicine}, 35(1), 115–129. \doi{10.1002/sim.6586}
#'
#' @export
adjMixBayes <- function(linked.data = NULL, priors = NULL){
 # Validate linked.data
 if (!is.null(linked.data)) {
  if (is.environment(linked.data) || is.list(linked.data)) {
   linked.data <- tryCatch(as.data.frame(linked.data), error = function(e) {
    stop("'linked.data' must be a data.frame or coercible to one.", call. = FALSE)
   })
  } else if (!is.data.frame(linked.data)) {
   stop("'linked.data' must be a data.frame, list, or environment.", call. = FALSE)
  }
 }

 # Validate priors
 if (!is.null(priors) && !is.list(priors)) {
  stop("'priors' must be a named list or NULL.", call. = FALSE)
 }

 # Construct and Return the S3 Object with Reference Semantics
 data_ref <- new.env(parent = emptyenv())
 data_ref$data <- linked.data

 out <- structure(
  list(
   data_ref = data_ref,
   priors = priors
  ),
  class = c("adjMixBayes", "adjustment")
 )

 return(out)
}
