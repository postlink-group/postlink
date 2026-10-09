#' Posterior Draws of Bayesian Mixture Fits as a Single Matrix
#'
#' Collects every stored posterior draw of a \code{glmMixBayes} or
#' \code{survMixBayes} fit into one draws x parameters matrix. Draws are in
#' chain order (after burn-in and thinning), so the matrix can be handed to
#' MCMC diagnostics such as \code{coda::effectiveSize()} or
#' \code{coda::traceplot()}. \code{coda::as.mcmc()} returns the same draws as a
#' \pkg{coda} object (the methods are registered when \pkg{coda} is loaded).
#'
#' @section Column names:
#' Each column is named after the element of \code{x$estimates} it comes from,
#' followed by the parameter name in brackets for blocks with several
#' parameters (the column names of the design matrix \code{X} or of the
#' linkage covariates \code{Z}; \code{glmMixBayes()} and
#' \code{survregMixBayes()} name unnamed columns \code{"(Intercept)"} (a
#' column of ones), \code{"X<j>"} or \code{"Z<j>"}, and in fits saved by
#' earlier versions columns without names are named by their positions). The
#' blocks appear in this order, when present in the fit:
#' \describe{
#'   \item{\code{coefficients[...]}}{regression coefficients of component 1
#'     (the correct matches; for Weibull fits with an intercept column, the
#'     intercept is the identified \code{(Intercept) + log(scale)}, as in
#'     \code{coefficients2[...]});}
#'   \item{\code{coefficients2[...]}}{regression coefficients of component 2
#'     (the mismatches);}
#'   \item{\code{m.coefficients[...]}}{coefficients of the mismatch-indicator
#'     model on the logit scale of the mismatch probability, as the
#'     \code{m.coefficients} of \code{\link{glmMixture}()} and
#'     \code{\link{coxphMixture}()} (\code{-gamma}, or
#'     \code{m.coefficients[(Intercept)]} \code{= qlogis(1 - theta)} without
#'     linkage covariates);}
#'   \item{\code{theta}}{the mixing weight of component 1 (without linkage
#'     covariates);}
#'   \item{\code{gamma[...]}}{the coefficients of the logistic regression of
#'     the match probability (with linkage covariates, for their original,
#'     uncentred values);}
#'   \item{\code{dispersion}, \code{dispersion2}}{the dispersions of components
#'     1 and 2 (Gaussian and Gamma GLMs);}
#'   \item{\code{shape}, \code{shape2}}{the shapes of components 1 and 2
#'     (survival models);}
#'   \item{\code{scale}, \code{scale2}}{the Weibull scale multipliers of
#'     components 1 and 2 (not identified separately from the intercepts when
#'     the design matrix has an intercept column);}
#'   \item{\code{lp}}{the log posterior of each reported draw, up to a
#'     constant, with the latent match indicators integrated out (for a draw
#'     whose labels were exchanged after sampling, that of the relabelled
#'     draw; see \code{diagnostics$lp} in \code{\link{glmMixBayes}}).}
#' }
#' With linkage covariates, \code{gamma[...]} and \code{m.coefficients[...]}
#' hold the same parameters with opposite signs, so the draws have a singular
#' covariance matrix: compute diagnostics that need to invert it per parameter
#' or on a subset of the columns (e.g.
#' \code{coda::gelman.diag(..., multivariate = FALSE)}).
#' Up to postlink 0.1.2 the parameters of component 2 were named
#' \code{m.coefficients}, \code{m.dispersion}, \code{m.shape} and
#' \code{m.scale}.
#'
#' @param x A fitted \code{glmMixBayes} or \code{survMixBayes} object.
#' @param ... Not used.
#' @return A numeric matrix (draws x parameters); for \code{as.mcmc()} a
#'   \code{coda::mcmc} object with \code{start}, \code{end} and \code{thin}
#'   attributes taken from the MCMC settings.
#'
#' @examples
#' set.seed(1)
#' X <- cbind(1, rnorm(120)); colnames(X) <- c("(Intercept)", "x")
#' y <- rnorm(120, X %*% c(1, 0.5))
#' fit <- glmMixBayes(X, y, family = "gaussian",
#'                    control = list(iterations = 400, burnin.iterations = 100, seed = 1))
#' draws <- posterior_draws(fit)
#' colnames(draws)
#' if (requireNamespace("coda", quietly = TRUE)) {
#'   coda::effectiveSize(coda::as.mcmc(fit))
#' }
#' @export
posterior_draws <- function(x, ...) {
 UseMethod("posterior_draws", x)
}

# Column-bind a block of draws with bracketed names.
.draw_block <- function(name, draws) {
 if (is.null(draws)) return(NULL)
 if (is.matrix(draws)) {
  cn <- colnames(draws)
  if (is.null(cn)) cn <- seq_len(ncol(draws))
  colnames(draws) <- paste0(name, "[", cn, "]")
  draws
 } else {
  m <- matrix(as.numeric(draws), ncol = 1L)
  colnames(m) <- name
  m
 }
}

.mixbayes_draws <- function(x) {
 est <- .mixbayes_compat(x)$estimates
 blocks <- list(
  .draw_block("coefficients", est$coefficients),
  .draw_block("coefficients2", est$coefficients2),
  .draw_block("m.coefficients", est$m.coefficients),
  .draw_block("theta", est$theta),
  .draw_block("gamma", est$gamma),
  .draw_block("dispersion", est$dispersion),
  .draw_block("dispersion2", est$dispersion2),
  .draw_block("shape", est$shape),
  .draw_block("shape2", est$shape2),
  .draw_block("scale", est$scale),
  .draw_block("scale2", est$scale2),
  .draw_block("lp", x$diagnostics$lp)
 )
 do.call(cbind, Filter(Negate(is.null), blocks))
}

#' @rdname posterior_draws
#' @export
posterior_draws.glmMixBayes <- function(x, ...) .mixbayes_draws(x)

#' @rdname posterior_draws
#' @export
posterior_draws.survMixBayes <- function(x, ...) .mixbayes_draws(x)

# as.mcmc() methods for coda's generic, registered lazily (delayed S3
# registration) so that the package does not need to import coda.
.mixbayes_as_mcmc <- function(x) {
 if (!requireNamespace("coda", quietly = TRUE)) {
  stop("Package 'coda' is required for as.mcmc(); install it with install.packages(\"coda\").",
       call. = FALSE)
 }
 draws <- .mixbayes_draws(x)
 s <- x$diagnostics$settings
 thin <- if (!is.null(s$thin)) as.integer(s$thin) else 1L
 start <- if (!is.null(s$burnin.iterations)) as.integer(s$burnin.iterations) + thin else 1L
 coda::mcmc(draws, start = start, thin = thin)
}

#' @rdname posterior_draws
#' @exportS3Method coda::as.mcmc
as.mcmc.glmMixBayes <- function(x, ...) .mixbayes_as_mcmc(x)

#' @rdname posterior_draws
#' @exportS3Method coda::as.mcmc
as.mcmc.survMixBayes <- function(x, ...) .mixbayes_as_mcmc(x)
