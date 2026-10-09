#' Methods for Bayesian Mixture GLM Fits
#'
#' @description
#' S3 methods for objects returned by \code{glmMixBayes()}, including printing,
#' summarizing fitted models, computing credible intervals and posterior
#' covariance matrices, generating predictions, and pooling regression results
#' across posterior match classifications.
#'
#' @name mixture_bayesglm_methods
#' @keywords internal
NULL

#' Pool Parameter Estimates Across Posterior Draws
#'
#' Generic function for pooling parameter estimates from Bayesian mixture
#' models using posterior draws of the latent component indicators.
#'
#' @param object A fitted Bayesian mixture model object.
#' @param ... Additional arguments passed to methods.
#'
#' @return A pooled model object combining parameter estimates across
#' posterior component-indicator draws.
#'
#' @examples
#' # mi_with() is a generic function for posterior allocation–based pooling.
#' # See ?mi_with.glmMixBayes for a complete example illustrating its use
#' # with Bayesian GLM mixture models.
#'
#' @export
mi_with <- function(object, ...) {
 UseMethod("mi_with", object)
}

#' Print a glmMixBayes Model Object
#'
#' Prints the call, the posterior means of the regression coefficients of
#' component 1 (the correct-match component), the posterior mean of the mixing
#' weight \code{theta} or, when the match probability was modelled with linkage
#' covariates, of the coefficients of the mismatch-indicator model
#' (\code{m.coefficients}), the posterior mean of the dispersion of component 1
#' (Gaussian and Gamma fits), and the numbers of records and of safe matches.
#'
#' @param x An object of class \code{glmMixBayes}.
#' @param digits Minimum number of significant digits to show.
#' @param ... Further arguments (unused).
#' @return The input \code{x}, invisibly.
#'
#' @examples
#' data(lifem)
#'
#' # lifem data preprocessing
#' # For computational efficiency in the example, we work with a subset of the lifem data.
#' lifem <- lifem[order(-(lifem$commf + lifem$comml)), ]
#' lifem_small <- rbind(
#'   head(subset(lifem, hndlnk == 1), 100),
#'   head(subset(lifem, hndlnk == 0), 20)
#' )
#'
#' # priors on the scale of the outcome (age in years) and of the slopes of the
#' # cubic polynomial: the defaults (normal(0, 10) for the intercepts,
#' # normal(0, 5) for the slopes) suit outcomes and covariates of order one
#' adj <- adjMixBayes(
#'   linked.data = lifem_small,
#'   priors = list(
#'     theta = "beta(2, 2)",
#'     intercept1 = "normal(60, 20)",
#'     intercept2 = "normal(60, 20)",
#'     beta1 = "normal(0, 100)",
#'     beta2 = "normal(0, 100)"
#'   )
#' )
#'
#' fit <- plglm(
#'   age_at_death ~ poly(unit_yob, 3, raw = TRUE),
#'   family = "gaussian",
#'   adjustment = adj,
#'   control = list(
#'     iterations = 200,
#'     burnin.iterations = 100,
#'     seed = 123
#'   )
#' )
#'
#' print(fit)
#'
#' @export
#' @method print glmMixBayes
print.glmMixBayes <- function(x, digits = max(3L, getOption("digits") - 3L),...){
  obj <- .mixbayes_compat(x)
  est <- obj$estimates
  n_safe <- if (is.null(obj$diagnostics$n_safe)) 0L else obj$diagnostics$n_safe
  among <- if (n_safe > 0) " among records not flagged as safe matches" else ""
  pm <- function(v) print(format(signif(v, digits)), print.gap = 2, quote = FALSE)

  cat("Call:\n")
  print(obj$call, quote = FALSE, digits = digits)
  cat("\n")

  cat("Coefficients (component 1 = correct-match; posterior means):", sep = "\n")
  pm(colMeans(.mixbayes_block_draws(obj, "coefficients")))
  cat("\n")

  if (isTRUE(obj$use_logistic) && !is.null(est$m.coefficients)) {
   cat(paste0("Mismatch Model Coefficients (logit of the mismatch probability", among,
              "; posterior means):"), sep = "\n")
   pm(colMeans(.mixbayes_block_draws(obj, "m.coefficients")))
   cat("\n")
  } else if (!is.null(est$theta)) {
   cat(paste0("Match probability theta (mixing weight of component 1", among, "; posterior mean): ",
              format(signif(mean(est$theta), digits))), sep = "\n")
  }
  if (!is.null(est$dispersion)) {
   disp_label <- if (identical(obj$family, "gaussian")) "residual variance sigma^2" else "1 / shape"
   cat(paste0("Dispersion (", disp_label, ", component 1; posterior mean): ",
              format(signif(mean(est$dispersion), digits))), sep = "\n")
  }
  if (is.matrix(obj$m_samples)) {
   cat(sprintf("Records: %d (safe matches: %d)\n", ncol(obj$m_samples), as.integer(n_safe)))
  }
  cat("\n")

  invisible(x)
}

#' Summary Method for glmMixBayes Models
#'
#' Posterior summaries of every parameter block of a \code{glmMixBayes} fit.
#' Each table has one row per parameter and the columns \code{"Estimate"} (the
#' posterior mean), \code{"Std. Error"} (the posterior standard deviation) and
#' \code{"2.5 \%"}, \code{"97.5 \%"} (the posterior quantiles bounding the
#' central 95% credible interval), the column names used by
#' \code{\link{glmMixture}()} and \code{stats::confint()}. The tables of the
#' coefficients and the dispersion of component 1 and of the match-probability
#' model (\code{theta}, \code{m.coefficients}, \code{gamma}) add the columns
#' \code{"MCSE"} (the Monte Carlo standard error of the posterior mean: the
#' posterior standard deviation divided by the square root of the effective
#' sample size), \code{"ESS"} (the effective sample size, as in
#' \code{diagnostics$ess}) and \code{"Rhat"} (the single-chain split R-hat,
#' as in \code{diagnostics$rhat}; see \emph{Sampling algorithm} in
#' \code{\link{glmMixBayes}}); for fits that do not store these diagnostics
#' they are computed from the draws, for the tables and for the notes on
#' them, and so are the values of design columns that share a name (e.g.
#' given duplicated names, or unnamed columns in fits saved by earlier
#' versions).
#'
#' @param object An object of class \code{glmMixBayes}.
#' @param ... Not used.
#' @return An object of class \code{"summary.glmMixBayes"}, which is printed
#'   with a custom method. It contains \code{call}, \code{family}, and the
#'   tables of the regression coefficients of component 1
#'   (\code{coefficients}, the correct matches) and of component 2
#'   (\code{coefficients2}, the mismatches); for Gaussian and Gamma fits the
#'   tables of the dispersion of both components (\code{dispersion},
#'   \code{dispersion2}: \eqn{\sigma^2}{sigma^2} or \eqn{1/\phi}{1/phi}); the
#'   table of the coefficients of the mismatch-indicator model
#'   (\code{m.coefficients}, the logit of the mismatch probability, as in
#'   \code{\link{glmMixture}()}); the table of the mixing weight \code{theta}
#'   or, when the match probability was modelled with linkage covariates, of
#'   its logistic regression coefficients \code{gamma} (\code{= -m.coefficients}),
#'   for the original linkage covariates, together with \code{z_center}, the
#'   centre of the linkage covariates of the fit, and
#'   \code{gamma_intercept_prior}, the \code{gamma_intercept} prior used, which
#'   applies to a record whose linkage covariates equal \code{z_center} (both
#'   printed under the \code{gamma} table; see \code{\link{glmMixBayes}});
#'   \code{match.prob}, the posterior probability that each record is a correct
#'   match; \code{match.rate}, a list with the average of \code{match.prob}
#'   over all records (safe matches counted as 1) and over the records that
#'   are not flagged as safe matches (\code{avg}), and the numbers of these
#'   records (\code{n}); \code{use_logistic}; \code{n_safe}, the number of
#'   known correct matches (with safe matches, \code{theta}, \code{gamma} and
#'   \code{m.coefficients} describe the records that are not flagged as safe);
#'   and the diagnostics printed as notes: \code{low_ess} (effective sample
#'   sizes below 100), \code{high_rhat} (split R-hat values above 1.05; the
#'   note says that R-hat is noisy when fewer than 200 draws were stored, their
#'   number being \code{n_draws}),
#'   \code{joint_short} (when the joint moves were off because the burn-in was
#'   too short: the burn-in used and the one that would suffice) and
#'   \code{mixed_labels}. The printed tables show the effective sample size
#'   as a whole number and the split R-hat with three decimals.
#'   Up to postlink 0.1.2 the coefficients of component 2 were stored as
#'   \code{m.coefficients}, the first column of the coefficient tables was
#'   named \code{"Estimates"}, and the dispersion tables had no interval
#'   columns; summary objects saved with these versions are printed with the
#'   component-2 tables under component 2.
#'
#' @examples
#' data(lifem)
#'
#' # lifem data preprocessing
#' # For computational efficiency in the example, we work with a subset of the lifem data.
#' lifem <- lifem[order(-(lifem$commf + lifem$comml)), ]
#' lifem_small <- rbind(
#'   head(subset(lifem, hndlnk == 1), 100),
#'   head(subset(lifem, hndlnk == 0), 20)
#' )
#'
#' # priors on the scale of the outcome (age in years) and of the slopes of the
#' # cubic polynomial: the defaults (normal(0, 10) for the intercepts,
#' # normal(0, 5) for the slopes) suit outcomes and covariates of order one
#' adj <- adjMixBayes(
#'   linked.data = lifem_small,
#'   priors = list(
#'     theta = "beta(2, 2)",
#'     intercept1 = "normal(60, 20)",
#'     intercept2 = "normal(60, 20)",
#'     beta1 = "normal(0, 100)",
#'     beta2 = "normal(0, 100)"
#'   )
#' )
#'
#' fit <- plglm(
#'   age_at_death ~ poly(unit_yob, 3, raw = TRUE),
#'   family = "gaussian",
#'   adjustment = adj,
#'   control = list(
#'     iterations = 200,
#'     burnin.iterations = 100,
#'     seed = 123
#'   )
#' )
#'
#' summary(fit)
#'
#' @export
summary.glmMixBayes <- function(object, ...) {
 object <- .mixbayes_compat(object)
 est <- object$estimates
 dg <- object$diagnostics
 # effective sample sizes and split R-hat (computed from the draws for fits
 # saved without them, so that the tables and the notes agree)
 dg <- .mixbayes_mcmc_diagnostics(object, dg)
 # Monte Carlo columns (MCSE, ESS, Rhat) of the tables of component 1 and of
 # the match-probability model, named by their block in posterior_draws()
 mc <- function(block) list(block = block, ess = dg$ess, rhat = dg$rhat)

 # posterior mean, SD and central 95% interval of every parameter
 out <- list(
  call          = object$call,
  family        = object$family,
  coefficients  = .posterior_table(est$coefficients, "(Intercept)", mcmc = mc("coefficients")),
  coefficients2 = .posterior_table(est$coefficients2, "(Intercept)")
 )
 if (!is.null(est$dispersion)) {
  out$dispersion <- .posterior_table(est$dispersion, "dispersion", mcmc = mc("dispersion"))
 }
 if (!is.null(est$dispersion2)) out$dispersion2 <- .posterior_table(est$dispersion2, "dispersion2")

 # Mismatch-indicator model (as the m.coefficients of glmMixture()), and the
 # match-probability model it is derived from: the constant mixing weight
 # theta, or the logistic regression coefficients gamma when covariates were
 # supplied via m.formula
 if (!is.null(est$m.coefficients)) {
  out$m.coefficients <- .posterior_table(est$m.coefficients, "(Intercept)", mcmc = mc("m.coefficients"))
 }
 out$use_logistic <- isTRUE(object$use_logistic)
 if (out$use_logistic && !is.null(est$gamma)) {
  out$gamma <- .posterior_table(est$gamma, "(Intercept)", mcmc = mc("gamma"))
  # the gamma_intercept prior applies at the centre of the linkage covariates
  out$z_center <- object$z_center
  out$gamma_intercept_prior <- .gamma_intercept_prior(object$priors)
 } else if (!is.null(est$theta)) {
  out$theta <- .posterior_table(est$theta, "theta", mcmc = mc("theta"))
 }
 out$n_safe <- if (is.null(dg$n_safe)) 0L else dg$n_safe
 out$match.prob <- object$match.prob
 if (!is.null(object$match.prob)) out$match.rate <- .match_rate_summary(object$match.prob, out$n_safe)
 out$low_ess <- .low_ess(dg$ess)
 out$high_rhat <- .high_rhat(dg$rhat)
 out$n_draws <- NROW(est$coefficients)
 out$joint_short <- .joint_moves_short(dg)
 out$mixed_labels <- .mixed_labels(dg)

 class(out) <- "summary.glmMixBayes"
 out
}

#' @noRd
#' @export
print.summary.glmMixBayes <- function(x, digits = max(3L, getOption("digits") - 3L),
                                      signif.stars = getOption("show.signif.stars"),...){
  out <- x
  # summaries saved by postlink <= 0.1.2 (and by the development versions
  # before the names were aligned with adjMixture()) stored the component-2
  # tables as m.coefficients and m.dispersion and had no mismatch-model table
  if (is.null(x$coefficients2) && !is.null(x$m.coefficients)) {
   x$coefficients2 <- x$m.coefficients
   x$dispersion2 <- x$m.dispersion
   x$m.coefficients <- NULL
   x$m.dispersion <- NULL
  }

  cat("Call:", sep="\n")
  print(x$call,quote=F)
  cat(" ", sep="\n")
  cat("Family:", x$family, " ", sep="\n")

  # significant digits per column, so that small values do not print as 0
  # (the effective sample size as a whole number, the split R-hat with three
  # decimals)
  print_mat <- function(mat) .print_posterior_table(mat, digits)

  disp_label <- if (identical(x$family, "gaussian")) "Dispersion (residual variance sigma^2):" else
   "Dispersion (1 / shape):"

  cat("(Component 1 = Correct-match):", sep = "\n")

  cat("Outcome Model Coefficients:", sep="\n")
  print_mat(x$coefficients)
  cat("\n")

  if (x$family %in% c("gamma", "gaussian")){
   cat(disp_label, sep="\n")
   print_mat(x$dispersion)
    cat("\n")
  }

  cat("(Component 2 = Incorrect-match):", sep = "\n")

  cat("Outcome Model Coefficients:", sep="\n")
  print_mat(x$coefficients2)
  cat("\n")

  if (x$family %in% c("gamma", "gaussian")){
   cat(disp_label, sep="\n")
   print_mat(x$dispersion2)
   cat("\n")
  }

  heads <- .match_headings(x$n_safe > 0)
  if (!is.null(x$m.coefficients)) {
   cat(paste0(heads[["m.coefficients"]], ":"), sep = "\n")
   print_mat(x$m.coefficients)
   cat("\n")
  }

  if (isTRUE(x$use_logistic) && !is.null(x$gamma)) {
   cat(paste0(heads[["gamma"]], ":"), sep = "\n")
   print_mat(x$gamma)
   .print_z_center_note(x$z_center, x$gamma_intercept_prior, digits)
   cat("\n")
  } else if (!is.null(x$theta)) {
   cat(paste0(heads[["theta"]], ":"), sep = "\n")
   print_mat(x$theta)
   cat("\n")
  }
  if (!is.null(x$match.rate)) {
   .print_match_rate(x$match.rate, digits)
   cat("\n")
  }
  cat("(Estimate: posterior mean; Std. Error: posterior standard deviation; 2.5 % and 97.5 %: posterior quantiles)\n")
  if ("MCSE" %in% colnames(x$coefficients)) {
   cat("(MCSE: Monte Carlo standard error of the estimate; ESS: effective sample size; Rhat: split R-hat of the single chain)\n")
  }
  .print_low_ess(x$low_ess)
  .print_high_rhat(x$high_rhat, x$n_draws)
  .print_joint_short(x$joint_short)
  .print_mixed_labels(x$mixed_labels, "glmMixBayes")

  invisible(out)
}

#' Posterior Covariance Matrix for glmMixBayes Coefficients
#'
#' @param object A \code{glmMixBayes} model object.
#' @param block Which parameter block: \code{"coefficients"} (component 1, the
#'   correct-match component; default), \code{"coefficients2"} (component 2),
#'   \code{"m.coefficients"} (the mismatch-indicator model, as in
#'   \code{\link{glmMixture}()}), \code{"theta"}, \code{"gamma"},
#'   \code{"dispersion"} or \code{"dispersion2"}; see
#'   \code{\link{confint.glmMixBayes}}.
#' @param ... Not used; other arguments raise an error.
#' @return Posterior covariance matrix of the parameters of the selected block
#'   (by default the regression coefficients of component 1), i.e. the
#'   covariance of their posterior draws.
#'
#' @examples
#' data(lifem)
#'
#' # lifem data preprocessing
#' # For computational efficiency in the example, we work with a subset of the lifem data.
#' lifem <- lifem[order(-(lifem$commf + lifem$comml)), ]
#' lifem_small <- rbind(
#'   head(subset(lifem, hndlnk == 1), 100),
#'   head(subset(lifem, hndlnk == 0), 20)
#' )
#'
#' # priors on the scale of the outcome (age in years) and of the slopes of the
#' # cubic polynomial: the defaults (normal(0, 10) for the intercepts,
#' # normal(0, 5) for the slopes) suit outcomes and covariates of order one
#' adj <- adjMixBayes(
#'   linked.data = lifem_small,
#'   priors = list(
#'     theta = "beta(2, 2)",
#'     intercept1 = "normal(60, 20)",
#'     intercept2 = "normal(60, 20)",
#'     beta1 = "normal(0, 100)",
#'     beta2 = "normal(0, 100)"
#'   )
#' )
#'
#' fit <- plglm(
#'   age_at_death ~ poly(unit_yob, 3, raw = TRUE),
#'   family = "gaussian",
#'   adjustment = adj,
#'   control = list(
#'     iterations = 200,
#'     burnin.iterations = 100,
#'     seed = 123
#'   )
#' )
#'
#' vcov(fit)
#'
#' # posterior covariance of the mixing weight
#' vcov(fit, block = "theta")
#'
#' @export
vcov.glmMixBayes <- function(object,
                             block = c("coefficients", "coefficients2", "m.coefficients", "theta",
                                       "gamma", "dispersion", "dispersion2"),
                             ...) {
 .check_unused_args(list(...), "vcov.glmMixBayes", "Select a parameter block with `block`.")
 block <- match.arg(block)
 stats::cov(.mixbayes_block_draws(.mixbayes_compat(object), block))
}

#' Credible Intervals for Parameters from a glmMixBayes Fit
#'
#' Computes posterior credible intervals for the parameters of a fitted
#' \code{glmMixBayes} model. By default, intervals are returned for all
#' regression coefficients in component 1 (the correct-match component) of the
#' mixture model; \code{block} selects another parameter block and a subset of
#' its parameters can be selected using \code{parm}.
#'
#' @param object A \code{glmMixBayes} model object.
#' @param parm Optional. Names or numeric indices of the parameters within the
#'   selected \code{block} for which credible intervals should be returned. If
#'   \code{NULL}, intervals are returned for all parameters of the block. When
#'   \code{block} is not given, a single block name (any of the names accepted
#'   by \code{block}, e.g. \code{parm = "theta"}, \code{"gamma"} or
#'   \code{"m.coefficients"}) selects that block (as in
#'   \code{\link{confint.survMixBayes}}), unless it is also the name of a
#'   coefficient of component 1.
#' @param level Probability level for the credible intervals. Defaults to
#'   \code{0.95}.
#' @param block Which parameter block to summarize: \code{"coefficients"}
#'   (component 1, the correct matches; default), \code{"coefficients2"}
#'   (component 2, the mismatches), \code{"m.coefficients"} (the
#'   mismatch-indicator model on the logit scale of the mismatch probability, as
#'   the \code{m.coefficients} of \code{\link{glmMixture}()}), \code{"theta"}
#'   (the mixing weight, i.e. the probability of a correct match, when
#'   \code{m.formula} has no covariates), \code{"gamma"} (the logistic
#'   regression coefficients of the match probability, \code{-m.coefficients},
#'   when \code{m.formula} has covariates), \code{"dispersion"} or
#'   \code{"dispersion2"} (Gaussian and Gamma fits). Up to postlink 0.1.2,
#'   \code{"m.coefficients"} selected the coefficients of component 2.
#' @param ... Not used; other arguments raise an error.
#' @return
#' A matrix with one row per parameter and two columns giving the lower and
#' upper credible interval bounds (the posterior quantiles), labelled as by
#' \code{stats::confint()} (\code{"2.5 \%"} and \code{"97.5 \%"} for
#' \code{level = 0.95}). Row names correspond to parameter names.
#'
#' @examples
#' data(lifem)
#'
#' # lifem data preprocessing
#' # For computational efficiency in the example, we work with a subset of the lifem data.
#' lifem <- lifem[order(-(lifem$commf + lifem$comml)), ]
#' lifem_small <- rbind(
#'   head(subset(lifem, hndlnk == 1), 100),
#'   head(subset(lifem, hndlnk == 0), 20)
#' )
#'
#' # priors on the scale of the outcome (age in years) and of the slopes of the
#' # cubic polynomial: the defaults (normal(0, 10) for the intercepts,
#' # normal(0, 5) for the slopes) suit outcomes and covariates of order one
#' adj <- adjMixBayes(
#'   linked.data = lifem_small,
#'   priors = list(
#'     theta = "beta(2, 2)",
#'     intercept1 = "normal(60, 20)",
#'     intercept2 = "normal(60, 20)",
#'     beta1 = "normal(0, 100)",
#'     beta2 = "normal(0, 100)"
#'   )
#' )
#'
#' fit <- plglm(
#'   age_at_death ~ poly(unit_yob, 3, raw = TRUE),
#'   family = "gaussian",
#'   adjustment = adj,
#'   control = list(
#'     iterations = 200,
#'     burnin.iterations = 100,
#'     seed = 123
#'   )
#' )
#'
#' confint(fit)
#'
#' # credible interval of the probability of a correct match
#' confint(fit, block = "theta")
#'
#' # the same on the logit scale of the mismatch probability, as the
#' # m.coefficients of adjMixture() fits
#' confint(fit, block = "m.coefficients")
#'
#' @export
#' @method confint glmMixBayes
confint.glmMixBayes <- function(object, parm = NULL, level = 0.95,
                                block = c("coefficients", "coefficients2", "m.coefficients", "theta",
                                          "gamma", "dispersion", "dispersion2"),
                                ...) {
 .check_unused_args(list(...), "confint.glmMixBayes", "Select a parameter block with `block`.")
 object <- .mixbayes_compat(object)
 # confint(fit, parm = "gamma") selects a block, as in confint.survMixBayes()
 # (coefficient names of component 1 take precedence over block names)
 if (missing(block) && is.character(parm) && length(parm) == 1L &&
     parm %in% .mixbayes_blocks(object) &&
     !parm %in% colnames(object$estimates$coefficients)) {
  block <- parm
  parm <- NULL
 }
 block <- match.arg(block)
 if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1) {
  stop("`level` must be a single number strictly between 0 and 1.", call. = FALSE)
 }
 draws <- .mixbayes_block_draws(object, block)
 alpha <- 1 - level
 .select_parm(.draw_quantiles(draws, c(alpha / 2, 1 - alpha / 2)), parm)
}


#' Predictions from a glmMixBayes Model
#'
#' Predictions of the correct-match component (component 1). New data are
#' given as a data frame in \code{newdata}, as for \code{\link{glmMixture}()}
#' fits (this needs the model terms stored by \code{\link{plglm}()}), or as a
#' model matrix in \code{newx} (for fits from \code{glmMixBayes()}).
#'
#' @param object A \code{glmMixBayes} model object.
#' @param newdata Optional data frame in which to look for the variables of the
#'   outcome model; its model matrix is built from the terms of a fit from
#'   \code{plglm()} (stored with \code{model = TRUE}, the default). A numeric
#'   matrix is taken to be a model matrix, as \code{newx}, so that calls of
#'   postlink 0.1.2 such as \code{predict(fit, X, "response")} keep working.
#' @param type Either \code{"link"} or \code{"response"}, indicating the scale of predictions.
#' @param se.fit Logical; if \code{TRUE}, also return posterior SD of predictions.
#' @param interval Either \code{"none"} or \code{"credible"}, indicating whether to compute a credible interval.
#' @param level Probability level for the credible interval (default 0.95).
#' @param na.action Function determining what to do with missing values in
#'   \code{newdata} (default \code{stats::na.pass}, which gives \code{NA}
#'   predictions).
#' @param newx Optional numeric matrix of new observations (n_new x K) with
#'   columns aligned to the design matrix \code{X} used for fitting (matched by
#'   name when the coefficient names are non-empty and unique and each names
#'   exactly one column, otherwise by position, with a warning when a column
#'   bears the name of a coefficient at another position). It follows the
#'   arguments of \code{predict.glmMixture()}, so give it by name.
#' @param ... Not used; other arguments (e.g. a misspelled \code{newdata})
#'   raise an error.
#'
#' @details When neither \code{newdata} nor \code{newx} is given, predictions
#'   are made for the analysed records of the fit (in the order of
#'   \code{match.prob}): from the design matrix stored by
#'   \code{plglm(..., x = TRUE)}, or else from the stored model frame. Records
#'   dropped at fit time because their linkage covariates are missing are
#'   left out.
#'
#' @return If \code{se.fit = FALSE} and \code{interval = "none"}, a numeric vector of predicted values.
#'   Otherwise, a matrix with columns for the fit, (optional) \code{se.fit}, and (optional)
#'   credible interval bounds labelled as by \code{stats::confint()}
#'   (\code{"2.5 \%"}, \code{"97.5 \%"}). The fitted value is the prediction at the
#'   posterior mean of the component-1 coefficients (for \code{type = "response"}
#'   the inverse link of \eqn{x^\top E[\beta]}{x' E[beta]}), whereas \code{se.fit}
#'   and the credible interval summarise the posterior draws of the prediction.
#'   Predictions are named after the rows of the model matrix (the row names
#'   of \code{newdata}, the row names of \code{newx} if it has any, or the
#'   analysed records), as by \code{\link[stats]{predict.glm}()}; rows with
#'   missing covariates give \code{NA}.
#'
#' @examples
#' data(lifem)
#'
#' # lifem data preprocessing
#' # For computational efficiency in the example, we work with a subset of the lifem data.
#' lifem <- lifem[order(-(lifem$commf + lifem$comml)), ]
#' lifem_small <- rbind(
#'   head(subset(lifem, hndlnk == 1), 100),
#'   head(subset(lifem, hndlnk == 0), 20)
#' )
#'
#' # priors on the scale of the outcome (age in years) and of the slopes of the
#' # cubic polynomial: the defaults (normal(0, 10) for the intercepts,
#' # normal(0, 5) for the slopes) suit outcomes and covariates of order one
#' adj <- adjMixBayes(
#'   linked.data = lifem_small,
#'   priors = list(
#'     theta = "beta(2, 2)",
#'     intercept1 = "normal(60, 20)",
#'     intercept2 = "normal(60, 20)",
#'     beta1 = "normal(0, 100)",
#'     beta2 = "normal(0, 100)"
#'   )
#' )
#'
#' fit <- plglm(
#'   age_at_death ~ poly(unit_yob, 3, raw = TRUE),
#'   family = "gaussian",
#'   adjustment = adj,
#'   control = list(
#'     iterations = 200,
#'     burnin.iterations = 100,
#'     seed = 123
#'   )
#' )
#'
#' # new data as a data frame (variables of the outcome formula) ...
#' predict(fit, newdata = data.frame(unit_yob = c(0.2, 0.5, 0.8)), type = "response")
#'
#' # ... or as a model matrix
#' newx <- cbind(1, poly(c(0.2, 0.5, 0.8), 3, raw = TRUE))
#' predict(fit, newx = newx, type = "response")
#'
#' @export
predict.glmMixBayes <- function(object, newdata = NULL,
                                type = c("link", "response"),
                                se.fit = FALSE,
                                interval = c("none", "credible"),
                                level = 0.95,
                                na.action = stats::na.pass,
                                newx = NULL,
                                ...) {
  .check_unused_args(list(...), "predict.glmMixBayes",
                     "Pass new data as a data frame in `newdata` or as a model matrix in `newx`.")
  type <- match.arg(type)
  interval <- match.arg(interval)

  if (!is.null(newdata) && !is.null(newx)) {
    stop("Supply either `newdata` or `newx`, not both.", call. = FALSE)
  }
  if (!is.null(newdata)) {
    if (is.matrix(newdata) && is.numeric(newdata)) {
      newx <- newdata                   # a model matrix, as `newx`
    } else if (is.data.frame(newdata)) {
      newx <- .mixbayes_newdata_matrix(object, newdata, na.action)
    } else {
      stop("`newdata` must be a data frame (or a numeric model matrix).", call. = FALSE)
    }
  }
  # no new data: the records of the fit, from the stored design matrix or
  # model frame (restricted to the analysed records)
  if (is.null(newx)) newx <- .mixbayes_fit_design(object)
  if (is.null(newx)) {
    stop("No new data were given and the fit stores neither its design matrix nor its model frame: ",
         "pass the model matrix as `newx` (or refit with plglm(..., x = TRUE)).", call. = FALSE)
  }
  if (!is.matrix(newx) || !is.numeric(newx)) {
    stop("`newx` must be a numeric matrix with columns aligned to the fitted design matrix.",
         call. = FALSE)
  }
  newx <- .align_newx(newx, object$estimates$coefficients, "newx")
  if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1) {
    stop("`level` must be a single number strictly between 0 and 1.", call. = FALSE)
  }

  mean_coef <- apply(object$estimates$coefficients, 2, mean)

  # point predictions, named after the rows of the model matrix (the rows of
  # newdata or the analysed records), as by stats::predict.glm()
  if (type == "link") {
    predictions <- as.vector(newx %*% mean_coef)
  } else {
    eta <- as.vector(newx %*% mean_coef)
    predictions <- switch(object$family,
      gaussian = eta,
      gamma    = exp(eta),
      poisson  = exp(eta),
      binomial = stats::plogis(eta),
      stop("Unknown family in object.")
    )
  }
  names(predictions) <- rownames(newx)

  if (!se.fit && interval == "none") {
    return(predictions)
  }

  # posterior predictive draws: rows are new data points, columns are draws
  if (type == "link") {
    all_predictions <- newx %*% t(object$estimates$coefficients)
  } else {
    eta_all <- newx %*% t(object$estimates$coefficients)
    all_predictions <- switch(object$family,
      gaussian = eta_all,
      gamma    = exp(eta_all),
      poisson  = exp(eta_all),
      binomial = stats::plogis(eta_all),
      stop("Unknown family in object.")
    )
  }

  se.predictions <- apply(all_predictions, 1, stats::sd)

  if (interval == "none") {
    vals <- cbind(fit = predictions, se.fit = se.predictions)
    colnames(vals) <- c("fit", "se.fit")
    return(vals)
  }

  alpha <- 1 - level
  ci <- t(apply(all_predictions, 1, stats::quantile, probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE))
  lower <- ci[, 1]
  upper <- ci[, 2]

  if (!se.fit && interval == "credible") {
    vals <- cbind(fit = predictions, lower = lower, upper = upper)
    colnames(vals) <- c("fit", .ci_labels(c(alpha / 2, 1 - alpha / 2)))
    return(vals)
  }

  vals <- cbind(fit = predictions, se.fit = se.predictions, lower = lower, upper = upper)
  colnames(vals) <- c("fit", "se.fit", .ci_labels(c(alpha / 2, 1 - alpha / 2)))
  vals
}


#' Pooling Regression Fits Across Posterior Draws of Correct-Match Classifications
#'
#' @description
#' Use posterior draws of the latent match indicators from \code{glmMixBayes()}
#' to repeatedly identify which records are treated as correct matches, refit the
#' requested regression model on those records, and pool the resulting estimates.
#'
#' Each retained posterior draw defines one subset of records classified as
#' correct matches. The function fits the specified \code{lm()} or \code{glm()}
#' model to that subset, extracts the estimated coefficients and their covariance
#' matrix, and combines the results across draws using multiple-imputation
#' pooling rules. At least two draws must be usable, since the between-draw
#' variance cannot be estimated from one; fewer raise an error.
#'
#' @param object A \code{glmMixBayes} model object containing posterior draws of
#'   the latent match indicators.
#' @param data A data.frame with the records used in the model. For fits from
#'   \code{plglm()} / \code{plsurvreg()} this is normally the \code{linked.data}
#'   of the adjustment object: its rows are matched to the analysed records by
#'   row name, so records dropped by \code{subset} or \code{na.action} are
#'   handled correctly. Otherwise \code{data} must contain exactly the analysed
#'   records, in the order used by the model.
#' @param formula Model formula for refitting on each draw. If omitted, the
#'   outcome formula of a \code{plglm()} fit is used (taken from the stored
#'   model frame, or from a formula written directly in the call).
#'   Data-dependent terms such as \code{poly()}, \code{scale()} or
#'   \code{splines::ns()} keep in every refit the basis computed once from
#'   all analysed records in \code{data} (for the stored outcome formula, the
#'   basis of the fit, computed from the records of its model frame), so that
#'   the pooled coefficients refer to a single basis, as those of a fit to all
#'   records do, rather than to a basis recomputed from the records of each
#'   draw. Only the terms whose prediction call R records (the
#'   \code{"predvars"} of the terms, as used by \code{predict()}) are fixed:
#'   top-level calls such as \code{poly()}, \code{scale()},
#'   \code{splines::ns()} or \code{splines::bs()}. A data-dependent call
#'   nested in another one, such as \code{I(scale(x)^2)} or
#'   \code{log(scale(x) + 3)}, is evaluated from the records of each draw, as
#'   \code{predict()} would evaluate it; compute such variables in
#'   \code{data} beforehand.
#' @param family A \code{stats::family()} object (or a family function, or its
#'   name) for the refitted model. If not supplied, the family and link of the
#'   mixture components are used: \code{gaussian()}, \code{poisson()},
#'   \code{binomial()} or \code{Gamma(link = "log")}. The name \code{"gamma"}
#'   (in any case, e.g. \code{"Gamma"}) gives \code{Gamma(link = "log")}, as in
#'   \code{\link{glmMixture}()}. A family with a different link (e.g. the
#'   family function \code{Gamma}, whose default link is the inverse) puts the
#'   pooled coefficients on a different scale from
#'   \code{object$estimates$coefficients}. Gaussian models with the identity
#'   link are refitted with \code{lm()}, all others with \code{glm()}.
#' @param min_n Minimum number of records required to fit the model for a given
#'   posterior draw. The default is \code{p + 1}, where \code{p} is the number
#'   of columns in the model matrix.
#' @param quietly If \code{TRUE}, draws that lead to fitting errors are skipped
#'   without printing the full error message. A warning reports how many draws
#'   were not used.
#' @param ... Additional arguments passed through (currently unused).
#' @return An object of class \code{c("mi_link_pool_glm", "mi_link_pool")}
#'   containing pooled coefficient estimates (\code{coef}), their covariance
#'   matrix (\code{vcov}), standard errors (\code{se}), 95% confidence
#'   intervals (\code{ci95}), the within- and between-draw covariances
#'   (\code{Ubar}, \code{B}), the fraction of missing information
#'   (\code{lambda}), the degrees of freedom (\code{df}), the number and indices
#'   of the draws used (\code{m}, \code{kept_draws}), and the refitted model:
#'   \code{refit} (\code{"lm"} or \code{"glm"}), \code{family}, \code{link} and
#'   the \code{call} of the mixture fit.
#'
#' @examples
#' data(lifem)
#'
#' # lifem data preprocessing
#' # For computational efficiency in the example, we work with a subset of the lifem data.
#' lifem <- lifem[order(-(lifem$commf + lifem$comml)), ]
#' lifem_small <- rbind(
#'   head(subset(lifem, hndlnk == 1), 100),
#'   head(subset(lifem, hndlnk == 0), 20)
#' )
#'
#' # priors on the scale of the outcome (age in years) and of the slopes of the
#' # cubic polynomial: the defaults (normal(0, 10) for the intercepts,
#' # normal(0, 5) for the slopes) suit outcomes and covariates of order one
#' adj <- adjMixBayes(
#'   linked.data = lifem_small,
#'   priors = list(
#'     theta = "beta(2, 2)",
#'     intercept1 = "normal(60, 20)",
#'     intercept2 = "normal(60, 20)",
#'     beta1 = "normal(0, 100)",
#'     beta2 = "normal(0, 100)"
#'   )
#' )
#'
#' fit <- plglm(
#'   age_at_death ~ poly(unit_yob, 3, raw = TRUE),
#'   family = "gaussian",
#'   adjustment = adj,
#'   control = list(
#'     iterations = 200,
#'     burnin.iterations = 100,
#'     seed = 123
#'   )
#' )
#'
#' pooled_fit <- mi_with(
#'   object = fit,
#'   data = lifem_small,
#'   formula = age_at_death ~ poly(unit_yob, 3, raw = TRUE),
#'   family = gaussian()
#' )
#'
#' print(pooled_fit)
#'
#' @export
#' @importFrom stats coef vcov cov lm glm qt model.frame model.matrix gaussian poisson binomial Gamma
mi_with.glmMixBayes <- function(object, data, formula,
                            family = NULL, min_n = NULL, quietly = TRUE, ...) {

  # Resolve formula: the outcome model of the fit when not provided
  if (missing(formula) || is.null(formula)) {
   formula <- .mixbayes_mi_formula(object, parent.frame())
  }

  # Ensure we ended up with a valid formula
  if (!inherits(formula, "formula")) {
   stop("`formula` must be a valid formula object.", call. = FALSE)
  }

  # Rows of `data` aligned with the records used in the fit
  data <- .mixbayes_mi_data(object, data)

  if (is.null(family)) {
    family <- switch(object$family,
      gaussian = stats::gaussian(),
      poisson  = stats::poisson(),
      binomial = stats::binomial(),
      gamma    = stats::Gamma(link = "log"),
      stats::gaussian()
    )
  }
  if (is.character(family)) {
    if (length(family) != 1L || is.na(family)) {
      stop("`family` must be a single family name, a family function or a family object.", call. = FALSE)
    }
    if (tolower(family) == "gamma") {
      # the name used by these fits; as in glmMixture(), the log link of the
      # mixture components (base::gamma() is not a family)
      family <- stats::Gamma(link = "log")
    } else {
      family_fun <- tryCatch(get(family, mode = "function", envir = parent.frame()),
                             error = function(e) NULL)
      if (is.null(family_fun)) stop("Unknown family '", family, "'.", call. = FALSE)
      family <- family_fun
    }
  }
  if (is.function(family)) family <- family()
  if (!inherits(family, "family")) {
    stop("`family` must be a family object such as gaussian(), a family function or its name.",
         call. = FALSE)
  }
  # lm() for Gaussian models with the identity link, glm() otherwise
  refit <- if (identical(family$family, "gaussian") && identical(family$link, "identity")) "lm" else "glm"

  # Extract the matrix of component allocations (S x N)
  z_samples <- object$m_samples

  collapse_z_g <- function(z_samples, g) {
    if (is.list(z_samples)) {
      parts <- lapply(z_samples, function(Z) {
        if (is.null(dim(Z))) {
          stop("Each list element of `z_samples` must be a matrix-like object.")
        }
        Z[g, ]
      })
      as.integer(do.call(c, parts))
    } else {
      if (is.null(dim(z_samples)) || length(dim(z_samples)) != 2L) {
        stop("`z_samples` must be an S x N matrix, or a list of such matrices.")
      }
      as.integer(z_samples[g, ])
    }
  }

  # Terms of the refitted model, resolved once from all analysed records: their
  # "predvars" keep the bases of data-dependent terms (poly(), scale(),
  # splines::ns(), ...) fixed in every refit, so that the pooled coefficients
  # refer to one basis (the stored terms of the fit when no formula is given)
  mf <- stats::model.frame(formula, data = data)
  tt <- attr(mf, "terms")

  fit_once <- function(df) {
    if (refit == "lm") {
      fit <- stats::lm(tt, data = df)
    } else {
      fit <- stats::glm(tt, data = df, family = family)
    }
    list(coef = stats::coef(fit), vcov = stats::vcov(fit))
  }

  # Determine p and default min_n
  Xtmp <- stats::model.matrix(tt, mf)
  p <- ncol(Xtmp)
  if (is.null(min_n)) min_n <- p + 1L

  # Number of posterior draws S
  S <- if (is.list(z_samples)) nrow(z_samples[[1]]) else nrow(z_samples)
  if (length(S) == 0L || is.null(S) || !is.finite(S)) {
    stop("Unable to infer the number of posterior draws (S).")
  }

  coefs_list <- list()
  vcovs_list <- list()
  kept <- logical(S)
  n_small <- 0L
  first_error <- NULL

  for (s in seq_len(S)) {
    zs <- collapse_z_g(z_samples, s)
    if (!all(zs %in% c(1L, 2L))) {
      stop("`m_samples` must contain only 1/2 indicators.")
    }

    if (sum(zs == 1L) < min_n) {
      n_small <- n_small + 1L
      next
    }
    dat_s <- data[which(zs == 1L), , drop = FALSE]

    res <- try(fit_once(dat_s), silent = quietly)
    if (!inherits(res, "try-error")) {
      coefs_list[[length(coefs_list) + 1L]] <- res$coef
      vcovs_list[[length(vcovs_list) + 1L]] <- res$vcov
      kept[s] <- TRUE
    } else if (is.null(first_error)) {
      first_error <- conditionMessage(attr(res, "condition"))
    }
  }

  m <- length(coefs_list)
  if (m == 0L) .mi_no_valid_error(S, n_small, min_n, first_error)
  .mi_check_m(m, S, n_small)
  .mi_skipped_warning(m, S)

  parnames <- names(coefs_list[[1L]])
  coefs_mat <- do.call(rbind, lapply(coefs_list, function(b) b[parnames]))

  Ubar <- Reduce("+", vcovs_list) / m
  B <- stats::cov(coefs_mat)
  Tmat <- Ubar + (1 + 1 / m) * B

  qbar <- colMeans(coefs_mat)
  se <- sqrt(diag(Tmat))
  lambda <- diag((1 + 1 / m) * B) / diag(Tmat)

  r <- diag((1 + 1 / m) * B) / diag(Ubar)
  dfold <- (m - 1) * (1 + 1 / r)^2
  df <- pmax(dfold, 3)

  tcrit <- stats::qt(0.975, df = df)
  lwr <- qbar - tcrit * se
  upr <- qbar + tcrit * se
  ci95 <- cbind(lwr = lwr, upr = upr)
  rownames(ci95) <- names(qbar)

  out <- list(
    m          = m,
    coef       = qbar,
    vcov       = Tmat,
    se         = se,
    ci95       = ci95,
    Ubar       = Ubar,
    B          = B,
    lambda     = lambda,
    df         = df,
    kept_draws = which(kept),
    refit      = refit,
    family     = family$family,
    link       = family$link,
    call       = object$call
  )
  class(out) <- c("mi_link_pool_glm", "mi_link_pool")
  out
}

#' Print Pooled Regression Results
#'
#' @param x An object of class \code{mi_link_pool_glm}, typically returned by
#'   \code{mi_with()} for a \code{glmMixBayes} fit.
#' @param digits the number of significant digits to print.
#' @param ... further arguments (unused).
#' @return The input \code{x}, invisibly.
#'
#' @examples
#' data(lifem)
#'
#' # lifem data preprocessing
#' # For computational efficiency in the example, we work with a subset of the lifem data.
#' lifem <- lifem[order(-(lifem$commf + lifem$comml)), ]
#' lifem_small <- rbind(
#'   head(subset(lifem, hndlnk == 1), 100),
#'   head(subset(lifem, hndlnk == 0), 20)
#' )
#'
#' # priors on the scale of the outcome (age in years) and of the slopes of the
#' # cubic polynomial: the defaults (normal(0, 10) for the intercepts,
#' # normal(0, 5) for the slopes) suit outcomes and covariates of order one
#' adj <- adjMixBayes(
#'   linked.data = lifem_small,
#'   priors = list(
#'     theta = "beta(2, 2)",
#'     intercept1 = "normal(60, 20)",
#'     intercept2 = "normal(60, 20)",
#'     beta1 = "normal(0, 100)",
#'     beta2 = "normal(0, 100)"
#'   )
#' )
#'
#' fit <- plglm(
#'   age_at_death ~ poly(unit_yob, 3, raw = TRUE),
#'   family = "gaussian",
#'   adjustment = adj,
#'   control = list(
#'     iterations = 200,
#'     burnin.iterations = 100,
#'     seed = 123
#'   )
#' )
#'
#' pooled_fit <- mi_with(
#'   object = fit,
#'   data = lifem_small,
#'   formula = age_at_death ~ poly(unit_yob, 3, raw = TRUE),
#'   family = gaussian()
#' )
#'
#' print(pooled_fit, digits = 4)
#'
#' @export
print.mi_link_pool_glm <- function(x, digits = max(3L, getOption("digits") - 2L), ...) {
  cat("Pooled regression results across posterior match classifications:\n")
  cat("  Retained imputations (m):", x$m, "\n")
  if (!is.null(x$refit)) {
    cat("  Refit model: ", x$refit, " (family ", x$family, ", link ", x$link, ")\n", sep = "")
  }
  cat("\n")

  stats::printCoefmat(.mi_pool_table(x), digits = digits, has.Pvalue = FALSE, ...)
  invisible(x)
}
