#' Methods for Bayesian Mixture Survival Regression Fits
#'
#' @description
#' S3 methods for objects returned by \code{survregMixBayes()}, including
#' printing, summarizing fitted models, computing credible intervals and
#' posterior covariance matrices, generating predictions, and pooling Cox
#' regression fits across posterior match classifications.
#'
#' @name mixture_bayessurvreg_methods
#' @keywords internal
NULL

#' Print a survMixBayes Model Object
#'
#' Prints the model call, the posterior means of the regression coefficients
#' of the first mixture component of the fitted survival model, the posterior
#' mean of the mixing weight \code{theta} or, when the match probability was
#' modelled with linkage covariates, of the coefficients of the
#' mismatch-indicator model (\code{m.coefficients}), the posterior mean of the
#' shape of component 1, and the numbers of records and of safe matches. In
#' this package, component 1 is interpreted as the correct-match component and
#' component 2 as the incorrect-match component. For Weibull fits whose design
#' matrix has an intercept column, the intercept is the identified
#' \code{(Intercept) + log(scale)} (see \code{\link{survregMixBayes}}).
#'
#' @param x An object of class \code{survMixBayes}.
#' @param digits Minimum number of significant digits to show.
#' @param ... Further arguments (unused).
#'
#' @return The input \code{x}, invisibly.
#'
#' @examples
#' set.seed(301)
#' n <- 150
#' trt <- rbinom(n, 1, 0.5)
#'
#' # Simulate Weibull AFT data
#' true_time <- rweibull(n, shape = 1.5, scale = exp(1 + 0.8 * trt))
#' cens_time <- rexp(n, rate = 0.1)
#' true_obs_time <- pmin(true_time, cens_time)
#' true_status <- as.integer(true_time <= cens_time)
#'
#' # Induce linkage mismatch errors in approximately 20% of records
#' is_mismatch <- rbinom(n, 1, 0.2)
#' obs_time <- true_obs_time
#' obs_status <- true_status
#' mismatch_idx <- which(is_mismatch == 1)
#'
#' shuffled <- sample(mismatch_idx)
#' obs_time[mismatch_idx] <- obs_time[shuffled]
#' obs_status[mismatch_idx] <- obs_status[shuffled]
#'
#' linked_df <- data.frame(time = obs_time, status = obs_status, trt = trt)
#' adj <- adjMixBayes(linked.data = linked_df)
#'
#' fit <- plsurvreg(
#'   survival::Surv(time, status) ~ trt,
#'   dist = "weibull",
#'   adjustment = adj,
#'   control = list(iterations = 200, burnin.iterations = 100, seed = 123)
#' )
#'
#' print(fit)
#'
#' @export
#' @method print survMixBayes
print.survMixBayes <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  obj <- .mixbayes_compat(x)
  est <- obj$estimates
  n_safe <- if (is.null(obj$diagnostics$n_safe)) 0L else obj$diagnostics$n_safe
  among <- if (n_safe > 0) " among records not flagged as safe matches" else ""
  # significant digits, so that small coefficients do not print as 0
  pm <- function(v) print(format(signif(v, digits)), print.gap = 2, quote = FALSE)

  cat("Bayesian two-component mixture survival regression\n")
  cat("Family:", obj$dist, "\n\n")
  if (!is.null(obj$call)) {
    cat("Call:\n")
    print(obj$call)
    cat("\n")
  }
  b <- est$coefficients
  if (is.matrix(b)) {
    cat("Coefficients (posterior means, component 1 = correct-match):\n")
    pm(colMeans(b))
    if (isTRUE(obj$intercept_includes_logscale)) {
      cat("(the intercept is the identified Weibull intercept (Intercept) + log(scale))\n")
    }
    cat("\n")
  }
  if (isTRUE(obj$use_logistic) && !is.null(est$m.coefficients)) {
    cat(paste0("Mismatch Model Coefficients (logit of the mismatch probability", among,
               "; posterior means):\n"))
    pm(colMeans(.mixbayes_block_draws(obj, "m.coefficients")))
    cat("\n")
  } else if (!is.null(est$theta)) {
    cat(paste0("Match probability theta (mixing weight of component 1", among, "; posterior mean): ",
               format(signif(mean(est$theta), digits)), "\n"))
  }
  if (!is.null(est$shape)) {
    cat(paste0("Shape (component 1; posterior mean): ", format(signif(mean(est$shape), digits)), "\n"))
  }
  if (is.matrix(obj$m_samples)) {
    cat(sprintf("Records: %d (safe matches: %d)\n", ncol(obj$m_samples), as.integer(n_safe)))
  }
  invisible(x)
}

#' Summary Method for survMixBayes Models
#'
#' Computes posterior summaries for the regression coefficients, mixing weight,
#' and component-specific distribution parameters in a fitted
#' \code{survMixBayes} model. Throughout, component 1 is interpreted as the
#' correct-match component and component 2 as the incorrect-match component.
#'
#' @param object An object of class \code{survMixBayes}.
#' @param probs Numeric vector of probabilities of the posterior quantiles
#'   reported for every parameter. The default, \code{c(0.025, 0.975)}, gives
#'   the central 95% credible interval (columns \code{"2.5 \%"} and
#'   \code{"97.5 \%"}); add \code{0.5} for the posterior median.
#' @param ... Further arguments (unused).
#'
#' @return An object of class \code{summary.survMixBayes}. Every parameter
#'   block is a table with one row per parameter and the columns
#'   \code{"Estimate"} (the posterior mean), \code{"Std. Error"} (the posterior
#'   standard deviation) and the posterior quantiles given by \code{probs}
#'   (\code{"2.5 \%"} and \code{"97.5 \%"} by default), the columns of
#'   \code{\link{summary.glmMixBayes}}. The tables of the coefficients and the
#'   shape (and \code{1 / shape}) of component 1 and of the match-probability
#'   model add the columns \code{"MCSE"} (the Monte Carlo standard error of the
#'   posterior mean), \code{"ESS"} (the effective sample size) and
#'   \code{"Rhat"} (the single-chain split R-hat; see \emph{Sampling
#'   algorithm} in \code{\link{survregMixBayes}}), computed from the draws
#'   for fits that do not store these diagnostics (also for the notes on
#'   them). The blocks are the regression
#'   coefficients of both mixture components (\code{coef1}, \code{coef2}); the
#'   coefficients of the mismatch-indicator model (\code{m.coefficients}, the
#'   logit of the mismatch probability, as in \code{\link{coxphMixture}()}); the
#'   mixing weight (\code{theta}) or, when the match probability was modelled
#'   with linkage covariates, its logistic regression coefficients
#'   (\code{gamma = -m.coefficients}, for the original linkage covariates,
#'   with \code{z_center}, the centre of the linkage covariates of the fit, and
#'   \code{gamma_intercept_prior}, the \code{gamma_intercept} prior used, which
#'   applies to a record whose linkage covariates equal \code{z_center}; both
#'   are printed under the \code{gamma} table); and the family-specific distribution
#'   parameters (\code{shape1}, \code{shape2} and, for the Weibull
#'   distribution, the scale multipliers \code{scale1}, \code{scale2}, which are
#'   not the \code{Scale} of \code{survival::survreg()}). For the Weibull
#'   distribution, \code{inv_shape1} summarises \code{1 / shape} of component
#'   1, which is comparable to the \code{Scale} reported by
#'   \code{survival::survreg()}; when the design matrix has an intercept
#'   column (\code{intercept_includes_logscale = TRUE}), the intercepts in
#'   \code{coef1} and \code{coef2} are the identified intercepts
#'   \code{(Intercept) + log(scale)} (comparable to the \code{(Intercept)} of
#'   \code{survival::survreg()}), and \code{scale1}, \code{scale2} are not
#'   identified separately from them (their posteriors reflect their priors).
#'   The summary also
#'   contains \code{match.prob}, the posterior probability that each record is
#'   a correct match, \code{match.rate}, a list with the average of
#'   \code{match.prob} over all records (safe matches counted as 1) and over
#'   the records not flagged as safe matches (\code{avg}) and the numbers of
#'   these records (\code{n}), \code{use_logistic}, \code{n_safe}, the
#'   number of known correct matches (with safe matches, \code{theta},
#'   \code{gamma} and \code{m.coefficients} describe the records that are not
#'   flagged as safe), the number of stored draws \code{n_draws}, and the
#'   diagnostics printed as notes (\code{low_ess}, \code{high_rhat},
#'   \code{joint_short} and \code{mixed_labels}; see
#'   \code{\link{summary.glmMixBayes}}). Component 1 corresponds to the correct-match component and
#'   component 2 to the incorrect-match component. Up to postlink 0.1.2 the
#'   blocks were matrices of posterior quantiles with one column per parameter
#'   (by default including the median).
#'
#' @examples
#' set.seed(301)
#' n <- 150
#' trt <- rbinom(n, 1, 0.5)
#'
#' # Simulate Weibull AFT data
#' true_time <- rweibull(n, shape = 1.5, scale = exp(1 + 0.8 * trt))
#' cens_time <- rexp(n, rate = 0.1)
#' true_obs_time <- pmin(true_time, cens_time)
#' true_status <- as.integer(true_time <= cens_time)
#'
#' # Induce linkage mismatch errors in approximately 20% of records
#' is_mismatch <- rbinom(n, 1, 0.2)
#' obs_time <- true_obs_time
#' obs_status <- true_status
#' mismatch_idx <- which(is_mismatch == 1)
#'
#' shuffled <- sample(mismatch_idx)
#' obs_time[mismatch_idx] <- obs_time[shuffled]
#' obs_status[mismatch_idx] <- obs_status[shuffled]
#'
#' linked_df <- data.frame(time = obs_time, status = obs_status, trt = trt)
#' adj <- adjMixBayes(linked.data = linked_df)
#'
#' fit <- plsurvreg(
#'   survival::Surv(time, status) ~ trt,
#'   dist = "weibull",
#'   adjustment = adj,
#'   control = list(iterations = 200, burnin.iterations = 100, seed = 123)
#' )
#'
#' fit_summary <- summary(fit)
#' print(fit_summary)
#'
#' # posterior medians in addition
#' summary(fit, probs = c(0.025, 0.5, 0.975))$coef1
#'
#' @export
#' @method summary survMixBayes
summary.survMixBayes <- function(object, probs = c(0.025, 0.975), ...) {
  .check_probs(probs)
  object <- .mixbayes_compat(object)
  est <- object$estimates
  dg <- object$diagnostics
  # effective sample sizes and split R-hat (computed from the draws for fits
  # saved without them, so that the tables and the notes agree)
  dg <- .mixbayes_mcmc_diagnostics(object, dg)
  # Monte Carlo columns (MCSE, ESS, Rhat) for the tables of component 1 and of
  # the match-probability model, named by their block in posterior_draws()
  tab <- function(draws, rowname, block = NULL) {
    .posterior_table(draws, rowname, probs = probs,
                     mcmc = if (!is.null(block)) list(block = block, ess = dg$ess, rhat = dg$rhat))
  }
  b1 <- est$coefficients
  b2 <- est$coefficients2
  if (!is.matrix(b1) || !is.matrix(b2)) stop("Expected coefficient draws as matrices.", call. = FALSE)

  s <- list(
    call = object$call,
    family = object$dist,
    coef1 = tab(b1, "(Intercept)", "coefficients"),
    coef2 = tab(b2, "(Intercept)")
  )

  # Mismatch-indicator model (as the m.coefficients of coxphMixture()), and
  # the match-probability model it is derived from: the constant mixing weight
  # theta, or the logistic regression coefficients gamma when covariates were
  # supplied via m.formula
  if (!is.null(est$m.coefficients)) s$m.coefficients <- tab(est$m.coefficients, "(Intercept)", "m.coefficients")
  use_logistic <- isTRUE(object$use_logistic)
  if (use_logistic && !is.null(est$gamma)) {
    s$gamma <- tab(est$gamma, "(Intercept)", "gamma")
    # the gamma_intercept prior applies at the centre of the linkage covariates
    s$z_center <- object$z_center
    s$gamma_intercept_prior <- .gamma_intercept_prior(object$priors)
  } else if (!is.null(est$theta)) {
    s$theta <- tab(est$theta, "theta", "theta")
  }

  if (!is.null(est$shape))  s$shape1 <- tab(est$shape, "shape", "shape")
  if (!is.null(est$shape2)) s$shape2 <- tab(est$shape2, "shape2")
  if (!is.null(est$scale))  s$scale1 <- tab(est$scale, "scale")
  if (!is.null(est$scale2)) s$scale2 <- tab(est$scale2, "scale2")
  # Weibull: 1 / shape is the Scale reported by survival::survreg() (its
  # Monte Carlo columns are computed from the draws of 1 / shape)
  if ((identical(object$dist, "weibull") || !is.null(est$scale)) && !is.null(est$shape)) {
    s$inv_shape1 <- tab(1 / est$shape, "1/shape", "1/shape")
  }
  s$use_logistic <- use_logistic
  # Weibull with an intercept column: the intercepts of coef1 and coef2 are
  # the identified (Intercept) + log(scale)
  s$intercept_includes_logscale <- isTRUE(object$intercept_includes_logscale)
  s$n_safe <- if (is.null(object$diagnostics$n_safe)) 0L else object$diagnostics$n_safe
  s$match.prob <- object$match.prob
  if (!is.null(object$match.prob)) s$match.rate <- .match_rate_summary(object$match.prob, s$n_safe)
  s$low_ess <- .low_ess(dg$ess)
  s$high_rhat <- .high_rhat(dg$rhat)
  s$n_draws <- NROW(b1)
  s$joint_short <- .joint_moves_short(dg)
  s$mixed_labels <- .mixed_labels(dg)

  class(s) <- "summary.survMixBayes"
  s
}

#' @noRd
#' @export
print.summary.survMixBayes <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("Summary of Bayesian mixture survival regression\n")
  if (!is.null(x$call)) {
    cat("\nCall:\n")
    print(x$call)
  }
  cat("\nFamily:", x$family, "\n")

  fmt <- function(tab, nm) {
    if (is.null(tab)) return(invisible(NULL))
    cat("\n", nm, ":\n", sep = "")
    # summaries of postlink <= 0.1.2: posterior quantiles with one column per parameter
    if (!is.matrix(tab)) tab <- matrix(tab, nrow = 1L, dimnames = list("", names(tab)))
    else if (!"Estimate" %in% colnames(tab)) tab <- t(tab)
    # significant digits, so that small values do not print as 0 (the
    # effective sample size as a whole number, the split R-hat with three
    # decimals)
    .print_posterior_table(tab, digits)
  }
  heads <- .match_headings(x$n_safe > 0)

  fmt(x$coef1, "Coefficients (component 1 = correct-match)")
  fmt(x$coef2, "Coefficients (component 2 = incorrect-match)")
  if (isTRUE(x$intercept_includes_logscale)) {
    cat("(Weibull intercepts: the identified (Intercept) + log(scale), comparable to the (Intercept) of survival::survreg().)\n")
  }

  if (!is.null(x$m.coefficients)) {
    fmt(x$m.coefficients, heads[["m.coefficients"]])
  }
  if (isTRUE(x$use_logistic) && !is.null(x$gamma)) {
    fmt(x$gamma, heads[["gamma"]])
    .print_z_center_note(x$z_center, x$gamma_intercept_prior, digits)
  } else if (!is.null(x$theta)) {
    fmt(x$theta, heads[["theta"]])
  }

  fmt(x$shape1, "Shape (component 1 = correct-match)")
  fmt(x$shape2, "Shape (component 2 = incorrect-match)")
  fmt(x$inv_shape1, "1 / shape (component 1 = correct-match; comparable to the Scale of survival::survreg())")
  scale_lab <- if (isTRUE(x$intercept_includes_logscale)) {
    "Scale multiplier of exp(X beta) (component %d = %s; not identified separately from the intercept: reflects its prior)"
  } else {
    "Scale multiplier of exp(X beta) (component %d = %s)"
  }
  fmt(x$scale1, sprintf(scale_lab, 1L, "correct-match"))
  fmt(x$scale2, sprintf(scale_lab, 2L, "incorrect-match"))
  if (!is.null(x$scale1)) {
    cat("(The Weibull scale multiplies exp(X beta); it is not survreg()'s Scale, which is 1 / shape.)\n")
  }
  if (!is.null(x$match.rate)) {
    cat("\n")
    .print_match_rate(x$match.rate, digits)
  }
  if (!is.null(x$coef1) && "Std. Error" %in% colnames(x$coef1)) {
    if ("MCSE" %in% colnames(x$coef1)) {
      cat(paste0("\n(Estimate: posterior mean; Std. Error: posterior standard deviation; columns labelled %: ",
                 "posterior quantiles;\n MCSE: Monte Carlo standard error of the estimate; ESS: effective sample ",
                 "size; Rhat: split R-hat of the single chain)\n"))
    } else {
      cat("\n(Estimate: posterior mean; Std. Error: posterior standard deviation; other columns: posterior quantiles)\n")
    }
  }
  if (length(x$low_ess) > 0L || length(x$high_rhat) > 0L || !is.null(x$joint_short) ||
      !is.null(x$mixed_labels)) cat("\n")
  .print_low_ess(x$low_ess)
  .print_high_rhat(x$high_rhat, x$n_draws)
  .print_joint_short(x$joint_short)
  .print_mixed_labels(x$mixed_labels, "survregMixBayes")

  invisible(x)
}

# Names of the parameter blocks of a survMixBayes fit in the list returned by
# confint() (and in summary()), named by the blocks of `estimates`.
.surv_list_names <- c(coefficients = "coef1", coefficients2 = "coef2", m.coefficients = "m.coefficients",
                      theta = "theta", gamma = "gamma", shape = "shape1", shape2 = "shape2",
                      scale = "scale1", scale2 = "scale2")

# Block of `estimates` for a block name of confint.survMixBayes() or
# vcov.survMixBayes(): the names of `estimates` and the list names coef1,
# coef2, shape1, scale1 are accepted; NA for unknown names.
.surv_block <- function(block) {
 if (block %in% names(.surv_list_names)) return(block)
 names(.surv_list_names)[match(block, .surv_list_names)]
}

#' Credible Intervals for Parameters from a survMixBayes Fit
#'
#' Computes posterior credible intervals for the regression coefficients,
#' mixing weight, and family-specific distribution parameters from a fitted
#' \code{survMixBayes} model. By default the intervals of all parameter blocks
#' are returned as a named list; \code{block} returns one block as a matrix,
#' and \code{parm} with names or indices of coefficients of component 1
#' returns those rows as a matrix, as \code{\link{confint.glmMixBayes}} does.
#' Component 1 corresponds to the correct-match component and component 2 to
#' the incorrect-match component.
#'
#' @param object An object of class \code{survMixBayes}.
#' @param parm Optional. Either a character vector of block names selecting
#'   elements of the returned list: \code{"coef1"}, \code{"coef2"},
#'   \code{"m.coefficients"} (the mismatch-indicator model, as in
#'   \code{\link{coxphMixture}()}), \code{"theta"} (constant match probability)
#'   or \code{"gamma"} (match probability modelled with linkage covariates),
#'   \code{"shape1"}, \code{"shape2"} and, for Weibull fits, \code{"scale1"},
#'   \code{"scale2"} (for example, \code{"theta"} returns only the credible
#'   interval for the mixing weight); the block names accepted by
#'   \code{block} (e.g. \code{"coefficients"}, \code{"shape"}) select the
#'   corresponding list elements too. Or names or indices of coefficients of
#'   component 1, which return the matrix of those rows (unnamed columns of
#'   \code{X} are named \code{"(Intercept)"} or \code{"X<j>"} by
#'   \code{survregMixBayes()}; in fits saved by earlier versions, coefficients
#'   without names are named by their indices). With \code{block},
#'   names or indices of parameters within that block. If \code{NULL},
#'   credible intervals are returned for all available parameter blocks;
#'   unknown names raise an error.
#' @param level Probability level for the credible intervals. Defaults to
#'   \code{0.95}.
#' @param block Optional name of one parameter block, returned as a matrix:
#'   \code{"coefficients"} (component 1), \code{"coefficients2"} (component 2),
#'   \code{"m.coefficients"}, \code{"theta"}, \code{"gamma"}, \code{"shape"},
#'   \code{"shape2"}, \code{"scale"} or \code{"scale2"} (the names of
#'   \code{object$estimates}, as for \code{\link{confint.glmMixBayes}}); the
#'   list names \code{"coef1"}, \code{"coef2"}, \code{"shape1"} and
#'   \code{"scale1"} are accepted too.
#' @param ... Not used; other arguments raise an error.
#'
#' @return Without \code{block}, and unless \code{parm} names coefficients of
#'   component 1, a named list of credible intervals. Elements \code{coef1} and
#'   \code{coef2} are matrices with one row per regression coefficient and two
#'   columns giving the lower and upper interval bounds for components 1 and 2,
#'   respectively, where component 1 is the correct-match component and
#'   component 2 is the incorrect-match component. Elements such as \code{theta}, \code{shape1},
#'   \code{shape2}, \code{scale1}, and \code{scale2} are numeric vectors of
#'   length 2 containing the lower and upper credible interval bounds for the
#'   corresponding scalar parameters. When the match probability is modeled
#'   with covariates (\code{m.formula} in \code{\link{adjMixBayes}}), a
#'   \code{gamma} matrix (one row per logistic regression coefficient)
#'   replaces \code{theta}. The last element, \code{m.coefficients}, is the
#'   matrix of the coefficients of the mismatch-indicator model (on the logit
#'   scale of the mismatch probability; \code{-gamma}, or
#'   \code{qlogis(1 - theta)} without covariates). With \code{block}, or with
#'   coefficient names or indices in \code{parm}, a matrix with one row per
#'   parameter. Interval bounds are labelled as by \code{stats::confint()}
#'   (\code{"2.5 \%"} and \code{"97.5 \%"} for \code{level = 0.95}). For
#'   Weibull fits whose design matrix has an intercept column, the intercept
#'   rows of \code{coef1} and \code{coef2} refer to the identified intercepts
#'   \code{(Intercept) + log(scale)} (see \code{\link{survregMixBayes}}).
#'
#' @examples
#' set.seed(301)
#' n <- 150
#' trt <- rbinom(n, 1, 0.5)
#'
#' # Simulate Weibull AFT data
#' true_time <- rweibull(n, shape = 1.5, scale = exp(1 + 0.8 * trt))
#' cens_time <- rexp(n, rate = 0.1)
#' true_obs_time <- pmin(true_time, cens_time)
#' true_status <- as.integer(true_time <= cens_time)
#'
#' # Induce linkage mismatch errors in approximately 20% of records
#' is_mismatch <- rbinom(n, 1, 0.2)
#' obs_time <- true_obs_time
#' obs_status <- true_status
#' mismatch_idx <- which(is_mismatch == 1)
#'
#' shuffled <- sample(mismatch_idx)
#' obs_time[mismatch_idx] <- obs_time[shuffled]
#' obs_status[mismatch_idx] <- obs_status[shuffled]
#'
#' linked_df <- data.frame(time = obs_time, status = obs_status, trt = trt)
#'
#' adj <- adjMixBayes(linked.data = linked_df)
#'
#' fit <- plsurvreg(
#'   survival::Surv(time, status) ~ trt,
#'   dist = "weibull",
#'   adjustment = adj,
#'   control = list(iterations = 200, burnin.iterations = 100, seed = 123)
#' )
#'
#' # Calculate 95% credible intervals for all parameters
#' confint(fit, level = 0.95)
#'
#' # Extract credible intervals specifically for the mixing weight
#' confint(fit, parm = "theta", level = 0.90)
#'
#' # One block as a matrix, e.g. the mismatch-indicator model (logit of the
#' # mismatch probability, as the m.coefficients of adjMixture() fits)
#' confint(fit, block = "m.coefficients")
#'
#' # Coefficients of component 1 by name, as a matrix
#' confint(fit, parm = "trt")
#'
#' @export
#' @method confint survMixBayes
confint.survMixBayes <- function(object, parm = NULL, level = 0.95, block = NULL, ...) {
 .check_unused_args(list(...), "confint.survMixBayes", "Select interval blocks with `parm` or `block`.")
 if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1) {
  stop("`level` must be a single number strictly between 0 and 1.", call. = FALSE)
 }
 object <- .mixbayes_compat(object)
 alpha <- (1 - level) / 2
 probs <- c(alpha, 1 - alpha)
 q <- function(b) .draw_quantiles(.mixbayes_block_draws(object, b), probs)
 # names of the blocks in the returned list
 list_name <- .surv_list_names

 # block = "<name>": that block as a matrix (parm then selects within it)
 if (!is.null(block)) {
  if (!is.character(block) || length(block) != 1L || is.na(block)) {
   stop("`block` must be a single block name.", call. = FALSE)
  }
  b <- .surv_block(block)
  if (is.na(b)) {
   stop("Unknown block '", block, "'. Available: ", paste(names(list_name), collapse = ", "), ".",
        call. = FALSE)
  }
  return(.select_parm(q(b), parm))
 }

 # parm with coefficient names or indices of component 1: those rows as a
 # matrix, as in confint.glmMixBayes() (survregMixBayes() names unnamed
 # columns of X; coefficients without names, in fits saved by earlier
 # versions, are named by their indices)
 cn <- colnames(.mixbayes_block_draws(object, "coefficients"))
 if (!is.null(parm) && (is.numeric(parm) ||
                        (is.character(parm) && !any(parm %in% list_name) && all(parm %in% cn)))) {
  return(.select_parm(q("coefficients"), parm))
 }

 # all blocks as a named list (the layout of postlink 0.1.2, with the
 # mismatch-indicator model added at the end)
 scalar <- function(b) {
  v <- q(b)
  stats::setNames(as.numeric(v), colnames(v))
 }
 est <- object$estimates
 out <- list()
 out$coef1 <- q("coefficients")
 out$coef2 <- q("coefficients2")
 if (isTRUE(object$use_logistic) && !is.null(est$gamma)) {
  out$gamma <- q("gamma")
 } else if (!is.null(est$theta)) {
  out$theta <- scalar("theta")
 }
 if (!is.null(est$shape))  out$shape1 <- scalar("shape")
 if (!is.null(est$shape2)) out$shape2 <- scalar("shape2")
 if (!is.null(est$scale))  out$scale1 <- scalar("scale")
 if (!is.null(est$scale2)) out$scale2 <- scalar("scale2")
 if (!is.null(est$m.coefficients)) out$m.coefficients <- q("m.coefficients")

 # Optional filtering: select interval blocks by name
 if (!is.null(parm)) {
  if (!is.character(parm)) {
   stop("`parm` must be block names or names or indices of coefficients of component 1.",
        call. = FALSE)
  }
  # the names of `estimates` accepted by `block` (e.g. "coefficients",
  # "shape") select the list elements of those blocks
  given <- parm
  est_name <- parm %in% names(list_name)
  parm[est_name] <- unname(list_name[parm[est_name]])
  if ("theta" %in% parm && is.null(out$theta) && !is.null(out$gamma)) {
   stop("No 'theta' draws: the match probability was modeled with covariates. ",
        "Use parm = \"gamma\" (or \"m.coefficients\").", call. = FALSE)
  }
  if ("gamma" %in% parm && is.null(out$gamma) && !is.null(out$theta)) {
   stop("No 'gamma' draws: the match probability was not modeled with covariates. ",
        "Use parm = \"theta\" (or \"m.coefficients\").", call. = FALSE)
  }
  bad <- !parm %in% names(out)
  if (any(bad) && all(given[bad] %in% cn)) {
   stop("`parm` combines coefficient names of component 1 (", paste(given[bad], collapse = ", "),
        ") with block names; select coefficients and blocks in separate calls.", call. = FALSE)
  }
  if (any(bad)) {
   stop("Unknown block(s) in `parm`: ", paste(given[bad], collapse = ", "), ". Available: ",
        paste(names(out), collapse = ", "), " (or names or indices of the coefficients of ",
        "component 1: ", paste(cn, collapse = ", "), ").", call. = FALSE)
  }
  out <- out[names(out) %in% parm]
 }

 out
}

#' Posterior Covariance Matrix for survMixBayes Coefficients
#'
#' Returns the empirical posterior covariance matrix of the regression
#' coefficients for component 1 of a fitted \code{survMixBayes} model, or of
#' another parameter block. In this package, component 1 is interpreted as the
#' correct-match component.
#'
#' @param object A \code{survMixBayes} model object.
#' @param block Which parameter block: \code{"coefficients"} (component 1;
#'   default), \code{"coefficients2"} (component 2), \code{"m.coefficients"}
#'   (the mismatch-indicator model), \code{"theta"}, \code{"gamma"},
#'   \code{"shape"}, \code{"shape2"}, \code{"scale"} or \code{"scale2"}; the
#'   names of the list returned by \code{\link{confint.survMixBayes}}
#'   (\code{"coef1"}, \code{"coef2"}, \code{"shape1"}, \code{"scale1"}) are
#'   accepted too, as by its \code{block} argument.
#' @param ... Not used; other arguments raise an error.
#'
#' @return Posterior covariance matrix of the parameters of the selected block
#'   (by default the regression coefficients of component 1, interpreted as
#'   the correct-match component), i.e. the covariance of their posterior
#'   draws. For Weibull fits whose design matrix has an intercept column, the
#'   intercept is the identified \code{(Intercept) + log(scale)} (see
#'   \code{\link{survregMixBayes}}).
#'
#' @examples
#' set.seed(301)
#' n <- 150
#' trt <- rbinom(n, 1, 0.5)
#'
#' # Simulate Weibull AFT data
#' true_time <- rweibull(n, shape = 1.5, scale = exp(1 + 0.8 * trt))
#' cens_time <- rexp(n, rate = 0.1)
#' true_obs_time <- pmin(true_time, cens_time)
#' true_status <- as.integer(true_time <= cens_time)
#'
#' # Induce linkage mismatch errors in approximately 20% of records
#' is_mismatch <- rbinom(n, 1, 0.2)
#' obs_time <- true_obs_time
#' obs_status <- true_status
#' mismatch_idx <- which(is_mismatch == 1)
#'
#' shuffled <- sample(mismatch_idx)
#' obs_time[mismatch_idx] <- obs_time[shuffled]
#' obs_status[mismatch_idx] <- obs_status[shuffled]
#'
#' linked_df <- data.frame(time = obs_time, status = obs_status, trt = trt)
#'
#' adj <- adjMixBayes(linked.data = linked_df)
#'
#' fit <- plsurvreg(
#'   survival::Surv(time, status) ~ trt,
#'   dist = "weibull",
#'   adjustment = adj,
#'   control = list(iterations = 200, burnin.iterations = 100, seed = 123)
#' )
#'
#' # Extract the empirical posterior covariance matrix for component 1
#' vcov_mat <- vcov(fit)
#' print(vcov_mat)
#'
#' # posterior covariance of the coefficients of the mismatch-indicator model
#' vcov(fit, block = "m.coefficients")
#'
#' @export
#' @method vcov survMixBayes
vcov.survMixBayes <- function(object,
                              block = c("coefficients", "coefficients2", "m.coefficients", "theta",
                                        "gamma", "shape", "shape2", "scale", "scale2"),
                              ...) {
  .check_unused_args(list(...), "vcov.survMixBayes", "Select a parameter block with `block`.")
  # the list names of confint() (coef1, coef2, shape1, scale1) are accepted too
  if (is.character(block) && length(block) == 1L && !is.na(block) && !is.na(.surv_block(block))) {
   block <- .surv_block(block)
  }
  block <- match.arg(block)
  stats::cov(.mixbayes_block_draws(.mixbayes_compat(object), block))
}

#' Predictions from a survMixBayes Model
#'
#' Computes posterior predictions for each latent component of a
#' \code{survMixBayes} model. By default, predictions are returned on the
#' linear predictor scale for both components.
#'
#' For the gamma model the linear predictor is \eqn{X\beta}{X beta}, the log of
#' the mean survival time. For the Weibull model it is
#' \eqn{\log(\mathrm{scale}) + X\beta}{log(scale) + X beta}, the log of each
#' record's Weibull scale parameter (the scale of
#' \code{predict(survival::survreg(...), type = "lp")}). When the design
#' matrix has an intercept column, the reported intercepts already include
#' \eqn{\log(\mathrm{scale})}{log(scale)} (the identified intercepts; see
#' \code{\link{survregMixBayes}}), so the linear predictor is \eqn{X\beta}{X beta}
#' with these coefficients; otherwise \eqn{\log(\mathrm{scale})}{log(scale)}
#' is added to \eqn{X\beta}{X beta}.
#'
#' Component 1 is interpreted as the correct-match component and
#' component 2 as the incorrect-match component (after label-switching
#' correction).
#'
#' @param object A \code{survMixBayes} model object.
#' @param newdata Optional new data: a data frame in which to look for the
#'   covariates of the outcome model, whose model matrix is built from the
#'   terms of a fit from \code{plsurvreg()} (stored with \code{model = TRUE},
#'   the default), as by \code{\link{predict.glmMixBayes}}; or a numeric
#'   matrix of new observations (\eqn{n_{new} \times K}{n_new x K}) with
#'   columns aligned to the design matrix used for fitting (matched by name
#'   when the coefficient names are non-empty and unique and each names
#'   exactly one column, otherwise by position, with a warning when a column
#'   bears the name of a coefficient at another position), which fits from
#'   \code{survregMixBayes()} need. If \code{NULL}, predictions are made for
#'   the analysed records of the fit, from the design matrix stored by
#'   \code{plsurvreg(..., x = TRUE)} or else from the stored model frame
#'   (records dropped at fit time because their linkage covariates are missing
#'   are left out).
#' @param se.fit Logical; if \code{TRUE}, also return posterior SD of predictions.
#' @param interval Either \code{"none"} or \code{"credible"}, indicating whether
#'   to compute credible intervals.
#' @param level Probability level for the credible interval (default 0.95).
#' @param na.action Function determining what to do with missing values in a
#'   data frame \code{newdata} (default \code{stats::na.pass}, which gives
#'   \code{NA} predictions).
#' @param ... Not used; other arguments (e.g. \code{type}) raise an error.
#'
#' @return A list with two components, \code{component1} and \code{component2},
#'   corresponding to the two latent mixture components.
#'   If \code{se.fit = FALSE} and \code{interval = "none"}, each element is a
#'   numeric vector of posterior mean linear predictors.
#'   Otherwise, each element is a matrix containing the fitted values and,
#'   optionally, posterior SDs and credible interval bounds. Predictions are
#'   named after the rows of the design matrix (the row names of a data frame
#'   in \code{newdata}, or the analysed records), as by
#'   \code{\link[stats]{predict.glm}()}; rows with missing covariates give
#'   \code{NA}.
#'
#' @examples
#' set.seed(301)
#' n <- 150
#' trt <- rbinom(n, 1, 0.5)
#'
#' # Simulate Weibull AFT data
#' true_time <- rweibull(n, shape = 1.5, scale = exp(1 + 0.8 * trt))
#' cens_time <- rexp(n, rate = 0.1)
#' true_obs_time <- pmin(true_time, cens_time)
#' true_status <- as.integer(true_time <= cens_time)
#'
#' # Induce linkage mismatch errors in approximately 20% of records
#' is_mismatch <- rbinom(n, 1, 0.2)
#' obs_time <- true_obs_time
#' obs_status <- true_status
#' mismatch_idx <- which(is_mismatch == 1)
#'
#' shuffled <- sample(mismatch_idx)
#' obs_time[mismatch_idx] <- obs_time[shuffled]
#' obs_status[mismatch_idx] <- obs_status[shuffled]
#'
#' linked_df <- data.frame(time = obs_time, status = obs_status, trt = trt)
#' adj <- adjMixBayes(linked.data = linked_df)
#'
#' fit <- plsurvreg(
#'   survival::Surv(time, status) ~ trt,
#'   dist = "weibull",
#'   adjustment = adj,
#'   control = list(
#'     iterations = 200,
#'     burnin.iterations = 100,
#'     seed = 123
#'   )
#' )
#'
#' # Predict posterior mean linear predictors for each latent component
#' preds <- predict(fit, newdata = data.frame(trt = c(0, 1)), se.fit = TRUE,
#'                  interval = "credible")
#' print(preds$component1)
#' print(preds$component2)
#'
#' # the same with a design matrix
#' newx <- stats::model.matrix(~ trt, data = data.frame(trt = c(0, 1)))
#' predict(fit, newdata = newx)$component1
#'
#' @export
#' @method predict survMixBayes
predict.survMixBayes <- function(object, newdata = NULL,
                                 se.fit = FALSE,
                                 interval = c("none", "credible"),
                                 level = 0.95,
                                 na.action = stats::na.pass,
                                 ...) {
 interval <- match.arg(interval)
 extra <- list(...)
 if (length(extra) > 0L) {
  nm <- names(extra)
  if (is.null(nm)) nm <- rep("", length(extra))
  nm[!nzchar(nm)] <- "(unnamed)"
  stop("Unused argument(s) in predict.survMixBayes(): ", paste(nm, collapse = ", "),
       ". Predictions are returned on the linear-predictor scale; pass new data as `newdata`.",
       call. = FALSE)
 }

 if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1) {
  stop("`level` must be a single number strictly between 0 and 1.", call. = FALSE)
 }

 X <- newdata
 # a data frame: model matrix from the stored terms, as in predict.glmMixBayes()
 if (is.data.frame(X)) X <- .mixbayes_newdata_matrix(object, X, na.action, matrix_arg = "newdata")
 # no new data: the analysed records of the fit, from the stored design
 # matrix or model frame
 if (is.null(X)) X <- .mixbayes_fit_design(object)

 if (is.null(X)) {
  stop("No new data were given and the fit stores neither its design matrix nor its model frame: ",
       "pass the model matrix as `newdata` (or refit with plsurvreg(..., x = TRUE)).", call. = FALSE)
 }

 if (!is.matrix(X) || !is.numeric(X)) {
  stop(
   "`newdata` must be a data frame or a numeric matrix with columns aligned to the design matrix ",
   "used for fitting.",
   call. = FALSE
  )
 }

 object <- .mixbayes_compat(object)
 b1 <- object$estimates$coefficients
 b2 <- object$estimates$coefficients2

 if (!is.matrix(b1) || !is.matrix(b2)) {
  stop(
   "Posterior coefficient draws are not stored in the expected matrix format.",
   call. = FALSE
  )
 }

 X <- .align_newx(X, b1, "newdata")

 pred1 <- X %*% t(b1)
 pred2 <- X %*% t(b2)
 # Weibull: the scale parameter multiplies exp(X beta), so the identified
 # linear predictor is log(scale) + X beta, the log of each record's Weibull
 # scale. When X has an intercept column the reported intercepts already
 # include log(scale) (intercept_includes_logscale); otherwise, and for fits
 # stored before the intercepts included it, log(scale) is added here.
 if (!isTRUE(object$intercept_includes_logscale)) {
  if (!is.null(object$estimates$scale))  pred1 <- sweep(pred1, 2, log(object$estimates$scale), "+")
  if (!is.null(object$estimates$scale2)) pred2 <- sweep(pred2, 2, log(object$estimates$scale2), "+")
 }

 # predictions are named after the rows of the design matrix (the rows of a
 # data frame in newdata, or the analysed records), as by stats::predict.glm()
 summarize_component <- function(all_predictions, se.fit, interval, level) {
  fit <- rowMeans(all_predictions)
  names(fit) <- rownames(X)

  if (!se.fit && interval == "none") {
   return(fit)
  }

  out <- cbind(fit = fit)

  if (se.fit) {
   out <- cbind(out, se.fit = apply(all_predictions, 1, stats::sd))
  }

  if (interval == "credible") {
   alpha <- 1 - level
   # rows of new data with missing covariates (na.action = na.pass) are NA
   ci <- t(apply(
    all_predictions, 1,
    stats::quantile,
    probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE
   ))
   out <- cbind(out, lower = ci[, 1], upper = ci[, 2])

   labels <- .ci_labels(c(alpha / 2, 1 - alpha / 2))
   colnames(out) <- if (se.fit) c("fit", "se.fit", labels) else c("fit", labels)
  } else {
   colnames(out) <- c("fit", "se.fit")
  }

  out
 }

 list(
  component1 = summarize_component(pred1, se.fit, interval, level),
  component2 = summarize_component(pred2, se.fit, interval, level)
 )
}

#' Pool Regression Fits Across Posterior Draws of Correct-Match Classifications
#'
#' @description
#' Use posterior draws of the latent match indicators from \code{survregMixBayes()}
#' to repeatedly identify which records are treated as correct matches, refit a
#' Cox proportional hazards model on those records, and pool the resulting
#' estimates using multiple-imputation pooling rules.
#'
#' Each retained posterior draw defines one subset of records classified as
#' correct matches. The function fits the specified \code{survival::coxph()}
#' model to that subset, extracts the estimated coefficients and covariance
#' matrix, and combines the results across draws using Rubin's rules. At
#' least two draws must be usable, since the between-draw variance cannot be
#' estimated from one; fewer raise an error.
#'
#' The pooled estimates are therefore Cox log hazard ratios, on the scale of
#' \code{\link{coxphMixture}()}, and not the accelerated failure time
#' coefficients returned by \code{coef()} and \code{summary()} of the mixture
#' fit. For \code{dist = "weibull"} the two describe the same effects with
#' opposite signs: the log hazard ratios are about \code{-shape * beta}
#' (compare with \code{colMeans(-fit$estimates$shape * fit$estimates$coefficients)}).
#' For \code{dist = "gamma"} the components are not proportional hazards
#' models, so there is no exact correspondence.
#'
#' @param object A \code{survMixBayes} model object containing posterior draws of
#'   the latent match indicators.
#' @param data A data.frame with the records used in the model. For fits from
#'   \code{plglm()} / \code{plsurvreg()} this is normally the \code{linked.data}
#'   of the adjustment object: its rows are matched to the analysed records by
#'   row name, so records dropped by \code{subset} or \code{na.action} are
#'   handled correctly. Otherwise \code{data} must contain exactly the analysed
#'   records, in the order used by the model.
#' @param formula Model formula for refitting on each draw, typically of the
#'   form \code{survival::Surv(time, event) ~ ...}. If omitted, the outcome
#'   formula of a \code{plsurvreg()} fit is used (taken from the stored model
#'   frame, or from a formula written directly in the call). Data-dependent
#'   terms such as \code{poly()}, \code{scale()} or \code{splines::ns()} keep
#'   in every refit the basis computed once from all analysed records in
#'   \code{data} (for the stored outcome formula, the basis of the fit), as in
#'   \code{\link{mi_with.glmMixBayes}}. Only top-level terms of this kind are
#'   fixed (those whose prediction call R records): a data-dependent call
#'   nested in another one, such as \code{I(scale(x)^2)}, is evaluated from
#'   the records of each draw, so compute such variables in \code{data}
#'   beforehand. Penalised \code{pspline()} terms keep
#'   the boundary knots and number of knots of all records, and are refitted
#'   in each draw with the penalty arguments as written (\code{df},
#'   \code{theta}, \code{method}, ...).
#' @param min_n Minimum number of records required to fit the model for a given
#'   posterior draw. The default is \code{p + 2}, where \code{p} is the number
#'   of non-intercept columns in the model matrix.
#' @param quietly If \code{TRUE} (a warning then reports how many draws were
#'   not used), draws that lead to fitting errors are skipped
#'   without printing the full error message.
#' @param ties Method for handling tied event times in \code{survival::coxph()}.
#'   Default is \code{"efron"}.
#' @param ... Additional arguments passed to \code{survival::coxph()}.
#'
#' @return An object of class \code{c("mi_link_pool_survreg", "mi_link_pool")}
#'   containing pooled coefficient estimates, standard errors, confidence
#'   intervals, and related summary information.
#'
#' @examples
#' set.seed(301)
#' n <- 150
#' trt <- rbinom(n, 1, 0.5)
#'
#' # Simulate Weibull AFT data
#' true_time <- rweibull(n, shape = 1.5, scale = exp(1 + 0.8 * trt))
#' cens_time <- rexp(n, rate = 0.1)
#' true_obs_time <- pmin(true_time, cens_time)
#' true_status <- as.integer(true_time <= cens_time)
#'
#' # Induce linkage mismatch errors in approximately 20% of records
#' is_mismatch <- rbinom(n, 1, 0.2)
#' obs_time <- true_obs_time
#' obs_status <- true_status
#' mismatch_idx <- which(is_mismatch == 1)
#'
#' shuffled <- sample(mismatch_idx)
#' obs_time[mismatch_idx] <- obs_time[shuffled]
#' obs_status[mismatch_idx] <- obs_status[shuffled]
#'
#' linked_df <- data.frame(time = obs_time, status = obs_status, trt = trt)
#' adj <- adjMixBayes(linked.data = linked_df)
#'
#' fit <- plsurvreg(
#'   survival::Surv(time, status) ~ trt,
#'   dist = "weibull",
#'   adjustment = adj,
#'   control = list(iterations = 200, burnin.iterations = 100, seed = 123)
#' )
#'
#' pooled_obj <- mi_with(
#'   object = fit,
#'   data = linked_df,
#'   formula = survival::Surv(time, status) ~ trt
#' )
#'
#' print(pooled_obj)
#'
#' @export
#' @method mi_with survMixBayes
#' @importFrom stats coef vcov cov qt model.frame model.matrix
#' @importFrom survival coxph Surv
mi_with.survMixBayes <- function(object, data, formula,
                                 min_n = NULL, quietly = TRUE,
                                 ties = "efron", ...) {

 if (missing(formula) || is.null(formula)) {
  formula <- .mixbayes_mi_formula(object, parent.frame())
 }

 if (!inherits(formula, "formula")) {
  stop("`formula` must be a valid formula object.", call. = FALSE)
 }

 # Rows of `data` aligned with the records used in the fit
 data <- .mixbayes_mi_data(object, data)

 z_samples <- object$m_samples

 collapse_z_g <- function(z_samples, g) {
  if (is.list(z_samples)) {
   parts <- lapply(z_samples, function(Z) {
    if (is.null(dim(Z))) {
     stop("Each list element of `z_samples` must be a matrix-like object.", call. = FALSE)
    }
    Z[g, ]
   })
   as.integer(do.call(c, parts))
  } else {
   if (is.null(dim(z_samples)) || length(dim(z_samples)) != 2L) {
    stop("`z_samples` must be an S x N matrix, or a list of such matrices.", call. = FALSE)
   }
   as.integer(z_samples[g, ])
  }
 }

 # Terms of the refitted model, resolved once from all analysed records (the
 # stored terms of the fit when no formula is given): their "predvars" keep
 # the bases of data-dependent terms (poly(), scale(), splines::ns(), ...)
 # fixed in every refit. coxph() does not accept a terms object, so these
 # terms are evaluated once and stored in the data, or, for penalised terms
 # such as pspline(), kept as written in the formula with the knots of all
 # records added (see .mi_fixed_basis()).
 mf <- stats::model.frame(formula, data = data)
 tt <- attr(mf, "terms")
 fixed <- .mi_fixed_basis(tt, data)

 fit_once <- function(df) {
  fit <- survival::coxph(formula = fixed$formula, data = df, ties = ties, ...)
  cf <- stats::coef(fit)
  V <- stats::vcov(fit)
  names(cf) <- .mi_restore_names(names(cf), fixed$labels)
  if (!is.null(dimnames(V))) dimnames(V) <- lapply(dimnames(V), .mi_restore_names, labels = fixed$labels)
  list(coef = cf, vcov = V)
 }

 Xtmp <- stats::model.matrix(tt, mf)
 if ("(Intercept)" %in% colnames(Xtmp)) {
  Xtmp <- Xtmp[, colnames(Xtmp) != "(Intercept)", drop = FALSE]
 }
 p <- ncol(Xtmp)
 if (is.null(min_n)) min_n <- p + 2L

 S <- if (is.list(z_samples)) nrow(z_samples[[1]]) else nrow(z_samples)
 if (length(S) == 0L || is.null(S) || !is.finite(S)) {
  stop("Unable to infer the number of posterior draws (S).", call. = FALSE)
 }

 coefs_list <- list()
 vcovs_list <- list()
 kept <- logical(S)
 n_small <- 0L
 first_error <- NULL

 for (s in seq_len(S)) {
  zs <- collapse_z_g(z_samples, s)

  if (!all(zs %in% c(1L, 2L))) {
   stop("`m_samples` must contain only 1/2 indicators.", call. = FALSE)
  }

  idx <- which(zs == 1L)
  if (length(idx) < min_n) {
   n_small <- n_small + 1L
   next
  }

  dat_s <- fixed$data[idx, , drop = FALSE]

  res <- try(fit_once(dat_s), silent = quietly)
  if (!inherits(res, "try-error") &&
      is.numeric(res$coef) &&
      is.matrix(res$vcov) &&
      length(res$coef) > 0L) {
   coefs_list[[length(coefs_list) + 1L]] <- res$coef
   vcovs_list[[length(vcovs_list) + 1L]] <- res$vcov
   kept[s] <- TRUE
  } else if (is.null(first_error)) {
   first_error <- if (inherits(res, "try-error")) conditionMessage(attr(res, "condition")) else
    "the refitted Cox model returned no coefficients"
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
 r[!is.finite(r)] <- 0
 dfold <- (m - 1) * (1 + 1 / pmax(r, .Machine$double.eps))^2
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
  dist       = object$dist,
  refit      = "coxph",
  call       = object$call
 )

 class(out) <- c("mi_link_pool_survreg", "mi_link_pool")
 out
}

#' Print Pooled Cox Regression Results
#'
#' @param x An object of class \code{mi_link_pool_survreg}, typically returned by
#'   \code{mi_with()} for a \code{survMixBayes} fit.
#' @param digits the number of significant digits to print.
#' @param ... further arguments (unused).
#'
#' @return The input \code{x}, invisibly.
#'
#' @examples
#' set.seed(301)
#' n <- 150
#' trt <- rbinom(n, 1, 0.5)
#'
#' # Simulate Weibull AFT data
#' true_time <- rweibull(n, shape = 1.5, scale = exp(1 + 0.8 * trt))
#' cens_time <- rexp(n, rate = 0.1)
#' true_obs_time <- pmin(true_time, cens_time)
#' true_status <- as.integer(true_time <= cens_time)
#'
#' # Induce linkage mismatch errors in approximately 20% of records
#' is_mismatch <- rbinom(n, 1, 0.2)
#' obs_time <- true_obs_time
#' obs_status <- true_status
#' mismatch_idx <- which(is_mismatch == 1)
#'
#' shuffled <- sample(mismatch_idx)
#' obs_time[mismatch_idx] <- obs_time[shuffled]
#' obs_status[mismatch_idx] <- obs_status[shuffled]
#'
#' linked_df <- data.frame(time = obs_time, status = obs_status, trt = trt)
#' adj <- adjMixBayes(linked.data = linked_df)
#'
#' fit <- plsurvreg(
#'   survival::Surv(time, status) ~ trt,
#'   dist = "weibull",
#'   adjustment = adj,
#'   control = list(iterations = 200, burnin.iterations = 100, seed = 123)
#' )
#'
#' pooled_obj <- mi_with(
#'   object = fit,
#'   data = linked_df,
#'   formula = survival::Surv(time, status) ~ trt
#' )
#'
#' print(pooled_obj, digits = 4)
#'
#' @export
#' @method print mi_link_pool_survreg
print.mi_link_pool_survreg <- function(x,
                                       digits = max(3L, getOption("digits") - 2L),
                                       ...) {
 cat("Pooled Cox regression results across posterior match classifications:\n")
 cat("  Retained imputations (m):", x$m, "\n")
 cat("  Mixture model distribution:", x$dist, "\n")
 cat("  Refit model: coxph\n\n")

 stats::printCoefmat(.mi_pool_table(x), digits = digits, has.Pvalue = FALSE, ...)
 invisible(x)
}
