#' Bayesian Mixture Helpers (Priors, Engine Glue, Label Alignment)
#'
#' Internal helper functions for the Bayesian mixture workers.
#' These are used by glmMixBayes() / survregMixBayes() and the corresponding
#' fit*.adjMixBayes() dispatch methods.
#'
#' @noRd

################################################################################
# Default priors
################################################################################

fill_defaults <- function(priors = list(), p_family, model_type = 'glm') {
  # Null prior means all defaults need to be filled. Create empty list
  if (is.null(priors)) {
    priors <- list()
  }

  # very simple validation of prior format
  if (!is.list(priors)) {
    stop("Invalid argument: 'priors' must be a named list.")
  }

  # Create list of defaults based on family and model
  # Intercept priors (intercept1/intercept2) are decoupled from slope priors
  # (beta1/beta2) to allow separate regularization.
  if (model_type == 'glm') {
    defaults <- switch(
      p_family,
      "gaussian" = list(
        intercept1 = "normal(0,10)",
        intercept2 = "normal(0,10)",
        beta1 = "normal(0,5)",
        sigma1 = "cauchy(0,2.5)",
        beta2  = "normal(0,5)",
        sigma2 = "cauchy(0,2.5)",
        theta = "beta(1,1)"
      ),
      "poisson" = list(
        intercept1 = "normal(0,10)",
        intercept2 = "normal(0,10)",
        beta1 = "normal(0,5)",
        beta2 = "normal(0,5)",
        theta = "beta(1,1)"
      ),
      "binomial" = list(
        intercept1 = "normal(0,10)",
        intercept2 = "normal(0,10)",
        beta1 = "normal(0,2.5)",
        beta2 = "normal(0,5)",
        theta = "beta(1,1)"
      ),
      "gamma" = list(
        intercept1 = "normal(0,10)",
        intercept2 = "normal(0,10)",
        beta1 = "normal(0,5)",
        beta2 = "normal(0,5)",
        phi1 = "gamma(2,0.1)",
        phi2 = "gamma(2,0.1)",
        theta = "beta(1,1)"
      ),
      stop("`p_family` must be one of 'gaussian','poisson','binomial','gamma'",
           " for model_type == 'glm'")
    )
  } else if (model_type == 'survival') {
    defaults <- switch(
      p_family,
      "gamma" = list(
        # priors for intercept and regression coefficients
        intercept1 = "normal(0,10)",
        intercept2 = "normal(0,10)",
        beta1 = "normal(0,5)",
        beta2 = "normal(0,5)",
        theta = "beta(1,1)",
        # priors for shape parameters
        phi1 = "exponential(1)",
        phi2 = "exponential(1)"
      ),
      "weibull" = list(
        intercept1 = "normal(0,10)",
        intercept2 = "normal(0,10)",
        beta1 = "normal(0,2)",
        beta2 = "normal(0,2)",
        shape1 = "gamma(2,1)",
        shape2 = "gamma(2,1)",
        scale1 = "gamma(2,1)",
        scale2 = "gamma(2,1)",
        theta = "beta(1,1)"
      ),
      stop("`p_family` must be 'gamma' or 'weibull' for model_type 'survival'.")
    )
  } else {
    stop("Unknown model_type in fill_defaults")
  }

  # merge user priors and appropriate defaults
  utils::modifyList(defaults, priors)
}

################################################################################
# Prior strings
################################################################################

# Distributions accepted for each kind of parameter. Strings follow the Stan
# parameterization: normal(mu, sd), student_t(nu, mu, sd), cauchy(loc, scale),
# gamma(shape, rate), exponential(rate), lognormal(meanlog, sdlog),
# inv_gamma(shape, scale), beta(a, b). Distributions with support on the whole
# real line are truncated at zero when used for a positive parameter.
.prior_menu <- list(
 coefficient = "normal",
 positive    = c("normal", "cauchy", "student_t", "gamma", "exponential", "lognormal", "inv_gamma"),
 # the Weibull scale is sampled jointly with the coefficients on the log scale,
 # which requires a prior that is log-concave in log(scale)
 scale       = c("gamma", "exponential", "lognormal"),
 theta       = "beta"
)

# Number of arguments and positivity constraints of each distribution.
.prior_spec <- list(
 normal      = list(n = 2L, positive = 2L, names = c("mu", "sd")),
 cauchy      = list(n = 2L, positive = 2L, names = c("loc", "scale")),
 student_t   = list(n = 3L, positive = c(1L, 3L), names = c("nu", "mu", "sd")),
 gamma       = list(n = 2L, positive = c(1L, 2L), names = c("shape", "rate")),
 exponential = list(n = 1L, positive = 1L, names = "rate"),
 lognormal   = list(n = 2L, positive = 2L, names = c("meanlog", "sdlog")),
 inv_gamma   = list(n = 2L, positive = c(1L, 2L), names = c("shape", "scale")),
 beta        = list(n = 2L, positive = c(1L, 2L), names = c("a", "b"))
)

# Regex parser to safely extract numbers from strings (e.g., "normal(0, 5)")
parse_prior_string <- function(prior_str) {
 if (!is.character(prior_str) || length(prior_str) != 1L) {
  stop("Prior specifications must be single character strings such as \"normal(0, 5)\".", call. = FALSE)
 }
 prior_str <- trimws(prior_str)

 pattern <- "^([a-zA-Z0-9_]+)\\s*\\(([^)]*)\\)$"
 if (!grepl(pattern, prior_str)) {
  stop(sprintf("Prior string '%s' must be in format 'dist(arg1, arg2)'.", prior_str), call. = FALSE)
 }

 match <- regexec(pattern, prior_str)
 parts <- regmatches(prior_str, match)[[1]]

 args_str <- parts[3]
 args_num <- suppressWarnings(as.numeric(strsplit(args_str, ",")[[1]]))

 if (anyNA(args_num)) {
  stop(sprintf("Could not parse numeric arguments from '%s'.", prior_str), call. = FALSE)
 }

 list(dist = tolower(parts[2]), args = args_num)
}

# Parse and validate one prior string against a menu of distributions.
# Returns list(dist, args, written); exponential(rate) is reported as the
# equivalent gamma(1, rate) so that downstream code only needs one code path,
# with written = "exponential" (the distribution name as written otherwise),
# so that fit$priors can keep the spelling.
.parse_prior <- function(key, prior_str, allowed) {
 parsed <- parse_prior_string(prior_str)
 dist <- parsed$dist
 if (!(dist %in% allowed)) {
  stop(sprintf("Prior for '%s' uses '%s', which is not supported here. Supported: %s.",
               key, dist, paste(allowed, collapse = ", ")), call. = FALSE)
 }
 spec <- .prior_spec[[dist]]
 if (length(parsed$args) != spec$n) {
  stop(sprintf("Prior for '%s' expected %d argument%s: %s(%s).", key, spec$n,
               if (spec$n > 1L) "s" else "", dist, paste(spec$names, collapse = ", ")), call. = FALSE)
 }
 if (any(!is.finite(parsed$args))) stop(sprintf("Prior for '%s' has non-finite arguments.", key), call. = FALSE)
 bad <- spec$positive[parsed$args[spec$positive] <= 0]
 if (length(bad) > 0L) {
  stop(sprintf("Prior for '%s': argument '%s' of %s() must be positive.", key,
               spec$names[bad[1L]], dist), call. = FALSE)
 }
 if (dist == "exponential") return(list(dist = "gamma", args = c(1, parsed$args[1L]), written = dist))
 list(dist = dist, args = parsed$args, written = dist)
}

# Exact moments of logit(X) for X ~ Beta(a, b)
logit_beta_moments <- function(a, b) {
 c(mu = digamma(a) - digamma(b), sd = sqrt(trigamma(a) + trigamma(b)))
}

# Convert a beta prior on theta to a normal prior on the logit scale
#
# Used when the user specifies \code{theta = "beta(a, b)"} but the logistic
# regression parameterization is active (Path B).  The conversion uses exact
# moment matching for the logit-beta distribution:
# \itemize{
#   \item \code{mean = digamma(a) - digamma(b)}
#   \item \code{sd   = sqrt(trigamma(a) + trigamma(b))}
# }
#
# @param prior_str Character string of the form \code{"beta(a, b)"}.
# @return Character string of the form \code{"normal(mu, sd)"}.
beta_to_logit_normal <- function(prior_str) {
  parsed <- parse_prior_string(prior_str)
  if (parsed$dist != "beta") {
    stop("beta_to_logit_normal() expects a beta prior, got '", parsed$dist, "'.",
         call. = FALSE)
  }
  m <- logit_beta_moments(parsed$args[1], parsed$args[2])
  sprintf("normal(%s, %s)", round(m[["mu"]], 6), round(m[["sd"]], 6))
}

# Convert (m.rate, sigma_theta) to a beta prior on theta via direct moment matching
#
# Given the expected mismatch rate \code{m.rate} and a desired probability-scale
# standard deviation \code{sigma_theta} of \eqn{\theta}, find Beta(a, b) such that
# \code{E[theta] = 1 - m.rate} and \code{Var[theta] = sigma_theta^2}:
# \deqn{a + b = \frac{\mathrm{m.rate}\,(1-\mathrm{m.rate})}{\sigma_\theta^2} - 1,}
# which requires \eqn{\sigma_\theta^2 < \mathrm{m.rate}\,(1-\mathrm{m.rate})}.
#
# @param m.rate Numeric in (0, 1): expected mismatch rate.
# @param sigma_theta Numeric > 0: desired SD of \eqn{\theta} on the probability scale.
# @return Named list with \code{alpha}, \code{beta}.
mrate_to_beta <- function(m.rate, sigma_theta = 0.1) {
  v <- m.rate * (1 - m.rate)
  if (sigma_theta^2 >= v) {
    stop(sprintf(
      "sigma_theta^2 (= %g) must be strictly less than m.rate*(1-m.rate) (= %g). ",
      sigma_theta^2, v),
      "Choose a smaller sigma_theta or note that this variance exceeds the ",
      "maximum possible for a Beta with mean = 1 - m.rate.",
      call. = FALSE)
  }
  s <- v / sigma_theta^2 - 1   # concentration a + b
  list(alpha = (1 - m.rate) * s, beta = m.rate * s)
}

# Validate the (m.rate, m.rate.sd) pair used to derive the match-rate prior:
# each value on its own always, and whether the pair defines a beta prior
# (m.rate.sd^2 < m.rate * (1 - m.rate)) only when `check_pair` is TRUE, i.e.
# when m.rate is actually used (no explicit theta or, with linkage
# covariates, gamma_intercept prior takes precedence).
.validate_mrate <- function(m.rate, m.rate.sd, check_pair = TRUE) {
 if (!is.null(m.rate)) {
  if (!is.numeric(m.rate) || length(m.rate) != 1L || is.na(m.rate) || m.rate <= 0 || m.rate >= 1) {
   stop("'m.rate' must be a single numeric value strictly between 0 and 1.", call. = FALSE)
  }
  if (is.null(m.rate.sd)) {
   stop("'m.rate.sd' must be a single positive numeric value when 'm.rate' is supplied.", call. = FALSE)
  }
 }
 if (!is.null(m.rate.sd)) {
  if (!is.numeric(m.rate.sd) || length(m.rate.sd) != 1L || is.na(m.rate.sd) || m.rate.sd <= 0) {
   stop("'m.rate.sd' must be a single positive numeric value.", call. = FALSE)
  }
  if (isTRUE(check_pair) && !is.null(m.rate) && m.rate.sd^2 >= m.rate * (1 - m.rate)) {
   stop(sprintf(paste0(
    "With m.rate = %g the prior SD of the match probability must be strictly less than ",
    "sqrt(m.rate * (1 - m.rate)) = %.4g, but m.rate.sd = %g (the default is 0.1). ",
    "Supply a smaller 'm.rate.sd'."),
    m.rate, sqrt(m.rate * (1 - m.rate)), m.rate.sd), call. = FALSE)
  }
 }
 invisible(TRUE)
}

# Bound on m.rate.sd for the beta prior implied by (m.rate, m.rate.sd): at the
# bound the smaller shape parameter equals 1 (for m.rate < 0.5,
# beta((1 - m.rate) / m.rate, 1), whose density still rises to its maximum at
# theta = 1), and below it both shape parameters exceed 1, so that the density
# is bounded (it is unbounded at 0 or 1 when a shape parameter is below 1).
.mrate_sd_unimodal <- function(m.rate) {
 sqrt(min(m.rate^2 * (1 - m.rate) / (1 + m.rate), m.rate * (1 - m.rate)^2 / (2 - m.rate)))
}

# A positive number rounded down to `digits` significant digits, so that a
# bound printed in a message is never above the exact value (format() and
# signif() round to the nearest value, e.g. 0.047559 to 0.0476).
.floor_signif <- function(x, digits = 4L) {
 if (!is.finite(x) || x <= 0) return(x)
 m <- 10^(digits - 1L - floor(log10(x)))
 floor(x * m * (1 + 1e-12)) / m
}

# The bound of .mrate_sd_unimodal() for messages: rounded down to 4
# significant digits.
.mrate_sd_bound_text <- function(m.rate) {
 format(.floor_signif(.mrate_sd_unimodal(m.rate), 4L), digits = 4L, scientific = FALSE)
}

# Which shape parameters of a beta distribution are below 1, allowing for
# rounding error: at m.rate.sd = .mrate_sd_unimodal(m.rate) the smaller shape
# parameter is 1 up to rounding (e.g. 1 - 1e-15), which is not reported as
# below 1.
.shape_below_one <- function(ab) {
 as.numeric(ab) < 1 - 1e-8
}

# Stems of the component-specific prior entries, which need the component
# number (intercept1, intercept2, beta1, ...).
.component_prior_stems <- c("intercept", "beta", "sigma", "phi", "shape", "scale")

# Distribution menu of a prior key (NULL for an unrecognised key, including a
# component-specific stem without its component number, e.g. "intercept").
.key_menu <- function(key) {
 if (key %in% c("gamma_intercept", "gamma_slope", "gamma")) return(.prior_menu$coefficient)
 if (identical(key, "theta")) return(.prior_menu$theta)
 if (!grepl(paste0("^(", paste(.component_prior_stems, collapse = "|"), ")[12]$"), key)) return(NULL)
 switch(sub("[12]$", "", key),
        intercept = , beta = .prior_menu$coefficient,
        sigma = , phi = , shape = .prior_menu$positive,
        scale = .prior_menu$scale)
}

# Hint added to the warning on unrecognised prior entries when some of them
# are component-specific stems without the component number.
.bare_stem_hint <- function(unknown) {
 bare <- intersect(unknown, .component_prior_stems)
 if (length(bare) == 0L) return("")
 others <- bare[-1L]
 sprintf(" Component-specific priors need the component number: write %s1 and/or %s2%s.",
         bare[1L], bare[1L],
         if (length(others) > 0L) sprintf(" (likewise for %s)", paste(others, collapse = ", ")) else "")
}

# Validate a user `priors` list without knowing the family: it must be a
# named list of prior strings, and every recognised entry must use a
# distribution of its menu with valid arguments. Unrecognised names are
# reported with a warning and returned (invisibly). Used by adjMixBayes() so
# that mistakes surface when the adjustment object is created.
.check_prior_list <- function(priors) {
 if (is.null(priors)) return(invisible(NULL))
 if (!is.list(priors) || (length(priors) > 0L && (is.null(names(priors)) || anyNA(names(priors)) ||
                                                   any(!nzchar(names(priors)))))) {
  stop("'priors' must be a named list or NULL.", call. = FALSE)
 }
 unknown <- character(0)
 for (key in names(priors)) {
  menu <- .key_menu(key)
  if (is.null(menu)) { unknown <- c(unknown, key); next }
  .parse_prior(key, priors[[key]], menu)
 }
 if (length(unknown) > 0L) {
  warning("Ignoring unrecognised prior entries: ", paste(unknown, collapse = ", "),
          ". Recognised entries: intercept1, intercept2, beta1, beta2, theta, sigma1, sigma2, ",
          "phi1, phi2, shape1, shape2, scale1, scale2, gamma_intercept, gamma_slope.",
          .bare_stem_hint(unknown), call. = FALSE)
 }
 invisible(unknown)
}

# Keys of the `priors` list that a given model understands.
.prior_keys <- function(family, model_type, use_logistic) {
 keys <- c("intercept1", "intercept2", "beta1", "beta2", "theta")
 if (model_type == "glm" && family == "gaussian") keys <- c(keys, "sigma1", "sigma2")
 if (family == "gamma") keys <- c(keys, "phi1", "phi2")
 if (family == "weibull") keys <- c(keys, "shape1", "shape2", "scale1", "scale2")
 if (use_logistic) keys <- c(keys, "gamma_intercept", "gamma_slope", "gamma")
 keys
}

# Warning issued when the beta prior implied by (m.rate, m.rate.sd) has a
# shape parameter below 1 (`theta_ab`, from mrate_to_beta()). Without linkage
# covariates this beta is the prior on theta: its density is unbounded at 0 or
# 1. With linkage covariates (`use_logistic`) the intercept of the
# match-probability model receives a normal prior with the logit-scale moments
# of this beta, which the message describes on the probability scale (median,
# central 95% interval and standard deviation of the match probability of a
# record with average linkage covariates, which can differ markedly from
# m.rate.sd (it is much larger for small mismatch rates), and the
# standard deviation at the bound on m.rate.sd). The bound on m.rate.sd is
# rounded down (see .mrate_sd_bound_text()).
.mrate_dispersed_message <- function(m.rate, m.rate.sd, theta_ab, use_logistic) {
 bound <- .mrate_sd_bound_text(m.rate)
 num <- function(v) format(signif(v, 3))
 theta_ab <- unname(as.numeric(theta_ab))
 if (isTRUE(use_logistic)) {
  mom <- logit_beta_moments(theta_ab[1], theta_ab[2])
  q <- stats::plogis(mom[["mu"]] + c(0, -1, 1) * stats::qnorm(0.975) * mom[["sd"]])
  psd <- .logitnormal_moments(mom[["mu"]], mom[["sd"]])[["sd"]]
  return(sprintf(paste0(
   "m.rate = %g with m.rate.sd = %g gives the intercept of the match-probability model the prior ",
   "normal(%s, %s) on the logit scale, so that the match probability of a record with average ",
   "linkage covariates has prior median %s, 95%% interval (%s, %s) and standard deviation %s ",
   "(rather than m.rate.sd = %g). With m.rate.sd = %s this standard deviation would be %s; ",
   "smaller values give a more concentrated prior."),
   m.rate, m.rate.sd, num(mom[["mu"]]), num(mom[["sd"]]), num(q[1]), num(q[2]), num(q[3]), num(psd),
   m.rate.sd, bound, num(.mrate_logit_sd(m.rate, as.numeric(bound)))))
 }
 low <- .shape_below_one(theta_ab)
 where <- if (all(low)) "0 and 1" else if (low[2]) "1" else "0"
 sprintf(paste0(
  "m.rate = %g with m.rate.sd = %g implies a beta(%s, %s) prior on the match probability, whose ",
  "density is unbounded at %s because a shape parameter is below 1. With m.rate.sd at most %s both ",
  "shape parameters exceed 1 and the density is bounded."),
  m.rate, m.rate.sd, num(theta_ab[1]), num(theta_ab[2]), where, bound)
}

# Probability-scale standard deviation of the normal prior on the logit
# scale that (m.rate, m.rate.sd) give the intercept of the match-probability
# model (with linkage covariates).
.mrate_logit_sd <- function(m.rate, m.rate.sd) {
 bp <- mrate_to_beta(m.rate, m.rate.sd)
 mom <- logit_beta_moments(bp$alpha, bp$beta)
 .logitnormal_moments(mom[["mu"]], mom[["sd"]])[["sd"]]
}

# Mean and standard deviation of plogis(Z) for Z ~ normal(mu, sd), the
# probability-scale moments of a normal prior on the logit scale (numerical
# integration over mu +- 12 sd).
.logitnormal_moments <- function(mu, sd) {
 mom <- function(k) {
  stats::integrate(function(z) stats::plogis(z)^k * stats::dnorm(z, mu, sd),
                   mu - 12 * sd, mu + 12 * sd, rel.tol = 1e-8)$value
 }
 m1 <- mom(1)
 c(mean = m1, sd = sqrt(max(mom(2) - m1^2, 0)))
}

# Parse the `priors` list into a flat list of numeric hyperparameters.
#
# @param priors Named list of prior strings (may be NULL); missing entries are
#   filled by fill_defaults().
# @param family Likelihood family.
# @param model_type "glm" or "survival".
# @param use_logistic TRUE for the logistic match-probability model (Path B).
# @param m.rate,m.rate.sd Optional expected mismatch rate and prior SD of theta
#   on the probability scale. They inform the theta prior (Path A: Beta by
#   moment matching) or the intercept of gamma (Path B: exact logit moments of
#   the same Beta) when the user gives no explicit prior.
prepare_mixbayes_priors <- function(priors, family, model_type,
                                    use_logistic = FALSE, m.rate = NULL,
                                    m.rate.sd = 0.1) {
 if (!is.null(priors) && !is.list(priors)) stop("'priors' must be a named list or NULL.", call. = FALSE)
 raw <- if (is.null(priors)) list() else priors
 if (length(raw) > 0L && (is.null(names(raw)) || any(!nzchar(names(raw))))) {
  stop("'priors' must be a named list.", call. = FALSE)
 }
 known <- .prior_keys(family, model_type, use_logistic = TRUE)
 unknown <- setdiff(names(raw), known)
 if (length(unknown) > 0L) {
  warning("Ignoring unrecognised prior entries: ", paste(unknown, collapse = ", "),
          ". Recognised entries for this model: ",
          paste(.prior_keys(family, model_type, use_logistic), collapse = ", "), ".",
          .bare_stem_hint(unknown), call. = FALSE)
 }
 gamma_keys <- intersect(names(raw), c("gamma_intercept", "gamma_slope", "gamma"))
 if (!use_logistic && length(gamma_keys) > 0L) {
  warning("Priors ", paste(gamma_keys, collapse = ", "),
          " apply to the logistic match-probability model only and are ignored ",
          "because no linkage covariates were supplied (Path A).", call. = FALSE)
 }
 user_theta <- !is.null(raw[["theta"]])
 user_gi <- use_logistic && !is.null(raw[["gamma_intercept"]])
 # m.rate is used only without an explicit theta (or, with linkage
 # covariates, gamma_intercept) prior; the pair (m.rate, m.rate.sd) must then
 # define a beta prior
 mrate_used <- !is.null(m.rate) && !user_theta && !user_gi
 .validate_mrate(m.rate, m.rate.sd, check_pair = mrate_used)
 if (!is.null(m.rate) && (user_theta || user_gi)) {
  message("'m.rate' is not used because an explicit '",
          if (user_gi) "gamma_intercept" else "theta", "' prior was supplied.")
 }
 if (user_theta && user_gi) {
  message("The 'theta' prior is not used: with linkage covariates the explicit 'gamma_intercept' ",
          "prior applies to the intercept of the match-probability model.")
 }
 if (use_logistic && !is.null(raw[["gamma_slope"]]) && !is.null(raw[["gamma"]])) {
  message("The 'gamma' prior (an alias of 'gamma_slope') is not used because a 'gamma_slope' ",
          "prior was supplied.")
 }

 priors <- fill_defaults(priors, family, model_type)
 coef_prior <- function(key) .parse_prior(key, priors[[key]], .prior_menu$coefficient)$args

 # Core outcome-model priors (shared between Path A and Path B)
 i1 <- coef_prior("intercept1"); i2 <- coef_prior("intercept2")
 b1 <- coef_prior("beta1"); b2 <- coef_prior("beta2")
 out <- list(
  prior_intercept1_mu = i1[1], prior_intercept1_sd = i1[2],
  prior_intercept2_mu = i2[1], prior_intercept2_sd = i2[2],
  prior_beta1_mu = b1[1], prior_beta1_sd = b1[2],
  prior_beta2_mu = b2[1], prior_beta2_sd = b2[2]
 )

 # Beta prior on theta: explicit user prior > (m.rate, m.rate.sd) > default
 theta_ab <- .parse_prior("theta", priors[["theta"]], .prior_menu$theta)$args
 if (mrate_used) {
  bp <- mrate_to_beta(m.rate, m.rate.sd)
  theta_ab <- c(bp$alpha, bp$beta)
  if (any(.shape_below_one(theta_ab))) {
   warning(.mrate_dispersed_message(m.rate, m.rate.sd, theta_ab, use_logistic), call. = FALSE)
  }
 }

 if (use_logistic) {
  # Path B: gamma = (intercept, slopes) of the logistic match model.
  # Intercept: explicit gamma_intercept > logit-scale moments of the theta prior
  # (user theta, or the Beta derived from m.rate/m.rate.sd, or the default).
  if (!is.null(raw[["gamma_intercept"]])) {
   gi <- .parse_prior("gamma_intercept", raw[["gamma_intercept"]], .prior_menu$coefficient)$args
  } else {
   gi <- unname(logit_beta_moments(theta_ab[1], theta_ab[2]))
  }
  # Slopes: explicit gamma_slope > legacy 'gamma' > normal(0, 2.5)
  if (!is.null(raw[["gamma_slope"]])) {
   gs <- .parse_prior("gamma_slope", raw[["gamma_slope"]], .prior_menu$coefficient)$args
  } else if (!is.null(raw[["gamma"]])) {
   gs <- .parse_prior("gamma", raw[["gamma"]], .prior_menu$coefficient)$args
  } else {
   gs <- c(0, 2.5)
  }
  out$prior_gamma_intercept_mu <- gi[1]; out$prior_gamma_intercept_sd <- gi[2]
  out$prior_gamma_slope_mu <- gs[1]; out$prior_gamma_slope_sd <- gs[2]
 } else {
  out$prior_theta_alpha <- theta_ab[1]
  out$prior_theta_beta  <- theta_ab[2]
 }

 # Family-specific dispersion / shape priors: distribution + arguments (an
 # exponential prior is used as gamma(1, rate); "written" keeps its name for
 # fit$priors)
 add_scalar <- function(key, allowed = .prior_menu$positive) {
  p <- .parse_prior(key, priors[[key]], allowed)
  out[[paste0("prior_", key, "_dist")]] <<- p$dist
  out[[paste0("prior_", key, "_args")]] <<- p$args
  if (!identical(p$written, p$dist)) out[[paste0("prior_", key, "_written")]] <<- p$written
 }
 if (model_type == "glm" && family == "gaussian") {
  add_scalar("sigma1"); add_scalar("sigma2")
 } else if (family == "gamma") {
  add_scalar("phi1"); add_scalar("phi2")
 } else if (family == "weibull") {
  add_scalar("shape1"); add_scalar("shape2")
  add_scalar("scale1", allowed = .prior_menu$scale); add_scalar("scale2", allowed = .prior_menu$scale)
 }

 out
}

################################################################################
# Engine glue: translate the hyperparameters into the vectors used by the C++
# Gibbs sampler, run it with a reproducible seed, and align labels.
################################################################################

# Numeric code of a positive-scalar prior distribution understood by the sampler.
.scalar_prior_code <- function(dist) {
 code <- c(normal = 1, cauchy = 2, gamma = 3, exponential = 4, lognormal = 5,
           student_t = 6, inv_gamma = 7)[dist]
 if (is.na(code)) {
  stop(sprintf(paste0("Prior distribution '%s' is not supported for a positive scalar ",
                      "parameter. Use one of: %s."), dist,
               paste(.prior_menu$positive, collapse = ", ")), call. = FALSE)
 }
 unname(code)
}

# Build the prior list consumed by mixbayes_gibbs_cpp().
#
# @param flat Output of prepare_mixbayes_priors().
# @param family Likelihood family ("gaussian", "poisson", "binomial", "gamma",
#   "weibull").
# @param model_type "glm" or "survival".
# @param K Number of columns of X.
# @param M Number of columns of Z (0 in Path A; the first column, a column of
#   ones, receives the gamma intercept prior).
# @param intercept Whether the first column of X is an intercept (a column of
#   ones). It then receives the intercept prior and the other columns the
#   slope prior; otherwise every column receives the slope prior.
build_engine_priors <- function(flat, family, model_type, K, M = 0L, intercept = TRUE) {
 expand <- function(int_val, slope_val) {
  if (isTRUE(intercept)) c(int_val, rep(slope_val, max(K - 1L, 0L))) else rep(slope_val, K)
 }
 out <- list(
  beta1_mean = expand(flat$prior_intercept1_mu, flat$prior_beta1_mu),
  beta1_sd   = expand(flat$prior_intercept1_sd, flat$prior_beta1_sd),
  beta2_mean = expand(flat$prior_intercept2_mu, flat$prior_beta2_mu),
  beta2_sd   = expand(flat$prior_intercept2_sd, flat$prior_beta2_sd)
 )
 if (M > 0L) {
  out$gamma_mean <- c(flat$prior_gamma_intercept_mu, rep(flat$prior_gamma_slope_mu, M - 1L))
  out$gamma_sd   <- c(flat$prior_gamma_intercept_sd, rep(flat$prior_gamma_slope_sd, M - 1L))
 } else {
  out$theta <- c(flat$prior_theta_alpha, flat$prior_theta_beta)
 }
 for (nm in grep("_sd$", names(out), value = TRUE)) {
  if (!all(is.finite(out[[nm]])) || any(out[[nm]] <= 0)) {
   stop("Prior standard deviations must be positive (", nm, ").", call. = FALSE)
  }
 }
 if (!is.null(out$theta) && (!all(is.finite(out$theta)) || any(out$theta <= 0))) {
  stop("The beta prior on theta must have positive shape parameters.", call. = FALSE)
 }

 scalar <- function(key) {
  args <- flat[[paste0("prior_", key, "_args")]]
  c(.scalar_prior_code(flat[[paste0("prior_", key, "_dist")]]), args, rep(0, 3L - length(args)))
 }
 if (model_type == "glm" && family == "gaussian") {
  out$disp1 <- scalar("sigma1"); out$disp2 <- scalar("sigma2")
 } else if (family == "gamma") {
  out$disp1 <- scalar("phi1"); out$disp2 <- scalar("phi2")
 } else if (family == "weibull") {
  out$disp1 <- scalar("shape1"); out$disp2 <- scalar("shape2")
  out$scale1 <- scalar("scale1"); out$scale2 <- scalar("scale2")
 }
 # the sampler uses the intercept column for the starting values of the
 # mismatch component
 out$intercept <- isTRUE(intercept)
 out
}

# TRUE when the first column of a design matrix is an intercept (all ones).
.has_intercept_column <- function(X) {
 ncol(X) >= 1L && nrow(X) >= 1L && all(X[, 1L] == 1)
}

# Without an intercept column in the design matrix every column receives the
# beta1 / beta2 prior (see build_engine_priors()), so intercept1 and
# intercept2 entries of `priors` are not used: warn when they were supplied
# and drop them, so that the fit neither reports them in fit$priors nor
# treats them as component-specific priors. With the design matrix `X`, the
# warning names a column of ones that is not first (named "(Intercept)" when
# it had no name, see .normalize_design_names()), which also receives the
# beta1 / beta2 prior.
.drop_unused_intercept_priors <- function(priors, intercept, X = NULL) {
 if (isTRUE(intercept) || !is.list(priors) || length(priors) == 0L) return(priors)
 given <- intersect(c("intercept1", "intercept2"), names(priors))
 given <- given[!vapply(given, function(k) is.null(priors[[k]]), logical(1))]
 if (length(given) > 0L) {
  j <- if (is.matrix(X)) .intercept_column(X) else 0L
  why <- if (j > 0L) {
   nm <- colnames(X)[j]
   sprintf(paste0("the intercept prior applies only to a first column of ones, and the column of ones of ",
                  "the design matrix (%s) is not first, so every column, this one included, receives the ",
                  "beta1 / beta2 prior. Put the column of ones first to give it the intercept prior."),
           if (!is.null(nm) && !is.na(nm) && nzchar(nm)) sprintf("column %d, \"%s\"", j, nm) else
            sprintf("column %d", j))
  } else {
   "the design matrix has no intercept column (a first column of ones), so every column receives the beta1 / beta2 prior."
  }
  warning(sprintf("The %s %s not used: %s",
                  paste(given, collapse = " and "), if (length(given) > 1L) "priors are" else "prior is", why),
          call. = FALSE)
  priors[given] <- NULL
 }
 priors
}

# Column names of a design matrix of the engines (X, or Z with
# prefix = "Z") with every name filled in, so that the posterior draws, the
# messages and predict() can refer to the columns by name: an unnamed (empty
# or NA) column of ones becomes "(Intercept)" (the first such column, unless
# another column already has this name), any other unnamed column
# "<prefix><j>", j being its position (e.g. "X2" for X = cbind(1, x)); a
# filled-in name that is already taken gets a suffix ".1", ".2", ... Names
# that are given are kept, also when duplicated.
.normalize_design_names <- function(M, prefix = "X") {
 K <- ncol(M)
 cn <- colnames(M)
 if (is.null(cn)) cn <- rep("", K)
 unnamed <- is.na(cn) | !nzchar(cn)
 if (!any(unnamed)) return(M)
 ones <- if (nrow(M) > 0L) colSums(M != 1) == 0 else rep(FALSE, K)
 j <- which(unnamed & ones)
 if (length(j) > 0L && !"(Intercept)" %in% cn[!unnamed]) {
  cn[j[1L]] <- "(Intercept)"
  unnamed[j[1L]] <- FALSE
 }
 for (k in which(unnamed)) {
  nm <- base <- paste0(prefix, k)
  i <- 0L
  while (nm %in% cn[-k]) {
   i <- i + 1L
   nm <- paste0(base, ".", i)
  }
  cn[k] <- nm
 }
 colnames(M) <- cn
 M
}

# Effective prior specification as prior strings (for storage in fitted
# objects and for printing), from the output of prepare_mixbayes_priors().
# Without an intercept column in the design matrix (`intercept = FALSE`)
# every column receives the beta1 / beta2 prior, so intercept1 and intercept2
# are left out. An exponential prior keeps its name (it is used as the same
# distribution gamma(1, rate)).
.prior_strings <- function(flat, family, model_type, use_logistic, intercept = TRUE) {
 num <- function(v) vapply(v, function(x) format(signif(x, 4), trim = TRUE), "")
 str <- function(dist, args) paste0(dist, "(", paste(num(args), collapse = ", "), ")")
 nrm <- function(key) str("normal", c(flat[[paste0("prior_", key, "_mu")]], flat[[paste0("prior_", key, "_sd")]]))
 out <- c(intercept1 = nrm("intercept1"), beta1 = nrm("beta1"),
          intercept2 = nrm("intercept2"), beta2 = nrm("beta2"))
 if (!isTRUE(intercept)) out <- out[c("beta1", "beta2")]
 scal <- switch(if (model_type == "glm" && family == "gaussian") "gaussian" else family,
                gaussian = c("sigma1", "sigma2"), gamma = c("phi1", "phi2"),
                weibull = c("shape1", "shape2", "scale1", "scale2"), character(0))
 for (k in scal) {
  args <- flat[[paste0("prior_", k, "_args")]]
  out[k] <- if (identical(flat[[paste0("prior_", k, "_written")]], "exponential")) {
   str("exponential", args[2L])
  } else {
   str(flat[[paste0("prior_", k, "_dist")]], args)
  }
 }
 if (use_logistic) {
  out["gamma_intercept"] <- nrm("gamma_intercept")
  out["gamma_slope"] <- nrm("gamma_slope")
 } else {
  out["theta"] <- str("beta", c(flat$prior_theta_alpha, flat$prior_theta_beta))
 }
 out
}

# Starting state of the main chain as reported by the sampler, for the
# warnings on the label orientation (see align_mixture_labels()): the share of
# the records not flagged as safe matches allocated to component 1 before any
# orientation (on average over the last sweeps of the selected pilot chain,
# else in the initial allocation), whether the labels were exchanged
# (`pre_oriented`) or the
# exchange was refused because the component-specific priors make the
# labelling with component 1 as the majority much less probable (`refused`,
# with the estimated log ratio `lp_change` and the result of the check chains
# `check`), whether the user supplied starting values (`init`), which are never
# reoriented, and whether the component-specific priors are the defaults of
# the family (`default_priors`, see .default_component_priors()).
.orientation_start <- function(posterior, mcmc, default_priors = FALSE) {
 list(share1 = posterior$start_share1,
      pre_oriented = isTRUE(posterior$pre_oriented),
      init = length(mcmc$init) > 0L,
      refused = isTRUE(posterior$orient_refused),
      lp_change = posterior$orient_lp_change,
      check = posterior$orient_check,
      default_priors = isTRUE(default_priors))
}

# Whether every component-specific prior (intercepts, slopes, sigma / phi /
# shape, Weibull scale) is the default of the family, from the output of
# prepare_mixbayes_priors(); used to word the warning issued when the
# orientation of the starting state was refused (only the "binomial" defaults
# differ between the components).
.default_component_priors <- function(prior_flat, family, model_type) {
 def <- prepare_mixbayes_priors(NULL, family, model_type)
 # the distributions as used (an exponential prior is gamma(1, rate), however written)
 keys <- grep("^prior_(intercept|beta|sigma|phi|shape|scale)[12]_(mu|sd|dist|args)$", names(def), value = TRUE)
 all(vapply(keys, function(k) identical(as.vector(def[[k]]), as.vector(prior_flat[[k]])), logical(1)))
}

# Backstop for the orientation of the starting state under the majority
# convention with component-specific priors that differ (`exchangeable`
# FALSE; for exchangeable components the exchange leaves the posterior
# unchanged): the sampler exchanges the labels only when the exchanged
# labelling is not much less probable, so the main chain is expected to stay
# near the log posterior of the selected pilot chain, plus the change of level
# when negative: the difference between the check chains started with the
# labels exchanged and kept when they were run (see mixbayes_gibbs_cpp()),
# otherwise the estimated log ratio of the two labellings. Warn when the
# average log posterior of the stored draws as sampled (`posterior$lp`, before
# any relabelling) lies below that level by more than max(`min_margin`, 2 SD
# of the log posterior of the stored draws): the chain is then probably
# trapped in a minor mode. The SD term allows for the Monte Carlo noise of
# both levels in large models, whose log posterior varies (and drifts) by
# several units. Under the match-rate prior (`orientation`) the sampler
# exchanges the labels only towards the more probable labelling, which is not
# checked. Returns TRUE when it warned.
.check_orientation_lp <- function(posterior, exchangeable, orientation = "majority component",
                                  min_margin = 5) {
 if (!isTRUE(posterior$pre_oriented) || isTRUE(exchangeable) ||
     !identical(orientation, "majority component")) {
  return(invisible(FALSE))
 }
 pilot <- posterior$pilot_lp[is.finite(posterior$pilot_lp)]
 lp <- posterior$lp[is.finite(posterior$lp)]
 if (length(pilot) == 0L || length(lp) == 0L) return(invisible(FALSE))
 chk <- posterior$orient_check
 change <- if (all(c("exchanged", "kept") %in% names(chk)) && all(is.finite(chk[c("exchanged", "kept")]))) {
  chk[["exchanged"]] - chk[["kept"]]
 } else {
  posterior$orient_lp_change
 }
 expected <- max(pilot) + if (length(change) == 1L && is.finite(change)) min(0, change) else 0
 margin <- max(min_margin, if (length(lp) > 1L) 2 * stats::sd(lp) else 0)
 if (mean(lp) < expected - margin) {
  warning(sprintf(paste0(
   "The sampler exchanged the labels of its starting state so that component 1 is the majority ",
   "component, but the average log posterior of the stored draws (%.1f) lies %.1f below that ",
   "of the selected pilot chain (%.1f), more than the margin of %.1f: the chain is probably ",
   "trapped in a minor mode of the labelling with component 1 as the majority. Give the ",
   "expected mismatch rate with 'm.rate', supply safe matches or starting values ",
   "(control$init, which are never reoriented), or check the component-specific priors."),
   mean(lp), max(pilot) - mean(lp), max(pilot), margin), call. = FALSE)
  return(invisible(TRUE))
 }
 invisible(FALSE)
}

# TRUE when a one-sided formula has no covariates (e.g. `~ 1`), i.e. the match
# probability is constant (Path A). Comparing formulas with identical() is not
# reliable because it also compares their environments.
.is_intercept_only_formula <- function(f) {
 if (is.null(f)) return(TRUE)
 length(attr(stats::terms(f), "term.labels")) == 0L
}

# Refuse the modelling arguments that the Bayesian engines do not support
# (weights, offset, ...), given their names.
.refuse_modelling_args <- function(arg_names) {
 modelling <- intersect(arg_names, c("weights", "offset", "start", "etastart", "mustart"))
 if (length(modelling) > 0L) {
  stop("The Bayesian mixture models do not support the argument(s) ",
       paste(modelling, collapse = ", "), ".", call. = FALSE)
 }
 invisible(NULL)
}

# Names of the arguments in `...` without evaluating them (base::...names()
# needs R >= 4.1), so that an unsupported argument whose value cannot be
# evaluated, e.g. offset = log(e) with `e` a column of the linked data, is
# refused by name instead of failing with "object 'e' not found".
.dots_names <- function(...) {
 nm <- names(substitute(list(...)))[-1L]
 if (is.null(nm)) character(0) else nm
}

# Resolve an MCMC control value: `...` overrides `control`, which overrides the
# default.
.mixbayes_control <- function(name, dots, control, default) {
 if (name %in% names(dots)) return(dots[[name]])
 if (is.list(control) && name %in% names(control)) return(control[[name]])
 default
}

# Collect the MCMC settings of glmMixBayes()/survregMixBayes() from `control`
# and `...`, validate them and return them normalised (integers as integer).
# Unknown settings are reported (the Stan-based engine accepted `cores`,
# `adapt_delta` and `max_treedepth`; only `cores` is still tolerated, silently).
# Modelling arguments that the Bayesian engines do not support (weights,
# offset, ...) are refused instead of being ignored.
.mixbayes_settings <- function(dcontrols, control) {
 if (!is.null(control) && !is.list(control)) {
  stop("`control` must be a named list of MCMC settings.", call. = FALSE)
 }
 if (length(control) > 0L && (is.null(names(control)) || anyNA(names(control)) ||
                              any(!nzchar(names(control))))) {
  stop("`control` must be a named list of MCMC settings (e.g. list(iterations = 4000)).", call. = FALSE)
 }
 given <- unique(c(names(dcontrols), names(control)))
 given <- given[!is.na(given) & nzchar(given)]
 .refuse_modelling_args(given)
 if ("data" %in% given) {
  warning("`data` is ignored: the Bayesian mixture models take the data from the adjustment ",
          "object (adjMixBayes(linked.data = ...)).", call. = FALSE)
 }
 known <- c("iterations", "burnin.iterations", "thin", "seed", "init", "pilots",
            "pilot.iterations", "collapse", "verbose", "cores", "priors")
 unknown <- setdiff(given, c(known, "data"))
 if (length(unknown) > 0L) {
  warning("Ignoring unknown MCMC control setting(s): ", paste(unknown, collapse = ", "),
          ". Supported: ", paste(setdiff(known, c("cores", "priors")), collapse = ", "), ".",
          call. = FALSE)
 }
 iterations <- .mixbayes_control("iterations", dcontrols, control, 1e4)
 # without an explicit burn-in, half of the iterations are discarded, at most
 # 1000, so that the number of stored draws grows with `iterations`
 burnin_default <- if (is.numeric(iterations) && length(iterations) == 1L && is.finite(iterations)) {
  min(1e3, floor(iterations / 2))
 } else {
  1e3
 }
 .normalize_mcmc_settings(list(
  iterations        = iterations,
  burnin.iterations = .mixbayes_control("burnin.iterations", dcontrols, control, burnin_default),
  thin              = .mixbayes_control("thin", dcontrols, control, 1L),
  seed              = .mixbayes_control("seed", dcontrols, control, sample.int(.Machine$integer.max, 1)),
  init              = .mixbayes_control("init", dcontrols, control, NULL),
  pilots            = .mixbayes_control("pilots", dcontrols, control, 5L),
  pilot.iterations  = .mixbayes_control("pilot.iterations", dcontrols, control, 200L),
  collapse          = .mixbayes_control("collapse", dcontrols, control, TRUE),
  verbose           = .mixbayes_control("verbose", dcontrols, control, FALSE)
 ))
}

# Validate the MCMC settings and store integer settings as integers. When
# starting values are supplied no pilot chains are run, which is recorded as
# pilots = 0.
.normalize_mcmc_settings <- function(s) {
 int_setting <- function(name, min) {
  v <- s[[name]]
  if (is.null(v) || !is.numeric(v) || length(v) != 1L || !is.finite(v) || v < min ||
      v != round(v) || v > .Machine$integer.max) {
   stop(sprintf("`%s` must be a single integer >= %d.", name, min), call. = FALSE)
  }
  as.integer(v)
 }
 s$iterations <- int_setting("iterations", 2L)
 s$burnin.iterations <- int_setting("burnin.iterations", 0L)
 if (s$burnin.iterations >= s$iterations) {
  stop("`burnin.iterations` must be smaller than `iterations`.", call. = FALSE)
 }
 s$thin <- int_setting("thin", 1L)
 if ((s$iterations - s$burnin.iterations) %/% s$thin < 1L) {
  stop("No draws would be stored: increase `iterations` or reduce `thin`.", call. = FALSE)
 }
 s$pilots <- int_setting("pilots", 0L)
 if (length(s$init) > 0L) s$pilots <- 0L
 s$pilot.iterations <- if (s$pilots > 0L) int_setting("pilot.iterations", 1L) else 0L
 for (nm in c("collapse", "verbose")) {
  if (is.null(s[[nm]])) s[[nm]] <- identical(nm, "collapse")
  if (!(isTRUE(s[[nm]]) || isFALSE(s[[nm]]))) {
   stop(sprintf("`%s` must be TRUE or FALSE.", nm), call. = FALSE)
  }
 }
 seed <- s$seed
 if (!is.null(seed)) {
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) || seed != round(seed) ||
      abs(seed) > .Machine$integer.max) {
   stop("`seed` must be a single integer.", call. = FALSE)
  }
  s$seed <- as.integer(seed)
 }
 s
}

# Evaluate `expr` with R's RNG seeded by `seed`, restoring the previous RNG
# state afterwards so that model fitting does not disturb the caller's stream.
.with_rng_seed <- function(seed, expr) {
 if (!is.null(seed)) {
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = globalenv(), inherits = FALSE) else NULL
  on.exit({
   if (had_seed) {
    assign(".Random.seed", old_seed, envir = globalenv())
   } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    rm(".Random.seed", envir = globalenv())
   }
  }, add = TRUE)
  set.seed(seed)
 }
 expr
}

# Linkage information carried by an adjMixBayes object, aligned with the rows
# of the outcome model through `idx_map` (row positions in `full_data`):
# the design matrix Z of the match-probability model (NULL for a constant
# match probability), the 0/1 indicator of known correct matches (NULL when
# none were given) and the prior mismatch rate with its probability-scale SD
# (NULL when no mismatch rate was given). Records whose linkage covariates are
# missing are dropped with a warning (as in the frequentist adjMixture()
# framework); `keep` flags the records of `idx_map` that remain.
.mixbayes_link_inputs <- function(adjustment, full_data, idx_map) {
 m_formula <- adjustment$m.formula
 if (is.null(m_formula)) m_formula <- ~1
 keep <- rep(TRUE, length(idx_map))

 Z <- NULL
 if (!.is_intercept_only_formula(m_formula)) {
  tt <- stats::terms(m_formula)
  if (attr(tt, "intercept") == 0L) {
   stop("'m.formula' must include an intercept (the baseline match probability); ",
        "remove '- 1' / '+ 0'.", call. = FALSE)
  }
  if (!is.null(attr(tt, "offset"))) {
   stop("offset() terms are not supported in 'm.formula'.", call. = FALSE)
  }
  data_subset <- full_data[idx_map, , drop = FALSE]
  Z <- tryCatch({
   Z_frame <- stats::model.frame(m_formula, data = data_subset, na.action = stats::na.pass,
                                 drop.unused.levels = TRUE)
   stats::model.matrix(m_formula, Z_frame)
  }, error = function(e) {
   stop("Could not build the match-probability design matrix from 'm.formula' (",
        deparse(m_formula), "): ", conditionMessage(e), call. = FALSE)
  })
  if (any(is.infinite(Z))) {
   stop("Infinite values found in the linkage covariates of 'm.formula' ",
        "(check transformations such as log() of zero).", call. = FALSE)
  }
  keep <- stats::complete.cases(Z)
  if (!any(keep)) {
   stop("No record has complete linkage covariates for 'm.formula'.", call. = FALSE)
  }
  if (!all(keep)) {
   warning(sprintf(paste0("Dropped %d observation(s) due to missing values in the linkage ",
                          "covariates of 'm.formula'."), sum(!keep)), call. = FALSE)
   Z <- Z[keep, , drop = FALSE]
  }
 }

 safe <- NULL
 if (!is.null(adjustment$safe.matches)) {
  if (length(adjustment$safe.matches) != nrow(full_data)) {
   stop("'safe.matches' must have the same length as 'nrow(linked.data)' (",
        nrow(full_data), ").", call. = FALSE)
  }
  safe <- as.integer(adjustment$safe.matches[idx_map][keep])
 }

 m_rate_sd <- NULL
 if (!is.null(adjustment$m.rate)) {
  m_rate_sd <- adjustment$m.rate.sd
  if (is.null(m_rate_sd)) m_rate_sd <- 0.1
 }

 list(Z = Z, safe = safe, m.rate = adjustment$m.rate, m.rate.sd = m_rate_sd, keep = keep)
}

# Shared front end of fitglm.adjMixBayes() and fitsurvreg.adjMixBayes():
# check the linked data, map the rows of the outcome model to the linked data
# through the row names, resolve the prior list (function argument > `...` >
# control$priors > adjustment$priors) and build the linkage inputs. Rows
# whose linkage covariates are missing are removed from `x` and `y`.
.mixbayes_fit_inputs <- function(x, y, adjustment, control, priors, dots) {
 full_data <- adjustment$data_ref$data
 if (is.null(full_data)) {
  stop("The 'adjustment' object does not contain linked data. ",
       "Please recreate the object with 'linked.data' provided.", call. = FALSE)
 }
 if (nrow(x) == 0L) {
  stop("No observations are left after applying 'subset' and 'na.action'.", call. = FALSE)
 }

 subset_names <- rownames(x)
 if (is.null(subset_names)) {
  if (nrow(x) != nrow(full_data)) {
   stop("Row mismatch: Model matrix 'x' has no row names and its length (", nrow(x),
        ") differs from the adjustment data (", nrow(full_data), "). ",
        "Ensure 'linked.data' matches the data passed to the upstream wrapper.", call. = FALSE)
  }
  idx_map <- seq_len(nrow(full_data))
 } else {
  idx_map <- match(subset_names, rownames(full_data))
  if (anyNA(idx_map)) {
   stop("Row mismatch: Some observations in the model matrix could not be matched ",
        "to the adjustment data. This usually happens if the upstream 'data' ",
        "differs from the data used to create the adjustment object.", call. = FALSE)
  }
 }

 if (anyNA(x) || anyNA(y)) {
  stop("NA values found in x or y. Upstream wrapper should remove missingness.", call. = FALSE)
 }

 # the linkage information comes from the adjustment object
 linkage <- intersect(names(dots), c("m.formula", "m.rate", "m.rate.sd", "safe.matches", "safe", "Z"))
 if (length(linkage) > 0L) {
  stop("Linkage information (", paste(linkage, collapse = ", "), ") is given in adjMixBayes(), ",
       "not in plglm() or plsurvreg().", call. = FALSE)
 }

 final_priors <- priors
 if (is.null(final_priors) && "priors" %in% names(dots)) final_priors <- dots$priors
 dots$priors <- NULL
 if (is.null(final_priors) && is.list(control) && "priors" %in% names(control)) {
  final_priors <- control$priors
 }
 if (is.null(final_priors)) {
  # entries of adjustment$priors that adjMixBayes() reported as not
  # recognised when it created the object (e.g. "intercept" without its
  # component number; stored in the attribute "unrecognised", see
  # .check_prior_list()) are dropped, so that the fit does not repeat the
  # warning. Other unrecognised entries (added to the object after its
  # creation) are reported when the priors are parsed.
  final_priors <- adjustment$priors
  reported <- attr(final_priors, "unrecognised")
  if (is.list(final_priors) && length(reported) > 0L && !is.null(names(final_priors))) {
   final_priors <- final_priors[!names(final_priors) %in% reported]
  }
 }

 link <- .mixbayes_link_inputs(adjustment, full_data, idx_map)
 if (!all(link$keep)) {
  x <- x[link$keep, , drop = FALSE]
  y <- if (is.matrix(y)) y[link$keep, , drop = FALSE] else y[link$keep]
 }
 obs_names <- rownames(x)

 list(x = x, y = y, priors = final_priors, link = link, dots = dots, obs_names = obs_names)
}

# Align the `data` argument of mi_with() with the records used in the fit:
# rows are matched by row name when the fit stores the names of the analysed
# records (fits from plglm()/plsurvreg()) and `data` contains all of them;
# otherwise `data` must hold exactly the analysed records in model order.
.mixbayes_mi_data <- function(object, data) {
 if (missing(data) || is.null(data) || !is.data.frame(data)) {
  stop("`data` must be provided as a data.frame.", call. = FALSE)
 }
 zs <- object$m_samples
 N <- if (is.list(zs)) sum(vapply(zs, ncol, 1L)) else ncol(zs)
 on <- object$obs_names
 if (!is.null(on) && length(on) == N && !is.null(rownames(data)) && all(on %in% rownames(data))) {
  return(data[match(on, rownames(data)), , drop = FALSE])
 }
 if (nrow(data) == N) return(data)
 stop(sprintf(paste0(
  "`data` has %d rows but the model was fitted to %d records. Supply the data frame used to ",
  "create the adjustment object (rows are matched by row name) or the %d analysed records ",
  "in the order used by the model."), nrow(data), N, N), call. = FALSE)
}

# Table printed for the pooled results of mi_with(): estimate, standard error,
# 95% confidence interval and degrees of freedom, with the column names of
# summary() and confint().
.mi_pool_table <- function(x) {
 cbind(Estimate = x$coef, `Std. Error` = x$se,
       `2.5 %` = x$ci95[, "lwr"], `97.5 %` = x$ci95[, "upr"], df = x$df)
}

# Stop when mi_with() can use fewer than two posterior draws: the
# between-draw variance, and with it the pooled standard errors, need at least
# two (`m` draws used of `S` stored, `n_small` of which had fewer than `min_n`
# correct-match records).
.mi_check_m <- function(m, S, n_small) {
 if (m >= 2L) return(invisible(TRUE))
 stop(sprintf(paste0(
  "Only %d of %d posterior draws could be used, but at least 2 are needed to estimate the ",
  "between-draw variance of the pooled estimates. Store more draws (increase `iterations`, or ",
  "reduce `burnin.iterations` or `thin`)%s."),
  m, S, if (n_small > 0L) " or lower `min_n`" else ""), call. = FALSE)
}

# Warn when mi_with() could not use every posterior draw.
.mi_skipped_warning <- function(m, S) {
 if (m < S) {
  warning(sprintf(paste0("%d of %d posterior draws were not used: the refit failed or the draw ",
                         "had fewer than `min_n` correct-match records."), S - m, S), call. = FALSE)
 }
 invisible(NULL)
}

# Columns of new data for predict(): matched to the coefficients by name when
# the coefficient names are all non-empty and unique and each of them is the
# name of exactly one column of `newx` (so a reordered model matrix is
# handled), otherwise by position. Fits store non-empty names (see
# .normalize_design_names()); fits saved by earlier versions may hold empty
# or duplicated names (e.g. for X = cbind(1, x)), which are matched by
# position. A warning is given when columns are matched by position although
# some of them bear the name of a coefficient (one that no other coefficient
# has) at another position, e.g. a reordered model matrix for a fit saved
# with the coefficient names c("", "x").
.align_newx <- function(newx, draws, arg) {
 cn <- colnames(draws)
 nx <- colnames(newx)
 by_name <- !is.null(cn) && !anyNA(cn) && all(nzchar(cn)) && !anyDuplicated(cn) &&
  !is.null(nx) && all(cn %in% nx) && !anyDuplicated(nx[nx %in% cn])
 if (by_name) {
  return(newx[, match(cn, nx), drop = FALSE])
 }
 if (ncol(newx) != ncol(draws)) {
  stop(sprintf("`%s` has %d column(s) but the model has %d coefficient(s).",
               arg, ncol(newx), ncol(draws)), call. = FALSE)
 }
 if (!is.null(cn) && !is.null(nx)) {
  named <- !is.na(cn) & nzchar(cn) & !cn %in% cn[duplicated(cn)]
  nx_ok <- !is.na(nx) & nzchar(nx) & !nx %in% nx[duplicated(nx)]
  pos <- match(cn, ifelse(nx_ok, nx, NA_character_))
  moved <- named & !is.na(pos) & pos != seq_along(cn)
  if (any(moved)) {
   many <- sum(moved) > 1L
   warning(sprintf(paste0(
    "The columns of `%s` are matched to the coefficients by position (the coefficient names are not ",
    "all non-empty and unique, or not all found in `%s`), but %s %s %s than the %s of %s name: ",
    "check the order of the columns."),
    arg, arg, if (many) "the columns" else "the column",
    paste(sprintf("\"%s\"", cn[moved]), collapse = ", "),
    if (many) "are at other positions" else "is at another position",
    if (many) "coefficients" else "coefficient", if (many) "the same" else "this"), call. = FALSE)
  }
 }
 newx
}

# Terms of the outcome model of a fit from plglm() / plsurvreg(): stored as
# `terms` or with the model frame (model = TRUE, the default); NULL for fits
# from the engines, which only see a model matrix.
.mixbayes_terms <- function(object) {
 if (!is.null(object$terms)) return(object$terms)
 if (is.data.frame(object$model)) return(attr(object$model, "terms"))
 NULL
}

# Model matrix of a data frame of new records for predict(), built from the
# stored terms of the outcome model (as predict.glmMixture() does), with the
# factor levels of the fitted model frame. `matrix_arg` names the argument
# that takes a model matrix instead (newx for glmMixBayes fits, a numeric
# matrix in newdata for survMixBayes fits).
.mixbayes_newdata_matrix <- function(object, newdata, na.action = stats::na.pass,
                                     matrix_arg = "newx") {
 tt <- .mixbayes_terms(object)
 if (is.null(tt)) {
  engines <- if (inherits(object, "survMixBayes")) c("survregMixBayes()", "plsurvreg") else
   c("glmMixBayes()", "plglm")
  stop("A data frame in `newdata` needs the terms of the outcome model, which this fit does not store ",
       "(fits from ", engines[1], ", or from ", engines[2], "(..., model = FALSE)): pass the model ",
       "matrix as `", matrix_arg, "`.", call. = FALSE)
 }
 xlev <- if (is.data.frame(object$model)) stats::.getXlevels(tt, object$model) else NULL
 tt <- stats::delete.response(tt)
 mf <- stats::model.frame(tt, data = newdata, na.action = na.action, xlev = xlev)
 stats::model.matrix(tt, mf)
}

# Design matrix of the analysed records of a fit, for predict() without new
# data: the matrix stored by plglm()/plsurvreg() with x = TRUE, or else the
# model matrix of the stored model frame (model = TRUE). Both hold every
# record of the outcome model frame, including records dropped at fit time
# because their linkage covariates are missing, so the rows are restricted to
# the analysed records (obs_names). NULL when the fit stores neither.
.mixbayes_fit_design <- function(object) {
 X <- object$x
 if (is.null(X) && is.data.frame(object$model)) {
  tt <- .mixbayes_terms(object)
  if (!is.null(tt)) X <- stats::model.matrix(tt, object$model)
 }
 if (is.null(X)) return(NULL)
 on <- object$obs_names
 if (is.matrix(X) && !is.null(on) && !is.null(rownames(X))) {
  idx <- match(on, rownames(X))
  if (!anyNA(idx)) X <- X[idx, , drop = FALSE]
 }
 X
}

# Outcome formula used by mi_with() when none is given: the terms stored with
# the model frame of a plglm() / plsurvreg() fit (model = TRUE, the default),
# returned as a terms object so that their "predvars" (the bases of
# data-dependent terms such as poly(), computed from the records of the
# outcome model frame) are kept, or else the formula written literally in the
# call. A formula passed through a variable is not re-evaluated, because the
# variable may have been changed or removed since the fit.
.mixbayes_mi_formula <- function(object, env) {
 tt <- if (is.data.frame(object$model)) attr(object$model, "terms") else NULL
 if (!is.null(tt)) return(tt)
 f <- if (is.call(object$call)) object$call$formula else NULL
 if (is.call(f) && identical(f[[1L]], as.name("~"))) return(eval(f, envir = env))
 stop("The outcome formula cannot be recovered from the fit: supply `formula` ",
      "(e.g. mi_with(fit, data, formula = y ~ x)).", call. = FALSE)
}

# Formula and data for refitting with survival::coxph(), which does not
# accept a terms object, while keeping the bases of data-dependent terms
# fixed across the posterior draws (see mi_with.survMixBayes()): every
# variable of the terms `tt` whose prediction call ("predvars", e.g.
# poly(x, 2, coefs = ...) for poly(x, 2)) differs from the variable as
# written is evaluated once on all records of `data`, stored in `data` under
# a placeholder name and replaced by it in the formula, so that each refit
# uses rows of this column instead of recomputing the basis from its own
# records. `labels` maps each placeholder to the variable as written, for the
# coefficient names (see .mi_restore_names()). Penalised terms (class
# "coxph.penalty", e.g. survival::pspline()) cannot be passed as a data
# column, because coxph() needs the call to fit the penalty: the call as
# written stays in the formula, with the arguments that fix the basis added
# (see .mi_penalty_call()). Variables written as calls are evaluated even
# when their prediction call is the same, to find penalised terms whose
# prediction call is not built (makepredictcall() recognises pspline() but
# not survival::pspline()).
.mi_fixed_basis <- function(tt, data) {
 f <- stats::formula(tt)
 vars <- as.list(attr(tt, "variables"))[-1L]
 pv <- attr(tt, "predvars")
 out <- list(formula = f, data = data, labels = character(0))
 if (is.null(pv) || length(pv) != length(vars) + 1L) return(out)
 pv <- as.list(pv)[-1L]
 env <- environment(tt)
 if (is.null(env)) env <- parent.frame()
 from <- list()
 to <- list()
 for (i in seq_along(vars)) {
  same <- identical(pv[[i]], vars[[i]])
  if (same && !is.call(vars[[i]])) next
  value <- if (same) tryCatch(eval(vars[[i]], data, env), error = function(e) NULL) else
   eval(pv[[i]], data, env)
  if (inherits(value, "coxph.penalty")) {
   cl <- .mi_penalty_call(vars[[i]], value, env)
   pen <- if (identical(cl, vars[[i]])) value else eval(cl, data, env)
   if (inherits(pen, "coxph.penalty")) {
    if (!identical(cl, vars[[i]])) {
     from[[length(from) + 1L]] <- vars[[i]]
     to[[length(to) + 1L]] <- cl
    }
    next
   }
   # an unpenalised basis (pspline(..., penalty = FALSE)), with the knots of
   # all records: stored as a data column
   value <- pen
  } else if (same) {
   next
  }
  ph <- sprintf(".mi_fixed_%d_", i)
  out$data[[ph]] <- value
  from[[length(from) + 1L]] <- vars[[i]]
  to[[length(to) + 1L]] <- as.name(ph)
  out$labels[ph] <- paste(deparse(vars[[i]], width.cutoff = 500L), collapse = " ")
 }
 if (length(from) == 0L) return(out)
 # replace the variables by their placeholders (or penalised calls); the other
 # variables are kept as written (and not searched, e.g. I(scale(x)^2) when
 # scale(x) is replaced)
 swap <- function(e) {
  for (k in seq_along(from)) if (identical(e, from[[k]])) return(to[[k]])
  for (v in vars) if (identical(e, v)) return(e)
  if (is.call(e) && length(e) > 1L) {
   for (k in 2:length(e)) e[[k]] <- swap(e[[k]])
  }
  e
 }
 for (k in 2:length(f)) f[[k]] <- swap(f[[k]])
 out$formula <- f
 out
}

# Call of a penalised term for the coxph() refits of mi_with(), given its
# `value` on all records: for survival::pspline(), the call as `written`, so
# that its penalty arguments (df, theta, method, eps, ...) are used in every
# refit, with the arguments that fix the basis (number of knots, boundary
# knots, degree, intercept and combine) set to the attributes of `value`, as
# makepredictcall() does. Its prediction call cannot be used itself, because
# it keeps only these arguments and df (pspline(x, theta = 0.9) would be
# refitted with df = 4, pspline(x, method = "aic") with method "df"). The
# written call is first matched to the arguments of pspline(), so that an
# argument written by position is not given twice. Other penalised terms
# (e.g. frailty()) are returned as written.
.mi_penalty_call <- function(written, value, env) {
 fun <- tryCatch(eval(written[[1L]], env), error = function(e) NULL)
 if (!identical(fun, survival::pspline)) return(written)
 at <- attributes(value)
 basis <- intersect(c("nterm", "intercept", "Boundary.knots", "combine", "degree"), names(at))
 cl <- tryCatch(match.call(fun, written), error = function(e) written)
 for (a in basis) cl[[a]] <- at[[a]]
 cl
}

# Coefficient names of a refit with the formula of .mi_fixed_basis(): the
# placeholders are replaced by the variables as written, which gives the
# names of a fit with the original formula (e.g. "poly(x, 2)1").
.mi_restore_names <- function(nm, labels) {
 if (length(labels) == 0L || is.null(nm)) return(nm)
 for (ph in names(labels)) nm <- gsub(ph, labels[[ph]], nm, fixed = TRUE)
 nm
}

# Error of mi_with() when no posterior draw could be used: say whether the
# draws had too few correct-match records or the refit itself failed (and how).
.mi_no_valid_error <- function(S, n_small, min_n, first_error) {
 if (n_small == S) {
  stop(sprintf(paste0("No valid imputations: every posterior draw has fewer than `min_n` = %d ",
                      "records allocated to the correct-match component."), min_n), call. = FALSE)
 }
 stop(sprintf("No valid imputations: the refit failed for %d of %d posterior draws. First error: %s",
              S - n_small, S, first_error), call. = FALSE)
}

# Labels of interval bounds, formatted as stats::confint() does: the
# percentages are formatted together with 3 significant digits, so that both
# bounds get the decimals the more precise one needs ("0.05 %", "99.95 %" for
# level = 0.999, where formatting 99.95 on its own gives "100"). Trailing
# zeros are then dropped from labels that remain exact without them, so that
# the median of summary(probs = c(0.025, 0.5, 0.975)) is labelled "50 %", not
# "50.0 %". The two bounds of an interval, c(a, 1 - a), never have such zeros
# (a rounded bound such as "51.0" for 51.05 keeps them), so their labels are
# those of stats::confint().
.ci_labels <- function(probs) {
 pct <- 100 * probs
 lab <- format(pct, trim = TRUE, scientific = FALSE, digits = 3)
 short <- sub("\\.?0+$", "", lab)
 exact <- grepl(".", lab, fixed = TRUE) &
  abs(suppressWarnings(as.numeric(short)) - pct) <= 1e-8 * pmax(1, abs(pct))
 exact <- exact %in% TRUE
 lab[exact] <- short[exact]
 paste(lab, "%")
}

# Posterior quantiles of each column of a matrix of draws: one row per column
# of `draws`, one column per probability (labelled as by stats::confint()).
.draw_quantiles <- function(draws, probs) {
 q <- vapply(seq_len(ncol(draws)),
             function(j) stats::quantile(draws[, j], probs = probs, names = FALSE, na.rm = TRUE),
             numeric(length(probs)))
 q <- t(matrix(q, nrow = length(probs)))
 dimnames(q) <- list(colnames(draws), .ci_labels(probs))
 q
}

# Validate the probabilities of posterior quantiles.
.check_probs <- function(probs, arg = "probs") {
 if (!is.numeric(probs) || length(probs) < 1L || anyNA(probs) || any(probs < 0 | probs > 1)) {
  stop(sprintf("`%s` must be a numeric vector of probabilities between 0 and 1.", arg), call. = FALSE)
 }
 invisible(TRUE)
}

# Posterior table of a vector of draws (one row) or of each column of a matrix
# of draws: posterior mean ("Estimate"), posterior standard deviation
# ("Std. Error") and posterior quantiles (by default the central 95% interval,
# "2.5 %" and "97.5 %"). With `mcmc = list(block, ess, rhat)` the Monte Carlo
# columns "MCSE", "ESS" and "Rhat" follow (see .mcmc_columns(); `block` names
# the block of posterior_draws(), `ess` and `rhat` are the stored diagnostics).
.posterior_table <- function(draws, rowname = "", probs = c(0.025, 0.975), mcmc = NULL) {
 scalar <- !is.matrix(draws)
 if (scalar) {
  draws <- matrix(as.numeric(draws), ncol = 1L, dimnames = list(NULL, rowname))
 }
 if (is.null(colnames(draws))) colnames(draws) <- as.character(seq_len(ncol(draws)))
 tab <- cbind(
  Estimate     = colMeans(draws),
  `Std. Error` = apply(draws, 2, stats::sd),
  .draw_quantiles(draws, probs)
 )
 if (!is.null(mcmc)) {
  tab <- cbind(tab, .mcmc_columns(draws, mcmc$block, scalar = scalar, ess = mcmc$ess, rhat = mcmc$rhat))
 }
 rownames(tab) <- colnames(draws)
 tab
}

# Print a posterior table of the summary() methods. The columns keep `digits`
# significant digits, except the Monte Carlo columns: the effective sample
# size ("ESS") is printed as a whole number and the split R-hat ("Rhat") with
# three decimals, as rstan prints it (with significant digits, values such as
# 1.0002 print as 1). Tables without these columns print as before.
.print_posterior_table <- function(tab, digits) {
 if (!is.matrix(tab) || !is.numeric(tab) || !all(c("ESS", "Rhat") %in% colnames(tab))) {
  print(tab, digits = digits)
  return(invisible(tab))
 }
 out <- matrix("", nrow(tab), ncol(tab), dimnames = dimnames(tab))
 for (j in seq_len(ncol(tab))) {
  v <- tab[, j]
  out[, j] <- switch(colnames(tab)[j],
                     ESS  = sprintf("%.0f", v),
                     Rhat = sprintf("%.3f", v),
                     format(v, digits = digits))
 }
 print(out, quote = FALSE, right = TRUE)
 invisible(tab)
}

# Posterior draws of the mismatch-indicator model on the scale and with the
# sign of the m.coefficients of glmMixture() / coxphMixture(), i.e.
# P(mismatch) = plogis(Z %*% m.coefficients): -gamma with linkage covariates
# (Path B), and a single "(Intercept)" column logit(1 - theta) with a constant
# match probability (Path A). Like theta and gamma, these refer to the records
# that are not flagged as safe matches.
.mismatch_coef_draws <- function(theta = NULL, gamma = NULL) {
 if (!is.null(gamma)) {
  if (!is.matrix(gamma)) {
   gamma <- matrix(as.numeric(gamma), ncol = 1L, dimnames = list(NULL, "(Intercept)"))
  }
  return(-gamma)
 }
 if (is.null(theta)) return(NULL)
 # -logit(theta) = logit(1 - theta), without forming 1 - theta
 matrix(-stats::qlogis(as.numeric(theta)), ncol = 1L, dimnames = list(NULL, "(Intercept)"))
}

# Posterior probability that each record is a correct match: the proportion of
# stored draws allocating it to component 1. The allocations are 1/2 integers,
# so this is 2 - colMeans(z), which needs no S x N temporary.
.match_prob <- function(z) {
 2 - colMeans(z)
}

# Bring a fit stored by an earlier version to the current layout of the
# estimates. In postlink <= 0.1.2 (and in the development versions before the
# names were aligned with adjMixture()), estimates$m.coefficients,
# m.dispersion, m.shape and m.scale held the parameters of component 2; they
# are now coefficients2, dispersion2, shape2 and scale2, and m.coefficients
# holds the mismatch-indicator model. Fits without match.prob get it from
# their allocation draws.
.mixbayes_compat <- function(object) {
 est <- object$estimates
 if (is.list(est) && is.null(est$coefficients2) && !is.null(est$m.coefficients)) {
  renamed <- c(m.coefficients = "coefficients2", m.dispersion = "dispersion2",
               m.shape = "shape2", m.scale = "scale2")
  for (old in names(renamed)) {
   if (!is.null(est[[old]])) {
    est[[renamed[[old]]]] <- est[[old]]
    est[[old]] <- NULL
   }
  }
  est$m.coefficients <- .mismatch_coef_draws(est$theta, est$gamma)
  object$estimates <- est
 }
 if (is.null(object$match.prob) && is.matrix(object$m_samples)) {
  object$match.prob <- .match_prob(object$m_samples)
 }
 object
}

# Parameter blocks of the Bayesian fits, as named in `estimates` (and in
# confint(), vcov() and posterior_draws()).
.mixbayes_blocks <- function(object) {
 if (inherits(object, "survMixBayes")) {
  c("coefficients", "coefficients2", "m.coefficients", "theta", "gamma",
    "shape", "shape2", "scale", "scale2")
 } else {
  c("coefficients", "coefficients2", "m.coefficients", "theta", "gamma",
    "dispersion", "dispersion2")
 }
}

# Posterior draws of one parameter block as a named draws x parameters matrix
# (scalar parameters give one column named after the block); informative
# errors for blocks that the fit does not have. `object` must have the current
# layout (see .mixbayes_compat()).
.mixbayes_block_draws <- function(object, block) {
 draws <- object$estimates[[block]]
 if (is.null(draws)) {
  if (block == "gamma") {
   stop("No 'gamma' draws: the match probability was not modeled with covariates ",
        "(see 'm.formula' in adjMixBayes()). Use block = \"theta\" (or \"m.coefficients\").",
        call. = FALSE)
  }
  if (block == "theta") {
   if (isTRUE(object$use_logistic)) {
    stop("No 'theta' draws: the match probability was modeled with covariates. ",
         "Use block = \"gamma\" (or \"m.coefficients\").", call. = FALSE)
   }
   stop("This fit does not store posterior draws of theta (fits created with postlink 0.1.2 or ",
        "earlier do not); refit the model with the current version.", call. = FALSE)
  }
  if (block == "m.coefficients") {
   stop("This fit stores neither theta nor gamma draws, so the mismatch model is not available ",
        "(fits created with postlink 0.1.2 or earlier); refit the model with the current version.",
        call. = FALSE)
  }
  if (block %in% c("dispersion", "dispersion2")) {
   stop(sprintf("No '%s' draws: the %s family has no dispersion parameter.", block, object$family),
        call. = FALSE)
  }
  if (block %in% c("scale", "scale2")) {
   stop(sprintf("No '%s' draws: the scale parameters apply to dist = \"weibull\" only.", block),
        call. = FALSE)
  }
  stop("No posterior draws found for block '", block, "'.", call. = FALSE)
 }
 if (!is.matrix(draws)) draws <- matrix(as.numeric(draws), ncol = 1L, dimnames = list(NULL, block))
 if (is.null(colnames(draws))) colnames(draws) <- as.character(seq_len(ncol(draws)))
 draws
}

# Rows of an interval matrix selected by `parm` (names or indices); all rows
# for parm = NULL.
.select_parm <- function(vals, parm) {
 if (is.null(parm)) return(vals)
 if (is.character(parm)) {
  bad <- setdiff(parm, rownames(vals))
  if (length(bad) > 0L) {
   stop("Unknown parameter(s) in `parm`: ", paste(bad, collapse = ", "), ". Available: ",
        paste(rownames(vals), collapse = ", "), ".", call. = FALSE)
  }
 } else if (!is.numeric(parm) || anyNA(parm) || any(parm < 1 | parm > nrow(vals) | parm != round(parm))) {
  stop("`parm` must be parameter names or indices between 1 and ", nrow(vals), ".", call. = FALSE)
 }
 vals[parm, , drop = FALSE]
}

# Refuse arguments that a method does not use (they would otherwise be
# swallowed by `...` without effect).
.check_unused_args <- function(extra, fun, hint = "") {
 if (length(extra) == 0L) return(invisible(TRUE))
 nm <- names(extra)
 if (is.null(nm)) nm <- rep("", length(extra))
 nm[!nzchar(nm)] <- "(unnamed)"
 stop("Unused argument(s) in ", fun, "(): ", paste(nm, collapse = ", "), ".",
      if (nzchar(hint)) paste0(" ", hint) else "", call. = FALSE)
}

# Average posterior probability of a correct match over all records (safe
# matches count as 1) and over the records that are not flagged as safe
# matches, with the numbers of records. Safe matches have probability exactly
# 1, so the second average is (sum - n_safe) / (n - n_safe).
.match_rate_summary <- function(match_prob, n_safe) {
 n <- length(match_prob)
 if (is.null(n_safe)) n_safe <- 0L
 n_other <- n - n_safe
 list(
  avg = c(all = if (n > 0L) mean(match_prob) else NA_real_,
          not.safe = if (n_other > 0L) (sum(match_prob) - n_safe) / n_other else NA_real_),
  n = c(all = n, not.safe = n_other)
 )
}

# Print the two average correct-match probabilities of a summary.
.print_match_rate <- function(mr, digits) {
 if (is.null(mr)) return(invisible(NULL))
 fmtv <- function(v) if (is.finite(v)) format(signif(v, digits)) else "NA (no such records)"
 lab <- c(sprintf("all records (n = %d; safe matches counted as 1):", as.integer(mr$n[["all"]])),
          sprintf("records not flagged as safe matches (n = %d):", as.integer(mr$n[["not.safe"]])))
 lab <- formatC(lab, width = -max(nchar(lab)))
 cat("Average Correct Match Probability (posterior):\n")
 cat("  ", lab[1], " ", fmtv(mr$avg[["all"]]), "\n", sep = "")
 cat("  ", lab[2], " ", fmtv(mr$avg[["not.safe"]]), "\n", sep = "")
 invisible(NULL)
}

# Headings of the mismatch-model, gamma and theta tables printed by the
# summary() methods of both classes (without the final colon). With safe
# matches (`safe` TRUE) these parameters describe the records not flagged as
# safe matches, which every heading says inside its parentheses.
.match_headings <- function(safe = FALSE) {
 among <- if (isTRUE(safe)) ", among records not flagged as safe matches" else ""
 c(m.coefficients = paste0("Mismatch Model Coefficients (logit of the mismatch probability", among, ")"),
   gamma = paste0("Match Probability Model gamma (logistic regression on the logit scale of the ",
                  "correct-match probability; gamma = -m.coefficients", among, ")"),
   theta = paste0("Match Probability theta (mixing weight of component 1 = correct-match", among, ")"))
}

# The gamma_intercept prior string of fit$priors, or NULL (fits saved before
# the element existed).
.gamma_intercept_prior <- function(priors) {
 if (!is.character(priors) || !"gamma_intercept" %in% names(priors)) return(NULL)
 unname(priors[["gamma_intercept"]])
}

# Note printed under the gamma table of summary() for fits with linkage
# covariates: the sampler centres them at z_center, so the gamma_intercept
# prior (`prior`, the string stored in fit$priors) applies to a record whose
# linkage covariates equal z_center, while the table refers to the original
# covariates. Nothing is printed without centre (intercept-only Z, or fits
# saved before the centring was introduced).
.print_z_center_note <- function(z_center, prior, digits) {
 if (length(z_center) == 0L) return(invisible(NULL))
 nm <- names(z_center)
 if (is.null(nm) || anyNA(nm) || any(!nzchar(nm))) nm <- paste0("Z[, ", seq_along(z_center) + 1L, "]")
 vals <- paste0(nm, " = ", vapply(unname(z_center), format, "", digits = digits), collapse = ", ")
 pr <- if (is.character(prior) && length(prior) == 1L && !is.na(prior)) paste0(" ", prior) else ""
 cat("(the gamma_intercept prior", pr, " applies at the average linkage covariates\n",
     " z_center: ", vals, "; the table refers to the original covariates)\n", sep = "")
 invisible(NULL)
}

# Validate the design matrix X of the engines.
.check_design <- function(X) {
 if (!is.matrix(X) || !is.numeric(X)) stop("`X` must be a numeric matrix.", call. = FALSE)
 if (ncol(X) < 1L) stop("`X` must have at least one column.", call. = FALSE)
 if (nrow(X) < 1L) stop("`X` must have at least one row.", call. = FALSE)
 if (!all(is.finite(X))) {
  stop("`X` must not contain missing or infinite values (check transformations such as log() of zero).",
       call. = FALSE)
 }
 # unnamed columns are named as in the fit (see .normalize_design_names())
 aliased <- .aliased_columns(.normalize_design_names(X))
 if (length(aliased) > 0L) {
  warning("The design matrix has linearly dependent columns (", paste(aliased, collapse = ", "),
          "): their coefficients are identified by the priors only (glm() would report NA).",
          call. = FALSE)
 }
 invisible(TRUE)
}

# Names of the columns of a matrix that are linear combinations of the other
# columns (pivoted QR decomposition, as glm() uses to report NA
# coefficients), "column <j>" for unnamed ones; empty when it has full
# column rank.
.aliased_columns <- function(M) {
 q <- qr(M)
 if (q$rank >= ncol(M)) return(character(0))
 j <- q$pivot[(q$rank + 1L):ncol(M)]
 nm <- colnames(M)
 if (is.null(nm)) nm <- rep("", ncol(M))
 lab <- ifelse(is.na(nm) | !nzchar(nm), paste0("column ", seq_along(nm)), nm)
 lab[j]
}

# Validate the design matrix Z of the match-probability model (Path B), and
# warn when its columns are linearly dependent over the records that inform
# the match-probability model (those not flagged as safe matches, `safe` == 0;
# all records when every record is safe or `safe` is NULL): the aliased
# coefficients are then identified by the gamma_slope prior only.
.check_Z <- function(Z, N, safe = NULL) {
 if (!is.matrix(Z) || !is.numeric(Z)) stop("`Z` must be a numeric matrix.", call. = FALSE)
 if (nrow(Z) != N) stop("`Z` must have the same number of rows as `X`.", call. = FALSE)
 if (ncol(Z) < 1L) stop("`Z` must have at least one column.", call. = FALSE)
 if (!all(is.finite(Z))) stop("`Z` must not contain missing or infinite values.", call. = FALSE)
 if (!all(Z[, 1L] == 1)) {
  stop("The first column of `Z` must be the intercept of the match-probability model (a column of ones).",
       call. = FALSE)
 }
 other <- !is.null(safe) && length(safe) == N && any(safe == 0L) && any(safe != 0L)
 rows <- if (other) safe == 0L else rep(TRUE, N)
 Zs <- Z[rows, , drop = FALSE]
 # unnamed columns are named as in the fit (see .normalize_design_names())
 colnames(Zs) <- colnames(.normalize_design_names(Z, prefix = "Z"))
 aliased <- .aliased_columns(Zs)
 if (length(aliased) > 0L) {
  warning("The design matrix of the match-probability model (Z, from 'm.formula') has linearly dependent ",
          "columns", if (other) " over the records not flagged as safe matches, which inform it" else "",
          " (", paste(aliased, collapse = ", "), "): their coefficients are identified by the ",
          "gamma_slope prior only (glm() would report NA).", call. = FALSE)
 }
 invisible(TRUE)
}

# Centring of the linkage covariates (Path B). The sampler works with the
# non-intercept columns of Z centred at their mean over the records not
# flagged as safe matches (over all records when every record is safe), so
# that the intercept of the centred model, which receives the gamma_intercept
# prior (derived from m.rate, from the theta prior, the default or given by
# the user), refers to a record with average linkage covariates: the point at
# which adjMixture() applies m.rate. The likelihood does not depend on the
# centring and the slopes are unchanged; .uncenter_gamma() maps the draws back
# to the original covariates.
#
# .z_center() returns the centre (a named vector of length ncol(Z) - 1, empty
# for an intercept-only Z); `safe` is the 0/1 vector of safe matches.
.z_center <- function(Z, safe) {
 if (ncol(Z) < 2L) return(stats::setNames(numeric(0), character(0)))
 rows <- if (any(safe == 0L)) safe == 0L else rep(TRUE, nrow(Z))
 colMeans(Z[rows, -1L, drop = FALSE])
}

# Z with its non-intercept columns centred at `center`.
.center_Z <- function(Z, center) {
 if (length(center) == 0L) return(Z)
 Z[, -1L] <- sweep(Z[, -1L, drop = FALSE], 2L, center)
 Z
}

# Draws of gamma on the original covariates from draws of the centred model:
# the intercept is gamma_c[, 1] - gamma_c[, -1] %*% center; the slopes are
# unchanged.
.uncenter_gamma <- function(gamma, center) {
 if (length(center) == 0L || !is.matrix(gamma)) return(gamma)
 gamma[, 1L] <- gamma[, 1L] - drop(gamma[, -1L, drop = FALSE] %*% center)
 gamma
}

# Starting values control$init$gamma are given on the original covariates;
# the sampler needs them for the centred model (intercept + slopes' centre).
# Invalid values are left for .validate_init() to report.
.center_init_gamma <- function(init, center) {
 if (length(center) == 0L || !is.list(init)) return(init)
 g <- init$gamma
 if (!is.numeric(g) || length(g) != length(center) + 1L || !all(is.finite(g))) return(init)
 init$gamma <- c(g[1L] + sum(g[-1L] * center), g[-1L])
 init
}

# Mean linkage covariates of the records of an adjMixBayes object that are not
# flagged as safe matches (all records when every record is safe), as used by
# the fit (.z_center()), for print.adjMixBayes(). NULL when the object has no
# linked data or the design matrix of m.formula cannot be built. The fit uses
# the analysed records only (e.g. after `subset` or missing outcomes), so its
# centre (fit$z_center) can differ slightly. The attribute "all_safe" is TRUE
# when every record used is a safe match (the mean is then over all of them).
.adj_z_center <- function(x) {
 dat <- if (is.environment(x$data_ref)) x$data_ref$data else NULL
 if (is.null(dat) || .is_intercept_only_formula(x$m.formula)) return(NULL)
 tryCatch({
  mf <- stats::model.frame(x$m.formula, data = dat, na.action = stats::na.pass)
  Z <- stats::model.matrix(x$m.formula, mf)
  safe <- if (is.null(x$safe.matches)) rep(0L, nrow(dat)) else as.integer(x$safe.matches)
  keep <- stats::complete.cases(Z) & apply(is.finite(Z), 1L, all)
  if (!any(keep) || length(safe) != nrow(Z)) return(NULL)
  zc <- .z_center(Z[keep, , drop = FALSE], safe[keep])
  attr(zc, "all_safe") <- all(safe[keep] == 1L)
  zc
 }, error = function(e) NULL)
}

# Position of the intercept column of a design matrix (the first column whose
# entries are all 1, whatever its name), or 0 when there is none.
.intercept_column <- function(X) {
 if (nrow(X) < 1L || ncol(X) < 1L) return(0L)
 j <- which(colSums(X != 1) == 0)
 if (length(j) == 0L) 0L else as.integer(j[1L])
}

# Validate the `safe.matches` argument of the engines: NULL, logical or 0/1
# numeric of length N; returns an integer 0/1 vector.
.normalize_safe <- function(safe, N) {
 if (is.null(safe)) return(rep(0L, N))
 if (length(safe) != N) stop("`safe.matches` must have length nrow(X).", call. = FALSE)
 if (anyNA(safe)) stop("`safe.matches` must not contain NA values.", call. = FALSE)
 if (is.logical(safe)) return(as.integer(safe))
 if (is.numeric(safe) && all(safe %in% c(0, 1))) return(as.integer(safe))
 stop("`safe.matches` must be a logical vector or a 0/1 vector.", call. = FALSE)
}

# Whether the prior on the match probability allows correct matches to be the
# majority (it does not say that they are the minority); used for the warning
# on atypical safe matches. With linkage covariates (Path B) the sampler
# centres them (see .z_center()), so the intercept prior refers to a record
# with average linkage covariates: its mean is the prior mean linear predictor
# at that record, also for slope priors with a non-zero mean.
.prior_majority <- function(prior_flat, use_logistic) {
 if (use_logistic) prior_flat$prior_gamma_intercept_mu >= 0 else
  prior_flat$prior_theta_alpha >= prior_flat$prior_theta_beta
}

# Which rule identifies component 1 as the correct-match component (section
# 'Label switching' of ?glmMixBayes), decided before sampling from the safe
# matches and the priors used by the sampler (`engine_priors`, the output of
# build_engine_priors(); with linkage covariates the gamma intercept refers to
# the centred covariates):
#   * has_safe: some records are flagged as safe matches;
#   * match_prior_identifies: the prior on the match probability is not
#     symmetric under theta -> 1 - theta (Path A: a beta(a, b) prior with
#     a != b) or gamma -> -gamma (Path B: a gamma_intercept or gamma_slope
#     prior with a non-zero mean);
#   * component_priors_identical: every component-specific prior (intercepts,
#     slopes, sigma / phi / shape, Weibull scale) is the same for both
#     components;
#   * exchangeable: none of the above identifies the labels, so exchanging
#     them leaves the posterior unchanged.
# Returns these flags with `orientation` ("safe matches", "match-rate prior"
# or "majority component") and `pre_orient`, how the sampler orients its
# starting state: "none" (safe matches), "mass" (towards the labelling with
# the larger posterior probability, under the match-rate prior) or
# "majority" (so that component 1 is the majority component; it refuses when
# the component-specific priors make that labelling much less probable, see
# .ORIENT_TOL). Under the majority convention the draws are relabelled after
# sampling (see align_mixture_labels()).
.label_rule <- function(engine_priors, safe) {
 has_safe <- any(safe != 0L)
 match_prior_identifies <- if (!is.null(engine_priors$gamma_mean)) {
  any(engine_priors$gamma_mean != 0)
 } else {
  engine_priors$theta[1L] != engine_priors$theta[2L]
 }
 same <- function(a, b) length(a) == length(b) && all(a == b)
 keys <- list(c("beta1_mean", "beta2_mean"), c("beta1_sd", "beta2_sd"),
              c("disp1", "disp2"), c("scale1", "scale2"))
 component_priors_identical <- all(vapply(keys, function(k) {
  same(engine_priors[[k[1L]]], engine_priors[[k[2L]]])
 }, logical(1)))
 orientation <- if (has_safe) "safe matches" else if (match_prior_identifies) "match-rate prior" else
  "majority component"
 list(orientation = orientation,
      pre_orient = switch(orientation, "safe matches" = "none", "match-rate prior" = "mass", "majority"),
      exchangeable = !has_safe && !match_prior_identifies && component_priors_identical,
      has_safe = has_safe,
      match_prior_identifies = match_prior_identifies,
      component_priors_identical = component_priors_identical)
}

# Log of the mean of exp(x) over the values of x that are not NA, computed
# stably; NA when there are none.
.log_mean_exp <- function(x) {
 x <- x[!is.na(x)]
 if (length(x) == 0L) return(NA_real_)
 m <- max(x)
 if (!is.finite(m)) return(m)
 m + log(mean(exp(x - m)))
}

# Log posterior and exchange change of the reported draws. The sampler
# reports, for every stored draw, its log posterior `lp` and the change in it
# when the labels of the draw are exchanged (`exchange_lp`, NULL or NA with
# safe matches, see mixbayes_gibbs_cpp()). A draw whose labels were exchanged
# after sampling (`swapped`) has log posterior lp + exchange_lp and exchange
# change -exchange_lp; both are unchanged when the components are
# exchangeable (exchange change 0).
.relabel_lp <- function(lp, exchange_lp, swapped) {
 if (length(exchange_lp) == length(lp) && any(swapped)) {
  lp[swapped] <- lp[swapped] + exchange_lp[swapped]
  exchange_lp[swapped] <- -exchange_lp[swapped]
 }
 list(lp = lp, exchange_lp = exchange_lp)
}

# Estimated share of the posterior probability held by the labelling with the
# two components exchanged, relative to the two labellings, from the exchange
# changes of the reported draws (see .relabel_lp()): the posterior
# probability of the mirror image of a region is the posterior expectation of
# exp(exchange change) over the region, so log mean exp(exchange_lp)
# estimates the log ratio of the posterior probabilities of the other and the
# reported labelling (importance sampling; it assumes that the draws stay in
# one labelling). 0.5 when the components are exchangeable; NA with safe
# matches (no exchange changes), which identify the labels.
.other_labelling_share <- function(exchange_lp) {
 lr <- .log_mean_exp(exchange_lp)
 if (is.na(lr)) NA_real_ else stats::plogis(lr)
}

# Warning on the labelling under the match-rate prior (rule 2 of 'Label
# switching'): the prior on the match probability identifies the labels only
# as far as it makes one labelling more probable than the other. The sampler
# orients its starting state towards the more probable labelling, and the
# draws are reported as sampled; warn when the other labelling holds more than
# `threshold` of the posterior probability (`share`, see
# .other_labelling_share()): the labels are then identified only weakly, or,
# above one half, the chain sat in the less probable labelling. Not checked
# when the draws mix both labellings (`ecr_share` above `ecr_warn`), which
# align_mixture_labels() reports. Returns TRUE when it warned.
.check_label_mass <- function(share, orientation, ecr_share = 0, help = "glmMixBayes",
                              threshold = 0.05, ecr_warn = 0.1) {
 if (!identical(orientation, "match-rate prior") || length(share) != 1L || !is.finite(share) ||
     share <= threshold || isTRUE(ecr_share > ecr_warn)) {
  return(invisible(FALSE))
 }
 fix <- sprintf(paste0("give a more informative prior on the match probability ('m.rate' with a smaller ",
                       "'m.rate.sd', or a theta or gamma_intercept prior), or supply safe matches (see ",
                       "'Label switching' in ?%s)."), help)
 if (share > 0.5) {
  warning(sprintf(paste0(
   "The reported draws describe the less probable of the two labellings of the mixture: the ",
   "labelling with the two components exchanged holds about %.0f%% of the posterior probability ",
   "(estimated from the stored draws; diagnostics$other_labelling), so component 1 %s ",
   "describes the mismatches. The chain did not reach that labelling: refit with another seed, more ",
   "pilot chains or more iterations, or, to identify the correct-match component more firmly, %s"),
   100 * share, if (share > 0.9) "probably" else "may", fix), call. = FALSE)
 } else {
  warning(sprintf(paste0(
   "The prior on the match probability identifies the components only weakly: the labelling with ",
   "the two components exchanged holds about %.0f%% of the posterior probability (estimated from ",
   "the stored draws; diagnostics$other_labelling), while the reported draws describe only the more ",
   "probable labelling, so component 1 may describe the mismatches. To identify the correct-match ",
   "component, %s"), 100 * share, fix), call. = FALSE)
 }
 invisible(TRUE)
}

# Tolerance of the orientation of the starting state under the majority
# convention (section 'Label switching' of ?glmMixBayes). Exchanging all labels
# leaves the likelihood with the indicators integrated out and a symmetric
# prior on the match probability unchanged, so it changes the posterior only
# through the component-specific priors. The sampler estimates, over the last
# sweeps of the selected pilot chain, the log ratio of the posterior
# probabilities of the labelling with component 1 as the majority and of the
# labelling of its starting state; it is exactly 0 for identical component
# priors, and the labels are exchanged when it is at least -.ORIENT_TOL (the
# exchanged labelling is at most exp(2), about 7.4, times less probable).
# Below, the estimate refers to the mirror image of the region visited by the
# pilot chain, while a chain started there can move on to another mode of the
# same labelling, so two check chains decide (see mixbayes_gibbs_cpp()): the
# exchange is refused only when the chain started with the labels exchanged
# returns to the original labelling, or stays below the chain started without
# exchange by more than max(.ORIENT_TOL, 2 se), se being the Monte Carlo
# standard error of the difference of their average log posteriors (about 0
# when the chains mix well). For the binomial defaults (beta1 = normal(0,
# 2.5), beta2 = normal(0, 5)) moving a slope b from component 2 to component 1
# changes the estimate by -0.5 * b^2 * (1 / 2.5^2 - 1 / 5^2) = -0.06 * b^2, so
# it falls below -2 once the squared slopes of the majority exceed those of
# the minority by more than about 33 (many covariates, or large slopes); the
# check chains then usually find a mode of comparable probability with
# component 1 as the majority, and the labels are exchanged. A tight prior on
# component 2 that the majority of the records contradict (e.g. beta2 =
# normal(0, 0.5) or tighter when the majority has a slope of 2) is refused.
.ORIENT_TOL <- 2

# Warn when every record is flagged as a known correct match: the mismatch
# component and the match probability are then informed by the priors only.
.check_all_safe <- function(safe) {
 if (length(safe) > 0L && all(safe == 1L)) {
  warning("All records are flagged as safe matches: the mismatch component and the match ",
          "probability are not informed by the data (their posteriors are the priors).",
          call. = FALSE)
 }
 invisible(NULL)
}

# Effective sample size of one chain: n var(x) / S(0), with the spectral
# density at frequency zero S(0) estimated from an autoregression whose order
# is chosen by AIC (the estimator of coda::effectiveSize()).
.ess_ar <- function(x) {
 x <- as.numeric(x)
 n <- length(x)
 if (n < 10L || !all(is.finite(x))) return(NA_real_)
 v <- stats::var(x)
 if (!is.finite(v) || v <= 0) return(NA_real_)
 fit <- tryCatch(stats::ar(x, aic = TRUE), error = function(e) NULL)
 if (is.null(fit)) return(NA_real_)
 s0 <- fit$var.pred / (1 - sum(fit$ar))^2
 if (!is.finite(s0) || s0 <= 0) return(NA_real_)
 n * v / s0
}

# Effective sample sizes of all stored parameters (named like the columns of
# posterior_draws()).
.mixbayes_ess <- function(fit) {
 d <- .mixbayes_draws(fit)
 stats::setNames(vapply(seq_len(ncol(d)), function(j) .ess_ar(d[, j]), numeric(1)), colnames(d))
}

# Effective sample sizes below `threshold` among the reported parameters (the
# correct-match coefficients and theta or, with linkage covariates, the
# mismatch-model coefficients m.coefficients = -gamma, whose effective sample
# sizes are those of gamma; named as in .mixbayes_rhat()), smallest first.
.low_ess <- function(ess, threshold = 100) {
 if (is.null(ess)) return(numeric(0))
 logistic <- any(grepl("^gamma\\[", names(ess)))
 key <- ess[grepl(if (logistic) "^coefficients\\[|^m\\.coefficients\\[" else "^coefficients\\[|^theta$",
                  names(ess))]
 key <- key[is.finite(key)]
 sort(key[key < threshold])
}

# Note printed by the summary() methods when effective sample sizes are low.
.print_low_ess <- function(low) {
 if (length(low) == 0L) return(invisible(NULL))
 shown <- utils::head(low, 3L)
 cat(sprintf(paste0("Note: low effective sample size for %s%s; these summaries are imprecise. ",
                    "Run more iterations (see fit$diagnostics$ess).\n"),
             paste0(names(shown), " (", round(shown), ")", collapse = ", "),
             if (length(low) > 3L) sprintf(" and %d more", length(low) - 3L) else ""))
 invisible(NULL)
}

# Share of the stored draws in which the chain visited the other labelling of
# the mixture (diagnostics$ecr_share), when the draws were reported as sampled
# and it exceeds `threshold` (the level of the warning issued at fit time);
# NULL otherwise. Used by the summary() methods.
.mixed_labels <- function(diagnostics, threshold = 0.1) {
 share <- diagnostics$ecr_share
 if (is.null(share) || isTRUE(diagnostics$relabelled) || !isTRUE(share > threshold)) return(NULL)
 share
}

# Note printed by the summary() methods when the reported draws mix both
# labellings of the mixture; `help` names the help page of the fitting
# function ("glmMixBayes" or "survregMixBayes").
.print_mixed_labels <- function(share, help = "glmMixBayes") {
 if (is.null(share)) return(invisible(NULL))
 cat(sprintf(paste0("Note: the chain visited both labellings of the mixture (%.1f%% of the stored ",
                    "draws; fit$diagnostics$ecr_share), so these summaries average over the two ",
                    "components and should not be interpreted (see 'Label switching' in ",
                    "?%s).\n"), 100 * share, help))
 invisible(NULL)
}

# Warn when the chain is too short for the reported parameters: the smallest
# effective sample size of the correct-match coefficients and of theta / gamma
# is below `threshold` (rstan applied this threshold to every parameter of the
# Stan-based 0.1.2 fits). Runs with fewer than `min_sweeps` iterations after
# the burn-in (`n_sweeps`, whatever the thinning) are treated as exploratory:
# they are not warned about, and summary() notes low effective sample sizes
# instead. Thinning does not silence the warning, since it lowers the number
# of stored draws (`n_draws`, reported in the message) but not the length of
# the run. The effective sample size cannot much exceed the number of stored
# draws, so a thinned run that stores fewer than twice `threshold` draws is
# also advised to reduce `thin`.
.check_ess <- function(ess, n_draws, n_sweeps = n_draws, threshold = 100, min_sweeps = 1000) {
 low <- .low_ess(ess, threshold)
 if (n_sweeps < min_sweeps || length(low) == 0L) return(invisible(FALSE))
 shown <- utils::head(low, 3L)
 thin_advice <- if (n_draws < n_sweeps && n_draws < 2 * threshold) {
  ", or reduce `thin` to store more draws (the effective sample size cannot much exceed their number)"
 } else ""
 warning(sprintf(paste0(
  "Low effective sample size: %s (out of %d stored draws from %d iterations after the burn-in). ",
  "Posterior means and intervals of these parameters are imprecise; increase `iterations`%s (see ",
  "coda::effectiveSize(coda::as.mcmc(fit)) and fit$diagnostics$ess)."),
  paste0(names(shown), " ", round(shown), collapse = ", "), n_draws, n_sweeps, thin_advice), call. = FALSE)
 invisible(TRUE)
}

# Rank-normalised split R-hat of a single chain (Vehtari et al., 2021), as
# rstan::Rhat() (rstan 2.32) computes it for one chain: the chain is split
# into two halves (the middle draw is dropped when the number of draws is
# odd), the draws of both halves are replaced by the normal scores of their
# joint ranks, and the potential scale reduction factor of the two halves is
# computed; R-hat is the larger of this "bulk" value and the same statistic
# for the draws folded around their median ("tail"). As in rstan, a version
# is undefined when its normal scores are constant (constancy is judged after
# the rank transformation, whatever the scale of the draws). NA for fewer
# than 4 draws, non-finite draws or when neither version is defined; when
# only one is defined, its value, where rstan returns NA: in particular Inf
# when each half of the chain is constant but the halves differ (a chain
# stuck at a different value in each half, whose folded draws are
# constant), so that this case is reported as not converged. A single chain
# cannot show a mode that it never visited: R-hat only detects drift within
# the run.
.rhat_split <- function(x) {
 x <- as.numeric(x)
 n <- length(x)
 if (n < 4L || !all(is.finite(x))) return(NA_real_)
 half <- n %/% 2L
 psrf <- function(v) {
  s <- c(v[seq_len(half)], v[(n - half + 1L):n])
  z <- matrix(stats::qnorm((rank(s, ties.method = "average") - 0.5) / length(s)), ncol = 2L)
  if (max(z) - min(z) < .Machine$double.eps) return(NA_real_)
  w <- mean(c(stats::var(z[, 1L]), stats::var(z[, 2L])))
  b <- half * stats::var(colMeans(z))
  sqrt((b / w + half - 1) / half)
 }
 r <- c(psrf(x), psrf(abs(x - stats::median(x))))
 if (all(is.na(r))) NA_real_ else max(r, na.rm = TRUE)
}

# Single-chain split R-hat (see .rhat_split()) of the parameters checked after
# sampling: the correct-match coefficients, theta (constant match probability)
# or the mismatch-model coefficients (with linkage covariates; their R-hat
# equals that of gamma), the component-1 dispersion (GLMs) or shape (survival
# models), and the log posterior lp. Named as the columns of
# posterior_draws().
.mixbayes_rhat <- function(fit) {
 d <- .mixbayes_draws(fit)
 key <- grepl("^coefficients\\[|^theta$|^dispersion$|^shape$|^lp$", colnames(d))
 if (isTRUE(fit$use_logistic)) key <- key | grepl("^m\\.coefficients\\[", colnames(d))
 d <- d[, key, drop = FALSE]
 stats::setNames(vapply(seq_len(ncol(d)), function(j) .rhat_split(d[, j]), numeric(1)), colnames(d))
}

# R-hat values above `threshold`, largest first (Inf included: see
# .rhat_split()).
.high_rhat <- function(rhat, threshold = 1.05) {
 if (is.null(rhat)) return(numeric(0))
 rhat <- rhat[!is.na(rhat)]
 sort(rhat[rhat > threshold], decreasing = TRUE)
}

# Warn when the single-chain split R-hat of a checked parameter exceeds
# `threshold` (the level rstan used for the Stan-based 0.1.2 fits). As for the
# effective sample size, runs with fewer than `min_sweeps` iterations after the
# burn-in are exploratory: summary() notes high values instead. So are runs
# with fewer than `min_draws` stored draws (heavily thinned runs): the split
# R-hat of nearly independent draws then exceeds 1.05 by chance too often (for
# 100 draws, in about 3% of the parameters, against 0.2% for 200 draws).
.check_rhat <- function(rhat, n_draws, n_sweeps, threshold = 1.05, min_sweeps = 1000, min_draws = 200) {
 high <- .high_rhat(rhat, threshold)
 if (n_sweeps < min_sweeps || n_draws < min_draws || length(high) == 0L) return(invisible(FALSE))
 shown <- utils::head(high, 3L)
 warning(sprintf(paste0(
  "Split R-hat above %s: %s%s (single chain of %d stored draws). The two halves of the chain ",
  "disagree, so it has not converged: increase `iterations` and `burnin.iterations`, or use more ",
  "informative priors (see fit$diagnostics$rhat). A single-chain R-hat detects drift within the ",
  "run but not a mode that the chain never visited: compare fits with several seeds."),
  format(threshold), paste0(names(shown), " ", sprintf("%.3f", shown), collapse = ", "),
  if (length(high) > 3L) sprintf(" and %d more", length(high) - 3L) else "", n_draws), call. = FALSE)
 invisible(TRUE)
}

# Note printed by the summary() methods when the split R-hat of a checked
# parameter is high. With fewer than `min_draws` stored draws (`n_draws`; NULL
# for summaries saved without it) the note says that R-hat is noisy, since
# nearly independent draws then exceed 1.05 by chance more often (the reason
# why .check_rhat() does not warn about such runs).
.print_high_rhat <- function(high, n_draws = NULL, min_draws = 200) {
 if (length(high) == 0L) return(invisible(NULL))
 shown <- utils::head(high, 3L)
 noisy <- length(n_draws) == 1L && isTRUE(n_draws < min_draws)
 verdict <- if (noisy) {
  sprintf(paste0("the chain may not have converged within the run, although the split R-hat of fewer than ",
                 "%d stored draws (here %d) is noisy and exceeds 1.05 by chance more often. Run more ",
                 "iterations (storing at least %d draws)"), min_draws, as.integer(n_draws), min_draws)
 } else {
  "the chain has not converged within the run. Run more iterations"
 }
 cat(sprintf(paste0("Note: split R-hat above 1.05 for %s%s; %s and compare fits with several seeds ",
                    "(see fit$diagnostics$rhat).\n"),
             paste0(names(shown), " (", sprintf("%.3f", shown), ")", collapse = ", "),
             if (length(high) > 3L) sprintf(" and %d more", length(high) - 3L) else "", verdict))
 invisible(NULL)
}

# The diagnostics of a fit with diagnostics$ess and diagnostics$rhat computed
# from the draws when the fit was saved without them (by an earlier version),
# for the summary() methods; left unchanged when the draws cannot be
# collected.
.mixbayes_mcmc_diagnostics <- function(object, dg = object$diagnostics) {
 if (is.null(dg)) dg <- list()
 if (is.null(dg$ess)) dg$ess <- tryCatch(.mixbayes_ess(object), error = function(e) NULL)
 if (is.null(dg$rhat)) dg$rhat <- tryCatch(.mixbayes_rhat(object), error = function(e) NULL)
 dg
}

# Monte Carlo columns of the posterior tables of summary(): the Monte Carlo
# standard error of the posterior mean ("MCSE", the posterior standard
# deviation divided by the square root of the effective sample size), the
# effective sample size ("ESS", rounded) and the single-chain split R-hat
# ("Rhat") of every column of `draws`. Values stored in the fit
# (diagnostics$ess and diagnostics$rhat, named as the columns of
# posterior_draws()) are used when present, the others are computed from the
# draws. `block` names the parameter block in posterior_draws(): a column is
# looked up as block[name], or as block for a parameter stored as a vector
# (`scalar = TRUE`, e.g. theta). A stored value is used only when its name
# occurs once: design columns without names (e.g. cbind(1, x, z)) or with the
# same name give several parameters the same name, whose values are then
# computed from their own column of draws.
.mcmc_columns <- function(draws, block, scalar = FALSE, ess = NULL, rhat = NULL) {
 keys <- if (scalar) rep(block, ncol(draws)) else paste0(block, "[", colnames(draws), "]")
 value <- function(stored, j, fun) {
  if (!is.null(stored) && sum(names(stored) == keys[j]) == 1L && sum(keys == keys[j]) == 1L) {
   as.numeric(stored[[keys[j]]])
  } else {
   fun(draws[, j])
  }
 }
 e <- vapply(seq_len(ncol(draws)), function(j) value(ess, j, .ess_ar), numeric(1))
 r <- vapply(seq_len(ncol(draws)), function(j) value(rhat, j, .rhat_split), numeric(1))
 sdv <- apply(draws, 2, stats::sd)
 cbind(MCSE = sdv / sqrt(e), ESS = round(e), Rhat = r)
}

# Smallest burn-in that collects `needed` draws for the proposal of the joint
# moves, which are collected after the first quarter of the burn-in.
.joint_min_burnin <- function(needed) {
 b <- max(0, floor(4 * needed / 3) - 4)
 while (b - b %/% 4 < needed) b <- b + 1
 b
}

# When the joint moves with the indicators integrated out were off because the
# burn-in was too short to estimate their proposal (diagnostics$accept_joint is
# NA and fewer burn-in draws were collected than needed; see
# diagnostics$joint_burnin): the burn-in used and the smallest one that would
# have sufficed. NULL otherwise (and for fits saved before joint_burnin was
# stored).
.joint_moves_short <- function(diagnostics) {
 aj <- unlist(diagnostics$accept_joint)
 jb <- diagnostics$joint_burnin
 if (length(aj) == 0L || !all(is.na(aj)) || !all(c("collected", "needed") %in% names(jb)) ||
     !isTRUE(jb[["collected"]] < jb[["needed"]])) {
  return(NULL)
 }
 burnin <- diagnostics$settings$burnin.iterations
 c(burnin = if (is.null(burnin)) NA_real_ else as.numeric(burnin),
   needed = .joint_min_burnin(jb[["needed"]]))
}

# Note printed by the summary() methods when the joint moves were off because
# the burn-in was too short (see .joint_moves_short()).
.print_joint_short <- function(js) {
 if (is.null(js)) return(invisible(NULL))
 cat(sprintf(paste0("Note: the joint moves with the match indicators integrated out were off ",
                    "(fit$diagnostics$accept_joint is NA) because the burn-in%s was too short to estimate ",
                    "their proposal; it needs at least %d iterations for this model. The chain may mix ",
                    "more slowly: increase `burnin.iterations`.\n"),
             if (is.finite(js[["burnin"]])) sprintf(" of %d iterations", as.integer(js[["burnin"]])) else "",
             as.integer(js[["needed"]])))
 invisible(NULL)
}

# Names of the coefficients for messages: the column names of the design
# matrix, with an unnamed intercept column (the first column of a design with
# an intercept) called "(Intercept)" and any other unnamed column "X<j>", its
# position (e.g. for X = cbind(1, x)).
.coef_labels <- function(cn, K, intercept = FALSE) {
 if (is.null(cn)) cn <- rep("", K)
 unnamed <- is.na(cn) | !nzchar(cn)
 cn[unnamed] <- paste0("X", which(unnamed))
 if (isTRUE(intercept) && K >= 1L && unnamed[1L]) cn[1L] <- "(Intercept)"
 cn
}

# Which coefficients of component 1 use a default prior, i.e. a prior that the
# user did not supply: the intercept column (the first column of a design with
# an intercept) unless `intercept1` was given, and the other columns unless
# `beta1` was given. `user_priors` is the `priors` list as supplied (NULL for
# none); an entry that is present but NULL counts as not supplied, as in
# prepare_mixbayes_priors(), which then uses the default.
.default_prior_columns <- function(user_priors, K, intercept = FALSE) {
 supplied <- function(key) is.list(user_priors) && !is.null(user_priors[[key]])
 out <- rep(!supplied("beta1"), K)
 if (isTRUE(intercept) && K >= 1L) out[1L] <- !supplied("intercept1")
 out
}

# Estimates of the outcome-model coefficients without the mixture and their
# standard errors, for the prior-scale check of .check_prior_data_conflict():
# a fit to the safe matches when there are at least p + 5 of them (they are
# correct matches), otherwise, or when that fit fails or warns, a fit to all
# records that ignores the linkage errors (the mismatches attenuate the
# estimates but leave their scale). GLMs are fitted by stats::glm.fit() with
# the links of the mixture components (Gaussian models by stats::lm.fit()),
# with the standard errors of summary.glm() (residual variance for the
# Gaussian family, Pearson dispersion for the gamma family); survival models
# by survival::survreg() with dist = "weibull", whose slopes are on the
# accelerated failure time scale of both distributions (its intercept
# includes log(scale), and for the gamma distribution it differs from the log
# mean by a term of order one, both immaterial at the scale of the check).
# A fit to an outcome without variation (e.g. safe matches that all have
# y = 1 in a logistic model, whose intercept then diverges without a warning)
# or, for survival models, without events counts as failed. Returns
# list(est, se, source), with NA for coefficients or standard errors that
# cannot be estimated, or NULL when both fits fail or warn (no convergence,
# separation, ...), so that no conclusion is drawn from them.
.naive_coefficients <- function(X, y, family, event = NULL, safe = NULL) {
 p <- ncol(X)
 # standard errors of a (weighted) least squares fit from its QR decomposition
 qr_se <- function(f, disp) {
  se <- rep(NA_real_, p)
  r <- f$rank
  if (r >= 1L && is.finite(disp)) {
   se[f$qr$pivot[seq_len(r)]] <- sqrt(diag(chol2inv(f$qr$qr[seq_len(r), seq_len(r), drop = FALSE])) * disp)
  }
  se
 }
 fit_rows <- function(rows) {
  Xs <- X[rows, , drop = FALSE]
  if (length(unique(y[rows])) < 2L) return(NULL)
  if (is.null(event)) {
   # least squares for the Gaussian family (glm.fit() computes the deviance,
   # which overflows for outcomes in extreme units)
   if (family == "gaussian") {
    f <- stats::lm.fit(Xs, y[rows])
    disp <- if (f$df.residual > 0L) sum(f$residuals^2) / f$df.residual else NA_real_
    return(list(est = f$coefficients, se = qr_se(f, disp)))
   }
   fam <- switch(family, poisson = stats::poisson(), binomial = stats::binomial(),
                 gamma = stats::Gamma(link = "log"), NULL)
   if (is.null(fam)) return(NULL)
   f <- stats::glm.fit(Xs, y[rows], family = fam)
   if (!isTRUE(f$converged)) return(NULL)
   disp <- if (family != "gamma") 1 else if (f$df.residual > 0L) {
    sum(f$weights * f$residuals^2) / f$df.residual
   } else NA_real_
   return(list(est = f$coefficients, se = qr_se(f, disp)))
  }
  tm <- y[rows]
  ev <- event[rows]
  if (!any(ev == 1L)) return(NULL)
  f <- survival::survreg(survival::Surv(tm, ev) ~ Xs - 1, dist = "weibull")
  v <- diag(stats::vcov(f))
  # the last row and column of vcov() belong to Log(scale)
  list(est = stats::coef(f), se = if (length(v) == p + 1L) sqrt(abs(v[seq_len(p)])) else NULL)
 }
 try_rows <- function(rows) {
  r <- tryCatch(fit_rows(rows), warning = function(w) NULL, error = function(e) NULL)
  if (is.null(r) || length(r$est) != p) return(NULL)
  est <- unname(as.numeric(r$est))
  est[!is.finite(est)] <- NA_real_
  se <- if (length(r$se) == p) unname(as.numeric(r$se)) else rep(NA_real_, p)
  se[!is.finite(se) | is.na(est)] <- NA_real_
  list(est = est, se = se)
 }
 if (!is.null(safe) && sum(safe == 1L) >= p + 5L) {
  r <- try_rows(which(safe == 1L))
  if (!is.null(r)) return(c(r, list(source = "the safe matches")))
 }
 r <- try_rows(seq_len(nrow(X)))
 if (is.null(r)) NULL else c(r, list(source = "all records"))
}

# Prior-scale and prior-data conflict check of the coefficients of the
# correct-match component, issued after sampling as one warning per fit:
#   * before sampling, the fitting functions estimate the coefficients
#     without the mixture (`naive`, see .naive_coefficients(): a fit to the
#     safe matches, or to all records ignoring the linkage errors); a
#     coefficient with a default prior (`default`, see
#     .default_prior_columns()) is flagged when this estimate lies more than
#     `k` prior standard deviations from the prior mean: the default priors are
#     not rescaled to the data and suit outcomes and covariates of order one.
#     When the data carry little information relative to such a prior (e.g. a
#     Gaussian outcome in very large units, whose residual SD absorbs the
#     outcome, or a covariate in very small units), the posterior collapses
#     onto the prior, which the posterior check below cannot see; the message
#     then says that the prior dominates the estimate (posterior standard
#     deviation above `dominance` times the prior standard deviation): the
#     posterior stays close to the prior, or, when its mean lies more than `k`
#     prior standard deviations from the prior mean as well, its standard
#     deviation stays close to the prior standard deviation;
#   * after sampling, a coefficient with a supplied prior is flagged when the
#     prior dominates its estimate although the naive estimate lies more than
#     `k` prior standard deviations from the prior mean and, allowing for its
#     standard error `se`, more than `k` times sqrt(prior SD^2 + se^2) from it
#     (the data then contradict the prior; a supplied informative prior may
#     dominate data that carry little information without a warning);
#   * after sampling, any coefficient (default or supplied prior) whose
#     posterior mean lies more than `k` prior standard deviations from its
#     prior mean is flagged: the prior is far from the data.
# The message describes the coefficient farthest from its prior and counts the
# others. `b1` holds the draws as sampled (for Weibull fits the sampled
# intercept, on which its prior is placed), or is NULL when sampling stopped
# (only the first check is then made, and `coef_names` names the
# coefficients); `engine_priors` is the output of build_engine_priors(),
# `default = NULL` treats every prior as a default, and `help` names the help
# page of the fitting function. `reported` (Weibull fits with an intercept
# column) holds the posterior means as summary() reports them, i.e. with
# log(scale) added to the intercept; the message then quotes both when it
# describes the intercept. Returns TRUE when it warned.
.check_prior_data_conflict <- function(b1, engine_priors, k = 3, help = "glmMixBayes",
                                       default = NULL, naive = NULL, dominance = 0.9,
                                       coef_names = colnames(b1), reported = NULL) {
 mu <- engine_priors$beta1_mean
 sdp <- engine_priors$beta1_sd
 K <- if (is.null(b1)) length(mu) else ncol(b1)
 if (is.null(default)) default <- rep(TRUE, K)
 pm <- if (is.null(b1)) rep(NA_real_, K) else colMeans(b1)
 psd <- if (!is.null(b1) && nrow(b1) > 1L) apply(b1, 2, stats::sd) else rep(NA_real_, K)
 z_post <- (pm - mu) / sdp
 z_naive <- if (is.null(naive)) rep(NA_real_, K) else (naive$est - mu) / sdp
 # distance of the naive estimate from the prior mean allowing for its
 # standard error (NA when it is not available)
 se_naive <- if (is.null(naive$se)) rep(NA_real_, K) else naive$se
 z_pred <- if (is.null(naive)) rep(NA_real_, K) else (naive$est - mu) / sqrt(sdp^2 + se_naive^2)
 far <- is.finite(z_naive) & abs(z_naive) > k
 dominated <- far & is.finite(psd) & psd > dominance * sdp
 scale_flag <- (default & far) | (!default & dominated & is.finite(z_pred) & abs(z_pred) > k)
 conflict <- is.finite(z_post) & abs(z_post) > k
 flagged <- scale_flag | conflict
 if (!any(flagged)) return(invisible(FALSE))
 score <- pmax(ifelse(scale_flag, abs(z_naive), 0), ifelse(conflict, abs(z_post), 0))
 j <- which(flagged)[which.max(score[flagged])]
 cn <- .coef_labels(coef_names, K, isTRUE(engine_priors$intercept))
 num <- function(v) format(signif(v, 4))
 # distances in prior standard deviations (huge ones in exponent notation)
 nsd <- function(v) if (abs(v) < 1e5) sprintf("%.1f", abs(v)) else format(signif(abs(v), 3))

 what <- sprintf("The %s prior of '%s' in the correct-match component, normal(%s, %s), ",
                 if (default[j]) "default" else "supplied", cn[j], num(mu[j]), num(sdp[j]))
 why <- if (scale_flag[j]) {
  # the safe matches are correct matches: only a fit to all records ignores
  # the linkage errors
  fit_desc <- if (identical(naive$source, "the safe matches")) {
   "a fit of the outcome model to the safe matches"
  } else {
   sprintf("a fit of the outcome model to %s that ignores the linkage errors", naive$source)
  }
  # a default prior is flagged for its scale, a supplied one (possibly a
  # tight prior on the right scale) for contradicting the data
  paste0(sprintf("%s: %s gives %s, %s prior standard deviations from the prior mean",
                 if (default[j]) "does not suit the scale of the data" else "conflicts with the data",
                 fit_desc, num(naive$est[j]), nsd(z_naive[j])),
         if (dominated[j] && !conflict[j]) {
          sprintf(paste0(", while the posterior (mean %s, standard deviation %s) stays close to the ",
                         "prior: the prior dominates the estimate"), num(pm[j]), num(psd[j]))
         } else if (dominated[j]) {
          sprintf(paste0(", and the posterior mean (%s) lies %s prior standard deviations from it, while ",
                         "the posterior standard deviation (%s) stays close to the prior standard deviation: ",
                         "the prior still dominates the estimate"),
                  num(pm[j]), nsd(z_post[j]), num(psd[j]))
         } else if (conflict[j]) {
          sprintf(", and the posterior mean (%s) lies %s prior standard deviations from it",
                  num(pm[j]), nsd(z_post[j]))
         } else if (is.finite(pm[j])) {
          sprintf(", while the posterior mean (%s) lies closer to the prior mean: the prior may distort the estimate",
                  num(pm[j]))
         } else "")
 } else {
  sprintf(paste0("is far from the data: the posterior mean (%s) lies %s prior standard deviations ",
                 "from the prior mean, and a prior on another scale than the data can distort the fit"),
          num(pm[j]), nsd(z_post[j]))
 }
 # Weibull fits: the prior is placed on the sampled intercept, whereas
 # summary() reports the intercept with log(scale) added
 if (!is.null(reported) && length(reported) == K && is.finite(pm[j]) && is.finite(reported[j]) &&
     reported[j] != pm[j]) {
  why <- paste0(why, sprintf(paste0(" (the prior and the posterior quoted here refer to the sampled intercept; ",
                                    "summary() reports the intercept with log(scale) added, whose posterior ",
                                    "mean is %s)"), num(reported[j])))
 }
 others <- setdiff(which(flagged), j)
 count <- if (length(others) > 0L) {
  sprintf(" (%d coefficients are affected; also %s%s)", sum(flagged),
          paste0("'", utils::head(cn[others], 3L), "'", collapse = ", "),
          if (length(others) > 3L) ", ..." else "")
 } else ""
 advice <- if (any(default[flagged])) {
  paste0("The default priors are not rescaled to the data: they suit outcomes and covariates of order ",
         "one. Supply priors that suit the scale of the outcome and covariates, or standardise them")
 } else {
  paste0("Check that the supplied priors suit the scale of the outcome and covariates and are not more ",
         "informative than intended, or standardise the outcome and covariates")
 }
 warning(paste0(what, why, count, ". ", advice, sprintf(" (see 'Prior distributions' in ?%s).", help)),
         call. = FALSE)
 invisible(TRUE)
}

# Validate user starting values (control$init) for a given engine family,
# number of coefficients K and number of match-model coefficients M (0 in
# Path A). Only the elements that apply to the model are accepted: beta1,
# beta2 (length K); theta (Path A) or gamma (length M, Path B); disp1, disp2
# (gaussian residual SD, gamma shape, survival shape); scale1, scale2 (Weibull).
.validate_init <- function(init, family, K, M) {
 if (is.null(init)) return(list())
 if (!is.list(init)) stop("`control$init` must be a named list.", call. = FALSE)
 if (length(init) == 0L) return(list())
 if (is.null(names(init)) || anyNA(names(init)) || any(!nzchar(names(init)))) {
  stop("`control$init` must be a named list.", call. = FALSE)
 }
 all_names <- c("beta1", "beta2", "theta", "gamma", "disp1", "disp2", "scale1", "scale2")
 unknown <- setdiff(names(init), all_names)
 if (length(unknown) > 0L) {
  stop("Unknown element(s) in `control$init`: ", paste(unknown, collapse = ", "),
       ". Allowed: ", paste(all_names, collapse = ", "), ".", call. = FALSE)
 }
 has_disp <- family %in% c("gaussian", "gamma", "surv_gamma", "surv_weibull")
 expected <- list(beta1 = K, beta2 = K)
 if (M > 0L) expected$gamma <- M else expected$theta <- 1L
 if (has_disp) expected$disp1 <- expected$disp2 <- 1L
 if (identical(family, "surv_weibull")) expected$scale1 <- expected$scale2 <- 1L
 why <- c(theta = "the match probability is modelled with covariates (Z / m.formula): use `gamma`",
          gamma = "no linkage covariates (Z / m.formula) were supplied: use `theta`",
          disp1 = "this family has no dispersion or shape parameter",
          disp2 = "this family has no dispersion or shape parameter",
          scale1 = "the scale parameters apply to dist = \"weibull\" only",
          scale2 = "the scale parameters apply to dist = \"weibull\" only")
 bad <- setdiff(names(init), names(expected))
 if (length(bad) > 0L) {
  stop(paste0("`control$init$", bad, "` does not apply to this model: ", why[bad], ".",
              collapse = " "), call. = FALSE)
 }
 for (nm in names(init)) {
  v <- init[[nm]]
  if (!is.numeric(v) || length(v) != expected[[nm]] || anyNA(v) || any(!is.finite(v))) {
   stop(sprintf("`control$init$%s` must be a finite numeric vector of length %d.", nm, expected[[nm]]),
        call. = FALSE)
  }
  if (nm == "theta" && (v <= 0 || v >= 1)) {
   stop("`control$init$theta` must lie strictly between 0 and 1.", call. = FALSE)
  }
  if (nm %in% c("disp1", "disp2", "scale1", "scale2") && v <= 0) {
   stop(sprintf("`control$init$%s` must be positive.", nm), call. = FALSE)
  }
  init[[nm]] <- as.double(v)
 }
 init
}

# Run the Gibbs sampler.
#
# @param family Engine family: "gaussian", "poisson", "binomial", "gamma",
#   "surv_gamma" or "surv_weibull".
# @param X Design matrix. @param y Outcome (survival time for survival models).
# @param event Integer event indicator (survival) or NULL.
# @param Z Design matrix of the match-probability model (Path B) or NULL.
# @param safe Integer vector of known-correct-match flags (0/1).
# @param priors Output of build_engine_priors().
# @param control Evaluated MCMC settings (see .mixbayes_settings()).
# @param pre_orient How the sampler orients the starting state of the main
#   chain (see .label_rule(); ignored with starting values in control$init
#   and with safe matches): "none"; "majority", exchanging its labels when
#   component 2 holds the majority of the records (TRUE is accepted for it,
#   FALSE for "none"); "mass", exchanging them when the labelling with the
#   components exchanged has the larger estimated posterior probability.
# @param orient_tol Under "majority": the labels of the starting state are
#   exchanged when the estimated log ratio of the posterior probabilities of
#   the exchanged and the original labelling is at least -orient_tol; below,
#   check chains decide, with a margin of at least orient_tol (`Inf`
#   exchanges the labels whatever the component-specific priors, without
#   check chains; see .ORIENT_TOL).
# @param on_lp_stop Optional function without arguments, called before the
#   fit stops because no stored draw has a finite log posterior (see
#   .check_lp()); the fitting functions report the prior-scale check there.
run_mixbayes_engine <- function(family, X, y, event, Z, safe, priors, control, pre_orient = FALSE,
                                orient_tol = .ORIENT_TOL, on_lp_stop = NULL) {
 storage.mode(X) <- "double"
 if (is.null(Z)) Z <- matrix(0, nrow(X), 0L) else storage.mode(Z) <- "double"
 if (is.null(event)) event <- integer(0)
 init <- .validate_init(control$init, family, ncol(X), ncol(Z))
 control <- .normalize_mcmc_settings(control)

 # the stored allocations take S x N integers; say so before sampling when that is large
 S <- (control$iterations - control$burnin.iterations) %/% control$thin
 gb <- as.double(S) * nrow(X) * 4 / 1e9
 if (gb > 1) {
  warning(sprintf(paste0(
   "Storing the match indicators of %d draws x %d records needs about %.1f GB of memory ",
   "(and about as much again while the labels are aligned); use `thin` to store fewer draws."),
   S, nrow(X), gb), call. = FALSE)
 }

 orient_code <- if (is.logical(pre_orient)) {
  as.integer(isTRUE(pre_orient))
 } else {
  match(match.arg(pre_orient, c("none", "majority", "mass")), c("none", "majority", "mass")) - 1L
 }
 post <- .with_rng_seed(
  control$seed,
  mixbayes_gibbs_cpp(family, X, as.double(y), as.integer(event), Z, as.integer(safe),
                     priors, init, control$iterations, control$burnin.iterations, control$thin,
                     control$pilots, max(1L, control$pilot.iterations),
                     isTRUE(control$collapse), isTRUE(control$verbose), orient_code,
                     as.double(orient_tol))
 )
 # scalar parameters come back as S x 1 matrices; return plain vectors
 for (nm in intersect(names(post), c("theta", "disp1", "disp2", "scale1", "scale2", "lp", "exchange_lp"))) {
  post[[nm]] <- as.numeric(post[[nm]])
 }
 .check_lp(post$lp, post$pilot_lp, on_stop = on_lp_stop)
 post
}

# Checks of the log posterior reported by the sampler, whether or not pilot
# chains were run (they are skipped with control$init and with pilots = 0):
# stop when no stored draw has a finite log posterior (the draws are then not
# valid), warn with the share of the stored draws whose log posterior is not
# finite, and warn when no pilot chain reached a finite log posterior (the
# main chain then started from the initial values). Each points to numerical
# problems, typically an outcome or covariates on extreme scales. `on_stop`
# (a function without arguments, or NULL) is called before stopping, e.g. to
# report the prior-scale check, which such scales usually trigger as well.
.check_lp <- function(lp, pilot_lp = numeric(0), on_stop = NULL) {
 advice <- paste0("This points to numerical problems, typically an outcome or covariates on extreme ",
                  "scales: consider centring and scaling them (or changing their units).")
 if (length(lp) > 0L && !any(is.finite(lp))) {
  if (is.function(on_stop)) on_stop()
  stop("The log posterior was not finite in any stored draw, so the posterior draws are not valid. ",
       advice, call. = FALSE)
 }
 if (length(pilot_lp) > 0L && !any(is.finite(pilot_lp))) {
  warning("The log posterior was not finite in any pilot chain, so the main chain started from the ",
          "initial values. ", advice, call. = FALSE)
 }
 bad <- !is.finite(lp)
 if (any(bad)) {
  warning(sprintf(paste0("The log posterior was not finite in %d of %d stored draws (%.1f%%): the ",
                         "results may not be valid. %s"),
                  sum(bad), length(lp), 100 * mean(bad), advice), call. = FALSE)
 }
 invisible(TRUE)
}

# Acceptance rates of the coefficient blocks under the reported labels: the
# sampler reports them under its raw labels, which a global relabelling exchanges.
.relabel_accept <- function(accept, flipped) {
 if (isTRUE(flipped)) accept[c("beta1", "beta2")] <- accept[c("beta2", "beta1")]
 accept
}

# Warn when a Metropolis-Hastings block was rarely accepted: the Laplace
# proposal then approximates the conditional posterior poorly (typically a
# weakly identified component under diffuse priors). `n_sweeps` is the number
# of sweeps the rates refer to (those after the burn-in, all sweeps without
# burn-in: iterations - burnin.iterations in both cases); with fewer than
# `min_sweeps` the rates are too imprecise for a warning (with one sweep they
# are 0 or 1).
.check_acceptance <- function(accept, n_sweeps = Inf, threshold = 0.2, min_sweeps = 50) {
 acc <- unlist(accept)
 if (n_sweeps < min_sweeps) return(invisible(acc))
 low <- names(acc)[!is.na(acc) & acc < threshold]
 if (length(low) > 0L) {
  role <- c(beta1 = "the coefficients of component 1 (correct matches)",
            beta2 = "the coefficients of component 2 (mismatches)",
            gamma = "the match-probability coefficients (gamma)")
  lab <- ifelse(low %in% names(role), role[low], low)
  warning(sprintf(paste0(
   "Low Metropolis-Hastings acceptance rate for %s (%s). The posterior of this ",
   "block may be poorly explored; consider more informative priors, centring and scaling ",
   "the covariates, or more iterations."),
   paste(lab, collapse = "; "),
   paste(sprintf("%.2f", acc[low]), collapse = ", ")), call. = FALSE)
 }
 invisible(acc)
}

################################################################################
# Label alignment across MCMC draws
################################################################################

# ECR-ITERATIVE-1 relabelling (Papastamoulis & Iliopoulos, 2010; Rodriguez &
# Walker, 2014) specialised to two components. Equivalent to
# label.switching::ecr.iterative.1() with K = 2: each draw is either kept or
# swapped so as to minimise the number of disagreements with a pivot
# allocation, and the pivot is iteratively re-estimated as the modal allocation.
#
# @param z S x N integer matrix of allocations (values 1 or 2).
# @param cols Optional column indices: only these records are used (e.g. the
#   records not flagged as safe matches); all records when NULL.
# @return Logical vector of length S; TRUE where the two labels of that draw
#   must be exchanged.
#
# The iterations run in C++ (ecr_two_cpp() in src/mixbayes_labels.cpp), which
# needs no S x N temporaries (and no copy of the selected columns).
ecr_iterative_two <- function(z, maxiter = 100L, threshold = 1e-6, cols = NULL) {
 if (!is.matrix(z)) stop("`z` must be an S x N matrix.", call. = FALSE)
 ecr_two_cpp(.as_integer_matrix(z), as.integer(maxiter), as.double(threshold),
             if (is.null(cols)) integer(0) else as.integer(cols))
}

# z as an integer matrix (no copy when it already is one).
.as_integer_matrix <- function(z) {
 if (is.integer(z)) return(z)
 matrix(as.integer(z), nrow(z), ncol(z), dimnames = dimnames(z))
}

# Warning issued under the majority convention when component 1 holds a
# minority of the records on average: in the fits, when the
# component-specific priors make the labelling with component 1 as the
# majority much less probable, so that neither the starting state
# (`start$refused`) nor the aligned draws (`start$flip_refused`, with the
# estimated log ratio `start$flip_lp_change`; see align_mixture_labels())
# were exchanged (the names avoid partial matching by `$`). Worded after the starting state of the main chain (`start`,
# see .orientation_start(); NULL when unknown): without a refusal, its share
# of component 1 after any exchange of its labels decides whether the chain
# started with component 1 as the majority. `help` names the help page of the
# fitting function.
.minority_message <- function(share1, start = NULL, help = "glmMixBayes") {
 head <- sprintf(paste0(
  "Only %.0f%% of the records are allocated to component 1 on average, although under the ",
  "majority convention component 1 is meant to be the majority component. "),
  100 * share1)
 # share of component 1 in the starting state of the main chain
 share_start <- if (is.numeric(start$share1) && length(start$share1) == 1L && is.finite(start$share1)) {
  if (isTRUE(start$pre_oriented)) 1 - start$share1 else start$share1
 } else NA_real_
 fix <- "give the expected mismatch rate with 'm.rate', or supply safe matches."
 lower_by <- function(change) {
  if (length(change) == 1L && is.finite(change)) {
   sprintf("lower by about %.1f (more than %g)", -change, .ORIENT_TOL)
  } else {
   "infinitely lower (that labelling has posterior probability zero)"
  }
 }
 why <- if (isTRUE(start$refused) || isTRUE(start$flip_refused)) {
  intro <- if (isTRUE(start$refused)) {
   # result of the check chains (see mixbayes_gibbs_cpp()), when they were run
   chk <- start$check
   checked <- if (!all(c("exchanged", "kept", "margin", "share1") %in% names(chk)) ||
                  !is.finite(chk[["share1"]])) {
    ""
   } else if (chk[["share1"]] < 0.5) {
    paste0(", and a check chain started with the labels exchanged returned to the labelling in ",
           "which component 1 is the minority")
   } else if (is.finite(chk[["exchanged"]])) {
    sprintf(paste0(", and a check chain started with the labels exchanged stayed %.1f below one ",
                   "started without them exchanged (more than the margin of %.1f)"),
            chk[["kept"]] - chk[["exchanged"]], chk[["margin"]])
   } else {
    ", and a check chain started with the labels exchanged did not reach a finite log posterior"
   }
   sprintf(paste0(
    "The component-specific priors favour the labelling in which component 1 is the minority: ",
    "the estimated log posterior probability of the labelling with component 1 as the majority ",
    "is %s%s, so the sampler did not exchange the labels of its starting state%s ",
    "(diagnostics$orient_refused; see 'Label switching' in ?%s). "), lower_by(start$lp_change), checked,
    if (isTRUE(start$flip_refused)) ", and the draws were not exchanged after sampling either" else "", help)
  } else {
   sprintf(paste0(
    "The component-specific priors favour the labelling in which component 1 is the minority: ",
    "the log posterior probability of the labelling with component 1 as the majority, estimated ",
    "from the stored draws, is %s, so the draws were not exchanged after sampling (see 'Label ",
    "switching' in ?%s). "), lower_by(start$flip_lp_change), help)
  }
  if (isTRUE(start$default_priors)) {
   paste0(intro,
    "These are the default priors, which differ between the components (for \"binomial\", ",
    "beta1 = normal(0, 2.5) and beta2 = normal(0, 5), so that large slopes are shrunk less in ",
    "component 2) without being meant to identify them: component 1 may describe the mismatches. ",
    "Give the expected mismatch rate with 'm.rate', or supply safe matches, to identify the ",
    "correct-match component.")
  } else {
   paste0(intro,
    "If these priors describe the two components as intended, component 1 is the correct-match ",
    "component; giving the expected mismatch rate with 'm.rate', or safe matches, identifies the ",
    "components without this warning. Otherwise check the component-specific priors.")
  }
 } else if (isTRUE(start$init) && isTRUE(start$share1 < 0.5)) {
  paste0("The starting values in control$init put component 1 on the minority of the records, ",
         "and starting values are never reoriented. Change the starting values, ", fix)
 } else if (isTRUE(share_start >= 0.5)) {
  paste0(sprintf(paste0("The sampler started with component 1 as the majority (%.0f%% of the records in ",
                        "its starting state%s), but the chain moved to the labelling in which it is ",
                        "the minority. Check the component-specific priors, "),
                 100 * share_start,
                 if (isTRUE(start$pre_oriented)) ", after the exchange of its labels" else ""), fix)
 } else {
  paste0("Check the component-specific priors, ", fix)
 }
 paste0(head, why)
}

# Exchange the two components of a pair of draws for the selected iterations.
# Draws can be numeric vectors (length S) or matrices (S x p).
.swap_pair <- function(a1, a2, swap) {
 if (is.matrix(a1)) {
  tmp <- a1[swap, , drop = FALSE]
  a1[swap, ] <- a2[swap, , drop = FALSE]
  a2[swap, ] <- tmp
 } else {
  tmp <- a1[swap]
  a1[swap] <- a2[swap]
  a2[swap] <- tmp
 }
 list(a1, a2)
}

# Align component labels across posterior draws (see 'Label switching' in
# ?glmMixBayes and .label_rule()).
#
# 1. ECR-ITERATIVE-1 (Papastamoulis & Iliopoulos, 2010), implemented natively
#    for two components, is run on the records not flagged as safe matches in
#    every case. The share of the draws that it exchanges, or would exchange
#    (`ecr_share`), is a diagnostic: a large share means that the chain
#    visited both labellings of the mixture.
# 2. Relabelling after sampling (`relabel = TRUE`) is applied under the
#    majority convention (no safe matches and a prior on the match
#    probability that is symmetric under exchanging the labels). The draws
#    are aligned by the ECR step and oriented globally so that component 1 is
#    the majority component, unless the component-specific priors make that
#    labelling much less probable: the log ratio of the posterior
#    probabilities of the exchanged and the aligned labelling, estimated from
#    the exchange changes of the aligned draws (`exchange_lp`, see
#    .relabel_lp() and .other_labelling_share()), is below -`orient_tol`, the
#    tolerance of the orientation of the starting state (see .ORIENT_TOL).
#    Component 1 then stays the minority component (`flip_refused`; the
#    priors identify the labels, and the warning on a minority component 1
#    says so). Relabelling renames the components of each draw without
#    changing the fit it describes, so the relabelled draws are posterior
#    draws of the parameters of the majority component (of the components as
#    aligned when the global exchange is refused), also when the
#    component-specific priors differ and the chain visited both labellings;
#    only their log posterior changes (see .relabel_lp()), unless the priors
#    are identical.
# 3. Otherwise (safe matches or an informative prior on the match
#    probability) the draws are returned exactly as sampled, and a warning is
#    issued when the ECR step would have exchanged more than a share
#    `ecr_warn` of them.
#
# @param z S x N allocation matrix.
# @param pairs Named list; each element is list(a1, a2) holding the draws of a
#   component-specific parameter (vector of length S or S x p matrix).
# @param theta Optional vector of mixing-weight draws (Path A).
# @param gamma Optional vector/matrix of match-model coefficient draws (Path B).
# @param safe Optional 0/1 vector of known correct matches.
# @param orientation What identifies component 1: "safe matches", "match-rate
#   prior" or "majority component" (see .label_rule()); by default "safe
#   matches" with safe matches and "majority component" otherwise.
# @param relabel Whether the draws are relabelled after sampling; by default
#   TRUE for "majority component". Only under the majority convention, never
#   with safe matches.
# @param expect_majority With safe matches: whether the prior on the match
#   probability allows correct matches to be the majority; a warning is then
#   issued when component 1 (anchored by the safe matches) holds a minority of
#   the other records.
# @param verbose Whether to report a global label swap with a message.
# @param ecr_warn Share of draws above which the ECR diagnostic warns when no
#   relabelling is applied.
# @param start Optional description of the starting state of the main chain
#   (see .orientation_start()), used to word the warning issued under the
#   majority convention when component 1 holds a minority of the records:
#   list(share1, pre_oriented, init, refused, lp_change, check,
#   default_priors).
# @param help Help page named by the warnings: "glmMixBayes" or
#   "survregMixBayes".
# @param exchange_lp Optional vector of length S: the change in the log
#   posterior of each draw as sampled when its labels are exchanged (see
#   mixbayes_gibbs_cpp()); with relabelling, it decides whether the global
#   orientation may exchange the labels (treated as 0, i.e. exchangeable
#   components, when not given).
# @param orient_tol Tolerance of that decision (see .ORIENT_TOL).
# @return list(z, pairs, theta, gamma, swapped, orientation, flipped,
#   flip_refused, flip_lp_change, relabelled, ecr_share, share1) with aligned
#   draws: `swapped` flags the draws whose labels were exchanged after
#   sampling, `flipped` whether the global orientation exchanged the labels,
#   `flip_refused` whether it was refused because of the component-specific
#   priors, with the estimated log ratio `flip_lp_change` (NA when no global
#   exchange was considered), `relabelled` whether relabelling after
#   sampling was applied, `ecr_share` the share of draws that the ECR step
#   exchanges (or would exchange), and `share1` the average share of the
#   records not flagged as safe matches that are allocated to component 1.
align_mixture_labels <- function(z, pairs, theta = NULL, gamma = NULL, safe = NULL,
                                 orientation = NULL, relabel = NULL,
                                 expect_majority = FALSE, verbose = FALSE, ecr_warn = 0.1,
                                 start = NULL, help = "glmMixBayes",
                                 exchange_lp = NULL, orient_tol = .ORIENT_TOL) {
 if (!is.matrix(z)) stop("`z` must be an S x N matrix.", call. = FALSE)
 z <- .as_integer_matrix(z)
 S <- nrow(z)
 if (is.null(safe)) safe <- rep(0L, ncol(z))
 has_safe <- any(safe != 0)
 if (is.null(orientation)) orientation <- if (has_safe) "safe matches" else "majority component"
 if (is.null(relabel)) relabel <- orientation == "majority component"
 if (isTRUE(relabel) && (has_safe || orientation != "majority component")) {
  stop("Draws can be relabelled after sampling only under the majority convention.",
       call. = FALSE)
 }
 other <- which(safe == 0)               # records not flagged as safe matches
 flipped <- FALSE
 flip_refused <- FALSE
 flip_lp_change <- NA_real_

 # --- ECR-ITERATIVE-1 on the records not flagged as safe matches ------------
 ecr <- if (length(other) == 0L) logical(S) else
  ecr_iterative_two(z, cols = if (has_safe) other else NULL)
 ecr_share <- mean(ecr)

 if (isTRUE(relabel)) {
  # --- majority convention: align, then orient by the majority label -------
  swap <- ecr
  if (count_label2_cpp(z, swap, integer(0)) > as.double(S) * ncol(z) / 2) {
   # exchanging all labels changes the posterior through the
   # component-specific priors only: estimated log ratio of the posterior
   # probabilities of the exchanged and the aligned labelling (0 for
   # identical priors; NA, which leaves the majority convention in force,
   # when it cannot be estimated)
   flip_lp_change <- if (length(exchange_lp) == S) .log_mean_exp(ifelse(swap, -exchange_lp, exchange_lp)) else 0
   if (is.na(flip_lp_change) || flip_lp_change >= -orient_tol) {
    flipped <- TRUE
    if (isTRUE(verbose)) message("Global label swap performed: label 2 dominates label 1.")
    swap <- !swap
   } else {
    flip_refused <- TRUE
   }
  }
 } else {
  # --- labels identified: the draws are reported as sampled ------------------
  swap <- logical(S)
  if (ecr_share > ecr_warn) {
   why <- switch(orientation,
    "safe matches" = "the safe matches identify component 1",
    "match-rate prior" = "the prior on the match probability identifies component 1",
    paste0("the component-specific priors identify component 1 (the sampler refused to exchange ",
           "the labels of its starting state; diagnostics$orient_refused)"))
   warning(sprintf(paste0(
    "The chain visited both labellings of the mixture: %.1f%% of the stored draws allocate ",
    "the records not flagged as safe matches in the opposite way to the modal allocation, so ",
    "the posterior is multimodal or the components are weakly separated. No relabelling was ",
    "applied because %s (see 'Label switching' in ?%s). Consider more informative ",
    "priors, safe matches, more iterations, or starting values."),
    100 * ecr_share, why, help), call. = FALSE)
  }
 }

 if (any(swap)) {
  z <- swap_rows_cpp(z, swap)
  pairs <- lapply(pairs, function(p) .swap_pair(p[[1L]], p[[2L]], swap))
  if (!is.null(theta)) theta[swap] <- 1 - theta[swap]
  if (!is.null(gamma)) {                 # inv_logit(-x) = 1 - inv_logit(x)
   if (is.matrix(gamma)) gamma[swap, ] <- -gamma[swap, , drop = FALSE] else gamma[swap] <- -gamma[swap]
  }
 }

 # average share of the records not flagged as safe allocated to component 1
 share1 <- if (length(other) > 0L) {
  1 - count_label2_cpp(z, logical(S), as.integer(other)) / (as.double(S) * length(other))
 } else NA_real_
 if (is.finite(share1) && share1 < 0.5) {
  if (has_safe && isTRUE(expect_majority)) {
   warning(sprintf(paste0(
    "Only %.0f%% of the records not flagged as safe matches are allocated to component 1, ",
    "which the safe matches identify as the correct-match component, although the prior on ",
    "the match probability does not make correct matches the minority. The safe matches may ",
    "be atypical of the correct matches (e.g. hand-checked pairs with unusual outcomes), so ",
    "that component 1 describes the mismatches. Check the safe matches, or give the expected ",
    "mismatch rate with 'm.rate'."), 100 * share1), call. = FALSE)
  } else if (orientation == "majority component") {
   if (flip_refused) start <- c(start, list(flip_refused = TRUE, flip_lp_change = flip_lp_change))
   warning(.minority_message(share1, start, help), call. = FALSE)
  }
 }

 list(z = z, pairs = pairs, theta = theta, gamma = gamma, swapped = swap,
      orientation = orientation, flipped = flipped, flip_refused = flip_refused,
      flip_lp_change = flip_lp_change, relabelled = isTRUE(relabel),
      ecr_share = ecr_share, share1 = share1)
}
