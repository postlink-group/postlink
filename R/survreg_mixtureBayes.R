#' Bayesian Two-Component Survival Mixture Model (Gamma or Weibull)
#'
#' Fit a two-component Bayesian parametric survival regression model in Stan
#' with either Gamma or Weibull component distributions. Supports right-censored data.
#' This function builds the Stan code with specified prior distributions, runs MCMC sampling,
#' and corrects label switching in the posterior draws.
#'
#' Note: The response variable must be a \code{Surv} object (from the \pkg{survival} package), and any data with missing values will be rejected by this function.
#'
#' @param formula Survival formula, e.g. \code{Surv(time, status) ~ x1 + x2}. The response must be a \code{Surv} object.
#' @param data A data frame (or named list) containing the variables in the formula.
#'   Must have no missing values in predictors or survival outcome.
#' @param family Distribution family for components: either \code{"gamma"} or \code{"weibull"}.
#'   (Both mixture components use the same family.)
#' @param priors A named list specifying the prior distributions for parameters (optional).
#'   Names should correspond to parameters:
#'   \itemize{
#'     \item For \code{family="gamma"}: \code{beta1}, \code{beta2} (coefficients);
#'           \code{phi1}, \code{phi2} (shape parameters of the two components); \code{theta} (mixing weight).
#'     \item For \code{family="weibull"}: \code{beta1}, \code{beta2} (coefficients);
#'           \code{shape1}, \code{shape2} (Weibull shape parameters);
#'           \code{scale1}, \code{scale2} (Weibull scale parameters); \code{theta} (mixing weight).
#'   }
#'   Prior values should be strings of valid Stan distribution right-hand-sides
#'   (e.g. \code{"normal(0,5)"}).
#'   Omit or set \code{priors=NULL} to use default weakly-informative priors (which are filled in automatically).
#' @param control A list of MCMC control parameters. Supported entries:
#'   \code{iterations} (total iterations per chain, default 10000),
#'   \code{burnin.iterations} (warm-up iterations, default 1000),
#'   \code{seed} (random seed, default a random integer),
#'   \code{cores} (parallel cores for Stan, default uses \code{getOption("mc.cores",1)}).
#'
#' @param ... Additional arguments to override elements of \code{control}. For example,
#'   \code{iterations=2000}, \code{seed=123}.
#'
#' @details
#' The model assumes two latent groups. For each observation, if $z_i=1$, the outcome follows the specified survival distribution with parameters $(\beta_1, \phi_1)$ or $(\beta_1, \text{shape}_1,\text{scale}_1)$; if $z_i=2$, it follows the distribution with $(\beta_2, \phi_2)$ or $(\beta_2, \text{shape}_2,\text{scale}_2)$. The mixture weight $\theta = P(z_i=1)$ is given a prior (default $\text{Beta}(1,1)$). Covariates enter the model through a linear predictor: for Gamma, $\mu_{j}(i) = \exp(X_i \beta_j)$ is the mean for component $j$ and shape $\phi_j$ defines variability; for Weibull, the scale parameter for component $j$ is modulated by $\exp(-X_i \beta_j)$ (accelerated failure time model form).
#'
#' Internally, the function will:
#' \enumerate{
#'   \item Validate inputs and fill in default priors via \code{fill_defaults(priors, p_family=family, model_type="survival")}.
#'   \item Construct the design matrix \code{X} and outcome vectors from the formula and data.
#'   \item Generate Stan code with \code{generate_stan_surv(components, priors)} and compile it.
#'   \item Run MCMC sampling using \code{rstan::sampling} (with a single chain by default to avoid label-switching issues).
#'   \item Extract posterior draws and perform label-switching correction: a global swap (if needed) followed by the ECR algorithm for per-iteration label alignment (using \code{label.switching::label.switching}).
#' }
#'
#' @return An object of class \code{"surv_mixtureBayes"}, which is a list with elements:
#' \describe{
#'   \item{\code{m_samples}}{Matrix of dimension (S x N) of latent class labels for each draw (after label alignment), where S is the number of post-warmup draws and N is the number of observations.}
#'   \item{\code{estimates}}{A list of posterior draws for each parameter:
#'         \code{coefficients} (S x p matrix for component 1 coefficients);
#'         \code{m.coefficients} (S x p matrix for component 2 coefficients);
#'         \code{phi1}/\code{phi2} (Gamma shape parameters) or \code{shape1}/\code{shape2} (Weibull shape parameters);
#'         \code{scale1}/\code{scale2} (Weibull scale parameters);
#'         \code{theta} (mixing weight).}
#'   \item{\code{family}}{The distribution family (\code{"gamma"} or \code{"weibull"}).}
#'   \item{\code{call}}{The matched call.}
#' }
#'
#' @examples
#' \dontrun{
#' library(survival)
#' set.seed(123)
#' n <- 200
#' x1 <- rnorm(n); x2 <- rbinom(n,1,0.5)
#' # Simulate survival times from two latent subpopulations
#' z <- rbinom(n, 1, 0.4) + 1  # latent group indicators (1 or 2)
#' T1 <- rgamma(n, shape=2, rate=0.1)        # component 1 event times
#' T2 <- rweibull(n, shape=3, scale=2)       # component 2 event times
#' Tobs <- ifelse(z==1, T1, T2)
#' Censoring <- rexp(n, rate=0.02)
#' time <- pmin(Tobs, Censoring)            # observed time is min(event, censor)
#' status <- as.integer(Tobs <= Censoring)  # status=1 if event observed
#' dat <- data.frame(time, status, x1, x2)
#'
#' fit <- survreg_mixtureBayes(Surv(time, status) ~ x1 + x2, data=dat, family="weibull",
#'                             control=list(iterations=2000, burnin.iterations=1000, seed=101))
#' }
#'
#' @export
#' @importFrom stats model.frame model.response model.matrix
#' @importFrom rstan stan_model sampling
#' @import label.switching

survreg_mixtureBayes <- function(formula, data, family = "gamma", priors = NULL,
                                 control = list(iterations = 1e4,
                                                burnin.iterations = 1e3,
                                                seed = sample.int(.Machine$integer.max, 1),
                                                cores = getOption("mc.cores", 1L)),
                                 ...) {

  dcontrols <- list(...)
if ("iterations" %in% names(dcontrols)){
  iterations <- dcontrols$iterations
 } else {
  iterations <- ifelse("iterations" %in% names(control),
                       control$iterations, 1e4)
 }
 if ("burnin.iterations" %in% names(dcontrols)){
  burnin.iterations <- dcontrols$burnin.iterations
 } else {
  burnin.iterations <- ifelse("burnin.iterations" %in% names(control),
                              control$burnin.iterations, 1e3)
 }

 if ("seed" %in% names(dcontrols)){
  seed <- dcontrols$seed
 } else {
  seed <- ifelse("seed" %in% names(control),
                 control$seed, sample.int(.Machine$integer.max, 1))
 }

 if ("cores" %in% names(dcontrols)){
  cores <- dcontrols$cores
 } else {
  cores <- ifelse("cores" %in% names(control),
                  control$cores, getOption("mc.cores", 1))
 }

if(missing(formula)){ stop("Error: a formula for the outcome
                            model is required")}
 if(!inherits(formula, "formula")){ stop("Error: formula should be
                                         a formula object")}

 if(!missing(data) && (!is.data.frame(data) && !is.list(data))){
  stop("Error: data should be a data.frame or list")}
  if(!(family %in% c("gamma", "weibull"))) {
    stop("Error: family must be either 'gamma' or 'weibull'.")}

    components <- rep(family, 2) # e.g. c("weibull","weibull")

  # Fill priors with defaults for survival model if needed
  priors <- fill_defaults(priors, p_family = family, model_type = "survival")

  # Prepare design matrix and response
  if (!missing(data) && ncol(data) <= 1) {
    stop("Data should contain at least one predictor in addition to the survival outcome.")
  }
  if(missing(data)) data <- list()  # model.frame will look in environment if data not provided

  model_frame <- model.frame(formula, data = data)  # will drop NAs automatically
  # Check that response is a survival object
  if (!survival::is.Surv(model_frame[[1]])) {
    stop("Error: the response must be a 'Surv' survival object (e.g. Surv(time, status)).")
  }
  # Only right-censoring supported
  if(attr(model_frame[[1]], "type") != "right") {
    stop("Error: only right-censoring (type='right') is currently supported.")
  }
  # Extract survival response data
  y <- model_frame[[1]][, "time"]    # event or censoring times
  status <- model_frame[[1]][, "status"]  # event status (1=event, 0=censored)
  # Construct model matrix of covariates (include intercept)
  X <- cbind(1, model_frame[ , -1, drop = FALSE])  # intercept + covariates
  X <- unique(X, MARGIN = 2)             # remove duplicate columns (e.g. if intercept was already present)

  if(anyNA(X) || anyNA(y) || anyNA(status)) {
    stop("Error: NA values found in the data (missing values are not allowed).")
  }

  # Prepare data list for Stan
  stan_data <- list(N = nrow(X), K = ncol(X), X = X, y = as.vector(y), status = as.integer(status))

  # Generate Stan program code for the specified mixture
  model_code <- generate_stan_surv(components, priors = priors)
  # Compile Stan model
  stan_model <- rstan::stan_model(model_code = model_code)


  # Run MCMC sampling (single chain by default, user can increase 'cores' for parallel chains)
  fit <- rstan::sampling(stan_model, data = stan_data, iter = iterations,
                         warmup = burnin.iterations, chains = 1,
                         seed = seed, cores = cores)
  posterior <- rstan::extract(fit)

  # Extract MCMC draws of parameters
  z_samples <- posterior$z         # S x N matrix of latent labels (Stan stores z as 1/2 values)
  beta1.p   <- posterior$beta1    # S x K matrix of beta1 draws
  beta2.p   <- posterior$beta2    # S x K matrix of beta2 draws

  # Family-specific parameter draws
  if(family == "gamma") {
    phi1.p <- posterior$phi1      # S-length vector of phi1 draws
    phi2.p <- posterior$phi2      # S-length vector of phi2 draws
  } else if(family == "weibull") {
    shape1.p <- posterior$shape1  # S-length vector of shape1 draws
    shape2.p <- posterior$shape2  # S-length vector of shape2 draws
    scale1.p <- posterior$scale1  # S-length vector of scale1 draws
    scale2.p <- posterior$scale2  # S-length vector of scale2 draws
  }

  ##### 1) Global label swap if component 2 dominates overall #####
  # Count how many total labels == 2 across all draws and observations:
  total_labels <- length(z_samples)  # S * N
  count_label2 <- sum(z_samples == 2, na.rm = TRUE)
  if(is.finite(count_label2) && count_label2 > total_labels/2) {
    message("Global label swap performed: component 2 was more frequent than component 1 overall.")
    # Flip all z labels: (1 -> 2, 2 -> 1)
    z_samples <- 3L - z_samples  # since labels are 1/2, 3 - current swaps them
    # Swap component-specific parameter draws for all iterations
    temp <- beta1.p;  beta1.p <- beta2.p;  beta2.p <- temp
    if(family == "gamma") {
      temp <- phi1.p;  phi1.p <- phi2.p;  phi2.p <- temp
    }
    if(family == "weibull") {
      temp <- shape1.p;  shape1.p <- shape2.p;  shape2.p <- temp
      temp <- scale1.p;  scale1.p <- scale2.p;  scale2.p <- temp
    }
  }

  ##### 2) Iterative ECR algorithm for per-iteration label alignment #####
  ls_out <- label.switching::label.switching(method = "ECR-ITERATIVE-1",
                                             z = z_samples, K = 2)
  perm <- ls_out$permutations[["ECR-ITERATIVE-1"]]  # S x 2 matrix of optimal label mappings for each draw
  # Apply the label permutation to each draw's labels and parameters
  S <- nrow(perm)
  for(i in seq_len(S)) {
    # perm[i,] is like c(1,2) (identity) or c(2,1) (swap); if first element is 2, that implies a swap of labels 1 and 2
    if(perm[i, 1] == 2) {
      # Swap labels in this draw's z
      z_samples[i, ] <- 3L - z_samples[i, ]
      # Swap this draw's component-specific parameter values
      # (since beta1.p, beta2.p are matrices S x K, we swap the i-th rows)
      beta1_row <- beta1.p[i, , drop = FALSE]
      beta1.p[i, ] <- beta2.p[i, ]
      beta2.p[i, ] <- beta1_row
      if(family == "gamma") {
        # swap phi values in draw i
        temp_phi <- phi1.p[i]
        phi1.p[i] <- phi2.p[i]
        phi2.p[i] <- temp_phi
      }
      if(family == "weibull") {
        temp_shape <- shape1.p[i]
        shape1.p[i] <- shape2.p[i]
        shape2.p[i] <- temp_shape
        temp_scale <- scale1.p[i]
        scale1.p[i] <- scale2.p[i]
        scale2.p[i] <- temp_scale
      }
    }
  }

  # Construct return object
  out <- list()
  out$m_samples <- z_samples  # aligned membership matrix (S x N)
  # Posterior draws organized in list "estimates"
  out$estimates <- list(
    coefficients    = beta1.p,
    m.coefficients  = beta2.p
  )
  if(family == "gamma") {
    out$estimates$shape    = phi1.p    # phi (shape) draws for component 1
    out$estimates$m.shape  = phi2.p    # phi (shape) draws for component 2
  }
  if(family == "weibull") {
    out$estimates$shape    = shape1.p  # shape parameter draws (comp 1)
    out$estimates$m.shape  = shape2.p  # shape parameter draws (comp 2)
    out$estimates$scale    = scale1.p  # scale parameter draws (comp 1)
    out$estimates$m.scale  = scale2.p  # scale parameter draws (comp 2)
  }
  out$family <- family
  out$call <- match.call()
  class(out) <- "surv_mixtureBayes"
  return(out)
}
