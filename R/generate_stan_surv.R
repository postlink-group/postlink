#' Generate Stan code for a two-component survival mixture 
#'
#' Creates a complete Stan program as a single string for a two-component mixture survival model 
#' with both components of the same family (Gamma or Weibull). The Stan code will include user-defined 
#' priors and any hyperparameters.
#'
#' @param components A length-2 character vector specifying the component distributions. 
#'   Only \code{c("gamma","gamma")} or \code{c("weibull","weibull")} are supported.
#' @param priors A named list of prior distribution strings (as in \code{survreg_mixtureBayes}). 
#'   Defaults will be used for any missing entries.
#' @return A single string containing the Stan program.
#' @keywords internal
#' @noRd
generate_stan_surv <- function(components, priors = list()) {
  # Get any function and data block definitions needed for custom priors
  defs <- get_stan_definitions(priors)
  function_definitions  <- defs$function_defs
  variable_definitions  <- defs$variable_defs
  
  if (identical(components, c("gamma", "gamma"))) {
    # Stan code for Gamma-Gamma mixture survival model
    stan_code <- paste(
      "functions {",
      function_definitions,
      "}",
      "data {",
      "  int<lower=1> N;              // number of observations",
      "  int<lower=1> K;              // number of predictors",
      "  matrix[N, K] X;              // design matrix (including intercept as first column)",
      "  vector<lower=0>[N] y;        // observed survival times (non-negative)",
      "}",
      "transformed data {",
      variable_definitions,
      "}",
      "parameters {",
      "  real<lower=0, upper=1> theta; // mixture weight (Pr(component 1))",
      "  vector[K] beta1;               // regression coefficients for component 1 (linear predictor)",
      "  vector[K] beta2;               // regression coefficients for component 2",
      "  real<lower=0> phi1;            // shape parameter (Gamma) for component 1",
      "  real<lower=0> phi2;            // shape parameter for component 2",
      "}",
      "model {",
      "  vector[N] mu1 = exp(X * beta1);   // mean of component 1 (mu = exp(Xb))",
      "  vector[N] mu2 = exp(X * beta2);   // mean of component 2",
      "  // Log-likelihood contributions",
      "  for (n in 1:N) {",
      "    // Each component's log-likelihood for observation n",
      "    real log_lik1 = gamma_lpdf(y[n] | phi1, phi1 / mu1[n]);",
      "    real log_lik2 = gamma_lpdf(y[n] | phi2, phi2 / mu2[n]);",
      "    // Mixture log-likelihood via log_sum_exp",
      "    target += log_mix(theta, log_lik1, log_lik2);",
      "  }",
      "  // Priors:",
      "  beta1 ~ ", priors$beta1, ";",
      "  beta2 ~ ", priors$beta2, ";",
      "  theta ~ ", priors$theta, ";",
      "  // Priors for shape parameters",
      "  phi1 ~ ", priors$phi1, ";",
      "  phi2 ~ ", priors$phi2, ";",
      "}",
      "generated quantities {",
      "  array[N] int<lower=1, upper=2> z;  // latent mixture labels",
      "  // Draw component labels for each observation",
      "  vector[N] mu1_gen = exp(X * beta1);",
      "  vector[N] mu2_gen = exp(X * beta2);",
      "  for (n in 1:N) {",
      "    real log_prob1 = log(theta) + gamma_lpdf(y[n] | phi1, phi1 / mu1_gen[n]);",
      "    real log_prob2 = log1m(theta) + gamma_lpdf(y[n] | phi2, phi2 / mu2_gen[n]);",
      "    // Sample label using unnormalized log-probs",
      "    z[n] = categorical_rng( softmax([log_prob1, log_prob2]') );",
      "  }",
      "}",
      sep = "\n"
    )
    return(stan_code)
    
  } else if (identical(components, c("weibull", "weibull"))) {
    # Stan code for Weibull-Weibull mixture survival model
    stan_code <- paste(
      "functions {",
      function_definitions,
      "}",
      "data {",
      "  int<lower=1> N; ",
      "  int<lower=1> K; ",
      "  matrix[N, K] X; ",
      "  vector<lower=0>[N] y; ",
      "  array[N] int<lower=0, upper=1> status; // 1 if event observed, 0 if right-censored",
      "}",
      "transformed data {",
      variable_definitions,
      "}",
      "parameters {",
      "  vector[K] beta1; ",
      "  vector[K] beta2; ",
      "  real<lower=0> shape1; ",
      "  real<lower=0> shape2; ",
      "  real<lower=0> scale1; ",
      "  real<lower=0> scale2; ",
      "  real<lower=0, upper=1> theta; ",
      "}",
      "model {",
      "  for (n in 1:N) {",
      "    // Compute linear predictors and scale factors for each component",
      "    real linpred1 = X[n] * beta1; ",
      "    real linpred2 = X[n] * beta2; ",
      "    real scale_factor1 = scale1 * exp(-linpred1); ",
      "    real scale_factor2 = scale2 * exp(-linpred2); ",
      "    // Contribution to log likelihood depends on status (event or censored)",
      "    real log_p1 = log(theta); ",
      "    real log_p2 = log1m(theta); ",
      "    if(status[n] == 1) {",
      "      log_p1 += weibull_lpdf(y[n] | shape1, scale_factor1); ",
      "      log_p2 += weibull_lpdf(y[n] | shape2, scale_factor2); ",
      "    } else {",
      "      log_p1 += weibull_lccdf(y[n] | shape1, scale_factor1); ",
      "      log_p2 += weibull_lccdf(y[n] | shape2, scale_factor2); ",
      "    }",
      "    target += log_sum_exp(log_p1, log_p2); ",
      "  }",
      "  // Priors:",
      "  beta1 ~ ", priors$beta1, ";",
      "  beta2 ~ ", priors$beta2, ";",
      "  shape1 ~ ", priors$shape1, ";",
      "  shape2 ~ ", priors$shape2, ";",
      "  scale1 ~ ", priors$scale1, ";",
      "  scale2 ~ ", priors$scale2, ";",
      "  theta  ~ ", priors$theta, ";",
      "}",
      "generated quantities {",
      "  array[N] int<lower=1, upper=2> z; ",
      "  vector[N] log_lik; ",
      "  // Draw latent labels and compute log-likelihood for WAIC/LOO if needed",
      "  for(n in 1:N) {",
      "    vector[2] log_p; ",
      "    real linpred1 = X[n] * beta1; ",
      "    real linpred2 = X[n] * beta2; ",
      "    real scale_factor1 = scale1 * exp(-linpred1); ",
      "    real scale_factor2 = scale2 * exp(-linpred2); ",
      "    // Use event likelihood for both (for label sampling, assume event distribution)",
      "    log_p[1] = log(theta) + weibull_lpdf(y[n] | shape1, scale_factor1); ",
      "    log_p[2] = log1m(theta) + weibull_lpdf(y[n] | shape2, scale_factor2); ",
      "    z[n] = categorical_rng( softmax(log_p) ); ",
      "    log_lik[n] = log_sum_exp(log_p[1], log_p[2]); ",
      "  }",
      "}",
      sep = "\n"
    )
    return(stan_code)
  } else {
    stop("Unsupported components; must be two of 'gamma' or two of 'weibull'.")
  }
}
