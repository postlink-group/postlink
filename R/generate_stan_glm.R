#' Generate Stan code for two-component GLM mixtures
#'
#' Build a complete Stan program (as a single character string) for a
#' two-component mixture model with one of the supported GLM families:
#' \code{"gaussian"}, \code{"poisson"}, \code{"gamma"}, or \code{"binomial"}.
#' The function stitches together family-specific Stan blocks and injects
#' user-specified prior strings and hyperparameter definitions.
#'
#' @param components A character vector of length 2 giving the component
#'   likelihoods; each element must be one of
#'   \code{c("gaussian","poisson","gamma","binomial")}.
#'   Current implementation supports only the symmetric pairs
#'   \code{c("gaussian","gaussian")}, \code{c("poisson","poisson")},
#'   \code{c("gamma","gamma")}, and \code{c("binomial","binomial")}.
#' @param priors A named \code{list} of prior strings and (optionally)
#'   hyperparameter definitions that appear in the Stan code. Prior strings
#'   should be valid Stan sampling statements' right-hand sides (e.g.,
#'   \code{"normal(0,1)"} or \code{"multi_normal(mu, Sigma)"}). Expected names
#'   depend on the family:
#'   \itemize{
#'     \item Gaussian: \code{beta1}, \code{beta2}, \code{sigma1}, \code{sigma2}, \code{theta}
#'     \item Poisson:  \code{beta1}, \code{beta2}, \code{theta}
#'     \item Gamma:    \code{beta1}, \code{beta2}, \code{phi1}, \code{phi2}, \code{theta}
#'     \item Binomial: \code{beta1}, \code{beta2}, \code{theta}
#'   }
#'   Any symbols referenced inside these strings (e.g., hyperparameters like
#'   \code{mu}, \code{Sigma}, \code{phi_loc}, \code{phi_scale}) must be defined
#'   via \code{get_stan_definitions(priors)} in the \code{transformed data}
#'   block (see Details).
#'
#' @details
#' Internally, the function calls \code{get_stan_definitions(priors)} (see
#' \code{priors_helpers.R}) to obtain two text fragments:
#' \itemize{
#'   \item \strong{\code{function_defs}}: any user-defined Stan functions placed in the \code{functions} block;
#'   \item \strong{\code{variable_defs}}: declarations/assignments for hyperparameters placed in the \code{transformed data} block.
#' }
#' The returned Stan program contains:
#' \enumerate{
#'   \item a \code{functions} block (if provided by the helper);
#'   \item \code{data} and \code{parameters} blocks specific to each family;
#'   \item a \code{model} block that encodes the mixture via \code{log_sum_exp} or
#'         \code{log_mix} and applies the injected prior strings;
#'   \item a \code{generated quantities} block that samples per-observation
#'         mixture memberships \code{z[n]} using \code{categorical_rng(softmax(lw))}.
#' }
#' For Poisson and Binomial families, the linear predictor is \code{eta = X * beta}
#' with canonical links (log for Poisson, logit for Binomial). For Gamma, the
#' parameterization uses shape \code{phi} and rate \code{phi/exp(eta)} so that
#' \code{E[y] = exp(eta)}.
#'
#' @return A single character string containing a complete Stan program tailored
#'   to the requested \code{components} and \code{priors}.
#'
#' @section Expected data shapes in Stan:
#' The generated programs expect:
#' \itemize{
#'   \item \code{N}: number of observations (\code{int});
#'   \item \code{K}: number of predictors (\code{int});
#'   \item \code{X}: predictor matrix (\code{matrix[N, K]});
#'   \item \code{y}: response (\code{vector[N]} for Gaussian/Gamma,
#'         \code{int[N]} nonnegative counts for Poisson, \code{int[N]} in \{0,1\} for Binomial).
#' }
#'
#' @keywords internal
#' @noRd
generate_stan <- function(components, priors = list()) {
  
  # helper defs from priors_helpers.R
  defs <- get_stan_definitions(priors)
  function_definitions <- defs$function_defs
  variable_definitions <- defs$variable_defs

  # checks that inputs are as expected
  if (identical(components, c("gaussian","gaussian"))) {
    # Fixed Stan code for linear-linear mixture
    stan_code <- paste(
      "functions {",
      function_definitions,
      "}",
      "data {",
      "  int<lower=1> N;             // Number of data points",
      "  int<lower=1> K;             // Number of predictors",
      "  matrix[N, K] X;             // Predictor matrix",
      "  vector[N] y;                // Response vector",
      "}",
      "transformed data {",
      variable_definitions, # defines any hyperparameter variables used in prior string
      "}",
      "parameters {",
      "  real<lower=0, upper=1> theta; // Mixture weight for the first component",
      "  real<lower=0> sigma1;         // Standard deviation of the first component",
      "  real<lower=0> sigma2;         // Standard deviation of the second component",
      "  vector[K] beta1;              // Regression coefficients for the first component",
      "  vector[K] beta2;              // Regression coefficients for the second component",
      "}",
      "model {",
      "  // priors",
      "  sigma1 ~ ", priors$sigma1, ";",
      "  sigma2 ~ ", priors$sigma2, ";",
      "  beta1 ~ ", priors$beta1, ";",
      "  beta2 ~ ", priors$beta2, ";",
      "  theta ~ ", priors$theta, ";",
      "  // Mixture model likelihood",
      "  for (n in 1:N) {",
      "    target += log_sum_exp(",
      "      log(theta) + normal_lpdf(y[n] | dot_product(X[n], beta1), sigma1),",
      "      log1m(theta) + normal_lpdf(y[n] | dot_product(X[n], beta2), sigma2)",
      "    );",
      "  }",
      "}",
      "generated quantities {",
      "  int<lower=1, upper=2> z[N];    // Mixture membership",
      "  for (n in 1:N) {",
      "    // Calculate unnormalized log probabilities for each component",
      "    vector[2] lw;",
      "    lw[1] = log(theta) + normal_lpdf(y[n] | dot_product(X[n], beta1), sigma1);",
      "    lw[2] = log1m(theta) + normal_lpdf(y[n] | dot_product(X[n], beta2), sigma2);",
      "    ",
      "    // Normalize probabilities using softmax",
      "    ",
      "    // Sample z[n] based on the posterior probabilities",
      "    z[n] = categorical_rng(softmax(lw));",
      "  }",
      "}",
      sep = "\n"
    )

    return(stan_code)
  
 } else if (identical(components, c("poisson", "poisson"))) {
  # Poisson-Poisson mixture code
    stan_code <- paste(
      "functions {",
      function_definitions,
      "}",
      "data {",
      "  int<lower=1> N;             // Number of observations",
      "  int<lower=0> y[N];          // Poisson response variable (counts)",
      "  int<lower=1> K;             // Number of predictors",
      "  matrix[N, K] X;             // Predictor matrix",
      "}",
      "",
      "transformed data {",
      variable_definitions,
      "}",
      "parameters {",
      "  real<lower=0, upper=1> theta;           // Mixing proportions (constrained to sum to 1)",
      "  vector[K] beta1;            // Regression coefficients for component 1",
      "  vector[K] beta2;            // Regression coefficients for component 2",
      "}",
      "",
      "model {",
      "  vector[N] log_lik1;         // Log-likelihood for component 1",
      "  vector[N] log_lik2;         // Log-likelihood for component 2",
      "  ",
      "  // Linear predictors for each component",
      "  vector[N] eta1 = X * beta1; // Linear predictor for component 1",
      "  vector[N] eta2 = X * beta2; // Linear predictor for component 2",
      "",
      "  // Calculate log-likelihoods for each component",
      "  for (n in 1:N) {",
      "    log_lik1[n] = poisson_log_lpmf(y[n] | eta1[n]); // Component 1",
      "    log_lik2[n] = poisson_log_lpmf(y[n] | eta2[n]); // Component 2",
      "    target += log_mix(theta, log_lik1[n], log_lik2[n]);",
      "  }",
      "",
      "  // Priors for regression coefficients",
      "  beta1 ~ ", priors$beta1, ";",
      "  beta2 ~ ", priors$beta2, ";",
      "  theta ~ ", priors$theta, ";",
      "}",
      "",
      "generated quantities {",
      "  int<lower=1, upper=2> z[N];    // Mixture membership",
      "  for (n in 1:N) {",
      "    // Calculate unnormalized log probabilities for each component",
      "    vector[2] lw;",
      "    lw[1] = log(theta) + poisson_log_lpmf(y[n] | dot_product(X[n], beta1));",
      "    lw[2] = log1m(theta) + poisson_log_lpmf(y[n] | dot_product(X[n], beta2));",
      "    // Normalize probabilities using softmax",
      "    // Sample z[n] based on the posterior probabilities",
      "    z[n] = categorical_rng(softmax(lw));",
      "  }",
      "}",
      sep = "\n"
    )
    return(stan_code)
  
  } else if (identical(components, c("gamma","gamma"))) {
    stan_code <- paste(
      "functions {",
      function_definitions,
      "}",
      "data {",
      "  int<lower=1> N;               // Number of observations",
      "  int<lower=1> K;               // Number of predictors",
      "  vector<lower=0>[N] y;         // Response variable (positive values)",
      "  matrix[N, K] X;               // Predictor matrix",
      "}",
      "",
      "transformed data {",
      variable_definitions,
      "}",
      "parameters {",
      "  real<lower=0, upper=1> theta; // Mixing proportions (must sum to 1)",
      "  vector[K] beta1;              // Regression coefficients for component 1",
      "  vector[K] beta2;              // Regression coefficients for component 2",
      "  real<lower=0> phi1;           // Shape parameter for component 1",
      "  real<lower=0> phi2;           // Shape parameter for component 2",
      "}",
      "",
      "model {",
      "  vector[N] eta1 = X * beta1;  // eta for component 1",
      "  vector[N] eta2 = X * beta2;  // eta for component 2",
      "  vector[N] log_lik1;",
      "  vector[N] log_lik2;",
      "",
      "  // Calculate log-likelihoods for each component",
      "  // likelihood: shape=phi, rate=phi/exp(eta)",
      "  for (n in 1:N) {",
      "    log_lik1[n] = gamma_lpdf(y[n] | phi1, phi1 * exp(-eta1[n]));",
      "    log_lik2[n] = gamma_lpdf(y[n] | phi2, phi2 * exp(-eta2[n]));",
      "    target += log_mix(theta, log_lik1[n], log_lik2[n]);",
      "  }",
      "",
      "  // Priors for regression coefficients and mix proportion",
      "  beta1 ~ ", priors$beta1, ";",
      "  beta2 ~ ", priors$beta2, ";",
      "  // Priors for shape parameters",
      "  phi1 ~ ", priors$phi1, ";",
      "  phi2 ~ ", priors$phi2, ";",
      "  theta ~ ", priors$theta, ";",
      "}",
      "generated quantities {",
      "  int<lower=1, upper=2> z[N];    // Mixture membership",
      "  vector[N] eta1 = X * beta1;  // Recompute eta1 for generated quantities",
      "  vector[N] eta2 = X * beta2;  // Recompute eta2 for generated quantities",
      "  for (n in 1:N) {",
      "    // Calculate unnormalized log probabilities for each component",
      "    vector[2] lw;",
      "    lw[1] = log(theta) + gamma_lpdf(y[n] | phi1, phi1 * exp(-eta1[n]));",
      "    lw[2] = log1m(theta) + gamma_lpdf(y[n] | phi2, phi2 * exp(-eta2[n]));",
      "    // Normalize probabilities using softmax",
      "    // Sample z[n] based on the posterior probabilities",
      "    z[n] = categorical_rng(softmax(lw));",
      "  }",
      "}",
      sep = "\n"
    )
    return(stan_code)

    } else if (identical(components, c("binomial","binomial"))) {
      # Logistic-logistic mixture Stan code
    stan_code <- paste(
      "functions {",
      function_definitions,
      "}",
      "data {",
      "  int<lower=1> N;           // Number of observations",
      "  int<lower=1> K;           // Number of predictors",
      "  matrix[N, K] X;           // Predictor matrix",
      "  int<lower=0, upper=1> y[N]; // Binary outcome",
      "}",
      "",
      "transformed data {",
      variable_definitions,
      "}",
      "parameters {",
      "  real<lower=0, upper=1> theta;     // Mixing proportion (for component 1)",
      "  vector[K] beta1;                 // Regression coefficients for component 1",
      "  vector[K] beta2;                 // Regression coefficients for component 2",
      "}",
      "",
      "model {",
      "  vector[N] eta1 = X * beta1;  // eta1 from component 1",
      "  vector[N] eta2 = X * beta2;  // eta2 from component 2",
      "  vector[N] log_lik1;  // Log-likelihood contributions from component 1",
      "  vector[N] log_lik2;  // Log-likelihood contributions from component 2",
      "  // priors",
      "  beta1 ~ ", priors$beta1, ";",
      "  beta2 ~ ", priors$beta2, ";",
      "  theta ~ ", priors$theta, ";",
      "",
      "  // Mixture model likelihood",
      "  for (n in 1:N) {",
      "    log_lik1[n] = bernoulli_logit_lpmf(y[n] | eta1[n]);",
      "    log_lik2[n] = bernoulli_logit_lpmf(y[n] | eta2[n]);",
      "    target += log_mix(theta, log_lik1[n], log_lik2[n]);",
      "  }",
      "}",
      "generated quantities {",
      "  int<lower=1, upper=2> z[N];      // Mixture membership for each observation",
      "  for (n in 1:N) {",
      "    vector[2] lw;",
      "    lw[1] = log(theta) + bernoulli_logit_lpmf(y[n] | dot_product(X[n], beta1));",
      "    lw[2] = log1m(theta) + bernoulli_logit_lpmf(y[n] | dot_product(X[n], beta2));",
      "    z[n] = categorical_rng(softmax(lw)); // Sample membership",
      "  }",
      "}",
      sep = "\n"
    )
    return(stan_code)
    
  } else {
    stop("Invalid mixture inputs. Must be: gaussian, poisson, gamma, or binomial")
  }
}
