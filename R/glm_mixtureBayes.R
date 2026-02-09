#' Bayesian Two-Component Mixture GLM (with label-switching adjustment)
#'
#' Fit a two-component Bayesian mixture generalized linear model (GLM) in Stan
#' with families \code{"gaussian"}, \code{"poisson"}, \code{"binomial"}, or \code{"gamma"}.
#' The function builds Stan code, compiles, samples, and then applies a
#' label-switching correction (global majority swap + ECR-ITERATIVE-1) to align
#' component labels across MCMC draws.
#'
#' @param formula A model formula (e.g., \code{y ~ x1 + x2}).
#' @param data A \code{data.frame} (or a named \code{list}) containing the variables
#'   in \code{formula}. Missing values are not allowed.
#' @param family One of \code{"gaussian"}, \code{"poisson"}, \code{"binomial"}, \code{"gamma"}.
#'   Controls the component-specific likelihood.
#' @param priors A named \code{list} (or \code{NULL}) of prior specifications used
#'   by the Stan generator. Any missing entries are filled by
#'   \code{fill_defaults(priors, p_family = family, model_type = "glm")}.
#' @param control A named \code{list} of tuning parameters with defaults:
#'   \itemize{
#'     \item \code{iterations} (default \code{1e4}) total iterations per chain;
#'     \item \code{burnin.iterations} (default \code{1e3}) warm-up iterations;
#'     \item \code{seed} (default random integer);
#'     \item \code{cores} (default \code{getOption("mc.cores", 1L)}).
#'   }
#'   Values in \code{...} override \code{control}.
#' @param ... Optional overrides for elements in \code{control}, e.g.
#'   \code{iterations = 4000}, \code{burnin.iterations = 1000}, \code{seed = 123},
#'   \code{cores = 2}.
#'
#' @details
#' Internally, the function:
#' \enumerate{
#'   \item Validates inputs and merges user priors with weakly-informative defaults;
#'   \item Constructs \code{X}, \code{y} via \code{model.frame} / \code{model.matrix};
#'   \item Generates Stan code using \code{generate_stan(components, priors = priors)};
#'   \item Calls \code{rstan::stan_model()} and \code{rstan::sampling()} (using a single MCMC chain by default to avoid label-switching issues);
#'   \item Extracts posterior draws and performs label-switching adjustment:
#'     a global majority swap if label 2 dominates, then
#'     \code{label.switching::label.switching(method = "ECR-ITERATIVE-1")}.
#' }
#'
#' @return An object of class \code{"glm_mixtureBayes"} containing (at least):
#' \describe{
#'   \item{\code{m_samples}}{Aligned \eqn{z} label matrix (S x N).}
#'   \item{\code{estimates$coefficients}}{Component 1 coefficient draws (S x p).}
#'   \item{\code{estimates$m.coefficients}}{Component 2 coefficient draws (S x p).}
#'   \item{\code{estimates$dispersion}}{Component 1 dispersion (family-specific).}
#'   \item{\code{estimates$m.dispersion}}{Component 2 dispersion (family-specific).}
#'   \item{\code{family}}{The GLM family string.}
#'   \item{\code{call}}{The matched call.}
#' }
#'
#' @section Label switching:
#' We first perform an optional global swap \eqn{(1 \leftrightarrow 2)} if label 2
#' is more frequent overall, then align per-draw labels using
#' \code{ECR-ITERATIVE-1} permutations. Component-specific parameters are permuted
#' accordingly (e.g., \code{beta1}/\code{beta2}, and \code{sigma}/\code{phi}).
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' n  <- 200
#' x1 <- rnorm(n)
#' x2 <- rbinom(n, 1, 0.5)
#' y  <- 0.5 + 1.2*x1 - 0.8*x2 + rnorm(n, sd = 0.7)
#' dat <- data.frame(y, x1, x2)
#'
#' fit <- glm_mixtureBayes(
#'   formula = y ~ x1 + x2,
#'   data    = dat,
#'   family  = "gaussian",
#'   control = list(iterations = 2000, burnin.iterations = 1000, seed = 123, cores = 1)
#' )
#'
#' }
#'
#' @export
#' @importFrom rstan stan_model sampling
#' @importFrom stats model.frame model.response model.matrix
#' @import label.switching


glm_mixtureBayes <- function(formula, data, family = "gaussian", priors = NULL,
                             control = list(iterations = 1e4,
                                            burnin.iterations = 1e3,
                                            seed = sample.int(.Machine$integer.max, 1),
                                            cores = getOption("mc.cores", 1L)),
                             ...){

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

 if(!(family %in% c("gaussian", "poisson", "binomial", "gamma"))){
  stop("Error: the family should be gaussian, poisson, binomial, or gamma")}

  components <- rep(family, 2)

  # build priors (default+user)
  priors <- fill_defaults(priors, p_family = family, model_type = "glm")

  # Prepare the data
  if (!missing(data) && ncol(data) <= 1) {
    stop("Data should have at least one predictor and one response variable.")
  }
  if(missing(data)){data <- list()}

 # Handle the formula & prepare the data
  model_frame <- stats::model.frame(formula, data)
  y <- stats::model.response(model_frame)
  X <- stats::model.matrix(formula, model_frame)
  if(anyNA(X) || anyNA(y)) stop("NA values found in data")

  # stan data
  stan_data <- list(N = nrow(X), K = ncol(X), X = X, y = y)

  # compile & sample
  model_code <- generate_stan(components, priors = priors)
  sm <- rstan::stan_model(model_code = model_code)
  fit <- rstan::sampling(sm, data = stan_data, iter = iterations,
                         warmup = burnin.iterations, chains = 1,
                         seed = seed, cores = cores)  # , refresh = 0

  posterior <- rstan::extract(fit)
  z_samples <- posterior$z            # get z-samples from fit

  ##### Label switching adjustment #####

  # --- 0) Optional global pre-alignment by majority label -----------------------
  # If label 2 is globally more frequent than label 1 across all draws and units,
  # flip all z labels (1 <-> 2) AND swap the corresponding component-specific
  # posterior draws in `posterior` to keep parameters consistent before ECR.
  #
  # This improves stability and reduces the work for the iterative ECR step.
  # Note: we flip only if count(2) > count(1). Since z \in {1,2}, that is
  # equivalent to: sum(z==2) > (S*N)/2.
  count_label2 <- sum(z_samples == 2L, na.rm = TRUE)
  total_labels <- length(z_samples)  # S * N
  if (is.finite(count_label2) && count_label2 > total_labels / 2) {
    message("Global label swap performed: label 2 dominates label 1.")

    # Flip z: 1 -> 2, 2 -> 1
    z_samples <- structure(3L - z_samples, dim = dim(z_samples))

    # Swap component-specific parameters inside `posterior` if they exist
    swap_if_present <- function(lst, a, b) {
      if (all(c(a, b) %in% names(lst))) {
        tmp <- lst[[a]]
        lst[[a]] <- lst[[b]]
        lst[[b]] <- tmp
      }
      lst
    }
    posterior <- swap_if_present(posterior, "beta1", "beta2")
    posterior <- swap_if_present(posterior, "sigma1", "sigma2")  # gaussian
    posterior <- swap_if_present(posterior, "phi1",   "phi2")    # gamma
  }

  # --- 1) Iterative ECR alignment of labels across MCMC draws -------------------
  ls_out <- label.switching::label.switching(
    method = "ECR-ITERATIVE-1",
    z      = z_samples,  # S x N matrix with labels in {1,2}
    K      = 2
  )

  # S x 2 permutation per draw (row i is the permutation for draw i)
  perm <- ls_out$permutations[["ECR-ITERATIVE-1"]]

  # Map z by the inverse permutation for each draw
  map_z <- function(z, perm) {
    if (!is.matrix(z)) stop("`z` must be an S x N matrix.")
    S <- nrow(z); K <- ncol(perm)
    if (K != 2L || nrow(perm) != S) {
      stop("`perm` must be an S x 2 matrix of permutations.")
    }
    for (i in seq_len(S)) {
      inv <- integer(K)
      inv[ perm[i, ] ] <- seq_len(K)   # inverse permutation for draw i
      z[i, ] <- inv[ z[i, ] ]
    }
    z
  }
  z_samples <- map_z(z_samples, perm)

  # --- 2) Helper to permute paired component-specific parameters ----------------
  perm_pair <- function(a1, a2, perm){
    # a1, a2: numeric vector (length S) OR numeric matrix (S x p)
    if (!is.numeric(a1) || !is.numeric(a2)) stop("Inputs must be numeric.")

    # Normalize to S x p matrices
    if (length(dim(a1)) == 1L) {  # vector size
      S <- length(a1); p <- 1L
      A1 <- matrix(a1, ncol = 1L)
      A2 <- matrix(a2, ncol = 1L)
    } else {
      S <- nrow(a1); p <- ncol(a1)
      A1 <- a1; A2 <- a2
    }
    if (!identical(dim(A1), dim(A2))) {
      stop("Shapes differ: ", paste(dim(A1), collapse = "x"), " vs ",
           paste(dim(A2), collapse = "x"))
    }
    if (!is.matrix(perm) || nrow(perm) != S || ncol(perm) != 2L) {
      stop("`perm` must be an S x 2 matrix (one permutation per draw).")
    }

    # Stack to S x 2 x p, then apply permutations per draw
    arr <- array(NA_real_, dim = c(S, 2L, p))
    arr[, 1L, ] <- A1
    arr[, 2L, ] <- A2
    arrp <- label.switching::permute.mcmc(arr, permutations = perm)[[1]]

    # If p == 1, permute.mcmc may return S x 2 (drop dims). Handle both cases.
    if (length(dim(arrp)) == 2L) {  # S x 2
      out1 <- as.numeric(arrp[, 1L])
      out2 <- as.numeric(arrp[, 2L])
      return(list(`1` = out1, `2` = out2))
    }

    # Otherwise S x 2 x p
    out1 <- array(arrp[, 1L, , drop = FALSE], dim = c(S, p))
    out2 <- array(arrp[, 2L, , drop = FALSE], dim = c(S, p))
    if (p == 1L) { out1 <- as.numeric(out1); out2 <- as.numeric(out2) }
    list(`1` = out1, `2` = out2)
  }

  # --- 3) Apply permutations to component-specific parameters -------------------
  # Extract then permute to match the ECR-aligned labels (per-draw permutations)
  beta1.p <- posterior$beta1          # S x p
  beta2.p <- posterior$beta2          # S x p

  tmp <- perm_pair(beta1.p, beta2.p, perm)
  beta1.p <- tmp[[1]]
  beta2.p <- tmp[[2]]

  if (family == "gamma") {
    phi1.p <- posterior$phi1          # S (or S x p if modeled that way)
    phi2.p <- posterior$phi2
    tmp <- perm_pair(phi1.p, phi2.p, perm)
    phi1.p <- tmp[[1]]
    phi2.p <- tmp[[2]]
  }

  if (family == "gaussian") {
    sigma1.p <- posterior$sigma1      # S (or S x p if modeled that way)
    sigma2.p <- posterior$sigma2
    tmp <- perm_pair(sigma1.p, sigma2.p, perm)
    sigma1.p <- tmp[[1]]
    sigma2.p <- tmp[[2]]
  }


  if (is.null(posterior$z) || is.null(posterior$beta1) || is.null(posterior$beta2)) {
    stop("Error: Missing expected parameters in posterior")
  }

  # Calculate statistics --------
  terms <- all.vars(formula[-2])
  if(ncol(beta1.p) > length(terms)){
    colnames(beta1.p) <- c("(Intercept)", terms)
    colnames(beta2.p) <- c("(Intercept)", terms)
  } else{
    colnames(beta1.p) <- terms
    colnames(beta2.p) <- terms
  }

  x <- list(m_samples = z_samples,
            estimates = list(coefficients = beta1.p,
                             m.coefficients = beta2.p))
  if (family == "gamma") {
    x$estimates$dispersion <- 1/phi1.p
    x$estimates$m.dispersion <- 1/phi2.p
  }
  if (family == "gaussian") {
    x$estimates$dispersion <- sigma1.p^2
    x$estimates$m.dispersion <- sigma2.p^2
  }

  x <- c(x, family = family, call = match.call())

  class(x) <- "glm_mixtureBayes"
  return(x)
}
