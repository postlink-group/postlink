# NOTE:
# This file depends on helpers defined in:
#   - generate_stan_glm.R    (generate_stan)
#   - priors_helpers_bayes.R (fill_defaults, get_stan_definitions, ...)
# These helpers must be collated/loaded before this file.

#' Bayesian Two-Component Mixture GLM (with label-switching adjustment)
#'
#' Fit a two-component Bayesian mixture generalized linear model (GLM) in Stan
#' with families \code{"gaussian"}, \code{"poisson"}, \code{"binomial"}, or \code{"gamma"}.
#' This implementation function assumes upstream wrappers have already handled the
#' formula/data interface and produced a design matrix \code{X} and response \code{y}.
#'
#' The function builds Stan code, compiles, samples, and then applies a
#' label-switching correction (global majority swap + ECR-ITERATIVE-1) to align
#' component labels across MCMC draws.
#'
#' @param X A numeric design matrix (N x K), typically from \code{model.matrix()}
#'   upstream. Missing values are not allowed.
#' @param y A response vector of length N. For \code{"gaussian"} and \code{"gamma"},
#'   \code{y} should be numeric; for \code{"poisson"} and \code{"binomial"}, \code{y}
#'   should be integer-valued (or coercible without loss).
#' @param family One of \code{"gaussian"}, \code{"poisson"}, \code{"binomial"}, or \code{"gamma"}.
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
#' @return An object of class \code{"glmMixBayes"} containing (at least):
#' \describe{
#'   \item{\code{m_samples}}{Aligned \eqn{z} label matrix (S x N).}
#'   \item{\code{estimates$coefficients}}{Component 1 coefficient draws (S x K).}
#'   \item{\code{estimates$m.coefficients}}{Component 2 coefficient draws (S x K).}
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
#' @keywords internal
#' @export
#' @importFrom rstan stan_model sampling
#' @import label.switching
glmMixBayes <- function(X, y, family = "gaussian", priors = NULL,
                        control = list(iterations = 1e4,
                                       burnin.iterations = 1e3,
                                       seed = sample.int(.Machine$integer.max, 1),
                                       cores = getOption("mc.cores", 1L)),
                        ...) {

  dcontrols <- list(...)
  iterations <- if ("iterations" %in% names(dcontrols)) {
    dcontrols$iterations
  } else if ("iterations" %in% names(control)) {
    control$iterations
  } else {
    1e4
  }

  burnin.iterations <- if ("burnin.iterations" %in% names(dcontrols)) {
    dcontrols$burnin.iterations
  } else if ("burnin.iterations" %in% names(control)) {
    control$burnin.iterations
  } else {
    1e3
  }

  # normalize family input
  if (inherits(family, "family")) family <- family$family
  family <- tolower(trimws(as.character(family)[1]))

  seed <- if ("seed" %in% names(dcontrols)) {
    dcontrols$seed
  } else if ("seed" %in% names(control)) {
    control$seed
  } else {
    sample.int(.Machine$integer.max, 1)
  }

  cores <- if ("cores" %in% names(dcontrols)) {
    dcontrols$cores
  } else if ("cores" %in% names(control)) {
    control$cores
  } else {
    getOption("mc.cores", 1L)
  }

  if (!(family %in% c("gaussian", "poisson", "binomial", "gamma"))) {
   stop("Error: the family should be gaussian, poisson, binomial, or gamma", call. = FALSE)
  }

  # Minimal input checks (formula/data checks handled upstream)
  if (!is.matrix(X) || !is.numeric(X)) stop("`X` must be a numeric matrix.", call. = FALSE)
  if (!is.numeric(y) || length(y) != nrow(X)) stop("`y` must be numeric with length nrow(X).", call. = FALSE)
  if (anyNA(X) || anyNA(y)) stop("NA values found in X or y", call. = FALSE)

  components <- rep(family, 2)

  # build priors (default + user overrides)
  priors <- fill_defaults(priors, p_family = family, model_type = "glm")

  # stan data (coerce for count/binary families if needed)
  stan_y <- y
  if (family %in% c("poisson", "binomial")) {
    stan_y <- as.integer(y)
    if (anyNA(stan_y)) stop("`y` could not be safely coerced to integer for family = '", family, "'.", call. = FALSE)
  }

  stan_data <- list(N = nrow(X), K = ncol(X), X = X, y = stan_y)

  # compile & sample
  model_code <- generate_stan(components, priors = priors)
  sm <- rstan::stan_model(model_code = model_code)
  fit <- rstan::sampling(
    sm,
    data   = stan_data,
    iter   = iterations,
    warmup = burnin.iterations,
    chains = 1,
    seed   = seed,
    cores  = cores
  )

  posterior <- rstan::extract(fit)
  z_samples <- posterior$z

  ##### Label switching adjustment #####

  # --- 0) Optional global pre-alignment by majority label -----------------------
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
    z      = z_samples,
    K      = 2
  )
  perm <- ls_out$permutations[["ECR-ITERATIVE-1"]]

  # Map z by the inverse permutation for each draw
  map_z <- function(z, perm) {
    if (!is.matrix(z)) stop("`z` must be an S x N matrix.", call. = FALSE)
    S <- nrow(z); K <- ncol(perm)
    if (K != 2L || nrow(perm) != S) stop("`perm` must be an S x 2 matrix of permutations.", call. = FALSE)
    for (i in seq_len(S)) {
      inv <- integer(K)
      inv[perm[i, ]] <- seq_len(K)
      z[i, ] <- inv[z[i, ]]
    }
    z
  }
  z_samples <- map_z(z_samples, perm)

  # Helper to permute paired component-specific parameters
  perm_pair <- function(a1, a2, perm) {
   if (!is.numeric(a1) || !is.numeric(a2)) stop("Inputs must be numeric.", call. = FALSE)
   if (!is.matrix(perm) || ncol(perm) != 2L) stop("`perm` must be an S x 2 matrix.", call. = FALSE)

   # Coerce component params to matrices S x p
   to_Sp <- function(a) {
    d <- dim(a)
    if (is.null(d)) {
     # vector length S
     S <- length(a)
     A <- matrix(as.numeric(a), ncol = 1L)
     return(list(S = S, A = A))
    }
    if (length(d) == 1L) {
     # 1D array, treat as vector
     S <- d[1]
     A <- matrix(as.numeric(a), ncol = 1L)
     return(list(S = S, A = A))
    }
    if (length(d) == 2L) {
     # matrix S x p
     S <- d[1]
     A <- matrix(as.numeric(a), nrow = d[1], ncol = d[2])
     return(list(S = S, A = A))
    }
    stop("Unsupported parameter shape; expected vector or matrix.", call. = FALSE)
   }

   x1 <- to_Sp(a1)
   x2 <- to_Sp(a2)

   if (x1$S != x2$S) stop("Different number of draws between components.", call. = FALSE)
   if (!identical(dim(x1$A), dim(x2$A))) stop("Shapes differ between component parameters.", call. = FALSE)

   S <- x1$S
   p <- ncol(x1$A)

   if (nrow(perm) != S) stop("`perm` must have one permutation per draw (nrow(perm) == S).", call. = FALSE)

   arr <- array(NA_real_, dim = c(S, 2L, p))
   arr[, 1L, ] <- x1$A
   arr[, 2L, ] <- x2$A

   arrp <- label.switching::permute.mcmc(arr, permutations = perm)[[1]]

   out1 <- arrp[, 1L, , drop = FALSE]
   out2 <- arrp[, 2L, , drop = FALSE]

   if (p == 1L) {
    return(list(`1` = as.numeric(out1[, 1L, 1L]), `2` = as.numeric(out2[, 1L, 1L])))
   }

   list(`1` = matrix(out1[, 1L, ], nrow = S, ncol = p),
        `2` = matrix(out2[, 1L, ], nrow = S, ncol = p))
  }

  if (is.null(posterior$z) || is.null(posterior$beta1) || is.null(posterior$beta2)) {
    stop("Missing expected parameters in posterior", call. = FALSE)
  }

  beta1.p <- posterior$beta1
  beta2.p <- posterior$beta2
  tmp <- perm_pair(beta1.p, beta2.p, perm)
  beta1.p <- tmp[[1]]
  beta2.p <- tmp[[2]]

  if (family == "gamma") {
    tmp <- perm_pair(posterior$phi1, posterior$phi2, perm)
    phi1.p <- tmp[[1]]
    phi2.p <- tmp[[2]]
  }

  if (family == "gaussian") {
    tmp <- perm_pair(posterior$sigma1, posterior$sigma2, perm)
    sigma1.p <- tmp[[1]]
    sigma2.p <- tmp[[2]]
  }

  # Set coefficient names from X if available
  cn <- colnames(X)
  if (!is.null(cn) && is.matrix(beta1.p) && ncol(beta1.p) == length(cn)) {
    colnames(beta1.p) <- cn
    colnames(beta2.p) <- cn
  }

  out <- list(
    m_samples = z_samples,
    estimates = list(
      coefficients   = beta1.p,
      m.coefficients = beta2.p
    ),
    family = family,
    call = match.call()
  )

  if (family == "gamma") {
    out$estimates$dispersion   <- 1 / phi1.p
    out$estimates$m.dispersion <- 1 / phi2.p
  }
  if (family == "gaussian") {
    out$estimates$dispersion   <- sigma1.p^2
    out$estimates$m.dispersion <- sigma2.p^2
  }

  class(out) <- "glmMixBayes"
  out
}

# ------------------------------------------------------------------------------
# S3 dispatch method for the postlink internal generic fitglm()
# ------------------------------------------------------------------------------
#' @keywords internal
#' @noRd
fitglm.adjMixBayes <- function(x, y, family, adjustment, control, ...) {

  # -------------------------------------------------------------------------
  # 1. Validation and Data Retrieval (align with mixture_glm.R pattern)
  # -------------------------------------------------------------------------
  full_data <- adjustment$data_ref$data
  if (is.null(full_data)) {
    stop("The 'adjustment' object does not contain linked data. ",
         "Please recreate the object with 'linked.data' provided.", call. = FALSE)
  }

  # -------------------------------------------------------------------------
  # 2. Align Adjustment Data to Outcome Model (X, Y) via row names
  # -------------------------------------------------------------------------
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

  # -------------------------------------------------------------------------
  # 3. Ensure no missingness in x/y (upstream should handle)
  # -------------------------------------------------------------------------
  if (anyNA(x) || anyNA(y)) {
    stop("NA values found in x or y. Upstream wrapper should remove missingness.", call. = FALSE)
  }

  # -------------------------------------------------------------------------
  # 4. Extract priors from dots/control (so users can supply priors upstream)
  # -------------------------------------------------------------------------
  dots <- list(...)
  priors <- NULL
  if ("priors" %in% names(dots)) {
    priors <- dots$priors
    dots$priors <- NULL
  } else if (!is.null(control) && is.list(control) && "priors" %in% names(control)) {
    priors <- control$priors
  }

  # -------------------------------------------------------------------------
  # 5. Dispatch to computational engine
  # -------------------------------------------------------------------------
  fit <- do.call(
    glmMixBayes,
    c(
      list(X = x, y = y, family = family, priors = priors, control = control),
      dots
    )
  )

  # -------------------------------------------------------------------------
  # 6. Post-processing (match mixture_glm.R conventions)
  # -------------------------------------------------------------------------
  fit$adjustment <- adjustment
  fit$call <- match.call()

  # Store observation names for downstream predict/plot helpers (if needed)
  if (!is.null(subset_names)) {
    fit$obs_names <- subset_names
  }

  fit
}
