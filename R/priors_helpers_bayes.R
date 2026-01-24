#' Fill default priors for Bayesian mixture models
#'
#' @description
#' Fill in weakly-informative default priors for a supported likelihood family and
#' then override with any user-supplied prior strings. Prior strings are inserted
#' directly into dynamically generated Stan code.
#'
#' @param priors a named list of user-supplied prior strings (or hyperparameter
#' definitions), e.g. \code{list(beta1 = "normal(0,5)", theta = "beta(2,2)")}.
#' If \code{NULL}, all defaults are used.
#' @param p_family the likelihood family. For GLM mixtures, one of
#' \code{"gaussian"}, \code{"poisson"}, \code{"binomial"}, or \code{"gamma"}.
#' For survival mixtures, one of \code{"gamma"} or \code{"weibull"}.
#' @param model_type the model type: \code{"glm"} (generalized linear model; default)
#' or \code{"survival"} (parametric survival model).
#'
#' @returns a named list of prior strings/hyperparameters containing one entry for
#' every parameter required by the requested family and model type.
#'
#' @examples
#' \dontrun{
#' fill_defaults(
#'   priors = list(beta1 = "normal(0, 2)", theta = "beta(2, 2)"),
#'   p_family = "gaussian",
#'   model_type = "glm"
#' )
#' }
#'
#' @keywords internal

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
  if (model_type == 'glm') {
    defaults <- switch(
      p_family,
      "gaussian" = list(
        beta1 = "normal(0,5)",
        sigma1 = "cauchy(0,2.5)",
        beta2  = "normal(0,5)",
        sigma2 = "cauchy(0,2.5)",
        theta = "beta(1,1)"
      ),
      "poisson" = list(
        beta1 = "normal(0,5)",
        beta2 = "normal(0,5)",
        theta = "beta(1,1)"
      ),
      "binomial" = list(
        beta1 = "normal(0,2.5)",
        beta2 = "normal(0,5)",
        theta = "beta(1,1)"
      ),
      "gamma" = list(
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
        # priors for regression coefficients and mix proportion
        beta1 = "normal(0,5)",
        beta2 = "normal(0,5)",
        theta = "beta(1,1)",
        # priors for shape parameters
        phi1 = "exponential(1)",
        phi2 = "exponential(1)"
      ),
      "weibull" = list(
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

#' Mark a string as a Stan function definition for prior injection
#'
#' @description
#' Wrap a character string to indicate it should be injected verbatim into the
#' Stan \code{functions} block during dynamic code generation. This is intended
#' for advanced use cases where a prior distribution is implemented via a custom
#' Stan function. The string is lightly validated for basic structure.
#'
#' @param stan_function_str a single character string containing a Stan function
#' definition.
#'
#' @returns the input string with class \code{"stan_function_string"} so that it
#' can be detected and routed to the Stan \code{functions} block.
#'
#' @examples
#' \dontrun{
#' f <- stan_func("
#'   real my_prior_lpdf(real x) {
#'     return normal_lpdf(x | 0, 1);
#'   }
#' ")
#' is_func(f)
#' }
#'
#' @keywords internal

stan_func <- function(stan_function_str) {
  if (!is.character(stan_function_str)) {
    stop("`stan_function_str` passed to stan_func() must be a single string.", call. = FALSE)
  } else if (!is_valid_stan_func(stan_function_str)) {
    stop("stan_func() called on invalid stan function")
  }
  structure(stan_function_str, class = c("stan_function_string", "character"))
}

#' Lightly validate a Stan function string
#'
#' @description
#' Helper used by \code{stan_func()} to perform a minimal structural check that a
#' string resembles a Stan function definition and does not contain R assignment
#' syntax (\code{<-}).
#'
#' @param str a character string.
#'
#' @returns \code{TRUE} if the string appears to be a Stan function definition and
#' \code{FALSE} otherwise.
#'
#' @keywords internal

is_valid_stan_func <- function(str) {
  txt <- gsub("[\r\n]+", " ", str)
  grepl("\\([^)]*\\).*\\{.*return.*\\}", txt, perl = TRUE) &&
    !grepl("<-", txt, fixed = TRUE)
}

#' Validate a user-supplied priors list for Stan code generation
#'
#' @description
#' Perform basic validation of the \code{priors} list used for dynamic Stan code
#' generation. Checks include: supported likelihood family; \code{priors} is a
#' named list; each prior string follows a simple \code{dist(arg1, arg2, ...)}
#' pattern; and any non-numeric hyperparameters referenced in a prior string are
#' also defined as keys in \code{priors} (so they can be declared/assigned in Stan).
#'
#' @param priors a named list of Stan-style prior strings and/or hyperparameter
#' definitions.
#' @param p_family the likelihood family.
#'
#' @returns invisibly \code{TRUE} if all checks pass (warnings may still be issued).
#'
#' @keywords internal

validate_args <- function(priors, p_family) {
  # p_family must be one of the supported families
  supported_families <- c("gaussian","poisson","binomial","gamma", "weibull")
  if (!(p_family %in% supported_families)) {
    stop("`p_family` must be one of: ", paste(supported_families, collapse = ", "))
  }

  # If priors == NULL, no need to check its validity beyond that
  if (is.null(priors)) {
    return(invisible(TRUE))
  }

  # priors list must be a named list or null (in which case all defaults are used)
  if (!is.list(priors) || is.null(names(priors))) {
    stop("`priors` must be a named list of prior strings. Yours is:\n",
         priors)
  }

  for (elt in priors) {
    if (is.character(elt) && !is_func(elt)) { # only looking at prior strings
      trimmed_elt <- trimws(elt) # remove whitespace

      # verify basic structure <dist>(<args>)
      pattern <- "^([^(\\s]+)\\s*\\(([^)]*)\\)\\s*$"
      match <- regexec(pattern, trimmed_elt, perl = TRUE)
      parts <- regmatches(trimmed_elt, match)[[1]]
      dist = parts[2]
      hyperparam = parts[3]

      if (length(parts) == 0
          || (hyperparam == "" && all(!vapply(priors, is_func, logical(1))))
          || dist == "") {
        stop("Prior defining string '",
             elt,
             "' must be of form '<dist>(<comma-separated args>)'")
      }

      # send warning about untested distributions
      if (!(dist %in% c("normal", "cauchy", "beta", "gamma", "exponential", "multi_normal"))) {
        warning("An untested distribution was defined as prior. Tested ones",
                " include but are not limited to normal, cauchy, beta, gamma,",
                " exponential, and multi_normal. Correct output can be expected",
                " regardless if the distribution is a valid stan distribution,",
                " or defined by a stan function injection with stan_func()---",
                " see documentation for details.")
      }

      # validate args
      if (nzchar(trimws(hyperparam))) { # if args has non-whitespace
        # get individual args
        raw_hyperparam_list <- strsplit(hyperparam, ",", fixed = TRUE)[[1]]
        hyperparam_list <- trimws(raw_hyperparam_list) # trim white space

        if (any(hyperparam_list == "") || grepl(",\\s*$", hyperparam)) {
          stop("Missing values in list of hyperparameters of prior string")
        }

        for (hyperparam in hyperparam_list) {
          # if hyperparameter is not a number and not a key defined in the priors
          if (is.na(suppressWarnings(as.numeric(hyperparam)))
              && !(hyperparam %in% names(priors))) {
            stop("A variable used in a prior string is not defined.",
                 " If a string defining a prior in priors list is normal(x,y),",
                 " x and y must be their own keys in list 'priors' to be used.")
          }
        }
      }
    }
  }
  invisible(TRUE)
  }

#' Detect Stan function injections in a priors list
#'
#' @description
#' Returns \code{TRUE} if an element was created by \code{stan_func()} and should
#' be injected into the Stan \code{functions} block.
#'
#' @param stan_function an object to test.
#'
#' @returns \code{TRUE} if \code{stan_function} inherits from
#' \code{"stan_function_string"}, \code{FALSE} otherwise.
#'
#' @keywords internal

is_func <- function(stan_function) {
  inherits(stan_function, "stan_function_string")
}

#' Build Stan code fragments for prior hyperparameters and functions
#'
#' @description
#' Convert a \code{priors} list into two character strings: (i) declarations and
#' assignments for any hyperparameters referenced by prior strings, and (ii) any
#' Stan functions supplied via \code{stan_func()}.
#'
#' @param priors a named list of priors and hyperparameter definitions.
#'
#' @returns a list with components:
#' \item{variable_defs}{character string containing Stan declarations/assignments
#' to be placed in \code{transformed data}.}
#' \item{function_defs}{character string containing Stan function definitions to
#' be placed in \code{functions}.}
#'
#' @keywords internal

get_stan_definitions <- function(priors) {
  # generate from list of priors the necessary variable definition Stan strings
  # to concatenate before the prior definition
  # ex. 'vector[2] mu;' and 'mu = [1, 2];' from list item mu = c(1,2)
  variable_declarations <- "" # ex. cov_matrix[3] beta1_sigma;
  variable_definitions <- "" # ex. beta1_sigma = ...; or vector[3] vec = [1,2,3]'
  function_definitions <- "" # for injecting stan into function blocks
  for (item_key in names(priors)) {
    # concatenate stan code for variables in dynamic stan generation
    # generated from processing key-value pair in priors list

    processed_vars <- process_variable(priors[[item_key]], item_key)
    variable_declarations <- paste0(
      variable_declarations, processed_vars[["declaration"]])
    variable_definitions <- paste0(
      variable_definitions, processed_vars[["definition"]])
    function_definitions <- paste0(
      function_definitions, processed_vars[["stan_func"]]
    )
    }
  # combine separate declarations with definitions to create a single variable
  # holding a string of complete definitions of prior hyperparameters
  # These were generated separately because declarations must come before definitions
  variable_definitions <- paste0(variable_declarations, variable_definitions)

  list(
    variable_defs = variable_definitions,
    function_defs = function_definitions
  )
}

#' Process a single priors-list entry into Stan code
#'
#' @description
#' Helper for \code{get_stan_definitions()}. Converts one key-value pair from the
#' \code{priors} list into Stan code fragments for (optional) declaration,
#' definition, and function injection.
#'
#' @param value the value associated with \code{key} in the priors list. Can be a
#' character prior string, a scalar/vector numeric hyperparameter, a symmetric
#' matrix hyperparameter, or a \code{stan_function_string}.
#' @param key the name of the priors-list entry (used as the Stan variable name).
#'
#' @returns a list with elements \code{declaration}, \code{definition}, and
#' \code{stan_func}.
#'
#' @keywords internal

process_variable <- function(value, key) {
  if (is.character(value)) { # is string
    if (is_func(value)) { # the string is a function i.e. stan_func("...stan function...")
      warning("Detected stan function definition. Note that stan function
                injection does not fully validate stan syntax.")
      return(list(declaration="", definition="", stan_func=as.character(value)))
    } else { # the string is the prior definition itself ex. normal(x,y)
      return(list(declaration="", definition="", stan_func=""))
  }
  } else if (is.matrix(value)) {
    if (nrow(value) == ncol(value) && isSymmetric(value)) {
      dim = nrow(value)
    } else {
      stop("Non-square or unsymmetric matrix found in list 'priors'")
    }
    dec = paste0("cov_matrix[", dim, "] ", key, ";") # declaration
    component_defs = ""
    for (i in 1:dim) {
      for (j in 1:dim) {
        component <- paste0(key, "[", i, ", ", j, "] = ", value[i, j], ";" )
        component_defs <- paste0(component_defs, component)
      }
    }
    return(list(declaration=dec, definition=component_defs, stan_func=""))
  } else if (is.numeric(value)) { # non-matrix numeric
    if (length(value) == 1L) { # not a vector; scalar
      def = paste0("real ", key, " = ", value, ";")
      return(list(declaration="", definition=def, stan_func=""))
    } else { # value is a vector
      len = length(value)
      # convert vector to stan list as a string:
      stan_vector = paste0("[", paste(value, collapse=", "), "]'")
      def = paste0("vector[", len, "] ", key, " = ", stan_vector, ";")
      return(list(declaration="", definition=def, stan_func=""))
    }
  } else {
    stop("Element of invalid type or form found in list 'priors'")
  }
}
