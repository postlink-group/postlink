# =============================================================================
# ai_helpers.R
# explain() — plain-language interpretation of fitted postlink models
#
# Provider and API key are passed directly as function parameters.
# If not specified, provider is auto-detected from available environment
# variables. No helper functions required.
#
# Supported providers:
#   "anthropic" — ANTHROPIC_API_KEY  (default if key present)
#   "openai"    — OPENAI_API_KEY
#   "gemini"    — GEMINI_API_KEY
#   "ollama"    — local, no key needed
# =============================================================================


# ── Internal: provider detection ──────────────────────────────────────────────

#' Detect provider from env vars or explicit argument
#' @noRd
.detect_provider <- function(provider = NULL) {
  if (!is.null(provider)) return(provider)
  if (nchar(Sys.getenv("ANTHROPIC_API_KEY")) > 0) return("anthropic")
  if (nchar(Sys.getenv("OPENAI_API_KEY"))    > 0) return("openai")
  if (nchar(Sys.getenv("GEMINI_API_KEY"))    > 0) return("gemini")
  "anthropic"  # fallback — will error with a clear message if key missing
}


#' Resolve the API key: explicit argument > environment variable
#' @noRd
.resolve_key <- function(api_key, provider) {
  if (!is.null(api_key) && nchar(api_key) > 0) return(api_key)

  key_name <- switch(provider,
    anthropic = "ANTHROPIC_API_KEY",
    openai    = "OPENAI_API_KEY",
    gemini    = "GEMINI_API_KEY",
    ollama    = NULL,
    stop("Unknown provider '", provider, "'. ",
         "Choose one of: anthropic, openai, gemini, ollama.", call. = FALSE)
  )

  if (is.null(key_name)) return(NULL)   # ollama needs no key

  key <- Sys.getenv(key_name)
  if (nchar(key) == 0) {
    stop(
      "No API key found for provider '", provider, "'.\n",
      "Set ", key_name, " in your environment, for example:\n",
      "  Sys.setenv(", key_name, " = 'your-key-here')\n",
      "Or pass the key directly: explain(fit, api_key = 'your-key-here')",
      call. = FALSE
    )
  }
  key
}


#' Check httr2 is available
#' @noRd
.check_httr2 <- function() {
  if (!requireNamespace("httr2", quietly = TRUE))
    stop("Package 'httr2' is required for AI features.\n",
         "Install it with: install.packages('httr2')", call. = FALSE)
}


#' Resolve context: inline text or file path
#' @noRd
.resolve_context <- function(context) {
  if (is.null(context)) return(NULL)
  # If it looks like a file path and the file exists, read it
  if (nchar(context) < 500 && file.exists(context)) {
    return(paste(readLines(context, warn = FALSE), collapse = "\n"))
  }
  context
}


# ── Internal: LLM dispatch ────────────────────────────────────────────────────

#' Unified LLM call — dispatches to correct provider
#' Returns a list: list(text = "...", metadata = list(...))
#' @noRd
.call_llm <- function(user_msg, system_msg, provider, api_key,
                       model = NULL, max_tokens = 2048,
                       temperature = NULL, seed = NULL) {
  .check_httr2()
  switch(provider,
    anthropic = .call_anthropic(user_msg, system_msg, api_key, model,
                                max_tokens, temperature, seed),
    openai    = .call_openai(user_msg, system_msg, api_key, model,
                             max_tokens, temperature, seed),
    gemini    = .call_gemini(user_msg, system_msg, api_key, model,
                             max(max_tokens, 8192L), temperature),
    ollama    = .call_ollama(user_msg, system_msg, model, max_tokens, temperature),
    stop("Unknown provider '", provider, "'.", call. = FALSE)
  )
}


#' @noRd
.call_anthropic <- function(user_msg, system_msg, api_key,
                              model = NULL, max_tokens = 2048,
                              temperature = NULL, seed = NULL) {
  model <- model %||% "claude-sonnet-4-6"
  body  <- list(model = model, max_tokens = max_tokens, system = system_msg,
                messages = list(list(role = "user", content = user_msg)))
  if (!is.null(temperature)) body$temperature <- temperature
  # Note: Anthropic does not support a seed parameter; reproducibility is via temperature=0

  resp <- httr2::request("https://api.anthropic.com/v1/messages") |>
    httr2::req_headers("x-api-key" = api_key,
                       "anthropic-version" = "2023-06-01",
                       "content-type" = "application/json") |>
    httr2::req_body_json(body) |>
    httr2::req_timeout(60) |>
    httr2::req_error(is_error = function(r) FALSE) |>
    httr2::req_perform()

  if (httr2::resp_status(resp) != 200L) {
    err <- tryCatch(httr2::resp_body_json(resp)$error$message, error=function(e) "")
    stop("Anthropic API error: ", err, call. = FALSE)
  }
  parsed <- httr2::resp_body_json(resp)
  list(
    text = parsed$content[[1]]$text,
    metadata = list(
      provider         = "anthropic",
      model            = parsed$model,
      input_tokens     = parsed$usage$input_tokens,
      output_tokens    = parsed$usage$output_tokens,
      stop_reason      = parsed$stop_reason
    )
  )
}


#' @noRd
.call_openai <- function(user_msg, system_msg, api_key,
                          model = NULL, max_tokens = 2048,
                          temperature = NULL, seed = NULL) {
  model <- model %||% "gpt-4o"
  body  <- list(model = model, max_tokens = max_tokens,
                messages = list(list(role = "system",  content = system_msg),
                                list(role = "user",    content = user_msg)))
  if (!is.null(temperature)) body$temperature <- temperature
  if (!is.null(seed))        body$seed        <- as.integer(seed)

  resp <- httr2::request("https://api.openai.com/v1/chat/completions") |>
    httr2::req_headers("Authorization" = paste("Bearer", api_key)) |>
    httr2::req_body_json(body, auto_unbox = TRUE) |>
    httr2::req_timeout(60) |>
    httr2::req_error(is_error = function(r) FALSE) |>
    httr2::req_perform()

  if (httr2::resp_status(resp) != 200L) {
    err <- tryCatch(httr2::resp_body_json(resp)$error$message, error=function(e) "")
    stop("OpenAI API error: ", err, call. = FALSE)
  }
  parsed <- httr2::resp_body_json(resp)
  list(
    text = parsed$choices[[1]]$message$content,
    metadata = list(
      provider           = "openai",
      model              = parsed$model,
      input_tokens       = parsed$usage$prompt_tokens,
      output_tokens      = parsed$usage$completion_tokens,
      system_fingerprint = parsed$system_fingerprint,
      finish_reason      = parsed$choices[[1]]$finish_reason
    )
  )
}


#' @noRd
.call_gemini <- function(user_msg, system_msg, api_key,
                          model = NULL, max_tokens = 8192,
                          temperature = NULL) {
  model <- model %||% "gemini-3.8-flash"

  # ── Thinking-model guard ──────────────────────────────────────────────────
  # Gemini models whose names contain "thinking" output internal chain-of-thought
  # rather than a clean response, producing garbled or cut-off text.
  is_thinking <- grepl("thinking", model, ignore.case = TRUE)

  if (is_thinking) {
    fallback <- "gemini-3.8-flash"
    warning(
      "The model '", model, "' is a Gemini 'thinking' model that outputs ",
      "internal reasoning steps rather than a direct response. This produces ",
      "garbled or incomplete output for postlink AI features.\n",
      "Switching to '", fallback, "' instead.\n",
      "To suppress this warning, set model = '", fallback, "' explicitly.",
      call. = FALSE
    )
    model <- fallback
  }

  combined      <- paste0(system_msg, "\n\n", user_msg)
  gen_config    <- list(maxOutputTokens = max_tokens)
  if (!is.null(temperature)) gen_config$temperature <- temperature

  body <- list(
    contents        = list(list(parts = list(list(text = combined)))),
    generationConfig = gen_config
  )

  url  <- paste0("https://generativelanguage.googleapis.com/v1beta/models/",
                 model, ":generateContent")
  resp <- httr2::request(url) |>
    httr2::req_headers("content-type"    = "application/json",
                       "x-goog-api-key"  = api_key) |>
    httr2::req_body_json(body) |>
    httr2::req_timeout(60) |>
    httr2::req_error(is_error = function(r) FALSE) |>
    httr2::req_perform()

  if (httr2::resp_status(resp) != 200L) {
    err <- tryCatch(httr2::resp_body_json(resp)$error$message, error=function(e) "")
    if (grepl("429|RESOURCE_EXHAUSTED|quota|rate", err, ignore.case = TRUE)) {
      stop(
        "Gemini API rate limit reached (free tier).\n",
        "Options:\n",
        "  1. Wait a few minutes and try again.\n",
        "  2. Upgrade to a paid Gemini plan for higher rate limits.\n",
        "  3. Switch provider: explain(fit, provider = 'anthropic')",
        call. = FALSE
      )
    }
    stop("Gemini API error: ", err, call. = FALSE)
  }

  parsed     <- httr2::resp_body_json(resp)
  candidates <- parsed$candidates
  if (is.null(candidates) || length(candidates) == 0)
    stop("Gemini returned no candidates. Model may be unavailable or request blocked.",
         call. = FALSE)

  parts <- candidates[[1]]$content$parts
  if (is.null(parts) || length(parts) == 0)
    stop("Gemini returned an empty response. Check your API quota or try a different model.",
         call. = FALSE)

  texts <- Filter(Negate(is.null), lapply(parts, `[[`, "text"))
  text  <- paste(texts, collapse = "")

  usage <- parsed$usageMetadata
  list(
    text = text,
    metadata = list(
      provider      = "gemini",
      model         = model,
      input_tokens  = usage$promptTokenCount,
      output_tokens = usage$candidatesTokenCount,
      finish_reason = candidates[[1]]$finishReason
    )
  )
}


#' @noRd
.call_ollama <- function(user_msg, system_msg, model = NULL, max_tokens = 2048,
                         temperature = NULL) {
  host  <- Sys.getenv("OLLAMA_HOST", "http://localhost:11434")
  model <- model %||% "llama3"
  opts  <- list(num_predict = max_tokens)
  if (!is.null(temperature)) opts$temperature <- temperature
  body  <- list(model = model, stream = FALSE,
                prompt = paste0(system_msg, "\n\n", user_msg),
                options = opts)

  resp <- tryCatch(
    httr2::request(paste0(host, "/api/generate")) |>
      httr2::req_headers("content-type" = "application/json") |>
      httr2::req_body_json(body) |>
      httr2::req_timeout(120) |>
      httr2::req_error(is_error = function(r) FALSE) |>
      httr2::req_perform(),
    error = function(e)
      stop("Cannot connect to Ollama at ", host, ".\n",
           "Make sure Ollama is running: ollama serve\n",
           "Set OLLAMA_HOST if using a non-default address.",
           call. = FALSE)
  )
  if (httr2::resp_status(resp) != 200L) {
    err <- tryCatch(httr2::resp_body_json(resp)$error, error=function(e) "")
    stop("Ollama error: ", err, call. = FALSE)
  }
  parsed <- httr2::resp_body_json(resp)
  list(
    text = parsed$response,
    metadata = list(
      provider      = "ollama",
      model         = model,
      input_tokens  = NULL,
      output_tokens = NULL
    )
  )
}


# ── Internal: null-coalescing ─────────────────────────────────────────────────
`%||%` <- function(a, b) if (!is.null(a)) a else b


# ── Internal: method description ─────────────────────────────────────────────

#' @noRd
.describe_method <- function(fit) {

  if (inherits(fit, c("glmMixture","coxphMixture","ctableMixture"))) {
    mp         <- fit$match.prob
    match_info <- if (!is.null(mp))
      sprintf("mean posterior correct match probability = %.3f (implied mismatch rate ~%.1f%%)",
              mean(mp), 100*(1-mean(mp))) else NULL
    mrate_used <- tryCatch(as.character(fit$adjustment$m.rate), error=function(e) NULL)
    safe_used  <- !is.null(tryCatch(fit$adjustment$safe.matches, error=function(e) NULL))
    mformula   <- tryCatch(deparse(fit$adjustment$m.formula), error=function(e) "~1")

    extras <- Filter(Negate(is.null), c(
      if (!is.null(mrate_used)) paste0("global mismatch rate constraint: ", mrate_used),
      if (!mformula %in% c("~1","~ 1")) paste0("mismatch model covariates: ", mformula),
      if (safe_used) "safe matches locked at match probability 1",
      match_info,
      if (!is.null(fit$converged))
        if (isTRUE(fit$converged)) "EM converged" else "WARNING: EM did not converge"
    ))

    return(list(
      method_short = "EM Mixture Model (adjMixture)",
      method_description = paste0(
        "EM-based mixture modeling (Slawski et al., 2025). Each record pair is ",
        "modelled as correct match (C=1) or false link (C=0). EM iterates between ",
        "estimating posterior match probabilities (E-step) and updating regression ",
        "parameters for true matches (M-step). Coefficients are for correctly matched ",
        "records only — already linkage-error-corrected."),
      estimator_note = NULL,
      extras = extras
    ))
  }

  if (inherits(fit, c("glmELE","coxphELE"))) {
    coefs <- fit$coefficients
    est   <- if (is.matrix(coefs)) rownames(coefs) else "ratio"
    n_est <- if (is.matrix(coefs)) nrow(coefs) else 1L
    return(list(
      method_short = "ELE Weighting (adjELE)",
      method_description = paste0(
        "Bias-adjusted estimating equations under the Exchangeable Linkage Error ",
        "(ELE) model (Chambers, 2009). Within each block, mismatch errors are assumed ",
        "exchangeable. Estimating equations modified using block-level mismatch rates. ",
        "Already linkage-error-corrected estimates."),
      estimator_note = if (n_est > 1)
        paste0("ELE computed ", n_est, " estimator(s): ",
               paste(rownames(coefs), collapse=", "),
               ". Reporting first (", rownames(coefs)[1], ").") else NULL,
      n_estimators = n_est,
      estimator_names = if (is.matrix(coefs)) rownames(coefs) else NULL,
      extras = paste0("weight matrix estimator(s): ", paste(est, collapse=", "))
    ))
  }

  if (inherits(fit, c("glmMixBayes","survMixBayes"))) {
    return(list(
      method_short = "Bayesian Mixture Model (adjMixBayes)",
      method_description = paste0(
        "Bayesian mixture modeling with MCMC posterior inference (Gutman et al., 2016). ",
        "Posterior samples via Stan (No-U-Turn Sampler). Reported coefficients are ",
        "posterior estimates for correctly matched records — already corrected."),
      estimator_note = NULL,
      n_estimators = 1L,
      estimator_names = NULL,
      extras = character(0)
    ))
  }

  list(method_short = "No adjustment (Naive)",
       method_description = "No linkage error correction applied.",
       estimator_note = NULL,
       n_estimators = 1L,
       estimator_names = NULL,
       extras = character(0))
}


# ── Internal: coefficient extraction ─────────────────────────────────────────

#' Extract outcome model coefficients as a data frame with significance
#' Uses summary() to get accurate p-values rather than rule-of-thumb |z|>2
#' @noRd
.get_coef_df <- function(fit) {
  postlink_classes <- c("glmMixture","glmELE","glmMixBayes",
                        "coxphMixture","coxphELE","survMixBayes","ctableMixture")

  # Standard glm/lm — use summary() for exact p-values
  if (!any(sapply(postlink_classes, function(cl) inherits(fit, cl)))) {
    cf <- tryCatch(coef(summary(fit)), error = function(e) NULL)
    if (!is.null(cf) && ncol(cf) >= 4) {
      return(data.frame(
        term      = rownames(cf),
        estimate  = round(as.numeric(cf[, 1]), 4),
        se        = round(as.numeric(cf[, 2]), 4),
        statistic = round(as.numeric(cf[, 3]), 3),
        p_value   = round(as.numeric(cf[, 4]), 4),
        row.names = NULL, stringsAsFactors = FALSE
      ))
    }
  }

  # glmELE: use summary()$coefficients[[1]] (first estimator matrix)
  if (inherits(fit, "glmELE")) {
    cf <- tryCatch(summary(fit)$coefficients[[1]], error = function(e) NULL)
    if (!is.null(cf) && is.matrix(cf) && ncol(cf) >= 4) {
      return(data.frame(
        term      = rownames(cf),
        estimate  = round(as.numeric(cf[, 1]), 4),
        se        = round(as.numeric(cf[, 2]), 4),
        statistic = round(as.numeric(cf[, 3]), 3),
        p_value   = round(as.numeric(cf[, 4]), 4),
        row.names = NULL, stringsAsFactors = FALSE
      ))
    }
  }

  # glmMixBayes / survMixBayes: fit$coefficients and fit$var don't exist;
  # use posterior mean and SD from fit$estimates$coefficients
  if (inherits(fit, c("glmMixBayes", "survMixBayes"))) {
    post_coef <- tryCatch(fit$estimates$coefficients, error = function(e) NULL)
    if (!is.null(post_coef) && is.matrix(post_coef)) {
      est        <- colMeans(post_coef)
      se         <- apply(post_coef, 2, sd)
      term_names <- colnames(post_coef)
    } else {
      est        <- numeric(0)
      se         <- numeric(0)
      term_names <- character(0)
    }
    z_val <- as.numeric(est) / as.numeric(se)
    return(data.frame(
      term      = term_names,
      estimate  = round(as.numeric(est), 4),
      se        = round(as.numeric(se),  4),
      statistic = round(z_val, 3),
      p_value   = round(2 * pnorm(-abs(z_val)), 4),
      row.names = NULL, stringsAsFactors = FALSE
    ))
  }

  # Other postlink objects (glmMixture, coxphMixture, coxphELE, ctableMixture):
  # use summary() where available, fall back to fit$coefficients + fit$var
  cf <- tryCatch(coef(summary(fit)), error = function(e) NULL)
  if (!is.null(cf) && is.matrix(cf) && ncol(cf) >= 4) {
    return(data.frame(
      term      = rownames(cf),
      estimate  = round(as.numeric(cf[, 1]), 4),
      se        = round(as.numeric(cf[, 2]), 4),
      statistic = round(as.numeric(cf[, 3]), 3),
      p_value   = round(as.numeric(cf[, 4]), 4),
      row.names = NULL, stringsAsFactors = FALSE
    ))
  }

  # Last-resort fallback: manual extraction from fit$coefficients + fit$var
  if (is.matrix(fit$coefficients)) {
    est        <- fit$coefficients[1, ]
    term_names <- colnames(fit$coefficients)
  } else {
    est        <- fit$coefficients
    term_names <- names(fit$coefficients)
  }

  if (!is.null(fit$var)) {
    var_mat <- if (is.list(fit$var)) fit$var[[1]] else fit$var
    se      <- sqrt(diag(var_mat)[seq_len(length(est))])
  } else {
    se <- rep(NA_real_, length(est))
  }

  z_val <- as.numeric(est) / as.numeric(se)
  data.frame(
    term      = term_names,
    estimate  = round(as.numeric(est), 4),
    se        = round(as.numeric(se),  4),
    statistic = round(z_val, 3),
    p_value   = round(2 * pnorm(-abs(z_val)), 4),
    row.names = NULL, stringsAsFactors = FALSE
  )
}


#' Format coefficient data frame as text for LLM prompt
#' @noRd
.format_coef_lines <- function(coef_df) {
  has_p <- "p_value" %in% names(coef_df)
  paste(apply(coef_df, 1, function(r) {
    sig <- if (has_p) {
      p <- as.numeric(r["p_value"])
      if (!is.na(p)) {
        if (p < 0.001) " ***" else if (p < 0.01) " **" else if (p < 0.05) " *" else " (ns)"
      } else ""
    } else ""
    sprintf("  %s: estimate=%s, SE=%s, z=%s, p=%s%s",
            r["term"], r["estimate"], r["se"],
            r["statistic"],
            if (has_p) r["p_value"] else "NA",
            sig)
  }), collapse = "\n")
}


#' Extract mismatch model coefficients (m.coefficients) as text
#' @noRd
.get_mismatch_coef_text <- function(fit) {
  mc <- fit$m.coefficients
  if (is.null(mc)) return(NULL)
  lines <- paste(names(mc), round(as.numeric(mc), 4), sep = " = ", collapse = ", ")
  paste0("Mismatch indicator model coefficients (log-odds of false link): ", lines)
}


# ── Public: explain() ─────────────────────────────────────────────────────────

#' Explain One or More Fitted postlink Models in Plain Language
#'
#' Sends coefficient estimates and standard errors from one or more fitted
#' models to an LLM and returns a plain-language interpretation. When a named
#' list of models is provided, the LLM compares results across all models.
#'
#' @param fit A single fitted postlink object (e.g. from \code{\link{plglm}},
#'   \code{\link{plcoxph}}), or a named list of fitted objects (postlink or
#'   standard \code{glm}/\code{lm}). When a list is provided, comparison mode
#'   is triggered automatically.
#' @param context Optional character string or file path. Study context passed
#'   to the LLM to improve interpretation. Can be a short description or a
#'   path to a \code{.txt} file containing longer context (e.g. a variable
#'   dictionary or study background).
#' @param provider Character string. LLM provider to use. One of
#'   \code{"anthropic"} (default if \code{ANTHROPIC_API_KEY} is set),
#'   \code{"openai"}, \code{"gemini"}, or \code{"ollama"} (local, no key
#'   needed). If \code{NULL}, auto-detected from available environment
#'   variables.
#' @param api_key Character string. API key for the chosen provider. If
#'   \code{NULL} (default), the key is read from the appropriate environment
#'   variable (\code{ANTHROPIC_API_KEY}, \code{OPENAI_API_KEY}, or
#'   \code{GEMINI_API_KEY}). Set with
#'   \code{Sys.setenv(ANTHROPIC_API_KEY = "sk-ant-...")} or add to
#'   \code{~/.Renviron}.
#' @param model Character string. Model name for the chosen provider.
#'   Defaults: Anthropic \eqn{\to} \code{"claude-sonnet-4-6"}, OpenAI
#'   \eqn{\to} \code{"gpt-4o"}, Gemini \eqn{\to} \code{"gemini-3.8-flash"},
#'   Ollama \eqn{\to} \code{"llama3"}.
#' @param max_tokens Integer. Maximum number of tokens in the LLM response.
#'   Default \code{2048}. Increase for longer or more complex outputs.
#' @param temperature Numeric in \code{[0, 1]}. Controls randomness of the LLM
#'   output. \code{0} gives the most deterministic (reproducible) output;
#'   higher values increase creativity/variability. Default is the provider's
#'   own default (usually around 1.0). Set \code{temperature = 0} for
#'   reproducible results.
#' @param seed Integer. Random seed for reproducibility (supported by OpenAI;
#'   ignored by other providers). When set alongside \code{temperature = 0},
#'   OpenAI outputs are highly reproducible.
#' @param verbose Logical. Print to console. Default \code{TRUE}.
#'
#' @return Invisibly returns the explanation as a character string with a
#'   \code{"postlink_metadata"} attribute containing provider, model, and
#'   token usage information.
#'
#' @details
#' For postlink objects, the LLM is told the coefficients are already
#' linkage-error-corrected and will not warn about bias that has been addressed.
#'
#' For Mixture and MixBayes models, the mismatch indicator model
#' (\code{m.coefficients}) is also summarised and sent to the LLM alongside
#' the outcome model coefficients.
#'
#' For ELE models with multiple weighting methods, only the first estimator is
#' interpreted; a warning is issued if additional estimators are present.
#'
#' Requires the \code{httr2} package (\code{install.packages("httr2")}).
#'
#' @examples
#' \dontrun{
#' data(lifem)
#' adj <- adjMixture(linked.data = lifem, m.formula = ~ commf + comml,
#'                   m.rate = 0.05, safe.matches = hndlnk)
#' fit_adj <- plglm(age_at_death ~ poly(unit_yob, 3, raw = TRUE),
#'                  family = "gaussian", adjustment = adj)
#'
#' # Anthropic (key from environment variable)
#' Sys.setenv(ANTHROPIC_API_KEY = "sk-ant-...")
#' result <- explain(fit_adj)
#'
#' # Reproducible output
#' result <- explain(fit_adj, temperature = 0)
#'
#' # Check token usage
#' attr(result, "postlink_metadata")
#'
#' # Pass key directly
#' explain(fit_adj, api_key = "sk-ant-...")
#'
#' # OpenAI with reproducible seed
#' explain(fit_adj, provider = "openai", api_key = "sk-...",
#'         temperature = 0, seed = 42)
#'
#' # Local Ollama (no key needed)
#' explain(fit_adj, provider = "ollama", model = "llama3")
#'
#' # Context from a text file
#' explain(fit_adj, context = "path/to/study_background.txt")
#'
#' # Multi-model comparison
#' fit_naive <- glm(age_at_death ~ poly(unit_yob, 3, raw = TRUE),
#'                  data = lifem, family = "gaussian")
#' explain(list(naive = fit_naive, adjusted = fit_adj))
#' }
#'
#' @export
explain <- function(fit,
                    context     = NULL,
                    provider    = NULL,
                    api_key     = NULL,
                    model       = NULL,
                    max_tokens  = 2048,
                    temperature = NULL,
                    seed        = NULL,
                    verbose     = TRUE) {

  provider <- .detect_provider(provider)
  api_key  <- .resolve_key(api_key, provider)
  context  <- .resolve_context(context)

  is_fit_list <- is.list(fit) && !inherits(fit, c(
    "glm","lm","glmMixture","glmELE","glmMixBayes",
    "coxphMixture","coxphELE","survMixBayes","ctableMixture"
  ))

  if (is_fit_list) {
    if (is.null(names(fit)) || any(names(fit) == ""))
      stop("When 'fit' is a list, all elements must be named.\n",
           "Example: list(naive = fit_glm, adjusted = fit_plglm)", call. = FALSE)
    return(.explain_multiple(fit, context=context, provider=provider,
                              api_key=api_key, model=model, max_tokens=max_tokens,
                              temperature=temperature, seed=seed, verbose=verbose))
  }

  .explain_single(fit, context=context, provider=provider,
                  api_key=api_key, model=model, max_tokens=max_tokens,
                  temperature=temperature, seed=seed, verbose=verbose)
}


#' @noRd
.explain_single <- function(fit, context=NULL, provider, api_key, model,
                              max_tokens=2048, temperature=NULL, seed=NULL, verbose=TRUE) {

  coef_df     <- .get_coef_df(fit)
  method_info <- .describe_method(fit)
  family_name <- tryCatch(fit$family$family, error=function(e) "unknown")
  family_label <- switch(family_name,
    gaussian = "Gaussian (linear) regression",
    binomial = "logistic regression",
    poisson  = "Poisson regression",
    Gamma    = "Gamma regression",
    paste0(family_name, " regression"))

  # ELE: warn if multiple estimators present (only first is interpreted)
  if (inherits(fit, c("glmELE","coxphELE")) &&
      !is.null(method_info$n_estimators) && method_info$n_estimators > 1) {
    warning(
      "ELE fit contains ", method_info$n_estimators, " estimators (",
      paste(method_info$estimator_names, collapse=", "), "). ",
      "Only the first ('", method_info$estimator_names[1], "') is interpreted. ",
      "Run explain() separately for each estimator if needed.",
      call. = FALSE
    )
  }

  coef_lines   <- .format_coef_lines(coef_df)
  mismatch_txt <- .get_mismatch_coef_text(fit)

  extras_block <- if (length(method_info$extras) > 0)
    paste0("Adjustment details:\n",
           paste(paste0("  - ", method_info$extras), collapse="\n"), "\n\n") else ""
  estimator_note <- if (!is.null(method_info$estimator_note))
    paste0(method_info$estimator_note, "\n\n") else ""
  context_block  <- if (!is.null(context))
    paste0("Study context:\n", context, "\n\n") else ""
  mismatch_block <- if (!is.null(mismatch_txt))
    paste0("\n\nMismatch indicator model:\n", mismatch_txt,
           "\n(These predict whether a record pair is a false link, ",
           "not the substantive outcome.)") else ""

  system_msg <- paste0(
    "You are a statistician explaining post-linkage data analysis results. ",
    "CRITICAL: The outcome model coefficients are ALREADY ADJUSTED for record ",
    "linkage errors — do NOT warn about mismatch bias, it has been corrected. ",
    "Tasks: (1) Briefly name the adjustment method and what it corrected for. ",
    "(2) Interpret each corrected coefficient in plain language, noting direction ",
    "and size. (3) Comment on statistical significance using the p-values provided ",
    "(* p<0.05, ** p<0.01, *** p<0.001). ",
    "(4) If mismatch model coefficients are provided, briefly note what they tell us ",
    "about which record pairs were more likely to be false links. ",
    "(5) No causation unless context supports it. 4-6 sentences. No headers.")

  user_msg <- paste0(
    context_block,
    "Model: ", family_label, "\n\n",
    "Adjustment: ", method_info$method_short, "\n",
    "What it does: ", method_info$method_description, "\n\n",
    extras_block, estimator_note,
    "Outcome model coefficients (already linkage-error-corrected):\n",
    coef_lines,
    mismatch_block, "\n\n",
    "Interpret these results.")

  message("Querying ", provider, " API...")
  result <- .call_llm(user_msg, system_msg, provider, api_key, model,
                       max_tokens = max_tokens, temperature = temperature, seed = seed)

  if (verbose) {
    cat("\n--- postlink: explain() ---\n\n")
    cat(result$text)
    cat("\n\n[AI-generated. Verify before use.]\n")
  }

  out <- invisible(result$text)
  attr(out, "postlink_metadata") <- result$metadata
  out
}


#' @noRd
.explain_multiple <- function(fit_list, context=NULL, provider, api_key, model,
                               max_tokens=2048, temperature=NULL, seed=NULL, verbose=TRUE) {

  context_block <- if (!is.null(context))
    paste0("Study context:\n", context, "\n\n") else ""

  model_blocks <- lapply(names(fit_list), function(nm) {
    m       <- fit_list[[nm]]
    coef_df <- .get_coef_df(m)
    method  <- .describe_method(m)
    n_obs   <- tryCatch(nrow(m$model),
                        error=function(e) tryCatch(length(m$residuals),
                                                   error=function(e) NA))
    coef_lines <- .format_coef_lines(coef_df)
    paste0("Model: ", nm, "\nMethod: ", method$method_short,
           if (!is.na(n_obs)) paste0("\nn = ", format(n_obs, big.mark=",")) else "",
           "\nCoefficients:\n", coef_lines)
  })

  system_msg <- paste0(
    "You are a statistician comparing regression models from a post-linkage analysis. ",
    "CRITICAL: Coefficients from postlink-adjusted models are ALREADY CORRECTED for ",
    "linkage errors — do NOT warn about bias in adjusted models. ",
    "(1) Briefly describe what each model does. ",
    "(2) Compare coefficients — note key differences in magnitude and direction. ",
    "(3) Explain what the differences reveal about linkage error correction. ",
    "(4) Comment on statistical significance where p-values are provided. ",
    "(5) No causation unless context supports it. 4-6 sentences. No headers.")

  user_msg <- paste0(
    context_block,
    length(fit_list), " models were fitted:\n\n",
    paste(model_blocks, collapse="\n\n---\n\n"), "\n\n",
    "Compare these models and explain what the differences reveal.")

  message("Querying ", provider, " API (comparing ", length(fit_list), " models)...")
  result <- .call_llm(user_msg, system_msg, provider, api_key, model,
                       max_tokens = max_tokens, temperature = temperature, seed = seed)

  if (verbose) {
    cat("\n--- postlink: explain() — Model Comparison ---\n\n")
    cat(result$text)
    cat("\n\n[AI-generated. Verify before use.]\n")
  }

  out <- invisible(result$text)
  attr(out, "postlink_metadata") <- result$metadata
  out
}
