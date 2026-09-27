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


# ── Internal: LLM dispatch ────────────────────────────────────────────────────

#' Unified LLM call — dispatches to correct provider
#' @noRd
.call_llm <- function(user_msg, system_msg, provider, api_key,
                       model = NULL, max_tokens = 600) {
  .check_httr2()
  switch(provider,
    anthropic = .call_anthropic(user_msg, system_msg, api_key, model, max_tokens),
    openai    = .call_openai(user_msg, system_msg, api_key, model, max_tokens),
    gemini    = .call_gemini(user_msg, system_msg, api_key, model, max_tokens),
    ollama    = .call_ollama(user_msg, system_msg, model, max_tokens),
    stop("Unknown provider '", provider, "'.", call. = FALSE)
  )
}


#' @noRd
.call_anthropic <- function(user_msg, system_msg, api_key,
                              model = NULL, max_tokens = 600) {
  model <- model %||% "claude-sonnet-4-6"
  body  <- list(model = model, max_tokens = max_tokens, system = system_msg,
                messages = list(list(role = "user", content = user_msg)))

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
  httr2::resp_body_json(resp)$content[[1]]$text
}


#' @noRd
.call_openai <- function(user_msg, system_msg, api_key,
                          model = NULL, max_tokens = 600) {
  model <- model %||% "gpt-4o"
  body  <- list(model = model, max_tokens = max_tokens,
                messages = list(list(role = "system",  content = system_msg),
                                list(role = "user",    content = user_msg)))

  resp <- httr2::request("https://api.openai.com/v1/chat/completions") |>
    httr2::req_headers("Authorization" = paste("Bearer", api_key),
                       "content-type"  = "application/json") |>
    httr2::req_body_json(body) |>
    httr2::req_timeout(60) |>
    httr2::req_error(is_error = function(r) FALSE) |>
    httr2::req_perform()

  if (httr2::resp_status(resp) != 200L) {
    err <- tryCatch(httr2::resp_body_json(resp)$error$message, error=function(e) "")
    stop("OpenAI API error: ", err, call. = FALSE)
  }
  httr2::resp_body_json(resp)$choices[[1]]$message$content
}


#' @noRd
.call_gemini <- function(user_msg, system_msg, api_key,
                          model = NULL, max_tokens = 600) {
  model <- model %||% "gemini-1.5-flash"
  combined <- paste0(system_msg, "\n\n", user_msg)
  body  <- list(contents = list(list(parts = list(list(text = combined)))),
                generationConfig = list(maxOutputTokens = max_tokens))

  url  <- paste0("https://generativelanguage.googleapis.com/v1beta/models/",
                 model, ":generateContent?key=", api_key)
  resp <- httr2::request(url) |>
    httr2::req_headers("content-type" = "application/json") |>
    httr2::req_body_json(body) |>
    httr2::req_timeout(60) |>
    httr2::req_error(is_error = function(r) FALSE) |>
    httr2::req_perform()

  if (httr2::resp_status(resp) != 200L) {
    err <- tryCatch(httr2::resp_body_json(resp)$error$message, error=function(e) "")
    stop("Gemini API error: ", err, call. = FALSE)
  }
  httr2::resp_body_json(resp)$candidates[[1]]$content$parts[[1]]$text
}


#' @noRd
.call_ollama <- function(user_msg, system_msg, model = NULL, max_tokens = 600) {
  host  <- Sys.getenv("OLLAMA_HOST", "http://localhost:11434")
  model <- model %||% "llama3"
  body  <- list(model = model, stream = FALSE,
                prompt = paste0(system_msg, "\n\n", user_msg),
                options = list(num_predict = max_tokens))

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
  httr2::resp_body_json(resp)$response
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
    mrate_used <- tryCatch(as.character(fit$call$adjustment$m.rate), error=function(e) NULL)
    safe_used  <- !is.null(tryCatch(fit$call$adjustment$safe.matches, error=function(e) NULL))
    mformula   <- tryCatch(deparse(fit$call$adjustment$m.formula), error=function(e) "~1")

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
    return(list(
      method_short = "ELE Weighting (adjELE)",
      method_description = paste0(
        "Bias-adjusted estimating equations under the Exchangeable Linkage Error ",
        "(ELE) model (Chambers, 2009). Within each block, mismatch errors are assumed ",
        "exchangeable. Estimating equations modified using block-level mismatch rates. ",
        "Already linkage-error-corrected estimates."),
      estimator_note = if (is.matrix(coefs) && nrow(coefs) > 1)
        paste0("ELE computed ", nrow(coefs), " estimator(s): ",
               paste(rownames(coefs), collapse=", "),
               ". Reporting first (", rownames(coefs)[1], ").") else NULL,
      extras = paste0("weight matrix estimator(s): ", paste(est, collapse=", "))
    ))
  }

  if (inherits(fit, c("glmMixBayes","survregMixBayes"))) {
    return(list(
      method_short = "Bayesian Mixture Model (adjMixBayes)",
      method_description = paste0(
        "Bayesian mixture modeling with MCMC posterior inference (Gutman et al., 2016). ",
        "Posterior samples via Stan (No-U-Turn Sampler). Reported coefficients are ",
        "posterior estimates for correctly matched records — already corrected."),
      estimator_note = NULL,
      extras = character(0)
    ))
  }

  list(method_short = "No adjustment (Naive)",
       method_description = "No linkage error correction applied.",
       estimator_note = NULL,
       extras = character(0))
}


# ── Internal: coefficient extraction ─────────────────────────────────────────

#' @noRd
.get_coef_df <- function(fit) {
  postlink_classes <- c("glmMixture","glmELE","glmMixBayes",
                        "coxphMixture","coxphELE","survregMixBayes","ctableMixture")

  if (is.matrix(fit$coefficients)) {
    est        <- fit$coefficients[1, ]
    term_names <- colnames(fit$coefficients)
  } else if (!any(sapply(postlink_classes, function(cl) inherits(fit, cl)))) {
    cf <- coef(summary(fit))
    return(data.frame(term=rownames(cf), estimate=round(as.numeric(cf[,1]),4),
                      se=round(as.numeric(cf[,2]),4),
                      row.names=NULL, stringsAsFactors=FALSE))
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

  data.frame(term=term_names, estimate=round(as.numeric(est),4),
             se=round(as.numeric(se),4), row.names=NULL, stringsAsFactors=FALSE)
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
#'   to the LLM to make the interpretation more relevant. Can be a short
#'   description or a path to a \code{.txt} file containing longer context
#'   (e.g. a variable dictionary).
#' @param provider Character string. LLM provider to use. One of
#'   \code{"anthropic"} (default if \code{ANTHROPIC_API_KEY} is set),
#'   \code{"openai"}, \code{"gemini"}, or \code{"ollama"} (local, no key
#'   needed). If \code{NULL}, auto-detected from available environment
#'   variables.
#' @param api_key Character string. API key for the chosen provider. If
#'   \code{NULL} (default), the key is read from the appropriate environment
#'   variable (\code{ANTHROPIC_API_KEY}, \code{OPENAI_API_KEY}, or
#'   \code{GEMINI_API_KEY}). Set the variable with
#'   \code{Sys.setenv(ANTHROPIC_API_KEY = "sk-ant-...")} or add it to
#'   \code{~/.Renviron}.
#' @param model Character string. Model name for the chosen provider.
#'   Defaults: Anthropic → \code{"claude-sonnet-4-6"}, OpenAI →
#'   \code{"gpt-4o"}, Gemini → \code{"gemini-1.5-flash"},
#'   Ollama → \code{"llama3"}.
#' @param verbose Logical. Print to console. Default \code{TRUE}.
#'
#' @return Invisibly returns the explanation as a character string.
#'
#' @details
#' For postlink objects, the LLM is told the coefficients are already
#' linkage-error-corrected and will not warn about bias that has been addressed.
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
#' explain(fit_adj)
#'
#' # Pass key directly
#' explain(fit_adj, api_key = "sk-ant-...")
#'
#' # OpenAI
#' explain(fit_adj, provider = "openai", api_key = "sk-...")
#'
#' # Local Ollama (no key needed)
#' explain(fit_adj, provider = "ollama", model = "llama3")
#'
#' # With context
#' explain(fit_adj,
#'   context = "LIFE-M linked Ohio birth and death certificates, 1883-1906.")
#'
#' # Multi-model comparison
#' fit_naive <- glm(age_at_death ~ poly(unit_yob, 3, raw = TRUE),
#'                  data = lifem, family = "gaussian")
#' explain(list(naive = fit_naive, adjusted = fit_adj))
#' }
#'
#' @export
explain <- function(fit,
                    context  = NULL,
                    provider = NULL,
                    api_key  = NULL,
                    model    = NULL,
                    verbose  = TRUE) {

  provider <- .detect_provider(provider)
  api_key  <- .resolve_key(api_key, provider)

  # Resolve context: file path or inline text
  if (!is.null(context) && file.exists(context)) {
    context <- paste(readLines(context, warn = FALSE), collapse = "\n")
  }

  is_fit_list <- is.list(fit) && !inherits(fit, c(
    "glm","lm","glmMixture","glmELE","glmMixBayes",
    "coxphMixture","coxphELE","survregMixBayes","ctableMixture"
  ))

  if (is_fit_list) {
    if (is.null(names(fit)) || any(names(fit) == ""))
      stop("When 'fit' is a list, all elements must be named.\n",
           "Example: list(naive = fit_glm, adjusted = fit_plglm)", call. = FALSE)
    return(.explain_multiple(fit, context=context, provider=provider,
                              api_key=api_key, model=model, verbose=verbose))
  }

  .explain_single(fit, context=context, provider=provider,
                  api_key=api_key, model=model, verbose=verbose)
}


#' @noRd
.explain_single <- function(fit, context=NULL, provider, api_key, model, verbose) {

  coef_df     <- .get_coef_df(fit)
  method_info <- .describe_method(fit)
  family_name <- tryCatch(fit$family$family, error=function(e) "unknown")
  family_label <- switch(family_name,
    gaussian="Gaussian (linear) regression", binomial="logistic regression",
    poisson="Poisson regression", Gamma="Gamma regression",
    paste0(family_name," regression"))

  coef_lines <- paste(apply(coef_df,1,function(r)
    sprintf("  %s:  estimate=%s, SE=%s (|z|=%.2f)",
            r["term"],r["estimate"],r["se"],
            abs(as.numeric(r["estimate"])/as.numeric(r["se"])))),
    collapse="\n")

  extras_block <- if (length(method_info$extras)>0)
    paste0("Adjustment details:\n",
           paste(paste0("  - ",method_info$extras),collapse="\n"),"\n\n") else ""
  estimator_note <- if (!is.null(method_info$estimator_note))
    paste0(method_info$estimator_note,"\n\n") else ""
  context_line <- if (!is.null(context))
    paste0("Study context: ",context,"\n\n") else ""

  system_msg <- paste0(
    "You are a statistician explaining post-linkage data analysis results. ",
    "CRITICAL: The coefficients are ALREADY ADJUSTED for record linkage errors. ",
    "Do NOT warn about mismatch bias — it has been corrected by postlink. ",
    "(1) Note which method was used and what it corrected for. ",
    "(2) Interpret the corrected coefficients in plain language. ",
    "(3) Note significance using |estimate/SE| > 2. ",
    "(4) No causation unless context supports it. (5) 3-5 sentences. No headers.")

  user_msg <- paste0(context_line,
    "Model: ",family_label,"\n\n",
    "Adjustment: ",method_info$method_short,"\n",
    "What it does: ",method_info$method_description,"\n\n",
    extras_block, estimator_note,
    "Corrected coefficients:\n",coef_lines,"\n\n",
    "Interpret these already-adjusted results.")

  message("Querying ", provider, " API...")
  explanation <- .call_llm(user_msg, system_msg, provider, api_key, model)

  if (verbose) {
    cat("\n--- postlink: explain() ---\n\n")
    cat(explanation)
    cat("\n\n[AI-generated. Verify before use.]\n")
  }
  invisible(explanation)
}


#' @noRd
.explain_multiple <- function(fit_list, context=NULL, provider, api_key, model, verbose) {

  context_line <- if (!is.null(context))
    paste0("Study context: ",context,"\n\n") else ""

  model_blocks <- lapply(names(fit_list), function(nm) {
    m       <- fit_list[[nm]]
    coef_df <- .get_coef_df(m)
    method  <- .describe_method(m)
    n_obs   <- tryCatch(nrow(m$model),
                        error=function(e) tryCatch(length(m$residuals),
                                                   error=function(e) NA))
    coef_lines <- paste(apply(coef_df,1,function(r)
      sprintf("  %s: estimate=%s, SE=%s",r["term"],r["estimate"],r["se"])),
      collapse="\n")
    paste0("Model: ",nm,"\nMethod: ",method$method_short,
           if (!is.na(n_obs)) paste0("\nn = ",format(n_obs,big.mark=",")) else "",
           "\nCoefficients:\n",coef_lines)
  })

  system_msg <- paste0(
    "You are a statistician comparing regression models from a post-linkage analysis. ",
    "CRITICAL: Coefficients from postlink-adjusted models are ALREADY CORRECTED for ",
    "linkage errors. Do NOT warn about bias in adjusted models. ",
    "(1) Briefly describe what each model does. ",
    "(2) Compare coefficients — note key differences. ",
    "(3) Explain what the differences show about linkage error correction. ",
    "(4) No causation unless context supports it. (5) 4-6 sentences. No headers.")

  user_msg <- paste0(context_line,
    length(fit_list)," models were fitted:\n\n",
    paste(model_blocks,collapse="\n\n---\n\n"),"\n\n",
    "Compare these models and explain what the differences reveal.")

  message("Querying ",provider," API (comparing ",length(fit_list)," models)...")
  explanation <- .call_llm(user_msg, system_msg, provider, api_key, model,
                            max_tokens=700)

  if (verbose) {
    cat("\n--- postlink: explain() — Model Comparison ---\n\n")
    cat(explanation)
    cat("\n\n[AI-generated. Verify before use.]\n")
  }
  invisible(explanation)
}
