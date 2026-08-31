# =============================================================================
# report_gen.R
# Option 8: generate_report()
#
# fit parameter now accepts either:
#   - a single postlink fit object (shows that model only)
#   - a named list of fit objects (shows comparison table; the last postlink
#     object in the list is used as the primary adjusted model)
#
# comparison parameter removed — list handling built into fit.
#
# Output formats: "rmd" (default), "qmd", "html", "text"
# Requires: httr2 (optional — Suggests dependency)
# =============================================================================


#' Generate an Analysis Report from One or More Fitted Models
#'
#' Extracts numerical results from fitted model objects and produces a
#' structured analysis report with four sections: Data, Methods, Results,
#' and References.
#'
#' When \code{fit} is a named list, a comparison table is automatically
#' included in the Results section, and the last postlink object in the list
#' is used as the primary adjusted model for the Methods section and detailed
#' coefficient table.
#'
#' @param fit Either:
#'   \itemize{
#'     \item A single fitted postlink object (e.g. from \code{\link{plglm}}).
#'     \item A named list of fitted objects (mix of postlink and standard
#'       \code{glm}/\code{lm}). The last postlink object in the list is treated
#'       as the primary adjusted model. List names are used as column headers.
#'       Example: \code{list("Naive" = fit_glm, "Adjusted" = fit_plglm)}.
#'   }
#' @param output_file File path for the output including extension.
#'   Defaults to \code{"postlink_report.Rmd"}.
#' @param output_format One of \code{"rmd"} (default), \code{"qmd"},
#'   \code{"html"}, or \code{"text"}.
#' @param context An optional named list with study context. All fields are
#'   optional: \code{study}, \code{dataset}, \code{n_total},
#'   \code{outcome_label}, \code{predictor_label},
#'   \code{linkage_description}.
#' @param render Logical. Render after writing. Default \code{FALSE}.
#'
#' @return Invisibly returns the file path, or the report text if
#'   \code{output_format = "text"}.
#'
#' @examples
#' \dontrun{
#' postlink_set_key("sk-ant-...")
#'
#' data(lifem)
#' adj <- adjMixture(linked.data = lifem, m.formula = ~ commf + comml,
#'                   m.rate = 0.05, safe.matches = hndlnk)
#' fit_adj   <- plglm(age_at_death ~ poly(unit_yob, 3, raw = TRUE),
#'                    family = "gaussian", adjustment = adj)
#' fit_naive <- glm(age_at_death ~ poly(unit_yob, 3, raw = TRUE),
#'                  data = lifem, family = "gaussian")
#' fit_hl    <- glm(age_at_death ~ poly(unit_yob, 3, raw = TRUE),
#'                  data = lifem[lifem$hndlnk, ], family = "gaussian")
#'
#' # Single model
#' generate_report(fit_adj)
#'
#' # List of models — automatic comparison (matches vignette structure)
#' generate_report(
#'   list("Naive" = fit_naive, "Hand-Linked Only" = fit_hl, "Adjusted" = fit_adj),
#'   output_file = "lifem_report.Rmd",
#'   context = list(
#'     study         = "This study examines year of birth and age at death.",
#'     dataset       = "LIFE-M Ohio birth and death certificates, 1883-1906",
#'     n_total       = 156453,
#'     outcome_label = "age at death (years)",
#'     predictor_label = "year of birth"
#'   )
#' )
#' }
#'
#' @export
generate_report <- function(fit,
                             output_file   = NULL,
                             output_format = c("rmd", "qmd", "html", "text"),
                             context       = NULL,
                             provider      = NULL,
                             api_key       = NULL,
                             model         = NULL,
                             render        = FALSE) {

  output_format <- match.arg(output_format)

  if (!is.null(context) && !is.list(context))
    stop("'context' must be a named list or NULL.", call. = FALSE)

  # ── Detect single vs list ──────────────────────────────────────────────────
  is_fit_list <- is.list(fit) && !inherits(fit, c(
    "glm","lm","glmMixture","glmELE","glmMixBayes",
    "coxphMixture","coxphELE","survregMixBayes","ctableMixture"
  ))

  if (is_fit_list) {
    if (is.null(names(fit)) || any(names(fit) == ""))
      stop("When 'fit' is a list, all elements must be named.\n",
           "Example: list('Naive' = fit_glm, 'Adjusted' = fit_plglm)",
           call. = FALSE)

    # Primary model = last postlink object in the list
    postlink_classes <- c("glmMixture","glmELE","glmMixBayes",
                          "coxphMixture","coxphELE","survregMixBayes","ctableMixture")
    is_postlink <- sapply(fit, function(m) any(sapply(postlink_classes, function(cl) inherits(m, cl))))

    if (!any(is_postlink))
      stop("When 'fit' is a list, at least one element must be a fitted postlink object.",
           call. = FALSE)

    primary_idx  <- max(which(is_postlink))
    primary_fit  <- fit[[primary_idx]]
    fit_list     <- fit

  } else {
    primary_fit <- fit
    fit_list    <- NULL
  }

  # ── Extract info and build report ──────────────────────────────────────────
  info        <- .extract_fit_info(primary_fit, context, fit_list)
  provider <- if (!is.null(provider)) provider else .detect_provider(NULL)
  api_key  <- .resolve_key(api_key, provider)
  message("Querying ", provider, " API for plain-language summary...")
  ai_summary  <- .summarise_with_ai(info, context, fit_list, provider, api_key, model)
  report_text <- .build_report(info, ai_summary, output_format)

  if (output_format == "text") {
    cat(report_text)
    return(invisible(report_text))
  }

  if (is.null(output_file)) {
    output_file <- switch(output_format,
      rmd  = "postlink_report.Rmd",
      qmd  = "postlink_report.qmd",
      html = "postlink_report.html"
    )
  }

  out_dir <- dirname(output_file)
  if (!dir.exists(out_dir) && out_dir != ".") dir.create(out_dir, recursive = TRUE)
  writeLines(report_text, output_file)
  cat("\nReport written to:", normalizePath(output_file, mustWork = FALSE), "\n")

  if (render) {
    if (output_format == "rmd") {
      if (!requireNamespace("rmarkdown", quietly = TRUE)) {
        warning("'rmarkdown' not installed. Skipping render.", call. = FALSE)
      } else {
        message("Rendering...")
        rmarkdown::render(output_file, quiet = TRUE)
        cat("Rendered successfully.\n")
      }
    } else if (output_format == "qmd") {
      if (nchar(Sys.which("quarto")) == 0) {
        warning("Quarto CLI not found. Install from https://quarto.org", call. = FALSE)
      } else {
        message("Rendering...")
        system2("quarto", c("render", output_file))
      }
    }
  }

  cat("[Review all AI-generated text before submitting or publishing.]\n")
  invisible(output_file)
}


# ── Null-coalescing helper ────────────────────────────────────────────────────
`%||%` <- function(a, b) if (!is.null(a)) a else b


# ── Internal: extract info from primary fit ───────────────────────────────────
.extract_fit_info <- function(fit, context = NULL, fit_list = NULL) {

  is_mixture <- inherits(fit, c("glmMixture","coxphMixture","ctableMixture"))
  is_ele     <- inherits(fit, c("glmELE","coxphELE"))
  is_bayes   <- inherits(fit, c("glmMixBayes","survregMixBayes"))

  method_short <- if (is_mixture) "EM Mixture Model"
                  else if (is_ele) "ELE Weighting"
                  else if (is_bayes) "Bayesian Mixture Model"
                  else "postlink adjustment"

  ref <- if (is_mixture) {
    paste0("Slawski, M., West, B. T., Bukke, P., Wang, Z., Diao, G., & Ben-David, E. (2025). ",
           "A general framework for regression with mismatched data based on mixture modelling. ",
           "*Journal of the Royal Statistical Society Series A*, 188(3), 896-919.")
  } else if (is_ele) {
    paste0("Chambers, R. (2009). Regression analysis of probability-linked data. ",
           "*Official Statistics Research Series*, 4, 1-15.")
  } else {
    paste0("Gutman, R., Sammartino, C. J., Green, T. C., & Montague, B. T. (2016). ",
           "Error adjustments for file linking methods. ",
           "*Statistics in Medicine*, 35(1), 115-129.")
  }

  family_name  <- tryCatch(fit$family$family, error = function(e) "unknown")
  family_label <- switch(family_name,
    gaussian = "Gaussian (linear) regression",
    binomial = "logistic regression",
    poisson  = "Poisson regression",
    Gamma    = "Gamma regression",
    paste0(family_name, " regression")
  )

  outcome_label   <- context$outcome_label %||%
    tryCatch(as.character(fit$call$formula[[2]]), error = function(e) "outcome")
  predictor_label <- context$predictor_label %||% NULL
  formula_str     <- tryCatch(deparse(fit$call$formula), error = function(e) "")

  # Coefficient table — primary adjusted model only
  coef_df <- .get_coef_df(fit)
  coef_df$z_stat  <- round(coef_df$estimate / coef_df$se, 3)
  coef_df$p_value <- round(2 * pnorm(-abs(coef_df$z_stat)), 4)

  ci <- tryCatch(confint(fit), error = function(e) NULL)
  mcoef_df <- NULL
  if (!is.null(ci)) {
    coef_rows <- grepl("^coef ", rownames(ci))
    ci_coef   <- if (any(coef_rows)) ci[coef_rows, , drop=FALSE]
                 else ci[seq_len(nrow(coef_df)), , drop=FALSE]
    coef_df$ci_lower <- round(ci_coef[, 1], 4)
    coef_df$ci_upper <- round(ci_coef[, 2], 4)

    mcoef_rows <- grepl("^m\\.coef ", rownames(ci))
    if (any(mcoef_rows)) {
      ci_m <- ci[mcoef_rows, , drop=FALSE]
      mcoef_df <- data.frame(
        term     = sub("^m\\.coef ", "", rownames(ci_m)),
        ci_lower = round(ci_m[, 1], 4),
        ci_upper = round(ci_m[, 2], 4),
        row.names = NULL, stringsAsFactors = FALSE
      )
    }
  }

  n_fit   <- tryCatch(nrow(fit$model), error = function(e) NA_integer_)
  n_total <- context$n_total %||% n_fit

  match_summary <- NULL
  if (!is.null(fit$match.prob)) {
    mp <- fit$match.prob
    match_summary <- list(
      mean          = round(mean(mp), 3),
      pct_above_90  = round(100 * mean(mp > 0.9), 1),
      implied_mrate = round(100 * (1 - mean(mp)), 1)
    )
  }

  ele_info <- NULL
  if (is_ele && is.matrix(fit$coefficients))
    ele_info <- list(estimators   = rownames(fit$coefficients),
                     n_estimators = nrow(fit$coefficients))

  mixture_paradata <- NULL
  if (is_mixture) {
    mformula <- tryCatch(deparse(fit$call$adjustment$m.formula), error=function(e) NULL)
    mrate    <- tryCatch(as.character(fit$call$adjustment$m.rate), error=function(e) NULL)
    has_safe <- !is.null(tryCatch(fit$call$adjustment$safe.matches, error=function(e) NULL))
    mixture_paradata <- list(
      mformula       = mformula,
      mrate          = mrate,
      has_safe       = has_safe,
      intercept_only = (!is.null(mformula) && mformula %in% c("~1","~ 1"))
    )
  }

  # Comparison data from fit_list
  comparison_data <- NULL
  if (!is.null(fit_list)) {
    comparison_data <- lapply(names(fit_list), function(nm) {
      m   <- fit_list[[nm]]
      df  <- .get_coef_df(m)
      n_m <- tryCatch(nrow(m$model),
                      error = function(e) tryCatch(length(m$residuals),
                                                   error = function(e) NA_integer_))
      list(name = nm, terms = df$term, est = df$estimate, se = df$se, n = n_m)
    })
  }

  list(
    method_short     = method_short,
    is_mixture       = is_mixture,
    is_ele           = is_ele,
    is_bayes         = is_bayes,
    ref              = ref,
    family_name      = family_name,
    family_label     = family_label,
    outcome_label    = outcome_label,
    predictor_label  = predictor_label,
    formula_str      = formula_str,
    coef_df          = coef_df,
    mcoef_df         = mcoef_df,
    n_fit            = n_fit,
    n_total          = n_total,
    converged        = fit$converged,
    match_summary    = match_summary,
    ele_info         = ele_info,
    mixture_paradata = mixture_paradata,
    comparison_data  = comparison_data,
    has_context      = !is.null(context)
  )
}


# ── AI summary ────────────────────────────────────────────────────────────────
.summarise_with_ai <- function(info, context = NULL, fit_list = NULL, provider = NULL, api_key = NULL, model = NULL) {

  coef_lines <- paste(apply(info$coef_df, 1, function(r) {
    ci_part <- if (all(c("ci_lower","ci_upper") %in% names(r)))
      sprintf(", 95%% CI [%s, %s]", r["ci_lower"], r["ci_upper"]) else ""
    sprintf("  %s:  estimate = %s,  SE = %s,  p = %s%s",
            r["term"], r["estimate"], r["se"], r["p_value"], ci_part)
  }), collapse="\n")

  ctx_block <- ""
  if (info$has_context) {
    parts <- c(
      if (!is.null(context$study))               paste0("Study: ", context$study),
      if (!is.null(context$dataset))             paste0("Dataset: ", context$dataset),
      if (!is.null(context$linkage_description)) paste0("Linkage: ", context$linkage_description)
    )
    if (length(parts) > 0) ctx_block <- paste0(paste(parts, collapse="\n"), "\n\n")
  }

  method_context <- paste0(
    "Primary adjustment method: ", info$method_short, ". ",
    if (info$is_mixture)
      "The EM mixture model down-weights false links — adjusted coefficients reflect correctly matched records."
    else if (info$is_ele)
      "ELE weighting corrects for mismatch bias using block-level mismatch rates."
    else
      "Bayesian mixture model — posterior estimates for correctly matched records."
  )

  comp_block <- ""
  if (!is.null(info$comparison_data) && length(info$comparison_data) > 1) {
    comp_lines <- sapply(info$comparison_data, function(m) {
      ests <- paste(paste0(m$terms, "=", m$est), collapse=", ")
      paste0("  ", m$name, " (n=",
             if (!is.na(m$n)) format(m$n, big.mark=",") else "?", "): ", ests)
    })
    comp_block <- paste0(
      "For context, results from ", length(info$comparison_data), " models:\n",
      paste(comp_lines, collapse="\n"),
      "\nThe primary adjusted model is the last postlink object. ",
      "Please note key differences between naive/unadjusted and adjusted estimates.\n\n"
    )
  }

  system_msg <- paste0(
    "You are a statistician writing the Results section of an academic paper. ",
    "CRITICAL: The primary model coefficients are ALREADY ADJUSTED for record ",
    "linkage errors. Do NOT say estimates may be biased due to linkage — corrected. ",
    "If comparison models are provided, note how adjusted estimates differ from ",
    "naive or restricted estimates — this demonstrates the value of the correction. ",
    "Rules: (1) Only reference exact numbers provided. (2) No causation unless ",
    "context supports it. (3) Note p < 0.05 significance. (4) 3-5 sentences. ",
    "(5) No headers. (6) No disclaimer."
  )

  user_msg <- paste0(
    ctx_block,
    method_context, "\n\n",
    comp_block,
    "Primary adjusted model — corrected coefficients:\n", coef_lines, "\n\n",
    "Write a 3-5 sentence Results summary.",
    if (!is.null(info$comparison_data) && length(info$comparison_data) > 1)
      " Mention how adjusted estimates compare to naive." else ""
  )

  .call_llm(user_msg, system_msg, provider, api_key, model, max_tokens = 500)
}


# ── Build document ────────────────────────────────────────────────────────────
.build_report <- function(info, ai_summary, output_format) {

  date_str <- format(Sys.Date(), "%B %d, %Y")
  is_md    <- output_format %in% c("rmd","qmd","text")

  header <- switch(output_format,
    rmd = paste0("---\ntitle: \"postlink Analysis Report\"\ndate: \"", date_str, "\"\n",
                 "output:\n  html_document:\n    toc: true\n    toc_float: true\n",
                 "  pdf_document:\n    toc: true\n---\n\n"),
    qmd = paste0("---\ntitle: \"postlink Analysis Report\"\ndate: \"", date_str, "\"\n",
                 "format:\n  html:\n    toc: true\n  pdf:\n    toc: true\n---\n\n"),
    html = paste0(
      "<!DOCTYPE html><html><head><title>postlink Analysis Report</title><style>",
      "body{font-family:Georgia,serif;max-width:820px;margin:2.5em auto;",
      "line-height:1.7;color:#222}",
      "h1{font-size:1.8em;border-bottom:2px solid #333;padding-bottom:.3em}",
      "h2{font-size:1.3em;color:#444;margin-top:2em}h3{font-size:1.1em;color:#555}",
      "table{border-collapse:collapse;width:100%;margin:1em 0}",
      "th,td{border:1px solid #ccc;padding:8px 12px;text-align:left}",
      "th{background:#f5f5f5}",
      "blockquote{border-left:4px solid #f0ad4e;padding:.5em 1em;background:#fffbe6}",
      "code{background:#f4f4f4;padding:2px 5px;border-radius:3px}",
      "</style></head><body>\n<h1>postlink Analysis Report</h1>\n",
      "<p><em>", date_str, "</em></p>\n\n"
    ),
    text = paste0("postlink Analysis Report\n", strrep("=",40), "\n", date_str, "\n\n")
  )

  sec    <- function(title, body)
    if (is_md) paste0("## ", title, "\n\n", body, "\n")
    else       paste0("<h2>", title, "</h2>\n\n", body, "\n")
  subsec <- function(title)
    if (is_md) paste0("### ", title, "\n\n")
    else       paste0("<h3>", title, "</h3>\n")
  bold   <- function(x) if (is_md) paste0("**",x,"**") else paste0("<strong>",x,"</strong>")
  code   <- function(x) if (is_md) paste0("`",x,"`")   else paste0("<code>",x,"</code>")
  note   <- function(x)
    if (is_md) paste0("> ", x)
    else       paste0("<blockquote>",x,"</blockquote>")
  ref    <- function(x) if (is_md) paste0("- ",x) else paste0("<p>",x,"</p>")

  # Data
  n_str <- if (!is.na(info$n_total))
    paste0(format(info$n_total, big.mark=","), " records") else "N not available"

  if (info$has_context) {
    data_body <- paste0(
      if (!is.na(info$n_total)) paste0("The linked dataset contains ", bold(n_str), ". ") else "",
      if (!is.null(info$predictor_label))
        paste0("The outcome variable is ", bold(info$outcome_label),
               " and the main predictor is ", bold(info$predictor_label), ". ")
      else paste0("The outcome variable is ", bold(info$outcome_label), ". ")
    )
  } else {
    data_body <- paste0(
      "The analysis used ", n_str, ". ",
      "The outcome variable is ", bold(info$outcome_label), ". ",
      "Model formula: ", code(info$formula_str), "."
    )
  }

  # Methods
  methods_body <- .build_methods_text(info, bold, code, subsec, is_md)

  # Results
  disclaimer <- note(paste0(
    "The narrative below was drafted by an AI language model using exact numerical ",
    "outputs from the model. Review and edit before submitting or publishing."
  ))

  comp_block <- ""
  if (!is.null(info$comparison_data) && length(info$comparison_data) > 0) {
    comp_block <- paste0(
      subsec("Model Comparison"),
      "The table below compares coefficient estimates (SE) across all models. ",
      "The primary adjusted model corrects for record linkage error; differences ",
      "from the naive model reflect the impact of this correction.\n\n",
      .make_comparison_table(info$comparison_data, output_format),
      "\n\n"
    )
  }

  coef_tbl <- .make_coef_table(info$coef_df, output_format)
  sig_note  <- if (is_md) "_\\* p < 0.05_\n" else "<p><em>* p &lt; 0.05</em></p>"

  results_body <- paste0(
    disclaimer, "\n\n", ai_summary, "\n\n",
    comp_block,
    subsec("Adjusted Model — Coefficient Table"),
    coef_tbl, "\n\n", sig_note
  )

  # References
  refs_body <- paste0(
    ref(paste0("postlink R package (v0.1.0). Bukke, P., Kamat, G., Cui, J., ",
               "Gutman, R., & Slawski, M. https://postlink-group.github.io/postlink/")),
    "\n", ref(info$ref)
  )

  footer <- if (output_format == "html") "\n</body></html>" else ""

  paste0(header,
    sec("Data",       data_body),
    sec("Methods",    methods_body),
    sec("Results",    results_body),
    sec("References", refs_body),
    footer)
}


# ── Methods text (deterministic) ─────────────────────────────────────────────
.build_methods_text <- function(info, bold, code, subsec, is_md) {

  opening <- paste0(
    "Linkage error was accounted for using ", bold("postlink"),
    " (R package, v0.1.0; Bukke et al., 2026). ",
    "The outcome (", info$outcome_label, ") was modelled using ",
    info$family_label, " with formula ", code(info$formula_str), ". "
  )

  if (info$is_mixture) {
    pd <- info$mixture_paradata
    mp <- info$match_summary

    mismatch_desc <- if (!is.null(pd)) {
      if (isTRUE(pd$intercept_only))
        "The mismatch indicator model used an intercept-only specification (constant global mismatch rate). "
      else if (!is.null(pd$mformula))
        paste0("The probability of correct matching was modelled via logistic regression with covariates: ",
               code(pd$mformula), ". ")
      else ""
    } else ""

    mrate_desc <- if (!is.null(pd) && !is.null(pd$mrate) && pd$mrate != "NULL")
      paste0("A global mismatch rate constraint of ", pd$mrate, " was applied. ") else ""

    safe_desc <- if (!is.null(pd) && isTRUE(pd$has_safe))
      "Records designated as safe matches were fixed at posterior match probability 1. " else ""

    conv_desc <- if (!is.null(info$converged))
      if (isTRUE(info$converged)) "The EM algorithm converged. "
      else "Warning: the EM algorithm did not converge. " else ""

    mp_desc <- if (!is.null(mp))
      paste0("The mean posterior correct match probability was ", mp$mean,
             " (", mp$pct_above_90, "% of records exceeded 0.9), ",
             "implying an estimated mismatch rate of approximately ", mp$implied_mrate, "%. ")
    else ""

    method_para <- paste0(
      "Linkage error correction used an EM-based two-component mixture model ",
      "(Slawski et al., 2025). Each record pair is modelled as either a correct ",
      "match (C = 1) or a false link (C = 0). The EM algorithm alternates between ",
      "estimating posterior correct match probabilities (E-step) and updating the ",
      "regression parameters for true matches (M-step). ",
      mismatch_desc, mrate_desc, safe_desc, conv_desc, mp_desc
    )

    mcoef_block <- ""
    if (!is.null(info$mcoef_df)) {
      mcoef_block <- paste0(
        "\n\n", subsec("Mismatch Model Coefficients"),
        "95% confidence intervals for the mismatch model coefficients are shown ",
        "below. These describe how linkage quality covariates predict correct match ",
        "probability and are reported as supplementary information.\n\n",
        .make_mcoef_table(info$mcoef_df, if (is_md) "rmd" else "html")
      )
    }
    return(paste0(opening, method_para, mcoef_block))
  }

  if (info$is_ele) {
    ele <- info$ele_info
    est_desc <- if (!is.null(ele) && ele$n_estimators > 1)
      paste0(ele$n_estimators, " weighting estimators were computed (",
             paste(ele$estimators, collapse=", "),
             "). Results reported for the first (", ele$estimators[1], "). ")
    else if (!is.null(ele))
      paste0("Weighting estimator used: ", paste(ele$estimators, collapse=", "), ". ")
    else ""

    method_para <- paste0(
      "Linkage error correction used bias-adjusted estimating equations under the ",
      "Exchangeable Linkage Error (ELE) model (Chambers, 2009). Within each block, ",
      "mismatch errors are assumed exchangeable. Estimating equations are modified ",
      "using block-level mismatch rates; standard errors use a sandwich-type variance ",
      "estimator. ", est_desc
    )
    return(paste0(opening, method_para))
  }

  if (info$is_bayes) {
    method_para <- paste0(
      "Linkage error correction used a Bayesian mixture model (Gutman et al., 2016) ",
      "with MCMC posterior inference via Stan (No-U-Turn Sampler). Prior distributions ",
      "are specified for regression parameters among true links, false links, and the ",
      "marginal correct match probability. Reported estimates are posterior estimates ",
      "for correctly matched records. "
    )
    return(paste0(opening, method_para))
  }

  paste0(opening, "The postlink adjustment method was applied.")
}


# ── Comparison table ──────────────────────────────────────────────────────────
.make_comparison_table <- function(comparison_data, output_format) {

  is_md     <- output_format %in% c("rmd","qmd","text")
  all_terms <- comparison_data[[1]]$terms
  n_models  <- length(comparison_data)
  mod_names <- sapply(comparison_data, `[[`, "name")

  if (is_md) {
    hdr <- paste0("| Term | ", paste(mod_names, collapse=" | "), " |")
    sep <- paste0("|------|", paste(rep("------|", n_models), collapse=""))
    rows <- sapply(seq_along(all_terms), function(ti) {
      term  <- all_terms[ti]
      cells <- sapply(comparison_data, function(m) {
        idx <- which(m$terms == term)
        if (length(idx)==0) return("—")
        sprintf("%.4f (%.4f)", m$est[idx], m$se[idx])
      })
      paste0("| ", term, " | ", paste(cells, collapse=" | "), " |")
    })
    n_row <- paste0("| n | ",
                    paste(sapply(comparison_data, function(m)
                      if (!is.na(m$n)) format(m$n, big.mark=",") else "—"),
                      collapse=" | "), " |")
    paste(c(hdr, sep, rows, n_row), collapse="\n")
  } else {
    th <- function(x) paste0("<th>",x,"</th>")
    td <- function(x) paste0("<td>",x,"</td>")
    hcells <- paste0(th("Term"), paste(sapply(mod_names,th),collapse=""))
    rows_html <- sapply(seq_along(all_terms), function(ti) {
      term  <- all_terms[ti]
      cells <- sapply(comparison_data, function(m) {
        idx <- which(m$terms == term)
        if (length(idx)==0) return(td("—"))
        td(sprintf("%.4f (%.4f)", m$est[idx], m$se[idx]))
      })
      paste0("<tr>", td(term), paste(cells,collapse=""), "</tr>")
    })
    n_row_html <- paste0("<tr>", td("n"),
                         paste(sapply(comparison_data, function(m)
                           td(if (!is.na(m$n)) format(m$n,big.mark=",") else "—")),
                           collapse=""), "</tr>")
    paste0("<table>\n<thead><tr>",hcells,"</tr></thead>\n<tbody>\n",
           paste(c(rows_html,n_row_html),collapse="\n"),"\n</tbody>\n</table>")
  }
}


# ── Coefficient tables ────────────────────────────────────────────────────────
.make_coef_table <- function(coef_df, output_format) {
  has_ci <- all(c("ci_lower","ci_upper") %in% names(coef_df))
  if (output_format == "html") {
    th <- function(x) paste0("<th>",x,"</th>")
    td <- function(x) paste0("<td>",x,"</td>")
    hcells <- paste0(th("Term"),th("Estimate"),th("SE"),th("z"),th("p-value"),
                     if(has_ci) paste0(th("95% CI Lower"),th("95% CI Upper")) else "")
    rows <- paste(apply(coef_df,1,function(r){
      sig <- if(!is.na(as.numeric(r["p_value"]))&&as.numeric(r["p_value"])<0.05)" *" else ""
      paste0("<tr>",td(r["term"]),td(r["estimate"]),td(r["se"]),td(r["z_stat"]),
             td(paste0(r["p_value"],sig)),
             if(has_ci) paste0(td(r["ci_lower"]),td(r["ci_upper"])) else "","</tr>")
    }),collapse="\n")
    paste0("<table>\n<thead><tr>",hcells,"</tr></thead>\n<tbody>\n",rows,"\n</tbody>\n</table>")
  } else {
    if(has_ci){hdr<-"| Term | Estimate | SE | z | p-value | 95% CI Lower | 95% CI Upper |"
               sep<-"|------|----------|----|---|---------|--------------|--------------|"
    }else{hdr<-"| Term | Estimate | SE | z | p-value |";sep<-"|------|----------|----|---|---------|"}
    rows<-apply(coef_df,1,function(r){
      sig<-if(!is.na(as.numeric(r["p_value"]))&&as.numeric(r["p_value"])<0.05)" *" else ""
      if(has_ci) sprintf("| %s | %s | %s | %s | %s%s | %s | %s |",
                         r["term"],r["estimate"],r["se"],r["z_stat"],r["p_value"],sig,r["ci_lower"],r["ci_upper"])
      else sprintf("| %s | %s | %s | %s | %s%s |",r["term"],r["estimate"],r["se"],r["z_stat"],r["p_value"],sig)
    })
    paste(c(hdr,sep,rows),collapse="\n")
  }
}

.make_mcoef_table <- function(mcoef_df, output_format) {
  if (output_format == "html") {
    th<-function(x) paste0("<th>",x,"</th>"); td<-function(x) paste0("<td>",x,"</td>")
    hcells<-paste0(th("Term"),th("95% CI Lower"),th("95% CI Upper"))
    rows<-paste(apply(mcoef_df,1,function(r)
      paste0("<tr>",td(r["term"]),td(r["ci_lower"]),td(r["ci_upper"]),"</tr>")),collapse="\n")
    paste0("<table>\n<thead><tr>",hcells,"</tr></thead>\n<tbody>\n",rows,"\n</tbody>\n</table>")
  } else {
    hdr<-"| Term | 95% CI Lower | 95% CI Upper |";sep<-"|------|--------------|--------------|"
    rows<-apply(mcoef_df,1,function(r) sprintf("| %s | %s | %s |",r["term"],r["ci_lower"],r["ci_upper"]))
    paste(c(hdr,sep,rows),collapse="\n")
  }
}
