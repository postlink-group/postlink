# mixbayes_survreg_helpers.R
# Internal helpers for the Bayesian survival mixture worker.

# Validate the `dist` string of survregMixBayes().
.validate_survreg_dist <- function(dist) {
 if (missing(dist) || !is.character(dist) || length(dist) != 1L || is.na(dist)) {
  stop("`dist` must be a single character string.", call. = FALSE)
 }
 tolower(trimws(dist))
}

# Normalise a right-censored survival response to list(time, event).
# Accepts a Surv object of type "right", a two-column matrix or data frame
# (time, event) or a list with elements `time` and `event`.
.normalize_surv_y <- function(y) {
 if (inherits(y, "Surv")) {
  type <- attr(y, "type")
  if (!identical(type, "right")) {
   stop(sprintf(paste0("Only right-censored survival responses Surv(time, event) are supported ",
                       "(got Surv type '%s')."), type), call. = FALSE)
  }
  time <- as.numeric(unclass(y)[, 1L])
  event <- unclass(y)[, 2L]
 } else if (is.list(y) && all(c("time", "event") %in% names(y))) {
  time <- as.numeric(y$time)
  event <- y$event
 } else if (is.matrix(y) || is.data.frame(y)) {
  if (ncol(y) != 2L) {
   stop("`y` must be a 2-column matrix (time, event), a Surv object or a list with time/event.",
        call. = FALSE)
  }
  time <- as.numeric(y[, 1L])
  event <- y[, 2L]
 } else {
  stop("`y` must be a 2-column matrix (time, event), a Surv object or a list with time/event.",
       call. = FALSE)
 }

 if (is.logical(event)) event <- as.integer(event)
 if (!is.numeric(event) || length(event) != length(time)) {
  stop("The event indicator must be numeric (0/1) or logical, with one value per survival time.",
       call. = FALSE)
 }
 if (anyNA(time) || anyNA(event)) stop("NA values found in y.", call. = FALSE)
 if (any(!is.finite(time)) || any(time <= 0)) stop("Survival times must be positive and finite.", call. = FALSE)
 if (!all(event %in% c(0, 1))) {
  stop("The event indicator must be coded 0 (right-censored) or 1 (event).", call. = FALSE)
 }
 list(time = time, event = as.integer(event))
}
