#!/usr/bin/env Rscript

# A script to automatically enforce Title Case in Roxygen2 headers
# while respecting stylized package names, classes, and acronyms.

format_title <- function(text) {
 # Apply base R's strict title casing
 cased_text <- tools::toTitleCase(text)

 # Define a dictionary of terms that must retain exact casing
 exception_terms <- c(
  "postlink", "survreg",
  "BRFSS", "ELE", "GLM", "CoxPH",
  "adjELE", "adjMixBayes", "adjMixture",
  "plcoxph", "coxphELE", "coxphMixture",
  "plglm", "glmELE", "glmMixBayes", "glmMixture",
  "plsurvreg", "survMixBayes", "survregMixBayes",
  "plctable", "ctableMixture",
  "mi_with"
 )

 # Loop through and force the exact casing for each protected term
 for (term in exception_terms) {
  # \\b ensures we only match whole words
  regex_pattern <- paste0("\\b", term, "\\b")
  cased_text <- gsub(regex_pattern, term, cased_text, ignore.case = TRUE)
 }

 return(cased_text)
}

# Scan all R files in the R/ directory
r_files <- list.files("R", pattern = "\\.[Rr]$", full.names = TRUE)

for (file in r_files) {
 lines <- readLines(file)
 in_roxygen <- FALSE
 file_modified <- FALSE

 for (i in seq_along(lines)) {
  is_roxygen <- grepl("^#'", lines[i])

  # Check for implicit titles (the very first line of a Roxygen block without a tag)
  if (is_roxygen && !in_roxygen) {
   if (grepl("^#'\\s+[a-zA-Z]", lines[i]) && !grepl("^#'\\s+@\\w+", lines[i])) {
    original_text <- sub("^#'\\s+", "", lines[i])
    new_text <- format_title(original_text)
    if (original_text != new_text) {
     lines[i] <- paste0("#' ", new_text)
     file_modified <- TRUE
    }
   }
   in_roxygen <- TRUE
  } else if (!is_roxygen) {
   in_roxygen <- FALSE
  }

  # Check for @title tags
  if (grepl("^#'\\s+@title\\s+", lines[i])) {
   original_text <- sub("^#'\\s+@title\\s+", "", lines[i])
   new_text <- format_title(original_text)
   if (original_text != new_text) {
    lines[i] <- paste0("#' @title ", new_text)
    file_modified <- TRUE
   }
  }
 }

 # Only write to the file if a change was actually made
 if (file_modified) {
  writeLines(lines, file)
  message(paste("Reformatted titles in:", file))
 }
}
