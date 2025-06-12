#!/usr/bin/env Rscript

# Dependency Analysis Script for fmriparametry
# This script analyzes function calls to identify potentially dead code

library(stringr)
library(purrr)

# Function to extract function definitions from R code
extract_function_defs <- function(file_content) {
  # Split into lines for better matching
  lines <- str_split(file_content, "\n")[[1]]
  
  function_names <- c()
  
  # Match function definitions with various patterns
  # Pattern 1: function_name <- function(
  # Pattern 2: .function_name <- function(
  # Pattern 3: function_name = function(
  patterns <- c(
    "^\\s*([\\._[:alnum:]]+)\\s*<-\\s*function\\s*\\(",
    "^\\s*([\\._[:alnum:]]+)\\s*=\\s*function\\s*\\(",
    "^([\\._[:alnum:]]+)\\s*<-\\s*function\\s*\\(",
    "^([\\._[:alnum:]]+)\\s*=\\s*function\\s*\\("
  )
  
  for (line in lines) {
    for (pattern in patterns) {
      if (str_detect(line, pattern)) {
        match <- str_match(line, pattern)
        if (!is.na(match[1])) {
          func_name <- match[2]
          if (!is.na(func_name) && nchar(func_name) > 0) {
            function_names <- c(function_names, func_name)
          }
        }
      }
    }
  }
  
  return(unique(function_names))
}

# Function to extract function calls from R code
extract_function_calls <- function(file_content, all_functions) {
  calls <- c()
  
  # Remove comments (multiline aware)
  lines <- str_split(file_content, "\n")[[1]]
  lines <- sapply(lines, function(line) {
    str_replace(line, "#.*$", "")
  })
  file_content <- paste(lines, collapse = "\n")
  
  for (func in all_functions) {
    # Escape special regex characters in function names
    escaped_func <- str_replace_all(func, "([\\.|\\+|\\*|\\?|\\[|\\]|\\{|\\}|\\(|\\)|\\^|\\$])", "\\\\\\1")
    
    # Look for function calls with various patterns
    patterns <- c(
      paste0("\\b", escaped_func, "\\s*\\("),          # func(
      paste0("\\s", escaped_func, "\\s*\\("),          # space before func(
      paste0("\\(", escaped_func, "\\s*\\("),          # (func(
      paste0("!", escaped_func, "\\s*\\("),            # !func(
      paste0("\\$", escaped_func, "\\s*\\("),          # $func(
      paste0("::", escaped_func, "\\s*\\("),           # ::func(
      paste0(":::", escaped_func, "\\s*\\(")           # :::func(
    )
    
    # Check if this function is called (not defined)
    for (pattern in patterns) {
      if (str_detect(file_content, pattern)) {
        # Make sure it's not a definition
        def_pattern <- paste0(escaped_func, "\\s*<-\\s*function")
        if (!str_detect(file_content, def_pattern)) {
          calls <- c(calls, func)
          break
        }
      }
    }
  }
  
  return(unique(calls))
}

# Get all R files
r_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
test_files <- list.files("tests/testthat", pattern = "^test-.*\\.R$", full.names = TRUE)
all_files <- c(r_files, test_files)

# First pass: collect all function definitions
all_functions <- c()
functions_by_file <- list()

cat("=== Phase 1: Collecting all function definitions ===\n\n")

for (file in r_files) {
  content <- paste(readLines(file, warn = FALSE), collapse = "\n")
  functions <- extract_function_defs(content)
  
  if (length(functions) > 0) {
    cat(sprintf("File: %s\n", basename(file)))
    cat(sprintf("  Functions defined: %s\n\n", paste(functions, collapse = ", ")))
    
    all_functions <- c(all_functions, functions)
    functions_by_file[[file]] <- functions
  }
}

all_functions <- unique(all_functions)
cat(sprintf("\nTotal unique functions found: %d\n", length(all_functions)))

# Second pass: find where each function is called
cat("\n=== Phase 2: Analyzing function usage ===\n\n")

function_calls <- list()
for (func in all_functions) {
  function_calls[[func]] <- list(called_in = c())
}

for (file in all_files) {
  content <- paste(readLines(file, warn = FALSE), collapse = "\n")
  calls <- extract_function_calls(content, all_functions)
  
  for (call in calls) {
    function_calls[[call]]$called_in <- c(function_calls[[call]]$called_in, file)
  }
}

# Analyze results
cat("\n=== Phase 3: Analysis Results ===\n\n")

# Find exported functions
namespace_file <- "NAMESPACE"
exported_functions <- c()
if (file.exists(namespace_file)) {
  namespace_content <- readLines(namespace_file, warn = FALSE)
  export_lines <- grep("^export\\(", namespace_content, value = TRUE)
  exported_functions <- str_extract(export_lines, "(?<=export\\()[^)]+")
}

# Find main entry points
main_functions <- c("estimate_parametric_hrf", exported_functions)

# Functions never called
never_called <- names(function_calls)[sapply(function_calls, function(x) length(x$called_in) == 0)]

# Functions only called in tests
only_in_tests <- names(function_calls)[sapply(function_calls, function(x) {
  all(str_detect(x$called_in, "^tests/"))
})]

# Functions with suspicious patterns
suspicious_patterns <- c("test_", "debug_", "old_", "_backup", "_deprecated", "_legacy")
suspicious_functions <- all_functions[sapply(all_functions, function(f) {
  any(str_detect(f, suspicious_patterns))
})]

# Output report
cat("=== POTENTIALLY DEAD CODE REPORT ===\n\n")

cat("1. NEVER CALLED (not found in any file):\n")
if (length(never_called) > 0) {
  for (func in never_called) {
    # Find which file defines it
    def_file <- names(functions_by_file)[sapply(functions_by_file, function(x) func %in% x)]
    cat(sprintf("   - %s (defined in: %s)\n", func, basename(def_file)))
  }
} else {
  cat("   None found\n")
}

cat("\n2. ONLY CALLED IN TESTS:\n")
test_only <- setdiff(only_in_tests, never_called)
if (length(test_only) > 0) {
  for (func in test_only) {
    def_file <- names(functions_by_file)[sapply(functions_by_file, function(x) func %in% x)]
    cat(sprintf("   - %s (defined in: %s)\n", func, basename(def_file)))
  }
} else {
  cat("   None found\n")
}

cat("\n3. SUSPICIOUS FUNCTION NAMES:\n")
if (length(suspicious_functions) > 0) {
  for (func in suspicious_functions) {
    def_file <- names(functions_by_file)[sapply(functions_by_file, function(x) func %in% x)]
    n_calls <- length(function_calls[[func]]$called_in)
    cat(sprintf("   - %s (defined in: %s, called %d times)\n", 
                func, basename(def_file), n_calls))
  }
} else {
  cat("   None found\n")
}

# Analyze specific files known to be problematic
problematic_files <- c(
  "performance_enhancements.R",
  "performance_optimizations.R", 
  "smart_performance_dispatcher.R",
  "rock-solid-memory.R",
  "rock-solid-numerical.R",
  "rock-solid-recovery.R",
  "rock-solid-validation.R"
)

cat("\n4. ANALYSIS OF SUSPECT FILES:\n")
for (file in problematic_files) {
  full_path <- file.path("R", file)
  if (file.exists(full_path)) {
    funcs <- functions_by_file[[full_path]]
    if (!is.null(funcs)) {
      cat(sprintf("\n   %s:\n", file))
      for (func in funcs) {
        n_calls <- length(function_calls[[func]]$called_in)
        if (n_calls == 0) {
          cat(sprintf("     - %s: NEVER CALLED\n", func))
        } else {
          calling_files <- unique(basename(function_calls[[func]]$called_in))
          cat(sprintf("     - %s: called %d times in: %s\n", 
                      func, n_calls, paste(calling_files, collapse = ", ")))
        }
      }
    }
  }
}

# Save detailed results
cat("\n\nSaving detailed results to analysis/dead_code_analysis.md\n")
dir.create("analysis", showWarnings = FALSE)

# Write detailed markdown report
report <- c(
  "# Dead Code Analysis Report",
  sprintf("Generated: %s", Sys.Date()),
  "",
  "## Summary",
  sprintf("- Total functions analyzed: %d", length(all_functions)),
  sprintf("- Never called: %d", length(never_called)),
  sprintf("- Only in tests: %d", length(test_only)),
  sprintf("- Suspicious names: %d", length(suspicious_functions)),
  "",
  "## Detailed Results",
  ""
)

# Add detailed function list
report <- c(report, "### Complete Function Analysis")
for (func in sort(all_functions)) {
  def_file <- names(functions_by_file)[sapply(functions_by_file, function(x) func %in% x)]
  n_calls <- length(function_calls[[func]]$called_in)
  
  status <- if (n_calls == 0) "**NEVER CALLED**" 
            else if (func %in% only_in_tests) "Test only"
            else if (func %in% exported_functions) "Exported"
            else "Internal"
  
  report <- c(report, sprintf("- `%s` (%s) - %s, %d calls", 
                              func, basename(def_file), status, n_calls))
}

writeLines(report, "analysis/dead_code_analysis.md")

cat("\nAnalysis complete!\n")