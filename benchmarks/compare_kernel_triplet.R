#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(argv) {
  out <- list()
  if (length(argv) == 0) {
    return(out)
  }
  for (a in argv) {
    if (startsWith(a, "--")) {
      kv <- strsplit(sub("^--", "", a), "=", fixed = TRUE)[[1]]
      if (length(kv) == 2) {
        out[[kv[1]]] <- kv[2]
      } else if (length(kv) == 1) {
        out[[kv[1]]] <- TRUE
      }
    }
  }
  out
}

get_opt <- function(opts, key, default = NULL, required = FALSE) {
  val <- opts[[key]]
  if (is.null(val) || identical(val, "")) {
    if (required) {
      stop("Missing required option --", key, call. = FALSE)
    }
    return(default)
  }
  val
}

fmt <- function(x, digits = 6) {
  ifelse(is.finite(x), formatC(x, format = "f", digits = digits), "NA")
}

opts <- parse_args(args)

baseline_path <- get_opt(opts, "baseline", required = TRUE)
candidate_path <- get_opt(opts, "candidate", required = TRUE)
max_regression_pct <- as.numeric(get_opt(opts, "max-regression-pct", "20"))
max_speedup_drop_pct <- as.numeric(get_opt(opts, "max-speedup-drop-pct", "15"))
max_error_tol <- as.numeric(get_opt(opts, "max-error-tol", "1e-6"))
summary_file <- get_opt(opts, "summary-file", default = NULL)

if (!file.exists(baseline_path)) {
  stop("Baseline CSV not found: ", baseline_path, call. = FALSE)
}
if (!file.exists(candidate_path)) {
  stop("Candidate CSV not found: ", candidate_path, call. = FALSE)
}

baseline <- utils::read.csv(baseline_path, stringsAsFactors = FALSE, check.names = FALSE)
candidate <- utils::read.csv(candidate_path, stringsAsFactors = FALSE, check.names = FALSE)

required_cols <- c("kernel", "impl", "median_sec", "speedup_vs_reference", "max_abs_error")
missing_base <- setdiff(required_cols, colnames(baseline))
missing_cand <- setdiff(required_cols, colnames(candidate))
if (length(missing_base) > 0) {
  stop("Baseline CSV missing columns: ", paste(missing_base, collapse = ", "), call. = FALSE)
}
if (length(missing_cand) > 0) {
  stop("Candidate CSV missing columns: ", paste(missing_cand, collapse = ", "), call. = FALSE)
}

base_cur <- baseline[baseline$impl == "current", c("kernel", "median_sec", "speedup_vs_reference", "max_abs_error")]
cand_cur <- candidate[candidate$impl == "current", c("kernel", "median_sec", "speedup_vs_reference", "max_abs_error")]

colnames(base_cur) <- c("kernel", "baseline_median_sec", "baseline_speedup", "baseline_error")
colnames(cand_cur) <- c("kernel", "candidate_median_sec", "candidate_speedup", "candidate_error")

merged <- merge(base_cur, cand_cur, by = "kernel", all = FALSE)
if (nrow(merged) == 0) {
  stop("No overlapping 'current' kernels between baseline and candidate", call. = FALSE)
}

merged$median_change_pct <- 100 * (merged$candidate_median_sec / merged$baseline_median_sec - 1)
merged$speedup_change_pct <- 100 * (merged$candidate_speedup / merged$baseline_speedup - 1)

merged$fail_median <- is.finite(merged$median_change_pct) &
  merged$median_change_pct > max_regression_pct
merged$fail_speedup <- is.finite(merged$speedup_change_pct) &
  merged$speedup_change_pct < -max_speedup_drop_pct
merged$fail_error <- is.finite(merged$candidate_error) &
  merged$candidate_error > max_error_tol

merged$status <- ifelse(
  merged$fail_median | merged$fail_speedup | merged$fail_error,
  "FAIL",
  "OK"
)

has_fail <- any(merged$status == "FAIL")

cat("Kernel Triplet Comparison\n")
cat("=========================\n")
cat("baseline:  ", baseline_path, "\n", sep = "")
cat("candidate: ", candidate_path, "\n", sep = "")
cat(sprintf("thresholds: median_regression<=%.2f%%, speedup_drop<=%.2f%%, max_error<=%.3e\n\n",
            max_regression_pct, max_speedup_drop_pct, max_error_tol))

print(
  merged[, c(
    "kernel", "baseline_median_sec", "candidate_median_sec",
    "median_change_pct", "baseline_speedup", "candidate_speedup",
    "speedup_change_pct", "candidate_error", "status"
  )],
  row.names = FALSE
)

if (!is.null(summary_file)) {
  lines <- c(
    "# Kernel Triplet Comparison",
    "",
    sprintf("- baseline: `%s`", baseline_path),
    sprintf("- candidate: `%s`", candidate_path),
    sprintf(
      "- thresholds: median regression <= %.2f%%, speedup drop <= %.2f%%, max error <= %.3e",
      max_regression_pct, max_speedup_drop_pct, max_error_tol
    ),
    ""
  )

  header <- "| kernel | baseline_med_s | candidate_med_s | median_change_% | baseline_speedup | candidate_speedup | speedup_change_% | candidate_error | status |"
  sep <- "|---|---:|---:|---:|---:|---:|---:|---:|---|"
  rows <- vapply(seq_len(nrow(merged)), function(i) {
    sprintf(
      "| %s | %s | %s | %s | %s | %s | %s | %.3e | %s |",
      merged$kernel[i],
      fmt(merged$baseline_median_sec[i], 6),
      fmt(merged$candidate_median_sec[i], 6),
      fmt(merged$median_change_pct[i], 2),
      fmt(merged$baseline_speedup[i], 3),
      fmt(merged$candidate_speedup[i], 3),
      fmt(merged$speedup_change_pct[i], 2),
      merged$candidate_error[i],
      merged$status[i]
    )
  }, character(1))
  lines <- c(lines, header, sep, rows, "")
  if (has_fail) {
    lines <- c(lines, "**Result:** FAIL")
  } else {
    lines <- c(lines, "**Result:** PASS")
  }
  writeLines(lines, con = summary_file)
}

if (has_fail) {
  quit(status = 1L, save = "no")
}
