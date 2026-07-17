#!/usr/bin/env bash
set -eo pipefail

# ============================================================
# Check R package installation + library() health in Conda env
# Optimized version for MTD
# ============================================================

ENV_NAME="MTD"
MODE="isolated"
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
MTD_ROOT="$(cd -- "$SCRIPT_DIR/.." && pwd -P)"
PACKAGE_HEALTH_DIR="$MTD_ROOT/temp/package_health"
OUT_FILE="$PACKAGE_HEALTH_DIR/MTD_R_package_health.tsv"
SHOW_MODE="all"
MSG_WIDTH="80"
STRICT_WARNINGS="false"
IGNORE_BUILT_WARNINGS="false"

usage() {
  cat <<EOF
Usage:
  $0 [options]

Options:
  --env NAME                 Conda environment name [default: MTD]
  --isolated                 Test each package in a fresh R session [default]
  --fast                     Test all packages in the same R session
  --out FILE                 Save TSV report [default: MTD/temp/package_health/MTD_R_package_health.tsv]

  --show all                 Show all packages [default]
  --show warnings            Show only BUILD_WARN, LOAD_WARN, LOAD_FAIL, NOT_INSTALLED
  --show problems            Show only LOAD_FAIL and NOT_INSTALLED
  --show summary             Show only summary

  --width N                  Max message width in terminal table [default: 80]

  --ignore-built-warnings    Treat "built under R version" warnings as OK
  --strict-warnings          Exit with code 2 if BUILD_WARN or LOAD_WARN exists

  -h, --help                 Show this help

Examples:
  $0
  $0 --show warnings
  $0 --ignore-built-warnings
  $0 --fast --show problems
  $0 --env MTD --isolated --out report.tsv
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --env)
      ENV_NAME="${2:-}"
      shift 2
      ;;
    --isolated)
      MODE="isolated"
      shift
      ;;
    --fast)
      MODE="fast"
      shift
      ;;
    --out)
      OUT_FILE="${2:-}"
      shift 2
      ;;
    --show)
      SHOW_MODE="${2:-}"
      shift 2
      ;;
    --width)
      MSG_WIDTH="${2:-80}"
      shift 2
      ;;
    --ignore-built-warnings)
      IGNORE_BUILT_WARNINGS="true"
      shift
      ;;
    --strict-warnings)
      STRICT_WARNINGS="true"
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "ERROR: Unknown option: $1"
      usage
      exit 1
      ;;
  esac
done

case "$MODE" in
  isolated|fast) ;;
  *)
    echo "ERROR: invalid mode: $MODE"
    exit 1
    ;;
esac

case "$SHOW_MODE" in
  all|warnings|problems|summary) ;;
  *)
    echo "ERROR: invalid --show value: $SHOW_MODE"
    echo "Allowed: all, warnings, problems, summary"
    exit 1
    ;;
esac

if ! [[ "$MSG_WIDTH" =~ ^[0-9]+$ ]]; then
  echo "ERROR: --width must be an integer"
  exit 1
fi

# -------------------------
# Activate Conda env
# -------------------------
CONDA_SH=""

if [[ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]]; then
  CONDA_SH="$HOME/miniconda3/etc/profile.d/conda.sh"
elif [[ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]]; then
  CONDA_SH="$HOME/anaconda3/etc/profile.d/conda.sh"
else
  echo "ERROR: Could not find conda.sh."
  echo "Checked:"
  echo "  $HOME/miniconda3/etc/profile.d/conda.sh"
  echo "  $HOME/anaconda3/etc/profile.d/conda.sh"
  exit 1
fi

source "$CONDA_SH"
conda activate "$ENV_NAME"

if ! command -v Rscript >/dev/null 2>&1; then
  echo "ERROR: Rscript was not found after activating Conda env: $ENV_NAME"
  exit 1
fi

# -------------------------
# Main R script
# -------------------------
mkdir -p "$(dirname -- "$OUT_FILE")"

TMP_R="$(mktemp "${TMPDIR:-/tmp}/check_R_pkg_MTD.XXXXXX.R")"
trap 'rm -f "$TMP_R"' EXIT

cat > "$TMP_R" <<'RSCRIPT'
mode <- Sys.getenv("CHECK_MODE", "isolated")
out_file <- Sys.getenv("CHECK_OUT", "MTD_R_package_health.tsv")
env_name <- Sys.getenv("CONDA_DEFAULT_ENV", "unknown")
show_mode <- Sys.getenv("CHECK_SHOW", "all")
msg_width <- as.integer(Sys.getenv("CHECK_MSG_WIDTH", "80"))
strict_warnings <- identical(Sys.getenv("CHECK_STRICT_WARNINGS", "false"), "true")
ignore_built_warnings <- identical(Sys.getenv("CHECK_IGNORE_BUILT_WARNINGS", "false"), "true")

pkgs <- c(
  "textshaping",
  "ragg",
  "tidyverse",
  "car",
  "rstatix",
  "ggpubr",
  "plyr",
  "BiocManager",
  "rlang",
  "vctrs",
  "BiocGenerics",
  "S4Vectors",
  "IRanges",
  "UCSC.utils",
  "GenomeInfoDbData",
  "GenomeInfoDb",
  "matrixStats",
  "formatR",
  "lambda.r",
  "futile.options",
  "futile.logger",
  "RColorBrewer"
)

duplicates <- sort(unique(pkgs[duplicated(pkgs)]))
pkgs <- unique(pkgs)

installed <- installed.packages()
lib_paths <- .libPaths()

clip_msg <- function(x, n = 160) {
  x <- paste(x, collapse = " | ")
  x <- gsub("[\r\n\t]+", " ", x)
  x <- trimws(x)

  if (length(x) == 0 || is.na(x) || x == "") {
    return("")
  }

  if (nchar(x) > n) {
    paste0(substr(x, 1, n - 3), "...")
  } else {
    x
  }
}

is_built_warning <- function(warns) {
  grepl("was built under R version", warns, fixed = TRUE)
}

classify_load <- function(ok, warns, err_msg = "") {
  if (!ok) {
    return(list(
      status = "LOAD_FAIL",
      message = clip_msg(err_msg)
    ))
  }

  if (length(warns) == 0) {
    return(list(
      status = "OK",
      message = ""
    ))
  }

  built <- is_built_warning(warns)

  if (all(built)) {
    if (ignore_built_warnings) {
      return(list(
        status = "OK",
        message = ""
      ))
    }

    return(list(
      status = "BUILD_WARN",
      message = clip_msg(warns)
    ))
  }

  list(
    status = "LOAD_WARN",
    message = clip_msg(warns)
  )
}

test_load_fast <- function(pkg) {
  warns <- character()

  result <- tryCatch(
    withCallingHandlers(
      {
        suppressPackageStartupMessages(
          library(pkg, character.only = TRUE)
        )
        TRUE
      },
      warning = function(w) {
        warns <<- c(warns, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )

  if (identical(result, TRUE)) {
    classify_load(TRUE, warns)
  } else {
    classify_load(FALSE, warns, conditionMessage(result))
  }
}

parse_child_output <- function(output) {
  status_line <- grep("^__CHECK_STATUS__=", output, value = TRUE)
  message_line <- grep("^__CHECK_MESSAGE__=", output, value = TRUE)

  status <- if (length(status_line) > 0) {
    sub("^__CHECK_STATUS__=", "", status_line[[length(status_line)]])
  } else {
    "LOAD_FAIL"
  }

  message <- if (length(message_line) > 0) {
    sub("^__CHECK_MESSAGE__=", "", message_line[[length(message_line)]])
  } else {
    clip_msg(output)
  }

  list(status = status, message = message)
}

test_load_isolated <- function(pkg) {
  rscript_bin <- file.path(R.home("bin"), "Rscript")
  tmp_test <- tempfile(pattern = paste0("check_", pkg, "_"), fileext = ".R")
  on.exit(unlink(tmp_test), add = TRUE)

  code <- c(
    "args <- commandArgs(trailingOnly = TRUE)",
    "pkg <- args[[1]]",
    "ignore_built_warnings <- identical(Sys.getenv('CHECK_IGNORE_BUILT_WARNINGS', 'false'), 'true')",
    "",
    "clip_msg <- function(x, n = 160) {",
    "  x <- paste(x, collapse = ' | ')",
    "  x <- gsub('[\\r\\n\\t]+', ' ', x)",
    "  x <- trimws(x)",
    "  if (length(x) == 0 || is.na(x) || x == '') return('')",
    "  if (nchar(x) > n) paste0(substr(x, 1, n - 3), '...') else x",
    "}",
    "",
    "is_built_warning <- function(warns) {",
    "  grepl('was built under R version', warns, fixed = TRUE)",
    "}",
    "",
    "classify_load <- function(ok, warns, err_msg = '') {",
    "  if (!ok) return(list(status = 'LOAD_FAIL', message = clip_msg(err_msg)))",
    "  if (length(warns) == 0) return(list(status = 'OK', message = ''))",
    "  built <- is_built_warning(warns)",
    "  if (all(built)) {",
    "    if (ignore_built_warnings) return(list(status = 'OK', message = ''))",
    "    return(list(status = 'BUILD_WARN', message = clip_msg(warns)))",
    "  }",
    "  list(status = 'LOAD_WARN', message = clip_msg(warns))",
    "}",
    "",
    "warns <- character()",
    "",
    "res <- tryCatch(",
    "  withCallingHandlers(",
    "    {",
    "      suppressPackageStartupMessages(",
    "        library(pkg, character.only = TRUE)",
    "      )",
    "      TRUE",
    "    },",
    "    warning = function(w) {",
    "      warns <<- c(warns, conditionMessage(w))",
    "      invokeRestart('muffleWarning')",
    "    }",
    "  ),",
    "  error = function(e) e",
    ")",
    "",
    "ans <- if (identical(res, TRUE)) {",
    "  classify_load(TRUE, warns)",
    "} else {",
    "  classify_load(FALSE, warns, conditionMessage(res))",
    "}",
    "",
    "cat('__CHECK_STATUS__=', ans$status, '\\n', sep = '')",
    "cat('__CHECK_MESSAGE__=', ans$message, '\\n', sep = '')",
    "quit(save = 'no', status = 0, runLast = FALSE)"
  )

  writeLines(code, tmp_test)

  output <- suppressWarnings(
    tryCatch(
      system2(
        rscript_bin,
        args = c("--vanilla", tmp_test, pkg),
        stdout = TRUE,
        stderr = TRUE
      ),
      error = function(e) {
        structure(conditionMessage(e), status = 1)
      }
    )
  )

  status_code <- attr(output, "status")

  if (!is.null(status_code) && status_code != 0) {
    return(list(
      status = "LOAD_FAIL",
      message = clip_msg(output)
    ))
  }

  parsed <- parse_child_output(output)

  list(
    status = parsed$status,
    message = clip_msg(parsed$message)
  )
}

check_one <- function(pkg) {
  is_installed <- pkg %in% rownames(installed)

  if (!is_installed) {
    return(data.frame(
      Package = pkg,
      Installed = "no",
      Version = "-",
      Status = "NOT_INSTALLED",
      Message = "",
      stringsAsFactors = FALSE
    ))
  }

  version <- installed[pkg, "Version"]

  load_result <- if (identical(mode, "fast")) {
    test_load_fast(pkg)
  } else {
    test_load_isolated(pkg)
  }

  data.frame(
    Package = pkg,
    Installed = "yes",
    Version = version,
    Status = load_result$status,
    Message = load_result$message,
    stringsAsFactors = FALSE
  )
}

rows <- lapply(pkgs, check_one)
report <- do.call(rbind, rows)

write.table(
  report,
  file = out_file,
  sep = "\t",
  quote = TRUE,
  row.names = FALSE
)

summary_counts <- table(report$Status)

get_count <- function(name) {
  if (name %in% names(summary_counts)) {
    as.integer(summary_counts[[name]])
  } else {
    0L
  }
}

cat("\n")
cat("============================================================\n")
cat("MTD R package health check\n")
cat("============================================================\n")
cat("R version:              ", R.version.string, "\n", sep = "")
cat("Conda env:              ", env_name, "\n", sep = "")
cat("Check mode:             ", mode, "\n", sep = "")
cat("Packages listed:        ", length(pkgs), "\n", sep = "")
cat("Ignore built warnings:  ", ignore_built_warnings, "\n", sep = "")

if (length(duplicates) > 0) {
  cat("Duplicates removed:     ", paste(duplicates, collapse = ", "), "\n", sep = "")
}

cat("\nLibrary paths:\n")
for (p in lib_paths) {
  cat("  - ", p, "\n", sep = "")
}

cat("\nSummary:\n")
cat("  OK:             ", get_count("OK"), "\n", sep = "")
cat("  BUILD_WARN:     ", get_count("BUILD_WARN"), "\n", sep = "")
cat("  LOAD_WARN:      ", get_count("LOAD_WARN"), "\n", sep = "")
cat("  LOAD_FAIL:      ", get_count("LOAD_FAIL"), "\n", sep = "")
cat("  NOT_INSTALLED:  ", get_count("NOT_INSTALLED"), "\n", sep = "")

cat("\nReport saved to:\n")
cat("  ", normalizePath(out_file, mustWork = FALSE), "\n", sep = "")

cat("\nMeaning:\n")
cat("  OK             = installed and library() works\n")
cat("  BUILD_WARN     = package loads, only warning is 'built under R version'\n")
cat("  LOAD_WARN      = package loads, but produced another warning\n")
cat("  LOAD_FAIL      = installed, but library() failed\n")
cat("  NOT_INSTALLED  = package is missing\n")

display <- report

if (identical(show_mode, "warnings")) {
  display <- display[display$Status %in% c("BUILD_WARN", "LOAD_WARN", "LOAD_FAIL", "NOT_INSTALLED"), , drop = FALSE]
} else if (identical(show_mode, "problems")) {
  display <- display[display$Status %in% c("LOAD_FAIL", "NOT_INSTALLED"), , drop = FALSE]
} else if (identical(show_mode, "summary")) {
  display <- display[0, , drop = FALSE]
}

if (!identical(show_mode, "summary")) {
  cat("\n")

  if (nrow(display) == 0) {
    cat("No packages to show for --show ", show_mode, ".\n", sep = "")
  } else {
    display$Message <- vapply(display$Message, clip_msg, character(1), n = msg_width)

    cols <- c("Package", "Installed", "Version", "Status", "Message")
    display <- display[, cols, drop = FALSE]

    widths <- vapply(
      cols,
      function(col) {
        max(nchar(c(col, as.character(display[[col]]))), na.rm = TRUE)
      },
      integer(1)
    )

    make_line <- function() {
      paste0(
        "+",
        paste(
          vapply(
            widths,
            function(w) paste(rep("-", w + 2), collapse = ""),
            character(1)
          ),
          collapse = "+"
        ),
        "+"
      )
    }

    make_row <- function(values) {
      values <- as.character(values)
      paste0(
        "| ",
        paste(
          mapply(
            function(value, width) sprintf(paste0("%-", width, "s"), value),
            values,
            widths,
            USE.NAMES = FALSE
          ),
          collapse = " | "
        ),
        " |"
      )
    }

    cat(make_line(), "\n", sep = "")
    cat(make_row(cols), "\n", sep = "")
    cat(make_line(), "\n", sep = "")

    for (i in seq_len(nrow(display))) {
      cat(make_row(display[i, ]), "\n", sep = "")
    }

    cat(make_line(), "\n", sep = "")
  }
}

bad <- report$Status %in% c("LOAD_FAIL", "NOT_INSTALLED")
warn <- report$Status %in% c("BUILD_WARN", "LOAD_WARN")

if (any(bad)) {
  cat("\nProblematic packages detected.\n")
  quit(save = "no", status = 1, runLast = FALSE)
}

if (any(warn)) {
  if (strict_warnings) {
    cat("\nAll packages are installed and load, but warnings were detected.\n")
    quit(save = "no", status = 2, runLast = FALSE)
  } else {
    cat("\nAll packages are installed and load. Warnings were detected, but this is not treated as a failure.\n")
    quit(save = "no", status = 0, runLast = FALSE)
  }
}

cat("\nAll packages are installed and load correctly.\n")
quit(save = "no", status = 0, runLast = FALSE)
RSCRIPT

CHECK_MODE="$MODE" \
CHECK_OUT="$OUT_FILE" \
CHECK_SHOW="$SHOW_MODE" \
CHECK_MSG_WIDTH="$MSG_WIDTH" \
CHECK_STRICT_WARNINGS="$STRICT_WARNINGS" \
CHECK_IGNORE_BUILT_WARNINGS="$IGNORE_BUILT_WARNINGS" \
Rscript --vanilla "$TMP_R"
