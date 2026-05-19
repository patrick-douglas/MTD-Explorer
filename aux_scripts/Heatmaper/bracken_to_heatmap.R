#!/usr/bin/env Rscript

# bracken_to_heatmap_table_baseR.R
#
# Convert Bracken/MTD/DEG tables into a heatmap-ready matrix:
#   UNIQID    NAME    sample1    sample2    sample3 ...
#
# Base R only: no external R packages required.
# Designed as an Rscript replacement for bracken_to_heatmap_table_nopandas_EN.py.

ID_CANDIDATES <- c(
  "UNIQID", "uniqid", "id", "ID", "feature_id", "FeatureID",
  "taxid", "tax_id", "taxonomy_id", "taxonomyID", "NCBI_tax_id"
)

NAME_CANDIDATES <- c(
  "", "NAME", "Name", "name", "taxon", "Taxon", "taxonomy", "species",
  "Species", "clade_name", "Unnamed: 0", "#NAME", "# Name"
)

BRACKEN_VALUE_CANDIDATES <- c(
  "new_est_reads",
  "fraction_total_reads",
  "kraken_assigned_reads",
  "added_reads",
  "reads",
  "read_count",
  "count",
  "abundance"
)

DEFAULT_EXCLUDE_REGEX <- paste0(
  "(^$|baseMean|log2FoldChange|lfcSE|(^|[._-])stat([._-]|$)|",
  "pvalue|padj|qvalue|FDR|PValue|P.Value|",
  "significant|regulation|comparison|contrast|",
  "taxonomy_lvl|taxonomy_level|taxonomic_level|",
  "kraken_assigned_reads|added_reads|new_est_reads|fraction_total_reads)"
)

msg_info <- function(x) cat(sprintf("[INFO] %s\n", x), file = stderr())
msg_warn <- function(x) cat(sprintf("[WARNING] %s\n", x), file = stderr())
msg_error <- function(x, exit_code = 1) {
  cat(sprintf("[ERROR] %s\n", x), file = stderr())
  quit(status = exit_code, save = "no")
}

help_text <- function() {
  cat("bracken_to_heatmap_table_baseR.R\n\n")
  cat("Convert Bracken/MTD/DEG tables into a heatmap-ready matrix:\n")
  cat("  UNIQID    NAME    sample1    sample2    sample3 ...\n\n")
  cat("Base R only: no external R packages required.\n\n")
  cat("Usage:\n")
  cat("  Rscript bracken_to_heatmap_table_baseR.R --input FILE --output heatmap_table.txt\n")
  cat("  Rscript bracken_to_heatmap_table_baseR.R -i FILE -o heatmap_table.txt\n")
  cat("  Rscript bracken_to_heatmap_table_baseR.R -i sample1.bracken sample2.bracken -o heatmap.txt --value-col new_est_reads\n\n")
  cat("Examples:\n")
  cat("  Rscript bracken_to_heatmap_table_baseR.R \\\n    --input bracken_species_all_normalized.csv \\\n    --output heatmap_table.txt\n\n")
  cat("  Rscript bracken_to_heatmap_table_baseR.R \\\n    --input bracken_species_all_DEG.csv \\\n    --output heatmap_table_normtrans.txt \\\n    --prefer normtrans\n\n")
  cat("  Rscript bracken_to_heatmap_table_baseR.R \\\n    --input sample1.bracken sample2.bracken sample3.bracken \\\n    --output heatmap_from_reports.txt \\\n    --value-col new_est_reads\n\n")
  cat("  Rscript bracken_to_heatmap_table_baseR.R \\\n    --input bracken_species_all_DEG.csv \\\n    --output top50_heatmap.txt \\\n    --prefer normtrans \\\n    --top 50 \\\n    --transform row_zscore\n\n")
  cat("Options:\n")
  cat("  -i, --input              Input file(s). One wide matrix or multiple single-sample Bracken reports. Required.\n")
  cat("  -o, --output             Output file. Default: heatmap_table.txt\n")
  cat("      --sep                Input separator: auto, ',', ';', or '\\t'. Default: auto\n")
  cat("      --out-sep            Output separator. Use '\\t' for tab. Default: \\t\n")
  cat("      --id-col             ID column name. If missing, NAME is copied into UNIQID.\n")
  cat("      --name-col           Name/taxon column name. If missing, UNIQID is copied into NAME.\n")
  cat("      --value-col          For Bracken reports: abundance/value column. Default: auto\n")
  cat("      --prefer             For DEG/wide tables: auto, raw, norm, normtrans. Default: auto\n")
  cat("      --sample-regex       Regex to keep only selected sample columns. Example: 'LIVER|TEL'.\n")
  cat("      --exclude-regex      Regex for columns excluded automatically.\n")
  cat("      --fill-na            Value used to replace NA/non-numeric values. Default: 0\n")
  cat("      --transform          none, log2p1, row_zscore. Default: none\n")
  cat("      --min-sum            Remove features whose sample-value sum is lower than this.\n")
  cat("      --top                Keep only top N features by mean absolute value.\n")
  cat("      --sample-name-mode   keep or condition. condition renames LIVER/TEL to Liver/Telencephalon. Default: keep\n")
  cat("      --keep-suffix        Keep .norm and .normtrans suffixes in sample column names.\n")
  cat("  -h, --help               Show this help.\n")
}

is_option <- function(x) grepl("^-", x)

parse_args <- function(argv) {
  args <- list(
    input = character(0),
    output = "heatmap_table.txt",
    sep = "auto",
    out_sep = "\\t",
    id_col = NULL,
    name_col = NULL,
    value_col = "auto",
    prefer = "auto",
    sample_regex = NULL,
    exclude_regex = DEFAULT_EXCLUDE_REGEX,
    fill_na = 0,
    transform = "none",
    min_sum = NULL,
    top = NULL,
    sample_name_mode = "keep",
    keep_suffix = FALSE
  )

  if (length(argv) == 0) {
    help_text()
    quit(status = 0, save = "no")
  }

  i <- 1
  while (i <= length(argv)) {
    key <- argv[[i]]

    if (key %in% c("-h", "--help")) {
      help_text()
      quit(status = 0, save = "no")
    }

    if (key %in% c("-i", "--input")) {
      i <- i + 1
      vals <- character(0)
      while (i <= length(argv) && !is_option(argv[[i]])) {
        vals <- c(vals, argv[[i]])
        i <- i + 1
      }
      if (length(vals) == 0) msg_error("--input requires at least one file.")
      args$input <- vals
      next
    }

    if (key %in% c("--keep-suffix")) {
      args$keep_suffix <- TRUE
      i <- i + 1
      next
    }

    value_options <- c(
      "-o", "--output", "--sep", "--out-sep", "--id-col", "--name-col",
      "--value-col", "--prefer", "--sample-regex", "--exclude-regex",
      "--fill-na", "--transform", "--min-sum", "--top", "--sample-name-mode"
    )

    if (key %in% value_options) {
      if (i + 1 > length(argv)) msg_error(sprintf("%s requires a value.", key))
      val <- argv[[i + 1]]

      if (key %in% c("-o", "--output")) args$output <- val
      else if (key == "--sep") args$sep <- val
      else if (key == "--out-sep") args$out_sep <- val
      else if (key == "--id-col") args$id_col <- val
      else if (key == "--name-col") args$name_col <- val
      else if (key == "--value-col") args$value_col <- val
      else if (key == "--prefer") args$prefer <- val
      else if (key == "--sample-regex") args$sample_regex <- val
      else if (key == "--exclude-regex") args$exclude_regex <- val
      else if (key == "--fill-na") args$fill_na <- suppressWarnings(as.numeric(val))
      else if (key == "--transform") args$transform <- val
      else if (key == "--min-sum") args$min_sum <- suppressWarnings(as.numeric(val))
      else if (key == "--top") args$top <- suppressWarnings(as.integer(val))
      else if (key == "--sample-name-mode") args$sample_name_mode <- val

      i <- i + 2
      next
    }

    msg_error(sprintf("Unknown option: %s", key))
  }

  if (length(args$input) == 0) msg_error("Missing required --input.")
  if (!args$prefer %in% c("auto", "raw", "norm", "normtrans")) {
    msg_error("--prefer must be one of: auto, raw, norm, normtrans.")
  }
  if (!args$transform %in% c("none", "log2p1", "row_zscore")) {
    msg_error("--transform must be one of: none, log2p1, row_zscore.")
  }
  if (!args$sample_name_mode %in% c("keep", "condition")) {
    msg_error("--sample-name-mode must be one of: keep, condition.")
  }
  if (is.na(args$fill_na)) msg_error("--fill-na must be numeric.")
  if (!is.null(args$min_sum) && is.na(args$min_sum)) msg_error("--min-sum must be numeric.")
  if (!is.null(args$top) && is.na(args$top)) msg_error("--top must be an integer.")

  args
}

normalize_sep <- function(x) {
  if (identical(x, "\\t") || identical(x, "tab")) return("\t")
  x
}

detect_sep <- function(path, user_sep = "auto") {
  if (!identical(user_sep, "auto")) return(normalize_sep(user_sep))

  lines <- readLines(path, n = 50, warn = FALSE, encoding = "UTF-8")
  lines <- lines[nzchar(lines)]
  if (length(lines) == 0) msg_error(sprintf("Empty file: %s", path))

  first_line <- lines[[1]]
  seps <- c(",", "\t", ";")
  counts <- vapply(seps, function(s) lengths(regmatches(first_line, gregexpr(s, first_line, fixed = TRUE))), numeric(1))
  seps[[which.max(counts)]]
}

dedupe_headers <- function(headers) {
  headers <- trimws(as.character(headers))
  make.unique(headers, sep = ".")
}

read_any_table <- function(path, sep = "auto") {
  real_sep <- detect_sep(path, sep)

  raw <- tryCatch(
    utils::read.table(
      file = path,
      sep = real_sep,
      header = FALSE,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      quote = "\"",
      comment.char = "",
      fill = TRUE,
      blank.lines.skip = TRUE,
      colClasses = "character",
      na.strings = character(0)
    ),
    error = function(e) msg_error(sprintf("Could not read %s: %s", path, conditionMessage(e)))
  )

  if (nrow(raw) < 1) msg_error(sprintf("Empty file: %s", path))

  header <- dedupe_headers(as.character(raw[1, , drop = TRUE]))
  n <- length(header)

  if (nrow(raw) == 1) {
    df <- as.data.frame(matrix(character(0), nrow = 0, ncol = n), stringsAsFactors = FALSE)
  } else {
    df <- raw[-1, seq_len(n), drop = FALSE]
  }

  names(df) <- header

  if (nrow(df) > 0) {
    keep <- apply(df, 1, function(z) any(nzchar(trimws(as.character(z)))))
    df <- df[keep, , drop = FALSE]
  }

  list(columns = header, rows = df, sep = real_sep)
}

parse_float_vec <- function(value) {
  s <- trimws(as.character(value))
  s <- gsub('^"|"$', "", s)
  s <- gsub("^'|'$", "", s)

  bad <- s == "" | tolower(s) %in% c("na", "nan", "none", "null", "inf", "-inf")

  comma_decimal <- grepl("^-?\\d+,\\d+([eE][+-]?\\d+)?$", s)
  s[comma_decimal] <- gsub(",", ".", s[comma_decimal], fixed = TRUE)

  x <- suppressWarnings(as.numeric(s))
  x[bad | !is.finite(x)] <- NA_real_
  x
}

parse_float_one <- function(value) {
  x <- parse_float_vec(value)
  if (length(x) == 0 || is.na(x[[1]])) return(NA_real_)
  x[[1]]
}

format_number <- function(x) {
  x <- parse_float_vec(x)
  x[is.na(x)] <- 0
  x[abs(x) < 1e-12] <- 0
  out <- format(signif(x, 10), scientific = FALSE, trim = TRUE)
  out
}

find_column <- function(columns, candidates) {
  low <- tolower(columns)
  for (cand in candidates) {
    idx <- which(columns == cand)
    if (length(idx) > 0) return(columns[[idx[[1]]]])
    idx <- which(low == tolower(cand))
    if (length(idx) > 0) return(columns[[idx[[1]]]])
  }
  NULL
}

column_numeric_ratio <- function(rows, col) {
  vals <- rows[[col]]
  non_empty <- nzchar(trimws(as.character(vals)))
  n_non_empty <- sum(non_empty)
  if (n_non_empty == 0) return(c(non_empty = 0, numeric = 0, ratio = 0))
  numeric_ok <- !is.na(parse_float_vec(vals[non_empty]))
  n_numeric <- sum(numeric_ok)
  c(non_empty = n_non_empty, numeric = n_numeric, ratio = n_numeric / n_non_empty)
}

find_first_non_numeric_col <- function(columns, rows) {
  for (col in columns) {
    r <- column_numeric_ratio(rows, col)
    if (r[["non_empty"]] > 0 && r[["ratio"]] < 0.5) return(col)
  }
  NULL
}

choose_name_id_cols <- function(columns, rows, id_col = NULL, name_col = NULL) {
  if (!is.null(name_col)) {
    if (!name_col %in% columns) {
      msg_error(sprintf("--name-col '%s' does not exist. Available columns: %s", name_col, paste(columns, collapse = ", ")))
    }
    detected_name <- name_col
  } else {
    detected_name <- find_column(columns, NAME_CANDIDATES)
    if (is.null(detected_name)) detected_name <- find_first_non_numeric_col(columns, rows)
  }

  if (!is.null(id_col)) {
    if (!id_col %in% columns) {
      msg_error(sprintf("--id-col '%s' does not exist. Available columns: %s", id_col, paste(columns, collapse = ", ")))
    }
    detected_id <- id_col
  } else {
    detected_id <- find_column(columns, ID_CANDIDATES)
  }

  list(id = detected_id, name = detected_name)
}

clean_sample_name <- function(name, strip_suffix = TRUE) {
  n <- as.character(name)
  if (strip_suffix) {
    n <- sub("\\.normtrans$", "", n)
    n <- sub("\\.norm$", "", n)
  }
  n
}

condition_name_from_sample <- function(name) {
  n <- as.character(name)
  upper <- toupper(n)
  if (grepl("LIVER", upper)) return("Liver")
  if (grepl("(^|[_\\-.])TEL([_\\-.]|$)", upper) || grepl("TELENCEPHALON", upper)) return("Telencephalon")
  clean_sample_name(n, strip_suffix = TRUE)
}

detect_value_col <- function(columns, value_col = "auto") {
  low <- tolower(columns)
  if (!identical(value_col, "auto")) {
    idx <- which(columns == value_col)
    if (length(idx) > 0) return(columns[[idx[[1]]]])
    idx <- which(low == tolower(value_col))
    if (length(idx) > 0) return(columns[[idx[[1]]]])
    msg_error(sprintf("--value-col '%s' does not exist. Columns: %s", value_col, paste(columns, collapse = ", ")))
  }

  for (cand in BRACKEN_VALUE_CANDIDATES) {
    idx <- which(low == tolower(cand))
    if (length(idx) > 0) return(columns[[idx[[1]]]])
  }

  msg_error(sprintf(
    "Could not find a value column in the Bracken report. Try --value-col. Candidate columns: %s",
    paste(BRACKEN_VALUE_CANDIDATES, collapse = ", ")
  ))
}

is_bracken_single_report <- function(columns, rows, name_col, value_col) {
  low_cols <- tolower(columns)
  has_name <- !is.null(name_col) || any(tolower(NAME_CANDIDATES) %in% low_cols)

  if (!identical(value_col, "auto")) {
    has_value <- value_col %in% columns || tolower(value_col) %in% low_cols
  } else {
    has_value <- any(tolower(BRACKEN_VALUE_CANDIDATES) %in% low_cols)
  }

  numeric_cols <- 0
  for (col in columns) {
    r <- column_numeric_ratio(rows, col)
    if (r[["numeric"]] > 0) numeric_cols <- numeric_cols + 1
  }

  isTRUE(has_name && has_value && numeric_cols <= 8)
}

sample_name_from_path <- function(path) {
  base <- basename(path)
  base <- sub("\\.(bracken|tsv|txt|csv|report)$", "", base, ignore.case = TRUE)
  base <- sub("([._-]?bracken.*)$", "", base, ignore.case = TRUE)
  base <- sub("([._-]?kraken.*)$", "", base, ignore.case = TRUE)
  if (!nzchar(base)) basename(path) else base
}

numeric_candidate_columns <- function(columns, rows, ignore_cols = character(0), exclude_regex = DEFAULT_EXCLUDE_REGEX, sample_regex = NULL) {
  ignore_cols <- ignore_cols[!vapply(ignore_cols, is.null, logical(1))]
  candidates <- character(0)

  for (col in columns) {
    if (col %in% ignore_cols) next
    if (!is.null(exclude_regex) && nzchar(exclude_regex) && grepl(exclude_regex, col, ignore.case = TRUE, perl = TRUE)) next
    if (!is.null(sample_regex) && nzchar(sample_regex) && !grepl(sample_regex, col, perl = TRUE)) next

    r <- column_numeric_ratio(rows, col)
    if (r[["non_empty"]] > 0 && r[["ratio"]] >= 0.8) candidates <- c(candidates, col)
  }

  candidates
}

prefer_columns <- function(cols, prefer = "auto") {
  normtrans <- cols[grepl("\\.normtrans$", cols)]
  norm <- cols[grepl("\\.norm$", cols)]
  raw <- cols[!grepl("\\.norm(trans)?$", cols)]

  if (prefer == "normtrans") return(normtrans)
  if (prefer == "norm") return(norm)
  if (prefer == "raw") return(raw)

  if (length(normtrans) > 0) return(normtrans)
  if (length(norm) > 0) return(norm)
  raw
}

build_from_wide_matrix <- function(path, sep, id_col, name_col, prefer, sample_regex,
                                   exclude_regex, fill_na, strip_suffix, sample_name_mode) {
  tab <- read_any_table(path, sep)
  columns <- tab$columns
  rows <- tab$rows

  detected <- choose_name_id_cols(columns, rows, id_col = id_col, name_col = name_col)
  detected_id <- detected$id
  detected_name <- detected$name

  if (is.null(detected_name) && is.null(detected_id)) {
    msg_warn("No ID/name column was found. UNIQID/NAME will be created from row numbers.")
  }

  candidates <- numeric_candidate_columns(
    columns,
    rows,
    ignore_cols = c(detected_id, detected_name),
    exclude_regex = exclude_regex,
    sample_regex = sample_regex
  )

  chosen <- prefer_columns(candidates, prefer = prefer)
  if (length(chosen) == 0) {
    msg_error("Could not find numeric sample columns. Try --sample-regex, --prefer raw/norm/normtrans, or adjust --exclude-regex.")
  }

  if (sample_name_mode == "condition") {
    sample_cols <- make.unique(vapply(chosen, condition_name_from_sample, character(1)), sep = ".")
  } else {
    sample_cols <- make.unique(clean_sample_name(chosen, strip_suffix = strip_suffix), sep = ".")
  }

  n <- nrow(rows)
  if (!is.null(detected_id)) {
    uniqid <- trimws(as.character(rows[[detected_id]]))
  } else if (!is.null(detected_name)) {
    uniqid <- trimws(as.character(rows[[detected_name]]))
  } else {
    uniqid <- sprintf("feature_%d", seq_len(n))
  }

  if (!is.null(detected_name)) {
    name <- trimws(as.character(rows[[detected_name]]))
  } else if (!is.null(detected_id)) {
    name <- trimws(as.character(rows[[detected_id]]))
  } else {
    name <- uniqid
  }

  out <- data.frame(UNIQID = uniqid, NAME = name, stringsAsFactors = FALSE, check.names = FALSE)

  for (j in seq_along(chosen)) {
    vals <- parse_float_vec(rows[[chosen[[j]]]])
    vals[is.na(vals)] <- fill_na
    out[[sample_cols[[j]]]] <- vals
  }

  keep <- nzchar(trimws(out$UNIQID)) & nzchar(trimws(out$NAME))
  out[keep, , drop = FALSE]
}

build_from_single_bracken_reports <- function(paths, sep, id_col, name_col, value_col, fill_na) {
  sample_names <- make.unique(vapply(paths, sample_name_from_path, character(1)), sep = ".")

  key_order <- character(0)
  uniqid_by_key <- list()
  name_by_key <- list()
  values_by_key <- list()

  for (pidx in seq_along(paths)) {
    path <- paths[[pidx]]
    sample <- sample_names[[pidx]]

    tab <- read_any_table(path, sep)
    columns <- tab$columns
    rows <- tab$rows

    detected <- choose_name_id_cols(columns, rows, id_col = id_col, name_col = name_col)
    detected_id <- detected$id
    detected_name <- detected$name

    if (is.null(detected_name)) msg_error(sprintf("Could not find a name/taxon column in %s. Use --name-col.", path))
    val_col <- detect_value_col(columns, value_col)

    for (i in seq_len(nrow(rows))) {
      name <- trimws(as.character(rows[[detected_name]][[i]]))
      if (!nzchar(name)) next

      if (!is.null(detected_id) && !identical(detected_id, detected_name)) {
        uniqid <- trimws(as.character(rows[[detected_id]][[i]]))
        if (!nzchar(uniqid)) uniqid <- name
      } else {
        uniqid <- name
      }

      key <- paste(uniqid, name, sep = "\r")
      if (!key %in% key_order) {
        key_order <- c(key_order, key)
        uniqid_by_key[[key]] <- uniqid
        name_by_key[[key]] <- name
        values_by_key[[key]] <- setNames(rep(fill_na, length(sample_names)), sample_names)
      }

      x <- parse_float_one(rows[[val_col]][[i]])
      values_by_key[[key]][[sample]] <- ifelse(is.na(x), fill_na, x)
    }
  }

  out <- data.frame(
    UNIQID = vapply(key_order, function(k) uniqid_by_key[[k]], character(1)),
    NAME = vapply(key_order, function(k) name_by_key[[k]], character(1)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  for (sample in sample_names) {
    out[[sample]] <- vapply(key_order, function(k) as.numeric(values_by_key[[k]][[sample]]), numeric(1))
  }

  out
}

transform_matrix <- function(df, transform = "none") {
  if (transform == "none") return(df)
  if (ncol(df) <= 2) return(df)

  mat <- as.matrix(df[, -(1:2), drop = FALSE])
  storage.mode(mat) <- "numeric"
  mat[is.na(mat)] <- 0

  if (transform == "log2p1") {
    mat <- log2(pmax(mat, 0) + 1)
  } else if (transform == "row_zscore") {
    row_mean <- rowMeans(mat)
    if (ncol(mat) > 1) {
      row_sd <- apply(mat, 1, stats::sd)
    } else {
      row_sd <- rep(0, nrow(mat))
    }
    row_sd[is.na(row_sd) | row_sd == 0] <- NA_real_
    mat <- sweep(mat, 1, row_mean, "-")
    mat <- sweep(mat, 1, row_sd, "/")
    mat[is.na(mat)] <- 0
  } else {
    msg_error(sprintf("Unknown transformation: %s", transform))
  }

  df[, -(1:2)] <- as.data.frame(mat, stringsAsFactors = FALSE, check.names = FALSE)
  df
}

filter_and_sort <- function(df, min_sum = NULL, top = NULL) {
  if (ncol(df) <= 2) return(df)
  mat <- as.matrix(df[, -(1:2), drop = FALSE])
  storage.mode(mat) <- "numeric"
  mat[is.na(mat)] <- 0

  if (!is.null(min_sum)) {
    keep <- rowSums(mat) >= min_sum
    df <- df[keep, , drop = FALSE]
    mat <- mat[keep, , drop = FALSE]
  }

  if (!is.null(top) && top > 0 && nrow(df) > top) {
    score <- rowMeans(abs(mat))
    idx <- order(score, decreasing = TRUE)[seq_len(top)]
    df <- df[idx, , drop = FALSE]
  }

  df
}

write_heatmap_table <- function(path, df, out_sep = "\t") {
  out_sep <- normalize_sep(out_sep)
  if (ncol(df) > 2) {
    for (j in 3:ncol(df)) df[[j]] <- format_number(df[[j]])
  }
  utils::write.table(
    df,
    file = path,
    sep = out_sep,
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = "0"
  )
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))

  missing <- args$input[!file.exists(args$input)]
  if (length(missing) > 0) msg_error(sprintf("File(s) not found: %s", paste(missing, collapse = ", ")))

  if (length(args$input) > 1) {
    first_tab <- read_any_table(args$input[[1]], args$sep)
    first_detected <- choose_name_id_cols(first_tab$columns, first_tab$rows, id_col = args$id_col, name_col = args$name_col)
    use_bracken_merge <- is_bracken_single_report(first_tab$columns, first_tab$rows, first_detected$name, args$value_col)
  } else {
    use_bracken_merge <- FALSE
  }

  if (use_bracken_merge) {
    msg_info("Detected mode: multiple single-sample Bracken reports; merging by UNIQID/NAME.")
    df <- build_from_single_bracken_reports(
      paths = args$input,
      sep = args$sep,
      id_col = args$id_col,
      name_col = args$name_col,
      value_col = args$value_col,
      fill_na = args$fill_na
    )
  } else {
    if (length(args$input) > 1) {
      msg_warn("More than one input was provided, but they do not look like single-sample Bracken reports. Only the first file will be used.")
    }
    msg_info("Detected mode: wide/DEG matrix; selecting numeric sample columns.")
    df <- build_from_wide_matrix(
      path = args$input[[1]],
      sep = args$sep,
      id_col = args$id_col,
      name_col = args$name_col,
      prefer = args$prefer,
      sample_regex = args$sample_regex,
      exclude_regex = args$exclude_regex,
      fill_na = args$fill_na,
      strip_suffix = !args$keep_suffix,
      sample_name_mode = args$sample_name_mode
    )
  }

  df <- transform_matrix(df, args$transform)
  df <- filter_and_sort(df, min_sum = args$min_sum, top = args$top)
  write_heatmap_table(args$output, df, args$out_sep)

  sample_cols <- if (ncol(df) > 2) names(df)[-(1:2)] else character(0)
  msg_info(sprintf("Saved file: %s", args$output))
  msg_info(sprintf("Features/rows: %d", nrow(df)))
  msg_info(sprintf("Samples/numeric columns: %d", length(sample_cols)))
  if (length(sample_cols) > 0) {
    first <- paste(head(sample_cols, 8), collapse = ", ")
    if (length(sample_cols) > 8) first <- paste0(first, " ...")
    msg_info(sprintf("First sample columns: %s", first))
  }
}

main()

