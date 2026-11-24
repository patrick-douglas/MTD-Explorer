#!/usr/bin/env Rscript

# ============================================================
# MTD exploratory: detected microbiome ranked tables
#
# Input:
#   prevalence_vs_abundance_<rank>_relative.tsv
#
# Taxon Prioritization Index (TPI); legacy importance_score is retained:
#   TPI =
#     0.45 * prevalence_score +
#     0.30 * mean_abundance_score +
#     0.15 * max_abundance_score +
#     0.10 * total_abundance_score
#
# Outputs:
#   detected_microbiome_<rank>_ranked_by_importance.tsv
#   top20_detected_microbiome_<rank>_by_importance.tsv
#   top50_detected_microbiome_<rank>_by_importance.tsv
#   detected_microbiome_<rank>_ranked_by_prevalence.tsv
#   detected_microbiome_<rank>_ranked_by_mean_abundance.tsv
#   detected_microbiome_<rank>_ranked_by_max_abundance.tsv
#   detected_microbiome_<rank>_ranked_by_total_abundance.tsv
#   detected_microbiome_<rank>_core_prevalence_75.tsv
#   detected_microbiome_<rank>_strict_core_prevalence_100.tsv
# ============================================================


args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  hit <- which(args == flag)
  if (length(hit) == 0) return(default)
  if (hit == length(args)) return(default)
  args[hit + 1]
}

has_flag <- function(flag) {
  flag %in% args
}

print_help <- function() {
  cat("
============================================================
MTD exploratory: detected microbiome ranked tables
============================================================

Usage:

  Rscript MTD_exploratory_detected_microbiome_rank.R \\
    --input prevalence_vs_abundance_species_relative.tsv \\
    --output MTD_res/exploratory/taxonomy/detected_microbiome/species_relative \\
    --rank species

Required:
  --input FILE
      Table produced by MTD_exploratory_prevalence_abundance.R

  --output DIR
      Output directory for ranked detected microbiome tables.

Optional:
  --rank STRING
      Taxonomic rank label. Default: taxa

  --abundance_input FILE
      Original combined Bracken abundance table.
      Example: bracken_species_all

  --samplesheet FILE
      Sample metadata table used to associate samples with
      groups, tissues or experimental conditions.

  --mode relative|absolute
      Abundance representation in the per-sample outputs.
      relative = percentage within each sample.
      absolute = estimated read counts.
      Default: relative

  --presence_threshold NUMERIC
      Minimum abundance required to consider a taxon present.
      Default: 0

  --w_prevalence NUMERIC
      Weight for prevalence score. Default: 0.45

  --w_mean NUMERIC
      Weight for mean abundance score. Default: 0.30

  --w_max NUMERIC
      Weight for maximum abundance score. Default: 0.15

  --w_total NUMERIC
      Weight for total abundance score. Default: 0.10

  --core_threshold NUMERIC
      Prevalence percentage used to define core taxa. Default: 75

  --plot_top INTEGER
      Number of top taxa to show in the TPI figures.
      Default: 20

  --width NUMERIC
      Plot width in inches.
      Default: 9

  --height NUMERIC
      Plot height in inches.
      Default: 7

  --dpi INTEGER
      PNG resolution.
      Default: 300

  --no_pdf
      Do not save PDF versions of the figures.

  --help
      Show help.

============================================================
\n")
}

if (has_flag("--help") || length(args) == 0) {
  print_help()
  quit(save = "no", status = 0)
}

input_file <- get_arg("--input")
output_dir <- get_arg("--output")
rank_label <- get_arg("--rank", "taxa")

# Optional original abundance table and sample metadata
abundance_input <- get_arg("--abundance_input", NULL)
samplesheet_file <- get_arg("--samplesheet", NULL)

# Abundance mode used in the per-sample output
mode <- get_arg("--mode", "relative")

# Minimum abundance required to consider a taxon present
presence_threshold <- as.numeric(
  get_arg("--presence_threshold", "0")
)

w_prevalence <- as.numeric(get_arg("--w_prevalence", "0.45"))
w_mean <- as.numeric(get_arg("--w_mean", "0.30"))
w_max <- as.numeric(get_arg("--w_max", "0.15"))
w_total <- as.numeric(get_arg("--w_total", "0.10"))

core_threshold <- as.numeric(get_arg("--core_threshold", "75"))
plot_top <- as.integer(get_arg("--plot_top", "20"))
plot_width <- as.numeric(get_arg("--width", "9"))
plot_height <- as.numeric(get_arg("--height", "7"))
plot_dpi <- as.integer(get_arg("--dpi", "300"))
save_pdf <- !has_flag("--no_pdf")
if (is.null(input_file)) stop("[ERROR] Missing --input")
if (is.null(output_dir)) stop("[ERROR] Missing --output")
if (!file.exists(input_file)) stop(paste0("[ERROR] Input not found: ", input_file))

if (!is.null(abundance_input) && !file.exists(abundance_input)) {
  stop(
    paste0(
      "[ERROR] Abundance input not found: ",
      abundance_input
    )
  )
}

if (!is.null(samplesheet_file) && !file.exists(samplesheet_file)) {
  stop(
    paste0(
      "[ERROR] Samplesheet not found: ",
      samplesheet_file
    )
  )
}

if (!mode %in% c("relative", "absolute")) {
  stop("[ERROR] --mode must be either 'relative' or 'absolute'")
}

if (!is.finite(presence_threshold) || presence_threshold < 0) {
  stop("[ERROR] --presence_threshold must be a numeric value >= 0")
}

weight_sum <- w_prevalence + w_mean + w_max + w_total

if (any(!is.finite(c(w_prevalence, w_mean, w_max, w_total))) ||
    any(c(w_prevalence, w_mean, w_max, w_total) < 0) ||
    !is.finite(weight_sum) || weight_sum <= 0) {
  stop("[ERROR] Invalid weights. Sum must be > 0.")
}

# Normalize weights in case user changes them and they do not sum to 1.
w_prevalence <- w_prevalence / weight_sum
w_mean <- w_mean / weight_sum
w_max <- w_max / weight_sum
w_total <- w_total / weight_sum

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
required_packages <- c("ggplot2")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      "[ERROR] Required R package not installed: ", pkg,
      "\nInstall with: install.packages('", pkg, "')"
    )
  }
}

suppressPackageStartupMessages({
  library(ggplot2)
})

has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)
cat("============================================================\n")
cat("MTD exploratory: detected microbiome ranked tables\n")
cat("Input :", input_file, "\n")
cat("Output:", output_dir, "\n")
cat("Rank  :", rank_label, "\n")
cat(
  "Abundance input:",
  ifelse(is.null(abundance_input), "not provided", abundance_input),
  "\n"
)
cat(
  "Samplesheet:",
  ifelse(is.null(samplesheet_file), "not provided", samplesheet_file),
  "\n"
)
cat("Per-sample abundance mode:", mode, "\n")
cat("Presence threshold:", presence_threshold, "\n")
cat("Weights:\n")
cat("  prevalence :", w_prevalence, "\n")
cat("  mean       :", w_mean, "\n")
cat("  max        :", w_max, "\n")
cat("  total      :", w_total, "\n")
cat("Core threshold:", core_threshold, "%\n")
cat("============================================================\n")

df <- read.table(
  input_file,
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required <- c(
  "taxon",
  "rank",
  "n_samples",
  "prevalence_n",
  "prevalence_percent",
  "mean_abundance",
  "median_abundance",
  "max_abundance",
  "total_abundance",
  "prevalence_class"
)

missing <- setdiff(required, colnames(df))

if (length(missing) > 0) {
  stop("[ERROR] Missing required columns: ", paste(missing, collapse = ", "))
}

to_num <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

# ------------------------------------------------------------
# Functions for per-sample abundance outputs
# ------------------------------------------------------------

detect_separator <- function(file) {
  first_line <- readLines(
    file,
    n = 1,
    warn = FALSE
  )

  n_tab <- lengths(
    regmatches(
      first_line,
      gregexpr("\t", first_line)
    )
  )

  n_comma <- lengths(
    regmatches(
      first_line,
      gregexpr(",", first_line)
    )
  )

  n_semicolon <- lengths(
    regmatches(
      first_line,
      gregexpr(";", first_line)
    )
  )

  if (
    n_tab >= n_comma &&
    n_tab >= n_semicolon &&
    n_tab > 0
  ) {
    return("\t")
  }

  if (n_semicolon > n_comma) {
    return(";")
  }

  return(",")
}


clean_input_column_names <- function(x) {
  x <- as.character(x)
  x <- gsub("\\s+", "_", x)
  x <- gsub("\\.+", "_", x)
  x
}


clean_sample_names <- function(x) {
  x <- as.character(x)

  # Remove Bracken abundance suffixes
  x <- gsub(
    "_(num|frac)$",
    "",
    x,
    perl = TRUE
  )

  # Remove report prefix
  x <- gsub(
    "^Report_",
    "",
    x,
    perl = TRUE
  )

  # Remove rank and Bracken suffix
  x <- gsub(
    "_(species|genus|family|order|class|phylum)_bracken$",
    "",
    x,
    perl = TRUE
  )

  # Remove possible dot-style Bracken suffix
  x <- gsub(
    "\\.(species|genus|family|order|class|phylum)\\.bracken$",
    "",
    x,
    perl = TRUE
  )

  # Remove common FASTQ/file extensions
  x <- gsub(
    "\\.(fastq|fq)(\\.gz)?$",
    "",
    x,
    ignore.case = TRUE,
    perl = TRUE
  )

  # Remove R-added X before names starting with a number
  x <- gsub(
    "^X(?=[0-9])",
    "",
    x,
    perl = TRUE
  )

  x
}


find_first_column <- function(candidates, available_columns) {
  candidate_lower <- tolower(candidates)
  available_lower <- tolower(available_columns)

  hit <- match(
    candidate_lower,
    available_lower
  )

  hit <- hit[!is.na(hit)]

  if (length(hit) == 0) {
    return(NULL)
  }

  available_columns[hit[1]]
}


find_taxon_column <- function(columns) {
  candidates <- c(
    "taxon",
    "taxa",
    "name",
    "scientific_name",
    "species",
    "genus",
    "family",
    "feature",
    "id"
  )

  detected <- find_first_column(
    candidates,
    columns
  )

  if (is.null(detected)) {
    detected <- columns[1]
  }

  detected
}


read_abundance_matrix <- function(
  abundance_file,
  abundance_mode = "relative"
) {
  separator <- detect_separator(
    abundance_file
  )

  abundance_df <- read.table(
    abundance_file,
    header = TRUE,
    sep = separator,
    quote = "\"",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  colnames(abundance_df) <- clean_input_column_names(
    colnames(abundance_df)
  )

  taxon_column <- find_taxon_column(
    colnames(abundance_df)
  )

  cat("[INFO] Abundance taxon column:", taxon_column, "\n")

  num_columns <- grep(
    "_num$",
    colnames(abundance_df),
    value = TRUE,
    perl = TRUE
  )

  frac_columns <- grep(
    "_frac$",
    colnames(abundance_df),
    value = TRUE,
    perl = TRUE
  )

  # Prefer _num columns because they contain estimated read counts.
  # Relative abundance is calculated below from those counts.
  if (length(num_columns) > 0) {
    selected_columns <- num_columns
    detected_type <- "Bracken _num columns"
  } else if (length(frac_columns) > 0) {
    selected_columns <- frac_columns
    detected_type <- "Bracken _frac columns"
  } else {
    excluded_columns <- c(
      taxon_column,
      "taxonomy_id",
      "taxid",
      "tax_id",
      "taxonomy_lvl",
      "rank",
      "level",
      "kraken_assigned_reads",
      "added_reads",
      "new_est_reads",
      "fraction_total_reads"
    )

    candidate_columns <- setdiff(
      colnames(abundance_df),
      excluded_columns
    )

    numeric_like <- vapply(
      abundance_df[
        ,
        candidate_columns,
        drop = FALSE
      ],
      function(z) {
        numeric_z <- suppressWarnings(
          as.numeric(
            as.character(z)
          )
        )

        mean(!is.na(numeric_z)) >= 0.8
      },
      logical(1)
    )

    selected_columns <- candidate_columns[
      numeric_like
    ]

    detected_type <- "generic numeric columns"
  }

  if (length(selected_columns) == 0) {
    stop(
      paste0(
        "[ERROR] No abundance columns detected in: ",
        abundance_file
      )
    )
  }

  sample_names <- clean_sample_names(
    selected_columns
  )

  if (anyDuplicated(sample_names)) {
    stop(
      paste0(
        "[ERROR] Duplicate biological sample names after cleaning: ",
        paste(
          unique(
            sample_names[
              duplicated(sample_names)
            ]
          ),
          collapse = ", "
        )
      )
    )
  }

  abundance_matrix <- as.matrix(
    data.frame(
      lapply(
        abundance_df[
          ,
          selected_columns,
          drop = FALSE
        ],
        to_num
      ),
      check.names = FALSE
    )
  )

  abundance_matrix[
    is.na(abundance_matrix)
  ] <- 0

  abundance_matrix[
    abundance_matrix < 0
  ] <- 0

  taxa <- as.character(
    abundance_df[[taxon_column]]
  )

  taxa[
    is.na(taxa) |
    taxa == ""
  ] <- "Unclassified_or_empty_name"

  rownames(abundance_matrix) <- taxa
  colnames(abundance_matrix) <- sample_names

  # Aggregate duplicate taxon names if they exist.
  if (anyDuplicated(rownames(abundance_matrix))) {
    abundance_matrix <- rowsum(
      abundance_matrix,
      group = rownames(abundance_matrix),
      reorder = FALSE
    )
  }

  abundance_matrix <- abundance_matrix[
    rowSums(abundance_matrix) > 0,
    ,
    drop = FALSE
  ]

  if (abundance_mode == "relative") {
    column_totals <- colSums(
      abundance_matrix,
      na.rm = TRUE
    )

    column_totals[
      column_totals == 0
    ] <- NA

    abundance_matrix <- sweep(
      abundance_matrix,
      2,
      column_totals,
      "/"
    ) * 100

    abundance_matrix[
      is.na(abundance_matrix)
    ] <- 0
  }

  cat("[INFO] Abundance input type:", detected_type, "\n")
  cat("[INFO] Biological samples:", ncol(abundance_matrix), "\n")
  cat(
    "[INFO] Sample names:",
    paste(
      colnames(abundance_matrix),
      collapse = ", "
    ),
    "\n"
  )

  abundance_matrix
}


read_sample_groups <- function(
  samplesheet,
  matrix_sample_names
) {
  output <- data.frame(
    sample = matrix_sample_names,
    group = "group_not_available",
    stringsAsFactors = FALSE
  )

  if (is.null(samplesheet)) {
    cat("[WARNING] No samplesheet provided. Group classification skipped.\n")
    return(output)
  }

  separator <- detect_separator(
    samplesheet
  )

  sample_df <- read.table(
    samplesheet,
    header = TRUE,
    sep = separator,
    quote = "\"",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  if (ncol(sample_df) == 0) {
    cat("[WARNING] Samplesheet contains no columns.\n")
    return(output)
  }

  sample_column <- find_first_column(
    c(
      "sample",
      "sample_id",
      "sampleid",
      "sample_name",
      "samplename",
      "library",
      "library_id",
      "name"
    ),
    colnames(sample_df)
  )

  group_column <- find_first_column(
    c(
      "group",
      "condition",
      "tissue",
      "organ",
      "treatment",
      "experimental_group",
      "sample_group"
    ),
    colnames(sample_df)
  )

  # Fallback: first column = sample, second column = group
  if (is.null(sample_column)) {
    sample_column <- colnames(sample_df)[1]
  }

  if (
    is.null(group_column) &&
    ncol(sample_df) >= 2
  ) {
    group_column <- colnames(sample_df)[2]
  }

  if (is.null(group_column)) {
    cat(
      "[WARNING] Group column could not be detected in samplesheet.\n"
    )
    return(output)
  }

  sample_ids <- clean_sample_names(
    sample_df[[sample_column]]
  )

  sample_groups <- as.character(
    sample_df[[group_column]]
  )

  group_map <- setNames(
    sample_groups,
    sample_ids
  )

  matched_groups <- unname(
    group_map[output$sample]
  )

  valid_match <- !is.na(matched_groups) &
    matched_groups != ""

  output$group[valid_match] <- matched_groups[
    valid_match
  ]

  cat("[INFO] Samplesheet sample column:", sample_column, "\n")
  cat("[INFO] Samplesheet group column:", group_column, "\n")
  cat(
    "[INFO] Samples matched to groups:",
    sum(valid_match),
    "of",
    nrow(output),
    "\n"
  )

  if (any(!valid_match)) {
    cat(
      "[WARNING] Samples without group match:",
      paste(
        output$sample[!valid_match],
        collapse = ", "
      ),
      "\n"
    )
  }

  output
}

num_cols <- c(
  "n_samples",
  "prevalence_n",
  "prevalence_percent",
  "mean_abundance",
  "median_abundance",
  "max_abundance",
  "total_abundance"
)

for (cc in num_cols) {
  df[[cc]] <- to_num(df[[cc]])
  df[[cc]][is.na(df[[cc]])] <- 0
}

minmax01 <- function(x) {
  x <- as.numeric(x)
  x[is.na(x)] <- 0

  if (length(x) == 0) {
    return(x)
  }

  xmin <- min(x, na.rm = TRUE)
  xmax <- max(x, na.rm = TRUE)

  if (!is.finite(xmin) || !is.finite(xmax) || xmax == xmin) {
    return(rep(0, length(x)))
  }

  (x - xmin) / (xmax - xmin)
}

# ------------------------------------------------------------
# Scores
# ------------------------------------------------------------

df$prevalence_score <- df$prevalence_percent / 100

# Log transform abundance before min-max scaling.
# This avoids one extreme taxon dominating the whole score.
df$mean_abundance_score <- minmax01(log10(df$mean_abundance + 1))
df$max_abundance_score <- minmax01(log10(df$max_abundance + 1))
df$total_abundance_score <- minmax01(log10(df$total_abundance + 1))

df$importance_score <-
  w_prevalence * df$prevalence_score +
  w_mean * df$mean_abundance_score +
  w_max * df$max_abundance_score +
  w_total * df$total_abundance_score

# Retain unrounded weighted contributions to reproduce the legacy rounded TPI.
df$tpi_prevalence_contribution <- w_prevalence * df$prevalence_score
df$tpi_mean_contribution <- w_mean * df$mean_abundance_score
df$tpi_max_contribution <- w_max * df$max_abundance_score
df$tpi_total_contribution <- w_total * df$total_abundance_score
df$importance_score <- round(df$importance_score, 6)
df$prevalence_score <- round(df$prevalence_score, 6)
df$mean_abundance_score <- round(df$mean_abundance_score, 6)
df$max_abundance_score <- round(df$max_abundance_score, 6)
df$total_abundance_score <- round(df$total_abundance_score, 6)

df$importance_interpretation <- ifelse(
  df$prevalence_percent >= core_threshold & df$mean_abundance_score >= 0.50,
  "consistent_and_abundant",
  ifelse(
    df$prevalence_percent >= core_threshold,
    "consistent_low_abundance",
    ifelse(
      df$prevalence_percent < 50 & df$max_abundance_score >= 0.75,
      "sample_specific_high_abundance",
      ifelse(
        df$prevalence_percent < 25 & df$mean_abundance_score < 0.25,
        "rare_low_abundance",
        "intermediate_priority"
      )
    )
  )
)

write_table <- function(x, filename) {
  write.table(
    x,
    file = file.path(output_dir, filename),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

make_ranked <- function(x, order_cols) {
  out <- x[do.call(order, order_cols), , drop = FALSE]
  out$rank_position <- seq_len(nrow(out))
  out
}

final_cols <- c(
  "rank_position",
  "taxon",
  "rank",
  "n_samples",
  "prevalence_n",
  "prevalence_percent",
  "mean_abundance",
  "median_abundance",
  "max_abundance",
  "total_abundance",
  "prevalence_class",
  "importance_interpretation",
  "importance_score",
  "prevalence_score",
  "mean_abundance_score",
  "max_abundance_score",
  "total_abundance_score"
)

# ------------------------------------------------------------
# Main rankings
# ------------------------------------------------------------

ranked_importance <- df[order(
  -df$importance_score,
  -df$prevalence_percent,
  -df$mean_abundance,
  -df$max_abundance,
  df$taxon
), , drop = FALSE]

ranked_importance$rank_position <- seq_len(nrow(ranked_importance))
ranked_importance <- ranked_importance[, final_cols]
# Canonical TPI exports. The legacy score and ranking are unchanged.
tpi_ranked <- ranked_importance
tpi_ranked$tpi_rank <- tpi_ranked$rank_position
tpi_ranked$tpi_score <- tpi_ranked$importance_score
tpi_ranked$abundance_mode <- rep(mode, nrow(tpi_ranked))
tpi_ranked$tpi_interpretation <- tpi_ranked$importance_interpretation
idx_tpi <- match(tpi_ranked$taxon, df$taxon)
if (anyNA(idx_tpi) || anyDuplicated(df$taxon)) stop("TPI requires unique matching taxon identifiers")
for (nm in c("tpi_prevalence_contribution", "tpi_mean_contribution", "tpi_max_contribution", "tpi_total_contribution")) {
  tpi_ranked[[nm]] <- df[[nm]][idx_tpi]
}
tpi_component_cols <- c("tpi_prevalence_contribution", "tpi_mean_contribution", "tpi_max_contribution", "tpi_total_contribution")
if (any(abs(rowSums(tpi_ranked[, tpi_component_cols, drop = FALSE]) - tpi_ranked$tpi_score) > 5.1e-7)) {
  stop("TPI components do not sum to the rounded legacy score")
}
write_table(tpi_ranked, paste0("detected_microbiome_", rank_label, "_ranked_by_TPI.tsv"))
for (n in c(20, 50)) write_table(head(tpi_ranked, n), paste0("top", n, "_detected_microbiome_", rank_label, "_by_TPI.tsv"))
write_table(tpi_ranked[, c("tpi_rank", "taxon", "abundance_mode", "tpi_score", tpi_component_cols)],
  paste0("detected_microbiome_", rank_label, "_TPI_components.tsv"))
write_table(data.frame(component = c("prevalence", "mean", "maximum", "total"),
  weight = c(w_prevalence, w_mean, w_max, w_total)),
  paste0("detected_microbiome_", rank_label, "_TPI_weights.tsv"))


write_table(
  ranked_importance,
  paste0("detected_microbiome_", rank_label, "_ranked_by_importance.tsv")
)

write_table(
  head(ranked_importance, 20),
  paste0("top20_detected_microbiome_", rank_label, "_by_importance.tsv")
)

write_table(
  head(ranked_importance, 50),
  paste0("top50_detected_microbiome_", rank_label, "_by_importance.tsv")
)

ranked_prevalence <- df[order(
  -df$prevalence_percent,
  -df$mean_abundance,
  -df$importance_score,
  df$taxon
), , drop = FALSE]

ranked_prevalence$rank_position <- seq_len(nrow(ranked_prevalence))
ranked_prevalence <- ranked_prevalence[, final_cols]

write_table(
  ranked_prevalence,
  paste0("detected_microbiome_", rank_label, "_ranked_by_prevalence.tsv")
)

write_table(
  head(ranked_prevalence, 20),
  paste0("top20_detected_microbiome_", rank_label, "_by_prevalence.tsv")
)

ranked_mean <- df[order(
  -df$mean_abundance,
  -df$prevalence_percent,
  -df$importance_score,
  df$taxon
), , drop = FALSE]

ranked_mean$rank_position <- seq_len(nrow(ranked_mean))
ranked_mean <- ranked_mean[, final_cols]

write_table(
  ranked_mean,
  paste0("detected_microbiome_", rank_label, "_ranked_by_mean_abundance.tsv")
)

write_table(
  head(ranked_mean, 20),
  paste0("top20_detected_microbiome_", rank_label, "_by_mean_abundance.tsv")
)

ranked_max <- df[order(
  -df$max_abundance,
  -df$prevalence_percent,
  -df$importance_score,
  df$taxon
), , drop = FALSE]

ranked_max$rank_position <- seq_len(nrow(ranked_max))
ranked_max <- ranked_max[, final_cols]

write_table(
  ranked_max,
  paste0("detected_microbiome_", rank_label, "_ranked_by_max_abundance.tsv")
)

write_table(
  head(ranked_max, 20),
  paste0("top20_detected_microbiome_", rank_label, "_by_max_abundance.tsv")
)

ranked_total <- df[order(
  -df$total_abundance,
  -df$prevalence_percent,
  -df$importance_score,
  df$taxon
), , drop = FALSE]

ranked_total$rank_position <- seq_len(nrow(ranked_total))
ranked_total <- ranked_total[, final_cols]

write_table(
  ranked_total,
  paste0("detected_microbiome_", rank_label, "_ranked_by_total_abundance.tsv")
)

write_table(
  head(ranked_total, 20),
  paste0("top20_detected_microbiome_", rank_label, "_by_total_abundance.tsv")
)

core <- ranked_importance[
  ranked_importance$prevalence_percent >= core_threshold,
  ,
  drop = FALSE
]

write_table(
  core,
  paste0("detected_microbiome_", rank_label, "_core_prevalence_", core_threshold, ".tsv")
)

strict_core <- ranked_importance[
  ranked_importance$prevalence_percent == 100,
  ,
  drop = FALSE
]

write_table(
  strict_core,
  paste0("detected_microbiome_", rank_label, "_strict_core_prevalence_100.tsv")
)

# ------------------------------------------------------------
# Per-sample abundance and group distribution
# ------------------------------------------------------------

if (!is.null(abundance_input)) {

  cat("============================================================\n")
  cat("[INFO] Generating per-sample abundance outputs\n")
  cat("Abundance input:", abundance_input, "\n")
  cat("Abundance mode :", mode, "\n")
  cat("Presence threshold:", presence_threshold, "\n")
  cat("============================================================\n")

  sample_matrix <- read_abundance_matrix(
    abundance_file = abundance_input,
    abundance_mode = mode
  )

  sample_groups <- read_sample_groups(
    samplesheet = samplesheet_file,
    matrix_sample_names = colnames(sample_matrix)
  )

  # Keep taxa contained in the main ranked table.
  ranked_taxa <- ranked_importance$taxon

  taxa_found <- ranked_taxa[
    ranked_taxa %in% rownames(sample_matrix)
  ]

  taxa_not_found <- setdiff(
    ranked_taxa,
    rownames(sample_matrix)
  )

  cat(
    "[INFO] Ranked taxa matched to abundance matrix:",
    length(taxa_found),
    "of",
    length(ranked_taxa),
    "\n"
  )

  if (length(taxa_not_found) > 0) {
    cat(
      "[WARNING] Ranked taxa not matched:",
      length(taxa_not_found),
      "\n"
    )
  }

  sample_matrix_ranked <- sample_matrix[
    taxa_found,
    ,
    drop = FALSE
  ]

  # ----------------------------------------------------------
  # Wide per-sample abundance table
  # ----------------------------------------------------------

  wide_abundance <- data.frame(
    taxon = rownames(sample_matrix_ranked),
    sample_matrix_ranked,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  wide_file <- paste0(
    "detected_microbiome_",
    rank_label,
    "_sample_abundance_wide.tsv"
  )

  write_table(
    wide_abundance,
    wide_file
  )

  # ----------------------------------------------------------
  # Long per-sample abundance table
  # ----------------------------------------------------------

  long_abundance <- data.frame(
    taxon = rep(
      rownames(sample_matrix_ranked),
      times = ncol(sample_matrix_ranked)
    ),
    sample = rep(
      colnames(sample_matrix_ranked),
      each = nrow(sample_matrix_ranked)
    ),
    abundance = as.vector(sample_matrix_ranked),
    stringsAsFactors = FALSE
  )

  group_map <- setNames(
    sample_groups$group,
    sample_groups$sample
  )

  long_abundance$group <- unname(
    group_map[long_abundance$sample]
  )

  long_abundance$group[
    is.na(long_abundance$group) |
    long_abundance$group == ""
  ] <- "group_not_available"

  long_abundance$present <- (
    long_abundance$abundance >
    presence_threshold
  )

  long_abundance <- long_abundance[
    ,
    c(
      "taxon",
      "sample",
      "group",
      "abundance",
      "present"
    )
  ]

  long_file <- paste0(
    "detected_microbiome_",
    rank_label,
    "_sample_abundance_long.tsv"
  )

  write_table(
    long_abundance,
    long_file
  )

  # ----------------------------------------------------------
  # Presence/absence matrix
  # ----------------------------------------------------------

  presence_matrix <- (
    sample_matrix_ranked >
    presence_threshold
  ) * 1

  presence_wide <- data.frame(
    taxon = rownames(presence_matrix),
    presence_matrix,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  presence_file <- paste0(
    "detected_microbiome_",
    rank_label,
    "_presence_absence.tsv"
  )

  write_table(
    presence_wide,
    presence_file
  )

  # ----------------------------------------------------------
  # Detection summary by samples and groups
  # ----------------------------------------------------------

  detection_summary <- lapply(
    rownames(sample_matrix_ranked),
    function(current_taxon) {

      taxon_values <- sample_matrix_ranked[
        current_taxon,
        ,
        drop = TRUE
      ]

      detected_samples <- names(taxon_values)[
        taxon_values > presence_threshold
      ]

      detected_groups <- unique(
        unname(
          group_map[detected_samples]
        )
      )

      detected_groups <- detected_groups[
        !is.na(detected_groups) &
        detected_groups != "" &
        detected_groups != "group_not_available"
      ]

      samples_text <- if (
        length(detected_samples) == 0
      ) {
        "not_detected"
      } else {
        paste(
          detected_samples,
          collapse = ";"
        )
      }

      groups_text <- if (
        length(detected_groups) == 0
      ) {
        "group_not_available"
      } else {
        paste(
          detected_groups,
          collapse = ";"
        )
      }

      distribution_scope <- if (
        length(detected_samples) == 0
      ) {
        "not_detected"
      } else if (
        length(detected_groups) == 0
      ) {
        "group_not_available"
      } else if (
        length(detected_groups) == 1
      ) {
        paste0(
          detected_groups[1],
          "_only"
        )
      } else {
        "shared_between_groups"
      }

      data.frame(
        taxon = current_taxon,
        samples_present = samples_text,
        groups_present = groups_text,
        n_samples_present = length(detected_samples),
        n_groups_present = length(detected_groups),
        distribution_scope = distribution_scope,
        stringsAsFactors = FALSE
      )
    }
  )

  detection_summary <- do.call(
    rbind,
    detection_summary
  )

  # ----------------------------------------------------------
  # Merge ranking, detection summary and sample abundances
  # ----------------------------------------------------------

  ranked_with_samples <- merge(
    ranked_importance,
    detection_summary,
    by = "taxon",
    all.x = TRUE,
    sort = FALSE
  )

  ranked_with_samples <- merge(
    ranked_with_samples,
    wide_abundance,
    by = "taxon",
    all.x = TRUE,
    sort = FALSE
  )

  # Restore the original importance ranking order.
  ranked_with_samples <- ranked_with_samples[
    match(
      ranked_importance$taxon,
      ranked_with_samples$taxon
    ),
    ,
    drop = FALSE
  ]

  ranked_with_samples_file <- paste0(
    "detected_microbiome_",
    rank_label,
    "_ranked_with_samples.tsv"
  )

  write_table(
    ranked_with_samples,
    ranked_with_samples_file
  )

  tpi_with_samples <- ranked_with_samples
  tpi_match <- match(tpi_with_samples$taxon, tpi_ranked$taxon)
  for (nm in setdiff(names(tpi_ranked), names(tpi_with_samples))) {
    tpi_with_samples[[nm]] <- tpi_ranked[[nm]][tpi_match]
  }
  write_table(tpi_with_samples, paste0("detected_microbiome_", rank_label, "_TPI_ranked_with_samples.tsv"))
  cat("[OK] Ranked table with sample abundances:\n")
  cat(
    file.path(
      output_dir,
      ranked_with_samples_file
    ),
    "\n"
  )

  cat("[OK] Wide sample abundance table:\n")
  cat(
    file.path(
      output_dir,
      wide_file
    ),
    "\n"
  )

  cat("[OK] Long sample abundance table:\n")
  cat(
    file.path(
      output_dir,
      long_file
    ),
    "\n"
  )

  cat("[OK] Presence/absence matrix:\n")
  cat(
    file.path(
      output_dir,
      presence_file
    ),
    "\n"
  )
}

# ------------------------------------------------------------
# README / methods note
# ------------------------------------------------------------

readme_file <- file.path(
  output_dir,
  paste0("README_detected_microbiome_", rank_label, ".txt")
)

cat(
  "Detected microbiome ranked tables\n",
  "==================================\n\n",
  "Input table:\n",
  input_file, "\n\n",
  "Taxon Prioritization Index (TPI); legacy importance_score is identical.\n",
  "Abundance mode: ", mode, "\n",
  "TPI = ",
  round(w_prevalence, 4), " * prevalence_score + ",
  round(w_mean, 4), " * mean_abundance_score + ",
  round(w_max, 4), " * max_abundance_score + ",
  round(w_total, 4), " * total_abundance_score\n\n",
  "Metric definitions:\n",
  "prevalence_score = prevalence_percent / 100\n",
  "mean_abundance_score = min-max scaled log10(mean_abundance + 1)\n",
  "max_abundance_score = min-max scaled log10(max_abundance + 1)\n",
  "total_abundance_score = min-max scaled log10(total_abundance + 1)\n\n",
  "Interpretation:\n",
  "High importance taxa are those that are both consistently detected across samples and quantitatively abundant.\n",
  "High prevalence but low abundance taxa may represent stable low-abundance community members.\n",
  "Low prevalence but high maximum abundance taxa may represent sample-specific blooms or possible contaminants.\n",
  "Low prevalence and low abundance taxa should be interpreted cautiously.\n\n",
  "TPI is a descriptive, dataset- and mode-specific prioritization score, not a significance test,\n",
  "a measure of pathogenicity, or proof of biological interactions. Mean and total are dependent.\n",
  "Relative abundance is within-library Bracken composition, not absolute microbial burden.\n",
  "Weighted components are exported before rounding; TPI is rounded to six decimals.\n\n",
  "Generated files (legacy names retained):\n",
  "- detected_microbiome_", rank_label, "_ranked_by_TPI.tsv\n",
  "- detected_microbiome_", rank_label, "_TPI_components.tsv\n",
  "- detected_microbiome_", rank_label, "_TPI_weights.tsv\n",
  "- top20/top50_detected_microbiome_", rank_label, "_by_TPI.tsv\n",
  "- TPI lollipop, components, landscape and overview figures (if plot_top > 0).\n",
  "- detected_microbiome_", rank_label, "_ranked_by_importance.tsv\n",
  "- top20_detected_microbiome_", rank_label, "_by_importance.tsv\n",
  "- top50_detected_microbiome_", rank_label, "_by_importance.tsv\n",
  "- detected_microbiome_", rank_label, "_ranked_by_prevalence.tsv\n",
  "- detected_microbiome_", rank_label, "_ranked_by_mean_abundance.tsv\n",
  "- detected_microbiome_", rank_label, "_ranked_by_max_abundance.tsv\n",
  "- detected_microbiome_", rank_label, "_ranked_by_total_abundance.tsv\n",
  "- detected_microbiome_", rank_label, "_core_prevalence_", core_threshold, ".tsv\n",
  "- detected_microbiome_", rank_label, "_strict_core_prevalence_100.tsv\n",
  "- detected_microbiome_", rank_label, "_core_prevalence_", core_threshold, ".tsv\n",
  "- detected_microbiome_", rank_label, "_strict_core_prevalence_100.tsv\n",
  "- detected_microbiome_", rank_label, "_ranked_with_samples.tsv\n",
  "- detected_microbiome_", rank_label, "_sample_abundance_wide.tsv\n",
  "- detected_microbiome_", rank_label, "_sample_abundance_long.tsv\n",
  "- detected_microbiome_", rank_label, "_presence_absence.tsv\n",
  sep = "",
  file = readme_file
)

# ------------------------------------------------------------
# Figures
# ------------------------------------------------------------

# TPI plots: same scores and ranking as the legacy output.
if (plot_top > 0 && nrow(tpi_ranked) > 0) {
  top_plot <- head(tpi_ranked, plot_top)
  n_top <- nrow(top_plot)
  levels_top <- rev(top_plot$taxon)
  top_plot$taxon <- factor(top_plot$taxon, levels = levels_top)
  taxon_labels <- setNames(paste0(top_plot$tpi_rank, "  ",
    vapply(as.character(top_plot$taxon), function(z) {
      paste(strwrap(z, width = 52), collapse = "\n")
    }, character(1))), as.character(top_plot$taxon))
  mode_label <- if (mode == "relative") "Relative abundance" else "Estimated read counts"
  mean_label <- if (mode == "relative") "Mean relative abundance (%)" else "Mean estimated read count"
  max_label <- if (mode == "relative") "Maximum relative\nabundance (%)" else "Maximum estimated\nread count"
  prefix <- paste0("detected_microbiome_", rank_label, "_TPI")
  top_prefix <- paste0("top", plot_top, "_detected_microbiome_", rank_label)
  tpi_theme <- theme_bw(base_size = 11) + theme(
    plot.title = element_text(face = "bold", color = "#203642", size = 13),
    plot.subtitle = element_text(color = "#536777", size = 10),
    panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "#e5ebef"),
    panel.border = element_blank(), axis.line = element_line(color = "#9aabb4"),
    legend.key = element_blank(), plot.margin = margin(8, 12, 8, 8))
  score_scale <- function() scale_x_continuous(limits = c(0, 1), breaks = c(0, .5, 1),
    expand = expansion(mult = c(0, .04)))
  taxon_scale <- function() scale_y_discrete(limits = levels_top, drop = FALSE,
    expand = expansion(add = .6))
  prevalence_scale <- function() scale_size_area(max_size = 4, limits = c(0, 100),
    breaks = c(25, 50, 75, 100), name = "Detection prevalence (%)")
  p_rank <- ggplot(top_plot, aes(x = tpi_score, y = taxon)) +
    geom_segment(aes(x = 0, xend = tpi_score, yend = taxon), color = "#8ba6af", size = .45) +
    geom_point(aes(size = prevalence_percent), color = "#086b7d") +
    score_scale() + prevalence_scale() +
    labs(title = paste0("Top ", n_top, " taxa by TPI"),
      subtitle = paste0(rank_label, "-level classifications | ", mode_label),
      x = "Taxon Prioritization Index (TPI)", y = NULL) + tpi_theme +
    scale_y_discrete(limits = levels_top, labels = taxon_labels, drop = FALSE,
      expand = expansion(add = .6)) +
    theme(legend.position = "bottom", axis.text.y = element_text(size = 9),
      legend.title = element_text(size = 9), legend.text = element_text(size = 9))
  component_cols <- c("tpi_prevalence_contribution", "tpi_mean_contribution",
    "tpi_max_contribution", "tpi_total_contribution")
  component_labels <- sprintf(c("Prevalence (%.2f)", "Mean (%.2f)",
    "Maximum (%.2f)", "Total (%.2f)"), c(w_prevalence,w_mean,w_max,w_total))
  components_plot <- do.call(rbind, lapply(seq_along(component_cols), function(j) {
    data.frame(taxon = top_plot$taxon, component = component_labels[j],
      contribution = top_plot[[component_cols[j]]], stringsAsFactors = FALSE)
  }))
  components_plot$taxon <- factor(components_plot$taxon, levels = levels_top)
  components_plot$component <- factor(components_plot$component, levels = component_labels)
  p_components <- ggplot(components_plot, aes(x = contribution, y = taxon, fill = component)) +
    geom_col(width = .62, position = position_stack(reverse = TRUE)) +
    score_scale() + taxon_scale() +
    scale_fill_manual(values = setNames(c("#086b7d", "#448dbb", "#df9938", "#a6b4c0"), component_labels), drop = FALSE) +
    labs(title = "TPI components", subtitle = "Weighted contributions", x = "Weighted contribution", y = NULL, fill = NULL) +
    tpi_theme + theme(legend.position = "bottom", legend.text = element_text(size = 9)) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE))
  # The standalone component panel keeps the rank/taxon key.
  p_components_single <- p_components + scale_y_discrete(limits = levels_top, labels = taxon_labels,
    drop = FALSE, expand = expansion(add = .6)) + theme(axis.text.y = element_text(size = 9))
  p_components <- p_components + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  # All positive taxa appear in the log-scale landscape; zero means cannot be plotted on log10.
  landscape_data <- tpi_ranked[is.finite(tpi_ranked$mean_abundance) & tpi_ranked$mean_abundance > 0, , drop = FALSE]
  label_df <- landscape_data[landscape_data$tpi_rank <= n_top, , drop = FALSE]
  sample_counts <- unique(tpi_ranked$n_samples)
  sample_text <- if (length(sample_counts) == 1) as.character(sample_counts) else "available"
  landscape_note <- paste0("Numbers identify the TPI ranks in A and B; ", mode_label, ".")
  if (nrow(landscape_data) < nrow(tpi_ranked)) landscape_note <- paste0(landscape_note,
    " Zero/non-finite means omitted from log axis: ", nrow(tpi_ranked) - nrow(landscape_data), ".")
  p_landscape <- ggplot(landscape_data, aes(x = prevalence_percent, y = mean_abundance)) +
    geom_point(aes(size = max_abundance, color = tpi_score), alpha = .78) +
    scale_y_log10(labels = function(x) format(x, scientific = FALSE, trim = TRUE)) +
    scale_x_continuous(breaks = seq(0, 100, by = 25), expand = expansion(mult = c(.03,.10))) +
    coord_cartesian(xlim = c(0, 100), clip = "off") +
    scale_color_gradient(low = "#c9dce5", high = "#075566", limits = c(0,1), name = "TPI") +
    scale_size_area(max_size = 6, name = max_label) +
    labs(title = paste0("Prevalence and abundance across ", nrow(tpi_ranked), " taxa"),
      subtitle = mode_label, x = paste0("Prevalence across ", sample_text, " libraries (%)"),
      y = mean_label, caption = landscape_note) + tpi_theme +
    theme(legend.position = "right", plot.caption = element_text(hjust = 0, color = "#536777", size = 9))
  if (has_ggrepel) {
    p_landscape <- p_landscape + ggrepel::geom_text_repel(data = label_df,
      aes(label = tpi_rank), seed = 7955, size = 3.2, max.overlaps = Inf,
      box.padding = .35, point.padding = .25, min.segment.length = 0, segment.color = "#8ba6af",
      color = "#203642", show.legend = FALSE)
  } else {
    # No new required package; labels may overlap if ggrepel is unavailable.
    warning("ggrepel unavailable: TPI landscape rank labels may overlap.")
    p_landscape <- p_landscape + geom_text(data = label_df, aes(label = tpi_rank),
      size = 3, nudge_x = 1.4, vjust = -.6, color = "#203642", show.legend = FALSE)
  }
  # grid and gtable are already supplied by R and ggplot2. Align the two upper panels.
  overview <- function() {
    a <- ggplotGrob(p_rank + labs(tag = "A") + theme(plot.tag = element_text(face = "bold")))
    b <- ggplotGrob(p_components + labs(tag = "B") + theme(plot.tag = element_text(face = "bold")))
    if (length(a$heights) == length(b$heights)) {
      h <- grid::unit.pmax(a$heights, b$heights)
      a$heights <- h; b$heights <- h
    }
    cplot <- ggplotGrob(p_landscape + labs(tag = "C") + theme(plot.tag = element_text(face = "bold")))
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 2,
      widths = grid::unit(c(0.60, 0.40), "null"),
      heights = grid::unit(c(max(4.8, 0.25*label_lines + 1.8), 4.5), "null"))))
    draw_at <- function(g, row, col) {
      grid::pushViewport(grid::viewport(layout.pos.row = row, layout.pos.col = col))
      grid::grid.draw(g); grid::popViewport()
    }
    draw_at(a, 1, 1); draw_at(b, 1, 2); draw_at(cplot, 2, 1:2)
    grid::popViewport()
  }
  max_label_chars <- max(nchar(as.character(top_plot$taxon)))
  label_lines <- sum(vapply(as.character(top_plot$taxon), function(z) max(1L, length(strwrap(z, 52))), integer(1)))
  top_height <- max(plot_height, 2.2 + .27 * label_lines)
  top_width <- max(plot_width, 8.5, 6 + min(max_label_chars, 52) * .075)
  overview_width <- max(12, top_width + 3.3)
  overview_height <- max(10, top_height + 4.5)
  save_plot <- function(plot, stem, width, height) {
    ggsave(file.path(output_dir, paste0(stem, ".png")), plot = plot, width = width, height = height,
      dpi = plot_dpi, bg = "white", limitsize = FALSE)
    if (save_pdf) ggsave(file.path(output_dir, paste0(stem, ".pdf")), plot = plot,
      width = width, height = height, bg = "white", limitsize = FALSE)
  }
  save_plot(p_rank, paste0(top_prefix, "_TPI_lollipop"), top_width, top_height)
  save_plot(p_components_single, paste0(prefix, "_components"), top_width, top_height)
  save_plot(p_landscape, paste0(prefix, "_landscape"), max(plot_width, 10), max(plot_height, 6))
  grDevices::png(file.path(output_dir, paste0(prefix, "_overview.png")), width = overview_width,
    height = overview_height, units = "in", res = plot_dpi, bg = "white")
  tryCatch(overview(), finally = grDevices::dev.off())
  if (save_pdf) {
    grDevices::pdf(file.path(output_dir, paste0(prefix, "_overview.pdf")), width = overview_width,
      height = overview_height, useDingbats = FALSE)
    tryCatch(overview(), finally = grDevices::dev.off())
  }
  # Legacy paths remain available to downstream tools and documentation.
  for (ext in c("png", if (save_pdf) "pdf")) {
    if (!file.copy(file.path(output_dir, paste0(top_prefix, "_TPI_lollipop.", ext)),
      file.path(output_dir, paste0(top_prefix, "_importance_lollipop.", ext)), overwrite = TRUE)) stop("Cannot write legacy TPI lollipop alias")
    if (!file.copy(file.path(output_dir, paste0(prefix, "_landscape.", ext)),
      file.path(output_dir, paste0("detected_microbiome_",rank_label,"_prevalence_abundance_importance_landscape.",ext)), overwrite = TRUE)) stop("Cannot write legacy TPI landscape alias")
  }
  cat("[OK] TPI lollipop, components, landscape and overview saved.\n")
}

cat("[OK] Main ranked table:\n")
cat(file.path(output_dir, paste0("detected_microbiome_", rank_label, "_ranked_by_importance.tsv")), "\n")

cat("[OK] Top 20 by importance:\n")
cat(file.path(output_dir, paste0("top20_detected_microbiome_", rank_label, "_by_importance.tsv")), "\n")

cat("[OK] README:\n")
cat(readme_file, "\n")

cat("Taxa ranked:", nrow(ranked_importance), "\n")
cat("Core taxa at >=", core_threshold, "% prevalence:", nrow(core), "\n")
cat("Strict core taxa at 100% prevalence:", nrow(strict_core), "\n")
cat("============================================================\n")
