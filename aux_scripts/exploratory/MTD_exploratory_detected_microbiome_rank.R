#!/usr/bin/env Rscript

# ============================================================
# MTD exploratory: detected microbiome ranked tables
#
# Input:
#   prevalence_vs_abundance_<rank>_relative.tsv
#
# Main score:
#   importance_score =
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
      Number of top taxa to show in the importance score figure.
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

if (!is.finite(weight_sum) || weight_sum <= 0) {
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
  "Main ranking formula:\n",
  "importance_score = ",
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
  "Generated files:\n",
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

if (plot_top > 0) {

  top_plot <- head(ranked_importance, plot_top)

  top_plot$taxon <- factor(
    top_plot$taxon,
    levels = rev(top_plot$taxon)
  )

  # ----------------------------------------------------------
  # Figure 1: horizontal barplot / lollipop-style ranking
  # ----------------------------------------------------------

  p_rank <- ggplot(
    top_plot,
    aes(
      x = importance_score,
      y = taxon
    )
  ) +
    geom_segment(
      aes(
        x = 0,
        xend = importance_score,
        y = taxon,
        yend = taxon
      ),
      linewidth = 0.7,
      alpha = 0.7
    ) +
    geom_point(
      aes(size = prevalence_percent),
      alpha = 0.9
    ) +
    labs(
      title = paste0("Top ", plot_top, " detected ", rank_label, " by exploratory importance score"),
      subtitle = paste0(
        "Score = ",
        round(w_prevalence, 2), "×prevalence + ",
        round(w_mean, 2), "×mean abundance + ",
        round(w_max, 2), "×max abundance + ",
        round(w_total, 2), "×total abundance"
      ),
      x = "Exploratory importance score",
      y = rank_label,
      size = "Prevalence (%)"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.y = element_text(size = 9),
      panel.grid.minor = element_blank()
    )

  rank_png <- file.path(
    output_dir,
    paste0("top", plot_top, "_detected_microbiome_", rank_label, "_importance_lollipop.png")
  )

  ggsave(
    filename = rank_png,
    plot = p_rank,
    width = plot_width,
    height = plot_height,
    dpi = plot_dpi
  )

  cat("[OK] Importance ranking figure saved:\n")
  cat(rank_png, "\n")

  if (save_pdf) {
    rank_pdf <- file.path(
      output_dir,
      paste0("top", plot_top, "_detected_microbiome_", rank_label, "_importance_lollipop.pdf")
    )

    ggsave(
      filename = rank_pdf,
      plot = p_rank,
      width = plot_width,
      height = plot_height
    )

    cat("[OK] Importance ranking PDF saved:\n")
    cat(rank_pdf, "\n")
  }

  # ----------------------------------------------------------
  # Figure 2: prevalence-abundance landscape colored by score
  # ----------------------------------------------------------

  label_df <- top_plot

  p_landscape <- ggplot(
    ranked_importance,
    aes(
      x = prevalence_percent,
      y = mean_abundance,
      size = max_abundance,
      color = importance_score
    )
  ) +
    geom_point(alpha = 0.75) +
    scale_y_continuous(trans = "log10") +
    scale_x_continuous(
      limits = c(0, 100),
      breaks = seq(0, 100, by = 25)
    ) +
    labs(
      title = paste0("Detected ", rank_label, " prevalence-abundance landscape"),
      subtitle = paste0("Color indicates exploratory importance score; labels show top ", plot_top, " taxa"),
      x = "Prevalence across samples (%)",
      y = "Mean relative abundance (%)",
      size = "Max abundance",
      color = "Importance score"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    )

  if (has_ggrepel) {
    p_landscape <- p_landscape +
      ggrepel::geom_text_repel(
        data = label_df,
        aes(label = taxon),
        size = 3,
        max.overlaps = Inf,
        box.padding = 0.4,
        point.padding = 0.25,
        min.segment.length = 0
      )
  } else {
    p_landscape <- p_landscape +
      geom_text(
        data = label_df,
        aes(label = taxon),
        size = 3,
        check_overlap = TRUE,
        vjust = -0.8
      )
  }

  landscape_png <- file.path(
    output_dir,
    paste0("detected_microbiome_", rank_label, "_prevalence_abundance_importance_landscape.png")
  )

  ggsave(
    filename = landscape_png,
    plot = p_landscape,
    width = plot_width,
    height = plot_height,
    dpi = plot_dpi
  )

  cat("[OK] Prevalence-abundance importance landscape saved:\n")
  cat(landscape_png, "\n")

  if (save_pdf) {
    landscape_pdf <- file.path(
      output_dir,
      paste0("detected_microbiome_", rank_label, "_prevalence_abundance_importance_landscape.pdf")
    )

    ggsave(
      filename = landscape_pdf,
      plot = p_landscape,
      width = plot_width,
      height = plot_height
    )

    cat("[OK] Prevalence-abundance importance landscape PDF saved:\n")
    cat(landscape_pdf, "\n")
  }
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
