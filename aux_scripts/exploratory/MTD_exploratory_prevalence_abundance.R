#!/usr/bin/env Rscript

# ============================================================
# MTD exploratory plot: Prevalence vs Abundance
#
# Purpose:
#   Generate a prevalence x abundance scatter plot from a
#   Bracken/Kraken-like abundance table.
#
# Handles:
#   1. Wide Bracken combined tables
#   2. Wide generic abundance matrices
#   3. Long-format tables with sample/taxon/abundance columns
#
# Important Bracken fix:
#   Combined Bracken outputs often contain two numeric columns per sample:
#     sample_num
#     sample_frac
#
#   For prevalence calculations, this script keeps only one column per
#   biological sample. By default it keeps *_num columns because relative
#   abundance can be recalculated from estimated abundances.
#
# Outputs:
#   prevalence_vs_abundance_<rank>_<mode>.png
#   prevalence_vs_abundance_<rank>_<mode>.pdf
#   prevalence_vs_abundance_<rank>_<mode>.tsv
#
# ============================================================

# -----------------------------
# Basic argument parser
# -----------------------------

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
MTD exploratory: Prevalence vs Abundance
============================================================

Usage:

  Rscript MTD_exploratory_prevalence_abundance.R \\
    --input bracken_species_all \\
    --output MTD_res/exploratory/taxonomy/prevalence_abundance/species_relative \\
    --rank species \\
    --mode relative

Required:
  --input FILE
      Input abundance table.

Optional:
  --output DIR
      Output directory.
      Default: current directory.

  --rank STRING
      Taxonomic rank label used only in output filenames/titles.
      Example: species, genus, family, phylum.
      Default: taxa

  --mode relative|absolute
      relative = abundance as percentage per sample.
      absolute = raw abundance/count values.
      Default: relative

  --taxon_col COLUMN
      Name of the taxon column.
      If not provided, the script tries to detect it automatically.

  --sample_col COLUMN
      For long-format input: sample column name.

  --abundance_col COLUMN
      For long-format input: abundance/count column name.

  --presence_threshold NUMERIC
      Minimum abundance to consider a taxon present in a sample.
      Default: 0

  --min_prevalence INTEGER
      Keep only taxa detected in at least this number of samples.
      Default: 1

  --top_labels INTEGER
      Number of taxa to label.
      Labels are chosen by highest mean abundance.
      Default: 15

  --width NUMERIC
      Plot width in inches.
      Default: 8

  --height NUMERIC
      Plot height in inches.
      Default: 6

  --dpi INTEGER
      PNG resolution.
      Default: 300

  --title STRING
      Custom plot title.
      Default: automatic title.

  --no_pdf
      Do not save PDF.

  --help
      Show this help message.

============================================================
\n")
}

if (has_flag("--help") || length(args) == 0) {
  print_help()
  quit(save = "no", status = 0)
}

# -----------------------------
# Parameters
# -----------------------------

input_file <- get_arg("--input")
output_dir <- get_arg("--output", ".")
rank_label <- get_arg("--rank", "taxa")
mode <- get_arg("--mode", "relative")
taxon_col_user <- get_arg("--taxon_col", NULL)
sample_col_user <- get_arg("--sample_col", NULL)
abundance_col_user <- get_arg("--abundance_col", NULL)
presence_threshold_user <- get_arg("--presence_threshold", NULL)
min_prevalence <- as.integer(get_arg("--min_prevalence", "1"))
top_labels <- as.integer(get_arg("--top_labels", "15"))
plot_width <- as.numeric(get_arg("--width", "8"))
plot_height <- as.numeric(get_arg("--height", "6"))
plot_dpi <- as.integer(get_arg("--dpi", "300"))
custom_title <- get_arg("--title", NULL)
save_pdf <- !has_flag("--no_pdf")

if (is.null(input_file)) {
  stop("[ERROR] Missing required argument: --input")
}

if (!file.exists(input_file)) {
  stop(paste0("[ERROR] Input file not found: ", input_file))
}

if (!mode %in% c("relative", "absolute")) {
  stop("[ERROR] --mode must be either 'relative' or 'absolute'")
}

if (is.null(presence_threshold_user)) {
  presence_threshold <- 0
} else {
  presence_threshold <- as.numeric(presence_threshold_user)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Package loading
# -----------------------------

required_packages <- c("ggplot2")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste0(
      "[ERROR] Required R package not installed: ", pkg, "\n",
      "Install it with:\n",
      "install.packages('", pkg, "')\n"
    ))
  }
}

suppressPackageStartupMessages({
  library(ggplot2)
})

has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)

# -----------------------------
# Helper functions
# -----------------------------

detect_sep <- function(file) {
  first_line <- readLines(file, n = 1, warn = FALSE)

  n_tab <- lengths(regmatches(first_line, gregexpr("\t", first_line)))
  n_comma <- lengths(regmatches(first_line, gregexpr(",", first_line)))
  n_semicolon <- lengths(regmatches(first_line, gregexpr(";", first_line)))

  if (n_tab >= n_comma && n_tab >= n_semicolon && n_tab > 0) {
    "\t"
  } else if (n_semicolon > n_comma) {
    ";"
  } else {
    ","
  }
}

clean_column_names <- function(x) {
  x <- gsub("^X$", "taxon", x)
  x <- gsub("\\s+", "_", x)
  x <- gsub("\\.+", "_", x)
  x
}

to_numeric_safe <- function(x) {
  x <- as.character(x)
  x <- gsub(",", ".", x)
  suppressWarnings(as.numeric(x))
}

find_first_existing <- function(candidates, cols) {
  candidates_lower <- tolower(candidates)
  cols_lower <- tolower(cols)
  hit <- match(candidates_lower, cols_lower)
  hit <- hit[!is.na(hit)]
  if (length(hit) == 0) return(NULL)
  cols[hit[1]]
}

is_numeric_like <- function(x) {
  y <- to_numeric_safe(x)
  mean(!is.na(y)) >= 0.8
}

make_unique_taxa <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == ""] <- "Unclassified_or_empty_name"
  make.unique(x, sep = "_dup")
}

safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9_\\-]+", "_", x)
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}

classify_prevalence <- function(prev_fraction) {
  ifelse(
    prev_fraction >= 0.75, "core_high_prevalence",
    ifelse(
      prev_fraction >= 0.50, "common",
      ifelse(prev_fraction >= 0.25, "intermediate", "rare_low_prevalence")
    )
  )
}

clean_bracken_sample_names <- function(sample_cols) {
  clean <- sample_cols

  clean <- gsub("^Report_", "", clean)
  clean <- gsub("_num$", "", clean, ignore.case = TRUE)
  clean <- gsub("_frac$", "", clean, ignore.case = TRUE)
  clean <- gsub("\\.num$", "", clean, ignore.case = TRUE)
  clean <- gsub("\\.frac$", "", clean, ignore.case = TRUE)

  # Remove common Bracken rank/file suffixes.
  clean <- gsub("_phylum.*$", "", clean, ignore.case = TRUE)
  clean <- gsub("_genus.*$", "", clean, ignore.case = TRUE)
  clean <- gsub("_species.*$", "", clean, ignore.case = TRUE)

  clean <- gsub("_bracken$", "", clean, ignore.case = TRUE)
  clean <- gsub("\\.bracken$", "", clean, ignore.case = TRUE)

  clean <- gsub("_+$", "", clean)
  clean
}

select_one_bracken_column_per_sample <- function(numeric_sample_cols) {
  num_cols <- numeric_sample_cols[
    grepl("(_num$|\\.num$)", numeric_sample_cols, ignore.case = TRUE)
  ]

  frac_cols <- numeric_sample_cols[
    grepl("(_frac$|\\.frac$)", numeric_sample_cols, ignore.case = TRUE)
  ]

  if (length(num_cols) >= 2) {
    cat("[INFO] Bracken-style _num/_frac columns detected.\n")
    cat("[INFO] Keeping only _num columns for prevalence vs abundance.\n")
    sample_cols <- num_cols
  } else if (length(frac_cols) >= 2) {
    cat("[INFO] Bracken-style _frac columns detected, but no _num columns were found.\n")
    cat("[INFO] Keeping _frac columns for prevalence vs abundance.\n")
    sample_cols <- frac_cols
  } else {
    sample_cols <- numeric_sample_cols
  }

  names(sample_cols) <- clean_bracken_sample_names(sample_cols)

  sample_cols
}

# -----------------------------
# Read input
# -----------------------------

sep <- detect_sep(input_file)

df <- read.table(
  input_file,
  header = TRUE,
  sep = sep,
  quote = "\"",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

original_colnames <- colnames(df)
colnames(df) <- clean_column_names(colnames(df))

cat("============================================================\n")
cat("MTD exploratory: prevalence vs abundance\n")
cat("Input file  :", input_file, "\n")
cat("Separator   :", ifelse(sep == "\t", "TAB", sep), "\n")
cat("Rows        :", nrow(df), "\n")
cat("Columns     :", ncol(df), "\n")
cat("Rank        :", rank_label, "\n")
cat("Mode        :", mode, "\n")
cat("Output dir  :", output_dir, "\n")
cat("============================================================\n")

# -----------------------------
# Detect format
# -----------------------------

cols <- colnames(df)

taxon_candidates <- c(
  "taxon", "taxa", "name", "Name", "scientific_name",
  "species", "genus", "family", "phylum",
  "clade_name", "feature", "Feature", "ID", "id"
)

sample_candidates <- c(
  "sample", "Sample", "sample_id", "SampleID", "Sample_ID",
  "sample_name", "SampleName", "library", "Library"
)

abundance_candidates <- c(
  "abundance", "Abundance",
  "new_est_reads", "fraction_total_reads",
  "reads", "read_count", "count", "counts",
  "relative_abundance", "rel_abundance",
  "value", "Value"
)

known_metadata_cols <- c(
  "taxonomy_id", "taxid", "tax_id", "NCBI_tax_id",
  "taxonomy_lvl", "rank", "level",
  "kraken_assigned_reads", "added_reads",
  "new_est_reads", "fraction_total_reads",
  "score", "p_value", "pvalue", "padj", "qvalue", "q_value",
  "log2FoldChange", "logFC", "baseMean", "lfcSE", "stat"
)

if (!is.null(taxon_col_user)) {
  taxon_col <- taxon_col_user
} else {
  taxon_col <- find_first_existing(taxon_candidates, cols)
}

if (!is.null(sample_col_user)) {
  sample_col <- sample_col_user
} else {
  sample_col <- find_first_existing(sample_candidates, cols)
}

if (!is.null(abundance_col_user)) {
  abundance_col <- abundance_col_user
} else {
  abundance_col <- find_first_existing(abundance_candidates, cols)
}

if (!is.null(taxon_col) && !taxon_col %in% cols) {
  stop(paste0("[ERROR] Taxon column not found: ", taxon_col))
}

if (!is.null(sample_col) && !sample_col %in% cols) {
  stop(paste0("[ERROR] Sample column not found: ", sample_col))
}

if (!is.null(abundance_col) && !abundance_col %in% cols) {
  stop(paste0("[ERROR] Abundance column not found: ", abundance_col))
}

# -----------------------------
# Convert input to abundance matrix
# taxa as rows, samples as columns
# -----------------------------

abundance_matrix <- NULL

long_format_detected <- !is.null(taxon_col) &&
  !is.null(sample_col) &&
  !is.null(abundance_col) &&
  taxon_col != sample_col

if (long_format_detected) {

  cat("[INFO] Detected long-format table.\n")
  cat("[INFO] Taxon column     :", taxon_col, "\n")
  cat("[INFO] Sample column    :", sample_col, "\n")
  cat("[INFO] Abundance column :", abundance_col, "\n")

  long_df <- df[, c(taxon_col, sample_col, abundance_col)]
  colnames(long_df) <- c("taxon", "sample", "abundance")

  long_df$taxon <- as.character(long_df$taxon)
  long_df$sample <- as.character(long_df$sample)
  long_df$abundance <- to_numeric_safe(long_df$abundance)

  long_df <- long_df[!is.na(long_df$taxon) & !is.na(long_df$sample), ]
  long_df$abundance[is.na(long_df$abundance)] <- 0

  wide <- aggregate(
    abundance ~ taxon + sample,
    data = long_df,
    FUN = sum
  )

  abundance_matrix <- xtabs(abundance ~ taxon + sample, data = wide)
  abundance_matrix <- as.matrix(abundance_matrix)

} else {

  cat("[INFO] Detected wide-format table or single Bracken-like table.\n")

  if (is.null(taxon_col)) {
    non_numeric_cols <- cols[!sapply(df, is_numeric_like)]

    if (length(non_numeric_cols) > 0) {
      taxon_col <- non_numeric_cols[1]
    } else {
      taxon_col <- cols[1]
    }
  }

  cat("[INFO] Taxon column:", taxon_col, "\n")

  possible_sample_cols <- setdiff(cols, c(taxon_col, known_metadata_cols))

  numeric_sample_cols <- possible_sample_cols[
    sapply(df[, possible_sample_cols, drop = FALSE], is_numeric_like)
  ]

  if (length(numeric_sample_cols) == 0) {
    stop("[ERROR] Could not detect numeric sample columns.")
  }

  sample_cols <- select_one_bracken_column_per_sample(numeric_sample_cols)

  # Special case: individual Bracken file with one abundance column
  if (length(sample_cols) == 1) {
    sample_name <- tools::file_path_sans_ext(basename(input_file))
    cat("[WARNING] Only one numeric abundance column detected.\n")
    cat("[WARNING] Prevalence across samples is not very informative with one sample.\n")
    cat("[INFO] Using sample name:", sample_name, "\n")
  }

  cat("[INFO] Sample/abundance columns detected:", length(sample_cols), "\n")
  cat(paste0("       raw:   ", paste(sample_cols, collapse = ", "), "\n"))
  cat(paste0("       clean: ", paste(names(sample_cols), collapse = ", "), "\n"))

  taxa <- make_unique_taxa(df[[taxon_col]])

  abundance_matrix <- as.matrix(
    data.frame(
      lapply(df[, sample_cols, drop = FALSE], to_numeric_safe),
      check.names = FALSE
    )
  )

  rownames(abundance_matrix) <- taxa

  if (length(sample_cols) == 1) {
    colnames(abundance_matrix) <- sample_name
  } else if (!is.null(names(sample_cols)) &&
             all(!is.na(names(sample_cols))) &&
             all(names(sample_cols) != "")) {
    colnames(abundance_matrix) <- names(sample_cols)
  } else {
    colnames(abundance_matrix) <- sample_cols
  }
}

abundance_matrix[is.na(abundance_matrix)] <- 0
abundance_matrix[abundance_matrix < 0] <- 0

# Remove empty taxa
abundance_matrix <- abundance_matrix[rowSums(abundance_matrix) > 0, , drop = FALSE]

if (nrow(abundance_matrix) == 0) {
  stop("[ERROR] No taxa with abundance > 0 after filtering.")
}

# -----------------------------
# Relative or absolute conversion
# -----------------------------

plot_matrix <- abundance_matrix

if (mode == "relative") {

  max_val <- max(plot_matrix, na.rm = TRUE)
  col_sums <- colSums(plot_matrix, na.rm = TRUE)

  if (max_val <= 1) {
    cat("[INFO] Input appears to be fractional relative abundance. Multiplying by 100.\n")
    plot_matrix <- plot_matrix * 100
  } else if (max_val <= 100 && median(col_sums, na.rm = TRUE) <= 101) {
    cat("[INFO] Input appears to already be percentage relative abundance.\n")
  } else {
    cat("[INFO] Input appears to be counts. Converting each sample to percentage.\n")
    col_sums[col_sums == 0] <- NA
    plot_matrix <- sweep(plot_matrix, 2, col_sums, "/") * 100
    plot_matrix[is.na(plot_matrix)] <- 0
  }

} else {
  cat("[INFO] Absolute mode selected. Using raw abundance values.\n")
}

# -----------------------------
# Calculate prevalence and abundance
# -----------------------------

present_absent <- plot_matrix > presence_threshold

prevalence_n <- rowSums(present_absent, na.rm = TRUE)
n_samples <- ncol(plot_matrix)
prevalence_fraction <- prevalence_n / n_samples
prevalence_percent <- prevalence_fraction * 100

mean_abundance <- rowMeans(plot_matrix, na.rm = TRUE)
median_abundance <- apply(plot_matrix, 1, median, na.rm = TRUE)
max_abundance <- apply(plot_matrix, 1, max, na.rm = TRUE)
total_abundance <- rowSums(plot_matrix, na.rm = TRUE)

result <- data.frame(
  taxon = rownames(plot_matrix),
  rank = rank_label,
  n_samples = n_samples,
  prevalence_n = prevalence_n,
  prevalence_percent = prevalence_percent,
  mean_abundance = mean_abundance,
  median_abundance = median_abundance,
  max_abundance = max_abundance,
  total_abundance = total_abundance,
  prevalence_class = classify_prevalence(prevalence_fraction),
  stringsAsFactors = FALSE
)

result <- result[result$prevalence_n >= min_prevalence, , drop = FALSE]

if (nrow(result) == 0) {
  stop("[ERROR] No taxa left after min_prevalence filtering.")
}

result <- result[order(-result$mean_abundance, -result$prevalence_n), ]

label_taxa <- head(result$taxon, top_labels)
result$label <- ifelse(result$taxon %in% label_taxa, result$taxon, "")

# -----------------------------
# Save table
# -----------------------------

safe_rank <- safe_name(rank_label)
safe_mode <- safe_name(mode)

out_prefix <- file.path(
  output_dir,
  paste0("prevalence_vs_abundance_", safe_rank, "_", safe_mode)
)

table_file <- paste0(out_prefix, ".tsv")

write.table(
  result,
  file = table_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("[OK] Table saved:", table_file, "\n")

# -----------------------------
# Plot
# -----------------------------

y_label <- ifelse(
  mode == "relative",
  "Mean relative abundance (%)",
  "Mean absolute abundance"
)

x_label <- "Prevalence across samples (%)"

if (is.null(custom_title)) {
  plot_title <- paste0(
    "Prevalence vs abundance at ", rank_label, " level"
  )
} else {
  plot_title <- custom_title
}

subtitle_text <- paste0(
  "Presence threshold: > ", presence_threshold,
  " | Samples: ", n_samples,
  " | Taxa shown: ", nrow(result)
)

p <- ggplot(
  result,
  aes(
    x = prevalence_percent,
    y = mean_abundance,
    size = max_abundance,
    color = prevalence_class
  )
) +
  geom_point(alpha = 0.75) +
  scale_y_continuous(trans = "log10") +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 25)) +
  labs(
    title = plot_title,
    subtitle = subtitle_text,
    x = x_label,
    y = y_label,
    color = "Prevalence class",
    size = "Max abundance"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

if (top_labels > 0) {
  if (has_ggrepel) {
    p <- p +
      ggrepel::geom_text_repel(
        aes(label = label),
        size = 3,
        max.overlaps = Inf,
        box.padding = 0.4,
        point.padding = 0.25,
        min.segment.length = 0
      )
  } else {
    cat("[WARNING] Package 'ggrepel' not installed. Using geom_text with check_overlap.\n")
    cat("[WARNING] For better labels, install it with: install.packages('ggrepel')\n")

    p <- p +
      geom_text(
        aes(label = label),
        size = 3,
        check_overlap = TRUE,
        vjust = -0.8
      )
  }
}

png_file <- paste0(out_prefix, ".png")

ggsave(
  filename = png_file,
  plot = p,
  width = plot_width,
  height = plot_height,
  dpi = plot_dpi
)

cat("[OK] PNG saved:", png_file, "\n")

if (save_pdf) {
  pdf_file <- paste0(out_prefix, ".pdf")

  ggsave(
    filename = pdf_file,
    plot = p,
    width = plot_width,
    height = plot_height
  )

  cat("[OK] PDF saved:", pdf_file, "\n")
}

# -----------------------------
# Final summary
# -----------------------------

cat("============================================================\n")
cat("[DONE] Prevalence vs abundance plot completed.\n")
cat("Taxa analyzed :", nrow(result), "\n")
cat("Samples       :", n_samples, "\n")
cat("Mode          :", mode, "\n")
cat("Rank          :", rank_label, "\n")
cat("============================================================\n")
