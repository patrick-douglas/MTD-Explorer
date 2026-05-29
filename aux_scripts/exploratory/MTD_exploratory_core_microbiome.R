#!/usr/bin/env Rscript

# ============================================================
# MTD exploratory plot: Core microbiome
#
# Purpose:
#   Estimate and plot the number of taxa retained as "core"
#   across different prevalence thresholds.
#
# Input:
#   Wide Bracken/Kraken-like table:
#
#   name/sample        sample1   sample2   sample3
#   Escherichia coli   10        0         5
#   Fusarium poae      0         100       3
#
#   Also tries to handle long-format tables with:
#   sample, taxon/name, abundance/count/new_est_reads
#
# Outputs:
#   1. core_microbiome_<rank>_<mode>.png
#   2. core_microbiome_<rank>_<mode>.pdf
#   3. core_microbiome_<rank>_<mode>_summary.tsv
#   4. core_microbiome_<rank>_<mode>_all_taxa.tsv
#   5. core_microbiome_<rank>_<mode>_coreXX.tsv
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
  return(args[hit + 1])
}

has_flag <- function(flag) {
  flag %in% args
}

print_help <- function() {
  cat("
============================================================
MTD exploratory: Core microbiome
============================================================

Usage:

  Rscript MTD_exploratory_core_microbiome.R \\
    --input bracken_genus_all \\
    --output MTD_res/exploratory/taxonomy/core_microbiome/genus_relative \\
    --rank genus \\
    --mode relative

Required:
  --input FILE
      Input abundance table.

Optional:
  --output DIR
      Output directory.
      Default: current directory.

  --rank STRING
      Taxonomic rank label used in filenames and plot titles.
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
      For relative mode, this is percentage.
      For absolute mode, this is count/read abundance.
      Default: 0

  --thresholds STRING
      Comma-separated prevalence thresholds in percent.
      Default: 25,50,75,90,100

  --core_threshold NUMERIC
      Main threshold used to export a dedicated core table.
      Default: 75

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

Examples:

  Rscript MTD_exploratory_core_microbiome.R \\
    --input bracken_genus_all \\
    --output MTD_res/exploratory/taxonomy/core_microbiome/genus_relative \\
    --rank genus \\
    --mode relative

  Rscript MTD_exploratory_core_microbiome.R \\
    --input bracken_species_all \\
    --output MTD_res/exploratory/taxonomy/core_microbiome/species_absolute \\
    --rank species \\
    --mode absolute \\
    --presence_threshold 10

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
presence_threshold <- as.numeric(get_arg("--presence_threshold", "0"))
thresholds_string <- get_arg("--thresholds", "25,50,75,90,100")
core_threshold <- as.numeric(get_arg("--core_threshold", "75"))
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

thresholds <- as.numeric(strsplit(thresholds_string, ",")[[1]])
thresholds <- thresholds[!is.na(thresholds)]
thresholds <- sort(unique(thresholds))

if (length(thresholds) == 0) {
  stop("[ERROR] No valid thresholds were provided.")
}

if (any(thresholds <= 0 | thresholds > 100)) {
  stop("[ERROR] --thresholds must contain values > 0 and <= 100.")
}

if (is.na(core_threshold) || core_threshold <= 0 || core_threshold > 100) {
  stop("[ERROR] --core_threshold must be > 0 and <= 100.")
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


# -----------------------------
# Helper functions
# -----------------------------

detect_sep <- function(file) {
  first_line <- readLines(file, n = 1, warn = FALSE)

  n_tab <- lengths(regmatches(first_line, gregexpr("\t", first_line)))
  n_comma <- lengths(regmatches(first_line, gregexpr(",", first_line)))
  n_semicolon <- lengths(regmatches(first_line, gregexpr(";", first_line)))

  if (n_tab >= n_comma && n_tab >= n_semicolon && n_tab > 0) {
    return("\t")
  } else if (n_semicolon > n_comma) {
    return(";")
  } else {
    return(",")
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
  x <- gsub("^_|_$", "", x)
  x
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

colnames(df) <- clean_column_names(colnames(df))

cat("============================================================\n")
cat("MTD exploratory: core microbiome\n")
cat("Input file        :", input_file, "\n")
cat("Separator         :", ifelse(sep == "\t", "TAB", sep), "\n")
cat("Rows              :", nrow(df), "\n")
cat("Columns           :", ncol(df), "\n")
cat("Rank              :", rank_label, "\n")
cat("Mode              :", mode, "\n")
cat("Presence threshold:", presence_threshold, "\n")
cat("Core thresholds   :", paste(thresholds, collapse = ", "), "\n")
cat("Main core threshold:", core_threshold, "\n")
cat("Output dir        :", output_dir, "\n")
cat("============================================================\n")


# -----------------------------
# Detect table format
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

  sample_cols <- numeric_sample_cols

  cat("[INFO] Sample/abundance columns detected:", length(sample_cols), "\n")
  cat(paste0("       ", paste(sample_cols, collapse = ", "), "\n"))

  taxa <- make_unique_taxa(df[[taxon_col]])

  abundance_matrix <- as.matrix(
    data.frame(
      lapply(df[, sample_cols, drop = FALSE], to_numeric_safe),
      check.names = FALSE
    )
  )

  rownames(abundance_matrix) <- taxa
  colnames(abundance_matrix) <- sample_cols
}

abundance_matrix[is.na(abundance_matrix)] <- 0
abundance_matrix[abundance_matrix < 0] <- 0

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
# Calculate prevalence
# -----------------------------

present_absent <- plot_matrix > presence_threshold

n_samples <- ncol(plot_matrix)
prevalence_n <- rowSums(present_absent, na.rm = TRUE)
prevalence_percent <- (prevalence_n / n_samples) * 100

mean_abundance <- rowMeans(plot_matrix, na.rm = TRUE)
median_abundance <- apply(plot_matrix, 1, median, na.rm = TRUE)
max_abundance <- apply(plot_matrix, 1, max, na.rm = TRUE)
total_abundance <- rowSums(plot_matrix, na.rm = TRUE)

all_taxa <- data.frame(
  taxon = rownames(plot_matrix),
  rank = rank_label,
  n_samples = n_samples,
  prevalence_n = prevalence_n,
  prevalence_percent = prevalence_percent,
  mean_abundance = mean_abundance,
  median_abundance = median_abundance,
  max_abundance = max_abundance,
  total_abundance = total_abundance,
  stringsAsFactors = FALSE
)

for (thr in thresholds) {
  colname <- paste0("core_", safe_name(thr), "pct")
  all_taxa[[colname]] <- all_taxa$prevalence_percent >= thr
}

all_taxa <- all_taxa[order(-all_taxa$prevalence_percent, -all_taxa$mean_abundance), ]


# -----------------------------
# Summary by threshold
# -----------------------------

summary_df <- data.frame(
  threshold_percent = thresholds,
  min_samples_required = ceiling((thresholds / 100) * n_samples),
  core_taxa_n = NA_integer_,
  mean_abundance_sum = NA_real_,
  max_abundance_sum = NA_real_,
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(summary_df))) {
  thr <- summary_df$threshold_percent[i]
  core_idx <- all_taxa$prevalence_percent >= thr

  summary_df$core_taxa_n[i] <- sum(core_idx)
  summary_df$mean_abundance_sum[i] <- sum(all_taxa$mean_abundance[core_idx], na.rm = TRUE)
  summary_df$max_abundance_sum[i] <- sum(all_taxa$max_abundance[core_idx], na.rm = TRUE)
}


# -----------------------------
# Save tables
# -----------------------------

safe_rank <- safe_name(rank_label)
safe_mode <- safe_name(mode)
safe_core <- safe_name(core_threshold)

out_prefix <- file.path(
  output_dir,
  paste0("core_microbiome_", safe_rank, "_", safe_mode)
)

summary_file <- paste0(out_prefix, "_summary.tsv")
all_taxa_file <- paste0(out_prefix, "_all_taxa.tsv")
main_core_file <- paste0(out_prefix, "_core", safe_core, ".tsv")

write.table(
  summary_df,
  file = summary_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  all_taxa,
  file = all_taxa_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

main_core <- all_taxa[all_taxa$prevalence_percent >= core_threshold, , drop = FALSE]

write.table(
  main_core,
  file = main_core_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("[OK] Summary table saved:", summary_file, "\n")
cat("[OK] All taxa table saved:", all_taxa_file, "\n")
cat("[OK] Main core taxa table saved:", main_core_file, "\n")


# -----------------------------
# Plot
# -----------------------------

if (is.null(custom_title)) {
  plot_title <- paste0("Core microbiome at ", rank_label, " level")
} else {
  plot_title <- custom_title
}

subtitle_text <- paste0(
  "Presence threshold: > ", presence_threshold,
  " | Samples: ", n_samples,
  " | Total taxa: ", nrow(all_taxa),
  " | Main core threshold: ", core_threshold, "%"
)

p <- ggplot(
  summary_df,
  aes(
    x = threshold_percent,
    y = core_taxa_n
  )
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 3) +
  geom_text(
    aes(label = core_taxa_n),
    vjust = -0.8,
    size = 3.8
  ) +
  scale_x_continuous(
    limits = c(0, 100),
    breaks = sort(unique(c(0, thresholds, 100)))
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.15))
  ) +
  labs(
    title = plot_title,
    subtitle = subtitle_text,
    x = "Minimum prevalence required across samples (%)",
    y = paste0("Number of core ", rank_label)
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

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
cat("[DONE] Core microbiome plot completed.\n")
cat("Taxa analyzed       :", nrow(all_taxa), "\n")
cat("Samples             :", n_samples, "\n")
cat("Mode                :", mode, "\n")
cat("Rank                :", rank_label, "\n")
cat("Presence threshold  :", presence_threshold, "\n")
cat("Main core threshold :", core_threshold, "%\n")
cat("Core taxa at main threshold:", nrow(main_core), "\n")
cat("============================================================\n")
