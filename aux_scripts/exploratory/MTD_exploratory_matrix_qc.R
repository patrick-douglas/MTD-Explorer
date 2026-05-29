#!/usr/bin/env Rscript

# ============================================================
# MTD exploratory plot: Matrix QC
#
# Works for:
#   1) Host featureCounts matrix
#   2) Bracken/taxonomic abundance matrix
#   3) Generic feature x sample matrix
#
# Generates:
#   - PCA
#   - Top variable features heatmap
#   - Sample correlation heatmap
#   - Detected features per sample
#
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
MTD exploratory: Matrix QC
============================================================

Usage:

  Rscript MTD_exploratory_matrix_qc.R \\
    --input MATRIX.tsv \\
    --samplesheet samplesheet.csv \\
    --output OUTDIR \\
    --label host_expression \\
    --matrix_type featurecounts \\
    --normalization logcpm

Required:
  --input FILE
      Input matrix.

Optional:
  --samplesheet FILE
      Samplesheet CSV. If provided, group information is used.

  --output DIR
      Output directory. Default: current directory.

  --label STRING
      Label used in filenames and plot titles.
      Examples: host_expression, microbiome_genus.
      Default: matrix_qc

  --matrix_type auto|featurecounts|bracken|generic
      auto = try to detect format.
      featurecounts = featureCounts host matrix.
      bracken = combined Bracken output.
      generic = first non-numeric column is feature ID.
      Default: auto

  --normalization auto|logcpm|relative|none
      auto = logcpm for featureCounts, relative for Bracken/generic.
      logcpm = log2(counts per million + 1).
      relative = convert each sample to percent and log10(percent + 1).
      none = log10(raw value + 1).
      Default: auto

  --feature_col COLUMN
      Feature/taxon/gene column name.
      If omitted, auto-detected.

  --sample_id_col COLUMN
      Sample ID column in samplesheet. Default: first column.

  --group_col COLUMN
      Group column in samplesheet. Default: second column.

  --top_variable INTEGER
      Number of top variable features for heatmap.
      Default: 50

  --pca_top_features INTEGER
      Number of top variable features used for PCA.
      Default: 5000

  --detected_threshold NUMERIC
      Raw abundance/count threshold for detected feature.
      Default: 0

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
      Do not save PDF.

  --help
      Show this help message.

Outputs:
  <label>_normalized_matrix.tsv
  <label>_sample_qc.tsv
  <label>_top_variable_features.tsv
  <label>_pca_scores.tsv
  <label>_sample_correlation.tsv

  <label>_pca.png/pdf
  <label>_top_variable_heatmap.png/pdf
  <label>_sample_correlation_heatmap.png/pdf
  <label>_detected_features_per_sample.png/pdf

============================================================
\n")
}

if (has_flag("--help") || length(args) == 0) {
  print_help()
  quit(save = "no", status = 0)
}

input_file <- get_arg("--input")
samplesheet_file <- get_arg("--samplesheet", NULL)
output_dir <- get_arg("--output", ".")
label <- get_arg("--label", "matrix_qc")
matrix_type <- get_arg("--matrix_type", "auto")
normalization <- get_arg("--normalization", "auto")
feature_col_user <- get_arg("--feature_col", NULL)
sample_id_col_user <- get_arg("--sample_id_col", NULL)
group_col_user <- get_arg("--group_col", NULL)
top_variable <- as.integer(get_arg("--top_variable", "50"))
pca_top_features <- as.integer(get_arg("--pca_top_features", "5000"))
detected_threshold <- as.numeric(get_arg("--detected_threshold", "0"))
plot_width <- as.numeric(get_arg("--width", "9"))
plot_height <- as.numeric(get_arg("--height", "7"))
plot_dpi <- as.integer(get_arg("--dpi", "300"))
save_pdf <- !has_flag("--no_pdf")

if (is.null(input_file)) stop("[ERROR] Missing --input")
if (!file.exists(input_file)) stop(paste0("[ERROR] Input file not found: ", input_file))

if (!matrix_type %in% c("auto", "featurecounts", "bracken", "generic")) {
  stop("[ERROR] --matrix_type must be auto, featurecounts, bracken or generic")
}

if (!normalization %in% c("auto", "logcpm", "relative", "none")) {
  stop("[ERROR] --normalization must be auto, logcpm, relative or none")
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

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

to_numeric_safe <- function(x) {
  x <- as.character(x)
  x <- gsub(",", ".", x)
  suppressWarnings(as.numeric(x))
}

is_numeric_like <- function(x) {
  y <- to_numeric_safe(x)
  mean(!is.na(y)) >= 0.8
}

find_first_existing <- function(candidates, cols) {
  hit <- match(tolower(candidates), tolower(cols))
  hit <- hit[!is.na(hit)]
  if (length(hit) == 0) return(NULL)
  cols[hit[1]]
}

safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9_\\-]+", "_", x)
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}

make_unique_features <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == ""] <- "Unclassified_or_empty_name"
  make.unique(x, sep = "_dup")
}

row_zscore <- function(mat) {
  out <- t(scale(t(mat)))
  out[is.na(out)] <- 0
  out
}

read_samplesheet <- function(samplesheet_file, sample_id_col_user, group_col_user) {
  if (is.null(samplesheet_file) || !file.exists(samplesheet_file)) return(NULL)

  sep <- detect_sep(samplesheet_file)

  ss <- read.table(
    samplesheet_file,
    header = TRUE,
    sep = sep,
    quote = "\"",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  if (ncol(ss) < 2) return(NULL)

  sample_col <- if (!is.null(sample_id_col_user)) sample_id_col_user else colnames(ss)[1]
  group_col <- if (!is.null(group_col_user)) group_col_user else colnames(ss)[2]

  if (!sample_col %in% colnames(ss)) sample_col <- colnames(ss)[1]
  if (!group_col %in% colnames(ss)) group_col <- colnames(ss)[2]

  out <- data.frame(
    sample = as.character(ss[[sample_col]]),
    group = as.character(ss[[group_col]]),
    stringsAsFactors = FALSE
  )

  out <- out[!is.na(out$sample) & out$sample != "", , drop = FALSE]
  out
}

matrix_to_long <- function(mat, value_name = "value") {
  df <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
  colnames(df) <- c("feature", "sample", value_name)
  df
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

cat("============================================================\n")
cat("MTD exploratory: matrix QC\n")
cat("Input file     :", input_file, "\n")
cat("Separator      :", ifelse(sep == "\t", "TAB", sep), "\n")
cat("Rows           :", nrow(df), "\n")
cat("Columns        :", ncol(df), "\n")
cat("Label          :", label, "\n")
cat("Matrix type    :", matrix_type, "\n")
cat("Normalization  :", normalization, "\n")
cat("Output dir     :", output_dir, "\n")
cat("============================================================\n")

cols <- colnames(df)

# -----------------------------
# Detect matrix type
# -----------------------------

if (matrix_type == "auto") {
  if (all(c("Geneid", "Chr", "Start", "End", "Strand", "Length") %in% cols)) {
    matrix_type <- "featurecounts"
  } else if (any(tolower(cols) %in% c("taxonomy_id", "taxonomy_lvl", "new_est_reads", "fraction_total_reads"))) {
    matrix_type <- "bracken"
  } else {
    matrix_type <- "generic"
  }

  cat("[INFO] Auto-detected matrix type:", matrix_type, "\n")
}

# -----------------------------
# Detect feature and sample columns
# -----------------------------

feature_col <- NULL
metadata_cols <- character(0)

if (matrix_type == "featurecounts") {
  feature_col <- if (!is.null(feature_col_user)) feature_col_user else "Geneid"

  metadata_cols <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")

  if (!feature_col %in% cols) {
    feature_col <- cols[1]
  }

} else if (matrix_type == "bracken") {

  taxon_candidates <- c(
    "name", "Name", "taxon", "taxa", "scientific_name",
    "species", "genus", "family", "phylum",
    "clade_name", "feature", "Feature", "ID", "id"
  )

  feature_col <- if (!is.null(feature_col_user)) {
    feature_col_user
  } else {
    find_first_existing(taxon_candidates, cols)
  }

  if (is.null(feature_col)) {
    non_numeric_cols <- cols[!sapply(df, is_numeric_like)]
    if (length(non_numeric_cols) > 0) {
      feature_col <- non_numeric_cols[1]
    } else {
      feature_col <- cols[1]
    }
  }

  metadata_cols <- c(
    feature_col,
    "taxonomy_id", "taxid", "tax_id", "NCBI_tax_id",
    "taxonomy_lvl", "rank", "level",
    "kraken_assigned_reads", "added_reads",
    "new_est_reads", "fraction_total_reads"
  )

} else {
  feature_col <- if (!is.null(feature_col_user)) feature_col_user else NULL

  if (is.null(feature_col)) {
    non_numeric_cols <- cols[!sapply(df, is_numeric_like)]
    if (length(non_numeric_cols) > 0) {
      feature_col <- non_numeric_cols[1]
    } else {
      feature_col <- cols[1]
    }
  }

  metadata_cols <- c(feature_col)
}

if (!feature_col %in% cols) {
  stop(paste0("[ERROR] Feature column not found: ", feature_col))
}

possible_sample_cols <- setdiff(cols, metadata_cols)

sample_cols <- possible_sample_cols[
  sapply(df[, possible_sample_cols, drop = FALSE], is_numeric_like)
]

# ------------------------------------------------------------
# Bracken combined tables often contain two numeric columns per sample:
#   sample_num  = estimated reads / abundance
#   sample_frac = fraction of total reads
#
# For matrix QC, keep only one measurement per biological sample.
# Default: keep *_num columns, because relative normalization can be
# computed downstream from counts/abundances.
# ------------------------------------------------------------

if (matrix_type == "bracken") {
  num_cols <- sample_cols[grepl("(_num$|\\.num$)", sample_cols, ignore.case = TRUE)]
  frac_cols <- sample_cols[grepl("(_frac$|\\.frac$)", sample_cols, ignore.case = TRUE)]

  if (length(num_cols) >= 2) {
    cat("[INFO] Bracken table detected with _num/_frac columns.\n")
    cat("[INFO] Keeping only _num columns for matrix QC.\n")
    sample_cols <- num_cols
  } else if (length(frac_cols) >= 2) {
    cat("[INFO] Bracken table detected with _frac columns only.\n")
    cat("[INFO] Keeping _frac columns for matrix QC.\n")
    sample_cols <- frac_cols
  }

  clean_sample_names <- sample_cols

  clean_sample_names <- gsub("^Report_", "", clean_sample_names)
  clean_sample_names <- gsub("_bracken$", "", clean_sample_names, ignore.case = TRUE)
  clean_sample_names <- gsub("_num$", "", clean_sample_names, ignore.case = TRUE)
  clean_sample_names <- gsub("_frac$", "", clean_sample_names, ignore.case = TRUE)
  clean_sample_names <- gsub("\\.num$", "", clean_sample_names, ignore.case = TRUE)
  clean_sample_names <- gsub("\\.frac$", "", clean_sample_names, ignore.case = TRUE)
  clean_sample_names <- gsub("\\.bracken$", "", clean_sample_names, ignore.case = TRUE)

  # Remove common rank/file suffixes left by Bracken combined outputs
  clean_sample_names <- gsub("_phylum.*$", "", clean_sample_names, ignore.case = TRUE)
  clean_sample_names <- gsub("_genus.*$", "", clean_sample_names, ignore.case = TRUE)
  clean_sample_names <- gsub("_species.*$", "", clean_sample_names, ignore.case = TRUE)

  names(sample_cols) <- clean_sample_names
}

if (length(sample_cols) < 2) {
  stop("[ERROR] At least 2 numeric sample columns are required.")
}

cat("[INFO] Feature column:", feature_col, "\n")
cat("[INFO] Sample columns detected:", length(sample_cols), "\n")
cat(paste0("       ", paste(sample_cols, collapse = ", "), "\n"))

features <- make_unique_features(df[[feature_col]])

raw_matrix <- as.matrix(
  data.frame(
    lapply(df[, sample_cols, drop = FALSE], to_numeric_safe),
    check.names = FALSE
  )
)

rownames(raw_matrix) <- features
if (matrix_type == "bracken" && !is.null(names(sample_cols)) && all(names(sample_cols) != "")) {
  colnames(raw_matrix) <- names(sample_cols)
} else {
  colnames(raw_matrix) <- sample_cols
}

raw_matrix[is.na(raw_matrix)] <- 0
raw_matrix[raw_matrix < 0] <- 0

raw_matrix <- raw_matrix[rowSums(raw_matrix) > 0, , drop = FALSE]

if (nrow(raw_matrix) < 2) stop("[ERROR] At least 2 non-zero features are required.")

# -----------------------------
# Normalization / transformation
# -----------------------------

if (normalization == "auto") {
  if (matrix_type == "featurecounts") {
    normalization <- "logcpm"
  } else {
    normalization <- "relative"
  }

  cat("[INFO] Auto-selected normalization:", normalization, "\n")
}

norm_matrix <- raw_matrix

if (normalization == "logcpm") {
  lib_sizes <- colSums(raw_matrix, na.rm = TRUE)
  lib_sizes[lib_sizes == 0] <- NA

  cpm <- sweep(raw_matrix, 2, lib_sizes, "/") * 1e6
  cpm[is.na(cpm)] <- 0

  norm_matrix <- log2(cpm + 1)

} else if (normalization == "relative") {
  col_sums <- colSums(raw_matrix, na.rm = TRUE)
  max_val <- max(raw_matrix, na.rm = TRUE)

  if (max_val <= 1) {
    rel <- raw_matrix * 100
  } else if (max_val <= 100 && median(col_sums, na.rm = TRUE) <= 101) {
    rel <- raw_matrix
  } else {
    col_sums[col_sums == 0] <- NA
    rel <- sweep(raw_matrix, 2, col_sums, "/") * 100
    rel[is.na(rel)] <- 0
  }

  norm_matrix <- log10(rel + 1)

} else if (normalization == "none") {
  norm_matrix <- log10(raw_matrix + 1)
}

norm_matrix[is.na(norm_matrix)] <- 0
norm_matrix[is.infinite(norm_matrix)] <- 0

# Remove zero-variance features for PCA/correlation
feature_var <- apply(norm_matrix, 1, var, na.rm = TRUE)
var_matrix <- norm_matrix[feature_var > 0, , drop = FALSE]

if (nrow(var_matrix) < 2) {
  stop("[ERROR] Fewer than 2 variable features after normalization.")
}

# -----------------------------
# Metadata / sample QC
# -----------------------------

ss <- read_samplesheet(samplesheet_file, sample_id_col_user, group_col_user)

sample_qc <- data.frame(
  sample = colnames(raw_matrix),
  group = "All_samples",
  total_raw_abundance = colSums(raw_matrix, na.rm = TRUE),
  detected_features = colSums(raw_matrix > detected_threshold, na.rm = TRUE),
  zero_fraction = colMeans(raw_matrix <= detected_threshold, na.rm = TRUE),
  stringsAsFactors = FALSE
)

if (!is.null(ss)) {
  sample_qc <- merge(sample_qc, ss, by = "sample", all.x = TRUE, suffixes = c("", "_samplesheet"))
  sample_qc$group <- sample_qc$group_samplesheet
  sample_qc$group[is.na(sample_qc$group) | sample_qc$group == ""] <- "Ungrouped"
  sample_qc$group_samplesheet <- NULL
}

sample_qc <- sample_qc[match(colnames(raw_matrix), sample_qc$sample), , drop = FALSE]

safe_label <- safe_name(label)

normalized_file <- file.path(output_dir, paste0(safe_label, "_normalized_matrix.tsv"))
sample_qc_file <- file.path(output_dir, paste0(safe_label, "_sample_qc.tsv"))

write.table(
  data.frame(feature = rownames(norm_matrix), norm_matrix, check.names = FALSE),
  normalized_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  sample_qc,
  sample_qc_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("[OK] Normalized matrix saved:", normalized_file, "\n")
cat("[OK] Sample QC table saved:", sample_qc_file, "\n")

# -----------------------------
# Top variable features
# -----------------------------

feature_var_df <- data.frame(
  feature = rownames(var_matrix),
  variance = apply(var_matrix, 1, var, na.rm = TRUE),
  mean_normalized = rowMeans(var_matrix, na.rm = TRUE),
  total_raw_abundance = rowSums(raw_matrix[rownames(var_matrix), , drop = FALSE], na.rm = TRUE),
  stringsAsFactors = FALSE
)

feature_var_df <- feature_var_df[order(-feature_var_df$variance), ]

top_variable <- min(top_variable, nrow(feature_var_df))
pca_top_features <- min(pca_top_features, nrow(feature_var_df))

top_features <- head(feature_var_df$feature, top_variable)
pca_features <- head(feature_var_df$feature, pca_top_features)

top_var_file <- file.path(output_dir, paste0(safe_label, "_top_variable_features.tsv"))

write.table(
  feature_var_df,
  top_var_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("[OK] Variable features table saved:", top_var_file, "\n")

# -----------------------------
# PCA
# -----------------------------

pca_matrix <- t(var_matrix[pca_features, , drop = FALSE])

pca <- prcomp(pca_matrix, center = TRUE, scale. = FALSE)

pca_var <- (pca$sdev^2) / sum(pca$sdev^2) * 100

pca_df <- data.frame(
  sample = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = if (ncol(pca$x) >= 2) pca$x[, 2] else 0,
  stringsAsFactors = FALSE
)

pca_df <- merge(pca_df, sample_qc[, c("sample", "group")], by = "sample", all.x = TRUE)
pca_df <- pca_df[match(rownames(pca$x), pca_df$sample), ]

pca_file <- file.path(output_dir, paste0(safe_label, "_pca_scores.tsv"))

write.table(
  pca_df,
  pca_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("[OK] PCA scores saved:", pca_file, "\n")

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(shape = group), size = 3, alpha = 0.9) +
  labs(
    title = paste0(label, " PCA"),
    subtitle = paste0(
      "Normalization: ", normalization,
      " | PCA features: ", length(pca_features)
    ),
    x = paste0("PC1 (", round(pca_var[1], 1), "%)"),
    y = paste0("PC2 (", round(pca_var[2], 1), "%)"),
    shape = "Group"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

if (has_ggrepel) {
  p_pca <- p_pca + ggrepel::geom_text_repel(aes(label = sample), size = 3, max.overlaps = Inf)
} else {
  p_pca <- p_pca + geom_text(aes(label = sample), size = 3, vjust = -0.8, check_overlap = TRUE)
}

pca_png <- file.path(output_dir, paste0(safe_label, "_pca.png"))

ggsave(
  pca_png,
  p_pca,
  width = plot_width,
  height = plot_height,
  dpi = plot_dpi
)

cat("[OK] PCA PNG saved:", pca_png, "\n")

if (save_pdf) {
  pca_pdf <- file.path(output_dir, paste0(safe_label, "_pca.pdf"))
  ggsave(pca_pdf, p_pca, width = plot_width, height = plot_height)
  cat("[OK] PCA PDF saved:", pca_pdf, "\n")
}

# -----------------------------
# Top variable features heatmap
# -----------------------------

heat_matrix <- norm_matrix[top_features, , drop = FALSE]
heat_z <- row_zscore(heat_matrix)

# Cluster features/samples when possible
feature_order <- rownames(heat_z)
sample_order <- colnames(heat_z)

if (nrow(heat_z) >= 2) {
  feature_order <- rownames(heat_z)[hclust(dist(heat_z))$order]
}

if (ncol(heat_z) >= 2) {
  sample_order <- colnames(heat_z)[hclust(dist(t(heat_z)))$order]
}

heat_z <- heat_z[feature_order, sample_order, drop = FALSE]

heat_long <- matrix_to_long(heat_z, "zscore")

heat_long$feature <- factor(heat_long$feature, levels = rev(feature_order))
heat_long$sample <- factor(heat_long$sample, levels = sample_order)

p_heat <- ggplot(heat_long, aes(x = sample, y = feature, fill = zscore)) +
  geom_tile() +
  labs(
    title = paste0(label, " top variable features heatmap"),
    subtitle = paste0("Top variable features: ", top_variable, " | Values: row z-score"),
    x = "Sample",
    y = "Feature",
    fill = "Z-score"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

heat_png <- file.path(output_dir, paste0(safe_label, "_top_variable_heatmap.png"))

ggsave(
  heat_png,
  p_heat,
  width = plot_width,
  height = max(plot_height, min(14, 0.18 * top_variable + 3)),
  dpi = plot_dpi
)

cat("[OK] Top variable heatmap PNG saved:", heat_png, "\n")

if (save_pdf) {
  heat_pdf <- file.path(output_dir, paste0(safe_label, "_top_variable_heatmap.pdf"))
  ggsave(
    heat_pdf,
    p_heat,
    width = plot_width,
    height = max(plot_height, min(14, 0.18 * top_variable + 3))
  )
  cat("[OK] Top variable heatmap PDF saved:", heat_pdf, "\n")
}

# -----------------------------
# Sample correlation heatmap
# -----------------------------

cor_matrix <- cor(var_matrix, method = "pearson", use = "pairwise.complete.obs")
cor_matrix[is.na(cor_matrix)] <- 0

cor_order <- colnames(cor_matrix)

if (ncol(cor_matrix) >= 2) {
  cor_order <- colnames(cor_matrix)[hclust(as.dist(1 - cor_matrix))$order]
}

cor_matrix <- cor_matrix[cor_order, cor_order, drop = FALSE]

cor_file <- file.path(output_dir, paste0(safe_label, "_sample_correlation.tsv"))

write.table(
  data.frame(sample = rownames(cor_matrix), cor_matrix, check.names = FALSE),
  cor_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("[OK] Sample correlation table saved:", cor_file, "\n")

cor_long <- as.data.frame(as.table(cor_matrix), stringsAsFactors = FALSE)
colnames(cor_long) <- c("sample_x", "sample_y", "correlation")
cor_long$sample_x <- factor(cor_long$sample_x, levels = cor_order)
cor_long$sample_y <- factor(cor_long$sample_y, levels = rev(cor_order))

p_cor <- ggplot(cor_long, aes(x = sample_x, y = sample_y, fill = correlation)) +
  geom_tile() +
  labs(
    title = paste0(label, " sample correlation heatmap"),
    subtitle = paste0("Pearson correlation | Normalization: ", normalization),
    x = "Sample",
    y = "Sample",
    fill = "r"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

cor_png <- file.path(output_dir, paste0(safe_label, "_sample_correlation_heatmap.png"))

ggsave(
  cor_png,
  p_cor,
  width = plot_width,
  height = plot_height,
  dpi = plot_dpi
)

cat("[OK] Sample correlation heatmap PNG saved:", cor_png, "\n")

if (save_pdf) {
  cor_pdf <- file.path(output_dir, paste0(safe_label, "_sample_correlation_heatmap.pdf"))
  ggsave(cor_pdf, p_cor, width = plot_width, height = plot_height)
  cat("[OK] Sample correlation heatmap PDF saved:", cor_pdf, "\n")
}

# -----------------------------
# Detected features per sample
# -----------------------------

sample_qc$sample <- factor(
  sample_qc$sample,
  levels = sample_qc$sample[order(sample_qc$detected_features, decreasing = TRUE)]
)

p_detected <- ggplot(sample_qc, aes(x = sample, y = detected_features)) +
  geom_col(aes(fill = group), width = 0.75) +
  labs(
    title = paste0(label, " detected features per sample"),
    subtitle = paste0("Detected threshold: > ", detected_threshold),
    x = "Sample",
    y = "Detected features",
    fill = "Group"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

detected_png <- file.path(output_dir, paste0(safe_label, "_detected_features_per_sample.png"))

ggsave(
  detected_png,
  p_detected,
  width = plot_width,
  height = plot_height,
  dpi = plot_dpi
)

cat("[OK] Detected features PNG saved:", detected_png, "\n")

if (save_pdf) {
  detected_pdf <- file.path(output_dir, paste0(safe_label, "_detected_features_per_sample.pdf"))
  ggsave(detected_pdf, p_detected, width = plot_width, height = plot_height)
  cat("[OK] Detected features PDF saved:", detected_pdf, "\n")
}

cat("============================================================\n")
cat("[DONE] Matrix QC completed.\n")
cat("Label          :", label, "\n")
cat("Matrix type    :", matrix_type, "\n")
cat("Normalization  :", normalization, "\n")
cat("Features       :", nrow(raw_matrix), "\n")
cat("Samples        :", ncol(raw_matrix), "\n")
cat("Output         :", output_dir, "\n")
cat("============================================================\n")
