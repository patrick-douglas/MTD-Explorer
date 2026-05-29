#!/usr/bin/env Rscript

# ============================================================
# MTD exploratory plot: Read composition QC
#
# Generates:
#   A) Sankey-like read flow summary
#   B) Per-sample percentage stacked barplot
#   C) Outlier plot by classification stage
#
# Input:
#   kraken_global_read_composition_final.tsv
#
# Expected columns:
#   sample
#   total_reads
#   host           e.g. 12345 (80.12%)
#   microbiome     e.g. 123 (0.80%)
#   unclassified   e.g. 2934 (19.08%)
#   check_pct_sum
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
MTD exploratory: Read composition QC
============================================================

Usage:

  Rscript MTD_exploratory_read_composition_qc.R \\
    --input MTD_res/kraken/kraken_global_read_composition_final.tsv \\
    --output MTD_res/exploratory/pipeline_qc/read_composition \\
    --label final

Required:
  --input FILE
      Kraken global read composition table.

Optional:
  --output DIR
      Output directory.
      Default: current directory.

  --label STRING
      Label used in filenames and titles.
      Example: raw, final.
      Default: final

  --outlier_method mad|iqr
      Method used to flag outliers.
      Default: mad

  --outlier_cutoff NUMERIC
      Cutoff for modified z-score when --outlier_method mad.
      Default: 3.5

  --iqr_multiplier NUMERIC
      IQR multiplier when --outlier_method iqr.
      Default: 1.5

  --width NUMERIC
      Plot width in inches.
      Default: 10

  --height NUMERIC
      Plot height in inches.
      Default: 6

  --dpi INTEGER
      PNG resolution.
      Default: 300

  --no_pdf
      Do not save PDF.

  --help
      Show this help message.

Outputs:
  read_composition_qc_<label>_clean.tsv
  read_composition_qc_<label>_long.tsv
  read_composition_qc_<label>_summary.tsv
  read_composition_qc_<label>_outliers.tsv

  read_composition_flow_<label>.png/pdf
  read_composition_percent_barplot_<label>.png/pdf
  read_composition_outlier_plot_<label>.png/pdf

============================================================
\n")
}

if (has_flag("--help") || length(args) == 0) {
  print_help()
  quit(save = "no", status = 0)
}

input_file <- get_arg("--input")
output_dir <- get_arg("--output", ".")
label <- get_arg("--label", "final")
outlier_method <- get_arg("--outlier_method", "mad")
outlier_cutoff <- as.numeric(get_arg("--outlier_cutoff", "3.5"))
iqr_multiplier <- as.numeric(get_arg("--iqr_multiplier", "1.5"))
plot_width <- as.numeric(get_arg("--width", "10"))
plot_height <- as.numeric(get_arg("--height", "6"))
plot_dpi <- as.integer(get_arg("--dpi", "300"))
save_pdf <- !has_flag("--no_pdf")

if (is.null(input_file)) {
  stop("[ERROR] Missing required argument: --input")
}

if (!file.exists(input_file)) {
  stop(paste0("[ERROR] Input file not found: ", input_file))
}

if (!outlier_method %in% c("mad", "iqr")) {
  stop("[ERROR] --outlier_method must be 'mad' or 'iqr'")
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

# -----------------------------
# Helpers
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

safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9_\\-]+", "_", x)
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}

to_numeric_safe <- function(x) {
  x <- as.character(x)
  x <- gsub(",", ".", x)
  suppressWarnings(as.numeric(x))
}

extract_count <- function(x) {
  x <- as.character(x)
  out <- sub("^\\s*([0-9.]+).*", "\\1", x)
  to_numeric_safe(out)
}

extract_percent <- function(x) {
  x <- as.character(x)

  has_paren <- grepl("\\(([0-9.]+)%\\)", x)

  out <- rep(NA_real_, length(x))

  out[has_paren] <- to_numeric_safe(
    sub(".*\\(([0-9.]+)%\\).*", "\\1", x[has_paren])
  )

  no_paren <- !has_paren
  out[no_paren] <- to_numeric_safe(x[no_paren])

  out
}

fmt_count <- function(x) {
  format(round(x), big.mark = ",", scientific = FALSE, trim = TRUE)
}

fmt_pct <- function(x) {
  sprintf("%.2f%%", x)
}

stage_label <- function(x) {
  d <- c(
    host = "Host",
    microbiome = "Microbiome",
    unclassified = "Unclassified"
  )
  unname(d[x])
}

modified_z <- function(x) {
  med <- median(x, na.rm = TRUE)
  mad_val <- mad(x, constant = 1, na.rm = TRUE)

  if (is.na(mad_val) || mad_val == 0) {
    return(rep(0, length(x)))
  }

  0.6745 * (x - med) / mad_val
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
cat("MTD exploratory: read composition QC\n")
cat("Input file     :", input_file, "\n")
cat("Separator      :", ifelse(sep == "\t", "TAB", sep), "\n")
cat("Rows           :", nrow(df), "\n")
cat("Columns        :", ncol(df), "\n")
cat("Label          :", label, "\n")
cat("Outlier method :", outlier_method, "\n")
cat("Output dir     :", output_dir, "\n")
cat("============================================================\n")

required_cols <- c("sample", "total_reads", "host", "microbiome", "unclassified")
missing_cols <- setdiff(required_cols, colnames(df))

if (length(missing_cols) > 0) {
  stop(paste0(
    "[ERROR] Missing required columns: ",
    paste(missing_cols, collapse = ", ")
  ))
}

# -----------------------------
# Parse composition table
# -----------------------------

clean_df <- data.frame(
  sample = as.character(df$sample),
  total_reads = to_numeric_safe(df$total_reads),

  host_reads = extract_count(df$host),
  microbiome_reads = extract_count(df$microbiome),
  unclassified_reads = extract_count(df$unclassified),

  host_pct = extract_percent(df$host),
  microbiome_pct = extract_percent(df$microbiome),
  unclassified_pct = extract_percent(df$unclassified),

  stringsAsFactors = FALSE
)

# Recalculate percentages if needed
for (stage in c("host", "microbiome", "unclassified")) {
  pct_col <- paste0(stage, "_pct")
  reads_col <- paste0(stage, "_reads")

  missing_pct <- is.na(clean_df[[pct_col]])

  clean_df[[pct_col]][missing_pct] <- (
    clean_df[[reads_col]][missing_pct] / clean_df$total_reads[missing_pct]
  ) * 100
}

clean_df$check_pct_sum_recalculated <- clean_df$host_pct +
  clean_df$microbiome_pct +
  clean_df$unclassified_pct

clean_df$check_reads_sum <- clean_df$host_reads +
  clean_df$microbiome_reads +
  clean_df$unclassified_reads

clean_df$check_reads_difference <- clean_df$total_reads - clean_df$check_reads_sum

clean_df <- clean_df[!is.na(clean_df$sample) & clean_df$sample != "", , drop = FALSE]

if (nrow(clean_df) == 0) {
  stop("[ERROR] No valid samples found after parsing input.")
}

# -----------------------------
# Long table
# -----------------------------

long_df <- rbind(
  data.frame(
    sample = clean_df$sample,
    stage = "host",
    stage_label = "Host",
    reads = clean_df$host_reads,
    percent = clean_df$host_pct,
    total_reads = clean_df$total_reads,
    stringsAsFactors = FALSE
  ),
  data.frame(
    sample = clean_df$sample,
    stage = "microbiome",
    stage_label = "Microbiome",
    reads = clean_df$microbiome_reads,
    percent = clean_df$microbiome_pct,
    total_reads = clean_df$total_reads,
    stringsAsFactors = FALSE
  ),
  data.frame(
    sample = clean_df$sample,
    stage = "unclassified",
    stage_label = "Unclassified",
    reads = clean_df$unclassified_reads,
    percent = clean_df$unclassified_pct,
    total_reads = clean_df$total_reads,
    stringsAsFactors = FALSE
  )
)

long_df$stage <- factor(
  long_df$stage,
  levels = c("host", "microbiome", "unclassified")
)

long_df$stage_label <- factor(
  long_df$stage_label,
  levels = c("Host", "Microbiome", "Unclassified")
)

# -----------------------------
# Summary table
# -----------------------------

summary_df <- aggregate(
  cbind(reads, percent) ~ stage + stage_label,
  data = long_df,
  FUN = function(x) c(
    mean = mean(x, na.rm = TRUE),
    median = median(x, na.rm = TRUE),
    min = min(x, na.rm = TRUE),
    max = max(x, na.rm = TRUE),
    sd = sd(x, na.rm = TRUE)
  )
)

summary_out <- data.frame(
  stage = summary_df$stage,
  stage_label = summary_df$stage_label,

  mean_reads = summary_df$reads[, "mean"],
  median_reads = summary_df$reads[, "median"],
  min_reads = summary_df$reads[, "min"],
  max_reads = summary_df$reads[, "max"],
  sd_reads = summary_df$reads[, "sd"],

  mean_percent = summary_df$percent[, "mean"],
  median_percent = summary_df$percent[, "median"],
  min_percent = summary_df$percent[, "min"],
  max_percent = summary_df$percent[, "max"],
  sd_percent = summary_df$percent[, "sd"],

  stringsAsFactors = FALSE
)

# -----------------------------
# Outlier detection
# -----------------------------

outlier_rows <- list()

for (st in levels(long_df$stage)) {
  idx <- which(long_df$stage == st)
  x <- long_df$percent[idx]

  if (outlier_method == "mad") {
    z <- modified_z(x)
    is_outlier <- abs(z) >= outlier_cutoff

    tmp <- data.frame(
      sample = long_df$sample[idx],
      stage = as.character(long_df$stage[idx]),
      stage_label = as.character(long_df$stage_label[idx]),
      reads = long_df$reads[idx],
      percent = long_df$percent[idx],
      outlier_method = "mad",
      outlier_score = z,
      outlier_cutoff = outlier_cutoff,
      is_outlier = is_outlier,
      stringsAsFactors = FALSE
    )

  } else {
    q1 <- quantile(x, 0.25, na.rm = TRUE)
    q3 <- quantile(x, 0.75, na.rm = TRUE)
    iqr <- q3 - q1

    lower <- q1 - iqr_multiplier * iqr
    upper <- q3 + iqr_multiplier * iqr

    is_outlier <- x < lower | x > upper

    tmp <- data.frame(
      sample = long_df$sample[idx],
      stage = as.character(long_df$stage[idx]),
      stage_label = as.character(long_df$stage_label[idx]),
      reads = long_df$reads[idx],
      percent = long_df$percent[idx],
      outlier_method = "iqr",
      outlier_score = NA_real_,
      outlier_cutoff = NA_real_,
      iqr_lower = lower,
      iqr_upper = upper,
      is_outlier = is_outlier,
      stringsAsFactors = FALSE
    )
  }

  outlier_rows[[st]] <- tmp
}

outlier_df <- do.call(rbind, outlier_rows)
long_df <- merge(
  long_df,
  outlier_df[, c("sample", "stage", "is_outlier")],
  by = c("sample", "stage"),
  all.x = TRUE
)

long_df$is_outlier[is.na(long_df$is_outlier)] <- FALSE

# -----------------------------
# Save tables
# -----------------------------

safe_label <- safe_name(label)

clean_file <- file.path(output_dir, paste0("read_composition_qc_", safe_label, "_clean.tsv"))
long_file <- file.path(output_dir, paste0("read_composition_qc_", safe_label, "_long.tsv"))
summary_file <- file.path(output_dir, paste0("read_composition_qc_", safe_label, "_summary.tsv"))
outlier_file <- file.path(output_dir, paste0("read_composition_qc_", safe_label, "_outliers.tsv"))

write.table(clean_df, clean_file, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(long_df, long_file, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(summary_out, summary_file, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(outlier_df[outlier_df$is_outlier, , drop = FALSE], outlier_file, sep = "\t", quote = FALSE, row.names = FALSE)

cat("[OK] Clean table saved:", clean_file, "\n")
cat("[OK] Long table saved:", long_file, "\n")
cat("[OK] Summary table saved:", summary_file, "\n")
cat("[OK] Outlier table saved:", outlier_file, "\n")

# -----------------------------
# A) Sankey-like / read flow plot
# -----------------------------

flow_summary <- data.frame(
  stage = factor(
    c("Total reads", "Host", "Microbiome", "Unclassified"),
    levels = c("Total reads", "Host", "Microbiome", "Unclassified")
  ),
  reads = c(
    sum(clean_df$total_reads, na.rm = TRUE),
    sum(clean_df$host_reads, na.rm = TRUE),
    sum(clean_df$microbiome_reads, na.rm = TRUE),
    sum(clean_df$unclassified_reads, na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)

flow_summary$percent_of_total <- (flow_summary$reads / flow_summary$reads[flow_summary$stage == "Total reads"]) * 100
flow_summary$x <- c(1, 2, 3, 3)
flow_summary$y <- c(0.5, 0.5, 0.70, 0.30)
flow_summary$label <- paste0(
  flow_summary$stage,
  "\n",
  fmt_count(flow_summary$reads),
  " reads\n",
  fmt_pct(flow_summary$percent_of_total)
)

flow_edges <- data.frame(
  x = c(1.25, 2.25, 2.25),
  xend = c(1.75, 2.75, 2.75),
  y = c(0.5, 0.5, 0.5),
  yend = c(0.5, 0.70, 0.30),
  stage = c("Host", "Microbiome", "Unclassified"),
  reads = c(
    sum(clean_df$host_reads, na.rm = TRUE),
    sum(clean_df$microbiome_reads, na.rm = TRUE),
    sum(clean_df$unclassified_reads, na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)

max_edge <- max(flow_edges$reads, na.rm = TRUE)
flow_edges$edge_size <- ifelse(max_edge > 0, 2 + 10 * flow_edges$reads / max_edge, 2)

p_flow <- ggplot() +
  geom_curve(
    data = flow_edges,
    aes(x = x, y = y, xend = xend, yend = yend, linewidth = edge_size),
    curvature = 0.20,
    alpha = 0.45,
    lineend = "round"
  ) +
  geom_label(
    data = flow_summary,
    aes(x = x, y = y, label = label),
    size = 4,
    label.size = 0.25,
    label.padding = unit(0.35, "lines")
  ) +
  scale_linewidth_identity() +
  coord_cartesian(xlim = c(0.6, 3.4), ylim = c(0.05, 0.95), clip = "off") +
  labs(
    title = paste0("Read classification flow - ", label),
    subtitle = paste0("Aggregated across ", nrow(clean_df), " samples"),
    x = NULL,
    y = NULL
  ) +
  theme_void(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(hjust = 0)
  )

flow_png <- file.path(output_dir, paste0("read_composition_flow_", safe_label, ".png"))
ggsave(flow_png, p_flow, width = plot_width, height = plot_height, dpi = plot_dpi)
cat("[OK] Flow PNG saved:", flow_png, "\n")

if (save_pdf) {
  flow_pdf <- file.path(output_dir, paste0("read_composition_flow_", safe_label, ".pdf"))
  ggsave(flow_pdf, p_flow, width = plot_width, height = plot_height)
  cat("[OK] Flow PDF saved:", flow_pdf, "\n")
}

# -----------------------------
# B) Per-sample percentage barplot
# -----------------------------

sample_order <- clean_df$sample[order(clean_df$host_pct, decreasing = TRUE)]
long_df$sample <- factor(long_df$sample, levels = sample_order)

p_bar <- ggplot(
  long_df,
  aes(x = sample, y = percent, fill = stage_label)
) +
  geom_col(width = 0.8) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, by = 20),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    title = paste0("Read classification percentage per sample - ", label),
    subtitle = "Samples are ordered by host percentage",
    x = "Sample",
    y = "Reads (%)",
    fill = "Classification"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

bar_png <- file.path(output_dir, paste0("read_composition_percent_barplot_", safe_label, ".png"))
ggsave(bar_png, p_bar, width = plot_width, height = plot_height, dpi = plot_dpi)
cat("[OK] Percentage barplot PNG saved:", bar_png, "\n")

if (save_pdf) {
  bar_pdf <- file.path(output_dir, paste0("read_composition_percent_barplot_", safe_label, ".pdf"))
  ggsave(bar_pdf, p_bar, width = plot_width, height = plot_height)
  cat("[OK] Percentage barplot PDF saved:", bar_pdf, "\n")
}

# -----------------------------
# C) Outlier plot
# -----------------------------

p_outlier <- ggplot(
  long_df,
  aes(x = stage_label, y = percent)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.55) +
  geom_jitter(
    aes(shape = is_outlier),
    width = 0.15,
    size = 2.5,
    alpha = 0.85
  ) +
  geom_text(
    data = long_df[long_df$is_outlier, , drop = FALSE],
    aes(label = sample),
    vjust = -0.8,
    size = 3,
    check_overlap = TRUE
  ) +
  labs(
    title = paste0("Read classification outlier check - ", label),
    subtitle = paste0(
      "Outlier method: ",
      outlier_method,
      ifelse(outlier_method == "mad",
             paste0(" | modified z cutoff: ", outlier_cutoff),
             paste0(" | IQR multiplier: ", iqr_multiplier))
    ),
    x = "Classification",
    y = "Reads (%)",
    shape = "Flagged outlier"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

outlier_png <- file.path(output_dir, paste0("read_composition_outlier_plot_", safe_label, ".png"))
ggsave(outlier_png, p_outlier, width = plot_width, height = plot_height, dpi = plot_dpi)
cat("[OK] Outlier plot PNG saved:", outlier_png, "\n")

if (save_pdf) {
  outlier_pdf <- file.path(output_dir, paste0("read_composition_outlier_plot_", safe_label, ".pdf"))
  ggsave(outlier_pdf, p_outlier, width = plot_width, height = plot_height)
  cat("[OK] Outlier plot PDF saved:", outlier_pdf, "\n")
}

cat("============================================================\n")
cat("[DONE] Read composition QC completed.\n")
cat("Samples analyzed:", nrow(clean_df), "\n")
cat("Total reads     :", fmt_count(sum(clean_df$total_reads, na.rm = TRUE)), "\n")
cat("Outliers found  :", sum(outlier_df$is_outlier, na.rm = TRUE), "\n")
cat("Output          :", output_dir, "\n")
cat("============================================================\n")
