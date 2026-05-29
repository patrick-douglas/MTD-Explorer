#!/usr/bin/env Rscript

# ============================================================
# MTD exploratory plot: Alpha diversity
#
# Metrics:
#   Observed taxa
#   Shannon
#   Simpson
#   Inverse Simpson
#   Pielou evenness
#
# Outputs:
#   alpha_diversity_<rank>_<mode>.tsv
#   alpha_diversity_<rank>_<mode>_faceted.png/pdf
#   alpha_diversity_<rank>_<mode>_<metric>.png/pdf
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  hit <- which(args == flag)
  if (length(hit) == 0) return(default)
  if (hit == length(args)) return(default)
  args[hit + 1]
}

has_flag <- function(flag) flag %in% args

print_help <- function() {
  cat("
============================================================
MTD exploratory: Alpha diversity
============================================================

Usage:

  Rscript MTD_exploratory_alpha_diversity.R \\
    --input bracken_genus_all \\
    --samplesheet samplesheet.csv \\
    --output MTD_res/exploratory/taxonomy/alpha_diversity/genus_relative \\
    --rank genus \\
    --mode relative

Required:
  --input FILE
      Combined Bracken abundance table.

Optional:
  --samplesheet FILE
      Samplesheet CSV. If provided, the script uses it for sample grouping.

  --output DIR
      Output directory. Default: current directory.

  --rank STRING
      Taxonomic rank label. Default: taxa

  --mode relative|absolute
      relative = convert each sample to percentage.
      absolute = use raw values.
      Default: relative

  --taxon_col COLUMN
      Taxon column name. Auto-detected if not provided.

  --sample_id_col COLUMN
      Sample ID column in samplesheet. Default: first column.

  --group_col COLUMN
      Group column in samplesheet. Default: second column.

  --presence_threshold NUMERIC
      Minimum abundance to count a taxon as present.
      Default: 0

  --width NUMERIC
      Plot width in inches. Default: 9

  --height NUMERIC
      Plot height in inches. Default: 6

  --dpi INTEGER
      PNG resolution. Default: 300

  --no_pdf
      Do not save PDF.

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
samplesheet_file <- get_arg("--samplesheet", NULL)
output_dir <- get_arg("--output", ".")
rank_label <- get_arg("--rank", "taxa")
mode <- get_arg("--mode", "relative")
taxon_col_user <- get_arg("--taxon_col", NULL)
sample_id_col_user <- get_arg("--sample_id_col", NULL)
group_col_user <- get_arg("--group_col", NULL)
presence_threshold <- as.numeric(get_arg("--presence_threshold", "0"))
plot_width <- as.numeric(get_arg("--width", "9"))
plot_height <- as.numeric(get_arg("--height", "6"))
plot_dpi <- as.integer(get_arg("--dpi", "300"))
save_pdf <- !has_flag("--no_pdf")

if (is.null(input_file)) stop("[ERROR] Missing --input")
if (!file.exists(input_file)) stop(paste0("[ERROR] Input not found: ", input_file))
if (!mode %in% c("relative", "absolute")) stop("[ERROR] --mode must be relative or absolute")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

required_packages <- c("ggplot2")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste0("[ERROR] Required R package not installed: ", pkg,
                "\nInstall with: install.packages('", pkg, "')"))
  }
}

suppressPackageStartupMessages(library(ggplot2))

detect_sep <- function(file) {
  first_line <- readLines(file, n = 1, warn = FALSE)
  n_tab <- lengths(regmatches(first_line, gregexpr("\t", first_line)))
  n_comma <- lengths(regmatches(first_line, gregexpr(",", first_line)))
  n_semicolon <- lengths(regmatches(first_line, gregexpr(";", first_line)))

  if (n_tab >= n_comma && n_tab >= n_semicolon && n_tab > 0) "\t"
  else if (n_semicolon > n_comma) ";"
  else ","
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

  colnames(ss) <- clean_column_names(colnames(ss))

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

alpha_metrics <- function(v, presence_threshold = 0) {
  v[is.na(v)] <- 0
  v[v < 0] <- 0

  observed <- sum(v > presence_threshold)
  total <- sum(v)

  if (total <= 0 || observed == 0) {
    return(c(
      observed_taxa = observed,
      shannon = NA,
      simpson = NA,
      inverse_simpson = NA,
      pielou_evenness = NA,
      total_abundance = total
    ))
  }

  p <- v[v > 0] / total

  shannon <- -sum(p * log(p))
  simpson <- 1 - sum(p^2)
  inverse_simpson <- 1 / sum(p^2)

  if (observed > 1) {
    pielou <- shannon / log(observed)
  } else {
    pielou <- NA
  }

  c(
    observed_taxa = observed,
    shannon = shannon,
    simpson = simpson,
    inverse_simpson = inverse_simpson,
    pielou_evenness = pielou,
    total_abundance = total
  )
}

# -----------------------------
# Read abundance table
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
cat("MTD exploratory: alpha diversity\n")
cat("Input file        :", input_file, "\n")
cat("Samplesheet       :", ifelse(is.null(samplesheet_file), "none", samplesheet_file), "\n")
cat("Rank              :", rank_label, "\n")
cat("Mode              :", mode, "\n")
cat("Presence threshold:", presence_threshold, "\n")
cat("Output dir        :", output_dir, "\n")
cat("============================================================\n")

cols <- colnames(df)

taxon_candidates <- c(
  "taxon", "taxa", "name", "scientific_name",
  "species", "genus", "family", "phylum",
  "clade_name", "feature", "ID", "id"
)

known_metadata_cols <- c(
  "taxonomy_id", "taxid", "tax_id", "NCBI_tax_id",
  "taxonomy_lvl", "rank", "level",
  "kraken_assigned_reads", "added_reads",
  "new_est_reads", "fraction_total_reads"
)

if (!is.null(taxon_col_user)) {
  taxon_col <- taxon_col_user
} else {
  taxon_col <- find_first_existing(taxon_candidates, cols)
}

if (is.null(taxon_col)) {
  non_numeric_cols <- cols[!sapply(df, is_numeric_like)]
  if (length(non_numeric_cols) > 0) taxon_col <- non_numeric_cols[1] else taxon_col <- cols[1]
}

if (!taxon_col %in% cols) stop(paste0("[ERROR] Taxon column not found: ", taxon_col))

possible_sample_cols <- setdiff(cols, c(taxon_col, known_metadata_cols))

sample_cols <- possible_sample_cols[
  sapply(df[, possible_sample_cols, drop = FALSE], is_numeric_like)
]

# ------------------------------------------------------------
# Bracken combined tables often contain two numeric columns per sample:
#   sample_num  = estimated reads / abundance
#   sample_frac = fraction of total reads
#
# For alpha diversity, keep only one column per biological sample.
# We keep *_num by default because diversity metrics can be computed
# from counts/estimated abundances.
# ------------------------------------------------------------

num_cols <- sample_cols[grepl("(_num$|\\.num$)", sample_cols, ignore.case = TRUE)]
frac_cols <- sample_cols[grepl("(_frac$|\\.frac$)", sample_cols, ignore.case = TRUE)]

if (length(num_cols) >= 2) {
  cat("[INFO] Bracken-style _num/_frac columns detected.\n")
  cat("[INFO] Keeping only _num columns for alpha diversity.\n")
  sample_cols <- num_cols
} else if (length(frac_cols) >= 2) {
  cat("[INFO] Bracken-style _frac columns detected, but no _num columns were found.\n")
  cat("[INFO] Keeping _frac columns for alpha diversity.\n")
  sample_cols <- frac_cols
}

clean_sample_names <- sample_cols

clean_sample_names <- gsub("^Report_", "", clean_sample_names)
clean_sample_names <- gsub("_num$", "", clean_sample_names, ignore.case = TRUE)
clean_sample_names <- gsub("_frac$", "", clean_sample_names, ignore.case = TRUE)
clean_sample_names <- gsub("\\.num$", "", clean_sample_names, ignore.case = TRUE)
clean_sample_names <- gsub("\\.frac$", "", clean_sample_names, ignore.case = TRUE)

# Remove common Bracken rank/file suffixes
clean_sample_names <- gsub("_phylum.*$", "", clean_sample_names, ignore.case = TRUE)
clean_sample_names <- gsub("_genus.*$", "", clean_sample_names, ignore.case = TRUE)
clean_sample_names <- gsub("_species.*$", "", clean_sample_names, ignore.case = TRUE)

clean_sample_names <- gsub("_bracken$", "", clean_sample_names, ignore.case = TRUE)
clean_sample_names <- gsub("\\.bracken$", "", clean_sample_names, ignore.case = TRUE)

names(sample_cols) <- clean_sample_names

if (length(sample_cols) == 0) stop("[ERROR] Could not detect numeric sample columns.")

taxa <- make_unique_taxa(df[[taxon_col]])

abundance_matrix <- as.matrix(
  data.frame(
    lapply(df[, sample_cols, drop = FALSE], to_numeric_safe),
    check.names = FALSE
  )
)

rownames(abundance_matrix) <- taxa
if (!is.null(names(sample_cols)) && all(!is.na(names(sample_cols))) && all(names(sample_cols) != "")) {
  colnames(abundance_matrix) <- names(sample_cols)
} else {
  colnames(abundance_matrix) <- sample_cols
}
abundance_matrix[is.na(abundance_matrix)] <- 0
abundance_matrix[abundance_matrix < 0] <- 0
abundance_matrix <- abundance_matrix[rowSums(abundance_matrix) > 0, , drop = FALSE]

if (nrow(abundance_matrix) == 0) stop("[ERROR] No taxa with abundance > 0.")

plot_matrix <- abundance_matrix

if (mode == "relative") {
  max_val <- max(plot_matrix, na.rm = TRUE)
  col_sums <- colSums(plot_matrix, na.rm = TRUE)

  if (max_val <= 1) {
    cat("[INFO] Input appears fractional. Multiplying by 100.\n")
    plot_matrix <- plot_matrix * 100
  } else if (max_val <= 100 && median(col_sums, na.rm = TRUE) <= 101) {
    cat("[INFO] Input appears already percentage relative abundance.\n")
  } else {
    cat("[INFO] Input appears counts. Converting samples to percentage.\n")
    col_sums[col_sums == 0] <- NA
    plot_matrix <- sweep(plot_matrix, 2, col_sums, "/") * 100
    plot_matrix[is.na(plot_matrix)] <- 0
  }
}

# -----------------------------
# Alpha calculation
# -----------------------------

alpha_mat <- t(apply(plot_matrix, 2, alpha_metrics, presence_threshold = presence_threshold))
alpha_df <- as.data.frame(alpha_mat)
alpha_df$sample <- rownames(alpha_df)
alpha_df$rank <- rank_label
alpha_df$mode <- mode

alpha_df <- alpha_df[, c(
  "sample", "rank", "mode",
  "observed_taxa", "shannon", "simpson",
  "inverse_simpson", "pielou_evenness", "total_abundance"
)]

ss <- read_samplesheet(samplesheet_file, sample_id_col_user, group_col_user)

if (!is.null(ss)) {
  alpha_df <- merge(alpha_df, ss, by = "sample", all.x = TRUE)
  alpha_df$group[is.na(alpha_df$group) | alpha_df$group == ""] <- "Ungrouped"
} else {
  alpha_df$group <- "All_samples"
}

safe_rank <- safe_name(rank_label)
safe_mode <- safe_name(mode)

out_prefix <- file.path(output_dir, paste0("alpha_diversity_", safe_rank, "_", safe_mode))
table_file <- paste0(out_prefix, ".tsv")

write.table(alpha_df, table_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat("[OK] Alpha diversity table saved:", table_file, "\n")

# -----------------------------
# Long format for plots
# -----------------------------

metrics <- c("observed_taxa", "shannon", "simpson", "inverse_simpson", "pielou_evenness")

long_df <- do.call(
  rbind,
  lapply(metrics, function(m) {
    data.frame(
      sample = alpha_df$sample,
      group = alpha_df$group,
      metric = m,
      value = alpha_df[[m]],
      stringsAsFactors = FALSE
    )
  })
)

metric_labels <- c(
  observed_taxa = "Observed taxa",
  shannon = "Shannon",
  simpson = "Simpson",
  inverse_simpson = "Inverse Simpson",
  pielou_evenness = "Pielou evenness"
)

long_df$metric_label <- metric_labels[long_df$metric]

has_groups <- length(unique(long_df$group)) > 1

# Faceted plot

if (has_groups) {
  p_faceted <- ggplot(long_df, aes(x = group, y = value)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.65) +
    geom_jitter(width = 0.15, size = 2, alpha = 0.85) +
    facet_wrap(~ metric_label, scales = "free_y") +
    labs(
      title = paste0("Alpha diversity at ", rank_label, " level"),
      subtitle = paste0("Mode: ", mode, " | Presence threshold: > ", presence_threshold),
      x = "Group",
      y = "Alpha diversity value"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 35, hjust = 1),
      panel.grid.minor = element_blank()
    )
} else {
  p_faceted <- ggplot(long_df, aes(x = sample, y = value)) +
    geom_point(size = 2.5, alpha = 0.9) +
    facet_wrap(~ metric_label, scales = "free_y") +
    labs(
      title = paste0("Alpha diversity at ", rank_label, " level"),
      subtitle = paste0("Mode: ", mode, " | Presence threshold: > ", presence_threshold),
      x = "Sample",
      y = "Alpha diversity value"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank()
    )
}

faceted_png <- paste0(out_prefix, "_faceted.png")
ggsave(faceted_png, p_faceted, width = plot_width, height = plot_height, dpi = plot_dpi)
cat("[OK] Faceted PNG saved:", faceted_png, "\n")

if (save_pdf) {
  faceted_pdf <- paste0(out_prefix, "_faceted.pdf")
  ggsave(faceted_pdf, p_faceted, width = plot_width, height = plot_height)
  cat("[OK] Faceted PDF saved:", faceted_pdf, "\n")
}

# Individual plots

for (m in metrics) {
  dfm <- long_df[long_df$metric == m, , drop = FALSE]
  metric_name <- metric_labels[[m]]

  if (has_groups) {
    p <- ggplot(dfm, aes(x = group, y = value)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.65) +
      geom_jitter(width = 0.15, size = 2.3, alpha = 0.85) +
      labs(
        title = paste0(metric_name, " at ", rank_label, " level"),
        subtitle = paste0("Mode: ", mode, " | Presence threshold: > ", presence_threshold),
        x = "Group",
        y = metric_name
      ) +
      theme_bw(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 35, hjust = 1),
        panel.grid.minor = element_blank()
      )
  } else {
    p <- ggplot(dfm, aes(x = sample, y = value)) +
      geom_point(size = 2.5, alpha = 0.9) +
      labs(
        title = paste0(metric_name, " at ", rank_label, " level"),
        subtitle = paste0("Mode: ", mode, " | Presence threshold: > ", presence_threshold),
        x = "Sample",
        y = metric_name
      ) +
      theme_bw(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank()
      )
  }

  png_file <- paste0(out_prefix, "_", safe_name(m), ".png")
  ggsave(png_file, p, width = plot_width, height = plot_height, dpi = plot_dpi)
  cat("[OK] PNG saved:", png_file, "\n")

  if (save_pdf) {
    pdf_file <- paste0(out_prefix, "_", safe_name(m), ".pdf")
    ggsave(pdf_file, p, width = plot_width, height = plot_height)
    cat("[OK] PDF saved:", pdf_file, "\n")
  }
}

cat("============================================================\n")
cat("[DONE] Alpha diversity completed.\n")
cat("Samples:", ncol(plot_matrix), "\n")
cat("Taxa   :", nrow(plot_matrix), "\n")
cat("Output :", output_dir, "\n")
cat("============================================================\n")
