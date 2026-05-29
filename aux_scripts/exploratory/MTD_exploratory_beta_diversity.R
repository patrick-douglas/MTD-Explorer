#!/usr/bin/env Rscript

# ============================================================
# MTD exploratory plot: Beta diversity / ordination
#
# Methods:
#   PCoA using Bray-Curtis distance
#   NMDS using Bray-Curtis distance
#   Optional PERMANOVA if samplesheet has >=2 groups
#
# Outputs:
#   beta_diversity_<rank>_<mode>_distance_matrix.tsv
#   beta_diversity_<rank>_<mode>_pcoa_scores.tsv
#   beta_diversity_<rank>_<mode>_pcoa_bray.png/pdf
#   beta_diversity_<rank>_<mode>_nmds_scores.tsv
#   beta_diversity_<rank>_<mode>_nmds_bray.png/pdf
#   beta_diversity_<rank>_<mode>_permanova.tsv
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
MTD exploratory: Beta diversity / ordination
============================================================

Usage:

  Rscript MTD_exploratory_beta_diversity.R \\
    --input bracken_genus_all \\
    --samplesheet samplesheet.csv \\
    --output MTD_res/exploratory/taxonomy/beta_diversity/genus_relative \\
    --rank genus \\
    --mode relative

Required:
  --input FILE
      Combined Bracken abundance table.

Optional:
  --samplesheet FILE
      Samplesheet CSV. If provided, groups are used for plotting and PERMANOVA.

  --output DIR
      Output directory. Default: current directory.

  --rank STRING
      Taxonomic rank label. Default: taxa

  --mode relative|absolute
      relative = convert each sample to percentage.
      absolute = use raw values.
      Default: relative

  --distance bray
      Distance metric for vegan::vegdist.
      Default: bray

  --transform none|sqrt|log10
      Optional transform after relative/absolute conversion.
      Default: none

  --taxon_col COLUMN
      Taxon column name. Auto-detected if not provided.

  --sample_id_col COLUMN
      Sample ID column in samplesheet. Default: first column.

  --group_col COLUMN
      Group column in samplesheet. Default: second column.

  --label_samples yes|no
      Label sample names in ordination plots.
      Default: yes

  --ellipse yes|no
      Add group ellipses if possible.
      Default: no

  --width NUMERIC
      Plot width in inches. Default: 8

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
distance_method <- get_arg("--distance", "bray")
transform_method <- get_arg("--transform", "none")
taxon_col_user <- get_arg("--taxon_col", NULL)
sample_id_col_user <- get_arg("--sample_id_col", NULL)
group_col_user <- get_arg("--group_col", NULL)
label_samples <- get_arg("--label_samples", "yes")
ellipse <- get_arg("--ellipse", "no")
plot_width <- as.numeric(get_arg("--width", "8"))
plot_height <- as.numeric(get_arg("--height", "6"))
plot_dpi <- as.integer(get_arg("--dpi", "300"))
save_pdf <- !has_flag("--no_pdf")

if (is.null(input_file)) stop("[ERROR] Missing --input")
if (!file.exists(input_file)) stop(paste0("[ERROR] Input not found: ", input_file))
if (!mode %in% c("relative", "absolute")) stop("[ERROR] --mode must be relative or absolute")
if (!transform_method %in% c("none", "sqrt", "log10")) stop("[ERROR] --transform must be none, sqrt or log10")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

required_packages <- c("ggplot2", "vegan")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste0("[ERROR] Required R package not installed: ", pkg,
                "\nInstall with: install.packages('", pkg, "')"))
  }
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(vegan)
})

has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)

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
cat("MTD exploratory: beta diversity / ordination\n")
cat("Input file  :", input_file, "\n")
cat("Samplesheet :", ifelse(is.null(samplesheet_file), "none", samplesheet_file), "\n")
cat("Rank        :", rank_label, "\n")
cat("Mode        :", mode, "\n")
cat("Distance    :", distance_method, "\n")
cat("Transform   :", transform_method, "\n")
cat("Output dir  :", output_dir, "\n")
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
# For beta diversity, keep only one column per biological sample.
# We keep *_num by default because Bray-Curtis can be computed from
# counts/estimated abundances, and relative normalization is done below.
# ------------------------------------------------------------

num_cols <- sample_cols[grepl("(_num$|\\.num$)", sample_cols, ignore.case = TRUE)]
frac_cols <- sample_cols[grepl("(_frac$|\\.frac$)", sample_cols, ignore.case = TRUE)]

if (length(num_cols) >= 2) {
  cat("[INFO] Bracken-style _num/_frac columns detected.\n")
  cat("[INFO] Keeping only _num columns for beta diversity.\n")
  sample_cols <- num_cols
} else if (length(frac_cols) >= 2) {
  cat("[INFO] Bracken-style _frac columns detected, but no _num columns were found.\n")
  cat("[INFO] Keeping _frac columns for beta diversity.\n")
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

if (length(sample_cols) < 2) {
  stop("[ERROR] Beta diversity requires at least 2 sample columns.")
}

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

if (transform_method == "sqrt") {
  plot_matrix <- sqrt(plot_matrix)
} else if (transform_method == "log10") {
  plot_matrix <- log10(plot_matrix + 1)
}

# samples x taxa
sample_taxa <- t(plot_matrix)

# Remove empty samples
sample_taxa <- sample_taxa[rowSums(sample_taxa) > 0, , drop = FALSE]

# Remove taxa that became zero after sample filtering
sample_taxa <- sample_taxa[, colSums(sample_taxa) > 0, drop = FALSE]

if (nrow(sample_taxa) < 2) stop("[ERROR] At least 2 non-empty samples are required.")
if (ncol(sample_taxa) < 1) stop("[ERROR] No non-zero taxa available for ordination.")

ss <- read_samplesheet(samplesheet_file, sample_id_col_user, group_col_user)

meta_df <- data.frame(
  sample = rownames(sample_taxa),
  group = "All_samples",
  stringsAsFactors = FALSE
)

if (!is.null(ss)) {
  meta_df <- merge(meta_df[, "sample", drop = FALSE], ss, by = "sample", all.x = TRUE)
  meta_df$group[is.na(meta_df$group) | meta_df$group == ""] <- "Ungrouped"
}

rownames(meta_df) <- meta_df$sample
meta_df <- meta_df[rownames(sample_taxa), , drop = FALSE]

safe_rank <- safe_name(rank_label)
safe_mode <- safe_name(mode)
safe_distance <- safe_name(distance_method)

out_prefix <- file.path(
  output_dir,
  paste0("beta_diversity_", safe_rank, "_", safe_mode)
)

# -----------------------------
# Distance matrix
# -----------------------------

dist_obj <- vegan::vegdist(sample_taxa, method = distance_method)
dist_mat <- as.matrix(dist_obj)

dist_file <- paste0(out_prefix, "_distance_matrix.tsv")
write.table(
  data.frame(sample = rownames(dist_mat), dist_mat, check.names = FALSE),
  dist_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("[OK] Distance matrix saved:", dist_file, "\n")

# -----------------------------
# PCoA
# -----------------------------

pcoa <- cmdscale(dist_obj, k = 2, eig = TRUE)

positive_eig <- pcoa$eig[pcoa$eig > 0]
var_explained <- pcoa$eig / sum(positive_eig) * 100

pcoa_df <- data.frame(
  sample = rownames(pcoa$points),
  PCoA1 = pcoa$points[, 1],
  PCoA2 = pcoa$points[, 2],
  stringsAsFactors = FALSE
)

pcoa_df <- merge(pcoa_df, meta_df, by = "sample", all.x = TRUE)

pcoa_file <- paste0(out_prefix, "_pcoa_scores.tsv")
write.table(pcoa_df, pcoa_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat("[OK] PCoA scores saved:", pcoa_file, "\n")

p_pcoa <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2)) +
  geom_point(aes(shape = group), size = 3, alpha = 0.9) +
  labs(
    title = paste0("PCoA - ", rank_label, " level"),
    subtitle = paste0(
      "Distance: ", distance_method,
      " | Mode: ", mode,
      " | Transform: ", transform_method
    ),
    x = paste0("PCoA1 (", round(var_explained[1], 1), "%)"),
    y = paste0("PCoA2 (", round(var_explained[2], 1), "%)"),
    shape = "Group"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

if (ellipse == "yes" && length(unique(pcoa_df$group)) > 1) {
  group_counts <- table(pcoa_df$group)
  if (all(group_counts >= 3)) {
    p_pcoa <- p_pcoa + stat_ellipse(aes(group = group), linetype = 2)
  }
}

if (label_samples == "yes") {
  if (has_ggrepel) {
    p_pcoa <- p_pcoa + ggrepel::geom_text_repel(aes(label = sample), size = 3, max.overlaps = Inf)
  } else {
    p_pcoa <- p_pcoa + geom_text(aes(label = sample), size = 3, vjust = -0.8, check_overlap = TRUE)
  }
}

pcoa_png <- paste0(out_prefix, "_pcoa_", safe_distance, ".png")
ggsave(pcoa_png, p_pcoa, width = plot_width, height = plot_height, dpi = plot_dpi)
cat("[OK] PCoA PNG saved:", pcoa_png, "\n")

if (save_pdf) {
  pcoa_pdf <- paste0(out_prefix, "_pcoa_", safe_distance, ".pdf")
  ggsave(pcoa_pdf, p_pcoa, width = plot_width, height = plot_height)
  cat("[OK] PCoA PDF saved:", pcoa_pdf, "\n")
}

# -----------------------------
# NMDS
# -----------------------------

if (nrow(sample_taxa) >= 3) {
  set.seed(123)

  nmds <- tryCatch(
    vegan::metaMDS(
      sample_taxa,
      distance = distance_method,
      k = 2,
      trymax = 100,
      autotransform = FALSE,
      trace = FALSE
    ),
    error = function(e) e
  )

  if (inherits(nmds, "error")) {
    cat("[WARNING] NMDS failed:\n")
    cat(nmds$message, "\n")
  } else {
    nmds_points <- as.data.frame(scores(nmds, display = "sites"))
    nmds_points$sample <- rownames(nmds_points)

    if (!"NMDS1" %in% colnames(nmds_points)) {
      colnames(nmds_points)[1:2] <- c("NMDS1", "NMDS2")
    }

    nmds_df <- merge(nmds_points, meta_df, by = "sample", all.x = TRUE)

    nmds_file <- paste0(out_prefix, "_nmds_scores.tsv")
    write.table(nmds_df, nmds_file, sep = "\t", quote = FALSE, row.names = FALSE)
    cat("[OK] NMDS scores saved:", nmds_file, "\n")

    p_nmds <- ggplot(nmds_df, aes(x = NMDS1, y = NMDS2)) +
      geom_point(aes(shape = group), size = 3, alpha = 0.9) +
      labs(
        title = paste0("NMDS - ", rank_label, " level"),
        subtitle = paste0(
          "Distance: ", distance_method,
          " | Mode: ", mode,
          " | Stress: ", round(nmds$stress, 4)
        ),
        x = "NMDS1",
        y = "NMDS2",
        shape = "Group"
      ) +
      theme_bw(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank()
      )

    if (ellipse == "yes" && length(unique(nmds_df$group)) > 1) {
      group_counts <- table(nmds_df$group)
      if (all(group_counts >= 3)) {
        p_nmds <- p_nmds + stat_ellipse(aes(group = group), linetype = 2)
      }
    }

    if (label_samples == "yes") {
      if (has_ggrepel) {
        p_nmds <- p_nmds + ggrepel::geom_text_repel(aes(label = sample), size = 3, max.overlaps = Inf)
      } else {
        p_nmds <- p_nmds + geom_text(aes(label = sample), size = 3, vjust = -0.8, check_overlap = TRUE)
      }
    }

    nmds_png <- paste0(out_prefix, "_nmds_", safe_distance, ".png")
    ggsave(nmds_png, p_nmds, width = plot_width, height = plot_height, dpi = plot_dpi)
    cat("[OK] NMDS PNG saved:", nmds_png, "\n")

    if (save_pdf) {
      nmds_pdf <- paste0(out_prefix, "_nmds_", safe_distance, ".pdf")
      ggsave(nmds_pdf, p_nmds, width = plot_width, height = plot_height)
      cat("[OK] NMDS PDF saved:", nmds_pdf, "\n")
    }
  }
} else {
  cat("[WARNING] NMDS skipped because fewer than 3 samples are available.\n")
}

# -----------------------------
# PERMANOVA
# -----------------------------

if (length(unique(meta_df$group)) > 1) {
  permanova <- tryCatch(
    vegan::adonis2(dist_obj ~ group, data = meta_df, permutations = 999),
    error = function(e) e
  )

  if (inherits(permanova, "error")) {
    cat("[WARNING] PERMANOVA failed:\n")
    cat(permanova$message, "\n")
  } else {
    permanova_df <- as.data.frame(permanova)
    permanova_df$term <- rownames(permanova_df)
    permanova_df <- permanova_df[, c("term", setdiff(colnames(permanova_df), "term"))]

    permanova_file <- paste0(out_prefix, "_permanova.tsv")
    write.table(permanova_df, permanova_file, sep = "\t", quote = FALSE, row.names = FALSE)
    cat("[OK] PERMANOVA saved:", permanova_file, "\n")
  }
} else {
  cat("[INFO] PERMANOVA skipped because only one group is available.\n")
}

cat("============================================================\n")
cat("[DONE] Beta diversity completed.\n")
cat("Samples:", nrow(sample_taxa), "\n")
cat("Taxa   :", ncol(sample_taxa), "\n")
cat("Output :", output_dir, "\n")
cat("============================================================\n")
