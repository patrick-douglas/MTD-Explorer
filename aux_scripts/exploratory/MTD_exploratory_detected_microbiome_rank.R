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
