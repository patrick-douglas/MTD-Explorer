#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

args2 <- c()
for (a in args) {
  if (grepl("^--[^=]+=", a)) {
    args2 <- c(args2, sub("=.*$", "", a), sub("^--[^=]+=", "", a))
  } else {
    args2 <- c(args2, a)
  }
}
args <- args2

get_arg <- function(flag, default = NA_character_) {
  idx <- which(args == flag)
  if (length(idx) == 0) return(default)
  if (idx[1] == length(args)) return(default)
  args[idx[1] + 1]
}

input <- get_arg("--input")
samplesheet <- get_arg("--samplesheet", "")

rank_label <- get_arg("--rank", "Taxon")
top_n <- as.integer(get_arg("--top", "12"))
mode <- tolower(get_arg("--mode", "relative"))
title <- get_arg("--title", paste0(rank_label, "-level taxonomic profile"))

output_prefix_user <- get_arg("--output_prefix", NA_character_)
output_dir <- get_arg("--output_dir", ".")

width_arg <- get_arg("--width", "auto")
height_arg <- get_arg("--height", "auto")
dpi <- as.integer(get_arg("--dpi", "300"))

bar_width <- as.numeric(get_arg("--bar_width", "0.52"))

legend_position_user <- tolower(get_arg("--legend_position", "auto"))
legend_ncol_arg <- get_arg("--legend_ncol", "auto")
legend_max_rows <- as.integer(get_arg("--legend_max_rows", "12"))

legend_text_size_user <- get_arg("--legend_text_size", "auto")
legend_title_size_user <- get_arg("--legend_title_size", "auto")
legend_key_size_user <- get_arg("--legend_key_size", "auto")

axis_text_size_default <- get_arg("--axis_text_size", "9")
x_text_size <- as.numeric(get_arg("--x_text_size", axis_text_size_default))
y_text_size <- as.numeric(get_arg("--y_text_size", axis_text_size_default))
axis_title_size <- as.numeric(get_arg("--axis_title_size", "11"))
plot_title_size <- as.numeric(get_arg("--plot_title_size", "13"))
x_angle <- as.numeric(get_arg("--x_angle", "45"))

taxon_format <- tolower(get_arg("--taxon_format", "scientific"))

if (!taxon_format %in% c("scientific", "plain", "none")) {
  stop("--taxon_format must be one of: scientific, plain, none")
}

if (is.na(input) || !file.exists(input)) {
  stop("Missing or invalid --input file.")
}

if (!mode %in% c("relative", "absolute", "fraction")) {
  stop("--mode must be one of: relative, absolute, fraction")
}

if (!legend_position_user %in% c("auto", "right", "bottom", "none")) {
  stop("--legend_position must be one of: auto, right, bottom, none")
}

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package ggplot2 is required.")
}

library(ggplot2)
library(grid)

make_safe_name <- function(x) {
  x <- as.character(x)
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

auto_name <- paste0(
  make_safe_name(rank_label),
  "_",
  make_safe_name(mode),
  "_top",
  top_n
)

if (is.na(output_prefix_user) || output_prefix_user == "" || output_prefix_user == "auto") {
  output_prefix <- file.path(output_dir, auto_name)
} else {
  output_prefix <- output_prefix_user
}

message("[INFO] Input: ", input)
message("[INFO] Rank: ", rank_label)
message("[INFO] Top taxa: ", top_n)
message("[INFO] Mode: ", mode)
message("[INFO] Taxon format: ", taxon_format)
message("[INFO] Output prefix: ", output_prefix)

read_table_auto <- function(path) {
  x <- tryCatch(
    read.delim(
      path,
      header = TRUE,
      sep = "\t",
      check.names = FALSE,
      quote = "",
      comment.char = "",
      stringsAsFactors = FALSE
    ),
    error = function(e) NULL
  )

  if (is.null(x) || ncol(x) < 2) {
    x <- read.csv(
      path,
      header = TRUE,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }

  x
}

clean_sample_name <- function(x) {
  x <- basename(x)

  x <- sub("^Report_", "", x)

  x <- sub("\\.phylum\\.bracken(_num|_frac)?$", "", x)
  x <- sub("\\.genus\\.bracken(_num|_frac)?$", "", x)
  x <- sub("\\.species\\.bracken(_num|_frac)?$", "", x)

  x <- sub("_phylum\\.bracken(_num|_frac)?$", "", x)
  x <- sub("_genus\\.bracken(_num|_frac)?$", "", x)
  x <- sub("_species\\.bracken(_num|_frac)?$", "", x)

  x <- sub("\\.bracken(_num|_frac)?$", "", x)
  x <- sub("_bracken(_num|_frac)?$", "", x)
  x <- sub("-bracken(_num|_frac)?$", "", x)

  x <- sub("(_num|_frac)$", "", x)
  x <- sub("(\\.num|\\.frac)$", "", x)

  x
}

clean_taxon_name <- function(x) {
  x <- as.character(x)

  if (taxon_format == "none") {
    x[x == "" | is.na(x)] <- "Unassigned"
    return(x)
  }

  x <- gsub("_", " ", x)
  x <- gsub("\\s+", " ", x)
  x <- trimws(x)
  x[x == "" | is.na(x)] <- "Unassigned"

  x
}

to_numeric_clean <- function(x) {
  x <- as.character(x)
  x <- gsub(",", "", x)
  x <- gsub("%", "", x)
  x <- trimws(x)
  suppressWarnings(as.numeric(x))
}

collapse_duplicate_sample_cols <- function(mat) {
  cn <- colnames(mat)
  u <- unique(cn)

  out <- lapply(u, function(s) {
    rowSums(mat[, cn == s, drop = FALSE], na.rm = TRUE)
  })

  out <- as.data.frame(out, check.names = FALSE, stringsAsFactors = FALSE)
  colnames(out) <- u

  out
}

cap_first_lower_rest <- function(x) {
  x <- as.character(x)
  x <- trimws(x)

  if (is.na(x) || x == "") return(x)

  paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))))
}

is_placeholder_taxon <- function(x) {
  lx <- tolower(trimws(x))
  lx %in% c(
    "other", "unassigned", "unclassified", "unknown", "root",
    "na", "nan", "none", ""
  )
}

format_genus_plain <- function(label) {
  label <- clean_taxon_name(label)

  if (is_placeholder_taxon(label)) return(label)

  words <- unlist(strsplit(label, "\\s+"))
  if (length(words) == 0) return(label)

  if (tolower(words[1]) == "candidatus" && length(words) >= 2) {
    genus <- cap_first_lower_rest(words[2])
    rest <- if (length(words) > 2) paste(words[3:length(words)], collapse = " ") else ""
    return(trimws(paste("Candidatus", genus, rest)))
  }

  if (tolower(words[1]) %in% c("uncultured", "unclassified", "candidate", "environmental") && length(words) >= 2) {
    prefix <- tolower(words[1])
    genus <- cap_first_lower_rest(words[2])
    rest <- if (length(words) > 2) paste(words[3:length(words)], collapse = " ") else ""
    return(trimws(paste(prefix, genus, rest)))
  }

  genus <- cap_first_lower_rest(words[1])
  rest <- if (length(words) > 1) paste(words[2:length(words)], collapse = " ") else ""

  trimws(paste(genus, rest))
}

format_species_plain <- function(label) {
  label <- clean_taxon_name(label)

  if (is_placeholder_taxon(label)) return(label)

  words <- unlist(strsplit(label, "\\s+"))
  if (length(words) == 0) return(label)

  if (tolower(words[1]) == "candidatus" && length(words) >= 3) {
    genus <- cap_first_lower_rest(words[2])
    species <- tolower(words[3])
    rest <- if (length(words) > 3) paste(words[4:length(words)], collapse = " ") else ""
    return(trimws(paste("Candidatus", genus, species, rest)))
  }

  if (tolower(words[1]) %in% c("uncultured", "unclassified", "candidate", "environmental") && length(words) >= 2) {
    prefix <- tolower(words[1])
    genus <- cap_first_lower_rest(words[2])

    if (length(words) >= 3 && tolower(words[3]) %in% c("sp.", "sp", "spp.", "spp")) {
      rest <- if (length(words) > 3) paste(words[4:length(words)], collapse = " ") else ""
      return(trimws(paste(prefix, genus, "sp.", rest)))
    }

    if (length(words) >= 3) {
      species <- tolower(words[3])
      rest <- if (length(words) > 3) paste(words[4:length(words)], collapse = " ") else ""
      return(trimws(paste(prefix, genus, species, rest)))
    }

    return(trimws(paste(prefix, genus)))
  }

  genus <- cap_first_lower_rest(words[1])

  if (length(words) >= 2 && tolower(words[2]) %in% c("sp.", "sp", "spp.", "spp")) {
    rest <- if (length(words) > 2) paste(words[3:length(words)], collapse = " ") else ""
    return(trimws(paste(genus, "sp.", rest)))
  }

  if (length(words) >= 2) {
    species <- tolower(words[2])
    rest <- if (length(words) > 2) paste(words[3:length(words)], collapse = " ") else ""
    return(trimws(paste(genus, species, rest)))
  }

  genus
}

format_high_rank_plain <- function(label) {
  label <- clean_taxon_name(label)

  if (is_placeholder_taxon(label)) return(label)

  words <- unlist(strsplit(label, "\\s+"))

  words <- sapply(words, function(w) {
    if (nchar(w) <= 1) return(toupper(w))
    cap_first_lower_rest(w)
  })

  paste(words, collapse = " ")
}

format_taxon_plain <- function(label, rank_label) {
  r <- tolower(rank_label)

  if (taxon_format == "none") return(label)
  if (taxon_format == "plain") return(clean_taxon_name(label))

  if (grepl("species", r)) return(format_species_plain(label))
  if (grepl("genus", r)) return(format_genus_plain(label))

  format_high_rank_plain(label)
}

pm_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("\\\\", "\\\\\\\\", x)
  x <- gsub('"', '\\"', x)
  x
}

pm_quote <- function(x) {
  paste0('"', pm_escape(x), '"')
}

pm_italic <- function(x) {
  paste0("italic(", pm_quote(x), ")")
}

pm_paste <- function(parts) {
  parts <- parts[!is.na(parts) & parts != ""]
  if (length(parts) == 1) return(parts)
  paste0("paste(", paste(parts, collapse = ", "), ")")
}

format_genus_plotmath <- function(label) {
  plain <- format_genus_plain(label)

  if (is_placeholder_taxon(plain)) return(pm_quote(plain))

  words <- unlist(strsplit(plain, "\\s+"))
  if (length(words) == 0) return(pm_quote(plain))

  if (tolower(words[1]) == "candidatus" && length(words) >= 2) {
    genus <- words[2]
    rest <- if (length(words) > 2) paste(words[3:length(words)], collapse = " ") else ""

    return(pm_paste(c(
      pm_quote("Candidatus "),
      pm_italic(genus),
      if (rest != "") pm_quote(paste0(" ", rest)) else ""
    )))
  }

  if (tolower(words[1]) %in% c("uncultured", "unclassified", "candidate", "environmental") && length(words) >= 2) {
    prefix <- words[1]
    genus <- words[2]
    rest <- if (length(words) > 2) paste(words[3:length(words)], collapse = " ") else ""

    return(pm_paste(c(
      pm_quote(paste0(prefix, " ")),
      pm_italic(genus),
      if (rest != "") pm_quote(paste0(" ", rest)) else ""
    )))
  }

  genus <- words[1]
  rest <- if (length(words) > 1) paste(words[2:length(words)], collapse = " ") else ""

  pm_paste(c(
    pm_italic(genus),
    if (rest != "") pm_quote(paste0(" ", rest)) else ""
  ))
}

format_species_plotmath <- function(label) {
  plain <- format_species_plain(label)

  if (is_placeholder_taxon(plain)) return(pm_quote(plain))

  words <- unlist(strsplit(plain, "\\s+"))
  if (length(words) == 0) return(pm_quote(plain))

  if (tolower(words[1]) == "candidatus" && length(words) >= 3) {
    genus_species <- paste(words[2], words[3])
    rest <- if (length(words) > 3) paste(words[4:length(words)], collapse = " ") else ""

    return(pm_paste(c(
      pm_quote("Candidatus "),
      pm_italic(genus_species),
      if (rest != "") pm_quote(paste0(" ", rest)) else ""
    )))
  }

  if (tolower(words[1]) %in% c("uncultured", "unclassified", "candidate", "environmental") && length(words) >= 2) {
    prefix <- words[1]
    genus <- words[2]

    if (length(words) >= 3 && tolower(words[3]) %in% c("sp.", "sp", "spp.", "spp")) {
      rest <- if (length(words) > 3) paste(words[4:length(words)], collapse = " ") else ""

      return(pm_paste(c(
        pm_quote(paste0(prefix, " ")),
        pm_italic(genus),
        pm_quote(" sp."),
        if (rest != "") pm_quote(paste0(" ", rest)) else ""
      )))
    }

    if (length(words) >= 3) {
      genus_species <- paste(genus, words[3])
      rest <- if (length(words) > 3) paste(words[4:length(words)], collapse = " ") else ""

      return(pm_paste(c(
        pm_quote(paste0(prefix, " ")),
        pm_italic(genus_species),
        if (rest != "") pm_quote(paste0(" ", rest)) else ""
      )))
    }

    return(pm_paste(c(
      pm_quote(paste0(prefix, " ")),
      pm_italic(genus)
    )))
  }

  genus <- words[1]

  if (length(words) >= 2 && tolower(words[2]) %in% c("sp.", "sp", "spp.", "spp")) {
    rest <- if (length(words) > 2) paste(words[3:length(words)], collapse = " ") else ""

    return(pm_paste(c(
      pm_italic(genus),
      pm_quote(" sp."),
      if (rest != "") pm_quote(paste0(" ", rest)) else ""
    )))
  }

  if (length(words) >= 2) {
    genus_species <- paste(genus, words[2])
    rest <- if (length(words) > 2) paste(words[3:length(words)], collapse = " ") else ""

    return(pm_paste(c(
      pm_italic(genus_species),
      if (rest != "") pm_quote(paste0(" ", rest)) else ""
    )))
  }

  pm_italic(genus)
}

format_taxon_plotmath <- function(label, rank_label) {
  r <- tolower(rank_label)

  if (taxon_format == "none") return(pm_quote(label))
  if (taxon_format == "plain") return(pm_quote(clean_taxon_name(label)))

  if (grepl("species", r)) return(format_species_plotmath(label))
  if (grepl("genus", r)) return(format_genus_plotmath(label))

  pm_quote(format_high_rank_plain(label))
}

df <- read_table_auto(input)

if (nrow(df) == 0 || ncol(df) < 2) {
  stop("Input table appears empty or has fewer than 2 columns.")
}

lower_names <- tolower(names(df))

taxon_col <- NULL

for (candidate in c("name", "taxon", "taxonomy", "classification", "#classification")) {
  hit <- which(lower_names == candidate)

  if (length(hit) > 0) {
    taxon_col <- names(df)[hit[1]]
    break
  }
}

if (is.null(taxon_col)) {
  taxon_col <- names(df)[1]
}

message("[INFO] Taxon column detected: ", taxon_col)

candidate_cols <- setdiff(names(df), taxon_col)
candidate_lower <- tolower(candidate_cols)

num_cols <- candidate_cols[grepl("(_num|\\.num)$", candidate_lower)]
frac_cols <- candidate_cols[grepl("(_frac|\\.frac)$", candidate_lower)]

if (mode %in% c("relative", "absolute")) {
  value_cols <- num_cols
  value_source <- "Bracken count columns ending in _num"
} else {
  value_cols <- frac_cols
  value_source <- "Bracken fraction columns ending in _frac"
}

if (length(value_cols) == 0) {
  metadata_patterns <- paste(
    c(
      "name", "taxon", "taxonomy", "classification", "taxid", "tax_id",
      "taxonomy_id", "taxonomy_lvl", "level", "rank", "kraken",
      "assigned", "added", "new_est", "fraction", "percent", "percentage"
    ),
    collapse = "|"
  )

  value_cols <- c()

  for (cn in candidate_cols) {
    lname <- tolower(cn)

    if (grepl(metadata_patterns, lname)) {
      next
    }

    vals <- to_numeric_clean(df[[cn]])
    ok <- sum(!is.na(vals))

    if (ok > 0 && ok >= 0.5 * nrow(df)) {
      value_cols <- c(value_cols, cn)
    }
  }

  value_source <- "fallback numeric sample columns"
}

if (length(value_cols) == 0) {
  stop(
    "Could not detect abundance columns. Check your input table format.\n",
    "Detected columns were:\n",
    paste(names(df), collapse = ", ")
  )
}

message("[INFO] Value source: ", value_source)
message("[INFO] Value columns detected:")
message(paste("  -", value_cols, collapse = "\n"))

taxa_raw <- clean_taxon_name(df[[taxon_col]])
taxa_plain <- vapply(taxa_raw, format_taxon_plain, character(1), rank_label = rank_label)

mat <- as.data.frame(
  lapply(df[value_cols], to_numeric_clean),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

mat[is.na(mat)] <- 0
mat[mat < 0] <- 0

colnames(mat) <- clean_sample_name(colnames(mat))
mat <- collapse_duplicate_sample_cols(mat)

mat$Taxon <- taxa_plain
mat <- aggregate(. ~ Taxon, data = mat, FUN = sum)

taxa <- mat$Taxon
mat$Taxon <- NULL
rownames(mat) <- taxa

total_by_taxon <- rowSums(mat)
keep <- total_by_taxon > 0
mat <- mat[keep, , drop = FALSE]

if (nrow(mat) == 0) {
  stop("No taxa remained after filtering.")
}

if (!is.na(samplesheet) && samplesheet != "" && file.exists(samplesheet)) {
  ss <- read.csv(samplesheet, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

  if (ncol(ss) >= 1) {
    wanted_order <- as.character(ss[[1]])
    matched <- wanted_order[wanted_order %in% colnames(mat)]
    remaining <- setdiff(colnames(mat), matched)
    mat <- mat[, c(matched, remaining), drop = FALSE]
  }
}

if (mode == "relative") {
  sample_sums <- colSums(mat)
  sample_sums[sample_sums == 0] <- 1
  plot_mat <- sweep(mat, 2, sample_sums, "/") * 100
  y_lab <- "Relative abundance (%)"
  y_limits <- c(0, 100)
  y_breaks <- seq(0, 100, 25)
} else if (mode == "fraction") {
  plot_mat <- mat

  if (max(plot_mat, na.rm = TRUE) <= 1) {
    plot_mat <- plot_mat * 100
  }

  y_lab <- "Relative abundance (%)"
  y_limits <- c(0, 100)
  y_breaks <- seq(0, 100, 25)
} else {
  plot_mat <- mat
  y_lab <- "Estimated reads assigned by Bracken"
  y_limits <- NULL
  y_breaks <- waiver()
}

taxon_totals <- rowSums(mat)
top_taxa <- names(sort(taxon_totals, decreasing = TRUE))[seq_len(min(top_n, length(taxon_totals)))]

plot_top <- plot_mat[top_taxa, , drop = FALSE]

if (nrow(plot_mat) > length(top_taxa)) {
  other_values <- colSums(plot_mat[setdiff(rownames(plot_mat), top_taxa), , drop = FALSE])
  plot_final <- rbind(plot_top, Other = other_values)
} else {
  plot_final <- plot_top
}

plot_df <- data.frame()

for (s in colnames(plot_final)) {
  tmp <- data.frame(
    Sample = s,
    Taxon = rownames(plot_final),
    Abundance = as.numeric(plot_final[, s]),
    stringsAsFactors = FALSE
  )

  plot_df <- rbind(plot_df, tmp)
}

taxon_levels <- rownames(plot_final)
taxon_levels <- c(setdiff(taxon_levels, "Other"), intersect(taxon_levels, "Other"))

plot_df$Sample <- factor(plot_df$Sample, levels = colnames(plot_final))
plot_df$Taxon <- factor(plot_df$Taxon, levels = taxon_levels)

n_taxa <- length(levels(plot_df$Taxon))
n_samples <- length(levels(plot_df$Sample))

legend_label_text <- vapply(
  levels(plot_df$Taxon),
  format_taxon_plotmath,
  character(1),
  rank_label = rank_label
)

legend_labels <- parse(text = legend_label_text)
names(legend_labels) <- levels(plot_df$Taxon)

if (legend_position_user == "auto") {
  if (n_taxa <= 30) {
    legend_position <- "right"
  } else {
    legend_position <- "bottom"
  }
} else {
  legend_position <- legend_position_user
}

if (legend_ncol_arg == "auto") {
  if (legend_position == "right") {
    legend_ncol <- max(1, ceiling(n_taxa / legend_max_rows))
  } else if (legend_position == "bottom") {
    legend_ncol <- max(1, ceiling(n_taxa / legend_max_rows))
    legend_ncol <- min(legend_ncol, 10)
  } else {
    legend_ncol <- 1
  }
} else {
  legend_ncol <- as.integer(legend_ncol_arg)
}

if (is.na(legend_ncol) || legend_ncol < 1) {
  legend_ncol <- 1
}

legend_rows <- ceiling(n_taxa / legend_ncol)

legend_text_size_auto <- if (n_taxa <= 18) {
  8.2
} else if (n_taxa <= 35) {
  7.4
} else {
  6.8
}

legend_title_size_auto <- if (n_taxa <= 18) {
  10.5
} else if (n_taxa <= 35) {
  10
} else {
  9.5
}

legend_key_size_auto <- if (n_taxa <= 18) {
  0.42
} else if (n_taxa <= 35) {
  0.36
} else {
  0.32
}

legend_text_size <- if (legend_text_size_user == "auto") {
  legend_text_size_auto
} else {
  as.numeric(legend_text_size_user)
}

legend_title_size <- if (legend_title_size_user == "auto") {
  legend_title_size_auto
} else {
  as.numeric(legend_title_size_user)
}

legend_key_size <- if (legend_key_size_user == "auto") {
  legend_key_size_auto
} else {
  as.numeric(legend_key_size_user)
}

if (is.na(legend_text_size)) legend_text_size <- legend_text_size_auto
if (is.na(legend_title_size)) legend_title_size <- legend_title_size_auto
if (is.na(legend_key_size)) legend_key_size <- legend_key_size_auto

if (tolower(width_arg) == "auto") {
  if (legend_position == "right") {
    width <- max(9, 4.5 + n_samples * 0.65 + legend_ncol * 2.9)
  } else if (legend_position == "bottom") {
    width <- max(12, 5.0 + n_samples * 0.75 + legend_ncol * 1.15)
  } else {
    width <- max(8, 4.5 + n_samples * 0.7)
  }
} else {
  width <- as.numeric(width_arg)
}

if (tolower(height_arg) == "auto") {
  if (legend_position == "right") {
    height <- max(5.5, 3.2 + legend_rows * 0.30 * (legend_text_size / 7))
  } else if (legend_position == "bottom") {
    height <- max(6.2, 4.9 + legend_rows * 0.34 * (legend_text_size / 7))
  } else {
    height <- 5.5
  }
} else {
  height <- as.numeric(height_arg)
}

message("[INFO] Legend position: ", legend_position)
message("[INFO] Legend columns: ", legend_ncol)
message("[INFO] Legend rows: ", legend_rows)
message("[INFO] Legend text size: ", legend_text_size)
message("[INFO] Legend title size: ", legend_title_size)
message("[INFO] Legend key size: ", legend_key_size)
message("[INFO] X-axis text size: ", x_text_size)
message("[INFO] Y-axis text size: ", y_text_size)
message("[INFO] Figure width: ", width)
message("[INFO] Figure height: ", height)

n_taxa_palette <- length(levels(plot_df$Taxon))
pal <- grDevices::hcl.colors(n_taxa_palette, palette = "Dark 3")
names(pal) <- levels(plot_df$Taxon)

if ("Other" %in% names(pal)) {
  pal["Other"] <- "grey82"
}

p <- ggplot(plot_df, aes(x = Sample, y = Abundance, fill = Taxon)) +
  geom_bar(stat = "identity", width = bar_width, color = "grey25", linewidth = 0.10) +
  scale_x_discrete(expand = expansion(mult = c(0.08, 0.08))) +
  scale_fill_manual(
    values = pal,
    breaks = levels(plot_df$Taxon),
    labels = legend_labels,
    drop = FALSE
  ) +
  labs(
    title = title,
    x = NULL,
    y = y_lab,
    fill = rank_label
  ) +
guides(
  fill = guide_legend(
    ncol = legend_ncol,
    byrow = TRUE,
    title.position = "top",
    title.hjust = 0.5,
    label.position = "right",
    label.hjust = 0,
    label.vjust = 0.5,
    keyheight = unit(legend_key_size, "cm"),
    keywidth = unit(legend_key_size, "cm")
  )
) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = plot_title_size, hjust = 0),
    axis.text.x = element_text(angle = x_angle, hjust = 1, vjust = 1, size = x_text_size),
    axis.text.y = element_text(size = y_text_size),
    axis.title.y = element_text(face = "bold", size = axis_title_size),
    legend.title = element_text(face = "bold", size = legend_title_size),
    legend.text = element_text(size = legend_text_size, hjust = 0),
    legend.position = legend_position,
    legend.box = "vertical",
    legend.box.just = "left",
    legend.justification = "left",
    legend.spacing.y = unit(0.05, "cm"),
    legend.spacing.x = unit(0.15, "cm"),
    legend.margin = margin(2, 2, 2, 2),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 12, 10, 10)
  )

if (mode %in% c("relative", "fraction")) {
  p <- p +
    scale_y_continuous(
      limits = y_limits,
      breaks = y_breaks,
      expand = c(0, 0)
    )
} else {
  p <- p +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.05)),
      labels = function(x) format(round(x), big.mark = ",", scientific = FALSE)
    )
}

if (legend_position == "none") {
  p <- p + theme(legend.position = "none")
}

png_out <- paste0(output_prefix, ".png")
pdf_out <- paste0(output_prefix, ".pdf")
csv_out <- paste0(output_prefix, "_plot_table.csv")
matrix_out <- paste0(output_prefix, "_matrix.csv")

out_dir <- dirname(output_prefix)

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

ggsave(png_out, p, width = width, height = height, dpi = dpi, bg = "white")
ggsave(pdf_out, p, width = width, height = height, bg = "white")

write.csv(plot_df, csv_out, row.names = FALSE)
write.csv(plot_final, matrix_out, row.names = TRUE)

message("[OK] Saved:")
message("  ", png_out)
message("  ", pdf_out)
message("  ", csv_out)
message("  ", matrix_out)
