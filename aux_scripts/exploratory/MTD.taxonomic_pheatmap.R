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
rank_label <- get_arg("--rank", "Species")
top_n <- as.integer(get_arg("--top", "30"))
mode <- tolower(get_arg("--mode", "relative"))
transform <- tolower(get_arg("--transform", "log10"))

cluster_samples <- tolower(get_arg("--cluster_samples", "yes"))
cluster_taxa <- tolower(get_arg("--cluster_taxa", "no"))
cluster_distance <- tolower(get_arg("--cluster_distance", "correlation"))
cluster_method <- get_arg("--cluster_method", "complete")

output_dir <- get_arg("--output_dir", ".")
output_prefix <- get_arg("--output_prefix", "auto")

width <- as.numeric(get_arg("--width", "8"))
height <- as.numeric(get_arg("--height", "8"))
fontsize_row <- as.numeric(get_arg("--fontsize_row", "8"))
fontsize_col <- as.numeric(get_arg("--fontsize_col", "10"))
fontsize <- as.numeric(get_arg("--fontsize", "10"))

title <- get_arg("--title", paste0(rank_label, "-level abundance heatmap"))

if (is.na(input) || !file.exists(input)) {
  stop("Missing or invalid --input file.")
}

if (!requireNamespace("pheatmap", quietly = TRUE)) {
  stop("Package pheatmap is required. Install with: install.packages('pheatmap')")
}

library(pheatmap)

make_safe_name <- function(x) {
  x <- as.character(x)
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

if (output_prefix == "auto" || output_prefix == "") {
  output_prefix <- file.path(
    output_dir,
    paste0(
      make_safe_name(rank_label),
      "_pheatmap_",
      make_safe_name(mode),
      "_",
      make_safe_name(transform),
      "_top",
      top_n,
      ifelse(cluster_samples == "yes", "_sampleDendrogram", ""),
      ifelse(cluster_taxa == "yes", "_taxaDendrogram", "")
    )
  )
}

dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)

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
    x <- read.csv(path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
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
  x <- gsub("_", " ", x)
  x <- gsub("\\s+", " ", x)
  trimws(x)
}

cap_first_lower_rest <- function(x) {
  paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))))
}

format_taxon <- function(x, rank_label) {
  x <- clean_taxon_name(x)
  r <- tolower(rank_label)

  if (tolower(x) %in% c("other", "unknown", "unclassified", "unassigned")) {
    return(x)
  }

  words <- unlist(strsplit(x, "\\s+"))

  if (length(words) == 0) return(x)

  if (grepl("species", r)) {
    if (tolower(words[1]) == "candidatus" && length(words) >= 3) {
      genus <- cap_first_lower_rest(words[2])
      species <- tolower(words[3])
      rest <- if (length(words) > 3) paste(words[4:length(words)], collapse = " ") else ""
      return(trimws(paste("Candidatus", genus, species, rest)))
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

    return(genus)
  }

  if (grepl("genus", r)) {
    return(cap_first_lower_rest(words[1]))
  }

  paste(sapply(words, cap_first_lower_rest), collapse = " ")
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

make_distance <- function(x, target = "samples") {
  x <- as.matrix(x)
  x[!is.finite(x)] <- 0

  if (target == "samples") {
    if (cluster_distance == "euclidean") return(dist(t(x)))

    cm <- suppressWarnings(cor(x, use = "pairwise.complete.obs"))
    dm <- 1 - cm
    dm[!is.finite(dm)] <- 1
    diag(dm) <- 0
    return(as.dist(dm))
  }

  if (cluster_distance == "euclidean") return(dist(x))

  cm <- suppressWarnings(cor(t(x), use = "pairwise.complete.obs"))
  dm <- 1 - cm
  dm[!is.finite(dm)] <- 1
  diag(dm) <- 0
  as.dist(dm)
}

df <- read_table_auto(input)

lower_names <- tolower(names(df))

taxon_col <- NULL
for (candidate in c("name", "taxon", "taxonomy", "classification", "#classification")) {
  hit <- which(lower_names == candidate)
  if (length(hit) > 0) {
    taxon_col <- names(df)[hit[1]]
    break
  }
}

if (is.null(taxon_col)) taxon_col <- names(df)[1]

candidate_cols <- setdiff(names(df), taxon_col)
candidate_lower <- tolower(candidate_cols)

num_cols <- candidate_cols[grepl("(_num|\\.num)$", candidate_lower)]
frac_cols <- candidate_cols[grepl("(_frac|\\.frac)$", candidate_lower)]

if (mode %in% c("relative", "absolute")) {
  value_cols <- num_cols
} else {
  value_cols <- frac_cols
}

if (length(value_cols) == 0) {
  stop("Could not detect abundance columns ending in _num or _frac.")
}

taxa <- vapply(df[[taxon_col]], format_taxon, character(1), rank_label = rank_label)

mat <- as.data.frame(lapply(df[value_cols], to_numeric_clean), check.names = FALSE)
mat[is.na(mat)] <- 0
mat[mat < 0] <- 0

colnames(mat) <- clean_sample_name(colnames(mat))
mat <- collapse_duplicate_sample_cols(mat)

mat$Taxon <- taxa
mat <- aggregate(. ~ Taxon, data = mat, FUN = sum)

rownames(mat) <- mat$Taxon
mat$Taxon <- NULL

mat <- mat[rowSums(mat) > 0, , drop = FALSE]

if (!is.na(samplesheet) && samplesheet != "" && file.exists(samplesheet)) {
  ss <- read.csv(samplesheet, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  wanted_order <- as.character(ss[[1]])
  matched <- wanted_order[wanted_order %in% colnames(mat)]
  remaining <- setdiff(colnames(mat), matched)
  mat <- mat[, c(matched, remaining), drop = FALSE]
}

taxon_totals <- rowSums(mat)
top_taxa <- names(sort(taxon_totals, decreasing = TRUE))[seq_len(min(top_n, length(taxon_totals)))]
mat <- mat[top_taxa, , drop = FALSE]

if (mode == "relative") {
  sample_sums <- colSums(mat)
  sample_sums[sample_sums == 0] <- 1
  mat_plot <- sweep(mat, 2, sample_sums, "/") * 100
  legend_title <- "Relative abundance (%)"
} else if (mode == "fraction") {
  mat_plot <- mat
  if (max(mat_plot, na.rm = TRUE) <= 1) mat_plot <- mat_plot * 100
  legend_title <- "Relative abundance (%)"
} else {
  mat_plot <- mat
  legend_title <- "Estimated reads"
}

if (transform == "log10") {
  pseudocount <- if (mode == "absolute") 1 else 0.001
  mat_plot <- log10(mat_plot + pseudocount)
  legend_title <- paste0("log10(", legend_title, " + ", pseudocount, ")")
} else if (transform == "zscore") {
  pseudocount <- if (mode == "absolute") 1 else 0.001
  tmp <- log10(mat_plot + pseudocount)
  mat_plot <- t(apply(tmp, 1, function(z) {
    s <- sd(z, na.rm = TRUE)
    if (is.na(s) || s == 0) return(rep(0, length(z)))
    (z - mean(z, na.rm = TRUE)) / s
  }))
  colnames(mat_plot) <- colnames(mat)
  rownames(mat_plot) <- rownames(mat)
  legend_title <- "Row z-score"
}

cluster_cols <- cluster_samples == "yes"
cluster_rows <- cluster_taxa == "yes"

clustering_distance_cols <- if (cluster_cols) make_distance(mat_plot, "samples") else "euclidean"
clustering_distance_rows <- if (cluster_rows) make_distance(mat_plot, "taxa") else "euclidean"

colors <- colorRampPalette(grDevices::hcl.colors(100, "YlOrRd", rev = FALSE))(100)

if (transform == "zscore") {
  colors <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)
}

png_out <- paste0(output_prefix, ".png")
pdf_out <- paste0(output_prefix, ".pdf")
matrix_out <- paste0(output_prefix, "_matrix.csv")

png(png_out, width = width, height = height, units = "in", res = 300)
pheatmap::pheatmap(
  mat_plot,
  color = colors,
  cluster_cols = cluster_cols,
  cluster_rows = cluster_rows,
  clustering_distance_cols = clustering_distance_cols,
  clustering_distance_rows = clustering_distance_rows,
  clustering_method = cluster_method,
  fontsize = fontsize,
  fontsize_row = fontsize_row,
  fontsize_col = fontsize_col,
  main = title,
  angle_col = 45,
  border_color = "grey85",
  legend = TRUE,
  legend_breaks = NA,
  silent = FALSE
)
dev.off()

pdf(pdf_out, width = width, height = height)
pheatmap::pheatmap(
  mat_plot,
  color = colors,
  cluster_cols = cluster_cols,
  cluster_rows = cluster_rows,
  clustering_distance_cols = clustering_distance_cols,
  clustering_distance_rows = clustering_distance_rows,
  clustering_method = cluster_method,
  fontsize = fontsize,
  fontsize_row = fontsize_row,
  fontsize_col = fontsize_col,
  main = title,
  angle_col = 45,
  border_color = "grey85",
  legend = TRUE,
  silent = FALSE
)
dev.off()

write.csv(mat_plot, matrix_out, row.names = TRUE)

message("[OK] Saved:")
message("  ", png_out)
message("  ", pdf_out)
message("  ", matrix_out)
