#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  ix <- which(args == flag)
  if (length(ix) == 0) return(default)
  if (ix[1] == length(args)) stop("Missing value for argument: ", flag)
  args[ix[1] + 1]
}

scores_file <- get_arg("--scores")
samplesheet_file <- get_arg("--samplesheet")
outdir <- get_arg("--outdir")
top_var <- as.integer(get_arg("--top-var", "50"))
top_diff <- as.integer(get_arg("--top-diff", "20"))
pca_top_var <- as.integer(get_arg("--pca-top-var", "500"))

if (is.null(scores_file) || is.null(samplesheet_file) || is.null(outdir)) {
  stop(
    "Usage:\n",
    "Rscript plot_ssgsea_results.R ",
    "--scores ssgsea-results-scores.gct ",
    "--samplesheet samplesheet.csv ",
    "--outdir ssGSEA/plots ",
    "[--top-var 50] [--top-diff 20] [--pca-top-var 500]\n"
  )
}

if (!file.exists(scores_file)) {
  stop("Scores file not found: ", scores_file)
}

if (!file.exists(samplesheet_file)) {
  stop("Samplesheet not found: ", samplesheet_file)
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

message("[INFO] ssGSEA scores: ", scores_file)
message("[INFO] samplesheet: ", samplesheet_file)
message("[INFO] output directory: ", outdir)

read_gct_scores <- function(file) {
  header_lines <- readLines(file, n = 2)

  if (!grepl("^#1\\.", header_lines[1])) {
    stop("Input file does not look like a GCT file: ", file)
  }

  dims <- strsplit(header_lines[2], "[\t ]+")[[1]]
  dims <- dims[dims != ""]

  if (length(dims) < 2) {
    stop("Could not parse GCT dimensions from second line: ", header_lines[2])
  }

  n_samples <- as.integer(dims[2])

  tab <- read.delim(
    file,
    skip = 2,
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )

  if (ncol(tab) < n_samples + 1) {
    stop("GCT file has fewer columns than expected.")
  }

  term_ids <- as.character(tab[[1]])
  sample_cols <- tail(colnames(tab), n_samples)

  mat <- as.matrix(tab[, sample_cols, drop = FALSE])
  suppressWarnings(storage.mode(mat) <- "numeric")

  rownames(mat) <- term_ids

  desc_col <- NULL
  if (ncol(tab) >= 2) {
    desc_col <- as.character(tab[[2]])
  } else {
    desc_col <- term_ids
  }

  desc <- data.frame(
    term = term_ids,
    description = desc_col,
    stringsAsFactors = FALSE
  )

  list(mat = mat, desc = desc, sample_cols = sample_cols)
}

read_samplesheet <- function(file, score_samples) {
  ss <- read.csv(
    file,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  if (ncol(ss) < 2) {
    stop("Samplesheet must have at least two columns: sample and group.")
  }

  cn <- colnames(ss)
  cn_lower <- tolower(cn)

  sample_col <- cn[1]

  group_candidates <- c("group", "condition", "treatment", "grupo")
  group_col <- NULL

  for (candidate in group_candidates) {
    hit <- which(cn_lower == candidate)
    if (length(hit) > 0) {
      group_col <- cn[hit[1]]
      break
    }
  }

  if (is.null(group_col)) {
    group_col <- cn[2]
  }

  meta <- data.frame(
    sample = as.character(ss[[sample_col]]),
    group = as.character(ss[[group_col]]),
    stringsAsFactors = FALSE
  )

  meta$sample <- trimws(meta$sample)
  meta$group <- trimws(meta$group)

  meta <- meta[meta$sample != "" & meta$group != "", , drop = FALSE]
  meta <- meta[!duplicated(meta$sample), , drop = FALSE]

  missing_in_scores <- setdiff(meta$sample, score_samples)
  missing_in_sheet <- setdiff(score_samples, meta$sample)

  if (length(missing_in_scores) > 0) {
    warning(
      "Samples present in samplesheet but absent from ssGSEA scores: ",
      paste(missing_in_scores, collapse = ", ")
    )
  }

  if (length(missing_in_sheet) > 0) {
    warning(
      "Samples present in ssGSEA scores but absent from samplesheet: ",
      paste(missing_in_sheet, collapse = ", ")
    )
  }

  meta <- meta[meta$sample %in% score_samples, , drop = FALSE]

  if (nrow(meta) < 2) {
    stop("Fewer than two samples matched between samplesheet and ssGSEA scores.")
  }

  meta$group <- factor(meta$group, levels = unique(meta$group))

  meta
}

clean_matrix <- function(mat) {
  keep <- rowSums(is.finite(mat)) >= 2
  mat <- mat[keep, , drop = FALSE]

  for (i in seq_len(nrow(mat))) {
    bad <- !is.finite(mat[i, ])
    if (any(bad)) {
      mat[i, bad] <- mean(mat[i, !bad], na.rm = TRUE)
    }
  }

  vars <- apply(mat, 1, var, na.rm = TRUE)
  mat <- mat[is.finite(vars) & vars > 0, , drop = FALSE]

  mat
}

row_zscore <- function(mat) {
  z <- t(scale(t(mat)))
  z[!is.finite(z)] <- 0
  z
}

save_plot <- function(plot, filename_base, width = 10, height = 8) {
  png_file <- paste0(filename_base, ".png")
  pdf_file <- paste0(filename_base, ".pdf")

  ggsave(png_file, plot = plot, width = width, height = height, dpi = 180)
  ggsave(pdf_file, plot = plot, width = width, height = height)

  message("[OK] Saved: ", png_file)
  message("[OK] Saved: ", pdf_file)
}

plot_heatmap <- function(mat, meta, outfile_base, title) {
  z <- row_zscore(mat)

  sample_order <- meta$sample[order(meta$group, meta$sample)]
  z <- z[, sample_order, drop = FALSE]

  term_order <- rownames(z)

  df <- as.data.frame(as.table(z), stringsAsFactors = FALSE)
  colnames(df) <- c("Term", "Sample", "Z_score")

  df$Sample <- factor(df$Sample, levels = sample_order)
  df$Term <- factor(df$Term, levels = rev(term_order))

  p <- ggplot(df, aes(x = Sample, y = Term, fill = Z_score)) +
    geom_tile() +
    theme_bw() +
    labs(
      title = title,
      x = "Sample",
      y = "GO term",
      fill = "Row z-score"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = 6),
      plot.title = element_text(hjust = 0.5)
    )

  h <- max(7, 0.18 * nrow(z) + 3)
  save_plot(p, outfile_base, width = 11, height = h)
}

plot_sample_correlation <- function(mat, meta, outfile_base) {
  sample_order <- meta$sample[order(meta$group, meta$sample)]
  mat <- mat[, sample_order, drop = FALSE]

  cmat <- cor(mat, use = "pairwise.complete.obs", method = "spearman")
  df <- as.data.frame(as.table(cmat), stringsAsFactors = FALSE)
  colnames(df) <- c("Sample1", "Sample2", "Correlation")

  df$Sample1 <- factor(df$Sample1, levels = sample_order)
  df$Sample2 <- factor(df$Sample2, levels = rev(sample_order))

  p <- ggplot(df, aes(x = Sample1, y = Sample2, fill = Correlation)) +
    geom_tile() +
    theme_bw() +
    labs(
      title = "Sample correlation based on ssGSEA scores",
      x = "",
      y = "",
      fill = "Spearman r"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      plot.title = element_text(hjust = 0.5)
    )

  save_plot(p, outfile_base, width = 9, height = 8)
}

plot_pca <- function(mat, meta, outfile_base, pca_top_var = 500) {
  vars <- apply(mat, 1, var, na.rm = TRUE)
  selected <- names(sort(vars, decreasing = TRUE))[seq_len(min(pca_top_var, length(vars)))]

  pca_mat <- mat[selected, meta$sample, drop = FALSE]
  pca <- prcomp(t(pca_mat), center = TRUE, scale. = TRUE)

  var_exp <- summary(pca)$importance[2, ] * 100

  df <- data.frame(
    sample = rownames(pca$x),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    stringsAsFactors = FALSE
  )

  df <- merge(df, meta, by = "sample", all.x = TRUE)

  p <- ggplot(df, aes(x = PC1, y = PC2, shape = group)) +
    geom_point(size = 4) +
    geom_text(aes(label = sample), vjust = -0.8, size = 3) +
    theme_bw() +
    labs(
      title = "PCA based on ssGSEA scores",
      x = paste0("PC1 (", round(var_exp[1], 1), "%)"),
      y = paste0("PC2 (", round(var_exp[2], 1), "%)"),
      shape = "Group"
    ) +
    theme(plot.title = element_text(hjust = 0.5))

  save_plot(p, outfile_base, width = 9, height = 7)
}

run_differential <- function(mat, meta, desc, outdir, top_diff = 20) {
  groups <- levels(meta$group)

  if (length(groups) != 2) {
    msg <- paste0(
      "Differential ssGSEA boxplots skipped because the samplesheet has ",
      length(groups), " groups. This script currently makes differential plots only for two groups."
    )

    writeLines(msg, file.path(outdir, "ssGSEA_differential_skipped.txt"))
    message("[INFO] ", msg)
    return(invisible(NULL))
  }

  g1 <- groups[1]
  g2 <- groups[2]

  x1_idx <- meta$sample[meta$group == g1]
  x2_idx <- meta$sample[meta$group == g2]

  res_list <- lapply(rownames(mat), function(term) {
    x1 <- as.numeric(mat[term, x1_idx])
    x2 <- as.numeric(mat[term, x2_idx])

    p_t <- tryCatch(
      t.test(x1, x2)$p.value,
      error = function(e) NA_real_
    )

    p_w <- tryCatch(
      wilcox.test(x1, x2, exact = FALSE)$p.value,
      error = function(e) NA_real_
    )

    data.frame(
      term = term,
      mean_group1 = mean(x1, na.rm = TRUE),
      mean_group2 = mean(x2, na.rm = TRUE),
      diff_group1_minus_group2 = mean(x1, na.rm = TRUE) - mean(x2, na.rm = TRUE),
      pvalue_ttest = p_t,
      pvalue_wilcox = p_w,
      group1 = g1,
      group2 = g2,
      stringsAsFactors = FALSE
    )
  })

  res <- do.call(rbind, res_list)
  res$padj_ttest <- p.adjust(res$pvalue_ttest, method = "BH")
  res$padj_wilcox <- p.adjust(res$pvalue_wilcox, method = "BH")

  res <- merge(res, desc, by.x = "term", by.y = "term", all.x = TRUE)
  res <- res[order(res$padj_ttest, -abs(res$diff_group1_minus_group2)), ]

  out_tsv <- file.path(outdir, "ssGSEA_differential_scores.tsv")
  write.table(res, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

  message("[OK] Saved: ", out_tsv)

  top_terms <- res$term[is.finite(res$padj_ttest)]
  top_terms <- head(top_terms, min(top_diff, length(top_terms)))

  if (length(top_terms) == 0) {
    writeLines("No valid differential ssGSEA terms found.", file.path(outdir, "ssGSEA_differential_no_valid_terms.txt"))
    return(invisible(res))
  }

  long_list <- lapply(top_terms, function(term) {
    data.frame(
      term = term,
      sample = meta$sample,
      group = meta$group,
      score = as.numeric(mat[term, meta$sample]),
      stringsAsFactors = FALSE
    )
  })

  df <- do.call(rbind, long_list)

  df$term <- factor(df$term, levels = rev(top_terms))
  df$group <- factor(df$group, levels = groups)

  p <- ggplot(df, aes(x = group, y = score)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.12, size = 1.8, alpha = 0.8) +
    facet_wrap(~ term, scales = "free_y", ncol = 4) +
    theme_bw() +
    labs(
      title = paste0("Top differential ssGSEA GO scores: ", g1, " vs ", g2),
      x = "Group",
      y = "ssGSEA score"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 7),
      plot.title = element_text(hjust = 0.5)
    )

  h <- max(7, ceiling(length(top_terms) / 4) * 2.5)
  save_plot(p, file.path(outdir, "ssGSEA_top_differential_boxplots"), width = 13, height = h)

  invisible(res)
}

gct <- read_gct_scores(scores_file)
mat <- clean_matrix(gct$mat)

meta <- read_samplesheet(samplesheet_file, colnames(mat))
mat <- mat[, meta$sample, drop = FALSE]

desc <- gct$desc[gct$desc$term %in% rownames(mat), , drop = FALSE]

summary_file <- file.path(outdir, "ssGSEA_plot_summary.txt")
summary_lines <- c(
  paste0("scores_file\t", scores_file),
  paste0("samplesheet_file\t", samplesheet_file),
  paste0("n_terms_after_filter\t", nrow(mat)),
  paste0("n_samples\t", ncol(mat)),
  paste0("samples\t", paste(colnames(mat), collapse = ",")),
  paste0("groups\t", paste(levels(meta$group), collapse = ",")),
  paste0("top_var\t", top_var),
  paste0("top_diff\t", top_diff),
  paste0("pca_top_var\t", pca_top_var)
)
writeLines(summary_lines, summary_file)
message("[OK] Saved: ", summary_file)

vars <- apply(mat, 1, var, na.rm = TRUE)
top_var_terms <- names(sort(vars, decreasing = TRUE))[seq_len(min(top_var, length(vars)))]

plot_heatmap(
  mat[top_var_terms, , drop = FALSE],
  meta,
  file.path(outdir, "ssGSEA_top_variable_heatmap"),
  paste0("Top ", length(top_var_terms), " most variable ssGSEA GO scores")
)

plot_sample_correlation(
  mat,
  meta,
  file.path(outdir, "ssGSEA_sample_correlation_heatmap")
)

plot_pca(
  mat,
  meta,
  file.path(outdir, "ssGSEA_PCA_samples"),
  pca_top_var = pca_top_var
)

run_differential(
  mat,
  meta,
  desc,
  outdir,
  top_diff = top_diff
)

message("[OK] ssGSEA plotting step finished.")
