#!/usr/bin/env Rscript

# ============================================================
# MTD species overlap Venn + Euler diagrams
# ------------------------------------------------------------
# Correct source for MTD:
#   bracken_species_all + samplesheet.csv
#
# This avoids counting species directly from Krona files.
# Presence is calculated from the abundance matrix:
#   species is present in a group if abundance > threshold
#   in at least one sample from that group.
# ============================================================

print_help <- function(exit_status = 0) {
  cat("
Generate Venn/Euler diagrams for detected microbiome species using
the MTD Bracken species matrix and samplesheet.

Required:
  --input FILE              Combined Bracken species table, usually:
                            bracken_species_all

  --samplesheet FILE        samplesheet.csv used by MTD

  --group1 NAME             First group label, e.g. Liver

  --group2 NAME             Second group label, e.g. Telencephalon

Optional:
  --sample_col COL          Sample column in samplesheet.
                            Can be column name or column number.
                            Default: 1

  --group_col COL           Group column in samplesheet.
                            Can be column name or column number.
                            Default: 2

  --taxon_col COL           Taxon/species column in Bracken table.
                            Default: auto

  --presence_threshold NUM  Species is considered present if abundance > threshold.
                            Default: 0

  --title TEXT              Plot title.
                            Default: Total species detected

  --subtitle TEXT           Plot subtitle.
                            Default: automatic summary

  --output_prefix PREFIX    Output file prefix.
                            Default: species_overlap_<group1>_vs_<group2>

  --plot_type TYPE          both, venn, or euler.
                            Default: both

  --label_type TYPE         count, percent, or both.
                            Default: count

  --width NUM               Figure width in inches.
                            Default: 10

  --height NUM              Figure height in inches.
                            Default: 8

  --write_lists             Write species list files.

  --replace_underscores     Replace underscores with spaces in taxon names.

  --no_install              Do not try to install missing R packages.

Examples:
  Rscript MTD_species_overlap_venn_euler.R \\
    --input /path/to/bracken_species_all \\
    --samplesheet /path/to/samplesheet.csv \\
    --group1 Liver \\
    --group2 Telencephalon \\
    --output_prefix Liver_vs_Telencephalon_species_overlap \\
    --write_lists
")
  quit(status = exit_status)
}

parse_args <- function(args) {
  if (length(args) == 0 || any(args %in% c("-h", "--help"))) {
    print_help(0)
  }

  opt <- list(
    input = NULL,
    samplesheet = NULL,
    group1 = NULL,
    group2 = NULL,
    sample_col = "1",
    group_col = "2",
    taxon_col = "auto",
    presence_threshold = 0,
    title = "Total species detected",
    subtitle = NULL,
    output_prefix = NULL,
    plot_type = "both",
    label_type = "count",
    width = 10,
    height = 8,
    write_lists = FALSE,
    replace_underscores = FALSE,
    auto_install = TRUE
  )

  known <- c(
    "input",
    "samplesheet",
    "group1",
    "group2",
    "sample_col",
    "group_col",
    "taxon_col",
    "presence_threshold",
    "title",
    "subtitle",
    "output_prefix",
    "plot_type",
    "label_type",
    "width",
    "height",
    "write_lists",
    "replace_underscores",
    "no_install"
  )

  i <- 1

  while (i <= length(args)) {
    key <- args[i]

    if (!startsWith(key, "--")) {
      stop("Unexpected argument: ", key, call. = FALSE)
    }

    name <- sub("^--", "", key)

    if (!name %in% known) {
      stop("Unknown option: --", name, call. = FALSE)
    }

    if (name == "write_lists") {
      opt$write_lists <- TRUE
      i <- i + 1
      next
    }

    if (name == "replace_underscores") {
      opt$replace_underscores <- TRUE
      i <- i + 1
      next
    }

    if (name == "no_install") {
      opt$auto_install <- FALSE
      i <- i + 1
      next
    }

    i <- i + 1
    values <- character()

    while (i <= length(args) && !startsWith(args[i], "--")) {
      values <- c(values, args[i])
      i <- i + 1
    }

    if (length(values) == 0) {
      stop("Missing value for --", name, call. = FALSE)
    }

    if (name %in% c("width", "height", "presence_threshold")) {
      opt[[name]] <- as.numeric(values[1])
    } else {
      opt[[name]] <- paste(values, collapse = " ")
    }
  }

  required <- c("input", "samplesheet", "group1", "group2")

  for (x in required) {
    if (is.null(opt[[x]]) || trimws(opt[[x]]) == "") {
      stop("Missing required argument --", x, call. = FALSE)
    }
  }

  if (!file.exists(opt$input)) {
    stop("Input file does not exist: ", opt$input, call. = FALSE)
  }

  if (!file.exists(opt$samplesheet)) {
    stop("Samplesheet file does not exist: ", opt$samplesheet, call. = FALSE)
  }

  if (!opt$plot_type %in% c("both", "venn", "euler")) {
    stop("--plot_type must be one of: both, venn, euler", call. = FALSE)
  }

  if (!opt$label_type %in% c("count", "percent", "both")) {
    stop("--label_type must be one of: count, percent, both", call. = FALSE)
  }

  if (is.na(opt$presence_threshold)) {
    stop("--presence_threshold must be numeric.", call. = FALSE)
  }

  opt
}

setup_user_library <- function() {
  user_lib <- Sys.getenv("R_LIBS_USER")

  if (user_lib == "") {
    user_lib <- file.path(Sys.getenv("HOME"), "R", "library")
    Sys.setenv(R_LIBS_USER = user_lib)
  }

  if (!dir.exists(user_lib)) {
    dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
  }

  .libPaths(unique(c(user_lib, .libPaths())))

  invisible(user_lib)
}

check_and_install_packages <- function(pkgs, auto_install = TRUE) {
  missing <- pkgs[
    !vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
  ]

  if (length(missing) == 0) {
    return(invisible(TRUE))
  }

  message("[INFO] Missing R package(s): ", paste(missing, collapse = ", "))

  cmd <- paste0(
    "Rscript -e 'install.packages(c(",
    paste(sprintf("\"%s\"", missing), collapse = ", "),
    "), repos=\"https://cloud.r-project.org\")'"
  )

  if (!auto_install) {
    stop(
      "Missing package(s): ", paste(missing, collapse = ", "),
      "\nInstall manually with:\n  ", cmd,
      call. = FALSE
    )
  }

  setup_user_library()

  install.packages(
    missing,
    repos = "https://cloud.r-project.org",
    dependencies = TRUE
  )

  still_missing <- missing[
    !vapply(missing, requireNamespace, logical(1), quietly = TRUE)
  ]

  if (length(still_missing) > 0) {
    stop(
      "Could not install/load package(s): ",
      paste(still_missing, collapse = ", "),
      "\nInstall manually with:\n  ", cmd,
      call. = FALSE
    )
  }

  invisible(TRUE)
}

sanitize_name <- function(x) {
  x <- gsub("[[:space:]]+", "_", x)
  x <- gsub("[^A-Za-z0-9_.-]", "_", x)
  x
}

clean_taxon_name_one <- function(x, replace_underscores = FALSE) {
  x <- trimws(as.character(x))

  if (length(x) == 0 || is.na(x) || x == "") {
    return("")
  }

  # If taxonomy path is accidentally provided, keep the last species part.
  if (grepl("[|;]", x)) {
    parts <- unlist(strsplit(x, "[|;]"))
    parts <- trimws(parts)

    species_parts <- parts[grepl("^s__", parts)]

    if (length(species_parts) > 0) {
      x <- species_parts[length(species_parts)]
    } else {
      x <- parts[length(parts)]
    }
  }

  x <- sub("^s__", "", x)
  x <- sub("^.*s__", "", x)

  if (replace_underscores) {
    x <- gsub("_", " ", x)
  }

  x <- gsub("[[:space:]]+", " ", x)
  trimws(x)
}

clean_taxon_name <- function(x, replace_underscores = FALSE) {
  vapply(
    x,
    function(z) clean_taxon_name_one(
      z,
      replace_underscores = replace_underscores
    ),
    character(1)
  )
}

get_col <- function(df, col_spec, label) {
  nms <- colnames(df)

  if (grepl("^[0-9]+$", col_spec)) {
    idx <- as.integer(col_spec)

    if (idx < 1 || idx > ncol(df)) {
      stop(label, " column index is out of range: ", col_spec, call. = FALSE)
    }

    return(nms[idx])
  }

  if (!col_spec %in% nms) {
    stop(
      label, " column not found: ", col_spec,
      "\nAvailable columns:\n  ",
      paste(nms, collapse = ", "),
      call. = FALSE
    )
  }

  col_spec
}

guess_taxon_col <- function(df, taxon_col = "auto") {
  nms <- colnames(df)

  if (taxon_col != "auto") {
    if (grepl("^[0-9]+$", taxon_col)) {
      idx <- as.integer(taxon_col)

      if (idx < 1 || idx > ncol(df)) {
        stop("--taxon_col index is out of range: ", taxon_col, call. = FALSE)
      }

      return(nms[idx])
    }

    if (!taxon_col %in% nms) {
      stop("--taxon_col not found: ", taxon_col, call. = FALSE)
    }

    return(taxon_col)
  }

  candidates <- c(
    "name",
    "taxon",
    "Taxon",
    "taxonomy",
    "classification",
    "#Classification",
    "ID"
  )

  hit <- candidates[candidates %in% nms]

  if (length(hit) > 0) {
    return(hit[1])
  }

  nms[1]
}

to_numeric_safe <- function(x) {
  x <- as.character(x)
  x <- gsub(",", "", x)
  x <- gsub("%", "", x)
  x[x %in% c("", "NA", "NaN", "nan", "null", "NULL")] <- "0"

  out <- suppressWarnings(as.numeric(x))
  out[is.na(out)] <- 0
  out
}

match_sample_columns <- function(sample_names, matrix_cols, metadata_cols) {
  usable_cols <- setdiff(matrix_cols, metadata_cols)

  out <- character()
  missing <- character()

  for (s in sample_names) {
    if (s %in% usable_cols) {
      out <- c(out, s)
      next
    }

    hit <- usable_cols[grepl(s, usable_cols, fixed = TRUE)]

    if (length(hit) == 1) {
      out <- c(out, hit)
      next
    }

    if (length(hit) > 1) {
      hit <- hit[order(nchar(hit))]
      out <- c(out, hit[1])
      next
    }

    missing <- c(missing, s)
  }

  if (length(missing) > 0) {
    stop(
      "Could not match the following samples from samplesheet to columns in the Bracken table:\n  ",
      paste(missing, collapse = ", "),
      "\n\nAvailable abundance columns detected:\n  ",
      paste(usable_cols, collapse = ", "),
      call. = FALSE
    )
  }

  names(out) <- sample_names
  out
}

write_species_file <- function(x, file) {
  writeLines(sort(unique(x)), con = file)
}

safe_metric <- function(x, field) {
  tryCatch(
    {
      if (!is.null(x[[field]])) {
        as.numeric(x[[field]])[1]
      } else {
        NA_real_
      }
    },
    error = function(e) NA_real_
  )
}

save_ggplot_all <- function(plot, prefix, width, height) {
  png_file <- paste0(prefix, ".png")
  pdf_file <- paste0(prefix, ".pdf")
  svg_file <- paste0(prefix, ".svg")

  dir.create(dirname(png_file), recursive = TRUE, showWarnings = FALSE)

  grDevices::png(
    png_file,
    width = width,
    height = height,
    units = "in",
    res = 300,
    bg = "white"
  )
  print(plot)
  grDevices::dev.off()

  ggplot2::ggsave(
    filename = pdf_file,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    bg = "white"
  )

  grDevices::svg(svg_file, width = width, height = height, bg = "white")
  print(plot)
  grDevices::dev.off()

  message("[OK] Files generated:")
  message("  ", png_file)
  message("  ", pdf_file)
  message("  ", svg_file)
}

make_euler_grob <- function(fit, quantities_arg) {
  grid::grid.grabExpr({
    euler_plot <- plot(
      fit,
      fills = list(
        fill = c("#E9C9CF", "#BFE5B8"),
        alpha = 0.65
      ),
      edges = list(
        col = "black",
        lwd = 2
      ),
      labels = list(
        cex = 0.95,
        font = 2
      ),
      quantities = quantities_arg,
      main = NULL
    )

    if (!is.null(euler_plot)) {
      try(grid::grid.draw(euler_plot), silent = TRUE)
    }
  })
}

save_euler_all <- function(fit,
                           prefix,
                           width,
                           height,
                           title,
                           subtitle,
                           caption,
                           label_type) {
  png_file <- paste0(prefix, ".png")
  pdf_file <- paste0(prefix, ".pdf")
  svg_file <- paste0(prefix, ".svg")

  dir.create(dirname(png_file), recursive = TRUE, showWarnings = FALSE)

  if (label_type == "percent") {
    quantities_arg <- list(type = "percent", cex = 1.1, font = 2)
  } else if (label_type == "both") {
    quantities_arg <- list(type = c("counts", "percent"), cex = 1.0, font = 2)
  } else {
    quantities_arg <- list(type = "counts", cex = 1.2, font = 2)
  }

  euler_grob <- make_euler_grob(
    fit = fit,
    quantities_arg = quantities_arg
  )

  draw_one <- function() {
    grid::grid.newpage()

    grid::pushViewport(
      grid::viewport(
        x = 0.5,
        y = 0.50,
        width = 0.88,
        height = 0.70
      )
    )
    grid::grid.draw(euler_grob)
    grid::popViewport()

    grid::grid.text(
      title,
      x = 0.5,
      y = 0.965,
      gp = grid::gpar(fontsize = 18, fontface = "bold")
    )

    grid::grid.text(
      subtitle,
      x = 0.5,
      y = 0.925,
      gp = grid::gpar(fontsize = 10)
    )

    grid::grid.text(
      caption,
      x = 0.5,
      y = 0.035,
      gp = grid::gpar(fontsize = 8, col = "grey35")
    )
  }

  grDevices::png(
    png_file,
    width = width,
    height = height,
    units = "in",
    res = 300,
    bg = "white"
  )
  draw_one()
  grDevices::dev.off()

  grDevices::pdf(pdf_file, width = width, height = height, bg = "white")
  draw_one()
  grDevices::dev.off()

  grDevices::svg(svg_file, width = width, height = height, bg = "white")
  draw_one()
  grDevices::dev.off()

  message("[OK] Files generated:")
  message("  ", png_file)
  message("  ", pdf_file)
  message("  ", svg_file)
}

# ============================================================
# Main
# ============================================================

args <- commandArgs(trailingOnly = TRUE)
opt <- parse_args(args)

required_packages <- c("ggplot2", "ggVennDiagram", "eulerr")
check_and_install_packages(required_packages, auto_install = opt$auto_install)

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggVennDiagram)
  library(eulerr)
  library(grid)
})

bracken <- read.delim(
  opt$input,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  quote = "",
  comment.char = "",
  stringsAsFactors = FALSE
)

samplesheet <- read.csv(
  opt$samplesheet,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

sample_col <- get_col(samplesheet, opt$sample_col, "Sample")
group_col <- get_col(samplesheet, opt$group_col, "Group")
taxon_col <- guess_taxon_col(bracken, opt$taxon_col)

samplesheet[[sample_col]] <- trimws(as.character(samplesheet[[sample_col]]))
samplesheet[[group_col]] <- trimws(as.character(samplesheet[[group_col]]))

group1_samples <- samplesheet[[sample_col]][samplesheet[[group_col]] == opt$group1]
group2_samples <- samplesheet[[sample_col]][samplesheet[[group_col]] == opt$group2]

if (length(group1_samples) == 0) {
  stop("No samples found for --group1: ", opt$group1, call. = FALSE)
}

if (length(group2_samples) == 0) {
  stop("No samples found for --group2: ", opt$group2, call. = FALSE)
}

metadata_cols <- c(
  taxon_col,
  "taxonomy_id",
  "taxonomy_lvl",
  "rank",
  "taxid",
  "new_est_reads",
  "fraction_total_reads"
)

metadata_cols <- metadata_cols[metadata_cols %in% colnames(bracken)]

group1_cols <- match_sample_columns(
  sample_names = group1_samples,
  matrix_cols = colnames(bracken),
  metadata_cols = metadata_cols
)

group2_cols <- match_sample_columns(
  sample_names = group2_samples,
  matrix_cols = colnames(bracken),
  metadata_cols = metadata_cols
)

all_selected_cols <- unique(c(group1_cols, group2_cols))

taxa <- clean_taxon_name(
  bracken[[taxon_col]],
  replace_underscores = opt$replace_underscores
)

valid_taxa <- taxa != "" &
  !grepl(
    "^(NA|NaN|null|unknown|unclassified|uncultured)$",
    taxa,
    ignore.case = TRUE
  )

bracken <- bracken[valid_taxa, , drop = FALSE]
taxa <- taxa[valid_taxa]

abundance_matrix <- as.data.frame(
  lapply(bracken[, all_selected_cols, drop = FALSE], to_numeric_safe),
  check.names = FALSE
)

get_present_taxa <- function(cols) {
  mat <- as.data.frame(
    lapply(bracken[, cols, drop = FALSE], to_numeric_safe),
    check.names = FALSE
  )

  present <- rowSums(mat > opt$presence_threshold, na.rm = TRUE) > 0
  sort(unique(taxa[present]))
}

species1 <- get_present_taxa(group1_cols)
species2 <- get_present_taxa(group2_cols)

present_any_selected <- rowSums(
  abundance_matrix > opt$presence_threshold,
  na.rm = TRUE
) > 0

detected_selected <- sort(unique(taxa[present_any_selected]))

unique1 <- setdiff(species1, species2)
unique2 <- setdiff(species2, species1)
shared <- intersect(species1, species2)
union_species <- union(species1, species2)

n1 <- length(species1)
n2 <- length(species2)
n_unique1 <- length(unique1)
n_unique2 <- length(unique2)
n_shared <- length(shared)
n_union <- length(union_species)
n_detected_selected <- length(detected_selected)

if (n_union != n_detected_selected) {
  warning(
    "Union of group1/group2 species is not equal to all detected species ",
    "across selected samples. This usually means samples outside group1/group2 ",
    "exist in the samplesheet or sample matching needs checking."
  )
}

label1 <- opt$group1
label2 <- opt$group2

if (is.null(opt$output_prefix)) {
  output_prefix <- paste0(
    "species_overlap_",
    sanitize_name(label1),
    "_vs_",
    sanitize_name(label2)
  )
} else {
  output_prefix <- opt$output_prefix
}

dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)

message("============================================================")
message("MTD species overlap diagrams")
message("Input matrix: ", opt$input)
message("Samplesheet: ", opt$samplesheet)
message("Taxon column: ", taxon_col)
message("Group 1: ", label1, " | samples: ", paste(group1_samples, collapse = ", "))
message("Group 2: ", label2, " | samples: ", paste(group2_samples, collapse = ", "))
message("Presence threshold: > ", opt$presence_threshold)
message("Output prefix: ", output_prefix)
message("============================================================")

message("")
message("Summary:")
message("  ", label1, " total species: ", n1)
message("  ", label2, " total species: ", n2)
message("  Unique to ", label1, ": ", n_unique1)
message("  Unique to ", label2, ": ", n_unique2)
message("  Shared species: ", n_shared)
message("  Total union: ", n_union)
message("  Total detected across selected samples: ", n_detected_selected)
message("")

stopifnot(n_unique1 + n_shared == n1)
stopifnot(n_unique2 + n_shared == n2)
stopifnot(n_unique1 + n_unique2 + n_shared == n_union)

summary_df <- data.frame(
  comparison = paste(label1, "vs", label2),
  group1 = label1,
  group2 = label2,
  group1_samples = paste(group1_samples, collapse = ";"),
  group2_samples = paste(group2_samples, collapse = ";"),
  group1_total_species = n1,
  group2_total_species = n2,
  unique_to_group1 = n_unique1,
  unique_to_group2 = n_unique2,
  shared_species = n_shared,
  union_species = n_union,
  detected_species_across_selected_samples = n_detected_selected,
  presence_threshold = opt$presence_threshold,
  input_file = opt$input,
  stringsAsFactors = FALSE
)

summary_file <- paste0(output_prefix, "_summary.tsv")

write.table(
  summary_df,
  file = summary_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("[OK] Summary table written:")
message("  ", summary_file)

if (opt$write_lists) {
  write_species_file(species1, paste0(output_prefix, "_group1_species.txt"))
  write_species_file(species2, paste0(output_prefix, "_group2_species.txt"))
  write_species_file(unique1, paste0(output_prefix, "_unique_group1.txt"))
  write_species_file(unique2, paste0(output_prefix, "_unique_group2.txt"))
  write_species_file(shared, paste0(output_prefix, "_shared_species.txt"))
  write_species_file(union_species, paste0(output_prefix, "_union_species.txt"))

  message("[OK] Species TXT files written.")
}

if (is.null(opt$subtitle)) {
  subtitle_text <- paste0(
    label1, ": ", n1,
    " species | ",
    label2, ": ", n2,
    " species | Shared: ", n_shared,
    " | Union: ", n_union
  )
} else {
  subtitle_text <- opt$subtitle
}

caption_text <- paste0(
  "Presence defined as abundance > ",
  opt$presence_threshold,
  " in at least one sample from each group. ",
  "Numbers inside regions represent disjoint species counts."
)

# ============================================================
# Venn diagram
# ============================================================

if (opt$plot_type %in% c("both", "venn")) {
  message("[PLOT] Generating Venn diagram...")

  venn_input <- list()
  venn_input[[label1]] <- species1
  venn_input[[label2]] <- species2

  p_venn <- suppressMessages({
    ggVennDiagram::ggVennDiagram(
      venn_input,
      label = opt$label_type,
      label_alpha = 0,
      label_geom = "label",
      label_size = 5,
      edge_size = 1.1,
      set_size = 5.2
    ) +
      ggplot2::scale_fill_gradient(
        low = "#FFF8F9",
        high = "#BFE5B8",
        name = "Region\nspecies\ncount"
      ) +
      ggplot2::labs(
        title = opt$title,
        subtitle = subtitle_text,
        caption = caption_text
      ) +
      ggplot2::scale_x_continuous(
        expand = ggplot2::expansion(mult = c(0.10, 0.08))
      ) +
      ggplot2::scale_y_continuous(
        expand = ggplot2::expansion(mult = c(0.10, 0.10))
      ) +
      ggplot2::coord_fixed(clip = "off") +
      ggplot2::theme_void(base_size = 16) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(
          hjust = 0.5,
          face = "bold",
          size = 22,
          margin = ggplot2::margin(b = 8)
        ),
        plot.subtitle = ggplot2::element_text(
          hjust = 0.5,
          size = 13,
          margin = ggplot2::margin(b = 20)
        ),
        plot.caption = ggplot2::element_text(
          hjust = 0.5,
          size = 10,
          color = "grey35",
          margin = ggplot2::margin(t = 18)
        ),
        legend.position = "right",
        legend.title = ggplot2::element_text(size = 11, face = "bold"),
        legend.text = ggplot2::element_text(size = 10),
        plot.margin = ggplot2::margin(20, 25, 20, 55)
      )
  })

  save_ggplot_all(
    plot = p_venn,
    prefix = paste0(output_prefix, "_venn"),
    width = opt$width,
    height = opt$height
  )
}

# ============================================================
# Euler diagram
# ============================================================

if (opt$plot_type %in% c("both", "euler")) {
  message("[PLOT] Generating proportional Euler diagram...")

  euler_counts <- numeric(3)

  names(euler_counts) <- c(
    label1,
    label2,
    paste(label1, label2, sep = "&")
  )

  euler_counts[label1] <- n_unique1
  euler_counts[label2] <- n_unique2
  euler_counts[paste(label1, label2, sep = "&")] <- n_shared

  euler_fit <- eulerr::euler(euler_counts)

  euler_stress <- safe_metric(euler_fit, "stress")
  euler_diag_error <- safe_metric(euler_fit, "diagError")

  euler_caption <- paste0(
    "Euler diagram: circle areas are fitted to approximate disjoint species counts. ",
    "Stress: ", signif(euler_stress, 4),
    "; diagnostic error: ", signif(euler_diag_error, 4), "."
  )

  save_euler_all(
    fit = euler_fit,
    prefix = paste0(output_prefix, "_euler"),
    width = opt$width,
    height = opt$height,
    title = opt$title,
    subtitle = subtitle_text,
    caption = euler_caption,
    label_type = opt$label_type
  )
}

message("")
message("[DONE] Finished successfully.")
message("")
