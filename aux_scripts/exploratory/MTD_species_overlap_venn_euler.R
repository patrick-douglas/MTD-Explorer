#!/usr/bin/env Rscript

# ============================================================
# MTD Explorer
# MTD species overlap Venn + Euler diagrams
# VERSION: 3.0.0 MULTI-GROUP
# ------------------------------------------------------------
# Source:
#   combined Bracken abundance table + mandatory samplesheet.csv
#
# Presence rule:
#   a species is present in a group when abundance is greater
#   than the selected threshold in at least one sample belonging
#   to that group.
#
# Group behavior:
#   - default: use every non-empty group in samplesheet column 2
#   - optional: --groups GroupA,GroupB,GroupC
#   - backward compatibility: --group1 and --group2
#
# Outputs:
#   <prefix>_venn.png/.pdf/.svg
#   <prefix>_euler.png/.pdf/.svg
#   <prefix>_summary.tsv
#   optional species lists per group plus shared/union lists
# ============================================================

SCRIPT_VERSION <- "3.0.0"

print_help <- function(exit_status = 0L) {
  cat(
"MTD species overlap Venn + Euler diagrams — multi-group version

Required:
  --input FILE
      Combined Bracken table, usually bracken_species_all.

  --samplesheet FILE
      Mandatory MTD samplesheet.csv.

Group selection:
  --groups LIST
      Comma-separated group names to include.
      Default: all non-empty groups found in the samplesheet.

  --group1 NAME
  --group2 NAME
      Backward-compatible two-group selection.
      Use both options together. Ignored when --groups is supplied.

Optional:
  --sample_col COL
      Sample column name or number in samplesheet. Default: 1

  --group_col COL
      Group column name or number in samplesheet. Default: 2

  --taxon_col COL
      Species/taxon column in Bracken table. Default: auto

  --presence_threshold NUM
      Species is present when abundance is greater than this value
      in at least one sample from the corresponding group. Default: 0

  --title TEXT
      Plot title. Default: Total species detected

  --subtitle TEXT
      Plot subtitle. Default: automatic group/union summary

  --output_prefix PREFIX
      Output prefix. Default:
        species_overlap_<group1>_vs_<group2> for two groups
        species_overlap_all_<N>_groups for more than two groups

  --plot_type TYPE
      both, venn, or euler. Default: both

  --label_type TYPE
      count, percent, or both. Default: count

  --width NUM
      Figure width in inches. Default: 10

  --height NUM
      Figure height in inches. Default: 8

  --write_lists
      Write one species list per group plus shared-by-all and union lists.

  --replace_underscores
      Replace underscores with spaces in displayed/stored species names.

  --no_install
      Do not attempt installation of missing R packages.

  -v, --version
      Show version.

  -h, --help
      Show help.

Examples:
  # All groups in samplesheet:
  Rscript MTD_species_overlap_venn_euler.R \\
    --input /path/to/bracken_species_all \\
    --samplesheet /path/to/samplesheet.csv \\
    --output_prefix species_overlap_all_groups \\
    --write_lists

  # Selected groups:
  Rscript MTD_species_overlap_venn_euler.R \\
    --input /path/to/bracken_species_all \\
    --samplesheet /path/to/samplesheet.csv \\
    --groups Liver,Telencephalon,Kidney \\
    --output_prefix three_tissues_species_overlap
"
  )

  quit(status = exit_status)
}

parse_args <- function(args) {
  if (any(args %in% c("-v", "--version"))) {
    cat(
      "MTD_species_overlap_venn_euler.R version ",
      SCRIPT_VERSION,
      " — multi-group Venn + Euler\n",
      sep = ""
    )
    quit(status = 0L)
  }

  if (length(args) == 0L || any(args %in% c("-h", "--help"))) {
    print_help(0L)
  }

  opt <- list(
    input = NULL,
    samplesheet = NULL,
    groups = NULL,
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

  value_options <- c(
    "input",
    "samplesheet",
    "groups",
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
    "height"
  )

  flag_options <- c(
    "write_lists",
    "replace_underscores",
    "no_install"
  )

  known <- c(value_options, flag_options)

  i <- 1L

  while (i <= length(args)) {
    key <- args[i]

    if (!startsWith(key, "--")) {
      stop("Unexpected argument: ", key, call. = FALSE)
    }

    name <- sub("^--", "", key)

    if (!name %in% known) {
      stop("Unknown option: --", name, call. = FALSE)
    }

    if (name %in% flag_options) {
      if (name == "write_lists") {
        opt$write_lists <- TRUE
      } else if (name == "replace_underscores") {
        opt$replace_underscores <- TRUE
      } else if (name == "no_install") {
        opt$auto_install <- FALSE
      }

      i <- i + 1L
      next
    }

    i <- i + 1L
    values <- character()

    while (i <= length(args) && !startsWith(args[i], "--")) {
      values <- c(values, args[i])
      i <- i + 1L
    }

    if (length(values) == 0L) {
      stop("Missing value for --", name, call. = FALSE)
    }

    value <- paste(values, collapse = " ")

    if (name %in% c("presence_threshold", "width", "height")) {
      value <- suppressWarnings(as.numeric(value))

      if (is.na(value)) {
        stop("--", name, " must be numeric.", call. = FALSE)
      }
    }

    opt[[name]] <- value
  }

  for (required in c("input", "samplesheet")) {
    if (is.null(opt[[required]]) || trimws(opt[[required]]) == "") {
      stop("Missing required argument --", required, call. = FALSE)
    }
  }

  if (!file.exists(opt$input)) {
    stop("Input file does not exist: ", opt$input, call. = FALSE)
  }

  if (!file.exists(opt$samplesheet)) {
    stop("Samplesheet file does not exist: ", opt$samplesheet, call. = FALSE)
  }

  if (!opt$plot_type %in% c("both", "venn", "euler")) {
    stop("--plot_type must be both, venn, or euler.", call. = FALSE)
  }

  if (!opt$label_type %in% c("count", "percent", "both")) {
    stop("--label_type must be count, percent, or both.", call. = FALSE)
  }

  if (opt$width <= 0 || opt$height <= 0) {
    stop("--width and --height must be greater than zero.", call. = FALSE)
  }

  if (xor(is.null(opt$group1), is.null(opt$group2))) {
    stop(
      "--group1 and --group2 must be supplied together.",
      call. = FALSE
    )
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
    dir.create(
      user_lib,
      recursive = TRUE,
      showWarnings = FALSE
    )
  }

  .libPaths(unique(c(user_lib, .libPaths())))
  invisible(user_lib)
}

check_and_install_packages <- function(
    packages,
    auto_install = TRUE
) {
  missing <- packages[
    !vapply(
      packages,
      requireNamespace,
      logical(1),
      quietly = TRUE
    )
  ]

  if (length(missing) == 0L) {
    return(invisible(TRUE))
  }

  message(
    "[INFO] Missing R package(s): ",
    paste(missing, collapse = ", ")
  )

  if (!auto_install) {
    stop(
      "Missing package(s): ",
      paste(missing, collapse = ", "),
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
    !vapply(
      missing,
      requireNamespace,
      logical(1),
      quietly = TRUE
    )
  ]

  if (length(still_missing) > 0L) {
    stop(
      "Could not install/load package(s): ",
      paste(still_missing, collapse = ", "),
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

clean_taxon_name_one <- function(
    x,
    replace_underscores = FALSE
) {
  x <- trimws(as.character(x))

  if (length(x) == 0L || is.na(x) || x == "") {
    return("")
  }

  if (grepl("[|;]", x)) {
    parts <- unlist(strsplit(x, "[|;]"))
    parts <- trimws(parts)

    species_parts <- parts[
      grepl("^s__", parts)
    ]

    if (length(species_parts) > 0L) {
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

clean_taxon_name <- function(
    x,
    replace_underscores = FALSE
) {
  vapply(
    x,
    function(value) {
      clean_taxon_name_one(
        value,
        replace_underscores = replace_underscores
      )
    },
    character(1)
  )
}

get_col <- function(
    data,
    col_spec,
    label
) {
  names_available <- colnames(data)

  if (grepl("^[0-9]+$", col_spec)) {
    index <- as.integer(col_spec)

    if (index < 1L || index > ncol(data)) {
      stop(
        label,
        " column index is out of range: ",
        col_spec,
        call. = FALSE
      )
    }

    return(names_available[index])
  }

  if (!col_spec %in% names_available) {
    stop(
      label,
      " column not found: ",
      col_spec,
      "\nAvailable columns:\n  ",
      paste(names_available, collapse = ", "),
      call. = FALSE
    )
  }

  col_spec
}

guess_taxon_col <- function(
    data,
    taxon_col = "auto"
) {
  names_available <- colnames(data)

  if (taxon_col != "auto") {
    return(
      get_col(
        data,
        taxon_col,
        "Taxon"
      )
    )
  }

  candidates <- c(
    "name",
    "Name",
    "taxon",
    "Taxon",
    "taxonomy",
    "classification",
    "#Classification",
    "ID"
  )

  hits <- candidates[
    candidates %in% names_available
  ]

  if (length(hits) > 0L) {
    return(hits[1])
  }

  names_available[1]
}

to_numeric_safe <- function(x) {
  x <- as.character(x)
  x <- gsub(",", "", x, fixed = TRUE)
  x <- gsub("%", "", x, fixed = TRUE)
  x[
    x %in% c(
      "",
      "NA",
      "NaN",
      "nan",
      "null",
      "NULL"
    )
  ] <- "0"

  output <- suppressWarnings(as.numeric(x))
  output[is.na(output)] <- 0
  output
}

match_sample_columns <- function(
    sample_names,
    matrix_cols,
    metadata_cols
) {
  usable_cols <- setdiff(
    matrix_cols,
    metadata_cols
  )

  output <- character()
  missing <- character()

  for (sample_name in sample_names) {
    if (sample_name %in% usable_cols) {
      output <- c(output, sample_name)
      next
    }

    hits <- usable_cols[
      grepl(
        sample_name,
        usable_cols,
        fixed = TRUE
      )
    ]

    if (length(hits) == 1L) {
      output <- c(output, hits)
      next
    }

    if (length(hits) > 1L) {
      hits <- hits[
        order(nchar(hits))
      ]

      output <- c(output, hits[1])
      next
    }

    missing <- c(missing, sample_name)
  }

  if (length(missing) > 0L) {
    stop(
      "Could not match these samples to Bracken abundance columns:\n  ",
      paste(missing, collapse = ", "),
      "\nAvailable abundance columns:\n  ",
      paste(usable_cols, collapse = ", "),
      call. = FALSE
    )
  }

  names(output) <- sample_names
  output
}

write_species_file <- function(
    species,
    file
) {
  writeLines(
    sort(unique(species)),
    con = file
  )
}

safe_metric <- function(
    fit,
    field
) {
  tryCatch(
    {
      if (!is.null(fit[[field]])) {
        as.numeric(fit[[field]])[1]
      } else {
        NA_real_
      }
    },
    error = function(error) NA_real_
  )
}

safe_max_metric <- function(
    fit,
    field
) {
  tryCatch(
    {
      values <- as.numeric(fit[[field]])
      values <- values[is.finite(values)]

      if (length(values) == 0L) {
        NA_real_
      } else {
        max(abs(values))
      }
    },
    error = function(error) NA_real_
  )
}

group_palette <- function(number_of_sets) {
  base_colors <- c(
    "#E9C9CF",
    "#BFE5B8",
    "#BDD7EE",
    "#F7D9A6",
    "#D8C4E8",
    "#C4E3E8",
    "#E8D7B7"
  )

  rep(
    base_colors,
    length.out = number_of_sets
  )
}

save_ggplot_all <- function(
    plot,
    prefix,
    width,
    height
) {
  png_file <- paste0(prefix, ".png")
  pdf_file <- paste0(prefix, ".pdf")
  svg_file <- paste0(prefix, ".svg")

  dir.create(
    dirname(png_file),
    recursive = TRUE,
    showWarnings = FALSE
  )

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
    bg = "white",
    limitsize = FALSE
  )

  grDevices::svg(
    svg_file,
    width = width,
    height = height,
    bg = "white"
  )
  print(plot)
  grDevices::dev.off()

  message("[OK] Files generated:")
  message("  ", png_file)
  message("  ", pdf_file)
  message("  ", svg_file)
}

make_euler_grob <- function(
    fit,
    quantities_arg,
    number_of_sets
) {
  label_cex <- if (number_of_sets <= 3L) {
    0.95
  } else if (number_of_sets <= 5L) {
    0.78
  } else {
    0.64
  }

  edge_width <- if (number_of_sets <= 4L) {
    2
  } else {
    1.4
  }

  grid::grid.grabExpr({
    euler_plot <- plot(
      fit,
      fills = list(
        fill = group_palette(number_of_sets),
        alpha = 0.65
      ),
      edges = list(
        col = "black",
        lwd = edge_width
      ),
      labels = list(
        cex = label_cex,
        font = 2
      ),
      quantities = quantities_arg,
      main = NULL
    )

    if (!is.null(euler_plot)) {
      try(
        grid::grid.draw(euler_plot),
        silent = TRUE
      )
    }
  })
}

save_euler_all <- function(
    fit,
    prefix,
    width,
    height,
    title,
    subtitle,
    caption,
    label_type,
    number_of_sets
) {
  png_file <- paste0(prefix, ".png")
  pdf_file <- paste0(prefix, ".pdf")
  svg_file <- paste0(prefix, ".svg")

  dir.create(
    dirname(png_file),
    recursive = TRUE,
    showWarnings = FALSE
  )

  quantity_cex <- if (number_of_sets <= 3L) {
    1.05
  } else if (number_of_sets <= 5L) {
    0.82
  } else {
    0.64
  }

  if (label_type == "percent") {
    quantities_arg <- list(
      type = "percent",
      cex = quantity_cex,
      font = 2
    )
  } else if (label_type == "both") {
    quantities_arg <- list(
      type = c("counts", "percent"),
      cex = quantity_cex,
      font = 2
    )
  } else {
    quantities_arg <- list(
      type = "counts",
      cex = quantity_cex,
      font = 2
    )
  }

  euler_grob <- make_euler_grob(
    fit = fit,
    quantities_arg = quantities_arg,
    number_of_sets = number_of_sets
  )

  diagram_width <- if (number_of_sets <= 3L) {
    0.88
  } else if (number_of_sets <= 5L) {
    0.92
  } else {
    0.96
  }

  diagram_height <- if (number_of_sets <= 3L) {
    0.70
  } else {
    0.76
  }

  draw_one <- function() {
    grid::grid.newpage()

    grid::pushViewport(
      grid::viewport(
        x = 0.5,
        y = 0.50,
        width = diagram_width,
        height = diagram_height
      )
    )

    grid::grid.draw(euler_grob)
    grid::popViewport()

    grid::grid.text(
      title,
      x = 0.5,
      y = 0.965,
      gp = grid::gpar(
        fontsize = 18,
        fontface = "bold"
      )
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
      gp = grid::gpar(
        fontsize = 8,
        col = "grey35"
      )
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

  grDevices::pdf(
    pdf_file,
    width = width,
    height = height,
    bg = "white"
  )
  draw_one()
  grDevices::dev.off()

  grDevices::svg(
    svg_file,
    width = width,
    height = height,
    bg = "white"
  )
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

opt <- parse_args(
  commandArgs(trailingOnly = TRUE)
)

required_packages <- c(
  "ggplot2",
  "ggVennDiagram",
  "eulerr"
)

check_and_install_packages(
  required_packages,
  auto_install = opt$auto_install
)

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
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)

sample_col <- get_col(
  samplesheet,
  opt$sample_col,
  "Sample"
)

group_col <- get_col(
  samplesheet,
  opt$group_col,
  "Group"
)

taxon_col <- guess_taxon_col(
  bracken,
  opt$taxon_col
)

samplesheet[[sample_col]] <- trimws(
  as.character(samplesheet[[sample_col]])
)

samplesheet[[group_col]] <- trimws(
  as.character(samplesheet[[group_col]])
)

valid_rows <- (
  !is.na(samplesheet[[sample_col]]) &
  samplesheet[[sample_col]] != "" &
  !is.na(samplesheet[[group_col]]) &
  samplesheet[[group_col]] != ""
)

samplesheet <- samplesheet[
  valid_rows,
  ,
  drop = FALSE
]

if (nrow(samplesheet) == 0L) {
  stop(
    "No valid sample/group rows were found in samplesheet.",
    call. = FALSE
  )
}

duplicate_samples <- unique(
  samplesheet[[sample_col]][
    duplicated(samplesheet[[sample_col]])
  ]
)

if (length(duplicate_samples) > 0L) {
  stop(
    "Duplicated samples in samplesheet: ",
    paste(duplicate_samples, collapse = ", "),
    call. = FALSE
  )
}

available_groups <- unique(
  samplesheet[[group_col]]
)

if (!is.null(opt$groups)) {
  selected_groups <- trimws(
    unlist(
      strsplit(
        opt$groups,
        ",",
        fixed = TRUE
      )
    )
  )

  selected_groups <- selected_groups[
    selected_groups != ""
  ]
} else if (
  !is.null(opt$group1) &&
  !is.null(opt$group2)
) {
  selected_groups <- c(
    trimws(opt$group1),
    trimws(opt$group2)
  )
} else {
  selected_groups <- available_groups
}

selected_groups <- unique(selected_groups)

missing_groups <- setdiff(
  selected_groups,
  available_groups
)

if (length(missing_groups) > 0L) {
  stop(
    "Selected group(s) not found in samplesheet: ",
    paste(missing_groups, collapse = ", "),
    "\nAvailable groups: ",
    paste(available_groups, collapse = ", "),
    call. = FALSE
  )
}

if (length(selected_groups) < 2L) {
  stop(
    "At least two groups are required. Selected groups: ",
    paste(selected_groups, collapse = ", "),
    call. = FALSE
  )
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

metadata_cols <- metadata_cols[
  metadata_cols %in% colnames(bracken)
]

group_samples <- setNames(
  lapply(
    selected_groups,
    function(group_name) {
      samplesheet[[sample_col]][
        samplesheet[[group_col]] == group_name
      ]
    }
  ),
  selected_groups
)

empty_group_samples <- names(group_samples)[
  lengths(group_samples) == 0L
]

if (length(empty_group_samples) > 0L) {
  stop(
    "No samples found for group(s): ",
    paste(empty_group_samples, collapse = ", "),
    call. = FALSE
  )
}

all_selected_samples <- unique(
  unlist(
    group_samples,
    use.names = FALSE
  )
)

sample_columns <- match_sample_columns(
  sample_names = all_selected_samples,
  matrix_cols = colnames(bracken),
  metadata_cols = metadata_cols
)

group_columns <- setNames(
  lapply(
    group_samples,
    function(samples) {
      unname(sample_columns[samples])
    }
  ),
  selected_groups
)

taxa <- clean_taxon_name(
  bracken[[taxon_col]],
  replace_underscores = opt$replace_underscores
)

valid_taxa <- (
  taxa != "" &
  !grepl(
    "^(NA|NaN|null|unknown|unclassified|uncultured)$",
    taxa,
    ignore.case = TRUE
  )
)

bracken <- bracken[
  valid_taxa,
  ,
  drop = FALSE
]

taxa <- taxa[valid_taxa]

get_present_taxa <- function(columns) {
  abundance <- as.data.frame(
    lapply(
      bracken[
        ,
        columns,
        drop = FALSE
      ],
      to_numeric_safe
    ),
    check.names = FALSE
  )

  present <- rowSums(
    abundance > opt$presence_threshold,
    na.rm = TRUE
  ) > 0

  sort(unique(taxa[present]))
}

species_sets <- setNames(
  lapply(
    group_columns,
    get_present_taxa
  ),
  selected_groups
)

group_counts <- lengths(species_sets)
union_species <- sort(
  unique(
    unlist(
      species_sets,
      use.names = FALSE
    )
  )
)

shared_by_all <- sort(
  Reduce(
    intersect,
    species_sets
  )
)

all_selected_columns <- unique(
  unlist(
    group_columns,
    use.names = FALSE
  )
)

all_abundance <- as.data.frame(
  lapply(
    bracken[
      ,
      all_selected_columns,
      drop = FALSE
    ],
    to_numeric_safe
  ),
  check.names = FALSE
)

present_any_selected <- rowSums(
  all_abundance > opt$presence_threshold,
  na.rm = TRUE
) > 0

detected_selected <- sort(
  unique(taxa[present_any_selected])
)

if (!setequal(union_species, detected_selected)) {
  stop(
    "Internal validation failed: union of group species does not ",
    "match species detected across all selected samples.",
    call. = FALSE
  )
}

number_of_groups <- length(selected_groups)

if (is.null(opt$output_prefix)) {
  if (number_of_groups == 2L) {
    output_prefix <- paste0(
      "species_overlap_",
      sanitize_name(selected_groups[1]),
      "_vs_",
      sanitize_name(selected_groups[2])
    )
  } else {
    output_prefix <- paste0(
      "species_overlap_all_",
      number_of_groups,
      "_groups"
    )
  }
} else {
  output_prefix <- opt$output_prefix
}

dir.create(
  dirname(output_prefix),
  recursive = TRUE,
  showWarnings = FALSE
)

group_summary <- paste0(
  selected_groups,
  ": ",
  as.integer(group_counts[selected_groups])
)

if (is.null(opt$subtitle)) {
  shared_label <- if (number_of_groups == 2L) {
    "Shared"
  } else {
    "Shared by all"
  }

  subtitle_text <- paste(
    c(
      group_summary,
      paste0(
        shared_label,
        ": ",
        length(shared_by_all)
      ),
      paste0(
        "Union: ",
        length(union_species)
      )
    ),
    collapse = " | "
  )
} else {
  subtitle_text <- opt$subtitle
}

caption_text <- paste0(
  "Presence was defined as abundance > ",
  format(
    opt$presence_threshold,
    scientific = FALSE,
    trim = TRUE
  ),
  " in at least one sample from the corresponding group. ",
  "Numbers inside regions represent disjoint species counts."
)

message("============================================================")
message(
  "MTD species overlap diagrams | VERSION ",
  SCRIPT_VERSION,
  " | MULTI-GROUP"
)
message("Input matrix: ", opt$input)
message("Samplesheet: ", opt$samplesheet)
message("Taxon column: ", taxon_col)
message("Selected groups: ", paste(selected_groups, collapse = ", "))
message("Presence threshold: > ", opt$presence_threshold)
message("Output prefix: ", output_prefix)
message("============================================================")

for (group_name in selected_groups) {
  message(
    "[GROUP] ",
    group_name,
    " | samples=",
    length(group_samples[[group_name]]),
    " | species=",
    length(species_sets[[group_name]])
  )
  message(
    "        ",
    paste(
      group_samples[[group_name]],
      collapse = ", "
    )
  )
}

message(
  "[SUMMARY] Shared by all groups: ",
  length(shared_by_all)
)

message(
  "[SUMMARY] Union: ",
  length(union_species)
)

summary_rows <- do.call(
  rbind,
  lapply(
    selected_groups,
    function(group_name) {
      data.frame(
        record_type = "group",
        group = group_name,
        value = length(species_sets[[group_name]]),
        sample_count = length(group_samples[[group_name]]),
        samples = paste(
          group_samples[[group_name]],
          collapse = ";"
        ),
        presence_threshold = opt$presence_threshold,
        input_file = opt$input,
        stringsAsFactors = FALSE
      )
    }
  )
)

overall_rows <- data.frame(
  record_type = c(
    "shared_by_all_groups",
    "union",
    "detected_across_selected_samples"
  ),
  group = NA_character_,
  value = c(
    length(shared_by_all),
    length(union_species),
    length(detected_selected)
  ),
  sample_count = NA_integer_,
  samples = NA_character_,
  presence_threshold = opt$presence_threshold,
  input_file = opt$input,
  stringsAsFactors = FALSE
)

summary_df <- rbind(
  summary_rows,
  overall_rows
)

summary_file <- paste0(
  output_prefix,
  "_summary.tsv"
)

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
  for (group_name in selected_groups) {
    group_file <- paste0(
      output_prefix,
      "_",
      sanitize_name(group_name),
      "_species.txt"
    )

    write_species_file(
      species_sets[[group_name]],
      group_file
    )

    message(
      "[OK] Group species list: ",
      group_file
    )
  }

  shared_file <- paste0(
    output_prefix,
    "_shared_by_all_species.txt"
  )

  union_file <- paste0(
    output_prefix,
    "_union_species.txt"
  )

  write_species_file(
    shared_by_all,
    shared_file
  )

  write_species_file(
    union_species,
    union_file
  )

  message(
    "[OK] Shared-by-all species list: ",
    shared_file
  )

  message(
    "[OK] Union species list: ",
    union_file
  )
}

# ============================================================
# Venn diagram
# ============================================================

if (opt$plot_type %in% c("both", "venn")) {
  message(
    "[PLOT] Generating Venn diagram for ",
    number_of_groups,
    " groups..."
  )

  p_venn <- tryCatch(
    suppressMessages({
      ggVennDiagram::ggVennDiagram(
        species_sets,
        label = opt$label_type,
        label_alpha = 0,
        label_geom = "label",
        label_size = if (number_of_groups <= 4L) 5 else 3.5,
        edge_size = if (number_of_groups <= 4L) 1.1 else 0.8,
        set_size = if (number_of_groups <= 4L) 5.2 else 3.8
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
          expand = ggplot2::expansion(
            mult = c(0.10, 0.08)
          )
        ) +
        ggplot2::scale_y_continuous(
          expand = ggplot2::expansion(
            mult = c(0.10, 0.10)
          )
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
            size = 12,
            margin = ggplot2::margin(b = 20)
          ),
          plot.caption = ggplot2::element_text(
            hjust = 0.5,
            size = 9,
            color = "grey35",
            margin = ggplot2::margin(t = 18)
          ),
          legend.position = "right",
          legend.title = ggplot2::element_text(
            size = 11,
            face = "bold"
          ),
          legend.text = ggplot2::element_text(
            size = 10
          ),
          plot.margin = ggplot2::margin(
            20,
            25,
            20,
            55
          )
        )
    }),
    error = function(error) {
      stop(
        "Could not generate the multi-group Venn diagram: ",
        conditionMessage(error),
        call. = FALSE
      )
    }
  )

  save_ggplot_all(
    plot = p_venn,
    prefix = paste0(
      output_prefix,
      "_venn"
    ),
    width = opt$width,
    height = opt$height
  )
}

# ============================================================
# Euler diagram
# ============================================================

if (opt$plot_type %in% c("both", "euler")) {
  message(
    "[PLOT] Generating proportional Euler diagram for ",
    number_of_groups,
    " groups..."
  )

  euler_fit <- tryCatch(
    eulerr::euler(
      species_sets,
      shape = "ellipse"
    ),
    error = function(error) {
      stop(
        "Could not fit the multi-group Euler diagram: ",
        conditionMessage(error),
        call. = FALSE
      )
    }
  )

  euler_stress <- safe_metric(
    euler_fit,
    "stress"
  )

  euler_diag_error <- safe_metric(
    euler_fit,
    "diagError"
  )

  euler_max_region_error <- safe_max_metric(
    euler_fit,
    "regionError"
  )

  fit_warning <- if (
    is.finite(euler_diag_error) &&
    euler_diag_error > 0.05
  ) {
    paste0(
      " Warning: diagnostic error exceeds 0.05; ",
      "interpret ellipse areas cautiously."
    )
  } else {
    ""
  }

  euler_caption <- paste0(
    "Euler diagram fitted from all pairwise and higher-order ",
    "species intersections across ",
    number_of_groups,
    " groups. Stress: ",
    signif(euler_stress, 4),
    "; diagnostic error: ",
    signif(euler_diag_error, 4),
    "; maximum region error: ",
    signif(euler_max_region_error, 4),
    ". With three or more groups, an exact area-proportional ",
    "ellipse configuration may not exist; these metrics quantify ",
    "the approximation.",
    fit_warning
  )

  save_euler_all(
    fit = euler_fit,
    prefix = paste0(
      output_prefix,
      "_euler"
    ),
    width = opt$width,
    height = opt$height,
    title = opt$title,
    subtitle = subtitle_text,
    caption = euler_caption,
    label_type = opt$label_type,
    number_of_sets = number_of_groups
  )
}

message("")
message("[DONE] Finished successfully.")
message("")
