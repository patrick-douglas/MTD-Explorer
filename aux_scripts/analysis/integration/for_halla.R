#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (!length(args) %in% c(2L, 3L)) {
  stop(
    "Usage: Rscript for_halla.R ",
    "<ssGSEA_scores.gct> <samplesheet.csv> [metadata.csv]"
  )
}

scores_file <- normalizePath(
  args[1],
  mustWork = TRUE
)

samplesheet_file <- normalizePath(
  args[2],
  mustWork = TRUE
)

metadata_supplied <- (
  length(args) == 3L &&
    !is.na(args[3]) &&
    nzchar(trimws(args[3]))
)

metadata_file <- NULL

if (metadata_supplied) {
  metadata_file <- normalizePath(
    args[3],
    mustWork = TRUE
  )
}

output_root <- dirname(
  dirname(scores_file)
)

halla_directory <- file.path(
  output_root,
  "halla"
)

dir.create(
  halla_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# Read ssGSEA score matrix
# ------------------------------------------------------------

score <- read.table(
  scores_file,
  row.names = 1,
  sep = "\t",
  header = TRUE,
  skip = 2,
  check.names = FALSE,
  quote = "",
  comment.char = ""
)

score <- as.matrix(score)

suppressWarnings(
  storage.mode(score) <- "numeric"
)

# ------------------------------------------------------------
# Read samplesheet or extended metadata
# ------------------------------------------------------------

samplesheet <- read.csv(
  samplesheet_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)

if (ncol(samplesheet) < 2L) {
  stop(
    "Samplesheet must contain at least sample_name and group."
  )
}

if (metadata_supplied) {
  metadata <- read.csv(
    metadata_file,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("", "NA")
  )

  message(
    "[FOR_HALLA] Extended metadata supplied: ",
    metadata_file
  )
} else {
  metadata <- samplesheet[, 1:2, drop = FALSE]

  names(metadata)[1:2] <- c(
    "sample_name",
    "group"
  )

  message(
    "[FOR_HALLA] No extended metadata supplied; ",
    "ssGSEA scores will be preserved without batch removal."
  )
}

required_columns <- c(
  "sample_name",
  "group"
)

missing_columns <- setdiff(
  required_columns,
  names(metadata)
)

if (length(missing_columns) > 0L) {
  stop(
    "Metadata is missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

metadata$sample_name <- trimws(
  as.character(metadata$sample_name)
)

metadata$group <- factor(
  trimws(
    as.character(metadata$group)
  )
)

metadata <- metadata[
  !is.na(metadata$sample_name) &
    metadata$sample_name != "",
  ,
  drop = FALSE
]

metadata <- metadata[
  !duplicated(metadata$sample_name),
  ,
  drop = FALSE
]

missing_score_samples <- setdiff(
  metadata$sample_name,
  colnames(score)
)

if (length(missing_score_samples) > 0L) {
  warning(
    "Metadata samples absent from ssGSEA scores: ",
    paste(missing_score_samples, collapse = ", ")
  )
}

missing_metadata_samples <- setdiff(
  colnames(score),
  metadata$sample_name
)

if (length(missing_metadata_samples) > 0L) {
  warning(
    "ssGSEA samples absent from metadata: ",
    paste(missing_metadata_samples, collapse = ", ")
  )
}

common_samples <- colnames(score)[
  colnames(score) %in% metadata$sample_name
]

if (length(common_samples) < 2L) {
  stop(
    "Fewer than two samples matched between ",
    "ssGSEA scores and metadata."
  )
}

score <- score[
  ,
  common_samples,
  drop = FALSE
]

metadata <- metadata[
  match(
    common_samples,
    metadata$sample_name
  ),
  ,
  drop = FALSE
]

if (!identical(
  colnames(score),
  metadata$sample_name
)) {
  stop(
    "Could not align ssGSEA score columns with metadata rows."
  )
}

# ------------------------------------------------------------
# Preserve biological group and remove repeated-subject effect
# ------------------------------------------------------------

score_adjusted <- score

if (metadata_supplied) {

  biological_design <- model.matrix(
    ~ group,
    data = metadata
  )

  subject_candidates <- c(
    "animal",
    "animal_id",
    "subject",
    "subject_id",
    "individual",
    "individual_id",
    "participant",
    "participant_id",
    "patient",
    "patient_id",
    "pair",
    "pair_id"
  )

  metadata_names_lower <- tolower(
    names(metadata)
  )

  subject_column <- NULL

  for (candidate in subject_candidates) {

    candidate_index <- match(
      candidate,
      metadata_names_lower,
      nomatch = 0L
    )

    if (candidate_index == 0L) {
      next
    }

    candidate_name <- names(metadata)[candidate_index]

    candidate_values <- trimws(
      as.character(
        metadata[[candidate_name]]
      )
    )

    valid_values <- candidate_values[
      !is.na(candidate_values) &
        candidate_values != ""
    ]

    if (
      length(unique(valid_values)) > 1L &&
      anyDuplicated(valid_values) > 0L
    ) {
      subject_column <- candidate_name
      break
    }
  }

  batch_factor <- NULL

  if (!is.null(subject_column)) {
    batch_factor <- factor(
      metadata[[subject_column]]
    )
  }

  excluded_columns <- c(
    "sample_name",
    "group",
    subject_column
  )

  covariate_candidates <- setdiff(
    names(metadata),
    excluded_columns
  )

  numeric_covariates <- list()

  for (column_name in covariate_candidates) {

    raw_values <- metadata[[column_name]]

    if (is.numeric(raw_values) || is.integer(raw_values)) {
      numeric_values <- as.numeric(raw_values)
    } else {
      character_values <- trimws(
        as.character(raw_values)
      )

      numeric_values <- suppressWarnings(
        as.numeric(character_values)
      )

      valid_conversion <- (
        is.na(character_values) |
          character_values == "" |
          !is.na(numeric_values)
      )

      if (!all(valid_conversion)) {
        numeric_values <- NULL
      }
    }

    if (!is.null(numeric_values)) {
      observed_values <- numeric_values[
        !is.na(numeric_values)
      ]

      if (length(unique(observed_values)) > 1L) {
        numeric_covariates[[column_name]] <- numeric_values
      }
    }
  }

  covariate_matrix <- NULL

  if (length(numeric_covariates) > 0L) {
    covariate_matrix <- as.matrix(
      as.data.frame(
        numeric_covariates,
        check.names = FALSE
      )
    )
  }

  adjustment_arguments <- list(
    x = score,
    design = biological_design
  )

  if (!is.null(batch_factor)) {
    adjustment_arguments$batch <- batch_factor
  }

  if (!is.null(covariate_matrix)) {
    adjustment_arguments$covariates <- covariate_matrix
  }

  if (
    !is.null(batch_factor) ||
    !is.null(covariate_matrix)
  ) {
    score_adjusted <- do.call(
      limma::removeBatchEffect,
      adjustment_arguments
    )
  }

  message(
    "[FOR_HALLA] Biological group preserved: group"
  )

  message(
    "[FOR_HALLA] Repeated-measures block: ",
    if (is.null(subject_column)) {
      "none"
    } else {
      subject_column
    }
  )

  message(
    "[FOR_HALLA] Numeric covariates removed: ",
    if (length(numeric_covariates) == 0L) {
      "none"
    } else {
      paste(
        names(numeric_covariates),
        collapse = ", "
      )
    }
  )
}

output_file <- file.path(
  halla_directory,
  "Host_score.txt"
)

write.table(
  score_adjusted,
  output_file,
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

message(
  "[FOR_HALLA] Host pathway matrix saved: ",
  output_file
)
