#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
  stop(
    "Usage: Rscript rerun_GO_faceted_dotplot.R ",
    "/path/to/Host_DEG"
  )
}

host_deg <- normalizePath(
  args[1],
  mustWork = TRUE
)

setwd(host_deg)

input_file <- "biological_theme_comparison_GO_results.csv"

if (!file.exists(input_file)) {
  stop(
    "Input file not found: ",
    file.path(host_deg, input_file)
  )
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

## Number of terms selected within each ontology and direction.
TOP_N <- as.integer(
  Sys.getenv("GO_TOP_N", "5")
)

if (is.na(TOP_N) || TOP_N < 1) {
  stop("GO_TOP_N must be an integer greater than zero.")
}

cat("============================================================\n")
cat("[GO DOTPLOT] Faceted publication figure\n")
cat("Working directory: ", getwd(), "\n", sep = "")
cat("Input:             ", input_file, "\n", sep = "")
cat("Top terms/panel:   ", TOP_N, "\n", sep = "")
cat("R:                 ", R.version.string, "\n", sep = "")
cat("ggplot2:           ",
    as.character(packageVersion("ggplot2")), "\n", sep = "")
cat("============================================================\n")

## ------------------------------------------------------------
## Read results
## ------------------------------------------------------------

go_df <- read.csv(
  input_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_columns <- c(
  "Cluster",
  "ONTOLOGY",
  "ID",
  "Description",
  "GeneRatio",
  "p.adjust",
  "Count"
)

missing_columns <- setdiff(
  required_columns,
  colnames(go_df)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

## ------------------------------------------------------------
## Convert GeneRatio such as 25/2513 to a numeric proportion
## ------------------------------------------------------------

ratio_to_numeric <- function(x) {

  vapply(
    strsplit(
      as.character(x),
      split = "/",
      fixed = TRUE
    ),
    function(parts) {

      if (length(parts) != 2) {
        return(NA_real_)
      }

      numerator <- suppressWarnings(
        as.numeric(parts[1])
      )

      denominator <- suppressWarnings(
        as.numeric(parts[2])
      )

      if (
        is.na(numerator) ||
        is.na(denominator) ||
        denominator == 0
      ) {
        return(NA_real_)
      }

      numerator / denominator
    },
    FUN.VALUE = numeric(1)
  )
}

go_df$GeneRatio_numeric <- ratio_to_numeric(
  go_df$GeneRatio
)

go_df$p.adjust <- suppressWarnings(
  as.numeric(go_df$p.adjust)
)

go_df$Count <- suppressWarnings(
  as.numeric(go_df$Count)
)

go_df$ONTOLOGY <- toupper(
  trimws(
    as.character(go_df$ONTOLOGY)
  )
)

go_df$Cluster <- trimws(
  as.character(go_df$Cluster)
)

go_df$Description <- trimws(
  as.character(go_df$Description)
)

## ------------------------------------------------------------
## Remove invalid rows
## ------------------------------------------------------------

go_df <- go_df[
  go_df$ONTOLOGY %in% c("BP", "CC", "MF") &
  !is.na(go_df$Cluster) &
  nzchar(go_df$Cluster) &
  !is.na(go_df$Description) &
  nzchar(go_df$Description) &
  is.finite(go_df$GeneRatio_numeric) &
  go_df$GeneRatio_numeric > 0 &
  is.finite(go_df$p.adjust) &
  go_df$p.adjust >= 0 &
  is.finite(go_df$Count) &
  go_df$Count > 0,
  ,
  drop = FALSE
]

if (nrow(go_df) == 0) {
  stop("No valid GO rows remained after filtering.")
}

## Avoid Inf if an adjusted p-value is exactly zero.

positive_pvalues <- go_df$p.adjust[
  go_df$p.adjust > 0
]

pvalue_floor <- if (length(positive_pvalues) > 0) {
  min(positive_pvalues) / 10
} else {
  .Machine$double.xmin
}

go_df$neg_log10_padj <- -log10(
  pmax(
    go_df$p.adjust,
    pvalue_floor
  )
)

## ------------------------------------------------------------
## More readable labels for the comparison direction
## ------------------------------------------------------------

go_df$Direction <- dplyr::case_when(
  go_df$Cluster == "Liver_vs_Telencephalon_UP" ~ "UP",
  go_df$Cluster == "Liver_vs_Telencephalon_DOWN" ~ "DOWN",
  TRUE ~ go_df$Cluster
)

## If the contrast was defined as Liver minus Telencephalon,
## you may replace the labels above with:
##
## UP   = "Higher in liver"
## DOWN = "Higher in telencephalon"

go_df$Direction <- factor(
  go_df$Direction,
  levels = unique(
    c(
      "UP",
      "DOWN",
      go_df$Direction
    )
  )
)

go_df$ONTOLOGY_label <- factor(
  go_df$ONTOLOGY,
  levels = c("BP", "CC", "MF"),
  labels = c(
    "Biological Process",
    "Cellular Component",
    "Molecular Function"
  )
)

## ------------------------------------------------------------
## Select top N by adjusted p-value in each panel
## ------------------------------------------------------------

plot_df <- go_df %>%
  group_by(
    ONTOLOGY,
    ONTOLOGY_label,
    Direction
  ) %>%
  arrange(
    p.adjust,
    desc(GeneRatio_numeric),
    desc(Count),
    .by_group = TRUE
  ) %>%
  slice_head(
    n = TOP_N
  ) %>%
  ungroup()

if (nrow(plot_df) == 0) {
  stop("No GO terms were selected for plotting.")
}

## ------------------------------------------------------------
## Create a unique internal term key
##
## Description alone may be repeated in different ontologies or
## directions. The suffix avoids duplicated factor levels.
## ------------------------------------------------------------

plot_df$Term_key <- paste(
  plot_df$Description,
  plot_df$ID,
  plot_df$ONTOLOGY,
  plot_df$Direction,
  sep = "|||"
)

## Within each panel, place larger GeneRatio values nearer the top.

plot_df <- plot_df %>%
  arrange(
    ONTOLOGY_label,
    Direction,
    GeneRatio_numeric,
    neg_log10_padj
  )

plot_df$Term_key <- factor(
  plot_df$Term_key,
  levels = unique(plot_df$Term_key)
)

## Wrap long GO descriptions.

wrap_term <- function(x, width = 42) {

  descriptions <- sub(
    "\\|\\|\\|.*$",
    "",
    x
  )

  vapply(
    descriptions,
    function(term) {
      paste(
        strwrap(
          term,
          width = width
        ),
        collapse = "\n"
      )
    },
    FUN.VALUE = character(1)
  )
}

## ------------------------------------------------------------
## Build figure
## ------------------------------------------------------------

p <- ggplot(
  plot_df,
  aes(
    x = GeneRatio_numeric,
    y = Term_key
  )
) +
  geom_point(
    aes(
      size = Count,
      colour = neg_log10_padj
    ),
    alpha = 0.90
  ) +
  facet_grid(
    rows = vars(ONTOLOGY_label),
    cols = vars(Direction),
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  scale_y_discrete(
    labels = wrap_term
  ) +
  scale_x_continuous(
    labels = function(x) {
      sprintf("%.3f", x)
    },
    expand = expansion(
      mult = c(0.02, 0.10)
    )
  ) +
  scale_size(
    range = c(2.5, 8),
    name = "Gene count"
  ) +
  scale_colour_viridis_c(
    option = "C",
    direction = 1,
    name = expression(
      -log[10]("adjusted p-value")
    )
  ) +
  labs(
    title = "Gene Ontology enrichment by direction",
    subtitle = paste0(
      "Top ",
      TOP_N,
      " terms ranked by adjusted p-value ",
      "within each ontology and direction"
    ),
    x = "Gene ratio",
    y = NULL
  ) +
  theme_bw(
    base_size = 11
  ) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),

    panel.spacing.x = grid::unit(
      1.2,
      "lines"
    ),

    panel.spacing.y = grid::unit(
      0.9,
      "lines"
    ),

    strip.background = element_rect(
      fill = "grey95",
      colour = "grey50",
      linewidth = 0.4
    ),

    strip.text = element_text(
      face = "bold",
      size = 10
    ),

    strip.text.y.left = element_text(
      angle = 0
    ),

    axis.text.x = element_text(
      size = 9
    ),

    axis.text.y = element_text(
      size = 8,
      colour = "black"
    ),

    axis.title.x = element_text(
      face = "bold"
    ),

    legend.position = "right",

    plot.title = element_text(
      face = "bold",
      size = 15,
      hjust = 0.5
    ),

    plot.subtitle = element_text(
      size = 10,
      hjust = 0.5
    ),

    plot.margin = margin(
      10,
      15,
      10,
      10
    )
  )

## ------------------------------------------------------------
## Save figure and selected data
## ------------------------------------------------------------

pdf_file <- paste0(
  "biological_theme_comparison_GO_faceted_dotplot_top",
  TOP_N,
  ".pdf"
)

tiff_file <- paste0(
  "biological_theme_comparison_GO_faceted_dotplot_top",
  TOP_N,
  ".tiff"
)

data_file <- paste0(
  "biological_theme_comparison_GO_faceted_dotplot_top",
  TOP_N,
  "_data.csv"
)

ggsave(
  filename = pdf_file,
  plot = p,
  width = 14,
  height = 13,
  units = "in",
  limitsize = FALSE
)

ggsave(
  filename = tiff_file,
  plot = p,
  width = 14,
  height = 13,
  units = "in",
  dpi = 600,
  compression = "lzw",
  limitsize = FALSE
)

write.csv(
  plot_df[
    ,
    c(
      "Cluster",
      "Direction",
      "ONTOLOGY",
      "ID",
      "Description",
      "GeneRatio",
      "GeneRatio_numeric",
      "p.adjust",
      "neg_log10_padj",
      "Count"
    )
  ],
  data_file,
  row.names = FALSE
)

cat("\n============================================================\n")
cat("[OK] Faceted GO dotplot generated\n")
cat("[OK] PDF:  ", file.path(getwd(), pdf_file), "\n", sep = "")
cat("[OK] TIFF: ", file.path(getwd(), tiff_file), "\n", sep = "")
cat("[OK] Data: ", file.path(getwd(), data_file), "\n", sep = "")
cat("============================================================\n")
