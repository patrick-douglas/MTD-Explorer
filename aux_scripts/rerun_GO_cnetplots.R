#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
  stop(
    "Uso: Rscript /tmp/test_GO_cnetplot_by_ontology.R ",
    "/caminho/para/Host_DEG"
  )
}

host_deg <- normalizePath(args[1], mustWork = TRUE)
setwd(host_deg)

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(enrichplot)
  library(ggplot2)
})

input_file <- "biological_theme_comparison_GO_results.csv"

if (!file.exists(input_file)) {
  stop("Arquivo não encontrado: ", file.path(host_deg, input_file))
}

show_n <- as.integer(Sys.getenv("CNET_SHOW", "10"))

if (is.na(show_n) || show_n < 1) {
  stop("CNET_SHOW deve ser um número inteiro maior que zero.")
}

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
  "geneID",
  "Count"
)

missing_columns <- setdiff(required_columns, colnames(go_df))

if (length(missing_columns) > 0) {
  stop(
    "Colunas ausentes: ",
    paste(missing_columns, collapse = ", ")
  )
}

cat("============================================================\n")
cat("[TEST] GO cnetplot separated by ontology\n")
cat("Working directory: ", getwd(), "\n", sep = "")
cat("Input:             ", input_file, "\n", sep = "")
cat("showCategory:      ", show_n, "\n", sep = "")
cat("R:                 ", R.version.string, "\n", sep = "")
cat(
  "clusterProfiler:   ",
  as.character(packageVersion("clusterProfiler")),
  "\n",
  sep = ""
)
cat(
  "enrichplot:        ",
  as.character(packageVersion("enrichplot")),
  "\n",
  sep = ""
)
cat("============================================================\n")

available_ontologies <- unique(
  trimws(as.character(go_df$ONTOLOGY))
)

available_ontologies <- available_ontologies[
  !is.na(available_ontologies) &
  nzchar(available_ontologies)
]

cat(
  "\n[INFO] Ontologies found: ",
  paste(available_ontologies, collapse = ", "),
  "\n",
  sep = ""
)

preferred_order <- c("BP", "CC", "MF")

ontologies <- c(
  intersect(preferred_order, available_ontologies),
  setdiff(available_ontologies, preferred_order)
)

successful_outputs <- character()
failed_outputs <- character()
go_cnet_plots <- list()

for (ontology_name in ontologies) {

  cat("\n============================================================\n")
  cat("[ONTOLOGY] ", ontology_name, "\n", sep = "")
  cat("============================================================\n")

  go_sub <- go_df[
    trimws(as.character(go_df$ONTOLOGY)) == ontology_name,
    ,
    drop = FALSE
  ]

  go_sub <- go_sub[
    !is.na(go_sub$Cluster) &
    nzchar(trimws(as.character(go_sub$Cluster))) &
    !is.na(go_sub$Description) &
    nzchar(trimws(as.character(go_sub$Description))) &
    !is.na(go_sub$geneID) &
    nzchar(trimws(as.character(go_sub$geneID))),
    ,
    drop = FALSE
  ]

  if (nrow(go_sub) == 0) {
    cat("[SKIP] No valid rows for ontology ", ontology_name, "\n", sep = "")
    next
  }

  cluster_levels <- unique(
    as.character(go_sub$Cluster)
  )

  go_sub$Cluster <- factor(
    as.character(go_sub$Cluster),
    levels = cluster_levels
  )

  go_sub$Description <- as.character(go_sub$Description)
  go_sub$geneID <- as.character(go_sub$geneID)

  gene_clusters <- setNames(
    lapply(
      cluster_levels,
      function(cluster_name) {

        gene_strings <- go_sub$geneID[
          as.character(go_sub$Cluster) == cluster_name
        ]

        genes <- unlist(
          strsplit(
            gene_strings,
            split = "/",
            fixed = TRUE
          ),
          use.names = FALSE
        )

        genes <- trimws(genes)

        unique(
          genes[
            !is.na(genes) &
            nzchar(genes)
          ]
        )
      }
    ),
    cluster_levels
  )

  cgo_ontology <- methods::new(
    "compareClusterResult",
    compareClusterResult = go_sub,
    geneClusters = gene_clusters,
    fun = "enrichGO"
  )

  denominator <- sub(
    "^.*/",
    "",
    as.character(go_sub$GeneRatio)
  )

  denominator_summary <- unique(
    data.frame(
      Cluster = as.character(go_sub$Cluster),
      ONTOLOGY = ontology_name,
      GeneRatio_denominator = denominator,
      stringsAsFactors = FALSE
    )
  )

  cat("[INFO] Rows: ", nrow(go_sub), "\n", sep = "")
  cat("[INFO] GeneRatio denominators:\n")
  print(denominator_summary)

  cat(
    "[RUN] Testing fortify(showCategory = ",
    show_n,
    ")...\n",
    sep = ""
  )

  fortified <- tryCatch(
    {
      ggplot2::fortify(
        cgo_ontology,
        showCategory = show_n,
        includeAll = TRUE
      )
    },
    error = function(e) {
      cat(
        "[FAIL] fortify: ",
        conditionMessage(e),
        "\n",
        sep = ""
      )

      NULL
    }
  )

  if (is.null(fortified)) {
    failed_outputs <- c(failed_outputs, ontology_name)
    next
  }

  cluster_na <- sum(is.na(fortified$Cluster))

  cat(
    "[INFO] Rows after fortify: ",
    nrow(fortified),
    "\n",
    sep = ""
  )

  cat(
    "[INFO] NA in Cluster after fortify: ",
    cluster_na,
    "\n",
    sep = ""
  )

  if (cluster_na > 0) {

    diagnostic_file <- paste0(
      "GO_cnetplot_",
      ontology_name,
      "_fortify_NA_rows.csv"
    )

    write.csv(
      fortified[is.na(fortified$Cluster), , drop = FALSE],
      diagnostic_file,
      row.names = FALSE
    )

    cat(
      "[FAIL] fortify still generated NA values.\n",
      "[FAIL] Diagnostic: ",
      diagnostic_file,
      "\n",
      sep = ""
    )

    failed_outputs <- c(failed_outputs, ontology_name)
    next
  }

  output_file <- paste0(
    "biological_theme_comparison_GO_net_",
    ontology_name,
    ".pdf"
  )

  cat("[RUN] Creating cnetplot...\n")

  plot_result <- tryCatch(
    {
      cnetplot(
        cgo_ontology,
        showCategory = show_n,
        node_label = "category",
        cex_category = 1.10,
        cex_gene = 0.45,
        cex_label_category = 0.85
      ) +
        ggplot2::labs(
          title = switch(
            ontology_name,
            BP = "GO Biological Process (BP)",
            CC = "GO Cellular Component (CC)",
            MF = "GO Molecular Function (MF)",
            paste("GO", ontology_name)
          )
        ) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(
            hjust = 0.5,
            face = "bold",
            size = 14
          ),
          legend.position = "right"
        )
    },
    error = function(e) {
      cat(
        "[FAIL] cnetplot: ",
        conditionMessage(e),
        "\n",
        sep = ""
      )

      NULL
    }
  )

  if (is.null(plot_result)) {
    failed_outputs <- c(failed_outputs, ontology_name)
    next
  }

  cat("[RUN] Saving: ", output_file, "\n", sep = "")

  save_result <- tryCatch(
    {
      ggsave(
        filename = output_file,
        plot = plot_result,
        limitsize = FALSE,
        width = 14,
        height = 11
      )

      TRUE
    },
    error = function(e) {
      cat(
        "[FAIL] ggsave: ",
        conditionMessage(e),
        "\n",
        sep = ""
      )

      FALSE
    }
  )

  if (save_result && file.exists(output_file)) {

    successful_outputs <- c(
      successful_outputs,
      output_file
    )

    go_cnet_plots[[ontology_name]] <- plot_result

    cat("[OK] Generated: ", output_file, "\n", sep = "")

  } else {

    failed_outputs <- c(
      failed_outputs,
      ontology_name
    )
  }
}

## ------------------------------------------------------------
## Create the original expected multi-page PDF
## ------------------------------------------------------------

if (length(go_cnet_plots) > 0) {

  combined_file <- "biological_theme_comparison_GO_net.pdf"
  pdf_is_open <- FALSE

  combined_success <- tryCatch(
    {
      grDevices::pdf(
        file = combined_file,
        width = 14,
        height = 11,
        onefile = TRUE
      )

      pdf_is_open <- TRUE

      for (ontology_name in names(go_cnet_plots)) {
        print(go_cnet_plots[[ontology_name]])
      }

      grDevices::dev.off()
      pdf_is_open <- FALSE

      TRUE
    },
    error = function(e) {

      if (
        pdf_is_open &&
        grDevices::dev.cur() > 1
      ) {
        grDevices::dev.off()
      }

      cat(
        "[FAIL] Multi-page PDF: ",
        conditionMessage(e),
        "\n",
        sep = ""
      )

      FALSE
    }
  )

  if (
    combined_success &&
    file.exists(combined_file)
  ) {
    cat(
      "[OK] Generated multi-page PDF: ",
      combined_file,
      "\n",
      sep = ""
    )
  }
}

cat("\n============================================================\n")
cat("[SUMMARY]\n")
cat("============================================================\n")

if (length(successful_outputs) > 0) {
  cat("[OK] Generated files:\n")

  for (output_file in successful_outputs) {
    cat("  ", file.path(getwd(), output_file), "\n", sep = "")
  }
} else {
  cat("[WARNING] No PDF was generated.\n")
}

if (length(failed_outputs) > 0) {
  cat(
    "[WARNING] Failed ontologies: ",
    paste(unique(failed_outputs), collapse = ", "),
    "\n",
    sep = ""
  )
}

cat("============================================================\n")
