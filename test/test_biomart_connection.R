#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(biomaRt)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  stop(
    "Usage:\n",
    "Rscript test_biomart_connection.R <host_counts.txt> <HostSpecies.csv> <taxid> <output_cache.csv>\n\n",
    "Example:\n",
    "Rscript test_biomart_connection.R /home/me/projeto_morcego/MTD_res/host_counts.txt /home/me/MTD/HostSpecies.csv 59463 /home/me/projeto_morcego/MTD_res/Host_DEG/gene_ID_cache_online.csv"
  )
}

host_counts_file <- args[1]
host_species_file <- args[2]
taxid <- args[3]
output_cache <- args[4]

cat("------------------------------------------------------------\n")
cat("BioMart diagnostic and annotation fetcher\n")
cat("host_counts:", host_counts_file, "\n")
cat("HostSpecies:", host_species_file, "\n")
cat("TaxID:", taxid, "\n")
cat("Output cache:", output_cache, "\n")
cat("------------------------------------------------------------\n\n")

safe_run <- function(expr, label) {
  cat("\n====================\n")
  cat(label, "\n")
  cat("====================\n")

  tryCatch(
    {
      result <- eval.parent(substitute(expr))
      cat("[OK]", label, "\n")
      result
    },
    error = function(e) {
      cat("[FAIL]", label, "\n")
      cat("Reason:", conditionMessage(e), "\n")
      NULL
    }
  )
}

cat("R version:\n")
print(R.version.string)

cat("\nbiomaRt version:\n")
print(as.character(packageVersion("biomaRt")))

cat("\nSession info short:\n")
print(sessionInfo()[c("R.version", "platform", "locale")])

# ------------------------------------------------------------
# Read HostSpecies.csv and dataset
# ------------------------------------------------------------
host_sp <- read.csv(host_species_file, header = TRUE, check.names = FALSE)

if (!("Taxon_ID" %in% names(host_sp))) {
  stop("HostSpecies.csv does not contain column Taxon_ID")
}

row <- host_sp[as.character(host_sp$Taxon_ID) == as.character(taxid), ]

if (nrow(row) == 0) {
  stop("TaxID ", taxid, " not found in HostSpecies.csv")
}

dataset_use <- as.character(row[1, 2])

cat("\nDataset detected from HostSpecies.csv:\n")
print(dataset_use)

if (is.na(dataset_use) || dataset_use == "") {
  stop("Dataset for TaxID ", taxid, " is empty/NA in HostSpecies.csv")
}

# ------------------------------------------------------------
# Read genes from host_counts.txt
# ------------------------------------------------------------
cts <- read.table(
  host_counts_file,
  row.names = 1,
  sep = "\t",
  header = TRUE,
  quote = "",
  check.names = FALSE
)

gene_len <- cts["Length"]
colnames(gene_len) <- "gene_length"

# Same logic as DEG_Anno_Plot.R: drop first 5 featureCounts annotation columns
cts_counts <- cts[, -c(1:5), drop = FALSE]
cts_counts <- cts_counts[rowSums(cts_counts[-1]) > 0, , drop = FALSE]

genes <- rownames(cts_counts)

cat("\nGenes to query after count filtering:", length(genes), "\n")
cat("First genes:\n")
print(head(genes, 10))

# Use a small subset first to test getBM
test_genes <- head(genes, 20)

# ------------------------------------------------------------
# Basic system-level URL diagnostics
# ------------------------------------------------------------
urls <- c(
  "https://www.ensembl.org",
  "https://useast.ensembl.org",
  "https://asia.ensembl.org"
)

cat("\n------------------------------------------------------------\n")
cat("System curl diagnostics\n")
cat("------------------------------------------------------------\n")

for (u in urls) {
  cat("\nURL:", u, "\n")
  cmd <- paste("curl -I -L --connect-timeout 10 --max-time 20", shQuote(u), "2>&1 | head -n 20")
  out <- system(cmd, intern = TRUE)
  cat(paste(out, collapse = "\n"), "\n")
}

# ------------------------------------------------------------
# Try different connection strategies
# ------------------------------------------------------------
connectors <- list(
  list(
    name = "useEnsembl mirror www",
    fun = function() biomaRt::useEnsembl(
      biomart = "genes",
      dataset = dataset_use,
      mirror = "www"
    )
  ),
  list(
    name = "useEnsembl mirror useast",
    fun = function() biomaRt::useEnsembl(
      biomart = "genes",
      dataset = dataset_use,
      mirror = "useast"
    )
  ),
  list(
    name = "useEnsembl mirror asia",
    fun = function() biomaRt::useEnsembl(
      biomart = "genes",
      dataset = dataset_use,
      mirror = "asia"
    )
  ),
  list(
    name = "useMart host www.ensembl.org",
    fun = function() {
      mart <- biomaRt::useMart(
        biomart = "ENSEMBL_MART_ENSEMBL",
        host = "https://www.ensembl.org"
      )
      biomaRt::useDataset(dataset = dataset_use, mart = mart)
    }
  ),
  list(
    name = "useMart host useast.ensembl.org",
    fun = function() {
      mart <- biomaRt::useMart(
        biomart = "ENSEMBL_MART_ENSEMBL",
        host = "https://useast.ensembl.org"
      )
      biomaRt::useDataset(dataset = dataset_use, mart = mart)
    }
  ),
  list(
    name = "useMart host asia.ensembl.org",
    fun = function() {
      mart <- biomaRt::useMart(
        biomart = "ENSEMBL_MART_ENSEMBL",
        host = "https://asia.ensembl.org"
      )
      biomaRt::useDataset(dataset = dataset_use, mart = mart)
    }
  )
)

working_mart <- NULL
working_method <- NA_character_

for (connector in connectors) {
  cat("\n------------------------------------------------------------\n")
  cat("Trying:", connector$name, "\n")
  cat("------------------------------------------------------------\n")

  mart <- tryCatch(
    connector$fun(),
    error = function(e) {
      cat("[CONNECT FAIL]", connector$name, "\n")
      cat("Reason:", conditionMessage(e), "\n")
      NULL
    }
  )

  if (is.null(mart)) {
    next
  }

  cat("[CONNECT OK]", connector$name, "\n")

  # Test attributes
  attrs <- tryCatch(
    {
      a <- biomaRt::listAttributes(mart)
      cat("Attributes available:", nrow(a), "\n")
      a
    },
    error = function(e) {
      cat("[ATTRIBUTES FAIL]", conditionMessage(e), "\n")
      NULL
    }
  )

  if (is.null(attrs)) {
    next
  }

  required_attrs <- c(
    "external_gene_name",
    "ensembl_gene_id",
    "chromosome_name",
    "start_position",
    "end_position",
    "strand",
    "gene_biotype",
    "description"
  )

  missing_attrs <- setdiff(required_attrs, attrs$name)

  if (length(missing_attrs) > 0) {
    cat("[ATTRIBUTES MISSING]", paste(missing_attrs, collapse = ", "), "\n")
    next
  }

  # Test getBM on 20 genes
  test <- tryCatch(
    {
      biomaRt::getBM(
        filters = "ensembl_gene_id",
        attributes = required_attrs,
        values = test_genes,
        mart = mart
      )
    },
    error = function(e) {
      cat("[getBM TEST FAIL]", conditionMessage(e), "\n")
      NULL
    }
  )

  if (is.null(test)) {
    next
  }

  cat("[getBM TEST OK] Rows returned:", nrow(test), "\n")
  print(head(test))

  working_mart <- mart
  working_method <- connector$name
  break
}

if (is.null(working_mart)) {
  cat("\n============================================================\n")
  cat("No live BioMart connection worked.\n")
  cat("Recommendation: keep using local gene_ID_cache.csv for now.\n")
  cat("============================================================\n")
  quit(status = 2)
}

cat("\n============================================================\n")
cat("Working BioMart method:", working_method, "\n")
cat("Now fetching annotations for all genes...\n")
cat("============================================================\n")

attributes_use <- c(
  "external_gene_name",
  "ensembl_gene_id",
  "chromosome_name",
  "start_position",
  "end_position",
  "strand",
  "gene_biotype",
  "description"
)

# Query in chunks to reduce timeout risk
chunk_size <- 500
chunks <- split(genes, ceiling(seq_along(genes) / chunk_size))

all_gene_ID <- list()

for (i in seq_along(chunks)) {
  cat("Fetching chunk", i, "of", length(chunks), "genes:", length(chunks[[i]]), "\n")

  chunk_result <- tryCatch(
    {
      biomaRt::getBM(
        filters = "ensembl_gene_id",
        attributes = attributes_use,
        values = chunks[[i]],
        mart = working_mart
      )
    },
    error = function(e) {
      cat("[CHUNK FAIL]", i, conditionMessage(e), "\n")
      NULL
    }
  )

  if (!is.null(chunk_result) && nrow(chunk_result) > 0) {
    all_gene_ID[[length(all_gene_ID) + 1]] <- chunk_result
  }

  Sys.sleep(1)
}

if (length(all_gene_ID) == 0) {
  stop("No annotation rows were returned by BioMart.")
}

gene_ID <- do.call(rbind, all_gene_ID)
gene_ID <- unique(gene_ID)

names(gene_ID)[names(gene_ID) == "external_gene_name"] <- "gene_name"

gene_ID <- merge(
  gene_ID,
  gene_len,
  by.x = "ensembl_gene_id",
  by.y = "row.names",
  all.x = TRUE
)

gene_ID <- gene_ID[!is.na(gene_ID$ensembl_gene_id), ]
gene_ID <- gene_ID[!duplicated(gene_ID$ensembl_gene_id), ]

cat("\nAnnotations retrieved:", nrow(gene_ID), "\n")
cat("Coverage:", round(100 * nrow(gene_ID) / length(genes), 2), "%\n")

write.csv(
  gene_ID,
  output_cache,
  row.names = FALSE,
  quote = TRUE
)

cat("\nSaved online BioMart cache:\n")
cat(output_cache, "\n")

cat("\nDone.\n")
