#!/usr/bin/env Rscript

# ============================================================
# EV.volcano_with_ensembl_symbols.R
#
# One-step script:
#   1) Reads a DESeq2-like differential expression/abundance table.
#   2) Detects Ensembl gene IDs in the selected ID column.
#   3) Optionally queries Ensembl REST to convert IDs into gene symbols.
#   4) Writes a converted/proof table keeping the original ID column.
#   5) Creates a volcano plot using gene symbols when available.
#
# Required input columns for volcano:
#   log2FoldChange
#   padj
# ============================================================

# -----------------------------
# Help message
# -----------------------------
print_help <- function() {
    cat("\nEV volcano plot with optional Ensembl ID -> gene symbol conversion\n\n")
    cat("Usage:\n")
    cat("  Rscript EV.volcano_with_ensembl_symbols.R --de_results FILE.csv [options]\n\n")

    cat("Required:\n")
    cat("  --de_results FILE\n")
    cat("      Input table containing differential expression/abundance results.\n")
    cat("      Required columns: log2FoldChange and padj.\n\n")

    cat("ID / gene symbol options:\n")
    cat("  --id_column VALUE\n")
    cat("      Column containing feature IDs/names. Use a 1-based column number or a column name.\n")
    cat("      For convenience, 0 is also accepted and means the first column.\n")
    cat("      Default: 1\n\n")

    cat("  --convert auto|yes|no\n")
    cat("      auto: convert only values that look like Ensembl gene IDs.\n")
    cat("      yes : force conversion attempt for Ensembl-like IDs.\n")
    cat("      no  : do not query Ensembl; use the input labels directly.\n")
    cat("      Default: auto\n\n")

    cat("  --xref_fallback\n")
    cat("      For IDs without display_name from /lookup/id, try /xrefs/id as fallback.\n")
    cat("      This is slower, because it queries one ID at a time.\n\n")

    cat("  --label_unmapped_ensembl\n")
    cat("      By default, Ensembl IDs that could not be converted are kept in the output table\n")
    cat("      but are not used as labels in the volcano plot. Use this flag to label them too.\n\n")

    cat("  --keep_version\n")
    cat("      Keep version suffixes such as ENSG000001.5 during lookup.\n")
    cat("      Default behavior strips version suffixes before querying.\n\n")

    cat("  --cache FILE.csv\n")
    cat("      Local cache file for Ensembl ID -> symbol mappings.\n")
    cat("      Default: ensembl_symbol_cache.csv\n\n")

    cat("  --server URL\n")
    cat("      Ensembl REST server.\n")
    cat("      Default: https://rest.ensembl.org\n\n")

    cat("  --batch_size N\n")
    cat("      Number of IDs per POST /lookup/id request. Maximum recommended: 1000.\n")
    cat("      Default: 1000\n\n")

    cat("Output options:\n")
    cat("  --output FILE.pdf|FILE.png\n")
    cat("      Output volcano plot file. Extension controls the device.\n")
    cat("      Default: <input>_volcano.pdf\n\n")

    cat("  --output_table FILE.csv\n")
    cat("      Converted/proof table. Keeps the original IDs and adds mapping columns.\n")
    cat("      Default: <input>_gene_symbols_table.csv\n\n")

    cat("Input format options:\n")
    cat("  --sep auto|comma|tab|semicolon\n")
    cat("      Table separator. Default: auto\n\n")

    cat("Plot label options:\n")
    cat("  --labels\n")
    cat("      Show labels for all significant features.\n\n")

    cat("  --label_top N\n")
    cat("      Show labels only for the top N features with the lowest padj.\n")
    cat("      Overrides --labels.\n\n")

    cat("Threshold options:\n")
    cat("  --padj VALUE\n")
    cat("      Adjusted p-value / FDR cutoff. Default: 0.05\n\n")

    cat("  --logfc VALUE\n")
    cat("      Absolute log2 fold-change cutoff. Default: 2\n\n")

    cat("Title options:\n")
    cat("  --group_label1 TEXT\n")
    cat("      First group shown in the plot title. Default: Group1\n\n")

    cat("  --group_label2 TEXT\n")
    cat("      Second group shown in the plot title. Default: Group2\n\n")

    cat("  --subtitle TEXT\n")
    cat("      Plot subtitle. Default: Differential abundance/expression\n\n")

    cat("Plot size options:\n")
    cat("  --width VALUE\n")
    cat("      Plot width in inches. Default: 10\n\n")

    cat("  --height VALUE\n")
    cat("      Plot height in inches. Default: 8\n\n")

    cat("  --dpi VALUE\n")
    cat("      DPI used for PNG output. Default: 300\n\n")

    cat("Package options:\n")
    cat("  --no_install\n")
    cat("      Do not try to install missing R packages automatically.\n\n")

    cat("Examples:\n\n")
    cat("  1) Convert Ensembl IDs automatically, save table, plot top 30 labels:\n")
    cat("     Rscript EV.volcano_with_ensembl_symbols.R \\\n")
    cat("       --de_results host_counts_Liver_vs_Telencephalon.csv \\\n")
    cat("       --label_top 30 \\\n")
    cat("       --group_label1 Liver \\\n")
    cat("       --group_label2 Telencephalon\n\n")

    cat("  2) Same, but use the slower xrefs fallback for IDs without symbols:\n")
    cat("     Rscript EV.volcano_with_ensembl_symbols.R \\\n")
    cat("       --de_results host_counts_Liver_vs_Telencephalon.csv \\\n")
    cat("       --label_top 30 \\\n")
    cat("       --xref_fallback\n\n")

    cat("  3) Input already has gene symbols/species names; skip conversion:\n")
    cat("     Rscript EV.volcano_with_ensembl_symbols.R \\\n")
    cat("       --de_results microbiome_DE_results.csv \\\n")
    cat("       --convert no \\\n")
    cat("       --labels\n\n")

    cat("Output table columns added:\n")
    cat("  Feature_ID      = symbol used as main label when available; otherwise original ID/name\n")
    cat("  Original_ID     = original value from the selected ID column\n")
    cat("  Gene_symbol     = converted gene symbol, when found\n")
    cat("  Mapping_status  = converted, unmapped_ensembl, not_ensembl_format, or conversion_disabled\n")
    cat("  Label_in_plot   = TRUE/FALSE, whether the feature was eligible as a plot label\n\n")
}

# -----------------------------
# Argument parser
# -----------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0 || any(args %in% c("-h", "--help"))) {
    print_help()
    if (length(args) == 0) quit(status = 1) else quit(status = 0)
}

# Defaults
de_results <- NULL
output_file <- NULL
output_table <- NULL

id_column <- "1"
convert_mode <- "auto"
sep_option <- "auto"

show_labels <- FALSE
label_top <- NULL
label_unmapped_ensembl <- FALSE

padj_cutoff <- 0.05
logfc_cutoff <- 2

group_label1 <- "Group1"
group_label2 <- "Group2"
subtitle_text <- "Differential abundance/expression"

plot_width <- 10
plot_height <- 8
plot_dpi <- 300

cache_file <- "ensembl_symbol_cache.csv"
server_url <- "https://rest.ensembl.org"
batch_size <- 1000
http_timeout <- 60
http_retries <- 4
sleep_seconds <- 0.15
strip_version <- TRUE
xref_fallback <- FALSE
install_missing <- TRUE

get_next <- function(args, i, opt) {
    if (i + 1 > length(args)) {
        stop(paste0("Error: Missing value after ", opt), call. = FALSE)
    }
    args[i + 1]
}

i <- 1
while (i <= length(args)) {
    arg <- args[i]

    if (arg == "--de_results") {
        de_results <- get_next(args, i, arg)
        i <- i + 2
        next
    }

    if (arg == "--output") {
        output_file <- get_next(args, i, arg)
        i <- i + 2
        next
    }

    if (arg == "--output_table") {
        output_table <- get_next(args, i, arg)
        i <- i + 2
        next
    }

    if (arg == "--id_column") {
        id_column <- get_next(args, i, arg)
        i <- i + 2
        next
    }

    if (arg == "--convert") {
        convert_mode <- get_next(args, i, arg)
        i <- i + 2
        next
    }

    if (arg == "--sep") {
        sep_option <- get_next(args, i, arg)
        i <- i + 2
        next
    }

    if (arg == "--labels") {
        show_labels <- TRUE
        i <- i + 1
        next
    }

    if (arg == "--label_top") {
        label_top <- suppressWarnings(as.numeric(get_next(args, i, arg)))
        i <- i + 2
        next
    }

    if (arg == "--label_unmapped_ensembl") {
        label_unmapped_ensembl <- TRUE
        i <- i + 1
        next
    }

    if (arg == "--padj") {
        padj_cutoff <- suppressWarnings(as.numeric(get_next(args, i, arg)))
        i <- i + 2
        next
    }

    if (arg == "--logfc") {
        logfc_cutoff <- suppressWarnings(as.numeric(get_next(args, i, arg)))
        i <- i + 2
        next
    }

    if (arg == "--group_label1") {
        group_label1 <- get_next(args, i, arg)
        i <- i + 2
        next
    }

    if (arg == "--group_label2") {
        group_label2 <- get_next(args, i, arg)
        i <- i + 2
        next
    }

    if (arg == "--subtitle") {
        subtitle_text <- get_next(args, i, arg)
        i <- i + 2
        next
    }

    if (arg == "--width") {
        plot_width <- suppressWarnings(as.numeric(get_next(args, i, arg)))
        i <- i + 2
        next
    }

    if (arg == "--height") {
        plot_height <- suppressWarnings(as.numeric(get_next(args, i, arg)))
        i <- i + 2
        next
    }

    if (arg == "--dpi") {
        plot_dpi <- suppressWarnings(as.numeric(get_next(args, i, arg)))
        i <- i + 2
        next
    }

    if (arg == "--cache") {
        cache_file <- get_next(args, i, arg)
        i <- i + 2
        next
    }

    if (arg == "--server") {
        server_url <- get_next(args, i, arg)
        i <- i + 2
        next
    }

    if (arg == "--batch_size") {
        batch_size <- suppressWarnings(as.integer(get_next(args, i, arg)))
        i <- i + 2
        next
    }

    if (arg == "--timeout") {
        http_timeout <- suppressWarnings(as.numeric(get_next(args, i, arg)))
        i <- i + 2
        next
    }

    if (arg == "--retries") {
        http_retries <- suppressWarnings(as.integer(get_next(args, i, arg)))
        i <- i + 2
        next
    }

    if (arg == "--sleep") {
        sleep_seconds <- suppressWarnings(as.numeric(get_next(args, i, arg)))
        i <- i + 2
        next
    }

    if (arg == "--keep_version") {
        strip_version <- FALSE
        i <- i + 1
        next
    }

    if (arg == "--xref_fallback") {
        xref_fallback <- TRUE
        i <- i + 1
        next
    }

    if (arg == "--no_install") {
        install_missing <- FALSE
        i <- i + 1
        next
    }

    stop(
        paste0(
            "Error: Unknown option: ", arg, "\n\n",
            "Use:\n",
            "  Rscript EV.volcano_with_ensembl_symbols.R --help\n\n",
            "to see all available options."
        ),
        call. = FALSE
    )
}

# -----------------------------
# Validate arguments
# -----------------------------
if (is.null(de_results)) {
    cat("\nError: Missing required argument --de_results\n\n")
    print_help()
    quit(status = 1)
}

if (!file.exists(de_results)) {
    stop(paste0("Error: Input file not found: ", de_results), call. = FALSE)
}

if (!(convert_mode %in% c("auto", "yes", "no"))) {
    stop("Error: --convert must be one of: auto, yes, no.", call. = FALSE)
}

if (!(sep_option %in% c("auto", "comma", "tab", "semicolon"))) {
    stop("Error: --sep must be one of: auto, comma, tab, semicolon.", call. = FALSE)
}

if (!is.null(label_top) && (is.na(label_top) || label_top <= 0)) {
    stop("Error: --label_top must be a positive number.", call. = FALSE)
}

if (is.na(padj_cutoff) || padj_cutoff <= 0 || padj_cutoff >= 1) {
    stop("Error: --padj must be a number between 0 and 1.", call. = FALSE)
}

if (is.na(logfc_cutoff) || logfc_cutoff < 0) {
    stop("Error: --logfc must be a positive number or zero.", call. = FALSE)
}

if (is.na(plot_width) || plot_width <= 0) {
    stop("Error: --width must be a positive number.", call. = FALSE)
}

if (is.na(plot_height) || plot_height <= 0) {
    stop("Error: --height must be a positive number.", call. = FALSE)
}

if (is.na(plot_dpi) || plot_dpi <= 0) {
    stop("Error: --dpi must be a positive number.", call. = FALSE)
}

if (is.na(batch_size) || batch_size < 1 || batch_size > 1000) {
    stop("Error: --batch_size must be between 1 and 1000.", call. = FALSE)
}

if (is.na(http_timeout) || http_timeout <= 0) {
    stop("Error: --timeout must be a positive number.", call. = FALSE)
}

if (is.na(http_retries) || http_retries < 0) {
    stop("Error: --retries must be zero or a positive integer.", call. = FALSE)
}

if (is.na(sleep_seconds) || sleep_seconds < 0) {
    stop("Error: --sleep must be zero or a positive number.", call. = FALSE)
}

# -----------------------------
# Check/install packages
# -----------------------------
install_if_needed <- function(packages, install_missing = TRUE) {
    for (pkg in packages) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            if (!install_missing) {
                stop(
                    paste0(
                        "Missing R package: ", pkg, "\n",
                        "Install it or rerun without --no_install."
                    ),
                    call. = FALSE
                )
            }

            message("Installing missing package: ", pkg)

            if (pkg == "EnhancedVolcano") {
                if (!requireNamespace("BiocManager", quietly = TRUE)) {
                    install.packages("BiocManager", repos = "https://cloud.r-project.org/")
                }
                BiocManager::install(pkg, ask = FALSE, update = FALSE)
            } else {
                install.packages(pkg, dependencies = TRUE, repos = "https://cloud.r-project.org/")
            }
        }

        suppressPackageStartupMessages(
            library(pkg, character.only = TRUE)
        )
    }
}

necessary_packages <- c("ggplot2", "ggrepel", "EnhancedVolcano", "httr", "jsonlite")
install_if_needed(necessary_packages, install_missing = install_missing)

# -----------------------------
# Utility functions
# -----------------------------
detect_separator <- function(file) {
    first_lines <- readLines(file, n = 5, warn = FALSE)
    first_lines <- first_lines[nchar(first_lines) > 0]

    if (length(first_lines) == 0) {
        stop("Error: Input file appears to be empty.", call. = FALSE)
    }

    first_line <- first_lines[1]

    n_tabs <- lengths(regmatches(first_line, gregexpr("\t", first_line, fixed = TRUE)))
    n_commas <- lengths(regmatches(first_line, gregexpr(",", first_line, fixed = TRUE)))
    n_semicolons <- lengths(regmatches(first_line, gregexpr(";", first_line, fixed = TRUE)))

    if (n_tabs >= n_commas && n_tabs >= n_semicolons && n_tabs > 0) return("\t")
    if (n_semicolons >= n_commas && n_semicolons > 0) return(";")
    return(",")
}

resolve_column_index <- function(col_spec, col_names, ncols) {
    if (grepl("^[0-9]+$", col_spec)) {
        idx <- as.integer(col_spec)
        if (idx == 0) idx <- 1
        if (idx < 1 || idx > ncols) {
            stop(
                paste0("Error: --id_column ", col_spec, " is outside the column range 1..", ncols),
                call. = FALSE
            )
        }
        return(idx)
    }

    hit <- match(col_spec, col_names)
    if (is.na(hit)) {
        stop(
            paste0(
                "Error: Column '", col_spec, "' not found. Available columns: ",
                paste(col_names, collapse = ", ")
            ),
            call. = FALSE
        )
    }
    hit
}

normalize_id <- function(x, strip_version = TRUE) {
    x <- trimws(as.character(x))
    if (strip_version) {
        x <- sub("\\.[0-9]+$", "", x)
    }
    x
}

ensembl_gene_pattern <- "^ENS[A-Za-z0-9]*G[0-9]+(\\.[0-9]+)?$"

looks_like_ensembl_gene <- function(x) {
    grepl(ensembl_gene_pattern, x)
}

useful_symbol <- function(raw_id, candidate) {
    if (is.null(candidate) || length(candidate) == 0) return(NA_character_)
    candidate <- trimws(as.character(candidate[1]))
    if (is.na(candidate) || candidate == "") return(NA_character_)
    if (candidate == raw_id) return(NA_character_)
    if (looks_like_ensembl_gene(candidate)) return(NA_character_)
    candidate
}

make_default_output <- function(input_file) {
    base <- tools::file_path_sans_ext(input_file)
    paste0(base, "_volcano.pdf")
}

make_default_table <- function(input_file) {
    base <- tools::file_path_sans_ext(input_file)
    paste0(base, "_gene_symbols_table.csv")
}

# -----------------------------
# Cache functions
# -----------------------------
read_symbol_cache <- function(cache_file) {
    if (is.null(cache_file) || cache_file == "" || !file.exists(cache_file)) {
        out <- character(0)
        return(out)
    }

    cache_df <- tryCatch(
        read.csv(cache_file, stringsAsFactors = FALSE, check.names = FALSE),
        error = function(e) NULL
    )

    if (is.null(cache_df) || !("ensembl_id" %in% colnames(cache_df)) || !("gene_symbol" %in% colnames(cache_df))) {
        warning("Cache file exists but does not have columns ensembl_id and gene_symbol. Ignoring cache.")
        out <- character(0)
        return(out)
    }

    ids <- trimws(as.character(cache_df$ensembl_id))
    symbols <- trimws(as.character(cache_df$gene_symbol))
    symbols[symbols == ""] <- NA_character_

    keep <- !is.na(ids) & ids != ""
    ids <- ids[keep]
    symbols <- symbols[keep]

    out <- symbols
    names(out) <- ids
    out[!duplicated(names(out))]
}

write_symbol_cache <- function(cache_file, cache) {
    if (is.null(cache_file) || cache_file == "") return(invisible(NULL))
    if (length(cache) == 0) return(invisible(NULL))

    cache <- cache[order(names(cache))]
    cache_df <- data.frame(
        ensembl_id = names(cache),
        gene_symbol = ifelse(is.na(cache), "", as.character(cache)),
        stringsAsFactors = FALSE,
        check.names = FALSE
    )

    tmp <- paste0(cache_file, ".tmp")
    write.csv(cache_df, tmp, row.names = FALSE, quote = TRUE)
    if (file.exists(cache_file)) unlink(cache_file)
    ok <- file.rename(tmp, cache_file)
    if (!ok) stop(paste0("Could not write cache file: ", cache_file), call. = FALSE)
    invisible(NULL)
}

# -----------------------------
# Ensembl REST functions
# -----------------------------
http_json <- function(method, url, body = NULL, timeout = 60, retries = 4, verbose = FALSE) {
    last_error <- NULL
    last_status <- NA_integer_
    max_attempts <- retries + 1

    for (attempt in seq_len(max_attempts)) {
        resp <- tryCatch({
            if (toupper(method) == "POST") {
                httr::POST(
                    url,
                    httr::user_agent("EV-volcano-ensembl-symbols-R/1.0"),
                    httr::accept_json(),
                    httr::content_type_json(),
                    httr::timeout(timeout),
                    body = body,
                    encode = "json"
                )
            } else {
                httr::GET(
                    url,
                    httr::user_agent("EV-volcano-ensembl-symbols-R/1.0"),
                    httr::accept_json(),
                    httr::timeout(timeout)
                )
            }
        }, error = function(e) e)

        if (inherits(resp, "response")) {
            last_status <- httr::status_code(resp)

            if (!httr::http_error(resp)) {
                text <- httr::content(resp, as = "text", encoding = "UTF-8")
                if (is.null(text) || text == "") return(list())
                return(jsonlite::fromJSON(text, simplifyVector = FALSE))
            }

            last_error <- paste0("HTTP ", last_status, ": ", httr::content(resp, as = "text", encoding = "UTF-8"))

            retry_after <- httr::headers(resp)[["retry-after"]]
            if (!is.null(retry_after) && !is.na(suppressWarnings(as.numeric(retry_after)))) {
                wait <- as.numeric(retry_after)
            } else {
                wait <- min(2^attempt, 30)
            }
        } else {
            last_error <- conditionMessage(resp)
            wait <- min(2^attempt, 30)
        }

        if (attempt < max_attempts) {
            if (verbose) message("[WARN] Request failed for ", url, "; retrying in ", wait, " seconds.")
            Sys.sleep(wait)
        }
    }

    stop(
        paste0("Request failed after retries: ", url, "\n", last_error),
        call. = FALSE
    )
}

split_into_chunks <- function(x, size) {
    split(x, ceiling(seq_along(x) / size))
}

fetch_lookup_symbols <- function(ids, server, batch_size, timeout, retries, sleep_seconds) {
    if (length(ids) == 0) {
        out <- character(0)
        return(out)
    }

    endpoint <- paste0(sub("/$", "", server), "/lookup/id")
    results <- rep(NA_character_, length(ids))
    names(results) <- ids

    id_chunks <- split_into_chunks(ids, batch_size)

    for (batch_i in seq_along(id_chunks)) {
        batch <- id_chunks[[batch_i]]
        message("[INFO] Ensembl /lookup/id batch ", batch_i, "/", length(id_chunks), ": ", length(batch), " IDs")

        decoded <- http_json(
            method = "POST",
            url = endpoint,
            body = list(ids = as.list(batch)),
            timeout = timeout,
            retries = retries
        )

        for (gene_id in batch) {
            item <- decoded[[gene_id]]
            if (!is.null(item) && is.list(item) && !is.null(item$display_name)) {
                results[gene_id] <- useful_symbol(gene_id, item$display_name)
            } else {
                results[gene_id] <- NA_character_
            }
        }

        Sys.sleep(sleep_seconds)
    }

    results
}

score_xref_record <- function(rec) {
    get_chr <- function(x) {
        if (is.null(x) || length(x) == 0 || is.na(x[1])) return("")
        as.character(x[1])
    }

    dbname <- tolower(get_chr(rec$dbname))
    primary_id <- get_chr(rec$primary_id)
    display_id <- get_chr(rec$display_id)
    description <- tolower(get_chr(rec$description))

    score <- 0
    if (grepl("symbol", dbname, fixed = TRUE)) score <- score + 100
    if (dbname %in% c("hgnc", "hgnc symbol", "mgi", "mgi symbol", "rgd", "rgd symbol", "vgnc", "vgnc symbol")) score <- score + 80
    if (grepl("name", dbname, fixed = TRUE) || grepl("gene", dbname, fixed = TRUE)) score <- score + 20
    if (display_id != "" && !startsWith(display_id, "ENS")) score <- score + 10
    if (primary_id != "" && !startsWith(primary_id, "ENS")) score <- score + 5
    if (grepl("predicted", description, fixed = TRUE)) score <- score - 5
    score
}

fetch_xref_symbol_for_id <- function(gene_id, server, timeout, retries) {
    endpoint <- paste0(
        sub("/$", "", server),
        "/xrefs/id/",
        utils::URLencode(gene_id, reserved = TRUE),
        "?object_type=gene"
    )

    decoded <- http_json(
        method = "GET",
        url = endpoint,
        timeout = timeout,
        retries = retries
    )

    if (!is.list(decoded) || length(decoded) == 0) return(NA_character_)

    candidates <- data.frame(score = numeric(0), symbol = character(0), stringsAsFactors = FALSE)

    for (rec in decoded) {
        if (!is.list(rec)) next
        for (field in c("display_id", "primary_id")) {
            candidate <- useful_symbol(gene_id, rec[[field]])
            if (!is.na(candidate)) {
                candidates <- rbind(
                    candidates,
                    data.frame(score = score_xref_record(rec), symbol = candidate, stringsAsFactors = FALSE)
                )
            }
        }
    }

    if (nrow(candidates) == 0) return(NA_character_)
    candidates <- candidates[order(candidates$score, decreasing = TRUE), , drop = FALSE]
    candidates$symbol[1]
}

fetch_xref_fallback_symbols <- function(ids, server, timeout, retries, sleep_seconds) {
    if (length(ids) == 0) {
        out <- character(0)
        return(out)
    }

    results <- rep(NA_character_, length(ids))
    names(results) <- ids

    for (idx in seq_along(ids)) {
        gene_id <- ids[idx]
        if (idx == 1 || idx %% 100 == 0 || idx == length(ids)) {
            message("[INFO] Ensembl /xrefs/id fallback ", idx, "/", length(ids))
        }

        results[gene_id] <- tryCatch(
            fetch_xref_symbol_for_id(gene_id, server, timeout, retries),
            error = function(e) {
                warning("xref fallback failed for ", gene_id, ": ", conditionMessage(e))
                NA_character_
            }
        )

        Sys.sleep(sleep_seconds)
    }

    results
}

# -----------------------------
# Detect separator and load table
# -----------------------------
if (sep_option == "auto") {
    sep <- detect_separator(de_results)
} else if (sep_option == "comma") {
    sep <- ","
} else if (sep_option == "tab") {
    sep <- "\t"
} else if (sep_option == "semicolon") {
    sep <- ";"
}

message("Input file: ", de_results)
message("Detected/selected separator: ", ifelse(sep == "\t", "tab", sep))

res_raw <- read.table(
    de_results,
    header = TRUE,
    sep = sep,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "\"",
    comment.char = "",
    fill = TRUE
)

if (ncol(res_raw) < 3) {
    stop("Error: Input file has too few columns. Check the separator using --sep.", call. = FALSE)
}

id_col_idx <- resolve_column_index(id_column, colnames(res_raw), ncol(res_raw))
original_id_colname <- colnames(res_raw)[id_col_idx]
original_ids <- trimws(as.character(res_raw[[id_col_idx]]))
normalized_ids <- normalize_id(original_ids, strip_version = strip_version)
ensembl_like <- looks_like_ensembl_gene(normalized_ids)

message("Selected ID column: ", original_id_colname, " [column ", id_col_idx, "]")
message("Input rows: ", nrow(res_raw))
message("Ensembl-like gene IDs detected: ", sum(ensembl_like, na.rm = TRUE))

# -----------------------------
# Convert Ensembl IDs to symbols when needed
# -----------------------------
gene_symbol <- rep(NA_character_, nrow(res_raw))
mapping_status <- rep("not_ensembl_format", nrow(res_raw))

if (convert_mode == "no") {
    mapping_status[] <- "conversion_disabled"
    message("Conversion mode: no. No online lookup will be performed.")
} else if (!any(ensembl_like, na.rm = TRUE)) {
    message("Conversion mode: ", convert_mode, ". No Ensembl-like IDs found, so no online lookup is needed.")
} else {
    ids_to_convert <- sort(unique(normalized_ids[ensembl_like & !is.na(normalized_ids) & normalized_ids != ""]))

    message("Conversion mode: ", convert_mode)
    message("Unique Ensembl-like IDs to convert: ", length(ids_to_convert))

    cache <- read_symbol_cache(cache_file)
    found_in_cache <- ids_to_convert[ids_to_convert %in% names(cache)]
    missing_from_cache <- setdiff(ids_to_convert, names(cache))

    message("Found in cache: ", length(found_in_cache))
    message("To query online: ", length(missing_from_cache))

    if (length(missing_from_cache) > 0) {
        lookup_symbols <- fetch_lookup_symbols(
            ids = missing_from_cache,
            server = server_url,
            batch_size = batch_size,
            timeout = http_timeout,
            retries = http_retries,
            sleep_seconds = sleep_seconds
        )

        cache[names(lookup_symbols)] <- lookup_symbols
        write_symbol_cache(cache_file, cache)
    }

    if (xref_fallback) {
        still_missing <- ids_to_convert[is.na(cache[ids_to_convert]) | cache[ids_to_convert] == ""]
        message("IDs needing /xrefs/id fallback: ", length(still_missing))

        if (length(still_missing) > 0) {
            xref_symbols <- fetch_xref_fallback_symbols(
                ids = still_missing,
                server = server_url,
                timeout = http_timeout,
                retries = http_retries,
                sleep_seconds = sleep_seconds
            )

            has_xref <- !is.na(xref_symbols) & xref_symbols != ""
            cache[names(xref_symbols)[has_xref]] <- xref_symbols[has_xref]
            cache[names(xref_symbols)[!has_xref]] <- NA_character_
            write_symbol_cache(cache_file, cache)
        }
    }

    gene_symbol[ensembl_like] <- cache[normalized_ids[ensembl_like]]

    converted_rows <- ensembl_like & !is.na(gene_symbol) & gene_symbol != ""
    unmapped_rows <- ensembl_like & !converted_rows

    mapping_status[converted_rows] <- "converted"
    mapping_status[unmapped_rows] <- "unmapped_ensembl"

    message("Converted rows: ", sum(converted_rows, na.rm = TRUE))
    message("Unmapped Ensembl rows: ", sum(unmapped_rows, na.rm = TRUE))
    message("Cache file: ", cache_file)
}

feature_id <- ifelse(!is.na(gene_symbol) & gene_symbol != "", gene_symbol, original_ids)
feature_id <- trimws(as.character(feature_id))

plot_label <- feature_id
if (!label_unmapped_ensembl) {
    plot_label[mapping_status == "unmapped_ensembl"] <- ""
}
plot_label[is.na(plot_label)] <- ""

# -----------------------------
# Build converted/proof table
# -----------------------------
remaining_cols <- res_raw[, -id_col_idx, drop = FALSE]

res2 <- data.frame(
    Feature_ID = feature_id,
    Original_ID = original_ids,
    Gene_symbol = gene_symbol,
    Mapping_status = mapping_status,
    Label_in_plot = plot_label != "",
    remaining_cols,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

# Check required columns
required_cols <- c("log2FoldChange", "padj")
missing_cols <- setdiff(required_cols, colnames(res2))

if (length(missing_cols) > 0) {
    stop(
        paste("Error: Missing required column(s):", paste(missing_cols, collapse = ", ")),
        call. = FALSE
    )
}

# Convert important columns to numeric
res2$log2FoldChange <- suppressWarnings(as.numeric(res2$log2FoldChange))
res2$padj <- suppressWarnings(as.numeric(res2$padj))

# Output table path
if (is.null(output_table)) {
    output_table <- make_default_table(de_results)
}

write.csv(res2, output_table, row.names = FALSE, quote = TRUE)
message("Converted/proof table saved as: ", output_table)

# -----------------------------
# Clean rows for volcano only
# -----------------------------
n_before <- nrow(res2)

valid_rows <- !is.na(res2$log2FoldChange) &
    !is.na(res2$padj) &
    !is.na(res2$Feature_ID) &
    res2$Feature_ID != ""

res2 <- res2[valid_rows, , drop = FALSE]
plot_label <- plot_label[valid_rows]

n_after <- nrow(res2)

if (n_after == 0) {
    stop("Error: No valid rows left after filtering NA values in log2FoldChange, padj, or Feature_ID.", call. = FALSE)
}

if (n_after < n_before) {
    message("Removed ", n_before - n_after, " row(s) with missing/invalid values for plotting.")
}

# Avoid problems with padj = 0
res2$padj_plot <- ifelse(res2$padj <= 0, .Machine$double.xmin, res2$padj)

# -----------------------------
# Output plot file
# -----------------------------
if (is.null(output_file)) {
    output_file <- make_default_output(de_results)
}

if (!grepl("\\.(pdf|png)$", output_file, ignore.case = TRUE)) {
    output_file <- paste0(output_file, ".pdf")
}

message("Output volcano plot: ", output_file)

# -----------------------------
# Dynamic axis limits
# -----------------------------
x_min <- min(res2$log2FoldChange, na.rm = TRUE)
x_max <- max(res2$log2FoldChange, na.rm = TRUE)
x_abs <- max(abs(c(x_min, x_max)), na.rm = TRUE)
x_padding <- x_abs * 0.08
if (x_padding == 0 || is.na(x_padding)) x_padding <- 1
x_range <- c(x_min - x_padding, x_max + x_padding)

y_values <- -log10(res2$padj_plot)
y_max <- max(y_values, na.rm = TRUE)
if (is.infinite(y_max) || is.na(y_max)) y_max <- 1
y_range <- c(0, y_max * 1.12)

# -----------------------------
# Significant features
# -----------------------------
significant_features <- sum(
    res2$padj < padj_cutoff & abs(res2$log2FoldChange) >= logfc_cutoff,
    na.rm = TRUE
)

message("Total variables plotted: ", nrow(res2))
message("Significant variables: ", significant_features)

# -----------------------------
# Select labels
# -----------------------------
res2$Plot_Label <- plot_label
res_label_candidates <- res2[!is.na(res2$Plot_Label) & res2$Plot_Label != "", , drop = FALSE]

if (!is.null(label_top)) {
    ord <- order(res_label_candidates$padj, -abs(res_label_candidates$log2FoldChange), na.last = NA)
    labels_to_show <- unique(res_label_candidates$Plot_Label[ord])
    labels_to_show <- head(labels_to_show, label_top)
    message("Label mode: top ", label_top, " by lowest padj, excluding blank/suppressed labels.")
} else if (show_labels) {
    keep <- res_label_candidates$padj < padj_cutoff &
        abs(res_label_candidates$log2FoldChange) >= logfc_cutoff
    labels_to_show <- unique(res_label_candidates$Plot_Label[keep])
    message("Label mode: all significant features, excluding blank/suppressed labels.")
} else {
    labels_to_show <- character(0)
    message("Label mode: no labels.")
}

labels_to_show <- labels_to_show[!is.na(labels_to_show) & labels_to_show != ""]
message("Labels selected for plot: ", length(labels_to_show))

# -----------------------------
# Make volcano plot
# -----------------------------
volcano_plot <- EnhancedVolcano::EnhancedVolcano(
    res2,

    lab = res2$Plot_Label,
    selectLab = labels_to_show,

    x = "log2FoldChange",
    y = "padj_plot",

    title = bquote(.(group_label1)~italic(versus)~.(group_label2)),
    subtitle = subtitle_text,

    xlab = bquote(~Log[2]~"Fold Change"),
    ylab = bquote(-1*Log[10]~"FDR"),

    pCutoff = padj_cutoff,
    FCcutoff = logfc_cutoff,

    pointSize = 2.0,
    labSize = 3.0,
    labFace = "bold",
    boxedLabels = TRUE,

    colAlpha = 1,

    cutoffLineType = "blank",

    hline = c(padj_cutoff),
    hlineCol = c("grey75"),
    hlineType = "longdash",
    hlineWidth = 0.4,

    vline = c(-logfc_cutoff, logfc_cutoff),
    vlineCol = c("grey75"),
    vlineType = "longdash",
    vlineWidth = 0.4,

    drawConnectors = TRUE,
    colConnectors = "grey50",
    maxoverlapsConnectors = Inf,
    widthConnectors = 0.5,

    border = "full",
    borderColour = "black",
    borderWidth = 0.5,

    gridlines.major = TRUE,
    gridlines.minor = FALSE,

    legendPosition = "top",
    legendLabSize = 12,
    legendIconSize = 3.0,
    legendLabels = c(
        "NS",
        expression(Log[2]~FC),
        "FDR",
        expression(FDR~and~Log[2]~FC)
    )
) +
    ggplot2::coord_cartesian(
        xlim = x_range,
        ylim = y_range,
        clip = "off"
    ) +
    ggplot2::scale_x_continuous(
        breaks = seq(floor(x_range[1]), ceiling(x_range[2]), length.out = 5)
    ) +
    ggplot2::scale_y_continuous(
        breaks = seq(0, ceiling(y_range[2]), length.out = 5)
    ) +
    ggplot2::annotate(
        "text",
        x = logfc_cutoff + 1,
        y = max(y_range) * 0.95,
        label = paste0("+", logfc_cutoff),
        color = "grey75",
        hjust = 0,
        fontface = "italic",
        size = 4
    ) +
    ggplot2::annotate(
        "text",
        x = -logfc_cutoff - 1,
        y = max(y_range) * 0.95,
        label = paste0("-", logfc_cutoff),
        color = "grey75",
        hjust = 1,
        fontface = "italic",
        size = 4
    ) +
    ggplot2::annotate(
        "text",
        x = min(x_range),
        y = -log10(padj_cutoff) + 0.5,
        label = as.character(padj_cutoff),
        color = "grey75",
        hjust = 0,
        fontface = "italic",
        size = 4
    ) +
    ggplot2::labs(
        caption = paste(
            "Total =",
            nrow(res2),
            "variables | Significant =",
            significant_features
        )
    ) +
    ggplot2::theme(
        plot.margin = ggplot2::margin(10, 35, 10, 10)
    )

# ggsave automatically chooses the device from the file extension.
ggplot2::ggsave(
    filename = output_file,
    plot = volcano_plot,
    width = plot_width,
    height = plot_height,
    dpi = plot_dpi,
    units = "in"
)

message("Done.")
message("Volcano plot saved as: ", output_file)
message("Converted/proof table saved as: ", output_table)

