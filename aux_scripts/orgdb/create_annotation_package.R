#!/usr/bin/env Rscript

# ============================================================
# create_annotation_package.R
#
# Cria pacote OrgDb usando AnnotationForge::makeOrgPackageFromNCBI
# com suporte real a modo offline.
#
# Exemplo:
#
# Rscript $MTDIR/create_annotation_package.R \
#   --taxid $customized \
#   --offline /path/to/makeOrgPackageFromNCBI/ \
#   --copy \
#   --clean
# ============================================================

options(repos = c(CRAN = "https://cran.r-project.org"))
options(timeout = 600)

# ============================================================
# Pacotes
# ============================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

if (!requireNamespace("AnnotationForge", quietly = TRUE)) {
    BiocManager::install("AnnotationForge", force = TRUE, update = FALSE)
}

if (!requireNamespace("optparse", quietly = TRUE)) {
    install.packages("optparse")
}

if (!requireNamespace("readr", quietly = TRUE)) {
    install.packages("readr")
}

if (!requireNamespace("httr", quietly = TRUE)) {
    install.packages("httr")
}

suppressPackageStartupMessages({
    library(AnnotationForge)
    library(optparse)
    library(readr)
    library(httr)
})

# ============================================================
# Argumentos
# ============================================================

option_list <- list(
    make_option(
        c("-t", "--taxid"),
        type = "character",
        default = NULL,
        help = "NCBI Taxonomy ID of the species [required]"
    ),

    make_option(
        c("-d", "--dest_dir"),
        type = "character",
        default = getwd(),
        help = "Output directory for annotation package [default: current directory]"
    ),

    make_option(
        c("-o", "--offline"),
        type = "character",
        default = NULL,
        help = "Path to directory with pre-downloaded NCBI files"
    ),

    make_option(
        c("-c", "--copy"),
        action = "store_true",
        default = FALSE,
        help = "Copy offline NCBI files to ./NCBI before creating package"
    ),

    make_option(
        c("--clean"),
        action = "store_true",
        default = FALSE,
        help = "Remove existing ./NCBI before starting"
    ),

    make_option(
        c("--use_ensembl"),
        action = "store_true",
        default = FALSE,
        help = "Allow online Ensembl lookup. Default is FALSE for offline mode."
    ),

    make_option(
        c("--check_remote"),
        action = "store_true",
        default = FALSE,
        help = "Compare local NCBI file sizes with remote NCBI files. Requires internet."
    )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$taxid)) {
    stop("The taxid must be provided using -t or --taxid", call. = FALSE)
}

message("Output directory for annotation package: ", opt$dest_dir)

# ============================================================
# Arquivos NCBI necessários
# ============================================================

ncbi_files <- c(
    "gene2pubmed.gz",
    "gene2accession.gz",
    "gene2refseq.gz",
    "gene_info.gz",
    "gene2go.gz"
)

local_dir <- file.path(getwd(), "NCBI")

if (opt$clean && dir.exists(local_dir)) {
    message("Removing existing local NCBI directory: ", local_dir)
    unlink(local_dir, recursive = TRUE, force = TRUE)
}

dir.create(local_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# HostSpecies.csv
# ============================================================

csv_file_path <- file.path(getwd(), "HostSpecies.csv")

if (!file.exists(csv_file_path)) {
    stop(
        "The file HostSpecies.csv was not found in current directory: ",
        getwd(),
        call. = FALSE
    )
}

species_data <- read_csv(csv_file_path, show_col_types = FALSE)

if (!all(c("Taxon_ID", "Scientific_name") %in% colnames(species_data))) {
    stop(
        "The CSV file must contain 'Taxon_ID' and 'Scientific_name' columns.",
        call. = FALSE
    )
}

get_species_info_from_csv <- function(taxid, species_data) {
    matched <- species_data[species_data$Taxon_ID == as.numeric(taxid), ]

    if (nrow(matched) == 0) {
        stop("No species found for TaxID: ", taxid, call. = FALSE)
    }

    sp_name <- matched$Scientific_name[1]
    sp_split <- strsplit(sp_name, " ")[[1]]

    if (length(sp_split) < 2) {
        stop(
            "Scientific_name must be at least 'Genus species'. Found: ",
            sp_name,
            call. = FALSE
        )
    }

    genus <- sp_split[1]
    species <- sp_split[2]

    message("Taxid: ", taxid)
    message("Genus: ", genus)
    message("Species: ", species)

    list(genus = genus, species = species)
}

species_info <- get_species_info_from_csv(opt$taxid, species_data)

# ============================================================
# Copiar arquivos offline corretamente
# ============================================================

copy_ncbi_files_flat <- function(offline_dir, local_dir, ncbi_files) {
    if (!dir.exists(offline_dir)) {
        stop("Offline directory does not exist: ", offline_dir, call. = FALSE)
    }

    all_files <- list.files(
        offline_dir,
        recursive = TRUE,
        full.names = TRUE,
        include.dirs = FALSE
    )

    found <- setNames(rep(NA_character_, length(ncbi_files)), ncbi_files)

    for (f in ncbi_files) {
        hits <- all_files[basename(all_files) == f]

        if (length(hits) > 0) {
            found[f] <- hits[1]
        }
    }

    missing <- names(found)[is.na(found)]

    if (length(missing) > 0) {
        stop(
            "Missing required NCBI files in offline directory:\n  - ",
            paste(missing, collapse = "\n  - "),
            "\n\nSearched inside: ",
            offline_dir,
            call. = FALSE
        )
    }

    message("Copying required NCBI files to: ", local_dir)

    for (f in ncbi_files) {
        from <- found[f]
        to <- file.path(local_dir, f)

        ok <- file.copy(
            from = from,
            to = to,
            overwrite = TRUE
        )

        if (!ok) {
            stop("Failed to copy: ", from, " -> ", to, call. = FALSE)
        }

        message("Copied: ", f)
    }
}

if (!is.null(opt$offline)) {
    message("Offline mode enabled.")

    if (opt$copy) {
        copy_ncbi_files_flat(opt$offline, local_dir, ncbi_files)
    } else {
        message("You used --offline but not --copy.")
        message("Expecting files already present in: ", local_dir)
    }
}

# ============================================================
# Verificar arquivos locais
# ============================================================

check_local_ncbi_files <- function(local_dir, ncbi_files) {
    missing <- c()

    for (f in ncbi_files) {
        local_file <- file.path(local_dir, f)

        if (!file.exists(local_file)) {
            missing <- c(missing, f)
        } else {
            size_mb <- round(file.info(local_file)$size / 1024^2, 2)
            message("Local NCBI file OK: ", f, " (", size_mb, " MB)")
        }
    }

    if (length(missing) > 0) {
        stop(
            "Required NCBI files are missing from ",
            local_dir,
            ":\n  - ",
            paste(missing, collapse = "\n  - "),
            call. = FALSE
        )
    }
}

check_local_ncbi_files(local_dir, ncbi_files)

# ============================================================
# Checagem remota opcional
# Só use --check_remote se tiver internet/DNS funcionando.
# ============================================================

check_ncbi_file_consistency <- function(local_path, remote_url) {
    if (!file.exists(local_path)) {
        warning("Local file not found: ", local_path)
        return(FALSE)
    }

    local_size <- file.info(local_path)$size

    r <- tryCatch(
        httr::HEAD(remote_url),
        error = function(e) e
    )

    if (inherits(r, "error")) {
        warning("Could not contact remote NCBI URL: ", remote_url)
        return(FALSE)
    }

    remote_size <- as.numeric(httr::headers(r)[["content-length"]])

    if (is.na(remote_size)) {
        warning("Could not retrieve remote file size for: ", remote_url)
        return(FALSE)
    }

    if (local_size != remote_size) {
        warning(
            "File mismatch detected: ", basename(local_path), "\n",
            "  Local size: ", local_size, " bytes\n",
            "  Remote size: ", remote_size, " bytes\n",
            "Suggestion: re-download this file from NCBI."
        )
        return(FALSE)
    }

    message("Remote check OK: ", basename(local_path))
    TRUE
}

if (opt$check_remote) {
    message("Remote NCBI consistency check enabled.")

    ncbi_base_url <- "https://ftp.ncbi.nlm.nih.gov/gene/DATA/"

    for (f in ncbi_files) {
        local_file <- file.path(local_dir, f)
        remote_file <- paste0(ncbi_base_url, f)
        check_ncbi_file_consistency(local_file, remote_file)
    }
} else {
    message("Skipping remote NCBI consistency check.")
}

# ============================================================
# Patch 1: desativar consulta ao Ensembl em modo offline
# ============================================================

disable_ensembl_lookup <- function() {
    ns <- asNamespace("AnnotationForge")

    if (!exists("available.ensembl.datasets", envir = ns, inherits = FALSE)) {
        warning("Could not find AnnotationForge::available.ensembl.datasets to patch.")
        return(invisible(FALSE))
    }

    message("Disabling Ensembl lookup for offline mode.")

    replacement_fun <- function(...) {
        message("Skipping Ensembl dataset lookup.")
        stats::setNames(character(0), character(0))
    }

    was_locked <- bindingIsLocked("available.ensembl.datasets", ns)

    if (was_locked) {
        unlockBinding("available.ensembl.datasets", ns)
    }

    assign("available.ensembl.datasets", replacement_fun, envir = ns)

    if (was_locked) {
        lockBinding("available.ensembl.datasets", ns)
    }

    invisible(TRUE)
}

# ============================================================
# Patch 2: usar arquivos .gz locais ao reconstruir cache SQLite
#
# O AnnotationForge ainda imprime "starting download for",
# mas com este patch ele usa os arquivos locais e não baixa nada.
# ============================================================

enable_offline_ncbi_cache <- function(local_dir) {
    ns <- asNamespace("AnnotationForge")

    if (!exists(".tryDL", envir = ns, inherits = FALSE)) {
        warning("Could not find AnnotationForge::.tryDL to patch.")
        warning("If your AnnotationForge version uses another download function, this may still try internet.")
        return(invisible(FALSE))
    }

    message("Patching AnnotationForge download function for offline NCBI cache mode.")

    replacement_fun <- function(url, tmp, ...) {
        url_base <- basename(url)
        tmp_base <- basename(tmp)

        candidates <- unique(c(
            tmp,
            file.path(local_dir, tmp_base),
            file.path(local_dir, url_base)
        ))

        candidates <- candidates[!is.na(candidates)]

        existing <- candidates[file.exists(candidates)]

        if (length(existing) > 0) {
            src <- existing[1]

            if (!identical(normalizePath(src, mustWork = FALSE),
                           normalizePath(tmp, mustWork = FALSE))) {
                dir.create(dirname(tmp), showWarnings = FALSE, recursive = TRUE)
                ok <- file.copy(src, tmp, overwrite = TRUE)

                if (!ok) {
                    stop(
                        "Found offline file but failed to copy it:\n",
                        "  From: ", src, "\n",
                        "  To: ", tmp,
                        call. = FALSE
                    )
                }
            }

            message("Offline cache file found, skipping download: ", basename(src))
            return(invisible(0))
        }

        stop(
            "Offline mode is enabled, but required file was not found locally.\n",
            "Tried:\n  - ",
            paste(candidates, collapse = "\n  - "),
            "\nOriginal URL would have been:\n",
            url,
            call. = FALSE
        )
    }

    was_locked <- bindingIsLocked(".tryDL", ns)

    if (was_locked) {
        unlockBinding(".tryDL", ns)
    }

    assign(".tryDL", replacement_fun, envir = ns)

    if (was_locked) {
        lockBinding(".tryDL", ns)
    }

    invisible(TRUE)
}

# ============================================================
# Aplicar patches
# ============================================================

if (!opt$use_ensembl) {
    disable_ensembl_lookup()
} else {
    message("Ensembl lookup enabled. This requires DNS/internet access to www.ensembl.org.")
}

if (!is.null(opt$offline)) {
    enable_offline_ncbi_cache(local_dir)
}

# ============================================================
# Criar pacote
# ============================================================

message("Creating annotation package...")

httr::set_config(httr::config(ssl_verifypeer = FALSE))

make_args <- list(
    version = "0.1",
    author = "Generic Author",
    maintainer = "Generic Maintainer <maintainer@example.com>",
    outputDir = opt$dest_dir,
    tax_id = opt$taxid,
    genus = species_info$genus,
    species = species_info$species,
    NCBIFilesDir = local_dir
)

# Importante:
# rebuildCache precisa ser TRUE, senão o SQLite NCBI.sqlite pode não ser populado.
if ("rebuildCache" %in% names(formals(AnnotationForge::makeOrgPackageFromNCBI))) {
    make_args$rebuildCache <- TRUE
}

if ("verbose" %in% names(formals(AnnotationForge::makeOrgPackageFromNCBI))) {
    make_args$verbose <- TRUE
}

do.call(AnnotationForge::makeOrgPackageFromNCBI, make_args)

message("Annotation package created successfully.")
