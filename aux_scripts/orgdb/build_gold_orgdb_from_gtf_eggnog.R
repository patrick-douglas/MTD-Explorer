#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(AnnotationForge)
})

args <- commandArgs(trailingOnly = TRUE)

has_flag <- function(flag) flag %in% args

get_arg <- function(flag, default = NA_character_) {
  hit <- which(args == flag)
  if (length(hit) == 0) return(default)
  if (hit == length(args)) stop('Missing value after ', flag)
  args[hit + 1]
}

split_csv <- function(x, default) {
  if (is.na(x) || !nzchar(x)) return(default)
  trimws(unlist(strsplit(x, ',', fixed = TRUE)))
}

gtf_file <- get_arg('--gtf')
eggnog_file <- get_arg('--eggnog')
lib_dir <- get_arg('--lib', '/home/me/MTD/custom_R_libs')
build_dir <- get_arg('--build-dir', '/home/me/MTD/build/orgdb_gold')
version <- get_arg('--version', '0.1.0')
taxid <- get_arg('--taxid')
genus <- get_arg('--genus')
species <- get_arg('--species')
author <- get_arg('--author', 'MTD custom annotation builder')
maintainer <- get_arg('--maintainer', 'MTD maintainer <maintainer@example.com>')
symbol_mode <- get_arg('--symbol-mode', 'gene_id')
gene_attr_order <- split_csv(get_arg('--gene-attr-order', NA_character_), c('gene_id', 'ID', 'locus_tag', 'Name'))
gene_name_attr_order <- split_csv(get_arg('--gene-name-attr-order', NA_character_), c('gene_name', 'Name', 'gene', 'product'))
description_attr_order <- split_csv(get_arg('--description-attr-order', NA_character_), c('description', 'product', 'Note', 'gene_name', 'Name'))
query_gene_regex <- get_arg('--query-gene-regex', NA_character_)
force <- has_flag('--force')

if (is.na(gtf_file) || !file.exists(gtf_file)) stop('[ERROR] Missing or invalid --gtf: ', gtf_file)
if (is.na(eggnog_file) || !file.exists(eggnog_file)) stop('[ERROR] Missing or invalid --eggnog: ', eggnog_file)
if (is.na(taxid) || !nzchar(taxid)) stop('[ERROR] Missing --taxid')
if (is.na(genus) || !nzchar(genus)) stop('[ERROR] Missing --genus')
if (is.na(species) || !nzchar(species)) stop('[ERROR] Missing --species')
if (!symbol_mode %in% c('gene_id', 'gene_name_if_available')) stop('[ERROR] Invalid --symbol-mode')

dir.create(lib_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(build_dir, recursive = TRUE, showWarnings = FALSE)

message('[INFO] TaxID: ', taxid)
message('[INFO] Genus/species: ', genus, ' ', species)
message('[INFO] GTF/GFF: ', gtf_file)
message('[INFO] eggNOG annotations: ', eggnog_file)
message('[INFO] Build dir: ', build_dir)
message('[INFO] Install lib: ', lib_dir)
message('[INFO] Version: ', version)

# AnnotationForge creates package name from genus/species, usually org.Gspecies.eg.db
pkg_name_guess <- paste0('org.', substr(genus, 1, 1), species, '.eg.db')
pkg_install_dir_guess <- file.path(lib_dir, pkg_name_guess)
if (force && dir.exists(pkg_install_dir_guess)) {
  message('[INFO] Removing previous installed package guess: ', pkg_install_dir_guess)
  unlink(pkg_install_dir_guess, recursive = TRUE, force = TRUE)
}

open_text <- function(path) {
  if (grepl('\\.gz$', path)) gzfile(path, open = 'rt') else file(path, open = 'rt')
}

extract_gtf_attr <- function(x, key) {
  pattern <- paste0(key, ' "([^"]+)"')
  hit <- regexpr(pattern, x, perl = TRUE)
  out <- rep(NA_character_, length(x))
  ok <- hit > 0
  out[ok] <- sub(pattern, '\\1', regmatches(x, hit)[ok], perl = TRUE)
  out
}

extract_gff_attr <- function(x, key) {
  pattern <- paste0('(?:^|;)', key, '=([^;]+)')
  hit <- regexpr(pattern, x, perl = TRUE)
  out <- rep(NA_character_, length(x))
  ok <- hit > 0
  out[ok] <- sub(pattern, '\\1', regmatches(x, hit)[ok], perl = TRUE)
  out <- URLdecode(out)
  out
}

clean_gene_id <- function(x) {
  x <- as.character(x)
  x <- sub('^gene:', '', x)
  x <- sub('^transcript:', '', x)
  x <- sub('^ID=gene:', '', x)
  x <- sub('^ID=', '', x)
  x
}

extract_any_attr <- function(x, keys) {
  out <- rep(NA_character_, length(x))
  for (key in keys) {
    gtf_val <- extract_gtf_attr(x, key)
    gff_val <- extract_gff_attr(x, key)
    val <- ifelse(!is.na(gtf_val) & nzchar(gtf_val), gtf_val, gff_val)
    use <- (is.na(out) | !nzchar(out)) & !is.na(val) & nzchar(val)
    out[use] <- val[use]
  }
  clean_gene_id(out)
}

message('[INFO] Reading GTF/GFF...')
con <- open_text(gtf_file)
gtf <- read.delim(con, header = FALSE, sep = '\t', quote = '', comment.char = '#', stringsAsFactors = FALSE)
close(con)

if (ncol(gtf) < 9) stop('[ERROR] Annotation file has fewer than 9 columns.')
colnames(gtf)[1:9] <- c('seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute')

gene_rows <- gtf[gtf$feature == 'gene', , drop = FALSE]
if (nrow(gene_rows) == 0) {
  message('[WARNING] No feature == gene rows found. Falling back to all rows with gene identifiers.')
  gene_rows <- gtf
}

gene_rows$gene_id <- extract_any_attr(gene_rows$attribute, gene_attr_order)
gene_rows$gene_name <- extract_any_attr(gene_rows$attribute, gene_name_attr_order)
gene_rows$description <- extract_any_attr(gene_rows$attribute, description_attr_order)

gene_rows <- gene_rows[!is.na(gene_rows$gene_id) & nzchar(gene_rows$gene_id), , drop = FALSE]
if (nrow(gene_rows) == 0) stop('[ERROR] No gene IDs extracted. Check --gene-attr-order.')

gene_rows <- gene_rows[order(gene_rows$gene_id, gene_rows$start), , drop = FALSE]
gene_meta <- gene_rows[!duplicated(gene_rows$gene_id), , drop = FALSE]

gene_label <- gene_meta$description
missing_label <- is.na(gene_label) | !nzchar(gene_label)
gene_label[missing_label] <- gene_meta$gene_name[missing_label]
missing_label <- is.na(gene_label) | !nzchar(gene_label)
gene_label[missing_label] <- gene_meta$gene_id[missing_label]

symbol <- gene_meta$gene_id
if (symbol_mode == 'gene_name_if_available') {
  symbol <- gene_meta$gene_name
  missing_symbol <- is.na(symbol) | !nzchar(symbol)
  symbol[missing_symbol] <- gene_meta$gene_id[missing_symbol]
}

gene_info <- unique(data.frame(GID = gene_meta$gene_id, SYMBOL = symbol, GENENAME = gene_label, stringsAsFactors = FALSE))
chromosome <- unique(data.frame(GID = gene_meta$gene_id, CHR = gene_meta$seqname, stringsAsFactors = FALSE))

message('[INFO] Genes in gene_info: ', nrow(gene_info))
message('[INFO] First gene_info rows:')
print(head(gene_info, 5))

message('[INFO] Reading eggNOG annotations...')
eggnog_lines <- readLines(eggnog_file, warn = FALSE)
header_idx <- grep('^#query', eggnog_lines)
if (length(header_idx) == 0) stop('[ERROR] Could not find eggNOG header line starting with #query')
header <- strsplit(sub('^#', '', eggnog_lines[header_idx[1]]), '\t', fixed = TRUE)[[1]]
data_lines <- eggnog_lines[!grepl('^#', eggnog_lines) & nchar(eggnog_lines) > 0]
if (length(data_lines) == 0) stop('[ERROR] No eggNOG data lines found')

eggnog <- read.delim(text = paste(data_lines, collapse = '\n'), sep = '\t', header = FALSE, quote = '', comment.char = '', stringsAsFactors = FALSE, check.names = FALSE)
if (ncol(eggnog) != length(header)) stop('[ERROR] eggNOG header/data mismatch. Header: ', length(header), '; data: ', ncol(eggnog))
colnames(eggnog) <- header

query_col <- if ('query' %in% colnames(eggnog)) 'query' else colnames(eggnog)[1]
go_candidates <- c('GOs', 'GO_terms', 'go_terms', 'GO')
go_col <- go_candidates[go_candidates %in% colnames(eggnog)][1]
if (is.na(go_col)) stop('[ERROR] Could not find GO column in eggNOG file')

message('[INFO] eggNOG query column: ', query_col)
message('[INFO] eggNOG GO column: ', go_col)

query_ids <- as.character(eggnog[[query_col]])
go_raw <- as.character(eggnog[[go_col]])

if (!is.na(query_gene_regex) && nzchar(query_gene_regex)) {
  m <- regexpr(query_gene_regex, query_ids, perl = TRUE)
  gene_ids <- rep(NA_character_, length(query_ids))
  ok <- m > 0
  gene_ids[ok] <- regmatches(query_ids, m)[ok]
  # If regex has a capture group, sub() extracts it; if not, it returns full match.
  gene_ids[ok] <- sub(query_gene_regex, '\\1', query_ids[ok], perl = TRUE)
} else {
  gene_ids <- query_ids
}

gene_ids <- clean_gene_id(gene_ids)

# Genes whose eggNOG query IDs match the GTF/OrgDb GID namespace.
# For Myotis, this should be ENSMLUG...
valid_gene <- !is.na(gene_ids) &
  nzchar(gene_ids) &
  gene_ids %in% gene_info$GID

# GO-specific valid rows.
valid_go <- valid_gene &
  !is.na(go_raw) &
  nzchar(go_raw) &
  go_raw != '-'

eggnog_go <- data.frame(
  GID = gene_ids[valid_go],
  GO_raw = go_raw[valid_go],
  stringsAsFactors = FALSE
)

message('[INFO] eggNOG rows with matching gene ID and GO: ', nrow(eggnog_go))
message('[INFO] Unique genes with raw GO: ', length(unique(eggnog_go$GID)))

if (nrow(eggnog_go) == 0) stop('[ERROR] No eggNOG rows matched GTF gene IDs. Usually fix by running eggNOG on representative FASTA renamed to gene IDs.')

go_list <- vector('list', nrow(eggnog_go))
for (i in seq_len(nrow(eggnog_go))) {
  terms <- trimws(unlist(strsplit(eggnog_go$GO_raw[i], ',', fixed = TRUE)))
  terms <- terms[grepl('^GO:[0-9]{7}$', terms)]
  if (length(terms) > 0) {
    go_list[[i]] <- data.frame(GID = eggnog_go$GID[i], GO = terms, EVIDENCE = 'IEA', stringsAsFactors = FALSE)
  }
}

go <- unique(do.call(rbind, go_list))
if (is.null(go) || nrow(go) == 0) stop('[ERROR] No valid GO annotations extracted.')

message('[INFO] GO-gene pairs: ', nrow(go))
message('[INFO] Unique genes with GO: ', length(unique(go$GID)))
message('[INFO] Unique GO terms: ', length(unique(go$GO)))

# ------------------------------------------------------------
# eggNOG-derived KEGG tables
# These remain in the same ID namespace as the GTF/DEG matrix:
# GID = Ensembl-like gene ID, e.g. ENSMLUG...
# ------------------------------------------------------------

find_eggnog_col <- function(candidates) {
  hit <- candidates[candidates %in% colnames(eggnog)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

split_eggnog_terms <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  out <- trimws(unlist(strsplit(x, ",", fixed = TRUE)))
  out <- out[nzchar(out) & out != "-"]
  unique(out)
}

make_eggnog_term_table <- function(colname, out_col, valid_gene) {
  empty <- data.frame(GID = character(), TERM = character(), stringsAsFactors = FALSE)
  colnames(empty)[2] <- out_col

  if (is.na(colname) || !nzchar(colname) || !colname %in% colnames(eggnog)) {
    message("[WARNING] eggNOG column not found for ", out_col)
    return(empty)
  }

  raw <- as.character(eggnog[[colname]])
  valid <- valid_gene &
    !is.na(raw) &
    nzchar(raw) &
    raw != "-"

  if (!any(valid)) {
    message("[WARNING] No valid eggNOG entries for ", out_col, " using column ", colname)
    return(empty)
  }

  idx <- which(valid)
  rows <- vector("list", length(idx))
  j <- 0

  for (i in idx) {
    terms <- split_eggnog_terms(raw[i])

    if (length(terms) > 0) {
      j <- j + 1
      rows[[j]] <- data.frame(
        GID = gene_ids[i],
        TERM = terms,
        stringsAsFactors = FALSE
      )
    }
  }

  if (j == 0) return(empty)

  rows <- rows[seq_len(j)]
  out <- unique(do.call(rbind, rows))
  colnames(out)[2] <- out_col
  out
}

kegg_ko_col <- find_eggnog_col(c("KEGG_ko", "KEGG_KO", "KEGG KO"))
kegg_pathway_col <- find_eggnog_col(c("KEGG_Pathway", "KEGG_pathway", "KEGG_PATHWAY"))
kegg_module_col <- find_eggnog_col(c("KEGG_Module", "KEGG_module", "KEGG_MODULE"))
kegg_reaction_col <- find_eggnog_col(c("KEGG_Reaction", "KEGG_reaction", "KEGG_REACTION"))
kegg_brite_col <- find_eggnog_col(c("BRITE", "KEGG_BRITE", "KEGG_Brite", "KEGG_brite"))

message("[INFO] eggNOG KEGG KO column: ", kegg_ko_col)
message("[INFO] eggNOG KEGG Pathway column: ", kegg_pathway_col)
message("[INFO] eggNOG KEGG Module column: ", kegg_module_col)
message("[INFO] eggNOG KEGG Reaction column: ", kegg_reaction_col)
message("[INFO] eggNOG KEGG BRITE column: ", kegg_brite_col)

kegg_ko <- make_eggnog_term_table(
  colname = kegg_ko_col,
  out_col = "KEGG_KO",
  valid_gene = valid_gene
)

kegg_pathway <- make_eggnog_term_table(
  colname = kegg_pathway_col,
  out_col = "KEGG_PATHWAY",
  valid_gene = valid_gene
)

kegg_module <- make_eggnog_term_table(
  colname = kegg_module_col,
  out_col = "KEGG_MODULE",
  valid_gene = valid_gene
)

kegg_reaction <- make_eggnog_term_table(
  colname = kegg_reaction_col,
  out_col = "KEGG_REACTION",
  valid_gene = valid_gene
)

kegg_brite <- make_eggnog_term_table(
  colname = kegg_brite_col,
  out_col = "KEGG_BRITE",
  valid_gene = valid_gene
)

message("[INFO] KEGG KO gene-term pairs: ", nrow(kegg_ko))
message("[INFO] KEGG Pathway gene-term pairs: ", nrow(kegg_pathway))
message("[INFO] KEGG Module gene-term pairs: ", nrow(kegg_module))
message("[INFO] KEGG Reaction gene-term pairs: ", nrow(kegg_reaction))
message("[INFO] KEGG BRITE gene-term pairs: ", nrow(kegg_brite))

# Audit tables
if (force && dir.exists(file.path(build_dir, paste0('org.', substr(genus, 1, 1), species, '.eg.db')))) {
  unlink(file.path(build_dir, paste0('org.', substr(genus, 1, 1), species, '.eg.db')), recursive = TRUE, force = TRUE)
}
audit_dir <- file.path(build_dir, 'audit_tables')
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(gene_info, file.path(audit_dir, 'gene_info.csv'), row.names = FALSE)
write.csv(chromosome, file.path(audit_dir, 'chromosome.csv'), row.names = FALSE)
write.csv(go, file.path(audit_dir, 'go.csv'), row.names = FALSE)
write.csv(kegg_ko, file.path(audit_dir, 'kegg_ko.csv'), row.names = FALSE)
write.csv(kegg_pathway, file.path(audit_dir, 'kegg_pathway.csv'), row.names = FALSE)
write.csv(kegg_module, file.path(audit_dir, 'kegg_module.csv'), row.names = FALSE)
write.csv(kegg_reaction, file.path(audit_dir, 'kegg_reaction.csv'), row.names = FALSE)
write.csv(kegg_brite, file.path(audit_dir, 'kegg_brite.csv'), row.names = FALSE)

message('[INFO] Building OrgDb package...')

orgdb_tables <- list(
  gene_info = gene_info,
  chromosome = chromosome,
  go = go
)

# Only add non-empty KEGG-derived tables.
# AnnotationForge accepts extra data frames and creates matching columns/keytypes.
if (nrow(kegg_ko) > 0) {
  orgdb_tables$kegg_ko <- kegg_ko
}

if (nrow(kegg_pathway) > 0) {
  orgdb_tables$kegg_pathway <- kegg_pathway
}

if (nrow(kegg_module) > 0) {
  orgdb_tables$kegg_module <- kegg_module
}

if (nrow(kegg_reaction) > 0) {
  orgdb_tables$kegg_reaction <- kegg_reaction
}

if (nrow(kegg_brite) > 0) {
  orgdb_tables$kegg_brite <- kegg_brite
}

make_args <- c(
  orgdb_tables,
  list(
    version = version,
    maintainer = maintainer,
    author = author,
    outputDir = build_dir,
    tax_id = taxid,
    genus = genus,
    species = species,
    goTable = 'go',
    verbose = TRUE
  )
)

oldwd <- getwd()
setwd(build_dir)
pkg_path <- do.call(AnnotationForge::makeOrgPackage, make_args)
setwd(oldwd)

if (!dir.exists(pkg_path)) {
  dirs <- list.dirs(build_dir, recursive = FALSE, full.names = TRUE)
  org_dirs <- dirs[grepl('org\\..*\\.eg\\.db$', basename(dirs))]
  if (length(org_dirs) == 0) stop('[ERROR] Could not find generated OrgDb source directory')
  pkg_path <- org_dirs[which.max(file.info(org_dirs)$mtime)]
}

message('[INFO] Package source: ', pkg_path)
cmd_args <- c('CMD', 'INSTALL', '-l', lib_dir, pkg_path)
message('[RUN] R ', paste(shQuote(cmd_args), collapse = ' '))
status <- system2('R', args = cmd_args)
if (status != 0) stop('[ERROR] R CMD INSTALL failed')

.libPaths(c(lib_dir, .libPaths()))
pkg_name <- basename(pkg_path)
message('[OK] Installed package: ', pkg_name, ' into ', lib_dir)

suppressPackageStartupMessages({
  library(AnnotationDbi)
})

suppressPackageStartupMessages(library(pkg_name, character.only = TRUE))
db <- get(pkg_name)

message('[INFO] Loaded package from:')
print(find.package(pkg_name))
message('[INFO] Package version:')
print(packageVersion(pkg_name))

test_genes <- head(gene_info$GID, 10)
valid_symbol <- test_genes %in% keys(db, keytype = 'SYMBOL')
message('[INFO] Are first GTF genes valid SYMBOL keys?')
print(valid_symbol)
if (!all(valid_symbol)) stop('[ERROR] Installed OrgDb does not recognize test genes as SYMBOL')

available_cols <- AnnotationDbi::columns(db)

message('[INFO] Available OrgDb columns:')
print(available_cols)

kegg_cols <- intersect(
  c('KEGG_KO', 'KEGG_PATHWAY', 'KEGG_MODULE', 'KEGG_REACTION', 'KEGG_BRITE'),
  available_cols
)

message('[INFO] KEGG-derived columns found:')
print(kegg_cols)

select_cols <- intersect(
  c(
    'SYMBOL',
    'GID',
    'GENENAME',
    'GO',
    'ONTOLOGY',
    'EVIDENCE',
    'KEGG_KO',
    'KEGG_PATHWAY',
    'KEGG_MODULE',
    'KEGG_REACTION',
    'KEGG_BRITE'
  ),
  available_cols
)

tab <- AnnotationDbi::select(
  db,
  keys = test_genes,
  keytype = 'SYMBOL',
  columns = select_cols
)

message('[INFO] select() test:')
print(head(tab, 20))
message('[OK] Gold-standard custom OrgDb finished.')
message('[IMPORTANT] Use with:')
message('export R_LIBS_USER="', lib_dir, ':${R_LIBS_USER:-}"')
