#!/usr/bin/env Rscript

# ============================================================
# MTD detected species pie chart by taxonomic level
# ------------------------------------------------------------
# Input:
#   1) detected_microbiome_species_ranked_by_importance.tsv
#   2) Combined.mpa or graphlan_input.clean.tsv
#
# Supports:
#   --category_level auto
#   --category_level major
#   --category_level superkingdom
#   --category_level kingdom
#   --category_level phylum
#   --category_level class
#   --category_level order
#   --category_level family
#   --category_level genus
#   --category_level species
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  idx <- which(args == flag)

  if (length(idx) == 0) {
    return(default)
  }

  idx <- idx[1]

  if (idx == length(args)) {
    stop("Missing value for ", flag, call. = FALSE)
  }

  args[idx + 1]
}

has_flag <- function(flag) {
  flag %in% args
}

print_help <- function() {
  cat("
Usage:
  Rscript detected_species_pie_by_phylum.R \\
    --ranked FILE \\
    --taxonomy FILE \\
    --output_prefix PREFIX

Required:
  --ranked FILE
      detected_microbiome_species_ranked_by_importance.tsv

  --taxonomy FILE
      Combined.mpa or graphlan_input.clean.tsv

  --output_prefix PREFIX
      Output prefix without extension

Optional:
  --category_level LEVEL
      auto, major, superkingdom, kingdom, phylum, class, order,
      family, genus, or species.
      Default: phylum

  --top_n VALUE
      auto, 0/all, or positive integer.
      Default: auto

  --cumulative_percent NUM
      Used when --top_n auto.
      Default: 95

  --min_percent NUM
      Used when --top_n auto.
      Default: 2

  --min_slices NUM
      Used when --top_n auto.
      Default: 8

  --max_slices NUM
      Used when --top_n auto.
      Default: 15

  --auto_min_categories NUM
      Used when --category_level auto.
      Default: 4

  --auto_max_categories NUM
      Used when --category_level auto.
      Default: 18

  --auto_target_categories NUM
      Used when --category_level auto.
      Default: 12

  --title TEXT
      Plot title.
      Default: auto

  --legend_title TEXT
      Legend title.
      Default: auto

  --panel_label TEXT
      Panel label, e.g. A.
      Default: empty

  --width NUM
      Width in inches.
      Default: 9

  --height NUM
      Height in inches.
      Default: 7

Examples:
  Rscript detected_species_pie_by_phylum.R \\
    --ranked detected_microbiome_species_ranked_by_importance.tsv \\
    --taxonomy Combined.mpa \\
    --category_level phylum \\
    --output_prefix detected_species_by_phylum

  Rscript detected_species_pie_by_phylum.R \\
    --ranked detected_microbiome_species_ranked_by_importance.tsv \\
    --taxonomy Combined.mpa \\
    --category_level auto \\
    --output_prefix detected_species_by_auto_level
")
  quit(status = 0)
}

if (length(args) == 0 || has_flag("--help") || has_flag("-h")) {
  print_help()
}

ranked_file <- get_arg("--ranked")
taxonomy_file <- get_arg("--taxonomy")
output_prefix <- get_arg("--output_prefix")

if (is.null(ranked_file) || !file.exists(ranked_file)) {
  stop("Missing or invalid --ranked file.", call. = FALSE)
}

if (is.null(taxonomy_file) || !file.exists(taxonomy_file)) {
  stop("Missing or invalid --taxonomy file.", call. = FALSE)
}

if (is.null(output_prefix)) {
  stop("Missing --output_prefix.", call. = FALSE)
}

category_level <- tolower(get_arg("--category_level", "phylum"))
top_n_arg <- get_arg("--top_n", "auto")

cumulative_percent <- as.numeric(get_arg("--cumulative_percent", "95"))
min_percent <- as.numeric(get_arg("--min_percent", "2"))
min_slices <- as.integer(get_arg("--min_slices", "8"))
max_slices <- as.integer(get_arg("--max_slices", "15"))

auto_min_categories <- as.integer(get_arg("--auto_min_categories", "4"))
auto_max_categories <- as.integer(get_arg("--auto_max_categories", "18"))
auto_target_categories <- as.integer(get_arg("--auto_target_categories", "12"))

title_arg <- get_arg("--title", "auto")
legend_title_arg <- get_arg("--legend_title", "auto")
panel_label <- get_arg("--panel_label", "")
width <- as.numeric(get_arg("--width", "9"))
height <- as.numeric(get_arg("--height", "7"))

allowed_levels <- c(
  "auto",
  "major",
  "domain",
  "superkingdom",
  "kingdom",
  "phylum",
  "class",
  "order",
  "family",
  "genus",
  "species"
)

if (!category_level %in% allowed_levels) {
  stop(
    "--category_level must be one of: ",
    paste(allowed_levels, collapse = ", "),
    call. = FALSE
  )
}

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop(
    "Package ggplot2 is required. Install it with:\n",
    "install.packages('ggplot2')",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(ggplot2)
})

clean_name <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- sub("^[A-Za-z]__", "", x)
  x <- gsub("_", " ", x)
  x <- gsub("[[:space:]]+", " ", x)
  trimws(x)
}

clean_taxon_one <- function(x) {
  x <- trimws(as.character(x))

  if (is.na(x) || x == "") {
    return("")
  }

  if (grepl("[|;]", x)) {
    parts <- unlist(strsplit(x, "[|;]"))
    parts <- trimws(parts)

    sp <- parts[grepl("^s__", parts)]

    if (length(sp) > 0) {
      x <- sp[length(sp)]
    } else {
      x <- parts[length(parts)]
    }
  }

  clean_name(x)
}

clean_taxon <- function(x) {
  vapply(x, clean_taxon_one, character(1))
}

guess_col <- function(df, candidates, fallback = 1) {
  nms <- colnames(df)

  hit <- candidates[candidates %in% nms]

  if (length(hit) > 0) {
    return(hit[1])
  }

  nms[fallback]
}

get_rank_value <- function(parts, pattern, use_last = FALSE) {
  hit <- parts[grepl(pattern, parts)]

  if (length(hit) == 0) {
    return(NA_character_)
  }

  if (use_last) {
    return(clean_name(hit[length(hit)]))
  }

  clean_name(hit[1])
}

normalize_major <- function(superkingdom,
                            kingdom,
                            phylum,
                            class,
                            order,
                            family,
                            genus,
                            species) {
  ranks <- c(
    superkingdom = superkingdom,
    kingdom = kingdom,
    phylum = phylum,
    class = class,
    order = order,
    family = family,
    genus = genus,
    species = species
  )

  ranks <- as.character(ranks)
  ranks[is.na(ranks)] <- ""
  ranks_low <- tolower(ranks)

  any_rank_matches <- function(pattern) {
    any(grepl(pattern, ranks_low, ignore.case = TRUE))
  }

  any_rank_in <- function(values) {
    any(ranks_low %in% tolower(values))
  }

  # Already broad labels
  if (any_rank_matches("bacteria")) return("Bacteria")
  if (any_rank_matches("archaea")) return("Archaea")
  if (any_rank_matches("eukary|eukarya")) return("Eukaryota")
  if (any_rank_matches("virus|viruses|viral|virae")) return("Viruses")

  # Eukaryotic kingdoms, phyla, classes, and major clades
  if (any_rank_in(c(
    "Fungi", "Metazoa", "Viridiplantae", "Plantae",
    "Protozoa", "Chromista", "Alveolata", "Stramenopiles",
    "Rhizaria", "Amoebozoa", "Excavata",
    "Opisthokonta", "Holozoa", "Filasterea",
    "Cryptista", "Cryptophyta", "Cryptophyceae",
    "Metamonada", "Parabasalia", "Parabasaliae",
    "Ascomycota", "Basidiomycota", "Apicomplexa",
    "Nematoda", "Arthropoda", "Chordata", "Mollusca",
    "Annelida", "Platyhelminthes", "Cnidaria", "Porifera",
    "Streptophyta", "Chlorophyta", "Ciliophora",
    "Cercozoa", "Euglenozoa", "Oomycota", "Mucoromycota",
    "Microsporidia", "Fornicata",
    "Cryptomonadales", "Pyrenomonadales",
    "Cryptomonadaceae", "Geminigeraceae",
    "Capsaspora", "Cryptomonas", "Guillardia", "Trichomonas"
  ))) {
    return("Eukaryota")
  }

  # Bacterial high-level groups and phyla
  if (any_rank_in(c(
    "Pseudomonadati", "Bacillati", "Terrabacteria group",
    "FCB group", "PVC group", "Fusobacteriati",
    "Thermotogati", "Aquificati", "Spirochaetati",
    "Acidobacteriota", "Actinomycetota", "Aquificota",
    "Bacillota", "Bacteroidota", "Bdellovibrionota",
    "Campylobacterota", "Chlamydiota", "Chlorobiota",
    "Chloroflexota", "Cyanobacteriota", "Deinococcota",
    "Desulfobacterota", "Fibrobacterota", "Fusobacteriota",
    "Gemmatimonadota", "Mycoplasmatota", "Nitrospirota",
    "Planctomycetota", "Pseudomonadota", "Spirochaetota",
    "Synergistota", "Thermodesulfobacteriota", "Thermotogota",
    "Verrucomicrobiota", "Thermosulfidibacterota"
  ))) {
    return("Bacteria")
  }

  # Archaeal high-level groups and phyla
  if (any_rank_in(c(
    "Methanobacteriati", "Thermoproteati", "Halobacteriati",
    "Thermoplasmati", "Nanoarchaeota", "Asgard group",
    "Euryarchaeota", "Crenarchaeota", "Thermoproteota",
    "Halobacteriota", "Methanobacteriota", "Thermoplasmatota",
    "Asgardarchaeota", "Thaumarchaeota", "Korarchaeota"
  ))) {
    return("Archaea")
  }

  # Viral kingdom/realm/phylum-like labels
  if (any_rank_in(c(
    "Heunggongvirae", "Shotokuvirae", "Bamfordvirae",
    "Orthornavirae", "Pararnavirae", "Trapavirae",
    "Artverviricota", "Cossaviricota", "Duplornaviricota",
    "Hofneiviricota", "Kitrinoviricota", "Lenarviricota",
    "Negarnaviricota", "Nucleocytoviricota", "Peploviricota",
    "Pisuviricota", "Preplasmiviricota", "Uroviricota"
  ))) {
    return("Viruses")
  }

  if (any_rank_matches("virus|phage|viridae")) {
    return("Viruses")
  }

  "Unknown / unmapped"
}

extract_taxonomy_map <- function(taxonomy_file) {
  lines <- readLines(taxonomy_file, warn = FALSE)
  lines <- lines[trimws(lines) != ""]

  out <- list()
  k <- 1

  for (line in lines) {
    fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
    tax <- fields[1]

    if (tax %in% c("ID", "#Classification", "classification")) {
      next
    }

    parts <- unlist(strsplit(tax, "[|;]"))
    parts <- trimws(parts)

    species <- get_rank_value(parts, "^s__", use_last = TRUE)

    if (is.na(species) || species == "") {
      next
    }

    # d__ is true domain/superkingdom when available.
    # Do NOT replace missing d__ with k__, otherwise kingdom-level
    # groups such as Fungi, Metazoa, Bacillati, or Pseudomonadati
    # are incorrectly labeled as superkingdom.
    superkingdom <- get_rank_value(parts, "^d__", use_last = FALSE)
    kingdom <- get_rank_value(parts, "^k__", use_last = FALSE)

    phylum <- get_rank_value(parts, "^p__", use_last = FALSE)
    class <- get_rank_value(parts, "^c__", use_last = FALSE)
    order <- get_rank_value(parts, "^o__", use_last = FALSE)
    family <- get_rank_value(parts, "^f__", use_last = FALSE)
    genus <- get_rank_value(parts, "^g__", use_last = FALSE)

    major <- normalize_major(
      superkingdom = superkingdom,
      kingdom = kingdom,
      phylum = phylum,
      class = class,
      order = order,
      family = family,
      genus = genus,
      species = species
    )

    out[[k]] <- data.frame(
      taxon = species,
      major = major,
      domain = superkingdom,
      superkingdom = superkingdom,
      kingdom = kingdom,
      phylum = phylum,
      class = class,
      order = order,
      family = family,
      genus = genus,
      species = species,
      stringsAsFactors = FALSE
    )

    k <- k + 1
  }

  if (length(out) == 0) {
    stop(
      "No species-level taxonomy paths could be extracted from: ",
      taxonomy_file,
      call. = FALSE
    )
  }

  tax_map <- do.call(rbind, out)

  rank_cols <- c(
    "major", "domain", "superkingdom", "kingdom",
    "phylum", "class", "order", "family", "genus", "species"
  )

  tax_map$filled_ranks <- rowSums(
    !is.na(tax_map[, rank_cols, drop = FALSE]) &
      tax_map[, rank_cols, drop = FALSE] != ""
  )

  tax_map <- tax_map[order(tax_map$taxon, -tax_map$filled_ranks), ]
  tax_map <- tax_map[!duplicated(tax_map$taxon), , drop = FALSE]
  tax_map$filled_ranks <- NULL

  tax_map
}

get_category <- function(df, level) {
  if (level == "species") {
    out <- df$taxon
  } else {
    out <- df[[level]]
  }

  out <- as.character(out)
  missing <- is.na(out) | trimws(out) == ""

  # If the requested rank is missing, keep a biologically useful label
  # instead of generic Unknown/unmapped. Example for phylum:
  # Eukaryota; no phylum assignment
  if (any(missing) && level != "major" && "major" %in% colnames(df)) {
    major <- as.character(df$major)
    major[is.na(major) | trimws(major) == ""] <- "Unknown / unmapped"

    has_major <- major != "Unknown / unmapped"

    out[missing & has_major] <- paste0(
      major[missing & has_major],
      "; no ",
      level,
      " assignment"
    )

    out[missing & !has_major] <- "Unknown / unmapped"
  } else {
    out[missing] <- "Unknown / unmapped"
  }

  out
}

level_label <- function(level) {
  switch(
    level,
    major = "Major Taxonomic Group",
    domain = "Domain",
    superkingdom = "Superkingdom",
    kingdom = "Kingdom",
    phylum = "Phylum",
    class = "Class",
    order = "Order",
    family = "Family",
    genus = "Genus",
    species = "Species",
    auto = "Auto-selected Taxonomic Level",
    level
  )
}

choose_category_level_auto <- function(df,
                                       min_categories = 4,
                                       max_categories = 18,
                                       target_categories = 12) {
  candidate_levels <- c(
    "species",
    "genus",
    "family",
    "order",
    "class",
    "phylum",
    "kingdom",
    "superkingdom"
  )

  stats <- data.frame(
    level = character(),
    n_categories = integer(),
    n_unknown = integer(),
    specificity = integer(),
    stringsAsFactors = FALSE
  )

  for (i in seq_along(candidate_levels)) {
    lvl <- candidate_levels[i]
    cat <- get_category(df, lvl)

    n_cat <- length(unique(cat))
    n_unknown <- sum(cat == "Unknown / unmapped")

    stats <- rbind(
      stats,
      data.frame(
        level = lvl,
        n_categories = n_cat,
        n_unknown = n_unknown,
        specificity = length(candidate_levels) - i + 1,
        stringsAsFactors = FALSE
      )
    )
  }

  stats$distance_to_target <- abs(stats$n_categories - target_categories)

  good <- stats[
    stats$n_categories >= min_categories &
      stats$n_categories <= max_categories,
    ,
    drop = FALSE
  ]

  if (nrow(good) > 0) {
    good <- good[order(good$distance_to_target, -good$specificity), ]
    chosen <- good$level[1]
  } else {
    # Biological fallback: if nothing falls in the desired range,
    # phylum is usually the most interpretable compromise.
    if ("phylum" %in% stats$level && stats$n_categories[stats$level == "phylum"] > 1) {
      chosen <- "phylum"
    } else {
      stats <- stats[order(stats$distance_to_target, -stats$specificity), ]
      chosen <- stats$level[1]
    }
  }

  attr(chosen, "stats") <- stats
  chosen
}

choose_keep_n <- function(df,
                          top_n_arg = "auto",
                          cumulative_percent = 95,
                          min_percent = 2,
                          min_slices = 8,
                          max_slices = 15) {
  n_cat <- nrow(df)

  if (n_cat == 0) {
    return(0)
  }

  top_n_arg_lower <- tolower(as.character(top_n_arg))

  if (top_n_arg_lower %in% c("0", "all", "none", "false", "no")) {
    return(n_cat)
  }

  if (top_n_arg_lower != "auto") {
    forced_n <- suppressWarnings(as.integer(top_n_arg_lower))

    if (is.na(forced_n) || forced_n < 1) {
      stop("--top_n must be auto, 0/all, or a positive integer.", call. = FALSE)
    }

    return(min(forced_n, n_cat))
  }

  if (n_cat <= max_slices) {
    return(n_cat)
  }

  total <- sum(df$count, na.rm = TRUE)

  if (total <= 0) {
    return(min(max_slices, n_cat))
  }

  pct <- 100 * df$count / total
  cum_pct <- cumsum(pct)

  n_to_reach_cumulative <- which(cum_pct >= cumulative_percent)[1]

  if (is.na(n_to_reach_cumulative)) {
    n_to_reach_cumulative <- n_cat
  }

  n_above_min_percent <- sum(pct >= min_percent)

  keep_n <- max(
    min_slices,
    n_to_reach_cumulative,
    n_above_min_percent
  )

  keep_n <- min(keep_n, max_slices, n_cat)

  keep_n
}

make_labels <- function(df) {
  df$count <- as.numeric(df$count)
  total <- sum(df$count, na.rm = TRUE)

  if (total > 0) {
    df$percent <- 100 * df$count / total
  } else {
    df$percent <- 0
  }

  count_txt <- formatC(
    round(df$count),
    format = "d",
    big.mark = ","
  )

  df$legend_label <- paste0(
    df$category,
    " (",
    count_txt,
    " - ",
    sprintf("%.1f", df$percent),
    "%)"
  )

  df
}

# ============================================================
# Read input
# ============================================================

ranked <- read.delim(
  ranked_file,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  quote = "",
  comment.char = "",
  stringsAsFactors = FALSE
)

taxon_col <- guess_col(
  ranked,
  candidates = c("taxon", "Taxon", "name", "Name", "species", "Species"),
  fallback = 1
)

ranked_taxa <- unique(clean_taxon(ranked[[taxon_col]]))
ranked_taxa <- ranked_taxa[ranked_taxa != ""]

tax_map <- extract_taxonomy_map(taxonomy_file)

merged <- merge(
  data.frame(taxon = ranked_taxa, stringsAsFactors = FALSE),
  tax_map,
  by = "taxon",
  all.x = TRUE
)

rank_cols <- c(
  "major", "domain", "superkingdom", "kingdom",
  "phylum", "class", "order", "family", "genus", "species"
)

for (cc in rank_cols) {
  if (!cc %in% colnames(merged)) {
    merged[[cc]] <- NA_character_
  }
}

# Species always comes from the detected ranked table.
merged$species <- merged$taxon

if (category_level == "auto") {
  chosen <- choose_category_level_auto(
    df = merged,
    min_categories = auto_min_categories,
    max_categories = auto_max_categories,
    target_categories = auto_target_categories
  )

  auto_stats <- attr(chosen, "stats")
  category_level_final <- as.character(chosen)

  auto_stats_file <- paste0(output_prefix, "_auto_level_candidates.tsv")

  dir.create(dirname(auto_stats_file), recursive = TRUE, showWarnings = FALSE)

  write.table(
    auto_stats,
    file = auto_stats_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  cat("[INFO] Auto category-level selection:\n")
  cat("  selected:", category_level_final, "\n")
  cat("  candidates table:", auto_stats_file, "\n")
} else if (category_level == "domain") {
  category_level_final <- "domain"
} else {
  category_level_final <- category_level
}

merged$category <- get_category(merged, category_level_final)

unmapped <- merged[
  merged$category == "Unknown / unmapped" |
    is.na(merged$category) |
    merged$category == "",
  ,
  drop = FALSE
]

# ============================================================
# Summarize
# ============================================================

summary_full <- aggregate(
  taxon ~ category,
  data = merged,
  FUN = length
)

colnames(summary_full) <- c("category", "count")

summary_full <- summary_full[
  order(summary_full$count, decreasing = TRUE),
  ,
  drop = FALSE
]

summary_full <- make_labels(summary_full)

full_summary_file <- paste0(output_prefix, "_full_summary.tsv")
unmapped_file <- paste0(output_prefix, "_unmapped_species.tsv")

dir.create(dirname(full_summary_file), recursive = TRUE, showWarnings = FALSE)

write.table(
  summary_full,
  file = full_summary_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

if (nrow(unmapped) > 0) {
  write.table(
    unmapped[, c("taxon"), drop = FALSE],
    file = unmapped_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
} else {
  write.table(
    data.frame(taxon = character()),
    file = unmapped_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

plot_df <- summary_full[, c("category", "count"), drop = FALSE]

keep_n <- choose_keep_n(
  df = plot_df,
  top_n_arg = top_n_arg,
  cumulative_percent = cumulative_percent,
  min_percent = min_percent,
  min_slices = min_slices,
  max_slices = max_slices
)

cat("[INFO] Pie chart category display:\n")
cat("  requested category_level:", category_level, "\n")
cat("  final category_level:", category_level_final, "\n")
cat("  categories before collapsing:", nrow(plot_df), "\n")
cat("  categories kept:", keep_n, "\n")
cat("  top_n:", top_n_arg, "\n")

if (keep_n < nrow(plot_df)) {
  keep <- plot_df[seq_len(keep_n), , drop = FALSE]

  other <- data.frame(
    category = "Other",
    count = sum(plot_df$count[-seq_len(keep_n)], na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  plot_df <- rbind(keep, other)
}

plot_df <- make_labels(plot_df)

summary_file <- paste0(output_prefix, "_summary.tsv")

write.table(
  plot_df,
  file = summary_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

total_species <- sum(plot_df$count, na.rm = TRUE)

if (tolower(title_arg) == "auto") {
  title_text <- paste0(
    "Distribution of Detected Species by ",
    level_label(category_level_final)
  )
} else {
  title_text <- title_arg
}

if (tolower(legend_title_arg) == "auto") {
  legend_title <- level_label(category_level_final)
} else {
  legend_title <- legend_title_arg
}

subtitle_text <- paste0(
  "Detected species: ",
  formatC(total_species, format = "d", big.mark = ","),
  " | Taxonomic level: ",
  level_label(category_level_final)
)

caption_text <- paste0(
  "Detected species were grouped at the ",
  level_label(category_level_final),
  " level. ",
  "Low-frequency categories may be collapsed into Other."
)

plot_df$category <- factor(plot_df$category, levels = plot_df$category)
plot_df$legend_label <- factor(plot_df$legend_label, levels = plot_df$legend_label)

# ============================================================
# Plot
# ============================================================

p <- ggplot(
  plot_df,
  aes(
    x = "",
    y = count,
    fill = legend_label
  )
) +
  geom_col(
    width = 1,
    color = "white",
    linewidth = 0.35
  ) +
  coord_polar(theta = "y") +
  labs(
    title = title_text,
    subtitle = subtitle_text,
    caption = caption_text,
    fill = legend_title,
    tag = panel_label
  ) +
  theme_void(base_size = 14) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 18,
      margin = margin(b = 8)
    ),
    plot.subtitle = element_text(
      hjust = 0.5,
      size = 11,
      margin = margin(b = 12)
    ),
    plot.caption = element_text(
      hjust = 0.5,
      size = 9,
      color = "grey35",
      margin = margin(t = 12)
    ),
    plot.tag = element_text(
      face = "bold",
      size = 20
    ),
    plot.tag.position = c(0.02, 0.98),
    legend.position = "right",
    legend.title = element_text(
      face = "bold",
      size = 11
    ),
    legend.text = element_text(size = 9),
    plot.margin = margin(18, 18, 18, 18)
  )

png_file <- paste0(output_prefix, ".png")
pdf_file <- paste0(output_prefix, ".pdf")
svg_file <- paste0(output_prefix, ".svg")

dir.create(dirname(png_file), recursive = TRUE, showWarnings = FALSE)

png(
  png_file,
  width = width,
  height = height,
  units = "in",
  res = 300,
  bg = "white"
)
print(p)
dev.off()

pdf(
  pdf_file,
  width = width,
  height = height,
  bg = "white"
)
print(p)
dev.off()

svg(
  svg_file,
  width = width,
  height = height,
  bg = "white"
)
print(p)
dev.off()

cat("[OK] Files generated:\n")
cat("  ", png_file, "\n", sep = "")
cat("  ", pdf_file, "\n", sep = "")
cat("  ", svg_file, "\n", sep = "")
cat("  ", summary_file, "\n", sep = "")
cat("  ", full_summary_file, "\n", sep = "")
cat("  ", unmapped_file, "\n", sep = "")
