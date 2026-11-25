args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: Rscript humann_ID_translation_adjusted.R KEGG.tsv GO.tsv dictionary_dir")

translation <- function(tsv, dictionary, pattern) {
  id <- read.table(tsv, sep = "\t", comment.char = "", header = TRUE,
                   check.names = FALSE, stringsAsFactors = FALSE, quote = "")
  names(id)[1] <- "Gene_Family"
  names(id)[-1] <- sub("_Abundance-RPKs$", "", names(id)[-1])
  original <- as.character(id$Gene_Family)
  if (anyNA(original) || any(!nzchar(original)) || anyDuplicated(original)) {
    stop("Invalid or duplicate HUMAnN feature IDs: ", tsv)
  }
  # Extract canonical IDs from the unstratified portion only.
  base_id <- sub("\\|.*$", "", original)
  eligible <- grepl(pattern, base_id)
  key <- sub("^ko:", "", rownames(dictionary))
  if (anyDuplicated(key)) stop("Duplicate dictionary IDs: ", tsv)
  description <- as.character(dictionary[[1]])[match(base_id, key)]
  translated <- eligible & !is.na(description) & nzchar(description)
  # Keep unknown and special IDs intact. Include the original stratified ID
  # in translated labels so equal descriptions cannot hide distinct features.
  labels <- original
  labels[translated] <- paste0(description[translated], " [", original[translated], "]")
  if (anyDuplicated(labels)) stop("Duplicate translated feature labels: ", tsv)
  id$Gene_Family <- labels
  output <- sub("\\.tsv$", "_translated.tsv", tsv)
  if (identical(output, tsv)) stop("Expected .tsv input: ", tsv)
  write.table(id, output, sep = "\t", row.names = FALSE, quote = TRUE, qmethod = "double")
  mapping <- data.frame(Gene_Family = labels, original_id = original,
                        functional_id = base_id, eligible_for_functional_de = eligible,
                        dictionary_match = translated, stringsAsFactors = FALSE)
  write.table(mapping, sub("\\.tsv$", "_id_map.tsv", output),
              sep = "\t", row.names = FALSE, quote = TRUE, qmethod = "double")
  message("[HUMAnN translation] ", basename(tsv), ": ", sum(translated),
          " translated; ", sum(eligible & !translated), " canonical IDs without labels; ",
          sum(!eligible), " non-functional rows retained in complete table.")
}

terms <- read.table(file.path(args[3], "koterms.txt"), sep = "\t", quote = "", stringsAsFactors = FALSE)
translation(args[1], terms, "^K[0-9]{5}$")
terms <- read.table(file.path(args[3], "goterms.txt"), sep = "\t", quote = "", stringsAsFactors = FALSE)
translation(args[2], terms, "^GO:[0-9]{7}$")
