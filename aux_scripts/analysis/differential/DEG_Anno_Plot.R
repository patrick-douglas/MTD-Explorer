#!/usr/bin/Rscript
#options(echo = TRUE)
httr::set_config(httr::config(ssl_verifypeer = FALSE))
args = commandArgs(trailingOnly=TRUE) # Passing arguments to an R script from bash/shell command lines

# Make the custom OrgDb library visible even when this script is run
# directly, outside MTD_explorer.sh.
if (length(args) >= 4) {
  mtd_root <- dirname(
    normalizePath(
      args[4],
      mustWork = FALSE
    )
  )

  custom_r_lib <- file.path(
    mtd_root,
    "custom_R_libs"
  )

  if (dir.exists(custom_r_lib)) {
    .libPaths(
      unique(
        c(
          custom_r_lib,
          .libPaths()
        )
      )
    )

    message(
      "[INFO] Custom R library enabled: ",
      custom_r_lib
    )
  }

  message(
    "[INFO] Active R library paths: ",
    paste(.libPaths(), collapse = "; ")
  )
}


# MTD Explorer helper functions for robust Kraken host summaries

find_host_summary_path <- function(start_dir, max_depth = 6) {
    current_dir <- normalizePath(
        start_dir,
        mustWork = FALSE
    )

    for (depth in seq_len(max_depth)) {
        candidate <- file.path(
            current_dir,
            "kraken",
            "kraken_host_summary.tsv"
        )

        if (file.exists(candidate)) {
            return(candidate)
        }

        parent_dir <- dirname(current_dir)

        if (identical(parent_dir, current_dir)) {
            break
        }

        current_dir <- parent_dir
    }

    return(NA_character_)
}

read_host_report_counts <- function(report_file) {
    tab <- read.table(
        report_file,
        sep = "	",
        quote = "",
        fill = TRUE,
        comment.char = "",
        stringsAsFactors = FALSE
    )

    if (ncol(tab) < 5) {
        stop("Invalid Kraken2 host report: ", report_file)
    }

    unclassified <- suppressWarnings(
        as.numeric(tab[tab[[4]] == "U", 2][1])
    )

    classified_root <- suppressWarnings(
        as.numeric(tab[tab[[4]] == "R" & tab[[5]] == 1, 2][1])
    )

    if (is.na(unclassified)) {
        unclassified <- 0
    }

    if (is.na(classified_root)) {
        classified_root <- 0
    }

    total <- classified_root + unclassified

    if (total <= 0) {
        stop(
            "Kraken2 host report has zero total reads/read pairs: ",
            report_file
        )
    }

    list(
        classified = classified_root,
        unclassified = unclassified,
        total = total
    )
}

read_host_totals <- function(summary_path, report_files) {
    if (!is.na(summary_path) && file.exists(summary_path)) {
        summary <- read.table(
            summary_path,
            sep = "	",
            header = TRUE,
            quote = "",
            check.names = FALSE,
            stringsAsFactors = FALSE
        )

        required <- c(
            "sample",
            "host_classified_reads",
            "host_unclassified_reads"
        )

        missing <- setdiff(required, names(summary))

        if (length(missing) > 0) {
            stop(
                "Host Kraken summary is missing columns: ",
                paste(missing, collapse = ", ")
            )
        }

        totals <- suppressWarnings(
            as.numeric(summary$host_classified_reads) +
            as.numeric(summary$host_unclassified_reads)
        )

        names(totals) <- summary$sample

        if (any(is.na(totals)) || any(totals <= 0)) {
            stop(
                "Invalid total host read/read-pair counts in: ",
                summary_path
            )
        }

        return(totals)
    }

    if (length(report_files) == 0) {
        stop(
            "No host Kraken summary or Report_host_*.txt files were found."
        )
    }

    totals <- c()

    for (report_file in report_files) {
        counts <- read_host_report_counts(report_file)

        sample <- gsub(
            "^Report_host_|\\.txt$",
            "",
            basename(report_file)
        )

        totals <- c(
            totals,
            setNames(counts$total, sample)
        )
    }

    totals
}

read_host_ratios <- function(summary_path, report_files) {
    if (!is.na(summary_path) && file.exists(summary_path)) {
        summary <- read.table(
            summary_path,
            sep = "	",
            header = TRUE,
            quote = "",
            check.names = FALSE,
            stringsAsFactors = FALSE
        )

        required <- c(
            "sample",
            "host_classified_reads",
            "host_unclassified_reads"
        )

        missing <- setdiff(required, names(summary))

        if (length(missing) > 0) {
            stop(
                "Host Kraken summary is missing columns: ",
                paste(missing, collapse = ", ")
            )
        }

        classified <- suppressWarnings(
            as.numeric(summary$host_classified_reads)
        )

        unclassified <- suppressWarnings(
            as.numeric(summary$host_unclassified_reads)
        )

        total <- classified + unclassified

        unclassified_pct <- 100 * unclassified / total

        non_host_host_ratio <- ifelse(
            classified > 0,
            unclassified / classified,
            NA_real_
        )

        names(unclassified_pct) <- summary$sample
        names(non_host_host_ratio) <- summary$sample

        return(
            list(
                unclassified_reads_ratio_percent = unclassified_pct,
                non_host_host_reads_ratio = non_host_host_ratio
            )
        )
    }

    unclassified_pct <- c()
    non_host_host_ratio <- c()

    for (report_file in report_files) {
        counts <- read_host_report_counts(report_file)

        sample <- gsub(
            "^Report_host_|\\.txt$",
            "",
            basename(report_file)
        )

        pct <- 100 * counts$unclassified / counts$total

        ratio <- ifelse(
            counts$classified > 0,
            counts$unclassified / counts$classified,
            NA_real_
        )

        unclassified_pct <- c(
            unclassified_pct,
            setNames(pct, sample)
        )

        non_host_host_ratio <- c(
            non_host_host_ratio,
            setNames(ratio, sample)
        )
    }

    list(
        unclassified_reads_ratio_percent = unclassified_pct,
        non_host_host_reads_ratio = non_host_host_ratio
    )
}


# make folder DEG to store outputs
setwd(dirname(args[1]))
#system("mkdir -p DEG") # make new a folder DEG in current working directory

library(DESeq2)
library(tibble)
# Read & preprocess the input file cts before run deseq2 analysis
filename<-basename(args[1])
if (filename %in% c("humann_genefamilies_Abundance_go_translated.tsv","humann_genefamilies_Abundance_kegg_translated.tsv")){
  # read the data file
  #translated file has additional quote symbol "\""
  cts<-read.table(args[1],
                  row.names=1, sep="\t",header=T)
  #round the values before as.integer
  cts<-as.data.frame(lapply(lapply(cts,round),as.integer),row.names = row.names(cts))
  # make a folder for outputs
  if (filename == "humann_genefamilies_Abundance_go_translated.tsv"){
    dir.create("Nonhost_hmn_DEG/GO",recursive = T)
    setwd("Nonhost_hmn_DEG/GO")
  }
  if (filename == "humann_genefamilies_Abundance_kegg_translated.tsv"){
    dir.create("Nonhost_hmn_DEG/KEGG",recursive = T)
    setwd("Nonhost_hmn_DEG/KEGG")
  }
}
if (filename %in% c("bracken_species_all","bracken_phylum_all","bracken_genus_all")){
  # read the data file
  #bracken file (eg. bracken_species_all) without quote symbol "\""; mark empty quote as quote=""
  cts<-read.table(args[1],
                  row.names=1, sep="\t",header=T, quote="")
  # Extract columns with count numbers
  cts<-cts[,grepl("*_num", names(cts))]
  # match column name to sample name


#cts <- cts[, grepl("*_num", names(cts))]  # Filtra as colunas de cts
features <- cts  # Aqui, features pode ser definido como cts ou um subconjunto específico de cts




print(dim(cts))
print(dim(features))
  colnames(cts)<-gsub("^Report_|\\.species.bracken_num$|\\.phylum.bracken_num$|.genus.bracken_num$","",colnames(cts))
  # decontamination step
  # read a list of contamination organisms
  # conta_ls<-read.table("~/Dual-seq/conta_ls.txt",sep="\t",colClasses = c("character","NULL"))
  # conta_ls<-conta_ls[conta_ls != "",]
  # remove the entries matching the contamination organisms
  # for (c in conta_ls){
  #   cts<-cts[!grepl(c,row.names(cts)),]
  # }
  # make a folder for outputs
  dir.create("Nonhost_DEG", recursive = T)
  setwd("Nonhost_DEG")
  # save the decontaminated result for later reference
  # write.table(cts,paste0(filename,"_decontaminated"),sep="\t", quote = F, col.names = NA)
}
if (filename == "host_counts.txt"){
  # read the data file
  #featureCounts file (eg. host_counts.txt) without quote symbol "\""; mark empty quote as quote=""
  cts<-read.table(args[1],
                  row.names=1, sep="\t",header=T, quote="")
  # drop first 5 columns with information other than counts
  cts<-cts[,-c(1:5)]
  # drop rows of zero count
  cts<-cts[rowSums(cts[-1])>0,]
  # make a folder for outputs
  dir.create("Host_DEG", recursive = T)
  setwd("Host_DEG")
}

# read the original samplesheet
coldata0 <- read.csv(args[2], header = T, na.strings=c("","NA"))
# extract the contrast information for reference
coldata_vs <- coldata0[c("group1", "group2")]
coldata_vs <- coldata_vs[rowSums(is.na(coldata_vs)) == 0, ]
coldata_vs <- unique(coldata_vs)
# read samplesheet as factors (as.is = F) for Deseq2 statistical analysis

coldata_factor <- read.csv(args[2], header = T, as.is = F)
coldata_factor[]<-lapply(coldata_factor, factor)
coldata<-coldata_factor[,1:2]
# update the coldata if metadata is provided
if (length(args) == 5){
  coldata <- read.csv(args[5], header = T, as.is = F)
  coldata[]<-lapply(coldata, factor)
}

if (filename == "host_counts.txt"){
  # make cts(count matrix) has consistent order with samplesheet/metadata
  cts<-cts[coldata0$sample_name]
  # load the datastructure to DESeq
  dds <- DESeqDataSetFromMatrix(countData = cts,
                                colData = coldata,
                                design= ~ group)
}

if (filename %in% c("bracken_species_all","bracken_phylum_all","bracken_genus_all",
                    "humann_genefamilies_Abundance_go_translated.tsv","humann_genefamilies_Abundance_kegg_translated.tsv")){
  if (filename %in% c("bracken_species_all","bracken_phylum_all","bracken_genus_all")){
    files_h <- list.files(path=paste0(dirname(args[1]),"/temp"), pattern="^Report_host_.*\\.txt$", full.names=TRUE, recursive=FALSE)
  }
  if (filename %in% c("humann_genefamilies_Abundance_go_translated.tsv","humann_genefamilies_Abundance_kegg_translated.tsv")){
    humann_f <- gsub("/hmn_genefamily_abundance_files","",dirname(args[1]))
    files_h <- list.files(path=paste0(humann_f,"/temp"), pattern="^Report_host_.*\\.txt$", full.names=TRUE, recursive=FALSE)
  }
  # Robust host transcriptome size covariate.
  # Prefer the MTD_explorer.sh summary table, which is explicit and
  # works for both single-end reads and paired-end read pairs.
  # Fall back to Report_host_*.txt for backward compatibility.
  host_summary_path <- find_host_summary_path(
    dirname(args[1])
  )

  transcriptome_size <- read_host_totals(
    summary_path = host_summary_path,
    report_files = files_h
  )

  transcriptome_size <- log2(transcriptome_size)-mean(log2(transcriptome_size))
  coldata$order<-1:nrow(coldata)
  coldata<-merge(coldata,as.data.frame(transcriptome_size), by.x="sample_name",by.y="row.names")
  coldata<-coldata[order(coldata$order), ]
  coldata<-subset(coldata, select = -c(order))
  # make cts(count matrix) has consistent order with samplesheet/metadata
  cts<-cts[coldata0$sample_name]
  # load the datastructure to DESeq
  dds <- DESeqDataSetFromMatrix(countData = cts,
                                colData = coldata,
                                design= ~ group + transcriptome_size)
}

# adjust the design if metadata is provided
funNew <- function(x){
    as.formula(paste("~", paste(x, collapse = " + ")))
  }
if (length(args) == 5){
  design(dds)<-funNew(names(coldata)[2:ncol(coldata)])
}

# perform the DESeq analysis
dds <- DESeq(dds)  


# Data transformation for visualization (normalization included)
#rld<-rlog(dds,blind=F) # regularized log transformation (log2 based)
if (dim(results(dds))[1]  < 1000 || min(colSums(cts !=0)) < 1000){
  vsd<-varianceStabilizingTransformation(dds,blind=F) # vatiance stabilizing transformation
} else {
  vsd<-vst(dds,blind=F)
}
normtrans<-assay(vsd)

# save normalized & transformed data for visualization
write.csv(normtrans,file=paste0(sub(".tsv$|.txt$","",filename),"_normalized_transformed.csv"))

# save normalized (untransformed) data for reference
norm<-counts(dds,normalized=T)
write.csv(norm,file=paste0(sub(".tsv$|.txt$","",filename),"_normalized.csv"))

# merge and add suffixes; normalized and normalized&transformed
merge.nt<-merge(norm,normtrans,by="row.names", suffixes=c(".norm",".normtrans"))
if (filename == "host_counts.txt") {
  host_sp <- read.csv(args[4])

  names(merge.nt)[1] <- "GeneID"
  genes <- merge.nt$GeneID

  cache_file <- file.path(dirname(args[1]), "Host_DEG", "gene_ID_cache.csv")

  if (file.exists(cache_file)) {
    message("Using cached gene annotations: ", cache_file)

    gene_ID <- read.csv(cache_file, header = TRUE, check.names = FALSE)

    required_gene_id_cols <- c(
      "gene_name",
      "ensembl_gene_id",
      "chromosome_name",
      "start_position",
      "end_position",
      "strand",
      "gene_biotype",
      "description",
      "gene_length"
    )

    missing_gene_id_cols <- setdiff(required_gene_id_cols, names(gene_ID))
    if (length(missing_gene_id_cols) > 0) {
      stop("Cached gene_ID file is missing columns: ",
           paste(missing_gene_id_cols, collapse = ", "))
    }

    gene_ID <- gene_ID[gene_ID$ensembl_gene_id %in% genes, ]

  } else {
    message("No cached annotation found. Trying Ensembl online...")

    library("biomaRt")

    dataset_use <- as.character(host_sp[host_sp$Taxon_ID == args[3], 2])
    message("Using Ensembl dataset: ", dataset_use)

    ensembl <- NULL

    ensembl <- tryCatch({
      biomaRt::useEnsembl(
        biomart = "genes",
        dataset = dataset_use,
        mirror = "www"
      )
    }, error = function(e) {
      message("Failed to connect to Ensembl using useEnsembl(): ", conditionMessage(e))
      NULL
    })

    if (is.null(ensembl)) {
      stop("Failed to connect to Ensembl and no local cache was found: ", cache_file)
    }

    gene_ID <- biomaRt::getBM(
      filters = "ensembl_gene_id",
      attributes = c(
        "external_gene_name",
        "ensembl_gene_id",
        "chromosome_name",
        "start_position",
        "end_position",
        "strand",
        "gene_biotype",
        "description"
      ),
      values = genes,
      mart = ensembl
    )

    names(gene_ID)[names(gene_ID) == "external_gene_name"] <- "gene_name"

    gene_len <- read.table(args[1], row.names = 1, sep = "\t", header = TRUE, quote = "")
    gene_len <- gene_len["Length"]
    colnames(gene_len) <- "gene_length"

    gene_ID <- merge(gene_ID, gene_len, by.x = "ensembl_gene_id", by.y = "row.names")

    write.csv(gene_ID, cache_file, row.names = FALSE, quote = TRUE)
    message("Saved gene annotation cache: ", cache_file)
  }
}
# #to add entrezid separately; caution: may bring a issue of "duplicate" ENTREZID!!
# library(clusterProfiler)
# library(org.Mmu.eg.db)
# # may bring a issue of "duplicate" ENTREZID!!
# ENTREZID = bitr(gene_ID$external_gene_name,
#                 fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mmu.eg.db") #get entrezid with acceptable id duplicate
# gene_ID <- merge(gene_ID,ENTREZID,by.x="gene_name",by.y="SYMBOL",all.x = T)

# function for host
# function of contrast between groups, merge to annotation and count tables, and save to .csv files
comparison<-function(dds,coldata_vs,filename,gene_ID,cts,merge.nt){
  for (i in 1:nrow(coldata_vs)){
  group1<-coldata_vs$group1[i]
  group2<-coldata_vs$group2[i]
  comparison <- results(dds, contrast=c("group", group1, group2))
  comparison_f<-as.data.frame(comparison)
  comparison_f<-comparison_f[order(comparison_f$pvalue),]
  system(paste0("mkdir -p"," ",group1,"_vs_",group2)) #for output file structure
  setwd(paste0(group1,"_vs_",group2))
  write.csv(comparison_f,file=paste0(sub(".tsv$|.txt$","",filename),"_",group1,"_vs_",group2,".csv"))
  setwd("../")
  colnames(comparison_f)<-paste0(colnames(comparison_f),".",group1,"_vs_",group2)
  gene_ID<-merge(gene_ID,comparison_f,by.x="ensembl_gene_id", by.y="row.names") #merge each comparison to annotation
  }
  Anno.merge.all<-merge(cts,gene_ID,by.x="row.names",by.y="ensembl_gene_id")
  Anno.merge.all<-merge(Anno.merge.all,merge.nt,by.x="Row.names",by.y="GeneID")
  names(Anno.merge.all)[1]<-"gene_id" #rename the first column
  Anno.merge.all[Anno.merge.all==""]<-"-" #replace the empty cells to "-"
  Anno.merge.all<-Anno.merge.all[complete.cases(Anno.merge.all[,grep("pvalue|padj",names(Anno.merge.all))]),] # drop pvalue == NA rows; optional
  hybrid_name<-Anno.merge.all$gene_name # add a column of gene_id/name hybrid
  while ("-" %in% hybrid_name){
    hybrid_name[match("-",hybrid_name)] <-
      Anno.merge.all$gene_id[match("-",hybrid_name)]}
  Anno.merge.all<-add_column(Anno.merge.all, hybrid_name, .after = "gene_name")
  write.csv(Anno.merge.all,file=paste0(sub(".tsv$|.txt$","",filename),"_DEG.csv"),row.names=F)
}

#function for non-host
comparison_nonhost<-function(dds,coldata_vs,filename,cts,merge.nt){
  gene_ID_temp<-list()
  for (i in 1:nrow(coldata_vs)){
    group1<-coldata_vs$group1[i]
    group2<-coldata_vs$group2[i]
    comparison <- results(dds, contrast=c("group", group1, group2))
    comparison_f<-as.data.frame(comparison)
    comparison_f<-comparison_f[order(comparison_f$pvalue),]
    system(paste0("mkdir -p"," ",group1,"_vs_",group2)) #for output file structure
    setwd(paste0(group1,"_vs_",group2))
    write.csv(comparison_f,file=paste0(sub(".tsv$|.txt$","",filename),"_",group1,"_vs_",group2,".csv"))
    setwd("../")
    colnames(comparison_f)<-paste0(colnames(comparison_f),".",group1,"_vs_",group2)
    comparison_f <- tibble::rownames_to_column(comparison_f, "NAME")
    gene_ID_temp[[i]]<-comparison_f
  }
  library(dplyr)
  gene_ID <- purrr::reduce(gene_ID_temp,full_join, by = "NAME") #reduce multiple data.frames in a list to a single data.frame
  names(merge.nt)[1]<-"Name" #rename the first column
  Anno.merge.all<-merge(cts,gene_ID,by.x="row.names", by.y="NAME")
  Anno.merge.all<-merge(Anno.merge.all,merge.nt,by.x="Row.names",by.y="Name")
  names(Anno.merge.all)[1]<-"Name" #rename the first column
  Anno.merge.all[Anno.merge.all==""]<-"-" #replace the empty cells to "-"
  Anno.merge.all<-Anno.merge.all[complete.cases(Anno.merge.all[,grep("pvalue|padj",names(Anno.merge.all))]),] # drop pvalue == NA rows
  write.csv(Anno.merge.all,file=paste0(sub(".tsv$|.txt$","",filename),"_DEG.csv"),row.names=F)
}

# apply the comparison function
if (filename == "host_counts.txt"){
  comparison(dds, coldata_vs, filename,gene_ID,cts,merge.nt)
} else {comparison_nonhost(dds, coldata_vs, filename,cts,merge.nt)}

#to match with sample name from DESeq2, which is not allow "-" in the name
coldata[,1] <-gsub("-",".",coldata[,1])

## MaAsLin2 ##
library("Maaslin2")

dir.create(
  "MaAsLin2_results",
  showWarnings = FALSE,
  recursive = TRUE
)

features <- t(cts)

metadata <- coldata
rownames(metadata) <- as.character(metadata[, 1])
metadata <- metadata[, -1, drop = FALSE]

# ------------------------------------------------------------
# MaAsLin2 model construction
#
# Repeated-measures identifiers such as animal or subject must
# not be tested as ordinary fixed effects. They are represented
# as random intercepts.
# ------------------------------------------------------------

maaslin_subject_candidates <- c(
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

maaslin_subject_candidates <- intersect(
  maaslin_subject_candidates,
  names(metadata)
)

# A subject column is useful as a random effect only when
# identifiers occur in more than one sample.
maaslin_subject_candidates <- maaslin_subject_candidates[
  vapply(
    maaslin_subject_candidates,
    function(column) {
      values <- as.character(metadata[[column]])
      values <- values[
        !is.na(values) &
        values != ""
      ]

      length(unique(values)) > 1 &&
        anyDuplicated(values) > 0
    },
    logical(1)
  )
]

maaslin_random_effects <- NULL

if (length(maaslin_subject_candidates) > 0) {
  maaslin_random_effects <- maaslin_subject_candidates[1]
}

maaslin_fixed_effects <- setdiff(
  names(metadata),
  maaslin_random_effects
)

# Remove constant metadata columns because they cannot be fitted.
maaslin_fixed_effects <- maaslin_fixed_effects[
  vapply(
    maaslin_fixed_effects,
    function(column) {
      values <- metadata[[column]]

      length(
        unique(values[!is.na(values)])
      ) >= 2
    },
    logical(1)
  )
]

if (!"group" %in% maaslin_fixed_effects) {
  stop(
    "MaAsLin2 requires a non-constant 'group' ",
    "column in the metadata."
  )
}

# Build references for group and any additional categorical
# fixed effect. Subject identifiers are excluded because they
# are random effects.
build_maaslin_reference <- function(group_reference) {
  references <- paste0(
    "group,",
    group_reference
  )

  other_fixed_effects <- setdiff(
    maaslin_fixed_effects,
    "group"
  )

  for (effect in other_fixed_effects) {
    values <- metadata[[effect]]

    if (!is.numeric(values)) {
      effect_levels <- levels(
        factor(values)
      )

      if (length(effect_levels) >= 2) {
        references <- c(
          references,
          paste0(
            effect,
            ",",
            effect_levels[1]
          )
        )
      }
    }
  }

  paste(
    references,
    collapse = ";"
  )
}

run_mtd_maaslin2 <- function(
  input_features,
  output_dir,
  group_reference
) {
  dir.create(
    output_dir,
    showWarnings = FALSE,
    recursive = TRUE
  )

  reference_value <- build_maaslin_reference(
    group_reference
  )

  message(
    "[MaAsLin2] Fixed effects: ",
    paste(
      maaslin_fixed_effects,
      collapse = ", "
    )
  )

  message(
    "[MaAsLin2] Random effects: ",
    if (is.null(maaslin_random_effects)) {
      "none"
    } else {
      paste(
        maaslin_random_effects,
        collapse = ", "
      )
    }
  )

  message(
    "[MaAsLin2] Reference: ",
    reference_value
  )

  tryCatch(
    Maaslin2(
      input_features,
      metadata,
      output_dir,
      fixed_effects = maaslin_fixed_effects,
      random_effects = maaslin_random_effects,
      reference = reference_value,
      plot_heatmap = TRUE,
      plot_scatter = TRUE,
      min_abundance = 0.01,
      min_prevalence = 0.1,
      cores = 10
    ),
    error = function(e) {
      error_message <- paste0(
        "MaAsLin2 failed, but the remaining DEG analysis ",
        "will continue.\n",
        "Output: ",
        output_dir,
        "\nReason: ",
        conditionMessage(e),
        "\n"
      )

      warning(
        error_message,
        call. = FALSE
      )

      writeLines(
        error_message,
        file.path(
          output_dir,
          "MTD_MaAsLin2_ERROR.txt"
        )
      )

      NULL
    }
  )
}

group_references <- unique(
  as.character(coldata_vs[, 2])
)

group_references <- group_references[
  !is.na(group_references) &
  group_references != ""
]

for (group_reference in group_references) {
  fit_data <- run_mtd_maaslin2(
    input_features = features,
    output_dir = file.path(
      "MaAsLin2_results",
      paste0(
        "ref_",
        group_reference
      )
    ),
    group_reference = group_reference
  )
}

if (filename == "host_counts.txt") {
  # Repeat MaAsLin2 using non-duplicated gene symbols.
  dir.create(
    "MaAsLin2_results/gene_symbol",
    showWarnings = FALSE,
    recursive = TRUE
  )

  DEG <- read.csv(
    paste0(
      sub(
        ".tsv$|.txt$",
        "",
        filename
      ),
      "_DEG.csv"
    ),
    header = TRUE
  )

  df.DEG <- DEG[
    !duplicated(DEG$hybrid_name),
  ]

  df4features <- data.frame(
    df.DEG$hybrid_name,
    df.DEG[names(cts)],
    row.names = 1
  )

  features.symbol <- t(df4features)

  for (group_reference in group_references) {
    fit_data <- run_mtd_maaslin2(
      input_features = features.symbol,
      output_dir = file.path(
        "MaAsLin2_results",
        "gene_symbol",
        paste0(
          "ref_",
          group_reference
        )
      ),
      group_reference = group_reference
    )
  }

  write(
    paste(
      "The number of genes may be lower because",
      "duplicated gene symbols were removed."
    ),
    "MaAsLin2_results/gene_symbol/Readme.txt"
  )
}


if (filename == "host_counts.txt"){
  setwd(paste0(dirname(args[1]),"/Host_DEG")) # go back to the Host_DEG folder
}

### Plots ###
library(ggplot2)
library(ggrepel)
library(dplyr)
library(stringr) #for str_trunc function: restrict the showing character length
library(colorspace)
library(RColorBrewer)
# subset normtrans for visualization; consider all groups; to integrate in the comparison function?? DEG=Anno.merge.all
DEG<-read.csv(paste0(sub(".tsv$|.txt$","",filename),"_DEG.csv"),header = T)

# count gene TPM
if (filename == "host_counts.txt"){
  RPK <- DEG[,as.character(coldata[,1])]/DEG$gene_length
  TPM <- RPK*1000000/sum(RPK)
  TPM<-cbind(DEG[,c("gene_name","hybrid_name")],TPM)
  write.csv(TPM,file=paste0(sub(".tsv$|.txt$","",filename),"_TPM.csv"),row.names=F)
}

#### BEGIN BLOCK: bounded DEG heatmap selection ####
#
# Select genes using adjusted p-values whenever available and cap the
# number of displayed rows. This avoids multi-hundred-inch PDF files.
#
# The limit can be changed without editing the script:
#   export MTD_HEATMAP_MAX_GENES=200

heatmap_max_genes <- suppressWarnings(
  as.integer(
    Sys.getenv(
      "MTD_HEATMAP_MAX_GENES",
      "150"
    )
  )
)

if (
  is.na(heatmap_max_genes) ||
  heatmap_max_genes < 10
) {
  heatmap_max_genes <- 150L
}

heatmap_significant_rows <- list()
heatmap_ranked_rows <- list()

for (comparison_index in seq_len(nrow(coldata_vs))) {

  group1 <- as.character(
    coldata_vs$group1[comparison_index]
  )

  group2 <- as.character(
    coldata_vs$group2[comparison_index]
  )

  comparison_name <- paste0(
    group1,
    "_vs_",
    group2
  )

  adjusted_p_column <- paste0(
    "padj.",
    comparison_name
  )

  raw_p_column <- paste0(
    "pvalue.",
    comparison_name
  )

  fold_change_column <- paste0(
    "log2FoldChange.",
    comparison_name
  )

  if (!fold_change_column %in% names(DEG)) {
    next
  }

  if (adjusted_p_column %in% names(DEG)) {
    p_column <- adjusted_p_column
  } else if (raw_p_column %in% names(DEG)) {
    p_column <- raw_p_column
  } else {
    next
  }

  p_values <- suppressWarnings(
    as.numeric(DEG[[p_column]])
  )

  fold_changes <- suppressWarnings(
    as.numeric(DEG[[fold_change_column]])
  )

  complete_rows <- (
    !is.na(p_values) &
      is.finite(p_values) &
      !is.na(fold_changes) &
      is.finite(fold_changes)
  )

  ranked_rows <- DEG[
    complete_rows,
    ,
    drop = FALSE
  ]

  ranked_rows$.heatmap_p <- p_values[
    complete_rows
  ]

  ranked_rows$.heatmap_abs_fc <- abs(
    fold_changes[complete_rows]
  )

  heatmap_ranked_rows[[
    length(heatmap_ranked_rows) + 1L
  ]] <- ranked_rows

  significant_rows <- (
    complete_rows &
      p_values < 0.05 &
      abs(fold_changes) > 0.5
  )

  if (any(significant_rows)) {

    selected_rows <- DEG[
      significant_rows,
      ,
      drop = FALSE
    ]

    selected_rows$.heatmap_p <- p_values[
      significant_rows
    ]

    selected_rows$.heatmap_abs_fc <- abs(
      fold_changes[significant_rows]
    )

    heatmap_significant_rows[[
      length(heatmap_significant_rows) + 1L
    ]] <- selected_rows
  }
}

if (length(heatmap_significant_rows) > 0) {

  flt_groups_all <- do.call(
    rbind,
    heatmap_significant_rows
  )

  heatmap_selection_mode <- (
    "adjusted-p/significance-filtered"
  )

} else if (length(heatmap_ranked_rows) > 0) {

  flt_groups_all <- do.call(
    rbind,
    heatmap_ranked_rows
  )

  heatmap_selection_mode <- (
    "top-ranked fallback; no genes passed the heatmap cutoff"
  )

} else {

  stop(
    "No complete differential-expression rows were available ",
    "for the heatmap."
  )
}

flt_groups_all <- flt_groups_all[
  order(
    flt_groups_all$.heatmap_p,
    -flt_groups_all$.heatmap_abs_fc,
    na.last = TRUE
  ),
  ,
  drop = FALSE
]

heatmap_id_column <- if (
  "gene_id" %in% names(flt_groups_all)
) {
  "gene_id"
} else if (
  "Name" %in% names(flt_groups_all)
) {
  "Name"
} else if (
  "hybrid_name" %in% names(flt_groups_all)
) {
  "hybrid_name"
} else {
  NULL
}

if (!is.null(heatmap_id_column)) {
  flt_groups_all <- flt_groups_all[
    !duplicated(
      as.character(
        flt_groups_all[[heatmap_id_column]]
      )
    ),
    ,
    drop = FALSE
  ]
} else {
  flt_groups_all <- unique(
    flt_groups_all
  )
}

heatmap_rows_before_limit <- nrow(
  flt_groups_all
)

flt_groups_all <- head(
  flt_groups_all,
  heatmap_max_genes
)

flt_groups_all$.heatmap_p <- NULL
flt_groups_all$.heatmap_abs_fc <- NULL

message(
  "[HEATMAP] Selection mode: ",
  heatmap_selection_mode
)

message(
  "[HEATMAP] Genes available before limit: ",
  heatmap_rows_before_limit
)

message(
  "[HEATMAP] Genes plotted: ",
  nrow(flt_groups_all),
  " (maximum ",
  heatmap_max_genes,
  ")"
)
#### END BLOCK: bounded DEG heatmap selection ####

#subset the interested columnes
subset_normtrans<-flt_groups_all[,grep("normtrans|gene_name|gene_id|Name",names(flt_groups_all))]
if (filename == "host_counts.txt"){
  # pass undefined gene name
  while ("-" %in% subset_normtrans$gene_name){
    subset_normtrans$gene_name[match("-",subset_normtrans$gene_name)] <-
      subset_normtrans$gene_id[match("-",subset_normtrans$gene_name)]
  }
  ##check the duplicated rows
  #n_occur <- data.frame(table(subset_normtrans$gene_name))
  #n_occur[n_occur$Freq > 1,]
  ##drop the duplicated rows
  #flt_groups_all<-unique(flt_groups_all)
  # drop the duplicate rows
  subset_normtrans<-subset_normtrans[!duplicated(subset_normtrans$gene_name),]
  #convert column name (gene_name) to row name
  subset_normtrans2<-subset_normtrans[,-1:-2]
  rownames(subset_normtrans2)<-subset_normtrans[,2]
} else {
  # drop the duplicate rows
  subset_normtrans<-subset_normtrans[!duplicated(subset_normtrans$Name),]
  #convert column name (gene_name) to row name
  subset_normtrans2<-subset_normtrans[,-1]
  rownames(subset_normtrans2)<-subset_normtrans[,1]
}
#chop the sample column names
colnames(subset_normtrans2)<-gsub(".normtrans","",colnames(subset_normtrans2))
#heatmap require matrix
transdata <- as.matrix(subset_normtrans2)
# adjust the rowname length
for (i in 1:length(row.names(transdata))){
  l<-nchar(row.names(transdata)[i])
  if (l > 35){
    row.names(transdata)[i]<-str_trunc(row.names(transdata)[i],35)
  }
}
#read the column annotation information (derive from samplesheet.csv)
Anno_col<-as.data.frame(coldata[,2])
rownames(Anno_col)<-coldata[,1]
colnames(Anno_col)<-"Group"
# annotation color selection
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
n<-nrow(unique(Anno_col))
cols = gg_color_hue(n)
names(cols)<-as.character(unique(Anno_col)$Group)
anno_colors<-list(Group=cols)

# draw heat map
library(pheatmap)
hp_thumbnail<-pheatmap(transdata, cluster_cols = T, scale="row",
          annotation_col =Anno_col, annotation_colors=anno_colors,
          show_rownames=F,cex=1)
if (nrow(transdata)>100){
  hp<-pheatmap(transdata, cluster_cols = T, scale="row",
                         fontsize_row=5, annotation_col =Anno_col,
               annotation_colors=anno_colors)
} else {
  hp<-pheatmap(transdata, cluster_cols = T, scale="row",
               annotation_col =Anno_col,
               annotation_colors=anno_colors)
}
# hp<-pheatmap(transdata, cluster_cols = T, scale="row",
#              fontsize_row=5,cex=1, annotation_col =Anno_col, cex=0.9)

# save the plot
ggsave("heatmap_thumbnail.pdf",plot=hp_thumbnail)
if (nrow(transdata)>100){
  ggsave("heatmap.pdf",plot=hp,
                limitsize = F,
                height = 0.1*nrow(transdata),
                width = 0.8*ncol(transdata))
} else {
  ggsave("heatmap.pdf",plot=hp,height = 0.2*nrow(transdata))
}
## draw PCA
if (filename %in%
    c("humann_genefamilies_Abundance_go_translated.tsv",
      "humann_genefamilies_Abundance_kegg_translated.tsv",
      "host_counts.txt")){
  pcadata<-plotPCA(vsd,intgroup=c("group"), returnData=T)
  percentVar<-round(100*attr(pcadata,"percentVar"))

  ggplot(pcadata, aes(PC1, PC2, shape=group)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    coord_fixed() +
    theme_bw() +
    ggtitle("PCA")
  ggsave("PCA.pdf")

  ggplot(pcadata, aes(PC1, PC2, color=group)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    coord_fixed() +
    theme_bw() +
    ggtitle("PCA")
  ggsave("PCA_color.pdf")

  ggplot(pcadata, aes(PC1, PC2, shape=group)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    geom_text_repel(aes(label=pcadata$name), size=3, max.overlaps = Inf, segment.size = 0.5) +
    coord_fixed() +
    theme_bw() +
    ggtitle("PCA")
  ggsave("PCA_label.pdf")

  ggplot(pcadata, aes(PC1, PC2, color=group)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    geom_text_repel(aes(label=pcadata$name),size=3, max.overlaps = Inf) +
    coord_fixed() +
    theme_bw() +
    ggtitle("PCA")
  ggsave("PCA_label_color.pdf")
}
# draw PCoA, heatmap, diversity tests for microbiome
if (filename %in% c("bracken_species_all",
                    "bracken_phylum_all",
                    "bracken_genus_all")){
  library(phyloseq)
  library(vegan)
  #Imported original biom file into phyloseq object
  biomfilename=paste0(dirname(args[1]),"/temp/bracken_species_all0.biom")
  data<-import_biom(biomfilename,parseFunction = parse_taxonomy_default)
  colnames(tax_table(data)) <- c("Kingdom", "Phylum", "Class", "Order", "Family",  "Genus", "Species")
  #Imported host transcriptome size adjusted & normtrans biom file into phyloseq object
  biomfilename.adj=paste0(dirname(args[1]),"/bracken_species_all.biom")
  data.adj<-import_biom(biomfilename.adj,parseFunction = parse_taxonomy_default)
  colnames(tax_table(data.adj)) <- c("Kingdom", "Phylum", "Class", "Order", "Family",  "Genus", "Species")
  #Estimated and exported alpha diversity
  #The measures of diversity that aren't totally reliant on singletons, eg. Shannon/Simpson, are valid to use, and users can ignore the warning in phyloseq when calculating those measures.
data.alpha<-estimate_richness(data.adj, measures = c("Observed", "Shannon", "Simpson", "InvSimpson"))
  write.csv(data.alpha, file="alpha-diversity.csv")
  theme_set(theme_bw())
  p.alpha <- plot_richness(data,measures = c("Shannon","Simpson"))
  p.alpha
  ggsave("Alpha_diversity_sample.pdf")
  #Estimated and exported beta diversity (Bray-Curtis)
  #to transformation data by vst before beta diversity analysis
  if (length(args) == 5){
    samples4phyloseq<-as.data.frame(coldata[2:ncol(coldata)], row.names = as.character(coldata$sample_name))
    sampledata<-sample_data(samples4phyloseq)
    data<-merge_phyloseq(data,sampledata)
    dds1 <- phyloseq_to_deseq2(data, funNew(names(coldata)[2:ncol(coldata)]))
  } else {
    samples4phyloseq<-as.data.frame(coldata$group, row.names = as.character(coldata$sample_name))
    colnames(samples4phyloseq)<-"group"
    sampledata<-sample_data(samples4phyloseq)
    data<-merge_phyloseq(data,sampledata)
    dds1 <- phyloseq_to_deseq2(data, ~ group)
  }
  ## ANCOMBC ##
  library("ANCOMBC")
  pseq <- phyloseq::tax_glom(data, taxrank = "Species") # species taxid
  pseq1 <- microbiome::aggregate_taxa(data,"Species") # species name
  
  ancombc_out <- function(pseq,formula){
    ancombc(
      phyloseq = pseq, 
      formula = formula, 
      p_adj_method = "fdr", 
      zero_cut = 0.90, # by default prevalence filter of 10% is applied
      lib_cut = 0, 
      group = "group", 
      struc_zero = TRUE, 
      neg_lb = TRUE, 
      tol = 1e-5, 
      max_iter = 100, 
      conserve = TRUE, 
      alpha = 0.05, 
      global = TRUE
    )
  } 
  
  if (length(args) == 5){
    out <- ancombc_out(pseq,paste(names(coldata)[2:(ncol(coldata)-1)], collapse = " + "))
    out1 <- ancombc_out(pseq1,paste(names(coldata)[2:(ncol(coldata)-1)], collapse = " + "))
  } else {
    out <- ancombc_out(pseq,"group") # taxid
    out1 <- ancombc_out(pseq1,"group") # species name
  }

  res <- out$res # taxid
  res_global = out$res_global
  res1 <- out1$res # species name
  res_global1 = out1$res_global
  
  system("mkdir -p ANCOMBC_results")
  write.csv(res[["diff_abn"]], file="ANCOMBC_results/diff_abundance.csv")
  write.csv(res[["W"]], file="ANCOMBC_results/Test_statistics.csv")
  write.csv(res[["p_val"]], file="ANCOMBC_results/p_value.csv")
  write.csv(res[["q_val"]], file="ANCOMBC_results/q_value.csv")
  write.csv(res_global, file="ANCOMBC_results/Global_test.csv")
  
  system("mkdir -p ANCOMBC_results/with_species_names")
  write.csv(res1[["diff_abn"]], file="ANCOMBC_results/with_species_names/diff_abundance_name.csv")
  write.csv(res1[["W"]], file="ANCOMBC_results/with_species_names/Test_statistics_name.csv")
  write.csv(res1[["p_val"]], file="ANCOMBC_results/with_species_names/p_value_name.csv")
  write.csv(res1[["q_val"]], file="ANCOMBC_results/with_species_names/q_value_name.csv")
  write.csv(res_global1, file="ANCOMBC_results/with_species_names/Global_test_name.csv")
  
  ## heatmap for ANCOMBC results ##
  ANCOMBC_plot <- function(res){
    library(tidyr)
    heat.mw <- res[["W"]] %>% 
      rownames_to_column("taxid") %>% #preserve rownames
      gather(key, value, -taxid)
    heat.md <- res[["diff_abn"]] %>% 
      rownames_to_column("taxid") %>% #preserve rownames
      gather(key, value, -taxid)
    heat.m<-merge(heat.mw,heat.md,by="taxid")
    
    # Arrange the figure
    p <- ggplot(heat.m, aes(x = key.x, y = taxid, fill = value.x))
    p <- p + geom_tile() 
    p <- p + scale_fill_gradientn("value.x", name ="Test statistics",
                                  breaks = seq(from = -2, to = 2, by = 0.5), 
                                  colours = c("darkblue", "blue", "white", "red", "darkred"), 
                                  limits = c(-2,2)) 
    
    # Polish texts
    p <- p + theme(axis.text.x=element_text(angle = 90, hjust=1, face = "italic"),
                   axis.text.y=element_text(size = 8))
    p <- p + xlab("") + ylab("Taxonomy ID")
    
    # Mark the most significant cells with stars
    if (length(unique(heat.m$taxid))>100){
    p <- p + geom_text(data = subset(heat.m, value.y == "TRUE"),
                         aes(x = key.x, y = taxid, label = "+"))
    } else {
      p <- p + geom_text(data = subset(heat.m, value.y == "TRUE"), 
                         aes(x = key.x, y = taxid, label = "+"), col = "white", size = 3)
    }
    
    if (length(unique(heat.m$taxid))>100){
      ggsave("heatmap_ANCOMBC.pdf",plot=p,
             limitsize = F,
             height = 0.1*length(unique(heat.m$taxid)),
             width = 1.5+0.6*length(unique(heat.m$key.x)))
    } else {
      ggsave("heatmap_ANCOMBC.pdf",plot=p,height = 0.2*length(unique(heat.m$taxid)))
    }
  }
  setwd("ANCOMBC_results")
  ANCOMBC_plot(res)
  
  ANCOMBC_plot_sp <- function(res1){
    library(tidyr)
    heat.mw <- res1[["W"]] %>% 
      rownames_to_column("name") %>% #preserve rownames
      gather(key, value, -name)
    heat.md <- res1[["diff_abn"]] %>% 
      rownames_to_column("name") %>% #preserve rownames
      gather(key, value, -name)
    heat.m<-merge(heat.mw,heat.md,by="name")
    heat.m$name <- sub(".*s__", "", heat.m$name)
    heat.m$name <- stringr::str_trunc(heat.m$name, 31) # truncate sp_name
    # Arrange the figure
    p <- ggplot(heat.m, aes(x = key.x, y = name, fill = value.x))
    p <- p + geom_tile() 
    p <- p + scale_fill_gradientn("value.x", name ="Test statistics",
                                  breaks = seq(from = -2, to = 2, by = 0.5), 
                                  colours = c("darkblue", "blue", "white", "red", "darkred"), 
                                  limits = c(-2,2)) 
    
    # Polish texts
    p <- p + theme(axis.text.x=element_text(angle = 90, hjust=1, face = "italic"),
                   axis.text.y=element_text(size = 8))
    p <- p + xlab("") + ylab("Species name")
    
    # Mark the most significant cells with stars
    if (length(unique(heat.m$name))>100){
      p <- p + geom_text(data = subset(heat.m, value.y == "TRUE"), 
                         aes(x = key.x, y = name, label = "+"), col = "white", size = 2.5)
    } else {
      p <- p + geom_text(data = subset(heat.m, value.y == "TRUE"), 
                         aes(x = key.x, y = name, label = "+"), col = "white", size = 3)
    }
    
    if (length(unique(heat.m$name))>100){
      ggsave("heatmap_ANCOMBC.pdf",plot=p,
             limitsize = F,
             height = 0.1*length(unique(heat.m$name)),
             width = 2.5+0.8*length(unique(heat.m$key.x)))
    } else {
      ggsave("heatmap_ANCOMBC.pdf",plot=p,height = 0.2*length(unique(heat.m$name)))
    }
  }
  
  setwd("with_species_names")
  ANCOMBC_plot_sp(res1)
  setwd("../..")
  
  # continue for beta diversity
  data1<-data
  vsd1 <- varianceStabilizingTransformation(dds1)
  # normalized reads count with host transcriptome size and with avoiding removing variation associated with the other conditions
  if (length(args) == 5){
    mm <- model.matrix(funNew(names(coldata)[2:(ncol(coldata)-1)]), colData(vsd))
  } else {
    mm <- model.matrix(funNew(names(coldata)[2]), colData(vsd))
  }
  vsd1.df <- limma::removeBatchEffect(assay(vsd1), vsd$transcriptome_size, design=mm)
  vsd1.df[vsd1.df < 0.0] <- 0.0 #adjust negative values after vst
  otu_table(data1) <- otu_table(vsd1.df, taxa_are_rows = TRUE)

  braycurtis <- phyloseq::distance(data1, method = "bray")
  BCmat <- as.matrix(braycurtis)
  write.csv(BCmat, file = "braycurtis.csv")
  #Created and exported PCoA values
  braycurtis.pcoa <- ordinate(data1, method = "PCoA", distance = "bray")
  braycurtis.pcoa.export <- as.data.frame(braycurtis.pcoa$vectors, row.names = NULL, optional = FALSE, cut.names = FALSE, col.names = names(braycurtis.pcoa$vectors), fix.empty.names = TRUE, stringsAsFactors = default.stringsAsFactors())
  write.csv(braycurtis.pcoa.export, file="braycurtis-pcoa.csv")

  # add the group information in coldata
  coldata_order<-coldata # to keep the original order as the samplesheet
  coldata_order$order<-1:nrow(coldata_order)
  braycurtis.pcoa.export<-merge(braycurtis.pcoa.export, coldata_order,by.x="row.names",by.y = "sample_name")
  row.names(braycurtis.pcoa.export)<-braycurtis.pcoa.export[,1]
  braycurtis.pcoa.export<-braycurtis.pcoa.export[,-1]

  ## plot PCoA
  b<-braycurtis.pcoa[["values"]][["Relative_eig"]]
  #barplot(b[b>0],names.arg =colnames(braycurtis.pcoa.export))

  library(forcats)
  if (names(braycurtis.pcoa.export)[2]=="Axis.2"){
    # reorder the group column following the value of order column  
    braycurtis.pcoa.export %>%
      mutate(group = fct_reorder(group, order)) %>%
      ggplot(aes(Axis.1, Axis.2,color=group)) +
      geom_point(size=3) +
      xlab(paste0("PCoA1: ",round(100*b[1]),"% variance")) +
      ylab(paste0("PCoA2: ",round(100*b[2]),"% variance")) +
      geom_text_repel(aes(label=row.names(braycurtis.pcoa.export)), size=3, max.overlaps = Inf, segment.size = 0.5) +
      coord_fixed() +
      theme_bw() +
      ggtitle("Bray-Curtis Distances PCoA")
    ggsave("PCoA_label_color.pdf")
    
    braycurtis.pcoa.export %>%
      mutate(group = fct_reorder(group, order)) %>%
      ggplot(aes(Axis.1, Axis.2,color=group)) +
      geom_point(size=3) +
      xlab(paste0("PCoA1: ",round(100*b[1]),"% variance")) +
      ylab(paste0("PCoA2: ",round(100*b[2]),"% variance")) +
      coord_fixed() +
      theme_bw() +
      ggtitle("Bray-Curtis Distances PCoA")
    ggsave("PCoA_color.pdf")
  } else {
    write("No Axis.2 was found on PCoA","No_Axis2_on_PCoA.txt")
  }
  # anosim test
  data2<-t(as.matrix(data1@otu_table))
  pathotype.anosim <- anosim(data2, braycurtis.pcoa.export$group)
  # plot results
  pdf("ANOSIM.pdf")
  plot(pathotype.anosim,
       main="Analysis of Diversity in Groups",
       xlab="",
       ylab="")
  dev.off()
  # Start writing to an output file
  sink('ANOSIM-analysis-output.txt')
  summary(pathotype.anosim)
  # Stop writing to the file
  sink()

  # plot heatmap (Phyloseq style)
  sampledata<-as.data.frame(coldata[,2])
  row.names(sampledata)<-coldata[,1]
  colnames(sampledata)<-"Groups"
  sam  = sample_data(sampledata)
  sp = otu_table(norm, taxa_are_rows = TRUE)
  physeq<-phyloseq(sp, sam)
  if (nrow(norm)>100){
    ph<-plot_heatmap(physeq,max.label = nrow(norm)) +
      theme (axis.text.y = element_text(size=(2.2*225/nrow(norm))))
  } else {
    ph<-plot_heatmap(physeq)
  }
  ph$labels$y<-"Species"
  print(ph)
  ggsave("Heatmap_all.png")
  ggsave("Heatmap_all.pdf")

  # Add sample data
  tax  = tax_table(data.adj)
  otu  = otu_table(data.adj)
  data_sam <- phyloseq(otu,tax, sam)

  # Dynamic color scale for plot_bar
  # This avoids "Insufficient values in manual scale" when the number of phyla
  # is larger than a fixed palette.

  fix_tax_rank <- function(ps, rank = "Phylum", unknown = "Unclassified") {
    tt <- as(tax_table(ps), "matrix")

    if (!rank %in% colnames(tt)) {
      warning("[WARNING] Taxonomic rank not found in tax_table: ", rank)
      return(ps)
    }

    tt[, rank] <- as.character(tt[, rank])
    tt[is.na(tt[, rank]) | tt[, rank] == "" | tt[, rank] == "NA", rank] <- unknown

    tax_table(ps) <- tax_table(tt)
    ps
  }

  make_tax_palette <- function(ps, rank = "Phylum") {
    tt <- as(tax_table(ps), "matrix")

    if (!rank %in% colnames(tt)) {
      return(character(0))
    }

    taxa <- sort(unique(as.character(tt[, rank])))
    taxa <- taxa[!is.na(taxa) & taxa != "" & taxa != "NA"]

    n <- length(taxa)

    if (n == 0) {
      return(character(0))
    }

    pal <- grDevices::hcl.colors(n, palette = "Dark 3")
    names(pal) <- taxa

    pal
  }

  # Make sure missing phylum labels are explicit
  data <- fix_tax_rank(data, rank = "Phylum")
  data_sam <- fix_tax_rank(data_sam, rank = "Phylum")

  PhylaPalette_data <- make_tax_palette(data, rank = "Phylum")
  PhylaPalette_sam <- make_tax_palette(data_sam, rank = "Phylum")

  message("[INFO] Number of phyla in data: ", length(PhylaPalette_data))
  message("[INFO] Number of phyla in data_sam: ", length(PhylaPalette_sam))

  # phyloseq bar plot: all samples
  try({
    p_bar_phy <- plot_bar(data, fill = "Phylum") +
      scale_fill_manual(values = PhylaPalette_data, drop = FALSE)

    ggsave("Bar_phy.pdf", plot = p_bar_phy, limitsize = FALSE)
  }, silent = FALSE)

  # phyloseq bar plot: grouped/faceted
  try({
    p_bar_group_phy <- plot_bar(data_sam, fill = "Phylum", facet_grid = ~Groups) +
      scale_fill_manual(values = PhylaPalette_sam, drop = FALSE)

    ggsave("Bar_group_phy.pdf", plot = p_bar_group_phy, limitsize = FALSE)
  }, silent = FALSE)

  # relative abundance bar plot
  data_sam_relabund <- transform_sample_counts(data_sam, function(x) {
    if (sum(x) == 0) {
      return(x)
    } else {
      return(x / sum(x))
    }
  })

  data_sam_relabund <- fix_tax_rank(data_sam_relabund, rank = "Phylum")
  PhylaPalette_relabund <- make_tax_palette(data_sam_relabund, rank = "Phylum")

  theme_set(theme_grey())

  try({
    p_bar_relative_phy <- plot_bar(data_sam_relabund, fill = "Phylum") +
      geom_bar(stat = "identity", position = "stack") +
      labs(x = "", y = "Relative Abundance\n") +
      theme(panel.background = element_blank()) +
      scale_fill_manual(values = PhylaPalette_relabund, drop = FALSE)

    ggsave("Bar_relative_phy.pdf", plot = p_bar_relative_phy)
  }, silent = FALSE)

  # alpha diversity box plot
  theme_set(theme_bw())
  plot_richness(data_sam,"Groups",measures = c("Shannon","Simpson")) +
    geom_boxplot()
  ggsave("Alpha_diversity.pdf")
  theme_set(theme_grey())

  # alpha diversity comparisons
  data.alpha_g<-merge(data.alpha,coldata,by.x="row.names",by.y=1)
  library(ggpubr)
  my_comparisons<-list()
  for (i in 1:nrow(coldata_vs)){
    my_comparisons[[i]] <- c(coldata_vs$group1[i],coldata_vs$group2[i])
  }
  ggboxplot(data.alpha_g, x = "group", y = "Shannon",
                             color = "group", palette = "jco",
                             add = "jitter") +
    stat_compare_means(comparisons = my_comparisons,
                       method = "t.test") # Add pairwise comparisons p-value
  ggsave("Alpha_diversity_Shannon.pdf")

  ggboxplot(data.alpha_g, x = "group", y = "Simpson",
                             color = "group", palette = "jco",
                             add = "jitter") +
    stat_compare_means(comparisons = my_comparisons,
                       method = "t.test")
  ggsave("Alpha_diversity_Simpson.pdf")
}

## draw volcano and bar plot
#### BEGIN BLOCK: bounded DEG volcano and barplot ####

deg_volcano_label_per_direction <- suppressWarnings(
  as.integer(
    Sys.getenv(
      "MTD_VOLCANO_LABELS_PER_DIRECTION",
      "20"
    )
  )
)

if (
  is.na(deg_volcano_label_per_direction) ||
  deg_volcano_label_per_direction < 1
) {
  deg_volcano_label_per_direction <- 20L
}

deg_barplot_max_per_direction <- suppressWarnings(
  as.integer(
    Sys.getenv(
      "MTD_DEG_BARPLOT_MAX_PER_DIRECTION",
      "50"
    )
  )
)

if (
  is.na(deg_barplot_max_per_direction) ||
  deg_barplot_max_per_direction < 1
) {
  deg_barplot_max_per_direction <- 50L
}

for (i in 1:nrow(coldata_vs)){
  group1<-coldata_vs$group1[i]
  group2<-coldata_vs$group2[i]
  pvalue_name<-paste0("padj.",group1,"_vs_",group2) # get the adjusted p-value
  log2FoldChange_name<-paste0("log2FoldChange.",group1,"_vs_",group2)
  # nt4v<-cbind(DEG[,"gene_name"],DEG[,pvalue_name],
  #          DEG[,log2FoldChange_name])
  if (filename == "host_counts.txt"){
    nt4v<-cbind(DEG["gene_name"],DEG[pvalue_name],
                DEG[log2FoldChange_name])
    colnames(nt4v)<-c("gene_name","pvalue","log2FoldChange")
    # nt4v<-as.data.frame(nt4v) # volcano plot require data.frame
    # nt4v$pvalue<-as.numeric(nt4v$pvalue) #character to numeric
    # nt4v$log2FoldChange<-as.numeric(nt4v$log2FoldChange)
    # pass undefined gene name
    while ("-" %in% nt4v$gene_name){
      nt4v$gene_name[match("-",nt4v$gene_name)] <-
        DEG$gene_id[match("-",DEG$gene_name)]
    }
    # drop duplicated row
    nt4v<-nt4v[!duplicated(nt4v$gene_name),]
    # The significantly differentially expressed genes are the ones found in the upper-left and upper-right corners.
    # Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2FoldChange respectively positive or negative)
    # add a column of NAs
    nt4v$diffexpressed <- "NO"
    # if log2Foldchange > 0.5 and pvalue < 0.05, set as "UP"
    nt4v$diffexpressed[nt4v$log2FoldChange > 0.5 & nt4v$pvalue < 0.05] <- "UP"
    # if log2Foldchange < -0.5 and pvalue < 0.05, set as "DOWN"
    nt4v$diffexpressed[nt4v$log2FoldChange< -0.5 & nt4v$pvalue < 0.05] <- "DOWN"
    ## Volcano plotting
    # v<-ggplot(data=nt4v,
    #           aes(x=log2FoldChange,y=-log10(pvalue)),
    #           col=diffexpressed) +
    #   geom_point() +
    #   theme_minimal()
    # # Add vertical lines for log2FoldChange thresholds, and one horizontal line for the p-value threshold
    # v2<-v+geom_vline(xintercept = c(-0.5,0.5), col="red") +
    #   geom_hline(yintercept=-log10(0.05), col="red")
    # # Change point color
    # mycolors <- c("blue", "red", "black")
    # names(mycolors) <- c("DOWN", "UP", "NO")
    # v3 <- v2 + scale_colour_manual(values = mycolors)

    diffexpressed <- nt4v[
      nt4v$diffexpressed != "NO" &
        !is.na(nt4v$log2FoldChange) &
        is.finite(nt4v$log2FoldChange),
      ,
      drop = FALSE
    ]

    topFC_up <- head(
      diffexpressed[
        order(
          diffexpressed$log2FoldChange,
          decreasing = TRUE
        ),
        ,
        drop = FALSE
      ],
      deg_volcano_label_per_direction
    )

    topFC_down <- head(
      diffexpressed[
        order(
          diffexpressed$log2FoldChange,
          decreasing = FALSE
        ),
        ,
        drop = FALSE
      ],
      deg_volcano_label_per_direction
    )

    topFC <- unique(
      rbind(
        topFC_up,
        topFC_down
      )
    )

    topFC <- topFC[
      !duplicated(topFC$gene_name),
      ,
      drop = FALSE
    ]

    nt4v_label <- nt4v
    nt4v_label$label <- NA_character_

    label_rows <- nt4v_label$gene_name %in% topFC$gene_name

    nt4v_label$label[label_rows] <- as.character(
      nt4v_label$gene_name[label_rows]
    )
  } else {
    nt4v<-cbind(DEG["Name"],DEG[pvalue_name],
                DEG[log2FoldChange_name])
    colnames(nt4v)<-c("Name","pvalue","log2FoldChange")
    # pass undefined gene name
    while ("-" %in% nt4v$Name){
      nt4v$Name[match("-",nt4v$Name)] <-
        DEG$gene_id[match("-",DEG$Name)]
    }
    # drop duplicated row
    nt4v<-nt4v[!duplicated(nt4v$Name),]
    # The significantly differentially expressed genes are the ones found in the upper-left and upper-right corners.
    # Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2FoldChange respectively positive or negative)
    # add a column of NAs
    nt4v$diffexpressed <- "NO"
    # if log2Foldchange > 0.5 and pvalue < 0.05, set as "UP"
    nt4v$diffexpressed[nt4v$log2FoldChange > 0.5 & nt4v$pvalue < 0.05] <- "UP"
    # if log2Foldchange < -0.5 and pvalue < 0.05, set as "DOWN"
    nt4v$diffexpressed[nt4v$log2FoldChange< -0.5 & nt4v$pvalue < 0.05] <- "DOWN"

    diffexpressed <- nt4v[
      nt4v$diffexpressed != "NO" &
        !is.na(nt4v$log2FoldChange) &
        is.finite(nt4v$log2FoldChange),
      ,
      drop = FALSE
    ]

    topFC_up <- head(
      diffexpressed[
        order(
          diffexpressed$log2FoldChange,
          decreasing = TRUE
        ),
        ,
        drop = FALSE
      ],
      deg_volcano_label_per_direction
    )

    topFC_down <- head(
      diffexpressed[
        order(
          diffexpressed$log2FoldChange,
          decreasing = FALSE
        ),
        ,
        drop = FALSE
      ],
      deg_volcano_label_per_direction
    )

    topFC <- unique(
      rbind(
        topFC_up,
        topFC_down
      )
    )

    topFC <- topFC[
      !duplicated(topFC$Name),
      ,
      drop = FALSE
    ]

    nt4v_label <- nt4v
    nt4v_label$label <- NA_character_

    label_rows <- nt4v_label$Name %in% topFC$Name

    nt4v_label$label[label_rows] <- as.character(
      nt4v_label$Name[label_rows]
    )
  }

## Volcano plotting
  # plot adding up all layers we have seen so far
  ggplot(data=nt4v_label, aes(x=log2FoldChange, y=-log10(pvalue),
                              col=diffexpressed, label=label)) +
    geom_point(size = 1.5) +
    theme_minimal() +
    geom_text_repel(
      size = 3,
      max.overlaps = Inf,
      na.rm = TRUE
    ) +
    scale_color_manual(values=c("blue", "black", "red")) +
    geom_vline(xintercept=c(-0.5, 0.5), col="red", linetype="dashed", linewidth = 0.5) +
    geom_hline(yintercept=-log10(0.05), col="red", linetype="dashed", linewidth = 0.5) +
    ggtitle(paste0(group1,"_vs_",group2)) +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste0(group1,"_vs_",group2,"/Volcano_",group1,"_vs_",group2,".pdf"))

  ## bar plotting
  bar <- diffexpressed

  if (nrow(bar) != 0) {

    bar_up <- head(
      bar[
        order(
          bar$log2FoldChange,
          decreasing = TRUE
        ),
        ,
        drop = FALSE
      ],
      deg_barplot_max_per_direction
    )

    bar_down <- head(
      bar[
        order(
          bar$log2FoldChange,
          decreasing = FALSE
        ),
        ,
        drop = FALSE
      ],
      deg_barplot_max_per_direction
    )

    bar <- unique(
      rbind(
        bar_up,
        bar_down
      )
    )

    bar <- bar[
      order(
        bar$log2FoldChange,
        decreasing = TRUE
      ),
      ,
      drop = FALSE
    ]

    message(
      "[DEG BARPLOT] Significant genes available: ",
      nrow(diffexpressed),
      "; genes plotted: ",
      nrow(bar),
      " (maximum ",
      2L * deg_barplot_max_per_direction,
      ")"
    )

    # sort log2FoldChange top to down
    if (filename == "host_counts.txt"){
      bar$c <- with(bar,reorder(gene_name,log2FoldChange))
    } else {
      bar$c <- with(bar,reorder(Name,log2FoldChange))
    }
  
    bar_plot<-ggplot(bar,aes(x=log2FoldChange,y=c,fill=pvalue))+
      geom_bar(stat="identity",aes(fill=diffexpressed), width=0.5) +
      #scale_fill_brewer(palette="Blues")+
      scale_fill_manual(name="Expression",
                        labels = c("Down", "Up"),
                        values = c("DOWN"="#00ba38", "UP"="#f8766d")) +
      #scale_fill_gradient2(low=rgb(14,37,56,max=255),high=rgb(54,169,243,max=255), mid=rgb(52,109,157,max=255), midpoint=0.01)+
      #scale_fill_gradient(low=rgb(54,169,243,max=255), high=rgb(14,37,56,max=255))+
      #scale_fill_continuous_sequential(palette = "Blues3", begin=0.4)+
      labs(x="log2FoldChange",y=" ") + #could add ,fill="n=", title = "BMC_HvsBMC_L_DOWN")
      ggtitle(paste0(group1,"_vs_",group2)) +
      #geom_text(aes(x=log2FoldChange+0.6),label="*") +
      theme_bw() +
      scale_y_discrete(breaks=bar[,1],labels=str_trunc(bar[,1],40)) +
      guides(fill = guide_legend(reverse = TRUE)) #reverse the legend to put "Up" in an upper position
    deg_barplot_height <- max(
      6,
      min(
        30,
        0.24 * nrow(bar)
      )
    )

    ggsave(
      paste0(
        group1,
        "_vs_",
        group2,
        "/Barplot_",
        group1,
        "_vs_",
        group2,
        ".pdf"
      ),
      plot = bar_plot,
      height = deg_barplot_height,
      width = 10,
      limitsize = FALSE
    )
  }
}

#### END BLOCK: bounded DEG volcano and barplot ####

## Venn Diagram
if (length(unique(coldata$group))>=2){
  dv.list<-list()
  for (i in unique(coldata$group)){
    vd<-coldata[coldata$group==i,]
    sample_name<-as.character(vd[,1])
    sample_name.norm<-paste0(sample_name,".norm")
    if (filename == "host_counts.txt"){
      vd.df<-DEG[c("gene_name",sample_name.norm)]
    } else {
      vd.df<-DEG[c("Name",sample_name.norm)]
    }
    #Removing rows having all zeros, ignore 1st column with names
    vd.df<-vd.df[rowSums(vd.df[-1])>0,]
    #add to list
    dv.list[[i]]<-vd.df[,1]
  }

  #myCol <- brewer.pal(length(unique(coldata$group)), "Pastel2")

  library(VennDiagram)
  # Don't write log file for VennDiagram
  futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
  #venn.diagram(...)

  venn.diagram(
    x = dv.list,
    category.names = names(dv.list),
    filename = 'venn_diagramm.png',
    output=TRUE,

    # # Output features
     imagetype="png" ,
    # height = 480 ,
    # width = 480 ,
    # resolution = 300,
     compression = "lzw",

    # Circles
    size = 1,
    #lty = 'blank',
    col=cols,
    #fill = myCol,
    # keep same colors with heatmap and PCA
    fill = alpha(cols,0.6),

    # Numbers
    #cex = 1,
    fontface = "bold",
    fontfamily = "sans",

    # Set names
    #cat.cex = 1,
    cat.fontface = "bold",
    # cat.default.pos = "outer",
    # cat.pos = c(-27, 27, 135),
    # cat.dist = c(0.055, 0.055, 0.085),
    cat.fontfamily = "sans",
    #rotation = 1
  )
}


# non_host_host_reads_ratio comparison

if (filename %in% c("bracken_species_all","bracken_phylum_all","bracken_genus_all")){
  files_h <- list.files(paste0(dirname(args[1]),"/temp"), pattern="^Report_host_.*\\.txt$", full.names=TRUE, recursive=FALSE)
  host_summary_path <- find_host_summary_path(
    dirname(args[1])
  )

  host_ratios <- read_host_ratios(
    summary_path = host_summary_path,
    report_files = files_h
  )

  lh <- host_ratios$non_host_host_reads_ratio
  la <- host_ratios$unclassified_reads_ratio_percent

  lah<-cbind(la,lh)
  lah<-merge(lah,coldata,by.x="row.names",by.y="sample_name")
  colnames(lah)<-c("Sample","unclassified_reads_ratio_percent","non_host_host_reads_ratio","Groups")

  # to draw box plot with comparison of unclassified ratio
  my_comparisons<-list()
  for (i in 1:nrow(coldata_vs)){
    my_comparisons[[i]] <- c(coldata_vs$group1[i],coldata_vs$group2[i])
  }
  ggboxplot(lah, x = "Groups", y = "non_host_host_reads_ratio",
            color = "Groups", palette = "jco",
            add = "jitter") +
    stat_compare_means(comparisons = my_comparisons,
                       method = "t.test") # Add pairwise comparisons p-value
  ggsave("non-host_vs_host_reads_ratio.pdf")

  ggboxplot(lah, x = "Groups", y = "unclassified_reads_ratio_percent",
            color = "Groups", palette = "jco",
            add = "jitter") +
    stat_compare_means(comparisons = my_comparisons,
                       method = "t.test") # Add pairwise comparisons p-value
  ggsave("unclassified_reads_ratio.pdf")
}


# preprocess for graphlan
if (filename %in% c("bracken_species_all","bracken_phylum_all","bracken_genus_all")){
gph.anno<-read.table(paste0(dirname(args[1]),"/graphlan/annot.txt"),
                     header = F,blank.lines.skip=F, fill=TRUE,sep="\t",
                     comment.char = "",
                     col.names=c("V1","V2","V3","V4"))
gph.anno<-gph.anno[!(gph.anno[,1]=="" & gph.anno[,2]=="annotation" & gph.anno[,3]=="") ,] #clean the Unrecognized annotationline

gph.tree<-read.table(paste0(dirname(args[1]),"/graphlan/tree.txt"),
                     header = F, fill=TRUE)

library(stringi)
sp.name<-row.names(cts)
for ( i in 1:nrow(gph.tree)){
  if (length(strsplit(gph.tree[i,1],"\\.")[[1]])>=7){
    gph.sp<-tail(strsplit(gph.tree[i,1],"\\.")[[1]],1) #extract the last word of the species name from graphlan tree
    gph.sp1<-gsub("_","-",gph.sp) # To match the bracken style
    for ( n in 1:length(sp.name)){
      sp.name1<-tail(strsplit(sp.name," ")[[n]],1) #extract the last word of the species name from bracken file
      if (gph.sp1==sp.name1 || gph.sp==sp.name1){
        gph.sp2<-head(tail(strsplit(gph.tree[i,1],"\\.")[[1]],2),1) #extract the second last word of the species name from graphlan tree
        gph.sp6<-strsplit(gph.tree[i,1],"\\.")[[1]][6] #extract the sixth word of the species name from graphlan tree
        sp.name2<-head(strsplit(sp.name," ")[[n]],1) #extract the first word of the species name from bracken file
        if (grepl(sp.name2,gph.sp2) || grepl(sp.name2,gph.sp6)){
          sp.replace<-sp.name[n]
          sp.replace<-gsub(" ","_",sp.replace)
          #gph.tree[i,1]<-sub(paste0(gph.sp,"([^",gph.sp,"]*)$"),paste0(sp.replace,"\\1"), gph.tree[i,1])
          gph.tree[i,1]<-stri_replace_last_fixed(gph.tree[i,1],gph.sp,sp.replace)
        }
      }
    }
  }
}
write.table(gph.tree,paste0(dirname(args[1]),"/graphlan/tree.txt"),
            col.names = F, row.names = F, quote = F)

start_time <- Sys.time()
sp.replace<-gsub(" ","_",sp.name)
for (i in 31:nrow(gph.anno)){
  for ( n in 1:length(sp.name)){
    sp.name1<-tail(strsplit(sp.name," ")[[n]],1)
    sp.name1<-gsub("-","_",sp.name1) # To match the graphlan anno style
    sp.name2<-head(strsplit(sp.name," ")[[n]],1)
    if (length(grep(sp.name1,gph.anno[[1]])) > 1){
      # check if names match in both of species and genus level
      if (gph.anno[i,1]==sp.name1 & TRUE %in% grepl(sp.name2,gph.anno[(i-30):(i-1),1])){
        gph.anno[i,3]<-gsub(" ","_",gph.anno[i,3])
        if (grepl(gph.anno[i,1],gph.anno[i,3])){
          gph.anno[i,3]<-sub(gph.anno[i,1],sp.replace[n],gph.anno[i,3])
        }
        gph.anno[i,1]<-sp.replace[n]
      }
    } else {
      if (gph.anno[i,1]==sp.name1){
        if (grepl(gph.anno[i,1],gph.anno[i,3])){
          gph.anno[i,3]<-sub(gph.anno[i,1],sp.replace[n],gph.anno[i,3])
        }
        gph.anno[i,1]<-sp.replace[n]
      }
    }
  }
}
Sys.time() - start_time

# add rings to graphlan
ring.list<-list() # make a list for ring color of each group
for (i in unique(coldata$group)){
  ring.gp<-coldata[coldata$group==i,]
  sample_name<-as.character(ring.gp[,1])
  ring.df<-cts[sample_name]
  #Removing rows having all zeros
  ring.df<-ring.df[rowSums(ring.df)>0,]
  ring.df<-ring.df[complete.cases(ring.df),] # remove the NA rows
  ring.df<-gsub(" ","_",row.names(ring.df))
  ring.df<-gsub("\\[||\\]","",ring.df)
  #add to list
  ring.list[[i]]<-ring.df
}

group4gph<-list()
white.list<-data.frame()
for (i in 1:nrow(coldata_vs)){
  group1<-coldata_vs$group1[i]
  group2<-coldata_vs$group2[i]
  gp.vs<-paste0(group1,"_vs_",group2)
  gp.vs.p<-paste0("pvalue.",gp.vs)
  gp.vs.fc<-paste0("log2FoldChange.",gp.vs)
  group.p<-DEG[DEG[gp.vs.p]<0.05 & DEG[gp.vs.fc]>0.5,]
  group.gph<-group.p[["Name"]] # extract significant diff. expressed species
  group.gph<-gsub(" ","_",group.gph)
  group.gph<-gsub("\\[||\\]","",group.gph)
  level<-grep(group1,names(anno_colors$Group))
  ring.color<-ring.list[group1]
  ring.color<-ring.color[[1]]
  ring.color<-ring.color[which(!ring.color %in% group.gph)] #exclude the overlap color
  white.list<-append(white.list,data.frame(temp=group.gph)) # add sig. exp. in a white list with a temporary name
  names(white.list)[names(white.list)=="temp"]<-level #rename to level names
  group.gph<-rbind(c("ring_label", level, group1,""),
                   c("ring_label_color", level, anno_colors$Group[group1],""),
                   c("ring_internal_separator_thickness",level,"1",""),
                   c("ring_external_separator_thickness",level,"1",""),
                   cbind(group.gph,
                         rep("ring_color",length(group.gph)),
                         rep(level,length(group.gph)),
                         rep(anno_colors$Group[group1],length(group.gph))),
                   cbind(ring.color,
                         rep("ring_color",length(ring.color)),
                         rep(level,length(ring.color)),
                         rep(anno_colors$Group[group1],length(ring.color))),
                   cbind(ring.color,
                         rep("ring_alpha",length(ring.color)),
                         rep(level,length(ring.color)),
                         rep(0.2,length(ring.color))))

  group.p2<-DEG[DEG[gp.vs.p]<0.05 & DEG[gp.vs.fc]<(-0.5),]
  group.gph2<-group.p2[["Name"]]
  group.gph2<-gsub(" ","_",group.gph2)
  group.gph2<-gsub("\\[||\\]","",group.gph2)
  level2<-grep(group2,names(anno_colors$Group))
  ring.color2<-ring.list[group2]
  ring.color2<-ring.color2[[1]]
  ring.color2<-ring.color2[which(!ring.color2 %in% group.gph2)] #exclude the overlap color
  white.list<-append(white.list,data.frame(temp=group.gph2)) # add sig. exp. in a white list with a temporary name
  names(white.list)[names(white.list)=="temp"]<-level2 #rename to level names
  group.gph2<-rbind(c("ring_label", level2, group2,""),
                    c("ring_label_color", level2, anno_colors$Group[group2],""),
                    c("ring_internal_separator_thickness",level2,"1",""),
                    c("ring_external_separator_thickness",level2,"1",""),
                    cbind(group.gph2,
                          rep("ring_color",length(group.gph2)),
                          rep(level2,length(group.gph2)),
                          rep(anno_colors$Group[group2],length(group.gph2))),
                    cbind(ring.color2,
                          rep("ring_color",length(ring.color2)),
                          rep(level2,length(ring.color2)),
                          rep(anno_colors$Group[group2],length(ring.color2))),
                    cbind(ring.color2,
                          rep("ring_alpha",length(ring.color2)),
                          rep(level2,length(ring.color2)),
                          rep(0.2,length(ring.color2))))
  group4gph[[i]]<-rbind(unname(group.gph),unname(group.gph2))
}
group4gph.table<-do.call(rbind.data.frame,group4gph)
group4gph.table<-unique(group4gph.table)
names(group4gph.table)[1]<-"V1"

# remove "ring_alpha" from the sig. exp. white list
#group4gph.table1<-group4gph.table[!(group4gph.table[,1]==white.list & group4gph.table[,2]=="ring_alpha"),]
#group4gph.table<-group4gph.table[!((group4gph.table[,1] %in% white.list) & group4gph.table[,2]=="ring_alpha"),]
library(reshape2)
white.ls<-melt(white.list,na.rm=T)
for (g in unique(white.ls$L1)){
  group4gph.table<-group4gph.table[!((group4gph.table[,1] %in% white.ls[white.ls$L1==g,1]) &
                                       group4gph.table[,2]=="ring_alpha" &
                                       group4gph.table[,3]==g),]
}

# gatekeeper of species name format between graphlan and bracken
gph.match<-grep(paste(group4gph.table[[1]],collapse="|"), gph.anno[[1]], value=TRUE)
gph.nomatch<-unique(grep(paste(gph.match,collapse="|"),group4gph.table[[1]],value=TRUE, invert=T))
gph.nomatch<-gph.nomatch[which(!gph.nomatch %in%
                                 c("ring_label","ring_label_color",
                                   "ring_color","ring_alpha",
                                   "ring_internal_separator_thickness",
                                   "ring_external_separator_thickness"))]

for ( n in gph.nomatch){
  gph.nomatch1<-tail(strsplit(n,"_")[[1]],1) #extract the last word of the species name from bracken file
  gph.nomatch2<-head(strsplit(n,"_")[[1]],1) #extract the first word of the species name from bracken file
  count4del<-0
  for (m in gph.anno[[1]]){
    # replace only both 1st and last words matched terms
    if (grepl(gph.nomatch1,m) & grepl(gph.nomatch2,m)){
      group4gph.table[[1]]<-gsub(n, m, group4gph.table[[1]])
      count4del<-count4del+1
    }
  }
  if (count4del==0){
    group4gph.table<-group4gph.table[group4gph.table[,1]!=n,]
  }
}

gph.anno1<-rbind(gph.anno,c("total_plotted_degrees","330","",""),
                 c("start_rotation","270","",""),
                 group4gph.table)
write.table(gph.anno1,paste0(dirname(args[1]),"/graphlan/annot.txt"),
            col.names = F, row.names = F, quote = F, sep="\t", na = "")
}


# adjust covariance effect
#### BEGIN BLOCK: paired HAllA adjustment ####
#
# The experimental group is preserved in the expression matrix.
# Repeated measurements from the same animal/subject are handled as
# a categorical blocking factor. Only genuinely numeric metadata
# columns are passed as continuous covariates.
#
# This avoids converting sample names and categorical identifiers
# with as.numeric(), which previously generated NA values.

normtrans_adj <- normtrans

if (length(args) == 5) {

  halla_metadata <- as.data.frame(
    coldata,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  sample_column <- names(halla_metadata)[1]
  sample_ids <- as.character(
    halla_metadata[[sample_column]]
  )

  rownames(halla_metadata) <- sample_ids

  missing_halla_samples <- setdiff(
    colnames(normtrans),
    rownames(halla_metadata)
  )

  if (length(missing_halla_samples) > 0) {
    stop(
      "HAllA metadata is missing samples present in the ",
      "normalized host matrix: ",
      paste(missing_halla_samples, collapse = ", ")
    )
  }

  halla_metadata <- halla_metadata[
    colnames(normtrans),
    ,
    drop = FALSE
  ]

  if (!"group" %in% names(halla_metadata)) {
    stop(
      "The metadata must contain a group column for ",
      "HAllA preprocessing."
    )
  }

  halla_metadata$group <- factor(
    halla_metadata$group
  )

  # Preserve the biological comparison while removing nuisance effects.
  halla_design <- model.matrix(
    ~ group,
    data = halla_metadata
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
    names(halla_metadata)
  )

  subject_indices <- match(
    subject_candidates,
    metadata_names_lower,
    nomatch = 0L
  )

  subject_indices <- subject_indices[
    subject_indices > 0L
  ]

  halla_subject_column <- NULL

  for (subject_index in subject_indices) {
    candidate_name <- names(halla_metadata)[subject_index]

    candidate_values <- trimws(
      as.character(
        halla_metadata[[candidate_name]]
      )
    )

    candidate_values <- candidate_values[
      !is.na(candidate_values) &
        candidate_values != ""
    ]

    if (
      length(unique(candidate_values)) > 1 &&
      anyDuplicated(candidate_values) > 0
    ) {
      halla_subject_column <- candidate_name
      break
    }
  }

  halla_batch <- NULL

  if (!is.null(halla_subject_column)) {
    halla_batch <- factor(
      halla_metadata[[halla_subject_column]]
    )
  }

  excluded_columns <- c(
    sample_column,
    "group",
    halla_subject_column
  )

  nuisance_candidates <- setdiff(
    names(halla_metadata),
    excluded_columns
  )

  numeric_covariates <- list()
  ignored_categorical_covariates <- character()

  for (column_name in nuisance_candidates) {

    values <- halla_metadata[[column_name]]

    if (is.numeric(values) || is.integer(values)) {

      numeric_values <- as.numeric(values)

    } else {

      character_values <- trimws(
        as.character(values)
      )

      numeric_values <- suppressWarnings(
        as.numeric(character_values)
      )

      conversion_valid <- (
        is.na(character_values) |
          character_values == "" |
          !is.na(numeric_values)
      )

      if (!all(conversion_valid)) {
        numeric_values <- NULL
      }
    }

    if (!is.null(numeric_values)) {

      nonmissing_values <- numeric_values[
        !is.na(numeric_values)
      ]

      if (length(unique(nonmissing_values)) > 1) {
        numeric_covariates[[column_name]] <- numeric_values
      }

    } else {

      categorical_values <- as.character(values)
      categorical_values <- categorical_values[
        !is.na(categorical_values) &
          categorical_values != ""
      ]

      if (length(unique(categorical_values)) > 1) {
        ignored_categorical_covariates <- c(
          ignored_categorical_covariates,
          column_name
        )
      }
    }
  }

  halla_covariate_matrix <- NULL

  if (length(numeric_covariates) > 0) {
    halla_covariate_matrix <- as.matrix(
      as.data.frame(
        numeric_covariates,
        check.names = FALSE
      )
    )
  }

  remove_batch_arguments <- list(
    x = normtrans,
    design = halla_design
  )

  if (!is.null(halla_batch)) {
    remove_batch_arguments$batch <- halla_batch
  }

  if (!is.null(halla_covariate_matrix)) {
    remove_batch_arguments$covariates <- halla_covariate_matrix
  }

  if (
    !is.null(halla_batch) ||
    !is.null(halla_covariate_matrix)
  ) {
    normtrans_adj <- do.call(
      limma::removeBatchEffect,
      remove_batch_arguments
    )
  }

  message(
    "[HAllA] Biological group preserved: group"
  )

  message(
    "[HAllA] Repeated-measures block: ",
    if (is.null(halla_subject_column)) {
      "none"
    } else {
      halla_subject_column
    }
  )

  message(
    "[HAllA] Numeric covariates removed: ",
    if (length(numeric_covariates) == 0) {
      "none"
    } else {
      paste(
        names(numeric_covariates),
        collapse = ", "
      )
    }
  )

  if (length(ignored_categorical_covariates) > 0) {
    message(
      "[HAllA] Additional categorical metadata preserved: ",
      paste(
        ignored_categorical_covariates,
        collapse = ", "
      )
    )
  }

} else {

  # Without extended metadata there is no justified batch or
  # repeated-measures variable to remove. Preserve the transformed
  # matrix and retain the biological group effect.
  normtrans_adj <- normtrans

  message(
    "[HAllA] No extended metadata supplied; ",
    "no batch effect was removed and the biological group was preserved."
  )
}
#### END BLOCK: paired HAllA adjustment ####
# make a folder for outputs
dir.create("../halla",recursive = T)
setwd("../halla")
# save normalized & transformed & covariance corrected data for downstream correlation analysis (e.g. halla)
if (filename == "host_counts.txt"){
  DEG_name<-DEG[c("gene_id","hybrid_name")]
  normtrans_adj<-merge(DEG_name,normtrans_adj,by.x="gene_id",by.y="row.names")
  normtrans_adj<-normtrans_adj[!duplicated(normtrans_adj[2]),]
  #convert column name (gene_name) to row name
  normtrans_adj2<-normtrans_adj[,-1:-2]
  rownames(normtrans_adj2)<-normtrans_adj[,2]
  write.table(normtrans_adj2,"Host_gene.txt",sep="\t",quote=F,col.names=NA)
} else if (filename == "humann_genefamilies_Abundance_go_translated.tsv"){
  write.table(normtrans_adj,"Microbiomes_humann_go.txt",sep="\t",quote=F,col.names=NA)
} else if (filename == "humann_genefamilies_Abundance_kegg_translated.tsv"){
  write.table(normtrans_adj,"Microbiomes_humann_kegg.txt",sep="\t",quote=F,col.names=NA)
} else {
  write.table(normtrans_adj,"Microbiomes.txt",sep="\t",quote=F,col.names=NA)
}
setwd("../")


### Pathway enrichment for host genes ###
if (filename == "host_counts.txt"){
  setwd("Host_DEG")

  if (!"Taxon_ID" %in% names(host_sp)) {
    stop(
      "HostSpecies.csv does not contain a Taxon_ID column."
    )
  }

  if (!"OrgDb" %in% names(host_sp)) {
    stop(
      "HostSpecies.csv does not contain an OrgDb column."
    )
  }

  host_taxid <- trimws(
    as.character(args[3])
  )

  csv_taxid <- trimws(
    as.character(host_sp$Taxon_ID)
  )

  # which() discards NA comparisons instead of returning an
  # additional NA row in the selected OrgDb vector.
  host_rows <- which(
    !is.na(csv_taxid) &
      csv_taxid != "" &
      csv_taxid == host_taxid
  )

  do.db <- unique(
    trimws(
      as.character(
        host_sp$OrgDb[host_rows]
      )
    )
  )

  do.db <- do.db[
    !is.na(do.db) &
      do.db != "" &
      do.db != "-"
  ]

  if (length(do.db) == 0) {
    message(
      "[INFO] No OrgDb package is configured for host TaxID ",
      host_taxid,
      ". Host pathway enrichment will be skipped."
    )
  } else {

    if (length(do.db) > 1) {
      stop(
        "Multiple OrgDb packages are configured for TaxID ",
        host_taxid,
        ": ",
        paste(do.db, collapse = ", ")
      )
    }

    do.db <- do.db[[1]]

    message(
      "[INFO] Host OrgDb package: ",
      do.db
    )

    message(
      "[INFO] Active R library paths: ",
      paste(.libPaths(), collapse = "; ")
    )

    if (!requireNamespace(do.db, quietly = TRUE)) {
      stop(
        "The configured OrgDb package is not visible: ",
        do.db,
        "\nActive R library paths: ",
        paste(.libPaths(), collapse = "; ")
      )
    }

    suppressPackageStartupMessages({
      library(clusterProfiler)
      library(enrichplot)
      library(ggnewscale)
      library(do.db, character.only = TRUE)
      library(ggplot2)
    })

    do.db.obj <- get(
      do.db,
      envir = asNamespace(do.db)
    )
    
#### BEGIN BLOCK: OrgDb keytype auto-detection ####
get_orgdb_gene_keytype <- function(orgdb) {
  kt <- AnnotationDbi::keytypes(orgdb)
  cols <- AnnotationDbi::columns(orgdb)

  message("[INFO] OrgDb keytypes available:")
  print(kt)

  message("[INFO] OrgDb columns available:")
  print(cols)

  if ("ENTREZID" %in% kt && "ENTREZID" %in% cols) {
    return("ENTREZID")
  }

  if ("GID" %in% kt && "GID" %in% cols) {
    return("GID")
  }

  if ("ENSEMBL" %in% kt && "ENSEMBL" %in% cols) {
    return("ENSEMBL")
  }

  if ("SYMBOL" %in% kt && "SYMBOL" %in% cols) {
    return("SYMBOL")
  }

  stop(
    "No usable gene keytype found in OrgDb. Available keytypes: ",
    paste(kt, collapse = ", ")
  )
}

HOST_GENE_KEYTYPE <- get_orgdb_gene_keytype(do.db.obj)

HOST_ORGDB_COLUMNS <- AnnotationDbi::columns(do.db.obj)

HOST_HAS_ENTREZ <- "ENTREZID" %in% HOST_ORGDB_COLUMNS

HOST_HAS_KEGG_DERIVED <- any(
  c(
    "KEGG_KO",
    "KEGG_PATHWAY",
    "KEGG_MODULE",
    "KEGG_REACTION",
    "KEGG_BRITE"
  ) %in% HOST_ORGDB_COLUMNS
)

message("[INFO] Host OrgDb gene keytype selected: ", HOST_GENE_KEYTYPE)
message("[INFO] Host OrgDb has ENTREZID column: ", HOST_HAS_ENTREZ)
message("[INFO] Host OrgDb has eggNOG-derived KEGG columns: ", HOST_HAS_KEGG_DERIVED)
#### END BLOCK: OrgDb keytype auto-detection ####

    # function to make plots for GSEA results
plots4gsea <- function(edb, data, datax, edb0, genelist, group1, group2) {

  dir.create(edb, showWarnings = FALSE, recursive = TRUE)

  #### BEGIN BLOCK: legacy enrichplot compatibility ####
  #
  # enrichplot 1.14.x builds cnetplot category-node colours with
  # aes(color = I(...)). With newer ggplot2/scales combinations,
  # this AsIs mapping can reach scales::rescale() and fail.
  #
  # Convert only these legacy identity-colour mappings into fixed
  # layer parameters. Numeric fold-change colour mappings remain
  # untouched.

  fix_legacy_cnetplot_as_is <- function(plot_object) {

    if (
      is.null(plot_object) ||
      is.null(plot_object$layers) ||
      length(plot_object$layers) == 0
    ) {
      return(plot_object)
    }

    fixed_layers <- 0L
    identity_layer_number <- 0L

    default_identity_colours <- c(
      "#E5C494",
      "#B3B3B3"
    )

    for (layer_index in seq_along(plot_object$layers)) {

      layer <- plot_object$layers[[layer_index]]

      if (
        is.null(layer$mapping) ||
        length(layer$mapping) == 0
      ) {
        next
      }

      mapping_names <- names(layer$mapping)

      colour_mapping_name <- intersect(
        c("colour", "color"),
        mapping_names
      )

      if (length(colour_mapping_name) == 0) {
        next
      }

      colour_mapping_name <- colour_mapping_name[1]

      mapping_expression <- layer$mapping[[colour_mapping_name]]

      mapping_text <- tryCatch(
        {
          rlang::as_label(
            mapping_expression
          )
        },
        error = function(e) {
          paste(
            deparse(mapping_expression),
            collapse = ""
          )
        }
      )

      is_identity_mapping <- grepl(
        "I\\s*\\(",
        mapping_text
      )

      if (!is_identity_mapping) {
        next
      }

      identity_layer_number <- (
        identity_layer_number + 1L
      )

      fixed_colour <- default_identity_colours[
        min(
          identity_layer_number,
          length(default_identity_colours)
        )
      ]

      if (
        grepl(
          "#E5C494",
          mapping_text,
          fixed = TRUE
        )
      ) {
        fixed_colour <- "#E5C494"
      }

      if (
        grepl(
          "#B3B3B3",
          mapping_text,
          fixed = TRUE
        )
      ) {
        fixed_colour <- "#B3B3B3"
      }

      # Remove aes(colour = I(...)).
      layer$mapping[[colour_mapping_name]] <- NULL

      # Add the same colour outside aes().
      layer$aes_params$colour <- fixed_colour

      plot_object$layers[[layer_index]] <- layer

      fixed_layers <- fixed_layers + 1L
    }

    message(
      "[CNETPLOT] Legacy AsIs colour layers repaired: ",
      fixed_layers
    )

    plot_object
  }

  message(
    "[PLOT COMPAT] enrichplot version: ",
    as.character(
      utils::packageVersion("enrichplot")
    )
  )

  message(
    "[PLOT COMPAT] ggplot2 version: ",
    as.character(
      utils::packageVersion("ggplot2")
    )
  )
  #### END BLOCK: legacy enrichplot compatibility ####

  safe_filename <- function(x, max_chars = 120) {
    x <- as.character(x)
    if (length(x) == 0 || is.na(x) || x == "") {
      x <- "term"
    }

    x <- gsub("[/\\\\:*?\"<>|]", "_", x)
    x <- gsub("[[:cntrl:]]", "_", x)
    x <- gsub("\\s+", " ", x)
    x <- trimws(x)

    if (nchar(x) > max_chars) {
      x <- substr(x, 1, max_chars)
    }

    x
  }

  write_plot_skip <- function(plot_name, e) {
    msg <- paste0(plot_name, " skipped: ", conditionMessage(e))
    message(msg)

    write(
      msg,
      paste0(edb, "/GSEA_", edb0, "_", plot_name, "_skipped.txt")
    )
  }

  if (is.null(data) || nrow(data@result) == 0) {
    write(
      "No enrichment result was found.",
      paste0(edb, "/No_enrichment_result_found.txt")
    )
    return(invisible(NULL))
  }

#### BEGIN BLOCK: bounded enrichment graphics ####

  sanitize_enrichment_object <- function(object) {

    if (
      is.null(object) ||
      !methods::is(object, "S4") ||
      !"result" %in% methods::slotNames(object)
    ) {
      return(object)
    }

    result_table <- object@result

    for (column_name in names(result_table)) {
      if (inherits(result_table[[column_name]], "AsIs")) {
        result_table[[column_name]] <- as.vector(
          result_table[[column_name]]
        )
      }
    }

    character_columns <- intersect(
      c(
        "ID",
        "Description",
        "geneID",
        "core_enrichment",
        "setSize"
      ),
      names(result_table)
    )

    for (column_name in character_columns) {
      if (
        column_name != "setSize"
      ) {
        result_table[[column_name]] <- as.character(
          result_table[[column_name]]
        )
      }
    }

    if (
      "Description" %in% names(result_table) &&
      anyDuplicated(result_table$Description) > 0
    ) {
      result_table$Description <- make.unique(
        as.character(
          result_table$Description
        ),
        sep = " | "
      )
    }

    object@result <- result_table
    object
  }

  data <- sanitize_enrichment_object(data)
  datax <- sanitize_enrichment_object(datax)

  genelist_names <- names(genelist)
  genelist <- suppressWarnings(
    as.numeric(genelist)
  )
  names(genelist) <- genelist_names

  genelist <- genelist[
    !is.na(genelist) &
      is.finite(genelist) &
      !is.na(names(genelist)) &
      names(genelist) != ""
  ]

  max_plot_terms <- suppressWarnings(
    as.integer(
      Sys.getenv(
        "MTD_ENRICHMENT_MAX_PLOT_TERMS",
        "30"
      )
    )
  )

  if (
    is.na(max_plot_terms) ||
    max_plot_terms < 5
  ) {
    max_plot_terms <- 30L
  }

  nres <- nrow(data@result)
  plot_terms <- min(
    nres,
    max_plot_terms
  )

  message(
    "[ENRICHMENT PLOTS] Results available: ",
    nres,
    "; terms plotted: ",
    plot_terms
  )

#### END BLOCK: bounded enrichment graphics ####

  ## ----------------------------------------------------------
  ## Top GSEA plots
  ## ----------------------------------------------------------

  tryCatch({
top_n_gsea <- min(20, nres)

p.top <- gseaplot2(
  data,
  geneSetID = 1:top_n_gsea,
  rel_heights = c(2, 0.5, 1)
)

ggsave(
  paste0(edb, "/Top", top_n_gsea, "_GSEA_", edb0, ".pdf"),
  plot = p.top,
  height = 6 + 0.45 * top_n_gsea,
  width = 8
)
  }, error = function(e) {
    write_plot_skip("top_gseaplot", e)
  })

  ## ----------------------------------------------------------
  ## Individual GSEA plots
  ## ----------------------------------------------------------

  dir.create(
    paste0(edb, "/GSEA_all"),
    showWarnings = FALSE,
    recursive = TRUE
  )

  for (g in seq_len(nres)) {
    desc <- data@result$Description[g]
    desc <- safe_filename(desc)

    tryCatch({
      p.g <- gseaplot2(
        data,
        geneSetID = g,
        title = data@result$Description[g]
      )

      ggsave(
        paste0(edb, "/GSEA_all/", desc, "_GSEA_", edb0, ".pdf"),
        plot = p.g,
        height = 7,
        width = 7
      )
    }, error = function(e) {
      message("GSEA individual plot skipped for ", edb, " term ", g, ": ", conditionMessage(e))
      write(
        paste0("GSEA individual plot skipped for term ", g, ": ", conditionMessage(e)),
        paste0(edb, "/GSEA_all/term_", g, "_GSEA_", edb0, "_skipped.txt")
      )
    })
  }

  ## ----------------------------------------------------------
  ## Ridgeplot
  ## ----------------------------------------------------------

  tryCatch({
    p.ridge <- ridgeplot(data, showCategory = plot_terms)

    if (nres < 30) {
      ggsave(
        paste0(edb, "/GSEA_", edb0, "_ridgeplots.pdf"),
        plot = p.ridge,
        height = 2 + 0.48 * nres
      )
    } else {
      ggsave(
        paste0(edb, "/GSEA_", edb0, "_ridgeplots.pdf"),
        plot = p.ridge,
        height = 12
      )
    }
  }, error = function(e) {
    write_plot_skip("ridgeplot", e)
  })

  ## ----------------------------------------------------------
  ## Dotplot
  ## ----------------------------------------------------------

  tryCatch({
    p.dot <- dotplot(data, showCategory = plot_terms) + ggtitle("dotplot for GSEA")

    ggsave(
      paste0(edb, "/GSEA_", edb0, "_dotplot.pdf"),
      plot = p.dot,
      height = 4.8,
      width = 6
    )
  }, error = function(e) {
    write_plot_skip("dotplot", e)
  })

  ## ----------------------------------------------------------
  ## Cnetplot
  ## ----------------------------------------------------------

  tryCatch({
    p.net <- cnetplot(
      datax,
      foldChange = genelist,
      showCategory = plot_terms,
      node_label = "category",
      cex_label_category = 0.85,
      cex_label_gene = 0.6
    )

    p.net <- fix_legacy_cnetplot_as_is(
      p.net
    )

    if (nres < 25) {
      ggsave(
        paste0(edb, "/GSEA_", edb0, "_net.pdf"),
        plot = p.net
      )
    } else {
      ggsave(
        paste0(edb, "/GSEA_", edb0, "_net.pdf"),
        plot = p.net,
        scale = 2
      )
    }
  }, error = function(e) {
    write_plot_skip("net", e)
  })

  ## ----------------------------------------------------------
  ## Term similarity object
  ## ----------------------------------------------------------

  datax2 <- tryCatch({
    pairwise_termsim(datax)
  }, error = function(e) {
    message("pairwise_termsim skipped for ", edb, ": ", conditionMessage(e))
    write(
      paste0("pairwise_termsim skipped: ", conditionMessage(e)),
      paste0(edb, "/GSEA_", edb0, "_pairwise_termsim_skipped.txt")
    )
    NULL
  })

  ## ----------------------------------------------------------
  ## Term-similarity dendrogram
  ## ----------------------------------------------------------
  #
  # enrichplot::treeplot() from the Bioconductor 3.14 environment
  # can fail inside ggtree::geom_segment2() when combined with the
  # installed ggplot2 version:
  #
  #   'gpar' element 'lwd' must not be length 0
  #
  # Build the equivalent term-similarity hierarchy directly from
  # the pairwise_termsim() matrix. This avoids the incompatible
  # ggtree drawing layer while preserving the clustering result.

  if (!is.null(datax2)) {
    tryCatch({

      if (
        !"termsim" %in% methods::slotNames(datax2)
      ) {
        stop(
          "The enrichment object does not contain a termsim slot."
        )
      }

      tree_similarity <- methods::slot(
        datax2,
        "termsim"
      )

      tree_similarity <- as.matrix(
        tree_similarity
      )

      if (
        nrow(tree_similarity) < 2 ||
        ncol(tree_similarity) < 2
      ) {
        stop(
          "Fewer than two terms were available for clustering."
        )
      }

      if (
        is.null(rownames(tree_similarity)) ||
        is.null(colnames(tree_similarity))
      ) {
        stop(
          "The term-similarity matrix has no row or column names."
        )
      }

      tree_results <- as.data.frame(
        datax2@result,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )

      if (nrow(tree_results) < 2) {
        stop(
          "Fewer than two enrichment results were available."
        )
      }

      tree_similarity_names <- rownames(
        tree_similarity
      )

      id_matches <- 0L
      description_matches <- 0L

      if ("ID" %in% names(tree_results)) {
        id_matches <- sum(
          as.character(tree_results$ID) %in%
            tree_similarity_names
        )
      }

      if ("Description" %in% names(tree_results)) {
        description_matches <- sum(
          as.character(tree_results$Description) %in%
            tree_similarity_names
        )
      }

      if (
        id_matches == 0L &&
        description_matches == 0L
      ) {
        stop(
          "Could not match enrichment results to the ",
          "term-similarity matrix."
        )
      }

      tree_key_column <- if (
        id_matches >= description_matches
      ) {
        "ID"
      } else {
        "Description"
      }

      tree_results$.tree_key <- as.character(
        tree_results[[tree_key_column]]
      )

      tree_results <- tree_results[
        !is.na(tree_results$.tree_key) &
          tree_results$.tree_key != "" &
          tree_results$.tree_key %in%
            tree_similarity_names,
        ,
        drop = FALSE
      ]

      if (nrow(tree_results) < 2) {
        stop(
          "Fewer than two matched enrichment terms remained."
        )
      }

      tree_rank_column <- intersect(
        c(
          "p.adjust",
          "pvalue",
          "qvalue"
        ),
        names(tree_results)
      )

      if (length(tree_rank_column) > 0) {

        tree_rank_column <- tree_rank_column[1]

        tree_results$.tree_rank <- suppressWarnings(
          as.numeric(
            tree_results[[tree_rank_column]]
          )
        )

      } else {

        tree_results$.tree_rank <- seq_len(
          nrow(tree_results)
        )
      }

      tree_results$.tree_rank[
        is.na(tree_results$.tree_rank) |
          !is.finite(tree_results$.tree_rank)
      ] <- Inf

      tree_results <- tree_results[
        order(
          tree_results$.tree_rank,
          na.last = TRUE
        ),
        ,
        drop = FALSE
      ]

      tree_results <- tree_results[
        !duplicated(tree_results$.tree_key),
        ,
        drop = FALSE
      ]

      tree_results <- head(
        tree_results,
        plot_terms
      )

      selected_tree_keys <- as.character(
        tree_results$.tree_key
      )

      selected_tree_keys <- selected_tree_keys[
        selected_tree_keys %in%
          rownames(tree_similarity) &
          selected_tree_keys %in%
          colnames(tree_similarity)
      ]

      if (length(selected_tree_keys) < 2) {
        stop(
          "Fewer than two terms were available after ",
          "matching the similarity matrix."
        )
      }

      tree_similarity <- tree_similarity[
        selected_tree_keys,
        selected_tree_keys,
        drop = FALSE
      ]

      storage.mode(tree_similarity) <- "numeric"

      tree_similarity[
        !is.finite(tree_similarity)
      ] <- 0

      # Make the matrix symmetrical and constrain similarity values.
      tree_similarity <- (
        tree_similarity +
          t(tree_similarity)
      ) / 2

      tree_similarity <- pmax(
        0,
        pmin(
          1,
          tree_similarity
        )
      )

      diag(tree_similarity) <- 1

      tree_distance <- stats::as.dist(
        1 - tree_similarity
      )

      tree_hclust <- stats::hclust(
        tree_distance,
        method = "average"
      )

      number_of_leaves <- length(
        tree_hclust$labels
      )

      # Position of each original leaf after hclust ordering.
      leaf_x <- numeric(number_of_leaves)

      leaf_x[
        tree_hclust$order
      ] <- seq_len(number_of_leaves)

      node_x <- numeric(
        number_of_leaves - 1L
      )

      tree_segments <- vector(
        "list",
        number_of_leaves - 1L
      )

      get_tree_node_position <- function(
        merge_value,
        current_row
      ) {

        if (merge_value < 0) {

          leaf_index <- -merge_value

          return(
            list(
              x = leaf_x[leaf_index],
              height = 0
            )
          )
        }

        if (merge_value >= current_row) {
          stop(
            "Invalid hclust merge structure."
          )
        }

        list(
          x = node_x[merge_value],
          height = tree_hclust$height[merge_value]
        )
      }

      for (
        merge_row in seq_len(
          nrow(tree_hclust$merge)
        )
      ) {

        left_node <- get_tree_node_position(
          tree_hclust$merge[merge_row, 1],
          merge_row
        )

        right_node <- get_tree_node_position(
          tree_hclust$merge[merge_row, 2],
          merge_row
        )

        merge_height <- tree_hclust$height[
          merge_row
        ]

        node_x[merge_row] <- mean(
          c(
            left_node$x,
            right_node$x
          )
        )

        tree_segments[[merge_row]] <- data.frame(
          x = c(
            left_node$x,
            right_node$x,
            left_node$x
          ),
          y = c(
            left_node$height,
            right_node$height,
            merge_height
          ),
          xend = c(
            left_node$x,
            right_node$x,
            right_node$x
          ),
          yend = c(
            merge_height,
            merge_height,
            merge_height
          ),
          stringsAsFactors = FALSE
        )
      }

      tree_segments <- do.call(
        rbind,
        tree_segments
      )

      ordered_tree_keys <- tree_hclust$labels[
        tree_hclust$order
      ]

      tree_label_table <- tree_results[
        match(
          ordered_tree_keys,
          tree_results$.tree_key
        ),
        ,
        drop = FALSE
      ]

      if ("Description" %in% names(tree_label_table)) {

        tree_labels <- as.character(
          tree_label_table$Description
        )

      } else {

        tree_labels <- ordered_tree_keys
      }

      tree_labels[
        is.na(tree_labels) |
          tree_labels == ""
      ] <- ordered_tree_keys[
        is.na(tree_labels) |
          tree_labels == ""
      ]

      tree_labels <- make.unique(
        tree_labels,
        sep = " | "
      )

      tree_labels <- stringr::str_wrap(
        tree_labels,
        width = 35
      )

      tree_leaf_data <- data.frame(
        x = seq_len(number_of_leaves),
        y = 0,
        key = ordered_tree_keys,
        label = tree_labels,
        stringsAsFactors = FALSE
      )

      tree_color_column <- intersect(
        c(
          "p.adjust",
          "pvalue",
          "qvalue"
        ),
        names(tree_label_table)
      )

      if (length(tree_color_column) > 0) {

        tree_color_column <- tree_color_column[1]

        tree_leaf_data$significance <- suppressWarnings(
          as.numeric(
            tree_label_table[[tree_color_column]]
          )
        )

        positive_tree_p <- tree_leaf_data$significance[
          !is.na(tree_leaf_data$significance) &
            tree_leaf_data$significance > 0
        ]

        minimum_tree_p <- if (
          length(positive_tree_p) > 0
        ) {
          min(
            positive_tree_p,
            na.rm = TRUE
          )
        } else {
          .Machine$double.xmin
        }

        tree_leaf_data$significance <- pmax(
          tree_leaf_data$significance,
          minimum_tree_p * 0.1
        )

        tree_leaf_data$minus_log10_significance <- -log10(
          tree_leaf_data$significance
        )

      } else {

        tree_leaf_data$minus_log10_significance <- 1
      }

      write.csv(
        tree_leaf_data,
        paste0(
          edb,
          "/GSEA_",
          edb0,
          "_tree_selected_terms.csv"
        ),
        row.names = FALSE
      )

      write.csv(
        tree_similarity,
        paste0(
          edb,
          "/GSEA_",
          edb0,
          "_tree_similarity_matrix.csv"
        ),
        row.names = TRUE
      )

      p.tree <- ggplot2::ggplot() +
        ggplot2::geom_segment(
          data = tree_segments,
          ggplot2::aes(
            x = x,
            y = y,
            xend = xend,
            yend = yend
          ),
          linewidth = 0.6,
          lineend = "round"
        ) +
        ggplot2::geom_point(
          data = tree_leaf_data,
          ggplot2::aes(
            x = x,
            y = y,
            color = minus_log10_significance
          ),
          size = 3
        ) +
        ggplot2::scale_x_continuous(
          breaks = tree_leaf_data$x,
          labels = tree_leaf_data$label,
          expand = ggplot2::expansion(
            mult = c(0.02, 0.02)
          )
        ) +
        ggplot2::scale_color_gradient(
          low = "blue",
          high = "red",
          name = expression(
            -log[10]("p-value")
          )
        ) +
        ggplot2::labs(
          x = NULL,
          y = "1 - term similarity",
          title = paste0(
            "Term-similarity tree: ",
            edb0
          ),
          subtitle = paste0(
            "Average-linkage clustering of ",
            number_of_leaves,
            " enriched terms"
          )
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(
            angle = 60,
            hjust = 1,
            vjust = 1,
            size = 8
          ),
          plot.title = ggplot2::element_text(
            hjust = 0.5,
            face = "bold"
          ),
          plot.subtitle = ggplot2::element_text(
            hjust = 0.5
          ),
          panel.grid.minor = ggplot2::element_blank()
        )

      tree_plot_width <- max(
        12,
        min(
          24,
          0.55 * number_of_leaves
        )
      )

      ggsave(
        paste0(
          edb,
          "/GSEA_",
          edb0,
          "_tree.pdf"
        ),
        plot = p.tree,
        width = tree_plot_width,
        height = 8,
        limitsize = FALSE
      )

      message(
        "[TREE] Robust term-similarity dendrogram saved: ",
        paste0(
          edb,
          "/GSEA_",
          edb0,
          "_tree.pdf"
        )
      )

    }, error = function(e) {
      write_plot_skip(
        "tree",
        e
      )
    })
  }

  ## ----------------------------------------------------------
  ## Enrichment map
  ## ----------------------------------------------------------

  if (!is.null(datax2)) {
    tryCatch({
      p.emap <- emapplot(
        datax2,
        layout = "kk",
        showCategory = plot_terms
      )

      if (nres < 25) {
        ggsave(
          paste0(edb, "/GSEA_", edb0, "_map.pdf"),
          plot = p.emap
        )
      } else {
        ggsave(
          paste0(edb, "/GSEA_", edb0, "_map.pdf"),
          plot = p.emap,
          scale = 2
        )
      }
    }, error = function(e) {
      write_plot_skip("map", e)
    })
  }

  ## ----------------------------------------------------------
  ## Heatplot
  ## ----------------------------------------------------------

  if (!is.null(datax2)) {
    tryCatch({
      p.heat <- heatplot(
        datax2,
        foldChange = genelist,
        showCategory = plot_terms
      )

      plot_width <- 8
      if ("core_enrichment" %in% colnames(data@result)) {
        plot_width <- max(
          8,
          0.2 * round(max(nchar(data@result$core_enrichment), na.rm = TRUE) / 19)
        )
      }

      if (nres < 30) {
        ggsave(
          paste0(edb, "/GSEA_", edb0, "_heat.pdf"),
          plot = p.heat,
          height = 2 + 0.24 * nres,
          width = plot_width,
          limitsize = FALSE
        )
      } else {
        ggsave(
          paste0(edb, "/GSEA_", edb0, "_heat.pdf"),
          plot = p.heat,
          height = 6.5,
          width = plot_width,
          limitsize = FALSE
        )
      }
    }, error = function(e) {
      write_plot_skip("heat", e)
    })
  }

  ## ----------------------------------------------------------
  ## Upsetplot
  ## ----------------------------------------------------------

  tryCatch({
    if (nres < 30) {
      pdf(
        file = paste0(edb, "/GSEA_", edb0, "_upset.pdf"),
        height = max(4, 0.4 * nres),
        width = 12
      )

      p.upset <- upsetplot(data, n = plot_terms)

    } else {
      pdf(
        file = paste0(edb, "/GSEA_", edb0, "_upset_top30.pdf"),
        height = 12,
        width = 16
      )

      p.upset <- upsetplot(data, n = plot_terms)
    }

    print(p.upset)
    dev.off()

  }, error = function(e) {
    if (grDevices::dev.cur() != 1) {
      grDevices::dev.off()
    }

    write_plot_skip("upset", e)
  })

  ## ----------------------------------------------------------
  ## Barplot
  ## ----------------------------------------------------------

  tryCatch({
    bar <- data@result

    bar <- bar[
      !is.na(bar$Description) &
        !is.na(bar$pvalue) &
        !is.na(bar$enrichmentScore),
      ,
      drop = FALSE
    ]

    if (nrow(bar) > 0) {
      bar$pvalue_for_plot <- bar$pvalue

      positive_p <- bar$pvalue_for_plot[bar$pvalue_for_plot > 0]
      if (length(positive_p) > 0) {
        min_positive <- min(positive_p, na.rm = TRUE)
        bar$pvalue_for_plot[bar$pvalue_for_plot <= 0] <- min_positive * 0.1
      }

      bar <- bar[
        order(
          bar$pvalue_for_plot,
          -abs(bar$enrichmentScore),
          na.last = TRUE
        ),
        ,
        drop = FALSE
      ]

      bar <- head(
        bar,
        plot_terms
      )

      bar$c <- with(bar, reorder(Description, -log10(pvalue_for_plot)))

      p.bar <- ggplot(
        bar,
        aes(
          x = -log10(pvalue_for_plot),
          y = c,
          fill = enrichmentScore
        )
      ) +
        geom_bar(stat = "identity") +
        scale_fill_gradient2(low = "blue", high = "red") +
        labs(
          x = "-log10(p-value)",
          y = "Description",
          fill = "enrichmentScore"
        ) +
        ggtitle(paste0(group1, "_vs_", group2)) +
        theme_bw() +
        theme(text = element_text(size = 18)) +
        scale_y_discrete(
          breaks = bar$Description,
          labels = stringr::str_trunc(bar$Description, 40)
        )

      ggsave(
        paste0(edb, "/GSEA_", edb0, "_barplots.pdf"),
        plot = p.bar,
        width = 10,
        height = max(4, nrow(bar) / 12 * 5),
        limitsize = FALSE
      )
    }
  }, error = function(e) {
    write_plot_skip("barplot", e)
  })

  ## ----------------------------------------------------------
  ## Full result handling
  ## ----------------------------------------------------------

  if (nres > max_plot_terms) {
    writeLines(
      c(
        paste0(
          "The enrichment result contains ",
          nres,
          " terms."
        ),
        paste0(
          "Figures were limited to the top ",
          max_plot_terms,
          " terms to prevent oversized or unreadable files."
        ),
        "The complete results remain available in the CSV tables."
      ),
      paste0(
        edb,
        "/GSEA_",
        edb0,
        "_plot_limit_note.txt"
      )
    )
  }

  invisible(NULL)
}
#### BEGIN FUNCTION: pathview.p ####
# function for KEGG Pathview plots
pathview.p <- function(
  kk,
  ko.db,
  kegg_gene_list,
  dir = NULL,
  max_pathways = suppressWarnings(
    as.integer(
      Sys.getenv(
        "MTD_PATHVIEW_MAX_PATHWAYS",
        "20"
      )
    )
  )
) {
  if (is.null(kk) || nrow(kk@result) == 0) {
    write("No enrichment result was found on KEGG.", "No_KEGG_enrichment_result.txt")
    return(invisible(NULL))
  }

  ko.db <- trimws(as.character(ko.db))

  suppressPackageStartupMessages(library("pathview"))

  patch_pathview_mlf <- function() {
    old_fun <- pathview::kegg.species.code

    my_kegg_species_code <- function(species = "hsa", na.rm = FALSE, code.only = TRUE) {
      species_chr <- as.character(species)

      is_mlf <- species_chr %in% c(
        "mlf",
        "Myotis lucifugus",
        "59463"
      )

      if (all(is_mlf)) {
        if (code.only) {
          return(rep("mlf", length(species_chr)))
        } else {
          out <- matrix(
            rep(
              c(
                "T07795",
                "59463",
                "mlf",
                "Myotis lucifugus",
                "Little brown bat",
                "1",
                "102426605",
                "102426605",
                NA,
                NA
              ),
              length(species_chr)
            ),
            nrow = length(species_chr),
            byrow = TRUE
          )

          colnames(out) <- c(
            "ktax.id",
            "tax.id",
            "kegg.code",
            "scientific.name",
            "common.name",
            "entrez.gnodes",
            "kegg.geneid",
            "ncbi.geneid",
            "ncbi.proteinid",
            "uniprot"
          )

          return(out)
        }
      }

      old_fun(
        species = species,
        na.rm = na.rm,
        code.only = code.only
      )
    }

    unlockBinding("kegg.species.code", asNamespace("pathview"))
    assign(
      "kegg.species.code",
      my_kegg_species_code,
      envir = asNamespace("pathview")
    )
    lockBinding("kegg.species.code", asNamespace("pathview"))

    message("Pathview patched for Myotis lucifugus / mlf")
  }

  if (ko.db == "mlf") {
    patch_pathview_mlf()
  }

  gene_data <- kegg_gene_list
  gene_data <- gene_data[!is.na(gene_data)]
  gene_data <- gene_data[!is.na(names(gene_data))]
  gene_data <- gene_data[names(gene_data) != ""]

  # Para pathview com gene.idtype = "KEGG",
  # usar IDs numéricos KEGG/Entrez sem prefixo "mlf:".
  names(gene_data) <- sub("^.*:", "", names(gene_data))

  #### BEGIN BLOCK: bounded Pathview pathways ####

  if (
    is.na(max_pathways) ||
    max_pathways < 1
  ) {
    max_pathways <- 20L
  }

  pathview_results <- as.data.frame(
    kk@result,
    stringsAsFactors = FALSE
  )

  if ("p.adjust" %in% names(pathview_results)) {

    pathview_order <- order(
      suppressWarnings(
        as.numeric(pathview_results$p.adjust)
      ),
      na.last = TRUE
    )

  } else if ("pvalue" %in% names(pathview_results)) {

    pathview_order <- order(
      suppressWarnings(
        as.numeric(pathview_results$pvalue)
      ),
      na.last = TRUE
    )

  } else {

    pathview_order <- seq_len(
      nrow(pathview_results)
    )
  }

  pathview_results <- pathview_results[
    pathview_order,
    ,
    drop = FALSE
  ]

  total_pathview_results <- nrow(
    pathview_results
  )

  pathview_results <- head(
    pathview_results,
    max_pathways
  )

  write.csv(
    pathview_results,
    "Pathview_selected_pathways.csv",
    row.names = FALSE
  )

  message(
    "[PATHVIEW] Enriched pathways available: ",
    total_pathview_results,
    "; pathways selected: ",
    nrow(pathview_results),
    " (maximum ",
    max_pathways,
    ")"
  )

  for (g in seq_len(nrow(pathview_results))) {
    pathway_id <- as.character(pathview_results[g, 1])
    pathway_id <- sub("^path:", "", pathway_id)
    pathway_id <- sub(paste0("^", ko.db), "", pathway_id)

    message("Running pathview for pathway: ", pathway_id, " species: ", ko.db)

    tryCatch({
      pathview::pathview(
        gene.data = gene_data,
        pathway.id = pathway_id,
        species = ko.db,
        gene.idtype = "KEGG",
        limit = list(
          gene = max(abs(gene_data), na.rm = TRUE),
          cpd = 1
        ),
        kegg.native = TRUE
      )
    }, error = function(e) {
      message(
        "Pathview failed for pathway ",
        pathway_id,
        ": ",
        conditionMessage(e)
      )

      write(
        paste0(
          "Pathview failed for pathway ",
          pathway_id,
          ": ",
          conditionMessage(e)
        ),
        file = paste0("Pathview_failed_", pathway_id, ".txt")
      )
    })
  }

  invisible(NULL)
}
#### END FUNCTION: pathview.p ####

#### BEGIN FUNCTION: pathview.modules.p ####
# Convert KEGG Module IDs such as M00087 into linked KEGG pathway IDs,
# then run pathview on those pathways.
pathview.modules.p <- function(kk_module, ko.db, kegg_gene_list, max_modules = 20) {
  if (is.null(kk_module) || nrow(kk_module@result) == 0) {
    write(
      "No KEGG module enrichment result was found.",
      "No_KEGG_module_enrichment_result.txt"
    )
    return(invisible(NULL))
  }

  ko.db <- trimws(as.character(ko.db))

  suppressPackageStartupMessages(library("KEGGREST"))

  module_ids <- as.character(kk_module@result$ID)
  module_ids <- module_ids[grepl("^M[0-9]{5}$", module_ids)]

  if (length(module_ids) == 0) {
    write(
      "No valid KEGG Module IDs were found in kk_module@result.",
      "No_valid_KEGG_module_IDs_for_pathview.txt"
    )
    return(invisible(NULL))
  }

  module_ids <- head(unique(module_ids), max_modules)

  module_to_pathway <- list()

  for (mid in module_ids) {
    message("Resolving KEGG module to pathway: ", mid)

    kg <- tryCatch({
      KEGGREST::keggGet(mid)
    }, error = function(e) {
      message("KEGGREST::keggGet failed for ", mid, ": ", conditionMessage(e))
      NULL
    })

    if (is.null(kg) || length(kg) == 0) {
      next
    }

    pathways <- kg[[1]]$PATHWAY

    if (is.null(pathways) || length(pathways) == 0) {
      message("No linked pathways found for module: ", mid)
      next
    }

    pw_ids <- names(pathways)

    # Usually returns IDs like "map00071" or similar.
    # pathview wants the numeric pathway ID, e.g. "00071",
    # together with species = "mlf".
    pw_ids <- sub("^[a-z]{2,4}", "", pw_ids)
    pw_ids <- pw_ids[grepl("^[0-9]{5}$", pw_ids)]

    if (length(pw_ids) > 0) {
      module_to_pathway[[mid]] <- unique(pw_ids)
    }
  }

  if (length(module_to_pathway) == 0) {
    write(
      "No KEGG pathways could be resolved from enriched KEGG modules.",
      "No_pathways_resolved_from_KEGG_modules.txt"
    )
    return(invisible(NULL))
  }

  resolved_table <- data.frame(
    module_id = rep(names(module_to_pathway), lengths(module_to_pathway)),
    pathway_id = unlist(module_to_pathway),
    stringsAsFactors = FALSE
  )

  write.csv(
    resolved_table,
    "KEGG_modules_resolved_to_pathways.csv",
    row.names = FALSE
  )

  pathway_ids <- unique(resolved_table$pathway_id)

  message("Pathways resolved from KEGG modules:")
  print(pathway_ids)

  # Create a minimal object compatible with pathview.p:
  # pathview.p only needs kk@result and the first column/ID.
  fake_kk <- kk_module
  fake_kk@result <- data.frame(
    ID = pathway_ids,
    Description = paste0("Pathway resolved from KEGG module: ", pathway_ids),
    stringsAsFactors = FALSE
  )

  pathview.p(fake_kk, ko.db, kegg_gene_list)

  invisible(resolved_table)
}
#### END FUNCTION: pathview.modules.p ####
#### BEGIN FUNCTION: gseKEGG.safe ####
gseKEGG.safe <- function(geneList, ko.db, minGSSize = 8, maxGSSize = 500, pvalueCutoff = 0.05) {
  ko.db <- trimws(as.character(ko.db))

  geneList <- geneList[!is.na(geneList)]
  geneList <- geneList[!is.na(names(geneList))]
  geneList <- geneList[names(geneList) != ""]
  geneList <- geneList[!duplicated(names(geneList))]
  geneList <- sort(geneList, decreasing = TRUE)

  message("Trying gseKEGG with keyType = ncbi-geneid")

  kk <- tryCatch({
    clusterProfiler::gseKEGG(
      geneList = geneList,
      organism = ko.db,
      keyType = "ncbi-geneid",
      minGSSize = minGSSize,
      maxGSSize = maxGSSize,
      pvalueCutoff = pvalueCutoff,
      verbose = FALSE
    )
  }, error = function(e) {
    message("gseKEGG ncbi-geneid failed: ", conditionMessage(e))
    NULL
  })

  if (!is.null(kk) && nrow(kk@result) > 0) {
    message("gseKEGG ncbi-geneid worked.")
    return(kk)
  }

  message("Trying gseKEGG with keyType = kegg")

  kk <- tryCatch({
    clusterProfiler::gseKEGG(
      geneList = geneList,
      organism = ko.db,
      keyType = "kegg",
      minGSSize = minGSSize,
      maxGSSize = maxGSSize,
      pvalueCutoff = pvalueCutoff,
      verbose = FALSE
    )
  }, error = function(e) {
    message("gseKEGG kegg failed: ", conditionMessage(e))
    NULL
  })

  if (!is.null(kk) && nrow(kk@result) > 0) {
    message("gseKEGG kegg worked.")
    return(kk)
  }

  message("Falling back to KEGGREST::keggLink + clusterProfiler::GSEA")

  suppressPackageStartupMessages(library("KEGGREST"))

  links <- tryCatch({
    KEGGREST::keggLink("pathway", ko.db)
  }, error = function(e) {
    message("KEGGREST::keggLink failed: ", conditionMessage(e))
    NULL
  })

  if (is.null(links) || length(links) == 0) {
    message("No KEGG pathway links found for organism: ", ko.db)
    return(NULL)
  }

  df_links <- data.frame(
    from = names(links),
    to = unname(links),
    stringsAsFactors = FALSE
  )

  if (any(grepl(paste0("^", ko.db, ":"), df_links$from))) {
    gene_col <- df_links$from
    path_col <- df_links$to
  } else {
    gene_col <- df_links$to
    path_col <- df_links$from
  }

  TERM2GENE <- data.frame(
    term = sub("^path:", "", path_col),
    gene = sub(paste0("^", ko.db, ":"), "", gene_col),
    stringsAsFactors = FALSE
  )

  TERM2GENE <- TERM2GENE[
    grepl(paste0("^", ko.db, "[0-9]{5}$"), TERM2GENE$term) &
    TERM2GENE$gene != "",
  ]

  TERM2GENE <- unique(TERM2GENE)

  overlap <- length(intersect(names(geneList), TERM2GENE$gene))

  message("KEGG TERM2GENE rows: ", nrow(TERM2GENE))
  message("Unique KEGG pathway genes: ", length(unique(TERM2GENE$gene)))
  message("geneList genes overlapping KEGG pathways: ", overlap)

  if (overlap == 0) {
    message("No overlap between geneList and KEGG pathway genes.")
    return(NULL)
  }

  path_names <- tryCatch({
    KEGGREST::keggList("pathway", ko.db)
  }, error = function(e) {
    message("KEGGREST::keggList pathway failed: ", conditionMessage(e))
    NULL
  })

  if (is.null(path_names) || length(path_names) == 0) {
    TERM2NAME <- NULL
  } else {
    TERM2NAME <- data.frame(
      term = sub("^path:", "", names(path_names)),
      name = as.character(path_names),
      stringsAsFactors = FALSE
    )

    TERM2NAME$name <- sub(" - .*$", "", TERM2NAME$name)
  }

  kk <- tryCatch({
    clusterProfiler::GSEA(
      geneList = geneList,
      TERM2GENE = TERM2GENE,
      TERM2NAME = TERM2NAME,
      minGSSize = minGSSize,
      maxGSSize = maxGSSize,
      pvalueCutoff = pvalueCutoff,
      verbose = FALSE
    )
  }, error = function(e) {
    message("Manual KEGG pathway GSEA failed: ", conditionMessage(e))
    NULL
  })

  if (!is.null(kk)) {
    message("Manual KEGG pathway GSEA returned terms: ", nrow(kk@result))
  }

  return(kk)
}
#### END FUNCTION: gseKEGG.safe ####

#### BEGIN FUNCTION: dual KEGG helpers ####

dedup_ranked_gene_list <- function(x) {
  x <- x[!is.na(x)]
  x <- x[!is.na(names(x))]
  x <- x[names(x) != ""]
  x <- x[names(x) != "-"]

  if (length(x) == 0) {
    return(x)
  }

  ord <- order(abs(x), decreasing = TRUE)
  x <- x[ord]
  x <- x[!duplicated(names(x))]
  x <- sort(x, decreasing = TRUE)

  x
}

get_host_kegg_code <- function(host_sp, taxid) {
  ko.db <- trimws(as.character(
    host_sp[as.character(host_sp$Taxon_ID) == as.character(taxid), 6][1]
  ))

#  if ((is.na(ko.db) || ko.db == "") && as.character(taxid) == "59463") {
#    ko.db <- "mlf"
#  }

  ko.db
}

#### BEGIN FUNCTION: persistent official KEGG mapping ####

normalize_host_gene_id <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("\\.[0-9]+$", "", x)
  x
}

empty_gid_to_entrez_map <- function() {
  data.frame(
    ensembl_gene_id = character(),
    external_gene_name = character(),
    entrezgene_id = character(),
    stringsAsFactors = FALSE
  )
}

read_gid_to_entrez_map <- function(mapping_file, gids = NULL, source_label = "mapping") {
  if (is.na(mapping_file) || mapping_file == "" || !file.exists(mapping_file)) {
    return(empty_gid_to_entrez_map())
  }

  x <- tryCatch({
    read.delim(
      mapping_file,
      header = TRUE,
      sep = "\t",
      quote = "",
      comment.char = "",
      check.names = FALSE,
      colClasses = "character",
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    message(
      "[WARNING] Could not read ", source_label, " GID -> ENTREZ mapping: ",
      conditionMessage(e)
    )
    empty_gid_to_entrez_map()
  })

  required_columns <- c("ensembl_gene_id", "entrezgene_id")
  missing_columns <- setdiff(required_columns, colnames(x))

  if (length(missing_columns) > 0) {
    message(
      "[WARNING] Ignoring invalid ", source_label, " mapping: ", mapping_file,
      ". Missing columns: ", paste(missing_columns, collapse = ", ")
    )
    return(empty_gid_to_entrez_map())
  }

  if (!"external_gene_name" %in% colnames(x)) {
    x$external_gene_name <- ""
  }

  x <- x[, c(
    "ensembl_gene_id",
    "external_gene_name",
    "entrezgene_id"
  ), drop = FALSE]

  x$ensembl_gene_id <- normalize_host_gene_id(x$ensembl_gene_id)
  x$external_gene_name <- as.character(x$external_gene_name)
  x$entrezgene_id <- trimws(as.character(x$entrezgene_id))

  x <- x[
    !is.na(x$ensembl_gene_id) &
      x$ensembl_gene_id != "" &
      x$ensembl_gene_id != "-" &
      !is.na(x$entrezgene_id) &
      grepl("^[0-9]+$", x$entrezgene_id),
    ,
    drop = FALSE
  ]

  x <- unique(x)

  if (!is.null(gids)) {
    gids <- unique(normalize_host_gene_id(gids))
    gids <- gids[!is.na(gids) & gids != "" & gids != "-"]
    x <- x[x$ensembl_gene_id %in% gids, , drop = FALSE]
  }

  attr(x, "mapping_source") <- source_label
  attr(x, "mapping_file") <- normalizePath(mapping_file, mustWork = FALSE)

  message(
    "[INFO] Valid ", source_label, " GID -> ENTREZ mapping rows: ", nrow(x)
  )
  message("[INFO] Mapping file: ", normalizePath(mapping_file, mustWork = FALSE))

  x
}

write_gid_to_entrez_map <- function(x, out_file) {
  if (nrow(x) == 0) {
    return(invisible(FALSE))
  }

  dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)

  write.table(
    x[, c(
      "ensembl_gene_id",
      "external_gene_name",
      "entrezgene_id"
    ), drop = FALSE],
    file = out_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  invisible(TRUE)
}

build_gid_to_entrez_map <- function(gids, host_sp, taxid, out_file) {
  gids <- unique(normalize_host_gene_id(gids))
  gids <- gids[!is.na(gids) & gids != "" & gids != "-"]

  if (length(gids) == 0) {
    return(empty_gid_to_entrez_map())
  }

  # 1) Reuse the run-specific mapping first. This preserves exact rerun
  # reproducibility even if the host reference is later rebuilt.
  cached_map <- read_gid_to_entrez_map(
    mapping_file = out_file,
    gids = gids,
    source_label = "analysis_cache"
  )

  if (nrow(cached_map) > 0) {
    return(cached_map)
  }

  # 2) Use the persistent map generated by Create_custom_host.sh.
  mtd_root <- dirname(normalizePath(args[4], mustWork = FALSE))

  persistent_mapping_file <- file.path(
    mtd_root,
    paste0("ref_", taxid),
    "functional_annotation",
    paste0("GID_to_ENTREZ_taxid_", taxid, ".tsv")
  )

  persistent_map <- read_gid_to_entrez_map(
    mapping_file = persistent_mapping_file,
    gids = gids,
    source_label = "persistent_reference"
  )

  if (nrow(persistent_map) > 0) {
    write_gid_to_entrez_map(persistent_map, out_file)
    message("[INFO] Saved run-specific GID -> ENTREZ mapping: ", out_file)
    return(persistent_map)
  }

  # 3) BioMart remains a backward-compatible online fallback for old
  # references that were created before the persistent map feature.
  dataset_use <- as.character(
    host_sp[as.character(host_sp$Taxon_ID) == as.character(taxid), 2][1]
  )

  if (is.na(dataset_use) || dataset_use == "") {
    warning(
      "[WARNING] No persistent map or Ensembl/BioMart dataset found for TaxID ",
      taxid
    )
    return(empty_gid_to_entrez_map())
  }

  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    warning(
      "[WARNING] biomaRt is not installed and no persistent GID -> ENTREZ ",
      "mapping was available."
    )
    return(empty_gid_to_entrez_map())
  }

  message(
    "[INFO] Persistent GID -> ENTREZ mapping unavailable. ",
    "Trying BioMart fallback dataset: ", dataset_use
  )

  mart <- tryCatch({
    biomaRt::useEnsembl(
      biomart = "genes",
      dataset = dataset_use,
      mirror = "www"
    )
  }, error = function(e) {
    message("[WARNING] BioMart connection failed: ", conditionMessage(e))
    NULL
  })

  if (is.null(mart)) {
    return(empty_gid_to_entrez_map())
  }

  x <- tryCatch({
    biomaRt::getBM(
      attributes = c(
        "ensembl_gene_id",
        "external_gene_name",
        "entrezgene_id"
      ),
      filters = "ensembl_gene_id",
      values = gids,
      mart = mart
    )
  }, error = function(e) {
    message("[WARNING] BioMart getBM failed: ", conditionMessage(e))
    empty_gid_to_entrez_map()
  })

  if (nrow(x) == 0) {
    return(empty_gid_to_entrez_map())
  }

  temporary_mapping_file <- tempfile(fileext = ".tsv")
  write.table(
    x,
    file = temporary_mapping_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  x <- read_gid_to_entrez_map(
    mapping_file = temporary_mapping_file,
    gids = gids,
    source_label = "biomart_fallback"
  )
  unlink(temporary_mapping_file)

  if (nrow(x) > 0) {
    write_gid_to_entrez_map(x, out_file)
    message("[INFO] Saved BioMart fallback GID -> ENTREZ mapping: ", out_file)
  }

  x
}

prepare_official_kegg_gene_list <- function(
  genelist_gid_or_entrez,
  host_sp,
  taxid,
  do.db,
  host_gene_keytype,
  outdir
) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  if (length(genelist_gid_or_entrez) == 0) {
    return(genelist_gid_or_entrez)
  }

  # Case 1: geneList is already ENTREZ.
  if (host_gene_keytype == "ENTREZID") {
    out <- dedup_ranked_gene_list(genelist_gid_or_entrez)

    write.table(
      data.frame(ENTREZID = names(out), ranking_value = as.numeric(out)),
      file.path(outdir, "official_geneList_ENTREZ.tsv"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )

    return(out)
  }

  # Case 2: custom OrgDb uses GID/Ensembl IDs.
  mapping_file <- file.path(
    dirname(args[1]),
    "Host_DEG",
    paste0("GID_to_ENTREZ_taxid_", taxid, ".tsv")
  )

  map <- build_gid_to_entrez_map(
    gids = names(genelist_gid_or_entrez),
    host_sp = host_sp,
    taxid = taxid,
    out_file = mapping_file
  )

  if (nrow(map) == 0) {
    write(
      paste0(
        "Official KEGG was skipped because GID -> ENTREZ mapping was unavailable.\n",
        "TaxID: ", taxid, "\n",
        "Host gene keytype: ", host_gene_keytype, "\n",
        "Expected persistent map: ",
        file.path(
          dirname(normalizePath(args[4], mustWork = FALSE)),
          paste0("ref_", taxid),
          "functional_annotation",
          paste0("GID_to_ENTREZ_taxid_", taxid, ".tsv")
        ),
        "\n"
      ),
      file.path(outdir, "No_GID_to_ENTREZ_mapping_for_official_KEGG.txt")
    )

    return(numeric(0))
  }

  mapping_source <- attr(map, "mapping_source")
  mapping_source_file <- attr(map, "mapping_file")

  gid <- normalize_host_gene_id(names(genelist_gid_or_entrez))

  rank_df <- data.frame(
    GID = gid,
    ranking_value = as.numeric(genelist_gid_or_entrez),
    stringsAsFactors = FALSE
  )

  rank_df <- rank_df[
    !is.na(rank_df$GID) &
      rank_df$GID != "" &
      rank_df$GID != "-" &
      !is.na(rank_df$ranking_value) &
      is.finite(rank_df$ranking_value),
    ,
    drop = FALSE
  ]

  map_keep <- map[, c(
    "ensembl_gene_id",
    "external_gene_name",
    "entrezgene_id"
  ), drop = FALSE]

  mapped_df <- merge(
    rank_df,
    map_keep,
    by.x = "GID",
    by.y = "ensembl_gene_id",
    all.x = FALSE,
    all.y = FALSE
  )

  mapped_df <- mapped_df[
    !is.na(mapped_df$entrezgene_id) &
      grepl("^[0-9]+$", mapped_df$entrezgene_id) &
      !is.na(mapped_df$ranking_value) &
      is.finite(mapped_df$ranking_value),
    ,
    drop = FALSE
  ]

  mapped_df$entrezgene_id <- as.character(mapped_df$entrezgene_id)

  # Multiple Ensembl IDs can converge on one NCBI GeneID. Keep the entry
  # with the strongest absolute ranking value for GSEA.
  mapped_df <- mapped_df[
    order(abs(mapped_df$ranking_value), decreasing = TRUE),
    ,
    drop = FALSE
  ]

  mapped_unique <- mapped_df[
    !duplicated(mapped_df$entrezgene_id),
    ,
    drop = FALSE
  ]

  out <- mapped_unique$ranking_value
  names(out) <- mapped_unique$entrezgene_id
  out <- sort(out, decreasing = TRUE)

  out_table <- mapped_unique
  names(out_table)[names(out_table) == "entrezgene_id"] <- "ENTREZID"

  write.table(
    out_table,
    file.path(outdir, "official_GID_to_ENTREZ_geneList_used.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  input_unique_gids <- length(unique(rank_df$GID))
  mapped_unique_gids <- length(unique(mapped_df$GID))
  mapping_coverage <- if (input_unique_gids > 0) {
    100 * mapped_unique_gids / input_unique_gids
  } else {
    0
  }

  write(
    paste0(
      "Mapping source: ", ifelse(is.null(mapping_source), "unknown", mapping_source), "\n",
      "Mapping file: ", ifelse(is.null(mapping_source_file), "unknown", mapping_source_file), "\n",
      "Input ranked genes: ", length(genelist_gid_or_entrez), "\n",
      "Unique normalized input GIDs: ", input_unique_gids, "\n",
      "Mapped unique GIDs: ", mapped_unique_gids, "\n",
      "GID mapping coverage: ", sprintf("%.2f", mapping_coverage), "%\n",
      "Rows after GID -> ENTREZ merge: ", nrow(mapped_df), "\n",
      "Unique ENTREZ genes used for official KEGG: ", length(out), "\n",
      "Duplicated ENTREZ IDs removed: ", nrow(mapped_df) - length(out), "\n"
    ),
    file.path(outdir, "official_KEGG_mapping_summary.txt")
  )

  out
}

convert_gene_clusters_to_entrez <- function(
  gene_clusters,
  host_sp,
  taxid,
  host_gene_keytype,
  outdir
) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  clean_cluster <- function(ids) {
    ids <- unique(trimws(as.character(ids)))
    ids[!is.na(ids) & ids != "" & ids != "-"]
  }

  gene_clusters <- lapply(gene_clusters, clean_cluster)
  gene_clusters <- gene_clusters[lengths(gene_clusters) > 0]

  if (length(gene_clusters) == 0) {
    return(list())
  }

  if (host_gene_keytype == "ENTREZID") {
    gene_clusters <- lapply(
      gene_clusters,
      function(ids) ids[grepl("^[0-9]+$", ids)]
    )
    return(gene_clusters[lengths(gene_clusters) > 0])
  }

  all_gids <- unique(normalize_host_gene_id(unlist(gene_clusters, use.names = FALSE)))

  mapping_file <- file.path(
    dirname(args[1]),
    "Host_DEG",
    paste0("GID_to_ENTREZ_taxid_", taxid, ".tsv")
  )

  map <- build_gid_to_entrez_map(
    gids = all_gids,
    host_sp = host_sp,
    taxid = taxid,
    out_file = mapping_file
  )

  if (nrow(map) == 0) {
    return(list())
  }

  converted <- lapply(gene_clusters, function(ids) {
    ids <- normalize_host_gene_id(ids)
    entrez <- unique(map$entrezgene_id[map$ensembl_gene_id %in% ids])
    entrez <- entrez[
      !is.na(entrez) &
        entrez != "" &
        grepl("^[0-9]+$", entrez)
    ]
    as.character(entrez)
  })

  converted <- converted[lengths(converted) > 0]

  cluster_summary <- data.frame(
    Cluster = names(gene_clusters),
    Input_GIDs = lengths(gene_clusters),
    Mapped_ENTREZ = vapply(
      names(gene_clusters),
      function(cluster_name) {
        if (cluster_name %in% names(converted)) {
          length(converted[[cluster_name]])
        } else {
          0L
        }
      },
      integer(1)
    ),
    stringsAsFactors = FALSE
  )

  write.table(
    cluster_summary,
    file.path(outdir, "official_KEGG_compareCluster_mapping_summary.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  write.table(
    map,
    file.path(outdir, "official_KEGG_compareCluster_GID_to_ENTREZ.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  converted
}

#### END FUNCTION: persistent official KEGG mapping ####

make_kegg_derived_term2gene <- function(
  orgdb,
  column,
  prefix_regex = NULL,
  out_file = NULL
) {
  if (!column %in% AnnotationDbi::columns(orgdb)) {
    return(data.frame(term = character(), gene = character()))
  }

  all_gid <- AnnotationDbi::keys(orgdb, keytype = "GID")

  x <- AnnotationDbi::select(
    orgdb,
    keys = all_gid,
    keytype = "GID",
    columns = c("GID", column)
  )

  x <- x[
    !is.na(x[[column]]) &
      x[[column]] != "" &
      !is.na(x$GID) &
      x$GID != "",
    c(column, "GID"),
    drop = FALSE
  ]

  x <- unique(x)

  if (!is.null(prefix_regex)) {
    x <- x[grepl(prefix_regex, x[[column]]), , drop = FALSE]
  }

  term2gene <- data.frame(
    term = as.character(x[[column]]),
    gene = as.character(x$GID),
    stringsAsFactors = FALSE
  )

  term2gene <- unique(term2gene)

  if (!is.null(out_file)) {
    write.table(
      term2gene,
      file = out_file,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }

  term2gene
}

run_kegg_derived_gsea_and_enricher <- function(
  genelist_gid,
  sig_genes_gid,
  orgdb,
  column,
  prefix_regex,
  label,
  outdir,
  group1,
  group2
) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  term2gene_file <- file.path(outdir, paste0("TERM2GENE_", label, ".tsv"))

  term2gene <- make_kegg_derived_term2gene(
    orgdb = orgdb,
    column = column,
    prefix_regex = prefix_regex,
    out_file = term2gene_file
  )

  overlap <- length(intersect(names(genelist_gid), term2gene$gene))

  write(
    paste0(
      "Column: ", column, "\n",
      "Label: ", label, "\n",
      "TERM2GENE rows: ", nrow(term2gene), "\n",
      "Unique terms: ", length(unique(term2gene$term)), "\n",
      "Unique genes in TERM2GENE: ", length(unique(term2gene$gene)), "\n",
      "Ranked genes: ", length(genelist_gid), "\n",
      "Overlap ranked genes vs TERM2GENE: ", overlap, "\n",
      "Significant genes for enricher: ", length(sig_genes_gid), "\n"
    ),
    file.path(outdir, paste0(label, "_input_summary.txt"))
  )

  if (nrow(term2gene) == 0 || overlap < 10) {
    write(
      paste0(
        "No enough overlap for eggNOG-derived KEGG ", label, ".\n",
        "Overlap: ", overlap, "\n"
      ),
      file.path(outdir, paste0("No_", label, "_GSEA_result.txt"))
    )
    return(invisible(NULL))
  }

  term2name <- data.frame(
    term = unique(term2gene$term),
    name = unique(term2gene$term),
    stringsAsFactors = FALSE
  )

  kk_gsea <- tryCatch({
    clusterProfiler::GSEA(
      geneList = genelist_gid,
      TERM2GENE = term2gene,
      TERM2NAME = term2name,
      minGSSize = 10,
      maxGSSize = 500,
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      verbose = FALSE
    )
  }, error = function(e) {
    message("eggNOG-derived KEGG GSEA failed for ", label, ": ", conditionMessage(e))
    NULL
  })

  if (is.null(kk_gsea) || nrow(kk_gsea@result) == 0) {
    write(
      paste0("No eggNOG-derived KEGG GSEA result for ", label),
      file.path(outdir, paste0("No_", label, "_GSEA_result.txt"))
    )
  } else {
    write.csv(
      kk_gsea@result,
      file.path(outdir, paste0("GSEA_", label, "_results.csv")),
      row.names = FALSE
    )

    try(
      plots4gsea(
        edb = outdir,
        data = kk_gsea,
        datax = kk_gsea,
        edb0 = label,
        genelist = genelist_gid,
        group1 = group1,
        group2 = group2
      ),
      silent = TRUE
    )
  }

  sig_genes_gid <- unique(sig_genes_gid)
  sig_genes_gid <- sig_genes_gid[
    !is.na(sig_genes_gid) &
      sig_genes_gid != "" &
      sig_genes_gid != "-"
  ]

  sig_genes_gid <- intersect(sig_genes_gid, unique(term2gene$gene))

  if (length(sig_genes_gid) >= 5) {
    kk_enrich <- tryCatch({
      clusterProfiler::enricher(
        gene = sig_genes_gid,
        TERM2GENE = term2gene,
        TERM2NAME = term2name,
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        qvalueCutoff = 0.2,
        minGSSize = 10,
        maxGSSize = 500
      )
    }, error = function(e) {
      message("eggNOG-derived KEGG enricher failed for ", label, ": ", conditionMessage(e))
      NULL
    })

    if (!is.null(kk_enrich)) {
      write.csv(
        as.data.frame(kk_enrich),
        file.path(outdir, paste0("enricher_", label, "_results.csv")),
        row.names = FALSE
      )
    }
  } else {
    write(
      paste0("Too few significant genes for enricher: ", length(sig_genes_gid)),
      file.path(outdir, paste0("No_", label, "_enricher_result.txt"))
    )
  }

  invisible(kk_gsea)
}

#### END FUNCTION: dual KEGG helpers ####

    # function for pathway enrichment by using comparison results between groups
    enrichment <- function(coldata_vs,args1,do.db){
      for (i in 1:nrow(coldata_vs)){
        group1<-coldata_vs$group1[i]
        group2<-coldata_vs$group2[i]
        setwd(paste0(group1,"_vs_",group2)) # go into each comparison folder
dir.create("GO", showWarnings = FALSE, recursive = TRUE)

dir.create("KEGG", showWarnings = FALSE, recursive = TRUE)

dir.create("KEGG/official_gseKEGG", showWarnings = FALSE, recursive = TRUE)
dir.create("KEGG/official_gseKEGG/Modules", showWarnings = FALSE, recursive = TRUE)
dir.create("KEGG/official_gseKEGG/Pathview", showWarnings = FALSE, recursive = TRUE)
dir.create("KEGG/official_gseKEGG/Modules/Pathview", showWarnings = FALSE, recursive = TRUE)

dir.create("KEGG/eggNOG_KEGG_derived", showWarnings = FALSE, recursive = TRUE)
dir.create("KEGG/eggNOG_KEGG_derived/Pathway", showWarnings = FALSE, recursive = TRUE)
dir.create("KEGG/eggNOG_KEGG_derived/Module", showWarnings = FALSE, recursive = TRUE)
dir.create("KEGG/eggNOG_KEGG_derived/KO", showWarnings = FALSE, recursive = TRUE)
dir.create("KEGG/eggNOG_KEGG_derived/Reaction", showWarnings = FALSE, recursive = TRUE)
        # Prepare a complete ranked list for GSEA. Do not pre-filter by
        # significance: ORA/enricher uses the significant subset below.
        df.compare <- read.csv(
          paste0("host_counts_", group1, "_vs_", group2, ".csv"),
          header = TRUE,
          check.names = FALSE
        )

        #### BEGIN BLOCK: comparison gene-ID column normalization ####
        # write.csv() stores row names in the first column with an empty
        # header. With check.names = FALSE, that header remains empty rather
        # than being automatically converted to "X".
        #
        # Accept common gene-ID column names and normalize the selected column
        # to "X" for compatibility with the existing downstream code.
        if (!"X" %in% colnames(df.compare)) {

          comparison_id_aliases <- c(
            "",
            "...1",
            "row.names",
            "Row.names",
            "gene_id",
            "GeneID",
            "ensembl_gene_id",
            "GID"
          )

          comparison_id_candidates <- which(
            colnames(df.compare) %in% comparison_id_aliases
          )

          # The comparison CSV generated by this script always places the
          # row-name gene IDs in the first column. Use that position as a
          # guarded fallback when its header has an unexpected representation.
          if (
            length(comparison_id_candidates) == 0 &&
            ncol(df.compare) > 0
          ) {
            first_column_values <- as.character(df.compare[[1]])

            nonempty_first_values <- first_column_values[
              !is.na(first_column_values) &
                first_column_values != ""
            ]

            looks_like_gene_ids <- (
              length(nonempty_first_values) > 0 &&
              mean(
                grepl(
                  "^[[:alnum:]_.:-]+$",
                  nonempty_first_values
                )
              ) >= 0.90
            )

            if (looks_like_gene_ids) {
              comparison_id_candidates <- 1L
            }
          }

          if (length(comparison_id_candidates) != 1) {
            stop(
              "Could not identify exactly one gene-ID column in comparison table. ",
              "Detected columns: ",
              paste(
                paste0("[", colnames(df.compare), "]"),
                collapse = ", "
              )
            )
          }

          original_gene_id_column <- colnames(df.compare)[
            comparison_id_candidates
          ]

          colnames(df.compare)[comparison_id_candidates] <- "X"

          message(
            "[INFO] Comparison gene-ID column normalized to X from: [",
            original_gene_id_column,
            "]"
          )
        }

        df.compare$X <- normalize_host_gene_id(df.compare$X)

        required_de_columns <- c(
          "X",
          "log2FoldChange",
          "pvalue"
        )

        missing_de_columns <- setdiff(
          required_de_columns,
          colnames(df.compare)
        )

        if (length(missing_de_columns) > 0) {
          stop(
            "Comparison table is missing required columns: ",
            paste(missing_de_columns, collapse = ", "),
            ". Detected columns: ",
            paste(
              paste0("[", colnames(df.compare), "]"),
              collapse = ", "
            )
          )
        }

        if (
          all(is.na(df.compare$X)) ||
          all(df.compare$X == "")
        ) {
          stop(
            "The comparison gene-ID column was detected, ",
            "but it contains no usable gene IDs."
          )
        }

        message(
          "[INFO] Comparison genes loaded: ",
          length(unique(df.compare$X))
        )
        #### END BLOCK: comparison gene-ID column normalization ####

        df.compare.sig <- df.compare[
          !is.na(df.compare$pvalue) &
            df.compare$pvalue < 0.05 &
            !is.na(df.compare$log2FoldChange) &
            abs(df.compare$log2FoldChange) > 0.5,
          ,
          drop = FALSE
        ]

# Prepare geneList for GSEA using the keytype available in the OrgDb.
# log2FoldChange is retained as the ranking value because the same vector is
# also used by pathview for biologically interpretable pathway colouring.
deg_anno <- read.csv("../host_counts_DEG.csv", header = TRUE, check.names = FALSE)

df.compare <- df.compare[
  !is.na(df.compare$log2FoldChange) &
    is.finite(df.compare$log2FoldChange),
  ,
  drop = FALSE
]

if (HOST_GENE_KEYTYPE == "ENTREZID") {

  id_map <- deg_anno[, c("gene_id", "gene_name")]
  id_map <- id_map[
    !is.na(id_map$gene_name) &
      id_map$gene_name != "" &
      id_map$gene_name != "-",
  ]
  id_map <- id_map[!duplicated(id_map$gene_id), ]

  df.compare <- merge(df.compare, id_map, by.x = "X", by.y = "gene_id")

  genelist_symbol <- df.compare$log2FoldChange
  names(genelist_symbol) <- as.character(df.compare$gene_name)

  genelist_symbol <- genelist_symbol[!is.na(genelist_symbol)]
  genelist_symbol <- genelist_symbol[!is.na(names(genelist_symbol))]
  genelist_symbol <- genelist_symbol[names(genelist_symbol) != ""]
  genelist_symbol <- genelist_symbol[names(genelist_symbol) != "-"]

  ord <- order(abs(genelist_symbol), decreasing = TRUE)
  genelist_symbol <- genelist_symbol[ord]
  genelist_symbol <- genelist_symbol[!duplicated(names(genelist_symbol))]

  sym2ent <- AnnotationDbi::select(
    do.db,
    keys = unique(names(genelist_symbol)),
    keytype = "SYMBOL",
    columns = c("SYMBOL", "ENTREZID")
  )

  sym2ent <- sym2ent[!is.na(sym2ent$ENTREZID), ]
  sym2ent <- sym2ent[!duplicated(sym2ent$SYMBOL), ]

  genelist <- genelist_symbol[sym2ent$SYMBOL]
  names(genelist) <- sym2ent$ENTREZID

} else {

  # Para OrgDb custom, ex. org.Bglabrata.eg.db:
  # usar diretamente os IDs do featureCounts/host_counts, tipo BGLAX_...
  genelist <- df.compare$log2FoldChange
  names(genelist) <- as.character(df.compare$X)

  genelist <- genelist[!is.na(genelist)]
  genelist <- genelist[!is.na(names(genelist))]
  genelist <- genelist[names(genelist) != ""]
  genelist <- genelist[names(genelist) != "-"]

  valid_keys <- AnnotationDbi::keys(do.db, keytype = HOST_GENE_KEYTYPE)
  genelist <- genelist[names(genelist) %in% valid_keys]
}

ord <- order(abs(genelist), decreasing = TRUE)
genelist <- genelist[ord]
genelist <- genelist[!duplicated(names(genelist))]
genelist <- sort(genelist, decreasing = TRUE)

message("Genes sent to GO GSEA: ", length(genelist))
message("GO keyType used: ", HOST_GENE_KEYTYPE)
message("First geneList IDs for GSEA:")
print(head(names(genelist), 20))
        
        ## GSEA for GO ##
ego <- gseGO(geneList     = genelist,
             OrgDb        = do.db,
             keyType      = HOST_GENE_KEYTYPE,
             ont          = "ALL",
             minGSSize    = 10,
             maxGSSize    = 500,
             pvalueCutoff = 0.05,
             verbose      = FALSE)
        
        # save the full table of GSEA GO results
        write.csv(ego@result,"GO/GSEA_GO_results.csv")
        egox <- tryCatch({
          setReadable(ego, do.db, keyType = HOST_GENE_KEYTYPE)
        }, error = function(e) {
          message("setReadable failed for GO GSEA: ", conditionMessage(e))
          ego
        })
        write.csv(egox@result,"GO/GSEA_GO_results_symbol.csv")
        # draw plots for GO GSEA results
        try(plots4gsea("GO",ego,egox,"GO", genelist, group1, group2))
        
## Dual KEGG analysis ##
## 1) official_gseKEGG: uses ENTREZ / NCBI Gene IDs and KEGG organism code
## 2) eggNOG_KEGG_derived: uses GID / Ensembl IDs and custom TERM2GENE from OrgDb

ko.db <- get_host_kegg_code(host_sp, args[3])

message("KEGG organism code used for official KEGG: [", ko.db, "]")

# Keep a GID-based geneList for eggNOG-derived KEGG.
# For the custom OrgDb, genelist is already named with GID/ENSMLUG.
if (HOST_GENE_KEYTYPE == "GID") {
  genelist_gid <- genelist
} else {
  genelist_gid <- df.compare$log2FoldChange
  names(genelist_gid) <- as.character(df.compare$X)

  genelist_gid <- genelist_gid[!is.na(genelist_gid)]
  genelist_gid <- genelist_gid[!is.na(names(genelist_gid))]
  genelist_gid <- genelist_gid[names(genelist_gid) != ""]
  genelist_gid <- genelist_gid[names(genelist_gid) != "-"]

  if ("GID" %in% AnnotationDbi::keytypes(do.db)) {
    valid_gid <- AnnotationDbi::keys(do.db, keytype = "GID")
    genelist_gid <- genelist_gid[names(genelist_gid) %in% valid_gid]
  }

  genelist_gid <- dedup_ranked_gene_list(genelist_gid)
}

# Significant GID genes are reserved for ORA/enricher. GSEA above uses the
# complete ranked list.
sig_genes_gid <- unique(normalize_host_gene_id(df.compare.sig$X))
sig_genes_gid <- sig_genes_gid[
  !is.na(sig_genes_gid) &
    sig_genes_gid != "" &
    sig_genes_gid != "-"
]

if ("GID" %in% AnnotationDbi::keytypes(do.db)) {
  valid_gid_for_enricher <- AnnotationDbi::keys(do.db, keytype = "GID")
  valid_gid_for_enricher <- normalize_host_gene_id(valid_gid_for_enricher)
  sig_genes_gid <- intersect(sig_genes_gid, valid_gid_for_enricher)
}

write(
  paste0(
    "GID ranked genes for eggNOG-derived KEGG: ", length(genelist_gid), "\n",
    "Significant GID genes for enricher: ", length(sig_genes_gid), "\n",
    "HOST_GENE_KEYTYPE: ", HOST_GENE_KEYTYPE, "\n",
    "HOST_HAS_ENTREZ: ", HOST_HAS_ENTREZ, "\n",
    "HOST_HAS_KEGG_DERIVED: ", HOST_HAS_KEGG_DERIVED, "\n"
  ),
  "KEGG/KEGG_dual_input_summary.txt"
)

## ------------------------------------------------------------
## A) Official KEGG: gseKEGG / gseMKEGG
## ------------------------------------------------------------

if (is.na(ko.db) || ko.db == "") {

  write(
    "No KEGG organism code found for this host species. Official KEGG was skipped.",
    "KEGG/official_gseKEGG/No_KEGG_organism_code_found.txt"
  )

} else {

  official_gene_list <- prepare_official_kegg_gene_list(
    genelist_gid_or_entrez = genelist,
    host_sp = host_sp,
    taxid = args[3],
    do.db = do.db,
    host_gene_keytype = HOST_GENE_KEYTYPE,
    outdir = "KEGG/official_gseKEGG"
  )

  official_gene_list <- dedup_ranked_gene_list(official_gene_list)

  message("Genes sent to official KEGG: ", length(official_gene_list))
  message("First official KEGG gene IDs:")
  print(head(names(official_gene_list), 20))

  if (length(official_gene_list) < 10) {

    write(
      paste0(
        "Official KEGG was skipped because fewer than 10 genes could be mapped to ENTREZ.\n",
        "Mapped genes: ", length(official_gene_list), "\n",
        "Organism code: ", ko.db, "\n"
      ),
      "KEGG/official_gseKEGG/No_official_KEGG_geneList.txt"
    )

  } else {

    kk.p <- gseKEGG.safe(
      geneList = official_gene_list,
      ko.db = ko.db,
      minGSSize = 8,
      maxGSSize = 500,
      pvalueCutoff = 0.05
    )

    kk.f <- tryCatch({
      gseMKEGG(
        geneList = official_gene_list,
        organism = ko.db,
        keyType = "ncbi-geneid",
        minGSSize = 8,
        pvalueCutoff = 0.05
      )
    }, error = function(e) {
      message("Official KEGG module GSEA failed: ", conditionMessage(e))
      NULL
    })

    if (is.null(kk.p)) {

      write(
        paste0("Official KEGG pathway GSEA failed for organism code: ", ko.db),
        "KEGG/official_gseKEGG/No_KEGG_enrichment_result.txt"
      )

    } else {

      write.csv(
        kk.p@result,
        "KEGG/official_gseKEGG/GSEA_KEGG_official_results.csv",
        row.names = FALSE
      )

kkx.p <- kk.p

write(
  paste0(
    "setReadable was skipped for official KEGG because official KEGG uses ENTREZ/KEGG IDs, ",
    "while the custom OrgDb may use GID/Ensembl IDs.\n"
  ),
  "KEGG/official_gseKEGG/setReadable_skipped_official_KEGG.txt"
)

write.csv(
  kkx.p@result,
  "KEGG/official_gseKEGG/GSEA_KEGG_official_results_symbol.csv",
  row.names = FALSE
)

try(
  plots4gsea(
    "KEGG/official_gseKEGG",
    kk.p,
    kkx.p,
    "official_KEGG",
    official_gene_list,
    group1,
    group2
  ),
  silent = TRUE
)

      setwd("KEGG/official_gseKEGG/Pathview")
      try(pathview.p(kk.p, ko.db, official_gene_list), silent = TRUE)
      setwd(paste0(dirname(args[1]), "/Host_DEG/", group1, "_vs_", group2))
    }

    if (is.null(kk.f)) {

      write(
        paste0("Official KEGG module GSEA failed for organism code: ", ko.db),
        "KEGG/official_gseKEGG/Modules/No_KEGG_module_enrichment_result.txt"
      )

    } else {

      write.csv(
        kk.f@result,
        "KEGG/official_gseKEGG/Modules/GSEA_KEGG_modules_official_results.csv",
        row.names = FALSE
      )

kkx.f <- kk.f

write(
  paste0(
    "setReadable was skipped for official KEGG modules because official KEGG uses ENTREZ/KEGG IDs, ",
    "while the custom OrgDb may use GID/Ensembl IDs.\n"
  ),
  "KEGG/official_gseKEGG/Modules/setReadable_skipped_official_KEGG_modules.txt"
)
      write.csv(
        kkx.f@result,
        "KEGG/official_gseKEGG/Modules/GSEA_KEGG_modules_official_results_symbol.csv",
        row.names = FALSE
      )

      try(
        plots4gsea(
          "KEGG/official_gseKEGG/Modules",
          kk.f,
          kkx.f,
          "official_KEGG_modules",
          official_gene_list,
          group1,
          group2
        ),
        silent = TRUE
      )

      setwd("KEGG/official_gseKEGG/Modules/Pathview")
      try(
        pathview.modules.p(
          kk_module = kk.f,
          ko.db = ko.db,
          kegg_gene_list = official_gene_list
        ),
        silent = TRUE
      )
      setwd(paste0(dirname(args[1]), "/Host_DEG/", group1, "_vs_", group2))
    }
  }
}

## ------------------------------------------------------------
## B) eggNOG-derived KEGG: custom GSEA/enricher with GID IDs
## ------------------------------------------------------------

if (!HOST_HAS_KEGG_DERIVED) {

  write(
    paste0(
      "eggNOG-derived KEGG was skipped because this OrgDb does not provide KEGG-derived columns.\n",
      "Available columns: ",
      paste(AnnotationDbi::columns(do.db), collapse = ", "),
      "\n"
    ),
    "KEGG/eggNOG_KEGG_derived/No_KEGG_derived_columns_in_OrgDb.txt"
  )

} else {

  write(
    paste0(
      "Running eggNOG-derived KEGG using GID/Ensembl IDs.\n",
      "Ranked GID genes: ", length(genelist_gid), "\n",
      "Significant GID genes: ", length(sig_genes_gid), "\n"
    ),
    "KEGG/eggNOG_KEGG_derived/KEGG_derived_method.txt"
  )

  if ("KEGG_PATHWAY" %in% AnnotationDbi::columns(do.db)) {
    run_kegg_derived_gsea_and_enricher(
      genelist_gid = genelist_gid,
      sig_genes_gid = sig_genes_gid,
      orgdb = do.db,
      column = "KEGG_PATHWAY",
      prefix_regex = "^map[0-9]{5}$",
      label = "KEGG_PATHWAY_map",
      outdir = "KEGG/eggNOG_KEGG_derived/Pathway",
      group1 = group1,
      group2 = group2
    )
  }

  if ("KEGG_MODULE" %in% AnnotationDbi::columns(do.db)) {
    run_kegg_derived_gsea_and_enricher(
      genelist_gid = genelist_gid,
      sig_genes_gid = sig_genes_gid,
      orgdb = do.db,
      column = "KEGG_MODULE",
      prefix_regex = "^M[0-9]{5}$",
      label = "KEGG_MODULE",
      outdir = "KEGG/eggNOG_KEGG_derived/Module",
      group1 = group1,
      group2 = group2
    )
  }

  if ("KEGG_KO" %in% AnnotationDbi::columns(do.db)) {
    run_kegg_derived_gsea_and_enricher(
      genelist_gid = genelist_gid,
      sig_genes_gid = sig_genes_gid,
      orgdb = do.db,
      column = "KEGG_KO",
      prefix_regex = "^ko:K[0-9]+$",
      label = "KEGG_KO",
      outdir = "KEGG/eggNOG_KEGG_derived/KO",
      group1 = group1,
      group2 = group2
    )
  }

  if ("KEGG_REACTION" %in% AnnotationDbi::columns(do.db)) {
    run_kegg_derived_gsea_and_enricher(
      genelist_gid = genelist_gid,
      sig_genes_gid = sig_genes_gid,
      orgdb = do.db,
      column = "KEGG_REACTION",
      prefix_regex = "^R[0-9]+$",
      label = "KEGG_REACTION",
      outdir = "KEGG/eggNOG_KEGG_derived/Reaction",
      group1 = group1,
      group2 = group2
    )
  }
}
        
        setwd("../") # go back to the Host_DEG folder
      }
    }
    
    ## run GSEA enrichment analysis ##
    enrichment(coldata_vs,args[1],do.db.obj)
    
    # function for preparing genelist for biological theme comparison - compareCluster
BTC <- function(coldata_vs, do.db) {
  genelist.ct <- list()

  valid_keys <- AnnotationDbi::keys(do.db, keytype = HOST_GENE_KEYTYPE)

  for (i in 1:nrow(coldata_vs)) {
    group1 <- coldata_vs$group1[i]
    group2 <- coldata_vs$group2[i]

    df.btc <- read.csv(
      paste0(getwd(), "/", group1, "_vs_", group2, "/host_counts_", group1, "_vs_", group2, ".csv"),
      header = TRUE
    )

    flt_up <- df.btc[df.btc$log2FoldChange > 0.5 & df.btc$pvalue < 0.05, ]
    flt_down <- df.btc[df.btc$log2FoldChange < -0.5 & df.btc$pvalue < 0.05, ]

    get_ids_for_go <- function(x) {
      if (nrow(x) == 0) {
        return(character(0))
      }

      ids <- unique(as.character(x$X))
      ids <- ids[!is.na(ids) & ids != "" & ids != "-"]
      ids <- ids[ids %in% valid_keys]

      return(ids)
    }

    dedup_UP <- get_ids_for_go(flt_up)
    dedup_DOWN <- get_ids_for_go(flt_down)

    genelist.c <- list(dedup_UP, dedup_DOWN)
    names(genelist.c) <- c(
      paste0(group1, "_vs_", group2, "_UP"),
      paste0(group1, "_vs_", group2, "_DOWN")
    )

    genelist.ct <- c(genelist.ct, genelist.c)
  }

  return(genelist.ct)
}

    ## run biological theme comparison ##
    genelist.ct <- BTC(coldata_vs,do.db.obj)
# GO enrichment comparison
# User-adjustable settings
go_pvalue_cutoff <- 0.05
go_qvalue_cutoff <- 0.2
go_padjust_method <- "BH"
go_ont <- "ALL"

cgo <- tryCatch({
  compareCluster(
    genelist.ct,
    fun = enrichGO,
    OrgDb = do.db.obj,
    keyType = HOST_GENE_KEYTYPE,
    ont = go_ont,
    pAdjustMethod = go_padjust_method,
    pvalueCutoff = go_pvalue_cutoff,
    qvalueCutoff = go_qvalue_cutoff
  )
}, error = function(e) {
  warning(
    "GO compareCluster failed or found no enrichment. ",
    "Settings used: ont=", go_ont,
    ", pAdjustMethod=", go_padjust_method,
    ", pvalueCutoff=", go_pvalue_cutoff,
    ", qvalueCutoff=", go_qvalue_cutoff,
    ". Error: ", conditionMessage(e)
  )
  NULL
})

if (!is.null(cgo) && nrow(cgo@compareClusterResult) > 0) {

  cgo <- tryCatch({
    setReadable(
      cgo,
      OrgDb = do.db.obj,
      keyType = HOST_GENE_KEYTYPE
    )
  }, error = function(e) {
    warning("setReadable failed for GO compareCluster: ", conditionMessage(e))
    cgo
  })

  write.csv(
    cgo@compareClusterResult,
    "biological_theme_comparison_GO_results.csv",
    row.names = FALSE
  )

    ## ------------------------------------------------------------
  ## Combined GO dotplot: BP + CC + MF
  ## ------------------------------------------------------------

  p.cgo.dot <- dotplot(
    cgo,
    showCategory = nrow(cgo@compareClusterResult)
  ) +
    theme(
      axis.text.x = element_text(
        angle = 45,
        vjust = 1,
        hjust = 1
      )
    )

  ggsave(
    filename = "biological_theme_comparison_GO.pdf",
    plot = p.cgo.dot,
    limitsize = FALSE,
    height = 0.54 * nrow(cgo@compareClusterResult),
    width = 3 * length(
      unique(cgo@compareClusterResult$Cluster)
    )
  )
  ## ------------------------------------------------------------
  ## Additional publication-ready faceted GO dotplot
  ## ------------------------------------------------------------
  ## The original combined dotplot above is preserved.
  ##
  ## This additional script generates a faceted dotplot with:
  ##   rows    = BP, CC and MF;
  ##   columns = UP and DOWN;
  ##   x-axis  = GeneRatio;
  ##   size    = gene count;
  ##   colour  = -log10(adjusted p-value).
  ## ------------------------------------------------------------

  go_faceted_top_n <- 5L

  ## Locate the MTD installation directory from the path of the
  ## currently running DEG_Anno_Plot.R script.

  all_command_args <- commandArgs(
    trailingOnly = FALSE
  )

  current_script_argument <- grep(
    "^--file=",
    all_command_args,
    value = TRUE
  )

  if (length(current_script_argument) == 0L) {

    warning(
      "Could not determine the location of DEG_Anno_Plot.R. ",
      "The additional faceted GO dotplot will be skipped."
    )

    writeLines(
      "Could not determine the location of DEG_Anno_Plot.R.",
      "biological_theme_comparison_GO_faceted_dotplot_skipped.txt"
    )

  } else {

    current_script_path <- sub(
      "^--file=",
      "",
      current_script_argument[1]
    )

    current_script_path <- normalizePath(
      current_script_path,
      mustWork = TRUE
    )

    # HostSpecies.csv is stored in the MTD installation root.
    # Using args[4] avoids incorrectly treating the directory containing
    # DEG_Anno_Plot.R as the MTD root.
    mtd_installation_directory <- dirname(
      normalizePath(
        args[4],
        mustWork = TRUE
      )
    )

    go_faceted_script <- file.path(
      mtd_installation_directory,
      "aux_scripts",
      "GO",
      "GO_faceted_dotplot.R"
    )

    message(
      "[INFO] MTD installation directory: ",
      mtd_installation_directory
    )

    message(
      "[INFO] Faceted GO dotplot script path: ",
      go_faceted_script
    )

    go_host_deg_directory <- normalizePath(
      getwd(),
      mustWork = TRUE
    )

    go_faceted_log <- paste0(
      "biological_theme_comparison_GO_faceted_dotplot_top",
      go_faceted_top_n,
      ".log"
    )

    go_faceted_expected_pdf <- paste0(
      "biological_theme_comparison_GO_faceted_dotplot_top",
      go_faceted_top_n,
      ".pdf"
    )

    go_faceted_expected_tiff <- paste0(
      "biological_theme_comparison_GO_faceted_dotplot_top",
      go_faceted_top_n,
      ".tiff"
    )

    go_faceted_expected_data <- paste0(
      "biological_theme_comparison_GO_faceted_dotplot_top",
      go_faceted_top_n,
      "_data.csv"
    )

    if (!file.exists(go_faceted_script)) {

      warning(
        "Faceted GO dotplot script was not found: ",
        go_faceted_script
      )

      writeLines(
        c(
          "Faceted GO dotplot script was not found.",
          paste0(
            "Expected script: ",
            go_faceted_script
          )
        ),
        "biological_theme_comparison_GO_faceted_dotplot_skipped.txt"
      )

    } else {

      rscript_binary <- unname(
        Sys.which("Rscript")
      )

      if (!nzchar(rscript_binary)) {

        warning(
          "Rscript could not be found in the current environment."
        )

        writeLines(
          "Rscript could not be found in the current environment.",
          "biological_theme_comparison_GO_faceted_dotplot_skipped.txt"
        )

      } else {

        message(
          "============================================================"
        )

        message(
          "[GO DOTPLOT] Running additional faceted GO dotplot"
        )

        message(
          "[GO DOTPLOT] Script: ",
          go_faceted_script
        )

        message(
          "[GO DOTPLOT] Input: ",
          go_host_deg_directory
        )

        message(
          "[GO DOTPLOT] Top terms per ontology/direction: ",
          go_faceted_top_n
        )

        message(
          "============================================================"
        )

        go_faceted_execution <- tryCatch(
          {

            command_output <- suppressWarnings(
              system2(
                command = rscript_binary,
                args = c(
                  shQuote(go_faceted_script),
                  shQuote(go_host_deg_directory)
                ),
                env = paste0(
                  "GO_TOP_N=",
                  go_faceted_top_n
                ),
                stdout = TRUE,
                stderr = TRUE
              )
            )

            command_status <- attr(
              command_output,
              "status"
            )

            if (is.null(command_status)) {
              command_status <- 0L
            }

            list(
              output = command_output,
              status = as.integer(command_status),
              error = NULL
            )
          },
          error = function(e) {

            list(
              output = character(0),
              status = 1L,
              error = conditionMessage(e)
            )
          }
        )

        go_faceted_log_content <- c(
          paste0(
            "Script: ",
            go_faceted_script
          ),
          paste0(
            "Input directory: ",
            go_host_deg_directory
          ),
          paste0(
            "GO_TOP_N=",
            go_faceted_top_n
          ),
          paste0(
            "Exit status: ",
            go_faceted_execution$status
          ),
          "",
          go_faceted_execution$output
        )

        if (!is.null(go_faceted_execution$error)) {

          go_faceted_log_content <- c(
            go_faceted_log_content,
            "",
            paste0(
              "Execution error: ",
              go_faceted_execution$error
            )
          )
        }

        writeLines(
          go_faceted_log_content,
          go_faceted_log
        )

        if (length(go_faceted_execution$output) > 0) {

          message(
            paste(
              tail(
                go_faceted_execution$output,
                20
              ),
              collapse = "\n"
            )
          )
        }

        if (
          go_faceted_execution$status == 0L &&
          file.exists(go_faceted_expected_pdf)
        ) {

          unlink(
            "biological_theme_comparison_GO_faceted_dotplot_skipped.txt"
          )

          message(
            "[GO DOTPLOT] Faceted PDF generated: ",
            go_faceted_expected_pdf
          )

          if (file.exists(go_faceted_expected_tiff)) {

            message(
              "[GO DOTPLOT] Faceted TIFF generated: ",
              go_faceted_expected_tiff
            )
          }

          if (file.exists(go_faceted_expected_data)) {

            message(
              "[GO DOTPLOT] Selected-term table generated: ",
              go_faceted_expected_data
            )
          }

        } else {

          warning(
            "The additional faceted GO dotplot was not generated. ",
            "See log: ",
            go_faceted_log
          )

          writeLines(
            c(
              "The additional faceted GO dotplot was skipped.",
              paste0(
                "Exit status: ",
                go_faceted_execution$status
              ),
              paste0(
                "Expected PDF: ",
                go_faceted_expected_pdf
              ),
              paste0(
                "Log: ",
                go_faceted_log
              )
            ),
            "biological_theme_comparison_GO_faceted_dotplot_skipped.txt"
          )
        }
      }
    }
  }
  ## ------------------------------------------------------------
  ## GO cnetplots separated by ontology
  ## ------------------------------------------------------------
  ## enrichplot 1.14.2 can generate NA values in Cluster when
  ## BP, CC and MF are sent together to cnetplot().
  ##
  ## Therefore:
  ##   - the dotplot remains combined;
  ##   - cnetplots are generated separately for BP, CC and MF;
  ##   - only category labels are printed;
  ##   - gene nodes and edges remain in the network;
  ##   - individual and multi-page PDF files are generated.
  ## ------------------------------------------------------------

  go_cnet_show <- 10
  go_cnet_width <- 14
  go_cnet_height <- 11

  go_ontology_titles <- c(
    BP = "GO Biological Process (BP)",
    CC = "GO Cellular Component (CC)",
    MF = "GO Molecular Function (MF)"
  )

  build_go_cnetplot <- function(
    cgo_object,
    ontology_name,
    show_category = 10
  ) {

    go_df <- cgo_object@compareClusterResult

    if (!"ONTOLOGY" %in% colnames(go_df)) {
      stop(
        "The GO compareCluster result does not contain ",
        "an ONTOLOGY column."
      )
    }

    ontology_values <- toupper(
      trimws(
        as.character(go_df$ONTOLOGY)
      )
    )

    ontology_df <- go_df[
      !is.na(ontology_values) &
      ontology_values == ontology_name,
      ,
      drop = FALSE
    ]

    if (nrow(ontology_df) == 0) {
      stop(
        "No GO results were found for ontology ",
        ontology_name,
        "."
      )
    }

    required_columns <- c(
      "Cluster",
      "ID",
      "Description",
      "GeneRatio",
      "geneID",
      "Count"
    )

    missing_columns <- setdiff(
      required_columns,
      colnames(ontology_df)
    )

    if (length(missing_columns) > 0) {
      stop(
        "Missing columns in GO compareCluster result: ",
        paste(missing_columns, collapse = ", ")
      )
    }

    ## Remove rows that cannot be represented in the network.

    valid_rows <- (
      !is.na(ontology_df$Cluster) &
      nzchar(trimws(as.character(ontology_df$Cluster))) &
      !is.na(ontology_df$ID) &
      nzchar(trimws(as.character(ontology_df$ID))) &
      !is.na(ontology_df$Description) &
      nzchar(trimws(as.character(ontology_df$Description))) &
      !is.na(ontology_df$geneID) &
      nzchar(trimws(as.character(ontology_df$geneID))) &
      !is.na(ontology_df$Count) &
      ontology_df$Count > 0
    )

    ontology_df <- ontology_df[
      valid_rows,
      ,
      drop = FALSE
    ]

    if (nrow(ontology_df) == 0) {
      stop(
        "No valid GO rows remained for ontology ",
        ontology_name,
        "."
      )
    }

    ## Remove unused Cluster levels. This prevents fortify() from
    ## creating invalid cluster values.

    cluster_levels <- unique(
      as.character(ontology_df$Cluster)
    )

    ontology_df$Cluster <- factor(
      as.character(ontology_df$Cluster),
      levels = cluster_levels
    )

    ontology_df$ID <- as.character(
      ontology_df$ID
    )

    ontology_df$Description <- as.character(
      ontology_df$Description
    )

    ontology_df$geneID <- as.character(
      ontology_df$geneID
    )

    rownames(ontology_df) <- NULL

    ## Rebuild the geneClusters slot from the ontology subset.

    gene_clusters <- lapply(
      split(
        ontology_df$geneID,
        ontology_df$Cluster,
        drop = TRUE
      ),
      function(gene_strings) {

        genes <- unlist(
          strsplit(
            as.character(gene_strings),
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
    )

    ## Reconstruct a clean compareClusterResult containing only
    ## one ontology.

    cgo_ontology <- methods::new(
      "compareClusterResult",
      compareClusterResult = ontology_df,
      geneClusters = gene_clusters,
      fun = "enrichGO"
    )

    ## Validate the exact data that cnetplot() will consume.

    fortified <- ggplot2::fortify(
      cgo_ontology,
      showCategory = show_category,
      includeAll = TRUE
    )

    if (nrow(fortified) == 0) {
      stop(
        "fortify() returned no terms for ontology ",
        ontology_name,
        "."
      )
    }

    if (anyNA(fortified$Cluster)) {
      stop(
        "fortify() generated ",
        sum(is.na(fortified$Cluster)),
        " NA Cluster value(s) for ontology ",
        ontology_name,
        "."
      )
    }

    ontology_title <- if (
      ontology_name %in% names(go_ontology_titles)
    ) {
      unname(
        go_ontology_titles[ontology_name]
      )
    } else {
      paste("GO", ontology_name)
    }

    ## node_label = "category" hides gene-ID labels while keeping
    ## gene nodes and gene-category connections in the network.

    p <- cnetplot(
      cgo_ontology,
      showCategory = show_category,
      layout = "kk",
      node_label = "category",
      cex_category = 1.10,
      cex_gene = 0.45,
      cex_label_category = 0.85
    ) +
      ggplot2::labs(
        title = ontology_title
      ) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(
          hjust = 0.5,
          face = "bold",
          size = 14
        ),
        legend.position = "right"
      )

    return(p)
  }

  ## ------------------------------------------------------------
  ## Generate BP, CC and MF
  ## ------------------------------------------------------------

  available_ontologies <- unique(
    toupper(
      trimws(
        as.character(
          cgo@compareClusterResult$ONTOLOGY
        )
      )
    )
  )

  available_ontologies <- available_ontologies[
    !is.na(available_ontologies) &
    nzchar(available_ontologies)
  ]

  ontology_order <- c(
    "BP",
    "CC",
    "MF"
  )

  ontologies_to_plot <- intersect(
    ontology_order,
    available_ontologies
  )

  go_cnet_plots <- list()
  go_cnet_errors <- character()

  for (ontology_name in ontologies_to_plot) {

    message(
      "Creating GO compareCluster cnetplot for ontology: ",
      ontology_name
    )

    plot_attempt <- tryCatch(
      {

        plot_object <- build_go_cnetplot(
          cgo_object = cgo,
          ontology_name = ontology_name,
          show_category = go_cnet_show
        )

        list(
          plot = plot_object,
          error = NULL
        )
      },
      error = function(e) {

        list(
          plot = NULL,
          error = conditionMessage(e)
        )
      }
    )

    if (is.null(plot_attempt$plot)) {

      go_cnet_errors[ontology_name] <- plot_attempt$error

      message(
        "GO ",
        ontology_name,
        " cnetplot skipped: ",
        plot_attempt$error
      )

      next
    }

    ontology_output <- paste0(
      "biological_theme_comparison_GO_net_",
      ontology_name,
      ".pdf"
    )

    save_attempt <- tryCatch(
      {

        ggsave(
          filename = ontology_output,
          plot = plot_attempt$plot,
          limitsize = FALSE,
          width = go_cnet_width,
          height = go_cnet_height
        )

        list(
          success = TRUE,
          error = NULL
        )
      },
      error = function(e) {

        list(
          success = FALSE,
          error = conditionMessage(e)
        )
      }
    )

    if (save_attempt$success) {

      go_cnet_plots[[ontology_name]] <- plot_attempt$plot

      message(
        "Saved GO cnetplot: ",
        ontology_output
      )

    } else {

      go_cnet_errors[ontology_name] <- save_attempt$error

      message(
        "Failed to save GO ",
        ontology_name,
        " cnetplot: ",
        save_attempt$error
      )
    }
  }

  ## ------------------------------------------------------------
  ## Generate the original expected filename as a multi-page PDF
  ## ------------------------------------------------------------

  if (length(go_cnet_plots) > 0) {

    combined_pdf_open <- FALSE

    combined_attempt <- tryCatch(
      {

        grDevices::pdf(
          file = "biological_theme_comparison_GO_net.pdf",
          width = go_cnet_width,
          height = go_cnet_height,
          onefile = TRUE
        )

        combined_pdf_open <- TRUE

        for (ontology_name in names(go_cnet_plots)) {
          print(
            go_cnet_plots[[ontology_name]]
          )
        }

        grDevices::dev.off()
        combined_pdf_open <- FALSE

        list(
          success = TRUE,
          error = NULL
        )
      },
      error = function(e) {

        if (
          combined_pdf_open &&
          grDevices::dev.cur() > 1
        ) {
          grDevices::dev.off()
        }

        list(
          success = FALSE,
          error = conditionMessage(e)
        )
      }
    )

    if (combined_attempt$success) {

      unlink(
        "biological_theme_comparison_GO_net_skipped.txt"
      )

      message(
        "Saved multi-page GO cnetplot: ",
        "biological_theme_comparison_GO_net.pdf"
      )

    } else {

      message(
        "Failed to save multi-page GO cnetplot: ",
        combined_attempt$error
      )
    }

    if (length(go_cnet_errors) > 0) {

      writeLines(
        c(
          "Some GO ontology cnetplots were skipped:",
          paste0(
            names(go_cnet_errors),
            ": ",
            unname(go_cnet_errors)
          )
        ),
        "biological_theme_comparison_GO_net_partial_warnings.txt"
      )

    } else {

      unlink(
        "biological_theme_comparison_GO_net_partial_warnings.txt"
      )
    }

  } else {

    error_text <- if (length(go_cnet_errors) > 0) {

      paste0(
        names(go_cnet_errors),
        ": ",
        unname(go_cnet_errors),
        collapse = "\n"
      )

    } else {

      paste0(
        "No BP, CC or MF ontology was available in the ",
        "GO compareCluster result."
      )
    }

    message(
      "GO compareCluster cnetplots skipped: ",
      error_text
    )

    writeLines(
      c(
        "GO compareCluster cnetplots skipped.",
        error_text
      ),
      "biological_theme_comparison_GO_net_skipped.txt"
    )
  }

} else {

  warning(
    "No significant GO enrichment found by compareCluster. ",
    "GO compareCluster plots will be skipped."
  )

  write(
    paste0(
      "No significant GO enrichment found by compareCluster.\n",
      "This is not a pipeline failure.\n\n",
      "Settings used:\n",
      "  ont = ", go_ont, "\n",
      "  pAdjustMethod = ", go_padjust_method, "\n",
      "  pvalueCutoff = ", go_pvalue_cutoff, "\n",
      "  qvalueCutoff = ", go_qvalue_cutoff, "\n\n",
      "Input gene clusters were generated, but no GO term passed these cutoffs.\n",
      "GO compareCluster plots were skipped.\n"
    ),
    "No_GO_compareCluster_result.txt"
  )
}

#### BEGIN BLOCK: dual KEGG enrichment comparison ####

## A) eggNOG-derived KEGG compareCluster using GID clusters

if (HOST_HAS_KEGG_DERIVED && "KEGG_PATHWAY" %in% AnnotationDbi::columns(do.db.obj)) {

  term2gene_pathway <- make_kegg_derived_term2gene(
    orgdb = do.db.obj,
    column = "KEGG_PATHWAY",
    prefix_regex = "^map[0-9]{5}$",
    out_file = "TERM2GENE_KEGG_PATHWAY_map_compareCluster.tsv"
  )

  if (nrow(term2gene_pathway) > 0) {

    ck_derived <- tryCatch({
      compareCluster(
        genelist.ct,
        fun = "enricher",
        TERM2GENE = term2gene_pathway,
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2,
        minGSSize = 10,
        maxGSSize = 500
      )
    }, error = function(e) {
      message("eggNOG-derived KEGG compareCluster failed: ", conditionMessage(e))
      NULL
    })

    if (!is.null(ck_derived) && nrow(ck_derived@compareClusterResult) > 0) {

      write.csv(
        ck_derived@compareClusterResult,
        "biological_theme_comparison_KEGG_derived_PATHWAY_results.csv",
        row.names = FALSE
      )

      dotplot(
        ck_derived,
        showCategory = nrow(ck_derived@compareClusterResult)
      ) +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

      ggsave(
        "biological_theme_comparison_KEGG_derived_PATHWAY.pdf",
        limitsize = FALSE,
        height = 0.54 * nrow(ck_derived@compareClusterResult),
        width = 3 * length(unique(ck_derived@compareClusterResult$Cluster))
      )

      cnet_show <- max(1, round(nrow(ck_derived@compareClusterResult) / 6))

      try({
        ggsave(
          "biological_theme_comparison_KEGG_derived_PATHWAY_net.pdf",
          plot = cnetplot(
            ck_derived,
            cex_label_gene = 0.6,
            showCategory = cnet_show
          ),
          limitsize = FALSE
        )
      }, silent = TRUE)

    } else {

      write(
        paste0(
          "No significant eggNOG-derived KEGG pathway enrichment found by compareCluster.\n",
          "TERM2GENE rows: ", nrow(term2gene_pathway), "\n"
        ),
        "No_KEGG_derived_compareCluster_result.txt"
      )
    }
  }
}

## B) Official KEGG compareCluster
## Convert custom GID/Ensembl clusters to NCBI GeneID, obtain the
## organism-specific official KEGG pathway membership with KEGGREST,
## and run compareCluster through the generic enricher interface.
##
## This avoids the enrichKEGG ncbi-geneid conversion limitation observed
## for some non-model organisms, while retaining official KEGG pathway data.

#### BEGIN BLOCK: official KEGG TERM2GENE compareCluster ####

ko.db <- get_host_kegg_code(host_sp, args[3])

official_kegg_comparecluster_outputs <- c(
  "biological_theme_comparison_KEGG_official_results.csv",
  "biological_theme_comparison_KEGG_official.pdf",
  "No_KEGG_official_compareCluster_result.txt",
  paste0("TERM2GENE_KEGG_official_", ko.db, "_compareCluster.tsv"),
  paste0("TERM2NAME_KEGG_official_", ko.db, "_compareCluster.tsv"),
  paste0("KEGG_official_", ko.db, "_compareCluster_input_summary.tsv")
)

official_kegg_comparecluster_outputs <- official_kegg_comparecluster_outputs[
  !is.na(official_kegg_comparecluster_outputs) &
    !grepl("_NA_", official_kegg_comparecluster_outputs, fixed = TRUE)
]

unlink(
  official_kegg_comparecluster_outputs[
    file.exists(official_kegg_comparecluster_outputs)
  ]
)

if (!is.na(ko.db) && nzchar(ko.db)) {

  genelist.ct.entrez <- convert_gene_clusters_to_entrez(
    gene_clusters = genelist.ct,
    host_sp = host_sp,
    taxid = args[3],
    host_gene_keytype = HOST_GENE_KEYTYPE,
    outdir = getwd()
  )

  message(
    "[INFO] Official KEGG compareCluster groups with mapped ENTREZ IDs: ",
    length(genelist.ct.entrez)
  )

  if (length(genelist.ct.entrez) > 0) {

    suppressPackageStartupMessages(library("KEGGREST"))

    official_kegg_links_attempt <- tryCatch({
      list(
        value = KEGGREST::keggLink("pathway", ko.db),
        error = NULL
      )
    }, error = function(e) {
      list(
        value = NULL,
        error = conditionMessage(e)
      )
    })

    official_kegg_links <- official_kegg_links_attempt$value

    if (
      is.null(official_kegg_links) ||
      length(official_kegg_links) == 0
    ) {

      official_link_error <- official_kegg_links_attempt$error

      if (is.null(official_link_error)) {
        official_link_error <- "KEGG returned no pathway-gene links."
      }

      message(
        "[WARNING] Official KEGG compareCluster could not obtain pathway links: ",
        official_link_error
      )

      writeLines(
        c(
          "Official KEGG compareCluster was skipped because official KEGG",
          "pathway-to-gene links could not be obtained.",
          "",
          paste0("Organism code: ", ko.db),
          paste0("Mapped ENTREZ clusters: ", length(genelist.ct.entrez)),
          paste0("KEGGREST message: ", official_link_error)
        ),
        "No_KEGG_official_compareCluster_result.txt"
      )

    } else {

      official_kegg_links_df <- data.frame(
        from = names(official_kegg_links),
        to = unname(official_kegg_links),
        stringsAsFactors = FALSE
      )

      if (
        any(
          grepl(
            paste0("^", ko.db, ":"),
            official_kegg_links_df$from
          )
        )
      ) {
        official_gene_column <- official_kegg_links_df$from
        official_pathway_column <- official_kegg_links_df$to
      } else {
        official_gene_column <- official_kegg_links_df$to
        official_pathway_column <- official_kegg_links_df$from
      }

      official_term2gene <- data.frame(
        term = sub(
          "^path:",
          "",
          official_pathway_column
        ),
        gene = sub(
          paste0("^", ko.db, ":"),
          "",
          official_gene_column
        ),
        stringsAsFactors = FALSE
      )

      official_term2gene <- official_term2gene[
        grepl(
          paste0("^", ko.db, "[0-9]{5}$"),
          official_term2gene$term
        ) &
          !is.na(official_term2gene$gene) &
          official_term2gene$gene != "",
        ,
        drop = FALSE
      ]

      official_term2gene <- unique(official_term2gene)

      official_path_names_attempt <- tryCatch({
        list(
          value = KEGGREST::keggList("pathway", ko.db),
          error = NULL
        )
      }, error = function(e) {
        list(
          value = NULL,
          error = conditionMessage(e)
        )
      })

      official_path_names <- official_path_names_attempt$value

      if (
        is.null(official_path_names) ||
        length(official_path_names) == 0
      ) {
        official_term2name <- data.frame(
          term = unique(official_term2gene$term),
          name = unique(official_term2gene$term),
          stringsAsFactors = FALSE
        )

        message(
          "[WARNING] KEGG pathway names were unavailable; ",
          "pathway IDs will be used as names."
        )
      } else {
        official_term2name <- data.frame(
          term = sub(
            "^path:",
            "",
            names(official_path_names)
          ),
          name = as.character(official_path_names),
          stringsAsFactors = FALSE
        )

        official_term2name$name <- sub(
          " - .*$",
          "",
          official_term2name$name
        )

        official_term2name <- unique(official_term2name)

        official_term2name <- official_term2name[
          official_term2name$term %in% official_term2gene$term,
          ,
          drop = FALSE
        ]
      }

      official_term2gene_file <- paste0(
        "TERM2GENE_KEGG_official_",
        ko.db,
        "_compareCluster.tsv"
      )

      official_term2name_file <- paste0(
        "TERM2NAME_KEGG_official_",
        ko.db,
        "_compareCluster.tsv"
      )

      write.table(
        official_term2gene,
        official_term2gene_file,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
      )

      write.table(
        official_term2name,
        official_term2name_file,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
      )

      official_kegg_gene_ids <- unique(
        as.character(official_term2gene$gene)
      )

      cluster_names <- names(genelist.ct.entrez)

      if (
        is.null(cluster_names) ||
        any(is.na(cluster_names)) ||
        any(cluster_names == "")
      ) {
        cluster_names <- paste0(
          "Cluster_",
          seq_along(genelist.ct.entrez)
        )
        names(genelist.ct.entrez) <- cluster_names
      }

      official_cluster_input_counts <- lengths(genelist.ct.entrez)

      genelist.ct.kegg <- lapply(
        genelist.ct.entrez,
        function(ids) {
          ids <- unique(trimws(as.character(ids)))
          ids <- ids[
            !is.na(ids) &
              ids != "" &
              ids != "-"
          ]

          intersect(
            ids,
            official_kegg_gene_ids
          )
        }
      )

      official_cluster_overlap_counts <- lengths(genelist.ct.kegg)

      official_input_summary <- data.frame(
        Cluster = names(genelist.ct.entrez),
        mapped_ENTREZ_genes = as.integer(
          official_cluster_input_counts
        ),
        genes_in_official_KEGG_pathways = as.integer(
          official_cluster_overlap_counts
        ),
        overlap_percent = round(
          100 *
            official_cluster_overlap_counts /
            pmax(official_cluster_input_counts, 1),
          2
        ),
        stringsAsFactors = FALSE
      )

      official_input_summary_file <- paste0(
        "KEGG_official_",
        ko.db,
        "_compareCluster_input_summary.tsv"
      )

      write.table(
        official_input_summary,
        official_input_summary_file,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
      )

      message(
        "[INFO] Official KEGG TERM2GENE rows: ",
        nrow(official_term2gene)
      )

      message(
        "[INFO] Unique genes in official KEGG pathways: ",
        length(official_kegg_gene_ids)
      )

      for (summary_index in seq_len(nrow(official_input_summary))) {
        message(
          "[INFO] Official KEGG cluster ",
          official_input_summary$Cluster[summary_index],
          ": ",
          official_input_summary$genes_in_official_KEGG_pathways[
            summary_index
          ],
          "/",
          official_input_summary$mapped_ENTREZ_genes[
            summary_index
          ],
          " genes overlap KEGG pathways (",
          official_input_summary$overlap_percent[
            summary_index
          ],
          "%)"
        )
      }

      genelist.ct.kegg <- genelist.ct.kegg[
        lengths(genelist.ct.kegg) > 0
      ]

      message(
        "[INFO] Official KEGG compareCluster groups with pathway overlap: ",
        length(genelist.ct.kegg)
      )

      if (length(genelist.ct.kegg) > 0) {

        official_comparecluster_attempt <- tryCatch({
          list(
            result = compareCluster(
              genelist.ct.kegg,
              fun = "enricher",
              TERM2GENE = official_term2gene,
              TERM2NAME = official_term2name,
              pAdjustMethod = "BH",
              pvalueCutoff = 0.05,
              qvalueCutoff = 0.2,
              minGSSize = 10,
              maxGSSize = 500
            ),
            error = NULL
          )
        }, error = function(e) {
          list(
            result = NULL,
            error = conditionMessage(e)
          )
        })

        ck <- official_comparecluster_attempt$result
        ck_error <- official_comparecluster_attempt$error

        if (
          !is.null(ck) &&
          nrow(ck@compareClusterResult) > 0
        ) {

          official_kegg_result <- as.data.frame(
            ck@compareClusterResult,
            stringsAsFactors = FALSE
          )

          write.csv(
            official_kegg_result,
            "biological_theme_comparison_KEGG_official_results.csv",
            row.names = FALSE
          )

          message(
            "[OK] Official KEGG compareCluster returned terms: ",
            nrow(official_kegg_result)
          )

          #### BEGIN BLOCK: robust official KEGG compareCluster dotplot ####
          # enrichplot::dotplot() in the current Bioconductor environment can
          # fail when the same pathway description occurs in multiple clusters.
          # Build a faceted ggplot2 figure with an internal unique y-axis key.
          tryCatch({
            required_plot_columns <- c(
              "Cluster", "ID", "Description", "GeneRatio", "Count", "p.adjust"
            )

            missing_plot_columns <- setdiff(
              required_plot_columns,
              colnames(official_kegg_result)
            )

            if (length(missing_plot_columns) > 0) {
              stop(
                "Missing plot columns: ",
                paste(missing_plot_columns, collapse = ", ")
              )
            }

            parse_gene_ratio <- function(x) {
              vapply(
                strsplit(as.character(x), "/", fixed = TRUE),
                function(parts) {
                  if (length(parts) != 2) return(NA_real_)
                  a <- suppressWarnings(as.numeric(parts[1]))
                  b <- suppressWarnings(as.numeric(parts[2]))
                  if (is.na(a) || is.na(b) || b == 0) return(NA_real_)
                  a / b
                },
                numeric(1)
              )
            }

            plot_data <- official_kegg_result
            plot_data$GeneRatio_numeric <- parse_gene_ratio(plot_data$GeneRatio)
            plot_data$Count <- suppressWarnings(as.numeric(plot_data$Count))
            plot_data$p.adjust <- suppressWarnings(as.numeric(plot_data$p.adjust))

            plot_data <- plot_data[
              !is.na(plot_data$Cluster) & plot_data$Cluster != "" &
                !is.na(plot_data$Description) & plot_data$Description != "" &
                !is.na(plot_data$GeneRatio_numeric) &
                is.finite(plot_data$GeneRatio_numeric) &
                !is.na(plot_data$Count) & is.finite(plot_data$Count) &
                !is.na(plot_data$p.adjust) & is.finite(plot_data$p.adjust),
              ,
              drop = FALSE
            ]

            if (nrow(plot_data) == 0) {
              stop("No complete official KEGG rows were available for plotting.")
            }

            top_n_per_cluster <- 20L
            split_data <- split(plot_data, plot_data$Cluster, drop = TRUE)
            split_data <- lapply(split_data, function(x) {
              order_columns <- if ("pvalue" %in% colnames(x)) {
                order(x$p.adjust, x$pvalue, na.last = TRUE)
              } else {
                order(x$p.adjust, na.last = TRUE)
              }
              head(x[order_columns, , drop = FALSE], top_n_per_cluster)
            })

            plot_data <- do.call(rbind, split_data)
            rownames(plot_data) <- NULL

            positive_p <- plot_data$p.adjust[plot_data$p.adjust > 0]
            min_positive <- if (length(positive_p) > 0) {
              min(positive_p, na.rm = TRUE)
            } else {
              .Machine$double.xmin
            }

            plot_data$p_adjust_for_plot <- pmax(
              plot_data$p.adjust,
              min_positive * 0.1
            )
            plot_data$minus_log10_p_adjust <- -log10(plot_data$p_adjust_for_plot)
            plot_data$Description_display <- stringr::str_trunc(
              plot_data$Description,
              60
            )

            plot_data$plot_key <- paste(
              plot_data$Cluster,
              plot_data$ID,
              seq_len(nrow(plot_data)),
              sep = "___"
            )
            plot_data$plot_key <- factor(
              plot_data$plot_key,
              levels = rev(plot_data$plot_key)
            )

            label_map <- setNames(
              plot_data$Description_display,
              as.character(plot_data$plot_key)
            )

            write.csv(
              plot_data,
              "biological_theme_comparison_KEGG_official_plot_data.csv",
              row.names = FALSE
            )

            official_kegg_dotplot <- ggplot(
              plot_data,
              aes(
                x = GeneRatio_numeric,
                y = plot_key,
                size = Count,
                color = minus_log10_p_adjust
              )
            ) +
              geom_point(alpha = 0.85) +
              facet_grid(
                Cluster ~ .,
                scales = "free_y",
                space = "free_y",
                switch = "y"
              ) +
              scale_y_discrete(
                labels = function(keys) unname(label_map[as.character(keys)])
              ) +
              scale_color_gradient(low = "blue", high = "red") +
              labs(
                x = "Gene ratio",
                y = NULL,
                size = "Gene count",
                color = expression(-log[10]("adjusted p-value")),
                title = "Official KEGG pathway enrichment",
                subtitle = paste0("Top ", top_n_per_cluster, " pathways per cluster")
              ) +
              theme_bw() +
              theme(
                strip.placement = "outside",
                strip.text.y.left = element_text(angle = 0, face = "bold"),
                axis.text.y = element_text(size = 8),
                plot.title = element_text(hjust = 0.5),
                plot.subtitle = element_text(hjust = 0.5)
              )

            plot_height <- max(
              6,
              0.32 * nrow(plot_data) +
                1.2 * length(unique(plot_data$Cluster))
            )

            ggsave(
              "biological_theme_comparison_KEGG_official.pdf",
              plot = official_kegg_dotplot,
              limitsize = FALSE,
              height = plot_height,
              width = 11
            )

            message(
              "[OK] Official KEGG compareCluster plot saved: ",
              "biological_theme_comparison_KEGG_official.pdf"
            )

          }, error = function(e) {
            plot_error <- conditionMessage(e)
            message(
              "[WARNING] Official KEGG compareCluster enrichment succeeded, ",
              "but its plot was skipped: ",
              plot_error
            )
            writeLines(
              c(
                "Official KEGG compareCluster enrichment succeeded, but its plot failed.",
                paste0("Result rows saved in CSV: ", nrow(official_kegg_result)),
                paste0("Plot error: ", plot_error)
              ),
              "biological_theme_comparison_KEGG_official_plot_skipped.txt"
            )
          })
          #### END BLOCK: robust official KEGG compareCluster dotplot ####

        } else {

          if (is.null(ck_error)) {
            writeLines(
              c(
                "No significant official KEGG enrichment was found by compareCluster.",
                "This is not a pipeline failure.",
                "",
                paste0("Organism code: ", ko.db),
                paste0(
                  "Input clusters with KEGG pathway overlap: ",
                  length(genelist.ct.kegg)
                ),
                paste0(
                  "Official TERM2GENE rows: ",
                  nrow(official_term2gene)
                ),
                paste0(
                  "Unique official KEGG pathway genes: ",
                  length(official_kegg_gene_ids)
                ),
                "",
                "Settings:",
                "  method = compareCluster + enricher",
                "  pAdjustMethod = BH",
                "  pvalueCutoff = 0.05",
                "  qvalueCutoff = 0.2",
                "  minGSSize = 10",
                "  maxGSSize = 500"
              ),
              "No_KEGG_official_compareCluster_result.txt"
            )

            message(
              "[INFO] Official KEGG compareCluster completed, ",
              "but no term passed the significance cutoffs."
            )

          } else {

            writeLines(
              c(
                "Official KEGG compareCluster failed after building a valid",
                "official KEGG TERM2GENE table.",
                "",
                paste0("Organism code: ", ko.db),
                paste0(
                  "Input clusters with KEGG pathway overlap: ",
                  length(genelist.ct.kegg)
                ),
                paste0(
                  "Official TERM2GENE rows: ",
                  nrow(official_term2gene)
                ),
                paste0("Error: ", ck_error)
              ),
              "No_KEGG_official_compareCluster_result.txt"
            )

            message(
              "[WARNING] Official KEGG compareCluster failed: ",
              ck_error
            )
          }
        }

      } else {

        writeLines(
          c(
            "Official KEGG compareCluster was skipped because none of the",
            "mapped ENTREZ clusters overlapped organism-specific official",
            "KEGG pathway genes.",
            "",
            paste0("Organism code: ", ko.db),
            paste0(
              "Mapped ENTREZ clusters: ",
              length(genelist.ct.entrez)
            ),
            paste0(
              "Official TERM2GENE rows: ",
              nrow(official_term2gene)
            ),
            paste0(
              "Unique official KEGG pathway genes: ",
              length(official_kegg_gene_ids)
            )
          ),
          "No_KEGG_official_compareCluster_result.txt"
        )
      }
    }

  } else {

    writeLines(
      c(
        "Official KEGG compareCluster was skipped because no input",
        "cluster could be mapped to ENTREZ.",
        "",
        paste0("HOST_GENE_KEYTYPE: ", HOST_GENE_KEYTYPE),
        paste0("Organism code: ", ko.db)
      ),
      "No_KEGG_official_compareCluster_result.txt"
    )
  }

} else {

  writeLines(
    c(
      "Official KEGG compareCluster was skipped because no KEGG",
      "organism code was found.",
      "",
      paste0("HOST_GENE_KEYTYPE: ", HOST_GENE_KEYTYPE),
      paste0("Organism code: ", ko.db)
    ),
    "No_KEGG_official_compareCluster_result.txt"
  )
}

#### END BLOCK: official KEGG TERM2GENE compareCluster ####

#### END BLOCK: dual KEGG enrichment comparison ####
  }
}


