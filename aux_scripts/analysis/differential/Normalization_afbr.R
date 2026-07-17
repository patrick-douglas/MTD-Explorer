args = commandArgs(trailingOnly=TRUE) # Passing arguments to an R script from bash/shell command lines

# Firstly, to calculate the normalization of reads
library(DESeq2)
library(tibble)
# Read & preprocess the input file cts before run deseq2 analysis
filename<-basename(args[1])
if (filename %in% c("bracken_species_all","bracken_phylum_all","bracken_genus_all")){
  # read the data file
  #bracken file (eg. bracken_species_all) without quote symbol "\""; mark empty quote as quote=""
  cts<-read.table(args[1],
                  row.names=1, sep="\t",header=T, quote="")
  # Extract columns with count numbers
  cts<-cts[,grepl("*_num", names(cts))]
  # match column name to sample name
  colnames(cts)<-gsub("^Report_|\\.species.bracken_num$|\\.phylum.bracken_num$|.genus.bracken_num$","",colnames(cts))
}

# read the original samplesheet
coldata0 <- read.csv(args[2], header = T, na.strings=c("","NA"))
# extract the contrast information for reference
coldata_vs<- coldata0[c("group1","group2")]
coldata_vs<-coldata_vs[rowSums(is.na(coldata_vs)) == 0,] #remove the NA rows
# read samplesheet as factors (as.is = F) for Deseq2 statistical analysis
coldata_factor <- read.csv(args[2], header = T, as.is = F)
coldata_factor[]<-lapply(coldata_factor, factor)
coldata<-coldata_factor[,1:2]
# update the coldata if metadata is provided
if (length(args) == 4){
  coldata <- read.csv(args[4], header = T, as.is = F)
  coldata[]<-lapply(coldata, factor)
}

if (filename %in% c("bracken_species_all","bracken_phylum_all","bracken_genus_all")){
  files_h <- list.files(path=paste0(dirname(args[1]),"/temp"), pattern="^Report_host_.*\\.txt$", full.names=TRUE, recursive=FALSE)
  # to get a list for host
  transcriptome_size<-c() #generate a empty list
  for (i in files_h){
    t<-read.table(i,sep="\t", quote= "")
    total_reads<-t[1,2] + t[2,2] #get total reads abundance of a sample
    fn<-gsub("^Report_host_|\\.txt$","",basename(i)) #grep the sample name
    transcriptome_size<-c(transcriptome_size,setNames(total_reads,fn)) #add sample name with value to the list lh
  }
  transcriptome_size <- log2(transcriptome_size)-mean(log2(transcriptome_size))
  coldata$order<-1:nrow(coldata)
  coldata<-merge(coldata,as.data.frame(transcriptome_size), by.x="sample_name",by.y="row.names")
  coldata<-coldata[order(coldata$order), ]
  coldata<-subset(coldata, select = -c(order))
  # make cts(count matrix) has consistent order with samplesheet
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
if (length(args) == 4){
  design(dds)<-funNew(names(coldata)[2:ncol(coldata)])
}

# perform the DESeq analysis
dds <- DESeq(dds)

# normalized and transfored reads count
if (dim(results(dds))[1]  < 1000 || min(colSums(cts !=0)) < 1000){
  vsd<-varianceStabilizingTransformation(dds,blind=F) # vatiance stabilizing transformation
} else {
  vsd<-vst(dds,blind=F)
}
normtrans<-assay(vsd)

# norm<-counts(dds,normalized=T)

# normalized reads count with host transcriptome size and with avoiding removing variation associated with the other conditions
if (length(args) == 4){
  mm <- model.matrix(funNew(names(coldata)[2:(ncol(coldata)-1)]), colData(vsd))
} else {
  mm <- model.matrix(funNew(names(coldata)[2]), colData(vsd))
}
norm <- limma::removeBatchEffect(normtrans, vsd$transcriptome_size, design=mm)
# if (length(args) == 4){
#   coldata.n<-coldata
#   coldata.n[]<-lapply(coldata.n, as.numeric)
#   norm <- limma::removeBatchEffect(normtrans, covariates=coldata.n[,2:ncol(coldata.n)])
# }

# ------------------------------------------------------------
# Export transformed matrices
#
# Important:
# These matrices are intended for ordination, clustering,
# heatmaps and other multivariate visualizations.
#
# They must not be written back into Kraken/Bracken reports,
# because VST/batch-corrected values are not read counts and
# may contain negative or non-integer values.
# ------------------------------------------------------------

output_dir <- args[3]

dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

write_taxon_matrix <- function(matrix_object, output_file) {
    matrix_object <- as.matrix(matrix_object)

    if (is.null(rownames(matrix_object))) {
        stop(
            paste(
                "The matrix has no taxon row names:",
                output_file
            )
        )
    }

    if (any(!is.finite(matrix_object))) {
        stop(
            paste(
                "The matrix contains NA, NaN or infinite values:",
                output_file
            )
        )
    }

    output_table <- data.frame(
        taxon = rownames(matrix_object),
        matrix_object,
        check.names = FALSE
    )

    write.table(
        output_table,
        file = output_file,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        col.names = TRUE
    )
}


vst_output <- file.path(
    output_dir,
    paste0(filename, ".vst.tsv")
)

batch_corrected_output <- file.path(
    output_dir,
    paste0(filename, ".vst_batch_corrected.tsv")
)

normalized_counts_output <- file.path(
    output_dir,
    paste0(filename, ".deseq2_normalized_counts.tsv")
)


# Uncorrected VST matrix.
write_taxon_matrix(
    normtrans,
    vst_output
)


# VST matrix with the transcriptome-size/batch component removed.
#
# Negative values are valid here because this is a transformed
# visualization matrix, not a count table.
write_taxon_matrix(
    norm,
    batch_corrected_output
)


# DESeq2 size-factor-normalized counts.
#
# This is exported for descriptive/tabular use, but it is not
# converted into a Kraken hierarchy.
normalized_counts <- counts(
    dds,
    normalized = TRUE
)

write_taxon_matrix(
    normalized_counts,
    normalized_counts_output
)


# Preserve sample metadata used in the transformation.
metadata_output <- file.path(
    output_dir,
    paste0(filename, ".sample_metadata.tsv")
)

metadata_table <- data.frame(
    sample_name = rownames(as.data.frame(colData(vsd))),
    as.data.frame(colData(vsd)),
    check.names = FALSE
)

write.table(
    metadata_table,
    file = metadata_output,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
)


message("Transformation completed.")
message("VST matrix: ", vst_output)
message(
    "Batch-corrected VST matrix: ",
    batch_corrected_output
)
message(
    "DESeq2 normalized counts: ",
    normalized_counts_output
)
message("Sample metadata: ", metadata_output)
