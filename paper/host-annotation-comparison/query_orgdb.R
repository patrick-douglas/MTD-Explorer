args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2)
out <- args[1]
library(AnnotationDbi)
library(org.Bglabrata.eg.db)
db <- org.Bglabrata.eg.db
cat('Package version:', as.character(packageVersion('org.Bglabrata.eg.db')), '\n')
stopifnot(as.character(packageVersion('org.Bglabrata.eg.db'))=='1.0.0')
print(keytypes(db)); print(columns(db))
a <- AnnotationDbi::select(db, keys='BGLAX_026407', keytype='GID', columns=c('GO','EVIDENCE','ONTOLOGY'))
print(a)
stopifnot(nrow(a)==56L, all(a$GID=='BGLAX_026407'), all(!is.na(a$GO)), all(a$EVIDENCE=='IEA'), setequal(a$ONTOLOGY,c('BP','MF','CC')))
write.table(a, file.path(out,'orgdb_example_annotation.tsv'), sep='\t', quote=FALSE,row.names=FALSE)
# Evaluate the original catalog lookup only; do not call the BioMart server.
host_sp <- read.csv(args[2])
lookup <- host_sp[host_sp$Taxon_ID==6526,2]
cat('Original BioMart catalog lookup length:',length(lookup),'\n')
stopifnot(length(lookup)==0L)
writeLines(c('metric\tvalue','orgdb_version\t1.0.0','example_GO_associations\t56','original_biomart_catalog_lookup_length\t0','legacy_ssgsea_exit_status\t1'),file.path(out,'runtime_metrics.tsv'))
print(sessionInfo())
