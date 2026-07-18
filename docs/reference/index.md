---
hide:
  - navigation
---

# Academic references

This page lists academic and software references for tools used or supported by
[MTD Explorer][mtd-explorer].

The list is intended to help users cite the original methods, databases, and
software components used in an analysis.

For manuscripts, always check the citation policy of each tool and database,
and cite the exact software versions used when required.

## MTD and MTD Explorer

The original MTD pipeline should be cited whenever MTD Explorer is used, because
MTD Explorer is built upon and extends the original MTD workflow.

- Wu, F., Liu, Y.-Z., & Ling, B. (2022). **MTD: a unique pipeline for host and
  meta-transcriptome joint and integrative analyses of RNA-seq data**.
  *Briefings in Bioinformatics*, 23(3), bbac111.
  https://doi.org/10.1093/bib/bbac111

MTD Explorer itself should be cited according to the citation information
provided in the MTD Explorer repository, release, manuscript, or Zenodo archive
when available.

## Workflow execution and command-line utilities

### GNU Parallel

[GNU Parallel][gnu-parallel] is used to run jobs in parallel.

GNU Parallel explicitly asks users to cite it in publications. The most accurate
citation for the installed version can be obtained with:

```bash
parallel --citation
```

Recommended reference:

- Tange, O. (2011). **GNU Parallel: The Command-Line Power Tool**. *;login: The
  USENIX Magazine*, 36(1), 42–47.

Some GNU Parallel releases may recommend citing a version-specific Zenodo record.
Use `parallel --citation` to obtain the preferred citation for the version used
in a specific analysis.

### R

- R Core Team. (2024). **R: A language and environment for statistical
  computing**. R Foundation for Statistical Computing, Vienna, Austria.
  https://www.R-project.org/

### Python

- Van Rossum, G., & Drake, F. L. (2009). **Python 3 Reference Manual**.
  CreateSpace, Scotts Valley, CA.

## Read preprocessing and quality control

### fastp

- Chen, S., Zhou, Y., Chen, Y., & Gu, J. (2018). **fastp: an ultra-fast
  all-in-one FASTQ preprocessor**. *Bioinformatics*, 34(17), i884–i890.
  https://doi.org/10.1093/bioinformatics/bty560

## Taxonomic classification and abundance estimation

### Kraken 2

- Wood, D. E., Lu, J., & Langmead, B. (2019). **Improved metagenomic analysis
  with Kraken 2**. *Genome Biology*, 20, 257.
  https://doi.org/10.1186/s13059-019-1891-0

### Bracken

- Lu, J., Breitwieser, F. P., Thielen, P., & Salzberg, S. L. (2017).
  **Bracken: estimating species abundance in metagenomics data**.
  *PeerJ Computer Science*, 3, e104.
  https://doi.org/10.7717/peerj-cs.104

### KrakenTools

- Lu, J. (2022). **KrakenTools**. GitHub repository.
  https://github.com/jenniferlu717/KrakenTools

When KrakenTools is used for report manipulation or contaminant removal,
also cite [Kraken 2][kraken2] and [Bracken][bracken] when those outputs are
used downstream.

## Taxonomic visualization

### Krona

- Ondov, B. D., Bergman, N. H., & Phillippy, A. M. (2011). **Interactive
  metagenomic visualization in a Web browser**. *BMC Bioinformatics*, 12, 385.
  https://doi.org/10.1186/1471-2105-12-385

### GraPhlAn

- Asnicar, F., Weingart, G., Tickle, T. L., Huttenhower, C., & Segata, N.
  (2015). **Compact graphical representation of phylogenetic data and metadata
  with GraPhlAn**. *PeerJ*, 3, e1029.
  https://doi.org/10.7717/peerj.1029

## Host read alignment and gene counting

### HISAT2

- Kim, D., Paggi, J. M., Park, C., Bennett, C., & Salzberg, S. L. (2019).
  **Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype**.
  *Nature Biotechnology*, 37, 907–915.
  https://doi.org/10.1038/s41587-019-0201-4

### Magic-BLAST

- Boratyn, G. M., Thierry-Mieg, J., Thierry-Mieg, D., Busby, B., & Madden, T. L.
  (2019). **Magic-BLAST, an accurate RNA-seq aligner for long and short reads**.
  *BMC Bioinformatics*, 20, 405.
  https://doi.org/10.1186/s12859-019-2996-x

### featureCounts / Subread

- Liao, Y., Smyth, G. K., & Shi, W. (2014). **featureCounts: an efficient
  general purpose program for assigning sequence reads to genomic features**.
  *Bioinformatics*, 30(7), 923–930.
  https://doi.org/10.1093/bioinformatics/btt656

### SAMtools

- Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N.,
  Marth, G., Abecasis, G., & Durbin, R. (2009). **The Sequence
  Alignment/Map format and SAMtools**. *Bioinformatics*, 25(16), 2078–2079.
  https://doi.org/10.1093/bioinformatics/btp352

## Differential expression and visualization

### DESeq2

- Love, M. I., Huber, W., & Anders, S. (2014). **Moderated estimation of fold
  change and dispersion for RNA-seq data with DESeq2**. *Genome Biology*,
  15, 550. https://doi.org/10.1186/s13059-014-0550-8

### EnhancedVolcano

- Blighe, K., Rana, S., & Lewis, M. (2024). **EnhancedVolcano: Publication-ready
  volcano plots with enhanced colouring and labeling**. Bioconductor package.
  https://bioconductor.org/packages/EnhancedVolcano/

## Functional enrichment and gene-set activity

### Gene Ontology

- Ashburner, M., Ball, C. A., Blake, J. A., Botstein, D., Butler, H.,
  Cherry, J. M., Davis, A. P., Dolinski, K., Dwight, S. S., Eppig, J. T.,
  Harris, M. A., Hill, D. P., Issel-Tarver, L., Kasarskis, A., Lewis, S.,
  Matese, J. C., Richardson, J. E., Ringwald, M., Rubin, G. M., &
  Sherlock, G. (2000). **Gene Ontology: tool for the unification of biology**.
  *Nature Genetics*, 25, 25–29. https://doi.org/10.1038/75556

- The Gene Ontology Consortium. (2023). **The Gene Ontology knowledgebase in
  2023**. *Genetics*, 224(1), iyad031.
  https://doi.org/10.1093/genetics/iyad031

### KEGG

- Kanehisa, M., & Goto, S. (2000). **KEGG: Kyoto Encyclopedia of Genes and
  Genomes**. *Nucleic Acids Research*, 28(1), 27–30.
  https://doi.org/10.1093/nar/28.1.27

- Kanehisa, M., Furumichi, M., Sato, Y., Kawashima, M., & Ishiguro-Watanabe, M.
  (2023). **KEGG for taxonomy-based analysis of pathways and genomes**.
  *Nucleic Acids Research*, 51(D1), D587–D592.
  https://doi.org/10.1093/nar/gkac963

### clusterProfiler

- Yu, G., Wang, L.-G., Han, Y., & He, Q.-Y. (2012). **clusterProfiler: an R
  package for comparing biological themes among gene clusters**.
  *OMICS: A Journal of Integrative Biology*, 16(5), 284–287.
  https://doi.org/10.1089/omi.2011.0118

### GSEA and ssGSEA

- Subramanian, A., Tamayo, P., Mootha, V. K., Mukherjee, S., Ebert, B. L.,
  Gillette, M. A., Paulovich, A., Pomeroy, S. L., Golub, T. R.,
  Lander, E. S., & Mesirov, J. P. (2005). **Gene set enrichment analysis:
  a knowledge-based approach for interpreting genome-wide expression profiles**.
  *Proceedings of the National Academy of Sciences*, 102(43), 15545–15550.
  https://doi.org/10.1073/pnas.0506580102

- Barbie, D. A., Tamayo, P., Boehm, J. S., Kim, S. Y., Moody, S. E.,
  Dunn, I. F., Schinzel, A. C., Sandy, P., Meylan, E., Scholl, C.,
  Fröhling, S., Chan, E. M., Sos, M. L., Michel, K., Mermel, C.,
  Silver, S. J., Weir, B. A., Reiling, J. H., Sheng, Q., ... Gilliland, D. G.
  (2009). **Systematic RNA interference reveals that oncogenic KRAS-driven
  cancers require TBK1**. *Nature*, 462, 108–112.
  https://doi.org/10.1038/nature08460

## Microbial functional profiling

### HUMAnN / bioBakery

- Franzosa, E. A., McIver, L. J., Rahnavard, G., Thompson, L. R.,
  Schirmer, M., Weingart, G., Lipson, K. S., Knight, R., Caporaso, J. G.,
  Segata, N., & Huttenhower, C. (2018). **Species-level functional profiling
  of metagenomes and metatranscriptomes**. *Nature Methods*, 15, 962–968.
  https://doi.org/10.1038/s41592-018-0176-y

- Beghini, F., McIver, L. J., Blanco-Míguez, A., Dubois, L., Asnicar, F.,
  Maharjan, S., Mailyan, A., Manghi, P., Scholz, M., Thomas, A. M.,
  Valles-Colomer, M., Weingart, G., Zhang, Y., Zolfo, M., Huttenhower, C.,
  Franzosa, E. A., & Segata, N. (2021). **Integrating taxonomic, functional,
  and strain-level profiling of diverse microbial communities with bioBakery 3**.
  *eLife*, 10, e65088. https://doi.org/10.7554/eLife.65088

## Microbiome statistics and association models

### ANCOM-BC

- Lin, H., & Peddada, S. D. (2020). **Analysis of compositions of microbiomes
  with bias correction**. *Nature Communications*, 11, 3514.
  https://doi.org/10.1038/s41467-020-17041-7

### MaAsLin2

- Mallick, H., Rahnavard, A., McIver, L. J., Ma, S., Zhang, Y.,
  Nguyen, L. H., Tickle, T. L., Weingart, G., Ren, B., Schwager, E. H.,
  Chatterjee, S., Thompson, K. N., Wilkinson, J. E., Subramanian, A.,
  Lu, Y., Waldron, L., Paulson, J. N., Franzosa, E. A., Bravo, H. C., &
  Huttenhower, C. (2021). **Multivariable association discovery in
  population-scale meta-omics studies**. *PLOS Computational Biology*, 17(11),
  e1009442. https://doi.org/10.1371/journal.pcbi.1009442

### vegan / ANOSIM

- Oksanen, J., Simpson, G. L., Blanchet, F. G., Kindt, R., Legendre, P.,
  Minchin, P. R., O'Hara, R. B., Solymos, P., Stevens, M. H. H.,
  Szoecs, E., Wagner, H., Barbour, M., Bedward, M., Bolker, B.,
  Borcard, D., Carvalho, G., Chirico, M., De Caceres, M., Durand, S.,
  ... Weedon, J. (2024). **vegan: Community Ecology Package**. R package.
  https://CRAN.R-project.org/package=vegan

## Integration analysis

### HAllA

- Ghazi, A. R., Franzosa, E. A., Buhimschi, C. S., Annavajhala, M. K.,
  Hyman, R. W., Huh, J. W., Ravel, J., Buhimschi, I. A., & Huttenhower, C.
  (2022). **High-sensitivity pattern discovery in large, paired multi-omic
  datasets**. *Bioinformatics*, 38(Supplement_1), i378–i385.
  https://doi.org/10.1093/bioinformatics/btac234

## Plotting and R ecosystem

### ggplot2

- Wickham, H. (2016). **ggplot2: Elegant Graphics for Data Analysis**.
  Springer-Verlag New York. https://ggplot2.tidyverse.org

### pheatmap

- Kolde, R. (2019). **pheatmap: Pretty Heatmaps**. R package.
  https://CRAN.R-project.org/package=pheatmap

### RColorBrewer

- Neuwirth, E. (2022). **RColorBrewer: ColorBrewer Palettes**. R package.
  https://CRAN.R-project.org/package=RColorBrewer

## Notes for manuscripts

When preparing a manuscript based on MTD Explorer results, cite at minimum:

1. the original [MTD][mtd-explorer] paper;
2. MTD Explorer, if a DOI or manuscript is available;
3. the main tools used in the specific analysis;
4. the databases used for taxonomic, functional, host, GO, and KEGG annotation;
5. [GNU Parallel][gnu-parallel] when used to process data for publication.

For exact software-specific citations in R, users can also run:

```r
citation()
citation("DESeq2")
citation("ANCOMBC")
citation("Maaslin2")
citation("vegan")
citation("ggplot2")
```

[mtd-explorer]: https://github.com/patrick-douglas/MTD
[gnu-parallel]: https://www.gnu.org/software/parallel/
[kraken2]: https://ccb.jhu.edu/software/kraken2/index.shtml
[bracken]: https://github.com/jenniferlu717/Bracken
