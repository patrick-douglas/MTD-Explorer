# MTD Explorer manuscript supplements

This directory accompanies the joint host transcriptome and microbiome analysis of 15 *Biomphalaria glabrata* RNA-seq libraries from [BioProject PRJNA1306560](https://www.ncbi.nlm.nih.gov/bioproject/1306560). The local benchmark used [MTD Explorer revision 04d61d4](https://github.com/patrick-douglas/MTD-Explorer/tree/04d61d412955ac7b9142ca8ebcc4d2edeb6bac5d). The subsequent Taxon Prioritization Index (TPI) robustness analysis used the saved taxonomic profiles. Its computation time is outside the reported benchmark.

| Material | Files | Contents |
| --- | --- | --- |
| Figure 6 source data | [6A](figure-source-data/figure_6_panel_A_data.csv) · [6B](figure-source-data/figure_6_panel_B_data.csv) · [6C](figure-source-data/figure_6_panel_C_data.csv) | Infected versus uninfected libraries (five per group). Panels A and B list all 21,277 host genes and 392 microbial features plotted in the volcano plots, with log2 fold changes, BH-adjusted p values and label ranks. Panel C gives normalized ssGSEA scores for the four displayed GO sets in each of the ten libraries. |
| Figure 7 source data | [7A](figure-source-data/figure_7_panel_A_data.csv) · [7B](figure-source-data/figure_7_panel_B_data.csv) | Panel A gives Spearman correlations for all 945 displayed taxon–gene pairs across 15 libraries, with HAllA block membership and BH-adjusted p values from the full set of 8,839,575 tests. Panel B gives infected-versus-uninfected log2 fold changes and standard errors for the eight differentially expressed genes highlighted in the figure. |
| Supplementary Tables S1–S8 | [Excel workbook](supplementary-tables/MTD_Explorer_Supplementary_Tables.xlsx) and [CSV tables](supplementary-tables/CSV/) | Workflow comparison, environments, reference resources, study inputs, performance, output inventory, taxon rankings and host–microbiome associations. Each data table has an accompanying [legend](supplementary-tables/Supplementary_Table_Legends.txt). |
| Supplementary Table S9 | [Excel workbook](importance-index-validation/tables/Supplementary_Table_S9.xlsx) and [analysis results](importance-index-validation/results/) | Weight sensitivity, detection and taxon inclusion criteria, single-library omissions and stratified bootstrap results. |
| Supplementary Figure S1 | [Figure](supplementary-figure-S1/supplementary_figure_1_reference.svg) and [caption](supplementary-figure-S1/supplementary_figure_1_caption.txt) | How the host reference supplies the annotation database and the GO gene-set collection. Node and connection tables describe the diagram. The evidence table links to the code. |
| Supplementary Figure S2 | [Figure](importance-index-validation/figures/Supplementary_Figure_S2_importance_robustness.pdf) and [caption](importance-index-validation/figures/Supplementary_Figure_S2_caption.txt) | Sensitivity of taxon rankings to weights and library resampling. The same folder contains the numerical data for each panel. |
| Supplementary Figure S3 | [Figure](supplementary-figure-S3/supplementary_figure_3_performance.svg) and [caption](supplementary-figure-S3/supplementary_figure_3_caption.md) | Runtime by pipeline interval ([panel A](supplementary-figure-S3/supplementary_figure_3_panel_A_data.csv)) and sampled memory and CPU activity over time ([panels B and C](supplementary-figure-S3/supplementary_figure_3_panels_BC_data.csv)). |
| TPI analysis | [Instructions](importance-index-validation/README.md), [input matrices](importance-index-validation/data/), [protocol](importance-index-validation/protocol.json) and [scripts](importance-index-validation/src/) | Recompute and audit the rankings from the two observed abundance representations. Stored results include the bootstrap draws and replicate-level ranks. |
| Benchmark measurements | [Run records](benchmark-records/) | Sample sheet, stage timestamps and 39,020 resource observations underlying the reported run. |

The original sequencing reads are indexed by SRA study SRP609141. Supplementary Table S4C maps SRA run accessions to the library identifiers used here. The size reported in S6 is the logical size of the inventoried run outputs. This supplement contains the file inventory.

To verify the downloaded files on macOS or Linux, run the following command from this directory before regenerating outputs:

```sh
shasum -a 256 -c SHA256SUMS.txt
```

## Exploratory case study: Danio rerio

The [exploratory case study](exploratory-case-study/README.md) documents an independent six-library dataset (PRJNA1005695) using analytical revision `5be26aa`. It includes source tables, runtime and reference provenance, and repeatability checks. The original comparative materials above are retained unchanged; the exploratory run demonstrates software capabilities rather than biological inference. See the [online example](https://patrick-douglas.github.io/MTD-Explorer/user-guide/danio-rerio/).

## Comparative functional regression: Biomphalaria glabrata

The [comparative functional regression](comparative-functional-regression/README.md) documents the isolated GO/KEGG software test at revisions `af59f09` and `5be26aa`, including the original regression archive, before/after data, native plots, and audit. The original benchmark and historical supplements are preserved; corrected functional results must not be treated as numerically interchangeable with historical results.

## Reference-compatible host annotation comparison

The [host annotation comparison](host-annotation-comparison/README.md) provides the [source measurements for Figure 3](host-annotation-comparison/results/metrics.tsv) and documents the original MTD and MTD-Explorer source audit, fixed-input gene identifier compatibility, independently reproduced GO gene-set membership, the expected legacy ssGSEA no-overlap stop, and the custom OrgDb query. This comparison evaluates automation and identifier compatibility, not biological annotation accuracy.
