# MTD Explorer manuscript supplements

This directory accompanies the joint host transcriptome and microbiome analysis of 15 *Biomphalaria glabrata* RNA-seq libraries from [BioProject PRJNA1306560](https://www.ncbi.nlm.nih.gov/bioproject/1306560). The local benchmark used [MTD Explorer revision 04d61d4](https://github.com/patrick-douglas/MTD-Explorer/tree/04d61d412955ac7b9142ca8ebcc4d2edeb6bac5d). The subsequent Taxon Prioritization Index (TPI) robustness analysis used the saved taxonomic profiles. Its computation time is outside the reported benchmark.

| Material | Files | Contents |
| --- | --- | --- |
| Supplementary Tables S1–S8 | [Excel workbook](supplementary-tables/MTD_Explorer_Supplementary_Tables.xlsx) and [CSV tables](supplementary-tables/CSV/) | Workflow comparison, environments, reference resources, study inputs, performance, output inventory, taxon rankings and host–microbiome associations. Each data table has an accompanying [legend](supplementary-tables/Supplementary_Table_Legends.txt). |
| Supplementary Table S9 | [Excel workbook](importance-index-validation/tables/Supplementary_Table_S9.xlsx) and [analysis results](importance-index-validation/results/) | Weight sensitivity, detection and taxon inclusion criteria, single-library omissions and stratified bootstrap results. |
| Supplementary Figure S1 | [Figure](supplementary-figure-S1/supplementary_figure_1_reference.svg) and [caption](supplementary-figure-S1/supplementary_figure_1_caption.txt) | How the host reference supplies the annotation database and the GO gene-set collection. Node and connection tables describe the diagram. The evidence table links to the code. |
| Supplementary Figure S2 | [Figure](importance-index-validation/figures/Supplementary_Figure_S2_importance_robustness.pdf) and [caption](importance-index-validation/figures/Supplementary_Figure_S2_caption.txt) | Sensitivity of taxon rankings to weights and library resampling. The same folder contains the numerical data for each panel. |
| TPI analysis | [Instructions](importance-index-validation/README.md), [input matrices](importance-index-validation/data/), [protocol](importance-index-validation/protocol.json) and [scripts](importance-index-validation/src/) | Recompute and audit the rankings from the two observed abundance representations. Stored results include the bootstrap draws and replicate-level ranks. |
| Benchmark measurements | [Run records](benchmark-records/) | Sample sheet, stage timestamps and 39,020 resource observations underlying the reported run. |

The original sequencing reads are indexed by SRA study SRP609141. Supplementary Table S4C maps SRA run accessions to the library identifiers used here. The size reported in S6 is the logical size of the inventoried run outputs. This supplement contains the file inventory.

To verify the downloaded files on macOS or Linux, run the following command from this directory before regenerating outputs:

```sh
shasum -a 256 -c SHA256SUMS.txt
```
