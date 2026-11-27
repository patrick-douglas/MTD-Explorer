# Supplementary Tables S1–S8

These tables accompany the MTD Explorer analysis of 15 *Biomphalaria glabrata* RNA-seq libraries ([PRJNA1306560](https://www.ncbi.nlm.nih.gov/bioproject/1306560)). The local benchmark used [MTD Explorer revision 04d61d4](https://github.com/patrick-douglas/MTD-Explorer/tree/04d61d412955ac7b9142ca8ebcc4d2edeb6bac5d).

The [Excel workbook](MTD_Explorer_Supplementary_Tables.xlsx) and [21 CSV tables](CSV/) contain 12,465 data rows. Use the [index](Supplementary_Table_Index.csv) to locate a worksheet or CSV and the [legends](Supplementary_Table_Legends.txt) for definitions, units and interpretation.

| Table | Contents |
| --- | --- |
| S1 | Original MTD and MTD Explorer feature comparison at the inspected revisions |
| S2 | Software identities, versions and execution environments |
| S3A–C | Configured hosts, reference resources and host functional annotation coverage |
| S4A–E | Run configuration, parameter records, SRA run aliases, FASTQ identities and resource signatures |
| S5A–C | Durations for the 16 recorded pipeline intervals and computational resource measurements |
| S6A–C | Output families, file types and the complete inventory of 10,821 retained files |
| S7A–B | Taxon Prioritization Index (TPI) ranks and components for 425 features, in relative abundance and estimated reads |
| S8A–C | The 24 reported HAllA blocks, 49 member pairs and expression and annotation records for their 35 host genes |

For the TPI weight and resampling analyses, see [Supplementary Table S9](../importance-index-validation/tables/Supplementary_Table_S9.xlsx) and the [reproducibility files](../importance-index-validation/README.md).

## Reading and tracing the data

The CSV files are UTF-8 and retain the precision of the tabular export. Import accessions, gene identifiers, checksums, versions and timestamps as text when reading them into spreadsheet software. S5A includes workbook formulas for converting elapsed time to hours and percentages. Its CSV contains the numerical results.

The [companion records](Companion_records/) list the 281 retained directories and counts of completion markers for each of six stages. Paths in S6 are relative to the output directory on the benchmark machine. Placeholders such as `<INPUT_DIR>`, `<OUTPUT_DIR>` and `<MTD_ROOT>` in S4 represent the corresponding original run locations.

The [source catalog](Source_Catalog.csv) relates the source identifiers in the legends to repository revisions, public archive records and records from the reference run. Entries labelled “Reference-run record” or “Record or verification” identify source files used during preparation. They are historical identifiers, not paths relative to this folder. The [provenance record](Provenance_and_Verification.json) identifies the archived source data, row counts and checks performed on the original export. The workbook and CSVs use the current manuscript figure numbering. The S2 and S4A CSVs and the source catalog contain updated figure references or audit identifiers. The other 19 table CSVs, both companion inventories and the table index retain their original byte identities.
