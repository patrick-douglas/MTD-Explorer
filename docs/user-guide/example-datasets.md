# Example datasets

MTD Explorer is illustrated with two independent public RNA-seq datasets. The cases demonstrate complementary workflows and the generation of native software outputs, not a direct comparison between host species or a new biological investigation.

| Feature | [*Biomphalaria glabrata*](example-dataset.md) | [*Danio rerio*](danio-rerio.md) |
|---|---|---|
| BioProject | PRJNA1306560 | PRJNA1005695 |
| Host TaxID | 6526 | 7955 |
| Selected libraries | 15 | 6 |
| Library layout | Paired-end | Paired-end |
| Analytical groups | 3 | 1 |
| Analysis mode | Comparative | Exploratory |
| Host alignment | Magic-BLAST | HISAT2 |
| Bracken distribution | 75 bp | 150 bp |
| Taxonomic profiling and TPI | Yes | Yes |
| Host expression and functional profiling | Yes | Yes |
| Differential expression and enrichment | Yes | Not requested |
| ssGSEA and HAllA | Yes | Skipped |

## Comparative example

[Open the *Biomphalaria glabrata* case study](example-dataset.md). The existing guide describes the three-group design and the comparative output examples used throughout the manual. Its historical benchmark and figures retain their recorded provenance, including the separately regenerated TPI figures.

## Exploratory example

[Open the *Danio rerio* case study](danio-rerio.md). The single-group dataset demonstrates a full contrast-free run, automated custom-host reference preparation, an alternative Bracken read-length distribution, native TPI and descriptive expression/functional outputs. Source tables and reproduction settings are available on that page.

See [Pipeline benchmarking](benchmarking.md) for the independently measured performance of each case. The runs differ in dataset size, references, alignment backend, options and measurement methods; their durations should not be interpreted as a controlled speed comparison.
