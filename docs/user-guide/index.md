# User guide

This section explains how to prepare, run, inspect, and interpret an
[MTD Explorer][mtd-explorer] analysis.

Use this page as a map of the documentation.

## Start here

<div class="grid cards" markdown>

- **Prepare an analysis**

    Set up input FASTQ files, samplesheets, metadata, and analysis modes.

    [Input files](input-files.md) · [Analysis modes](analysis-modes.md)

- **Run from the command line**

    Review the main command-line options used to control host, microbiome,
    functional, and integration analyses.

    [Command-line reference](command-line.md)

- **Understand the output folder**

    Learn how the major result folders are organized before opening the more
    detailed output pages.

    [Output files](output-files.md)

- **Check reproducibility**

    Inspect run parameters, software paths, database settings, and analysis
    metadata.

    [Methods and reproducibility outputs](methods-reproducibility-outputs.md)

</div>

## Output guides

<div class="grid cards" markdown>

- **Taxonomic exploratory outputs**

    Exploratory taxonomic summaries, read composition, abundance landscapes,
    diversity plots, overlap diagrams, and matrix QC.

    [Open guide](taxonomic-exploratory-outputs.md)

- **Taxonomic visualizations**

    Interactive [Krona][krona] plots and [GraPhlAn][graphlan] cladograms for
    inspecting taxonomic composition.

    [Open guide](taxonomic-visualizations.md)

- **Host expression outputs**

    Host count matrices, PCA, heatmaps, differential expression, volcano plots,
    and host functional enrichment.

    [Open guide](host-expression-outputs.md)

- **Microbiome comparison outputs**

    Group comparisons for non-host taxa, including heatmaps, alpha diversity,
    ANOSIM, ANCOM-BC, and microbiome differential abundance.

    [Open guide](microbiome-comparison-outputs.md)

- **Functional profiling outputs**

    Non-host functional profiles derived from [HUMAnN][humann], including GO
    and KEGG heatmaps and PCA plots.

    [Open guide](functional-profiling-outputs.md)

- **ssGSEA outputs**

    Sample-level host gene-set activity summaries, PCA, correlation heatmaps,
    variable gene-set heatmaps, and differential boxplots.

    [Open guide](ssgsea-outputs.md)

- **HAllA integration outputs**

    Host–microbiome and host–function association summaries produced with
    [HAllA][halla].

    [Open guide](halla-integration-outputs.md)

- **Academic references**

    Citation-ready references for the original MTD pipeline and the main
    software packages used by MTD Explorer.

    [Open references](../reference/index.md)

</div>

## Recommended reading order

For a first complete analysis, read the documentation in this order:

```text
1. Installation
2. Verify installation
3. Input files
4. Analysis modes
5. Command-line reference
6. Output files
7. Taxonomic exploratory outputs
8. Taxonomic visualizations
9. Host expression outputs
10. Microbiome comparison outputs
11. Functional profiling outputs
12. ssGSEA outputs
13. HAllA integration outputs
14. Methods and reproducibility outputs
15. Academic references
```

## Result layers

MTD Explorer produces several result layers from the same RNA-seq dataset.

| Result layer | Main folder or file | Main guide |
| --- | --- | --- |
| Input setup | `samplesheet.txt`, FASTQ files, metadata | [Input files](input-files.md) |
| Analysis settings | command-line options | [Command-line reference](command-line.md) |
| Output structure | main output directory | [Output files](output-files.md) |
| Taxonomic exploration | `exploratory/taxonomy/` | [Taxonomic exploratory outputs](taxonomic-exploratory-outputs.md) |
| Taxonomic visualization | `krona/`, `graphlan/` | [Taxonomic visualizations](taxonomic-visualizations.md) |
| Host expression | `Host_DEG/` | [Host expression outputs](host-expression-outputs.md) |
| Microbiome comparison | `Nonhost_DEG/` | [Microbiome comparison outputs](microbiome-comparison-outputs.md) |
| Functional profiling | `hmn_genefamily_abundance_files/` | [Functional profiling outputs](functional-profiling-outputs.md) |
| Gene-set activity | `ssGSEA/` | [ssGSEA outputs](ssgsea-outputs.md) |
| Multi-layer integration | `halla/` | [HAllA integration outputs](halla-integration-outputs.md) |
| Reproducibility | `methods/` | [Methods and reproducibility outputs](methods-reproducibility-outputs.md) |

## Choosing the right page

| Question | Start with |
| --- | --- |
| How do I format my samplesheet? | [Input files](input-files.md) |
| Should I run comparison or exploratory mode? | [Analysis modes](analysis-modes.md) |
| Which command-line option controls trimming, Kraken2, Bracken, or read layout? | [Command-line reference](command-line.md) |
| Where are the main results saved? | [Output files](output-files.md) |
| How do I inspect detected microbiome composition? | [Taxonomic exploratory outputs](taxonomic-exploratory-outputs.md) |
| Where are Krona and GraPhlAn plots? | [Taxonomic visualizations](taxonomic-visualizations.md) |
| Where are host DEG plots and tables? | [Host expression outputs](host-expression-outputs.md) |
| Where are ANOSIM and ANCOM-BC outputs? | [Microbiome comparison outputs](microbiome-comparison-outputs.md) |
| Where are GO/KEGG functional profiles? | [Functional profiling outputs](functional-profiling-outputs.md) |
| Where are host gene-set activity plots? | [ssGSEA outputs](ssgsea-outputs.md) |
| Where are host–microbiome association outputs? | [HAllA integration outputs](halla-integration-outputs.md) |
| Which tools should I cite? | [Academic references](../reference/index.md) |

## Interpretation principles

MTD Explorer outputs are designed to be inspected together.

A taxonomic signal should be interpreted with read composition, abundance
tables, database settings, and sample metadata.

A host expression signal should be interpreted with the count matrix,
normalization, PCA, differential-expression tables, and biological context.

Functional and association outputs should be treated as exploratory summaries
unless supported by the corresponding tables, statistical results, and study
design.

## Related getting-started pages

- [System requirements](../getting-started/requirements.md)
- [Installation](../getting-started/installation.md)
- [Verify installation](../getting-started/verify-installation.md)
- [Quick start](../getting-started/quick-start.md)

[mtd-explorer]: https://github.com/patrick-douglas/MTD-Explorer
[krona]: https://github.com/marbl/Krona/wiki
[graphlan]: https://github.com/biobakery/graphlan
[humann]: https://github.com/biobakery/humann
[halla]: https://github.com/biobakery/halla
