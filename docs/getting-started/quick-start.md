# Quick start

This guide shows the minimal steps required to run MTD Explorer after a
successful installation.

Before running an analysis, complete:

1. [Installation](installation.md)
2. [Installation verification](verify-installation.md)

!!! important

    Run MTD Explorer only after the installation checker has completed
    successfully.

    ```bash
    bash MTD_check_installation.sh --mode full
    ```

## Basic command structure

The main MTD Explorer entry point is:

```text
MTD_explorer.sh
```

The minimal command structure is:

```bash
bash MTD_explorer.sh \
  --input samplesheet.csv \
  --output MTD_results \
  --hostid TAXID
```

Required arguments:

| Option | Description |
|---|---|
| `--input`, `-i` | Path to the samplesheet CSV file |
| `--output`, `-o` | Output directory |
| `--hostid`, `-h` | Host species NCBI Taxonomy ID used for host annotation and downstream analysis |

Display the current help message with:

```bash
bash MTD_explorer.sh --help
```

## Prepare the working directory

Create a project directory outside the MTD Explorer source-code folder.

Example:

```bash
mkdir -p ~/MTD_projects/test_run
cd ~/MTD_projects/test_run
```

A simple project layout is:

```text
test_run/
├── fastq/
│   ├── sample1.fastq.gz
│   ├── sample2.fastq.gz
│   └── sample3.fastq.gz
└── samplesheet.csv
```

For paired-end data:

```text
test_run/
├── fastq/
│   ├── sample1_R1.fastq.gz
│   ├── sample1_R2.fastq.gz
│   ├── sample2_R1.fastq.gz
│   ├── sample2_R2.fastq.gz
│   ├── sample3_R1.fastq.gz
│   └── sample3_R2.fastq.gz
└── samplesheet.csv
```

## Prepare the samplesheet

The samplesheet tells MTD Explorer which samples should be analyzed, which
experimental group each sample belongs to, and which group comparisons should
be performed.

Example samplesheets are available in the repository:

[Open the examples directory](https://github.com/patrick-douglas/MTD-Explorer/tree/main/examples)

After cloning the repository, these files are also available locally at:

```text
~/MTD/examples/
```

### Hypothetical samplesheet preview

A comparison samplesheet may look like this:

| sample_name | group | comparisons | group1 | vs | group2 |
|---|---|---|---|---|---|
| sample1 | Treatment1 |  | Treatment1 | vs | Control |
| sample2 | Treatment1 |  | Treatment2 | vs | Control |
| sample3 | Treatment1 |  |  |  |
| sample4 | Treatment2 |  |  |  |
| sample5 | Treatment2 |  |  |  |
| sample6 | Treatment2 |  |  |  |
| sample7 | Treatment2 |  |  |  |
| sample8 | Control |  |  |  |  |
| sample9 | Control |  |  |  |  |
| sample10 | Control |  |  |  |  |

The same file in CSV format would be:

```csv
sample_name,group,comparisons,group1,vs,group2
sample1,Treatment1,,Treatment1,vs,Control
sample2,Treatment1,,Treatment2,vs,Control
sample3,Treatment1,,,,
sample4,Treatment2,,,,
sample5,Treatment2,,,,
sample6,Treatment2,,,,
sample7,Treatment2,,,,
sample8,Control,,,,
sample9,Control,,,,
sample10,Control,,,,
```

### Samplesheet columns

| Column | Description |
|---|---|
| `sample_name` | Sample identifier used by MTD Explorer to find the corresponding FASTQ file or FASTQ pair |
| `group` | Experimental group assigned to each sample |
| `comparisons` | Header used to mark the comparison section of the samplesheet |
| `group1` | First group in the comparison |
| `vs` | Comparison separator |
| `group2` | Second group in the comparison |

In the example above, MTD Explorer will analyze the samples assigned to
`Treatment1`, `Treatment2`, and `Control`, and perform the following
comparisons:

```text
Treatment1 vs Control
Treatment2 vs Control
```

The group names in the comparison section must match the group names used in
the `group` column.

!!! important "Sample names and FASTQ files"

    The values in the `sample_name` column must match the FASTQ filenames.

    For example, the sample name:

    ```text
    sample1
    ```

    can match a single-end file such as:

    ```text
    sample1.fastq.gz
    ```

    or a paired-end pair such as:

    ```text
    sample1_R1.fastq.gz
    sample1_R2.fastq.gz
    ```

!!! note "Use the repository examples as templates"

    The `examples/` directory contains samplesheet files that can be copied
    and adapted for new analyses.

    Before running MTD Explorer, compare your samplesheet with the examples
    provided in the repository.

## FASTQ naming

MTD Explorer detects the sequencing layout from the FASTQ files.

### Supported single-end examples

```text
sample1.fastq.gz
sample1.fq.gz
Trimmed_sample1.fastq.gz
Trimmed_sample1.fq.gz
```

### Supported paired-end examples

```text
sample1_R1.fastq.gz
sample1_R2.fastq.gz
```

```text
sample1_1.fastq.gz
sample1_2.fastq.gz
```

```text
sample1_S1_L001_R1_001.fastq.gz
sample1_S1_L001_R2_001.fastq.gz
```

!!! warning "Do not mix single-end and paired-end samples"

    All samples in one MTD Explorer run must use the same read layout.

    Do not mix single-end and paired-end samples in the same samplesheet.

## Single-end example

From the project directory:

```bash
cd ~/MTD_projects/test_run
```

Run MTD Explorer using automatic read-layout detection:

```bash
bash ~/MTD/MTD_explorer.sh \
  --input samplesheet.csv \
  --output MTD_results_SE \
  --hostid 9606 \
  --read-layout se \
  --threads 16
```

In this example:

- `samplesheet.csv` describes the samples and groups;
- `MTD_results_SE` is the output directory;
- `9606` is the NCBI Taxonomy ID for the host species;
- `--read-layout se` forces single-end mode;
- `--threads 16` uses 16 CPU threads.

!!! tip

    If the FASTQ names are standard and all samples use the same layout,
    you can use:

    ```bash
    --read-layout auto
    ```

## Paired-end example

For paired-end data, use matching R1 and R2 FASTQ files for each sample.

Example:

```bash
bash ~/MTD/MTD_explorer.sh \
  --input samplesheet.csv \
  --output MTD_results_PE \
  --hostid 9606 \
  --read-layout pe \
  --threads 16
```

MTD Explorer validates paired-end inputs and checks that R1 and R2 are
resolved for each sample.

## Exploratory analysis

Use exploratory mode when you do not want to run group-comparison or
DEG-dependent steps.

```bash
bash ~/MTD/MTD_explorer.sh \
  --input samplesheet.csv \
  --output MTD_results_exploratory \
  --hostid 9606 \
  --analysis-mode exploratory \
  --read-layout auto \
  --threads 16
```

The aliases below are also accepted:

```bash
--exploratory
```

```bash
--no-comparison
```

Exploratory mode is useful for:

- initial data inspection;
- single-group datasets;
- datasets without a valid comparison design;
- checking host expression and microbiome profiles before formal comparison.

## Comparison analysis

Use comparison mode when the samplesheet contains experimental groups and
you want MTD Explorer to run group-dependent analysis steps.

```bash
bash ~/MTD/MTD_explorer.sh \
  --input samplesheet.csv \
  --output MTD_results_comparison \
  --hostid 9606 \
  --analysis-mode comparison \
  --read-layout auto \
  --threads 16
```

Comparison mode requires a valid group structure in the samplesheet.

## Automatic analysis mode

The default analysis mode is:

```text
auto
```

In automatic mode, MTD Explorer decides whether to run comparison or
exploratory outputs based on the samplesheet structure.

```bash
bash ~/MTD/MTD_explorer.sh \
  --input samplesheet.csv \
  --output MTD_results_auto \
  --hostid 9606 \
  --analysis-mode auto \
  --read-layout auto \
  --threads 16
```

## Trimming and alignment options

By default, MTD Explorer uses `fastp` for read preprocessing.

The default trimming profile is conservative. It is designed to remove very
low-quality terminal bases and problematic sequence features while avoiding
overly aggressive read loss.

!!! note "Trimming strategy"

    The default MTD Explorer trimming strategy is based on the rationale
    described by MacManes (2014), who evaluated quality trimming in
    high-throughput mRNA-seq data and showed that overly aggressive trimming
    can negatively affect transcriptome results.

    This is particularly relevant for RNA-seq and metatranscriptomic studies
    involving non-model organisms, where preserving informative reads may be
    important for host alignment, microbial classification, and functional
    profiling.

    Reference:

    MacManes MD. 2014. *On the optimal trimming of high-throughput mRNA
    sequence data*. Frontiers in Genetics 5:13.
    DOI: [10.3389/fgene.2014.00013](https://doi.org/10.3389/fgene.2014.00013)

To skip trimming:

```bash
bash ~/MTD/MTD_explorer.sh \
  --input samplesheet.csv \
  --output MTD_results_no_trim \
  --hostid 9606 \
  --no-trim \
  --read-layout auto \
  --threads 16
```

To use Magic-BLAST instead of HISAT2 for host alignment:

```bash
bash ~/MTD/MTD_explorer.sh \
  --input samplesheet.csv \
  --output MTD_results_blast \
  --hostid 9606 \
  --blast \
  --read-layout auto \
  --threads 16
```

## Custom Kraken2 databases
MTD Explorer can use custom Kraken2 databases for host filtering and
microbiome classification.

### Custom host-filtering database

```bash
bash ~/MTD/MTD_explorer.sh \
  --input samplesheet.csv \
  --output MTD_results_custom_host_kraken \
  --hostid 9606 \
  --kraken-host-db /path/to/kraken2_host_db \
  --read-layout auto \
  --threads 16
```

!!! note

    `--kraken-host-db` changes only the Kraken2 host-filtering database.

    It does not replace the host genome, GTF annotation, BLAST/HISAT2
    resources, featureCounts resources, or host downstream annotation
    associated with `--hostid`.

### Custom microbiome database

```bash
bash ~/MTD/MTD_explorer.sh \
  --input samplesheet.csv \
  --output MTD_results_custom_micro_db \
  --hostid 9606 \
  --kraken-micro-db /path/to/kraken2_microbiome_db \
  --read-layout auto \
  --threads 16
```

This can be used for targeted microbiome databases, such as a curated
taxonomic group database.

## ssGSEA GMT selection

The `--ssgsea-gmt` option controls the GMT file used for ssGSEA.

Available modes:

| Value | Meaning |
|---|---|
| `default` | Use the legacy MSigDB C2 symbols GMT |
| `auto` | Use the master eggNOG/GO GMT generated by `Create_custom_host.sh` for the selected `--hostid` |
| `FILE` | Use a user-provided GMT file |

Example using automatic host-specific GMT selection:

```bash
bash ~/MTD/MTD_explorer.sh \
  --input samplesheet.csv \
  --output MTD_results_ssgsea_auto \
  --hostid 9606 \
  --ssgsea-gmt auto \
  --read-layout auto \
  --threads 16
```

## Detected microbiome read extraction

Detected microbiome read extraction is disabled by default because it can
generate many files.

To enable extraction for the top detected taxa:

```bash
bash ~/MTD/MTD_explorer.sh \
  --input samplesheet.csv \
  --output MTD_results_extract_reads \
  --hostid 9606 \
  --extract-microbiome-reads \
  --extract-microbiome-reads-top-n 50 \
  --read-layout auto \
  --threads 16
```

Use:

```bash
--extract-microbiome-reads-top-n 0
```

to extract reads for all detected taxa.

## Output directory behavior

If the selected output directory already exists, MTD Explorer asks whether
it should delete and overwrite the directory.

!!! danger "Existing output directory"

    Confirming overwrite permanently removes the existing output directory.

    Use a new output directory name when you want to preserve a previous run.

## After the run starts

MTD Explorer prints progress messages and writes intermediate files inside
the selected output directory.

A successful run creates output sections for host analysis, Kraken2
classification, Bracken abundance estimation, HUMAnN outputs, functional
analysis, exploratory outputs, and integration analyses.

The exact output structure is described in the User Guide.

## Troubleshooting

If MTD Explorer stops with an error:

1. copy the complete command used;
2. save the terminal output;
3. record the current Git commit;
4. check whether the FASTQ filenames match the samplesheet;
5. check whether the selected `--hostid` is installed and supported;
6. check whether the installation checker passes.

Record the current commit with:

```bash
git log -1 --oneline
```

Run the installation checker again with:

```bash
bash ~/MTD/MTD_check_installation.sh \
  --mode full \
  --report-dir ./MTD_check_after_error
```

If the problem persists, open a GitHub issue:

[Open a GitHub issue](https://github.com/patrick-douglas/MTD-Explorer/issues)

## Next step

After completing this quick start, continue to the User Guide for detailed
information about input files, analysis modes, custom hosts, and output
interpretation.
