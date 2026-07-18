# Command-line reference

This page documents the main command-line options available in
`MTD_explorer.sh`.

For a practical first run, see the [Quick start guide](../getting-started/quick-start.md).

## Main command

The main MTD Explorer entry point is:

```bash
bash MTD_explorer.sh [options]
```

Display the current help message with:

```bash
bash MTD_explorer.sh --help
```

## Required arguments

Every MTD Explorer run requires an input samplesheet, an output directory,
and a host NCBI Taxonomy ID.

| Option | Argument | Description |
|---|---|---|
| `-i`, `--input` | `FILE` | Path to `samplesheet.csv` |
| `-o`, `--output` | `DIR` | Output directory |
| `-h`, `--hostid` | `TAXID` | Host species Taxonomy ID used for annotation and downstream host analysis |

Minimal command:

```bash
bash MTD_explorer.sh \
  --input samplesheet.csv \
  --output MTD_results \
  --hostid 9606
```

## Analysis mode options

MTD Explorer can run in automatic, comparison, or exploratory mode.

| Option | Argument | Description | Default |
|---|---|---|---|
| `--analysis-mode` | `auto`, `comparison`, `exploratory` | Selects how group-dependent analysis steps are handled | `auto` |
| `--exploratory` | — | Alias for `--analysis-mode exploratory` | — |
| `--no-comparison` | — | Alias for `--analysis-mode exploratory` | — |

### `auto`

In automatic mode, MTD Explorer inspects the samplesheet and decides whether
comparison analysis is possible.

```bash
bash MTD_explorer.sh \
  --input samplesheet.csv \
  --output MTD_results_auto \
  --hostid 9606 \
  --analysis-mode auto
```

### `comparison`

Comparison mode requires experimental groups and runs group-dependent
analysis steps.

```bash
bash MTD_explorer.sh \
  --input samplesheet.csv \
  --output MTD_results_comparison \
  --hostid 9606 \
  --analysis-mode comparison
```

### `exploratory`

Exploratory mode skips DEG-dependent comparison steps and generates
non-comparison outputs.

```bash
bash MTD_explorer.sh \
  --input samplesheet.csv \
  --output MTD_results_exploratory \
  --hostid 9606 \
  --analysis-mode exploratory
```

Equivalent aliases:

```bash
bash MTD_explorer.sh \
  --input samplesheet.csv \
  --output MTD_results_exploratory \
  --hostid 9606 \
  --exploratory
```

```bash
bash MTD_explorer.sh \
  --input samplesheet.csv \
  --output MTD_results_exploratory \
  --hostid 9606 \
  --no-comparison
```

## Input and metadata options

| Option | Argument | Description | Default |
|---|---|---|---|
| `-m`, `--metadata` | `FILE` | Optional metadata CSV file | none |
| `-p`, `--pdm` | `METHOD` | HAllA association metric | `spearman` |
| `--threads` | `INT` | Number of CPU threads | detected with `nproc` |

Supported HAllA metrics:

```text
spearman
pearson
mi
nmi
xicor
dcor
```

Example using metadata and 16 threads:

```bash
bash MTD_explorer.sh \
  --input samplesheet.csv \
  --metadata metadata.csv \
  --output MTD_results \
  --hostid 9606 \
  --threads 16
```

## Read preprocessing options

| Option | Argument | Description | Default |
|---|---|---|---|
| `-l`, `--trim-length` | `INT` | Minimum read length required by `fastp` | `35` |
| `-t`, `--no-trim` | — | Skip `fastp` trimming | trimming enabled |

By default, MTD Explorer uses `fastp` for read preprocessing.

The default trimming profile is conservative and is designed to remove very
low-quality terminal bases while avoiding overly aggressive read loss.

To skip trimming:

```bash
bash MTD_explorer.sh \
  --input samplesheet.csv \
  --output MTD_results_no_trim \
  --hostid 9606 \
  --no-trim
```

To change the minimum read length required after trimming:

```bash
bash MTD_explorer.sh \
  --input samplesheet.csv \
  --output MTD_results_len50 \
  --hostid 9606 \
  --trim-length 50
```

## Read layout options

| Option | Argument | Description | Default |
|---|---|---|---|
| `--read-layout` | `auto`, `se`, `pe` | Sequencing library layout | `auto` |

Available values:

| Value | Meaning |
|---|---|
| `auto` | Detect layout from FASTQ filenames |
| `se` | Force single-end input |
| `pe` | Force paired-end input |

All samples in one run must use the same layout.

Single-end example:

```bash
bash MTD_explorer.sh \
  --input samplesheet.csv \
  --output MTD_results_SE \
  --hostid 9606 \
  --read-layout se
```

Paired-end example:

```bash
bash MTD_explorer.sh \
  --input samplesheet.csv \
  --output MTD_results_PE \
  --hostid 9606 \
  --read-layout pe
```

Automatic layout detection:

```bash
bash MTD_explorer.sh \
  --input samplesheet.csv \
  --output MTD_results_auto_layout \
  --hostid 9606 \
  --read-layout auto
```

## Host processing options

| Option | Argument | Description | Default |
|---|---|---|---|
| `-b`, `--blast` | — | Use Magic-BLAST instead of HISAT2 for host alignment | HISAT2 |
| `--read-layout` | `auto`, `se`, `pe` | Select or detect sequencing layout | `auto` |

Use Magic-BLAST:

```bash
bash MTD_explorer.sh \
  --input samplesheet.csv \
  --output MTD_results_blast \
  --hostid 9606 \
  --blast
```

Use default HISAT2 host alignment:

```bash
bash MTD_explorer.sh \
  --input samplesheet.csv \
  --output MTD_results_hisat2 \
  --hostid 9606
```

## Kraken2 host filtering options

These options control the Kraken2 database and classification parameters
used during host read filtering.

| Option | Argument | Description | Default |
|---|---|---|---|
| `--kraken-host-db` | `DIR` | Optional custom Kraken2 host-filtering database | selected automatically from `--hostid` |
| `--kraken-host-confidence` | `FLOAT` | Kraken2 `--confidence` for host filtering | `0.05` |
| `--kraken-host-min-hit-groups` | `INT` | Kraken2 `--minimum-hit-groups` for host filtering | `3` |

Example:

```bash
bash MTD_explorer.sh \
  --input samplesheet.csv \
  --output MTD_results_custom_host_db \
  --hostid 59463 \
  --kraken-host-db /path/to/kraken2_host_db \
  --kraken-host-confidence 0.05 \
  --kraken-host-min-hit-groups 3
```

!!! important

    `--kraken-host-db` changes only the Kraken2 database used for host
    read filtering.

    It does not replace the host genome, GTF annotation, HISAT2 index,
    Magic-BLAST database, featureCounts resources, or host downstream
    annotation associated with `--hostid`.

## Kraken2 microbiome classification options

These options control the Kraken2 database and classification parameters
used for non-host or microbiome classification.

| Option | Argument | Description | Default |
|---|---|---|---|
| `--kraken-micro-db` | `DIR` | Optional custom Kraken2 microbiome or target database | `$MTDIR/kraken2DB_micro` |
| `--kraken-micro-confidence` | `FLOAT` | Kraken2 `--confidence` for microbiome classification | `0.10` |
| `--kraken-micro-min-hit-groups` | `INT` | Kraken2 `--minimum-hit-groups` for microbiome classification | `3` |

Example:

```bash
bash MTD_explorer.sh \
  --input samplesheet.csv \
  --output MTD_results_custom_micro_db \
  --hostid 59463 \
  --kraken-micro-db /path/to/kraken2_microbiome_db \
  --kraken-micro-confidence 0.10 \
  --kraken-micro-min-hit-groups 3
```

Custom microbiome databases can be useful for targeted analyses, such as
focused viral, bacterial, fungal, or parasite databases.

## Bracken options

| Option | Argument | Description | Default |
|---|---|---|---|
| `-r`, `--bracken-read-len` | `INT` | Bracken read length | `75` |
| `--bracken-threshold` | `INT` | Bracken `-t` minimum read threshold | `10` |

Example with a 100 bp Bracken read length:

```bash
bash MTD_explorer.sh \
  --input samplesheet.csv \
  --output MTD_results_bracken100 \
  --hostid 9606 \
  --bracken-read-len 100
```

!!! warning

    The selected Bracken read length must match a Bracken distribution file
    already built for the selected microbiome Kraken2 database.

    For example, `--bracken-read-len 75` requires a corresponding
    `database75mers.kmer_distrib` file inside the Kraken2 microbiome database.

## ssGSEA options

| Option | Argument | Description | Default |
|---|---|---|---|
| `--ssgsea-gmt` | `default`, `auto`, or `FILE` | Selects the GMT gene-set file used for ssGSEA | `default` |

Available values:

| Value | Meaning |
|---|---|
| `default` | Use the legacy MSigDB C2 symbols GMT |
| `auto` | Use the persistent master eggNOG/GO GMT generated by `Create_custom_host.sh` for the selected `--hostid` |
| `FILE` | Use an existing GMT file directly |

Example using default ssGSEA GMT:

```bash
bash MTD_explorer.sh \
  --input samplesheet.csv \
  --output MTD_results_ssgsea_default \
  --hostid 9606 \
  --ssgsea-gmt default
```

Example using host-specific eggNOG/GO GMT selection:

```bash
bash MTD_explorer.sh \
  --input samplesheet.csv \
  --output MTD_results_ssgsea_auto \
  --hostid 6526 \
  --ssgsea-gmt auto
```

Example using a custom GMT file:

```bash
bash MTD_explorer.sh \
  --input samplesheet.csv \
  --output MTD_results_ssgsea_custom \
  --hostid 9606 \
  --ssgsea-gmt /path/to/custom_gene_sets.gmt
```

## Detected microbiome read extraction options

Detected microbiome read extraction is disabled by default because it can
generate many files.

| Option | Argument | Description | Default |
|---|---|---|---|
| `--extract-microbiome-reads` | — | Extract Kraken-classified reads for detected microbiome taxa | disabled |
| `--extract-microbiome-reads-top-n` | `INT` | Number of top taxa to extract from the absolute detected microbiome ranking; use `0` for all detected taxa | `50` |
| `--extract-microbiome-reads-min-abundance` | `NUMERIC` | Minimum abundance required in a sample to extract reads | `0` |

Example:

```bash
bash MTD_explorer.sh \
  --input samplesheet.csv \
  --output MTD_results_extract_reads \
  --hostid 9606 \
  --extract-microbiome-reads \
  --extract-microbiome-reads-top-n 50 \
  --extract-microbiome-reads-min-abundance 0
```

Extract reads for all detected taxa:

```bash
bash MTD_explorer.sh \
  --input samplesheet.csv \
  --output MTD_results_extract_all_reads \
  --hostid 9606 \
  --extract-microbiome-reads \
  --extract-microbiome-reads-top-n 0
```

## Other options

| Option | Description |
|---|---|
| `--help` | Display the command-line help message and exit |

## Complete examples

### Default host database from `--hostid`

```bash
bash MTD_explorer.sh \
  --input samplesheet.csv \
  --output MTD_results_Myotis_auto \
  --hostid 59463 \
  --blast \
  --no-trim
```

### Custom Kraken2 host and microbiome databases

This example uses a custom Kraken2 host-filtering database while keeping
`--hostid` for host annotation and downstream host analysis.

```bash
bash MTD_explorer.sh \
  --input samplesheet.csv \
  --output MTD_results_Carollia_Myotis \
  --hostid 59463 \
  --blast \
  --no-trim \
  --kraken-host-db /home/me/MTD/kraken2DB_Carollia_Myotis/ \
  --kraken-micro-db /home/me/MTD/Kraken2DB_trematoda/ \
  --kraken-host-confidence 0.05 \
  --kraken-host-min-hit-groups 3 \
  --kraken-micro-confidence 0.10 \
  --kraken-micro-min-hit-groups 3 \
  --bracken-threshold 10
```

### Reference-matched eggNOG/GO GMT

This example uses the master eggNOG/GO GMT created by
`Create_custom_host.sh` for the selected host Taxonomy ID.

```bash
bash MTD_explorer.sh \
  --input samplesheet.csv \
  --output MTD_results_custom_ssGSEA \
  --hostid 6526 \
  --blast \
  --no-trim \
  --ssgsea-gmt auto
```

## Removed development options

The development-only `--custom-raw-path` option was removed before stable
documentation because it was used only for internal testing and is not part
of the public MTD Explorer command-line interface.

## Exit behavior

MTD Explorer exits with an error when required arguments are missing,
unsupported values are provided, input FASTQ files cannot be resolved, read
layouts are mixed, required databases are incomplete, or required downstream
files cannot be created.

For installation-level problems, run:

```bash
bash MTD_check_installation.sh --mode full
```

For usage examples, continue with the [Quick start guide](../getting-started/quick-start.md).
