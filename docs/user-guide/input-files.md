# Input files

This page describes the input files expected by [MTD Explorer][mtd-explorer],
including the samplesheet, FASTQ naming rules, sequencing layout detection,
and optional public accession handling.

For a minimal working example, see the
[Quick start guide](../getting-started/quick-start.md).

## Required input

A standard MTD Explorer run requires:

1. a samplesheet CSV file;
2. FASTQ files for each sample, or supported public run accessions;
3. a host Taxonomy ID provided with `--hostid`.

The samplesheet is passed with:

```bash
--input samplesheet.csv
```

The selected host is passed with:

```bash
--hostid TAXID
```

## Recommended project layout

Create one project directory per analysis, outside the MTD Explorer source
directory.

Example for single-end data:

```text
project/
├── fastq/
│   ├── sample1.fastq.gz
│   ├── sample2.fastq.gz
│   └── sample3.fastq.gz
└── samplesheet.csv
```

Example for paired-end data:

```text
project/
├── fastq/
│   ├── sample1_R1.fastq.gz
│   ├── sample1_R2.fastq.gz
│   ├── sample2_R1.fastq.gz
│   ├── sample2_R2.fastq.gz
│   ├── sample3_R1.fastq.gz
│   └── sample3_R2.fastq.gz
└── samplesheet.csv
```

Run MTD Explorer from the project directory or provide absolute paths.

Example:

```bash
cd ~/MTD_projects/my_analysis

bash ~/MTD/MTD_explorer.sh \
  --input samplesheet.csv \
  --output MTD_results \
  --hostid 9606 \
  --read-layout auto \
  --threads 16
```

!!! tip

    Absolute paths are safer when running analyses from another directory.

## Samplesheet format

The samplesheet tells MTD Explorer which samples should be analyzed and which
experimental group each sample belongs to.

A comparison samplesheet can also include a comparison section describing the
group contrasts to be tested.

### Example samplesheet

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

The same structure as a table:

| sample_name | group | comparisons | group1 | vs | group2 |
|---|---|---|---|---|---|
| sample1 | Treatment1 |  | Treatment1 | vs | Control |
| sample2 | Treatment1 |  | Treatment2 | vs | Control |
| sample3 | Treatment1 |  |  |  |  |
| sample4 | Treatment2 |  |  |  |  |
| sample5 | Treatment2 |  |  |  |  |
| sample6 | Treatment2 |  |  |  |  |
| sample7 | Treatment2 |  |  |  |  |
| sample8 | Control |  |  |  |  |
| sample9 | Control |  |  |  |  |
| sample10 | Control |  |  |  |  |

## Samplesheet columns

| Column | Required | Description |
|---|---|---|
| `sample_name` | yes | Sample identifier used to find the corresponding FASTQ file or public run accession |
| `group` | yes | Experimental group assigned to each sample |
| `comparisons` | recommended for comparison analyses | Header used to mark the comparison section |
| `group1` | recommended for comparison analyses | First group in the comparison |
| `vs` | recommended for comparison analyses | Comparison separator |
| `group2` | recommended for comparison analyses | Second group in the comparison |

The first column must contain the sample names or public run accessions.

The second column must contain the experimental group for each sample.

!!! important "Group names must match"

    Group names used in the comparison section must match the group names used
    in the `group` column.

    For example, if the group column contains `Treatment1`, the comparison
    section should also use `Treatment1`, not `treatment1`, `Treat1`, or
    `Treatment 1`.

## Sample names and FASTQ files

The values in the `sample_name` column must match the FASTQ filenames.

For example, this sample name:

```text
sample1
```

can match this single-end file:

```text
sample1.fastq.gz
```

or this paired-end pair:

```text
sample1_R1.fastq.gz
sample1_R2.fastq.gz
```

MTD Explorer resolves FASTQ files by sample name before preparing normalized
input files for downstream analysis.

## Supported single-end FASTQ names

For a sample named `sample1`, the following single-end names are supported:

```text
sample1.fastq
sample1.fq
sample1.fastq.gz
sample1.fq.gz
Trimmed_sample1.fastq
Trimmed_sample1.fq
Trimmed_sample1.fastq.gz
Trimmed_sample1.fq.gz
```

## Supported paired-end FASTQ names

For a sample named `sample1`, the following paired-end patterns are supported.

### `_R1` and `_R2`

```text
sample1_R1.fastq.gz
sample1_R2.fastq.gz
```

```text
sample1_R1.fq.gz
sample1_R2.fq.gz
```

### `_1` and `_2`

```text
sample1_1.fastq.gz
sample1_2.fastq.gz
```

```text
sample1_1.fq.gz
sample1_2.fq.gz
```

### Illumina-style names

```text
sample1_S1_L001_R1_001.fastq.gz
sample1_S1_L001_R2_001.fastq.gz
```

MTD Explorer also accepts related paired-end names containing `_R1_` and
`_R2_` after the sample name.

## Sequencing layout

MTD Explorer supports:

| Layout | Meaning |
|---|---|
| `se` | single-end |
| `pe` | paired-end |
| `auto` | detect layout from FASTQ filenames |

Use automatic detection:

```bash
--read-layout auto
```

Force single-end mode:

```bash
--read-layout se
```

Force paired-end mode:

```bash
--read-layout pe
```

!!! warning "Do not mix sequencing layouts"

    All samples in one MTD Explorer run must use the same sequencing layout.

    Do not mix single-end and paired-end samples in the same samplesheet.

## Public run accessions

MTD Explorer can also work with public run accessions when local FASTQ files
are not already present.

Supported accession prefixes are:

| Prefix | Archive |
|---|---|
| `SRR` | NCBI Sequence Read Archive |
| `ERR` | European Nucleotide Archive |
| `DRR` | DDBJ Sequence Read Archive |

Example samplesheet using public run accessions:

```csv
sample_name,group,comparisons,group1,vs,group2
SRR000001,Control,,Treatment,vs,Control
SRR000002,Control,,,,
SRR000003,Treatment,,,,
SRR000004,Treatment,,,,
```

When an accession is detected and no matching local FASTQ file exists, MTD
Explorer uses [SRA Toolkit][sra-toolkit] utilities to download and convert the
run into FASTQ files.

The generated FASTQ files are then validated with the same input-detection
logic used for local FASTQ files.

!!! note

    Public accession downloads depend on internet access, archive availability,
    and local [SRA Toolkit][sra-toolkit] configuration.

    For reproducible analyses, it is often useful to download FASTQ files once
    and keep them with the project.

## FASTQ preparation

After input detection, MTD Explorer creates normalized prepared FASTQ files in
the temporary run directory.

The normalized prepared filenames follow these patterns:

| Layout | Prepared files |
|---|---|
| single-end | `Trimmed_SAMPLE.fq.gz` |
| paired-end | `Trimmed_SAMPLE_R1.fq.gz` and `Trimmed_SAMPLE_R2.fq.gz` |

By default, MTD Explorer uses [fastp][fastp] for read preprocessing.

When `--no-trim` is used, MTD Explorer skips [fastp][fastp]. Compressed FASTQ
files are copied, while uncompressed FASTQ files are compressed before
downstream processing.

```bash
bash ~/MTD/MTD_explorer.sh \
  --input samplesheet.csv \
  --output MTD_results_no_trim \
  --hostid 9606 \
  --no-trim
```

## Paired-end validation

For paired-end data, MTD Explorer checks that each sample has exactly one R1
file and exactly one R2 file.

It also validates that paired FASTQ files contain matching record counts when
paired FASTQ validation is enabled.

A run will stop if:

- only R1 or only R2 is found;
- more than one R1 or R2 candidate is found for a sample;
- single-end and paired-end files are detected for the same sample;
- paired-end record counts do not match.

## Common input problems

### Sample name does not match FASTQ filename

If the samplesheet contains:

```text
sample1
```

but the file is named:

```text
Sample_1.fastq.gz
```

MTD Explorer will not treat them as the same sample.

Fix either the samplesheet or the FASTQ filename so the sample identifiers
match.

### Mixed single-end and paired-end inputs

A single run cannot contain both:

```text
sample1.fastq.gz
```

and:

```text
sample2_R1.fastq.gz
sample2_R2.fastq.gz
```

Split these into separate analyses.

### Ambiguous paired-end files

Avoid keeping duplicate files for the same sample in the input directory.

For example, do not keep all of these together for one sample:

```text
sample1_R1.fastq.gz
sample1_R1_001.fastq.gz
sample1_R2.fastq.gz
sample1_R2_001.fastq.gz
```

MTD Explorer expects one resolved R1 and one resolved R2 per sample.

### Group names are inconsistent

These are treated as different groups:

```text
Control
control
CTRL
Control_1
```

Use consistent group names throughout the samplesheet.

## Example input check before running

Before launching a full run, inspect your project directory:

```bash
ls -lh
ls -lh fastq 2>/dev/null || true
head samplesheet.csv
```

Check unique sample names and groups:

```bash
awk -F',' 'NR > 1 {print $1}' samplesheet.csv | sort -u

awk -F',' 'NR > 1 {print $2}' samplesheet.csv | sort -u
```

Then run MTD Explorer with automatic read-layout detection:

```bash
bash ~/MTD/MTD_explorer.sh \
  --input samplesheet.csv \
  --output MTD_results \
  --hostid 9606 \
  --read-layout auto \
  --threads 16
```

## Example files

Example samplesheets are available in the repository:

[Open the examples directory][examples]

After cloning the repository, they are also available locally at:

```text
~/MTD/examples/
```

Use these files as templates for new analyses.

## Next step

After preparing input files, continue to
[Analysis modes](analysis-modes.md) or review the
[Command-line reference](command-line.md).

[mtd-explorer]: https://github.com/patrick-douglas/MTD-Explorer-Explorer
[examples]: https://github.com/patrick-douglas/MTD-Explorer-Explorer/tree/main/examples
[fastp]: https://github.com/OpenGene/fastp
[sra-toolkit]: https://github.com/ncbi/sra-tools
