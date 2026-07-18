# Analysis modes

MTD Explorer can run in three analysis modes:

```text
auto
comparison
exploratory
```

The analysis mode controls whether MTD Explorer should run group-comparison
steps or produce only non-comparison exploratory outputs.

For details about input files and samplesheet formatting, see
[Input files](input-files.md).

## Overview

| Mode | Best used when | Main behavior |
|---|---|---|
| `auto` | You want MTD Explorer to decide from the samplesheet | Uses comparison mode when two or more groups are detected; otherwise uses exploratory mode |
| `comparison` | You have experimental groups and want formal group-dependent analyses | Requires at least two groups and runs comparison-dependent steps |
| `exploratory` | You have one group, no valid comparison, or want descriptive outputs only | Skips DEG-dependent comparison steps and produces non-comparison outputs |

The default mode is:

```bash
--analysis-mode auto
```

## What the analysis mode changes

The analysis mode affects group-dependent downstream analysis.

It does not change:

- input FASTQ detection;
- read layout detection;
- [fastp][fastp] preprocessing;
- [Kraken2][kraken2] host filtering;
- [Kraken2][kraken2] microbiome classification;
- [Bracken][bracken] abundance estimation;
- [HUMAnN][humann] functional profiling;
- the host selected with `--hostid`;
- custom Kraken2 database selection;
- Bracken read length settings.

The main difference is whether MTD Explorer should run analyses that depend on
experimental group comparisons, such as differential abundance or differential
expression steps using tools such as [DESeq2][deseq2].

## Automatic mode

Automatic mode is the default.

```bash
--analysis-mode auto
```

In this mode, MTD Explorer inspects the `group` column in the samplesheet.

If two or more groups are detected, MTD Explorer resolves the run as
comparison mode.

If fewer than two groups are detected, MTD Explorer resolves the run as
exploratory mode.

Example:

```bash
bash ~/MTD/MTD_explorer.sh \
  --input samplesheet.csv \
  --output MTD_results_auto \
  --hostid 9606 \
  --analysis-mode auto \
  --read-layout auto \
  --threads 16
```

### When to use automatic mode

Use automatic mode when:

- you want the safest default behavior;
- you are testing a new dataset;
- you are not sure whether the samplesheet contains a valid comparison design;
- you want MTD Explorer to avoid forcing comparison analyses for single-group data.

!!! tip

    For most users, `--analysis-mode auto` is the recommended starting point.

## Comparison mode

Comparison mode is used when the samplesheet contains experimental groups and
you want MTD Explorer to run group-dependent analysis steps.

```bash
--analysis-mode comparison
```

Example:

```bash
bash ~/MTD/MTD_explorer.sh \
  --input samplesheet.csv \
  --output MTD_results_comparison \
  --hostid 9606 \
  --analysis-mode comparison \
  --read-layout auto \
  --threads 16
```

Comparison mode requires at least two groups in the samplesheet.

If comparison mode is requested but fewer than two groups are detected, MTD
Explorer stops with an error.

### When to use comparison mode

Use comparison mode when:

- the samplesheet has at least two biological or experimental groups;
- the group labels are correct and consistent;
- you want comparison-dependent host or microbiome outputs;
- you want downstream association analysis with [HAllA][halla] where appropriate;
- the study design supports biological interpretation of group differences.

### Example comparison samplesheet

```csv
sample_name,group,comparisons,group1,vs,group2
sample1,Treatment,,Treatment,vs,Control
sample2,Treatment,,,,
sample3,Treatment,,,,
sample4,Control,,,,
sample5,Control,,,,
sample6,Control,,,,
```

In this example, MTD Explorer detects two groups:

```text
Treatment
Control
```

and the requested comparison is:

```text
Treatment vs Control
```

## Exploratory mode

Exploratory mode is used when you do not want formal group-comparison analysis.

```bash
--analysis-mode exploratory
```

The following aliases are also accepted:

```bash
--exploratory
```

```bash
--no-comparison
```

Example:

```bash
bash ~/MTD/MTD_explorer.sh \
  --input samplesheet.csv \
  --output MTD_results_exploratory \
  --hostid 9606 \
  --analysis-mode exploratory \
  --read-layout auto \
  --threads 16
```

Equivalent command using the alias:

```bash
bash ~/MTD/MTD_explorer.sh \
  --input samplesheet.csv \
  --output MTD_results_exploratory \
  --hostid 9606 \
  --no-comparison \
  --read-layout auto \
  --threads 16
```

### When to use exploratory mode

Use exploratory mode when:

- the dataset contains only one group;
- the samplesheet does not contain a valid comparison design;
- you are inspecting a dataset before defining contrasts;
- you want descriptive host, microbiome, functional, or quality-control outputs;
- you want to avoid DEG-dependent steps.

Exploratory mode is also useful for sparse microbiome or targeted database
analyses where some samples may have few or no detected taxa.

### Example exploratory samplesheet

```csv
sample_name,group,comparisons,group1,vs,group2
sample1,All,,,,
sample2,All,,,,
sample3,All,,,,
sample4,All,,,,
```

This samplesheet has one group:

```text
All
```

In automatic mode, MTD Explorer would resolve this run as exploratory.

## Choosing the right mode

| Situation | Recommended mode |
|---|---|
| First test of a dataset | `auto` |
| One-group dataset | `exploratory` |
| Dataset with no meaningful biological comparison | `exploratory` |
| Dataset with two or more valid groups | `comparison` or `auto` |
| You want MTD Explorer to decide safely | `auto` |
| You want to force comparison analysis | `comparison` |
| You want to skip DEG-dependent steps | `exploratory` |

## What happens with single-group data

Single-group data should not be forced into comparison mode.

For example:

```csv
sample_name,group,comparisons,group1,vs,group2
sample1,Control,,,,
sample2,Control,,,,
sample3,Control,,,,
```

This dataset has only one group.

Use:

```bash
--analysis-mode exploratory
```

or leave the default:

```bash
--analysis-mode auto
```

In automatic mode, MTD Explorer detects that fewer than two groups are present
and switches to exploratory behavior.

## What happens with two-group data

A two-group samplesheet can be used in comparison mode.

Example:

```csv
sample_name,group,comparisons,group1,vs,group2
sample1,Treated,,Treated,vs,Control
sample2,Treated,,,,
sample3,Treated,,,,
sample4,Control,,,,
sample5,Control,,,,
sample6,Control,,,,
```

Use:

```bash
--analysis-mode comparison
```

or:

```bash
--analysis-mode auto
```

In automatic mode, MTD Explorer detects two groups and resolves the run as a
comparison analysis.

## What happens with more than two groups

Datasets with more than two groups can still be used, but the samplesheet must
clearly define the group labels and intended comparisons.

Example:

```csv
sample_name,group,comparisons,group1,vs,group2
sample1,Treatment1,,Treatment1,vs,Control
sample2,Treatment1,,Treatment2,vs,Control
sample3,Treatment1,,,,
sample4,Treatment2,,,,
sample5,Treatment2,,,,
sample6,Treatment2,,,,
sample7,Control,,,,
sample8,Control,,,,
sample9,Control,,,,
```

This design contains three groups:

```text
Treatment1
Treatment2
Control
```

and two requested comparisons:

```text
Treatment1 vs Control
Treatment2 vs Control
```

!!! important

    Group names in the comparison section must match the group names in the
    `group` column exactly.

## Relationship with metadata

The analysis mode is controlled by the samplesheet groups, not by the optional
metadata file.

The metadata file can provide additional variables for downstream analysis or
interpretation, but it does not replace the required `group` column in the
samplesheet.

Example using metadata:

```bash
bash ~/MTD/MTD_explorer.sh \
  --input samplesheet.csv \
  --metadata metadata.csv \
  --output MTD_results_with_metadata \
  --hostid 9606 \
  --analysis-mode auto \
  --threads 16
```

## Common mistakes

### Forcing comparison mode with one group

This will fail:

```bash
bash ~/MTD/MTD_explorer.sh \
  --input one_group_samplesheet.csv \
  --output MTD_results \
  --hostid 9606 \
  --analysis-mode comparison
```

Use exploratory mode instead:

```bash
bash ~/MTD/MTD_explorer.sh \
  --input one_group_samplesheet.csv \
  --output MTD_results \
  --hostid 9606 \
  --analysis-mode exploratory
```

### Inconsistent group names

These are treated as different groups:

```text
Control
control
CTRL
Control_1
```

Use one consistent group name throughout the samplesheet.

### Confusing exploratory mode with incomplete analysis

Exploratory mode is not a failed or incomplete run.

It is a deliberate mode for non-comparison analyses. It still performs core
profiling steps but skips comparison-dependent steps that require a valid group
contrast.

## Recommended workflow

For a new dataset, start with:

```bash
--analysis-mode auto
```

After the first successful run, inspect the outputs and decide whether the
study design supports formal group comparisons.

For a single-group dataset, use:

```bash
--analysis-mode exploratory
```

For a validated two-group or multi-group design, use:

```bash
--analysis-mode comparison
```

## Next step

After choosing the analysis mode, review the
[Command-line reference](command-line.md) for all available
options.

[mtd-explorer]: https://github.com/patrick-douglas/MTD
[fastp]: https://github.com/OpenGene/fastp
[kraken2]: https://ccb.jhu.edu/software/kraken2/index.shtml
[bracken]: https://github.com/jenniferlu717/Bracken
[humann]: https://github.com/biobakery/humann
[deseq2]: https://bioconductor.org/packages/DESeq2/
[halla]: https://github.com/biobakery/halla
