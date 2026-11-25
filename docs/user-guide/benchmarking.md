# Pipeline benchmarking

MTD Explorer includes an end-to-end benchmark workflow for measuring complete
pipeline runs in either ordinary local mode or the optional [Slurm][slurm] HPC
mode.

The benchmark is designed for reproducible performance comparisons. It keeps
large scientific outputs in the normal MTD Explorer output directory and writes
a separate compact benchmark bundle containing commands, timings, resource
summaries, input records, output inventories, and diagnostics.

!!! important "Benchmark the complete pipeline"

    For local-versus-HPC comparisons, use the benchmark wrapper rather than
    timing individual commands manually. The wrapper records the complete
    user-facing wall time and, for HPC runs, collects Slurm-specific metrics
    separately.

## Main benchmark command

The pipeline benchmark is launched with:

```bash
bash benchmark/run_mtd_pipeline_benchmark.sh [options]
```

Display the current options with:

```bash
bash benchmark/run_mtd_pipeline_benchmark.sh --help
```

The benchmark requires:

| Option | Meaning |
|---|---|
| `--machine NAME` | Label for the machine or cluster being measured |
| `--input FILE` | Samplesheet used by the pipeline |
| `--output DIR` | Clean MTD Explorer scientific output directory |
| `--hostid TAXID` | Host NCBI Taxonomy ID |

Frequently used additional options include:

| Option | Default | Meaning |
|---|---:|---|
| `--dataset LABEL` | samplesheet parent directory | Human-readable dataset label |
| `--threads INT` | `nproc` | Controller thread setting passed to `MTD_explorer.sh` |
| `--run-number INT` | `1` | Repetition number |
| `--read-layout MODE` | `auto` | `auto`, `se`, or `pe` |
| `--analysis-mode MODE` | `auto` | `auto`, `comparison`, or `exploratory` |
| `--interval SECONDS` | `5` | Controller resource-sampling interval |
| `--benchmark-root DIR` | `$HOME/MTD_pipeline_benchmarks` | Root for compact benchmark bundles |
| `--execution-mode MODE` | `auto` | `auto`, `local`, or `hpc` |
| `--hpc-conf FILE` | none | Slurm configuration used for an HPC benchmark |
| `--extract` | disabled | Enable optional detected-microbiome read extraction |
| `--top-n INT` | `10` when extraction is enabled | Number of ranked taxa extracted; `0` means all detected taxa |
| `--hisat2` | [Magic-BLAST][magic-blast] enabled | Use [HISAT2][hisat2] instead of Magic-BLAST |
| `--yes` | disabled | Remove an existing benchmark pipeline output without prompting |

Additional MTD Explorer arguments can be passed after `--`:

```bash
bash benchmark/run_mtd_pipeline_benchmark.sh \
  --machine workstation \
  --input /path/to/samplesheet.csv \
  --output "$HOME/mtd_benchmark_output" \
  --hostid 9606 \
  -- \
  --metadata /path/to/metadata.csv \
  --pdm spearman
```

## Execution modes

The benchmark has three execution modes:

```text
auto
local
hpc
```

The resolution rules are:

```text
auto without --hpc-conf  -> local
auto with --hpc-conf     -> hpc
local                     -> rejects --hpc-conf
hpc                       -> requires --hpc-conf
```

Therefore, an existing benchmark command that does not specify
`--execution-mode` or `--hpc-conf` continues to run locally.

### Meaning of `--threads`

`--threads` always records and passes the controller thread setting to
`MTD_explorer.sh`.

In **local mode**, this is also the main pipeline thread count.

In **HPC mode**, distributed Slurm tasks use the resources allocated or detected
on the compute node according to the HPC configuration. The benchmark therefore
does not assume that `--threads 20`, for example, means that every heterogeneous
compute node receives 20 CPUs.

## Local benchmark

Example:

```bash
cd ~/MTD-Explorer

bash benchmark/run_mtd_pipeline_benchmark.sh \
  --execution-mode local \
  --machine workstation \
  --dataset example_dataset \
  --input /path/to/samplesheet.csv \
  --output "$HOME/MTD_benchmark_local" \
  --hostid 9606 \
  --threads 20 \
  --run-number 1
```

Because `auto` resolves to local mode when no HPC configuration is supplied,
`--execution-mode local` may be omitted when an explicit mode is not needed.

## Reference local benchmark: PRJNA1306560

A complete local reference benchmark was performed with the public
[PRJNA1306560][bioproject] *Biomphalaria glabrata* dataset described in the
[example dataset guide](example-dataset.md). The run used
all **15 paired-end RNA-seq samples** and completed successfully as an
**official clean** benchmark.

The benchmark measures the analysis workflow with the host reference,
databases, and software environment already installed. It does **not** include
installation, database construction, or custom-host reference construction.

### Benchmark configuration

| Setting | Reference run |
|---|---|
| Dataset | [PRJNA1306560][bioproject] |
| Host | [*Biomphalaria glabrata*][ncbi-taxonomy] |
| Host TaxID | `6526` |
| Samples | 15 paired-end RNA-seq samples |
| Controller / local threads | 20 |
| Read-layout request | `auto` (resolved to paired-end) |
| Analysis mode | `auto` |
| Host alignment | [Magic-BLAST][magic-blast] |
| Microbiome-read extraction | enabled, top 5 ranked taxa |
| Benchmark run kind | `official_clean` |
| Benchmark interval | 21–24 August 2026 (UTC) |
| `resume_heavy` | `0` |
| Git commit | `04d61d412955ac7b9142ca8ebcc4d2edeb6bac5d` |

The corresponding benchmark-wrapper configuration was equivalent to:

```bash
bash benchmark/run_mtd_pipeline_benchmark.sh \
  --machine local_linux \
  --dataset PRJNA1306560_Bglabrata_extract_top5 \
  --input /path/to/B.glabrata_fastq/samplesheet.csv \
  --output /path/to/clean/MTD_output \
  --hostid 6526 \
  --threads 20 \
  --read-layout auto \
  --analysis-mode auto \
  --extract \
  --top-n 5 \
  --run-number 1
```

### Benchmark hardware

| Resource | Reference system |
|---|---|
| Operating system | [Linux Mint][linux-mint] 21.1 |
| Kernel | Linux 6.8.0-134-generic |
| CPU | Intel Core i9-10900K @ 3.70 GHz |
| Logical CPUs | 20 |
| Physical memory | 128 GB class; Linux `MemTotal` 131,747,584 kB (~125.6 GiB) |
| Main/output filesystem | 4 TB-class XPG GAMMIX S70 BLADE NVMe SSD |
| Input-data filesystem | SATA HDD |

Storage location can materially affect runtime. The hardware table is therefore
part of the benchmark context rather than a general hardware recommendation.

### Reference results

| Metric | Observed value |
|---|---:|
| Pipeline status | **PASS** |
| Exit status | `0` |
| Completion marker | yes |
| Diagnostic hits | `0` |
| End-to-end wall time | **197,894.59 s (54 h 58 min 15 s)** |
| Peak process-tree RSS | **117.27 GiB** |
| Peak system memory used | **117.52 GiB** |
| Mean system CPU busy | **65.67%** |
| Peak process-tree CPU | **2022.18% of one core** |
| Controller physical disk read during run | **4,883.45 GiB** |
| Controller physical disk write during run | **2,871.87 GiB** |
| Pipeline output files | **10,821** |
| Pipeline output directories | **281** |
| Pipeline output size | **1.39 TB** (1,392,742,248,591 bytes) |
| Recorded pipeline stages | **16** |

The approximately 1.39 TB output footprint is specific to this run and its
enabled analyses, including top-5 microbiome-read extraction. It should not be
interpreted as a fixed storage requirement for every MTD Explorer analysis.

### Stage timing

The largest measured stages were:

| Stage | Duration | Share of wall time |
|---|---:|---:|
| [HUMAnN][humann] functional profiling | 30 h 39 min | 55.8% |
| Host-read processing and host expression | 11 h 53 min | 21.6% |
| ssGSEA pathway enrichment | 4 h 55 min | 9.0% |
| Preparing raw reads | 2 h 22 min | 4.3% |
| Optional contaminant removal | 1 h 56 min | 3.5% |
| [Kraken2][kraken2] / [Bracken][bracken] visualization preparation | 58 min | 1.8% |
| Host-read classification with Kraken2 | 38 min | 1.2% |
| Reclassification after contaminant removal | 28 min | 0.9% |
| Non-host classification with Kraken2 | 26 min | 0.8% |
| Host-gene association analysis with [HAllA][halla] | 18 min | 0.6% |
| Host-pathway association analysis with HAllA | 3 min | 0.1% |

Short bookkeeping and Bracken-combination stages account for the remaining
fraction. For this workload, HUMAnN plus host-read/host-expression processing
accounted for approximately **77.4%** of the complete wall time.

!!! important "How to interpret this benchmark"

    These values are a reproducible reference measurement, not a performance
    guarantee or a minimum system specification. Runtime, memory, I/O, and
    output size can change substantially with sequencing depth, sample count,
    host-genome complexity, databases, analysis options, storage performance,
    and available CPU resources.

    The benchmark bundle records `official_clean_benchmark=1`, `resume_heavy=0`,
    a clean Git state, and the exact MTD Explorer commit used for the run. The
    bundle does not independently record operating-system page-cache flushing,
    so the documented result is described as an **official clean pipeline run**
    rather than asserting a cold filesystem-cache state.

## Reference exploratory run: PRJNA1005695

A second reference measurement used six paired-end intestinal libraries of *Danio rerio* from [PRJNA1005695](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1005695). See the [exploratory case study](danio-rerio.md) for the exact input accessions, references, command and source tables.

| Setting | Exploratory reference run |
|---|---|
| Date | 25 November 2026 (UTC) |
| Analytical commit | `5be26aa2e03938abb12fe00a37a08b0c0cef3777` |
| Host TaxID | 7955 |
| Analysis mode | `exploratory` |
| Input | Six paired-end runs; 211,733,551 original read pairs |
| Host alignment | HISAT2 |
| Bracken read length | 150 bp |
| Optional contaminant filtering | Disabled for this case |
| Taxon-specific read extraction | Not requested |
| CPU | Intel Core i9-10900K, 20 logical CPUs |
| Threads | 20 |
| Physical memory | 128 GB class (~125 GiB available to Linux) |
| Measurement | GNU `/usr/bin/time -v` and pipeline stage markers |
| Run kind | Fresh output directory; `--resume-heavy` enabled but no checkpoints reused |

| Metric | Observed value |
|---|---:|
| Exit status | 0 |
| Wall time | 15,173 s (4 h 12 min 53 s) |
| Maximum reported RSS | 113,107,988 KiB (~107.87 GiB) |
| CPU utilization (GNU time) | 1089% |
| Swaps reported | 0 |
| Output files | 678 |
| Empty output files | 0 |
| Output footprint (`du -h`) | ~459 GiB |

Stage intervals below are elapsed times between progress markers, not exclusive per-tool measurements. Short intervals around contaminant handling, Bracken and visualization preparation are grouped. Initial and final overhead outside the marker interval is not included in the table.

| Stage interval | Approximate duration |
|---|---:|
| Preparing raw reads | 0 h 56 min 27 s |
| Kraken2 host classification | 0 h 14 min 51 s |
| Kraken2 non-host classification | 0 h 11 min 33 s |
| Read summary, optional filter handling, Bracken and visualization preparation | 0 h 01 min 21 s |
| HUMAnN and associated functional processing | 1 h 33 min 40 s |
| Host-read processing and host expression | 1 h 14 min 16 s |

The [stage timing TSV](../assets/data/example-dataset/danio-rerio/pipeline-stage-timing.tsv) provides the original timestamps.

### Separately measured reference preparation

| Process | Wall time | Maximum reported RSS (KiB) |
|---|---:|---:|
| Custom host reference construction | 1 h 30 min 05 s | 93,835,952 |
| Add Bracken 150-base distribution | 49 min 05 s | 115,804,492 |

These one-time preparation measurements are excluded from the analytical runtime; the microbial Kraken2 database already existed. The definitive run reused the host resources and new Bracken distribution.

!!! important "Do not interpret these cases as a controlled benchmark comparison"

    The *B. glabrata* run used the official instrumented benchmark wrapper, which records aggregate process-tree RSS and sampled system CPU use. The *D. rerio* run used GNU time, whose maximum RSS is not a sampled peak sum across the process tree; its CPU percentage is also a different metric from sampled mean system CPU busy. The cases differ in sample counts, sequencing depths, reference resources, alignment tools and enabled stages. The fresh output directory does not imply a cold OS cache. Neither run demonstrates a controlled speed advantage of one analysis mode over the other.

## HPC benchmark

First prepare and validate the HPC runtime as described in
[HPC / Slurm execution](hpc-slurm.md).

Then run the same scientific analysis with an HPC configuration:

```bash
cd ~/MTD-Explorer

bash benchmark/run_mtd_pipeline_benchmark.sh \
  --execution-mode hpc \
  --machine cluster_label \
  --dataset example_dataset \
  --input /path/to/samplesheet.csv \
  --output "$HOME/MTD_benchmark_hpc" \
  --hostid 9606 \
  --threads 20 \
  --hpc-conf "$HOME/MTD-Explorer/Installation/HPC/MTD_hpc_slurm.conf" \
  --run-number 1
```

The HPC collector recognizes the distributed stage names produced by the
current backend:

```text
fastp
kraken_host
kraken_micro_raw
kraken_micro_final
bracken
humann
magicblast
```

Magic-BLAST may contribute many task attempts because synchronized FASTQ input
can be split into multiple chunks before submission.

## Which time should be reported?

For both local and HPC scenarios, the primary comparison metric is the
**end-to-end pipeline wall time**.

For an HPC run, this duration includes:

- controller-side work;
- Slurm queueing;
- retries;
- input stage-in;
- distributed processing;
- output stage-out;
- final pipeline stages that remain local.

The benchmark also reports a **Slurm compute makespan**. This spans the earliest
start to the latest end among collected HPC task attempts. It is useful for
understanding the distributed part of the run, but it is not equivalent to the
complete pipeline wall time.

!!! warning "Do not compare summed parallel task time with local wall time"

    Multiple HPC tasks execute concurrently. The sum of their individual
    elapsed times is therefore not a valid replacement for end-to-end wall
    time when comparing local and distributed execution.

## Read-extraction benchmark behavior

The standard pipeline benchmark **does not extract detected microbiome reads by
default**.

To include the optional extraction feature:

```bash
--extract --top-n 10
```

When extraction is enabled and `--top-n` is not changed, the benchmark uses the
top 10 taxa from the absolute ranking.

Use:

```bash
--extract --top-n 0
```

only as an intentional stress test. Extracting every detected taxon can require
many repeated FASTQ scans, substantially more runtime, and additional storage.

## Benchmark bundle location

By default, compact benchmark bundles are stored under:

```text
~/MTD_pipeline_benchmarks/
```

The scientific pipeline output remains in the directory supplied with
`--output`.

This separation allows the benchmark bundle to remain small enough for
archiving and comparison while preserving the full analysis results elsewhere.

The raw `console.log` remains in the local benchmark directory but is excluded
from the compressed bundle. The bundle contains the compacted
`console_clean.log`, which collapses carriage-return progress updates instead of
archiving every terminal redraw. `bundle_manifest.tsv` records this policy.

## Common benchmark outputs

Every completed benchmark bundle can contain the following common files:

```text
benchmark_comparison_row.tsv
benchmark_run_mode.tsv
bundle_manifest.tsv
console_clean.log
diagnostic_hits.tsv
final_console_tail.txt
hardware.txt
input_files.tsv
input_samplesheet.csv
metadata.txt
output_extensions.tsv
output_inventory.tsv
pipeline_command.sh
pipeline_steps.tsv
pipeline_steps_raw.tsv
pipeline_summary.tsv
pipeline_summary.txt
resource_samples.csv
software.txt
summary.tsv
failure_report.txt          # only when incomplete or failed
```

Key files include:

| File | Purpose |
|---|---|
| `pipeline_command.sh` | Exact benchmarked pipeline command |
| `input_samplesheet.csv` | Copy of the benchmark samplesheet |
| `input_files.tsv` | Input-file inventory |
| `pipeline_steps.tsv` | Processed pipeline-stage timing table |
| `pipeline_summary.tsv` | Machine-readable benchmark summary |
| `pipeline_summary.txt` | Human-readable benchmark summary |
| `benchmark_run_mode.tsv` | Records `official_clean` versus development/resume execution state |
| `bundle_manifest.tsv` | Records which console logs are retained locally or included in the archive |
| `resource_samples.csv` | Time-series controller CPU, memory, I/O, network, and temperature samples |
| `hardware.txt` | Controller operating-system, CPU, memory, storage, and network inventory |
| `software.txt` | Benchmark wrapper and available controller software metadata |
| `output_inventory.tsv` | Inventory of produced scientific outputs |
| `benchmark_comparison_row.tsv` | One-row record intended for merging repetitions and comparing scenarios |
| `failure_report.txt` | Focused failure diagnostics when a run is incomplete or fails |

## Additional HPC benchmark outputs

An HPC benchmark additionally records files such as:

```text
hpc_configuration.conf
hpc_summary.tsv
hpc_jobs.tsv
hpc_tasks.tsv
hpc_stage_summary.tsv
hpc_node_summary.tsv
hpc_collection_warnings.txt   # only when collection is partial
```

The pipeline output also retains per-attempt HPC markers below the individual
stage `success/attempts/` directories. These preserve retry history, task timing,
and node hardware information instead of retaining only the final task result.

### `hpc_stage_summary.tsv`

This table summarizes each distributed stage, including:

- Slurm jobs and submitted task attempts;
- unique, successful, and failed tasks;
- retries;
- nodes used;
- earliest submission/start and latest end;
- parallel compute makespan;
- summed task elapsed time;
- allocated and reported CPU time when available;
- peak task RSS;
- CPUs per task.

### `hpc_node_summary.tsv`

This table reports **only compute nodes that actually executed task attempts in
that benchmark run**.

The node list is reconstructed from the per-attempt markers written by real MTD
Explorer jobs, with Slurm accounting used as a fallback. It is not generated
from a fixed node list, `sinfo`, or all nodes currently eligible in a partition.

Consequently:

- a node that remained `DRAIN`, `DOWN`, or otherwise unused does not appear;
- a node that executed a task remains part of the benchmark record even if its
  Slurm state changes later;
- heterogeneous node hardware can be reported according to the nodes that
  actually participated.

Node records can include CPU model, architecture, topology, memory, kernel,
Slurm state, executed stages, task counts, and accumulated processing metrics.

### `benchmark_comparison_row.tsv`

This one-row table is intended for merging benchmark repetitions. It records
fields such as:

- execution mode;
- machine or cluster label;
- controller threads;
- nodes used;
- end-to-end wall time;
- Slurm compute makespan when applicable;
- benchmark status;
- Git commit.

## Partial Slurm accounting

HPC task markers remain useful even when `sacct` is unavailable or returns
incomplete accounting information. The collector can still recover task-level
node, thread, memory, and elapsed-time information from the MTD Explorer HPC
records.

A successful run with incomplete Slurm accounting is labeled:

```text
PASS_WITH_PARTIAL_HPC_METRICS
```

rather than presenting partial metrics as though full accounting had been
collected.

Inspect `hpc_collection_warnings.txt` when that status is reported.

## Fair local-versus-HPC comparisons

Use the following rules when creating performance numbers for documentation,
reports, or manuscripts:

1. use the same MTD Explorer Git commit;
2. use identical FASTQ files;
3. use the same samplesheet and metadata;
4. use identical databases and reference versions;
5. use the same analysis options;
6. use separate clean scientific output directories;
7. preserve benchmark bundles from both successful and failed runs;
8. compare end-to-end wall time for local and HPC scenarios;
9. report Slurm compute makespan and actual node inventory separately for HPC;
10. do not compare summed parallel task durations directly against local wall
    time.

These rules are especially important for heterogeneous clusters, where the set
of compute nodes actually used can differ between repetitions as Slurm node
availability changes.

## Benchmark implementation files

The main pipeline benchmark is implemented by:

```text
benchmark/run_mtd_pipeline_benchmark.sh
benchmark/MTD_benchmark_install.sh
benchmark/MTD_pipeline_benchmark_report.py
benchmark/MTD_slurm_benchmark_report.py
```

`run_mtd_pipeline_benchmark.sh` validates the requested scenario, runs the
complete pipeline, and creates the final bundle.

`MTD_benchmark_install.sh` records controller-side wall time and resource use.

`MTD_pipeline_benchmark_report.py` creates the common local/HPC pipeline
summary, stage timings, output inventory, comparison row, and failure
information.

`MTD_slurm_benchmark_report.py` discovers MTD Explorer Slurm job histories,
combines `sacct` information with the pipeline's own task markers, and produces
stage- and node-level HPC summaries.

Normal pipeline runs are not instrumented with benchmark timestamps unless the
benchmark-specific environment variable used by the wrapper is set.

## Related pages

- [HPC / Slurm execution](hpc-slurm.md)
- [Command-line reference](command-line.md)
- [Methods and reproducibility outputs](methods-reproducibility-outputs.md)
- [Output files](output-files.md)
- [Troubleshooting](../troubleshooting/index.md)

[slurm]: https://slurm.schedmd.com/
[bioproject]: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1306560
[ncbi-taxonomy]: https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=6526
[linux-mint]: https://linuxmint.com/
[magic-blast]: https://ncbi.github.io/magicblast/
[kraken2]: https://ccb.jhu.edu/software/kraken2/index.shtml
[humann]: https://github.com/biobakery/humann
[halla]: https://github.com/biobakery/halla
[hisat2]: https://daehwankimlab.github.io/hisat2/
[bracken]: https://ccb.jhu.edu/software/bracken/
