# MTD Explorer pipeline benchmark

This benchmark measures a complete `MTD_explorer.sh` analysis in either local
or [Slurm](https://slurm.schedmd.com/) HPC mode. The large scientific output
remains separate from the compact benchmark bundle.

`MTD_benchmark_install.sh` remains the low-level monitor for the controller
machine. The pipeline-specific reporters add MTD stage timings, output
inventory, diagnostics and, in HPC mode, Slurm accounting and per-node hardware
records.

## Files

- `run_mtd_pipeline_benchmark.sh`
  - validates paths and benchmark options;
  - selects `local` or `hpc` execution;
  - removes only the previous analysis output after confirmation;
  - runs the complete MTD Explorer analysis;
  - creates the final benchmark bundle.
- `MTD_benchmark_install.sh`
  - records end-to-end wall time and resource use on the controller machine.
- `MTD_pipeline_benchmark_report.py`
  - creates the common local/HPC pipeline summary;
  - calculates durations between MTD progress stages;
  - inventories output files;
  - creates a one-row comparison table;
  - produces focused diagnostics after a failed run.
- `MTD_slurm_benchmark_report.py`
  - discovers all MTD Slurm stage job histories under the pipeline output;
  - queries `sacct` for jobs, array tasks, queue time, elapsed time, CPU and RSS;
  - combines accounting with MTD task markers;
  - summarizes performance by stage and by compute node.
- `MTD_explorer.sh`
  - writes progress timestamps only when
    `MTD_PIPELINE_BENCHMARK_STEPS_TSV` is set;
  - normal runs are unchanged.

## Execution modes

The default is:

```text
--execution-mode auto
```

Resolution rules:

```text
auto without --hpc-conf  -> local
auto with --hpc-conf     -> hpc
local                     -> rejects --hpc-conf
hpc                       -> requires --hpc-conf
```

`--threads` always describes the controller thread setting passed to
`MTD_explorer.sh`. In local mode this is also the pipeline thread count. In HPC
mode the Slurm jobs determine their allocated/detected resources independently.

## Local benchmark

```bash
cd ~/MTD-Explorer

bash benchmark/run_mtd_pipeline_benchmark.sh \
  --execution-mode local \
  --machine master \
  --dataset Bglabrata_PRJNA1306560 \
  --input /path/to/B.glabrata_fastq/samplesheet.csv \
  --output "$HOME/test_MTD_explorer_B.glabrata.local" \
  --hostid 6526 \
  --threads 20 \
  --run-number 1
```

The old command remains valid because `auto` selects local mode when no HPC
configuration is supplied:

```bash
bash benchmark/run_mtd_pipeline_benchmark.sh \
  --machine master \
  --dataset Bglabrata_PRJNA1306560 \
  --input /path/to/samplesheet.csv \
  --output "$HOME/test_MTD_explorer_B.glabrata.local" \
  --hostid 6526 \
  --threads 20 \
  --run-number 1
```

## Reference official-clean local run

The current documented local reference is the 15-sample paired-end
[PRJNA1306560](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1306560)
*Biomphalaria glabrata* analysis using host TaxID 6526, 20 local threads,
[Magic-BLAST](https://ncbi.github.io/magicblast/), automatic analysis/read-layout
selection, and top-5 microbiome-read extraction.

Benchmark identity:

```text
benchmark_run_kind       official_clean
resume_heavy             0
official_clean_benchmark 1
Git commit               04d61d412955ac7b9142ca8ebcc4d2edeb6bac5d
```

Observed results:

| Metric | Value |
|---|---:|
| Status | PASS |
| Benchmark interval | 21–24 August 2026 (UTC) |
| Wall time | 197,894.59 s (54 h 58 min 15 s) |
| Peak process-tree RSS | 117.27 GiB |
| Peak system memory used | 117.52 GiB |
| Mean system CPU busy | 65.67% |
| Pipeline output | 10,821 files, ~1.39 TB |
| Diagnostic hits | 0 |

The reference controller used [Linux Mint](https://linuxmint.com/) 21.1, an
Intel Core i9-10900K with 20 logical CPUs, 128 GB-class RAM, and an NVMe
main/output filesystem. Input data
were stored on a SATA HDD.

The largest stages were [HUMAnN](https://github.com/biobakery/humann)
functional profiling (30 h 39 min; 55.8%), host read/expression processing
(11 h 53 min; 21.6%), and ssGSEA (4 h 55 min; 9.0%).
These values describe one workload and are not minimum hardware or storage
requirements.

The public documentation contains the full benchmark configuration and
interpretation in `docs/user-guide/benchmarking.md`.

## Additional exploratory reference measurement

The six-library *Danio rerio* case (PRJNA1005695) was run locally at commit `5be26aa2e03938abb12fe00a37a08b0c0cef3777` with 20 threads, HISAT2, Bracken 150 bp and `--analysis-mode exploratory`. It completed on 25 November 2026 in 15,173 s (4 h 12 min 53 s), with GNU time maximum reported RSS of 113,107,988 KiB, 1089% CPU utilization and exit status 0. The output contained 678 files (~459 GiB).

This was a fresh-output execution timed with GNU `/usr/bin/time -v` and stage markers, not an `official_clean` benchmark-wrapper run. Reference preparation was measured separately. Do not compare maximum RSS or CPU percentage directly with process-tree or sampled-system metrics above, or infer a speed advantage across these distinct workloads. See the [exploratory case study](https://patrick-douglas.github.io/MTD-Explorer/user-guide/danio-rerio/) and [benchmark notes](https://patrick-douglas.github.io/MTD-Explorer/user-guide/benchmarking/#reference-exploratory-run-prjna1005695).

## HPC benchmark

```bash
cd ~/MTD-Explorer

bash benchmark/run_mtd_pipeline_benchmark.sh \
  --execution-mode hpc \
  --machine LBN_Slurm \
  --dataset Bglabrata_PRJNA1306560 \
  --input /path/to/B.glabrata_fastq/samplesheet.csv \
  --output "$HOME/test_MTD_explorer_B.glabrata.hpc" \
  --hostid 6526 \
  --threads 20 \
  --hpc-conf "$HOME/MTD-Explorer/Installation/HPC/MTD_hpc_slurm.conf" \
  --run-number 1
```

The end-to-end wall time remains the principal user-facing duration. It includes
controller work, Slurm queueing, retries, node-local staging, distributed
processing and final local analyses.

The Slurm compute makespan is reported separately. It spans the earliest start
to the latest end among the collected HPC task attempts and must not be confused
with the complete pipeline wall time.

When Slurm accounting is unavailable or incomplete, the task markers still
provide node, thread, memory and elapsed-time information. Such runs are labeled
`PASS_WITH_PARTIAL_HPC_METRICS` rather than silently presenting incomplete
accounting as a fully collected benchmark.

## Distributed stages

The collector recognizes the stage names written by the current HPC backend:

```text
fastp
kraken_host
kraken_micro_raw
kraken_micro_final
bracken
humann
magicblast
```

Magic-BLAST may produce many task attempts because samples are divided into
FASTQ chunks.

## Read extraction

The standard benchmark disables microbiome-read extraction. To benchmark the
optional extraction function:

```bash
--extract --top-n 10
```

Use `--extract --top-n 0` only as an explicit stress test because every detected
taxon can require repeated scans and substantial storage.

## Main outputs

Every run contains the common files:

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
failure_report.txt        # only when incomplete or failed
```

The raw `console.log` is retained in the local benchmark directory but excluded
from the compressed bundle. The compact `console_clean.log` is included instead,
and `bundle_manifest.tsv` records this policy.

An HPC run additionally contains:

```text
hpc_configuration.conf
hpc_summary.tsv
hpc_jobs.tsv
hpc_tasks.tsv
hpc_stage_summary.tsv
hpc_node_summary.tsv
hpc_collection_warnings.txt   # only when accounting is partial
```

The pipeline output also retains per-attempt task markers under each Slurm
stage `success/attempts/` directory. These preserve timings and node hardware
for retries instead of keeping only the final task outcome.

### `hpc_stage_summary.tsv`

Reports, for each distributed stage:

- Slurm jobs and submitted task attempts;
- unique, successful and failed tasks;
- retries and nodes used;
- earliest submit/start and latest end;
- parallel compute makespan;
- sum of task elapsed time;
- allocated and reported consumed CPU time;
- peak task RSS and CPUs per task.

### `hpc_node_summary.tsv`

Reports only the nodes that actually executed task attempts, including CPU
model, architecture, topology, memory, kernel, current Slurm state, stages,
task counts and accumulated processing metrics. The list is discovered from
the per-attempt markers written on the selected compute nodes, with `sacct`
used as a fallback. It is never built from a fixed node list, `sinfo`, or the
set of nodes currently eligible in the partition. Consequently, nodes in
`DRAIN`, `DOWN`, or otherwise unused during that run do not appear; a node that
executed a task remains in the report even if it is drained later. Hardware
values are recorded by the real MTD jobs and supplemented, when available,
with `scontrol show node` metadata.

### `benchmark_comparison_row.tsv`

This one-row file is intended for combining repetitions and comparing local and
HPC scenarios. It includes execution mode, machine/cluster label, controller
threads, nodes used, end-to-end wall time, Slurm compute makespan, status and Git
commit.

## Fair comparisons

For local-versus-HPC documentation:

- use the same Git commit;
- use identical FASTQ files, samplesheet, databases and analysis options;
- use separate clean output directories;
- preserve both successful and failed benchmark bundles;
- report the end-to-end wall time for both scenarios;
- report Slurm compute makespan and node inventory only for HPC;
- do not compare the sum of parallel task durations directly with local wall
  time.

The bundles are created under:

```text
~/MTD_pipeline_benchmarks/
```
