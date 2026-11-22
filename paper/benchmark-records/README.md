# Local benchmark records

The files in this folder are copies from the timed comparison-mode execution of MTD Explorer for PRJNA1306560. Supplementary Table S4 reports the software settings and identifiers, S5 the resource summaries, and S6 the final inventory.

| File | How to use it |
| --- | --- |
| `input_samplesheet.csv` | Identify the libraries, group labels and comparisons used in the run. See S4C for SRA run accessions. |
| `benchmark_run_mode.tsv` | Identify the clean-run and intermediate-reuse settings of the benchmark. |
| `pipeline_steps_raw.tsv` | Read the UTC progress events emitted during the run. |
| `pipeline_steps.tsv` | Read the 16 recorded intervals between progress events. Each duration includes all work in that interval. |
| `resource_samples.csv` | Inspect the 39,020 observations of CPU, memory, swap and other system measurements. The configured sampling interval was five seconds. |

Figure 7 uses system-wide CPU activity and the sampled sum of resident memory across monitored descendant processes. Multiple processes can count the same shared memory pages. The benchmark clock excludes software installation, reference preparation and the subsequent TPI robustness analysis.
