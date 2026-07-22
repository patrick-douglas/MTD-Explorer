# MTD Explorer installation benchmark suite

This package contains three scripts:

1. `MTD_benchmark_install.sh`  
   Runs the installer and measures the entire machine during the installation:
   wall time, CPU, RAM, process-tree RSS, disk I/O, network I/O, temperature,
   hardware, operating system, Git commit, hashes, logs and watched-directory
   sizes.

2. `MTD_make_instrumented_installer.sh`  
   Creates `Install_profiled.sh` without changing the original `Install.sh`.
   It records function-level timings in `steps.tsv`.

3. `MTD_benchmark_merge.py`  
   Combines runs from multiple machines into article-ready TSV tables.

## 1. Place the scripts in the MTD Explorer repository

```bash
cd ~/MTD-Explorer
cp /path/to/MTD_benchmark_suite/* .
chmod +x \
    benchmark/MTD_benchmark_install.sh \
    benchmark/MTD_make_instrumented_installer.sh \
    benchmark/MTD_fix_profiled_locale.sh \
    benchmark/run_mtd_clean_benchmark.sh
```

## 2. Create the profiled installer

```bash
bash ./benchmark/MTD_make_instrumented_installer.sh \
    --input ./Install.sh \
    --output ./Install_profiled.sh
```

Validate both scripts:

```bash
bash -n ./Install.sh
bash -n ./Install_profiled.sh
```

The original `Install.sh` is not modified.

## 3. Run a clean installation benchmark

Example for the master machine:

```bash
bash ./benchmark/MTD_benchmark_install.sh \
    --label master_cold_native_r1 \
    --interval 5 \
    --output-root "$HOME/MTD_benchmarks" \
    --watch-path "$HOME/miniconda3" \
    --watch-path "/media/me/MTD_install_cache" \
    --watch-path "$HOME/MTD-Explorer" \
    -- \
    bash ./Install_profiled.sh \
        -o "/media/me/MTD_install_cache"
```

Example for the independent lower-performance machine:

```bash
bash ./benchmark/MTD_benchmark_install.sh \
    --label secondary_cold_native_r1 \
    --interval 5 \
    --output-root "$HOME/MTD_benchmarks" \
    --watch-path "$HOME/miniconda3" \
    --watch-path "/path/to/MTD_install_cache" \
    --watch-path "$HOME/MTD-Explorer" \
    -- \
    bash ./Install_profiled.sh \
        -o "/path/to/MTD_install_cache"
```

Do not pass a password through `Install.sh -w` during a benchmark because the
executed command is recorded in the metadata. Let `sudo` prompt normally.

## 4. Output of each run

A run directory contains:

```text
console.log
git_state.txt
gnu_time.txt
hardware.txt
metadata.txt
resource_samples.csv
software.txt
steps.tsv
summary.tsv
watch_paths_before.tsv
watch_paths_after.tsv
```

`summary.tsv` contains whole-installation metrics.  
`resource_samples.csv` contains one resource sample every few seconds.  
`steps.tsv` contains function-level timings from the profiled installer.

Nested function durations overlap. For stage reporting, inspect
`parent_function`, `call_depth`, and the major installation functions rather
than summing every row blindly.

## 5. Merge results from both machines

Copy both benchmark root directories to one machine, then run:

```bash
python3 ./benchmark/MTD_benchmark_merge.py \
    --input "$HOME/MTD_benchmarks_master" \
    --input "$HOME/MTD_benchmarks_secondary" \
    --output "$HOME/MTD_benchmark_merged"
```

The merged directory contains:

```text
runs_wide.tsv
runs_article.tsv
steps_long.tsv
steps_summary.tsv
```

`runs_article.tsv` is the compact hardware/performance table intended for the
paper. `runs_wide.tsv` preserves every collected metric.

## 6. Recommended labels

```text
master_cold_native_r1
secondary_cold_native_r1

master_warm_native_r1
secondary_warm_native_r1

master_warm_t8_r1
master_warm_t8_r2
master_warm_t8_r3

secondary_warm_t8_r1
secondary_warm_t8_r2
secondary_warm_t8_r3
```

- `cold`: no Miniconda and no populated MTD cache.
- `warm`: downloaded inputs are already present in the persistent cache.
- `native`: the installer uses the machine's normal thread count.
- `t8`: both machines use eight threads after the installer gains a controlled
  thread option or `MTD_THREADS` environment variable.

## 7. Fair-run checklist

Before each timed run:

- use the same Git commit;
- record whether the cache is empty or populated;
- reboot the machine;
- stop unrelated analyses, backups and cloud synchronization;
- use the same network type;
- avoid running both large downloads simultaneously on the same connection;
- keep the machine connected to AC power;
- preserve every run directory, including failed runs;
- run `MTD_check_installation.sh --deep` after a successful installation.

The clean installation on each machine demonstrates reproducibility.
Repeated warm-cache component benchmarks provide the statistical performance
comparison.


## Important path rule

Do not quote a path beginning with a literal tilde:

```bash
# Wrong: the shell passes "~" literally
-o "~/MTD_install_cache"

# Correct
-o "$HOME/MTD_install_cache"

# Also correct
-o ~/MTD_install_cache
```

The benchmark wrapper now stops immediately if it detects a literal `~` path.

## Locale hotfix for profiler versions generated before this correction

Some locales format Bash `EPOCHREALTIME` with a decimal comma. If an older
`Install_profiled.sh` reports `value too great for base`, run:

```bash
bash ./benchmark/MTD_fix_profiled_locale.sh ./Install_profiled.sh
```

Alternatively, replace the suite scripts and regenerate `Install_profiled.sh`
from the unchanged original `Install.sh`.
