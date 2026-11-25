#!/usr/bin/env bash

set -euo pipefail

BASE="$HOME/MTD-Explorer_Exploratory_D.rerio"
REPO="$HOME/MTD-Explorer"

set +u
source "$HOME/miniconda3/etc/profile.d/conda.sh"

conda activate MTD
set -u

cd "$REPO"

mkdir -p "$BASE/logs" "$BASE/results" "$BASE/provenance"

export R_LIBS_USER="$REPO/custom_R_libs:${R_LIBS_USER:-}"

export MTD_PIPELINE_BENCHMARK_STEPS_TSV="$BASE/provenance/pipeline_stage_timing.tsv"

echo "Run started: $(date -Is)"

/usr/bin/time -v \
    bash MTD_explorer.sh \
    --input "$BASE/raw_fastq/samplesheet.csv" \
    --output "$BASE/results/Danio_rerio_exploratory" \
    --hostid 7955 \
    --analysis-mode exploratory \
    --read-layout pe \
    --threads 20 \
    --trim-length 40 \
    --bracken-read-len 150 \
    --bracken-threshold 10 \
    --kraken-host-confidence 0.05 \
    --kraken-host-min-hit-groups 3 \
    --kraken-micro-confidence 0.10 \
    --kraken-micro-min-hit-groups 3 \
    --ssgsea-gmt auto \
    --resume-heavy \
    2>&1 | tee -a "$BASE/logs/MTD_Danio_rerio_exploratory.log"

echo "Run completed: $(date -Is)"
