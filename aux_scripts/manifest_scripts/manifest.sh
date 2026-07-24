#!/usr/bin/env bash
set -euo pipefail

# Compatibility entry point retained for older MTD workflows.
# The current bacterial synchronizer always checks the live NCBI catalog.
offline_files_folder=""
offline_files_folder="${MTD_OFFLINE_FILES_FOLDER:-$offline_files_folder}"

if [[ -z "$offline_files_folder" ]]; then
    echo "[FAIL] offline_files_folder is not configured." >&2
    exit 1
fi

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
helper="${MTD_MANIFEST_HELPER:-$script_dir/sync_ncbi_cache.py}"

[[ -f "$helper" ]] || {
    echo "[FAIL] NCBI cache synchronization helper not found: $helper" >&2
    exit 1
}

args=(
    assemblies
    --library bacteria
    --summary-url "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt"
    --local-dir "$offline_files_folder/Kraken2DB_micro/library/bacteria/all"
    --metadata-dir "$offline_files_folder/Kraken2DB_micro/library/bacteria"
    --level "Complete Genome"
    --level "Chromosome"
    --min-count 1000
    --require-complete
)

[[ "${FULL_GZIP_CHECK:-0}" == "1" ]] && args+=(--full-gzip-check)

exec python3 "$helper" "${args[@]}"
