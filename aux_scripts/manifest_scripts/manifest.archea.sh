#!/usr/bin/env bash
set -euo pipefail

offline_files_folder="/media/me/4TB_BACKUP_LBN/Compressed/MTD"
offline_files_folder="${MTD_OFFLINE_FILES_FOLDER:-$offline_files_folder}"

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
helper="${MTD_MANIFEST_HELPER:-$script_dir/sync_ncbi_cache.py}"

[[ -f "$helper" ]] || {
    echo "[FAIL] NCBI cache synchronization helper not found: $helper" >&2
    exit 1
}

args=(
    assemblies
    --library archaea
    --summary-url "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/archaea/assembly_summary.txt"
    --local-dir "$offline_files_folder/Kraken2DB_micro/library/archaea/all"
    --metadata-dir "$offline_files_folder/Kraken2DB_micro/library/archaea"
    --level "Complete Genome"
    --level "Chromosome"
    --min-count 100
    --require-complete
)

[[ "${FULL_GZIP_CHECK:-0}" == "1" ]] && args+=(--full-gzip-check)

exec python3 "$helper" "${args[@]}"
