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

metadata_dir="$offline_files_folder/Kraken2DB_micro/library/viral"

args=(
    assemblies
    --library viral
    --summary-url "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/assembly_summary.txt"
    --local-dir "$metadata_dir/all"
    --metadata-dir "$metadata_dir"
    --min-count 1000
)

[[ "${REQUIRE_COMPLETE_COLLECTION:-0}" == "1" ]] && args+=(--require-complete)
[[ "${FULL_GZIP_CHECK:-0}" == "1" ]] && args+=(--full-gzip-check)

if [[ "${BUILD_COMBINED_FASTA:-0}" == "1" ]]; then
    args+=(--combined-fasta "$metadata_dir/all_viral_genomes.fna")
fi

[[ "${FORCE_COMBINED_FASTA:-0}" == "1" ]] && args+=(--force-combined)

exec python3 "$helper" "${args[@]}"
