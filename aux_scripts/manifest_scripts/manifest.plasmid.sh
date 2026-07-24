#!/usr/bin/env bash
set -euo pipefail

LOCAL_DIR="/media/me/18TB_BACKUP_LBN/drive.ifpa/LBN_RNA-Seq/Metatranscriptomics/MTD/MTD_Offline_Install_files/Kraken2DB_micro/library/plasmid"
LOCAL_DIR="${MTD_KRAKEN2_PLASMID_CACHE:-$LOCAL_DIR}"

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
helper="${MTD_MANIFEST_HELPER:-$script_dir/sync_ncbi_cache.py}"

[[ -f "$helper" ]] || {
    echo "[FAIL] NCBI cache synchronization helper not found: $helper" >&2
    exit 1
}

args=(
    release
    --library plasmid
    --base-url "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/plasmid"
    --release-number-url "https://ftp.ncbi.nlm.nih.gov/refseq/release/RELEASE_NUMBER"
    --pattern 'plasmid\.[0-9]+\.[0-9]+\.genomic\.fna\.gz'
    --local-dir "$LOCAL_DIR"
    --metadata-dir "$LOCAL_DIR"
    --min-count 10
    --require-complete
)

[[ "${FULL_GZIP_CHECK:-0}" == "1" ]] && args+=(--full-gzip-check)

exec python3 "$helper" "${args[@]}"
