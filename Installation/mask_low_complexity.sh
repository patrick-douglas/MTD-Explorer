#!/usr/bin/env bash
# MTD_KRAKEN2_PARALLEL_MASK_V2
#
# Drop-in replacement for Kraken2's legacy mask_low_complexity.sh.
# It preserves the legacy kraken2-build workflow while passing an
# explicit thread count to the multithreaded k2mask bundled with
# Kraken2 2.17.1.

set -euo pipefail

target="${1:-}"

if [[ -z "$target" ]]; then
    echo "[ERROR] mask_low_complexity.sh requires a file or directory." >&2
    exit 1
fi

masker_threads="${MTD_KRAKEN2_MASKER_THREADS:-${KRAKEN2_THREAD_CT:-1}}"

if ! [[ "$masker_threads" =~ ^[1-9][0-9]*$ ]]; then
    echo "[ERROR] Invalid Kraken2 masker thread count: $masker_threads" >&2
    exit 1
fi

protein_db="${KRAKEN2_PROTEIN_DB:-}"
masker="k2mask"

if [[ -n "$protein_db" ]]; then
    masker="segmasker"
fi

if ! command -v "$masker" >/dev/null 2>&1; then
    echo "[ERROR] Unable to find $masker in PATH." >&2
    exit 1
fi

echo "[MTD-KRAKEN2] Low-complexity masker: $masker" >&2

if [[ "$masker" == "k2mask" ]]; then
    echo "[MTD-KRAKEN2] k2mask threads: $masker_threads" >&2
fi

mask_one_file() {
    local file="$1"
    local temporary="${file}.tmp.$$"

    if [[ -e "${file}.masked" ]]; then
        return 0
    fi

    rm -f -- "$temporary"

    if [[ "$masker" == "k2mask" ]]; then
        if ! "$masker"             -in "$file"             -outfmt fasta             -threads "$masker_threads"             -r x             > "$temporary"; then

            rm -f -- "$temporary"
            echo "[ERROR] k2mask failed for: $file" >&2
            return 1
        fi
    else
        if ! "$masker"             -in "$file"             -outfmt fasta |
            sed -e '/^>/!s/[a-z]/x/g'             > "$temporary"; then

            rm -f -- "$temporary"
            echo "[ERROR] segmasker failed for: $file" >&2
            return 1
        fi
    fi

    if [[ ! -s "$temporary" ]]; then
        rm -f -- "$temporary"
        echo "[ERROR] Masked FASTA is empty: $file" >&2
        return 1
    fi

    mv -f -- "$temporary" "$file"
    touch -- "${file}.masked"
}

if [[ -d "$target" ]]; then
    found=0

    while IFS= read -r -d '' file; do
        found=1
        mask_one_file "$file"
    done < <(
        find "$target"             -type f             \( -name '*.fna' -o -name '*.faa' \)             -print0
    )

    if [[ "$found" -eq 0 ]]; then
        echo "[WARNING] No .fna or .faa files found under: $target" >&2
    fi
elif [[ -f "$target" ]]; then
    mask_one_file "$target"
else
    echo "[ERROR] Target must be a directory or regular file: $target" >&2
    exit 1
fi
