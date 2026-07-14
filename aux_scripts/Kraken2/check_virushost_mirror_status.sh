#!/usr/bin/env bash
set -euo pipefail

repository="${VIRUSHOST_MIRROR_REPOSITORY:-patrick-douglas/MTD}"
install_script="${INSTALL_SCRIPT:-Install.sh}"
official_url="${VIRUSHOST_OFFICIAL_DBREL_URL:-https://www.genome.jp/ftp/db/virushostdb/dbrel.txt}"

if [[ ! -s "$install_script" ]]; then
    echo "[FAIL] Install script not found: $install_script" >&2
    exit 1
fi

tag="$(
    sed -nE \
        's/^VIRUSHOST_MIRROR_TAG="([^"]+)"/\1/p' \
        "$install_script"
)"

if [[ -z "$tag" ]]; then
    echo "[FAIL] Could not read VIRUSHOST_MIRROR_TAG from $install_script" >&2
    exit 1
fi

mirror_url="https://github.com/${repository}/releases/download/${tag}/dbrel.txt"

work_dir="$(mktemp -d)"
trap 'rm -rf -- "$work_dir"' EXIT

official_file="$work_dir/official.dbrel.txt"
mirror_file="$work_dir/mirror.dbrel.txt"

download_file() {
    local label="$1"
    local url="$2"
    local destination="$3"

    echo "[INFO] Downloading $label..."

    curl \
        --fail \
        --silent \
        --show-error \
        --location \
        --retry 5 \
        --retry-delay 5 \
        --connect-timeout 30 \
        --output "$destination" \
        "$url"

    if [[ ! -s "$destination" ]]; then
        echo "[FAIL] Empty $label file downloaded." >&2
        exit 1
    fi

    sed -i 's/\r$//' "$destination"
}

extract_refseq_release() {
    sed -nE \
        's/.*RefSeq Release ([0-9]+).*/\1/p' \
        "$1" |
    head -n 1
}

extract_genbank_release() {
    sed -nE \
        's/.*GenBank Release ([0-9]+(\.[0-9]+)?).*/\1/p' \
        "$1" |
    head -n 1
}

version_is_greater() {
    local first="$1"
    local second="$2"
    local greatest

    [[ "$first" != "$second" ]] || return 1

    greatest="$(
        printf '%s\n%s\n' "$first" "$second" |
        sort -V |
        tail -n 1
    )"

    [[ "$greatest" == "$first" ]]
}

download_file \
    "official dbrel.txt" \
    "$official_url" \
    "$official_file"

download_file \
    "mirror dbrel.txt" \
    "$mirror_url" \
    "$mirror_file"

official_refseq="$(extract_refseq_release "$official_file")"
official_genbank="$(extract_genbank_release "$official_file")"
mirror_refseq="$(extract_refseq_release "$mirror_file")"
mirror_genbank="$(extract_genbank_release "$mirror_file")"

for value_name in \
    official_refseq \
    official_genbank \
    mirror_refseq \
    mirror_genbank
do
    value="${!value_name}"

    if [[ -z "$value" ]]; then
        echo "[FAIL] Could not parse release value: $value_name" >&2
        exit 1
    fi
done

echo
echo "Official Virus-Host DB:"
cat "$official_file"

echo
echo "Pinned mirror:"
cat "$mirror_file"

echo
echo "Mirror tag:"
echo "  $tag"
echo

if cmp -s "$official_file" "$mirror_file"; then
    echo "[PASS] MIRROR UP TO DATE"
    exit 0
fi

if [[ "$official_refseq" == "$mirror_refseq" &&
      "$official_genbank" == "$mirror_genbank" ]]
then
    echo "[WARNING] Release numbers match, but dbrel.txt metadata differs."
    exit 3
fi

official_newer=0
mirror_newer=0

if version_is_greater "$official_refseq" "$mirror_refseq"; then
    official_newer=1
elif version_is_greater "$mirror_refseq" "$official_refseq"; then
    mirror_newer=1
fi

if version_is_greater "$official_genbank" "$mirror_genbank"; then
    official_newer=1
elif version_is_greater "$mirror_genbank" "$official_genbank"; then
    mirror_newer=1
fi

if (( official_newer == 1 && mirror_newer == 0 )); then
    echo "[OUTDATED] MIRROR IS BEHIND THE OFFICIAL RELEASE"
    echo "[OUTDATED] Official: RefSeq $official_refseq; GenBank $official_genbank"
    echo "[OUTDATED] Mirror:   RefSeq $mirror_refseq; GenBank $mirror_genbank"
    exit 2
fi

if (( mirror_newer == 1 && official_newer == 0 )); then
    echo "[WARNING] Mirror appears newer than the current official dbrel.txt."
    exit 3
fi

echo "[WARNING] RefSeq and GenBank release comparisons are inconsistent."
exit 3
