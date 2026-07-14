#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# update_virushost_mirror.sh
#
# Maintainer-only workflow for mirroring a validated Virus-Host DB snapshot to
# a pinned GitHub Release. This script is not called by Install.sh.
#
# Optional overrides:
#   GITHUB_REPOSITORY=owner/repo
#   VIRUSHOST_BASE_URL=https://www.genome.jp/ftp/db/virushostdb
#   VIRUSHOST_MIRROR_WORKDIR=/path/to/workdir
#   VIRUSHOST_MIRROR_TAG=custom-tag
###############################################################################

REPOSITORY="${GITHUB_REPOSITORY:-patrick-douglas/MTD}"
SOURCE_BASE_URL="${VIRUSHOST_BASE_URL:-https://www.genome.jp/ftp/db/virushostdb}"
LOCAL_SOURCE_DIR="${VIRUSHOST_LOCAL_SOURCE_DIR:-}"

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
repo_dir="$(cd -- "$script_dir/../.." && pwd -P)"
work_root="${VIRUSHOST_MIRROR_WORKDIR:-$repo_dir/build/virushost_mirror}"

required_assets=(
    "virushostdb.genomic.fna.gz"
    "non-segmented_virus_list.tsv"
    "segmented_virus_list.tsv"
)

log_info() {
    printf '[INFO] %s\n' "$*"
}

log_ok() {
    printf '[OK] %s\n' "$*"
}

log_error() {
    printf '[ERROR] %s\n' "$*" >&2
}

require_command() {
    local command_name="$1"

    if ! command -v "$command_name" >/dev/null 2>&1; then
        log_error "Required command not found: $command_name"
        exit 127
    fi
}

for command_name in \
    awk \
    cmp \
    curl \
    date \
    gh \
    gzip \
    grep \
    mktemp \
    sed \
    sha256sum \
    stat \
    wc
do
    require_command "$command_name"
done

if ! gh auth status >/dev/null 2>&1; then
    log_error "GitHub CLI is not authenticated."
    log_error "Run: gh auth login"
    exit 1
fi

mkdir -p "$work_root"

stage_dir="$(mktemp -d "$work_root/.incoming.XXXXXX")"
cleanup_stage=1

cleanup() {
    local exit_status=$?

    if (( cleanup_stage == 1 )); then
        if (( exit_status == 0 )); then
            rm -rf -- "$stage_dir"
        else
            printf '%s\n' \
                "[WARNING] Mirror creation failed after downloading files." \
                "[WARNING] The staging directory was preserved for recovery:" \
                "[WARNING]   $stage_dir" \
                >&2
        fi
    fi
}
trap cleanup EXIT

download_file() {
    local remote_name="$1"
    local local_name="${2:-$remote_name}"
    local destination="$stage_dir/$local_name"
    local partial="${destination}.part"

    if [[ -n "$LOCAL_SOURCE_DIR" &&
          -s "$LOCAL_SOURCE_DIR/$remote_name" ]]
    then
        log_info "Using local validated source:"
        log_info "  $LOCAL_SOURCE_DIR/$remote_name"

        cp -f --             "$LOCAL_SOURCE_DIR/$remote_name"             "$destination"

        return 0
    fi

    log_info "Downloading $remote_name"

    rm -f -- "$partial"

    curl -4 \
        --fail \
        --location \
        --show-error \
        --retry 8 \
        --retry-all-errors \
        --retry-delay 5 \
        --connect-timeout 30 \
        --continue-at - \
        --output "$partial" \
        "$SOURCE_BASE_URL/$remote_name"

    if [[ ! -s "$partial" ]]; then
        log_error "Downloaded file is empty: $remote_name"
        return 1
    fi

    mv -f -- "$partial" "$destination"
}

###############################################################################
# 1. Download a consistent upstream snapshot
###############################################################################

download_file "dbrel.txt" "dbrel.before.txt"

for asset in "${required_assets[@]}"; do
    download_file "$asset"
done

download_file "dbrel.txt" "dbrel.after.txt"

if ! cmp -s \
    "$stage_dir/dbrel.before.txt" \
    "$stage_dir/dbrel.after.txt"
then
    log_error "Virus-Host DB release metadata changed during the download."
    diff -u \
        "$stage_dir/dbrel.before.txt" \
        "$stage_dir/dbrel.after.txt" \
        || true
    exit 1
fi

mv -f \
    "$stage_dir/dbrel.after.txt" \
    "$stage_dir/dbrel.txt"

rm -f "$stage_dir/dbrel.before.txt"

###############################################################################
# 2. Validate the downloaded snapshot
###############################################################################

gzip -t "$stage_dir/virushostdb.genomic.fna.gz"

for table in \
    non-segmented_virus_list.tsv \
    segmented_virus_list.tsv
do
    if ! awk -F $'\t' '
        NF >= 3 && $2 ~ /^[0-9]+$/ && $2 > 1 {
            valid_rows++
        }
        END {
            exit(valid_rows < 10)
        }
    ' "$stage_dir/$table"; then
        log_error "Virus-Host DB metadata table failed validation: $table"
        exit 1
    fi
done

fasta_records="$(
    gzip -dc "$stage_dir/virushostdb.genomic.fna.gz" |
    awk '/^>/ { records++ } END { print records + 0 }'
)"

if ! [[ "$fasta_records" =~ ^[0-9]+$ ]] ||
   (( fasta_records < 1000 )); then
    log_error "Unexpected Virus-Host DB FASTA record count: $fasta_records"
    exit 1
fi

non_segmented_rows="$(
    awk -F $'\t' '
        NF >= 3 && $2 ~ /^[0-9]+$/ && $2 > 1 {
            rows++
        }
        END {
            print rows + 0
        }
    ' "$stage_dir/non-segmented_virus_list.tsv"
)"

segmented_rows="$(
    awk -F $'\t' '
        NF >= 3 && $2 ~ /^[0-9]+$/ && $2 > 1 {
            rows++
        }
        END {
            print rows + 0
        }
    ' "$stage_dir/segmented_virus_list.tsv"
)"

###############################################################################
# 3. Build checksums, metadata, and an immutable tag
###############################################################################

(
    cd "$stage_dir"

    sha256sum \
        dbrel.txt \
        virushostdb.genomic.fna.gz \
        non-segmented_virus_list.tsv \
        segmented_virus_list.tsv \
        > SHA256SUMS
)

bundle_sha="$(
    sha256sum "$stage_dir/SHA256SUMS" |
    awk '{print $1}'
)"

refseq_release="$(
    sed -nE \
        's/.*RefSeq Release ([0-9]+).*/\1/p' \
        "$stage_dir/dbrel.txt" |
    head -n 1
)"

genbank_release="$(
    sed -nE \
        's/.*GenBank Release ([0-9.]+).*/\1/p' \
        "$stage_dir/dbrel.txt" |
    head -n 1
)"

: "${refseq_release:=unknown}"
: "${genbank_release:=unknown}"

default_tag="virushostdb-mirror-r${refseq_release}-g${genbank_release}-${bundle_sha:0:12}"
release_tag="${VIRUSHOST_MIRROR_TAG:-$default_tag}"
fetched_at="$(date -u '+%Y-%m-%dT%H:%M:%SZ')"

{
    printf 'field\tvalue\n'
    printf 'repository\t%s\n' "$REPOSITORY"
    printf 'release_tag\t%s\n' "$release_tag"
    printf 'source_base_url\t%s\n' "$SOURCE_BASE_URL"
    printf 'fetched_at_utc\t%s\n' "$fetched_at"
    printf 'refseq_release\t%s\n' "$refseq_release"
    printf 'genbank_release\t%s\n' "$genbank_release"
    printf 'fasta_records\t%s\n' "$fasta_records"
    printf 'non_segmented_rows\t%s\n' "$non_segmented_rows"
    printf 'segmented_rows\t%s\n' "$segmented_rows"
    printf 'bundle_sha256\t%s\n' "$bundle_sha"
    printf 'dbrel\t%s\n' "$(tr '\t\r\n' '   ' < "$stage_dir/dbrel.txt")"
} > "$stage_dir/MIRROR_METADATA.tsv"

cat > "$stage_dir/RELEASE_NOTES.md" <<EOF
Validated mirror of the official Virus-Host DB files used by MTD Explorer.

Upstream source:
$SOURCE_BASE_URL

Upstream release metadata:
$(cat "$stage_dir/dbrel.txt")

Validation:
- FASTA gzip integrity: PASS
- FASTA records: $fasta_records
- Non-segmented metadata rows: $non_segmented_rows
- Segmented metadata rows: $segmented_rows
- Bundle checksum: $bundle_sha

The SHA256SUMS asset contains checksums for all upstream files.
This release is intended to be pinned by Install.sh for reproducible installs.
EOF

###############################################################################
# 4. Create a draft GitHub Release
###############################################################################

if gh release view \
    "$release_tag" \
    --repo "$REPOSITORY" \
    >/dev/null 2>&1
then
    log_error "A GitHub Release already exists for tag:"
    log_error "  $release_tag"
    log_error "No asset was replaced."
    exit 1
fi

log_info "Creating draft GitHub Release: $release_tag"

gh release create \
    "$release_tag" \
    "$stage_dir/virushostdb.genomic.fna.gz" \
    "$stage_dir/non-segmented_virus_list.tsv" \
    "$stage_dir/segmented_virus_list.tsv" \
    "$stage_dir/dbrel.txt" \
    "$stage_dir/SHA256SUMS" \
    "$stage_dir/MIRROR_METADATA.tsv" \
    --repo "$REPOSITORY" \
    --target main \
    --title "Virus-Host DB mirror — RefSeq ${refseq_release}, GenBank ${genbank_release}" \
    --notes-file "$stage_dir/RELEASE_NOTES.md" \
    --draft

release_url="$(
    gh release view \
        "$release_tag" \
        --repo "$REPOSITORY" \
        --json url \
        --jq '.url'
)"

artifact_dir="$work_root/$release_tag"
rm -rf -- "$artifact_dir"
mv -- "$stage_dir" "$artifact_dir"
cleanup_stage=0

echo
log_ok "Validated mirror uploaded as a draft release."
printf 'Tag:       %s\n' "$release_tag"
printf 'Release:   %s\n' "$release_url"
printf 'Artifacts: %s\n' "$artifact_dir"
echo
printf 'Review assets:\n'
printf '  gh release view %q --repo %q --web\n' \
    "$release_tag" \
    "$REPOSITORY"
echo
printf 'Publish after review:\n'
printf 'Publish after review without marking it as latest:\n'
printf '  RELEASE_ID="$(gh api repos/%s/releases --jq '\''.[] | select(.tag_name == "%s") | .id'\'')"\n' \
    "$REPOSITORY" \
    "$release_tag"

printf '  gh api --method PATCH repos/%s/releases/"$RELEASE_ID" -F draft=false -f make_latest=false\n' \
    "$REPOSITORY"
echo
printf 'Pinned Install.sh base URL:\n'
printf '  https://github.com/%s/releases/download/%s\n' \
    "$REPOSITORY" \
    "$release_tag"
