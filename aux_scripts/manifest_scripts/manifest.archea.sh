#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Robust RefSeq Archaea downloader for Kraken2
###############################################################################

offline_files_folder="/media/me/4TB_BACKUP_LBN/Compressed/MTD"

LIBRARY="archaea"

new_download_dir="$offline_files_folder/Kraken2DB_micro/library/$LIBRARY/all"
home_download_dir="$HOME/MTD/kraken2DB_micro/library/$LIBRARY/all"

mkdir -p "$new_download_dir"
mkdir -p "$home_download_dir"
mkdir -p "$offline_files_folder/Kraken2DB_micro/library/$LIBRARY"

assembly_summary_file="$offline_files_folder/Kraken2DB_micro/library/$LIBRARY/assembly_summary_${LIBRARY}.txt"
manifest_list="$offline_files_folder/Kraken2DB_micro/library/$LIBRARY/manifest_${LIBRARY}.list.txt"
failed_downloads="$offline_files_folder/failed_downloads_${LIBRARY}.txt"
corrupted_list="$offline_files_folder/corrupted_${LIBRARY}.txt"
to_download_list="$offline_files_folder/to_download_${LIBRARY}.txt"
obsolete_list="$offline_files_folder/obsolete_local_${LIBRARY}.txt"

rm -f "$assembly_summary_file" "$manifest_list" "$failed_downloads" \
      "$corrupted_list" "$to_download_list" "$obsolete_list"

progress() {
    local current="$1"
    local total="$2"
    local msg="$3"
    msg="${msg:0:100}"
    printf "\r\033[2K[%d/%d] %s" "$current" "$total" "$msg"
}

###############################################################################
# 1. Download assembly_summary.txt and generate manifest
###############################################################################

echo "STEP 1: Downloading latest assembly_summary.txt for $LIBRARY..."

curl -4 -L --retry 20 --retry-delay 10 --retry-all-errors \
  --connect-timeout 30 -fsSL \
  -o "$assembly_summary_file" \
  "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/${LIBRARY}/assembly_summary.txt"

awk -F '\t' '
    /^#/ { next }

    # Column 12 = assembly_level
    # Column 20 = ftp_path
    ($12 == "Complete Genome" || $12 == "Chromosome") && $20 != "na" {
        ftp_path = $20

        sub(/^ftp:/, "https:", ftp_path)
        sub(/ftp\.ncbi\.nlm\.nih\.gov/, "ftp.ncbi.nlm.nih.gov", ftp_path)
        gsub(/\/+$/, "", ftp_path)

        n = split(ftp_path, a, "/")
        asm = a[n]

        sub(/_genomic\.fna\.gz$/, "", asm)
        sub(/\.fna\.gz$/, "", asm)
        sub(/_genomic$/, "", asm)

        if (asm != "") {
            print ftp_path "/" asm "_genomic.fna.gz"
        }
    }
' "$assembly_summary_file" > "$manifest_list"

echo "✅ $(wc -l < "$manifest_list") archaeal genomes listed on NCBI servers."

###############################################################################
# 2. Sync local folder against NCBI manifest
###############################################################################

echo
echo "STEP 2: Syncing local folder against NCBI manifest..."

mapfile -t remote_urls < "$manifest_list"
total_remote=${#remote_urls[@]}

declare -A remote_map

for url in "${remote_urls[@]}"; do
    fname="$(basename "$url")"
    remote_map["$fname"]="$url"
done

shopt -s nullglob
local_paths=("$new_download_dir"/*.gz)
shopt -u nullglob

: > "$obsolete_list"
: > "$to_download_list"
: > "$failed_downloads"

###############################################################################
# 2A. Remove obsolete local files
###############################################################################

echo
echo "STEP 2A: Checking LOCAL files for obsolete files..."

obsolete_count=0
local_total=${#local_paths[@]}
local_checked=0

if (( local_total == 0 )); then
    echo "No local .gz files found yet."
else
    for path in "${local_paths[@]}"; do
        ((local_checked+=1))
        fname="$(basename "$path")"

        progress "$local_checked" "$local_total" "[LOCAL] checking: $fname"

        if [[ -z "${remote_map[$fname]:-}" ]]; then
            echo
            echo "❌ Removing obsolete local file: $fname"
            echo "$fname" >> "$obsolete_list"
            rm -f "$path"
            rm -f "$home_download_dir/$fname"
            ((obsolete_count+=1))
        fi
    done
    echo
fi

###############################################################################
# Refresh local map after removals
###############################################################################

shopt -s nullglob
local_paths=("$new_download_dir"/*.gz)
shopt -u nullglob

declare -A local_map

for path in "${local_paths[@]}"; do
    fname="$(basename "$path")"
    local_map["$fname"]="$path"
done

###############################################################################
# Function to get remote file size
###############################################################################

get_remote_size() {
    local url="$1"

    curl -4 -L -fsSI \
      --retry 5 \
      --retry-delay 5 \
      --retry-all-errors \
      --connect-timeout 30 \
      "$url" 2>/dev/null \
      | awk 'BEGIN{IGNORECASE=1} /^Content-Length:/ {gsub("\r","",$2); print $2; exit}'
}

###############################################################################
# 2B. Detect missing and changed files
###############################################################################

echo
echo "STEP 2B: Checking REMOTE files against LOCAL files..."

missing_count=0
changed_count=0
size_check_failed_count=0
checked_count=0

for fname in "${!remote_map[@]}"; do
    ((checked_count+=1))
    url="${remote_map[$fname]}"

    progress "$checked_count" "$total_remote" "[REMOTE] checking: $fname"

    if [[ -z "${local_map[$fname]:-}" ]]; then
        echo "$url" >> "$to_download_list"
        ((missing_count+=1))
        continue
    fi

    local_size="$(stat -c%s "${local_map[$fname]}" 2>/dev/null || echo 0)"
    remote_size="$(get_remote_size "$url" || true)"

    if [[ -z "$remote_size" ]]; then
        ((size_check_failed_count+=1))
        continue
    fi

    if [[ "$remote_size" != "$local_size" ]]; then
        echo
        echo "Remote file changed, re-downloading: $fname"
        rm -f "${local_map[$fname]}"
        rm -f "$home_download_dir/$fname"
        echo "$url" >> "$to_download_list"
        ((changed_count+=1))
    fi
done

echo
echo "Finished checking remote/local files."

available_local=$(find "$new_download_dir" -maxdepth 1 -type f -name "*.gz" | wc -l)
to_download_count=$(wc -l < "$to_download_list" 2>/dev/null || echo 0)

echo
echo "Total remote             : $total_remote"
echo "Already local            : $available_local"
echo "Obsolete removed         : $obsolete_count"
echo "Missing files            : $missing_count"
echo "Changed files            : $changed_count"
echo "Remote size check failed : $size_check_failed_count"
echo "To download now          : $to_download_count"

###############################################################################
# 3. Download missing/changed files
###############################################################################

echo
echo "STEP 3: Downloading missing/changed files..."

cd "$new_download_dir"

download_one() {
    local url="$1"
    local file
    file="$(basename "$url")"

    for attempt in {1..5}; do
        aria2c --disable-ipv6=true \
               --continue=true \
               --auto-file-renaming=false \
               --max-tries=0 \
               --retry-wait=20 \
               --timeout=60 \
               --connect-timeout=30 \
               -x8 -s8 \
               -o "$file" \
               "$url" > /dev/null 2>&1 || true

        if [[ -f "$file" ]] && gzip -t "$file" 2>/dev/null; then
            [[ -f "$home_download_dir/$file" ]] || cp -p "$file" "$home_download_dir/"
            return 0
        else
            echo "❌ Attempt $attempt failed for $file"
            rm -f "$file"
            sleep 5
        fi
    done

    echo "❌ $file" >> "$failed_downloads"
}
export -f download_one
export failed_downloads
export home_download_dir

if (( to_download_count == 0 )); then
    echo "✅ Local folder is already synchronized with NCBI."
else
    echo "Downloading $to_download_count archaeal genome(s)..."
    mapfile -t download_urls < "$to_download_list"

    if command -v parallel >/dev/null 2>&1 && command -v aria2c >/dev/null 2>&1; then
        printf "%s\n" "${download_urls[@]}" | parallel --bar -j 4 download_one
    else
        echo "⚠️ aria2c/parallel not found – falling back to serial wget."

        download_checked=0

        for url in "${download_urls[@]}"; do
            ((download_checked+=1))
            file="$(basename "$url")"

            progress "$download_checked" "$to_download_count" "[DOWNLOAD] downloading: $file"

            wget -4 -q -c -O "$file" "$url" || {
                echo
                echo "❌ Failed: $file"
                echo "❌ Failed: $file" >> "$failed_downloads"
                rm -f "$file"
                continue
            }

            if gzip -t "$file" 2>/dev/null; then
                [[ -f "$home_download_dir/$file" ]] || cp -p "$file" "$home_download_dir/"
            else
                echo
                echo "❌ Corrupted after download: $file"
                echo "$file" >> "$corrupted_list"
                rm -f "$file"
            fi
        done

        echo
    fi

    echo "✅ Download stage completed."
fi

###############################################################################
# 4. Verify integrity and redownload corrupted files
###############################################################################

echo
echo "STEP 4: Verifying integrity of .gz files..."

: > "$corrupted_list"

shopt -s nullglob
gz_files=("$new_download_dir"/*.gz)
shopt -u nullglob

if (( ${#gz_files[@]} > 0 )); then
    total_gz=${#gz_files[@]}
    checked_gz=0

    for gz in "${gz_files[@]}"; do
        ((checked_gz+=1))
        fname="$(basename "$gz")"

        progress "$checked_gz" "$total_gz" "[GZIP TEST] testing: $fname"

        gzip -t "$gz" 2>/dev/null || echo "$fname" >> "$corrupted_list"
    done

    echo
else
    echo "⚠️ No .gz files found yet in $new_download_dir"
fi

corrupted_count=$(wc -l < "$corrupted_list")

if (( corrupted_count > 0 )); then
    echo "⚠️ $corrupted_count corrupted file(s) found. Re-downloading with wget, max 5 attempts..."

    corrupted_checked=0

    while IFS= read -r fname; do
        [[ -z "$fname" ]] && continue

        ((corrupted_checked+=1))
        url="${remote_map[$fname]:-}"

        if [[ -z "$url" ]]; then
            echo "❌ URL not found in manifest for $fname" >> "$failed_downloads"
            rm -f "$new_download_dir/$fname"
            rm -f "$home_download_dir/$fname"
            continue
        fi

        for attempt in {1..5}; do
            progress "$corrupted_checked" "$corrupted_count" "[REDOWNLOAD] $fname attempt $attempt"

            wget -4 -q -O "$new_download_dir/$fname" "$url" || true

            if gzip -t "$new_download_dir/$fname" 2>/dev/null; then
                echo
                echo "✅ Integrity OK after attempt $attempt: $fname"
                [[ -f "$home_download_dir/$fname" ]] || cp -p "$new_download_dir/$fname" "$home_download_dir/"
                break
            elif [[ $attempt -eq 5 ]]; then
                echo
                echo "❌ Still corrupted after 5 attempts: $fname"
                echo "❌ Failed after 5 attempts: $fname" >> "$failed_downloads"
                rm -f "$new_download_dir/$fname"
                rm -f "$home_download_dir/$fname"
            else
                rm -f "$new_download_dir/$fname"
                sleep 5
            fi
        done
    done < "$corrupted_list"
else
    echo "✅ All .gz files passed integrity check."
fi

rm -f "$corrupted_list"

###############################################################################
# 5. Final report
###############################################################################

echo
echo "STEP 5: Final report"
echo "Remote genomes listed    : $total_remote"
echo "Obsolete removed         : $obsolete_count"
echo "Missing detected         : $missing_count"
echo "Changed detected         : $changed_count"
echo "Remote size check failed : $size_check_failed_count"
echo "Downloaded this run      : $to_download_count"
echo "Local folder             : $new_download_dir"
echo "Home mirror folder       : $home_download_dir"

if [[ -s "$failed_downloads" ]]; then
    echo
    echo "⚠️ The following genomes could not be retrieved:"
    cat "$failed_downloads"
else
    echo
    echo "✅ All archaeal genomes synchronized and verified successfully!"
fi
