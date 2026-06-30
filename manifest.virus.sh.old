#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# 0. Diretórios e nomes de arquivos
###############################################################################
offline_files_folder="/media/me/18TB_BACKUP_LBN/lbn_workspace/RNA-Seq-LBN/viral-rna-seq/MTD/Compressed/MTD"

new_download_dir="$offline_files_folder/Kraken2DB_micro/library/viral/all"

mkdir -p "$new_download_dir"
mkdir -p "$offline_files_folder/Kraken2DB_micro/library/viral"

assembly_summary_file="$offline_files_folder/Kraken2DB_micro/library/viral/assembly_summary_viral.txt"
manifest_list="$offline_files_folder/Kraken2DB_micro/library/viral/manifest_viral.list.txt"
failed_downloads="$offline_files_folder/failed_downloads.txt"
corrupted_list="$offline_files_folder/corrupted_viral.txt"
to_download_list="$offline_files_folder/to_download_viral.txt"
obsolete_list="$offline_files_folder/obsolete_local_viral.txt"

rm -f "$assembly_summary_file" "$manifest_list" "$failed_downloads" \
      "$corrupted_list" "$to_download_list" "$obsolete_list"

progress() {
    local current="$1"
    local total="$2"
    local msg="$3"

    # corta mensagem muito longa para não quebrar a linha
    msg="${msg:0:100}"

    # \033[2K limpa a linha inteira antes de reescrever
    printf "\r\033[2K[%d/%d] %s" "$current" "$total" "$msg"
}

###############################################################################
# 1. Baixar assembly_summary.txt e gerar manifest
###############################################################################
echo "STEP 1: Downloading latest assembly_summary.txt viral..."

curl -4 --retry 5 --retry-delay 2 --connect-timeout 20 -fsSL \
  -o "$assembly_summary_file" \
  "https://ftp.ncbi.nih.gov/genomes/refseq/viral/assembly_summary.txt"

awk -F '\t' '
    /^#/ { next }
    $20 != "na" {
        ftp_path = $20
        sub(/^ftp:/, "https:", ftp_path)
        gsub(/\/+$/, "", ftp_path)

        n = split(ftp_path, a, "/")
        asm = a[n]

        # Segurança extra:
        # remove sufixos caso algum nome venha estranho
        sub(/_genomic\.fna\.gz$/, "", asm)
        sub(/\.fna\.gz$/, "", asm)
        sub(/_genomic$/, "", asm)

        if (asm != "") {
            print ftp_path "/" asm "_genomic.fna.gz"
        }
    }
' "$assembly_summary_file" > "$manifest_list"

echo "✅ $(wc -l < "$manifest_list") viral genomes listed on NCBI servers."

###############################################################################
# 2. Sincronizar pasta local com o NCBI
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
# 2A. Remover arquivos locais obsoletos
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
            ((obsolete_count+=1))
        fi
    done
    echo
fi

###############################################################################
# Atualizar lista local após remoções
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
# Função para obter tamanho remoto
###############################################################################
get_remote_size() {
    local url="$1"

    curl -4 -fsSI --connect-timeout 20 "$url" \
      | awk 'BEGIN{IGNORECASE=1} /^Content-Length:/ {gsub("\r","",$2); print $2; exit}'
}

###############################################################################
# 2B. Detectar faltantes e alterados
###############################################################################
echo
echo "STEP 2B: Checking REMOTE files against LOCAL files..."

missing_count=0
changed_count=0
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
    remote_size="$(get_remote_size "$url" || echo "")"

    if [[ -n "$remote_size" && "$remote_size" != "$local_size" ]]; then
        echo
        echo "Remote file changed, re-downloading: $fname"
        rm -f "${local_map[$fname]}"
        echo "$url" >> "$to_download_list"
        ((changed_count+=1))
    fi
done

echo
echo "Finished checking remote/local files."

available_local=$(find "$new_download_dir" -maxdepth 1 -type f -name "*.gz" | wc -l)
to_download_count=$(wc -l < "$to_download_list" 2>/dev/null || echo 0)

echo
echo "Total remote      : $total_remote"
echo "Already local     : $available_local"
echo "Obsolete removed  : $obsolete_count"
echo "Missing files     : $missing_count"
echo "Changed files     : $changed_count"
echo "To download now   : $to_download_count"

###############################################################################
# 3. Baixar arquivos faltantes/alterados
###############################################################################
echo
echo "STEP 3: Downloading missing/changed files..."

cd "$new_download_dir"

download_one() {
    local url="$1"
    local file
    file="$(basename "$url")"

    for attempt in {1..3}; do
        aria2c --disable-ipv6=true \
               --continue=true \
               --auto-file-renaming=false \
               -x16 -s16 \
               -o "$file" \
               "$url" > /dev/null 2>&1

        if [[ -f "$file" ]] && gzip -t "$file" 2>/dev/null; then
            return 0
        else
            echo "❌ Attempt $attempt failed for $file"
            rm -f "$file"
            sleep 1
        fi
    done

    echo "❌ $file" >> "$failed_downloads"
}
export -f download_one
export failed_downloads

if (( to_download_count == 0 )); then
    echo "✅ Local folder is already synchronized with NCBI."
else
    echo "Downloading $to_download_count genome(s)..."
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
            }
        done

        echo
    fi

    echo "✅ Download stage completed."
fi

###############################################################################
# 4. Checar integridade e re-baixar corrompidos
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
    echo "⚠️ $corrupted_count corrupted file(s) found. Re-downloading with wget, max 3 attempts..."

    corrupted_checked=0

    while IFS= read -r fname; do
        [[ -z "$fname" ]] && continue

        ((corrupted_checked+=1))
        url="${remote_map[$fname]:-}"

        if [[ -z "$url" ]]; then
            echo "❌ URL not found in manifest for $fname" >> "$failed_downloads"
            rm -f "$new_download_dir/$fname"
            continue
        fi

        for attempt in {1..3}; do
            progress "$corrupted_checked" "$corrupted_count" "[REDOWNLOAD] $fname attempt $attempt"

            wget -4 -q -O "$new_download_dir/$fname" "$url" || true

            if gzip -t "$new_download_dir/$fname" 2>/dev/null; then
                echo
                echo "✅ Integrity OK after attempt $attempt: $fname"
                break
            elif [[ $attempt -eq 3 ]]; then
                echo
                echo "❌ Still corrupted after 3 attempts: $fname"
                echo "❌ Failed after 3 attempts: $fname" >> "$failed_downloads"
                rm -f "$new_download_dir/$fname"
            else
                rm -f "$new_download_dir/$fname"
            fi
        done
    done < "$corrupted_list"
else
    echo "✅ All .gz files passed integrity check."
fi

rm -f "$corrupted_list"

###############################################################################
# 5. Descompactar e concatenar em um único FASTA
###############################################################################
echo
echo "STEP 5: Building combined FASTA..."

combined_fasta="$new_download_dir/all_viral_genomes.fna"
final_fasta="$offline_files_folder/Kraken2DB_micro/library/viral/all_viral_genomes.fna"

echo "Output FASTA: $final_fasta"

: > "$combined_fasta"

shopt -s nullglob
gz_files=("$new_download_dir"/*.gz)
shopt -u nullglob

if (( ${#gz_files[@]} == 0 )); then
    echo "❌ No .gz files available to concatenate."
    exit 1
fi

total_gz=${#gz_files[@]}
cat_count=0

for gz in "${gz_files[@]}"; do
    ((cat_count+=1))
    fname="$(basename "$gz")"

    progress "$cat_count" "$total_gz" "[FASTA] adding: $fname"

    zcat "$gz" >> "$combined_fasta"
done

echo

mv "$combined_fasta" "$final_fasta"

echo "✅ Combined FASTA created: $final_fasta"

###############################################################################
# 6. Relatório final
###############################################################################
echo
echo "STEP 6: Final report"
echo "Remote genomes listed : $total_remote"
echo "Obsolete removed      : $obsolete_count"
echo "Missing detected      : $missing_count"
echo "Changed detected      : $changed_count"
echo "Downloaded this run   : $to_download_count"

if [[ -s "$failed_downloads" ]]; then
    echo
    echo "⚠️ The following genomes could not be retrieved:"
    cat "$failed_downloads"
else
    echo
    echo "✅ All viral genomes synchronized, verified and concatenated successfully!"
fi
