#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# manifest.bacteria.local_sync.sh
#
# Sincronização otimizada dos genomas bacterianos RefSeq.
#
# Estratégia:
#   1. Faz apenas UM download do assembly_summary.txt.
#   2. Gera localmente a lista esperada de arquivos.
#   3. Compara nomes remotos e locais usando sort/comm.
#   4. NÃO executa curl -I individual para cada genoma.
#   5. Usa cache local de tamanho+mtime para aplicar gzip -t somente a arquivos
#      novos ou modificados desde a última execução.
#
# Para forçar auditoria completa de todos os .gz:
#   FULL_GZIP_CHECK=1 bash manifest.bacteria.local_sync.sh
#
# Para ajustar o número de testes gzip paralelos:
#   GZIP_CHECK_JOBS=4 bash manifest.bacteria.local_sync.sh
###############################################################################

###############################################################################
# 0. Diretórios e configurações
###############################################################################
offline_files_folder="/media/me/4TB_BACKUP_LBN/Compressed/MTD"
metadata_dir="$offline_files_folder/Kraken2DB_micro/library/bacteria"
new_download_dir="$metadata_dir/all"
home_download_dir="$HOME/MTD/kraken2DB_micro/library/bacteria/all"

FULL_GZIP_CHECK="${FULL_GZIP_CHECK:-0}"
GZIP_CHECK_JOBS="${GZIP_CHECK_JOBS:-4}"
DOWNLOAD_JOBS="${DOWNLOAD_JOBS:-4}"
ARIA_CONNECTIONS="${ARIA_CONNECTIONS:-16}"
DOWNLOAD_ATTEMPTS="${DOWNLOAD_ATTEMPTS:-3}"
FAILED_RETRY_ROUNDS="${FAILED_RETRY_ROUNDS:-2}"
RETRY_DELAY_SECONDS="${RETRY_DELAY_SECONDS:-5}"

if [[ "$FULL_GZIP_CHECK" != "0" && "$FULL_GZIP_CHECK" != "1" ]]; then
    echo "[FAIL] FULL_GZIP_CHECK must be 0 or 1." >&2
    exit 1
fi

if ! [[ "$GZIP_CHECK_JOBS" =~ ^[1-9][0-9]*$ ]]; then
    echo "[FAIL] GZIP_CHECK_JOBS must be a positive integer." >&2
    exit 1
fi

if ! [[ "$DOWNLOAD_JOBS" =~ ^[1-9][0-9]*$ ]]; then
    echo "[FAIL] DOWNLOAD_JOBS must be a positive integer." >&2
    exit 1
fi

if ! [[ "$ARIA_CONNECTIONS" =~ ^[1-9][0-9]*$ ]]; then
    echo "[FAIL] ARIA_CONNECTIONS must be a positive integer." >&2
    exit 1
fi

if ! [[ "$DOWNLOAD_ATTEMPTS" =~ ^[1-9][0-9]*$ ]]; then
    echo "[FAIL] DOWNLOAD_ATTEMPTS must be a positive integer." >&2
    exit 1
fi

if ! [[ "$FAILED_RETRY_ROUNDS" =~ ^[1-9][0-9]*$ ]]; then
    echo "[FAIL] FAILED_RETRY_ROUNDS must be a positive integer." >&2
    exit 1
fi

if ! [[ "$RETRY_DELAY_SECONDS" =~ ^[0-9]+$ ]]; then
    echo "[FAIL] RETRY_DELAY_SECONDS must be zero or a positive integer." >&2
    exit 1
fi

mkdir -p "$new_download_dir" "$home_download_dir" "$metadata_dir"

assembly_summary_file="$metadata_dir/assembly_summary_bacteria.txt"
manifest_list="$metadata_dir/manifest_bacteria.list.txt"
manifest_names="$metadata_dir/manifest_bacteria.names.txt"
integrity_cache="$metadata_dir/integrity_bacteria.stat.tsv"
failed_downloads="$offline_files_folder/failed_downloads_bacteria.txt"
corrupted_list="$offline_files_folder/corrupted_bacteria.txt"
to_download_list="$offline_files_folder/to_download_bacteria.txt"
obsolete_list="$offline_files_folder/obsolete_local_bacteria.txt"

work_dir="$(mktemp -d "$metadata_dir/.manifest_bacteria.tmp.XXXXXX")"
trap 'rm -rf "$work_dir"' EXIT

: > "$failed_downloads"
: > "$corrupted_list"
: > "$to_download_list"
: > "$obsolete_list"

progress() {
    local current="$1"
    local total="$2"
    local msg="$3"

    msg="${msg:0:100}"
    printf "\r\033[2K[%d/%d] %s" "$current" "$total" "$msg"
}

build_local_names() {
    find "$new_download_dir" \
        -maxdepth 1 \
        -type f \
        -name '*.gz' \
        -printf '%f\n' \
        | LC_ALL=C sort -u
}

build_local_state() {
    find "$new_download_dir" \
        -maxdepth 1 \
        -type f \
        -name '*.gz' \
        -printf '%f\t%s\t%T@\n' \
        | LC_ALL=C sort -t $'\t' -k1,1
}

###############################################################################
# 1. Baixar somente o catálogo do NCBI e gerar manifest local
###############################################################################
echo "STEP 1: Downloading latest bacterial assembly_summary.txt..."

tmp_summary="$work_dir/assembly_summary_bacteria.txt"
tmp_manifest="$work_dir/manifest_bacteria.list.txt"
tmp_remote_index="$work_dir/remote_name_url.tsv"
tmp_remote_names="$work_dir/remote_names.txt"

curl -4 \
    --retry 5 \
    --retry-delay 2 \
    --connect-timeout 20 \
    -fsSL \
    -o "$tmp_summary" \
    "https://ftp.ncbi.nih.gov/genomes/refseq/bacteria/assembly_summary.txt"

if [[ ! -s "$tmp_summary" ]]; then
    echo "[FAIL] Downloaded assembly_summary.txt is empty." >&2
    exit 1
fi

awk -F '\t' '
    /^#/ { next }
    ($12 == "Complete Genome" || $12 == "Chromosome") && $20 != "na" {
        ftp_path = $20
        sub(/^ftp:/, "https:", ftp_path)
        sub(/ftp\.ncbi\.nlm\.nih\.gov/, "ftp.ncbi.nih.gov", ftp_path)
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
' "$tmp_summary" | LC_ALL=C sort -u > "$tmp_manifest"

if [[ ! -s "$tmp_manifest" ]]; then
    echo "[FAIL] No genome URLs were generated from assembly_summary.txt." >&2
    exit 1
fi

# Índice local: filename<TAB>URL
awk -F/ 'NF { print $NF "\t" $0 }' "$tmp_manifest" \
    | LC_ALL=C sort -t $'\t' -k1,1 \
    > "$tmp_remote_index"

cut -f1 "$tmp_remote_index" > "$tmp_remote_names"

# Atualização atômica: arquivos antigos só são substituídos após sucesso.
mv -f "$tmp_summary" "$assembly_summary_file"
mv -f "$tmp_manifest" "$manifest_list"
cp -f "$tmp_remote_names" "$manifest_names"

total_remote="$(wc -l < "$manifest_names")"
echo "[PASS] $total_remote bacterial genomes listed in the NCBI catalog."
echo "[INFO] Remote per-file HEAD requests: 0"

###############################################################################
# 2. Comparação totalmente local: nomes esperados versus nomes presentes
###############################################################################
echo
echo "STEP 2: Comparing the NCBI catalog with local filenames..."

local_names_before="$work_dir/local_names.before.txt"
local_names_after="$work_dir/local_names.after.txt"
missing_names="$work_dir/missing_names.txt"

build_local_names > "$local_names_before"

###############################################################################
# 2A. Arquivos locais obsoletos
###############################################################################
echo
echo "STEP 2A: Detecting obsolete local files..."

# Presentes localmente, mas ausentes do catálogo atual do NCBI.
comm -23 "$local_names_before" "$manifest_names" > "$obsolete_list"
obsolete_count="$(wc -l < "$obsolete_list")"

if (( obsolete_count == 0 )); then
    echo "[PASS] No obsolete local files detected."
else
    removed=0
    while IFS= read -r fname; do
        [[ -z "$fname" ]] && continue
        ((removed+=1))
        progress "$removed" "$obsolete_count" "[REMOVE] $fname"
        rm -f -- "$new_download_dir/$fname"
        rm -f -- "$home_download_dir/$fname"
        rm -f -- "$new_download_dir/$fname.aria2"
    done < "$obsolete_list"
    echo
    echo "[PASS] Removed $obsolete_count obsolete local file(s)."
fi

# Reconstruir lista local após remoções.
build_local_names > "$local_names_after"

###############################################################################
# 2B. Arquivos ausentes
###############################################################################
echo
echo "STEP 2B: Detecting missing genomes locally..."

# Presentes no catálogo remoto, mas ausentes localmente.
comm -13 "$local_names_after" "$manifest_names" > "$missing_names"
missing_count="$(wc -l < "$missing_names")"

# Traduzir filename -> URL sem fazer qualquer consulta individual ao servidor.
awk -F '\t' '
    NR == FNR {
        wanted[$1] = 1
        next
    }
    ($1 in wanted) {
        print $2
    }
' "$missing_names" "$tmp_remote_index" > "$to_download_list"

to_download_count="$(wc -l < "$to_download_list")"
available_local="$(wc -l < "$local_names_after")"

echo "[PASS] Local comparison completed."
echo
echo "Remote catalog entries : $total_remote"
echo "Already local          : $available_local"
echo "Obsolete removed       : $obsolete_count"
echo "Missing files          : $missing_count"
echo "To download now        : $to_download_count"
echo "Remote HEAD requests   : 0"

retry_failed_downloads() {
    local phase_label="$1"
    local round current_failed retry_urls unmapped_names
    local total_failed remaining_count checked fname url attempt success

    [[ -s "$failed_downloads" ]] || return 0
    LC_ALL=C sort -u "$failed_downloads" -o "$failed_downloads"

    declare -A retry_url_map=()
    while IFS=$'\t' read -r fname url; do
        [[ -n "$fname" && -n "$url" ]] || continue
        retry_url_map["$fname"]="$url"
    done < "$tmp_remote_index"

    for ((round=1; round<=FAILED_RETRY_ROUNDS; round++)); do
        [[ -s "$failed_downloads" ]] || break

        current_failed="$work_dir/failed.${phase_label}.round${round}.txt"
        retry_urls="$work_dir/failed.${phase_label}.round${round}.urls.txt"
        unmapped_names="$work_dir/failed.${phase_label}.round${round}.unmapped.txt"

        cp -f "$failed_downloads" "$current_failed"
        : > "$retry_urls"
        : > "$unmapped_names"
        : > "$failed_downloads"

        total_failed="$(wc -l < "$current_failed")"
        echo
        echo "[WARN] Retry round $round/$FAILED_RETRY_ROUNDS for $total_failed failed genome(s) [$phase_label]..."

        while IFS= read -r fname; do
            [[ -z "$fname" ]] && continue
            url="${retry_url_map[$fname]:-}"

            if [[ -n "$url" ]]; then
                printf '%s\n' "$url" >> "$retry_urls"
            else
                echo "[WARN] URL not found in current manifest; skipping: $fname" >&2
                printf '%s\n' "$fname" >> "$unmapped_names"
            fi
        done < "$current_failed"

        if [[ -s "$retry_urls" ]]; then
            if command -v parallel >/dev/null 2>&1 \
                && command -v aria2c >/dev/null 2>&1; then
                parallel \
                    --bar \
                    --halt never \
                    -j "$DOWNLOAD_JOBS" \
                    download_one :::: "$retry_urls" || true
            elif command -v aria2c >/dev/null 2>&1; then
                while IFS= read -r url; do
                    [[ -z "$url" ]] && continue
                    download_one "$url" || true
                done < "$retry_urls"
            else
                checked=0
                while IFS= read -r url; do
                    [[ -z "$url" ]] && continue
                    ((checked+=1))
                    fname="$(basename "$url")"
                    success=0

                    for ((attempt=1; attempt<=DOWNLOAD_ATTEMPTS; attempt++)); do
                        progress "$checked" "$total_failed" \
                            "[RETRY $round/$FAILED_RETRY_ROUNDS] $fname attempt $attempt/$DOWNLOAD_ATTEMPTS"

                        rm -f -- "$new_download_dir/$fname"
                        if wget -4 -q -O "$new_download_dir/$fname" "$url" \
                            && gzip -t "$new_download_dir/$fname" 2>/dev/null; then
                            cp -p -- "$new_download_dir/$fname" "$home_download_dir/$fname"
                            success=1
                            break
                        fi

                        rm -f -- "$new_download_dir/$fname"
                        sleep "$RETRY_DELAY_SECONDS"
                    done

                    if (( success == 0 )); then
                        printf '%s\n' "$fname" >> "$failed_downloads"
                    fi
                done < "$retry_urls"
                echo
            fi
        fi

        if [[ -s "$unmapped_names" ]]; then
            cat "$unmapped_names" >> "$failed_downloads"
        fi

        if [[ -s "$failed_downloads" ]]; then
            LC_ALL=C sort -u "$failed_downloads" -o "$failed_downloads"
            remaining_count="$(wc -l < "$failed_downloads")"
            echo "[WARN] $remaining_count genome(s) still unavailable after retry round $round."
        else
            echo "[PASS] All previously failed genomes were recovered in retry round $round."
            break
        fi

        if (( round < FAILED_RETRY_ROUNDS )); then
            sleep "$RETRY_DELAY_SECONDS"
        fi
    done

    return 0
}

###############################################################################
# 3. Baixar apenas arquivos faltantes
###############################################################################
echo
echo "STEP 3: Downloading missing files..."

cd "$new_download_dir"

download_one() {
    local url="$1"
    local file
    local attempt

    file="$(basename "$url")"

    for ((attempt=1; attempt<=DOWNLOAD_ATTEMPTS; attempt++)); do
        rm -f -- "$file.aria2"

        if aria2c \
            --disable-ipv6=true \
            --continue=true \
            --auto-file-renaming=false \
            -x "$ARIA_CONNECTIONS" \
            -s "$ARIA_CONNECTIONS" \
            -o "$file" \
            "$url" > /dev/null 2>&1; then

            if [[ -f "$file" ]] && gzip -t "$file" 2>/dev/null; then
                cp -p -- "$file" "$home_download_dir/$file"
                return 0
            fi
        fi

        echo "[WARN] Attempt $attempt failed for $file" >&2
        rm -f -- "$file" "$file.aria2"
        sleep "$RETRY_DELAY_SECONDS"
    done

    printf '%s\n' "$file" >> "$failed_downloads"
    return 1
}

export -f download_one
export failed_downloads home_download_dir ARIA_CONNECTIONS DOWNLOAD_ATTEMPTS RETRY_DELAY_SECONDS

if (( to_download_count == 0 )); then
    echo "[PASS] Local folder already contains every genome in the catalog."
else
    echo "[INFO] Downloading $to_download_count genome(s)..."

    if command -v parallel >/dev/null 2>&1 && command -v aria2c >/dev/null 2>&1; then
        # --halt never: um download com falha não interrompe os demais.
        parallel \
            --bar \
            --halt never \
            -j "$DOWNLOAD_JOBS" \
            download_one :::: "$to_download_list" || true
    else
        echo "[WARN] aria2c and/or GNU parallel not found; using serial wget."

        download_checked=0
        while IFS= read -r url; do
            [[ -z "$url" ]] && continue

            ((download_checked+=1))
            file="$(basename "$url")"
            progress "$download_checked" "$to_download_count" "[DOWNLOAD] $file"

            success=0
            for ((attempt=1; attempt<=DOWNLOAD_ATTEMPTS; attempt++)); do
                rm -f -- "$file"

                if wget -4 -q -O "$file" "$url" && gzip -t "$file" 2>/dev/null; then
                    cp -p -- "$file" "$home_download_dir/$file"
                    success=1
                    break
                fi

                rm -f -- "$file"
                sleep "$RETRY_DELAY_SECONDS"
            done

            if (( success == 0 )); then
                printf '%s\n' "$file" >> "$failed_downloads"
            fi
        done < "$to_download_list"
        echo
    fi

    echo "[PASS] Initial download stage completed."
fi

# Nova rodada dedicada somente aos arquivos que falharam no lote inicial.
retry_failed_downloads "initial-downloads"

###############################################################################
# 4. Integridade local incremental
###############################################################################
echo
echo "STEP 4: Verifying local gzip integrity..."

current_state="$work_dir/current_local_state.tsv"
integrity_candidates="$work_dir/integrity_candidates.txt"
build_local_state > "$current_state"

if [[ "$FULL_GZIP_CHECK" == "1" ]]; then
    cut -f1 "$current_state" > "$integrity_candidates"
    integrity_mode="full (forced)"
elif [[ ! -s "$integrity_cache" ]]; then
    cut -f1 "$current_state" > "$integrity_candidates"
    integrity_mode="full (first cached run)"
else
    # Testar apenas arquivos novos ou cujo tamanho/mtime local mudou.
    awk -F '\t' '
        NR == FNR {
            previous[$1] = $2 FS $3
            next
        }
        {
            current = $2 FS $3
            if (!($1 in previous) || previous[$1] != current) {
                print $1
            }
        }
    ' "$integrity_cache" "$current_state" > "$integrity_candidates"
    integrity_mode="incremental"
fi

integrity_candidate_count="$(wc -l < "$integrity_candidates")"
: > "$corrupted_list"

echo "[INFO] Integrity mode : $integrity_mode"
echo "[INFO] Files to test  : $integrity_candidate_count"

if (( integrity_candidate_count == 0 )); then
    echo "[PASS] No new or locally modified gzip files require testing."
elif command -v parallel >/dev/null 2>&1; then
    export new_download_dir

    parallel \
        --bar \
        --keep-order \
        -j "$GZIP_CHECK_JOBS" \
        'gzip -t "$new_download_dir/{}" 2>/dev/null || printf "%s\n" "{}"' \
        :::: "$integrity_candidates" \
        > "$corrupted_list"
else
    checked_gz=0
    while IFS= read -r fname; do
        [[ -z "$fname" ]] && continue
        ((checked_gz+=1))
        progress "$checked_gz" "$integrity_candidate_count" "[GZIP TEST] $fname"

        gzip -t "$new_download_dir/$fname" 2>/dev/null \
            || printf '%s\n' "$fname" >> "$corrupted_list"
    done < "$integrity_candidates"
    echo
fi

corrupted_count="$(wc -l < "$corrupted_list")"

###############################################################################
# 4B. Rebaixar somente arquivos locais corrompidos
###############################################################################
if (( corrupted_count > 0 )); then
    echo "[WARN] $corrupted_count corrupted file(s) found. Re-downloading..."

    # Reconstruir mapa filename -> URL uma única vez em memória.
    declare -A remote_map=()
    while IFS=$'\t' read -r fname url; do
        [[ -n "$fname" && -n "$url" ]] || continue
        remote_map["$fname"]="$url"
    done < "$tmp_remote_index"

    corrupted_checked=0
    while IFS= read -r fname; do
        [[ -z "$fname" ]] && continue

        ((corrupted_checked+=1))
        url="${remote_map[$fname]:-}"

        if [[ -z "$url" ]]; then
            echo "[WARN] URL not found in manifest for $fname; skipping." >&2
            printf '%s\n' "$fname" >> "$failed_downloads"
            rm -f -- "$new_download_dir/$fname" "$home_download_dir/$fname"
            continue
        fi

        success=0
        for ((attempt=1; attempt<=DOWNLOAD_ATTEMPTS; attempt++)); do
            progress "$corrupted_checked" "$corrupted_count" \
                "[REDOWNLOAD] $fname attempt $attempt"

            rm -f -- "$new_download_dir/$fname" "$new_download_dir/$fname.aria2"

            if wget -4 -q -O "$new_download_dir/$fname" "$url" \
                && gzip -t "$new_download_dir/$fname" 2>/dev/null; then
                cp -p -- "$new_download_dir/$fname" "$home_download_dir/$fname"
                success=1
                break
            fi

            rm -f -- "$new_download_dir/$fname"
            sleep "$RETRY_DELAY_SECONDS"
        done

        if (( success == 0 )); then
            echo >&2
            echo "[WARN] Could not recover corrupted file after $DOWNLOAD_ATTEMPTS attempts: $fname" >&2
            printf '%s\n' "$fname" >> "$failed_downloads"
            rm -f -- "$new_download_dir/$fname" "$home_download_dir/$fname"
        fi
    done < "$corrupted_list"
    echo
else
    echo "[PASS] All tested gzip files passed integrity validation."
fi

# Remover duplicatas ocasionais do arquivo de falhas.
if [[ -s "$failed_downloads" ]]; then
    LC_ALL=C sort -u "$failed_downloads" -o "$failed_downloads"
fi

# Atualizar cache apenas depois de finalizar os testes e re-downloads.
final_state="$work_dir/final_local_state.tsv"
build_local_state > "$final_state"
mv -f "$final_state" "$integrity_cache"

rm -f "$corrupted_list"

###############################################################################
# 5. Relatório final
###############################################################################
final_local_count="$(find "$new_download_dir" -maxdepth 1 -type f -name '*.gz' | wc -l)"
failed_count="$(wc -l < "$failed_downloads")"

echo
echo "STEP 5: Final report"
echo "Remote genomes listed : $total_remote"
echo "Final local genomes    : $final_local_count"
echo "Obsolete removed       : $obsolete_count"
echo "Missing detected       : $missing_count"
echo "Requested downloads    : $to_download_count"
echo "Gzip files tested      : $integrity_candidate_count"
echo "Integrity mode         : $integrity_mode"
echo "Remote HEAD requests   : 0"
echo "Failed downloads       : $failed_count"
echo "Attempts per round     : $DOWNLOAD_ATTEMPTS"
echo "Extra retry rounds     : $FAILED_RETRY_ROUNDS"
echo "Local folder           : $new_download_dir"
echo "Home mirror folder     : $home_download_dir"
echo "Integrity cache        : $integrity_cache"

if [[ -s "$failed_downloads" ]]; then
    echo
    echo "[WARN] The following bacterial genomes could not be retrieved after all attempts:"
    cat "$failed_downloads"
    echo
    echo "[WARN] Synchronization completed with $failed_count unavailable genome(s)."
    echo "[WARN] These files were skipped. The script will finish with exit status 0."
else
    echo
    echo "[PASS] All bacterial genomes are synchronized and locally verified."
fi
