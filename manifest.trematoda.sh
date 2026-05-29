#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Manifest Trematoda / Schistosomatidae / Schistosoma for Kraken2 DB_micro
#
# This script:
#   1. Downloads/synchronizes RefSeq invertebrate assembly_summary.txt
#   2. Uses NCBI taxonomy names.dmp/nodes.dmp to find descendant taxids
#   3. Filters assemblies belonging to the selected target clades
#   4. Downloads genomic.fna.gz files
#   5. Verifies gzip integrity
#   6. Creates FASTA files with kraken:taxid headers
#   7. Writes add_to_library.list.txt for kraken2-build --add-to-library
###############################################################################

###############################################################################
# 0. User-editable settings
###############################################################################

offline_files_folder="/media/me/4TB_BACKUP_LBN/Compressed/MTD"

# RefSeq group to search.
# For Schistosoma and trematodes, this should normally be invertebrate.
NCBI_GROUP="${NCBI_GROUP:-invertebrate}"

# Target clades/species. The script finds all descendants of these names.
# You can override before running:
#   TARGET_NAMES="Schistosomatidae,Schistosoma" bash manifest.trematoda.sh
TARGET_NAMES="${TARGET_NAMES:-Trematoda,Schistosomatidae,Schistosoma}"

# Assembly levels to keep.
# Good compromise: Complete Genome, Chromosome, Scaffold.
# If too few results, include Contig:
#   ASSEMBLY_LEVELS="Complete Genome,Chromosome,Scaffold,Contig"
ASSEMBLY_LEVELS="${ASSEMBLY_LEVELS:-Complete Genome,Chromosome,Scaffold}"

# Number of parallel downloads if GNU parallel + aria2c are available.
DOWNLOAD_JOBS="${DOWNLOAD_JOBS:-4}"

###############################################################################
# 1. Paths
###############################################################################

base_dir="$offline_files_folder/Kraken2DB_micro/library/trematoda"
new_download_dir="$base_dir/all"
formatted_dir="$base_dir/formatted"

home_base_dir="$HOME/MTD/kraken2DB_micro/library/trematoda"
home_download_dir="$home_base_dir/all"

taxdump_dir="$offline_files_folder/Kraken2DB_micro/taxonomy"
fallback_taxdump_dir="$base_dir/taxdump"

assembly_summary_file="$base_dir/assembly_summary_${NCBI_GROUP}.txt"
manifest_tsv="$base_dir/manifest_trematoda.tsv"
manifest_list="$base_dir/manifest_trematoda.list.txt"

target_taxids_file="$base_dir/target_taxids.txt"
target_taxids_named_file="$base_dir/target_taxids_named.tsv"
filtered_assembly_summary="$base_dir/assembly_summary_trematoda_filtered.tsv"

failed_downloads="$offline_files_folder/failed_downloads_trematoda.txt"
corrupted_list="$offline_files_folder/corrupted_trematoda.txt"
to_download_list="$offline_files_folder/to_download_trematoda.txt"
obsolete_list="$offline_files_folder/obsolete_local_trematoda.txt"

add_to_library_list="$base_dir/add_to_library.list.txt"
format_log="$base_dir/format_fasta.log"

mkdir -p "$base_dir" "$new_download_dir" "$formatted_dir"
mkdir -p "$home_base_dir" "$home_download_dir"

rm -f "$manifest_tsv" "$manifest_list" "$target_taxids_file" \
      "$target_taxids_named_file" "$filtered_assembly_summary" \
      "$failed_downloads" "$corrupted_list" "$to_download_list" \
      "$obsolete_list" "$add_to_library_list" "$format_log"

progress() {
    local current="$1"
    local total="$2"
    local msg="$3"

    msg="${msg:0:110}"
    printf "\r\033[2K[%d/%d] %s" "$current" "$total" "$msg"
}

###############################################################################
# 2. Dependency checks
###############################################################################

echo "============================================================"
echo "MTD Trematoda manifest"
echo "offline_files_folder : $offline_files_folder"
echo "NCBI group           : $NCBI_GROUP"
echo "Target names         : $TARGET_NAMES"
echo "Assembly levels      : $ASSEMBLY_LEVELS"
echo "Base dir             : $base_dir"
echo "============================================================"

if ! command -v python3 >/dev/null 2>&1; then
    echo "[ERROR] python3 not found in PATH."
    exit 1
fi

if ! command -v curl >/dev/null 2>&1 && ! command -v wget >/dev/null 2>&1; then
    echo "[ERROR] Neither curl nor wget was found."
    exit 1
fi

###############################################################################
# 3. Get taxonomy files
###############################################################################

echo
echo "STEP 1: Locating NCBI taxonomy files..."

NAMES=""
NODES=""

candidate_tax_dirs=(
  "$offline_files_folder/Kraken2DB_micro/taxonomy"
  "$HOME/MTD/kraken2DB_micro/taxonomy"
  "$PWD/kraken2DB_micro/taxonomy"
  "$fallback_taxdump_dir"
)

for td in "${candidate_tax_dirs[@]}"; do
    if [[ -s "$td/names.dmp" && -s "$td/nodes.dmp" ]]; then
        NAMES="$td/names.dmp"
        NODES="$td/nodes.dmp"
        break
    fi
done

if [[ -z "$NAMES" || -z "$NODES" ]]; then
    echo "No local names.dmp/nodes.dmp found. Downloading taxdump..."
    mkdir -p "$fallback_taxdump_dir"

    taxdump_zip="$fallback_taxdump_dir/taxdmp.zip"

    if command -v curl >/dev/null 2>&1; then
        curl -4 --retry 5 --retry-delay 2 --connect-timeout 30 -fsSL \
          -o "$taxdump_zip" \
          "https://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip"
    else
        wget -4 -q -O "$taxdump_zip" \
          "https://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip"
    fi

    if command -v unzip >/dev/null 2>&1; then
        unzip -o -q "$taxdump_zip" -d "$fallback_taxdump_dir"
    else
        echo "[ERROR] unzip is required to extract taxdmp.zip."
        exit 1
    fi

    NAMES="$fallback_taxdump_dir/names.dmp"
    NODES="$fallback_taxdump_dir/nodes.dmp"
fi

if [[ ! -s "$NAMES" || ! -s "$NODES" ]]; then
    echo "[ERROR] Could not prepare names.dmp and nodes.dmp."
    exit 1
fi

echo "names.dmp: $NAMES"
echo "nodes.dmp: $NODES"

###############################################################################
# 4. Prepare/download assembly_summary
###############################################################################

echo
echo "STEP 2: Preparing RefSeq assembly_summary.txt for $NCBI_GROUP..."

assembly_url="https://ftp.ncbi.nih.gov/genomes/refseq/${NCBI_GROUP}/assembly_summary.txt"
tmp_assembly="${assembly_summary_file}.tmp.$$"
download_ok="no"

use_local_assembly_summary() {
    local candidates=(
        "$assembly_summary_file"
        "$base_dir/assembly_summary.txt"
        "$base_dir/assembly_summary_${NCBI_GROUP}.txt"
        "$offline_files_folder/Kraken2DB_micro/library/${NCBI_GROUP}/assembly_summary.txt"
        "$offline_files_folder/Kraken2DB_micro/library/${NCBI_GROUP}/assembly_summary_${NCBI_GROUP}.txt"
        "$offline_files_folder/Kraken2DB_micro/library/invertebrate/assembly_summary.txt"
        "$offline_files_folder/Kraken2DB_micro/library/invertebrate/assembly_summary_invertebrate.txt"
        "$HOME/MTD/kraken2DB_micro/library/${NCBI_GROUP}/assembly_summary.txt"
        "$HOME/MTD/kraken2DB_micro/library/${NCBI_GROUP}/assembly_summary_${NCBI_GROUP}.txt"
        "$HOME/MTD/kraken2DB_micro/library/invertebrate/assembly_summary.txt"
        "$HOME/MTD/kraken2DB_micro/library/invertebrate/assembly_summary_invertebrate.txt"
    )

    for cand in "${candidates[@]}"; do
        if [[ -s "$cand" ]]; then
            echo "[INFO] Using local assembly_summary:"
            echo "       $cand"

            mkdir -p "$(dirname "$assembly_summary_file")"

            if [[ "$(readlink -f "$cand" 2>/dev/null || echo "$cand")" != "$(readlink -f "$assembly_summary_file" 2>/dev/null || echo "$assembly_summary_file")" ]]; then
                cp -f "$cand" "$assembly_summary_file"
            fi

            return 0
        fi
    done

    return 1
}

if [[ "${OFFLINE:-auto}" == "1" || "${OFFLINE:-auto}" == "yes" ]]; then
    echo "[INFO] OFFLINE mode requested. Skipping download."
else
    echo "[INFO] Trying to download:"
    echo "       $assembly_url"

    rm -f "$tmp_assembly"

    if command -v curl >/dev/null 2>&1; then
        if curl -4 --retry 3 --retry-delay 2 --connect-timeout 30 -fsSL \
          -o "$tmp_assembly" \
          "$assembly_url"; then
            if [[ -s "$tmp_assembly" ]]; then
                mv -f "$tmp_assembly" "$assembly_summary_file"
                download_ok="yes"
            fi
        fi
    else
        if wget -4 -q -O "$tmp_assembly" "$assembly_url"; then
            if [[ -s "$tmp_assembly" ]]; then
                mv -f "$tmp_assembly" "$assembly_summary_file"
                download_ok="yes"
            fi
        fi
    fi

    rm -f "$tmp_assembly"
fi

if [[ "$download_ok" != "yes" ]]; then
    echo "[WARNING] Could not download assembly_summary from NCBI."
    echo "[WARNING] Trying to reuse a local/offline assembly_summary file..."

    if ! use_local_assembly_summary; then
        echo
        echo "[ERROR] No local assembly_summary file found for group: $NCBI_GROUP"
        echo
        echo "To run fully offline, place this file here:"
        echo "  $assembly_summary_file"
        echo
        echo "Expected source when online:"
        echo "  $assembly_url"
        echo
        echo "Example on a machine with internet:"
        echo "  wget -O assembly_summary_invertebrate.txt $assembly_url"
        echo
        echo "Then copy it to:"
        echo "  $assembly_summary_file"
        exit 1
    fi
fi

if [[ ! -s "$assembly_summary_file" ]]; then
    echo "[ERROR] assembly_summary file is empty or missing:"
    echo "        $assembly_summary_file"
    exit 1
fi

echo "assembly_summary: $assembly_summary_file"

###############################################################################
# 5. Resolve target taxids and filter assembly_summary
###############################################################################

echo
echo "STEP 3: Resolving target taxids and filtering assemblies..."

python3 - "$NAMES" "$NODES" "$assembly_summary_file" "$TARGET_NAMES" "$ASSEMBLY_LEVELS" \
  "$target_taxids_file" "$target_taxids_named_file" "$filtered_assembly_summary" "$manifest_tsv" "$manifest_list" << 'PY'
import csv
import sys
from collections import defaultdict, deque

names_dmp = sys.argv[1]
nodes_dmp = sys.argv[2]
assembly_summary = sys.argv[3]
target_names = [x.strip() for x in sys.argv[4].split(",") if x.strip()]
assembly_levels = set(x.strip() for x in sys.argv[5].split(",") if x.strip())

target_taxids_file = sys.argv[6]
target_taxids_named_file = sys.argv[7]
filtered_assembly_summary = sys.argv[8]
manifest_tsv = sys.argv[9]
manifest_list = sys.argv[10]

name_to_taxids = defaultdict(list)
taxid_to_name = {}

with open(names_dmp, "r", encoding="utf-8", errors="replace") as f:
    for line in f:
        parts = [p.strip() for p in line.split("|")]
        if len(parts) < 4:
            continue
        taxid, name_txt, unique_name, name_class = parts[:4]
        if name_class == "scientific name":
            name_to_taxids[name_txt.lower()].append(taxid)
            taxid_to_name[taxid] = name_txt

children = defaultdict(list)

with open(nodes_dmp, "r", encoding="utf-8", errors="replace") as f:
    for line in f:
        parts = [p.strip() for p in line.split("|")]
        if len(parts) < 2:
            continue
        taxid, parent_taxid = parts[0], parts[1]
        children[parent_taxid].append(taxid)

root_taxids = []

for name in target_names:
    hits = name_to_taxids.get(name.lower(), [])
    if not hits:
        sys.stderr.write(f"[WARNING] Target name not found in taxonomy: {name}\n")
    for h in hits:
        root_taxids.append(h)

root_taxids = sorted(set(root_taxids), key=lambda x: int(x) if x.isdigit() else x)

desc = set(root_taxids)
q = deque(root_taxids)

while q:
    current = q.popleft()
    for child in children.get(current, []):
        if child not in desc:
            desc.add(child)
            q.append(child)

with open(target_taxids_file, "w", encoding="utf-8") as out:
    for taxid in sorted(desc, key=lambda x: int(x) if x.isdigit() else x):
        out.write(taxid + "\n")

with open(target_taxids_named_file, "w", encoding="utf-8") as out:
    out.write("taxid\tname\n")
    for taxid in sorted(desc, key=lambda x: int(x) if x.isdigit() else x):
        out.write(f"{taxid}\t{taxid_to_name.get(taxid, '')}\n")

rows_kept = []
manifest_rows = []

with open(assembly_summary, "r", encoding="utf-8", errors="replace") as f:
    header = None

    for line in f:
        line = line.rstrip("\n")

        if line.startswith("#assembly_accession"):
            header = line.lstrip("#").split("\t")
            continue

        if line.startswith("#") or not line.strip():
            continue

        parts = line.split("\t")

        if len(parts) < 20:
            continue

        assembly_accession = parts[0]
        refseq_category = parts[4]
        taxid = parts[5]
        species_taxid = parts[6]
        organism_name = parts[7]
        version_status = parts[10]
        assembly_level = parts[11]
        ftp_path = parts[19]

        if ftp_path == "na":
            continue

        if version_status.lower() != "latest":
            continue

        if assembly_levels and assembly_level not in assembly_levels:
            continue

        if taxid not in desc and species_taxid not in desc:
            continue

        ftp_path = ftp_path.replace("ftp://", "https://")
        ftp_path = ftp_path.replace("ftp.ncbi.nlm.nih.gov", "ftp.ncbi.nih.gov")
        ftp_path = ftp_path.rstrip("/")

        asm = ftp_path.split("/")[-1]
        url = f"{ftp_path}/{asm}_genomic.fna.gz"
        fname = f"{asm}_genomic.fna.gz"

        rows_kept.append(parts)
        manifest_rows.append({
            "url": url,
            "filename": fname,
            "taxid": taxid,
            "species_taxid": species_taxid,
            "organism_name": organism_name,
            "assembly_accession": assembly_accession,
            "assembly_level": assembly_level,
            "refseq_category": refseq_category
        })

if header is None:
    header = [
        "assembly_accession", "bioproject", "biosample", "wgs_master",
        "refseq_category", "taxid", "species_taxid", "organism_name",
        "infraspecific_name", "isolate", "version_status", "assembly_level",
        "release_type", "genome_rep", "seq_rel_date", "asm_name",
        "submitter", "gbrs_paired_asm", "paired_asm_comp", "ftp_path",
        "excluded_from_refseq", "relation_to_type_material"
    ]

with open(filtered_assembly_summary, "w", encoding="utf-8", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow(header)
    for r in rows_kept:
        writer.writerow(r)

with open(manifest_tsv, "w", encoding="utf-8", newline="") as out:
    fieldnames = [
        "url", "filename", "taxid", "species_taxid", "organism_name",
        "assembly_accession", "assembly_level", "refseq_category"
    ]
    writer = csv.DictWriter(out, delimiter="\t", fieldnames=fieldnames)
    writer.writeheader()
    for r in manifest_rows:
        writer.writerow(r)

with open(manifest_list, "w", encoding="utf-8") as out:
    for r in manifest_rows:
        out.write(r["url"] + "\n")

sys.stderr.write(f"[INFO] Root target taxids: {', '.join(root_taxids) if root_taxids else 'none'}\n")
sys.stderr.write(f"[INFO] Descendant taxids: {len(desc)}\n")
sys.stderr.write(f"[INFO] Assemblies retained: {len(manifest_rows)}\n")
PY

echo "Target taxids:"
wc -l "$target_taxids_file"

echo "Assemblies retained:"
wc -l "$manifest_list"

if [[ ! -s "$manifest_list" ]]; then
    echo
    echo "[WARNING] No assemblies were retained."
    echo "Possible fixes:"
    echo "  1. Use broader target names:"
    echo "     TARGET_NAMES='Platyhelminthes,Trematoda,Schistosomatidae,Schistosoma'"
    echo
    echo "  2. Include Contig assemblies:"
    echo "     ASSEMBLY_LEVELS='Complete Genome,Chromosome,Scaffold,Contig'"
    echo
    echo "  3. Check files:"
    echo "     $target_taxids_named_file"
    echo "     $filtered_assembly_summary"
    exit 0
fi

###############################################################################
# 6. Sync local folder against manifest
###############################################################################

echo
echo "STEP 4: Syncing local Trematoda folder against NCBI manifest..."

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
# 6A. Remove obsolete local files
###############################################################################

echo
echo "STEP 4A: Checking LOCAL files for obsolete files..."

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
            echo "Removing obsolete local file: $fname"
            echo "$fname" >> "$obsolete_list"
            rm -f "$path"
            rm -f "$home_download_dir/$fname"
            ((obsolete_count+=1))
        fi
    done
    echo
fi

###############################################################################
# 6B. Build local map after obsolete removal
###############################################################################

shopt -s nullglob
local_paths=("$new_download_dir"/*.gz)
shopt -u nullglob

declare -A local_map

for path in "${local_paths[@]}"; do
    fname="$(basename "$path")"
    local_map["$fname"]="$path"
done

get_remote_size() {
    local url="$1"

    if command -v curl >/dev/null 2>&1; then
        curl -4 -fsSI \
          --retry 3 \
          --retry-delay 2 \
          --connect-timeout 20 \
          "$url" 2>/dev/null \
          | awk 'BEGIN{IGNORECASE=1} /^Content-Length:/ {gsub("\r","",$2); print $2; exit}'
    else
        wget --spider --server-response "$url" 2>&1 \
          | awk 'BEGIN{IGNORECASE=1} /^  Content-Length:/ {print $2; exit}'
    fi
}

###############################################################################
# 6C. Detect missing and changed files
###############################################################################

echo
echo "STEP 4B: Checking REMOTE files against LOCAL files..."

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
# 7. Download missing/changed files
###############################################################################

echo
echo "STEP 5: Downloading missing/changed files..."

cd "$new_download_dir"

download_one() {
    local url="$1"
    local file
    file="$(basename "$url")"

    for attempt in {1..3}; do
        if command -v aria2c >/dev/null 2>&1; then
            aria2c --disable-ipv6=true \
                   --continue=true \
                   --auto-file-renaming=false \
                   -x16 -s16 \
                   -o "$file" \
                   "$url" > /dev/null 2>&1 || true
        else
            wget -4 -q -c -O "$file" "$url" || true
        fi

        if [[ -f "$file" ]] && gzip -t "$file" 2>/dev/null; then
            [[ -f "$home_download_dir/$file" ]] || cp -p "$file" "$home_download_dir/"
            return 0
        else
            echo "Attempt $attempt failed for $file"
            rm -f "$file"
            sleep 1
        fi
    done

    echo "Failed after 3 attempts: $file" >> "$failed_downloads"
}
export -f download_one
export failed_downloads
export home_download_dir

if (( to_download_count == 0 )); then
    echo "Local Trematoda folder is already synchronized with NCBI."
else
    echo "Downloading $to_download_count genome(s)..."
    mapfile -t download_urls < "$to_download_list"

    if command -v parallel >/dev/null 2>&1; then
        printf "%s\n" "${download_urls[@]}" | parallel --bar -j "$DOWNLOAD_JOBS" download_one
    else
        echo "parallel not found; using serial download."

        download_checked=0

        for url in "${download_urls[@]}"; do
            ((download_checked+=1))
            file="$(basename "$url")"

            progress "$download_checked" "$to_download_count" "[DOWNLOAD] downloading: $file"

            download_one "$url"
        done

        echo
    fi

    echo "Download stage completed."
fi

###############################################################################
# 8. Verify integrity and redownload corrupted files
###############################################################################

echo
echo "STEP 6: Verifying integrity of .gz files..."

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
    echo "No .gz files found yet in $new_download_dir"
fi

corrupted_count=$(wc -l < "$corrupted_list" 2>/dev/null || echo 0)

if (( corrupted_count > 0 )); then
    echo "$corrupted_count corrupted file(s) found. Re-downloading..."

    corrupted_checked=0

    while IFS= read -r fname; do
        [[ -z "$fname" ]] && continue

        ((corrupted_checked+=1))
        url="${remote_map[$fname]:-}"

        if [[ -z "$url" ]]; then
            echo "URL not found in manifest for $fname" >> "$failed_downloads"
            rm -f "$new_download_dir/$fname"
            rm -f "$home_download_dir/$fname"
            continue
        fi

        for attempt in {1..3}; do
            progress "$corrupted_checked" "$corrupted_count" "[REDOWNLOAD] $fname attempt $attempt"

            if command -v wget >/dev/null 2>&1; then
                wget -4 -q -O "$new_download_dir/$fname" "$url" || true
            else
                curl -4 -fsSL -o "$new_download_dir/$fname" "$url" || true
            fi

            if gzip -t "$new_download_dir/$fname" 2>/dev/null; then
                echo
                echo "Integrity OK after attempt $attempt: $fname"
                [[ -f "$home_download_dir/$fname" ]] || cp -p "$new_download_dir/$fname" "$home_download_dir/"
                break
            elif [[ $attempt -eq 3 ]]; then
                echo
                echo "Still corrupted after 3 attempts: $fname"
                echo "Failed after 3 attempts: $fname" >> "$failed_downloads"
                rm -f "$new_download_dir/$fname"
                rm -f "$home_download_dir/$fname"
            else
                rm -f "$new_download_dir/$fname"
            fi
        done
    done < "$corrupted_list"
else
    echo "All .gz files passed integrity check."
fi

rm -f "$corrupted_list"

###############################################################################
# 9. Format FASTA headers for Kraken2 custom library
###############################################################################

echo
echo "STEP 7: Formatting FASTA headers with kraken:taxid..."

python3 - "$manifest_tsv" "$new_download_dir" "$formatted_dir" "$add_to_library_list" "$format_log" << 'PY'
import csv
import gzip
import os
import re
import sys
from pathlib import Path

manifest_tsv = Path(sys.argv[1])
download_dir = Path(sys.argv[2])
formatted_dir = Path(sys.argv[3])
add_to_library_list = Path(sys.argv[4])
format_log = Path(sys.argv[5])

formatted_dir.mkdir(parents=True, exist_ok=True)

records = {}

with open(manifest_tsv, "r", encoding="utf-8", errors="replace") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for r in reader:
        records[r["filename"]] = r

def safe_name(x):
    x = re.sub(r"[^A-Za-z0-9._-]+", "_", x)
    return x.strip("_")

formatted_files = []

with open(format_log, "w", encoding="utf-8") as log:
    for filename, rec in sorted(records.items()):
        in_path = download_dir / filename

        if not in_path.exists():
            log.write(f"[SKIP missing] {filename}\n")
            continue

        taxid = rec.get("taxid") or rec.get("species_taxid")
        organism = rec.get("organism_name", "")
        accession = rec.get("assembly_accession", "")

        out_name = safe_name(filename.replace(".fna.gz", "")) + ".kraken.fna"
        out_path = formatted_dir / out_name

        log.write(f"[FORMAT] {filename} -> {out_path.name} taxid={taxid} organism={organism}\n")

        with gzip.open(in_path, "rt", encoding="utf-8", errors="replace") as inp, \
             open(out_path, "w", encoding="utf-8") as out:

            for line in inp:
                if line.startswith(">"):
                    header = line[1:].rstrip("\n")
                    seqid = header.split()[0]

                    if "kraken:taxid|" in header:
                        out.write(">" + header + "\n")
                    else:
                        out.write(f">{seqid}|kraken:taxid|{taxid} {header} [{organism}; {accession}]\n")
                else:
                    out.write(line)

        formatted_files.append(str(out_path))

with open(add_to_library_list, "w", encoding="utf-8") as out:
    for f in formatted_files:
        out.write(f + "\n")

print(f"[INFO] Formatted FASTA files: {len(formatted_files)}")
print(f"[INFO] add_to_library list: {add_to_library_list}")
PY

formatted_count=$(wc -l < "$add_to_library_list" 2>/dev/null || echo 0)

###############################################################################
# 10. Final report
###############################################################################

echo
echo "STEP 8: Final report"
echo "Target names             : $TARGET_NAMES"
echo "Assembly levels          : $ASSEMBLY_LEVELS"
echo "Remote assemblies listed : $total_remote"
echo "Obsolete removed         : $obsolete_count"
echo "Missing detected         : $missing_count"
echo "Changed detected         : $changed_count"
echo "Remote size check failed : $size_check_failed_count"
echo "Downloaded this run      : $to_download_count"
echo "Formatted FASTA files    : $formatted_count"
echo "Local raw folder         : $new_download_dir"
echo "Home mirror folder       : $home_download_dir"
echo "Formatted folder         : $formatted_dir"
echo "Add-to-library list      : $add_to_library_list"

if [[ -s "$failed_downloads" ]]; then
    echo
    echo "The following genomes could not be retrieved:"
    cat "$failed_downloads"
else
    echo
    echo "All Trematoda/Schistosoma genomes synchronized and verified successfully!"
fi

echo
echo "Next step:"
echo "  while read -r f; do"
echo "    kraken2-build --add-to-library \"\$f\" --db kraken2DB_micro"
echo "  done < \"$add_to_library_list\""
echo
