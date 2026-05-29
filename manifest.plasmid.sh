#!/bin/bash
set -euo pipefail

###############################################################################
# BASIC CONFIGURATION
###############################################################################
BASE_URL="https://ftp.ncbi.nlm.nih.gov"
REMOTE_DIR="genomes/refseq/plasmid"
LOCAL_DIR="/media/me/18TB_BACKUP_LBN/drive.ifpa/LBN_RNA-Seq/Metatranscriptomics/MTD/MTD_Offline_Install_files/Kraken2DB_micro/library/plasmid"

# Temporary files and failed download list
TMP_DIR="$HOME"
mkdir -p "$TMP_DIR"

FAILED_DL="$TMP_DIR/failed_downloads.txt"
TMP_REMOTE_HTML="$TMP_DIR/ncbi_remote_plasmid.html"
TMP_REMOTE_NAMES="$TMP_DIR/ncbi_remote_names.txt"
TMP_LOCAL_LIST="$TMP_DIR/local_list.txt"
CORRUPTED_LIST="$TMP_DIR/corrupted_files.txt"

# Number of retries per file
MAX_ATTEMPTS=3

###############################################################################
# CREATE / CHECK DESTINATION FOLDER
###############################################################################
mkdir -p "$LOCAL_DIR"

###############################################################################
# FETCH REMOTE FILE LIST
###############################################################################
echo "Fetching remote file list..."

curl -4 -fsSL "$BASE_URL/$REMOTE_DIR/" -o "$TMP_REMOTE_HTML"

grep -oP 'href="\K[^"]+\.gz(?=")' "$TMP_REMOTE_HTML" | sort -u > "$TMP_REMOTE_NAMES"

if [[ ! -s "$TMP_REMOTE_NAMES" ]]; then
  echo "ERROR: Failed to retrieve remote .gz file list from $BASE_URL/$REMOTE_DIR/"
  rm -f "$TMP_REMOTE_HTML" "$TMP_REMOTE_NAMES"
  exit 1
fi

###############################################################################
# FETCH LOCAL FILE LIST
###############################################################################
find "$LOCAL_DIR" -maxdepth 1 -type f -printf "%f\n" | sort > "$TMP_LOCAL_LIST"

###############################################################################
# DEFINE FILES TO BE DOWNLOADED
###############################################################################
mapfile -t ALL_REMOTE < "$TMP_REMOTE_NAMES"

MISSING=()
for f in "${ALL_REMOTE[@]}"; do
  [[ -f "$LOCAL_DIR/$f" ]] || MISSING+=("$f")
done

TOTAL=${#ALL_REMOTE[@]}
NEEDED=${#MISSING[@]}
ALREADY_LOCAL=$((TOTAL - NEEDED))

echo -e "\n────────  SUMMARY  ────────"
echo "Total on remote  : $TOTAL"
echo "Already local    : $ALREADY_LOCAL"
echo "To be downloaded : $NEEDED"
echo -e "───────────────────────────"

###############################################################################
# FUNCTIONS
###############################################################################
draw_bar() {
  local done=$1
  local need=$2
  local current_file="${3:-}"
  local width=40
  local pct=0
  local fill=0
  local empty=0

  if [[ "$need" -gt 0 ]]; then
    pct=$((100 * done / need))
    fill=$((pct * width / 100))
  fi

  empty=$((width - fill))

  printf "\r\033[KProgress: [%s%s] %d/%d (%d%%)  %s" \
    "$(printf '=%.0s' $(seq 1 "$fill" 2>/dev/null || true))" \
    "$(printf ' %.0s' $(seq 1 "$empty" 2>/dev/null || true))" \
    "$done" "$need" "$pct" "$current_file"
}

download_and_verify() {
  local file="$1"
  local attempt
  local ok=0
  local url="$BASE_URL/$REMOTE_DIR/$file"

  for attempt in $(seq 1 "$MAX_ATTEMPTS"); do

    rm -f "$LOCAL_DIR/$file"

    aria2c \
      --disable-ipv6=true \
      --continue=false \
      --auto-file-renaming=false \
      --allow-overwrite=true \
      -x16 -s16 \
      -d "$LOCAL_DIR" \
      -o "$file" \
      "$url" \
      > /dev/null 2>&1 || true

    if [[ ! -s "$LOCAL_DIR/$file" ]]; then
      echo
      echo "  ❌ Attempt [$attempt/$MAX_ATTEMPTS] failed: file missing or empty after download: $file"
      rm -f "$LOCAL_DIR/$file"
      continue
    fi

    if gzip -t "$LOCAL_DIR/$file" 2>/dev/null; then
      ok=1
      break
    else
      echo
      echo "  ❌ Attempt [$attempt/$MAX_ATTEMPTS] failed: corrupted or partial file: $file"
      echo "     Removing and retrying..."
      rm -f "$LOCAL_DIR/$file"
    fi
  done

  if [[ "$ok" -ne 1 ]]; then
    echo "$file" >> "$FAILED_DL"
    return 1
  fi

  return 0
}

###############################################################################
# DOWNLOAD LOOP
###############################################################################
> "$FAILED_DL"

if [[ "$NEEDED" -eq 0 ]]; then
  echo "Nothing new to download."
else
  DL_DONE=0

  for file in "${MISSING[@]}"; do
    download_and_verify "$file" || true
    ((++DL_DONE))
    draw_bar "$DL_DONE" "$NEEDED" "$file"
  done

  echo
fi

###############################################################################
# FINAL EXTRA INTEGRITY CHECK
###############################################################################
echo -e "\nVerifying integrity of all local .gz files..."

> "$CORRUPTED_LIST"

shopt -s nullglob
for file in "$LOCAL_DIR"/*.gz; do
  if ! gzip -t "$file" 2>/dev/null; then
    echo "$(basename "$file")" >> "$CORRUPTED_LIST"
  fi
done
shopt -u nullglob

CORRUPTED_COUNT=$(wc -l < "$CORRUPTED_LIST")

if [[ "$CORRUPTED_COUNT" -gt 0 ]]; then
  echo "⚠️ Found $CORRUPTED_COUNT corrupted files after final verification."
  echo "Removing corrupted files..."

  while read -r file; do
    [[ -n "$file" ]] || continue
    rm -f "$LOCAL_DIR/$file"
    echo "$file" >> "$FAILED_DL"
  done < "$CORRUPTED_LIST"
else
  echo "✅ All .gz files passed final integrity check."
fi

rm -f "$CORRUPTED_LIST"

###############################################################################
# FINAL REPORT
###############################################################################
if [[ -s "$FAILED_DL" ]]; then
  echo -e "\n⚠️ Some files failed after $MAX_ATTEMPTS attempt(s):"
  sort -u "$FAILED_DL"
else
  echo -e "\n✅ All downloads completed successfully!"
fi

###############################################################################
# CLEANUP
###############################################################################
rm -f "$TMP_REMOTE_HTML" "$TMP_REMOTE_NAMES" "$TMP_LOCAL_LIST"
