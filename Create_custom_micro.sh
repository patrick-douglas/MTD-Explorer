#!/usr/bin/env bash

# Create_custom_micro.sh
# Build a custom Kraken 2 / Bracken microbiome database for MTD Explorer.
#
# Modes:
#   1) Public-library mode:
#      add Kraken 2 public libraries such as bacteria, viral, archaea, fungi, protozoa.
#
#   2) Custom-FASTA mode:
#      add one or more user FASTA files.
#
#   3) Mixed mode:
#      combine public libraries and user FASTA files in the same database.
#
# Notes:
#   - For custom FASTA files that do not have NCBI accession mapping, use
#     --add-fasta-with-taxid FILE:TAXID to rewrite headers with kraken:taxid|TAXID.
#   - Bracken read length must match the MTD Explorer --bracken-read-len value.
#
# This script intentionally does not use "set -e" or explicit "exit 1".

SCRIPT_VERSION="0.1.0"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
OUTPUT_ROOT="$SCRIPT_DIR"
DB_NAME=""
DB_PATH=""

# ------------------------------------------------------------
# MTD dedicated Kraken2/Bracken builder runtime
# ------------------------------------------------------------
KRAKEN2_ENV_NAME="${MTD_KRAKEN2_ENV_NAME:-MTD_kraken2}"
CONDA_ROOT="${MTD_CONDA_ROOT:-}"
KRAKEN2_ENV_DIR=""
KRAKEN2_BIN=""
KRAKEN2_BUILD_BIN=""
KRAKEN2_INSPECT_BIN=""
BRACKEN_BUILD_BIN=""
KRAKEN2_VERSION=""
BRACKEN_PACKAGE_VERSION=""
KRAKEN2_RUNTIME_READY=0

THREADS="${THREADS:-}"
if [[ -z "$THREADS" ]]; then
  if command -v nproc >/dev/null 2>&1; then
    THREADS="$(nproc)"
  else
    THREADS="4"
  fi
fi

KMER_LEN="35"
MINIMIZER_LEN="31"
MINIMIZER_SPACES="7"
BRACKEN_READ_LEN="75"

DOWNLOAD_TAXONOMY="1"
BUILD_KRAKEN="1"
BUILD_BRACKEN="1"
CLEAN_DB="0"
FORCE_CLEAN="0"
VALIDATE_ONLY="0"
REBUILD_BRACKEN_ONLY="0"
USE_FTP="0"

LIBRARIES=()
FASTA_FILES=()
FASTA_WITH_TAXID=()
FASTA_LISTS=()

print_help() {
  cat <<'HELP'
Create_custom_micro.sh - build a custom Kraken 2 / Bracken microbiome database

USAGE

  Public libraries only:

    bash Create_custom_micro.sh \
      --db-name kraken2DB_micro_custom \
      --libraries bacteria,viral,archaea,fungi,protozoa \
      --threads 20 \
      --bracken-read-len 75

  Custom FASTA with TaxID already in headers or accession-resolvable headers:

    bash Create_custom_micro.sh \
      --db-name kraken2DB_micro_custom \
      --add-fasta /path/to/sequences.fa \
      --threads 20

  Custom FASTA where all records belong to one TaxID:

    bash Create_custom_micro.sh \
      --db-name kraken2DB_micro_custom \
      --add-fasta-with-taxid /path/to/sequences.fa:10239 \
      --threads 20

  Mixed public + custom database:

    bash Create_custom_micro.sh \
      --db-name kraken2DB_micro_custom \
      --libraries bacteria,viral \
      --add-fasta-with-taxid /path/to/project_virus.fa:10239 \
      --threads 20

REQUIRED OUTPUT OPTION

  --db-name NAME
      Create the database under ./NAME or under --output-root.
      Example: --db-name kraken2DB_micro_custom

  --output-db PATH
      Full path to the output database directory.
      Example: --output-db /data/kraken2DB_micro_custom

INPUT OPTIONS

  --library LIB
      Add one Kraken 2 public library. Can be repeated.
      Examples: bacteria, viral, archaea, fungi, protozoa, plasmid, UniVec_Core

  --libraries LIB1,LIB2,LIB3
      Add multiple Kraken 2 public libraries using a comma-separated list.

  --add-fasta FILE
      Add a custom FASTA file directly with kraken2-build --add-to-library.
      Use this when headers already contain kraken:taxid|TAXID or are accession-resolvable.

  --add-fasta-with-taxid FILE:TAXID
      Add a custom FASTA file and assign the same NCBI TaxID to every sequence.
      The script creates a prepared FASTA with kraken:taxid|TAXID in each header.

  --fasta-list FILE
      Read FASTA inputs from a text file.
      Accepted line formats:
        /path/to/file.fa
        /path/to/file.fa<TAB>10239
        /path/to/file.fa,10239
      Lines starting with # are ignored.

BUILD OPTIONS

  --threads N
      Number of threads. Default: nproc when available.

  --kmer-len N
      Kraken 2 k-mer length. Default: 35.

  --minimizer-len N
      Kraken 2 minimizer length. Default: 31.

  --minimizer-spaces N
      Kraken 2 minimizer spaces. Default: 7.

  --bracken-read-len N
      Bracken read length. Default: 75.

  --output-root DIR
      Root directory used with --db-name. Default: directory containing this script.

  --download-taxonomy
      Download NCBI taxonomy before adding libraries. Default.

  --no-download-taxonomy
      Do not download taxonomy. Use only if taxonomy is already present in the DB.

  --use-ftp
      Pass --use-ftp to Kraken 2 taxonomy/library download commands.

  --skip-kraken-build
      Add/download files but do not run kraken2-build --build.

  --skip-bracken
      Do not run bracken-build.

  --validate-only
      Only validate an existing database and print a short report.

  --rebuild-bracken-only
      Rebuild only Bracken files for an existing Kraken 2 database.
      This is useful when the Kraken 2 DB already exists but the Bracken
      read length needs to be changed.
      Equivalent to:
        --no-download-taxonomy --skip-kraken-build
      No microbial input is required in this mode.

CLEANING OPTIONS

  --clean
      Remove the existing target database directory before building.
      The path is safety-checked before removal.

  --force
      Allow --clean even when the database path does not look like a Kraken/MTD DB path.

EXAMPLES

  Broad microbial database:

    bash Create_custom_micro.sh \
      --db-name kraken2DB_micro_custom \
      --libraries bacteria,viral,archaea,fungi,protozoa \
      --threads 20 \
      --bracken-read-len 75

  Virus-only targeted database:

    bash Create_custom_micro.sh \
      --db-name kraken2DB_viral_custom \
      --libraries viral \
      --threads 20 \
      --bracken-read-len 75

  Project FASTA where all sequences are viral root TaxID 10239:

    bash Create_custom_micro.sh \
      --db-name kraken2DB_project_viral \
      --add-fasta-with-taxid project_viral_sequences.fa:10239 \
      --threads 20

  Use in MTD Explorer:

    bash MTD_explorer.sh \
      -i samplesheet.csv \
      -o results \
      -h <host_taxon_id> \
      --kraken-micro-db /path/to/kraken2DB_micro_custom \
      --bracken-read-len 75

HELP
}

log() {
  echo "[CUSTOM_MICRO] $*"
}

warn() {
  echo "[WARN] $*" >&2
}

fail() {
  echo "[ERROR] $*" >&2
  return 1
}

append_library_csv() {
  local csv="$1"
  local old_ifs="$IFS"
  IFS=","
  for item in $csv; do
    item="$(echo "$item" | sed 's/^ *//; s/ *$//')"
    if [[ -n "$item" ]]; then
      LIBRARIES+=("$item")
    fi
  done
  IFS="$old_ifs"
}

read_fasta_list() {
  local list_file="$1"

  if [[ ! -s "$list_file" ]]; then
    warn "FASTA list not found or empty: $list_file"
    return 1
  fi

  while IFS= read -r line || [[ -n "$line" ]]; do
    line="$(echo "$line" | sed 's/\r$//')"

    [[ -z "$line" ]] && continue
    [[ "$line" =~ ^[[:space:]]*# ]] && continue

    if [[ "$line" == *$'\t'* ]]; then
      local file_part tax_part
      file_part="$(echo "$line" | cut -f1)"
      tax_part="$(echo "$line" | cut -f2)"
      if [[ -n "$file_part" && -n "$tax_part" ]]; then
        FASTA_WITH_TAXID+=("${file_part}:${tax_part}")
      fi
    elif [[ "$line" == *,* ]]; then
      local file_part tax_part
      file_part="${line%%,*}"
      tax_part="${line#*,}"
      if [[ -n "$file_part" && -n "$tax_part" ]]; then
        FASTA_WITH_TAXID+=("${file_part}:${tax_part}")
      fi
    else
      FASTA_FILES+=("$line")
    fi
  done < "$list_file"
}

parse_args() {
  while [[ $# -gt 0 ]]; do
    case "$1" in
      -h|--help)
        print_help
        return 2
        ;;
      --version)
        echo "$SCRIPT_VERSION"
        return 2
        ;;
      --db-name)
        DB_NAME="$2"
        shift 2
        ;;
      --output-db)
        DB_PATH="$2"
        shift 2
        ;;
      --output-root)
        OUTPUT_ROOT="$2"
        shift 2
        ;;
      --threads|-t)
        THREADS="$2"
        shift 2
        ;;
      --kmer-len)
        KMER_LEN="$2"
        shift 2
        ;;
      --minimizer-len)
        MINIMIZER_LEN="$2"
        shift 2
        ;;
      --minimizer-spaces)
        MINIMIZER_SPACES="$2"
        shift 2
        ;;
      --bracken-read-len|-r)
        BRACKEN_READ_LEN="$2"
        shift 2
        ;;
      --library)
        LIBRARIES+=("$2")
        shift 2
        ;;
      --libraries)
        append_library_csv "$2"
        shift 2
        ;;
      --add-fasta)
        FASTA_FILES+=("$2")
        shift 2
        ;;
      --add-fasta-with-taxid)
        FASTA_WITH_TAXID+=("$2")
        shift 2
        ;;
      --fasta-list)
        FASTA_LISTS+=("$2")
        shift 2
        ;;
      --download-taxonomy)
        DOWNLOAD_TAXONOMY="1"
        shift
        ;;
      --no-download-taxonomy)
        DOWNLOAD_TAXONOMY="0"
        shift
        ;;
      --use-ftp)
        USE_FTP="1"
        shift
        ;;
      --skip-kraken-build)
        BUILD_KRAKEN="0"
        shift
        ;;
      --skip-bracken)
        BUILD_BRACKEN="0"
        shift
        ;;
      --validate-only)
        VALIDATE_ONLY="1"
        DOWNLOAD_TAXONOMY="0"
        BUILD_KRAKEN="0"
        BUILD_BRACKEN="0"
        shift
        ;;
      --rebuild-bracken-only)
        REBUILD_BRACKEN_ONLY="1"
        DOWNLOAD_TAXONOMY="0"
        BUILD_KRAKEN="0"
        BUILD_BRACKEN="1"
        shift
        ;;
      --clean)
        CLEAN_DB="1"
        shift
        ;;
      --force)
        FORCE_CLEAN="1"
        shift
        ;;
      *)
        warn "Unknown option: $1"
        print_help
        return 1
        ;;
    esac
  done

  return 0
}

resolve_db_path() {
  if [[ -z "$DB_PATH" ]]; then
    if [[ -z "$DB_NAME" ]]; then
      fail "Use --db-name NAME or --output-db PATH."
      return 1
    fi
    DB_PATH="${OUTPUT_ROOT%/}/${DB_NAME}"
  fi

  DB_PATH="$(python3 - "$DB_PATH" <<'PY'
from pathlib import Path
import sys
print(Path(sys.argv[1]).expanduser().resolve())
PY
)"
}


resolve_dedicated_kraken_runtime() {
  if [[ "$KRAKEN2_RUNTIME_READY" == "1" ]]; then
    return 0
  fi

  if [[ -z "$CONDA_ROOT" && -s "$SCRIPT_DIR/condaPath" ]]; then
    IFS= read -r CONDA_ROOT < "$SCRIPT_DIR/condaPath"
    CONDA_ROOT="${CONDA_ROOT%$'\r'}"
  fi

  if [[ -z "$CONDA_ROOT" ]]; then
    CONDA_ROOT="$HOME/miniconda3"
  fi

  KRAKEN2_ENV_DIR="$CONDA_ROOT/envs/$KRAKEN2_ENV_NAME"
  KRAKEN2_BIN="$KRAKEN2_ENV_DIR/bin/kraken2"
  KRAKEN2_BUILD_BIN="$KRAKEN2_ENV_DIR/bin/kraken2-build"
  KRAKEN2_INSPECT_BIN="$KRAKEN2_ENV_DIR/bin/kraken2-inspect"
  BRACKEN_BUILD_BIN="$KRAKEN2_ENV_DIR/bin/bracken-build"

  local required
  for required in "$KRAKEN2_BIN" "$KRAKEN2_BUILD_BIN" "$KRAKEN2_INSPECT_BIN"; do
    if [[ ! -x "$required" ]]; then
      fail "Dedicated Kraken2 executable not found: $required"
      return 1
    fi
  done

  KRAKEN2_VERSION="$("$KRAKEN2_BIN" --version 2>/dev/null | awk 'NR == 1 { print $3 }')"
  if [[ "$KRAKEN2_VERSION" != "2.17.1" ]]; then
    fail "Expected Kraken2 2.17.1, observed: ${KRAKEN2_VERSION:-unknown}"
    return 1
  fi

  if [[ "$BUILD_BRACKEN" == "1" ]]; then
    if [[ ! -x "$BRACKEN_BUILD_BIN" ]]; then
      fail "Dedicated bracken-build not found: $BRACKEN_BUILD_BIN"
      return 1
    fi

    BRACKEN_PACKAGE_VERSION="$(find "$KRAKEN2_ENV_DIR/conda-meta" -maxdepth 1 -type f -name 'bracken-3.1p1-*.json' -printf '3.1p1\n' -quit 2>/dev/null)"
    if [[ "$BRACKEN_PACKAGE_VERSION" != "3.1p1" ]]; then
      fail "Expected Bracken Conda package 3.1p1 in $KRAKEN2_ENV_DIR"
      return 1
    fi
  fi

  KRAKEN2_RUNTIME_READY=1
  log "Dedicated Kraken2 runtime:"
  log "  environment: $KRAKEN2_ENV_NAME"
  log "  Kraken2:    $KRAKEN2_VERSION"
  log "  build:      $KRAKEN2_BUILD_BIN"
  log "  inspect:    $KRAKEN2_INSPECT_BIN"
  if [[ "$BUILD_BRACKEN" == "1" ]]; then
    log "  Bracken:    $BRACKEN_PACKAGE_VERSION"
    log "  build:      $BRACKEN_BUILD_BIN"
  fi
}

run_kraken2_build() {
  env PATH="$KRAKEN2_ENV_DIR/bin:$PATH" "$KRAKEN2_BUILD_BIN" "$@"
}

run_kraken2_inspect() {
  env PATH="$KRAKEN2_ENV_DIR/bin:$PATH" "$KRAKEN2_INSPECT_BIN" "$@"
}

run_bracken_build() {
  env PATH="$KRAKEN2_ENV_DIR/bin:$PATH" "$BRACKEN_BUILD_BIN" "$@"
}

check_dependencies() {
  local missing=0

  resolve_dedicated_kraken_runtime || return 1

  if ! command -v python3 >/dev/null 2>&1; then
    warn "Missing command: python3"
    missing=1
  fi

  if [[ "$missing" == "1" ]]; then
    fail "Missing required commands."
    return 1
  fi
}

safe_clean_db() {
  if [[ "$CLEAN_DB" != "1" ]]; then
    return 0
  fi

  if [[ -z "$DB_PATH" || "$DB_PATH" == "/" || "$DB_PATH" == "$HOME" ]]; then
    fail "Refusing to clean unsafe database path: '$DB_PATH'"
    return 1
  fi

  if [[ "$FORCE_CLEAN" != "1" ]]; then
    case "$DB_PATH" in
      *kraken*|*Kraken*|*KRAKEN*|*DB*|*db*)
        ;;
      *)
        fail "Refusing to clean path that does not look like a database path: $DB_PATH. Use --force only if you are sure."
        return 1
        ;;
    esac
  fi

  if [[ -d "$DB_PATH" ]]; then
    log "Cleaning existing database directory:"
    log "  $DB_PATH"
    rm -rf "$DB_PATH"
    local status=$?
    if [[ "$status" != "0" ]]; then
      fail "Could not remove existing database directory."
      return 1
    fi
  fi
}

prepare_fasta_with_taxid() {
  local input="$1"
  local taxid="$2"
  local out="$3"

  python3 - "$input" "$taxid" "$out" <<'PY'
from pathlib import Path
import gzip
import sys
import re

inp = Path(sys.argv[1]).expanduser()
taxid = sys.argv[2].strip()
out = Path(sys.argv[3]).expanduser()

if not inp.exists():
    print(f"[ERROR] FASTA not found: {inp}", file=sys.stderr)
    sys.exit(3)

if not re.fullmatch(r"[0-9]+", taxid):
    print(f"[ERROR] TaxID must be numeric: {taxid}", file=sys.stderr)
    sys.exit(4)

out.parent.mkdir(parents=True, exist_ok=True)

def open_in(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")

def open_out(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "wt", encoding="utf-8")
    return open(path, "wt", encoding="utf-8")

records = 0
with open_in(inp) as fin, open_out(out) as fout:
    for line in fin:
        if line.startswith(">"):
            records += 1
            header = line[1:].strip()
            if f"kraken:taxid|{taxid}" in header:
                fout.write(">" + header + "\n")
            elif "kraken:taxid|" in header:
                header = re.sub(r"\|?kraken:taxid\|[0-9]+", "", header).strip()
                fout.write(f">{header}|kraken:taxid|{taxid}\n")
            else:
                fout.write(f">{header}|kraken:taxid|{taxid}\n")
        else:
            fout.write(line)

print(f"[OK] Prepared {records} FASTA records with TaxID {taxid}: {out}")
PY
}

add_custom_fasta_files() {
  local prepared_dir="$DB_PATH/MTD_custom_fasta_prepared"
  mkdir -p "$prepared_dir"

  for fasta in "${FASTA_FILES[@]}"; do
    fasta="$(python3 - "$fasta" <<'PY'
from pathlib import Path
import sys
print(Path(sys.argv[1]).expanduser())
PY
)"
    if [[ ! -s "$fasta" ]]; then
      fail "Custom FASTA not found or empty: $fasta"
      return 1
    fi

    log "Adding custom FASTA:"
    log "  $fasta"
    run_kraken2_build --add-to-library "$fasta" --db "$DB_PATH"
    local status=$?
    if [[ "$status" != "0" ]]; then
      fail "kraken2-build --add-to-library failed for: $fasta"
      return 1
    fi
  done

  local idx=0
  for item in "${FASTA_WITH_TAXID[@]}"; do
    idx=$((idx + 1))

    local fasta taxid
    fasta="${item%:*}"
    taxid="${item##*:}"

    fasta="$(python3 - "$fasta" <<'PY'
from pathlib import Path
import sys
print(Path(sys.argv[1]).expanduser())
PY
)"

    if [[ ! -s "$fasta" ]]; then
      fail "Custom FASTA not found or empty: $fasta"
      return 1
    fi

    if [[ ! "$taxid" =~ ^[0-9]+$ ]]; then
      fail "TaxID must be numeric in --add-fasta-with-taxid FILE:TAXID: $item"
      return 1
    fi

    local base
    base="$(basename "$fasta")"
    base="${base%.gz}"
    local prepared="$prepared_dir/custom_${idx}_${taxid}_${base}"

    log "Preparing FASTA with TaxID:"
    log "  input:  $fasta"
    log "  taxid:  $taxid"
    log "  output: $prepared"

    prepare_fasta_with_taxid "$fasta" "$taxid" "$prepared"
    local prep_status=$?
    if [[ "$prep_status" != "0" ]]; then
      fail "Could not prepare FASTA with TaxID: $item"
      return 1
    fi

    log "Adding prepared FASTA:"
    log "  $prepared"
    run_kraken2_build --add-to-library "$prepared" --db "$DB_PATH"
    local add_status=$?
    if [[ "$add_status" != "0" ]]; then
      fail "kraken2-build --add-to-library failed for prepared FASTA: $prepared"
      return 1
    fi
  done
}

download_taxonomy() {
  if [[ "$DOWNLOAD_TAXONOMY" != "1" ]]; then
    log "Skipping taxonomy download."
    return 0
  fi

  local args=(--download-taxonomy --threads "$THREADS" --db "$DB_PATH")
  if [[ "$USE_FTP" == "1" ]]; then
    args+=(--use-ftp)
  fi

  log "Downloading NCBI taxonomy:"
  log "  $KRAKEN2_BUILD_BIN ${args[*]}"
  run_kraken2_build "${args[@]}"
  local status=$?
  if [[ "$status" != "0" ]]; then
    fail "Taxonomy download failed."
    return 1
  fi
}

download_libraries() {
  for lib in "${LIBRARIES[@]}"; do
    [[ -z "$lib" ]] && continue

    local args=(--download-library "$lib" --threads "$THREADS" --db "$DB_PATH")
    if [[ "$USE_FTP" == "1" ]]; then
      args+=(--use-ftp)
    fi

    log "Downloading Kraken 2 library: $lib"
    log "  $KRAKEN2_BUILD_BIN ${args[*]}"
    run_kraken2_build "${args[@]}"
    local status=$?
    if [[ "$status" != "0" ]]; then
      fail "Library download failed: $lib"
      return 1
    fi
  done
}

build_kraken_db() {
  if [[ "$BUILD_KRAKEN" != "1" ]]; then
    log "Skipping kraken2-build --build."
    return 0
  fi

  local args=(
    --build
    --threads "$THREADS"
    --db "$DB_PATH"
    --kmer-len "$KMER_LEN"
    --minimizer-len "$MINIMIZER_LEN"
    --minimizer-spaces "$MINIMIZER_SPACES"
  )

  log "Building Kraken 2 database:"
  log "  $KRAKEN2_BUILD_BIN ${args[*]}"
  run_kraken2_build "${args[@]}"
  local status=$?
  if [[ "$status" != "0" ]]; then
    fail "Kraken 2 database build failed."
    return 1
  fi
}

build_bracken_db() {
  if [[ "$BUILD_BRACKEN" != "1" ]]; then
    log "Skipping Bracken build."
    return 0
  fi

  log "Building Bracken database:"
  log "  $BRACKEN_BUILD_BIN -d $DB_PATH -t $THREADS -k $KMER_LEN -l $BRACKEN_READ_LEN"
  run_bracken_build \
    -d "$DB_PATH" \
    -t "$THREADS" \
    -k "$KMER_LEN" \
    -l "$BRACKEN_READ_LEN"

  local status=$?
  if [[ "$status" != "0" ]]; then
    fail "Bracken build failed."
    return 1
  fi
}

write_manifest() {
  local manifest="$DB_PATH/MTD_custom_micro_manifest.txt"

  mkdir -p "$DB_PATH"

  {
    echo "MTD Explorer custom microbiome database manifest"
    echo "================================================"
    echo "Date: $(date -Iseconds)"
    echo "Script: Create_custom_micro.sh"
    echo "Script_version: $SCRIPT_VERSION"
    echo "Database_path: $DB_PATH"
    echo "Threads: $THREADS"
    echo "Kraken2_kmer_len: $KMER_LEN"
    echo "Kraken2_minimizer_len: $MINIMIZER_LEN"
    echo "Kraken2_minimizer_spaces: $MINIMIZER_SPACES"
    echo "Bracken_read_len: $BRACKEN_READ_LEN"
    echo "Download_taxonomy: $DOWNLOAD_TAXONOMY"
    echo "Use_ftp: $USE_FTP"
    echo
    echo "Public_libraries:"
    if [[ "${#LIBRARIES[@]}" -eq 0 ]]; then
      echo "  none"
    else
      for lib in "${LIBRARIES[@]}"; do
        echo "  $lib"
      done
    fi
    echo
    echo "Custom_FASTA_direct:"
    if [[ "${#FASTA_FILES[@]}" -eq 0 ]]; then
      echo "  none"
    else
      for fasta in "${FASTA_FILES[@]}"; do
        echo "  $fasta"
      done
    fi
    echo
    echo "Custom_FASTA_with_taxid:"
    if [[ "${#FASTA_WITH_TAXID[@]}" -eq 0 ]]; then
      echo "  none"
    else
      for item in "${FASTA_WITH_TAXID[@]}"; do
        echo "  $item"
      done
    fi
  } > "$manifest"

  log "Manifest written:"
  log "  $manifest"
}

validate_database() {
  log "Validating database:"
  log "  $DB_PATH"

  local missing=0

  for f in hash.k2d opts.k2d taxo.k2d; do
    if [[ -s "$DB_PATH/$f" ]]; then
      log "Found $f"
    else
      warn "Missing $f"
      missing=1
    fi
  done

  local bracken_file="$DB_PATH/database${BRACKEN_READ_LEN}mers.kmer_distrib"
  if [[ -s "$bracken_file" ]]; then
    log "Found Bracken file: $(basename "$bracken_file")"
  else
    warn "Missing Bracken file: $bracken_file"
    if [[ "$BUILD_BRACKEN" == "1" ]]; then
      missing=1
    fi
  fi

  if [[ -x "$KRAKEN2_INSPECT_BIN" && -s "$DB_PATH/hash.k2d" ]]; then
    log "kraken2-inspect preview:"
    run_kraken2_inspect --db "$DB_PATH" | head -20
  else
    warn "kraken2-inspect unavailable or database hash missing; skipping taxonomy preview."
  fi

  if [[ "$missing" == "1" ]]; then
    warn "Validation found missing files."
    return 1
  fi

  log "Validation complete."
  return 0
}

main() {
  parse_args "$@"
  local parse_status=$?

  if [[ "$parse_status" == "2" ]]; then
    return 0
  elif [[ "$parse_status" != "0" ]]; then
    return 1
  fi

  for list_file in "${FASTA_LISTS[@]}"; do
    read_fasta_list "$list_file"
    local list_status=$?
    if [[ "$list_status" != "0" ]]; then
      return 1
    fi
  done

  resolve_db_path || return 1
  resolve_dedicated_kraken_runtime || return 1

  echo "============================================================"
  echo "Create_custom_micro.sh"
  echo "============================================================"
  echo "Database path:       $DB_PATH"
  echo "Threads:             $THREADS"
  echo "Public libraries:    ${LIBRARIES[*]:-none}"
  echo "Custom FASTA:        ${#FASTA_FILES[@]}"
  echo "Custom FASTA+TaxID:  ${#FASTA_WITH_TAXID[@]}"
  echo "Kraken k-mer length: $KMER_LEN"
  echo "Bracken read length: $BRACKEN_READ_LEN"
  echo "Rebuild Bracken only: $REBUILD_BRACKEN_ONLY"
  echo "============================================================"

  if [[ "$VALIDATE_ONLY" == "1" ]]; then
    validate_database
    return $?
  fi

  if [[ "${#LIBRARIES[@]}" -eq 0 && "${#FASTA_FILES[@]}" -eq 0 && "${#FASTA_WITH_TAXID[@]}" -eq 0 ]]; then
    if [[ "$REBUILD_BRACKEN_ONLY" == "1" ]]; then
      if [[ ! -s "$DB_PATH/hash.k2d" || ! -s "$DB_PATH/opts.k2d" || ! -s "$DB_PATH/taxo.k2d" ]]; then
        fail "--rebuild-bracken-only requires an existing Kraken 2 database with hash.k2d, opts.k2d, and taxo.k2d."
        return 1
      fi
      log "No microbial input provided, but --rebuild-bracken-only was requested."
      log "Using existing Kraken 2 database to rebuild Bracken files only."
    else
      fail "No microbial input provided. Use --libraries, --library, --add-fasta, --add-fasta-with-taxid, or --rebuild-bracken-only."
      return 1
    fi
  fi

  check_dependencies || return 1
  safe_clean_db || return 1

  mkdir -p "$DB_PATH"

  download_taxonomy || return 1
  download_libraries || return 1
  add_custom_fasta_files || return 1
  build_kraken_db || return 1
  build_bracken_db || return 1
  write_manifest
  validate_database
  local val_status=$?

  echo
  echo "============================================================"
  echo "MTD Explorer usage"
  echo "============================================================"
  echo "Use this database with:"
  echo
  echo "bash ~/MTD/MTD_explorer.sh \\"
  echo "  -i /path/to/samplesheet.csv \\"
  echo "  -o /path/to/output_directory \\"
  echo "  -h <host_taxon_id> \\"
  echo "  --kraken-micro-db \"$DB_PATH\" \\"
  echo "  --bracken-read-len $BRACKEN_READ_LEN"
  echo

  return $val_status
}

main "$@"
