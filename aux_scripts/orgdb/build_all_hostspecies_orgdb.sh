#!/usr/bin/env bash
# MTD Explorer — batch OrgDb builder
# Version 1.0.0
#
# Builds or validates one OrgDb package per HostSpecies.csv row by reusing:
#   aux_scripts/orgdb/install_gold_orgdb_from_hostspecies.sh
#
# Only the GTF and protein FASTA are required. The script runs eggNOG-mapper,
# builds/installs the OrgDb, detects the exact installed package name, validates
# it, and atomically updates the OrgDb column in HostSpecies.csv.
#
# It is intentionally serial: each eggNOG-mapper process uses --threads, and
# running several species at once would compete heavily for RAM, disk and the
# shared eggNOG database.

set -Eeuo pipefail
shopt -s nullglob

SCRIPT_VERSION="1.0.2"

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
MTDIR_DEFAULT="$(cd -- "$SCRIPT_DIR/../.." 2>/dev/null && pwd -P || true)"
if [[ ! -f "$MTDIR_DEFAULT/HostSpecies.csv" ]]; then
    MTDIR_DEFAULT="$(pwd -P)"
fi

MTDIR="${MTDIR:-$MTDIR_DEFAULT}"
HOSTSPECIES="$MTDIR/HostSpecies.csv"
SPECIES_ROOT="$HOME/org_species"
DOWNLOAD_REPORT=""
HELPER="$MTDIR/aux_scripts/orgdb/install_gold_orgdb_from_hostspecies.sh"
LIB="$MTDIR/custom_R_libs"
BUILD_ROOT="$MTDIR/build/orgdb_gold"
LOG_ROOT=""
REPORT=""
EGGNOG_DB=""
CONDA_R_ENV="${MTD_ORGDB_ENV:-MTD_orgdb}"
THREADS="$(nproc)"
ORGDB_VERSION="0.3.0"

FORCE=0
DRY_RUN=0
MAX_SPECIES=0
START_AT=""
STOP_AFTER=""
UPDATE_CSV=1
SKIP_IF_COMPATIBLE=1
DOWNLOAD_MISSING=1

declare -a SELECTED_TAXIDS=()

usage() {
    cat <<EOF

MTD Explorer — build all host OrgDb packages

Usage:
  bash $0 [options]

Main options:
  --hostspecies FILE       Curated HostSpecies.csv
                           [default: $HOSTSPECIES]
  --species-root DIR       Reference folders created by the downloader
                           [default: $SPECIES_ROOT]
  --download-report FILE   download_report.tsv used to recover exact folder names
                           [default: <species-root>/download_report.tsv]
  --threads N              Threads used by eggNOG-mapper [default: $THREADS]
  --taxid ID               Process only one TaxID; may be repeated
  --max-species N          Stop after N selected species (0 = no limit)
  --start-at ID            Begin when this TaxID is reached
  --stop-after ID          Stop after completing this TaxID

Build behavior:
  --force                  Rebuild even when an existing OrgDb is compatible
  --no-skip-if-compatible  Do not preserve a compatible existing OrgDb
  --no-download-missing    Fail instead of downloading a missing GTF/pep from
                           the curated GTF_URL/Pep_URL columns
  --orgdb-version VERSION  Generated package version [default: $ORGDB_VERSION]
  --conda-r-env NAME       R/AnnotationForge environment [default: $CONDA_R_ENV]
  --eggnog-db DIR          eggNOG-mapper database directory
                           [default: read from offlineCachePath]
  --lib DIR                Installed OrgDb library [default: $LIB]
  --build-root DIR         Per-TaxID build directories [default: $BUILD_ROOT]

Output:
  --log-root DIR           Per-species logs
                           [default: <species-root>/orgdb_logs]
  --report FILE            Batch report
                           [default: <species-root>/orgdb_report.tsv]
  --no-update-csv          Build/validate packages but do not change HostSpecies.csv
  --dry-run                Print the queue and commands without running them
  --version                Print version
  -h, --help               Show this help

Examples:

  Test one species:
    bash $0 --taxid 10090 --threads 20

  Test three species:
    bash $0 --taxid 10090 --taxid 8839 --taxid 6526 --threads 20

  Resume/process the complete CSV:
    bash $0 --threads 20

  Rebuild one package:
    bash $0 --taxid 6526 --threads 20 --force

EOF
}

die() {
    echo "[ERROR] $*" >&2
    exit 1
}

warn() {
    echo "[WARN] $*" >&2
}

info() {
    echo "[INFO] $*"
}

command_exists() {
    command -v "$1" >/dev/null 2>&1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --hostspecies)
            HOSTSPECIES="$2"
            shift 2
            ;;
        --species-root)
            SPECIES_ROOT="$2"
            shift 2
            ;;
        --download-report)
            DOWNLOAD_REPORT="$2"
            shift 2
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        --taxid)
            SELECTED_TAXIDS+=("$2")
            shift 2
            ;;
        --max-species)
            MAX_SPECIES="$2"
            shift 2
            ;;
        --start-at)
            START_AT="$2"
            shift 2
            ;;
        --stop-after)
            STOP_AFTER="$2"
            shift 2
            ;;
        --force)
            FORCE=1
            shift
            ;;
        --no-skip-if-compatible)
            SKIP_IF_COMPATIBLE=0
            shift
            ;;
        --no-download-missing)
            DOWNLOAD_MISSING=0
            shift
            ;;
        --orgdb-version)
            ORGDB_VERSION="$2"
            shift 2
            ;;
        --conda-r-env)
            CONDA_R_ENV="$2"
            shift 2
            ;;
        --eggnog-db)
            EGGNOG_DB="$2"
            shift 2
            ;;
        --lib)
            LIB="$2"
            shift 2
            ;;
        --build-root)
            BUILD_ROOT="$2"
            shift 2
            ;;
        --log-root)
            LOG_ROOT="$2"
            shift 2
            ;;
        --report)
            REPORT="$2"
            shift 2
            ;;
        --no-update-csv)
            UPDATE_CSV=0
            shift
            ;;
        --dry-run)
            DRY_RUN=1
            shift
            ;;
        --version)
            echo "$SCRIPT_VERSION"
            exit 0
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            die "Unknown option: $1"
            ;;
    esac
done

[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || die "--threads must be a positive integer"
[[ "$MAX_SPECIES" =~ ^[0-9]+$ ]] || die "--max-species must be zero or a positive integer"

HOSTSPECIES="$(readlink -f "$HOSTSPECIES")"
SPECIES_ROOT="$(readlink -m "$SPECIES_ROOT")"
HELPER="$(readlink -m "$HELPER")"
LIB="$(readlink -m "$LIB")"
BUILD_ROOT="$(readlink -m "$BUILD_ROOT")"

if [[ -z "$DOWNLOAD_REPORT" ]]; then
    DOWNLOAD_REPORT="$SPECIES_ROOT/download_report.tsv"
fi
DOWNLOAD_REPORT="$(readlink -m "$DOWNLOAD_REPORT")"

if [[ -z "$LOG_ROOT" ]]; then
    LOG_ROOT="$SPECIES_ROOT/orgdb_logs"
fi
LOG_ROOT="$(readlink -m "$LOG_ROOT")"

if [[ -z "$REPORT" ]]; then
    REPORT="$SPECIES_ROOT/orgdb_report.tsv"
fi
REPORT="$(readlink -m "$REPORT")"

[[ -s "$HOSTSPECIES" ]] || die "HostSpecies.csv not found or empty: $HOSTSPECIES"
[[ -s "$HELPER" ]] || die "OrgDb helper not found or empty: $HELPER"
command_exists python3 || die "python3 is required"
command_exists gzip || die "gzip is required"
command_exists flock || die "flock is required"

if [[ -z "$EGGNOG_DB" ]]; then
    OFFLINE_CACHE_FILE="$MTDIR/offlineCachePath"
    [[ -s "$OFFLINE_CACHE_FILE" ]] || die \
        "offlineCachePath was not found. Provide --eggnog-db explicitly."
    IFS= read -r OFFLINE_CACHE_ROOT < "$OFFLINE_CACHE_FILE"
    OFFLINE_CACHE_ROOT="${OFFLINE_CACHE_ROOT%$'\r'}"
    [[ -n "$OFFLINE_CACHE_ROOT" ]] || die "offlineCachePath is empty"
    EGGNOG_DB="$OFFLINE_CACHE_ROOT/eggNOG/emapperdb-5.0.2"
fi
EGGNOG_DB="$(readlink -m "$EGGNOG_DB")"

mkdir -p "$SPECIES_ROOT" "$LOG_ROOT" "$LIB" "$BUILD_ROOT" "$(dirname "$REPORT")"

# Prevent two batch runs from updating the CSV simultaneously.
LOCK_FILE="$SPECIES_ROOT/.build_all_hostspecies_orgdb.lock"
exec 9>"$LOCK_FILE"
if ! flock -n 9; then
    die "Another OrgDb batch appears to be running: $LOCK_FILE"
fi

# Validate required columns before any expensive work.
python3 - "$HOSTSPECIES" <<'PY'
import csv
import sys

path = sys.argv[1]
required = {
    "Taxon_ID",
    "Scientific_name",
    "OrgDb",
    "GTF_URL",
    "Pep_URL",
    "Reference_status",
}
with open(path, encoding="utf-8-sig", newline="") as handle:
    reader = csv.DictReader(handle)
    fields = set(reader.fieldnames or [])
missing = required - fields
if missing:
    raise SystemExit(
        "[ERROR] HostSpecies.csv is missing required column(s): "
        + ", ".join(sorted(missing))
    )
PY

# Validate selected TaxIDs early.
if ((${#SELECTED_TAXIDS[@]} > 0)); then
    python3 - "$HOSTSPECIES" "${SELECTED_TAXIDS[@]}" <<'PY'
import csv
import sys

path = sys.argv[1]
requested = sys.argv[2:]
with open(path, encoding="utf-8-sig", newline="") as handle:
    available = {
        str(row.get("Taxon_ID", "")).strip()
        for row in csv.DictReader(handle)
    }
missing = [taxid for taxid in requested if taxid not in available]
if missing:
    raise SystemExit(
        "[ERROR] Requested TaxID(s) not found in HostSpecies.csv: "
        + ", ".join(missing)
    )
PY
fi

# Check the shared eggNOG database now rather than failing after preparing a species.
for required_file in eggnog.db eggnog_proteins.dmnd eggnog.taxa.db; do
    [[ -s "$EGGNOG_DB/$required_file" ]] || die \
        "eggNOG database is incomplete; missing: $EGGNOG_DB/$required_file"
done

# Load conda so package validation can run after each helper call.
if ! command_exists conda; then
    if [[ -s "$MTDIR/condaPath" ]]; then
        CONDA_BASE="$(head -n 1 "$MTDIR/condaPath")"
        # shellcheck source=/dev/null
        source "$CONDA_BASE/etc/profile.d/conda.sh"
    elif [[ -s "$HOME/miniconda3/etc/profile.d/conda.sh" ]]; then
        # shellcheck source=/dev/null
        source "$HOME/miniconda3/etc/profile.d/conda.sh"
    fi
fi
command_exists conda || die "conda was not found"

if ! conda run -n "$CONDA_R_ENV" Rscript -e \
    'stopifnot(requireNamespace("AnnotationDbi", quietly=TRUE),
               requireNamespace("AnnotationForge", quietly=TRUE));
     cat("[OK] OrgDb R environment\n")' >/dev/null; then
    die "Conda environment '$CONDA_R_ENV' is unavailable or incomplete"
fi

# Create a single backup. Progress is then written atomically after each species,
# allowing the batch to resume safely after interruption.
BACKUP=""
if [[ "$UPDATE_CSV" == "1" && "$DRY_RUN" != "1" ]]; then
    TIMESTAMP="$(date +%Y%m%d-%H%M%S)"
    BACKUP="${HOSTSPECIES}.bak.orgdb-${TIMESTAMP}"
    cp -a -- "$HOSTSPECIES" "$BACKUP"
fi

# Create queue. The temporary package scientific name is kept binomial for
# AnnotationForge. Trinomials/hybrids combine all epithets after the genus,
# preventing package-name collisions.
QUEUE_FILE="$(mktemp)"
trap 'rm -f "$QUEUE_FILE"' EXIT

python3 - \
    "$HOSTSPECIES" \
    "$DOWNLOAD_REPORT" \
    "$QUEUE_FILE" \
    "${SELECTED_TAXIDS[@]}" <<'PY'
import csv
import os
import re
import sys

host_path, report_path, queue_path, *selected = sys.argv[1:]
selected_set = set(selected)

folders = {}
if os.path.isfile(report_path):
    with open(report_path, encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            taxid = str(row.get("taxid", "")).strip()
            folder = str(row.get("folder", "")).strip()
            if taxid and folder:
                folders[taxid] = folder

def fallback_folder(scientific_name):
    parts = scientific_name.split()
    if len(parts) < 2:
        return re.sub(r"[^A-Za-z0-9._-]+", "_", scientific_name)
    rest = "_".join(parts[1:])
    return f"{parts[0][0]}.{rest}"

def package_scientific_name(scientific_name):
    parts = scientific_name.split()
    if len(parts) < 2:
        raise ValueError(
            f"Scientific_name is not at least binomial: {scientific_name}"
        )
    if len(parts) == 2:
        return scientific_name

    def token(value):
        pieces = re.findall(r"[A-Za-z0-9]+", value)
        return "".join(
            piece[:1].upper() + piece[1:]
            for piece in pieces
        )

    combined = parts[1] + "".join(token(part) for part in parts[2:])
    return f"{parts[0]} {combined}"

with open(host_path, encoding="utf-8-sig", newline="") as handle:
    rows = list(csv.DictReader(handle))

with open(queue_path, "w", encoding="utf-8", newline="") as handle:
    # ASCII Unit Separator preserves empty fields when read by Bash.
    writer = csv.writer(handle, delimiter="\x1f", lineterminator="\n")
    writer.writerow([
        "taxid",
        "scientific_name",
        "package_scientific_name",
        "folder",
        "old_orgdb",
        "gtf_url",
        "pep_url",
        "reference_status",
    ])

    for row in rows:
        taxid = str(row.get("Taxon_ID", "")).strip()
        if selected_set and taxid not in selected_set:
            continue

        scientific_name = str(row.get("Scientific_name", "")).strip()
        writer.writerow([
            taxid,
            scientific_name,
            package_scientific_name(scientific_name),
            folders.get(taxid) or fallback_folder(scientific_name),
            str(row.get("OrgDb", "") or "").strip(),
            str(row.get("GTF_URL", "") or "").strip(),
            str(row.get("Pep_URL", "") or "").strip(),
            str(row.get("Reference_status", "") or "").strip(),
        ])
PY

validate_gzip_or_plain() {
    local path="$1"
    [[ -s "$path" ]] || return 1
    case "$path" in
        *.gz)
            gzip -t -- "$path" >/dev/null 2>&1
            ;;
        *)
            return 0
            ;;
    esac
}

validate_gtf() {
    local path="$1"
    validate_gzip_or_plain "$path" || return 1

    # The annotation only needs one valid non-comment record. Disable pipefail
    # inside this inspection pipeline because awk intentionally exits early;
    # otherwise gzip receives SIGPIPE and a valid file is reported as failed.
    (
        set +o pipefail
        case "$path" in
            *.gz)
                gzip -dc -- "$path"
                ;;
            *)
                cat -- "$path"
                ;;
        esac | awk -F'\t' '
            /^#/ || NF == 0 { next }
            NF >= 9 { found=1; exit }
            END { exit !found }
        '
    )
}

validate_protein_fasta() {
    local path="$1"
    validate_gzip_or_plain "$path" || return 1

    # grep -m1 exits as soon as it sees the first FASTA header. With global
    # pipefail that early exit makes gzip fail with SIGPIPE, so inspect in a
    # subshell where the grep result is the pipeline result.
    (
        set +o pipefail
        case "$path" in
            *.gz)
                gzip -dc -- "$path"
                ;;
            *)
                cat -- "$path"
                ;;
        esac | grep -qm1 '^>'
    )
}

download_file() {
    local url="$1"
    local destination="$2"
    local partial="${destination}.part"

    mkdir -p "$(dirname "$destination")"

    if [[ -s "$destination" ]] && validate_gzip_or_plain "$destination"; then
        return 0
    fi

    rm -f -- "$destination"

    info "Downloading: $url" >&2
    info "Destination: $destination" >&2

    if command_exists aria2c; then
        aria2c \
            --continue=true \
            --max-connection-per-server=8 \
            --split=8 \
            --min-split-size=10M \
            --retry-wait=20 \
            --max-tries=50 \
            --timeout=60 \
            --connect-timeout=30 \
            --allow-overwrite=true \
            --auto-file-renaming=false \
            --dir "$(dirname "$partial")" \
            --out "$(basename "$partial")" \
            "$url" >&2
    elif command_exists curl; then
        curl \
            --fail \
            --location \
            --connect-timeout 30 \
            --retry 20 \
            --retry-delay 20 \
            --retry-connrefused \
            --continue-at - \
            --output "$partial" \
            "$url" >&2
    elif command_exists wget; then
        wget \
            --continue \
            --tries=50 \
            --waitretry=20 \
            --timeout=60 \
            --read-timeout=60 \
            --retry-connrefused \
            --output-document="$partial" \
            "$url" >&2
    else
        die "aria2c, curl and wget are unavailable"
    fi

    validate_gzip_or_plain "$partial" || {
        rm -f -- "$partial"
        return 1
    }

    mv -f -- "$partial" "$destination"
}

resolve_reference() {
    local kind="$1"
    local stable_path="$2"
    local url="$3"
    local reference_dir="$4"

    if [[ -e "$stable_path" || -L "$stable_path" ]]; then
        local resolved
        resolved="$(readlink -f "$stable_path" 2>/dev/null || true)"
        if [[ -n "$resolved" && -s "$resolved" ]]; then
            case "$kind" in
                gtf)
                    if validate_gtf "$resolved"; then
                        printf '%s\n' "$resolved"
                        return 0
                    fi
                    ;;
                protein)
                    if validate_protein_fasta "$resolved"; then
                        printf '%s\n' "$resolved"
                        return 0
                    fi
                    ;;
            esac
        fi
        if [[ -n "$resolved" ]]; then
            warn "Existing $kind reference failed validation: $resolved"
        else
            warn "Broken $kind reference link: $stable_path"
        fi
        rm -f -- "$stable_path"
    fi

    [[ "$DOWNLOAD_MISSING" == "1" ]] || return 1
    [[ -n "$url" ]] || return 1

    local clean_url="${url%%\?*}"
    local filename
    filename="$(basename "$clean_url")"
    [[ -n "$filename" && "$filename" != "." ]] || return 1

    local destination="$reference_dir/$filename"
    download_file "$url" "$destination" || return 1

    case "$kind" in
        gtf)
            validate_gtf "$destination" || return 1
            ;;
        protein)
            validate_protein_fasta "$destination" || return 1
            ;;
    esac

    ln -sfn -- "$filename" "$stable_path"
    printf '%s\n' "$(readlink -f "$stable_path")"
}

make_package_csv() {
    local taxid="$1"
    local package_scientific_name="$2"
    local destination="$3"

    python3 - \
        "$HOSTSPECIES" \
        "$taxid" \
        "$package_scientific_name" \
        "$destination" <<'PY'
import csv
import os
import sys

source, taxid, package_name, destination = sys.argv[1:]
with open(source, encoding="utf-8-sig", newline="") as handle:
    reader = csv.DictReader(handle)
    fieldnames = reader.fieldnames
    rows = list(reader)

matches = 0
for row in rows:
    if str(row.get("Taxon_ID", "")).strip() == taxid:
        row["Scientific_name"] = package_name
        matches += 1

if matches != 1:
    raise SystemExit(
        f"[ERROR] Expected one CSV row for TaxID {taxid}, found {matches}"
    )

temporary = destination + ".tmp"
with open(temporary, "w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=fieldnames,
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(rows)
os.replace(temporary, destination)
PY
}

read_current_orgdb() {
    local taxid="$1"
    python3 - "$HOSTSPECIES" "$taxid" <<'PY'
import csv
import sys

path, taxid = sys.argv[1:]
with open(path, encoding="utf-8-sig", newline="") as handle:
    for row in csv.DictReader(handle):
        if str(row.get("Taxon_ID", "")).strip() == taxid:
            print(str(row.get("OrgDb", "") or "").strip())
            raise SystemExit(0)
raise SystemExit(f"[ERROR] TaxID {taxid} not found")
PY
}

update_orgdb_csv() {
    local taxid="$1"
    local package="$2"

    python3 - "$HOSTSPECIES" "$taxid" "$package" <<'PY'
import csv
import os
import stat
import sys

path, taxid, package = sys.argv[1:]

if not package.startswith("org.") or not package.endswith(".db"):
    raise SystemExit(f"[ERROR] Refusing invalid OrgDb package name: {package}")

with open(path, encoding="utf-8-sig", newline="") as handle:
    reader = csv.DictReader(handle)
    fieldnames = reader.fieldnames
    rows = list(reader)

matches = 0
old = ""
for row in rows:
    if str(row.get("Taxon_ID", "")).strip() == taxid:
        old = str(row.get("OrgDb", "") or "").strip()
        row["OrgDb"] = package
        matches += 1

if matches != 1:
    raise SystemExit(
        f"[ERROR] Expected one CSV row for TaxID {taxid}, found {matches}"
    )

temporary = path + ".tmp"
with open(temporary, "w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=fieldnames,
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(rows)

mode = stat.S_IMODE(os.stat(path).st_mode)
os.chmod(temporary, mode)
os.replace(temporary, path)

print(f"[CSV] TaxID {taxid}: {old or '<empty>'} -> {package}")
PY
}

detect_package_name() {
    local log_file="$1"
    local build_dir="$2"
    local old_orgdb="$3"
    local detected=""

    # Preferred source: exact message emitted by the R builder.
    detected="$(
        sed -n \
            's/.*\[OK\] Installed package: \([^[:space:]]*\) into .*/\1/p' \
            "$log_file" | tail -n 1
    )"

    # Fallback: source package directory inside this TaxID-specific build dir.
    if [[ -z "$detected" ]]; then
        detected="$(
            find "$build_dir" \
                -mindepth 1 \
                -maxdepth 1 \
                -type d \
                -name 'org.*.db' \
                -printf '%T@\t%f\n' 2>/dev/null \
            | sort -nr \
            | awk -F'\t' 'NR==1 {print $2}'
        )"
    fi

    # Compatible official/existing package: helper exits without rebuilding.
    if [[ -z "$detected" ]]; then
        detected="$old_orgdb"
    fi

    printf '%s\n' "$detected"
}

validate_installed_package() {
    local package="$1"
    local expected_taxid="$2"
    local validation_file="$3"

    ORGDB_PKG="$package" \
    ORGDB_LIB="$LIB" \
    EXPECTED_TAXID="$expected_taxid" \
    conda run -n "$CONDA_R_ENV" Rscript - <<'RS' >"$validation_file" 2>&1
pkg <- Sys.getenv("ORGDB_PKG")
lib <- Sys.getenv("ORGDB_LIB")
expected_taxid <- Sys.getenv("EXPECTED_TAXID")

.libPaths(unique(c(lib, .libPaths())))

if (!requireNamespace(pkg, quietly = TRUE)) {
  stop("OrgDb package is not loadable: ", pkg)
}

suppressPackageStartupMessages(library(AnnotationDbi))
suppressPackageStartupMessages(library(pkg, character.only = TRUE))

db <- get(pkg)
kt <- keytypes(db)
cols <- columns(db)
gid_keys <- tryCatch(keys(db, keytype = "GID"), error = function(e) character())
symbol_keys <- tryCatch(keys(db, keytype = "SYMBOL"), error = function(e) character())
md <- tryCatch(AnnotationDbi::metadata(db), error = function(e) data.frame())

tax_values <- character()
if (nrow(md) > 0 && all(c("name", "value") %in% names(md))) {
  tax_values <- as.character(
    md$value[toupper(md$name) %in% c("TAXID", "TAX_ID", "TAXONOMY ID")]
  )
  tax_values <- unique(tax_values[nzchar(tax_values)])
}

cat("package\t", pkg, "\n", sep = "")
cat("path\t", find.package(pkg), "\n", sep = "")
cat("version\t", as.character(packageVersion(pkg)), "\n", sep = "")
cat("gid_keys\t", length(gid_keys), "\n", sep = "")
cat("symbol_keys\t", length(symbol_keys), "\n", sep = "")
cat("taxid_metadata\t", paste(tax_values, collapse = ","), "\n", sep = "")
cat("keytypes\t", paste(kt, collapse = ","), "\n", sep = "")
cat("columns\t", paste(cols, collapse = ","), "\n", sep = "")

if (length(gid_keys) == 0 && length(symbol_keys) == 0) {
  stop("OrgDb contains neither GID nor SYMBOL keys")
}

if (length(tax_values) > 0 && !expected_taxid %in% tax_values) {
  warning(
    "Package taxonomy metadata does not contain expected TaxID ",
    expected_taxid,
    ": ",
    paste(tax_values, collapse = ",")
  )
}

cat("status\tPASS\n", sep = "")
RS
}

append_report() {
    local taxid="$1"
    local scientific_name="$2"
    local folder="$3"
    local old_orgdb="$4"
    local final_orgdb="$5"
    local status="$6"
    local elapsed="$7"
    local log_file="$8"
    local message="$9"

    python3 - \
        "$REPORT" \
        "$taxid" \
        "$scientific_name" \
        "$folder" \
        "$old_orgdb" \
        "$final_orgdb" \
        "$status" \
        "$elapsed" \
        "$log_file" \
        "$message" <<'PY'
import csv
import os
import sys

(
    path,
    taxid,
    scientific_name,
    folder,
    old_orgdb,
    final_orgdb,
    status,
    elapsed,
    log_file,
    message,
) = sys.argv[1:]

fields = [
    "taxid",
    "scientific_name",
    "folder",
    "old_orgdb",
    "final_orgdb",
    "status",
    "elapsed_seconds",
    "log_file",
    "message",
]

exists = os.path.isfile(path) and os.path.getsize(path) > 0
with open(path, "a", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n",
    )
    if not exists:
        writer.writeheader()
    writer.writerow({
        "taxid": taxid,
        "scientific_name": scientific_name,
        "folder": folder,
        "old_orgdb": old_orgdb,
        "final_orgdb": final_orgdb,
        "status": status,
        "elapsed_seconds": elapsed,
        "log_file": log_file,
        "message": message,
    })
PY
}

echo "============================================================"
echo "MTD Explorer — batch OrgDb builder"
echo "Version:          $SCRIPT_VERSION"
echo "HostSpecies.csv:  $HOSTSPECIES"
echo "Species root:     $SPECIES_ROOT"
echo "Download report: $DOWNLOAD_REPORT"
echo "eggNOG DB:        $EGGNOG_DB"
echo "OrgDb library:    $LIB"
echo "Build root:       $BUILD_ROOT"
echo "Threads/species:  $THREADS"
echo "R environment:    $CONDA_R_ENV"
echo "Update CSV:       $([[ "$UPDATE_CSV" == "1" ]] && echo yes || echo no)"
echo "Force:            $([[ "$FORCE" == "1" ]] && echo yes || echo no)"
echo "Dry run:          $([[ "$DRY_RUN" == "1" ]] && echo yes || echo no)"
if [[ -n "$BACKUP" ]]; then
    echo "CSV backup:       $BACKUP"
fi
echo "============================================================"

processed=0
completed=0
failed=0
skipped_before_start=0
start_reached=0

if [[ -z "$START_AT" ]]; then
    start_reached=1
fi

# Read the queue using ASCII Unit Separator. Unlike tab, this delimiter is
# not IFS whitespace, so an empty OrgDb field does not shift all later columns.
while IFS=$'\x1f' read -r \
    taxid scientific_name package_scientific_name folder old_orgdb \
    gtf_url pep_url reference_status; do

    [[ "$taxid" == "taxid" ]] && continue

    # Defensive validation against queue-column shifts or malformed rows.
    if [[ "$gtf_url" != http://* && "$gtf_url" != https://* ]]; then
        die "Malformed queue row for TaxID $taxid: invalid GTF_URL '$gtf_url'"
    fi
    if [[ "$pep_url" != http://* && "$pep_url" != https://* ]]; then
        die "Malformed queue row for TaxID $taxid: invalid Pep_URL '$pep_url'"
    fi
    if [[ "$reference_status" != COMPLETE* ]]; then
        warn "TaxID $taxid has Reference_status '$reference_status'"
    fi

    if [[ "$start_reached" != "1" ]]; then
        if [[ "$taxid" == "$START_AT" ]]; then
            start_reached=1
        else
            ((skipped_before_start += 1))
            continue
        fi
    fi

    if [[ "$MAX_SPECIES" -gt 0 && "$processed" -ge "$MAX_SPECIES" ]]; then
        break
    fi

    ((processed += 1))

    species_dir="$SPECIES_ROOT/$folder"
    reference_dir="$species_dir/references"
    orgdb_dir="$species_dir/orgdb"
    build_dir="$BUILD_ROOT/$taxid"
    log_file="$LOG_ROOT/${taxid}_${folder}.log"
    validation_file="$orgdb_dir/validation.tsv"
    package_file="$orgdb_dir/package_name.txt"
    status_file="$orgdb_dir/status.tsv"
    package_csv="$build_dir/HostSpecies.package_input.csv"

    mkdir -p "$reference_dir" "$orgdb_dir" "$build_dir" "$LOG_ROOT"

    echo
    echo "============================================================"
    echo "[$processed] TaxID $taxid — $scientific_name"
    echo "Folder:          $folder"
    echo "Package species: $package_scientific_name"
    echo "Current OrgDb:   ${old_orgdb:-<empty>}"
    echo "============================================================"

    if [[ "$DRY_RUN" == "1" ]]; then
        echo "[DRY-RUN] GTF link: $reference_dir/annotation.gtf.gz"
        echo "[DRY-RUN] PEP link: $reference_dir/proteins.pep.all.fa.gz"
        echo "[DRY-RUN] Build dir: $build_dir"
        echo "[DRY-RUN] Log: $log_file"
        echo "[DRY-RUN] Helper: $HELPER"
        if [[ "$taxid" == "$STOP_AFTER" ]]; then
            break
        fi
        continue
    fi

    start_epoch="$(date +%s)"
    final_orgdb=""
    species_status="FAILED"
    species_message=""

    # Resolve or fetch the two necessary inputs.
    if ! gtf_path="$(
        resolve_reference \
            gtf \
            "$reference_dir/annotation.gtf.gz" \
            "$gtf_url" \
            "$reference_dir"
    )"; then
        species_message="GTF missing or invalid"
        warn "$species_message"
        elapsed="$(( $(date +%s) - start_epoch ))"
        append_report \
            "$taxid" "$scientific_name" "$folder" "$old_orgdb" "" \
            "FAILED_INPUT" "$elapsed" "$log_file" "$species_message"
        printf "status\tFAILED_INPUT\nmessage\t%s\n" "$species_message" > "$status_file"
        ((failed += 1))
        [[ "$taxid" == "$STOP_AFTER" ]] && break
        continue
    fi

    if ! pep_path="$(
        resolve_reference \
            protein \
            "$reference_dir/proteins.pep.all.fa.gz" \
            "$pep_url" \
            "$reference_dir"
    )"; then
        species_message="Protein FASTA missing or invalid"
        warn "$species_message"
        elapsed="$(( $(date +%s) - start_epoch ))"
        append_report \
            "$taxid" "$scientific_name" "$folder" "$old_orgdb" "" \
            "FAILED_INPUT" "$elapsed" "$log_file" "$species_message"
        printf "status\tFAILED_INPUT\nmessage\t%s\n" "$species_message" > "$status_file"
        ((failed += 1))
        [[ "$taxid" == "$STOP_AFTER" ]] && break
        continue
    fi

    make_package_csv "$taxid" "$package_scientific_name" "$package_csv"

    helper_args=(
        --taxid "$taxid"
        --hostspecies "$package_csv"
        --gtf "$gtf_path"
        --protein-fasta "$pep_path"
        --eggnog-db "$EGGNOG_DB"
        --lib "$LIB"
        --build-dir "$build_dir"
        --threads "$THREADS"
        --version "$ORGDB_VERSION"
        --conda-r-env "$CONDA_R_ENV"
    )

    if [[ "$SKIP_IF_COMPATIBLE" == "1" ]]; then
        helper_args+=(--skip-if-compatible)
    fi
    if [[ "$FORCE" == "1" ]]; then
        helper_args+=(--force)
    fi

    {
        echo "[RUN] $(printf '%q ' bash "$HELPER" "${helper_args[@]}")"
        echo "[INFO] Started: $(date --iso-8601=seconds)"
    } > "$log_file"

    set +e
    bash "$HELPER" "${helper_args[@]}" 2>&1 | tee -a "$log_file"
    helper_status="${PIPESTATUS[0]}"
    set -e

    echo "[INFO] Finished: $(date --iso-8601=seconds)" >> "$log_file"
    echo "[INFO] Helper exit status: $helper_status" >> "$log_file"

    if [[ "$helper_status" -ne 0 ]]; then
        species_message="OrgDb helper exited with status $helper_status"
        warn "$species_message"
        elapsed="$(( $(date +%s) - start_epoch ))"
        append_report \
            "$taxid" "$scientific_name" "$folder" "$old_orgdb" "" \
            "FAILED_BUILD" "$elapsed" "$log_file" "$species_message"
        printf "status\tFAILED_BUILD\nmessage\t%s\n" "$species_message" > "$status_file"
        ((failed += 1))
        [[ "$taxid" == "$STOP_AFTER" ]] && break
        continue
    fi

    # The actual CSV may already have been updated by a previous completed row.
    current_orgdb="$(read_current_orgdb "$taxid")"
    final_orgdb="$(
        detect_package_name "$log_file" "$build_dir" "$current_orgdb"
    )"

    if [[ -z "$final_orgdb" ]]; then
        species_message="Could not determine the final OrgDb package name"
        warn "$species_message"
        elapsed="$(( $(date +%s) - start_epoch ))"
        append_report \
            "$taxid" "$scientific_name" "$folder" "$old_orgdb" "" \
            "FAILED_DETECTION" "$elapsed" "$log_file" "$species_message"
        printf "status\tFAILED_DETECTION\nmessage\t%s\n" "$species_message" > "$status_file"
        ((failed += 1))
        [[ "$taxid" == "$STOP_AFTER" ]] && break
        continue
    fi

    if ! validate_installed_package \
        "$final_orgdb" \
        "$taxid" \
        "$validation_file"; then
        species_message="Detected package failed R validation: $final_orgdb"
        warn "$species_message"
        elapsed="$(( $(date +%s) - start_epoch ))"
        append_report \
            "$taxid" "$scientific_name" "$folder" "$old_orgdb" "$final_orgdb" \
            "FAILED_VALIDATION" "$elapsed" "$log_file" "$species_message"
        printf "status\tFAILED_VALIDATION\nmessage\t%s\n" "$species_message" > "$status_file"
        ((failed += 1))
        [[ "$taxid" == "$STOP_AFTER" ]] && break
        continue
    fi

    printf '%s\n' "$final_orgdb" > "$package_file"

    if [[ "$UPDATE_CSV" == "1" ]]; then
        update_orgdb_csv "$taxid" "$final_orgdb"
    fi

    elapsed="$(( $(date +%s) - start_epoch ))"
    species_status="COMPLETE"
    species_message="OrgDb built or validated successfully"

    {
        printf "status\t%s\n" "$species_status"
        printf "taxid\t%s\n" "$taxid"
        printf "scientific_name\t%s\n" "$scientific_name"
        printf "package_scientific_name\t%s\n" "$package_scientific_name"
        printf "orgdb\t%s\n" "$final_orgdb"
        printf "elapsed_seconds\t%s\n" "$elapsed"
        printf "log\t%s\n" "$log_file"
        printf "validation\t%s\n" "$validation_file"
    } > "$status_file"

    append_report \
        "$taxid" "$scientific_name" "$folder" "$old_orgdb" "$final_orgdb" \
        "$species_status" "$elapsed" "$log_file" "$species_message"

    echo "[OK] TaxID $taxid -> $final_orgdb"
    echo "[OK] Elapsed: ${elapsed}s"
    ((completed += 1))

    if [[ "$taxid" == "$STOP_AFTER" ]]; then
        break
    fi
done < "$QUEUE_FILE"

echo
echo "============================================================"
echo "MTD Explorer — OrgDb batch summary"
echo "Selected/processed: $processed"
echo "Completed:          $completed"
echo "Failed:             $failed"
echo "Report:             $REPORT"
echo "Logs:               $LOG_ROOT"
if [[ -n "$BACKUP" ]]; then
    echo "Original CSV backup:$BACKUP"
fi
echo "============================================================"

if [[ "$failed" -gt 0 ]]; then
    exit 2
fi
