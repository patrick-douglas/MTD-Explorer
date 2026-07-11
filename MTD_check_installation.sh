#!/usr/bin/env bash

# ==============================================================================
# MTD Explorer installation checker
# Version: 2026.07.11-r8
#
# Aligned with the installer architecture in which:
#   - only the shared microbiome Kraken2 database is installed by default;
#   - predefined human, mouse, and rhesus host databases are NOT required;
#   - custom hosts are prepared separately with Create_custom_host.sh;
#   - OrgDb construction uses the dedicated MTD_orgdb environment when present.
#
# The checker is read-only, except for its report directory and temporary files.
# ==============================================================================

set -uo pipefail

CHECKER_VERSION="2026.07.11-r8"

MTD_DIR="$HOME/MTD"
CONDA_PATH=""
OFFLINE_DIR=""
READ_LEN=75
MODE="full"
REPORT_DIR=""
STRICT=0
KEEP_TEMP=0

PASS_COUNT=0
WARN_COUNT=0
FAIL_COUNT=0
SKIP_COUNT=0
TOTAL_COUNT=0

COLOR_RESET=""
COLOR_RED=""
COLOR_GREEN=""
COLOR_YELLOW=""
COLOR_BLUE=""

RESULTS_TSV=""
FULL_LOG=""
SUMMARY_TXT=""
TMP_WORK=""

usage() {
    cat <<'USAGE'
MTD Explorer installation checker

Usage:
  bash MTD_check_installation.sh [options]

Options:
  --mtd-dir PATH       MTD repository/installation directory
                       Default: $HOME/MTD

  --conda-path PATH    Miniconda directory
                       Default: value from MTD/condaPath, then $HOME/miniconda3

  --offline-dir PATH   Persistent installation cache
                       Default: value from MTD/offlineCachePath

  -r, --read-length N  Bracken read length
                       Default: 75

  --mode MODE          quick, full, or deep
                       quick: essential runtime and database checks
                       full:  source, environments, packages, cache, and databases
                       deep:  full plus archive/integrity and heavier runtime checks

  --report-dir PATH    Output directory for reports
  --strict             Treat warnings as final failure
  --keep-temp          Keep temporary checker files
  --version            Print checker version
  -h, --help           Show this help

Examples:
  bash MTD_check_installation.sh --mode full

  bash MTD_check_installation.sh \
      --mode deep \
      --offline-dir /path/to/installer-cache

Exit status:
  0  no failures; warnings allowed unless --strict is used
  1  one or more failures, or warnings with --strict
  2  invalid checker arguments
USAGE
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --mtd-dir)
            [[ $# -ge 2 ]] || { echo "ERROR: --mtd-dir requires a value." >&2; exit 2; }
            MTD_DIR="$2"
            shift 2
            ;;
        --conda-path)
            [[ $# -ge 2 ]] || { echo "ERROR: --conda-path requires a value." >&2; exit 2; }
            CONDA_PATH="$2"
            shift 2
            ;;
        --offline-dir|-o)
            [[ $# -ge 2 ]] || { echo "ERROR: --offline-dir requires a value." >&2; exit 2; }
            OFFLINE_DIR="$2"
            shift 2
            ;;
        -r|--read-length)
            [[ $# -ge 2 ]] || { echo "ERROR: --read-length requires a value." >&2; exit 2; }
            READ_LEN="$2"
            shift 2
            ;;
        --mode)
            [[ $# -ge 2 ]] || { echo "ERROR: --mode requires a value." >&2; exit 2; }
            MODE="$2"
            shift 2
            ;;
        --report-dir)
            [[ $# -ge 2 ]] || { echo "ERROR: --report-dir requires a value." >&2; exit 2; }
            REPORT_DIR="$2"
            shift 2
            ;;
        --strict)
            STRICT=1
            shift
            ;;
        --keep-temp)
            KEEP_TEMP=1
            shift
            ;;
        --version)
            printf '%s\n' "$CHECKER_VERSION"
            exit 0
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            printf 'ERROR: Unknown option: %s\n' "$1" >&2
            usage >&2
            exit 2
            ;;
    esac
done

case "$MODE" in
    quick|full|deep) ;;
    *)
        printf 'ERROR: --mode must be quick, full, or deep. Received: %s\n' "$MODE" >&2
        exit 2
        ;;
esac

if ! [[ "$READ_LEN" =~ ^[1-9][0-9]*$ ]]; then
    printf 'ERROR: --read-length must be a positive integer.\n' >&2
    exit 2
fi

expand_path() {
    local value="$1"
    if [[ "$value" == "~" ]]; then
        value="$HOME"
    elif [[ "$value" == "~/"* ]]; then
        value="$HOME/${value#~/}"
    fi
    readlink -m -- "$value" 2>/dev/null || printf '%s\n' "$value"
}

MTD_DIR="$(expand_path "$MTD_DIR")"

if [[ -z "$CONDA_PATH" && -s "$MTD_DIR/condaPath" ]]; then
    CONDA_PATH="$(head -n 1 "$MTD_DIR/condaPath" | tr -d '\r\n')"
fi
: "${CONDA_PATH:=$HOME/miniconda3}"
CONDA_PATH="$(expand_path "$CONDA_PATH")"

if [[ -z "$OFFLINE_DIR" && -s "$MTD_DIR/offlineCachePath" ]]; then
    OFFLINE_DIR="$(head -n 1 "$MTD_DIR/offlineCachePath" | tr -d '\r\n')"
fi
if [[ -n "$OFFLINE_DIR" ]]; then
    OFFLINE_DIR="$(expand_path "$OFFLINE_DIR")"
fi

if [[ -z "$REPORT_DIR" ]]; then
    REPORT_DIR="$MTD_DIR/installation_check_$(date +%Y%m%d_%H%M%S)"
fi
REPORT_DIR="$(expand_path "$REPORT_DIR")"

mkdir -p "$REPORT_DIR" || {
    printf 'ERROR: Could not create report directory: %s\n' "$REPORT_DIR" >&2
    exit 2
}

RESULTS_TSV="$REPORT_DIR/MTD_installation_check.tsv"
FULL_LOG="$REPORT_DIR/MTD_installation_check.log"
SUMMARY_TXT="$REPORT_DIR/MTD_installation_summary.txt"
TMP_WORK="$REPORT_DIR/.tmp"

mkdir -p "$TMP_WORK"
printf 'Status\tSection\tCheck\tDetails\n' > "$RESULTS_TSV"
: > "$FULL_LOG"

cleanup() {
    if [[ "$KEEP_TEMP" -eq 0 ]]; then
        rm -rf "$TMP_WORK"
    fi
}

on_signal() {
    cleanup
    exit 130
}

trap cleanup EXIT
trap on_signal INT TERM

init_colors() {
    if [[ -t 1 && -z "${NO_COLOR:-}" ]] && command -v tput >/dev/null 2>&1; then
        COLOR_RESET="$(tput sgr0 2>/dev/null || true)"
        COLOR_RED="$(tput setaf 1 2>/dev/null || true)"
        COLOR_GREEN="$(tput setaf 2 2>/dev/null || true)"
        COLOR_YELLOW="$(tput setaf 3 2>/dev/null || true)"
        COLOR_BLUE="$(tput setaf 4 2>/dev/null || true)"
    fi
}

sanitize_field() {
    printf '%s' "$1" |
        tr '\t\r\n' '   ' |
        sed -E 's/[[:space:]]+/ /g; s/^ //; s/ $//'
}

record() {
    local status="$1"
    local section="$2"
    local check="$3"
    local details="${4:-}"
    local color=""

    TOTAL_COUNT=$((TOTAL_COUNT + 1))

    case "$status" in
        PASS)
            PASS_COUNT=$((PASS_COUNT + 1))
            color="$COLOR_GREEN"
            ;;
        WARN)
            WARN_COUNT=$((WARN_COUNT + 1))
            color="$COLOR_YELLOW"
            ;;
        FAIL)
            FAIL_COUNT=$((FAIL_COUNT + 1))
            color="$COLOR_RED"
            ;;
        SKIP)
            SKIP_COUNT=$((SKIP_COUNT + 1))
            color="$COLOR_BLUE"
            ;;
    esac

    section="$(sanitize_field "$section")"
    check="$(sanitize_field "$check")"
    details="$(sanitize_field "$details")"

    printf '%s[%s]%s %-20s | %-43s | %s\n' \
        "$color" "$status" "$COLOR_RESET" "$section" "$check" "$details"

    printf '%s\t%s\t%s\t%s\n' \
        "$status" "$section" "$check" "$details" >> "$RESULTS_TSV"
}

file_size() {
    local path="$1"
    if [[ -e "$path" ]]; then
        du -h "$path" 2>/dev/null | awk 'NR==1 {print $1}'
    fi
}

capture_command() {
    local output_file="$1"
    shift

    "$@" > "$output_file" 2>&1
}

compact_output() {
    local path="$1"
    local max_lines="${2:-8}"

    [[ -s "$path" ]] || return 0

    sed -E \
        -e 's/\x1B\[[0-9;]*[[:alpha:]]//g' \
        -e 's/\r//g' \
        "$path" |
        tail -n "$max_lines" |
        tr '\n' ' ' |
        sed -E 's/[[:space:]]+/ /g; s/^ //; s/ $//'
}

check_required_file() {
    local section="$1"
    local label="$2"
    local path="$3"

    if [[ -s "$path" ]]; then
        record PASS "$section" "$label" "$path ($(file_size "$path"))"
    elif [[ -e "$path" ]]; then
        record FAIL "$section" "$label" "exists but is empty: $path"
    else
        record FAIL "$section" "$label" "missing: $path"
    fi
}

check_optional_file() {
    local section="$1"
    local label="$2"
    local path="$3"

    if [[ -s "$path" ]]; then
        record PASS "$section" "$label" "$path ($(file_size "$path"))"
    else
        record WARN "$section" "$label" "missing or empty: $path"
    fi
}

check_required_dir() {
    local section="$1"
    local label="$2"
    local path="$3"

    if [[ -d "$path" ]] && find "$path" -mindepth 1 -print -quit 2>/dev/null | grep -q .; then
        record PASS "$section" "$label" "$path"
    elif [[ -d "$path" ]]; then
        record FAIL "$section" "$label" "directory is empty: $path"
    else
        record FAIL "$section" "$label" "missing: $path"
    fi
}

check_command() {
    local command_name="$1"
    local resolved=""

    resolved="$(command -v "$command_name" 2>/dev/null || true)"
    if [[ -n "$resolved" ]]; then
        record PASS "System" "command: $command_name" "$resolved"
    else
        record FAIL "System" "command: $command_name" "not found in PATH"
    fi
}

check_shell_syntax() {
    local label="$1"
    local path="$2"
    local tmp="$TMP_WORK/syntax_$(basename "$path").txt"

    if [[ ! -f "$path" ]]; then
        record FAIL "Installer source" "syntax: $label" "missing: $path"
        return
    fi

    if capture_command "$tmp" bash -n "$path"; then
        record PASS "Installer source" "syntax: $label" "$(compact_output "$tmp")"
    else
        record FAIL "Installer source" "syntax: $label" "$(compact_output "$tmp" 20)"
    fi
}

check_perl_syntax() {
    local label="$1"
    local path="$2"
    local tmp="$TMP_WORK/perl_$(basename "$path").txt"

    if [[ ! -f "$path" ]]; then
        record FAIL "Installer source" "syntax: $label" "missing: $path"
        return
    fi

    if capture_command "$tmp" perl -c "$path"; then
        record PASS "Installer source" "syntax: $label" "$(compact_output "$tmp")"
    else
        record FAIL "Installer source" "syntax: $label" "$(compact_output "$tmp" 20)"
    fi
}

check_python_syntax() {
    local label="$1"
    local path="$2"
    local tmp="$TMP_WORK/python_$(basename "$path").txt"

    if [[ ! -f "$path" ]]; then
        record FAIL "Installer source" "compile: $label" "missing: $path"
        return
    fi

    if capture_command "$tmp" python3 -m py_compile "$path"; then
        record PASS "Installer source" "compile: $label" "$(compact_output "$tmp")"
    else
        record FAIL "Installer source" "compile: $label" "$(compact_output "$tmp" 20)"
    fi
}

conda_available() {
    [[ -x "$CONDA_PATH/bin/conda" ]]
}

env_prefix() {
    local env_name="$1"
    if [[ "$env_name" == "base" ]]; then
        printf '%s\n' "$CONDA_PATH"
    else
        printf '%s\n' "$CONDA_PATH/envs/$env_name"
    fi
}

env_exists() {
    local prefix
    prefix="$(env_prefix "$1")"
    [[ -d "$prefix/conda-meta" ]]
}

check_conda_env() {
    local env_name="$1"
    local prefix
    prefix="$(env_prefix "$env_name")"

    if [[ -d "$prefix/conda-meta" ]]; then
        record PASS "Conda" "environment: $env_name" "$prefix"
    else
        record FAIL "Conda" "environment: $env_name" "missing or invalid: $prefix"
    fi
}

check_env_command() {
    local env_name="$1"
    local command_name="$2"
    local tmp="$TMP_WORK/env_${env_name}_${command_name//[^A-Za-z0-9_.-]/_}.txt"

    if ! env_exists "$env_name"; then
        record FAIL "Tools/$env_name" "$command_name" "environment is unavailable"
        return
    fi

    if capture_command "$tmp" "$CONDA_PATH/bin/conda" run -n "$env_name" \
        bash -c "command -v '$command_name'"; then
        record PASS "Tools/$env_name" "$command_name" "$(compact_output "$tmp" 2)"
    else
        record FAIL "Tools/$env_name" "$command_name" "not found"
    fi
}

check_env_version() {
    local env_name="$1"
    local package_name="$2"
    local expected_prefix="${3:-}"
    local tmp="$TMP_WORK/version_${env_name}_${package_name}.txt"
    local version=""

    if ! env_exists "$env_name"; then
        record FAIL "Versions/$env_name" "$package_name" "environment is unavailable"
        return
    fi

    if capture_command "$tmp" "$CONDA_PATH/bin/conda" list -n "$env_name" \
        "$package_name" --json; then
        version="$(
            python3 - "$tmp" "$package_name" <<'PY'
import json
import sys
path, wanted = sys.argv[1:]
try:
    data = json.load(open(path))
except Exception:
    data = []
wanted = wanted.lower().replace("_", "-")
for item in data:
    name = str(item.get("name", "")).lower().replace("_", "-")
    if name == wanted:
        print(item.get("version", ""))
        break
PY
        )"
    fi

    if [[ -z "$version" ]]; then
        record FAIL "Versions/$env_name" "$package_name" "not installed"
    elif [[ -n "$expected_prefix" && "$version" != "$expected_prefix"* ]]; then
        record WARN "Versions/$env_name" "$package_name" \
            "$version (expected ${expected_prefix}*)"
    else
        if [[ -n "$expected_prefix" ]]; then
            record PASS "Versions/$env_name" "$package_name" \
                "$version (expected ${expected_prefix}*)"
        else
            record PASS "Versions/$env_name" "$package_name" "$version"
        fi
    fi
}

check_python_imports() {
    local env_name="$1"
    shift
    local modules=("$@")
    local joined
    local tmp="$TMP_WORK/imports_${env_name}.txt"

    joined="$(IFS=,; printf '%s' "${modules[*]}")"

    if capture_command "$tmp" "$CONDA_PATH/bin/conda" run -n "$env_name" \
        python - "$joined" <<'PY'
import importlib
import sys
mods = [x for x in sys.argv[1].split(",") if x]
failed = []
for mod in mods:
    try:
        importlib.import_module(mod)
    except Exception as exc:
        failed.append("%s: %s" % (mod, exc))
if failed:
    print(" | ".join(failed))
    raise SystemExit(1)
print("Imported: " + ", ".join(mods))
PY
    then
        record PASS "Python/$env_name" "module imports" "$(compact_output "$tmp" 6)"
    else
        record FAIL "Python/$env_name" "module imports" "$(compact_output "$tmp" 20)"
    fi
}

check_r_packages() {
    local env_name="$1"
    local label="$2"
    shift 2
    local packages=("$@")
    local package_string
    local tmp="$TMP_WORK/r_packages_${env_name}_${label//[^A-Za-z0-9]/_}.txt"

    package_string="$(IFS=,; printf '%s' "${packages[*]}")"

    if capture_command "$tmp" "$CONDA_PATH/bin/conda" run -n "$env_name" \
        Rscript - "$package_string" <<'RS'
args <- commandArgs(trailingOnly=TRUE)
pkgs <- strsplit(args[[1]], ",", fixed=TRUE)[[1]]
bad <- character()
for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly=TRUE)) {
        bad <- c(bad, paste0(pkg, "=NOT_INSTALLED"))
        next
    }
    result <- tryCatch({
        suppressPackageStartupMessages(
            library(pkg, character.only=TRUE, quietly=TRUE, warn.conflicts=FALSE)
        )
        paste0(pkg, "=", as.character(packageVersion(pkg)))
    }, error=function(e) paste0(pkg, "=LOAD_FAIL: ", conditionMessage(e)))
    cat(result, "\n")
    if (grepl("=LOAD_FAIL:", result, fixed=TRUE)) {
        bad <- c(bad, result)
    }
}
if (length(bad)) {
    cat("PROBLEM:", paste(bad, collapse=" | "), "\n")
    quit(status=1)
}
RS
    then
        record PASS "R/$env_name" "$label" \
            "${#packages[@]} package(s) installed and loadable"
    else
        record FAIL "R/$env_name" "$label" "$(compact_output "$tmp" 20)"
    fi
}

run_repo_r_checker() {
    local env_name="$1"
    local label="$2"
    local script_path="$3"
    local tmp="$TMP_WORK/repo_r_${env_name}_$(basename "$script_path").txt"

    if [[ ! -f "$script_path" ]]; then
        record WARN "R/$env_name" "$label" "checker script absent: $script_path"
        return
    fi

    if capture_command "$tmp" "$CONDA_PATH/bin/conda" run -n "$env_name" \
        bash "$script_path"; then
        record PASS "R/$env_name" "$label" "$(compact_output "$tmp" 8)"
    else
        record FAIL "R/$env_name" "$label" "$(compact_output "$tmp" 30)"
    fi
}

check_taxonomy_dir() {
    local section="$1"
    local label="$2"
    local taxonomy_dir="$3"
    local file

    for file in names.dmp nodes.dmp nucl_gb.accession2taxid nucl_wgs.accession2taxid; do
        check_required_file "$section" "$label/$file" "$taxonomy_dir/$file"
    done

    if [[ "$MODE" != "quick" ]]; then
        if [[ -s "$taxonomy_dir/nucl_gb.accession2taxid" ]] &&
           head -n 1 "$taxonomy_dir/nucl_gb.accession2taxid" |
           grep -q $'^accession\taccession.version\ttaxid'; then
            record PASS "$section" "$label nucl_gb header" "valid"
        else
            record FAIL "$section" "$label nucl_gb header" "invalid or unavailable"
        fi

        if [[ -s "$taxonomy_dir/nucl_wgs.accession2taxid" ]] &&
           head -n 1 "$taxonomy_dir/nucl_wgs.accession2taxid" |
           grep -q $'^accession\taccession.version\ttaxid'; then
            record PASS "$section" "$label nucl_wgs header" "valid"
        else
            record FAIL "$section" "$label nucl_wgs header" "invalid or unavailable"
        fi
    fi
}

check_kraken_core() {
    local database="$1"
    local label="$2"
    local file

    if [[ ! -d "$database" ]]; then
        record FAIL "Kraken DB" "$label directory" "missing: $database"
        return
    fi

    for file in hash.k2d opts.k2d taxo.k2d; do
        check_required_file "Kraken DB" "$label/$file" "$database/$file"
    done
}

check_library() {
    local database="$1"
    local library="$2"
    local severity="${3:-FAIL}"
    local library_dir="$database/library/$library"

    if [[ -d "$library_dir" ]]; then
        record PASS "Kraken library" "microbiome/$library directory" "$library_dir"
    else
        record "$severity" "Kraken library" "microbiome/$library directory" \
            "missing: $library_dir"
        return
    fi

local -a required_files=(
    library.fna
    prelim_map.txt
)

case "$library" in
    bacteria | archaea | protozoa | fungi | plasmid)
        required_files+=(manifest.txt)
        ;;

    UniVec | UniVec_Core)
        # UniVec libraries are distributed as a single source file.
        # Kraken2 creates library.fna and prelim_map.txt, but no manifest.txt.
        ;;
esac

for file in "${required_files[@]}"; do
    if [[ -s "$library_dir/$file" ]]; then
        record PASS "Kraken library" "microbiome/$library $file" \
            "$library_dir/$file ($(file_size "$library_dir/$file"))"
    else
        record "$severity" "Kraken library" "microbiome/$library $file" \
            "missing or empty: $library_dir/$file"
    fi
done

if [[ "$library" == "UniVec_Core" ]]; then
    if [[ -s "$library_dir/UniVec_Core" ]]; then
        record PASS \
            "Kraken library" \
            "microbiome/UniVec_Core source file" \
            "$library_dir/UniVec_Core ($(file_size "$library_dir/UniVec_Core"))"
    else
        record FAIL \
            "Kraken library" \
            "microbiome/UniVec_Core source file" \
            "missing or empty: $library_dir/UniVec_Core"
    fi

    local first_header=""
    local raw_count=0
    local library_count=0
    local map_count=0

    first_header="$(
        awk '/^>/{print; exit}' \
            "$library_dir/library.fna" \
            2>/dev/null || true
    )"

    if [[ "$first_header" == '>kraken:taxid|28384|'* ]]; then
        record PASS \
            "Kraken library" \
            "microbiome/UniVec_Core taxid" \
            "expected taxid 28384 detected"
    else
        record FAIL \
            "Kraken library" \
            "microbiome/UniVec_Core taxid" \
            "expected kraken:taxid|28384| was not detected"
    fi

    raw_count="$(
        grep -c '^>' "$library_dir/UniVec_Core" 2>/dev/null || true
    )"

    library_count="$(
        grep -c '^>' "$library_dir/library.fna" 2>/dev/null || true
    )"

    map_count="$(
        awk 'NF {n++} END {print n+0}' \
            "$library_dir/prelim_map.txt" \
            2>/dev/null
    )"

    if (( raw_count > 0 &&
          raw_count == library_count &&
          library_count == map_count )); then

        record PASS \
            "Kraken library" \
            "microbiome/UniVec_Core sequence alignment" \
            "$raw_count source sequences; $library_count library sequences; $map_count map entries"
    else
        record FAIL \
            "Kraken library" \
            "microbiome/UniVec_Core sequence alignment" \
            "source=$raw_count; library=$library_count; map=$map_count"
    fi

    if [[ -e "$library_dir/library.fna.masked" ]]; then
        record PASS \
            "Kraken library" \
            "microbiome/UniVec_Core masking marker" \
            "library.fna.masked exists; zero-byte marker is expected"
    else
        record WARN \
            "Kraken library" \
            "microbiome/UniVec_Core masking marker" \
            "library.fna.masked is absent"
    fi
fi
    if [[ "$MODE" == "deep" && -s "$library_dir/library.fna" ]]; then
        local header
        header="$(awk '/^>/{print; exit}' "$library_dir/library.fna" 2>/dev/null || true)"
        if [[ "$header" == *"kraken:taxid|"* ]]; then
            record PASS "Kraken library" "microbiome/$library FASTA header" \
                "Kraken taxid prefix detected"
        elif [[ "$library" == "plasmid" ]]; then
            record WARN "Kraken library" "microbiome/plasmid FASTA header" \
                "no kraken:taxid prefix in first header; prelim_map.txt is present"
        else
            record FAIL "Kraken library" "microbiome/$library FASTA header" \
                "expected kraken:taxid prefix not detected"
        fi
    fi
}

count_nonempty_lines() {
    local path="$1"
    if [[ -s "$path" ]]; then
        awk 'NF {n++} END {print n+0}' "$path"
    else
        printf '0\n'
    fi
}

check_manifest_alignment() {
    local library="$1"
    local installed="$MTD_DIR/kraken2DB_micro/library/$library/manifest.txt"
    local expected=""
    local installed_count=0
    local expected_count=0

    case "$library" in
        bacteria)
            expected="$OFFLINE_DIR/Kraken2DB_micro/library/bacteria/manifest_bacteria.list.txt"
            ;;
        archaea)
            expected="$OFFLINE_DIR/Kraken2DB_micro/library/archaea/manifest_archaea.list.txt"
            ;;
        plasmid)
            expected="$OFFLINE_DIR/Kraken2DB_micro/library/plasmid"
            ;;
    esac

    installed_count="$(count_nonempty_lines "$installed")"

    if [[ "$library" == "plasmid" ]]; then
        expected_count="$(
            find "$expected" -maxdepth 1 -type f -name '*.genomic.fna.gz' -size +0c \
                2>/dev/null | wc -l
        )"
    else
        expected_count="$(count_nonempty_lines "$expected")"
    fi

    if (( installed_count > 0 && expected_count > 0 && installed_count == expected_count )); then
        record PASS "Kraken library" "microbiome/$library installed/cache alignment" \
            "$installed_count installed manifest entries; $expected_count expected"
    elif (( installed_count > 0 && expected_count > 0 )); then
        record WARN "Kraken library" "microbiome/$library installed/cache alignment" \
            "$installed_count installed; $expected_count expected"
    else
        record WARN "Kraken library" "microbiome/$library installed/cache alignment" \
            "could not obtain both counts"
    fi
}

check_failed_download_record() {
    local name="$1"
    local path="$OFFLINE_DIR/failed_downloads_${name}.txt"

    if [[ ! -e "$path" || ! -s "$path" ]]; then
        record PASS "Installation cache" "failed-download record: $name" \
            "${path} is absent or empty"
    else
        local count
        count="$(count_nonempty_lines "$path")"
        record WARN "Installation cache" "failed-download record: $name" \
            "$count recorded item(s); inspect $path"
    fi
}

archive_integrity() {
    local label="$1"
    local path="$2"
    local tmp="$TMP_WORK/archive_$(basename "$path").txt"

    if [[ ! -s "$path" ]]; then
        record FAIL "Installation cache" "archive integrity: $label" "missing: $path"
        return
    fi

    case "$path" in
        *.tar.gz|*.tgz)
            if capture_command "$tmp" tar -tzf "$path"; then
                record PASS "Installation cache" "archive integrity: $label" \
                    "$(compact_output "$tmp" 4)"
            else
                record FAIL "Installation cache" "archive integrity: $label" \
                    "$(compact_output "$tmp" 20)"
            fi
            ;;
        *.gz)
            if capture_command "$tmp" gzip -t "$path"; then
                record PASS "Installation cache" "gzip integrity: $label" "valid"
            else
                record FAIL "Installation cache" "gzip integrity: $label" \
                    "$(compact_output "$tmp" 20)"
            fi
            ;;
    esac
}

init_colors

echo "============================================================"
echo "MTD Explorer installation check v$CHECKER_VERSION"
echo "============================================================"
echo "MTD directory : $MTD_DIR"
echo "Conda path    : $CONDA_PATH"
echo "Cache path    : ${OFFLINE_DIR:-not detected}"
echo "Mode          : $MODE"
echo "Bracken length: $READ_LEN"
echo "Report folder : $REPORT_DIR"
echo "============================================================"

# ==============================================================================
# Core paths and system
# ==============================================================================

if [[ -d "$MTD_DIR" ]]; then
    record PASS "Paths" "MTD directory" "$MTD_DIR"
else
    record FAIL "Paths" "MTD directory" "missing: $MTD_DIR"
fi

check_required_file "Paths" "current Install.sh" "$MTD_DIR/Install.sh"

if conda_available; then
    record PASS "Paths" "Conda executable" "$CONDA_PATH/bin/conda"
else
    record FAIL "Paths" "Conda executable" "missing: $CONDA_PATH/bin/conda"
fi

check_required_file "Paths" "Conda initialization" \
    "$CONDA_PATH/etc/profile.d/conda.sh"

if [[ -s "$MTD_DIR/condaPath" ]]; then
    stored_conda="$(head -n 1 "$MTD_DIR/condaPath" | tr -d '\r\n')"
    if [[ "$(expand_path "$stored_conda")" == "$CONDA_PATH" ]]; then
        record PASS "Paths" "MTD/condaPath" "$stored_conda"
    else
        record WARN "Paths" "MTD/condaPath" \
            "stored=$stored_conda; selected=$CONDA_PATH"
    fi
else
    record FAIL "Paths" "MTD/condaPath" "missing or empty"
fi

if [[ -n "$OFFLINE_DIR" && -d "$OFFLINE_DIR" ]]; then
    record PASS "Installation cache" "cache directory" "$OFFLINE_DIR"
else
    record FAIL "Installation cache" "cache directory" \
        "not detected or missing: ${OFFLINE_DIR:-unset}"
fi

for command_name in \
    bash awk sed grep find xargs gzip gunzip tar wget curl rsync pigz unpigz \
    perl python3 sha256sum timeout pkg-config aria2c parallel
do
    check_command "$command_name"
done

# ==============================================================================
# Installer contract: current architecture
# ==============================================================================

INSTALL_SH="$MTD_DIR/Install.sh"

if [[ -s "$INSTALL_SH" ]]; then
    check_shell_syntax "Install.sh" "$INSTALL_SH"

    for function_name in \
        parse_arguments validate_arguments ensure_sudo_credentials run_as_root \
        prepare_installation_cache prepare_shared_kraken2_taxonomy \
        install_shared_kraken2_taxonomy build_kraken2_database
    do
        if grep -Eq "^[[:space:]]*${function_name}[[:space:]]*\\(\\)" "$INSTALL_SH"; then
            record PASS "Installer contract" "function: $function_name" \
                "defined in Install.sh"
        else
            record FAIL "Installer contract" "function: $function_name" \
                "not defined"
        fi
    done

    if grep -q '^export CONDA_PLUGINS_AUTO_ACCEPT_TOS=yes' "$INSTALL_SH"; then
        record PASS "Installer contract" "Conda ToS automation" \
            "CONDA_PLUGINS_AUTO_ACCEPT_TOS=yes"
    else
        record WARN "Installer contract" "Conda ToS automation" \
            "expected export not detected"
    fi

    if grep -q 'sudo_with_pass' "$INSTALL_SH"; then
        record FAIL "Installer contract" "legacy sudo_with_pass" \
            "obsolete invocation remains"
    else
        record PASS "Installer contract" "legacy sudo_with_pass" \
            "not present"
    fi

    if grep -Eq '(^|[[:space:]])-w([[:space:]]|$)|sudo_password' "$INSTALL_SH"; then
        record WARN "Installer contract" "legacy sudo password option" \
            "-w or sudo_password reference detected"
    else
        record PASS "Installer contract" "legacy sudo password option" \
            "removed"
    fi

    if grep -Eq '(^|[[:space:]])-a([[:space:]]|$)|accept_required_conda_tos|accept_conda_tos' \
        "$INSTALL_SH"; then
        record WARN "Installer contract" "legacy ToS prompt option" \
            "old -a/prompt code detected"
    else
        record PASS "Installer contract" "legacy ToS prompt option" \
            "removed"
    fi

    if grep -q 'offlineCachePath' "$INSTALL_SH"; then
        record PASS "Installer contract" "writes offlineCachePath" \
            "referenced by Install.sh"
    else
        record FAIL "Installer contract" "writes offlineCachePath" \
            "reference not detected"
    fi

    if grep -q 'condaPath' "$INSTALL_SH"; then
        record PASS "Installer contract" "writes condaPath" \
            "referenced by Install.sh"
    else
        record FAIL "Installer contract" "writes condaPath" \
            "reference not detected"
    fi

    # The current default installation must not require predefined host databases.
    legacy_host_pattern='kraken2DB_(human|mice|rhesus)|hisat2_index_(human|mouse|rhesus)|human_blastdb|mouse_blastdb|rhesus_blastdb|build_(human|mouse|rhesus)_kraken_database|build_hisat2_host_index|build_magic_blast_databases'
    if grep -Eq "$legacy_host_pattern" "$INSTALL_SH"; then
        record WARN "Installer architecture" "predefined host database references" \
            "legacy host-database code is still referenced"
    else
        record PASS "Installer architecture" "predefined host databases removed" \
            "human/mouse/rhesus Kraken2, HISAT2, and BLAST are not required"
    fi

    if grep -q 'Customized_hosts' "$INSTALL_SH"; then
        record PASS "Installer architecture" "custom-host cache" \
            "Customized_hosts is prepared"
    else
        record WARN "Installer architecture" "custom-host cache" \
            "Customized_hosts reference not detected"
    fi

    if grep -q 'MTD_orgdb' "$INSTALL_SH"; then
        record PASS "Installer architecture" "MTD_orgdb integration" \
            "dedicated OrgDb environment referenced"
        EXPECT_MTD_ORGDB=1
    else
        EXPECT_MTD_ORGDB=0
        if [[ -d "$CONDA_PATH/envs/MTD_orgdb" ]]; then
            record WARN "Installer architecture" "MTD_orgdb integration" \
                "environment exists but reference was not detected in Install.sh"
            EXPECT_MTD_ORGDB=1
        else
            record SKIP "Installer architecture" "MTD_orgdb integration" \
                "not detected in this installer revision"
        fi
    fi
else
    EXPECT_MTD_ORGDB=0
fi

# ==============================================================================
# Source/configuration files
# ==============================================================================

if [[ "$MODE" != "quick" ]]; then
    for source_file in \
        Installation/MTD_fastp.yml \
        Installation/MTD.yml \
        Installation/MTD_R_additions.yml \
        Installation/py2.yml \
        Installation/halla0820.yml \
        Installation/R412.yml \
        Installation/pip.requirements \
        Installation/M33262_SIVMM239.fa \
        Installation/download_genomic_library.sh \
        Installation/download_genomic_library_plasmid.sh \
        Installation/rsync_from_ncbi.pl \
        Installation/rsync_from_ncbi_archaea.pl \
        Installation/rsync_from_ncbi_bacteria.pl \
        manifest.virus.sh \
        manifest.bacteria.sh \
        manifest.archea.sh \
        manifest.plasmid.sh \
        kraken2-build-download-taxonomy \
        Create_custom_host.sh \
        HostSpecies.csv \
        aux_scripts/orgdb/build_gold_orgdb_from_gtf_eggnog.R \
        aux_scripts/orgdb/install_gold_orgdb_from_hostspecies.sh \
        aux_scripts/orgdb/make_gene_representative_fasta.py
    do
        check_required_file "Installer source" "$source_file" "$MTD_DIR/$source_file"
    done

    for shell_script in \
        Install.sh \
        manifest.virus.sh \
        manifest.bacteria.sh \
        manifest.archea.sh \
        manifest.plasmid.sh \
        Create_custom_host.sh \
        aux_scripts/orgdb/install_gold_orgdb_from_hostspecies.sh \
        update_fix/Install.R.packages.MTD.sh \
        update_fix/check_R_pkg.MTD.sh \
        update_fix/Install.R.packages.R412_optimized.sh \
        update_fix/check_R_pkg.R412.sh \
        update_fix/check_R_pkg.halla0820.sh
    do
        if [[ -f "$MTD_DIR/$shell_script" ]]; then
            check_shell_syntax "$shell_script" "$MTD_DIR/$shell_script"
        fi
    done

    for perl_script in \
        Installation/rsync_from_ncbi.pl \
        Installation/rsync_from_ncbi_archaea.pl \
        Installation/rsync_from_ncbi_bacteria.pl
    do
        check_perl_syntax "$perl_script" "$MTD_DIR/$perl_script"
    done

    for python_script in \
        Installation/hisat2_extract_exons.py \
        Installation/hisat2_extract_splice_sites.py \
        aux_scripts/ssGSEA/resolve_ssgsea_go_terms.py \
        aux_scripts/orgdb/make_gene_representative_fasta.py
    do
        if [[ -f "$MTD_DIR/$python_script" ]]; then
            check_python_syntax "$python_script" "$MTD_DIR/$python_script"
        fi
    done
fi

# ==============================================================================
# Conda environments, tools, and package runtimes
# ==============================================================================

if conda_available; then
    for env_name in base MTD_fastp MTD py2 halla0820 R412; do
        check_conda_env "$env_name"
    done

    if (( EXPECT_MTD_ORGDB == 1 )); then
        check_conda_env MTD_orgdb
    fi

    check_env_version MTD_fastp fastp 1.3.3
    check_env_version MTD r-base 4.0.3
    check_env_version MTD kraken2 2.1.2
    check_env_version MTD bracken 2.6.0
    check_env_version MTD hisat2 2.2.1
    check_env_version MTD humann 3.1.1
    check_env_version py2 python 2.7
    check_env_version halla0820 python 3.10
    check_env_version halla0820 r-base 4.1.2
    check_env_version R412 r-base 4.1.2

    for command_name in fastp R Rscript; do
        check_env_command MTD_fastp "$command_name"
    done

    for command_name in \
        python R Rscript fastp kraken2 kraken2-build kraken2-inspect \
        bracken bracken-build hisat2 hisat2-build bowtie2 samtools \
        featureCounts makeblastdb blastdbcmd blastn magicblast humann \
        humann_config metaphlan diamond emapper.py datasets STAR \
        rsem-calculate-expression nextflow parallel
    do
        check_env_command MTD "$command_name"
    done

    for command_name in python hclust2.py; do
        check_env_command py2 "$command_name"
    done

    for command_name in python Rscript halla; do
        check_env_command halla0820 "$command_name"
    done

    for command_name in R Rscript kraken2; do
        check_env_command R412 "$command_name"
    done

    if (( EXPECT_MTD_ORGDB == 1 )); then
        for command_name in R Rscript jq yq; do
            check_env_command MTD_orgdb "$command_name"
        done
    fi

    if [[ "$MODE" != "quick" ]]; then
        check_python_imports MTD Bio biom numpy pandas scipy sklearn yaml rpy2
        check_python_imports py2 numpy pandas scipy matplotlib hclust2 biom
        check_python_imports halla0820 halla rpy2 numpy pandas scipy sklearn jinja2

        check_r_packages MTD_fastp "Venn/Euler packages" ggVennDiagram eulerr

        run_repo_r_checker MTD "required packages" \
            "$MTD_DIR/update_fix/check_R_pkg.MTD.sh"

        run_repo_r_checker R412 "required packages" \
            "$MTD_DIR/update_fix/check_R_pkg.R412.sh"

        run_repo_r_checker halla0820 "required packages" \
            "$MTD_DIR/update_fix/check_R_pkg.halla0820.sh"

        if (( EXPECT_MTD_ORGDB == 1 )); then
            check_r_packages MTD_orgdb "OrgDb construction packages" \
                AnnotationDbi AnnotationForge biomaRt GenomeInfoDb GO.db DBI RSQLite
        fi
    fi
else
    record FAIL "Conda" "environment checks" "Conda executable unavailable"
fi

# ==============================================================================
# Kraken2 helper restoration
# ==============================================================================

ACTIVE_LIBEXEC="$CONDA_PATH/envs/MTD/libexec"

if [[ "$MODE" != "quick" ]]; then
    if [[ -s "$MTD_DIR/Installation/rsync_from_ncbi.pl" &&
          -s "$ACTIVE_LIBEXEC/rsync_from_ncbi.pl" ]]; then
        if cmp -s "$MTD_DIR/Installation/rsync_from_ncbi.pl" \
            "$ACTIVE_LIBEXEC/rsync_from_ncbi.pl"; then
            record PASS "Kraken helpers" "rsync_from_ncbi.pl restored" \
                "active helper matches repository default"
        else
            record FAIL "Kraken helpers" "rsync_from_ncbi.pl restored" \
                "active helper differs from repository default"
        fi
    else
        record FAIL "Kraken helpers" "rsync_from_ncbi.pl restored" \
            "source or active helper is missing"
    fi

    if [[ -s "$MTD_DIR/Installation/download_genomic_library.sh" &&
          -s "$ACTIVE_LIBEXEC/download_genomic_library.sh" ]]; then
        if cmp -s "$MTD_DIR/Installation/download_genomic_library.sh" \
            "$ACTIVE_LIBEXEC/download_genomic_library.sh"; then
            record PASS "Kraken helpers" "download_genomic_library.sh restored" \
                "active helper matches repository default"
        else
            record FAIL "Kraken helpers" "download_genomic_library.sh restored" \
                "temporary plasmid helper may still be active"
        fi
    else
        record FAIL "Kraken helpers" "download_genomic_library.sh restored" \
            "source or active helper is missing"
    fi
fi

# ==============================================================================
# Default microbiome database only
# ==============================================================================

MICRO_DB="$MTD_DIR/kraken2DB_micro"

check_kraken_core "$MICRO_DB" "microbiome"
check_taxonomy_dir "Kraken taxonomy" "microbiome" "$MICRO_DB/taxonomy"

for library in bacteria archaea protozoa fungi plasmid UniVec_Core; do
    check_library "$MICRO_DB" "$library"
done

check_required_file "Kraken DB" "custom viral FASTA" \
    "$MTD_DIR/viruses4kraken.fa"

if [[ "$MODE" != "quick" && conda_available && -d "$MICRO_DB" ]]; then
    inspect_tmp="$TMP_WORK/kraken_inspect_microbiome.txt"
    if capture_command "$inspect_tmp" "$CONDA_PATH/bin/conda" run -n MTD \
        kraken2-inspect --db "$MICRO_DB"; then
        record PASS "Kraken DB" "inspect microbiome" \
            "$(head -n 6 "$inspect_tmp" | tr '\n' ' ' | sed -E 's/[[:space:]]+/ /g')"
    else
        record FAIL "Kraken DB" "inspect microbiome" \
            "$(compact_output "$inspect_tmp" 20)"
    fi
fi

if [[ -n "$OFFLINE_DIR" && "$MODE" != "quick" ]]; then
    check_manifest_alignment bacteria
    check_manifest_alignment archaea
    check_manifest_alignment plasmid
fi

# Predefined hosts are intentionally outside the default install.
record SKIP "Host databases" "predefined human database" \
    "not part of the current default installation"
record SKIP "Host databases" "predefined mouse database" \
    "not part of the current default installation"
record SKIP "Host databases" "predefined rhesus database" \
    "not part of the current default installation"
record SKIP "Host indexes" "predefined HISAT2/BLAST indexes" \
    "created later by the custom-host workflow when requested"

# ==============================================================================
# Bracken
# ==============================================================================

check_required_file "Bracken DB" "read length $READ_LEN distribution" \
    "$MICRO_DB/database${READ_LEN}mers.kmer_distrib"

check_required_file "Bracken DB" "database.kraken" \
    "$MICRO_DB/database.kraken"

# ==============================================================================
# HUMAnN
# ==============================================================================

HUMANN_ROOT="$MTD_DIR/HUMAnN/ref_database"
CHOCO_DIR="$HUMANN_ROOT/chocophlan"
UNIREF_DIR="$HUMANN_ROOT/uniref"
MAPPING_DIR="$HUMANN_ROOT/utility_mapping"

check_required_dir "HUMAnN" "ChocoPhlAn directory" "$CHOCO_DIR"
check_required_dir "HUMAnN" "UniRef directory" "$UNIREF_DIR"
check_required_dir "HUMAnN" "utility mapping directory" "$MAPPING_DIR"

choco_count="$(
    find "$CHOCO_DIR" -type f -name '*.ffn.gz' -size +0c 2>/dev/null | wc -l
)"
if (( choco_count > 0 )); then
    record PASS "HUMAnN" "ChocoPhlAn files" "$choco_count file(s)"
else
    record FAIL "HUMAnN" "ChocoPhlAn files" "none found"
fi

uniref_count="$(
    find "$UNIREF_DIR" -type f \
        \( -name '*.dmnd' -o -name '*.faa' -o -name '*.faa.gz' \) \
        -size +0c 2>/dev/null | wc -l
)"
if (( uniref_count > 0 )); then
    record PASS "HUMAnN" "UniRef database files" "$uniref_count candidate file(s)"
else
    record FAIL "HUMAnN" "UniRef database files" "no .dmnd/.faa files found"
fi

mapping_count="$(
    find "$MAPPING_DIR" -type f -size +0c 2>/dev/null | wc -l
)"
if (( mapping_count > 0 )); then
    record PASS "HUMAnN" "utility mapping files" "$mapping_count file(s)"
else
    record FAIL "HUMAnN" "utility mapping files" "none found"
fi

if conda_available && env_exists MTD; then
    humann_tmp="$TMP_WORK/humann_config.txt"
    if capture_command "$humann_tmp" "$CONDA_PATH/bin/conda" run -n MTD \
        humann_config; then
        humann_ok=1
        for expected_path in "$CHOCO_DIR" "$UNIREF_DIR" "$MAPPING_DIR"; do
            if ! grep -Fq "$expected_path" "$humann_tmp"; then
                humann_ok=0
            fi
        done

        if (( humann_ok == 1 )); then
            record PASS "HUMAnN" "configured database paths" \
                "all paths point inside $HUMANN_ROOT"
        else
            record FAIL "HUMAnN" "configured database paths" \
                "one or more expected paths are absent from humann_config"
        fi
    else
        record FAIL "HUMAnN" "configured database paths" \
            "$(compact_output "$humann_tmp" 20)"
    fi
fi

# ==============================================================================
# eggNOG database in persistent cache
# ==============================================================================

if [[ -n "$OFFLINE_DIR" ]]; then
    EGGNOG_DIR="$OFFLINE_DIR/eggNOG/emapperdb-5.0.2"

    check_required_dir "eggNOG DB" "operational database directory" "$EGGNOG_DIR"

    for eggnog_file in \
        eggnog.db \
        eggnog_proteins.dmnd \
        eggnog.taxa.db \
        eggnog.taxa.db.traverse.pkl
    do
        check_required_file "eggNOG DB" "$eggnog_file" "$EGGNOG_DIR/$eggnog_file"
    done

    if conda_available && env_exists MTD && [[ -s "$EGGNOG_DIR/eggnog_proteins.dmnd" ]]; then
        emapper_tmp="$TMP_WORK/emapper_version.txt"
        if capture_command "$emapper_tmp" "$CONDA_PATH/bin/conda" run -n MTD \
            emapper.py --version; then
            record PASS "eggNOG DB" "emapper runtime" \
                "$(compact_output "$emapper_tmp" 4)"
        else
            record FAIL "eggNOG DB" "emapper runtime" \
                "$(compact_output "$emapper_tmp" 20)"
        fi

        if [[ "$MODE" == "deep" ]]; then
            diamond_tmp="$TMP_WORK/diamond_dbinfo.txt"
            if capture_command "$diamond_tmp" "$CONDA_PATH/bin/conda" run -n MTD \
                diamond dbinfo -d "$EGGNOG_DIR/eggnog_proteins.dmnd"; then
                record PASS "eggNOG DB" "DIAMOND database information" \
                    "$(compact_output "$diamond_tmp" 8)"
            else
                record FAIL "eggNOG DB" "DIAMOND database information" \
                    "$(compact_output "$diamond_tmp" 20)"
            fi

            sqlite_tmp="$TMP_WORK/eggnog_sqlite.txt"
            if capture_command "$sqlite_tmp" python3 - "$EGGNOG_DIR/eggnog.db" <<'PY'
import sqlite3
import sys
db = sqlite3.connect("file:%s?mode=ro" % sys.argv[1], uri=True)
n = db.execute("select count(*) from sqlite_master where type='table'").fetchone()[0]
print("SQLite tables:", n)
db.close()
PY
            then
                record PASS "eggNOG DB" "eggnog.db read-only SQLite open" \
                    "$(compact_output "$sqlite_tmp" 4)"
            else
                record FAIL "eggNOG DB" "eggnog.db read-only SQLite open" \
                    "$(compact_output "$sqlite_tmp" 20)"
            fi
        fi
    fi
fi

# ==============================================================================
# Persistent cache
# ==============================================================================

if [[ -n "$OFFLINE_DIR" ]]; then
    check_required_file "Installation cache" \
        "Ref_genomes/MTD_virus/virushostdb.genomic.fna.gz" \
        "$OFFLINE_DIR/Ref_genomes/MTD_virus/virushostdb.genomic.fna.gz"

    for humann_archive in \
        full_mapping_v201901b.tar.gz \
        uniref90_annotated_v201901b_full.tar.gz \
        full_chocophlan.v201901_v31.tar.gz
    do
        check_required_file "Installation cache" \
            "HUMAnN/$humann_archive" \
            "$OFFLINE_DIR/HUMAnN/$humann_archive"
    done

    check_required_dir "Installation cache" "Kraken2 viral cache" \
        "$OFFLINE_DIR/Kraken2DB_micro/library/viral"
    check_required_dir "Installation cache" "Kraken2 bacteria cache" \
        "$OFFLINE_DIR/Kraken2DB_micro/library/bacteria"
    check_required_dir "Installation cache" "Kraken2 archaea cache" \
        "$OFFLINE_DIR/Kraken2DB_micro/library/archaea"
    check_required_dir "Installation cache" "Kraken2 plasmid cache" \
        "$OFFLINE_DIR/Kraken2DB_micro/library/plasmid"
    check_required_dir "Installation cache" "custom-host cache directory" \
        "$OFFLINE_DIR/Customized_hosts"

    check_failed_download_record bacteria
    check_failed_download_record archaea
    check_failed_download_record viral

    TAX_CACHE="$OFFLINE_DIR/Kraken2_taxonomy_cache"
    check_required_dir "Kraken taxonomy" "shared cache directory" "$TAX_CACHE"
    check_required_file "Kraken taxonomy" "shared cache completion marker" \
        "$TAX_CACHE/.mtd_taxonomy_complete"
    check_taxonomy_dir "Kraken taxonomy" "shared-cache" "$TAX_CACHE/taxonomy"

    if [[ "$MODE" == "deep" ]]; then
        archive_integrity \
            "Ref_genomes/MTD_virus/virushostdb.genomic.fna.gz" \
            "$OFFLINE_DIR/Ref_genomes/MTD_virus/virushostdb.genomic.fna.gz"

        archive_integrity \
            "HUMAnN/full_mapping_v201901b.tar.gz" \
            "$OFFLINE_DIR/HUMAnN/full_mapping_v201901b.tar.gz"

        archive_integrity \
            "HUMAnN/uniref90_annotated_v201901b_full.tar.gz" \
            "$OFFLINE_DIR/HUMAnN/uniref90_annotated_v201901b_full.tar.gz"

        archive_integrity \
            "HUMAnN/full_chocophlan.v201901_v31.tar.gz" \
            "$OFFLINE_DIR/HUMAnN/full_chocophlan.v201901_v31.tar.gz"
    fi
fi

# ==============================================================================
# Custom-host workflow
# ==============================================================================

if [[ "$MODE" != "quick" ]]; then
    if [[ -x "$MTD_DIR/Create_custom_host.sh" || -f "$MTD_DIR/Create_custom_host.sh" ]]; then
        record PASS "Custom host" "Create_custom_host.sh" \
            "$MTD_DIR/Create_custom_host.sh"
    else
        record FAIL "Custom host" "Create_custom_host.sh" "missing"
    fi

    if [[ -s "$MTD_DIR/HostSpecies.csv" ]]; then
        host_rows="$(
            awk -F',' 'NR>1 && NF {n++} END {print n+0}' "$MTD_DIR/HostSpecies.csv"
        )"
        record PASS "Custom host" "HostSpecies.csv" "$host_rows data row(s)"
    else
        record FAIL "Custom host" "HostSpecies.csv" "missing or empty"
    fi

    for path in \
        aux_scripts/orgdb/build_gold_orgdb_from_gtf_eggnog.R \
        aux_scripts/orgdb/install_gold_orgdb_from_hostspecies.sh \
        aux_scripts/orgdb/make_gene_representative_fasta.py
    do
        check_required_file "Custom host" "$path" "$MTD_DIR/$path"
    done

    if [[ -s "$MTD_DIR/Create_custom_host.sh" ]]; then
        if grep -q 'offlineCachePath' "$MTD_DIR/Create_custom_host.sh"; then
            record PASS "Custom host" "persistent cache integration" \
                "Create_custom_host.sh reads offlineCachePath"
        else
            record WARN "Custom host" "persistent cache integration" \
                "offlineCachePath reference not detected"
        fi

        if grep -q 'MTD_orgdb' "$MTD_DIR/Create_custom_host.sh"; then
            record PASS "Custom host" "MTD_orgdb usage" \
                "dedicated OrgDb environment referenced"
        elif (( EXPECT_MTD_ORGDB == 1 )); then
            record WARN "Custom host" "MTD_orgdb usage" \
                "environment expected but reference not detected in Create_custom_host.sh"
        else
            record SKIP "Custom host" "MTD_orgdb usage" \
                "not used by this revision"
        fi
    fi
fi

# ==============================================================================
# Summary
# ==============================================================================

FINAL_STATUS="PASS"
FINAL_EXIT=0

if (( FAIL_COUNT > 0 )); then
    FINAL_STATUS="FAIL"
    FINAL_EXIT=1
elif (( STRICT == 1 && WARN_COUNT > 0 )); then
    FINAL_STATUS="FAIL (strict warnings)"
    FINAL_EXIT=1
elif (( WARN_COUNT > 0 )); then
    FINAL_STATUS="PASS WITH WARNINGS"
fi

{
    echo "============================================================"
    echo "MTD Explorer installation check summary v$CHECKER_VERSION"
    echo "============================================================"
    echo "Date:        $(date '+%Y-%m-%d %H:%M:%S %z')"
    echo "Host:        $(hostname -f 2>/dev/null || hostname)"
    echo "MTD:         $MTD_DIR"
    echo "Conda:       $CONDA_PATH"
    echo "Mode:        $MODE"
    echo "Cache:       ${OFFLINE_DIR:-not detected}"
    echo "Checks:      $TOTAL_COUNT"
    echo "PASS:        $PASS_COUNT"
    echo "WARN:        $WARN_COUNT"
    echo "FAIL:        $FAIL_COUNT"
    echo "SKIP:        $SKIP_COUNT"
    echo "Results TSV: $RESULTS_TSV"
    echo "Full log:    $FULL_LOG"
    echo "============================================================"
    echo "FINAL STATUS: $FINAL_STATUS"
} | tee "$SUMMARY_TXT"

# Store the terminal-style results and summary in a compact full log.
{
    cat "$RESULTS_TSV"
    echo
    cat "$SUMMARY_TXT"
} > "$FULL_LOG"

exit "$FINAL_EXIT"
