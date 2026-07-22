#!/usr/bin/env bash

# ==============================================================================
# MTD Explorer installation checker
# Version: 2026.07.22-r9.3
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

CHECKER_VERSION="2026.07.22-r9.3"

MTD_DIR="$HOME/MTD"
CONDA_PATH=""
OFFLINE_DIR=""
READ_LEN=75
MODE="full"
REPORT_DIR=""
STRICT=0
KEEP_TEMP=0
HOST_TAXID=""

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

  --hostid TAXID       Check one installed custom-host reference
                       Default: automatically detect numeric host references

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
        --hostid)
            [[ $# -ge 2 ]] || { echo "ERROR: --hostid requires a value." >&2; exit 2; }
            HOST_TAXID="$2"
            shift 2
            ;;
        --hostid=*)
            HOST_TAXID="${1#*=}"
            shift
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

if [[ -n "$HOST_TAXID" ]] &&
   ! [[ "$HOST_TAXID" =~ ^[1-9][0-9]*$ ]]
then
    printf 'ERROR: --hostid must be a positive NCBI Taxonomy ID.\n' >&2
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
    local pycache_root="$TMP_WORK/python_bytecode"

    if [[ ! -f "$path" ]]; then
        record FAIL "Installer source" "compile: $label" "missing: $path"
        return
    fi

    mkdir -p "$pycache_root"

    if capture_command "$tmp" \
        env PYTHONPYCACHEPREFIX="$pycache_root" \
        python3 -m py_compile "$path"
    then
        record PASS "Installer source" "compile: $label" \
            "syntax valid; bytecode kept under report temporary directory"
    else
        record FAIL "Installer source" "compile: $label" \
            "$(compact_output "$tmp" 20)"
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
        *.tar)
            if capture_command "$tmp" tar -tf "$path"; then
                record PASS "Installation cache" "archive integrity: $label" \
                    "$(compact_output "$tmp" 4)"
            else
                record FAIL "Installation cache" "archive integrity: $label" \
                    "$(compact_output "$tmp" 20)"
            fi
            ;;
        *.bz2)
            if capture_command "$tmp" bzip2 -t "$path"; then
                record PASS "Installation cache" "bzip2 integrity: $label" "valid"
            else
                record FAIL "Installation cache" "bzip2 integrity: $label" \
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
        *)
            record SKIP "Installation cache" "archive integrity: $label" \
                "unsupported extension: $path"
            ;;
    esac
}
# ==============================================================================
# MTD CHECKER R9: dedicated HUMAnN and custom-host validation
# MTD CHECKER R9.1: repository-wide layout audit
# MTD CHECKER R9.2: interpreter-aware syntax validation
# MTD CHECKER R9.3: curated viral reference contract
# ==============================================================================

run_humann_isolated() {
    local prefix="$CONDA_PATH/envs/MTD_humann"

    env \
        -u PYTHONPATH \
        -u PYTHONHOME \
        PYTHONNOUSERSITE=1 \
        PATH="$prefix/bin:$PATH" \
        "$@"
}

check_humann_runtime() {
    local prefix="$CONDA_PATH/envs/MTD_humann"
    local tmp=""
    local command_name=""
    local -a required_commands=(
        python
        humann
        humann_config
        humann_test
        humann_join_tables
        humann_renorm_table
        humann_split_stratified_table
        humann_regroup_table
        metaphlan
        diamond
        bowtie2
        glpsol
        hclust2.py
    )

    if ! env_exists MTD_humann; then
        record FAIL "HUMAnN runtime" "dedicated environment" \
            "missing: $prefix"
        return
    fi

    for command_name in "${required_commands[@]}"; do
        if [[ -x "$prefix/bin/$command_name" ]]; then
            record PASS "Tools/MTD_humann" "$command_name" \
                "$prefix/bin/$command_name"
        else
            record FAIL "Tools/MTD_humann" "$command_name" \
                "missing or not executable: $prefix/bin/$command_name"
        fi
    done

    tmp="$TMP_WORK/humann39_version.txt"
    if capture_command "$tmp" run_humann_isolated \
        "$prefix/bin/humann" --version &&
       grep -Fq 'humann v3.9' "$tmp"
    then
        record PASS "Versions/MTD_humann" "HUMAnN" \
            "$(compact_output "$tmp" 4)"
    else
        record FAIL "Versions/MTD_humann" "HUMAnN" \
            "expected humann v3.9; observed: $(compact_output "$tmp" 8)"
    fi

    tmp="$TMP_WORK/metaphlan411_version.txt"
    if capture_command "$tmp" run_humann_isolated \
        "$prefix/bin/metaphlan" --version &&
       grep -Fq 'MetaPhlAn version 4.1.1' "$tmp"
    then
        record PASS "Versions/MTD_humann" "MetaPhlAn" \
            "$(compact_output "$tmp" 4)"
    else
        record FAIL "Versions/MTD_humann" "MetaPhlAn" \
            "expected 4.1.1; observed: $(compact_output "$tmp" 8)"
    fi

    tmp="$TMP_WORK/humann_diamond_version.txt"
    if capture_command "$tmp" run_humann_isolated \
        "$prefix/bin/diamond" version &&
       grep -Fq '2.0.15' "$tmp"
    then
        record PASS "Versions/MTD_humann" "DIAMOND" \
            "$(compact_output "$tmp" 4)"
    else
        record FAIL "Versions/MTD_humann" "DIAMOND" \
            "expected 2.0.15; observed: $(compact_output "$tmp" 8)"
    fi

    tmp="$TMP_WORK/humann_bowtie2_version.txt"
    if capture_command "$tmp" run_humann_isolated \
        "$prefix/bin/bowtie2" --version &&
       grep -Fq 'version 2.5.4' "$tmp"
    then
        record PASS "Versions/MTD_humann" "Bowtie2" \
            "$(compact_output "$tmp" 4)"
    else
        record FAIL "Versions/MTD_humann" "Bowtie2" \
            "expected 2.5.4; observed: $(compact_output "$tmp" 8)"
    fi

    tmp="$TMP_WORK/humann_python_isolation.txt"
    if capture_command "$tmp" run_humann_isolated \
        "$prefix/bin/python" - "$prefix" <<'PY_HUMANN_ISOLATION'
import os
import site
import sys
from pathlib import Path

import Cython
import humann
import numpy
import pysam
import simplejson

expected = Path(sys.argv[1]).resolve()
observed = Path(sys.prefix).resolve()

assert observed == expected, (observed, expected)
assert sys.version_info[:2] == (3, 10), sys.version
assert numpy.__version__ == "1.26.4", numpy.__version__
assert simplejson.__version__.split(".", 1)[0] == "3", simplejson.__version__
assert site.ENABLE_USER_SITE is False, site.ENABLE_USER_SITE
assert all("/.local/" not in entry for entry in sys.path), sys.path

for module in (Cython, humann, numpy, pysam, simplejson):
    module_path = Path(module.__file__).resolve()
    assert str(module_path).startswith(str(expected) + os.sep), module_path

print("Python", sys.version.split()[0])
print("NumPy", numpy.__version__)
print("simplejson", simplejson.__version__)
print("user-site disabled and module paths isolated")
PY_HUMANN_ISOLATION
    then
        record PASS "Python/MTD_humann" "isolation and pinned modules" \
            "$(compact_output "$tmp" 8)"
    else
        record FAIL "Python/MTD_humann" "isolation and pinned modules" \
            "$(compact_output "$tmp" 20)"
    fi
}

check_custom_host_reference() {
    local taxid="$1"
    local db="$MTD_DIR/kraken2DB_${taxid}"
    local clean_fasta="$db/genome_${taxid}.fa"
    local kraken_fasta="$db/genome_${taxid}.kraken.fa"
    local gtf="$MTD_DIR/ref_${taxid}/ref_${taxid}.gtf.gz"
    local hisat_prefix="$MTD_DIR/hisat2_index_${taxid}/genome_tran"
    local inspect_bin="$CONDA_PATH/envs/MTD/bin/hisat2-inspect"
    local extension=""
    local part=""
    local clean_headers=0
    local clean_tagged=0
    local kraken_headers=0
    local kraken_tagged=0
    local tmp=""

    check_kraken_core "$db" "host $taxid"
    check_required_file "Host $taxid" "clean alignment FASTA" "$clean_fasta"
    check_required_file "Host $taxid" "Kraken2-only FASTA" "$kraken_fasta"
    check_required_file "Host $taxid" "GTF annotation" "$gtf"

    if [[ -s "$clean_fasta" ]]; then
        clean_headers="$(grep -c '^>' "$clean_fasta" 2>/dev/null || true)"
        clean_tagged="$(grep -c '^>kraken:taxid|' "$clean_fasta" 2>/dev/null || true)"

        if (( clean_headers > 0 && clean_tagged == 0 )); then
            record PASS "Host $taxid" "clean FASTA identifiers" \
                "$clean_headers headers; no Kraken prefix"
        else
            record FAIL "Host $taxid" "clean FASTA identifiers" \
                "headers=$clean_headers; Kraken-prefixed=$clean_tagged"
        fi
    fi

    if [[ -s "$kraken_fasta" ]]; then
        kraken_headers="$(grep -c '^>' "$kraken_fasta" 2>/dev/null || true)"
        kraken_tagged="$(
            grep -c "^>kraken:taxid|${taxid}|" \
                "$kraken_fasta" 2>/dev/null || true
        )"

        if (( kraken_headers > 0 && kraken_headers == kraken_tagged )); then
            record PASS "Host $taxid" "Kraken FASTA identifiers" \
                "$kraken_tagged/$kraken_headers headers use TaxID $taxid"
        else
            record FAIL "Host $taxid" "Kraken FASTA identifiers" \
                "matching=$kraken_tagged; total=$kraken_headers"
        fi
    fi

    if (( clean_headers > 0 && kraken_headers > 0 )); then
        if (( clean_headers == kraken_headers )); then
            record PASS "Host $taxid" "clean/Kraken FASTA sequence count" \
                "$clean_headers sequences in each FASTA"
        else
            record FAIL "Host $taxid" "clean/Kraken FASTA sequence count" \
                "clean=$clean_headers; Kraken=$kraken_headers"
        fi
    fi

    if [[ -s "${hisat_prefix}.1.ht2" ]]; then
        extension="ht2"
    elif [[ -s "${hisat_prefix}.1.ht2l" ]]; then
        extension="ht2l"
    fi

    if [[ -n "$extension" ]]; then
        for part in 1 2 3 4 5 6 7 8; do
            check_required_file "Host $taxid" \
                "HISAT2 genome_tran.$part.$extension" \
                "${hisat_prefix}.${part}.${extension}"
        done
    else
        record FAIL "Host $taxid" "HISAT2 genome_tran index" \
            "no .ht2 or .ht2l index detected: $hisat_prefix"
    fi

    if [[ "$MODE" != "quick" && -s "$clean_fasta" && -s "$gtf" ]]; then
        tmp="$TMP_WORK/host_${taxid}_fasta_gtf.txt"

        if capture_command "$tmp" python3 - \
            "$clean_fasta" "$gtf" <<'PY_FASTA_GTF'
import gzip
import sys
from pathlib import Path

fasta_path = Path(sys.argv[1])
gtf_path = Path(sys.argv[2])

fasta_ids = set()
with fasta_path.open("rt", encoding="utf-8", errors="replace") as handle:
    for line in handle:
        if line.startswith(">"):
            identifier = line[1:].split(None, 1)[0].strip()
            if identifier:
                fasta_ids.add(identifier)

gtf_ids = set()
with gzip.open(gtf_path, "rt", encoding="utf-8", errors="replace") as handle:
    for line in handle:
        if not line or line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if fields and fields[0]:
            gtf_ids.add(fields[0])

assert fasta_ids, "no FASTA identifiers"
assert gtf_ids, "no GTF contig identifiers"
assert not any(x.startswith("kraken:taxid|") for x in fasta_ids)

missing = sorted(gtf_ids - fasta_ids)
if missing:
    print("GTF contigs absent from FASTA:", len(missing))
    for item in missing[:20]:
        print(item)
    raise SystemExit(1)

print("clean FASTA identifiers:", len(fasta_ids))
print("GTF contigs:", len(gtf_ids))
print("shared contigs:", len(fasta_ids & gtf_ids))
PY_FASTA_GTF
        then
            record PASS "Host $taxid" "clean FASTA/GTF compatibility" \
                "$(compact_output "$tmp" 8)"
        else
            record FAIL "Host $taxid" "clean FASTA/GTF compatibility" \
                "$(compact_output "$tmp" 20)"
        fi
    fi

    if [[ "$MODE" != "quick" && -x "$inspect_bin" && -n "$extension" ]]; then
        tmp="$TMP_WORK/host_${taxid}_hisat2_names.txt"

        if capture_command "$tmp" "$inspect_bin" -n "$hisat_prefix"; then
            if grep -q '^kraken:taxid|' "$tmp"; then
                record FAIL "Host $taxid" "HISAT2 clean identifiers" \
                    "Kraken prefixes remain in the HISAT2 index"
            else
                record PASS "Host $taxid" "HISAT2 clean identifiers" \
                    "no Kraken prefix detected"
            fi

            if [[ -s "$gtf" ]]; then
                local compare_tmp="$TMP_WORK/host_${taxid}_hisat2_gtf.txt"

                if capture_command "$compare_tmp" python3 - \
                    "$tmp" "$gtf" <<'PY_HISAT_GTF'
import gzip
import sys
from pathlib import Path

index_path = Path(sys.argv[1])
gtf_path = Path(sys.argv[2])

index_ids = {
    line.strip()
    for line in index_path.read_text(
        encoding="utf-8", errors="replace"
    ).splitlines()
    if line.strip()
}

gtf_ids = set()
with gzip.open(gtf_path, "rt", encoding="utf-8", errors="replace") as handle:
    for line in handle:
        if not line or line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if fields and fields[0]:
            gtf_ids.add(fields[0])

missing = sorted(gtf_ids - index_ids)
if missing:
    print("GTF contigs absent from HISAT2 index:", len(missing))
    for item in missing[:20]:
        print(item)
    raise SystemExit(1)

print("HISAT2 identifiers:", len(index_ids))
print("GTF contigs:", len(gtf_ids))
print("shared contigs:", len(index_ids & gtf_ids))
PY_HISAT_GTF
                then
                    record PASS "Host $taxid" "HISAT2/GTF compatibility" \
                        "$(compact_output "$compare_tmp" 8)"
                else
                    record FAIL "Host $taxid" "HISAT2/GTF compatibility" \
                        "$(compact_output "$compare_tmp" 20)"
                fi
            fi
        else
            record FAIL "Host $taxid" "HISAT2 index inspection" \
                "$(compact_output "$tmp" 20)"
        fi
    elif [[ "$MODE" != "quick" && ! -x "$inspect_bin" ]]; then
        record FAIL "Host $taxid" "hisat2-inspect" \
            "missing or not executable: $inspect_bin"
    fi
}

check_installed_custom_hosts() {
    local -a taxids=()
    local db=""
    local taxid=""

    if [[ -n "$HOST_TAXID" ]]; then
        taxids=("$HOST_TAXID")
    else
        while IFS= read -r db; do
            taxid="${db##*/kraken2DB_}"
            if [[ "$taxid" =~ ^[1-9][0-9]*$ ]]; then
                taxids+=("$taxid")
            fi
        done < <(
            find "$MTD_DIR" \
                -maxdepth 1 \
                -type d \
                -name 'kraken2DB_[0-9]*' \
                -print 2>/dev/null | sort
        )
    fi

    if (( ${#taxids[@]} == 0 )); then
        record SKIP "Custom hosts" "installed references" \
            "no numeric custom-host directory detected"
        return
    fi

    for taxid in "${taxids[@]}"; do
        check_custom_host_reference "$taxid"
    done
}

# ==============================================================================
# MTD CHECKER R9.3: curated viral reference contract
# ==============================================================================

viral_summary_metric() {
    local summary_file="$1"
    local metric="$2"

    awk -F $'\t' -v wanted="$metric" '
        $1 == wanted {
            print $2
            exit
        }
    ' "$summary_file" 2>/dev/null
}

check_curated_viral_installer_contract() {
    local installer="$MTD_DIR/Install.sh"
    local required_marker=""

    if [[ ! -s "$installer" ]]; then
        record FAIL "Viral contract" "Install.sh" \
            "missing or empty: $installer"
        return
    fi

    for required_marker in \
        prepare_virushost_release_cache \
        Ref_genomes/MTD_virus/official_current \
        viral_genomes_combined_nonredundant.fna \
        viral_genomes_combined_nonredundant.summary.tsv \
        viral_genomes_combined_nonredundant.details.tsv \
        records_without_taxid \
        records_without_accession
    do
        if grep -Fq "$required_marker" "$installer"; then
            record PASS "Viral contract" \
                "Install.sh marker: $required_marker" \
                "current curated viral workflow detected"
        else
            record FAIL "Viral contract" \
                "Install.sh marker: $required_marker" \
                "required marker is absent"
        fi
    done

    if grep -Fq 'viruses4kraken.fa' "$installer"; then
        record FAIL "Viral contract" "legacy installer FASTA" \
            "Install.sh still references viruses4kraken.fa"
    else
        record PASS "Viral contract" "legacy installer FASTA removed" \
            "Install.sh does not reference viruses4kraken.fa"
    fi

    if grep -Fq \
        'VIRUSHOST_MIRROR_SHA256SUMS_SHA256=' \
        "$installer"
    then
        record PASS "Viral contract" \
            "pinned Virus-Host checksum manifest" \
            "checksum-of-manifest constant is present"
    else
        record FAIL "Viral contract" \
            "pinned Virus-Host checksum manifest" \
            "constant is absent from Install.sh"
    fi
}

check_curated_viral_reference_cache() {
    local release_dir=""
    local viral_dir=""
    local viral_download_dir=""
    local final_fasta=""
    local final_summary=""
    local final_details=""
    local old_root_fasta=""
    local old_release_fasta=""
    local required_file=""
    local summary_records_seen=""
    local summary_records_written=""
    local summary_without_taxid=""
    local summary_without_accession=""
    local observed_manifest_hash=""
    local details_header=""
    local expected_details_header=""
    local tmp=""
    local header_limit=10000
    local expected_count=0
    local local_count=0
    local stale_marker=""

    check_curated_viral_installer_contract

    old_root_fasta="$MTD_DIR/viruses4kraken.fa"

    if [[ -e "$old_root_fasta" ]]; then
        record WARN "Viral legacy" "viruses4kraken.fa" \
            "obsolete file is still present but is no longer consumed: $old_root_fasta"
    else
        record SKIP "Viral legacy" "viruses4kraken.fa" \
            "obsolete root-level FASTA is correctly absent"
    fi

    if [[ -z "$OFFLINE_DIR" ]]; then
        record SKIP "Viral cache" "curated reference stack" \
            "persistent cache path was not detected"
        return
    fi

    release_dir="$OFFLINE_DIR/Ref_genomes/MTD_virus/official_current"
    viral_dir="$OFFLINE_DIR/Kraken2DB_micro/library/viral"
    viral_download_dir="$viral_dir/all"

    final_fasta="$viral_dir/viral_genomes_combined_nonredundant.fna"
    final_summary="$viral_dir/viral_genomes_combined_nonredundant.summary.tsv"
    final_details="$viral_dir/viral_genomes_combined_nonredundant.details.tsv"

    old_release_fasta="$OFFLINE_DIR/Ref_genomes/MTD_virus/virushostdb.genomic.fna.gz"
    stale_marker="$viral_dir/all_viral_genomes.fna.STALE"

    if [[ -e "$old_release_fasta" ]]; then
        record WARN "Viral legacy" "old Virus-Host DB cache path" \
            "obsolete location remains; current release is under official_current/: $old_release_fasta"
    else
        record SKIP "Viral legacy" "old Virus-Host DB cache path" \
            "obsolete cache location is correctly absent"
    fi

    check_required_dir "Viral cache" \
        "pinned Virus-Host DB release" \
        "$release_dir"

    check_required_dir "Viral cache" \
        "RefSeq viral genome directory" \
        "$viral_download_dir"

    # Immutable mirrored release and provenance.
    for required_file in \
        virushostdb.genomic.fna.gz \
        non-segmented_virus_list.tsv \
        segmented_virus_list.tsv \
        dbrel.txt \
        SHA256SUMS \
        MIRROR_METADATA.tsv
    do
        check_required_file "Virus-Host DB" \
            "official_current/$required_file" \
            "$release_dir/$required_file"
    done

    if [[ -s "$release_dir/SHA256SUMS" ]]; then
        observed_manifest_hash="$(
            sha256sum "$release_dir/SHA256SUMS" 2>/dev/null |
                awk '{print $1}'
        )"

        if [[ "$observed_manifest_hash" == \
              "a250b2e61d9f9365773205d04d019e0976a778ebd589553f3a0a0e6f159f4bec" ]]
        then
            record PASS "Virus-Host DB" \
                "SHA256SUMS manifest identity" \
                "$observed_manifest_hash"
        else
            record FAIL "Virus-Host DB" \
                "SHA256SUMS manifest identity" \
                "expected=a250b2e61d9f9365773205d04d019e0976a778ebd589553f3a0a0e6f159f4bec; observed=${observed_manifest_hash:-unavailable}"
        fi
    fi

    if [[ -s "$release_dir/virushostdb.genomic.fna.gz" ]]; then
        tmp="$TMP_WORK/virushost_gzip_r9_3.txt"

        if capture_command "$tmp" \
            gzip -t "$release_dir/virushostdb.genomic.fna.gz"
        then
            record PASS "Virus-Host DB" "compressed FASTA integrity" \
                "gzip stream is valid"
        else
            record FAIL "Virus-Host DB" "compressed FASTA integrity" \
                "$(compact_output "$tmp" 20)"
        fi
    fi

    if [[ "$MODE" != "quick" && -s "$release_dir/SHA256SUMS" ]]; then
        tmp="$TMP_WORK/virushost_checksums_r9_3.txt"

        if capture_command "$tmp" \
            bash -c '
                cd "$1" &&
                sha256sum -c SHA256SUMS
            ' _ "$release_dir"
        then
            record PASS "Virus-Host DB" \
                "mirrored release checksums" \
                "$(compact_output "$tmp" 12)"
        else
            record FAIL "Virus-Host DB" \
                "mirrored release checksums" \
                "$(compact_output "$tmp" 30)"
        fi
    fi

    # Derived release files and accession mappings.
    if [[ "$MODE" != "quick" ]]; then
        for required_file in \
            virushostdb.genomic.fna \
            virushostdb_accession2taxid.tsv \
            virushostdb_accession_conflicts.tsv
        do
            check_required_file "Virus-Host DB" \
                "derived/$required_file" \
                "$release_dir/$required_file"
        done
    fi

    # RefSeq viral synchronization products.
    # all_viral_genomes.fna is an intermediate RefSeq aggregate, not the final
    # Kraken2 input.
    for required_file in \
        assembly_summary_viral.txt \
        manifest_viral.list.txt \
        manifest_viral.names.txt \
        all_viral_genomes.fna
    do
        check_required_file "RefSeq viral cache" \
            "$required_file" \
            "$viral_dir/$required_file"
    done

    if [[ "$MODE" != "quick" ]]; then
        for required_file in \
            integrity_viral.stat.tsv \
            all_viral_genomes.source_state.tsv \
            all_viral_genomes.output_state.tsv \
            refseq_viral_accession2taxid.tsv \
            refseq_viral_accession_conflicts.tsv
        do
            check_required_file "RefSeq viral cache" \
                "$required_file" \
                "$viral_dir/$required_file"
        done

        check_required_file "RefSeq viral cache" \
            "copied manifest.virus.sh" \
            "$OFFLINE_DIR/Kraken2DB_micro/library/manifest.virus.sh"
    fi

    if [[ -e "$stale_marker" ]]; then
        if [[ "$MODE" == "quick" ]]; then
            record WARN "RefSeq viral cache" \
                "collection completeness marker" \
                "partial/stale marker exists: $stale_marker"
        else
            record FAIL "RefSeq viral cache" \
                "collection completeness marker" \
                "installer requires a complete collection but stale marker exists: $stale_marker"
        fi
    else
        record PASS "RefSeq viral cache" \
            "collection completeness marker" \
            "no stale marker detected"
    fi

    if [[ -s "$viral_dir/manifest_viral.names.txt" &&
          -d "$viral_download_dir" ]]
    then
        expected_count="$(
            awk 'NF {n++} END {print n+0}' \
                "$viral_dir/manifest_viral.names.txt"
        )"

        local_count="$(
            find "$viral_download_dir" \
                -maxdepth 1 \
                -type f \
                -name '*.gz' \
                -size +0c \
                2>/dev/null |
            awk 'END {print NR + 0}'
        )"

        if (( expected_count > 0 && local_count == expected_count )); then
            record PASS "RefSeq viral cache" \
                "catalog/local completeness" \
                "$local_count local genome(s); $expected_count expected"
        else
            record FAIL "RefSeq viral cache" \
                "catalog/local completeness" \
                "$local_count local genome(s); $expected_count expected"
        fi
    fi

    # Final curated collection consumed by kraken2-build.
    check_required_file "Curated viral library" \
        "nonredundant taxid-aware FASTA" \
        "$final_fasta"

    check_required_file "Curated viral library" \
        "deduplication summary" \
        "$final_summary"

    check_required_file "Curated viral library" \
        "deduplication details" \
        "$final_details"

    if [[ -s "$final_summary" ]]; then
        summary_records_seen="$(
            viral_summary_metric "$final_summary" records_seen
        )"
        summary_records_written="$(
            viral_summary_metric "$final_summary" records_written
        )"
        summary_without_taxid="$(
            viral_summary_metric "$final_summary" records_without_taxid
        )"
        summary_without_accession="$(
            viral_summary_metric "$final_summary" records_without_accession
        )"

        if [[ "$summary_records_seen" =~ ^[0-9]+$ &&
              "$summary_records_written" =~ ^[0-9]+$ &&
              "$summary_without_taxid" =~ ^[0-9]+$ &&
              "$summary_without_accession" =~ ^[0-9]+$ ]]
        then
            record PASS "Curated viral library" \
                "summary schema and numeric metrics" \
                "records_seen=$summary_records_seen; records_written=$summary_records_written"
        else
            record FAIL "Curated viral library" \
                "summary schema and numeric metrics" \
                "one or more required metrics are absent or non-numeric"
        fi

        if [[ "$summary_records_written" =~ ^[0-9]+$ ]] &&
           (( summary_records_written > 0 ))
        then
            record PASS "Curated viral library" \
                "records written" \
                "$summary_records_written sequence record(s)"
        else
            record FAIL "Curated viral library" \
                "records written" \
                "expected a positive records_written metric"
        fi

        if [[ "$summary_records_seen" =~ ^[0-9]+$ &&
              "$summary_records_written" =~ ^[0-9]+$ ]] &&
           (( summary_records_seen >= summary_records_written ))
        then
            record PASS "Curated viral library" \
                "summary count relationship" \
                "records_seen >= records_written"
        else
            record FAIL "Curated viral library" \
                "summary count relationship" \
                "invalid or inconsistent records_seen/records_written values"
        fi

        if [[ "$summary_without_taxid" == "0" ]]; then
            record PASS "Curated viral library" \
                "records without TaxID" \
                "0"
        else
            record FAIL "Curated viral library" \
                "records without TaxID" \
                "${summary_without_taxid:-metric missing}"
        fi

        if [[ "$summary_without_accession" == "0" ]]; then
            record PASS "Curated viral library" \
                "records without accession" \
                "0"
        else
            record FAIL "Curated viral library" \
                "records without accession" \
                "${summary_without_accession:-metric missing}"
        fi
    fi

    if [[ -s "$final_details" ]]; then
        details_header="$(
            head -n 1 "$final_details" 2>/dev/null || true
        )"

        expected_details_header="$(
            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' \
                source status length accession taxid canonical_sha256 \
                header matched_source matched_accession matched_taxid \
                matched_header
        )"

        if [[ "$details_header" == "$expected_details_header" ]]; then
            record PASS "Curated viral library" \
                "details table schema" \
                "expected 11-column header detected"
        else
            record FAIL "Curated viral library" \
                "details table schema" \
                "unexpected header: ${details_header:-empty}"
        fi
    fi

    if [[ -s "$final_fasta" ]]; then
        if [[ "$MODE" == "deep" ]]; then
            header_limit=0
        else
            header_limit=10000
        fi

        tmp="$TMP_WORK/viral_final_headers_r9_3.txt"

        if capture_command "$tmp" \
            awk -v limit="$header_limit" '
                /^>/ {
                    headers++

                    if ($0 !~ /^>kraken:taxid\|[1-9][0-9]*\|/) {
                        bad++
                    }

                    if (limit > 0 && headers >= limit) {
                        exit
                    }
                }

                END {
                    print "headers_checked=" (headers + 0)
                    print "invalid_headers=" (bad + 0)

                    if (headers == 0 || bad > 0) {
                        exit 1
                    }
                }
            ' "$final_fasta"
        then
            if (( header_limit == 0 )); then
                record PASS "Curated viral library" \
                    "all FASTA headers carry TaxIDs" \
                    "$(compact_output "$tmp" 4)"
            else
                record PASS "Curated viral library" \
                    "sampled FASTA headers carry TaxIDs" \
                    "$(compact_output "$tmp" 4)"
            fi
        else
            record FAIL "Curated viral library" \
                "FASTA TaxID header validation" \
                "$(compact_output "$tmp" 20)"
        fi
    fi
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
echo "Host TaxID    : ${HOST_TAXID:-auto-detect}"
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
    perl python3 sha256sum md5sum bzip2 timeout pkg-config aria2c parallel
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
        install_shared_kraken2_taxonomy build_kraken2_database \
        prepare_humann_metaphlan_cache validate_humann_metaphlan_cache \
        validate_humann_environment validate_installed_humann_databases \
        configure_humann_database_paths install_humann_databases
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

    installer_getopts="$(
    sed -nE \
        's/.*while[[:space:]]+getopts[[:space:]]+"([^"]+)".*/\1/p' \
        "$INSTALL_SH" |
    head -n 1
)"

legacy_sudo_password=0

if grep -Eq \
    '(^|[^[:alnum:]_])sudo_password([^[:alnum:]_]|$)|^[[:space:]]*w\)[[:space:]]' \
    "$INSTALL_SH"; then
    legacy_sudo_password=1
fi

if [[ "$installer_getopts" == *"w:"* ]]; then
    legacy_sudo_password=1
fi

if grep -Eq \
    '^[[:space:]]*-w[[:space:]]+(TEXT|PASSWORD|PASS)' \
    "$INSTALL_SH"; then
    legacy_sudo_password=1
fi

if (( legacy_sudo_password == 1 )); then
    record WARN \
        "Installer contract" \
        "legacy sudo password option" \
        "old -w or sudo_password implementation detected"
else
    record PASS \
        "Installer contract" \
        "legacy sudo password option" \
        "removed; Bash file-test operator -w is allowed"
fi

legacy_tos_prompt=0

if grep -Eq \
    'accept_required_conda_tos|accept_conda_tos|^[[:space:]]*a\)[[:space:]]' \
    "$INSTALL_SH"; then
    legacy_tos_prompt=1
fi

if [[ "$installer_getopts" == *a* ]]; then
    legacy_tos_prompt=1
fi

if grep -Eq \
    '^[[:space:]]*-a[[:space:]]+(Accept|AUTO|Auto|Anaconda|TOS|ToS)' \
    "$INSTALL_SH"; then
    legacy_tos_prompt=1
fi

if (( legacy_tos_prompt == 1 )); then
    record WARN \
        "Installer contract" \
        "legacy ToS prompt option" \
        "old -a option or interactive ToS implementation detected"
else
    record PASS \
        "Installer contract" \
        "legacy ToS prompt option" \
        "removed; unrelated command options such as cp -a are allowed"
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
# Repository layout and source/configuration files
# ==============================================================================
# This contract follows the current MTD Explorer tree. It distinguishes:
#   1. required runtime/source directories;
#   2. critical files consumed by the installer and pipeline;
#   3. optional documentation/development directories;
#   4. every file tracked by Git, when a checkout is available.
#
# The single-cell workflow is discontinued and is intentionally excluded from
# required runtime checks even if legacy files remain in the repository.
# ==============================================================================

if [[ "$MODE" != "quick" ]]; then
    required_repo_dirs=(
        Installation
        Tools
        Tools/KrakenTools
        Tools/export2graphlan
        Tools/graphlan
        Tools/ssGSEA2.0
        Tools/ssGSEA2.0/db/msigdb
        aux_scripts
        aux_scripts/EV
        aux_scripts/GO
        aux_scripts/Heatmaper
        aux_scripts/Kraken2
        aux_scripts/analysis
        aux_scripts/analysis/differential
        aux_scripts/analysis/humann
        aux_scripts/analysis/integration
        aux_scripts/exploratory
        aux_scripts/host_reference
        aux_scripts/manifest_scripts
        aux_scripts/orgdb
        aux_scripts/ssGSEA
        update_fix
        update_fix/pvr_pkg
    )

    for repo_dir in "${required_repo_dirs[@]}"; do
        check_required_dir \
            "Repository layout" \
            "$repo_dir" \
            "$MTD_DIR/$repo_dir"
    done

    optional_repo_dirs=(
        Tutorial
        benchmark
        docs
        examples
        old_scripts_files
        test
    )

    for repo_dir in "${optional_repo_dirs[@]}"; do
        if [[ -d "$MTD_DIR/$repo_dir" ]]; then
            record PASS "Repository extras" "$repo_dir" \
                "$MTD_DIR/$repo_dir"
        else
            record SKIP "Repository extras" "$repo_dir" \
                "optional documentation/development directory not installed"
        fi
    done

    if [[ -d "$MTD_DIR/aux_scripts/analysis/single_cell" ]]; then
        record SKIP "Repository scope" "single-cell analysis" \
            "legacy directory present but intentionally excluded from MTD Explorer scope"
    else
        record SKIP "Repository scope" "single-cell analysis" \
            "discontinued and not required"
    fi

    # -------------------------------------------------------------------------
    # Core scripts and configuration.
    # -------------------------------------------------------------------------
    required_source_files=(
        Install.sh
        MTD_explorer.sh
        MTD_check_installation.sh
        Create_custom_host.sh
        Create_custom_micro.sh
        HostSpecies.csv
        conta_ls.txt

        Installation/M33262_SIVMM239.fa
        Installation/MTD.yml
        Installation/MTD_R_additions.yml
        Installation/MTD_fastp.yml
        Installation/MTD_humann.yml
        Installation/MTD_orgdb.yml
        Installation/R412.yml
        Installation/R_packages_installation.R
        Installation/check_MTD_orgdb.R
        Installation/check_R_packages_installation.R
        Installation/download_genomic_library.sh
        Installation/download_genomic_library_plasmid.sh
        Installation/halla0820.yml
        Installation/hisat2_extract_exons.py
        Installation/hisat2_extract_splice_sites.py
        Installation/pip.requirements
        Installation/py2.yml
        Installation/repair_R_packages_installation.R
        Installation/rsync_from_ncbi.pl
        Installation/rsync_from_ncbi_2.pl
        Installation/rsync_from_ncbi_archaea.pl
        Installation/rsync_from_ncbi_bacteria.pl
        Installation/rsync_from_ncbi_offline.pl

        aux_scripts/manifest_scripts/manifest.sh
        aux_scripts/manifest_scripts/manifest.virus.sh
        aux_scripts/manifest_scripts/manifest.bacteria.sh
        aux_scripts/manifest_scripts/manifest.archea.sh
        aux_scripts/manifest_scripts/manifest.plasmid.sh

        aux_scripts/Kraken2/build_multi_genome_kraken2_db.sh
        aux_scripts/Kraken2/build_nonredundant_viral_fasta.py
        aux_scripts/Kraken2/build_refseq_viral_taxid_map.py
        aux_scripts/Kraken2/build_virushost_taxid_map.py
        aux_scripts/Kraken2/check_virushost_mirror_status.sh
        aux_scripts/Kraken2/download_kraken2_taxonomy_https.sh
        aux_scripts/Kraken2/download_ncbi_taxon_genomes_manifest.sh
        aux_scripts/Kraken2/download_ncbi_taxons_from_file.sh
        aux_scripts/Kraken2/kraken2-build-download-taxonomy
        aux_scripts/Kraken2/update_virushost_mirror.sh

        aux_scripts/host_reference/build_gid_to_entrez_from_ncbi.py
        aux_scripts/host_reference/ensure_hostspecies_entry.py
        aux_scripts/host_reference/resolve_host_reference_from_csv.py

        aux_scripts/orgdb/build_all_hostspecies_orgdb.sh
        aux_scripts/orgdb/build_curated_hostspecies_csv.py
        aux_scripts/orgdb/build_gold_orgdb_from_gtf_eggnog.R
        aux_scripts/orgdb/create_annotation_package.R
        aux_scripts/orgdb/host_reference_overrides.tsv
        aux_scripts/orgdb/install_gold_orgdb_from_hostspecies.sh
        aux_scripts/orgdb/make_gene_representative_fasta.py
        aux_scripts/orgdb/update_hostspecies_common_kegg.py

        aux_scripts/ssGSEA/build_master_gmt_from_eggnog.py
        aux_scripts/ssGSEA/filter_master_gmt_for_host_gct.py
        aux_scripts/ssGSEA/make_custom_ssgsea_gmt_from_taxid_auto.sh
        aux_scripts/ssGSEA/plot_ssgsea_results.R
        aux_scripts/ssGSEA/resolve_ssgsea_go_terms.py

        aux_scripts/analysis/differential/DEG_Anno_Plot.R
        aux_scripts/analysis/differential/Normalization_afbr.R
        aux_scripts/analysis/differential/venn_diagram.R
        aux_scripts/analysis/humann/goterms.txt
        aux_scripts/analysis/humann/humann_ID_translation.R
        aux_scripts/analysis/humann/humann_ID_translation_adjusted.R
        aux_scripts/analysis/humann/koterms.txt
        aux_scripts/analysis/integration/for_halla.R
        aux_scripts/analysis/integration/gct_making.R
        aux_scripts/analysis/integration/generate_halla_heatmap.py
        aux_scripts/analysis/integration/kmeans_clustering.py
        aux_scripts/analysis/integration/pls_da_analysis.py

        aux_scripts/exploratory/MTD.find_target_taxa.py
        aux_scripts/exploratory/MTD.taxonomic_pheatmap.R
        aux_scripts/exploratory/MTD.taxonomic_stacked_bar.R
        aux_scripts/exploratory/MTD_detected_species_pie_by_phylum.R
        aux_scripts/exploratory/MTD_exploratory_alpha_diversity.R
        aux_scripts/exploratory/MTD_exploratory_beta_diversity.R
        aux_scripts/exploratory/MTD_exploratory_core_microbiome.R
        aux_scripts/exploratory/MTD_exploratory_detected_microbiome_rank.R
        aux_scripts/exploratory/MTD_exploratory_matrix_qc.R
        aux_scripts/exploratory/MTD_exploratory_prevalence_abundance.R
        aux_scripts/exploratory/MTD_exploratory_read_composition_qc.R
        aux_scripts/exploratory/MTD_extract_reads_by_detected_microbiome.py
        aux_scripts/exploratory/MTD_species_overlap_venn_euler.R

        aux_scripts/EV/EV.volcano.R
        aux_scripts/GO/GO_faceted_dotplot.R
        aux_scripts/Heatmaper/bracken_to_heatmap.R

        Tools/combine_bracken_outputs.py
        Tools/KrakenTools/combine_kreports.py
        Tools/KrakenTools/combine_mpa.py
        Tools/KrakenTools/extract_kraken_reads.py
        Tools/KrakenTools/filter_bracken.out.py
        Tools/KrakenTools/fix_unmapped.py
        Tools/KrakenTools/kreport2krona.py
        Tools/KrakenTools/kreport2mpa.py
        Tools/KrakenTools/make_kreport.py
        Tools/KrakenTools/make_ktaxonomy.py
        Tools/export2graphlan/export2graphlan.py
        Tools/graphlan/graphlan.py
        Tools/graphlan/graphlan_annotate.py
        Tools/graphlan/verify_and_correct_annotations.py
        Tools/ssGSEA2.0/config.yaml
        Tools/ssGSEA2.0/db/msigdb/c2.all.v7.5.1.symbols.gmt
        Tools/ssGSEA2.0/ssgsea-cli.R

        update_fix/Install.R.AnnotPackages.base.sh
        update_fix/Install.R.packages.MTD.sh
        update_fix/Install.R.packages.R412.sh
        update_fix/Install.R.packages.R412_optimized.sh
        update_fix/check_R_pkg.MTD.sh
        update_fix/check_R_pkg.R412.sh
        update_fix/check_R_pkg.halla0820.sh
        update_fix/hclust2.py
        update_fix/patch_halla_matplotlib.py
        update_fix/verify_NCBI_makeOrgPackageFromNCBI.sh
    )

    for source_file in "${required_source_files[@]}"; do
        check_required_file \
            "Repository source" \
            "$source_file" \
            "$MTD_DIR/$source_file"
    done

    # The manual GO replacement table is an optional override. The runtime has
    # built-in and QuickGO fallbacks when it is absent.
    if [[ -s "$MTD_DIR/aux_scripts/ssGSEA/go_replacement_manual_map.tsv" ]]; then
        record PASS "Repository source" \
            "aux_scripts/ssGSEA/go_replacement_manual_map.tsv" \
            "optional manual GO replacement map is available"
    else
        record SKIP "Repository source" \
            "aux_scripts/ssGSEA/go_replacement_manual_map.tsv" \
            "optional; resolver falls back to built-in/QuickGO mappings"
    fi

    # -------------------------------------------------------------------------
    # Syntax checks at their current paths.
    # -------------------------------------------------------------------------
    shell_scripts=(
        Install.sh
        MTD_explorer.sh
        MTD_check_installation.sh
        Create_custom_host.sh
        Create_custom_micro.sh
        Installation/download_genomic_library.sh
        Installation/download_genomic_library_plasmid.sh
        aux_scripts/manifest_scripts/manifest.sh
        aux_scripts/manifest_scripts/manifest.virus.sh
        aux_scripts/manifest_scripts/manifest.bacteria.sh
        aux_scripts/manifest_scripts/manifest.archea.sh
        aux_scripts/manifest_scripts/manifest.plasmid.sh
        aux_scripts/Kraken2/build_multi_genome_kraken2_db.sh
        aux_scripts/Kraken2/check_virushost_mirror_status.sh
        aux_scripts/Kraken2/download_kraken2_taxonomy_https.sh
        aux_scripts/Kraken2/download_ncbi_taxon_genomes_manifest.sh
        aux_scripts/Kraken2/download_ncbi_taxons_from_file.sh
        aux_scripts/Kraken2/update_virushost_mirror.sh
        aux_scripts/orgdb/build_all_hostspecies_orgdb.sh
        aux_scripts/orgdb/install_gold_orgdb_from_hostspecies.sh
        aux_scripts/ssGSEA/make_custom_ssgsea_gmt_from_taxid_auto.sh
        update_fix/Install.R.AnnotPackages.base.sh
        update_fix/Install.R.packages.MTD.sh
        update_fix/Install.R.packages.R412.sh
        update_fix/Install.R.packages.R412_optimized.sh
        update_fix/check_R_pkg.MTD.sh
        update_fix/check_R_pkg.R412.sh
        update_fix/check_R_pkg.halla0820.sh
        update_fix/verify_NCBI_makeOrgPackageFromNCBI.sh
    )

    for shell_script in "${shell_scripts[@]}"; do
        check_shell_syntax "$shell_script" "$MTD_DIR/$shell_script"
    done

    # This Kraken2 helper intentionally has no .pl extension. Its shebang is
    # the source of truth for syntax validation.
    extensionless_script="aux_scripts/Kraken2/kraken2-build-download-taxonomy"
    extensionless_path="$MTD_DIR/$extensionless_script"

    if [[ ! -f "$extensionless_path" ]]; then
        record FAIL "Installer source" \
            "syntax by shebang: $extensionless_script" \
            "missing: $extensionless_path"
    else
        extensionless_shebang="$(
            head -n 1 "$extensionless_path" 2>/dev/null || true
        )"

        case "$extensionless_shebang" in
            *perl*)
                check_perl_syntax \
                    "$extensionless_script" \
                    "$extensionless_path"
                ;;
            *bash*|*"/sh"*|*" sh")
                check_shell_syntax \
                    "$extensionless_script" \
                    "$extensionless_path"
                ;;
            *python*)
                check_python_syntax \
                    "$extensionless_script" \
                    "$extensionless_path"
                ;;
            *)
                record FAIL "Installer source" \
                    "syntax by shebang: $extensionless_script" \
                    "unsupported or missing shebang: ${extensionless_shebang:-<empty>}"
                ;;
        esac
    fi

    perl_scripts=(
        Installation/rsync_from_ncbi.pl
        Installation/rsync_from_ncbi_2.pl
        Installation/rsync_from_ncbi_archaea.pl
        Installation/rsync_from_ncbi_bacteria.pl
        Installation/rsync_from_ncbi_offline.pl
    )

    for perl_script in "${perl_scripts[@]}"; do
        check_perl_syntax "$perl_script" "$MTD_DIR/$perl_script"
    done

    # Project-owned Python 3 scripts. Vendored legacy Python utilities are
    # checked for presence above and exercised through their runtime workflows.
    python_scripts=(
        Installation/hisat2_extract_exons.py
        Installation/hisat2_extract_splice_sites.py
        aux_scripts/Kraken2/build_nonredundant_viral_fasta.py
        aux_scripts/Kraken2/build_refseq_viral_taxid_map.py
        aux_scripts/Kraken2/build_virushost_taxid_map.py
        aux_scripts/host_reference/build_gid_to_entrez_from_ncbi.py
        aux_scripts/host_reference/ensure_hostspecies_entry.py
        aux_scripts/host_reference/resolve_host_reference_from_csv.py
        aux_scripts/orgdb/build_curated_hostspecies_csv.py
        aux_scripts/orgdb/make_gene_representative_fasta.py
        aux_scripts/orgdb/update_hostspecies_common_kegg.py
        aux_scripts/ssGSEA/build_master_gmt_from_eggnog.py
        aux_scripts/ssGSEA/filter_master_gmt_for_host_gct.py
        aux_scripts/ssGSEA/resolve_ssgsea_go_terms.py
        aux_scripts/analysis/integration/generate_halla_heatmap.py
        aux_scripts/analysis/integration/kmeans_clustering.py
        aux_scripts/analysis/integration/pls_da_analysis.py
        aux_scripts/exploratory/MTD.find_target_taxa.py
        aux_scripts/exploratory/MTD_extract_reads_by_detected_microbiome.py
        Tools/combine_bracken_outputs.py
        update_fix/patch_halla_matplotlib.py
    )

    for python_script in "${python_scripts[@]}"; do
        check_python_syntax "$python_script" "$MTD_DIR/$python_script"
    done

    # -------------------------------------------------------------------------
    # Static contracts between current core scripts and the current tree.
    # -------------------------------------------------------------------------
    if grep -Fq 'MANIFEST_SCRIPTS_DIR="$dir/aux_scripts/manifest_scripts"' \
        "$MTD_DIR/Install.sh"
    then
        record PASS "Source path contract" "manifest directory" \
            "Install.sh uses aux_scripts/manifest_scripts"
    else
        record FAIL "Source path contract" "manifest directory" \
            "Install.sh does not declare the current manifest directory"
    fi

    if grep -Fq 'KRAKEN_AUX_DIR="$dir/aux_scripts/Kraken2"' \
        "$MTD_DIR/Install.sh"
    then
        record PASS "Source path contract" "Kraken2 helper directory" \
            "Install.sh uses aux_scripts/Kraken2"
    else
        record FAIL "Source path contract" "Kraken2 helper directory" \
            "Install.sh does not declare the current Kraken2 helper directory"
    fi

    if grep -Fq 'ANALYSIS_SCRIPTS_DIR="$MTDIR/aux_scripts/analysis"' \
        "$MTD_DIR/MTD_explorer.sh"
    then
        record PASS "Source path contract" "analysis script directory" \
            "MTD_explorer.sh uses aux_scripts/analysis"
    else
        record FAIL "Source path contract" "analysis script directory" \
            "MTD_explorer.sh does not declare aux_scripts/analysis"
    fi

    for active_runtime_path in \
        Tools/KrakenTools/extract_kraken_reads.py \
        Tools/KrakenTools/kreport2krona.py \
        Tools/KrakenTools/kreport2mpa.py \
        Tools/KrakenTools/combine_mpa.py \
        Tools/combine_bracken_outputs.py \
        Tools/export2graphlan/export2graphlan.py \
        Tools/graphlan/verify_and_correct_annotations.py \
        Tools/graphlan/graphlan_annotate.py \
        Tools/graphlan/graphlan.py \
        Tools/ssGSEA2.0/ssgsea-cli.R \
        Tools/ssGSEA2.0/config.yaml \
        Tools/ssGSEA2.0/db/msigdb/c2.all.v7.5.1.symbols.gmt
    do
        if grep -Fq "$active_runtime_path" "$MTD_DIR/MTD_explorer.sh"; then
            record PASS "Runtime path contract" "$active_runtime_path" \
                "referenced by MTD_explorer.sh"
        else
            record FAIL "Runtime path contract" "$active_runtime_path" \
                "expected runtime reference not found in MTD_explorer.sh"
        fi
    done

    # -------------------------------------------------------------------------
    # Full tracked-tree audit. This automatically follows future committed
    # renames and checks every tracked path, including docs and tests, without
    # hard-coding them as runtime requirements.
    # -------------------------------------------------------------------------
    if command -v git >/dev/null 2>&1 &&
       git -C "$MTD_DIR" rev-parse --is-inside-work-tree \
           >/dev/null 2>&1
    then
        tracked_list="$TMP_WORK/git_tracked_files.txt"
        missing_list="$TMP_WORK/git_missing_tracked_files.txt"

        git -C "$MTD_DIR" ls-files > "$tracked_list"
        : > "$missing_list"

        tracked_count=0
        while IFS= read -r tracked_file; do
            [[ -n "$tracked_file" ]] || continue
            tracked_count=$((tracked_count + 1))

            if [[ ! -e "$MTD_DIR/$tracked_file" &&
                  ! -L "$MTD_DIR/$tracked_file" ]]
            then
                printf '%s
' "$tracked_file" >> "$missing_list"
            fi
        done < "$tracked_list"

        missing_count="$(wc -l < "$missing_list" | tr -d '[:space:]')"

        if [[ "$missing_count" -eq 0 ]]; then
            record PASS "Git repository" "tracked-file tree" \
                "$tracked_count tracked paths are present"
        else
            record FAIL "Git repository" "tracked-file tree" \
                "$missing_count tracked path(s) missing; see $missing_list"
        fi

        for top_level in \
            Installation Tools aux_scripts update_fix docs benchmark examples test
        do
            top_count="$((
                $(awk -F/ -v top="$top_level" '$1 == top {n++} END {print n+0}' \
                    "$tracked_list")
            ))"

            if (( top_count > 0 )); then
                record PASS "Git tree inventory" "$top_level" \
                    "$top_count tracked file(s)"
            else
                record SKIP "Git tree inventory" "$top_level" \
                    "no tracked files in this checkout"
            fi
        done
    else
        record SKIP "Git repository" "tracked-file tree" \
            ".git metadata unavailable; explicit runtime/source contract was used"
    fi
fi

# ==============================================================================
# Conda environments, tools, and package runtimes
# ==============================================================================

if conda_available; then
    for env_name in base MTD_fastp MTD MTD_humann py2 halla0820 R412; do
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
        featureCounts makeblastdb blastdbcmd blastn magicblast diamond \
        emapper.py datasets STAR \
        rsem-calculate-expression nextflow parallel
    do
        check_env_command MTD "$command_name"
    done

    check_humann_runtime

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
        check_python_imports MTD_humann Cython humann numpy pysam simplejson
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

# The previous viruses4kraken.fa root-level artifact was replaced by the
# persistent, taxid-aware nonredundant viral collection. Its absence is not
# an installation failure; the current contract is validated below.

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
# Installed custom-host references
# ==============================================================================

check_installed_custom_hosts

# ==============================================================================
# Bracken
# ==============================================================================

check_required_file "Bracken DB" "read length $READ_LEN distribution" \
    "$MICRO_DB/database${READ_LEN}mers.kmer_distrib"

check_required_file "Bracken DB" "database.kraken" \
    "$MICRO_DB/database.kraken"

# ==============================================================================
# HUMAnN 3.9 / MetaPhlAn 4.1.1 dedicated database stack
# ==============================================================================

HUMANN_ROOT="$MTD_DIR/HUMAnN/ref_database"
CHOCO_DIR="$HUMANN_ROOT/chocophlan"
UNIREF_DIR="$HUMANN_ROOT/uniref"
MAPPING_DIR="$HUMANN_ROOT/utility_mapping"
METAPHLAN_DIR="$HUMANN_ROOT/metaphlan"
METAPHLAN_INDEX="mpa_vJun23_CHOCOPhlAnSGB_202403"
HUMANN_ENV_PREFIX="$CONDA_PATH/envs/MTD_humann"

check_required_dir "HUMAnN" "ChocoPhlAn directory" "$CHOCO_DIR"
check_required_dir "HUMAnN" "UniRef directory" "$UNIREF_DIR"
check_required_dir "HUMAnN" "utility mapping directory" "$MAPPING_DIR"
check_required_dir "MetaPhlAn" "database directory" "$METAPHLAN_DIR"
check_required_file "HUMAnN" "installation completion marker" \
    "$HUMANN_ROOT/.mtd_humann_databases_complete"

choco_count="$(
    find "$CHOCO_DIR" \
        -type f \
        -name '*.ffn.gz' \
        -size +0c \
        2>/dev/null | awk 'END {print NR + 0}'
)"

if (( choco_count > 0 )); then
    record PASS "HUMAnN" "ChocoPhlAn files" "$choco_count file(s)"
else
    record FAIL "HUMAnN" "ChocoPhlAn files" "none found"
fi

check_required_file "HUMAnN" "UniRef90 DIAMOND database" \
    "$UNIREF_DIR/uniref90_201901b_full.dmnd"

mapping_count="$(
    find "$MAPPING_DIR" -type f -size +0c 2>/dev/null |
        awk 'END {print NR + 0}'
)"

if (( mapping_count >= 10 )); then
    record PASS "HUMAnN" "utility mapping files" \
        "$mapping_count non-empty file(s)"
else
    record FAIL "HUMAnN" "utility mapping files" \
        "$mapping_count file(s); expected at least 10"
fi

check_required_file "HUMAnN" "UniRef90 to KO mapping" \
    "$MAPPING_DIR/map_ko_uniref90.txt.gz"
check_required_file "HUMAnN" "UniRef90 to GO mapping" \
    "$MAPPING_DIR/map_go_uniref90.txt.gz"

for metaphlan_file in \
    "${METAPHLAN_INDEX}.pkl" \
    "${METAPHLAN_INDEX}.1.bt2l" \
    "${METAPHLAN_INDEX}.2.bt2l" \
    "${METAPHLAN_INDEX}.3.bt2l" \
    "${METAPHLAN_INDEX}.4.bt2l" \
    "${METAPHLAN_INDEX}.rev.1.bt2l" \
    "${METAPHLAN_INDEX}.rev.2.bt2l" \
    "${METAPHLAN_INDEX}.nwk" \
    "${METAPHLAN_INDEX}_marker_info.txt.bz2" \
    "${METAPHLAN_INDEX}_species.txt.bz2"
do
    check_required_file "MetaPhlAn" "$metaphlan_file" \
        "$METAPHLAN_DIR/$metaphlan_file"
done

if [[ "$MODE" != "quick" ]]; then
    for bz2_file in \
        "${METAPHLAN_INDEX}_marker_info.txt.bz2" \
        "${METAPHLAN_INDEX}_species.txt.bz2"
    do
        bz2_tmp="$TMP_WORK/installed_${bz2_file}.txt"
        if capture_command "$bz2_tmp" bzip2 -t "$METAPHLAN_DIR/$bz2_file"; then
            record PASS "MetaPhlAn" "installed bzip2 integrity: $bz2_file" "valid"
        else
            record FAIL "MetaPhlAn" "installed bzip2 integrity: $bz2_file" \
                "$(compact_output "$bz2_tmp" 20)"
        fi
    done
fi

if env_exists MTD_humann; then
    humann_config_tmp="$TMP_WORK/humann_config_r9.txt"

    if capture_command "$humann_config_tmp" run_humann_isolated \
        "$HUMANN_ENV_PREFIX/bin/humann_config" --print
    then
        humann_config_ok=1

        for expected_path in \
            "$CHOCO_DIR" \
            "$UNIREF_DIR" \
            "$MAPPING_DIR"
        do
            if ! grep -Fq "$expected_path" "$humann_config_tmp"; then
                humann_config_ok=0
            fi
        done

        if (( humann_config_ok == 1 )); then
            record PASS "HUMAnN" "configured database paths" \
                "all paths point inside $HUMANN_ROOT"
        else
            record FAIL "HUMAnN" "configured database paths" \
                "one or more expected paths are absent: $(compact_output "$humann_config_tmp" 20)"
        fi
    else
        record FAIL "HUMAnN" "configured database paths" \
            "$(compact_output "$humann_config_tmp" 20)"
    fi
fi

if [[ "$MODE" == "deep" && -x "$HUMANN_ENV_PREFIX/bin/humann_test" ]]; then
    humann_test_tmp="$TMP_WORK/humann_test_r9.txt"

    if capture_command "$humann_test_tmp" timeout 300 env \
        -u PYTHONPATH \
        -u PYTHONHOME \
        PYTHONNOUSERSITE=1 \
        PATH="$HUMANN_ENV_PREFIX/bin:$PATH" \
        "$HUMANN_ENV_PREFIX/bin/humann_test" &&
       grep -Eq '^Ran [0-9]+ tests' "$humann_test_tmp" &&
       grep -Eq '^OK$' "$humann_test_tmp"
    then
        record PASS "HUMAnN" "unit test suite" \
            "$(grep -E '^Ran [0-9]+ tests|^OK$' "$humann_test_tmp" | tr '
' ' ')"
    else
        record FAIL "HUMAnN" "unit test suite" \
            "$(compact_output "$humann_test_tmp" 30)"
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
    check_curated_viral_reference_cache

    for humann_archive in \
        full_mapping_v201901b.tar.gz \
        uniref90_annotated_v201901b_full.tar.gz \
        full_chocophlan.v201901_v31.tar.gz
    do
        check_required_file "Installation cache" \
            "HUMAnN/$humann_archive" \
            "$OFFLINE_DIR/HUMAnN/$humann_archive"
    done


    METAPHLAN_CACHE_INDEX="mpa_vJun23_CHOCOPhlAnSGB_202403"
    METAPHLAN_CACHE_DIR="$OFFLINE_DIR/HUMAnN/metaphlan_vJun23_202403_archives"

    for metaphlan_cache_file in \
        "${METAPHLAN_CACHE_INDEX}.tar" \
        "${METAPHLAN_CACHE_INDEX}.md5" \
        "${METAPHLAN_CACHE_INDEX}_bt2.tar" \
        "${METAPHLAN_CACHE_INDEX}_bt2.md5" \
        "${METAPHLAN_CACHE_INDEX}.nwk" \
        "${METAPHLAN_CACHE_INDEX}_marker_info.txt.bz2" \
        "${METAPHLAN_CACHE_INDEX}_species.txt.bz2"
    do
        check_required_file "Installation cache" \
            "MetaPhlAn/$metaphlan_cache_file" \
            "$METAPHLAN_CACHE_DIR/$metaphlan_cache_file"
    done

    check_required_file "Installation cache" \
        "MetaPhlAn completion marker" \
        "$METAPHLAN_CACHE_DIR/.metaphlan_vJun23_202403_cache_complete"

    if [[ "$MODE" != "quick" && \
          -s "$METAPHLAN_CACHE_DIR/${METAPHLAN_CACHE_INDEX}.nwk" ]]
    then
        metaphlan_newick_tail="$(
            tail -c 1024 \
                "$METAPHLAN_CACHE_DIR/${METAPHLAN_CACHE_INDEX}.nwk" \
                2>/dev/null | tr -d '\r\n[:space:]'
        )"

        if [[ -n "$metaphlan_newick_tail" && \
              "${metaphlan_newick_tail: -1}" == ";" ]]
        then
            record PASS "Installation cache" "MetaPhlAn taxonomy tree" \
                "Newick terminator detected"
        else
            record FAIL "Installation cache" "MetaPhlAn taxonomy tree" \
                "file does not end with ';'"
        fi
    fi

    if [[ "$MODE" == "deep" ]]; then
        archive_integrity "MetaPhlAn main archive" \
            "$METAPHLAN_CACHE_DIR/${METAPHLAN_CACHE_INDEX}.tar"
        archive_integrity "MetaPhlAn Bowtie2 archive" \
            "$METAPHLAN_CACHE_DIR/${METAPHLAN_CACHE_INDEX}_bt2.tar"
        archive_integrity "MetaPhlAn marker information" \
            "$METAPHLAN_CACHE_DIR/${METAPHLAN_CACHE_INDEX}_marker_info.txt.bz2"
        archive_integrity "MetaPhlAn species information" \
            "$METAPHLAN_CACHE_DIR/${METAPHLAN_CACHE_INDEX}_species.txt.bz2"

        main_md5="$(
            md5sum "$METAPHLAN_CACHE_DIR/${METAPHLAN_CACHE_INDEX}.tar" \
                2>/dev/null | awk '{print $1}'
        )"
        bt2_md5="$(
            md5sum "$METAPHLAN_CACHE_DIR/${METAPHLAN_CACHE_INDEX}_bt2.tar" \
                2>/dev/null | awk '{print $1}'
        )"

        if [[ "$main_md5" == "d985de75a217cd319e721863f68e7d33" ]]; then
            record PASS "Installation cache" "MetaPhlAn main archive MD5" \
                "$main_md5"
        else
            record FAIL "Installation cache" "MetaPhlAn main archive MD5" \
                "observed=${main_md5:-unavailable}"
        fi

        if [[ "$bt2_md5" == "8caae86b4d2931416cbdbb92f5985cef" ]]; then
            record PASS "Installation cache" "MetaPhlAn Bowtie2 archive MD5" \
                "$bt2_md5"
        else
            record FAIL "Installation cache" "MetaPhlAn Bowtie2 archive MD5" \
                "observed=${bt2_md5:-unavailable}"
        fi

        if (
            cd "$METAPHLAN_CACHE_DIR" &&
            md5sum -c "${METAPHLAN_CACHE_INDEX}.md5" >/dev/null 2>&1 &&
            md5sum -c "${METAPHLAN_CACHE_INDEX}_bt2.md5" >/dev/null 2>&1
        ); then
            record PASS "Installation cache" "MetaPhlAn official MD5 manifests" \
                "both manifests validated"
        else
            record FAIL "Installation cache" "MetaPhlAn official MD5 manifests" \
                "one or both manifests failed"
        fi
    fi

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
        # Virus-Host integrity is validated by the r9.3 curated stack.
        :

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
    echo "Host TaxID:  ${HOST_TAXID:-auto-detect}"
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
