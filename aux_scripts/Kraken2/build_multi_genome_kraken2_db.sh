#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# build_multi_genome_kraken2_db.sh
#
# Build a custom Kraken2 database from multiple genome FASTA files.
#
# Input genome_info file:
#   taxid<TAB>/path/to/genome.fasta
#
# Example:
#   40233   /home/me/MTD/genomes/Carollia_perspicillata.fa
#   89673   /home/me/MTD/genomes/Phyllostomus_discolor.fa
#   59463   /home/me/MTD/genomes/Myotis_lucifugus.fa
#   9430    /home/me/MTD/genomes/Desmodus_rotundus.fa
#
# The script:
#   1) Reads taxid + FASTA path from --genome_info
#   2) Adds kraken:taxid|TAXID| to each FASTA header
#   3) Creates a DB folder named kraken2DB_TAXID1_TAXID2_...
#   4) Copies taxonomy from --kraken_offline if provided
#   5) Otherwise downloads taxonomy using kraken2-build --download-taxonomy --use-ftp
#   6) Checks whether all supplied taxids exist in taxonomy/nodes.dmp
#   7) Retrieves scientific names from taxonomy/names.dmp
#   8) Shows species names beside taxids in the terminal
#   9) Adds all prepared FASTAs to the Kraken2 library
#  10) Builds the Kraken2 database
# ============================================================

# ============================================================
# Pretty terminal output
# ============================================================

if [[ -t 1 ]]; then
    BOLD="$(tput bold || true)"
    DIM="$(tput dim || true)"
    ITALIC="$(tput sitm || true)"
    RESET_ITALIC="$(tput ritm || true)"
    RESET="$(tput sgr0 || true)"

    GREEN="$(tput setaf 2 || true)"
    RED="$(tput setaf 1 || true)"
    YELLOW="$(tput setaf 3 || true)"
    BLUE="$(tput setaf 4 || true)"
    MAGENTA="$(tput setaf 5 || true)"
    CYAN="$(tput setaf 6 || true)"
else
    BOLD=""
    DIM=""
    ITALIC=""
    RESET_ITALIC=""
    RESET=""

    GREEN=""
    RED=""
    YELLOW=""
    BLUE=""
    MAGENTA=""
    CYAN=""
fi

info() {
    echo "${CYAN}${BOLD}[INFO]${RESET} $*"
}

ok() {
    echo "${GREEN}${BOLD}[OK]${RESET} $*"
}

warn() {
    echo "${YELLOW}${BOLD}[WARNING]${RESET} $*"
}

err() {
    echo "${RED}${BOLD}[ERROR]${RESET} $*" >&2
}

step() {
    echo
    echo "${BLUE}${BOLD}============================================================${RESET}"
    echo "${BLUE}${BOLD}$*${RESET}"
    echo "${BLUE}${BOLD}============================================================${RESET}"
}

species_name_fmt() {
    local name="$1"
    echo "${ITALIC}${name}${RESET}"
}

taxid_fmt() {
    local taxid="$1"
    echo "${GREEN}${BOLD}${taxid}${RESET}"
}

# ============================================================
# Help
# ============================================================

usage() {
    cat <<EOF

${BOLD}Usage:${RESET}
  $0 --genome_info genomes.txt [options]

${BOLD}Required:${RESET}
  --genome_info FILE          Tab-delimited file with two columns:
                              taxid<TAB>/path/to/fasta

${BOLD}Optional:${RESET}
  --outdir DIR                Output parent directory.
                              Default: current directory

  --threads INT               Number of threads.
                              Default: auto-detect with nproc minus 1

  --kraken_offline DIR        Existing Kraken2 DB folder or taxonomy folder.
                              If provided, taxonomy will be copied from it.
                              Accepted:
                                /path/to/old_db
                                /path/to/old_db/taxonomy

  --db_name NAME              Custom DB folder name.
                              Default: kraken2DB_TAXID1_TAXID2_...

  --min_len INT               Minimum sequence length to keep.
                              Default: 0

  --force                     Remove existing DB folder before rebuilding.

  --no_build                  Prepare files and add library, but do not run build.

  --ncbi_fallback             If a taxid is missing from local names.dmp,
                              try to retrieve the name using NCBI datasets CLI.
                              Requires internet and the datasets command.
                              Important: Kraken2 still needs the taxid inside
                              local taxonomy/nodes.dmp.
  --micro-db                  Mark this DB as a microbiome/target DB for MTD.
                              This automatically enables Bracken database build.

  --build-bracken             Build Bracken k-mer distribution file after
                              Kraken2 DB build.

  --bracken-read-len INT      Read length for Bracken build.
                              This creates databaseINTmers.kmer_distrib.
                              Default: 75

  --bracken-kmer-len INT      Kraken2 k-mer length used by Bracken.
                              Default: 35
  -h, --help                  Show this help message.

${BOLD}Example with offline taxonomy:${RESET}
  $0 \\
    --genome_info genomes.txt \\
    --kraken_offline /home/me/MTD/kraken2DB_host_carollia_mCarPer1.2 \\
    --outdir /home/me/MTD

${BOLD}Example with manual threads:${RESET}
  $0 \\
    --genome_info genomes.txt \\
    --threads 20 \\
    --kraken_offline /home/me/MTD/kraken2DB_host_carollia_mCarPer1.2 \\
    --outdir /home/me/MTD

${BOLD}Example downloading taxonomy with FTP:${RESET}
  $0 \\
    --genome_info genomes.txt \\
    --outdir /home/me/MTD
${BOLD}Example building a Trematoda DB for MTD micro/target step:${RESET}
  $0 \\
    --genome_info trematoda_genome_info.tsv \\
    --outdir /home/me/MTD \\
    --db_name kraken2DB_Trematoda \\
    --kraken_offline /home/me/MTD/kraken2DB_micro \\
    --threads 20 \\
    --micro-db \\
    --bracken-read-len 75
EOF
}

# ============================================================
# Defaults
# ============================================================

GENOME_INFO=""
OUTDIR="$(pwd)"

# Auto-detect threads if --threads is not provided
if command -v nproc >/dev/null 2>&1; then
    detected_threads="$(nproc)"
    if [[ "$detected_threads" -gt 1 ]]; then
        THREADS=$((detected_threads - 2))
    else
        THREADS=1
    fi
else
    THREADS=1
fi

KRAKEN_OFFLINE=""
DB_NAME=""
MIN_LEN=0
FORCE=0
NO_BUILD=0
NCBI_FALLBACK=0

# Optional Bracken build.
# Useful when this Kraken2 DB will be used as DB_micro / target DB in MTD.
BUILD_BRACKEN=0
DB_ROLE="custom"
BRACKEN_READ_LEN=75
BRACKEN_KMER_LEN=35

# ============================================================
# Parse arguments
# ============================================================

while [[ $# -gt 0 ]]; do
    case "$1" in
        --genome_info)
            GENOME_INFO="$2"
            shift 2
            ;;
        --outdir)
            OUTDIR="$2"
            shift 2
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        --kraken_offline)
            KRAKEN_OFFLINE="$2"
            shift 2
            ;;
        --db_name)
            DB_NAME="$2"
            shift 2
            ;;
        --min_len)
            MIN_LEN="$2"
            shift 2
            ;;
        --force)
            FORCE=1
            shift
            ;;
        --no_build)
            NO_BUILD=1
            shift
            ;;
        --ncbi_fallback)
            NCBI_FALLBACK=1
            shift
            ;;
        --micro-db|--micro_db|--db-role-micro)
            DB_ROLE="micro"
            BUILD_BRACKEN=1
            shift
            ;;
        --build-bracken|--build_bracken)
            BUILD_BRACKEN=1
            shift
            ;;
        --bracken-read-len|--bracken-read-length)
            BRACKEN_READ_LEN="$2"
            shift 2
            ;;
        --bracken-read-len=*|--bracken-read-length=*)
            BRACKEN_READ_LEN="${1#*=}"
            shift
            ;;
        --bracken-kmer-len|--bracken-kmer-length)
            BRACKEN_KMER_LEN="$2"
            shift 2
            ;;
        --bracken-kmer-len=*|--bracken-kmer-length=*)
            BRACKEN_KMER_LEN="${1#*=}"
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            err "Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done


# ============================================================
# MTD dedicated Kraken2/Bracken builder runtime
# ============================================================
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
MTD_ROOT="$(cd -- "$SCRIPT_DIR/../.." && pwd -P)"
KRAKEN2_ENV_NAME="${MTD_KRAKEN2_ENV_NAME:-MTD_kraken2}"
CONDA_ROOT="${MTD_CONDA_ROOT:-}"
KRAKEN2_ENV_DIR=""
KRAKEN2_BIN=""
KRAKEN2_BUILD_BIN=""
KRAKEN2_INSPECT_BIN=""
BRACKEN_BUILD_BIN=""
KRAKEN2_VERSION=""
BRACKEN_PACKAGE_VERSION=""

resolve_dedicated_kraken_runtime() {
    if [[ -z "$CONDA_ROOT" && -s "$MTD_ROOT/condaPath" ]]; then
        IFS= read -r CONDA_ROOT < "$MTD_ROOT/condaPath"
        CONDA_ROOT="${CONDA_ROOT%$'\r'}"
    fi
    [[ -n "$CONDA_ROOT" ]] || CONDA_ROOT="$HOME/miniconda3"

    KRAKEN2_ENV_DIR="$CONDA_ROOT/envs/$KRAKEN2_ENV_NAME"
    KRAKEN2_BIN="$KRAKEN2_ENV_DIR/bin/kraken2"
    KRAKEN2_BUILD_BIN="$KRAKEN2_ENV_DIR/bin/kraken2-build"
    KRAKEN2_INSPECT_BIN="$KRAKEN2_ENV_DIR/bin/kraken2-inspect"
    BRACKEN_BUILD_BIN="$KRAKEN2_ENV_DIR/bin/bracken-build"

    local required
    for required in "$KRAKEN2_BIN" "$KRAKEN2_BUILD_BIN" "$KRAKEN2_INSPECT_BIN"; do
        if [[ ! -x "$required" ]]; then
            err "Dedicated Kraken2 executable not found: $required"
            exit 1
        fi
    done

    KRAKEN2_VERSION="$("$KRAKEN2_BIN" --version 2>/dev/null | awk 'NR == 1 { print $3 }')"
    if [[ "$KRAKEN2_VERSION" != "2.17.1" ]]; then
        err "Expected Kraken2 2.17.1; observed: ${KRAKEN2_VERSION:-unknown}"
        exit 1
    fi

    if [[ "$BUILD_BRACKEN" -eq 1 ]]; then
        [[ -x "$BRACKEN_BUILD_BIN" ]] || { err "Dedicated bracken-build not found: $BRACKEN_BUILD_BIN"; exit 1; }
        BRACKEN_PACKAGE_VERSION="$(find "$KRAKEN2_ENV_DIR/conda-meta" -maxdepth 1 -type f -name 'bracken-3.1p1-*.json' -printf '3.1p1\n' -quit 2>/dev/null)"
        if [[ "$BRACKEN_PACKAGE_VERSION" != "3.1p1" ]]; then
            err "Expected Bracken Conda package 3.1p1 in $KRAKEN2_ENV_DIR"
            exit 1
        fi
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

# ============================================================
# Basic checks
# ============================================================

if [[ -z "$GENOME_INFO" ]]; then
    err "--genome_info is required."
    usage
    exit 1
fi

if [[ ! -s "$GENOME_INFO" ]]; then
    err "genome_info file not found or empty: $GENOME_INFO"
    exit 1
fi

resolve_dedicated_kraken_runtime

if ! [[ "$THREADS" =~ ^[0-9]+$ ]]; then
    err "--threads must be an integer. Got: $THREADS"
    exit 1
fi

if [[ "$THREADS" -lt 1 ]]; then
    err "--threads must be >= 1. Got: $THREADS"
    exit 1
fi

if ! [[ "$MIN_LEN" =~ ^[0-9]+$ ]]; then
    err "--min_len must be an integer. Got: $MIN_LEN"
    exit 1
fi
if ! [[ "$BRACKEN_READ_LEN" =~ ^[0-9]+$ ]] || [[ "$BRACKEN_READ_LEN" -lt 1 ]]; then
    err "--bracken-read-len must be a positive integer. Got: $BRACKEN_READ_LEN"
    exit 1
fi

if ! [[ "$BRACKEN_KMER_LEN" =~ ^[0-9]+$ ]] || [[ "$BRACKEN_KMER_LEN" -lt 1 ]]; then
    err "--bracken-kmer-len must be a positive integer. Got: $BRACKEN_KMER_LEN"
    exit 1
fi

mkdir -p "$OUTDIR"

# ============================================================
# Parse taxids and create default DB name
# ============================================================

taxid_list=$(
    awk '
        BEGIN{FS="[ \t]+"}
        NF >= 2 && $1 !~ /^#/ {
            print $1
        }
    ' "$GENOME_INFO" | awk '!seen[$0]++' | paste -sd "_" -
)

if [[ -z "$taxid_list" ]]; then
    err "No valid taxids found in genome_info file."
    exit 1
fi

if [[ -z "$DB_NAME" ]]; then
    DB_NAME="kraken2DB_${taxid_list}"
fi

DB_DIR="${OUTDIR%/}/${DB_NAME}"
PREP_DIR="$DB_DIR/prepared_fastas"
LOG_DIR="$DB_DIR/logs"

if [[ -d "$DB_DIR" && "$FORCE" -eq 1 ]]; then
    warn "--force used. Removing existing DB folder: $DB_DIR"
    rm -rf "$DB_DIR"
fi

if [[ -d "$DB_DIR" && "$FORCE" -eq 0 ]]; then
    err "DB folder already exists: $DB_DIR"
    echo "Use --force to remove and rebuild it, or use --db_name with another name."
    exit 1
fi

mkdir -p "$DB_DIR" "$PREP_DIR" "$LOG_DIR"

echo
echo "${GREEN}${BOLD}============================================================${RESET}"
echo "${GREEN}${BOLD} Multi-genome Kraken2 DB builder${RESET}"
echo "${GREEN}${BOLD}============================================================${RESET}"
echo "${BOLD}Genome info:${RESET}      $GENOME_INFO"
echo "${BOLD}Output dir:${RESET}       $OUTDIR"
echo "${BOLD}DB name:${RESET}          $DB_NAME"
echo "${BOLD}DB dir:${RESET}           $DB_DIR"
echo "${BOLD}Threads:${RESET}          $THREADS"
echo "${BOLD}Kraken2 runtime:${RESET}   $KRAKEN2_VERSION"
echo "${BOLD}Kraken2 build:${RESET}     $KRAKEN2_BUILD_BIN"
if [[ "$BUILD_BRACKEN" -eq 1 ]]; then
    echo "${BOLD}Bracken package:${RESET}  $BRACKEN_PACKAGE_VERSION"
    echo "${BOLD}Bracken build:${RESET}    $BRACKEN_BUILD_BIN"
fi
echo "${BOLD}Taxids:${RESET}           ${GREEN}${BOLD}${taxid_list}${RESET}"
echo "${BOLD}Min length:${RESET}       $MIN_LEN"
echo "${BOLD}Offline taxonomy:${RESET} ${KRAKEN_OFFLINE:-none}"
echo "${BOLD}NCBI fallback:${RESET}    $NCBI_FALLBACK"
echo "${BOLD}DB role:${RESET}          $DB_ROLE"
echo "${BOLD}Build Bracken:${RESET}    $BUILD_BRACKEN"
echo "${BOLD}Bracken read length:${RESET} $BRACKEN_READ_LEN"
echo "${BOLD}Bracken k-mer length:${RESET} $BRACKEN_KMER_LEN"
echo "${GREEN}${BOLD}============================================================${RESET}"
echo

# ============================================================
# Functions
# ============================================================

prepare_fasta() {
    local taxid="$1"
    local fasta="$2"
    local out_fasta="$3"
    local min_len="$4"

    info "Preparing FASTA for TaxID $(taxid_fmt "$taxid")"
    echo "       ${DIM}Input:${RESET}  $fasta"
    echo "       ${DIM}Output:${RESET} $out_fasta"

    if [[ ! -s "$fasta" ]]; then
        err "FASTA not found or empty: $fasta"
        exit 1
    fi

    local reader
    if [[ "$fasta" == *.gz ]]; then
        reader="gzip -cd"
    else
        reader="cat"
    fi

    $reader "$fasta" | awk -v taxid="$taxid" -v min_len="$min_len" '
        function flush_record() {
            if (header != "") {
                if (length(seq) >= min_len) {
                    n++

                    clean_header = header
                    sub(/^kraken:taxid\|[0-9]+\|/, "", clean_header)

                    print ">kraken:taxid|" taxid "|" clean_header
                    print seq
                }
            }
        }

        BEGIN {
            header = ""
            seq = ""
            n = 0
        }

        /^>/ {
            flush_record()
            header = substr($0, 2)
            gsub(/\r/, "", header)
            seq = ""
            next
        }

        {
            line = toupper($0)

            # Keep common IUPAC nucleotide symbols.
            # Replace unusual characters with N.
            gsub(/[^ACGTRYSWKMBDHVN]/, "N", line)

            seq = seq line
        }

        END {
            flush_record()
            print "[INFO] sequences_kept=" n > "/dev/stderr"
        }
    ' > "$out_fasta"
}

get_tax_name_local() {
    local taxid="$1"
    local names_dmp="$2"

    awk -v t="$taxid" '
        BEGIN { FS="\\|"; OFS="\t" }

        {
            id=$1
            name=$2
            class=$4

            gsub(/^[ \t]+|[ \t]+$/, "", id)
            gsub(/^[ \t]+|[ \t]+$/, "", name)
            gsub(/^[ \t]+|[ \t]+$/, "", class)

            if (id == t && class == "scientific name") {
                print name
                exit
            }
        }
    ' "$names_dmp"
}

get_tax_name_ncbi() {
    local taxid="$1"

    if ! command -v datasets >/dev/null 2>&1; then
        echo ""
        return 0
    fi

    # Requires internet.
    # Parse JSON simply without jq dependency.
    datasets summary taxonomy taxon "$taxid" 2>/dev/null \
        | grep -m 1 '"sci_name"' \
        | sed 's/.*"sci_name"[[:space:]]*:[[:space:]]*"//; s/".*//'
}

# ============================================================
# Prepare all FASTAs
# ============================================================

step "[STEP] Preparing FASTA files"

prepared_manifest="$DB_DIR/prepared_fastas.tsv"
: > "$prepared_manifest"

line_number=0

while IFS=$'\t ' read -r taxid fasta extra; do
    line_number=$((line_number + 1))

    # Skip blank lines and comments
    [[ -z "${taxid:-}" ]] && continue
    [[ "$taxid" =~ ^# ]] && continue

    if [[ -z "${fasta:-}" ]]; then
        err "Invalid line $line_number in $GENOME_INFO"
        echo "Expected:"
        echo "  taxid<TAB>/path/to/fasta"
        exit 1
    fi

    if [[ ! "$taxid" =~ ^[0-9]+$ ]]; then
        err "Invalid taxid on line $line_number: $taxid"
        exit 1
    fi

    if [[ ! -s "$fasta" ]]; then
        err "FASTA from genome_info does not exist or is empty:"
        echo "Line: $line_number"
        echo "TaxID: $taxid"
        echo "FASTA: $fasta"
        exit 1
    fi

    fasta_base=$(basename "$fasta")
    fasta_base=${fasta_base%.gz}
    fasta_base=${fasta_base%.fasta}
    fasta_base=${fasta_base%.fa}
    fasta_base=${fasta_base%.fna}
    fasta_base=$(echo "$fasta_base" | sed 's/[^A-Za-z0-9_.-]/_/g')

    out_fasta="$PREP_DIR/${taxid}_${fasta_base}.kraken.fa"

    prepare_fasta "$taxid" "$fasta" "$out_fasta" "$MIN_LEN" \
        2> "$LOG_DIR/${taxid}_${fasta_base}.prepare.log"

    if [[ ! -s "$out_fasta" ]]; then
        err "Prepared FASTA is empty: $out_fasta"
        echo "Check min length setting or input FASTA."
        exit 1
    fi

    echo -e "${taxid}\t${fasta}\t${out_fasta}" >> "$prepared_manifest"

done < "$GENOME_INFO"

if [[ ! -s "$prepared_manifest" ]]; then
    err "No FASTA files were prepared."
    exit 1
fi

echo
ok "Prepared FASTA manifest: $prepared_manifest"
column -s $'\t' -t "$prepared_manifest" || cat "$prepared_manifest"
echo

# ============================================================
# Taxonomy setup
# ============================================================

step "[STEP] Taxonomy setup"

if [[ -n "$KRAKEN_OFFLINE" ]]; then

    if [[ -d "$KRAKEN_OFFLINE/taxonomy" ]]; then
        TAX_SRC="$KRAKEN_OFFLINE/taxonomy"
    elif [[ -d "$KRAKEN_OFFLINE" && -s "$KRAKEN_OFFLINE/names.dmp" && -s "$KRAKEN_OFFLINE/nodes.dmp" ]]; then
        TAX_SRC="$KRAKEN_OFFLINE"
    else
        err "--kraken_offline must be either:"
        echo "  1) an existing Kraken2 DB folder containing taxonomy/"
        echo "  2) a taxonomy folder containing names.dmp and nodes.dmp"
        echo
        echo "Provided:"
        echo "$KRAKEN_OFFLINE"
        exit 1
    fi

    info "Copying taxonomy from: $TAX_SRC"
    mkdir -p "$DB_DIR/taxonomy"
    rsync -a "$TAX_SRC/" "$DB_DIR/taxonomy/"

else
    info "No --kraken_offline provided."
    info "Downloading taxonomy using kraken2-build --download-taxonomy --use-ftp"

    run_kraken2_build --download-taxonomy --use-ftp --db "$DB_DIR" \
        > "$LOG_DIR/download_taxonomy.log" 2>&1
fi

if [[ ! -s "$DB_DIR/taxonomy/names.dmp" || ! -s "$DB_DIR/taxonomy/nodes.dmp" ]]; then
    err "Taxonomy setup failed. Missing names.dmp or nodes.dmp."
    echo "Check:"
    echo "$DB_DIR/taxonomy"
    echo
    echo "If download failed, check:"
    echo "$LOG_DIR/download_taxonomy.log"
    exit 1
fi

ok "Taxonomy ready: $DB_DIR/taxonomy"

# ============================================================
# Check taxids and show species names
# ============================================================

step "[STEP] Checking taxids and retrieving scientific names"

tax_check="$DB_DIR/taxid_species_check.tsv"
echo -e "taxid\tstatus\tscientific_name\tfasta" > "$tax_check"

missing_taxids=0

while IFS=$'\t' read -r taxid original_fasta prepared_fasta; do
    status="OK"
    scientific_name=""

    pretty_taxid="$(taxid_fmt "$taxid")"

    if grep -wq "^${taxid}" "$DB_DIR/taxonomy/nodes.dmp"; then
        scientific_name=$(get_tax_name_local "$taxid" "$DB_DIR/taxonomy/names.dmp")

        if [[ -z "$scientific_name" ]]; then
            scientific_name="FOUND_IN_NODES_BUT_NAME_NOT_FOUND"
            warn "TaxID ${pretty_taxid} found in nodes.dmp, but scientific name was not found in names.dmp"
        else
            pretty_species="$(species_name_fmt "$scientific_name")"
            ok "TaxID ${pretty_taxid} = ${pretty_species}"
        fi
    else
        status="MISSING_IN_LOCAL_TAXONOMY"
        missing_taxids=$((missing_taxids + 1))

        if [[ "$NCBI_FALLBACK" -eq 1 ]]; then
            ncbi_name=$(get_tax_name_ncbi "$taxid")
            if [[ -n "$ncbi_name" ]]; then
                scientific_name="$ncbi_name"
                status="MISSING_LOCAL_FOUND_NCBI"

                pretty_species="$(species_name_fmt "$scientific_name")"
                warn "TaxID ${pretty_taxid} missing in local taxonomy, but NCBI fallback found: ${pretty_species}"
            else
                scientific_name="NOT_FOUND"
                err "TaxID ${pretty_taxid} missing in local taxonomy and not found by NCBI fallback"
            fi
        else
            scientific_name="NOT_FOUND"
            err "TaxID ${pretty_taxid} missing in local taxonomy"
        fi
    fi

    echo -e "${taxid}\t${status}\t${scientific_name}\t${original_fasta}" >> "$tax_check"

done < "$prepared_manifest"

echo
info "Taxid/species check table:"
column -s $'\t' -t "$tax_check" || cat "$tax_check"

echo
ok "Taxid/species check saved to: $tax_check"

if [[ "$missing_taxids" -gt 0 ]]; then
    echo
    err "One or more taxids are missing from local taxonomy/nodes.dmp."
    echo "This usually means:"
    echo "  1) the taxid is wrong, or"
    echo "  2) the copied taxonomy is old, or"
    echo "  3) the taxonomy folder is incomplete."
    echo
    echo "Fix the taxid or use a newer taxonomy folder."
    echo "You can also run with --ncbi_fallback to display possible names from NCBI,"
    echo "but Kraken2 still needs the taxid inside local taxonomy/nodes.dmp."
    exit 1
fi

echo
ok "All supplied taxids exist in local taxonomy."

# ============================================================
# Add FASTAs to Kraken2 library
# ============================================================

step "[STEP] Adding FASTAs to Kraken2 library"

while IFS=$'\t' read -r taxid original_fasta prepared_fasta; do
    scientific_name=$(awk -F '\t' -v t="$taxid" '$1==t {print $3; exit}' "$tax_check")

    pretty_taxid="$(taxid_fmt "$taxid")"
    pretty_species="$(species_name_fmt "${scientific_name:-unknown}")"

    echo "${MAGENTA}${BOLD}[ADD]${RESET} TaxID ${pretty_taxid} | ${pretty_species}"
    echo "      ${DIM}Original:${RESET} $original_fasta"
    echo "      ${DIM}Prepared:${RESET} $prepared_fasta"

    run_kraken2_build \
        --add-to-library "$prepared_fasta" \
        --db "$DB_DIR" \
        > "$LOG_DIR/add_${taxid}_$(basename "$prepared_fasta").log" 2>&1

done < "$prepared_manifest"

echo
ok "All FASTAs added to library."

# ============================================================
# Build database
# ============================================================

if [[ "$NO_BUILD" -eq 1 ]]; then
    warn "--no_build used. Skipping kraken2-build --build."
    info "DB folder prepared at: $DB_DIR"
    echo
    echo "${GREEN}${BOLD}============================================================${RESET}"
    echo "${GREEN}${BOLD}[DONE]${RESET}"
    echo "${BOLD}Kraken2 DB prepared:${RESET}"
    echo "$DB_DIR"
    echo "${GREEN}${BOLD}============================================================${RESET}"
    exit 0
fi

step "[STEP] Building Kraken2 DB"

run_kraken2_build \
    --build \
    --threads "$THREADS" \
    --db "$DB_DIR" \
    > "$LOG_DIR/build.log" 2>&1

ok "Kraken2 DB built successfully: $DB_DIR"

# ============================================================
# Optional Bracken build
# ============================================================

if [[ "$BUILD_BRACKEN" -eq 1 ]]; then
    step "[STEP] Building Bracken database"

    BRACKEN_DIST="$DB_DIR/database${BRACKEN_READ_LEN}mers.kmer_distrib"

    echo "${BOLD}Bracken DB:${RESET} $DB_DIR"
    echo "${BOLD}Bracken read length:${RESET} $BRACKEN_READ_LEN"
    echo "${BOLD}Bracken k-mer length:${RESET} $BRACKEN_KMER_LEN"
    echo "${BOLD}Expected distribution file:${RESET} $BRACKEN_DIST"

    if [[ -s "$BRACKEN_DIST" ]]; then
        ok "Bracken distribution file already exists: $BRACKEN_DIST"
    else
        run_bracken_build \
            -d "$DB_DIR" \
            -t "$THREADS" \
            -k "$BRACKEN_KMER_LEN" \
            -l "$BRACKEN_READ_LEN" \
            > "$LOG_DIR/bracken-build_L${BRACKEN_READ_LEN}_K${BRACKEN_KMER_LEN}.log" 2>&1

        if [[ ! -s "$BRACKEN_DIST" ]]; then
            err "Bracken build finished, but expected distribution file was not created."
            echo "Expected:"
            echo "  $BRACKEN_DIST"
            echo
            echo "Check log:"
            echo "  $LOG_DIR/bracken-build_L${BRACKEN_READ_LEN}_K${BRACKEN_KMER_LEN}.log"
            exit 1
        fi

        ok "Bracken distribution file created: $BRACKEN_DIST"
    fi
else
    info "Bracken build not requested. Use --micro-db or --build-bracken to enable it."
fi

# ============================================================
# Inspect database
# ============================================================

step "[STEP] Inspecting database"

run_kraken2_inspect --db "$DB_DIR" > "$LOG_DIR/kraken2_inspect.txt"

ok "kraken2-inspect output saved to: $LOG_DIR/kraken2_inspect.txt"

echo
info "Taxids/species found in inspect:"

while IFS=$'\t' read -r taxid status scientific_name fasta; do
    [[ "$taxid" == "taxid" ]] && continue

    pretty_taxid="$(taxid_fmt "$taxid")"
    pretty_species="$(species_name_fmt "$scientific_name")"

    echo "${BLUE}${BOLD}------------------------------------------------------------${RESET}"
    echo "${BOLD}TaxID:${RESET} $pretty_taxid"
    echo "${BOLD}Name:${RESET}  $pretty_species"

    grep -w "$taxid" "$LOG_DIR/kraken2_inspect.txt" | head || true
done < "$tax_check"

# ============================================================
# Final message
# ============================================================

echo
echo "${GREEN}${BOLD}============================================================${RESET}"
echo "${GREEN}${BOLD}[DONE] Multi-genome Kraken2 DB ready${RESET}"
echo "${GREEN}${BOLD}============================================================${RESET}"
echo "${BOLD}Kraken2 DB:${RESET}"
echo "$DB_DIR"
echo
echo "${BOLD}Taxid/species table:${RESET}"
echo "$tax_check"
echo
echo "${BOLD}Prepared FASTA manifest:${RESET}"
echo "$prepared_manifest"
echo
echo "${BOLD}Logs:${RESET}"
echo "$LOG_DIR"
echo
echo "${BOLD}Use in your MTD script as:${RESET}"

if [[ "$DB_ROLE" == "micro" ]]; then
    echo "DB_micro=\"$DB_DIR\""
    echo
    echo "${BOLD}Or with the new MTD_SE option:${RESET}"
    echo "  --kraken-micro-db \"$DB_DIR\""
else
    echo "DB_host=\"$DB_DIR\""
fi

if [[ "$BUILD_BRACKEN" -eq 1 ]]; then
    echo
    echo "${BOLD}Bracken distribution file:${RESET}"
    echo "$DB_DIR/database${BRACKEN_READ_LEN}mers.kmer_distrib"
fi
echo "${GREEN}${BOLD}============================================================${RESET}"
