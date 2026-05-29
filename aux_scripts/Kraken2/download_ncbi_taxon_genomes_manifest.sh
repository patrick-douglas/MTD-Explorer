#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# download_ncbi_taxon_genomes_manifest.sh
#
# Download all available NCBI genome assemblies under a taxon
# using NCBI Datasets CLI and create a genome_info TSV ready for
# build_multi_genome_kraken2_db.sh:
#
#   taxid<TAB>/path/to/genomic.fna
#
# Example:
#   bash download_ncbi_taxon_genomes_manifest.sh \
#     --taxid 6178 \
#     --outdir /media/me/18TB_BACKUP_LBN/drive.ifpa/LBN_RNA-Seq/Metatranscriptomics/MTD/MTD_Offline_Install_files/Ref_genomes/Trematoda \
#     --prefix trematoda \
#     --assembly-level complete,chromosome,scaffold,contig \
#     --assembly-source all \
#     --threads 20 \
#     --force
#
# Requirements:
#   datasets, dataformat, unzip, awk, sed, find
# ============================================================

# -----------------------------
# Pretty output
# -----------------------------
if [[ -t 1 ]]; then
    BOLD="$(tput bold || true)"
    RESET="$(tput sgr0 || true)"
    GREEN="$(tput setaf 2 || true)"
    RED="$(tput setaf 1 || true)"
    YELLOW="$(tput setaf 3 || true)"
    BLUE="$(tput setaf 4 || true)"
    CYAN="$(tput setaf 6 || true)"
else
    BOLD=""; RESET=""; GREEN=""; RED=""; YELLOW=""; BLUE=""; CYAN=""
fi

info() { echo "${CYAN}${BOLD}[INFO]${RESET} $*"; }
ok()   { echo "${GREEN}${BOLD}[OK]${RESET} $*"; }
warn() { echo "${YELLOW}${BOLD}[WARNING]${RESET} $*"; }
err()  { echo "${RED}${BOLD}[ERROR]${RESET} $*" >&2; }
step() {
    echo
    echo "${BLUE}${BOLD}============================================================${RESET}"
    echo "${BLUE}${BOLD}$*${RESET}"
    echo "${BLUE}${BOLD}============================================================${RESET}"
}

usage() {
cat <<EOF_USAGE

${BOLD}Usage:${RESET}
  $0 --taxid TAXID --outdir DIR [options]

${BOLD}Required:${RESET}
  --taxid TAXID                  Root NCBI TaxID to download.
                                 Example: Trematoda = 6178

  --outdir DIR                   Output directory.

${BOLD}Optional:${RESET}
  --prefix NAME                  Prefix for output files.
                                 Default: taxon_TAXID

  --assembly-level LIST          Assembly levels to include.
                                 Default: complete,chromosome,scaffold,contig
                                 NCBI values: complete,chromosome,scaffold,contig

  --assembly-source SOURCE       all, RefSeq, or GenBank.
                                 Default: all

  --include LIST                 Files to download from NCBI Datasets.
                                 Default: genome
                                 For Kraken2, genome is enough.
                                 Examples: genome,gff3,gtf or genome,cds,protein

  --api-key KEY                  Optional NCBI API key.

  --threads INT                  Threads used by datasets rehydrate.
                                 Default: auto-detect with nproc minus 2

  --keep-zip                     Keep the downloaded dehydrated zip.
                                 Default: remove zip after successful rehydrate

  --no-exclude-atypical          Do not pass --exclude-atypical.
                                 Default: exclude atypical assemblies

  --exclude-multi-isolate        Pass --exclude-multi-isolate.
                                 Default: not used

  --preview-only                 Only preview/list matching genomes; do not download.

  --force                        Remove previous output package dir/files.

  -h, --help                     Show this help.

${BOLD}Outputs:${RESET}
  PREFIX_assembly_metadata.tsv   Assembly metadata table
  PREFIX_accessions.txt          Assembly accession list
  PREFIX_genome_info.tsv         taxid<TAB>genomic_fasta path for Kraken2 builder
  PREFIX_missing_fastas.tsv      Assemblies with metadata but no genomic FASTA found

${BOLD}Typical Trematoda example:${RESET}
  bash $0 \
    --taxid 6178 \
    --outdir /media/me/18TB_BACKUP_LBN/drive.ifpa/LBN_RNA-Seq/Metatranscriptomics/MTD/MTD_Offline_Install_files/Ref_genomes/Trematoda \
    --prefix trematoda \
    --assembly-level complete,chromosome,scaffold,contig \
    --assembly-source all \
    --threads 20 \
    --force

EOF_USAGE
}

# -----------------------------
# Defaults
# -----------------------------
TAXID=""
OUTDIR=""
PREFIX=""
ASSEMBLY_LEVEL="complete,chromosome,scaffold,contig"
ASSEMBLY_SOURCE="all"
INCLUDE="genome"
API_KEY=""
KEEP_ZIP=0
EXCLUDE_ATYPICAL=1
EXCLUDE_MULTI_ISOLATE=0
PREVIEW_ONLY=0
FORCE=0

if command -v nproc >/dev/null 2>&1; then
    NPROC=$(nproc)
    if [[ "$NPROC" -gt 2 ]]; then THREADS=$((NPROC - 2)); else THREADS=1; fi
else
    THREADS=1
fi

# -----------------------------
# Parse args
# -----------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --taxid) TAXID="$2"; shift 2 ;;
        --outdir) OUTDIR="$2"; shift 2 ;;
        --prefix) PREFIX="$2"; shift 2 ;;
        --assembly-level) ASSEMBLY_LEVEL="$2"; shift 2 ;;
        --assembly-source) ASSEMBLY_SOURCE="$2"; shift 2 ;;
        --include) INCLUDE="$2"; shift 2 ;;
        --api-key) API_KEY="$2"; shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        --keep-zip) KEEP_ZIP=1; shift ;;
        --no-exclude-atypical) EXCLUDE_ATYPICAL=0; shift ;;
        --exclude-multi-isolate) EXCLUDE_MULTI_ISOLATE=1; shift ;;
        --preview-only) PREVIEW_ONLY=1; shift ;;
        --force) FORCE=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) err "Unknown option: $1"; usage; exit 1 ;;
    esac
done

# -----------------------------
# Checks
# -----------------------------
if [[ -z "$TAXID" ]]; then err "--taxid is required."; usage; exit 1; fi
if [[ -z "$OUTDIR" ]]; then err "--outdir is required."; usage; exit 1; fi
if [[ -z "$PREFIX" ]]; then PREFIX="taxon_${TAXID}"; fi
if ! [[ "$TAXID" =~ ^[0-9]+$ ]]; then err "--taxid must be numeric. Got: $TAXID"; exit 1; fi
if ! [[ "$THREADS" =~ ^[0-9]+$ ]] || [[ "$THREADS" -lt 1 ]]; then err "--threads must be >= 1. Got: $THREADS"; exit 1; fi

for cmd in datasets dataformat unzip awk sed find sort; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        err "Required command not found in PATH: $cmd"
        echo "Install/activate NCBI Datasets CLI first if missing: datasets and dataformat."
        exit 1
    fi
done

mkdir -p "$OUTDIR"
OUTDIR_ABS=$(cd "$OUTDIR" && pwd)

ZIP_FILE="$OUTDIR_ABS/${PREFIX}_ncbi_dataset_dehydrated.zip"
PACKAGE_DIR="$OUTDIR_ABS/${PREFIX}_ncbi_dataset"
LOG_DIR="$OUTDIR_ABS/${PREFIX}_logs"
METADATA_JSONL="$PACKAGE_DIR/ncbi_dataset/data/assembly_data_report.jsonl"
METADATA_TSV="$OUTDIR_ABS/${PREFIX}_assembly_metadata.tsv"
ACCESSIONS="$OUTDIR_ABS/${PREFIX}_accessions.txt"
GENOME_INFO="$OUTDIR_ABS/${PREFIX}_genome_info.tsv"
MISSING_FASTAS="$OUTDIR_ABS/${PREFIX}_missing_fastas.tsv"
PREVIEW_LOG="$LOG_DIR/${PREFIX}_preview.txt"

if [[ "$FORCE" -eq 1 ]]; then
    warn "--force used. Removing previous outputs for prefix: $PREFIX"
    rm -rf "$PACKAGE_DIR" "$LOG_DIR"
    rm -f "$ZIP_FILE" "$METADATA_TSV" "$ACCESSIONS" "$GENOME_INFO" "$MISSING_FASTAS"
fi

mkdir -p "$LOG_DIR"

# -----------------------------
# Build datasets args
# -----------------------------
DATASETS_ARGS=(genome taxon "$TAXID")
DATASETS_ARGS+=(--assembly-level "$ASSEMBLY_LEVEL")
DATASETS_ARGS+=(--assembly-source "$ASSEMBLY_SOURCE")
DATASETS_ARGS+=(--include "$INCLUDE")
DATASETS_ARGS+=(--no-progressbar)

if [[ -n "$API_KEY" ]]; then DATASETS_ARGS+=(--api-key "$API_KEY"); fi
if [[ "$EXCLUDE_ATYPICAL" -eq 1 ]]; then DATASETS_ARGS+=(--exclude-atypical); fi
if [[ "$EXCLUDE_MULTI_ISOLATE" -eq 1 ]]; then DATASETS_ARGS+=(--exclude-multi-isolate); fi

# -----------------------------
# Header
# -----------------------------
echo
echo "${GREEN}${BOLD}============================================================${RESET}"
echo "${GREEN}${BOLD} NCBI taxon genome downloader + Kraken2 manifest maker${RESET}"
echo "${GREEN}${BOLD}============================================================${RESET}"
echo "TaxID:              $TAXID"
echo "Outdir:             $OUTDIR_ABS"
echo "Prefix:             $PREFIX"
echo "Assembly levels:    $ASSEMBLY_LEVEL"
echo "Assembly source:    $ASSEMBLY_SOURCE"
echo "Include:            $INCLUDE"
echo "Exclude atypical:   $EXCLUDE_ATYPICAL"
echo "Exclude multi-iso:  $EXCLUDE_MULTI_ISOLATE"
echo "Threads:            $THREADS"
echo "Package dir:        $PACKAGE_DIR"
echo "${GREEN}${BOLD}============================================================${RESET}"
echo

# -----------------------------
# Preview
# -----------------------------
step "[STEP 1] Preview matching NCBI genome package"

set +e
datasets download "${DATASETS_ARGS[@]}" --preview > "$PREVIEW_LOG" 2>&1
preview_status=$?
set -e

if [[ "$preview_status" -ne 0 ]]; then
    warn "datasets preview exited with status $preview_status. Continuing, but check: $PREVIEW_LOG"
else
    ok "Preview saved to: $PREVIEW_LOG"
    echo
    sed -n '1,80p' "$PREVIEW_LOG" || true
fi

if [[ "$PREVIEW_ONLY" -eq 1 ]]; then
    warn "--preview-only used. Stopping before download."
    exit 0
fi

# -----------------------------
# Download dehydrated package
# -----------------------------
step "[STEP 2] Download dehydrated NCBI Datasets package"

if [[ -e "$ZIP_FILE" ]]; then
    err "Zip already exists: $ZIP_FILE"
    echo "Use --force to overwrite or choose another --prefix."
    exit 1
fi

# Use dehydrated download first because it is safer for large taxon-wide packages.
datasets download "${DATASETS_ARGS[@]}" \
    --dehydrated \
    --filename "$ZIP_FILE" \
    > "$LOG_DIR/${PREFIX}_download.log" 2>&1

ok "Dehydrated zip downloaded: $ZIP_FILE"

# -----------------------------
# Unzip
# -----------------------------
step "[STEP 3] Unzip package"

mkdir -p "$PACKAGE_DIR"
unzip -q "$ZIP_FILE" -d "$PACKAGE_DIR"
ok "Package extracted to: $PACKAGE_DIR"

# -----------------------------
# Rehydrate
# -----------------------------
step "[STEP 4] Rehydrate package: download genome FASTA files"

# Some datasets versions support --threads for rehydrate; if not, retry without it.
set +e
datasets rehydrate --directory "$PACKAGE_DIR" --threads "$THREADS" \
    > "$LOG_DIR/${PREFIX}_rehydrate.log" 2>&1
rehydrate_status=$?
set -e

if [[ "$rehydrate_status" -ne 0 ]]; then
    warn "rehydrate with --threads failed. Retrying without --threads."
    datasets rehydrate --directory "$PACKAGE_DIR" \
        > "$LOG_DIR/${PREFIX}_rehydrate_retry_no_threads.log" 2>&1
fi

ok "Rehydration complete."

if [[ ! -s "$METADATA_JSONL" ]]; then
    err "Metadata JSONL not found: $METADATA_JSONL"
    echo "Check logs in: $LOG_DIR"
    exit 1
fi

# -----------------------------
# Metadata table
# -----------------------------
step "[STEP 5] Create metadata table"

dataformat tsv genome \
    --inputfile "$METADATA_JSONL" \
    --fields accession,organism-tax-id,organism-name,assminfo-level,source_database,assminfo-name \
    > "$METADATA_TSV"

if [[ ! -s "$METADATA_TSV" ]]; then
    err "Metadata TSV was not created or is empty: $METADATA_TSV"
    exit 1
fi

awk 'NR > 1 {print $1}' "$METADATA_TSV" | sort -u > "$ACCESSIONS"

ok "Metadata TSV: $METADATA_TSV"
ok "Accessions:   $ACCESSIONS"

# -----------------------------
# Create genome_info.tsv for Kraken2 builder
# -----------------------------
step "[STEP 6] Create Kraken2 genome_info manifest"

: > "$GENOME_INFO"
echo -e "accession\ttaxid\torganism_name\tassembly_level\tsource_database\tassembly_name\tstatus" > "$MISSING_FASTAS"

DATA_DIR="$PACKAGE_DIR/ncbi_dataset/data"
found=0
missing=0

# Expected columns from dataformat:
# accession, organism-tax-id, organism-name, assminfo-level, source_database, assminfo-name
# Use tail -n +2 to skip header.
while IFS=$'\t' read -r accession taxid organism_name assembly_level source_database assembly_name rest; do
    [[ -z "${accession:-}" ]] && continue
    [[ -z "${taxid:-}" ]] && taxid="NA"

    asm_dir="$DATA_DIR/$accession"

    fasta=""
    if [[ -d "$asm_dir" ]]; then
        fasta=$(find "$asm_dir" -maxdepth 1 -type f \( -name "*_genomic.fna" -o -name "*_genomic.fna.gz" -o -name "*.fna" -o -name "*.fna.gz" \) | sort | head -n 1 || true)
    fi

    if [[ -n "$fasta" && -s "$fasta" && "$taxid" =~ ^[0-9]+$ ]]; then
        echo -e "${taxid}\t${fasta}" >> "$GENOME_INFO"
        found=$((found + 1))
    else
        echo -e "${accession}\t${taxid}\t${organism_name}\t${assembly_level}\t${source_database}\t${assembly_name}\tMISSING_FASTA_OR_TAXID" >> "$MISSING_FASTAS"
        missing=$((missing + 1))
    fi

done < <(tail -n +2 "$METADATA_TSV")

if [[ ! -s "$GENOME_INFO" ]]; then
    err "No genomic FASTA files were found. genome_info is empty: $GENOME_INFO"
    echo "Check package dir: $PACKAGE_DIR"
    echo "Check missing table: $MISSING_FASTAS"
    exit 1
fi

# Remove accidental duplicates while preserving useful sorted output.
sort -u "$GENOME_INFO" -o "$GENOME_INFO"

ok "Kraken2 genome_info manifest: $GENOME_INFO"

# -----------------------------
# Summary
# -----------------------------
step "[STEP 7] Summary"

n_accessions=$(wc -l < "$ACCESSIONS" | awk '{print $1}')
n_manifest=$(wc -l < "$GENOME_INFO" | awk '{print $1}')
n_taxids=$(cut -f1 "$GENOME_INFO" | sort -u | wc -l | awk '{print $1}')

# Count metadata species names, if available.
n_orgs=$(tail -n +2 "$METADATA_TSV" | cut -f3 | sort -u | wc -l | awk '{print $1}')

echo "Assembly accessions listed:     $n_accessions"
echo "Genome FASTAs in manifest:      $n_manifest"
echo "Unique organism taxids:         $n_taxids"
echo "Unique organism names in table: $n_orgs"
echo "Missing FASTA/taxid records:    $missing"
echo

echo "First lines of genome_info:"
head "$GENOME_INFO" || true

if [[ "$missing" -gt 0 ]]; then
    warn "Some metadata records had no FASTA/taxid. See: $MISSING_FASTAS"
fi

if [[ "$KEEP_ZIP" -eq 0 ]]; then
    rm -f "$ZIP_FILE"
    ok "Removed dehydrated zip. Use --keep-zip if you want to keep it next time."
fi

echo
echo "${GREEN}${BOLD}============================================================${RESET}"
echo "${GREEN}${BOLD}[DONE] Taxon genome download + manifest complete${RESET}"
echo "${GREEN}${BOLD}============================================================${RESET}"
echo "Genome info for Kraken2 builder:"
echo "$GENOME_INFO"
echo
echo "Next example:"
echo "bash build_multi_genome_kraken2_db.sh \\
  --genome_info \"$GENOME_INFO\" \\
  --outdir /path/to/MTD_Offline_Install_files \\
  --db_name Kraken2DB_trematoda \\
  --kraken_offline /path/to/existing/Kraken2DB_micro \\
  --threads 20 \\
  --force"
echo "${GREEN}${BOLD}============================================================${RESET}"

