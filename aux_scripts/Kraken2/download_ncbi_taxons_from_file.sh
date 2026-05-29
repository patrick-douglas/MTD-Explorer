#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# download_ncbi_taxons_from_file.sh
#
# Wrapper para rodar download_ncbi_taxon_genomes_manifest.sh
# usando vários TaxIDs em um arquivo taxons.txt.
#
# Ele gera:
#   PREFIX_ALL_genome_info.tsv
#
# Esse arquivo final pode ser usado no build_multi_genome_kraken2_db.sh
# ============================================================

usage() {
cat <<EOF

Usage:
  $0 --taxons-file taxons.txt --outdir DIR [options]

Required:
  --taxons-file FILE             Arquivo com TaxIDs, um por linha.
                                 Comentários com # são permitidos.

  --outdir DIR                   Diretório de saída.

Optional:
  --script FILE                  Script original que baixa um TaxID.
                                 Default: ./download_ncbi_taxon_genomes_manifest.sh

  --prefix NAME                  Prefixo geral dos outputs.
                                 Default: multi_taxon

  --assembly-level LIST          Default: complete,chromosome,scaffold,contig

  --assembly-source SOURCE       all, RefSeq, or GenBank.
                                 Default: all

  --include LIST                 Default: genome

  --api-key KEY                  Optional NCBI API key.

  --threads INT                  Default: auto-detect.

  --keep-zip                     Mantém os zips baixados.

  --no-exclude-atypical          Não exclui assemblies atypical.

  --exclude-multi-isolate        Exclui assemblies multi-isolate.

  --preview-only                 Só faz preview, não baixa.

  --force                        Remove outputs anteriores para cada TaxID.

Example:
  bash $0 \\
    --taxons-file taxons.txt \\
    --outdir /media/me/18TB_BACKUP_LBN/drive.ifpa/LBN_RNA-Seq/Metatranscriptomics/MTD/MTD_Offline_Install_files/Ref_genomes/Bglabrata_parasites \\
    --script ./download_ncbi_taxon_genomes_manifest.sh \\
    --prefix bglabrata_parasites \\
    --assembly-level complete,chromosome,scaffold,contig \\
    --assembly-source all \\
    --threads 20 \\
    --force

EOF
}

# -----------------------------
# Defaults
# -----------------------------
TAXONS_FILE=""
OUTDIR=""
SCRIPT="./download_ncbi_taxon_genomes_manifest.sh"
PREFIX="multi_taxon"
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
        --taxons-file) TAXONS_FILE="$2"; shift 2 ;;
        --outdir) OUTDIR="$2"; shift 2 ;;
        --script) SCRIPT="$2"; shift 2 ;;
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
        *) echo "[ERROR] Unknown option: $1" >&2; usage; exit 1 ;;
    esac
done

# -----------------------------
# Checks
# -----------------------------
if [[ -z "$TAXONS_FILE" ]]; then echo "[ERROR] --taxons-file is required." >&2; usage; exit 1; fi
if [[ -z "$OUTDIR" ]]; then echo "[ERROR] --outdir is required." >&2; usage; exit 1; fi
if [[ ! -s "$TAXONS_FILE" ]]; then echo "[ERROR] Taxons file not found or empty: $TAXONS_FILE" >&2; exit 1; fi
if [[ ! -x "$SCRIPT" ]]; then
    echo "[ERROR] Original script not found or not executable: $SCRIPT" >&2
    echo "Run: chmod +x $SCRIPT" >&2
    exit 1
fi
if ! [[ "$THREADS" =~ ^[0-9]+$ ]] || [[ "$THREADS" -lt 1 ]]; then
    echo "[ERROR] --threads must be >= 1. Got: $THREADS" >&2
    exit 1
fi

mkdir -p "$OUTDIR"
OUTDIR_ABS=$(cd "$OUTDIR" && pwd)

CLEAN_TAXIDS="$OUTDIR_ABS/${PREFIX}_taxids_clean.txt"
FINAL_GENOME_INFO="$OUTDIR_ABS/${PREFIX}_ALL_genome_info.tsv"
FINAL_METADATA="$OUTDIR_ABS/${PREFIX}_ALL_assembly_metadata.tsv"
FINAL_MISSING="$OUTDIR_ABS/${PREFIX}_ALL_missing_fastas.tsv"
RUN_LOG="$OUTDIR_ABS/${PREFIX}_run_summary.tsv"

# -----------------------------
# Prepare clean TaxID list
# -----------------------------
: > "$CLEAN_TAXIDS"

while IFS= read -r line || [[ -n "$line" ]]; do
    # remove comentários
    clean="${line%%#*}"

    # remove espaços extras
    taxid=$(echo "$clean" | awk '{print $1}')

    [[ -z "${taxid:-}" ]] && continue

    if [[ "$taxid" =~ ^[0-9]+$ ]]; then
        echo "$taxid" >> "$CLEAN_TAXIDS"
    else
        echo "[WARNING] Ignoring non-numeric line: $line"
    fi
done < "$TAXONS_FILE"

sort -u "$CLEAN_TAXIDS" -o "$CLEAN_TAXIDS"

if [[ ! -s "$CLEAN_TAXIDS" ]]; then
    echo "[ERROR] No valid numeric TaxIDs found in: $TAXONS_FILE" >&2
    exit 1
fi

echo
echo "============================================================"
echo " Multi-TaxID NCBI downloader for Kraken2"
echo "============================================================"
echo "Taxons file:       $TAXONS_FILE"
echo "Clean TaxIDs:      $CLEAN_TAXIDS"
echo "Outdir:            $OUTDIR_ABS"
echo "Original script:   $SCRIPT"
echo "Prefix:            $PREFIX"
echo "Assembly levels:   $ASSEMBLY_LEVEL"
echo "Assembly source:   $ASSEMBLY_SOURCE"
echo "Include:           $INCLUDE"
echo "Threads:           $THREADS"
echo "============================================================"
echo

echo "TaxIDs to process:"
cat "$CLEAN_TAXIDS"
echo

if [[ "$PREVIEW_ONLY" -eq 0 ]]; then
    : > "$FINAL_GENOME_INFO"
    : > "$FINAL_METADATA"
    echo -e "taxid\tstatus\tgenome_info_file\tmetadata_file\tmissing_file" > "$RUN_LOG"
fi

# -----------------------------
# Run original script for each TaxID
# -----------------------------
while IFS= read -r taxid; do
    [[ -z "$taxid" ]] && continue

    taxon_prefix="${PREFIX}_${taxid}"

    echo
    echo "============================================================"
    echo "[RUNNING] TaxID: $taxid"
    echo "Prefix: $taxon_prefix"
    echo "============================================================"

    cmd=(
        "$SCRIPT"
        --taxid "$taxid"
        --outdir "$OUTDIR_ABS"
        --prefix "$taxon_prefix"
        --assembly-level "$ASSEMBLY_LEVEL"
        --assembly-source "$ASSEMBLY_SOURCE"
        --include "$INCLUDE"
        --threads "$THREADS"
    )

    if [[ -n "$API_KEY" ]]; then cmd+=(--api-key "$API_KEY"); fi
    if [[ "$KEEP_ZIP" -eq 1 ]]; then cmd+=(--keep-zip); fi
    if [[ "$EXCLUDE_ATYPICAL" -eq 0 ]]; then cmd+=(--no-exclude-atypical); fi
    if [[ "$EXCLUDE_MULTI_ISOLATE" -eq 1 ]]; then cmd+=(--exclude-multi-isolate); fi
    if [[ "$PREVIEW_ONLY" -eq 1 ]]; then cmd+=(--preview-only); fi
    if [[ "$FORCE" -eq 1 ]]; then cmd+=(--force); fi

    set +e
    "${cmd[@]}"
    status=$?
    set -e

    genome_info="$OUTDIR_ABS/${taxon_prefix}_genome_info.tsv"
    metadata="$OUTDIR_ABS/${taxon_prefix}_assembly_metadata.tsv"
    missing="$OUTDIR_ABS/${taxon_prefix}_missing_fastas.tsv"

    if [[ "$status" -ne 0 ]]; then
        echo "[WARNING] TaxID $taxid failed with status $status"
        if [[ "$PREVIEW_ONLY" -eq 0 ]]; then
            echo -e "${taxid}\tFAILED\t${genome_info}\t${metadata}\t${missing}" >> "$RUN_LOG"
        fi
        continue
    fi

    if [[ "$PREVIEW_ONLY" -eq 1 ]]; then
        continue
    fi

    if [[ -s "$genome_info" ]]; then
        cat "$genome_info" >> "$FINAL_GENOME_INFO"
        echo -e "${taxid}\tOK\t${genome_info}\t${metadata}\t${missing}" >> "$RUN_LOG"
    else
        echo "[WARNING] No genome_info found for TaxID $taxid: $genome_info"
        echo -e "${taxid}\tNO_GENOME_INFO\t${genome_info}\t${metadata}\t${missing}" >> "$RUN_LOG"
    fi

    if [[ -s "$metadata" ]]; then
        if [[ ! -s "$FINAL_METADATA" ]]; then
            cat "$metadata" >> "$FINAL_METADATA"
        else
            tail -n +2 "$metadata" >> "$FINAL_METADATA"
        fi
    fi

    if [[ -s "$missing" ]]; then
        if [[ ! -s "$FINAL_MISSING" ]]; then
            cat "$missing" >> "$FINAL_MISSING"
        else
            tail -n +2 "$missing" >> "$FINAL_MISSING"
        fi
    fi

done < "$CLEAN_TAXIDS"

if [[ "$PREVIEW_ONLY" -eq 1 ]]; then
    echo
    echo "============================================================"
    echo "[DONE] Preview-only finished."
    echo "============================================================"
    exit 0
fi

# -----------------------------
# Final cleanup/deduplication
# -----------------------------
if [[ ! -s "$FINAL_GENOME_INFO" ]]; then
    echo "[ERROR] Final genome_info is empty: $FINAL_GENOME_INFO" >&2
    echo "Check run log: $RUN_LOG" >&2
    exit 1
fi

sort -u "$FINAL_GENOME_INFO" -o "$FINAL_GENOME_INFO"

n_fastas=$(wc -l < "$FINAL_GENOME_INFO" | awk '{print $1}')
n_taxids=$(cut -f1 "$FINAL_GENOME_INFO" | sort -u | wc -l | awk '{print $1}')

echo
echo "============================================================"
echo "[DONE] Multi-TaxID download complete"
echo "============================================================"
echo "Final Kraken2 genome_info:"
echo "$FINAL_GENOME_INFO"
echo
echo "Final metadata:"
echo "$FINAL_METADATA"
echo
echo "Final missing FASTAs:"
echo "$FINAL_MISSING"
echo
echo "Run summary:"
echo "$RUN_LOG"
echo
echo "Genome FASTAs in final manifest: $n_fastas"
echo "Unique organism TaxIDs:          $n_taxids"
echo
echo "First lines of final genome_info:"
head "$FINAL_GENOME_INFO"
echo "============================================================"
