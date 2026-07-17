#!/usr/bin/env bash
set -eo pipefail

# ============================================================
# download_kraken2_taxonomy_https.sh
#
# Substitui:
#   kraken2-build --download-taxonomy --db DBNAME
#
# Usando HTTPS em vez de rsync.
#
# Uso:
#   ./aux_scripts/Kraken2/download_kraken2_taxonomy_https.sh --db kraken2DB_59463 --threads 20
#
# Opções:
#   --db <folder>      Pasta do banco Kraken2
#   --threads <N>      Número de threads/conexões desejadas
#   --skip-maps        Não baixa accession2taxid
#   --protein          Baixa mapa protein: prot.accession2taxid
# ============================================================

DBNAME=""
THREADS="${threads:-4}"
SKIP_MAPS="false"
PROTEIN_DB="false"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --db)
            DBNAME="$2"
            shift 2
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        --skip-maps)
            SKIP_MAPS="true"
            shift
            ;;
        --protein)
            PROTEIN_DB="true"
            shift
            ;;
        -h|--help)
            echo "Usage:"
            echo "  $0 --db <kraken2_db_folder> [--threads N] [--skip-maps] [--protein]"
            exit 0
            ;;
        *)
            echo "Erro: argumento desconhecido: $1"
            exit 1
            ;;
    esac
done

if [[ -z "$DBNAME" ]]; then
    echo "Erro: use --db <kraken2_db_folder>"
    exit 1
fi

# ------------------------------------------------------------
# Validar threads
# ------------------------------------------------------------
if ! [[ "$THREADS" =~ ^[0-9]+$ ]]; then
    echo "Aviso: --threads '$THREADS' não é um número válido. Usando 4."
    THREADS=4
fi

if [[ "$THREADS" -lt 1 ]]; then
    THREADS=1
fi

# aria2c aceita no máximo 16 conexões por servidor
ARIA2_CONNECTIONS="$THREADS"

if [[ "$ARIA2_CONNECTIONS" -gt 16 ]]; then
    ARIA2_CONNECTIONS=16
fi

if [[ "$ARIA2_CONNECTIONS" -lt 1 ]]; then
    ARIA2_CONNECTIONS=1
fi

TAXONOMY_DIR="$DBNAME/taxonomy"
BASE_URL="https://ftp.ncbi.nlm.nih.gov/pub/taxonomy"

mkdir -p "$TAXONOMY_DIR"

echo "------------------------------------------------------------"
echo "Kraken2 taxonomy HTTPS downloader"
echo "DB: $DBNAME"
echo "Taxonomy dir: $TAXONOMY_DIR"
echo "Requested threads/connections: $THREADS"
echo "aria2 connections used: $ARIA2_CONNECTIONS"
echo "Skip accession maps: $SKIP_MAPS"
echo "Protein DB mode: $PROTEIN_DB"
echo "------------------------------------------------------------"

download_file() {
    local url="$1"
    local output="$2"

    if [[ -s "$output" ]]; then
        echo "[OK] Already exists: $output"
        return 0
    fi

    echo "[DOWNLOAD] $output"
    echo "URL: $url"

    if command -v aria2c >/dev/null 2>&1; then
        aria2c \
            --continue=true \
            --max-connection-per-server="$ARIA2_CONNECTIONS" \
            --split="$ARIA2_CONNECTIONS" \
            --min-split-size=1M \
            --retry-wait=5 \
            --max-tries=10 \
            --timeout=60 \
            --dir "$(dirname "$output")" \
            --out "$(basename "$output")" \
            "$url"

    elif command -v wget >/dev/null 2>&1; then
        wget \
            -c \
            --tries=10 \
            --timeout=60 \
            -O "$output" \
            "$url"

    elif command -v curl >/dev/null 2>&1; then
        curl \
            -L \
            --retry 10 \
            --retry-delay 5 \
            --connect-timeout 60 \
            -o "$output" \
            "$url"

    else
        echo "Erro: instale aria2c, wget ou curl."
        exit 1
    fi

    if [[ ! -s "$output" ]]; then
        echo "Erro: download falhou ou arquivo vazio: $output"
        exit 1
    fi
}

cd "$TAXONOMY_DIR" || exit 1

# ------------------------------------------------------------
# 1. Accession to TaxID maps
# ------------------------------------------------------------
if [[ "$SKIP_MAPS" == "false" && ! -e "accmap.dlflag" ]]; then

    if [[ "$PROTEIN_DB" == "true" ]]; then
        download_file \
            "$BASE_URL/accession2taxid/prot.accession2taxid.gz" \
            "prot.accession2taxid.gz"
    else
        download_file \
            "$BASE_URL/accession2taxid/nucl_gb.accession2taxid.gz" \
            "nucl_gb.accession2taxid.gz"

        download_file \
            "$BASE_URL/accession2taxid/nucl_wgs.accession2taxid.gz" \
            "nucl_wgs.accession2taxid.gz"
    fi

    touch accmap.dlflag
    echo "[OK] Accession maps downloaded."
else
    echo "[SKIP] Accession maps already downloaded or skipped."
fi

# ------------------------------------------------------------
# 2. Taxdump
# ------------------------------------------------------------
if [[ ! -e "taxdump.dlflag" ]]; then
    download_file \
        "$BASE_URL/taxdump.tar.gz" \
        "taxdump.tar.gz"

    touch taxdump.dlflag
    echo "[OK] taxdump.tar.gz downloaded."
else
    echo "[SKIP] taxdump.tar.gz already downloaded."
fi

# ------------------------------------------------------------
# 3. Uncompress accession maps
# ------------------------------------------------------------
shopt -s nullglob
gz_maps=( *.accession2taxid.gz )

if [[ ${#gz_maps[@]} -gt 0 ]]; then
    echo "[PROCESS] Uncompressing accession2taxid files..."

    for gz in "${gz_maps[@]}"; do
        out="${gz%.gz}"

        if [[ -s "$out" ]]; then
            echo "[OK] Already uncompressed: $out"
            rm -f "$gz"
        else
            echo "[PROCESS] gunzip $gz"
            gunzip -f "$gz"
            echo "[OK] Created: $out"
        fi
    done
else
    echo "[SKIP] No accession2taxid.gz files to uncompress."
fi

# ------------------------------------------------------------
# 4. Extract taxdump
# ------------------------------------------------------------
if [[ ! -e "taxdump.untarflag" ]]; then
    echo "[PROCESS] Extracting taxdump.tar.gz..."

    if [[ ! -s "taxdump.tar.gz" ]]; then
        echo "Erro: taxdump.tar.gz não existe ou está vazio."
        exit 1
    fi

    tar -xzf taxdump.tar.gz
    touch taxdump.untarflag
    echo "[OK] taxdump extracted."
else
    echo "[SKIP] taxdump already extracted."
fi

# ------------------------------------------------------------
# 5. Basic validation
# ------------------------------------------------------------
echo "------------------------------------------------------------"
echo "Validation"
echo "------------------------------------------------------------"

required_files=(
    "names.dmp"
    "nodes.dmp"
)

for f in "${required_files[@]}"; do
    if [[ -s "$f" ]]; then
        echo "[OK] $f"
    else
        echo "[ERROR] Missing or empty: $TAXONOMY_DIR/$f"
        exit 1
    fi
done

if [[ "$SKIP_MAPS" == "false" ]]; then
    if [[ "$PROTEIN_DB" == "true" ]]; then
        map_files=( "prot.accession2taxid" )
    else
        map_files=(
            "nucl_gb.accession2taxid"
            "nucl_wgs.accession2taxid"
        )
    fi

    for f in "${map_files[@]}"; do
        if [[ -s "$f" ]]; then
            echo "[OK] $f"
        else
            echo "[ERROR] Missing or empty: $TAXONOMY_DIR/$f"
            exit 1
        fi
    done
fi

# ------------------------------------------------------------
# 6. Optional: show species name if DBNAME ends in TaxID
# Example:
#   DBNAME=kraken2DB_59463
# ------------------------------------------------------------
taxid="${DBNAME##*_}"

if [[ "$taxid" =~ ^[0-9]+$ ]]; then
    scientific_name=$(awk -F'\t\\|\t' -v taxid="$taxid" '
        $1 == taxid && $4 ~ /^scientific name/ {
            gsub(/^[ \t]+|[ \t]+$/, "", $2)
            print $2
            exit
        }
    ' names.dmp)

    echo "------------------------------------------------------------"
    echo "DB information"
    echo "------------------------------------------------------------"
    echo "DBNAME: $DBNAME"
    echo "TaxID from DBNAME: $taxid"

    if [[ -n "$scientific_name" ]]; then
        echo "Scientific name: $scientific_name"
    else
        echo "Scientific name: not found in names.dmp"
    fi
else
    echo "------------------------------------------------------------"
    echo "DB information"
    echo "------------------------------------------------------------"
    echo "DBNAME: $DBNAME"
    echo "TaxID from DBNAME: not detected"
fi

echo "------------------------------------------------------------"
echo "Done."
echo "Taxonomy prepared at:"
echo "$TAXONOMY_DIR"
echo "------------------------------------------------------------"
