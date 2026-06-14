#!/bin/bash
# v5: gold OrgDb workflow + safe cleanup + shared Kraken taxonomy cache + fixed HostSpecies awk

MTDIR=~/MTD
condapath=~/miniconda3
threads=$(nproc)

protein_fasta=""
offline_files_folder=""
SKIP_ORGDB=0
FORCE_ORGDB=0
SKIP_IF_COMPATIBLE=1
ORGDB_VERSION="0.3.0"
CLEAN_PREVIOUS=1
CLEAN_ONLY=0
USE_KRAKEN_TAXONOMY_CACHE=1
REBUILD_KRAKEN_TAXONOMY_CACHE=0
KRAKEN_TAXONOMY_CACHE=""


# ------------------------------------------------------------
# Add a path to a colon-separated environment variable only once
# ------------------------------------------------------------
prepend_unique_path_var() {
    local var_name="$1"
    local new_path="$2"
    local current="${!var_name:-}"
    local cleaned=""
    local item

    # Remove duplicates and empty fields from the current variable.
    IFS=':' read -r -a _path_items <<< "$current"
    for item in "${_path_items[@]}"; do
        [[ -z "$item" ]] && continue
        [[ "$item" == "$new_path" ]] && continue
        case ":$cleaned:" in
            *":$item:"*) ;;
            *)
                if [[ -z "$cleaned" ]]; then cleaned="$item"; else cleaned="$cleaned:$item"; fi
                ;;
        esac
    done

    if [[ -z "$cleaned" ]]; then
        export "$var_name=$new_path"
    else
        export "$var_name=$new_path:$cleaned"
    fi
}

# Função para mostrar instruções de uso
usage() {
    echo ""
    echo "Usage:"
    echo "  bash $0 --genome <genome.fa.gz> --gtf-file <annotations.gtf.gz> --ncbi-taxon-id <TaxonID> --protein-fasta <proteins.fa.gz>"
    echo ""
    echo "Required options:"
    echo "  -d, --genome           Path to host genome FASTA (.fa.gz)"
    echo "  -g, --gtf-file         Path to GTF annotation file (.gtf.gz)"
    echo "  -c, --ncbi-taxon-id    NCBI Taxon ID of the host species"
    echo ""
    echo "Gold OrgDb options:"
    echo "  -p, --protein-fasta    Protein FASTA from the same annotation/release as the GTF. Required to create the gold OrgDb."
    echo "      --skip-orgdb       Skip gold OrgDb creation"
    echo "      --force-orgdb      Rebuild gold OrgDb even if an existing compatible OrgDb is detected"
    echo "      --no-skip-if-compatible"
    echo "                         Do not skip when an existing OrgDb appears compatible"
    echo "      --orgdb-version    Version for the generated OrgDb package [default: 0.3.0]"
    echo "      --no-clean         Do not remove previous files for this TaxID/reference before starting"
    echo "      --clean-only       Remove previous files for this TaxID/reference and exit"
    echo ""
    echo "Kraken taxonomy cache options:"
    echo "      --kraken-taxonomy-cache DIR"
    echo "                         Shared Kraken2 taxonomy cache directory [default: $MTDIR/kraken2_taxonomy_cache]"
    echo "      --no-kraken-taxonomy-cache"
    echo "                         Do not use shared cache; download taxonomy into the species DB as before"
    echo "      --rebuild-kraken-taxonomy-cache"
    echo "                         Delete and rebuild the shared Kraken2 taxonomy cache before using it"
    echo ""
    echo "Deprecated/ignored:"
    echo "  -o, --offline-folder   Kept for backwards compatibility with the old NCBI OrgDb workflow. Not used by the gold OrgDb workflow."
    echo ""
    echo "Other:"
    echo "      --help             Show this help message"
    echo ""
    echo "Example:"
    echo "  bash $0 \\\n    --genome /path/to/Biomphalaria.dna.toplevel.fa.gz \\\n    --gtf-file /path/to/Biomphalaria.62.gtf.gz \\\n    --protein-fasta /path/to/Biomphalaria.pep.all.fa.gz \\\n    --ncbi-taxon-id 6526"
    echo ""
    exit 1
}

# Define as opções curtas e longas
OPTIONS=c:d:g:o:p:
LONGOPTIONS=ncbi-taxon-id:,genome:,gtf-file:,offline-folder:,protein-fasta:,skip-orgdb,force-orgdb,no-skip-if-compatible,orgdb-version:,no-clean,clean-only,kraken-taxonomy-cache:,no-kraken-taxonomy-cache,rebuild-kraken-taxonomy-cache,help

# Processa os argumentos
PARSED=$(getopt --options=$OPTIONS --longoptions=$LONGOPTIONS --name "$0" -- "$@")
if [[ $? -ne 0 ]]; then
    usage
fi

eval set -- "$PARSED"

# Variáveis
while true; do
    case "$1" in
        -c|--ncbi-taxon-id)
            customized="$2"
            shift 2
            ;;
        -d|--genome)
            download="$2"
            shift 2
            ;;
        -g|--gtf-file)
            gtf="$2"
            shift 2
            ;;
        -o|--offline-folder)
            offline_files_folder="$2"
            shift 2
            ;;
        -p|--protein-fasta)
            protein_fasta="$2"
            shift 2
            ;;
        --skip-orgdb)
            SKIP_ORGDB=1
            shift
            ;;
        --force-orgdb)
            FORCE_ORGDB=1
            shift
            ;;
        --no-skip-if-compatible)
            SKIP_IF_COMPATIBLE=0
            shift
            ;;
        --orgdb-version)
            ORGDB_VERSION="$2"
            shift 2
            ;;
        --no-clean)
            CLEAN_PREVIOUS=0
            shift
            ;;
        --clean-only)
            CLEAN_ONLY=1
            shift
            ;;
        --kraken-taxonomy-cache)
            KRAKEN_TAXONOMY_CACHE="$2"
            shift 2
            ;;
        --no-kraken-taxonomy-cache)
            USE_KRAKEN_TAXONOMY_CACHE=0
            shift
            ;;
        --rebuild-kraken-taxonomy-cache)
            REBUILD_KRAKEN_TAXONOMY_CACHE=1
            shift
            ;;
        --help)
            usage
            ;;
        --)
            shift
            break
            ;;
        *)
            echo "Error: invalid arguments: $1"
            usage
            ;;
    esac
done

# Verificação básica de argumentos obrigatórios
if [[ -z "${download:-}" || -z "${gtf:-}" || -z "${customized:-}" ]]; then
    echo "Error: Missing required arguments"
    usage
fi

if [[ ! -s "$download" ]]; then
    echo "Error: genome FASTA not found or empty: $download"
    exit 1
fi

if [[ ! -s "$gtf" ]]; then
    echo "Error: GTF file not found or empty: $gtf"
    exit 1
fi

if [[ "$SKIP_ORGDB" != "1" && -z "${protein_fasta:-}" ]]; then
    echo "WARNING: --protein-fasta was not provided. Gold OrgDb creation will be skipped."
    echo "WARNING: To create a reference-matched OrgDb, provide the protein FASTA from the same annotation/release as the GTF."
    SKIP_ORGDB=1
fi

if [[ "$SKIP_ORGDB" != "1" && ! -s "$protein_fasta" ]]; then
    echo "Error: protein FASTA not found or empty: $protein_fasta"
    exit 1
fi

# get MTD folder place; same as Install.sh script file path (in the MTD folder)
dir=$(dirname $(readlink -f $0))
cd $dir # MTD folder place

# Normalize MTDIR after possible tilde expansion.
MTDIR=$(cd "$MTDIR" && pwd)

if [[ -z "$KRAKEN_TAXONOMY_CACHE" ]]; then
    KRAKEN_TAXONOMY_CACHE="$MTDIR/kraken2_taxonomy_cache"
fi

# ------------------------------------------------------------
# Clean previous files from the same TaxID/reference
# ------------------------------------------------------------
get_orgdb_from_hostspecies() {
    local taxid="$1"
    local host_csv="$2"

    if [[ ! -s "$host_csv" ]]; then
        return 0
    fi

    # HostSpecies.csv is expected to contain at least:
    # Taxon_ID,...,OrgDb,...
    # We remove Windows CR characters before awk to avoid CR-character regex portability issues.
    tr -d '\r' < "$host_csv" | awk -F',' -v id="$taxid" '
        NR == 1 {
            for (i = 1; i <= NF; i++) {
                gsub(/^[[:space:]]+|[[:space:]]+$/, "", $i)
                if ($i == "Taxon_ID") tax_col = i
                if ($i == "OrgDb") org_col = i
            }
            next
        }

        tax_col && org_col {
            tax_id = $tax_col
            orgdb = $org_col

            gsub(/^[[:space:]]+|[[:space:]]+$/, "", tax_id)
            gsub(/^[[:space:]]+|[[:space:]]+$/, "", orgdb)

            if (tax_id == id && orgdb != "" && orgdb != "NA") {
                print orgdb
                exit
            }
        }
    '
}

remove_path_if_exists() {
    local p="$1"
    if [[ -e "$p" || -L "$p" ]]; then
        echo "[CLEAN] Removing: $p"
        rm -rf --one-file-system "$p"
    else
        echo "[CLEAN] Not found, skipping: $p"
    fi
}

cleanup_previous_reference() {
    local taxid="$1"
    local mtdir="$2"
    local host_csv="$mtdir/HostSpecies.csv"
    local orgdb_pkg=""

    echo "------------------------------------------------------------"
    echo "Cleaning previous files for TaxID ${taxid}"
    echo "MTD directory: $mtdir"
    echo "------------------------------------------------------------"

    # Main MTD reference outputs for this TaxID.
    remove_path_if_exists "$mtdir/kraken2DB_${taxid}"
    remove_path_if_exists "$mtdir/ref_${taxid}"
    remove_path_if_exists "$mtdir/hisat2_index_${taxid}"
    remove_path_if_exists "$mtdir/blastdb_${taxid}"

    # Gold OrgDb build/cache for this TaxID.
    remove_path_if_exists "$mtdir/build/orgdb_gold/${taxid}"

    # Older experimental build folders used during development.
    remove_path_if_exists "$mtdir/build/orgdb_${taxid}"
    remove_path_if_exists "$mtdir/build/orgdb_taxid_${taxid}"

    # Remove the custom OrgDb package from the custom R library, but only
    # the package listed for this TaxID in HostSpecies.csv.
    orgdb_pkg=$(get_orgdb_from_hostspecies "$taxid" "$host_csv" || true)

    if [[ -n "$orgdb_pkg" ]]; then
        echo "[CLEAN] OrgDb listed for TaxID ${taxid}: $orgdb_pkg"
        remove_path_if_exists "$mtdir/custom_R_libs/$orgdb_pkg"
        remove_path_if_exists "$mtdir/custom_R_libs/00LOCK-$orgdb_pkg"

        # Remove package source/sqlite leftovers in MTD root only if they match this OrgDb.
        remove_path_if_exists "$mtdir/$orgdb_pkg"
        remove_path_if_exists "$mtdir/${orgdb_pkg%.db}.sqlite"
    else
        echo "[CLEAN] No OrgDb package listed for TaxID ${taxid} in HostSpecies.csv."
    fi

    echo "[CLEAN] Done."
    echo "------------------------------------------------------------"
}


# ------------------------------------------------------------
# Shared Kraken2 taxonomy cache
# These files are universal and very large, so we do not remove
# or re-download them for every species/reference.
# ------------------------------------------------------------
validate_kraken_taxonomy_dir() {
    local taxdir="$1"
    local required=(names.dmp nodes.dmp nucl_gb.accession2taxid nucl_wgs.accession2taxid)
    local missing=0

    for f in "${required[@]}"; do
        if [[ ! -s "$taxdir/$f" ]]; then
            missing=1
            echo "[KRAKEN-CACHE] Missing or empty: $taxdir/$f"
        fi
    done

    [[ "$missing" == "0" ]]
}

copy_kraken_taxonomy_cache_to_db() {
    local cache_taxdir="$1"
    local db_taxdir="$2"

    rm -rf "$db_taxdir"
    mkdir -p "$db_taxdir"

    # Prefer hardlinks to avoid duplicating ~8+ GB of taxonomy files.
    # If source and destination are on different filesystems, fall back to normal copy.
    if cp -al "$cache_taxdir/." "$db_taxdir/" 2>/dev/null; then
        echo "[KRAKEN-CACHE] Linked taxonomy cache into DB taxonomy folder."
    else
        echo "[KRAKEN-CACHE] Hardlink failed; copying taxonomy cache into DB taxonomy folder."
        rm -rf "$db_taxdir"
        mkdir -p "$db_taxdir"
        cp -a "$cache_taxdir/." "$db_taxdir/"
    fi
}

prepare_kraken_taxonomy_for_db() {
    local mtdir="$1"
    local dbname="$2"
    local threads="$3"
    local cache_dir="$4"
    local rebuild_cache="$5"
    local use_cache="$6"

    if [[ "$use_cache" != "1" ]]; then
        echo "[KRAKEN-CACHE] Shared cache disabled. Downloading taxonomy directly into: $dbname"
        "$mtdir/download_kraken2_taxonomy_https.sh" --db "$dbname" --threads "$threads"
        return
    fi

    local cache_taxdir="$cache_dir/taxonomy"
    local db_taxdir="$dbname/taxonomy"

    echo "------------------------------------------------------------"
    echo "Kraken2 shared taxonomy cache"
    echo "Cache dir: $cache_dir"
    echo "DB taxonomy dir: $db_taxdir"
    echo "------------------------------------------------------------"

    if [[ "$rebuild_cache" == "1" ]]; then
        echo "[KRAKEN-CACHE] --rebuild-kraken-taxonomy-cache was used. Removing cache: $cache_dir"
        rm -rf --one-file-system "$cache_dir"
    fi

    if validate_kraken_taxonomy_dir "$cache_taxdir"; then
        echo "[KRAKEN-CACHE] Existing shared taxonomy cache is valid. Reusing it."
    else
        echo "[KRAKEN-CACHE] Shared taxonomy cache missing/incomplete. Creating it once."
        rm -rf --one-file-system "$cache_dir"
        mkdir -p "$cache_dir"
        "$mtdir/download_kraken2_taxonomy_https.sh" --db "$cache_dir" --threads "$threads"
    fi

    if ! validate_kraken_taxonomy_dir "$cache_taxdir"; then
        echo "[ERROR] Kraken taxonomy cache is still incomplete after preparation: $cache_taxdir"
        exit 1
    fi

    copy_kraken_taxonomy_cache_to_db "$cache_taxdir" "$db_taxdir"

    if ! validate_kraken_taxonomy_dir "$db_taxdir"; then
        echo "[ERROR] Kraken taxonomy folder in DB is incomplete after cache copy/link: $db_taxdir"
        exit 1
    fi

    echo "[KRAKEN-CACHE] Taxonomy ready for DB: $dbname"
    echo "------------------------------------------------------------"
}

if [[ "$CLEAN_PREVIOUS" == "1" ]]; then
    cleanup_previous_reference "$customized" "$MTDIR"
else
    echo "[INFO] --no-clean was used. Keeping previous files for TaxID ${customized}."
fi

if [[ "$CLEAN_ONLY" == "1" ]]; then
    echo "[INFO] --clean-only was used. Cleanup finished; exiting before building reference."
    exit 0
fi

# get conda path
condapath=$(head -n 1 $MTDIR/condaPath)
# activate MTD conda environment
source $condapath/etc/profile.d/conda.sh

# Custom R library used by gold OrgDb packages generated by MTD
mkdir -p "$MTDIR/custom_R_libs"
prepend_unique_path_var R_LIBS_USER "$MTDIR/custom_R_libs"
echo "[INFO] R_LIBS_USER=$R_LIBS_USER"

conda activate MTD

# Kraken2 database building - Customized
DBNAME=kraken2DB_${customized}
rm -rf $DBNAME 
mkdir -p $DBNAME
cd $DBNAME
#cp /media/me/4TB_BACKUP_LBN/Compressed/MTD/Calidris_pugnax.ASM143184v1.dna.toplevel.fa.gz .
#cp /media/me/4TB_BACKUP_LBN/Compressed/MTD/Myotis_lucifugus/Myotis_lucifugus.Myoluc2.0.dna.toplevel.fa.gz .
#wget -c $download #http://ftp.ensembl.org/pub/release-104/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
cp -u $download .
#bash ~/MTD/Customized_host.sh -d http://ftp.ensembl.org/pub/release-111/fasta/myotis_lucifugus/dna/Myotis_lucifugus.Myoluc2.0.dna.toplevel.fa.gz -c 59463 -g http://ftp.ensembl.org/pub/release-111/gtf/myotis_lucifugus/Myotis_lucifugus.Myoluc2.0.111.gtf.gz

#Morcego offline
#bash ~/MTD/Customized_host_offline.sh -d /media/me/4TB_BACKUP_LBN/Compressed/MTD/Myotis_lucifugus/Myotis_lucifugus.Myoluc2.0.dna.toplevel.fa.gz -g /media/me/4TB_BACKUP_LBN/Compressed/MTD/Myotis_lucifugus/Myotis_lucifugus.Myoluc2.0.111.gtf.gz -c 59463

#Calidris pugnax
#bash  ~/MTD/Customized_host.sh -t 20 -d https://ftp.ensembl.org/pub/release-111/fasta/calidris_pugnax/dna/Calidris_pugnax.ASM143184v1.dna.toplevel.fa.gz -c 198806 -g https://ftp.ensembl.org/pub/release-111/gtf/calidris_pugnax/Calidris_pugnax.ASM143184v1.111.gtf.gz
#Rattus norvegicus
#bash ~/MTD/Customized_host_offline.sh -d /media/me/4TB_BACKUP_LBN/Compressed/MTD/Rattus_norvegicus/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa.gz -g /media/me/4TB_BACKUP_LBN/Compressed/MTD/-c 10116

#Gallus gallus
#bash ~/MTD/Customized_host_offline.sh -d /media/me/4TB_BACKUP_LBN/Compressed/MTD/Gallus_gallus/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa.gz -g /media/me/4TB_BACKUP_LBN/Compressed/MTD/Gallus_gallus/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.111.gtf.gz -c 9031

unpigz *.fa.gz
mv *.fa genome_${customized}.fa

# Extraia o nome científico da espécie baseado no Taxon_ID
species_name=$(awk -F, -v taxid="$customized" '$1 == taxid {print $3}' "$MTDIR/HostSpecies.csv")

# Verifique se o nome da espécie foi encontrado
if [ -z "$species_name" ]; then
  echo "Error: species name not found for Taxon_ID $customized."
  exit 1
fi

# Extraia o assembly_name do cabeçalho da sequência de entrada
assembly_name=$(grep -m 1 '^>' genome_${customized}.fa | sed -n 's/.*dna:primary_assembly \([^ ]*\).*/\1/p')
echo ''
echo -e "Selected host species:\e[3m $species_name\e[0m"
#echo "Selected host species:$species_name"
echo "Taxon ID: $customized"
echo ''

# Check the header format of the FASTA file
# Add Kraken taxid to FASTA headers and append scientific name from HostSpecies.csv

fa="genome_${customized}.fa"
taxid="${customized}"
host_csv="${MTDIR}/HostSpecies.csv"

if [[ ! -s "$fa" ]]; then
    echo "ERROR: FASTA file not found or empty: $fa"
    exit 1
fi

if [[ -z "${MTDIR:-}" ]]; then
    echo "ERROR: MTDIR variable is not defined."
    exit 1
fi

if [[ ! -s "$host_csv" ]]; then
    echo "ERROR: Host species CSV not found or empty:"
    echo "  $host_csv"
    exit 1
fi

# ------------------------------------------------------------
# Get scientific name from $MTDIR/HostSpecies.csv
# Expected columns:
# Taxon_ID,MartDatasets,Scientific_name,OrgDb,Common_name,kegg
# ------------------------------------------------------------

species_from_taxid=$(
    awk -F',' -v id="$taxid" '
    NR == 1 {
        for (i = 1; i <= NF; i++) {
            gsub(/\r/, "", $i)
            if ($i == "Taxon_ID") tax_col = i
            if ($i == "Scientific_name") sci_col = i
        }

        if (!tax_col || !sci_col) {
            exit 2
        }

        next
    }

    {
        gsub(/\r/, "", $0)

        tax_id = $tax_col
        sci_name = $sci_col

        gsub(/^[ \t]+|[ \t]+$/, "", tax_id)
        gsub(/^[ \t]+|[ \t]+$/, "", sci_name)

        if (tax_id == id) {
            print sci_name
            exit
        }
    }
    ' "$host_csv"
)

awk_status=$?

if [[ "$awk_status" -eq 2 ]]; then
    echo "ERROR: Could not find required columns Taxon_ID and Scientific_name in:"
    echo "  $host_csv"
    exit 1
fi

if [[ -z "$species_from_taxid" ]]; then
    echo "WARNING: TaxID ${taxid} was not found in:"
    echo "  $host_csv"

    if [[ -n "${species_name:-}" ]]; then
        species_from_taxid="$species_name"
        echo "Using existing species_name variable: $species_from_taxid"
    else
        species_from_taxid="TaxID_${taxid}"
        echo "Using fallback name: $species_from_taxid"
    fi
else
    echo "Detected scientific name from HostSpecies.csv:"
    echo "  TaxID ${taxid} -> ${species_from_taxid}"
fi

# ------------------------------------------------------------
# Add Kraken taxid while preserving original FASTA headers
# ------------------------------------------------------------

echo "Adding Kraken taxid ${taxid} to FASTA headers..."
echo "Appending organism name: $species_from_taxid"
echo "Original FASTA headers will be preserved."

tmp=$(mktemp)

awk -v taxid="$taxid" -v organism="$species_from_taxid" '
  /^>/ {
    h = $0

    # Remove initial ">"
    sub(/^>/, "", h)

    # Remove old Kraken taxid prefix if present
    sub(/^kraken:taxid\|[0-9]+\|/, "", h)

    # Remove previous organism tag at the end, if present
    sub(/[[:space:]]+\[organism=[^]]+\][[:space:]]*$/, "", h)

    print ">kraken:taxid|" taxid "|" h " [organism=" organism "]"
    next
  }

  { print }
' "$fa" > "$tmp" && mv "$tmp" "$fa"

# ------------------------------------------------------------
# Validate result
# ------------------------------------------------------------

total_headers=$(grep -c '^>' "$fa")
tagged_headers=$(grep -c "^>kraken:taxid|${taxid}|" "$fa")

echo "Total FASTA headers: $total_headers"
echo "Headers with Kraken taxid ${taxid}: $tagged_headers"

if [[ "$total_headers" -ne "$tagged_headers" ]]; then
    echo "ERROR: Not all FASTA headers received the Kraken taxid."
    exit 1
fi

echo "Done."

echo "Header modification complete."
rm -rf $MTDIR/blastdb_$customized
mkdir -p $MTDIR/blastdb_$customized
#cp genome_${customized}.fa $MTDIR/blastdb_$customized 

cd ..
# Prepare Kraken2 taxonomy using a shared universal cache when enabled.
# This avoids re-downloading huge accession2taxid files for every custom host.
prepare_kraken_taxonomy_for_db     "$MTDIR"     "$DBNAME"     "$threads"     "$KRAKEN_TAXONOMY_CACHE"     "$REBUILD_KRAKEN_TAXONOMY_CACHE"     "$USE_KRAKEN_TAXONOMY_CACHE"

kraken2-build --add-to-library $DBNAME/genome_${customized}.fa --threads $threads --db $DBNAME
kraken2-build --build --threads $threads --db $DBNAME
# download host GTF
#wget -c $gtf -P ref_${customized} -O ref_${customized}.gtf.gz
echo "Copying GTF file to ref_${customized}/ref_${customized}.gtf.gz"
rm -rf ref_${customized}
mkdir -p ref_${customized}
cp $gtf ref_${customized}
cd ref_${customized}
mv *.gtf.gz ref_${customized}.gtf.gz

#cp /media/me/4TB_BACKUP_LBN/Compressed/MTD/Calidris_pugnax.ASM143184v1.111.gtf.gz .
#cp /media/me/4TB_BACKUP_LBN/Compressed/MTD/Myotis_lucifugus/Myotis_lucifugus.Myoluc2.0.111.gtf.gz .
cd ..
echo "Building indexes for hisat2"
rm -rf hisat2_index_${customized}
mkdir -p hisat2_index_${customized}
cd hisat2_index_${customized}
cp ../ref_${customized}/ref_${customized}.gtf.gz .
gzip -d *.gtf.gz
mv *.gtf genome.gtf
python $dir/Installation/hisat2_extract_splice_sites.py genome.gtf > genome.ss
python $dir/Installation/hisat2_extract_exons.py genome.gtf > genome.exon
if [[ ! -s "../$DBNAME/genome_${customized}.fa" ]]; then
    echo "ERROR: Genome FASTA expected for HISAT2 was not found: ../$DBNAME/genome_${customized}.fa"
    exit 1
fi
cp ../$DBNAME/genome_${customized}.fa genome.fa
hisat2-build -p $threads --exon genome.exon --ss genome.ss genome.fa genome_tran
cd ..

echo "Creating blast databases for custom reference $customized"
mkdir -p $MTDIR/blastdb_$customized
cd $MTDIR/blastdb_$customized
#cp $DBNAME/genome_${customized}.fa .
cp $download . 
gunzip *.fa.gz 
mv *.fa blastdb_$customized

makeblastdb -in $MTDIR/blastdb_$customized/blastdb_$customized -dbtype nucl -out $MTDIR/blastdb_$customized/blastdb_$customized -parse_seqids

echo "Creating gold reference-matched OrgDb package"
echo -e "Selected host species:\e[3m $species_name\e[0m"
echo "Taxon ID: $customized"
echo ''

if [[ "$SKIP_ORGDB" == "1" ]]; then
    echo "[WARNING] Gold OrgDb creation skipped."
else
    GOLD_ORGDB_SCRIPT="$MTDIR/aux_scripts/orgdb/install_gold_orgdb_from_hostspecies.sh"

    if [[ ! -s "$GOLD_ORGDB_SCRIPT" ]]; then
        echo "ERROR: Gold OrgDb installer not found:"
        echo "  $GOLD_ORGDB_SCRIPT"
        echo "Install/copy the new orgdb helper scripts into $MTDIR/aux_scripts/orgdb first."
        exit 1
    fi

    GOLD_ORGDB_ARGS=(
        --taxid "$customized"
        --hostspecies "$MTDIR/HostSpecies.csv"
        --gtf "$MTDIR/ref_${customized}/ref_${customized}.gtf.gz"
        --genome "$MTDIR/kraken2DB_${customized}/genome_${customized}.fa"
        --protein-fasta "$protein_fasta"
        --lib "$MTDIR/custom_R_libs"
        --build-dir "$MTDIR/build/orgdb_gold/${customized}"
        --threads "$threads"
        --version "$ORGDB_VERSION"
    )

    if [[ "$SKIP_IF_COMPATIBLE" == "1" ]]; then
        GOLD_ORGDB_ARGS+=(--skip-if-compatible)
    fi

    if [[ "$FORCE_ORGDB" == "1" ]]; then
        GOLD_ORGDB_ARGS+=(--force)
    fi

    echo "[INFO] Running gold OrgDb installer:"
    printf '  %q' bash "$GOLD_ORGDB_SCRIPT" "${GOLD_ORGDB_ARGS[@]}"
    echo

    if ! bash "$GOLD_ORGDB_SCRIPT" "${GOLD_ORGDB_ARGS[@]}"; then
        echo "ERROR: Gold OrgDb installer failed. Stopping customized host reference build."
        exit 1
    fi

    prepend_unique_path_var R_LIBS_USER "$MTDIR/custom_R_libs"

    echo "[INFO] Gold OrgDb check:"
    conda deactivate
    conda activate R412
    Rscript - <<'RS'
cat("R_LIBS_USER=", Sys.getenv("R_LIBS_USER"), "\n", sep="")
cat(".libPaths():\n")
print(.libPaths())
libs <- strsplit(Sys.getenv("R_LIBS_USER"), .Platform$path.sep, fixed = TRUE)[[1]]
libs <- unique(libs[nzchar(libs) & dir.exists(libs)])
hits <- unique(unlist(lapply(libs, function(p) list.files(p, pattern = "^org[.].*[.]eg[.]db$", full.names = TRUE))))
cat("Gold/custom OrgDb packages found in R_LIBS_USER:\n")
print(hits)
RS
    conda deactivate
fi

echo "Customized host reference building is done"
