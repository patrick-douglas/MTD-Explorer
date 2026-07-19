#!/bin/bash
# v5: gold OrgDb workflow + safe cleanup + shared Kraken taxonomy cache + fixed HostSpecies awk

dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
MTDIR="$dir"

condapath=~/miniconda3
threads="$(nproc)"

protein_fasta=""
scientific_name_arg=""
REFERENCE_MODE=""
offline_files_folder=""

REFERENCE_RESOLVER="$MTDIR/aux_scripts/host_reference/resolve_host_reference_from_csv.py"
HOSTSPECIES_ENTRY_HELPER="$MTDIR/aux_scripts/host_reference/ensure_hostspecies_entry.py"

CSV_REFERENCE_FOUND=0
CSV_TAXID=""
CSV_SCIENTIFIC_NAME=""
CSV_REFERENCE_TAXID=""
CSV_REFERENCE_SCIENTIFIC_NAME=""
CSV_ENSEMBL_NAME=""
CSV_ENSEMBL_DIVISION=""
CSV_ENSEMBL_RELEASE=""
CSV_ASSEMBLY=""
CSV_GENOME_URL=""
CSV_GTF_URL=""
CSV_PEP_URL=""
CSV_REFERENCE_STATUS=""


EGGNOG_DB_CACHE=""
CUSTOM_HOST_CACHE=""
RESOLVED_REFERENCE=""

SKIP_ORGDB=0
FORCE_ORGDB=0
SKIP_IF_COMPATIBLE=1
ORGDB_VERSION="1.0.0"
CLEAN_PREVIOUS=1
CLEAN_ONLY=0
USE_KRAKEN_TAXONOMY_CACHE=1
REBUILD_KRAKEN_TAXONOMY_CACHE=0
KRAKEN_TAXONOMY_CACHE=""
MTD_ORGDB_ENV="${MTD_ORGDB_ENV:-MTD_orgdb}"
export MTD_ORGDB_ENV

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
# ------------------------------------------------------------
# Persistent installation cache
# ------------------------------------------------------------

load_persistent_installation_cache() {
    local cache_path_file="$MTDIR/offlineCachePath"

    # -o remains available as an explicit override.
    if [[ -z "$offline_files_folder" ]]; then
        if [[ ! -s "$cache_path_file" ]]; then
            echo "[ERROR] Persistent MTD cache path was not found."
            echo "[ERROR] Expected file:"
            echo "  $cache_path_file"
            echo
            echo "[ERROR] Run the updated MTD installer first or provide:"
            echo "  --offline-folder /path/to/cache"
            exit 1
        fi

        IFS= read -r offline_files_folder < "$cache_path_file"
        offline_files_folder="${offline_files_folder%$'\r'}"
    fi

    if [[ -z "$offline_files_folder" ]]; then
        echo "[ERROR] The persistent cache path is empty."
        exit 1
    fi

    if [[ -e "$offline_files_folder" && ! -d "$offline_files_folder" ]]; then
        echo "[ERROR] Cache path exists but is not a directory:"
        echo "  $offline_files_folder"
        exit 1
    fi

    if ! mkdir -p "$offline_files_folder"; then
        echo "[ERROR] Could not create or access the persistent cache:"
        echo "  $offline_files_folder"
        exit 1
    fi

    offline_files_folder="$(readlink -f "$offline_files_folder")"

    if [[ ! -w "$offline_files_folder" ]]; then
        echo "[ERROR] Persistent cache is not writable:"
        echo "  $offline_files_folder"
        exit 1
    fi

    if [[ -z "$KRAKEN_TAXONOMY_CACHE" ]]; then
        KRAKEN_TAXONOMY_CACHE="$offline_files_folder/Kraken2_taxonomy_cache"
    fi

    EGGNOG_DB_CACHE="$offline_files_folder/eggNOG/emapperdb-5.0.2"
    CUSTOM_HOST_CACHE="$offline_files_folder/Customized_hosts/$customized"

    mkdir -p \
        "$KRAKEN_TAXONOMY_CACHE" \
        "$EGGNOG_DB_CACHE" \
        "$CUSTOM_HOST_CACHE/genome" \
        "$CUSTOM_HOST_CACHE/annotation" \
        "$CUSTOM_HOST_CACHE/protein"

    echo "------------------------------------------------------------"
    echo "[CACHE] Persistent MTD installation cache"
    echo "[CACHE] Root:             $offline_files_folder"
    echo "[CACHE] Kraken taxonomy:  $KRAKEN_TAXONOMY_CACHE"
    echo "[CACHE] eggNOG database:  $EGGNOG_DB_CACHE"
    echo "[CACHE] Custom host:      $CUSTOM_HOST_CACHE"
    echo "------------------------------------------------------------"
}

# ------------------------------------------------------------
# Persistent NCBI Ensembl -> Entrez mapping
# ------------------------------------------------------------
prepare_ncbi_gid_to_entrez_mapping() {
    local requested_taxid="$1"
    local reference_taxid="$2"
    local gtf_genes="$3"
    local functional_dir="$4"

    local helper
    local cache_dir
    local cache_file
    local partial_file
    local output_file
    local summary_file
    local legacy_cache

    helper="$MTDIR/aux_scripts/host_reference/build_gid_to_entrez_from_ncbi.py"

    cache_dir="$offline_files_folder/NCBI_gene"
    cache_file="$cache_dir/gene2ensembl.gz"
    partial_file="${cache_file}.part"

    output_file="$functional_dir/GID_to_ENTREZ_taxid_${requested_taxid}.tsv"
    summary_file="$functional_dir/GID_to_ENTREZ_taxid_${requested_taxid}.summary.tsv"

    legacy_cache="$MTDIR/reference_cache/ncbi_gene/gene2ensembl.gz"

    GID_TO_ENTREZ_MAP=""
    GID_TO_ENTREZ_SUMMARY=""
    NCBI_GENE2ENSEMBL_CACHE="$cache_file"
    ENTREZ_REFERENCE_TAXID="$reference_taxid"

    if [[ ! -s "$helper" ]]; then
        echo "[WARNING] NCBI Ensembl-to-Entrez helper was not found:"
        echo "  $helper"
        echo "[WARNING] Official KEGG input mapping will be unavailable."
        return 0
    fi

    mkdir -p "$cache_dir" "$functional_dir"

    # Adopt the older local cache when available, avoiding another 274 MB
    # download after upgrading an existing MTD Explorer installation.
    if [[ ! -s "$cache_file" && -s "$legacy_cache" ]]; then
        echo "[NCBI GENE] Importing existing local cache:"
        echo "  Source:      $legacy_cache"
        echo "  Destination: $cache_file"

        if ! cp --reflink=auto -f "$legacy_cache" "$cache_file" 2>/dev/null; then
            cp -f "$legacy_cache" "$cache_file"
        fi
    fi

    if [[ -s "$cache_file" ]]; then
        if gzip -t "$cache_file" >/dev/null 2>&1; then
            echo "[NCBI GENE] Reusing persistent gene2ensembl cache:"
            echo "  $cache_file"
        else
            echo "[WARNING] Removing corrupt NCBI gene2ensembl cache:"
            echo "  $cache_file"
            rm -f "$cache_file"
        fi
    fi

    if [[ ! -s "$cache_file" ]]; then
        echo "[NCBI GENE] Downloading official gene2ensembl.gz"
        echo "  Destination: $cache_file"

        rm -f "$partial_file"

        if command -v aria2c >/dev/null 2>&1; then
            if ! aria2c \
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
                --dir "$cache_dir" \
                --out "$(basename "$partial_file")" \
                "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2ensembl.gz"
            then
                echo "[WARNING] NCBI gene2ensembl download failed."
                rm -f "$partial_file"
                return 0
            fi

        elif command -v curl >/dev/null 2>&1; then
            if ! curl \
                --fail \
                --location \
                --connect-timeout 30 \
                --retry 10 \
                --retry-delay 10 \
                --retry-connrefused \
                --output "$partial_file" \
                "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2ensembl.gz"
            then
                echo "[WARNING] NCBI gene2ensembl download failed."
                rm -f "$partial_file"
                return 0
            fi

        elif command -v wget >/dev/null 2>&1; then
            if ! wget \
                --tries=20 \
                --waitretry=10 \
                --timeout=60 \
                --output-document="$partial_file" \
                "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2ensembl.gz"
            then
                echo "[WARNING] NCBI gene2ensembl download failed."
                rm -f "$partial_file"
                return 0
            fi

        else
            echo "[WARNING] aria2c, curl and wget are unavailable."
            echo "[WARNING] NCBI Ensembl-to-Entrez mapping was skipped."
            return 0
        fi

        if ! gzip -t "$partial_file" >/dev/null 2>&1; then
            echo "[WARNING] Downloaded gene2ensembl.gz failed validation."
            rm -f "$partial_file"
            return 0
        fi

        mv -f "$partial_file" "$cache_file"

        echo "[OK] Persistent NCBI gene2ensembl cache created:"
        echo "  $cache_file"
    fi

    echo "============================================================"
    echo "[NCBI GENE MAPPING]"
    echo "Requested TaxID: $requested_taxid"
    echo "Reference TaxID: $reference_taxid"
    echo "GTF gene list:   $gtf_genes"
    echo "============================================================"

    if ! conda run \
        --no-capture-output \
        -n MTD \
        python3 "$helper" \
            --gene2ensembl "$cache_file" \
            --reference-taxid "$reference_taxid" \
            --requested-taxid "$requested_taxid" \
            --gtf-genes "$gtf_genes" \
            --output "$output_file" \
            --summary "$summary_file"
    then
        echo "[WARNING] No usable NCBI Ensembl-to-Entrez map was created."
        echo "[WARNING] Host reference creation will continue."
        echo "[WARNING] Official KEGG enrichment may be unavailable."
        rm -f "$output_file" "$summary_file"
        return 0
    fi

    if [[ -s "$output_file" ]]; then
        GID_TO_ENTREZ_MAP="$output_file"
    fi

    if [[ -s "$summary_file" ]]; then
        GID_TO_ENTREZ_SUMMARY="$summary_file"
    fi
}

is_url() {
    case "$1" in
        http://*|https://*|ftp://*)
            return 0
            ;;
        *)
            return 1
            ;;
    esac
}

validate_reference_filename() {
    local kind="$1"
    local filename="$2"

    case "$kind" in
        genome)
            case "$filename" in
                *.fa|*.fna|*.fasta|*.fa.gz|*.fna.gz|*.fasta.gz)
                    return 0
                    ;;
            esac
            ;;

        gtf)
            case "$filename" in
                *.gtf|*.gtf.gz)
                    return 0
                    ;;
            esac
            ;;

        protein)
            case "$filename" in
                *.fa|*.faa|*.fna|*.fasta|*.pep|\
                *.fa.gz|*.faa.gz|*.fna.gz|*.fasta.gz|*.pep.gz)
                    return 0
                    ;;
            esac
            ;;
    esac

    echo "[ERROR] Unsupported $kind filename:"
    echo "  $filename"
    return 1
}

validate_cached_reference() {
    local cached_file="$1"
    local original_filename="$2"

    if [[ ! -s "$cached_file" ]]; then
        echo "[ERROR] Cached reference is missing or empty:"
        echo "  $cached_file"
        return 1
    fi

    case "$original_filename" in
        *.gz)
            if ! gzip -t "$cached_file" >/dev/null 2>&1; then
                echo "[ERROR] Cached gzip file is corrupt or incomplete:"
                echo "  $cached_file"
                return 1
            fi
            ;;
    esac

    return 0
}

download_reference_to_cache() {
    local url="$1"
    local destination="$2"
    local original_filename="$3"
    local partial_file="${destination}.part"

    mkdir -p "$(dirname "$destination")"

    if [[ -s "$destination" ]]; then
        if validate_cached_reference "$destination" "$original_filename"; then
            echo "[CACHE] Reusing downloaded reference:"
            echo "  $destination"
            return 0
        fi

        echo "[WARNING] Removing invalid cached reference:"
        echo "  $destination"
        rm -f "$destination"
    fi

    echo "[DOWNLOAD] URL:"
    echo "  $url"
    echo "[DOWNLOAD] Cache destination:"
    echo "  $destination"

    if command -v aria2c >/dev/null 2>&1; then
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
            --dir "$(dirname "$partial_file")" \
            --out "$(basename "$partial_file")" \
            "$url"

    elif command -v curl >/dev/null 2>&1; then
        curl \
            --fail \
            --location \
            --connect-timeout 30 \
            --retry 10 \
            --retry-delay 20 \
            --retry-connrefused \
            --continue-at - \
            --output "$partial_file" \
            "$url"

    elif command -v wget >/dev/null 2>&1; then
        wget \
            --continue \
            --tries=50 \
            --waitretry=20 \
            --timeout=60 \
            --read-timeout=60 \
            --retry-connrefused \
            --output-document="$partial_file" \
            "$url"

    else
        echo "[ERROR] aria2c, curl and wget are unavailable."
        return 1
    fi

    if ! validate_cached_reference "$partial_file" "$original_filename"; then
        echo "[ERROR] Downloaded reference failed validation."
        rm -f "$partial_file"
        return 1
    fi

    mv -f "$partial_file" "$destination"

    echo "[OK] Reference stored in persistent cache:"
    echo "  $destination"
}

resolve_reference_input() {
    local label="$1"
    local kind="$2"
    local source_value="$3"
    local cache_subdirectory="$4"

    local clean_source
    local filename
    local destination
    local source_absolute

    clean_source="${source_value%%\?*}"
    filename="$(basename "$clean_source")"

    if [[ -z "$filename" || "$filename" == "." || "$filename" == "/" ]]; then
        echo "[ERROR] Could not determine filename for $label:"
        echo "  $source_value"
        exit 1
    fi

    if ! validate_reference_filename "$kind" "$filename"; then
        exit 1
    fi

    destination="$CUSTOM_HOST_CACHE/$cache_subdirectory/$filename"

    if is_url "$source_value"; then
        if ! download_reference_to_cache \
            "$source_value" \
            "$destination" \
            "$filename"; then

            echo "[ERROR] Failed to download $label."
            exit 1
        fi
    else
        if [[ ! -s "$source_value" ]]; then
            echo "[ERROR] Local $label was not found or is empty:"
            echo "  $source_value"
            exit 1
        fi

        source_absolute="$(readlink -f "$source_value")"

        if [[ "$source_absolute" != "$destination" ]]; then
            if [[ -s "$destination" ]] &&
               validate_cached_reference "$destination" "$filename"; then

                echo "[CACHE] Reusing cached $label:"
                echo "  $destination"
            else
                rm -f "$destination"

                echo "[CACHE] Copying local $label into persistent cache:"
                echo "  Source:      $source_absolute"
                echo "  Destination: $destination"

                if ! cp --reflink=auto -f "$source_absolute" "$destination" 2>/dev/null; then
                    cp -f "$source_absolute" "$destination"
                fi
            fi
        fi

        if ! validate_cached_reference "$destination" "$filename"; then
            echo "[ERROR] Cached copy of $label failed validation."
            exit 1
        fi
    fi

    RESOLVED_REFERENCE="$(readlink -f "$destination")"
}

materialize_fasta_reference() {
    local source_file="$1"
    local destination_file="$2"
    local temporary_file="${destination_file}.tmp"

    rm -f "$temporary_file"

    case "$source_file" in
        *.gz)
            if ! gzip -dc -- "$source_file" > "$temporary_file"; then
                echo "[ERROR] Could not decompress FASTA:"
                echo "  $source_file"
                rm -f "$temporary_file"
                return 1
            fi
            ;;

        *)
            if ! cp -f -- "$source_file" "$temporary_file"; then
                echo "[ERROR] Could not copy FASTA:"
                echo "  $source_file"
                rm -f "$temporary_file"
                return 1
            fi
            ;;
    esac

    if [[ ! -s "$temporary_file" ]] ||
       ! grep -qm1 '^>' "$temporary_file"; then

        echo "[ERROR] Materialized FASTA is empty or invalid:"
        echo "  $source_file"
        rm -f "$temporary_file"
        return 1
    fi

    mv -f "$temporary_file" "$destination_file"
}

prepare_cached_gtf_gz() {
    local source_file="$1"
    local destination_file="$2"
    local temporary_file="${destination_file}.tmp"

    rm -f "$temporary_file"

    case "$source_file" in
        *.gz)
            if ! gzip -t "$source_file" >/dev/null 2>&1; then
                echo "[ERROR] Invalid compressed GTF:"
                echo "  $source_file"
                return 1
            fi

            cp -f -- "$source_file" "$temporary_file"
            ;;

        *)
            if ! gzip -c -- "$source_file" > "$temporary_file"; then
                echo "[ERROR] Could not compress GTF:"
                echo "  $source_file"
                rm -f "$temporary_file"
                return 1
            fi
            ;;
    esac

    if ! gzip -t "$temporary_file" >/dev/null 2>&1; then
        echo "[ERROR] Prepared GTF failed gzip validation:"
        echo "  $temporary_file"
        rm -f "$temporary_file"
        return 1
    fi

    mv -f "$temporary_file" "$destination_file"
}

# Função para mostrar instruções de uso
usage() {
echo "Usage:"
    echo ""
    echo "Automatic mode for a species in HostSpecies.csv:"
    echo "  bash $0 --ncbi-taxon-id TAXID"
    echo ""
    echo "Manual mode for a species outside HostSpecies.csv:"
    echo "  bash $0 \\"
    echo "    --ncbi-taxon-id TAXID \\"
    echo "    --scientific-name \"Genus species\" \\"
    echo "    --genome genome.fa.gz \\"
    echo "    --gtf-file annotation.gtf.gz \\"
    echo "    --protein-fasta proteins.fa.gz"
    echo ""
    echo "Required options:"
    echo " -d, --genome           Manual genome path or URL."
    echo "                         Not required when TaxID exists in the curated CSV."
    echo " -g, --gtf-file         Manual GTF path or URL."
    echo "                         Not required when TaxID exists in the curated CSV."
    echo " -p, --protein-fasta    Manual protein FASTA path or URL."
    echo "                         Automatically resolved for curated species."
    echo "                 Supported: .gtf and .gtf.gz"
    echo "  -c, --ncbi-taxon-id    NCBI Taxon ID of the host species"
    echo ""
    echo "Gold OrgDb options:"
    echo "      --skip-orgdb       Skip gold OrgDb creation"
    echo "Skip only the R OrgDb package construction."
    echo "Representative proteins, eggNOG annotations and the master ssGSEA GMT are still generated."
    echo "      --force-orgdb      Rebuild gold OrgDb even if an existing compatible OrgDb is detected"
    echo "      --no-skip-if-compatible"
    echo "                         Do not skip when an existing OrgDb appears compatible"
    echo "      --orgdb-version    Version for the generated OrgDb package [default: 0.3.0]"
    echo "      --no-clean         Do not remove previous files for this TaxID/reference before starting"
    echo "      --clean-only       Remove previous files for this TaxID/reference and exit"
    echo " -s, --scientific-name NAME"
    echo "                         Scientific name for a TaxID absent from"
    echo "                         HostSpecies.csv. Required in manual mode"
    echo "                         when creating a new species entry."
    echo ""
    echo "Kraken taxonomy cache options:"
    echo "      --kraken-taxonomy-cache DIR"
    echo "                         Shared Kraken2 taxonomy cache directory [default: $MTDIR/kraken2_taxonomy_cache]"
    echo "      --no-kraken-taxonomy-cache"
    echo "                         Do not use shared cache; download taxonomy into the species DB as before"
    echo "      --rebuild-kraken-taxonomy-cache"
    echo "                         Delete and rebuild the shared Kraken2 taxonomy cache before using it"
    echo ""
    echo "Installation cache:"
    echo " -o, --offline-folder Override the persistent cache recorded by the MTD installer."
    echo "                       Normally this option is unnecessary because the path is"
    echo "                       read automatically from: $MTDIR/offlineCachePath"
    echo ""
    echo "Other:"
    echo "      --help             Show this help message"
    echo ""
echo "Examples:"
echo ""
echo "Local files:"
echo " bash $0 \\"
echo "   --genome /path/to/species.dna.toplevel.fa.gz \\"
echo "   --gtf-file /path/to/species.annotation.gtf.gz \\"
echo "   --protein-fasta /path/to/species.pep.all.fa.gz \\"
echo "   --ncbi-taxon-id 6526"
echo ""
echo "URLs:"
echo " bash $0 \\"
echo "   --genome https://server/species.dna.toplevel.fa.gz \\"
echo "   --gtf-file https://server/species.annotation.gtf.gz \\"
echo "   --protein-fasta https://server/species.pep.all.fa.gz \\"
echo "   --ncbi-taxon-id 6526"
echo ""    exit 1
}

# Define as opções curtas e longas
OPTIONS=c:d:g:o:p:s:
LONGOPTIONS=ncbi-taxon-id:,genome:,gtf-file:,offline-folder:,protein-fasta:,scientific-name:,skip-orgdb,force-orgdb,no-skip-if-compatible,orgdb-version:,no-clean,clean-only,kraken-taxonomy-cache:,no-kraken-taxonomy-cache,rebuild-kraken-taxonomy-cache,help

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
        -s|--scientific-name)
            scientific_name_arg="$2"
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

# ------------------------------------------------------------
# Resolve automatic or manual reference mode
# ------------------------------------------------------------

if [[ -z "${customized:-}" ]]; then
    echo "[ERROR] --ncbi-taxon-id is required."
    usage
fi

if ! [[ "$customized" =~ ^[1-9][0-9]*$ ]]; then
    echo "[ERROR] Invalid NCBI Taxon ID:"
    echo " $customized"
    exit 1
fi

HOSTSPECIES_CSV="$MTDIR/HostSpecies.csv"

if [[ ! -s "$HOSTSPECIES_CSV" ]]; then
    echo "[ERROR] HostSpecies.csv was not found or is empty:"
    echo " $HOSTSPECIES_CSV"
    exit 1
fi

if [[ ! -s "$REFERENCE_RESOLVER" ]]; then
    echo "[ERROR] Automatic reference resolver was not found:"
    echo " $REFERENCE_RESOLVER"
    exit 1
fi

if [[ ! -s "$HOSTSPECIES_ENTRY_HELPER" ]]; then
    echo "[ERROR] HostSpecies entry helper was not found:"
    echo " $HOSTSPECIES_ENTRY_HELPER"
    exit 1
fi

# Genome and GTF determine whether this is automatic or manual mode.
if [[ -z "${download:-}" && -z "${gtf:-}" ]]; then
    if [[ -n "${protein_fasta:-}" ]]; then
        echo "[ERROR] --protein-fasta was provided without --genome and --gtf-file."
        echo "[ERROR] Use either:"
        echo "  automatic mode: only --ncbi-taxon-id"
        echo "or:"
        echo "  manual mode: --genome + --gtf-file + --protein-fasta"
        exit 1
    fi

    REFERENCE_MODE="automatic"

elif [[ -n "${download:-}" && -n "${gtf:-}" ]]; then
    REFERENCE_MODE="manual"

else
    echo "[ERROR] Partial manual reference input detected."
    echo "[ERROR] --genome and --gtf-file must be supplied together."
    echo "[ERROR] Do not mix automatic and partial manual references."
    exit 1
fi

echo "------------------------------------------------------------"
echo "[REFERENCE] Input mode: $REFERENCE_MODE"
echo "[REFERENCE] Requested TaxID: $customized"
echo "------------------------------------------------------------"

if [[ "$REFERENCE_MODE" == "automatic" ]]; then
    reference_env="$(mktemp)"

    if python3 "$REFERENCE_RESOLVER" \
        --csv "$HOSTSPECIES_CSV" \
        --taxid "$customized" \
        > "$reference_env"; then

        # The resolver emits shell-escaped assignments.
        # shellcheck disable=SC1090
        source "$reference_env"
        rm -f "$reference_env"

    else
        resolver_status=$?
        rm -f "$reference_env"

        if [[ "$resolver_status" -eq 3 ]]; then
            echo "[ERROR] TaxID $customized is not present in the curated CSV."
            echo
            echo "Use manual mode:"
            echo
            echo "  bash $0 \\"
            echo "    --ncbi-taxon-id $customized \\"
            echo "    --scientific-name \"Genus species\" \\"
            echo "    --genome /path/genome.fa.gz \\"
            echo "    --gtf-file /path/annotation.gtf.gz \\"
            echo "    --protein-fasta /path/proteins.fa.gz"
            exit 1
        fi

        echo "[ERROR] Could not resolve the automatic reference."
        exit 1
    fi

    download="$CSV_GENOME_URL"
    gtf="$CSV_GTF_URL"
    protein_fasta="$CSV_PEP_URL"

    echo "[REFERENCE] Scientific name: $CSV_SCIENTIFIC_NAME"
    echo "[REFERENCE] Assembly:        $CSV_ASSEMBLY"
    echo "[REFERENCE] Ensembl release: $CSV_ENSEMBL_RELEASE"
    echo "[REFERENCE] Division:        $CSV_ENSEMBL_DIVISION"
    echo "[REFERENCE] Genome URL:      $download"
    echo "[REFERENCE] GTF URL:         $gtf"
    echo "[REFERENCE] Protein URL:     $protein_fasta"

    if [[ "$CSV_REFERENCE_TAXID" != "$customized" ]]; then
        echo "[WARNING] Requested TaxID and reference TaxID differ:"
        echo " Requested TaxID: $customized"
        echo " Reference TaxID: $CSV_REFERENCE_TAXID"
        echo " Reference name:  $CSV_REFERENCE_SCIENTIFIC_NAME"
    fi

else
    # Manual mode. Ensure the TaxID exists in HostSpecies.csv.
    entry_args=(
        --csv "$HOSTSPECIES_CSV"
        --taxid "$customized"
        --add-if-missing
    )

    if [[ -n "$scientific_name_arg" ]]; then
        entry_args+=(
            --scientific-name "$scientific_name_arg"
        )
    fi

    if ! python3 "$HOSTSPECIES_ENTRY_HELPER" \
        "${entry_args[@]}"; then

        echo
        echo "[ERROR] Manual species registration failed."
        echo "[ERROR] For a TaxID outside HostSpecies.csv, provide:"
        echo " --scientific-name \"Genus species\""
        exit 1
    fi

    echo "[REFERENCE] Manual genome:  $download"
    echo "[REFERENCE] Manual GTF:     $gtf"

    if [[ -n "${protein_fasta:-}" ]]; then
        echo "[REFERENCE] Manual protein: $protein_fasta"
    fi
fi

load_persistent_installation_cache

resolve_reference_input \
    "host genome" \
    "genome" \
    "$download" \
    "genome"

download="$RESOLVED_REFERENCE"

resolve_reference_input \
    "GTF annotation" \
    "gtf" \
    "$gtf" \
    "annotation"

gtf="$RESOLVED_REFERENCE"

if [[ -z "${protein_fasta:-}" ]]; then
    echo "[ERROR] A protein FASTA is required to create a complete"
    echo "[ERROR] MTD custom host reference."
    echo
    echo "[ERROR] The protein FASTA is used once by Create_custom_host.sh"
    echo "[ERROR] to generate:"
    echo "  - representative proteins per gene"
    echo "  - eggNOG annotations"
    echo "  - reference master GO GMT"
    echo "  - optional custom OrgDb"
    echo
    echo "[ERROR] Provide --protein-fasta using the same annotation"
    echo "[ERROR] and release as the GTF."
    exit 1
fi

resolve_reference_input \
    "protein FASTA" \
    "protein" \
    "$protein_fasta" \
    "protein"

protein_fasta="$RESOLVED_REFERENCE"

cd "$MTDIR" || {
    echo "[ERROR] Could not enter MTD directory:"
    echo "  $MTDIR"
    exit 1
}

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
        "$mtdir/aux_scripts/Kraken2/download_kraken2_taxonomy_https.sh" --db "$dbname" --threads "$threads"
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
        "$mtdir/aux_scripts/Kraken2/download_kraken2_taxonomy_https.sh" --db "$cache_dir" --threads "$threads"
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
# Load Conda
source "$condapath/etc/profile.d/conda.sh"

# ------------------------------------------------------------
# Validate dedicated OrgDb environment
# ------------------------------------------------------------

if [[ "$SKIP_ORGDB" != "1" ]]; then
    ORGDB_CHECK_SCRIPT="$MTDIR/Installation/check_MTD_orgdb.R"

    if [[ ! -s "$ORGDB_CHECK_SCRIPT" ]]; then
        echo "[ERROR] MTD_orgdb validation script was not found:"
        echo "  $ORGDB_CHECK_SCRIPT"
        exit 1
    fi

    echo "[INFO] Validating OrgDb Conda environment:"
    echo "  $MTD_ORGDB_ENV"

    if ! conda run -n "$MTD_ORGDB_ENV" \
        Rscript "$ORGDB_CHECK_SCRIPT"
    then
        echo "[ERROR] The dedicated OrgDb environment is unavailable or incomplete:"
        echo "  $MTD_ORGDB_ENV"
        echo
        echo "Reinstall MTD or create the environment with:"
        echo "  conda env create -f $MTDIR/Installation/MTD_orgdb.yml"
        exit 1
    fi
fi

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

echo "[INFO] Preparing host genome FASTA from cache:"
echo "  $download"

if ! materialize_fasta_reference \
    "$download" \
    "genome_${customized}.fa"; then

    exit 1
fi

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
echo "[INFO] Preparing GTF annotation from cache:"
echo "  $gtf"

rm -rf "$MTDIR/ref_${customized}"
mkdir -p "$MTDIR/ref_${customized}"

if ! prepare_cached_gtf_gz \
    "$gtf" \
    "$MTDIR/ref_${customized}/ref_${customized}.gtf.gz"; then

    exit 1
fi

echo "Building indexes for hisat2"
rm -rf hisat2_index_${customized}
mkdir -p hisat2_index_${customized}
cd hisat2_index_${customized}
if ! gzip -dc \
    "$MTDIR/ref_${customized}/ref_${customized}.gtf.gz" \
    > genome.gtf; then

    echo "[ERROR] Could not decompress GTF for HISAT2."
    exit 1
fi
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
rm -rf "$MTDIR/blastdb_${customized}"
mkdir -p "$MTDIR/blastdb_${customized}"
cd "$MTDIR/blastdb_${customized}" || exit 1

if ! materialize_fasta_reference \
    "$download" \
    "blastdb_${customized}"; then

    exit 1
fi

makeblastdb \
    -in "$MTDIR/blastdb_${customized}/blastdb_${customized}" \
    -dbtype nucl \
    -out "$MTDIR/blastdb_${customized}/blastdb_${customized}" \
    -parse_seqids

echo "============================================================"
echo "Creating reference-matched functional host resources"
echo "Host: $species_name"
echo "Taxon ID: $customized"
echo "============================================================"

GOLD_ORGDB_SCRIPT="$MTDIR/aux_scripts/orgdb/install_gold_orgdb_from_hostspecies.sh"

if [[ ! -s "$GOLD_ORGDB_SCRIPT" ]]; then
    echo "[ERROR] Gold OrgDb installer not found:"
    echo "  $GOLD_ORGDB_SCRIPT"
    exit 1
fi

GOLD_ORGDB_ARGS=(
    --taxid "$customized"
    --hostspecies "$MTDIR/HostSpecies.csv"
    --gtf "$MTDIR/ref_${customized}/ref_${customized}.gtf.gz"
    --genome "$MTDIR/kraken2DB_${customized}/genome_${customized}.fa"
    --protein-fasta "$protein_fasta"
    --eggnog-db "$EGGNOG_DB_CACHE"
    --lib "$MTDIR/custom_R_libs"
    --build-dir "$MTDIR/build/orgdb_gold/${customized}"
    --threads "$threads"
    --version "$ORGDB_VERSION"
    --conda-r-env "$MTD_ORGDB_ENV"
)

if [[ "$SKIP_IF_COMPATIBLE" == "1" ]]; then
    GOLD_ORGDB_ARGS+=(--skip-if-compatible)
fi

if [[ "$FORCE_ORGDB" == "1" ]]; then
    GOLD_ORGDB_ARGS+=(--force)
fi

if [[ "$SKIP_ORGDB" == "1" ]]; then
    GOLD_ORGDB_ARGS+=(--skip-orgdb-build)
fi

echo "[INFO] Running reference-matched annotation builder:"
printf ' %q' bash "$GOLD_ORGDB_SCRIPT" "${GOLD_ORGDB_ARGS[@]}"
echo

if ! bash "$GOLD_ORGDB_SCRIPT" "${GOLD_ORGDB_ARGS[@]}"; then
    echo "[ERROR] Reference-matched annotation builder failed."
    echo "[ERROR] Stopping customized host reference build."
    exit 1
fi

# ------------------------------------------------------------
# Persist eggNOG/GO resources with the custom reference
# ------------------------------------------------------------

FUNCTIONAL_DIR="$MTDIR/ref_${customized}/functional_annotation"

ORGDB_WORK_DIR="$MTDIR/build/orgdb_gold/${customized}/work_taxid_${customized}"

ORGDB_EGGNOG_ANNOT="$ORGDB_WORK_DIR/eggnog_results/taxid_${customized}_gold_orgdb.emapper.annotations"

ORGDB_GTF_GENES="$ORGDB_WORK_DIR/gtf_gene_ids.txt"

ORGDB_REP_REPORT="$ORGDB_WORK_DIR/taxid_${customized}_representative_proteins_report.tsv"

ORGDB_EGGNOG_COMPAT="$ORGDB_WORK_DIR/eggnog_gtf_compatibility.tsv"

ORGDB_INPUT_COMPAT="$ORGDB_WORK_DIR/input_compatibility_report.tsv"

MASTER_GMT_BUILDER="$MTDIR/aux_scripts/ssGSEA/build_master_gmt_from_eggnog.py"

if [[ ! -s "$ORGDB_EGGNOG_ANNOT" ]]; then
    echo "[ERROR] Reference-matched eggNOG annotation was not found:"
    echo "  $ORGDB_EGGNOG_ANNOT"
    exit 1
fi

if [[ ! -s "$ORGDB_GTF_GENES" ]]; then
    echo "[ERROR] Reference GTF gene list was not found:"
    echo "  $ORGDB_GTF_GENES"
    exit 1
fi

if [[ ! -s "$MASTER_GMT_BUILDER" ]]; then
    echo "[ERROR] Master GMT builder was not found:"
    echo "  $MASTER_GMT_BUILDER"
    exit 1
fi

mkdir -p "$FUNCTIONAL_DIR"

REFERENCE_TAXID_FOR_ENTREZ="${CSV_REFERENCE_TAXID:-$customized}"

if [[ -z "$REFERENCE_TAXID_FOR_ENTREZ" ]]; then
    REFERENCE_TAXID_FOR_ENTREZ="$customized"
fi

prepare_ncbi_gid_to_entrez_mapping \
    "$customized" \
    "$REFERENCE_TAXID_FOR_ENTREZ" \
    "$ORGDB_GTF_GENES" \
    "$FUNCTIONAL_DIR"

PERSISTENT_EGGNOG_ANNOT="$FUNCTIONAL_DIR/taxid_${customized}_gold_orgdb.emapper.annotations"

cp -f \
    "$ORGDB_EGGNOG_ANNOT" \
    "$PERSISTENT_EGGNOG_ANNOT"

if [[ -s "$ORGDB_REP_REPORT" ]]; then
    cp -f \
        "$ORGDB_REP_REPORT" \
        "$FUNCTIONAL_DIR/"
fi

if [[ -s "$ORGDB_EGGNOG_COMPAT" ]]; then
    cp -f \
        "$ORGDB_EGGNOG_COMPAT" \
        "$FUNCTIONAL_DIR/"
fi

if [[ -s "$ORGDB_INPUT_COMPAT" ]]; then
    cp -f \
        "$ORGDB_INPUT_COMPAT" \
        "$FUNCTIONAL_DIR/"
fi

echo "============================================================"
echo "[FUNCTIONAL ANNOTATION]"
echo "TaxID: $customized"
echo "eggNOG annotations:"
echo "  $PERSISTENT_EGGNOG_ANNOT"
echo "Reference GTF genes:"
echo "  $ORGDB_GTF_GENES"
echo "Output directory:"
echo "  $FUNCTIONAL_DIR"
echo "============================================================"

if ! conda run \
    --no-capture-output \
    -n MTD \
    python3 "$MASTER_GMT_BUILDER" \
        --taxid "$customized" \
        --eggnog "$PERSISTENT_EGGNOG_ANNOT" \
        --gtf-genes "$ORGDB_GTF_GENES" \
        --outdir "$FUNCTIONAL_DIR"
then
    echo "[ERROR] Persistent master GMT generation failed."
    exit 1
fi

MASTER_GMT="$FUNCTIONAL_DIR/custom_taxid_${customized}_eggNOG_GO_master.gmt"

if [[ ! -s "$MASTER_GMT" ]]; then
    echo "[ERROR] Master GMT was not created:"
    echo "  $MASTER_GMT"
    exit 1
fi

FUNCTIONAL_MANIFEST="$FUNCTIONAL_DIR/functional_annotation_manifest.tsv"

{
    printf 'parameter\tvalue\n'
    printf 'taxid\t%s\n' "$customized"
    printf 'scientific_name\t%s\n' "$species_name"
    printf 'gtf\t%s\n' "$MTDIR/ref_${customized}/ref_${customized}.gtf.gz"
    printf 'protein_fasta_cache\t%s\n' "$protein_fasta"
    printf 'eggnog_annotations\t%s\n' "$PERSISTENT_EGGNOG_ANNOT"
    printf 'master_gmt\t%s\n' "$MASTER_GMT"
    printf 'created_at\t%s\n' "$(date -Is)"
} > "$FUNCTIONAL_MANIFEST"

(
    cd "$FUNCTIONAL_DIR" || exit 1

    find . \
        -maxdepth 1 \
        -type f \
        ! -name 'checksums.sha256' \
        -printf '%P\n' |
    sort |
    xargs -r sha256sum \
        > checksums.sha256
)

echo "[OK] Persistent host functional resources:"
echo "  $FUNCTIONAL_DIR"

if [[ "$SKIP_ORGDB" != "1" ]]; then
    prepend_unique_path_var \
        R_LIBS_USER \
        "$MTDIR/custom_R_libs"

    echo "[INFO] Gold OrgDb check:"

    conda deactivate
    conda activate "$MTD_ORGDB_ENV"

    Rscript - <<'RS'
cat("R_LIBS_USER=", Sys.getenv("R_LIBS_USER"), "\n", sep = "")
cat(".libPaths():\n")
print(.libPaths())

libs <- strsplit(
    Sys.getenv("R_LIBS_USER"),
    .Platform$path.sep,
    fixed = TRUE
)[[1]]

libs <- unique(libs[nzchar(libs) & dir.exists(libs)])

hits <- unique(
    unlist(
        lapply(
            libs,
            function(path) {
                list.files(
                    path,
                    pattern = "^org[.].*[.]eg[.]db$",
                    full.names = TRUE
                )
            }
        )
    )
)

cat("Gold/custom OrgDb packages found in R_LIBS_USER:\n")
print(hits)
RS

    conda deactivate
else
    echo "[INFO] --skip-orgdb was used."
    echo "[INFO] eggNOG annotations and master GMT were still generated."
fi

echo "Customized host reference building is done"
