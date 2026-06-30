#!/usr/bin/env bash

# ==============================================================================
# MTD installer
# ==============================================================================
# This version reorganizes the original installer into functions while preserving
# the execution order and the local Kraken2 helper-script workflow.
#
# Deliberately not using `set -e` here. The original installer continues after
# some non-critical commands fail, and enabling it would change existing behavior.
# ==============================================================================

# ------------------------------------------------------------------------------
# Defaults and global paths
# ------------------------------------------------------------------------------

kmer=""                     # Kraken2 --kmer-len, used only with --build
min_l=""                    # Kraken2 --minimizer-len, used only with --build
min_s=""                    # Kraken2 --minimizer-spaces, used only with --build
read_len=75                 # Bracken read length
threads="$(nproc)"
condapath="${HOME}/miniconda3"
offline_files_folder=""
sudo_password=""

dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"

KRAKEN_ENV_LIBEXEC=""
KRAKEN_PKG_LIBEXEC=""
kraken_build_opts=()
ORIGINAL_CHANNEL_PRIORITY=""

# ------------------------------------------------------------------------------
# Terminal formatting and messages
# ------------------------------------------------------------------------------

init_colors() {
    if command -v tput >/dev/null 2>&1 && [[ -t 1 ]]; then
        w="$(tput sgr0)"
        r="$(tput setaf 1)"
        g="$(tput setaf 2)"
        y="$(tput setaf 3)"
        p="$(tput setaf 5)"
    else
        w=""
        r=""
        g=""
        y=""
        p=""
    fi
}

print_rule() {
    printf '%s\n' "============================================================"
}

log_info() {
    printf '%s[INFO]%s %s\n' "$g" "$w" "$*"
}

log_ok() {
    printf '%s[OK]%s %s\n' "$g" "$w" "$*"
}

log_warning() {
    printf '%s[WARNING]%s %s\n' "$y" "$w" "$*" >&2
}

log_error() {
    printf '%s[ERROR]%s %s\n' "$r" "$w" "$*" >&2
}

show_progress() {
    local bar="$1"
    local percent="$2"
    local message="$3"

    echo "${g}"
    echo "MTD installation progress:"
    printf '%-20s[%s]\n' "$bar" "$percent"
    echo "$message"
    echo "${w}"
}

restore_conda_channel_priority() {
    if [[ -n "$ORIGINAL_CHANNEL_PRIORITY" ]] && command -v conda >/dev/null 2>&1; then
        conda config --set channel_priority "$ORIGINAL_CHANNEL_PRIORITY" >/dev/null 2>&1 || true
    fi
}

on_exit() {
    local exit_status=$?
    restore_conda_channel_priority

    if (( exit_status != 0 && exit_status != 130 )); then
        echo
        log_error "MTD installation failed with exit status $exit_status."
    fi
}

on_interrupt() {
    echo
    log_error "Installation stopped by user."
    exit 130
}

trap on_exit EXIT
trap on_interrupt INT TERM

# ------------------------------------------------------------------------------
# Command-line arguments and validation
# ------------------------------------------------------------------------------

usage() {
    cat <<USAGE
Usage:
  $0 -p <condapath> -o <offline_files_folder> [options]

Required:
  -p PATH   Conda installation path
  -o PATH   Folder containing the offline installation files

Optional:
  -k INT    Kraken2 k-mer length used during --build
  -m INT    Kraken2 minimizer length used during --build
  -s INT    Kraken2 minimizer spaces used during --build
  -r INT    Bracken read length (default: 75)
  -w TEXT   Sudo password used by the existing expect helper
  -h        Show this help message
USAGE
}

parse_arguments() {
    if [[ $# -eq 1 && "$1" == "-h" ]]; then
        usage
        exit 0
    fi

    if [[ $# -lt 4 ]]; then
        usage
        exit 1
    fi

    while getopts ":p:o:k:m:s:r:w:h" option; do
        case "$option" in
            p) condapath="$OPTARG" ;;
            o) offline_files_folder="$OPTARG" ;;
            k) kmer="$OPTARG" ;;
            m) min_l="$OPTARG" ;;
            s) min_s="$OPTARG" ;;
            r) read_len="$OPTARG" ;;
            w) sudo_password="$OPTARG" ;;
            h)
                usage
                exit 0
                ;;
            :)
                log_error "Option -$OPTARG requires a value."
                usage
                exit 1
                ;;
            \?)
                log_error "Unknown option: -$OPTARG"
                usage
                exit 1
                ;;
        esac
    done
}

validate_arguments() {
    if [[ -z "$condapath" || -z "$offline_files_folder" ]]; then
        usage
        exit 1
    fi

    condapath="$(readlink -f "$condapath")"
    offline_files_folder="$(readlink -f "$offline_files_folder")"

    if [[ ! -f "$condapath/etc/profile.d/conda.sh" ]]; then
        log_error "Conda initialization script not found:"
        log_error "  $condapath/etc/profile.d/conda.sh"
        exit 1
    fi

    if [[ ! -d "$offline_files_folder" ]]; then
        log_error "Offline installation folder not found:"
        log_error "  $offline_files_folder"
        exit 1
    fi

    if ! [[ "$threads" =~ ^[1-9][0-9]*$ ]]; then
        log_error "Invalid CPU thread count: $threads"
        exit 1
    fi

    if ! [[ "$read_len" =~ ^[1-9][0-9]*$ ]]; then
        log_error "Invalid Bracken read length: $read_len"
        exit 1
    fi
}

configure_paths_and_options() {
    KRAKEN_ENV_LIBEXEC="$condapath/envs/MTD/libexec"
    KRAKEN_PKG_LIBEXEC="$condapath/pkgs/kraken2-2.1.2-pl5262h7d875b9_0/libexec"

    kraken_build_opts=()

    if [[ -n "$kmer" ]]; then
        kraken_build_opts+=(--kmer-len "$kmer")
    fi
    if [[ -n "$min_l" ]]; then
        kraken_build_opts+=(--minimizer-len "$min_l")
    fi
    if [[ -n "$min_s" ]]; then
        kraken_build_opts+=(--minimizer-spaces "$min_s")
    fi
}

# ------------------------------------------------------------------------------
# General helpers
# ------------------------------------------------------------------------------

sudo_with_pass() {
    local cmd="$1"

    if command -v expect >/dev/null 2>&1; then
        expect <<EXPECT_EOF
            set timeout -1
            spawn bash -c "$cmd"
            expect {
                "*password*" {
                    send "$sudo_password\r"
                    exp_continue
                }
                eof
            }
EXPECT_EOF
    elif [[ -n "$sudo_password" ]]; then
        log_warning "expect is unavailable; using sudo -S fallback."
        printf '%s\n' "$sudo_password" | sudo -S bash -c "${cmd#sudo }"
    else
        log_warning "expect is unavailable and no sudo password was supplied."
        bash -c "$cmd"
    fi
}

safe_conda_deactivate() {
    conda deactivate >/dev/null 2>&1 || true
}

run_required_command() {
    local description="$1"
    shift

    log_info "$description"
    "$@"
    local exit_status=$?

    if (( exit_status != 0 )); then
        log_error "$description failed with exit status $exit_status."
        printf '[ERROR] Command:' >&2
        printf ' %q' "$@" >&2
        printf '\n' >&2
        exit "$exit_status"
    fi
}

activate_required_env() {
    local env_name="$1"

    if ! conda activate "$env_name"; then
        log_error "Could not activate required Conda environment: $env_name"
        exit 1
    fi

    if [[ "${CONDA_DEFAULT_ENV:-}" != "$env_name" ]]; then
        log_error "Conda reported an unexpected active environment."
        log_error "Expected: $env_name"
        log_error "Active:   ${CONDA_DEFAULT_ENV:-none}"
        exit 1
    fi

    log_ok "Activated Conda environment: $env_name"
}

require_env_command() {
    local env_name="$1"
    local command_name="$2"

    if ! conda run -n "$env_name" bash -c "command -v '$command_name' >/dev/null 2>&1"; then
        log_error "Required command '$command_name' is missing from Conda environment '$env_name'."
        exit 1
    fi

    log_ok "Found '$command_name' in Conda environment '$env_name'."
}

copy_required_file() {
    local source_file="$1"
    local destination="$2"

    if [[ ! -f "$source_file" ]]; then
        log_error "Required file not found: $source_file"
        exit 1
    fi

    if ! cp -f "$source_file" "$destination"; then
        log_error "Could not copy required file:"
        log_error "  Source:      $source_file"
        log_error "  Destination: $destination"
        exit 1
    fi
}

run_required_script() {
    local script="$1"
    shift

    if [[ ! -f "$script" ]]; then
        log_error "Required script not found: $script"
        exit 1
    fi

    chmod +x "$script"
    if ! "$script" "$@"; then
        log_error "Required script failed: $script"
        exit 1
    fi
}

retry_until_success() {
    local description="$1"
    shift

    local attempt=1
    local retry_delay="${RETRY_INITIAL_DELAY:-20}"
    local max_retry_delay="${RETRY_MAX_DELAY:-300}"
    local exit_status=0

    while true; do
        print_rule
        echo "[RETRY] $description"
        echo "[RETRY] Attempt: $attempt"
        printf '[RETRY] Command:'
        printf ' %q' "$@"
        printf '\n'
        print_rule

        if "$@"; then
            print_rule
            log_ok "$description"
            print_rule
            return 0
        else
            exit_status=$?
        fi

        if (( exit_status == 126 || exit_status == 127 )); then
            print_rule
            log_error "Permanent command failure while running: $description"
            log_error "Exit status: $exit_status"
            log_error "The command is unavailable or cannot be executed; retrying will not help."
            print_rule
            return "$exit_status"
        fi

        print_rule
        log_warning "Command failed: $description"
        log_warning "Exit status: $exit_status"
        log_warning "Retrying in $retry_delay seconds. Press Ctrl+C to stop."
        print_rule

        sleep "$retry_delay"
        attempt=$((attempt + 1))

        if (( retry_delay < max_retry_delay )); then
            retry_delay=$((retry_delay * 2))
            if (( retry_delay > max_retry_delay )); then
                retry_delay=$max_retry_delay
            fi
        fi
    done
}

# ------------------------------------------------------------------------------
# Kraken2 helpers
# ------------------------------------------------------------------------------

download_kraken2_library_until_success() {
    local database="$1"
    local library="$2"
    shift 2

    if ! retry_until_success \
        "Kraken2 library '$library' for database '$database'" \
        kraken2-build \
        "$@" \
        --download-library "$library" \
        --threads "$threads" \
        --db "$database"; then
        log_error "Kraken2 library download cannot continue: $library"
        exit 1
    fi
}

download_kraken2_taxonomy_until_success() {
    local database="$1"
    local downloader="$2"
    shift 2

    if ! retry_until_success \
        "NCBI taxonomy for Kraken2 database '$database'" \
        "$downloader" \
        --download-taxonomy \
        "$@" \
        --threads "$threads" \
        --db "$database"; then
        log_error "Kraken2 taxonomy download cannot continue for database: $database"
        exit 1
    fi
}

build_kraken2_database() {
    local database="$1"

    kraken2-build \
        --build \
        --threads "$threads" \
        --db "$database" \
        "${kraken_build_opts[@]}"
}

install_kraken_helper() {
    local source_file="$1"
    local target_name="$2"

    if [[ ! -d "$KRAKEN_ENV_LIBEXEC" ]]; then
        log_error "Kraken2 environment libexec directory not found:"
        log_error "  $KRAKEN_ENV_LIBEXEC"
        exit 1
    fi

    copy_required_file "$source_file" "$KRAKEN_ENV_LIBEXEC/$target_name"
    chmod +x "$KRAKEN_ENV_LIBEXEC/$target_name"

    # Preserve the original package-cache patch when that exact Kraken2 package
    # directory exists. The active environment copy above is the required one.
    if [[ -d "$KRAKEN_PKG_LIBEXEC" ]]; then
        cp -f "$source_file" "$KRAKEN_PKG_LIBEXEC/$target_name"
        chmod +x "$KRAKEN_PKG_LIBEXEC/$target_name"
    else
        log_warning "Kraken2 package-cache libexec not found; skipped optional copy:"
        log_warning "  $KRAKEN_PKG_LIBEXEC"
    fi
}

restore_default_rsync_helper() {
    install_kraken_helper "$dir/Installation/rsync_from_ncbi.pl" "rsync_from_ncbi.pl"
}

restore_default_genomic_library_helper() {
    install_kraken_helper \
        "$dir/Installation/download_genomic_library.sh" \
        "download_genomic_library.sh"
}

patch_perl_local_download_dir() {
    local perl_script="$1"
    local local_directory="$2"

    sed -i \
        's|^[[:space:]]*my \$local_download_dir = .*;|my $local_download_dir = "'"$local_directory"'";|' \
        "$perl_script"

    if ! grep -Fq 'my $local_download_dir = "'"$local_directory"'";' "$perl_script"; then
        log_error "Could not patch Perl local_download_dir in: $perl_script"
        exit 1
    fi
}
patch_shell_local_download_dir() {
    local shell_script="$1"
    local local_directory="$2"

    sed -i \
        "s|^[[:space:]]*local_download_dir=.*|    local_download_dir=\"$local_directory\"|" \
        "$shell_script"

    if ! grep -Fq "local_download_dir=\"$local_directory\"" "$shell_script"; then
        log_error "Could not patch shell local_download_dir in: $shell_script"
        exit 1
    fi
}

copy_manifest_with_offline_folder() {
    local source_script="$1"
    local destination_script="$2"

    copy_required_file "$source_script" "$destination_script"
    sed -i \
        "s|^offline_files_folder=.*|offline_files_folder=$offline_files_folder|" \
        "$destination_script"
    run_required_script "$destination_script"
}

prepare_local_kraken_host_genome() {
    local database="$1"
    local source_gz="$2"
    local compressed_name="$3"
    local fasta_name="$4"

    mkdir -p "$database"

    if [[ -s "$database/$fasta_name" ]]; then
        log_info "Using existing host FASTA: $database/$fasta_name"
        return 0
    fi

    copy_required_file "$source_gz" "$database/$compressed_name"
    unpigz -f "$database/$compressed_name"
    mv -f "$database/${compressed_name%.gz}" "$database/$fasta_name"
}

# ------------------------------------------------------------------------------
# Installation stages
# ------------------------------------------------------------------------------

initialize_installation() {
    cd "$dir" || exit 1
    printf '%s\n' "$condapath" > "$dir/condaPath"
    # shellcheck disable=SC1090
    source "$condapath/etc/profile.d/conda.sh"

    ORIGINAL_CHANNEL_PRIORITY="$(conda config --show channel_priority 2>/dev/null | awk '{print $2}')"
    : "${ORIGINAL_CHANNEL_PRIORITY:=flexible}"

    log_info "Original Conda channel priority: $ORIGINAL_CHANNEL_PRIORITY"
    run_required_command \
        "Setting Conda channel priority to flexible for legacy MTD environments" \
        conda config --set channel_priority flexible
}

install_system_dependencies() {
    log_info "Installing system dependencies..."

    # Commands are intentionally kept in the same order as the working script.
    sudo_with_pass "sudo apt-get update"
    sudo_with_pass "sudo apt-get install libgeos-dev -y"
    sudo_with_pass "sudo apt install libharfbuzz-dev libfribidi-dev libfreetype6-dev -y"
    sudo_with_pass "sudo apt install libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev -y"
    sudo_with_pass "sudo apt install libharfbuzz-dev rsync libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev pigz -y"
}

create_conda_environments() {
    safe_conda_deactivate

    run_required_command \
        "Creating MTD-fastp environment" \
        conda env create -f "$dir/Installation/MTD_fastp.yml"

    run_required_command \
        "Creating MTD environment" \
        conda env create -f "$dir/Installation/MTD.yml"

    run_required_command \
        "Updating MTD environment with R additions" \
        conda env update -n MTD -f "$dir/Installation/MTD_R_additions.yml"

    run_required_command \
        "Installing R packages in the MTD environment" \
        conda run -n MTD bash "$dir/update_fix/Install.R.packages.MTD.sh"

    run_required_command \
        "Checking R packages in the MTD environment" \
        conda run -n MTD bash "$dir/update_fix/check_R_pkg.MTD.sh"

    log_info "Creating Python 2 and HAllA environments..."
    sed -i 's/^rpy2[>=<]/# &/' "$dir/Installation/pip.requirements"

    run_required_command \
        "Creating py2 environment" \
        conda env create -f "$dir/Installation/py2.yml"

    run_required_command \
        "Creating halla0820 environment" \
        conda env create -f "$dir/Installation/halla0820.yml"

    safe_conda_deactivate

    run_required_command \
        "Creating initial R412 environment" \
        conda env create -f "$dir/Installation/R412.yml"

    sed -i '/^# *rpy2/s/^# *//' "$dir/Installation/pip.requirements"

    chmod +x "$dir/aux_scripts/ssGSEA/resolve_ssgsea_go_terms.py"
    run_required_command \
        "Checking ssGSEA GO resolver syntax" \
        python3 -m py_compile "$dir/aux_scripts/ssGSEA/resolve_ssgsea_go_terms.py"

    require_env_command MTD kraken2-build
    require_env_command MTD bracken-build
    require_env_command MTD humann
    require_env_command MTD hisat2-build
    require_env_command MTD makeblastdb
}

install_halla_dependencies() {
    activate_required_env halla0820

    conda install -n halla0820 -y -c conda-forge pkg-config
    conda install -n halla0820 -y -c conda-forge ca-certificates openssl libcurl curl
    conda install -n halla0820 -y -c conda-forge libuv

    R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/lattice/lattice_0.22-7.tar.gz', repos=NULL, type='source', Ncpus=$threads)"
    R -e "install.packages('$dir/update_fix/pvr_pkg/Matrix_1.6-5.tar.gz', repos=NULL, type='source', Ncpus=$threads)"
    R -e "install.packages('$dir/update_fix/pvr_pkg/mnormt_2.1.0.tar.gz', repos=NULL, type='source', Ncpus=$threads)"
    R -e "install.packages('$dir/update_fix/pvr_pkg/nlme_3.1-167.tar.gz', repos=NULL, type='source', Ncpus=$threads)"
    R -e "install.packages('$dir/update_fix/pvr_pkg/GPArotation_2024.3-1.tar.gz', repos=NULL, type='source', Ncpus=$threads)"
    R -e "install.packages('$dir/update_fix/pvr_pkg/psych_2.5.3.tar.gz', repos=NULL, type='source', Ncpus=$threads)"
    R -e "install.packages('$dir/update_fix/pvr_pkg/foreign_0.8-89.tar.gz', repos=NULL, type='source', Ncpus=$threads)"
    R -e "install.packages('$dir/update_fix/pvr_pkg/R.methodsS3_1.8.2.tar.gz', repos=NULL, type='source', Ncpus=$threads)"
    R -e "install.packages('$dir/update_fix/pvr_pkg/R.oo_1.27.0.tar.gz', repos=NULL, type='source', Ncpus=$threads)"
    R -e "install.packages('$dir/update_fix/pvr_pkg/rtf_0.4-14.tar.gz', repos=NULL, type='source', Ncpus=$threads)"
    R -e "install.packages('$dir/update_fix/pvr_pkg/psychTools_2.4.3.tar.gz', repos=NULL, type='source', Ncpus=$threads)"
    R -e "install.packages('https://cran.r-project.org/src/contrib/XICOR_0.4.1.tar.gz', repos=NULL, type='source', Ncpus=$threads)"
    R -e "install.packages('https://cran.r-project.org/src/contrib/mclust_6.1.2.tar.gz', repos=NULL, type='source', Ncpus=$threads)"
    R -e 'install.packages("BiocManager", repos = "https://cloud.r-project.org")'
    R -e "install.packages('$dir/update_fix/pvr_pkg/MASS_7.3-60.tar.gz', repos=NULL, type='source')"
    R -e "install.packages('$dir/update_fix/pvr_pkg/preprocessCore_1.72.0.tar.gz', repos=NULL, type='source')"
    R -e 'install.packages("remotes", repos="https://cloud.r-project.org")'
    R -e 'remotes::install_url("https://cran.r-project.org/src/contrib/EnvStats_3.1.0.tar.gz", dependencies=TRUE)'
    R -e 'remotes::install_version("Hmisc", version = "4.8-0", repos = "https://cloud.r-project.org")'
    R -e "install.packages('https://cran.r-project.org/src/contrib/00Archive/eva/eva_0.2.6.tar.gz', repos=NULL, type='source')"

    conda run -n halla0820 "$dir/update_fix/check_R_pkg.halla0820.sh"
    safe_conda_deactivate
}

install_mtd_extra_tools() {
    activate_required_env MTD

    # Preserved from the original installer; rsync is also installed in the
    # initial system-dependency stage.
    sudo_with_pass "sudo apt-get update"
    sudo_with_pass "sudo apt-get install rsync -y"

    safe_conda_deactivate
    activate_required_env MTD
    run_required_command \
        "Installing pkg-config in MTD" \
        conda install -n MTD -y -c conda-forge pkg-config
    run_required_command \
        "Installing NCBI Datasets CLI in MTD" \
        conda install -n MTD -y -c conda-forge ncbi-datasets-cli
    run_required_command \
        "Installing eggNOG-mapper and DIAMOND in MTD" \
        conda install -n MTD -y -c conda-forge -c bioconda eggnog-mapper diamond

    require_env_command MTD kraken2-build
}

prepare_virome_files() {
    copy_required_file \
        "$offline_files_folder/Ref_genomes/MTD_virus/virushostdb.genomic.fna.gz" \
        "$dir/virushostdb.genomic.fna.gz"

    unpigz -f "$dir/virushostdb.genomic.fna.gz"
    cat \
        "$dir/Installation/M33262_SIVMM239.fa" \
        "$dir/virushostdb.genomic.fna" \
        > "$dir/viruses4kraken.fa"

    log_info "Preparing additional NCBI RefSeq viral files..."
    copy_manifest_with_offline_folder \
        "$dir/manifest.virus.sh" \
        "$offline_files_folder/Kraken2DB_micro/library/manifest.virus.sh"

    sed -i -E \
        's/^>([^ ]+) (.+)/>\1 [\1] \2./' \
        "$offline_files_folder/Kraken2DB_micro/library/viral/all_viral_genomes.fna"

    # The original installer intentionally leaves this append operation disabled:
    # cat "$offline_files_folder/Kraken2DB_micro/library/viral/all_viral_genomes.fna" >> "$dir/viruses4kraken.fa"
}

install_default_kraken_helpers() {
    restore_default_rsync_helper
    restore_default_genomic_library_helper
}

prepare_microbiome_manifests() {
    copy_manifest_with_offline_folder \
        "$dir/manifest.bacteria.sh" \
        "$offline_files_folder/Kraken2DB_micro/library/manifest.bacteria.sh"

    copy_required_file \
        "$dir/manifest.sh" \
        "$offline_files_folder/Kraken2DB_micro/library/manifest.sh"

    sed -i \
        "s|^offline_files_folder=.*|offline_files_folder=$offline_files_folder|" \
        "$offline_files_folder/Kraken2DB_micro/library/manifest.sh"
}

add_local_archaea_library() {
    local database="$1"
    local helper="$KRAKEN_ENV_LIBEXEC/rsync_from_ncbi.pl"

    copy_manifest_with_offline_folder \
        "$dir/manifest.archea.sh" \
        "$offline_files_folder/Kraken2DB_micro/library/manifest.archea.sh"

    install_kraken_helper \
        "$dir/Installation/rsync_from_ncbi_archaea.pl" \
        "rsync_from_ncbi.pl"

    patch_perl_local_download_dir \
        "$helper" \
        "$offline_files_folder/Kraken2DB_micro/library/archaea/all/"

    copy_required_file \
        "$offline_files_folder/Kraken2DB_micro/library/archaea/assembly_summary_archaea.txt" \
        "$offline_files_folder/Kraken2DB_micro/library/archaea/assembly_summary.txt"

    chmod +x "$helper"
    log_info "Adding local archaeal sequences to Kraken2 database..."
    download_kraken2_library_until_success "$database" "archaea" --use-ftp

    restore_default_rsync_helper
}

add_local_bacteria_library() {
    local database="$1"
    local helper="$KRAKEN_ENV_LIBEXEC/rsync_from_ncbi.pl"

    install_kraken_helper \
        "$dir/Installation/rsync_from_ncbi_bacteria.pl" \
        "rsync_from_ncbi.pl"

    patch_perl_local_download_dir \
        "$helper" \
        "$offline_files_folder/Kraken2DB_micro/library/bacteria/all/"

    copy_required_file \
        "$offline_files_folder/Kraken2DB_micro/library/bacteria/assembly_summary_bacteria.txt" \
        "$offline_files_folder/Kraken2DB_micro/library/bacteria/assembly_summary.txt"

    chmod +x "$helper"
    log_info "Adding local bacterial sequences to Kraken2 database..."
    download_kraken2_library_until_success "$database" "bacteria" --use-ftp

    restore_default_rsync_helper
}

add_local_plasmid_library() {
    local database="$1"
    local manifest_destination="$offline_files_folder/Kraken2DB_micro/library/manifest.plasmid.sh"
    local custom_helper="$dir/Installation/download_genomic_library_plasmid.sh"

    copy_required_file "$dir/manifest.plasmid.sh" "$manifest_destination"
    sed -i \
        "s|^LOCAL_DIR=.*|LOCAL_DIR=$offline_files_folder/Kraken2DB_micro/library/plasmid/|" \
        "$manifest_destination"
    run_required_script "$manifest_destination"

    patch_shell_local_download_dir \
        "$custom_helper" \
        "$offline_files_folder/Kraken2DB_micro/library/plasmid/"

    install_kraken_helper "$custom_helper" "download_genomic_library.sh"

    log_info "Adding local plasmid sequences to Kraken2 database..."
    kraken2-build \
        --download-library plasmid \
        --threads "$threads" \
        --db "$database"

    restore_default_genomic_library_helper
}

build_microbiome_kraken_database() {
    local database="$dir/kraken2DB_micro"

    prepare_microbiome_manifests

    chmod +x "$dir/kraken2-build-download-taxonomy"
    log_info "Downloading NCBI taxonomy database with Kraken2..."
    download_kraken2_taxonomy_until_success \
        "$database" \
        "$dir/kraken2-build-download-taxonomy"

    add_local_archaea_library "$database"
    add_local_bacteria_library "$database"

    log_info "Downloading RefSeq Protozoa library..."
    download_kraken2_library_until_success "$database" "protozoa" --use-ftp

    log_info "Downloading RefSeq Fungi library..."
    download_kraken2_library_until_success "$database" "fungi" --use-ftp

    add_local_plasmid_library "$database"

    log_info "Downloading UniVec_Core library..."
    download_kraken2_library_until_success "$database" "UniVec_Core" --use-ftp

    log_info "Adding custom viral sequences to Kraken2 database..."
    kraken2-build \
        --add-to-library "$dir/viruses4kraken.fa" \
        --threads "$threads" \
        --db "$database"

    log_info "Building final Kraken2 microbiome database..."
    build_kraken2_database "$database"
}

build_human_kraken_database() {
    local database="$dir/kraken2DB_human"

    download_kraken2_taxonomy_until_success \
        "$database" \
        kraken2-build \
        --use-ftp

    download_kraken2_library_until_success \
        "$database" \
        "human" \
        --use-ftp

    build_kraken2_database "$database"
}

build_mouse_kraken_database() {
    local database="$dir/kraken2DB_mice"
    local fasta_name="GCF_000001635.27_GRCm39_genomic.fa"

    prepare_local_kraken_host_genome \
        "$database" \
        "$offline_files_folder/GCF_000001635.27_GRCm39_genomic.fna.gz" \
        "GCF_000001635.27_GRCm39_genomic.fna.gz" \
        "$fasta_name"

    download_kraken2_taxonomy_until_success \
        "$database" \
        kraken2-build \
        --use-ftp

    kraken2-build \
        --add-to-library "$database/$fasta_name" \
        --threads "$threads" \
        --db "$database"

    build_kraken2_database "$database"
}

build_rhesus_kraken_database() {
    local database="$dir/kraken2DB_rhesus"
    local fasta_name="GCF_003339765.1_Mmul_10_genomic.fa"

    prepare_local_kraken_host_genome \
        "$database" \
        "$offline_files_folder/GCF_003339765.1_Mmul_10_genomic.fna.gz" \
        "GCF_003339765.1_Mmul_10_genomic.fna.gz" \
        "$fasta_name"

    download_kraken2_taxonomy_until_success \
        "$database" \
        kraken2-build \
        --use-ftp

    kraken2-build \
        --add-to-library "$database/$fasta_name" \
        --threads "$threads" \
        --db "$database"

    build_kraken2_database "$database"
}

build_bracken_database() {
    local database="$dir/kraken2DB_micro"

    if [[ -z "$kmer" ]]; then
        bracken-build \
            -d "$database" \
            -t "$threads" \
            -l "$read_len"
    else
        bracken-build \
            -d "$database" \
            -t "$threads" \
            -l "$read_len" \
            -k "$kmer"
    fi
}

install_humann_databases() {
    local humann_dir="$dir/HUMAnN/ref_database"

    mkdir -p "$humann_dir/chocophlan"
    mkdir -p "$humann_dir/uniref"
    mkdir -p "$humann_dir/utility_mapping"

    copy_required_file \
        "$offline_files_folder/HUMAnN/full_chocophlan.v201901_v31.tar.gz" \
        "$humann_dir/full_chocophlan.v201901_v31.tar.gz"

    copy_required_file \
        "$offline_files_folder/HUMAnN/uniref90_annotated_v201901b_full.tar.gz" \
        "$humann_dir/uniref90_annotated_v201901b_full.tar.gz"

    copy_required_file \
        "$offline_files_folder/HUMAnN/full_mapping_v201901b.tar.gz" \
        "$humann_dir/full_mapping_v201901b.tar.gz"

    tar xzvf \
        "$humann_dir/full_chocophlan.v201901_v31.tar.gz" \
        -C "$humann_dir/chocophlan/"

    tar xzvf \
        "$humann_dir/uniref90_annotated_v201901b_full.tar.gz" \
        -C "$humann_dir/uniref/"

    tar xzvf \
        "$humann_dir/full_mapping_v201901b.tar.gz" \
        -C "$humann_dir/utility_mapping/"

    humann_config --update database_folders nucleotide "$humann_dir/chocophlan"
    humann_config --update database_folders protein "$humann_dir/uniref"
    humann_config --update database_folders utility_mapping "$humann_dir/utility_mapping"
}

copy_host_reference_gtfs() {
    mkdir -p "$dir/ref_rhesus"
    mkdir -p "$dir/ref_human"
    mkdir -p "$dir/ref_mouse"

    copy_required_file \
        "$offline_files_folder/Ref_genomes/Macaca_mulatta/Macaca_mulatta.Mmul_10.104.gtf.gz" \
        "$dir/ref_rhesus/Macaca_mulatta.Mmul_10.104.gtf.gz"

    copy_required_file \
        "$offline_files_folder/Ref_genomes/Homo_sapiens/Homo_sapiens.GRCh38.104.gtf.gz" \
        "$dir/ref_human/Homo_sapiens.GRCh38.104.gtf.gz"

    copy_required_file \
        "$offline_files_folder/Ref_genomes/Mus_musculus/Mus_musculus.GRCm39.104.gtf.gz" \
        "$dir/ref_mouse/Mus_musculus.GRCm39.104.gtf.gz"
}

build_hisat2_host_index() {
    local index_directory="$1"
    local gtf_gz="$2"
    local genome_gz="$3"

    mkdir -p "$index_directory"

    gzip -dc "$gtf_gz" > "$index_directory/genome.gtf"
    python "$dir/Installation/hisat2_extract_splice_sites.py" \
        "$index_directory/genome.gtf" \
        > "$index_directory/genome.ss"
    python "$dir/Installation/hisat2_extract_exons.py" \
        "$index_directory/genome.gtf" \
        > "$index_directory/genome.exon"

    gzip -dc "$genome_gz" > "$index_directory/genome.fa"

    hisat2-build \
        --large-index \
        -p "$threads" \
        --exon "$index_directory/genome.exon" \
        --ss "$index_directory/genome.ss" \
        "$index_directory/genome.fa" \
        "$index_directory/genome_tran"
}

build_magic_blast_databases() {
    mkdir -p "$dir/human_blastdb"
    mkdir -p "$dir/mouse_blastdb"
    mkdir -p "$dir/rhesus_blastdb"

    makeblastdb \
        -in "$dir/hisat2_index_human/genome.fa" \
        -dbtype nucl \
        -parse_seqids \
        -out "$dir/human_blastdb/human_blastdb"

    makeblastdb \
        -in "$dir/hisat2_index_mouse/genome.fa" \
        -dbtype nucl \
        -parse_seqids \
        -out "$dir/mouse_blastdb/mouse_blastdb"

    makeblastdb \
        -in "$dir/hisat2_index_rhesus/genome.fa" \
        -dbtype nucl \
        -parse_seqids \
        -out "$dir/rhesus_blastdb/rhesus_blastdb"
}

install_r412_and_annotation_packages() {
    safe_conda_deactivate

    conda install -n py2 -y -c conda-forge pkg-config
    activate_required_env R412
    run_required_command \
        "Installing pkg-config in R412" \
        conda install -n R412 -y -c conda-forge pkg-config

    # Preserved libcurl troubleshooting step from the original installer.
    wget \
        -T 300 \
        -t 5 \
        -N \
        --no-if-modified-since \
        https://cran.r-project.org/src/contrib/Archive/curl/curl_4.3.2.tar.gz

    if [[ -f /usr/lib/x86_64-linux-gnu/pkgconfig/libcurl.pc ]]; then
        locate_lib=/usr/lib/x86_64-linux-gnu/pkgconfig
    else
        locate_lib="$(dirname "$(locate libcurl 2>/dev/null | grep '\.pc' | head -n 1)")"
    fi

    # `locate_lib` is intentionally retained for compatibility with the old
    # troubleshooting block, although it is not consumed later in this script.
    : "${locate_lib:=}"

    cd "$dir" || exit 1
    safe_conda_deactivate

    conda env remove -n R412 -y
    rm -rf "$condapath/envs/R412"
    rm -rf "$condapath/pkgs/r-base-4.1.2-hde4fec0_0"
    rm -f "$condapath"/pkgs/r-base-4.1.2-hde4fec0_0*.tar.bz2
    rm -f "$condapath"/pkgs/r-base-4.1.2-hde4fec0_0*.conda

    conda clean --packages --tarballs -y
    run_required_command \
        "Setting Conda channel priority to strict for R412 recreation" \
        conda config --set channel_priority strict
    run_required_command \
        "Recreating R412 environment" \
        conda env create -f "$dir/Installation/R412.yml"
    activate_required_env R412

    run_required_script "$dir/update_fix/Install.R.packages.R412_optimized.sh"

    Rscript -e \
        'install.packages("https://bioconductor.org/packages/3.19/bioc/src/contrib/UCSC.utils_1.0.0.tar.gz", repos=NULL, type="source", dependencies=FALSE); library(UCSC.utils); packageVersion("UCSC.utils")'

    run_required_script "$dir/update_fix/check_R_pkg.R412.sh"
    safe_conda_deactivate

    # Install annotation tools in the base environment, as in the original.
    R -e 'BiocManager::install("GenomeInfoDb")'
    bash "$dir/update_fix/Install.R.AnnotPackages.base.sh"
}

show_r_package_versions() {
    echo "${g}"
    echo "*********************************"
    echo "R packages version for conda envs"
    echo "*********************************"
    echo "${w}"

    conda run -n MTD "$dir/update_fix/check_R_pkg.MTD.sh"
    conda run -n R412 "$dir/update_fix/check_R_pkg.R412.sh"
    conda run -n halla0820 "$dir/update_fix/check_R_pkg.halla0820.sh"
}

# ------------------------------------------------------------------------------
# Main installation sequence
# ------------------------------------------------------------------------------

main() {
    init_colors
    parse_arguments "$@"
    validate_arguments
    configure_paths_and_options
    initialize_installation

    install_system_dependencies
    create_conda_environments

    show_progress ">>" "10%" "Installing HAllA dependencies..."
    install_halla_dependencies

    show_progress ">>>" "15%" "Preparing virome database and MTD tools..."
    install_mtd_extra_tools
    prepare_virome_files
    install_default_kraken_helpers

    show_progress ">>>>" "20%" \
        "Preparing microbiome (virus, bacteria, archaea, protozoa, fungi, plasmid, UniVec_Core) database..."
    build_microbiome_kraken_database

    show_progress ">>>>>>" "30%" "Preparing host (human) Kraken2 database..."
    build_human_kraken_database

    show_progress ">>>>>>>" "35%" "Preparing host (mouse) Kraken2 database..."
    build_mouse_kraken_database

    show_progress ">>>>>>>>" "40%" "Preparing host (rhesus monkey) Kraken2 database..."
    build_rhesus_kraken_database

    show_progress ">>>>>>>>>" "45%" "Building Bracken database..."
    build_bracken_database

    show_progress ">>>>>>>>>>>" "55%" "Installing HUMAnN3 databases..."
    install_humann_databases

    show_progress ">>>>>>>>>>>>>>" "70%" \
        "Fetching host (rhesus, human, mouse) references from local storage: $offline_files_folder"
    copy_host_reference_gtfs

    show_progress ">>>>>>>>>>>>>>>" "75%" \
        "Building host indexes (rhesus monkey) for HISAT2..."
    build_hisat2_host_index \
        "$dir/hisat2_index_rhesus" \
        "$dir/ref_rhesus/Macaca_mulatta.Mmul_10.104.gtf.gz" \
        "$offline_files_folder/Ref_genomes/Macaca_mulatta/Macaca_mulatta.Mmul_10.dna.toplevel.fa.gz"

    show_progress ">>>>>>>>>>>>>>>>" "80%" \
        "Building host indexes (mouse) for HISAT2..."
    build_hisat2_host_index \
        "$dir/hisat2_index_mouse" \
        "$dir/ref_mouse/Mus_musculus.GRCm39.104.gtf.gz" \
        "$offline_files_folder/Ref_genomes/Mus_musculus/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"

    show_progress ">>>>>>>>>>>>>>>>>" "85%" \
        "Building host indexes (human) for HISAT2..."
    build_hisat2_host_index \
        "$dir/hisat2_index_human" \
        "$dir/ref_human/Homo_sapiens.GRCh38.104.gtf.gz" \
        "$offline_files_folder/Ref_genomes/Homo_sapiens/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"

    log_info "Creating BLAST databases for Magic-BLAST..."
    build_magic_blast_databases

    show_progress ">>>>>>>>>>>>>>>>>>" "90%" "Installing R packages..."
    install_r412_and_annotation_packages
    show_r_package_versions

    echo "${g}"
    echo "*********************************"
    echo
    echo "MTD installation progress:"
    echo ">>>>>>>>>>>>>>>>>>>>[100%]"
    echo "MTD installation is finished"
    echo "${w}"
}

main "$@"
