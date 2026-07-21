#!/usr/bin/env bash

# Automatically accept the Anaconda Terms of Service required by
# the default Conda channels during unattended MTD installation.
export CONDA_PLUGINS_AUTO_ACCEPT_TOS=yes

# ==============================================================================
# MTD installer
# ==============================================================================
# This version reorganizes the original installer into functions while preserving
# the execution order and the local Kraken2 helper-script workflow.
#
# Deliberately not using `set -e` here. The original installer continues after
# some non-critical commands fail, and enabling it would change existing behavior.
# ==============================================================================
#
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

dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
MANIFEST_SCRIPTS_DIR="$dir/aux_scripts/manifest_scripts"
KRAKEN_AUX_DIR="$dir/aux_scripts/Kraken2"

KRAKEN_ENV_LIBEXEC=""
KRAKEN_PKG_LIBEXEC=""
kraken_build_opts=()
ORIGINAL_CHANNEL_PRIORITY=""
KRAKEN_TAXONOMY_CACHE=""
# Validated and immutable Virus-Host DB mirror used by MTD Explorer.
VIRUSHOST_MIRROR_REPOSITORY="patrick-douglas/MTD"
VIRUSHOST_MIRROR_TAG="virushostdb-mirror-r235-g271.0-a250b2e61d9f"

VIRUSHOST_MIRROR_BASE_URL="https://github.com/${VIRUSHOST_MIRROR_REPOSITORY}/releases/download/${VIRUSHOST_MIRROR_TAG}"

# SHA-256 of the SHA256SUMS asset from the validated mirror release.
VIRUSHOST_MIRROR_SHA256SUMS_SHA256="a250b2e61d9f9365773205d04d019e0976a778ebd589553f3a0a0e6f159f4bec"


# Pinned MetaPhlAn database used by HUMAnN 3.9 / MetaPhlAn 4.1.1.
HUMANN_METAPHLAN_INDEX="mpa_vJun23_CHOCOPhlAnSGB_202403"
HUMANN_METAPHLAN_CACHE_DIRNAME="metaphlan_vJun23_202403_archives"
HUMANN_METAPHLAN_BASE_URL="https://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases"
HUMANN_METAPHLAN_BT2_BASE_URL="${HUMANN_METAPHLAN_BASE_URL}/bowtie2_indexes"
HUMANN_METAPHLAN_MAIN_ARCHIVE_MD5="d985de75a217cd319e721863f68e7d33"
HUMANN_METAPHLAN_BT2_ARCHIVE_MD5="8caae86b4d2931416cbdbb92f5985cef"

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
        c="$(tput setaf 6)"
    else
        w=""
        r=""
        g=""
        y=""
        p=""
        c=""
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
    local percent="$1"
    local message="$2"
    local width=30
    local filled
    local empty
    local filled_bar
    local empty_bar

    if ! [[ "$percent" =~ ^[0-9]+$ ]] ||
       (( percent < 0 || percent > 100 )); then
        log_error "Invalid installation progress percentage: $percent"
        return 1
    fi

    filled=$((percent * width / 100))
    empty=$((width - filled))

    printf -v filled_bar '%*s' "$filled" ''
    printf -v empty_bar '%*s' "$empty" ''

    filled_bar="${filled_bar// /#}"
    empty_bar="${empty_bar// /-}"

    echo
    print_rule
    printf '%sMTD Explorer installation%s\n' "$p" "$w"
    printf '[%s%s%s%s] %s%3d%%%s\n' \
        "$g" "$filled_bar" \
        "$w" "$empty_bar" \
        "$c" "$percent" "$w"
    printf '%s%s%s\n' "$g" "$message" "$w"
    print_rule
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

  $0 -o <offline_files_folder> [options]

Required:

  -o PATH  Persistent installation cache; created and populated if absent

Optional:

  -p PATH  Miniconda installation path
           Default: \$HOME/miniconda3

  -k INT   Kraken2 k-mer length used during --build
  -m INT   Kraken2 minimizer length used during --build
  -s INT   Kraken2 minimizer spaces used during --build
  -r INT   Bracken read length (default: 75)
  -h       Show this help message

The installer automatically downloads and installs Miniconda.

WARNING:
If the selected Miniconda directory already exists, the installer will ask
for explicit confirmation before permanently deleting it.

USAGE
}

parse_arguments() {
    if [[ $# -eq 1 && "$1" == "-h" ]]; then
        usage
        exit 0
    fi

   if [[ $# -lt 2 ]]; then
        usage
        exit 1
    fi

    while getopts ":p:o:k:m:s:r:h" option; do
        case "$option" in
            p) condapath="$OPTARG" ;;
            o) offline_files_folder="$OPTARG" ;;
            k) kmer="$OPTARG" ;;
            m) min_l="$OPTARG" ;;
            s) min_s="$OPTARG" ;;
            r) read_len="$OPTARG" ;;
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
    if [[ -z "$offline_files_folder" ]]; then
        log_error "The persistent installation cache path is required."
        usage
        exit 1
    fi

    # Expand a literal leading ~, including when the user passed it quoted.
    if [[ "$condapath" == "~" ]]; then
        condapath="$HOME"
    elif [[ "$condapath" == "~/"* ]]; then
        condapath="$HOME/${condapath#~/}"
    fi

    if [[ "$offline_files_folder" == "~" ]]; then
        offline_files_folder="$HOME"
    elif [[ "$offline_files_folder" == "~/"* ]]; then
        offline_files_folder="$HOME/${offline_files_folder#~/}"
    fi

    # -m allows the final path to be resolved even before Miniconda exists.
    if ! condapath="$(readlink -m -- "$condapath")"; then
        log_error "Could not resolve the Miniconda installation path."
        exit 1
    fi

    if [[ -z "$condapath" || "$condapath" == "/" || "$condapath" == "$HOME" ]]; then
        log_error "Unsafe Miniconda installation path:"
        log_error "  $condapath"
        log_error "Refusing to continue because this path must never be removed."
        exit 1
    fi

    if [[ -e "$offline_files_folder" && ! -d "$offline_files_folder" ]]; then
        log_error "The cache path exists but is not a directory:"
        log_error "  $offline_files_folder"
        exit 1
    fi

    if ! mkdir -p "$offline_files_folder"; then
        log_error "Could not create the installation cache directory:"
        log_error "  $offline_files_folder"
        exit 1
    fi

    if ! offline_files_folder="$(readlink -f -- "$offline_files_folder")"; then
        log_error "Could not resolve the installation cache path."
        exit 1
    fi

    if [[ ! -w "$offline_files_folder" ]]; then
        log_error "Installation cache directory is not writable:"
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

    log_info "Miniconda installation path:"
    log_info "  $condapath"

    log_info "Persistent installation cache:"
    log_info "  $offline_files_folder"
}

configure_paths_and_options() {
    KRAKEN_ENV_LIBEXEC="$condapath/envs/MTD/libexec"
    KRAKEN_PKG_LIBEXEC="$condapath/pkgs/kraken2-2.1.2-pl5262h7d875b9_0/libexec"
    KRAKEN_TAXONOMY_CACHE="$offline_files_folder/Kraken2_taxonomy_cache"

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

ensure_sudo_credentials() {
    if (( EUID == 0 )); then
        return 0
    fi

    if ! command -v sudo >/dev/null 2>&1; then
        log_error "sudo is required to install system dependencies."
        return 127
    fi

    # Reuse an existing valid sudo authentication, when available.
    if sudo -n -v >/dev/null 2>&1; then
        return 0
    fi

    if [[ ! -t 0 ]]; then
        log_error "Interactive sudo authentication is required."
        log_error "Run the installer directly from an interactive terminal."
        return 1
    fi

    log_info "Administrator privileges are required to install system dependencies."

    if ! sudo -v; then
        log_error "Could not obtain administrator privileges."
        return 1
    fi
}

run_as_root() {
    if (( EUID == 0 )); then
        "$@"
        return $?
    fi

    ensure_sudo_credentials || return $?

    # Authentication was validated above, so this command must not prompt.
    sudo -n "$@"
}
confirm_miniconda_removal() {
    local confirmation=""

    print_rule
    log_warning "An existing Miniconda installation was found at:"
    log_warning "  $condapath"
    echo
    log_warning "Continuing will permanently delete this directory."
    log_warning "All Conda environments and packages stored inside it will be lost."
    log_warning "This operation cannot be undone."
    print_rule

    if [[ ! -t 0 ]]; then
        log_error "Interactive confirmation is required before deleting Miniconda."
        log_error "Run the installer directly from an interactive terminal."
        exit 1
    fi

    printf '%s' "${y}Remove the existing Miniconda installation? [y/N]: ${w}"
    read -r confirmation

    case "$confirmation" in
        y | Y)
            echo
            log_warning "Removal confirmed by the user."
            ;;
        *)
            echo
            log_warning "Miniconda removal was not confirmed."
            log_warning "Installation cancelled without deleting anything."
            exit 0
            ;;
    esac
}

detect_miniconda_architecture() {
    local machine_arch

    machine_arch="$(uname -m)"

    case "$machine_arch" in
        x86_64 | amd64)
            printf '%s\n' "x86_64"
            ;;

        aarch64 | arm64)
            printf '%s\n' "aarch64"
            ;;

        *)
            log_error "Unsupported CPU architecture for Miniconda:"
            log_error "  $machine_arch"
            return 1
            ;;
    esac
}

remove_existing_miniconda() {
    if [[ ! -e "$condapath" ]]; then
        return 0
    fi

    case "$condapath" in
        "" | "/" | "$HOME")
            log_error "Refusing to remove unsafe path:"
            log_error "  $condapath"
            exit 1
            ;;
    esac

    confirm_miniconda_removal

    log_info "Removing the previous Miniconda installation:"
    log_info "  $condapath"

    if ! rm -rf -- "$condapath"; then
        log_error "Could not remove the previous Miniconda installation:"
        log_error "  $condapath"
        exit 1
    fi

    if [[ -e "$condapath" ]]; then
        log_error "The previous Miniconda directory still exists:"
        log_error "  $condapath"
        exit 1
    fi

    log_ok "Previous Miniconda installation removed."
}

install_miniconda() {
    local miniconda_arch
    local miniconda_url
    local miniconda_installer

    miniconda_arch="$(detect_miniconda_architecture)" || exit 1

    miniconda_url="https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-${miniconda_arch}.sh"
    miniconda_installer="${TMPDIR:-/tmp}/mtd-miniconda-installer-${USER:-user}-$$.sh"

    remove_existing_miniconda

    log_info "Downloading the latest Miniconda installer..."
    log_info "Architecture: $miniconda_arch"
    log_info "URL: $miniconda_url"

    rm -f -- "$miniconda_installer"

    if ! curl \
        --fail \
        --location \
        --show-error \
        --retry 5 \
        --retry-delay 10 \
        --connect-timeout 30 \
        --output "$miniconda_installer" \
        "$miniconda_url"; then

        rm -f -- "$miniconda_installer"

        log_error "Could not download the Miniconda installer."
        exit 1
    fi

    if [[ ! -s "$miniconda_installer" ]]; then
        rm -f -- "$miniconda_installer"

        log_error "The downloaded Miniconda installer is empty."
        exit 1
    fi

    log_ok "Miniconda installer downloaded successfully."

    log_info "Installing Miniconda in batch mode:"
    log_info "  $condapath"

    if ! bash "$miniconda_installer" -b -p "$condapath"; then
        rm -f -- "$miniconda_installer"

        log_error "Miniconda installation failed."
        exit 1
    fi

    rm -f -- "$miniconda_installer"

    if [[ ! -x "$condapath/bin/conda" ]]; then
        log_error "The Conda executable was not created:"
        log_error "  $condapath/bin/conda"
        exit 1
    fi

    if [[ ! -f "$condapath/etc/profile.d/conda.sh" ]]; then
        log_error "The Conda initialization script was not created:"
        log_error "  $condapath/etc/profile.d/conda.sh"
        exit 1
    fi

    # Remove cached shell command paths inherited from an older installation.
    hash -r

    log_ok "Miniconda installed and validated successfully."
    log_info "Installed Conda version:"
    "$condapath/bin/conda" --version
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
# Persistent installation cache
# ------------------------------------------------------------------------------

validate_downloaded_cache_file() {
    local downloaded_file="$1"
    local expected_name="$2"
    local newick_tail=""

    if [[ ! -s "$downloaded_file" ]]; then
        log_warning "Downloaded cache file is empty: $downloaded_file"
        return 1
    fi

    case "$expected_name" in
        *.tar.gz)
            tar -tzf "$downloaded_file" >/dev/null 2>&1
            ;;
        *.tar)
            tar -tf "$downloaded_file" >/dev/null 2>&1
            ;;
        *.gz)
            gzip -t "$downloaded_file" >/dev/null 2>&1
            ;;
        *.bz2)
            bzip2 -t "$downloaded_file" >/dev/null 2>&1
            ;;
        *.md5)
            grep -Eq \
                '^[[:xdigit:]]{32}[[:space:]]+\*?[^[:space:]]+$' \
                "$downloaded_file"
            ;;
        *.nwk)
            newick_tail="$(
                tail -c 1024 "$downloaded_file" 2>/dev/null |
                    tr -d '\r\n[:space:]'
            )"
            [[ -n "$newick_tail" && "${newick_tail: -1}" == ";" ]]
            ;;
        *)
            return 0
            ;;
    esac
}

download_cache_file_once() {
    local description="$1"
    local url="$2"
    local destination="$3"
    local partial_file="${destination}.part"
    local exit_status=0

    if ! mkdir -p "$(dirname "$destination")"; then
        log_error "Could not create cache destination directory:"
        log_error "  $(dirname "$destination")"
        return 1
    fi

    log_info "$description"
    log_info "URL:         $url"
    log_info "Destination: $destination"

    if command -v curl >/dev/null 2>&1; then
        curl \
            --fail \
            --location \
            --connect-timeout 30 \
            --retry 5 \
            --retry-delay 10 \
            --continue-at - \
            --output "$partial_file" \
            "$url"
        exit_status=$?
    elif command -v wget >/dev/null 2>&1; then
        wget \
            --continue \
            --tries=5 \
            --timeout=60 \
            --read-timeout=60 \
            --output-document="$partial_file" \
            "$url"
        exit_status=$?
    else
        log_error "Neither curl nor wget is available for downloading cache files."
        return 127
    fi

    if (( exit_status != 0 )); then
        return "$exit_status"
    fi

    if ! validate_downloaded_cache_file "$partial_file" "$destination"; then
        log_warning "Downloaded file failed archive validation:"
        log_warning "  $destination"
        rm -f "$partial_file"
        return 1
    fi

    if ! mv -f "$partial_file" "$destination"; then
        log_error "Could not finalize cached file:"
        log_error "  $destination"
        return 1
    fi

    log_ok "Cached successfully: $destination"
}

ensure_cached_file() {
    local description="$1"
    local url="$2"
    local destination="$3"

    if [[ -s "$destination" ]]; then
        if validate_downloaded_cache_file "$destination" "$destination"; then
            log_ok "Using validated cached file: $destination"
            return 0
        fi

        log_warning "Cached file failed validation and will be downloaded again:"
        log_warning "  $destination"

        rm -f -- \
            "$destination" \
            "${destination}.part"

    elif [[ -e "$destination" ]]; then
        log_warning "Removing empty cache file: $destination"
        rm -f -- "$destination"
    fi

    if ! retry_until_success \
        "$description" \
        download_cache_file_once \
        "$description" \
        "$url" \
        "$destination"
    then
        log_error "Could not prepare required cache file:"
        log_error "  $destination"
        exit 1
    fi
}

download_humann_metaphlan_cache_assets() {
    local cache_dir="$offline_files_folder/HUMAnN/$HUMANN_METAPHLAN_CACHE_DIRNAME"
    local index="$HUMANN_METAPHLAN_INDEX"

    ensure_cached_file \
        "MetaPhlAn vJun23 database archive" \
        "$HUMANN_METAPHLAN_BASE_URL/${index}.tar" \
        "$cache_dir/${index}.tar"

    ensure_cached_file \
        "MetaPhlAn vJun23 database MD5 manifest" \
        "$HUMANN_METAPHLAN_BASE_URL/${index}.md5" \
        "$cache_dir/${index}.md5"

    ensure_cached_file \
        "MetaPhlAn vJun23 Bowtie2 indexes" \
        "$HUMANN_METAPHLAN_BT2_BASE_URL/${index}_bt2.tar" \
        "$cache_dir/${index}_bt2.tar"

    ensure_cached_file \
        "MetaPhlAn vJun23 Bowtie2 MD5 manifest" \
        "$HUMANN_METAPHLAN_BT2_BASE_URL/${index}_bt2.md5" \
        "$cache_dir/${index}_bt2.md5"

    ensure_cached_file \
        "MetaPhlAn vJun23 taxonomy tree" \
        "$HUMANN_METAPHLAN_BASE_URL/${index}.nwk" \
        "$cache_dir/${index}.nwk"

    ensure_cached_file \
        "MetaPhlAn vJun23 marker information" \
        "$HUMANN_METAPHLAN_BASE_URL/${index}_marker_info.txt.bz2" \
        "$cache_dir/${index}_marker_info.txt.bz2"

    ensure_cached_file \
        "MetaPhlAn vJun23 species information" \
        "$HUMANN_METAPHLAN_BASE_URL/${index}_species.txt.bz2" \
        "$cache_dir/${index}_species.txt.bz2"
}

validate_humann_metaphlan_cache() {
    local cache_dir="$offline_files_folder/HUMAnN/$HUMANN_METAPHLAN_CACHE_DIRNAME"
    local index="$HUMANN_METAPHLAN_INDEX"
    local main_archive="$cache_dir/${index}.tar"
    local bt2_archive="$cache_dir/${index}_bt2.tar"
    local observed_md5=""
    local required_file=""
    local main_listing=""
    local bt2_listing=""
    local newick_tail=""

    local -a required_files=(
        "${index}.tar"
        "${index}.md5"
        "${index}_bt2.tar"
        "${index}_bt2.md5"
        "${index}.nwk"
        "${index}_marker_info.txt.bz2"
        "${index}_species.txt.bz2"
    )

    local -a required_bt2_files=(
        "${index}.1.bt2l"
        "${index}.2.bt2l"
        "${index}.3.bt2l"
        "${index}.4.bt2l"
        "${index}.rev.1.bt2l"
        "${index}.rev.2.bt2l"
    )

    for required_file in "${required_files[@]}"; do
        if [[ ! -s "$cache_dir/$required_file" ]]; then
            log_warning "Missing or empty MetaPhlAn cache file:"
            log_warning "  $cache_dir/$required_file"
            return 1
        fi
    done

    observed_md5="$(md5sum "$main_archive" | awk '{print $1}')"

    if [[ "$observed_md5" != "$HUMANN_METAPHLAN_MAIN_ARCHIVE_MD5" ]]; then
        log_warning "MetaPhlAn database archive MD5 mismatch."
        log_warning "Expected: $HUMANN_METAPHLAN_MAIN_ARCHIVE_MD5"
        log_warning "Observed: $observed_md5"
        return 1
    fi

    observed_md5="$(md5sum "$bt2_archive" | awk '{print $1}')"

    if [[ "$observed_md5" != "$HUMANN_METAPHLAN_BT2_ARCHIVE_MD5" ]]; then
        log_warning "MetaPhlAn Bowtie2 archive MD5 mismatch."
        log_warning "Expected: $HUMANN_METAPHLAN_BT2_ARCHIVE_MD5"
        log_warning "Observed: $observed_md5"
        return 1
    fi

    if ! (
        cd "$cache_dir" &&
        md5sum -c "${index}.md5" >/dev/null 2>&1 &&
        md5sum -c "${index}_bt2.md5" >/dev/null 2>&1
    ); then
        log_warning "MetaPhlAn official MD5 manifest validation failed."
        return 1
    fi

    if ! bzip2 -t \
        "$cache_dir/${index}_marker_info.txt.bz2" \
        >/dev/null 2>&1
    then
        log_warning "MetaPhlAn marker information failed bzip2 validation."
        return 1
    fi

    if ! bzip2 -t \
        "$cache_dir/${index}_species.txt.bz2" \
        >/dev/null 2>&1
    then
        log_warning "MetaPhlAn species information failed bzip2 validation."
        return 1
    fi

    newick_tail="$(
        tail -c 1024 "$cache_dir/${index}.nwk" 2>/dev/null |
            tr -d '\r\n[:space:]'
    )"

    if [[ -z "$newick_tail" || "${newick_tail: -1}" != ";" ]]; then
        log_warning "MetaPhlAn taxonomy tree does not end with ';'."
        return 1
    fi

    main_listing="$(mktemp)" || return 1

    bt2_listing="$(mktemp)" || {
        rm -f -- "$main_listing"
        return 1
    }

    if ! tar -tf "$main_archive" > "$main_listing"; then
        rm -f -- "$main_listing" "$bt2_listing"
        log_warning "Could not read the MetaPhlAn database archive listing."
        return 1
    fi

    if ! tar -tf "$bt2_archive" > "$bt2_listing"; then
        rm -f -- "$main_listing" "$bt2_listing"
        log_warning "Could not read the MetaPhlAn Bowtie2 archive listing."
        return 1
    fi

    if ! grep -Fxq "${index}.pkl" "$main_listing"; then
        rm -f -- "$main_listing" "$bt2_listing"
        log_warning "MetaPhlAn archive does not contain ${index}.pkl."
        return 1
    fi

    for required_file in "${required_bt2_files[@]}"; do
        if ! grep -Fxq "$required_file" "$bt2_listing"; then
            rm -f -- "$main_listing" "$bt2_listing"
            log_warning "MetaPhlAn Bowtie2 archive is missing: $required_file"
            return 1
        fi
    done

    rm -f -- "$main_listing" "$bt2_listing"

    return 0
}

prepare_humann_metaphlan_cache() {
    local cache_dir="$offline_files_folder/HUMAnN/$HUMANN_METAPHLAN_CACHE_DIRNAME"
    local index="$HUMANN_METAPHLAN_INDEX"
    local complete_marker="$cache_dir/.metaphlan_vJun23_202403_cache_complete"

    if ! mkdir -p "$cache_dir"; then
        log_error "Could not create the MetaPhlAn cache directory:"
        log_error "  $cache_dir"
        exit 1
    fi

    download_humann_metaphlan_cache_assets

    log_info "Validating the complete MetaPhlAn vJun23 cache..."

    if ! validate_humann_metaphlan_cache; then
        log_warning "The MetaPhlAn cache is incomplete or invalid."
        log_warning "All pinned MetaPhlAn assets will be downloaded again."

        rm -f -- \
            "$cache_dir/${index}.tar" \
            "$cache_dir/${index}.md5" \
            "$cache_dir/${index}_bt2.tar" \
            "$cache_dir/${index}_bt2.md5" \
            "$cache_dir/${index}.nwk" \
            "$cache_dir/${index}_marker_info.txt.bz2" \
            "$cache_dir/${index}_species.txt.bz2" \
            "$complete_marker"

        download_humann_metaphlan_cache_assets

        if ! validate_humann_metaphlan_cache; then
            log_error "MetaPhlAn cache validation failed after re-downloading."
            exit 1
        fi
    fi

    {
        echo "status=complete"
        echo "database_index=$index"
        echo "main_archive_md5=$HUMANN_METAPHLAN_MAIN_ARCHIVE_MD5"
        echo "bowtie2_archive_md5=$HUMANN_METAPHLAN_BT2_ARCHIVE_MD5"
        echo "validated_at=$(date --iso-8601=seconds)"
        echo "archive_directory=$cache_dir"
    } > "$complete_marker"

    log_ok "MetaPhlAn vJun23 persistent cache is complete and valid."
    log_info "MetaPhlAn cache:"
    log_info "  $cache_dir"
}

prepare_virushost_release_cache() {
    local base_url="$VIRUSHOST_MIRROR_BASE_URL"
    local release_dir
    local staging_dir
    local current_release
    local staged_release
    local current_manifest_hash
    local staged_manifest_hash
    local required_file
    local cache_is_complete=1

    local -a mirrored_files=(
        virushostdb.genomic.fna.gz
        non-segmented_virus_list.tsv
        segmented_virus_list.tsv
        dbrel.txt
        SHA256SUMS
        MIRROR_METADATA.tsv
    )

    release_dir="$offline_files_folder/Ref_genomes/MTD_virus/official_current"
    current_release="$release_dir/dbrel.txt"

    if ! mkdir -p "$release_dir"; then
        log_error "Could not create the Virus-Host DB cache directory:"
        log_error "  $release_dir"
        exit 1
    fi

    staging_dir="$(mktemp -d "$release_dir/.incoming.XXXXXX")" || {
        log_error "Could not create the Virus-Host DB staging directory."
        exit 1
    }

    staged_release="$staging_dir/dbrel.txt"

    log_info "Checking the pinned Virus-Host DB mirror release:"
    log_info "  $VIRUSHOST_MIRROR_TAG"

    # Download the small provenance and checksum files first.
    for required_file in \
        dbrel.txt \
        SHA256SUMS \
        MIRROR_METADATA.tsv
    do
        if ! retry_until_success \
            "Virus-Host DB mirror file: $required_file" \
            download_cache_file_once \
            "Virus-Host DB mirror file: $required_file" \
            "$base_url/$required_file" \
            "$staging_dir/$required_file"
        then
            rm -rf -- "$staging_dir"

            log_error "Could not retrieve the Virus-Host DB mirror metadata."
            exit 1
        fi
    done

    staged_manifest_hash="$(
        sha256sum "$staging_dir/SHA256SUMS" |
        awk '{print $1}'
    )"

    if [[ "$staged_manifest_hash" != "$VIRUSHOST_MIRROR_SHA256SUMS_SHA256" ]]
    then
        rm -rf -- "$staging_dir"

        log_error "The downloaded Virus-Host DB checksum manifest is unexpected."
        log_error "Expected:"
        log_error "  $VIRUSHOST_MIRROR_SHA256SUMS_SHA256"
        log_error "Observed:"
        log_error "  $staged_manifest_hash"
        exit 1
    fi

    # Check whether the existing cache is already this exact validated release.
    for required_file in "${mirrored_files[@]}"; do
        if [[ ! -s "$release_dir/$required_file" ]]; then
            cache_is_complete=0
            break
        fi
    done

    if (( cache_is_complete == 1 )) &&
       cmp -s "$current_release" "$staged_release"
    then
        current_manifest_hash="$(
            sha256sum "$release_dir/SHA256SUMS" |
            awk '{print $1}'
        )"

        if [[ "$current_manifest_hash" == "$VIRUSHOST_MIRROR_SHA256SUMS_SHA256" ]] &&
           (
               cd "$release_dir" &&
               sha256sum -c SHA256SUMS
           ) &&
           gzip -t "$release_dir/virushostdb.genomic.fna.gz"
        then
            rm -rf -- "$staging_dir"

            log_ok "Using the validated cached Virus-Host DB mirror."
            log_info "Mirror tag:"
            log_info "  $VIRUSHOST_MIRROR_TAG"

            while IFS= read -r release_line; do
                [[ -n "$release_line" ]] &&
                    log_info "  $release_line"
            done < "$current_release"

            return 0
        fi

        log_warning "The cached Virus-Host DB release failed validation."
        log_warning "It will be downloaded again from the pinned mirror."
    fi

    log_info "Downloading the validated Virus-Host DB mirror..."
    log_info "Mirror URL:"
    log_info "  $base_url"

    for required_file in \
        virushostdb.genomic.fna.gz \
        non-segmented_virus_list.tsv \
        segmented_virus_list.tsv
    do
        if ! retry_until_success \
            "Virus-Host DB mirror file: $required_file" \
            download_cache_file_once \
            "Virus-Host DB mirror file: $required_file" \
            "$base_url/$required_file" \
            "$staging_dir/$required_file"
        then
            rm -rf -- "$staging_dir"

            log_error "Could not synchronize the Virus-Host DB mirror."
            exit 1
        fi
    done

    log_info "Validating Virus-Host DB mirror checksums..."

    if ! (
        cd "$staging_dir" &&
        sha256sum -c SHA256SUMS
    ); then
        rm -rf -- "$staging_dir"

        log_error "Virus-Host DB mirror checksum validation failed."
        exit 1
    fi

    if ! gzip -t "$staging_dir/virushostdb.genomic.fna.gz"; then
        rm -rf -- "$staging_dir"

        log_error "The Virus-Host DB FASTA failed gzip validation."
        exit 1
    fi

    # Promote the complete validated release atomically, file by file.
    # dbrel.txt is moved last and acts as the release marker.
    for required_file in \
        virushostdb.genomic.fna.gz \
        non-segmented_virus_list.tsv \
        segmented_virus_list.tsv \
        SHA256SUMS \
        MIRROR_METADATA.tsv
    do
        if ! mv -f -- \
            "$staging_dir/$required_file" \
            "$release_dir/$required_file"
        then
            rm -rf -- "$staging_dir"

            log_error "Could not promote Virus-Host DB mirror file:"
            log_error "  $required_file"
            exit 1
        fi
    done

    if ! mv -f -- \
        "$staged_release" \
        "$current_release"
    then
        rm -rf -- "$staging_dir"

        log_error "Could not update the Virus-Host DB release marker."
        exit 1
    fi

    rm -rf -- "$staging_dir"

    # Derived files must be regenerated when the source release changes.
    rm -f -- \
        "$release_dir/virushostdb.genomic.fna" \
        "$release_dir/virushostdb_accession2taxid.tsv" \
        "$release_dir/virushostdb_accession_conflicts.tsv"

    log_ok "Validated Virus-Host DB mirror synchronized successfully."
    log_info "Mirror tag:"
    log_info "  $VIRUSHOST_MIRROR_TAG"
}

prepare_installation_cache() {
    log_info "Preparing persistent MTD installation cache:"
    log_info "  $offline_files_folder"

    if ! mkdir -p \
        "$offline_files_folder/Ref_genomes/MTD_virus" \
        "$offline_files_folder/Ref_genomes/MTD_virus/official_current" \
        "$offline_files_folder/HUMAnN" \
        "$offline_files_folder/Kraken2_taxonomy_cache" \
        "$offline_files_folder/Kraken2DB_micro/library/viral/all" \
        "$offline_files_folder/Kraken2DB_micro/library/bacteria/all" \
        "$offline_files_folder/Kraken2DB_micro/library/archaea/all" \
        "$offline_files_folder/Kraken2DB_micro/library/plasmid" \
        "$offline_files_folder/eggNOG/emapperdb-5.0.2" \
        "$offline_files_folder/Customized_hosts"
    then
        log_error "Could not create the persistent cache structure:"
        log_error "  $offline_files_folder"
        exit 1
    fi

    prepare_virushost_release_cache

    ensure_cached_file \
    "HUMAnN utility mapping database" \
    "https://huttenhower.sph.harvard.edu/humann_data/full_mapping_v201901b.tar.gz" \
    "$offline_files_folder/HUMAnN/full_mapping_v201901b.tar.gz"

    ensure_cached_file \
        "HUMAnN UniRef90 annotated database" \
        "https://huttenhower.sph.harvard.edu/humann_data/uniprot/uniref_annotated/uniref90_annotated_v201901b_full.tar.gz" \
        "$offline_files_folder/HUMAnN/uniref90_annotated_v201901b_full.tar.gz"

    ensure_cached_file \
        "HUMAnN ChocoPhlAn database" \
        "https://huttenhower.sph.harvard.edu/humann_data/chocophlan/full_chocophlan.v201901_v31.tar.gz" \
        "$offline_files_folder/HUMAnN/full_chocophlan.v201901_v31.tar.gz"

    prepare_humann_metaphlan_cache

    log_ok "Persistent installation cache is ready."
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

validate_kraken2_taxonomy_dir() {
    local taxonomy_dir="$1"
    local required_file

    for required_file in \
        names.dmp \
        nodes.dmp \
        nucl_gb.accession2taxid \
        nucl_wgs.accession2taxid
    do
        if [[ ! -s "$taxonomy_dir/$required_file" ]]; then
            log_warning "Missing or empty taxonomy file:"
            log_warning "  $taxonomy_dir/$required_file"
            return 1
        fi
    done

    if ! head -n 1 "$taxonomy_dir/nucl_gb.accession2taxid" |
        grep -q $'^accession\taccession.version\ttaxid'; then
        log_warning "Invalid header in nucl_gb.accession2taxid."
        return 1
    fi

    if ! head -n 1 "$taxonomy_dir/nucl_wgs.accession2taxid" |
        grep -q $'^accession\taccession.version\ttaxid'; then
        log_warning "Invalid header in nucl_wgs.accession2taxid."
        return 1
    fi

    if ! head -n 10 "$taxonomy_dir/names.dmp" | grep -q $'\t|\t'; then
        log_warning "names.dmp does not appear to have the expected format."
        return 1
    fi

    if ! head -n 10 "$taxonomy_dir/nodes.dmp" | grep -q $'\t|\t'; then
        log_warning "nodes.dmp does not appear to have the expected format."
        return 1
    fi

    return 0
}

download_shared_kraken2_taxonomy_once() {
    local downloader="$KRAKEN_AUX_DIR/kraken2-build-download-taxonomy"
    local taxonomy_dir="$KRAKEN_TAXONOMY_CACHE/taxonomy"
    local complete_marker="$KRAKEN_TAXONOMY_CACHE/.mtd_taxonomy_complete"

    if [[ ! -f "$downloader" ]]; then
        log_error "Kraken2 taxonomy downloader not found:"
        log_error "  $downloader"
        return 127
    fi

    if ! chmod +x "$downloader"; then
        log_error "Could not make taxonomy downloader executable:"
        log_error "  $downloader"
        return 126
    fi

    log_info "Removing any incomplete taxonomy cache before download..."

    rm -rf "$taxonomy_dir"
    rm -f "$complete_marker"

    if ! mkdir -p "$KRAKEN_TAXONOMY_CACHE"; then
        log_error "Could not create taxonomy cache directory:"
        log_error "  $KRAKEN_TAXONOMY_CACHE"
        return 1
    fi

    log_info "Downloading the shared NCBI taxonomy cache..."

    if ! "$downloader" \
        --download-taxonomy \
        --threads "$threads" \
        --db "$KRAKEN_TAXONOMY_CACHE"; then

        log_warning "Shared taxonomy download failed."
        log_warning "Removing incomplete files before the next attempt."

        rm -rf "$taxonomy_dir"
        rm -f "$complete_marker"
        return 1
    fi

    if ! validate_kraken2_taxonomy_dir "$taxonomy_dir"; then
        log_warning "Downloaded taxonomy did not pass validation."
        log_warning "Removing incomplete or invalid taxonomy files."

        rm -rf "$taxonomy_dir"
        rm -f "$complete_marker"
        return 1
    fi

    date -u '+%Y-%m-%dT%H:%M:%SZ' > "$complete_marker"

    log_ok "Shared Kraken2 taxonomy downloaded and validated."
    log_info "Taxonomy cache:"
    log_info "  $KRAKEN_TAXONOMY_CACHE"

    return 0
}

prepare_shared_kraken2_taxonomy() {
    local taxonomy_dir="$KRAKEN_TAXONOMY_CACHE/taxonomy"
    local complete_marker="$KRAKEN_TAXONOMY_CACHE/.mtd_taxonomy_complete"

    if [[ -f "$complete_marker" ]] &&
       validate_kraken2_taxonomy_dir "$taxonomy_dir"; then

        log_ok "Using the existing validated Kraken2 taxonomy cache."
        log_info "Taxonomy cache:"
        log_info "  $KRAKEN_TAXONOMY_CACHE"
        return 0
    fi

    if [[ -e "$taxonomy_dir" || -e "$complete_marker" ]]; then
        log_warning "Existing taxonomy cache is incomplete or invalid."
        log_warning "It will be removed and downloaded again."
        rm -rf "$taxonomy_dir"
        rm -f "$complete_marker"
    else
        log_info "No valid shared Kraken2 taxonomy cache was found."
    fi

    if ! retry_until_success \
        "Shared NCBI taxonomy for all Kraken2 databases" \
        download_shared_kraken2_taxonomy_once; then

        log_error "The shared Kraken2 taxonomy could not be prepared."
        exit 1
    fi
}

install_shared_kraken2_taxonomy() {
    local database="$1"
    local source_taxonomy="$KRAKEN_TAXONOMY_CACHE/taxonomy"
    local destination_taxonomy="$database/taxonomy"

    prepare_shared_kraken2_taxonomy

    log_info "Installing cached NCBI taxonomy into:"
    log_info "  $database"

    if ! mkdir -p "$database"; then
        log_error "Could not create Kraken2 database directory:"
        log_error "  $database"
        exit 1
    fi

    rm -rf "$destination_taxonomy"

    if ! cp -a "$source_taxonomy" "$database/"; then
        log_error "Could not copy the shared taxonomy into:"
        log_error "  $database"
        exit 1
    fi

    if ! validate_kraken2_taxonomy_dir "$destination_taxonomy"; then
        log_error "Copied taxonomy failed validation:"
        log_error "  $destination_taxonomy"
        exit 1
    fi

    log_ok "Cached taxonomy installed and validated:"
    log_ok "  $destination_taxonomy"
}

build_kraken2_database() {
    local database="$1"

    kraken2-build \
        --build \
        --threads "$threads" \
        --db "$database" \
        "${kraken_build_opts[@]}"
}
validate_built_kraken2_database() {
    local database="$1"
    local required_file

    for required_file in \
        hash.k2d \
        opts.k2d \
        taxo.k2d
    do
        if [[ ! -s "$database/$required_file" ]]; then
            log_error "Required Kraken2 database file is missing or empty:"
            log_error "  $database/$required_file"
            exit 1
        fi
    done

    if [[ -s "$database/unmapped.txt" ]]; then
        log_error "Kraken2 left reference sequences without taxonomy mapping:"
        log_error "  $database/unmapped.txt"

        head -n 20 \
            "$database/unmapped.txt" \
            >&2

        exit 1
    fi

    log_ok "Kraken2 database outputs were created successfully:"
    log_ok "  $database"
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

copy_manifest_with_offline_folder() {
    local source_script="$1"
    local destination_script="$2"

    copy_required_file "$source_script" "$destination_script"
    sed -i \
    "s|^offline_files_folder=.*|offline_files_folder=\"$offline_files_folder\"|" \
    "$destination_script"
    run_required_script "$destination_script"
}


# ------------------------------------------------------------------------------
# Installation stages
# ------------------------------------------------------------------------------

initialize_installation() {
    cd "$dir" || exit 1
    printf '%s\n' "$condapath" > "$dir/condaPath"
    printf '%s\n' "$offline_files_folder" > "$dir/offlineCachePath"
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
    local -a system_packages=(
        libgeos-dev
        libharfbuzz-dev
        libfribidi-dev
        libfreetype6-dev
        libpng-dev
        libtiff5-dev
        libjpeg-dev
        rsync
        pigz
        curl
        wget
        aria2
        parallel
        ca-certificates
        build-essential
        pkg-config
        libssl-dev
    )

    log_info "Installing system dependencies..."

    run_required_command \
        "Updating APT package indexes" \
        run_as_root apt-get update

    run_required_command \
        "Installing required system packages" \
        run_as_root env \
            DEBIAN_FRONTEND=noninteractive \
            apt-get install -y "${system_packages[@]}"

    run_required_command \
        "Installing SRA Toolkit without optional Debian Med configuration" \
        run_as_root env \
            DEBIAN_FRONTEND=noninteractive \
            apt-get install -y --no-install-recommends sra-toolkit

    log_ok "System dependencies installed successfully."
}

validate_humann_environment() {
    local env_name="MTD_humann"
    local command_name=""
    local humann_version=""
    local metaphlan_version=""

    local -a required_commands=(
        python
        humann
        humann_config
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

    for command_name in "${required_commands[@]}"; do
        require_env_command "$env_name" "$command_name"
    done

    run_required_command \
        "Validating MTD_humann Python isolation and core modules" \
        conda run -n "$env_name" \
        python -c \
        'import os, site, sys, Cython, humann, numpy, pysam, simplejson; prefix=os.path.realpath(os.environ["CONDA_PREFIX"]); assert os.path.realpath(sys.prefix) == prefix; assert numpy.__version__ == "1.26.4"; assert simplejson.__version__.split(".", 1)[0] == "3"; assert site.ENABLE_USER_SITE is False; assert all("/.local/" not in entry for entry in sys.path); modules=(Cython, humann, numpy, pysam, simplejson); assert all(os.path.realpath(module.__file__).startswith(prefix + os.sep) for module in modules); print("MTD_humann Python isolation: OK")'

    humann_version="$(
        conda run -n "$env_name" humann --version 2>&1
    )"

    if [[ "$humann_version" != *"humann v3.9"* ]]; then
        log_error "Unexpected HUMAnN version in $env_name:"
        log_error "  $humann_version"
        exit 1
    fi

    metaphlan_version="$(
        conda run -n "$env_name" metaphlan --version 2>&1
    )"

    if [[ "$metaphlan_version" != *"MetaPhlAn version 4.1.1"* ]]; then
        log_error "Unexpected MetaPhlAn version in $env_name:"
        log_error "  $metaphlan_version"
        exit 1
    fi

    log_ok "Dedicated MTD_humann environment passed validation."
    log_info "  $humann_version"
    log_info "  $metaphlan_version"
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
        "Creating dedicated MTD_humann environment" \
        conda env create -f "$dir/Installation/MTD_humann.yml"

    validate_humann_environment

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

    # ----------------------------------------------------------
    # Dedicated environment for custom OrgDb construction
    # ----------------------------------------------------------

    run_required_command \
        "Creating dedicated MTD_orgdb environment" \
        conda env create -f "$dir/Installation/MTD_orgdb.yml"

    require_env_command MTD_orgdb Rscript
    require_env_command MTD_orgdb jq
    require_env_command MTD_orgdb yq

    run_required_command \
        "Validating MTD_orgdb R packages" \
        conda run -n MTD_orgdb \
        Rscript "$dir/Installation/check_MTD_orgdb.R"

    sed -i '/^# *rpy2/s/^# *//' "$dir/Installation/pip.requirements"

    chmod +x "$dir/aux_scripts/ssGSEA/resolve_ssgsea_go_terms.py"
    run_required_command \
        "Checking ssGSEA GO resolver syntax" \
        python3 -m py_compile "$dir/aux_scripts/ssGSEA/resolve_ssgsea_go_terms.py"

    require_env_command MTD kraken2-build
    require_env_command MTD bracken-build
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
    R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/mclust/mclust_6.1.2.tar.gz', repos=NULL, type='source', Ncpus=$threads)"
    R -e 'install.packages("BiocManager", repos = "https://cloud.r-project.org")'
    R -e "install.packages('$dir/update_fix/pvr_pkg/MASS_7.3-60.tar.gz', repos=NULL, type='source')"
    R -e "install.packages('$dir/update_fix/pvr_pkg/preprocessCore_1.72.0.tar.gz', repos=NULL, type='source')"
    R -e 'install.packages("remotes", repos="https://cloud.r-project.org")'
    R -e 'remotes::install_url("https://cran.r-project.org/src/contrib/EnvStats_3.1.0.tar.gz", dependencies=TRUE)'
    R -e 'remotes::install_version("Hmisc", version = "4.8-0", repos = "https://cloud.r-project.org")'
#    R -e "install.packages('https://cran.r-project.org/src/contrib/00Archive/eva/eva_0.2.6.tar.gz', repos=NULL, type='source')"
    R -e "install.packages('$dir/update_fix/pvr_pkg/eva_0.2.6.tar.gz', repos=NULL, type='source', Ncpus=$threads)"

    run_required_command \
        "Applying HAllA Matplotlib compatibility patch" \
        conda run -n halla0820 \
        python "$dir/update_fix/patch_halla_matplotlib.py"

    run_required_command \
        "Validating HAllA Matplotlib compatibility patch" \
        conda run -n halla0820 \
        python "$dir/update_fix/patch_halla_matplotlib.py" --check

    run_required_command \
        "Checking R packages in the halla0820 environment" \
        conda run -n halla0820 \
        bash "$dir/update_fix/check_R_pkg.halla0820.sh"

    safe_conda_deactivate
}

install_mtd_extra_tools() {
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
    local viral_cache_dir
    local viral_download_dir
    local manifest_destination
    local all_viral_fasta
    local siv_fasta

    local combined_fasta
    local combined_summary
    local combined_details

    local viral_builder
    local refseq_map_builder
    local virushost_map_builder
    local taxonomy_dir

    local virushost_release_dir
    local virushost_fasta_gz
    local virushost_fasta
    local virushost_release_file
    local virushost_non_segmented
    local virushost_segmented
    local virushost_taxid_map
    local virushost_conflicts

    local refseq_assembly_summary
    local refseq_taxid_map
    local refseq_conflicts

    local unresolved_taxids
    local records_without_accession
    local required_file
    local decompressed_tmp

    viral_cache_dir="$offline_files_folder/Kraken2DB_micro/library/viral"
    viral_download_dir="$viral_cache_dir/all"

    manifest_destination="$offline_files_folder/Kraken2DB_micro/library/manifest.virus.sh"

    all_viral_fasta="$viral_cache_dir/all_viral_genomes.fna"

    refseq_assembly_summary="$viral_cache_dir/assembly_summary_viral.txt"
    refseq_taxid_map="$viral_cache_dir/refseq_viral_accession2taxid.tsv"
    refseq_conflicts="$viral_cache_dir/refseq_viral_accession_conflicts.tsv"

    virushost_release_dir="$offline_files_folder/Ref_genomes/MTD_virus/official_current"

    virushost_fasta_gz="$virushost_release_dir/virushostdb.genomic.fna.gz"
    virushost_fasta="$virushost_release_dir/virushostdb.genomic.fna"
    virushost_release_file="$virushost_release_dir/dbrel.txt"

    virushost_non_segmented="$virushost_release_dir/non-segmented_virus_list.tsv"
    virushost_segmented="$virushost_release_dir/segmented_virus_list.tsv"

    virushost_taxid_map="$virushost_release_dir/virushostdb_accession2taxid.tsv"
    virushost_conflicts="$virushost_release_dir/virushostdb_accession_conflicts.tsv"

    siv_fasta="$dir/Installation/M33262_SIVMM239.fa"

    combined_fasta="$viral_cache_dir/viral_genomes_combined_nonredundant.fna"
    combined_summary="$viral_cache_dir/viral_genomes_combined_nonredundant.summary.tsv"
    combined_details="$viral_cache_dir/viral_genomes_combined_nonredundant.details.tsv"

    viral_builder="$dir/aux_scripts/Kraken2/build_nonredundant_viral_fasta.py"
    refseq_map_builder="$dir/aux_scripts/Kraken2/build_refseq_viral_taxid_map.py"
    virushost_map_builder="$dir/aux_scripts/Kraken2/build_virushost_taxid_map.py"

    taxonomy_dir="$KRAKEN_TAXONOMY_CACHE/taxonomy"

    for required_file in \
        "$virushost_fasta_gz" \
        "$virushost_release_file" \
        "$virushost_non_segmented" \
        "$virushost_segmented" \
        "$siv_fasta" \
        "$viral_builder" \
        "$refseq_map_builder" \
        "$virushost_map_builder"
    do
        if [[ ! -s "$required_file" ]]; then
            log_error "Required virome input is missing or empty:"
            log_error "  $required_file"
            exit 1
        fi
    done

    log_info "Using the official Virus-Host DB release:"

    while IFS= read -r release_line; do
        [[ -n "$release_line" ]] &&
            log_info "  $release_line"
    done < "$virushost_release_file"

    if [[ ! -s "$virushost_fasta" ||
          "$virushost_fasta_gz" -nt "$virushost_fasta" ]]
    then
        log_info "Decompressing the official Virus-Host DB FASTA..."

        decompressed_tmp="${virushost_fasta}.tmp.$$"
        rm -f -- "$decompressed_tmp"

        if ! gzip -dc -- "$virushost_fasta_gz" > "$decompressed_tmp"; then
            rm -f -- "$decompressed_tmp"
            log_error "Could not decompress the official Virus-Host DB FASTA."
            exit 1
        fi

        if [[ ! -s "$decompressed_tmp" ]]; then
            rm -f -- "$decompressed_tmp"
            log_error "The decompressed Virus-Host DB FASTA is empty."
            exit 1
        fi

        mv -f -- "$decompressed_tmp" "$virushost_fasta"
    else
        log_ok "Using the existing decompressed Virus-Host DB FASTA."
    fi

    log_info "Synchronizing the NCBI RefSeq viral collection..."

copy_required_file \
    "$MANIFEST_SCRIPTS_DIR/manifest.virus.sh" \
    "$manifest_destination"

chmod +x "$manifest_destination"

run_required_command \
    "Synchronizing NCBI RefSeq viral genomes" \
    env \
    MTD_OFFLINE_FILES_FOLDER="$offline_files_folder" \
    BUILD_COMBINED_FASTA=1 \
    REQUIRE_COMPLETE_COLLECTION=1 \
    "$manifest_destination"

    if [[ ! -s "$all_viral_fasta" ]]; then
        log_error "The RefSeq viral combined FASTA is missing or empty:"
        log_error "  $all_viral_fasta"
        exit 1
    fi

    if [[ ! -s "$refseq_assembly_summary" ]]; then
        log_error "The RefSeq viral assembly summary is missing or empty:"
        log_error "  $refseq_assembly_summary"
        exit 1
    fi

    log_info "Preparing taxonomy for viral deduplication..."

    prepare_shared_kraken2_taxonomy

    if ! validate_kraken2_taxonomy_dir "$taxonomy_dir"; then
        log_error "Kraken2 taxonomy is unavailable or invalid:"
        log_error "  $taxonomy_dir"
        exit 1
    fi

    run_required_command \
        "Building the RefSeq viral accession-to-TaxID map" \
        python3 "$refseq_map_builder" \
        --assembly-summary "$refseq_assembly_summary" \
        --fasta-dir "$viral_download_dir" \
        --output "$refseq_taxid_map" \
        --conflicts "$refseq_conflicts"

    run_required_command \
        "Building the Virus-Host DB accession-to-TaxID map" \
        python3 "$virushost_map_builder" \
        --non-segmented "$virushost_non_segmented" \
        --segmented "$virushost_segmented" \
        --output "$virushost_taxid_map" \
        --conflicts "$virushost_conflicts"

    run_required_command \
        "Combining and deduplicating viral reference sequences" \
        python3 "$viral_builder" \
        --primary "$all_viral_fasta" \
        --primary-taxid-map "$refseq_taxid_map" \
        --virushost "$virushost_fasta" \
        --virushost-taxid-map "$virushost_taxid_map" \
        --extra "$siv_fasta" \
        --taxonomy-dir "$taxonomy_dir" \
        --output "$combined_fasta" \
        --summary "$combined_summary" \
        --details "$combined_details"

    for required_file in \
        "$combined_fasta" \
        "$combined_summary" \
        "$combined_details"
    do
        if [[ ! -s "$required_file" ]]; then
            log_error "Required combined viral output is missing or empty:"
            log_error "  $required_file"
            exit 1
        fi
    done

    unresolved_taxids="$(
        awk -F $'\t' \
            '$1 == "records_without_taxid" { print $2; exit }' \
            "$combined_summary"
    )"

    records_without_accession="$(
        awk -F $'\t' \
            '$1 == "records_without_accession" { print $2; exit }' \
            "$combined_summary"
    )"

    if ! [[ "$unresolved_taxids" =~ ^[0-9]+$ ]]; then
        log_error "Could not read records_without_taxid from:"
        log_error "  $combined_summary"
        exit 1
    fi

    if ! [[ "$records_without_accession" =~ ^[0-9]+$ ]]; then
        log_error "Could not read records_without_accession from:"
        log_error "  $combined_summary"
        exit 1
    fi

    if (( unresolved_taxids != 0 ||
          records_without_accession != 0 ))
    then
        log_error "The viral collection is not safe to add to Kraken2."
        log_error "Records without TaxID:     $unresolved_taxids"
        log_error "Records without accession: $records_without_accession"
        log_error "Details:"
        log_error "  $combined_details"
        exit 1
    fi

    log_ok "Nonredundant viral collection completed with complete taxonomy:"
    log_ok "  $combined_fasta"

    log_info "Viral deduplication summary:"

    while IFS=$'\t' read -r metric value; do
        [[ "$metric" == "metric" ]] && continue
        printf '  %-55s %s\n' "$metric" "$value"
    done < "$combined_summary"
}

install_default_kraken_helpers() {
    restore_default_rsync_helper
    restore_default_genomic_library_helper
}

prepare_microbiome_manifests() {
    copy_manifest_with_offline_folder \
        "$MANIFEST_SCRIPTS_DIR/manifest.bacteria.sh" \
        "$offline_files_folder/Kraken2DB_micro/library/manifest.bacteria.sh"

    copy_required_file \
        "$MANIFEST_SCRIPTS_DIR/manifest.sh" \
        "$offline_files_folder/Kraken2DB_micro/library/manifest.sh"

    sed -i \
        "s|^offline_files_folder=.*|offline_files_folder=$offline_files_folder|" \
        "$offline_files_folder/Kraken2DB_micro/library/manifest.sh"
}

add_local_archaea_library() {
    local database="$1"
    local helper="$KRAKEN_ENV_LIBEXEC/rsync_from_ncbi.pl"

    copy_manifest_with_offline_folder \
        "$MANIFEST_SCRIPTS_DIR/manifest.archea.sh" \
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
    local manifest_destination
    local custom_helper
    local plasmid_cache
    local plasmid_cache_count
    local download_status=0

    manifest_destination="$offline_files_folder/Kraken2DB_micro/library/manifest.plasmid.sh"
    custom_helper="$dir/Installation/download_genomic_library_plasmid.sh"
    plasmid_cache="$offline_files_folder/Kraken2DB_micro/library/plasmid"

    log_info "Preparing cached plasmid genomic sequences..."

    copy_required_file \
    "$MANIFEST_SCRIPTS_DIR/manifest.plasmid.sh" \
    "$manifest_destination"

    sed -i \
        "s|^LOCAL_DIR=.*|LOCAL_DIR=$plasmid_cache/|" \
        "$manifest_destination"

    run_required_script "$manifest_destination"

    if [[ ! -d "$plasmid_cache" ]]; then
        log_error "Plasmid cache directory not found:"
        log_error "  $plasmid_cache"
        exit 1
    fi

    plasmid_cache_count="$(
        find "$plasmid_cache" \
            -maxdepth 1 \
            -type f \
            -name '*.genomic.fna.gz' \
            -size +0c |
        wc -l
    )"

    if (( plasmid_cache_count == 0 )); then
        log_error "No usable genomic plasmid FASTA archives were found:"
        log_error "  $plasmid_cache"
        exit 1
    fi

    log_ok "Usable genomic plasmid FASTA archives: $plasmid_cache_count"
    log_info "Plasmid cache:"
    log_info "  $plasmid_cache"

    export MTD_KRAKEN2_PLASMID_CACHE="$plasmid_cache"

    install_kraken_helper \
        "$custom_helper" \
        "download_genomic_library.sh"

    log_info "Adding local plasmid sequences to Kraken2 database..."

    retry_until_success \
        "Kraken2 plasmid library for database '$database'" \
        kraken2-build \
        --download-library plasmid \
        --threads "$threads" \
        --db "$database" ||
        download_status=$?

    log_info "Restoring the standard Kraken2 genomic-library helper..."

    restore_default_genomic_library_helper

    unset MTD_KRAKEN2_PLASMID_CACHE

    if (( download_status != 0 )); then
        log_error "The Kraken2 plasmid library could not be prepared."
        exit "$download_status"
    fi

    local plasmid_library_dir="$database/library/plasmid"
    local required_file

    for required_file in \
        library.fna \
        manifest.txt \
        prelim_map.txt
    do
        if [[ ! -s "$plasmid_library_dir/$required_file" ]]; then
            log_error "Required plasmid-library output is missing or empty:"
            log_error "  $plasmid_library_dir/$required_file"
            exit 1
        fi
    done

    log_ok "Plasmid library prepared and validated:"
    log_ok "  $plasmid_library_dir"
}

build_microbiome_kraken_database() {
    local database="$dir/kraken2DB_micro"
    local viral_library

    viral_library="$offline_files_folder/Kraken2DB_micro/library/viral/viral_genomes_combined_nonredundant.fna"

    show_progress 40 "Preparing microbiome manifests"
    prepare_microbiome_manifests

    chmod +x "$KRAKEN_AUX_DIR/kraken2-build-download-taxonomy"

    show_progress 42 "Preparing the shared NCBI taxonomy"
    install_shared_kraken2_taxonomy "$database"

    show_progress 45 "Adding archaeal genomes to the microbiome database"
    add_local_archaea_library "$database"

    show_progress 49 "Adding bacterial genomes to the microbiome database"
    add_local_bacteria_library "$database"

    show_progress 55 "Adding RefSeq protozoan genomes"
    download_kraken2_library_until_success "$database" "protozoa" --use-ftp

    show_progress 58 "Adding RefSeq fungal genomes"
    download_kraken2_library_until_success "$database" "fungi" --use-ftp

    show_progress 61 "Adding plasmid sequences"
    add_local_plasmid_library "$database"

    show_progress 64 "Adding the UniVec_Core library"
    download_kraken2_library_until_success "$database" "UniVec_Core" --use-ftp

    show_progress 66 "Adding the nonredundant viral collection"

    if [[ ! -s "$viral_library" ]]; then
        log_error "The nonredundant viral library is missing or empty:"
        log_error " $viral_library"
        exit 1
    fi

    run_required_command \
    "Adding the nonredundant viral collection to Kraken2" \
    kraken2-build \
    --add-to-library "$viral_library" \
    --threads "$threads" \
    --db "$database"

    show_progress 68 "Building the final Kraken2 microbiome database"

run_required_command \
    "Building the final Kraken2 microbiome database" \
    build_kraken2_database \
    "$database"

validate_built_kraken2_database \
    "$database"

show_progress 72 "Kraken2 microbiome database completed"
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

validate_installed_humann_databases() {
    local humann_dir="$1"
    local chocophlan_dir="$humann_dir/chocophlan"
    local uniref_dir="$humann_dir/uniref"
    local utility_dir="$humann_dir/utility_mapping"
    local metaphlan_dir="$humann_dir/metaphlan"
    local index="$HUMANN_METAPHLAN_INDEX"
    local required_file=""
    local utility_count=0

    if ! find "$chocophlan_dir" \
        -type f \
        -name '*.ffn.gz' \
        -size +0c \
        -print -quit 2>/dev/null | grep -q .
    then
        log_warning "No usable ChocoPhlAn .ffn.gz files were found:"
        log_warning "  $chocophlan_dir"
        return 1
    fi

    if [[ ! -s "$uniref_dir/uniref90_201901b_full.dmnd" ]]; then
        log_warning "The HUMAnN UniRef90 DIAMOND database is missing or empty:"
        log_warning "  $uniref_dir/uniref90_201901b_full.dmnd"
        return 1
    fi

    utility_count="$(
        find "$utility_dir" \
            -type f \
            -size +0c \
            2>/dev/null | awk 'END { print NR + 0 }'
    )"

    if ! [[ "$utility_count" =~ ^[0-9]+$ ]] ||
       (( utility_count < 10 ))
    then
        log_warning "The HUMAnN utility-mapping installation is incomplete:"
        log_warning "  $utility_dir"
        log_warning "Non-empty files detected: ${utility_count:-invalid}"
        return 1
    fi

    for required_file in \
        "$metaphlan_dir/${index}.pkl" \
        "$metaphlan_dir/${index}.1.bt2l" \
        "$metaphlan_dir/${index}.2.bt2l" \
        "$metaphlan_dir/${index}.3.bt2l" \
        "$metaphlan_dir/${index}.4.bt2l" \
        "$metaphlan_dir/${index}.rev.1.bt2l" \
        "$metaphlan_dir/${index}.rev.2.bt2l" \
        "$metaphlan_dir/${index}.nwk" \
        "$metaphlan_dir/${index}_marker_info.txt.bz2" \
        "$metaphlan_dir/${index}_species.txt.bz2"
    do
        if [[ ! -s "$required_file" ]]; then
            log_warning "Required MetaPhlAn database file is missing or empty:"
            log_warning "  $required_file"
            return 1
        fi
    done

    if ! bzip2 -t \
        "$metaphlan_dir/${index}_marker_info.txt.bz2" \
        >/dev/null 2>&1
    then
        log_warning "Installed MetaPhlAn marker information is invalid."
        return 1
    fi

    if ! bzip2 -t \
        "$metaphlan_dir/${index}_species.txt.bz2" \
        >/dev/null 2>&1
    then
        log_warning "Installed MetaPhlAn species information is invalid."
        return 1
    fi

    return 0
}

configure_humann_database_paths() {
    local humann_dir="$1"
    local env_name="MTD_humann"
    local config_output=""
    local expected_path=""

    run_required_command \
        "Configuring the HUMAnN nucleotide database" \
        conda run -n "$env_name" \
        humann_config --update database_folders nucleotide \
        "$humann_dir/chocophlan"

    run_required_command \
        "Configuring the HUMAnN protein database" \
        conda run -n "$env_name" \
        humann_config --update database_folders protein \
        "$humann_dir/uniref"

    run_required_command \
        "Configuring the HUMAnN utility-mapping database" \
        conda run -n "$env_name" \
        humann_config --update database_folders utility_mapping \
        "$humann_dir/utility_mapping"

    config_output="$(
        conda run -n "$env_name" humann_config --print 2>&1
    )"

    for expected_path in \
        "$humann_dir/chocophlan" \
        "$humann_dir/uniref" \
        "$humann_dir/utility_mapping"
    do
        if ! grep -Fq "$expected_path" <<< "$config_output"; then
            log_error "HUMAnN configuration does not contain the expected path:"
            log_error "  $expected_path"
            log_error "Current HUMAnN configuration:"
            printf '%s\n' "$config_output" >&2
            exit 1
        fi
    done

    log_ok "HUMAnN database paths were configured in MTD_humann."
}

install_humann_databases() {
    local humann_dir="$dir/HUMAnN/ref_database"
    local chocophlan_dir="$humann_dir/chocophlan"
    local uniref_dir="$humann_dir/uniref"
    local utility_dir="$humann_dir/utility_mapping"
    local metaphlan_dir="$humann_dir/metaphlan"
    local cache_dir="$offline_files_folder/HUMAnN"
    local metaphlan_cache="$cache_dir/$HUMANN_METAPHLAN_CACHE_DIRNAME"
    local index="$HUMANN_METAPHLAN_INDEX"
    local complete_marker="$humann_dir/.mtd_humann_databases_complete"
    local installed_is_valid=0

    validate_humann_environment

    if [[ -f "$complete_marker" ]] &&
       validate_installed_humann_databases "$humann_dir"
    then
        installed_is_valid=1
        log_ok "Using the existing validated HUMAnN and MetaPhlAn databases."
        log_info "  $humann_dir"
    fi

    if (( installed_is_valid == 0 )); then
        log_info "Installing HUMAnN and MetaPhlAn databases from the persistent cache..."

        if ! validate_humann_metaphlan_cache; then
            log_error "The persistent MetaPhlAn cache failed validation."
            exit 1
        fi

        for required_file in \
            "$cache_dir/full_chocophlan.v201901_v31.tar.gz" \
            "$cache_dir/uniref90_annotated_v201901b_full.tar.gz" \
            "$cache_dir/full_mapping_v201901b.tar.gz" \
            "$metaphlan_cache/${index}.tar" \
            "$metaphlan_cache/${index}_bt2.tar" \
            "$metaphlan_cache/${index}.nwk" \
            "$metaphlan_cache/${index}_marker_info.txt.bz2" \
            "$metaphlan_cache/${index}_species.txt.bz2"
        do
            if [[ ! -s "$required_file" ]]; then
                log_error "Required HUMAnN/MetaPhlAn cache file is missing or empty:"
                log_error "  $required_file"
                exit 1
            fi
        done

        rm -rf -- "$humann_dir"

        if ! mkdir -p \
            "$chocophlan_dir" \
            "$uniref_dir" \
            "$utility_dir" \
            "$metaphlan_dir"
        then
            log_error "Could not create the HUMAnN database directories:"
            log_error "  $humann_dir"
            exit 1
        fi

        run_required_command \
            "Extracting the HUMAnN ChocoPhlAn database" \
            tar -xzf \
            "$cache_dir/full_chocophlan.v201901_v31.tar.gz" \
            -C "$chocophlan_dir"

        run_required_command \
            "Extracting the HUMAnN UniRef90 database" \
            tar -xzf \
            "$cache_dir/uniref90_annotated_v201901b_full.tar.gz" \
            -C "$uniref_dir"

        run_required_command \
            "Extracting the HUMAnN utility-mapping database" \
            tar -xzf \
            "$cache_dir/full_mapping_v201901b.tar.gz" \
            -C "$utility_dir"

        run_required_command \
            "Extracting the MetaPhlAn database metadata" \
            tar -xf \
            "$metaphlan_cache/${index}.tar" \
            -C "$metaphlan_dir"

        run_required_command \
            "Extracting the MetaPhlAn Bowtie2 indexes" \
            tar -xf \
            "$metaphlan_cache/${index}_bt2.tar" \
            -C "$metaphlan_dir"

        copy_required_file \
            "$metaphlan_cache/${index}.nwk" \
            "$metaphlan_dir/${index}.nwk"

        copy_required_file \
            "$metaphlan_cache/${index}_marker_info.txt.bz2" \
            "$metaphlan_dir/${index}_marker_info.txt.bz2"

        copy_required_file \
            "$metaphlan_cache/${index}_species.txt.bz2" \
            "$metaphlan_dir/${index}_species.txt.bz2"

        if ! validate_installed_humann_databases "$humann_dir"; then
            log_error "The installed HUMAnN/MetaPhlAn databases failed validation."
            exit 1
        fi

        {
            echo "status=complete"
            echo "humann_environment=MTD_humann"
            echo "humann_version=3.9"
            echo "metaphlan_version=4.1.1"
            echo "metaphlan_index=$index"
            echo "nucleotide_database=$chocophlan_dir"
            echo "protein_database=$uniref_dir"
            echo "utility_mapping_database=$utility_dir"
            echo "metaphlan_database=$metaphlan_dir"
            echo "installed_at=$(date --iso-8601=seconds)"
        } > "$complete_marker"

        log_ok "HUMAnN and MetaPhlAn databases were installed and validated."
    fi

    configure_humann_database_paths "$humann_dir"

    if ! validate_installed_humann_databases "$humann_dir"; then
        log_error "Final HUMAnN/MetaPhlAn database validation failed."
        exit 1
    fi

    log_ok "HUMAnN 3.9 database installation is ready."
    log_info "HUMAnN database root:"
    log_info "  $humann_dir"
    log_info "MetaPhlAn database index:"
    log_info "  $HUMANN_METAPHLAN_INDEX"
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

    log_info "OrgDb annotation packages are managed by the dedicated MTD_orgdb environment."

    run_required_command \
        "Rechecking MTD_orgdb environment" \
        conda run -n MTD_orgdb \
        Rscript "$dir/Installation/check_MTD_orgdb.R"
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

    echo
    echo "${g}MTD_orgdb environment:${w}"

    run_required_command \
        "Reporting MTD_orgdb package versions" \
        conda run -n MTD_orgdb \
        Rscript "$dir/Installation/check_MTD_orgdb.R"
}

# ------------------------------------------------------------------------------
# Main installation sequence
# ------------------------------------------------------------------------------

main() {
    init_colors

    parse_arguments "$@"
    validate_arguments

    show_progress 2 "Installing operating-system dependencies"
    install_system_dependencies

    show_progress 5 "Installing Miniconda"
    install_miniconda

    configure_paths_and_options
    initialize_installation

    show_progress 8 "Preparing the persistent installation cache"
    prepare_installation_cache

    show_progress 12 "Creating the MTD Conda environments"
    create_conda_environments

    show_progress 27 "Installing HAllA dependencies"
    install_halla_dependencies

    show_progress 34 "Installing MTD tools and preparing virome files"
    install_mtd_extra_tools
    prepare_virome_files
    install_default_kraken_helpers

    build_microbiome_kraken_database

    show_progress 72 "Building the Bracken database"

run_required_command \
    "Building the Bracken database" \
    build_bracken_database

    show_progress 82 "Installing and configuring HUMAnN databases"
    install_humann_databases

    show_progress 90 "Installing R412 and annotation packages"
    install_r412_and_annotation_packages

    show_progress 98 "Validating installed R packages"
    show_r_package_versions

    show_progress 100 "MTD Explorer installation completed successfully"

    echo
    log_ok "MTD Explorer is ready."
    log_info "No default host database was installed."
    log_info "Create a host database before analyzing real data:"
    log_info "  $dir/Create_custom_host.sh"
    echo
}

main "$@"
