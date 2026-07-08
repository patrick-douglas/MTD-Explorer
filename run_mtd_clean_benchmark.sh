#!/usr/bin/env bash

set -Eeuo pipefail

PROGRAM_NAME="$(basename "$0")"
MACHINE_NAME=""
CACHE_INPUT=""

usage() {
    cat <<USAGE
Usage:
  bash $PROGRAM_NAME -n MACHINE -o CACHE_DIRECTORY

Required:
  -n NAME   Machine name used in the benchmark label, e.g. s5 or master
  -o PATH   Persistent MTD cache passed to Install.sh with -o

Examples:
  bash $PROGRAM_NAME -n s5 -o "$HOME/MTD_install_cache"
  bash $PROGRAM_NAME -n master -o /media/me/18TB_BACKUP_LBN/MTD_cache/installer

Behavior:
  - Detects whether the cache is cold or warm.
  - Creates the cache directory if it does not exist.
  - Clones and validates the current GitHub repository before deleting anything.
  - Preserves the cache while removing the previous MTD/Conda benchmark state.
  - Generates Install_profiled.sh and runs the installation benchmark.
  - Does not use the removed -a installer option.
USAGE
}

die() {
    printf '[STOP] %s\n' "$*" >&2
    exit 1
}

while getopts ":n:o:h" option; do
    case "$option" in
        n) MACHINE_NAME="$OPTARG" ;;
        o) CACHE_INPUT="$OPTARG" ;;
        h)
            usage
            exit 0
            ;;
        :)
            die "Option -$OPTARG requires a value."
            ;;
        \?)
            die "Unknown option: -$OPTARG"
            ;;
    esac
done

shift $((OPTIND - 1))

[[ $# -eq 0 ]] || die "Unexpected positional arguments: $*"
[[ -n "$MACHINE_NAME" ]] || die "Machine name is required with -n."
[[ -n "$CACHE_INPUT" ]] || die "Cache directory is required with -o."

case "$CACHE_INPUT" in
    "~") CACHE_INPUT="$HOME" ;;
    "~/"*) CACHE_INPUT="$HOME/${CACHE_INPUT#~/}" ;;
esac

MACHINE_SAFE="$(printf '%s' "$MACHINE_NAME" | tr -cs 'A-Za-z0-9._-' '_')"
MACHINE_SAFE="${MACHINE_SAFE##_}"
MACHINE_SAFE="${MACHINE_SAFE%%_}"
[[ -n "$MACHINE_SAFE" ]] || die "Machine name does not contain usable characters."

MTD_DIR="$HOME/MTD"
CONDA_DIR="$HOME/miniconda3"
BENCHMARK_ROOT="$HOME/MTD_benchmarks"
SCRIPT_PATH="$(readlink -f "$0")"

if [[ "$SCRIPT_PATH" == "$MTD_DIR" || "$SCRIPT_PATH" == "$MTD_DIR/"* ]]; then
    die "Place this runner outside $MTD_DIR, for example $HOME/$PROGRAM_NAME"
fi

if [[ -e "$CACHE_INPUT" && ! -d "$CACHE_INPUT" ]]; then
    die "Cache path exists but is not a directory: $CACHE_INPUT"
fi

if [[ ! -e "$CACHE_INPUT" ]]; then
    printf '[INFO] Cache does not exist; creating a cold cache:\n  %s\n' "$CACHE_INPUT"
    mkdir -p -- "$CACHE_INPUT"
fi

CACHE="$(readlink -f -- "$CACHE_INPUT")"

case "$CACHE" in
    ""|"/"|"$HOME")
        die "Unsafe cache path: $CACHE"
        ;;
esac

for removable_path in "$MTD_DIR" "$CONDA_DIR" "$BENCHMARK_ROOT" "$HOME/R"; do
    if [[ "$CACHE" == "$removable_path" || "$CACHE" == "$removable_path/"* ]]; then
        die "Cache is inside a directory scheduled for removal: $removable_path"
    fi
done

if find "$CACHE" -type f -size +0c -print -quit 2>/dev/null | grep -q .; then
    CACHE_MODE="warm"
else
    CACHE_MODE="cold"
fi

BENCHMARK_LABEL="${MACHINE_SAFE}_${CACHE_MODE}_cache_auto_tos_r1"

STAGING_ROOT="$(mktemp -d "$HOME/.mtd_clone_staging.XXXXXX")"
STAGED_REPO="$STAGING_ROOT/MTD"

cleanup_staging() {
    if [[ -n "${STAGING_ROOT:-}" && -d "$STAGING_ROOT" ]]; then
        rm -rf -- "$STAGING_ROOT"
    fi
}
trap cleanup_staging EXIT

printf '%s\n' "============================================================"
printf '%s\n' "MTD CLEAN BENCHMARK PRECHECK"
printf '%-18s %s\n' "Machine:" "$MACHINE_NAME"
printf '%-18s %s\n' "Safe label:" "$MACHINE_SAFE"
printf '%-18s %s\n' "Cache:" "$CACHE"
printf '%-18s %s\n' "Cache mode:" "$CACHE_MODE"
printf '%-18s %s\n' "Benchmark label:" "$BENCHMARK_LABEL"
printf '%s\n' "============================================================"
du -sh "$CACHE" 2>/dev/null || true

echo
printf '[INFO] Cloning the current repository into a temporary staging area...\n'
git clone --depth 1 https://github.com/patrick-douglas/MTD.git "$STAGED_REPO"

INSTALLER="$STAGED_REPO/Install.sh"
[[ -s "$INSTALLER" ]] || die "Cloned Install.sh is missing or empty."
bash -n "$INSTALLER" || die "Cloned Install.sh failed Bash syntax validation."

# The final installer must embed automatic Conda ToS acceptance and must not
# contain the retired -a option or the old custom interactive prompt.
if ! grep -Eq '^[[:space:]]*export[[:space:]]+CONDA_PLUGINS_AUTO_ACCEPT_TOS=(yes|1|true)[[:space:]]*$' "$INSTALLER"; then
    die "The cloned Install.sh does not embed CONDA_PLUGINS_AUTO_ACCEPT_TOS=yes. Nothing was deleted."
fi

if grep -qE 'accept_required_conda_tos|accept_conda_tos|Accept the Anaconda Terms of Service\?' "$INSTALLER"; then
    die "The cloned Install.sh still contains the old interactive ToS implementation. Nothing was deleted."
fi

HELP_OUTPUT="$(bash "$INSTALLER" -h 2>&1 || true)"
if printf '%s\n' "$HELP_OUTPUT" | grep -Eq '(^|[[:space:]])-a([[:space:]]|$)'; then
    die "The cloned Install.sh still advertises the removed -a option. Nothing was deleted."
fi

for required_script in \
    MTD_benchmark_install.sh \
    MTD_make_instrumented_installer.sh \
    MTD_benchmark_merge.py \
    MTD_fix_profiled_locale.sh; do
    [[ -s "$STAGED_REPO/$required_script" ]] || die "Missing benchmark component: $required_script"
done

CLONED_COMMIT="$(git -C "$STAGED_REPO" rev-parse HEAD)"
printf '[PASS] Repository validated before cleanup.\n'
printf '[INFO] Commit: %s\n' "$CLONED_COMMIT"
printf '[PASS] Embedded Conda ToS autoaccept detected.\n'
printf '[PASS] Retired -a option and interactive ToS prompt are absent.\n'

echo
printf '%s\n' "The following installation state will be removed:"
for path in \
    "$MTD_DIR" \
    "$CONDA_DIR" \
    "$BENCHMARK_ROOT" \
    "$HOME/.conda" \
    "$HOME/.continuum" \
    "$HOME/.config/conda" \
    "$HOME/.cache/conda" \
    "$HOME/.cache/pip" \
    "$HOME/R"; do
    if [[ -e "$path" || -L "$path" ]]; then
        printf '  [REMOVE] %s\n' "$path"
    else
        printf '  [ABSENT] %s\n' "$path"
    fi
done

printf '\n[PRESERVE] %s\n' "$CACHE"
read -r -p "Type DELETE to clean and start the benchmark: " confirmation
[[ "$confirmation" == "DELETE" ]] || die "Cancelled. Nothing was removed."

conda deactivate >/dev/null 2>&1 || true
unset \
    CONDA_PREFIX \
    CONDA_DEFAULT_ENV \
    CONDA_EXE \
    CONDA_PYTHON_EXE \
    CONDA_SHLVL \
    _CE_CONDA \
    _CE_M \
    PYTHONPATH \
    R_LIBS \
    R_LIBS_USER \
    LD_LIBRARY_PATH || true
unset -f conda >/dev/null 2>&1 || true
hash -r

rm -rf -- \
    "$MTD_DIR" \
    "$CONDA_DIR" \
    "$BENCHMARK_ROOT" \
    "$HOME/.conda" \
    "$HOME/.continuum" \
    "$HOME/.config/conda" \
    "$HOME/.cache/conda" \
    "$HOME/.cache/pip" \
    "$HOME/R"

rm -f -- \
    "$HOME/.condarc" \
    "$HOME/MTD_git_commit_${MACHINE_SAFE}.txt" \
    "$HOME/MTD_preinstallation_state_${MACHINE_SAFE}.txt" \
    "$HOME/mtd_install_steps.tsv"

find "$HOME" -maxdepth 1 -type f \
    \( -name 'MTD_install_*.log' \
       -o -name 'MTD_clean_install_*.log' \
       -o -name 'mtd_install_*.log' \
       -o -name 'plasmid_library_preparation.log' \) \
    -print -delete 2>/dev/null || true

if [[ -f "$HOME/.bashrc" ]] && grep -q '# >>> conda initialize >>>' "$HOME/.bashrc"; then
    cp "$HOME/.bashrc" "$HOME/.bashrc.backup_before_mtd_$(date +%Y%m%d_%H%M%S)"
    sed -i '/# >>> conda initialize >>>/,/# <<< conda initialize <<</d' "$HOME/.bashrc"
fi

[[ -d "$CACHE" ]] || die "Cache disappeared during cleanup: $CACHE"

mv "$STAGED_REPO" "$MTD_DIR"
rmdir "$STAGING_ROOT" 2>/dev/null || true
STAGING_ROOT=""

cd "$MTD_DIR"
printf '%s\n' "$CLONED_COMMIT" | tee "$HOME/MTD_git_commit_${MACHINE_SAFE}.txt"

chmod +x \
    Install.sh \
    MTD_benchmark_install.sh \
    MTD_make_instrumented_installer.sh \
    MTD_benchmark_merge.py \
    MTD_fix_profiled_locale.sh

bash -n Install.sh
bash -n MTD_benchmark_install.sh
bash -n MTD_make_instrumented_installer.sh
bash -n MTD_fix_profiled_locale.sh
python3 -m py_compile MTD_benchmark_merge.py

rm -f Install_profiled.sh
bash ./MTD_make_instrumented_installer.sh \
    --input ./Install.sh \
    --output ./Install_profiled.sh

bash ./MTD_fix_profiled_locale.sh ./Install_profiled.sh
bash -n ./Install_profiled.sh

grep -q 'BEGIN MTD FUNCTION PROFILER' ./Install_profiled.sh || \
    die "Function profiler marker was not found in Install_profiled.sh."

grep -Eq '^[[:space:]]*export[[:space:]]+CONDA_PLUGINS_AUTO_ACCEPT_TOS=(yes|1|true)[[:space:]]*$' ./Install_profiled.sh || \
    die "The instrumented installer lost the embedded Conda ToS autoaccept."

if grep -qE 'accept_required_conda_tos|accept_conda_tos|Accept the Anaconda Terms of Service\?' ./Install_profiled.sh; then
    die "The instrumented installer contains the old interactive ToS implementation."
fi

{
    echo "Date: $(date --iso-8601=seconds)"
    echo "Machine label: $MACHINE_NAME"
    echo "Hostname: $(hostname)"
    echo "Cache mode: $CACHE_MODE"
    echo "Cache: $CACHE"
    echo "Git commit: $CLONED_COMMIT"
    echo "Install.sh SHA-256: $(sha256sum Install.sh | awk '{print $1}')"
    echo "Kernel: $(uname -a)"
    echo
    echo "================ CPU ================"
    lscpu
    echo
    echo "=============== MEMORY =============="
    free -h
    echo
    echo "=============== STORAGE ============="
    lsblk -o NAME,SIZE,TYPE,FSTYPE,ROTA,TRAN,MODEL,MOUNTPOINTS
    echo
    echo "============= FILESYSTEMS ==========="
    df -hT
    echo
    echo "================ CACHE ==============="
    du -sh "$CACHE"
} | tee "$HOME/MTD_preinstallation_state_${MACHINE_SAFE}.txt"

sudo -v

echo
printf '%s\n' "============================================================"
printf '%s\n' "STARTING MTD BENCHMARK"
printf '%-18s %s\n' "Machine:" "$MACHINE_NAME"
printf '%-18s %s\n' "Cache mode:" "$CACHE_MODE"
printf '%-18s %s\n' "Label:" "$BENCHMARK_LABEL"
printf '%-18s %s\n' "Cache:" "$CACHE"
printf '%-18s %s\n' "Output root:" "$BENCHMARK_ROOT"
printf '%s\n' "============================================================"

set +e
bash ./MTD_benchmark_install.sh \
    --label "$BENCHMARK_LABEL" \
    --interval 5 \
    --output-root "$BENCHMARK_ROOT" \
    --watch-path "$CONDA_DIR" \
    --watch-path "$CACHE" \
    --watch-path "$MTD_DIR" \
    -- \
    bash ./Install_profiled.sh \
        -o "$CACHE"
benchmark_status=$?
set -e

LATEST_BENCHMARK="$(
    find "$BENCHMARK_ROOT" \
        -mindepth 1 \
        -maxdepth 1 \
        -type d \
        -name "${BENCHMARK_LABEL}__*" \
        -printf '%T@ %p\n' 2>/dev/null |
    sort -nr |
    head -n 1 |
    cut -d' ' -f2-
)"

echo
printf '%s\n' "============================================================"
printf '%s\n' "BENCHMARK FINISHED"
printf 'Exit status: %s\n' "$benchmark_status"

if [[ -n "$LATEST_BENCHMARK" && -d "$LATEST_BENCHMARK" ]]; then
    printf 'Results: %s\n' "$LATEST_BENCHMARK"
    find "$LATEST_BENCHMARK" -maxdepth 1 -type f -printf '  %f\t%s bytes\n' | sort

    if [[ -s "$LATEST_BENCHMARK/summary.tsv" ]]; then
        echo
        echo "Summary:"
        if command -v column >/dev/null 2>&1; then
            column -t -s $'\t' "$LATEST_BENCHMARK/summary.tsv"
        else
            cat "$LATEST_BENCHMARK/summary.tsv"
        fi
    fi
else
    echo "[WARN] Benchmark result directory was not located."
fi

printf '[PASS] Cache preserved: %s\n' "$CACHE"
du -sh "$CACHE" 2>/dev/null || true

exit "$benchmark_status"
