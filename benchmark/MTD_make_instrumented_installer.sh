#!/usr/bin/env bash
# ==============================================================================
# Create an instrumented copy of the MTD installer
# ==============================================================================
# The original Install.sh is never edited. This script inserts a lightweight
# function profiler before the installation workflow starts and writes a new
# Install_profiled.sh.
#
# Usage:
#   bash benchmark/MTD_make_instrumented_installer.sh \
#       --input ./Install.sh \
#       --output ./Install_profiled.sh
#
# Run the generated installer through MTD_benchmark_install.sh. The profiler
# will automatically write steps.tsv into the benchmark run directory.
# ==============================================================================

set -euo pipefail

input="./Install.sh"
output="./Install_profiled.sh"

usage() {
    cat <<'USAGE'
Usage:
  MTD_make_instrumented_installer.sh [options]

Options:
  --input PATH     Original installer. Default: ./Install.sh
  --output PATH    Instrumented copy. Default: ./Install_profiled.sh
  -h, --help       Show this help.

The source installer is not modified.
USAGE
}

die() {
    printf '[ERROR] %s\n' "$*" >&2
    exit 1
}

while (($# > 0)); do
    case "$1" in
        --input)
            (($# >= 2)) || die "--input requires a path."
            input="$2"
            shift 2
            ;;
        --output)
            (($# >= 2)) || die "--output requires a path."
            output="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            die "Unknown option: $1"
            ;;
    esac
done

[[ -f "$input" ]] || die "Installer not found: $input"
[[ "$input" != "$output" ]] || die "Input and output must be different files."
grep -q '^#!/usr/bin/env bash' "$input" ||
    die "The input does not look like the expected Bash installer."
grep -q 'install_miniconda[[:space:]]*()' "$input" ||
    die "The input does not contain the expected install_miniconda() function."
grep -q 'prepare_installation_cache[[:space:]]*()' "$input" ||
    die "The input does not contain prepare_installation_cache()."

if grep -q 'BEGIN MTD FUNCTION PROFILER' "$input"; then
    die "The input is already instrumented."
fi

mkdir -p "$(dirname "$output")"

# The current MTD installer defines its functions and then starts the workflow
# with a bare init_colors call. Use the final exact call, never the function
# definition. A parse_arguments fallback is retained for future reorganizations.
insertion_line="$(
    grep -nE '^[[:space:]]*init_colors[[:space:]]*$' "$input" |
    tail -n 1 |
    cut -d: -f1 || true
)"

if [[ -z "$insertion_line" ]]; then
    insertion_line="$(
        grep -nE '^[[:space:]]*parse_arguments[[:space:]]+"\$@"[[:space:]]*$' "$input" |
        tail -n 1 |
        cut -d: -f1 || true
    )"
fi

[[ -n "$insertion_line" ]] ||
    die "Could not identify where the installation workflow starts."

tmp="$(mktemp "${TMPDIR:-/tmp}/mtd-profiled-installer.XXXXXX")"
trap 'rm -f "$tmp"' EXIT

head -n "$((insertion_line - 1))" "$input" > "$tmp"

cat >> "$tmp" <<'PROFILER'
# ==============================================================================
# BEGIN MTD FUNCTION PROFILER
# Added automatically by MTD_make_instrumented_installer.sh.
# Remove this block by regenerating the file from the original Install.sh.
# ==============================================================================

MTD_BENCHMARK_STEPS_TSV="${MTD_BENCHMARK_STEPS_TSV:-${PWD}/mtd_install_steps.tsv}"
MTD_BENCHMARK_RUN_ID="${MTD_BENCHMARK_RUN_ID:-manual_$(date -u '+%Y%m%dT%H%M%SZ')}"
MTD_BENCHMARK_LABEL="${MTD_BENCHMARK_LABEL:-manual_profile}"
MTD_BENCHMARK_HOST="${MTD_BENCHMARK_HOST:-$(hostname -s 2>/dev/null || hostname)}"
MTD_BENCHMARK_MIN_SECONDS="${MTD_BENCHMARK_MIN_SECONDS:-0.05}"
MTD_BENCHMARK_APPEND="${MTD_BENCHMARK_APPEND:-0}"
MTD_BENCH_STACK=()

mtd_bench_epoch_us() {
    local realtime seconds fraction

    if [[ -n "${EPOCHREALTIME:-}" ]]; then
        # Bash may format EPOCHREALTIME with a comma under locales such as
        # pt_BR. Normalize it before splitting seconds and microseconds.
        realtime="${EPOCHREALTIME/,/.}"
        seconds="${realtime%%.*}"
        fraction="${realtime#*.}000000"
        fraction="${fraction:0:6}"
        printf '%s%s\n' "$seconds" "$fraction"
    else
        # GNU date fallback for older Bash versions.
        date +%s%6N
    fi
}

mtd_bench_seconds_to_us() {
    local value="$1"
    local whole fraction

    # Accept either decimal dot or decimal comma in user-supplied values.
    value="${value/,/.}"
    whole="${value%%.*}"
    if [[ "$value" == *.* ]]; then
        fraction="${value#*.}000000"
        fraction="${fraction:0:6}"
    else
        fraction="000000"
    fi

    printf '%s\n' "$((10#$whole * 1000000 + 10#$fraction))"
}

mtd_bench_format_elapsed() {
    local elapsed_us="$1"
    printf '%d.%06d' \
        "$((elapsed_us / 1000000))" \
        "$((elapsed_us % 1000000))"
}

mtd_bench_epoch_us_to_utc() {
    local epoch_us="$1"
    date -u -d "@$((epoch_us / 1000000))" '+%Y-%m-%dT%H:%M:%SZ'
}

mtd_bench_initialize() {
    local parent_dir
    parent_dir="$(dirname "$MTD_BENCHMARK_STEPS_TSV")"
    mkdir -p "$parent_dir"

    MTD_BENCHMARK_MIN_US="$(
        mtd_bench_seconds_to_us "$MTD_BENCHMARK_MIN_SECONDS"
    )"

    if [[ "$MTD_BENCHMARK_APPEND" != "1" || ! -s "$MTD_BENCHMARK_STEPS_TSV" ]]; then
        printf '%s\n' \
            $'run_id\tlabel\thost\tfunction\tparent_function\tcall_depth\tstart_utc\tend_utc\tstart_epoch_us\tend_epoch_us\telapsed_seconds\tstatus\tpid' \
            > "$MTD_BENCHMARK_STEPS_TSV"
    fi

    printf '[BENCHMARK] Function-level timing enabled.\n'
    printf '[BENCHMARK] Step metrics: %s\n' "$MTD_BENCHMARK_STEPS_TSV"
    printf '[BENCHMARK] Minimum recorded duration: %s seconds\n' \
        "$MTD_BENCHMARK_MIN_SECONDS"
}

mtd_bench_should_record() {
    local elapsed_us="$1"
    ((elapsed_us >= MTD_BENCHMARK_MIN_US))
}

mtd_bench_append_row() {
    local function_name="$1"
    local parent_function="$2"
    local depth="$3"
    local start_utc="$4"
    local end_utc="$5"
    local start_epoch_us="$6"
    local end_epoch_us="$7"
    local elapsed="$8"
    local status="$9"

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$MTD_BENCHMARK_RUN_ID" \
        "$MTD_BENCHMARK_LABEL" \
        "$MTD_BENCHMARK_HOST" \
        "$function_name" \
        "$parent_function" \
        "$depth" \
        "$start_utc" \
        "$end_utc" \
        "$start_epoch_us" \
        "$end_epoch_us" \
        "$elapsed" \
        "$status" \
        "$$" \
        >> "$MTD_BENCHMARK_STEPS_TSV"
}

mtd_bench_run() {
    local function_name="$1"
    local original_function="$2"
    shift 2

    local start_us end_us elapsed_us start_utc end_utc elapsed status
    local depth parent_function stack_index

    MTD_BENCH_STACK+=("$function_name")
    depth="${#MTD_BENCH_STACK[@]}"
    parent_function="TOP_LEVEL"
    if ((depth > 1)); then
        parent_function="${MTD_BENCH_STACK[$((depth - 2))]}"
    fi

    start_us="$(mtd_bench_epoch_us)"

    "$original_function" "$@"
    status=$?

    end_us="$(mtd_bench_epoch_us)"
    elapsed_us=$((end_us - start_us))

    if mtd_bench_should_record "$elapsed_us"; then
        start_utc="$(mtd_bench_epoch_us_to_utc "$start_us")"
        end_utc="$(mtd_bench_epoch_us_to_utc "$end_us")"
        elapsed="$(mtd_bench_format_elapsed "$elapsed_us")"

        mtd_bench_append_row \
            "$function_name" \
            "$parent_function" \
            "$depth" \
            "$start_utc" \
            "$end_utc" \
            "$start_us" \
            "$end_us" \
            "$elapsed" \
            "$status"
    fi

    stack_index=$((depth - 1))
    unset 'MTD_BENCH_STACK[stack_index]'
    return "$status"
}

mtd_bench_wrap_function() {
    local function_name="$1"
    local original_name definition

    case "$function_name" in
        mtd_bench_*|\
        log_info|log_ok|log_warning|log_error|\
        print_rule|show_progress|\
        on_exit|on_interrupt|restore_conda_channel_priority)
            return 0
            ;;
    esac

    declare -F "$function_name" >/dev/null 2>&1 || return 0

    original_name="__mtd_bench_original_${function_name}"
    declare -F "$original_name" >/dev/null 2>&1 && return 0

    definition="$(declare -f "$function_name")"
    definition="$(
        printf '%s\n' "$definition" |
        sed -E "1s/^${function_name}[[:space:]]*\\(\\)/${original_name} ()/"
    )"

    eval "$definition"
    eval "${function_name}() { mtd_bench_run '${function_name}' '${original_name}' \"\$@\"; }"
}

mtd_bench_initialize

# Snapshot function names before wrapping so the newly created original copies
# are not themselves wrapped.
mapfile -t __mtd_functions_to_wrap < <(
    declare -F |
    awk '{print $3}' |
    LC_ALL=C sort
)

for __mtd_function_name in "${__mtd_functions_to_wrap[@]}"; do
    mtd_bench_wrap_function "$__mtd_function_name"
done

unset __mtd_function_name
unset __mtd_functions_to_wrap

# ==============================================================================
# END MTD FUNCTION PROFILER
# ==============================================================================

PROFILER

tail -n "+$insertion_line" "$input" >> "$tmp"

bash -n "$tmp" || die "The generated installer did not pass bash -n."

cp -f "$tmp" "$output"
chmod --reference="$input" "$output" 2>/dev/null || chmod +x "$output"

printf '[OK] Instrumented installer created:\n'
printf '  %s\n' "$(readlink -f "$output")"
printf '[OK] Original installer was not modified:\n'
printf '  %s\n' "$(readlink -f "$input")"
printf '\nRun it through the system benchmark wrapper, for example:\n\n'
printf '  bash ./benchmark/MTD_benchmark_install.sh \\\n'
printf '    --label master_cold_native_r1 \\\n'
printf '    --watch-path "$HOME/miniconda3" \\\n'
printf '    --watch-path "/path/to/cache" \\\n'
printf '    -- bash ./Install_profiled.sh -o "/path/to/cache"\n'
