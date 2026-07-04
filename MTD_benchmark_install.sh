#!/usr/bin/env bash
# ==============================================================================
# MTD installation benchmark wrapper
# ==============================================================================
# Runs any installation command in the foreground while collecting:
#   - complete terminal log;
#   - GNU time statistics, when /usr/bin/time is available;
#   - CPU, RAM, process-tree RSS, disk I/O, network I/O and temperature samples;
#   - operating-system, CPU, memory, storage and network hardware information;
#   - Git commit and SHA-256 hashes;
#   - watched-directory sizes and filesystem space before and after installation.
#
# The monitored command keeps access to the terminal, so sudo prompts and other
# interactive confirmations continue to work.
#
# Example:
#   bash MTD_benchmark_install.sh \
#     --label master_cold_allthreads \
#     --interval 5 \
#     --output-root "$HOME/MTD_benchmarks" \
#     --watch-path "$HOME/miniconda3" \
#     --watch-path "/data/MTD_install_cache" \
#     -- bash ./Install.sh -o /data/MTD_install_cache
# ==============================================================================

set -uo pipefail

PROGRAM_NAME="$(basename "$0")"
LABEL=""
INTERVAL=5
OUTPUT_ROOT="${PWD}/MTD_benchmarks"
WATCH_PATHS=()
COMMAND=()

usage() {
    cat <<USAGE
Usage:
  $PROGRAM_NAME [benchmark options] -- <installation command> [arguments...]

Benchmark options:
  --label TEXT          Required run label, e.g. master_cold_allthreads
  --interval SECONDS    Sampling interval; default: 5
  --output-root PATH    Parent output directory; default: ./MTD_benchmarks
  --watch-path PATH     Directory to measure before/after; may be repeated
  -h, --help            Show this help

Example:
  bash $PROGRAM_NAME \\
    --label master_cold_allthreads \\
    --interval 5 \\
    --output-root "\$HOME/MTD_benchmarks" \\
    --watch-path "\$HOME/miniconda3" \\
    --watch-path "/data/MTD_install_cache" \\
    -- bash ./Install.sh -o /data/MTD_install_cache
USAGE
}

die() {
    printf '[ERROR] %s\n' "$*" >&2
    exit 1
}

while (($# > 0)); do
    case "$1" in
        --label)
            (($# >= 2)) || die "--label requires a value."
            LABEL="$2"
            shift 2
            ;;
        --interval)
            (($# >= 2)) || die "--interval requires a value."
            INTERVAL="$2"
            shift 2
            ;;
        --output-root)
            (($# >= 2)) || die "--output-root requires a value."
            OUTPUT_ROOT="$2"
            shift 2
            ;;
        --watch-path)
            (($# >= 2)) || die "--watch-path requires a value."
            WATCH_PATHS+=("$2")
            shift 2
            ;;
        --)
            shift
            COMMAND=("$@")
            break
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            die "Unknown benchmark option: $1"
            ;;
    esac
done

[[ -n "$LABEL" ]] || die "--label is required."
[[ "$INTERVAL" =~ ^[0-9]+([.][0-9]+)?$ ]] || die "--interval must be numeric."
awk -v x="$INTERVAL" 'BEGIN { exit !(x > 0) }' || die "--interval must be greater than zero."
((${#COMMAND[@]} > 0)) || die "No installation command was provided after --."
command -v python3 >/dev/null 2>&1 || die "python3 is required for monitoring."

safe_label="$(printf '%s' "$LABEL" | tr -cs 'A-Za-z0-9._-' '_')"
timestamp="$(date -u '+%Y%m%dT%H%M%SZ')"
host_short="$(hostname -s 2>/dev/null || hostname)"
run_id="${safe_label}__${host_short}__${timestamp}"
output_dir="${OUTPUT_ROOT%/}/${run_id}"

mkdir -p "$output_dir" || die "Could not create output directory: $output_dir"
output_dir="$(readlink -f "$output_dir")"

console_log="$output_dir/console.log"
time_log="$output_dir/gnu_time.txt"
samples_csv="$output_dir/resource_samples.csv"
summary_tsv="$output_dir/summary.tsv"
metadata_file="$output_dir/metadata.txt"
hardware_file="$output_dir/hardware.txt"
software_file="$output_dir/software.txt"
git_file="$output_dir/git_state.txt"
watch_before="$output_dir/watch_paths_before.tsv"
watch_after="$output_dir/watch_paths_after.tsv"
stop_file="$output_dir/.stop_monitor"
monitor_file="$output_dir/.mtd_resource_monitor.py"

rm -f "$stop_file"

# Expose the run directory to an instrumented installer. These variables do
# not change a normal Install.sh; they are used only by Install_profiled.sh.
export MTD_BENCHMARK_RUN_DIR="$output_dir"
export MTD_BENCHMARK_RUN_ID="$run_id"
export MTD_BENCHMARK_LABEL="$LABEL"
export MTD_BENCHMARK_HOST="$host_short"
export MTD_BENCHMARK_STEPS_TSV="${MTD_BENCHMARK_STEPS_TSV:-$output_dir/steps.tsv}"

# Preserve interactive stdin while logging stdout and stderr.
exec > >(tee -a "$console_log") 2>&1

quote_command() {
    printf '%q ' "${COMMAND[@]}"
    printf '\n'
}

nearest_existing_path() {
    local candidate="$1"
    if [[ "$candidate" != /* ]]; then
        candidate="$PWD/$candidate"
    fi
    while [[ ! -e "$candidate" && "$candidate" != "/" ]]; do
        candidate="$(dirname "$candidate")"
    done
    printf '%s\n' "$candidate"
}

path_snapshot() {
    local destination="$1"
    local phase="$2"
    local requested existing exists size_bytes source fstype avail_bytes total_bytes

    printf 'phase\trequested_path\texists\tsize_bytes\tfilesystem_source\tfilesystem_type\tfilesystem_total_bytes\tfilesystem_available_bytes\n' > "$destination"

    for requested in "${WATCH_PATHS[@]}"; do
        if [[ -e "$requested" ]]; then
            exists=1
            size_bytes="$(du -sx --block-size=1 "$requested" 2>/dev/null | awk '{print $1}' || true)"
        else
            exists=0
            size_bytes=0
        fi

        existing="$(nearest_existing_path "$requested")"
        source="$(findmnt -n -o SOURCE -T "$existing" 2>/dev/null || printf 'unknown')"
        fstype="$(findmnt -n -o FSTYPE -T "$existing" 2>/dev/null || printf 'unknown')"
        total_bytes="$(df -PB1 "$existing" 2>/dev/null | awk 'NR==2 {print $2}' || true)"
        avail_bytes="$(df -PB1 "$existing" 2>/dev/null | awk 'NR==2 {print $4}' || true)"

        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$phase" "$requested" "$exists" "${size_bytes:-NA}" \
            "$source" "$fstype" "${total_bytes:-NA}" "${avail_bytes:-NA}" \
            >> "$destination"
    done
}

collect_hardware() {
    {
        echo "=== Host identity ==="
        printf 'hostname: '; hostname 2>/dev/null || true
        printf 'hostnamectl: '; hostnamectl 2>/dev/null | head -n 8 || true
        printf 'kernel: '; uname -a
        printf 'architecture: '; uname -m
        printf 'logical_cpus: '; nproc
        echo

        echo "=== Operating system ==="
        cat /etc/os-release 2>/dev/null || true
        echo

        echo "=== CPU ==="
        lscpu 2>/dev/null || cat /proc/cpuinfo
        echo

        echo "=== Memory ==="
        free -b 2>/dev/null || true
        echo
        cat /proc/meminfo 2>/dev/null || true
        echo

        echo "=== DMI / machine ==="
        for f in \
            /sys/devices/virtual/dmi/id/sys_vendor \
            /sys/devices/virtual/dmi/id/product_name \
            /sys/devices/virtual/dmi/id/product_version \
            /sys/devices/virtual/dmi/id/board_vendor \
            /sys/devices/virtual/dmi/id/board_name
        do
            [[ -r "$f" ]] && printf '%s: %s\n' "$(basename "$f")" "$(cat "$f")"
        done
        echo

        echo "=== Block devices ==="
        lsblk -e 7 -b -o NAME,KNAME,TYPE,SIZE,ROTA,TRAN,MODEL,FSTYPE,MOUNTPOINTS 2>/dev/null || lsblk
        echo

        echo "=== Mounted filesystems ==="
        findmnt -o TARGET,SOURCE,FSTYPE,OPTIONS 2>/dev/null || true
        echo

        echo "=== Filesystem capacity ==="
        df -B1 -T
        echo

        echo "=== Network interfaces ==="
        ip -brief link 2>/dev/null || true
        ip -brief address 2>/dev/null || true
        echo

        echo "=== PCI storage/network controllers ==="
        if command -v lspci >/dev/null 2>&1; then
            lspci | grep -Ei 'ethernet|network|sata|raid|non-volatile|nvme|scsi|storage' || true
        fi
        echo

        echo "=== Thermal zones at baseline ==="
        for f in /sys/class/thermal/thermal_zone*/temp; do
            [[ -r "$f" ]] || continue
            type_file="${f%/temp}/type"
            printf '%s\t%s\t%s\n' \
                "$(basename "${f%/temp}")" \
                "$([[ -r "$type_file" ]] && cat "$type_file" || printf unknown)" \
                "$(cat "$f")"
        done
    } > "$hardware_file"
}

collect_software() {
    {
        echo "=== Benchmark wrapper ==="
        printf 'wrapper_path: %s\n' "$(readlink -f "$0")"
        sha256sum "$0" 2>/dev/null || true
        echo

        echo "=== Command ==="
        quote_command
        echo

        echo "=== Tool versions available before installation ==="
        for cmd in bash python3 git curl wget time awk sed grep tar gzip xz; do
            if command -v "$cmd" >/dev/null 2>&1; then
                printf '\n--- %s ---\n' "$cmd"
                case "$cmd" in
                    bash) bash --version | head -n 1 ;;
                    python3) python3 --version ;;
                    git) git --version ;;
                    curl) curl --version | head -n 1 ;;
                    wget) wget --version | head -n 1 ;;
                    time) /usr/bin/time --version 2>&1 | head -n 1 ;;
                    awk) awk --version 2>&1 | head -n 1 || true ;;
                    sed) sed --version 2>&1 | head -n 1 || true ;;
                    grep) grep --version 2>&1 | head -n 1 || true ;;
                    tar) tar --version 2>&1 | head -n 1 || true ;;
                    gzip) gzip --version 2>&1 | head -n 1 || true ;;
                    xz) xz --version 2>&1 | head -n 1 || true ;;
                esac
            fi
        done
    } > "$software_file"
}

collect_git_state() {
    {
        echo "=== Current directory ==="
        pwd
        echo

        echo "=== Git repository ==="
        if git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
            printf 'root: '; git rev-parse --show-toplevel
            printf 'commit: '; git rev-parse HEAD
            printf 'branch: '; git branch --show-current
            printf 'describe: '; git describe --always --dirty --tags 2>/dev/null || true
            echo
            echo "--- remotes ---"
            git remote -v || true
            echo
            echo "--- status --short ---"
            git status --short || true
        else
            echo "Not executed inside a Git work tree."
        fi
        echo

        echo "=== Existing command files and SHA-256 ==="
        for arg in "${COMMAND[@]}"; do
            if [[ -f "$arg" ]]; then
                sha256sum "$arg" || true
            fi
        done
    } > "$git_file"
}

cat > "$monitor_file" <<'PYMONITOR'
#!/usr/bin/env python3
import argparse
import csv
import datetime as dt
import os
import pathlib
import time

parser = argparse.ArgumentParser()
parser.add_argument("--root-pid", type=int, required=True)
parser.add_argument("--interval", type=float, required=True)
parser.add_argument("--output", required=True)
parser.add_argument("--stop-file", required=True)
args = parser.parse_args()

ROOT_PID = args.root_pid
SELF_PID = os.getpid()
CLK_TCK = os.sysconf(os.sysconf_names["SC_CLK_TCK"])
NCPU = os.cpu_count() or 1
SECTOR_BYTES = 512

def utc_now():
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")

def read_cpu_total():
    with open("/proc/stat", "r", encoding="utf-8") as handle:
        parts = handle.readline().split()
    values = [int(x) for x in parts[1:]]
    idle = values[3] + (values[4] if len(values) > 4 else 0)
    total = sum(values)
    return total, idle

def read_meminfo():
    result = {}
    with open("/proc/meminfo", "r", encoding="utf-8") as handle:
        for line in handle:
            key, value = line.split(":", 1)
            result[key] = int(value.strip().split()[0])  # kB
    return result

def read_load1():
    with open("/proc/loadavg", "r", encoding="utf-8") as handle:
        return float(handle.read().split()[0])

def read_proc_table():
    table = {}
    for entry in pathlib.Path("/proc").iterdir():
        if not entry.name.isdigit():
            continue
        pid = int(entry.name)
        try:
            raw = (entry / "stat").read_text(encoding="utf-8")
            right = raw.rfind(")")
            fields = raw[right + 2:].split()
            ppid = int(fields[1])
            utime = int(fields[11])
            stime = int(fields[12])

            rss_kb = 0
            swap_kb = 0
            threads = 0
            for line in (entry / "status").read_text(encoding="utf-8").splitlines():
                if line.startswith("VmRSS:"):
                    rss_kb = int(line.split()[1])
                elif line.startswith("VmSwap:"):
                    swap_kb = int(line.split()[1])
                elif line.startswith("Threads:"):
                    threads = int(line.split()[1])
            table[pid] = {
                "ppid": ppid,
                "ticks": utime + stime,
                "rss_kb": rss_kb,
                "swap_kb": swap_kb,
                "threads": threads,
            }
        except (FileNotFoundError, ProcessLookupError, PermissionError, ValueError, IndexError):
            continue
    return table

def descendants(table):
    children = {}
    for pid, info in table.items():
        children.setdefault(info["ppid"], []).append(pid)

    selected = set()
    stack = list(children.get(ROOT_PID, []))
    while stack:
        pid = stack.pop()
        if pid == SELF_PID or pid in selected:
            continue
        selected.add(pid)
        stack.extend(children.get(pid, []))
    return selected

def physical_disks():
    names = []
    sys_block = pathlib.Path("/sys/block")
    if not sys_block.exists():
        return names
    for entry in sys_block.iterdir():
        name = entry.name
        if name.startswith(("loop", "ram", "zram", "fd")):
            continue
        names.append(name)
    return names

DISKS = physical_disks()

def read_disk_bytes():
    sectors_read = 0
    sectors_written = 0
    for name in DISKS:
        try:
            fields = pathlib.Path(f"/sys/block/{name}/stat").read_text().split()
            sectors_read += int(fields[2])
            sectors_written += int(fields[6])
        except (FileNotFoundError, ValueError, IndexError, PermissionError):
            continue
    return sectors_read * SECTOR_BYTES, sectors_written * SECTOR_BYTES

def read_network_bytes():
    rx = 0
    tx = 0
    with open("/proc/net/dev", "r", encoding="utf-8") as handle:
        for line in handle.readlines()[2:]:
            iface, data = line.split(":", 1)
            if iface.strip() == "lo":
                continue
            fields = data.split()
            rx += int(fields[0])
            tx += int(fields[8])
    return rx, tx

def read_max_temp_c():
    values = []
    for path in pathlib.Path("/sys/class/thermal").glob("thermal_zone*/temp"):
        try:
            raw = float(path.read_text().strip())
            values.append(raw / 1000.0 if raw > 1000 else raw)
        except (OSError, ValueError):
            continue
    return max(values) if values else None

start_monotonic = time.monotonic()
prev_sample_time = start_monotonic
prev_cpu_total, prev_cpu_idle = read_cpu_total()
base_disk_read, base_disk_write = read_disk_bytes()
base_net_rx, base_net_tx = read_network_bytes()
previous_ticks = {}

fieldnames = [
    "timestamp_utc",
    "elapsed_seconds",
    "system_cpu_busy_percent",
    "process_tree_cpu_percent_one_core",
    "process_tree_cpu_percent_total_capacity",
    "system_memory_used_gib",
    "system_memory_available_gib",
    "process_tree_rss_gib",
    "process_tree_swap_gib",
    "load_1min",
    "process_count",
    "thread_count",
    "disk_read_gib_since_start",
    "disk_write_gib_since_start",
    "network_rx_gib_since_start",
    "network_tx_gib_since_start",
    "max_temperature_c",
]

with open(args.output, "w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=fieldnames)
    writer.writeheader()

    while True:
        now = time.monotonic()
        interval = max(now - prev_sample_time, 1e-9)

        cpu_total, cpu_idle = read_cpu_total()
        delta_total = max(cpu_total - prev_cpu_total, 1)
        delta_idle = max(cpu_idle - prev_cpu_idle, 0)
        system_cpu = 100.0 * (delta_total - delta_idle) / delta_total

        table = read_proc_table()
        selected = descendants(table)

        rss_kb = sum(table[pid]["rss_kb"] for pid in selected if pid in table)
        swap_kb = sum(table[pid]["swap_kb"] for pid in selected if pid in table)
        threads = sum(table[pid]["threads"] for pid in selected if pid in table)

        delta_ticks = 0
        current_ticks = {}
        for pid in selected:
            if pid not in table:
                continue
            ticks = table[pid]["ticks"]
            # A process first observed in this sample starts with its current
            # tick count as baseline; this prevents pre-existing helper processes
            # from creating an artificial CPU spike.
            old_ticks = previous_ticks.get(pid, ticks)
            if ticks >= old_ticks:
                delta_ticks += ticks - old_ticks
            current_ticks[pid] = ticks

        proc_cpu_one_core = 100.0 * (delta_ticks / CLK_TCK) / interval
        proc_cpu_total = proc_cpu_one_core / NCPU

        mem = read_meminfo()
        mem_total_kb = mem.get("MemTotal", 0)
        mem_available_kb = mem.get("MemAvailable", mem.get("MemFree", 0))
        mem_used_kb = max(mem_total_kb - mem_available_kb, 0)

        disk_read, disk_write = read_disk_bytes()
        net_rx, net_tx = read_network_bytes()
        temperature = read_max_temp_c()

        writer.writerow({
            "timestamp_utc": utc_now(),
            "elapsed_seconds": f"{now - start_monotonic:.3f}",
            "system_cpu_busy_percent": f"{system_cpu:.3f}",
            "process_tree_cpu_percent_one_core": f"{proc_cpu_one_core:.3f}",
            "process_tree_cpu_percent_total_capacity": f"{proc_cpu_total:.3f}",
            "system_memory_used_gib": f"{mem_used_kb / 1024 / 1024:.6f}",
            "system_memory_available_gib": f"{mem_available_kb / 1024 / 1024:.6f}",
            "process_tree_rss_gib": f"{rss_kb / 1024 / 1024:.6f}",
            "process_tree_swap_gib": f"{swap_kb / 1024 / 1024:.6f}",
            "load_1min": f"{read_load1():.3f}",
            "process_count": len(selected),
            "thread_count": threads,
            "disk_read_gib_since_start": f"{max(disk_read - base_disk_read, 0) / 1024**3:.6f}",
            "disk_write_gib_since_start": f"{max(disk_write - base_disk_write, 0) / 1024**3:.6f}",
            "network_rx_gib_since_start": f"{max(net_rx - base_net_rx, 0) / 1024**3:.6f}",
            "network_tx_gib_since_start": f"{max(net_tx - base_net_tx, 0) / 1024**3:.6f}",
            "max_temperature_c": "" if temperature is None else f"{temperature:.3f}",
        })
        handle.flush()

        previous_ticks = current_ticks
        prev_cpu_total, prev_cpu_idle = cpu_total, cpu_idle
        prev_sample_time = now

        if os.path.exists(args.stop_file):
            break
        time.sleep(args.interval)
PYMONITOR
chmod 700 "$monitor_file"

collect_hardware
collect_software
collect_git_state
path_snapshot "$watch_before" "before"

start_utc="$(date -u '+%Y-%m-%dT%H:%M:%SZ')"
start_epoch_ns="$(date +%s%N)"

{
    printf 'run_id\t%s\n' "$run_id"
    printf 'label\t%s\n' "$LABEL"
    printf 'host\t%s\n' "$host_short"
    printf 'start_utc\t%s\n' "$start_utc"
    printf 'sampling_interval_seconds\t%s\n' "$INTERVAL"
    printf 'output_directory\t%s\n' "$output_dir"
    printf 'working_directory\t%s\n' "$PWD"
    printf 'command\t'
    quote_command
} > "$metadata_file"

printf '\n============================================================\n'
printf '[BENCHMARK] Run ID: %s\n' "$run_id"
printf '[BENCHMARK] Output: %s\n' "$output_dir"
printf '[BENCHMARK] Command: '
quote_command
printf '============================================================\n\n'

python3 "$monitor_file" \
    --root-pid "$$" \
    --interval "$INTERVAL" \
    --output "$samples_csv" \
    --stop-file "$stop_file" &
monitor_pid=$!

cleanup_monitor() {
    touch "$stop_file" 2>/dev/null || true
    if [[ -n "${monitor_pid:-}" ]]; then
        wait "$monitor_pid" 2>/dev/null || true
    fi
}
trap cleanup_monitor EXIT

status=0
if [[ -x /usr/bin/time ]]; then
    /usr/bin/time -v -o "$time_log" -- "${COMMAND[@]}"
    status=$?
else
    printf '[WARNING] /usr/bin/time was not found; GNU time metrics will be absent.\n'
    "${COMMAND[@]}"
    status=$?
fi

end_epoch_ns="$(date +%s%N)"
end_utc="$(date -u '+%Y-%m-%dT%H:%M:%SZ')"

touch "$stop_file"
wait "$monitor_pid" 2>/dev/null || true
monitor_pid=""

path_snapshot "$watch_after" "after"

python3 - "$samples_csv" "$summary_tsv" "$time_log" "$watch_before" "$watch_after" \
    "$run_id" "$LABEL" "$host_short" "$start_utc" "$end_utc" \
    "$start_epoch_ns" "$end_epoch_ns" "$status" <<'PYSUMMARY'
import csv
import pathlib
import re
import sys

(
    samples_path,
    summary_path,
    time_path,
    watch_before_path,
    watch_after_path,
    run_id,
    label,
    host,
    start_utc,
    end_utc,
    start_ns,
    end_ns,
    exit_status,
) = sys.argv[1:]

rows = []
sample_file = pathlib.Path(samples_path)
if sample_file.exists():
    with sample_file.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

def floats(column):
    out = []
    for row in rows:
        value = row.get(column, "")
        if value not in ("", None):
            try:
                out.append(float(value))
            except ValueError:
                pass
    return out

def maximum(column):
    values = floats(column)
    return max(values) if values else None

def mean(column):
    values = floats(column)
    return sum(values) / len(values) if values else None

def last(column):
    values = floats(column)
    return values[-1] if values else None

gnu = {}
time_file = pathlib.Path(time_path)
if time_file.exists():
    for line in time_file.read_text(encoding="utf-8", errors="replace").splitlines():
        if ":" not in line:
            continue
        key, value = line.split(":", 1)
        gnu[key.strip()] = value.strip()

def parse_float(text):
    if text is None:
        return None
    match = re.search(r"-?[0-9]+(?:\.[0-9]+)?", text)
    return float(match.group(0)) if match else None

def watch_sizes(path):
    result = {}
    p = pathlib.Path(path)
    if not p.exists():
        return result
    with p.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            result[row["requested_path"]] = row
    return result

before = watch_sizes(watch_before_path)
after = watch_sizes(watch_after_path)

wall_seconds = (int(end_ns) - int(start_ns)) / 1_000_000_000

metrics = [
    ("run_id", run_id, "text"),
    ("label", label, "text"),
    ("host", host, "text"),
    ("start_utc", start_utc, "ISO-8601 UTC"),
    ("end_utc", end_utc, "ISO-8601 UTC"),
    ("exit_status", exit_status, "integer"),
    ("wall_time_seconds", f"{wall_seconds:.6f}", "seconds"),
    ("sample_count", str(len(rows)), "count"),
    ("mean_system_cpu_busy", mean("system_cpu_busy_percent"), "percent of total CPU capacity"),
    ("peak_system_cpu_busy", maximum("system_cpu_busy_percent"), "percent of total CPU capacity"),
    ("mean_process_tree_cpu_total_capacity", mean("process_tree_cpu_percent_total_capacity"), "percent of total CPU capacity"),
    ("peak_process_tree_cpu_total_capacity", maximum("process_tree_cpu_percent_total_capacity"), "percent of total CPU capacity"),
    ("peak_process_tree_cpu_one_core", maximum("process_tree_cpu_percent_one_core"), "percent of one logical CPU"),
    ("peak_process_tree_rss", maximum("process_tree_rss_gib"), "GiB; sampled sum across active descendants"),
    ("peak_process_tree_swap", maximum("process_tree_swap_gib"), "GiB; sampled sum across active descendants"),
    ("peak_system_memory_used", maximum("system_memory_used_gib"), "GiB; MemTotal minus MemAvailable"),
    ("minimum_system_memory_available", min(floats("system_memory_available_gib")) if floats("system_memory_available_gib") else None, "GiB"),
    ("peak_load_1min", maximum("load_1min"), "load average"),
    ("peak_process_count", maximum("process_count"), "count"),
    ("peak_thread_count", maximum("thread_count"), "count"),
    ("total_physical_disk_read", last("disk_read_gib_since_start"), "GiB; system-wide physical disks"),
    ("total_physical_disk_write", last("disk_write_gib_since_start"), "GiB; system-wide physical disks"),
    ("total_network_received", last("network_rx_gib_since_start"), "GiB; system-wide excluding loopback"),
    ("total_network_transmitted", last("network_tx_gib_since_start"), "GiB; system-wide excluding loopback"),
    ("peak_reported_temperature", maximum("max_temperature_c"), "degrees Celsius"),
    ("gnu_time_user_seconds", parse_float(gnu.get("User time (seconds)")), "seconds"),
    ("gnu_time_system_seconds", parse_float(gnu.get("System time (seconds)")), "seconds"),
    ("gnu_time_cpu_percent", parse_float(gnu.get("Percent of CPU this job got")), "percent"),
    ("gnu_time_max_rss", parse_float(gnu.get("Maximum resident set size (kbytes)")), "KiB; GNU time definition"),
    ("gnu_time_major_page_faults", parse_float(gnu.get("Major (requiring I/O) page faults")), "count"),
    ("gnu_time_minor_page_faults", parse_float(gnu.get("Minor (reclaiming a frame) page faults")), "count"),
    ("gnu_time_fs_inputs", parse_float(gnu.get("File system inputs")), "implementation-defined blocks"),
    ("gnu_time_fs_outputs", parse_float(gnu.get("File system outputs")), "implementation-defined blocks"),
]

for path in sorted(set(before) | set(after)):
    b = before.get(path, {})
    a = after.get(path, {})
    bsize = b.get("size_bytes")
    asize = a.get("size_bytes")
    metrics.extend([
        (f"watch_before_size_bytes::{path}", bsize, "bytes"),
        (f"watch_after_size_bytes::{path}", asize, "bytes"),
        (f"watch_filesystem_source::{path}", a.get("filesystem_source", b.get("filesystem_source")), "text"),
        (f"watch_filesystem_type::{path}", a.get("filesystem_type", b.get("filesystem_type")), "text"),
        (f"watch_after_available_bytes::{path}", a.get("filesystem_available_bytes"), "bytes"),
    ])

with open(summary_path, "w", newline="", encoding="utf-8") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["metric", "value", "unit_or_definition"])
    for name, value, unit in metrics:
        if isinstance(value, float):
            value = f"{value:.6f}"
        if value is None:
            value = "NA"
        writer.writerow([name, value, unit])
PYSUMMARY

{
    printf 'end_utc\t%s\n' "$end_utc"
    printf 'exit_status\t%s\n' "$status"
    printf 'wall_time_seconds\t'
    awk -v a="$start_epoch_ns" -v b="$end_epoch_ns" 'BEGIN { printf "%.6f\n", (b-a)/1000000000 }'
} >> "$metadata_file"

printf '\n============================================================\n'
if ((status == 0)); then
    printf '[BENCHMARK] Installation command completed successfully.\n'
else
    printf '[BENCHMARK] Installation command failed with exit status %s.\n' "$status"
fi
printf '[BENCHMARK] Results: %s\n' "$output_dir"
printf '[BENCHMARK] Main summary: %s\n' "$summary_tsv"
printf '============================================================\n'

exit "$status"
