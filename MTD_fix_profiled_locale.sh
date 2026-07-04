#!/usr/bin/env bash
# Fix locale-sensitive EPOCHREALTIME parsing in an existing Install_profiled.sh.
set -euo pipefail

target="${1:-./Install_profiled.sh}"

die() {
    printf '[ERROR] %s\n' "$*" >&2
    exit 1
}

[[ -f "$target" ]] || die "File not found: $target"
grep -q 'BEGIN MTD FUNCTION PROFILER' "$target" ||
    die "The file does not contain the MTD function profiler."

if grep -Fq 'realtime="${EPOCHREALTIME/,/.}"' "$target"; then
    printf '[OK] Locale correction is already present in %s\n' "$target"
    exit 0
fi

backup="${target}.bak.$(date -u '+%Y%m%dT%H%M%SZ')"
cp -a "$target" "$backup"

python3 - "$target" <<'PY'
from pathlib import Path
import sys

path = Path(sys.argv[1])
text = path.read_text(encoding="utf-8")

old = '        realtime="$EPOCHREALTIME"\n'
new = (
    '        # Normalize a locale decimal comma before parsing EPOCHREALTIME.\n'
    '        realtime="${EPOCHREALTIME/,/.}"\n'
)
if old not in text:
    raise SystemExit("[ERROR] Expected EPOCHREALTIME line was not found.")
text = text.replace(old, new, 1)

old2 = (
    'mtd_bench_seconds_to_us() {\n'
    '    local value="$1"\n'
    '    local whole fraction\n\n'
    '    whole="${value%%.*}"\n'
)
new2 = (
    'mtd_bench_seconds_to_us() {\n'
    '    local value="$1"\n'
    '    local whole fraction\n\n'
    '    value="${value/,/.}"\n'
    '    whole="${value%%.*}"\n'
)
if old2 in text:
    text = text.replace(old2, new2, 1)

path.write_text(text, encoding="utf-8")
PY

if ! bash -n "$target"; then
    cp -a "$backup" "$target"
    die "Syntax validation failed; the original file was restored from $backup"
fi

printf '[OK] Corrected profiler: %s\n' "$(readlink -f "$target")"
printf '[OK] Backup preserved:  %s\n' "$(readlink -f "$backup")"
