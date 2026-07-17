#!/usr/bin/env bash
# Download the minimum Ensembl inputs required for batch EggNOG/OrgDb creation:
#   1) gene annotation in GTF format
#   2) translated proteins in pep.all FASTA format
#
# The script reads HostSpecies.csv and creates one directory per species, e.g.:
#   ~/org_species/M.musculus/
#
# It uses only Python 3 standard-library modules plus the system gzip command.

set -Eeuo pipefail

SCRIPT_VERSION="1.1.0"

usage() {
    cat <<'EOF'
Usage:
  download_all_hostspecies_orgdb_inputs.sh [options]

Required input:
  -i, --hostspecies FILE   HostSpecies.csv file
                           Default: ./HostSpecies.csv

Output:
  -o, --output DIR         Root output directory
                           Default: ~/org_species

Selection:
  -t, --taxid TAXID        Download only this TaxID. May be repeated.
  -n, --limit N            Process only the first N selected rows.

Download behavior:
  -j, --jobs N             Number of species processed concurrently.
                           Default: 2
      --retries N          Download/discovery attempts. Default: 5
      --timeout SEC        Network timeout per request. Default: 90
      --force              Download files again even when valid copies exist.
      --dry-run            Resolve URLs and create a report, but download nothing.

Other:
      --keep-partials      Keep .part files after a final download failure.
      --version            Show version.
  -h, --help               Show this help.

Examples:
  # Safe first test with mouse only
  bash download_all_hostspecies_orgdb_inputs.sh \
      --hostspecies HostSpecies.csv \
      --output ~/org_species \
      --taxid 10090

  # Process every species in HostSpecies.csv
  bash download_all_hostspecies_orgdb_inputs.sh \
      --hostspecies HostSpecies.csv \
      --output ~/org_species \
      --jobs 3

Output for Mus musculus:
  ~/org_species/M.musculus/
  ├── references/
  │   ├── <original Ensembl GTF filename>.gtf.gz
  │   ├── <original Ensembl protein filename>.pep.all.fa.gz
  │   ├── annotation.gtf.gz -> <original GTF filename>
  │   └── proteins.pep.all.fa.gz -> <original protein filename>
  ├── metadata.tsv
  ├── checksums.sha256
  └── download.log

A global report is written to:
  <output>/download_report.tsv
EOF
}

HOSTSPECIES="./HostSpecies.csv"
OUTPUT_DIR="${HOME}/org_species"
JOBS=2
RETRIES=5
TIMEOUT=90
FORCE=0
DRY_RUN=0
KEEP_PARTIALS=0
LIMIT=0
TAXIDS=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--hostspecies)
            [[ $# -ge 2 ]] || { echo "[ERROR] Missing value after $1" >&2; exit 2; }
            HOSTSPECIES="$2"
            shift 2
            ;;
        -o|--output)
            [[ $# -ge 2 ]] || { echo "[ERROR] Missing value after $1" >&2; exit 2; }
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -t|--taxid)
            [[ $# -ge 2 ]] || { echo "[ERROR] Missing value after $1" >&2; exit 2; }
            TAXIDS+=("$2")
            shift 2
            ;;
        -n|--limit)
            [[ $# -ge 2 ]] || { echo "[ERROR] Missing value after $1" >&2; exit 2; }
            LIMIT="$2"
            shift 2
            ;;
        -j|--jobs)
            [[ $# -ge 2 ]] || { echo "[ERROR] Missing value after $1" >&2; exit 2; }
            JOBS="$2"
            shift 2
            ;;
        --retries)
            [[ $# -ge 2 ]] || { echo "[ERROR] Missing value after $1" >&2; exit 2; }
            RETRIES="$2"
            shift 2
            ;;
        --timeout)
            [[ $# -ge 2 ]] || { echo "[ERROR] Missing value after $1" >&2; exit 2; }
            TIMEOUT="$2"
            shift 2
            ;;
        --force)
            FORCE=1
            shift
            ;;
        --dry-run)
            DRY_RUN=1
            shift
            ;;
        --keep-partials)
            KEEP_PARTIALS=1
            shift
            ;;
        --version)
            echo "$SCRIPT_VERSION"
            exit 0
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        --)
            shift
            break
            ;;
        *)
            echo "[ERROR] Unknown option: $1" >&2
            usage >&2
            exit 2
            ;;
    esac
done

for value_name in JOBS RETRIES TIMEOUT LIMIT; do
    value="${!value_name}"
    if [[ ! "$value" =~ ^[0-9]+$ ]]; then
        echo "[ERROR] $value_name must be a non-negative integer: $value" >&2
        exit 2
    fi
done

if (( JOBS < 1 )); then
    echo "[ERROR] --jobs must be at least 1." >&2
    exit 2
fi
if (( RETRIES < 1 )); then
    echo "[ERROR] --retries must be at least 1." >&2
    exit 2
fi
if (( TIMEOUT < 1 )); then
    echo "[ERROR] --timeout must be at least 1 second." >&2
    exit 2
fi

command -v python3 >/dev/null 2>&1 || {
    echo "[ERROR] python3 is required." >&2
    exit 1
}

command -v gzip >/dev/null 2>&1 || {
    echo "[ERROR] gzip is required." >&2
    exit 1
}

if [[ ! -f "$HOSTSPECIES" ]]; then
    echo "[ERROR] HostSpecies.csv not found: $HOSTSPECIES" >&2
    exit 1
fi

mkdir -p "$OUTPUT_DIR"
HOSTSPECIES="$(cd -- "$(dirname -- "$HOSTSPECIES")" && pwd -P)/$(basename -- "$HOSTSPECIES")"
OUTPUT_DIR="$(cd -- "$OUTPUT_DIR" && pwd -P)"

TAXID_JOINED=""
if (( ${#TAXIDS[@]} > 0 )); then
    TAXID_JOINED="$(IFS=,; echo "${TAXIDS[*]}")"
fi

export MTD_DL_HOSTSPECIES="$HOSTSPECIES"
export MTD_DL_OUTPUT_DIR="$OUTPUT_DIR"
export MTD_DL_JOBS="$JOBS"
export MTD_DL_RETRIES="$RETRIES"
export MTD_DL_TIMEOUT="$TIMEOUT"
export MTD_DL_FORCE="$FORCE"
export MTD_DL_DRY_RUN="$DRY_RUN"
export MTD_DL_KEEP_PARTIALS="$KEEP_PARTIALS"
export MTD_DL_LIMIT="$LIMIT"
export MTD_DL_TAXIDS="$TAXID_JOINED"
export MTD_DL_SCRIPT_VERSION="$SCRIPT_VERSION"

python3 - <<'PY'
from __future__ import annotations

import concurrent.futures
import contextlib
import csv
import datetime as dt
import gzip
import hashlib
import html.parser
import json
import os
import re
import shutil
import socket
import subprocess
import sys
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

HOSTSPECIES = Path(os.environ["MTD_DL_HOSTSPECIES"])
OUTPUT_DIR = Path(os.environ["MTD_DL_OUTPUT_DIR"])
JOBS = int(os.environ["MTD_DL_JOBS"])
RETRIES = int(os.environ["MTD_DL_RETRIES"])
TIMEOUT = int(os.environ["MTD_DL_TIMEOUT"])
FORCE = os.environ["MTD_DL_FORCE"] == "1"
DRY_RUN = os.environ["MTD_DL_DRY_RUN"] == "1"
KEEP_PARTIALS = os.environ["MTD_DL_KEEP_PARTIALS"] == "1"
LIMIT = int(os.environ["MTD_DL_LIMIT"])
REQUESTED_TAXIDS = {
    value.strip()
    for value in os.environ.get("MTD_DL_TAXIDS", "").split(",")
    if value.strip()
}
SCRIPT_VERSION = os.environ["MTD_DL_SCRIPT_VERSION"]

USER_AGENT = f"MTD-Explorer-OrgDb-input-downloader/{SCRIPT_VERSION}"
PRINT_LOCK = threading.Lock()

ENSEMBL_REST = "https://rest.ensembl.org"
DIVISIONS = (
    "EnsemblVertebrates",
    "EnsemblMetazoa",
    "EnsemblFungi",
    "EnsemblPlants",
    "EnsemblProtists",
)

# Ensembl release 116 uses /pub/current/{gtf,fasta}; the older
# /pub/current_gtf and /pub/current_fasta aliases are not reliable.
MAIN_GTF_ROOT = "https://ftp.ensembl.org/pub/current/gtf"
MAIN_FASTA_ROOT = "https://ftp.ensembl.org/pub/current/fasta"
GENOMES_ROOT = "https://ftp.ensemblgenomes.ebi.ac.uk/pub"

REPORT_FIELDS = [
    "taxid",
    "scientific_name",
    "folder",
    "ensembl_name",
    "division",
    "assembly",
    "release",
    "source",
    "gtf_url",
    "pep_url",
    "gtf_file",
    "pep_file",
    "gtf_status",
    "pep_status",
    "status",
    "message",
]


class LinkParser(html.parser.HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.links: List[str] = []

    def handle_starttag(self, tag: str, attrs: Sequence[Tuple[str, Optional[str]]]) -> None:
        if tag.lower() != "a":
            return
        for key, value in attrs:
            if key.lower() == "href" and value:
                self.links.append(value)


def say(message: str) -> None:
    with PRINT_LOCK:
        print(message, flush=True)


def request_bytes(url: str, *, headers: Optional[Dict[str, str]] = None) -> bytes:
    merged = {"User-Agent": USER_AGENT, "Accept": "*/*"}
    if headers:
        merged.update(headers)

    last_error: Optional[BaseException] = None
    attempts_used = 0

    for attempt in range(1, RETRIES + 1):
        attempts_used = attempt
        try:
            request = urllib.request.Request(url, headers=merged)
            with urllib.request.urlopen(request, timeout=TIMEOUT) as response:
                return response.read()

        except urllib.error.HTTPError as exc:
            last_error = exc

            # Permanent client errors will not become valid after retrying.
            # In particular, avoiding five retries for every 404 makes species
            # discovery much faster when a source/division does not contain it.
            if exc.code in {400, 401, 403, 404, 405, 410}:
                break

            if attempt == RETRIES:
                break

            retry_after = exc.headers.get("Retry-After")
            if retry_after and retry_after.isdigit():
                delay = min(120, int(retry_after))
            else:
                delay = min(30, 2 ** (attempt - 1))
            time.sleep(delay)

        except (urllib.error.URLError, TimeoutError, socket.timeout) as exc:
            last_error = exc
            if attempt == RETRIES:
                break
            time.sleep(min(30, 2 ** (attempt - 1)))

    raise RuntimeError(
        f"Unable to read URL after {attempts_used} attempt(s): "
        f"{url}: {last_error}"
    )


def request_json(url: str) -> dict:
    payload = request_bytes(url, headers={"Content-Type": "application/json"})
    return json.loads(payload.decode("utf-8"))


def list_directory(url: str) -> List[str]:
    payload = request_bytes(url).decode("utf-8", errors="replace")
    parser = LinkParser()
    parser.feed(payload)

    names: List[str] = []
    for href in parser.links:
        parsed = urllib.parse.urlsplit(href)
        name = urllib.parse.unquote(Path(parsed.path).name)
        if name and name not in {".", ".."}:
            names.append(name)
    return sorted(set(names))


def safe_slug(text: str) -> str:
    value = text.strip().lower()
    value = re.sub(r"[^a-z0-9]+", "_", value)
    return value.strip("_")


def abbreviated_folder(scientific_name: str, taxid: str, used: Dict[str, str]) -> str:
    parts = [safe_slug(part) for part in scientific_name.split() if safe_slug(part)]
    if not parts:
        candidate = f"taxid_{taxid}"
    elif len(parts) == 1:
        candidate = parts[0]
    else:
        candidate = f"{parts[0][0].upper()}.{parts[1]}"
        if len(parts) > 2:
            candidate += "_" + "_".join(parts[2:])

    previous = used.get(candidate)
    if previous is not None and previous != taxid:
        candidate = f"{candidate}_{taxid}"
    used[candidate] = taxid
    return candidate


def load_rows() -> List[dict]:
    with HOSTSPECIES.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"Taxon_ID", "Scientific_name"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise RuntimeError(
                "HostSpecies.csv is missing required column(s): " + ", ".join(sorted(missing))
            )
        rows = list(reader)

    selected: List[dict] = []
    seen: set[str] = set()
    for row in rows:
        taxid = str(row.get("Taxon_ID", "")).strip()
        scientific_name = str(row.get("Scientific_name", "")).strip()
        if not taxid or not scientific_name:
            continue
        if not taxid.isdigit():
            say(f"[WARN] Skipping non-numeric TaxID: {taxid!r}")
            continue
        if taxid in seen:
            say(f"[WARN] Skipping duplicate TaxID in CSV: {taxid}")
            continue
        seen.add(taxid)
        if REQUESTED_TAXIDS and taxid not in REQUESTED_TAXIDS:
            continue
        selected.append(row)

    if REQUESTED_TAXIDS:
        missing_taxids = sorted(REQUESTED_TAXIDS.difference(seen))
        if missing_taxids:
            raise RuntimeError(
                "Requested TaxID(s) not found in HostSpecies.csv: " + ", ".join(missing_taxids)
            )

    if LIMIT > 0:
        selected = selected[:LIMIT]
    return selected


def fetch_species_catalog() -> Dict[str, List[dict]]:
    index: Dict[str, List[dict]] = {}
    errors: List[str] = []

    for division in DIVISIONS:
        query = urllib.parse.urlencode({"division": division, "content-type": "application/json"})
        url = f"{ENSEMBL_REST}/info/species?{query}"
        try:
            data = request_json(url)
        except Exception as exc:
            errors.append(f"{division}: {exc}")
            continue

        for item in data.get("species", []):
            taxid = str(item.get("taxon_id", "")).strip()
            if not taxid:
                continue
            enriched = dict(item)
            enriched.setdefault("division", division)
            index.setdefault(taxid, []).append(enriched)

    if not index:
        raise RuntimeError(
            "Could not retrieve the Ensembl species catalog. " + " | ".join(errors)
        )

    if errors:
        say("[WARN] Some Ensembl catalog divisions could not be read:")
        for error in errors:
            say(f"       {error}")
    return index


def choose_catalog_entry(row: dict, entries: Sequence[dict]) -> Optional[dict]:
    if not entries:
        return None

    scientific = safe_slug(str(row.get("Scientific_name", "")))
    dataset = str(row.get("MartDatasets", "")).strip().lower()

    def score(item: dict) -> Tuple[int, int, str]:
        name = safe_slug(str(item.get("name", "")))
        display = safe_slug(str(item.get("display_name", "")))
        aliases = [safe_slug(str(alias)) for alias in item.get("aliases", []) or []]
        value = 0
        if display == scientific:
            value += 100
        if name == scientific:
            value += 90
        if scientific in aliases:
            value += 80
        prefix = dataset.split("_gene_")[0]
        if prefix and prefix in name.replace("_", ""):
            value += 20
        if str(item.get("assembly", "")).strip():
            value += 5
        return (value, -len(name), name)

    return max(entries, key=score)


def source_candidates(entry: Optional[dict], scientific_name: str) -> List[Tuple[str, str, str]]:
    names: List[str] = []
    if entry:
        name = safe_slug(str(entry.get("name", "")))
        if name:
            names.append(name)
    derived = safe_slug(scientific_name)
    if derived and derived not in names:
        names.append(derived)

    division = str((entry or {}).get("division", ""))
    division_short = {
        "EnsemblMetazoa": "metazoa",
        "EnsemblFungi": "fungi",
        "EnsemblPlants": "plants",
        "EnsemblProtists": "protists",
        "EnsemblBacteria": "bacteria",
    }.get(division)

    candidates: List[Tuple[str, str, str]] = []
    for name in names:
        # The main Ensembl FTP is checked first because HostSpecies.csv is primarily
        # derived from the main Ensembl BioMart catalog.
        candidates.append(
            (
                "Ensembl",
                f"{MAIN_GTF_ROOT}/{name}/",
                f"{MAIN_FASTA_ROOT}/{name}/pep/",
            )
        )

        if division_short:
            candidates.append(
                (
                    f"EnsemblGenomes:{division_short}",
                    f"{GENOMES_ROOT}/{division_short}/current/gtf/{name}/",
                    f"{GENOMES_ROOT}/{division_short}/current/fasta/{name}/pep/",
                )
            )

    # Biomphalaria and other external/custom rows may be absent from the REST result.
    # Probe the non-bacterial Ensembl Genomes divisions as a fallback.
    for name in names:
        for short in ("metazoa", "fungi", "plants", "protists"):
            item = (
                f"EnsemblGenomes:{short}",
                f"{GENOMES_ROOT}/{short}/current/gtf/{name}/",
                f"{GENOMES_ROOT}/{short}/current/fasta/{name}/pep/",
            )
            if item not in candidates:
                candidates.append(item)
    return candidates


def select_gtf(names: Iterable[str]) -> Optional[str]:
    candidates = [name for name in names if name.lower().endswith(".gtf.gz")]
    if not candidates:
        return None

    preferred = [
        name
        for name in candidates
        if ".abinitio." not in name.lower()
        and ".chr." not in name.lower()
        and "chromosome" not in name.lower()
    ]
    if preferred:
        candidates = preferred

    # Standard whole-annotation GTF files generally have the shortest name.
    return sorted(candidates, key=lambda value: (len(value), value.lower()))[0]


def select_pep(names: Iterable[str]) -> Optional[str]:
    exact = [name for name in names if name.lower().endswith(".pep.all.fa.gz")]
    if exact:
        return sorted(exact, key=lambda value: (len(value), value.lower()))[0]

    fallback = [
        name
        for name in names
        if name.lower().endswith((".pep.fa.gz", ".pep.fasta.gz"))
        and "abinitio" not in name.lower()
    ]
    return sorted(fallback, key=lambda value: (len(value), value.lower()))[0] if fallback else None


def discover_files(
    entry: Optional[dict], scientific_name: str
) -> Tuple[str, str, str, str, str]:
    errors: List[str] = []
    for source, gtf_dir, pep_dir in source_candidates(entry, scientific_name):
        try:
            gtf_name = select_gtf(list_directory(gtf_dir))
            if not gtf_name:
                errors.append(f"{source}: no GTF in {gtf_dir}")
                continue
            pep_name = select_pep(list_directory(pep_dir))
            if not pep_name:
                errors.append(f"{source}: no pep.all FASTA in {pep_dir}")
                continue
            return (
                source,
                urllib.parse.urljoin(gtf_dir, urllib.parse.quote(gtf_name)),
                urllib.parse.urljoin(pep_dir, urllib.parse.quote(pep_name)),
                gtf_name,
                pep_name,
            )
        except Exception as exc:
            errors.append(f"{source}: {exc}")

    raise RuntimeError("Reference files could not be resolved. " + " | ".join(errors))


def gzip_ok(path: Path) -> bool:
    if not path.is_file() or path.stat().st_size == 0:
        return False
    result = subprocess.run(
        ["gzip", "-t", str(path)],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=False,
    )
    return result.returncode == 0


def gtf_content_ok(path: Path) -> bool:
    try:
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                return len(line.rstrip("\n").split("\t")) >= 9
    except OSError:
        return False
    return False


def pep_content_ok(path: Path) -> bool:
    try:
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if not line.strip():
                    continue
                return line.startswith(">")
    except OSError:
        return False
    return False


def validate(path: Path, kind: str) -> bool:
    if not gzip_ok(path):
        return False
    if kind == "gtf":
        return gtf_content_ok(path)
    if kind == "pep":
        return pep_content_ok(path)
    return True


def download(url: str, destination: Path, kind: str, log) -> str:
    if destination.exists() and not FORCE and validate(destination, kind):
        log.write(f"[SKIP] Valid existing file: {destination}\n")
        return "EXISTING_VALID"

    if DRY_RUN:
        log.write(f"[DRY-RUN] {url} -> {destination}\n")
        return "DRY_RUN"

    if destination.exists():
        invalid = destination.with_name(destination.name + ".invalid")
        with contextlib.suppress(FileNotFoundError):
            invalid.unlink()
        destination.replace(invalid)
        log.write(f"[WARN] Existing file moved to: {invalid}\n")

    partial = destination.with_name(destination.name + ".part")
    destination.parent.mkdir(parents=True, exist_ok=True)

    last_error: Optional[BaseException] = None
    for attempt in range(1, RETRIES + 1):
        try:
            offset = partial.stat().st_size if partial.exists() else 0
            headers = {"User-Agent": USER_AGENT, "Accept": "*/*"}
            if offset > 0:
                headers["Range"] = f"bytes={offset}-"

            request = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(request, timeout=TIMEOUT) as response:
                status = getattr(response, "status", response.getcode())
                append = offset > 0 and status == 206
                mode = "ab" if append else "wb"
                if offset > 0 and not append:
                    log.write("[INFO] Server did not honor Range; restarting file.\n")

                with partial.open(mode) as output:
                    while True:
                        chunk = response.read(8 * 1024 * 1024)
                        if not chunk:
                            break
                        output.write(chunk)

            if not validate(partial, kind):
                raise RuntimeError("downloaded gzip/content validation failed")

            os.replace(partial, destination)
            log.write(f"[OK] Downloaded and validated: {destination}\n")
            return "DOWNLOADED"

        except Exception as exc:
            last_error = exc
            log.write(f"[WARN] Attempt {attempt}/{RETRIES} failed: {exc}\n")
            log.flush()
            if attempt < RETRIES:
                time.sleep(min(60, 2 ** (attempt - 1)))

    if partial.exists() and not KEEP_PARTIALS:
        partial.unlink()
    raise RuntimeError(f"download failed after {RETRIES} attempts: {last_error}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(8 * 1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def replace_symlink(link: Path, target_name: str) -> None:
    if link.is_symlink() or link.exists():
        link.unlink()
    link.symlink_to(target_name)


def write_metadata(path: Path, values: Sequence[Tuple[str, str]]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["field", "value"])
        writer.writerows(values)
    os.replace(temporary, path)


def process_species(row: dict, entry: Optional[dict], folder_name: str) -> dict:
    taxid = str(row["Taxon_ID"]).strip()
    scientific_name = str(row["Scientific_name"]).strip()
    species_dir = OUTPUT_DIR / folder_name
    reference_dir = species_dir / "references"
    species_dir.mkdir(parents=True, exist_ok=True)
    reference_dir.mkdir(parents=True, exist_ok=True)
    log_path = species_dir / "download.log"

    result = {field: "" for field in REPORT_FIELDS}
    result.update(
        {
            "taxid": taxid,
            "scientific_name": scientific_name,
            "folder": folder_name,
            "ensembl_name": str((entry or {}).get("name", "")),
            "division": str((entry or {}).get("division", "")),
            "assembly": str((entry or {}).get("assembly", "")),
            "release": str((entry or {}).get("release", "")),
            "status": "FAILED",
        }
    )

    say(f"[INFO] [{taxid}] {scientific_name} -> {folder_name}")

    with log_path.open("a", encoding="utf-8") as log:
        started = dt.datetime.now(dt.timezone.utc).isoformat()
        log.write("=" * 72 + "\n")
        log.write(f"Started UTC: {started}\n")
        log.write(f"TaxID: {taxid}\n")
        log.write(f"Scientific name: {scientific_name}\n")

        try:
            source, gtf_url, pep_url, gtf_name, pep_name = discover_files(
                entry, scientific_name
            )
            result.update(
                {
                    "source": source,
                    "gtf_url": gtf_url,
                    "pep_url": pep_url,
                    "gtf_file": gtf_name,
                    "pep_file": pep_name,
                }
            )
            log.write(f"Source: {source}\n")
            log.write(f"GTF URL: {gtf_url}\n")
            log.write(f"PEP URL: {pep_url}\n")

            gtf_path = reference_dir / gtf_name
            pep_path = reference_dir / pep_name

            result["gtf_status"] = download(gtf_url, gtf_path, "gtf", log)
            result["pep_status"] = download(pep_url, pep_path, "pep", log)

            if not DRY_RUN:
                replace_symlink(reference_dir / "annotation.gtf.gz", gtf_name)
                replace_symlink(reference_dir / "proteins.pep.all.fa.gz", pep_name)

                checksums = species_dir / "checksums.sha256"
                with checksums.open("w", encoding="utf-8") as handle:
                    handle.write(f"{sha256(gtf_path)}  references/{gtf_name}\n")
                    handle.write(f"{sha256(pep_path)}  references/{pep_name}\n")

            metadata = [
                ("taxid", taxid),
                ("scientific_name", scientific_name),
                ("folder", folder_name),
                ("mart_dataset", str(row.get("MartDatasets", ""))),
                ("orgdb_previous", str(row.get("OrgDb", ""))),
                ("ensembl_name", result["ensembl_name"]),
                ("division", result["division"]),
                ("assembly", result["assembly"]),
                ("release", result["release"]),
                ("source", source),
                ("gtf_url", gtf_url),
                ("pep_url", pep_url),
                ("gtf_file", gtf_name),
                ("pep_file", pep_name),
                ("downloaded_utc", dt.datetime.now(dt.timezone.utc).isoformat()),
                ("status", "DRY_RUN" if DRY_RUN else "COMPLETE"),
            ]
            write_metadata(species_dir / "metadata.tsv", metadata)

            result["status"] = "DRY_RUN" if DRY_RUN else "COMPLETE"
            result["message"] = ""
            say(f"[OK]   [{taxid}] {scientific_name}: {result['status']}")

        except Exception as exc:
            result["message"] = str(exc).replace("\t", " ").replace("\n", " ")
            log.write(f"[ERROR] {exc}\n")
            say(f"[FAIL] [{taxid}] {scientific_name}: {exc}")

    return result


def write_report(rows: Sequence[dict]) -> Path:
    report = OUTPUT_DIR / "download_report.tsv"
    temporary = report.with_suffix(report.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=REPORT_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, report)
    return report


def main() -> int:
    rows = load_rows()
    if not rows:
        raise RuntimeError("No species were selected from HostSpecies.csv.")

    say("============================================================")
    say("MTD Explorer — HostSpecies OrgDb input downloader")
    say(f"HostSpecies.csv: {HOSTSPECIES}")
    say(f"Output:          {OUTPUT_DIR}")
    say(f"Species:         {len(rows)}")
    say(f"Parallel jobs:   {JOBS}")
    say(f"Dry run:         {'yes' if DRY_RUN else 'no'}")
    say("============================================================")
    say("[INFO] Reading the current Ensembl species catalog...")

    catalog = fetch_species_catalog()
    folder_map: Dict[str, str] = {}
    work: List[Tuple[dict, Optional[dict], str]] = []

    for row in rows:
        taxid = str(row["Taxon_ID"]).strip()
        scientific_name = str(row["Scientific_name"]).strip()
        entry = choose_catalog_entry(row, catalog.get(taxid, []))
        folder = abbreviated_folder(scientific_name, taxid, folder_map)
        work.append((row, entry, folder))

    results: List[dict] = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=JOBS) as executor:
        future_map = {
            executor.submit(process_species, row, entry, folder): str(row["Taxon_ID"])
            for row, entry, folder in work
        }
        for future in concurrent.futures.as_completed(future_map):
            results.append(future.result())

    order = {str(row["Taxon_ID"]).strip(): index for index, row in enumerate(rows)}
    results.sort(key=lambda item: order.get(item["taxid"], 10**9))
    report = write_report(results)

    complete = sum(result["status"] == "COMPLETE" for result in results)
    dry = sum(result["status"] == "DRY_RUN" for result in results)
    failed = sum(result["status"] == "FAILED" for result in results)

    say("============================================================")
    say(f"Complete: {complete}")
    say(f"Dry run:  {dry}")
    say(f"Failed:   {failed}")
    say(f"Report:   {report}")
    say("============================================================")

    return 1 if failed else 0


try:
    raise SystemExit(main())
except KeyboardInterrupt:
    say("\n[ERROR] Interrupted by user.")
    raise SystemExit(130)
except Exception as exc:
    say(f"[ERROR] {exc}")
    raise SystemExit(1)
PY
