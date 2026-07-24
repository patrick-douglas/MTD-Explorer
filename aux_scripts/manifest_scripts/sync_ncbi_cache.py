#!/usr/bin/env python3
from __future__ import annotations

import argparse
import concurrent.futures
import csv
import gzip
import hashlib
import html.parser
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
import urllib.error
import urllib.request
from dataclasses import dataclass
from email.utils import parsedate_to_datetime
from pathlib import Path
from typing import Iterable


USER_AGENT = "MTD-Explorer-NCBI-cache-sync/1.0"
DEFAULT_RETRIES = 5
DEFAULT_TIMEOUT = 60


@dataclass(frozen=True)
class RemoteEntry:
    name: str
    url: str
    size: int
    modified: str
    signature: str


class LinkParser(html.parser.HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.links: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag.lower() != "a":
            return
        for key, value in attrs:
            if key.lower() == "href" and value:
                self.links.append(value)


def log(level: str, message: str) -> None:
    print(f"[{level}] {message}", flush=True)


def fail(message: str) -> "NoReturn":
    log("FAIL", message)
    raise SystemExit(1)


def request_bytes(url: str, *, method: str = "GET", retries: int = DEFAULT_RETRIES) -> tuple[bytes, dict[str, str]]:
    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        request = urllib.request.Request(
            url,
            method=method,
            headers={"User-Agent": USER_AGENT},
        )
        try:
            with urllib.request.urlopen(request, timeout=DEFAULT_TIMEOUT) as response:
                data = response.read() if method != "HEAD" else b""
                headers = {key.lower(): value for key, value in response.headers.items()}
                return data, headers
        except (urllib.error.URLError, TimeoutError, OSError) as error:
            last_error = error
            if attempt < retries:
                log("WARN", f"Network attempt {attempt}/{retries} failed for {url}: {error}")
                time.sleep(min(5 * attempt, 20))
    raise RuntimeError(f"Could not retrieve {url}: {last_error}")


def atomic_write(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temp_name = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    temp = Path(temp_name)
    try:
        with os.fdopen(fd, "wb") as handle:
            handle.write(data)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temp, path)
    finally:
        temp.unlink(missing_ok=True)


def atomic_write_text(path: Path, text: str) -> None:
    atomic_write(path, text.encode("utf-8"))


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def gzip_ok(path: Path) -> bool:
    if not path.is_file() or path.stat().st_size <= 0:
        return False
    try:
        with gzip.open(path, "rb") as handle:
            while handle.read(8 * 1024 * 1024):
                pass
        return True
    except (OSError, EOFError):
        return False


def load_catalog(path: Path) -> dict[str, RemoteEntry]:
    entries: dict[str, RemoteEntry] = {}
    if not path.is_file():
        return entries
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"name", "url", "size", "modified", "signature"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            return {}
        for row in reader:
            name = row["name"].strip()
            if not name:
                continue
            try:
                size = int(row["size"] or 0)
            except ValueError:
                size = 0
            entries[name] = RemoteEntry(
                name=name,
                url=row["url"].strip(),
                size=size,
                modified=row["modified"].strip(),
                signature=row["signature"].strip(),
            )
    return entries


def catalog_text(entries: Iterable[RemoteEntry]) -> str:
    lines = ["name\turl\tsize\tmodified\tsignature"]
    for entry in entries:
        lines.append(
            "\t".join(
                [
                    entry.name,
                    entry.url,
                    str(entry.size),
                    entry.modified,
                    entry.signature,
                ]
            )
        )
    return "\n".join(lines) + "\n"


def load_integrity_state(path: Path) -> dict[str, tuple[int, int]]:
    state: dict[str, tuple[int, int]] = {}
    if not path.is_file():
        return state
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            try:
                state[fields[0]] = (int(fields[1]), int(float(fields[2])))
            except ValueError:
                continue
    return state


def write_integrity_state(path: Path, local_dir: Path, names: Iterable[str]) -> None:
    rows: list[str] = []
    for name in sorted(set(names)):
        file_path = local_dir / name
        if not file_path.is_file():
            continue
        stat_result = file_path.stat()
        rows.append(f"{name}\t{stat_result.st_size}\t{stat_result.st_mtime_ns}")
    atomic_write_text(path, "\n".join(rows) + ("\n" if rows else ""))


def remove_paths(paths: Iterable[Path]) -> None:
    for path in paths:
        path.unlink(missing_ok=True)
        Path(f"{path}.aria2").unlink(missing_ok=True)
        Path(f"{path}.part").unlink(missing_ok=True)
        Path(f"{path}.part.aria2").unlink(missing_ok=True)


def download_one(entry: RemoteEntry, local_dir: Path, retries: int) -> tuple[str, str | None]:
    destination = local_dir / entry.name
    partial = local_dir / f"{entry.name}.part"
    aria_sidecar = Path(f"{partial}.aria2")
    local_dir.mkdir(parents=True, exist_ok=True)

    for attempt in range(1, retries + 1):
        partial.unlink(missing_ok=True)
        aria_sidecar.unlink(missing_ok=True)

        try:
            if shutil.which("aria2c"):
                command = [
                    "aria2c",
                    "--disable-ipv6=true",
                    "--continue=true",
                    "--auto-file-renaming=false",
                    "--allow-overwrite=true",
                    "--max-tries=3",
                    "--retry-wait=5",
                    "--timeout=60",
                    "--connect-timeout=30",
                    "-x4",
                    "-s4",
                    "-d",
                    str(local_dir),
                    "-o",
                    partial.name,
                    entry.url,
                ]
                completed = subprocess.run(
                    command,
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.PIPE,
                    text=True,
                )
                if completed.returncode != 0:
                    raise RuntimeError(completed.stderr.strip() or f"aria2c exit {completed.returncode}")
            else:
                request = urllib.request.Request(entry.url, headers={"User-Agent": USER_AGENT})
                with urllib.request.urlopen(request, timeout=DEFAULT_TIMEOUT) as response, partial.open("wb") as output:
                    shutil.copyfileobj(response, output, length=8 * 1024 * 1024)

            if entry.size > 0 and partial.stat().st_size != entry.size:
                raise RuntimeError(
                    f"size mismatch: expected {entry.size}, observed {partial.stat().st_size}"
                )

            if not gzip_ok(partial):
                raise RuntimeError("gzip integrity validation failed")

            os.replace(partial, destination)
            aria_sidecar.unlink(missing_ok=True)
            return entry.name, None
        except Exception as error:
            partial.unlink(missing_ok=True)
            aria_sidecar.unlink(missing_ok=True)
            if attempt < retries:
                time.sleep(min(5 * attempt, 20))
            else:
                return entry.name, str(error)

    return entry.name, "unexpected download failure"


def download_entries(
    entries: list[RemoteEntry],
    local_dir: Path,
    jobs: int,
    retries: int,
) -> list[tuple[str, str]]:
    failures: list[tuple[str, str]] = []
    if not entries:
        return failures

    log("INFO", f"Downloading {len(entries)} missing or changed file(s) with {jobs} job(s).")
    with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:
        future_map = {
            executor.submit(download_one, entry, local_dir, retries): entry
            for entry in entries
        }
        completed_count = 0
        for future in concurrent.futures.as_completed(future_map):
            completed_count += 1
            entry = future_map[future]
            try:
                name, error = future.result()
            except Exception as exc:
                name, error = entry.name, str(exc)
            if error:
                failures.append((name, error))
                log("WARN", f"[{completed_count}/{len(entries)}] Failed: {name}: {error}")
            else:
                log("OK", f"[{completed_count}/{len(entries)}] Downloaded and validated: {name}")
    return failures


# MTD_NCBI_ASSEMBLY_HEADER_PARSER_V2
def parse_assembly_summary(
    data: bytes,
    summary_url: str,
    allowed_levels: set[str],
) -> list[RemoteEntry]:
    # utf-8-sig removes a possible UTF-8 BOM without affecting normal UTF-8.
    text = data.decode("utf-8-sig", errors="replace")
    header: list[str] | None = None
    rows: list[list[str]] = []

    for raw_line in text.splitlines():
        # NCBI files normally use "# assembly_accession", but some mirrors,
        # releases, and generated subsets may use "#assembly_accession",
        # leading whitespace, or an UTF-8 BOM. Normalize the comment prefix
        # before deciding whether this is the tab-delimited header.
        line = raw_line.lstrip("\ufeff")
        stripped = line.lstrip()

        if not stripped:
            continue

        if stripped.startswith("#"):
            candidate = stripped[1:].lstrip()
            fields = [field.strip() for field in candidate.split("\t")]

            if fields and fields[0].lower() == "assembly_accession":
                header = fields

            continue

        rows.append(line.rstrip("\r\n").split("\t"))

    if header is None:
        nonempty = [
            line.strip()
            for line in text.splitlines()
            if line.strip()
        ]
        preview = " | ".join(nonempty[:5])
        preview = preview[:800]

        lowered = text.lstrip().lower()
        if lowered.startswith("<!doctype html") or lowered.startswith("<html"):
            fail(
                "NCBI returned HTML instead of an assembly summary: "
                f"{summary_url}. Response preview: {preview}"
            )

        fail(
            "The NCBI assembly summary has no recognized tab-delimited "
            f"assembly_accession header: {summary_url}. "
            f"Response preview: {preview}"
        )

    # Remove any accidental leading comment marker or whitespace from every
    # column name while preserving the official NCBI names.
    header = [
        field.lstrip("\ufeff#").strip()
        for field in header
    ]

    index = {
        name: position
        for position, name in enumerate(header)
    }

    required = {
        "assembly_accession",
        "assembly_level",
        "version_status",
        "seq_rel_date",
        "asm_name",
        "ftp_path",
    }
    missing = required - set(index)

    if missing:
        fail(
            "Assembly summary is missing required columns: "
            + ", ".join(sorted(missing))
        )

    entries: list[RemoteEntry] = []
    seen: set[str] = set()

    for fields in rows:
        if len(fields) < len(header):
            continue

        ftp_path = fields[index["ftp_path"]].strip()
        assembly_level = fields[index["assembly_level"]].strip()

        if ftp_path == "na" or not ftp_path:
            continue

        if allowed_levels and assembly_level not in allowed_levels:
            continue

        ftp_path = re.sub(r"^ftp:", "https:", ftp_path)
        ftp_path = ftp_path.rstrip("/")
        assembly_dir = ftp_path.rsplit("/", 1)[-1]

        if not re.match(r"^GC[AF]_\d+\.\d+_", assembly_dir):
            continue

        name = f"{assembly_dir}_genomic.fna.gz"

        if name in seen:
            continue

        seen.add(name)
        url = f"{ftp_path}/{name}"

        signature_source = "\t".join(
            [
                fields[index["assembly_accession"]].strip(),
                assembly_level,
                fields[index["version_status"]].strip(),
                fields[index["seq_rel_date"]].strip(),
                fields[index["asm_name"]].strip(),
                ftp_path,
            ]
        )
        signature = hashlib.sha256(
            signature_source.encode("utf-8")
        ).hexdigest()

        entries.append(
            RemoteEntry(
                name=name,
                url=url,
                size=0,
                modified=fields[index["seq_rel_date"]].strip(),
                signature=signature,
            )
        )

    entries.sort(key=lambda item: item.name)
    return entries

def remote_head_entry(name: str, url: str, release: str) -> RemoteEntry:
    _, headers = request_bytes(url, method="HEAD")
    try:
        size = int(headers.get("content-length", "0"))
    except ValueError:
        size = 0
    modified = headers.get("last-modified", "")
    signature_source = f"{release}\t{name}\t{size}\t{modified}\t{url}"
    signature = hashlib.sha256(signature_source.encode("utf-8")).hexdigest()
    return RemoteEntry(name=name, url=url, size=size, modified=modified, signature=signature)


def parse_release_directory(
    base_url: str,
    release_number_url: str,
    pattern: re.Pattern[str],
    jobs: int,
) -> tuple[str, list[RemoteEntry]]:
    release_data, _ = request_bytes(release_number_url)
    release = release_data.decode("utf-8", errors="replace").strip()
    if not release:
        fail(f"Empty NCBI RefSeq release number from {release_number_url}")

    html_data, _ = request_bytes(base_url.rstrip("/") + "/")
    parser = LinkParser()
    parser.feed(html_data.decode("utf-8", errors="replace"))
    names = sorted({Path(link).name for link in parser.links if pattern.fullmatch(Path(link).name)})
    if not names:
        fail(f"No files matching {pattern.pattern!r} were found at {base_url}")

    entries: list[RemoteEntry] = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=min(jobs, len(names))) as executor:
        futures = {
            executor.submit(
                remote_head_entry,
                name,
                f"{base_url.rstrip('/')}/{name}",
                release,
            ): name
            for name in names
        }
        for future in concurrent.futures.as_completed(futures):
            name = futures[future]
            try:
                entries.append(future.result())
            except Exception as error:
                fail(f"Could not read remote metadata for {name}: {error}")

    entries.sort(key=lambda item: item.name)
    return release, entries


def validate_catalog_size(current_count: int, previous_count: int, minimum_count: int) -> None:
    if current_count < minimum_count:
        fail(
            f"Remote catalog is unexpectedly small: {current_count} entries; "
            f"minimum accepted is {minimum_count}. Existing cache was not modified."
        )
    if previous_count > 0:
        minimum_relative = max(minimum_count, previous_count * 80 // 100)
        if current_count < minimum_relative:
            fail(
                f"Remote catalog shrank unexpectedly: previous={previous_count}, "
                f"current={current_count}, minimum accepted={minimum_relative}. "
                "Existing cache was not modified."
            )


def synchronize(
    *,
    library: str,
    entries: list[RemoteEntry],
    local_dir: Path,
    metadata_dir: Path,
    summary_data: bytes | None,
    minimum_count: int,
    require_complete: bool,
    download_jobs: int,
    gzip_jobs: int,
    retries: int,
    full_gzip_check: bool,
    trust_existing_without_catalog: bool,
    combined_fasta: Path | None = None,
    force_combined: bool = False,
    check_only: bool = False,
) -> None:
    metadata_dir.mkdir(parents=True, exist_ok=True)
    local_dir.mkdir(parents=True, exist_ok=True)

    catalog_path = metadata_dir / f"remote_catalog_{library}.tsv"
    names_path = metadata_dir / f"manifest_{library}.names.txt"
    urls_path = metadata_dir / f"manifest_{library}.list.txt"
    integrity_path = metadata_dir / f"integrity_{library}.stat.tsv"
    failed_path = metadata_dir / f"failed_downloads_{library}.txt"
    obsolete_path = metadata_dir / f"obsolete_local_{library}.txt"
    catalog_hash_path = metadata_dir / f"remote_catalog_{library}.sha256"
    complete_marker = metadata_dir / f".mtd_{library}_cache_complete"

    previous = load_catalog(catalog_path)
    current = {entry.name: entry for entry in entries}
    validate_catalog_size(len(current), len(previous), minimum_count)

    current_names = set(current)
    local_names = {
        path.name
        for path in local_dir.glob("*.gz")
        if path.is_file()
    }

    obsolete = sorted(local_names - current_names)
    changed = sorted(
        name
        for name in current_names & set(previous)
        if previous[name].signature != current[name].signature
    )

    if not previous and not trust_existing_without_catalog:
        changed = sorted(current_names & local_names)

    integrity_state = load_integrity_state(integrity_path)
    integrity_candidates: list[str] = []
    for name in sorted(current_names & local_names):
        path = local_dir / name
        stat_result = path.stat()
        current_state = (stat_result.st_size, stat_result.st_mtime_ns)
        if full_gzip_check or integrity_state.get(name) != current_state:
            integrity_candidates.append(name)

    invalid: list[str] = []
    if integrity_candidates:
        log("INFO", f"Testing gzip integrity for {len(integrity_candidates)} cached file(s).")
        with concurrent.futures.ThreadPoolExecutor(max_workers=gzip_jobs) as executor:
            checks = {executor.submit(gzip_ok, local_dir / name): name for name in integrity_candidates}
            for future in concurrent.futures.as_completed(checks):
                name = checks[future]
                if not future.result():
                    invalid.append(name)

    missing = sorted(current_names - local_names)
    to_download_names = sorted(set(missing) | set(changed) | set(invalid))

    log("INFO", f"{library}: remote catalog entries: {len(current_names)}")
    log("INFO", f"{library}: existing local files: {len(local_names)}")
    log("INFO", f"{library}: obsolete local files: {len(obsolete)}")
    log("INFO", f"{library}: catalog-changed files: {len(changed)}")
    log("INFO", f"{library}: invalid gzip files: {len(invalid)}")
    log("INFO", f"{library}: files to download: {len(to_download_names)}")

    if check_only:
        if obsolete or to_download_names:
            log(
                "INFO",
                f"{library}: check-only found cache changes; no files were modified.",
            )
        else:
            log("OK", f"{library}: check-only found the cache current.")
        return

    # The remote catalog has already passed sanity checks, so stale local files
    # may now be removed safely.
    atomic_write_text(obsolete_path, "\n".join(obsolete) + ("\n" if obsolete else ""))
    remove_paths(local_dir / name for name in obsolete)
    remove_paths(local_dir / name for name in set(changed) | set(invalid))

    failures = download_entries(
        [current[name] for name in to_download_names],
        local_dir,
        download_jobs,
        retries,
    )
    if failures:
        failure_text = "\n".join(f"{name}\t{error}" for name, error in failures) + "\n"
        atomic_write_text(failed_path, failure_text)
    else:
        atomic_write_text(failed_path, "")

    final_missing = sorted(
        name
        for name in current_names
        if not (local_dir / name).is_file() or (local_dir / name).stat().st_size <= 0
    )

    if final_missing:
        log("WARN", f"{library}: {len(final_missing)} remote files are still missing.")
        if require_complete:
            fail(
                f"{library} cache is incomplete. First missing file: {final_missing[0]}. "
                f"See {failed_path}"
            )
    else:
        log("OK", f"{library}: local cache matches the current NCBI catalog.")

    catalog_data = catalog_text(entries)
    catalog_hash = sha256_bytes(catalog_data.encode("utf-8"))
    previous_hash = catalog_hash_path.read_text(encoding="utf-8").strip() if catalog_hash_path.is_file() else ""

    atomic_write_text(catalog_path, catalog_data)
    atomic_write_text(names_path, "\n".join(entry.name for entry in entries) + "\n")
    atomic_write_text(urls_path, "\n".join(entry.url for entry in entries) + "\n")
    atomic_write_text(catalog_hash_path, catalog_hash + "\n")
    if summary_data is not None:
        atomic_write(metadata_dir / f"assembly_summary_{library}.txt", summary_data)

    write_integrity_state(integrity_path, local_dir, current_names)

    if combined_fasta is not None:
        should_build = (
            force_combined
            or not combined_fasta.is_file()
            or bool(to_download_names)
            or bool(obsolete)
            or previous_hash != catalog_hash
        )
        if should_build:
            if final_missing:
                fail("Cannot build the combined FASTA from an incomplete collection.")
            combined_fasta.parent.mkdir(parents=True, exist_ok=True)
            temp = combined_fasta.with_name(f".{combined_fasta.name}.tmp.{os.getpid()}")
            temp.unlink(missing_ok=True)
            log("INFO", f"Building combined FASTA: {combined_fasta}")
            try:
                with temp.open("wb") as output:
                    for index, entry in enumerate(entries, start=1):
                        source = local_dir / entry.name
                        log("INFO", f"[{index}/{len(entries)}] Combining {entry.name}")
                        with gzip.open(source, "rb") as input_handle:
                            shutil.copyfileobj(input_handle, output, length=8 * 1024 * 1024)
                    output.flush()
                    os.fsync(output.fileno())
                if temp.stat().st_size <= 0:
                    fail("Combined FASTA is empty.")
                with temp.open("rb") as handle:
                    if handle.read(1) != b">":
                        fail("Combined FASTA does not begin with a FASTA header.")
                os.replace(temp, combined_fasta)
            finally:
                temp.unlink(missing_ok=True)
            log("OK", f"Combined FASTA rebuilt: {combined_fasta}")
        else:
            log("OK", f"Combined FASTA is current: {combined_fasta}")

    if not final_missing:
        marker = (
            f"status=complete\n"
            f"library={library}\n"
            f"remote_entries={len(entries)}\n"
            f"catalog_sha256={catalog_hash}\n"
            f"validated_at={time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}\n"
        )
        atomic_write_text(complete_marker, marker)


def command_assemblies(args: argparse.Namespace) -> None:
    summary_data, _ = request_bytes(args.summary_url)
    levels = set(args.level or [])
    entries = parse_assembly_summary(summary_data, args.summary_url, levels)
    synchronize(
        library=args.library,
        entries=entries,
        local_dir=Path(args.local_dir),
        metadata_dir=Path(args.metadata_dir),
        summary_data=summary_data,
        minimum_count=args.min_count,
        require_complete=args.require_complete,
        download_jobs=args.download_jobs,
        gzip_jobs=args.gzip_jobs,
        retries=args.retries,
        full_gzip_check=args.full_gzip_check,
        trust_existing_without_catalog=True,
        combined_fasta=Path(args.combined_fasta) if args.combined_fasta else None,
        force_combined=args.force_combined,
        check_only=args.check_only,
    )


def command_release(args: argparse.Namespace) -> None:
    pattern = re.compile(args.pattern)
    release, entries = parse_release_directory(
        args.base_url,
        args.release_number_url,
        pattern,
        args.metadata_jobs,
    )
    metadata_dir = Path(args.metadata_dir or args.local_dir)
    synchronize(
        library=args.library,
        entries=entries,
        local_dir=Path(args.local_dir),
        metadata_dir=metadata_dir,
        summary_data=None,
        minimum_count=args.min_count,
        require_complete=args.require_complete,
        download_jobs=args.download_jobs,
        gzip_jobs=args.gzip_jobs,
        retries=args.retries,
        full_gzip_check=args.full_gzip_check,
        # Fixed release-volume filenames are reused between RefSeq releases.
        # Without a trusted prior catalog, all existing files must be refreshed.
        trust_existing_without_catalog=False,
        check_only=args.check_only,
    )
    if not args.check_only:
        atomic_write_text(metadata_dir / f"remote_release_{args.library}.txt", release + "\n")
    log("OK", f"{args.library}: synchronized against RefSeq release {release}.")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Synchronize MTD Explorer NCBI caches against current remote catalogs."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--library", required=True)
    common.add_argument("--local-dir", required=True)
    common.add_argument("--metadata-dir")
    common.add_argument("--min-count", type=int, default=1)
    common.add_argument("--download-jobs", type=int, default=int(os.environ.get("DOWNLOAD_JOBS", "4")))
    common.add_argument("--gzip-jobs", type=int, default=int(os.environ.get("GZIP_CHECK_JOBS", "4")))
    common.add_argument("--retries", type=int, default=int(os.environ.get("DOWNLOAD_ATTEMPTS", "5")))
    common.add_argument("--full-gzip-check", action="store_true")
    common.add_argument("--require-complete", action="store_true")
    common.add_argument(
        "--check-only",
        action="store_true",
        help="Compare the live remote catalog with the cache without modifying files.",
    )

    assemblies = subparsers.add_parser("assemblies", parents=[common])
    assemblies.add_argument("--summary-url", required=True)
    assemblies.add_argument("--level", action="append")
    assemblies.add_argument("--combined-fasta")
    assemblies.add_argument("--force-combined", action="store_true")
    assemblies.set_defaults(func=command_assemblies)

    release = subparsers.add_parser("release", parents=[common])
    release.add_argument("--base-url", required=True)
    release.add_argument("--release-number-url", required=True)
    release.add_argument("--pattern", required=True)
    release.add_argument("--metadata-jobs", type=int, default=8)
    release.set_defaults(func=command_release)

    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    for name in ("download_jobs", "gzip_jobs", "retries"):
        value = getattr(args, name)
        if value < 1:
            fail(f"--{name.replace('_', '-')} must be a positive integer")
    if args.min_count < 1:
        fail("--min-count must be a positive integer")

    if args.metadata_dir is None:
        args.metadata_dir = args.local_dir

    args.func(args)


if __name__ == "__main__":
    main()
