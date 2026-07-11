#!/usr/bin/env python3
"""
Update the Common_name and kegg columns of MTD Explorer HostSpecies.csv.

Data sources
------------
Common_name:
    Ensembl REST API /info/species, matched by NCBI Taxonomy ID.

kegg:
    KEGG REST API /list/genome/<taxid>, with an exact scientific-name
    fallback against /list/genome when the TaxID is not a supported rank.

Safety
------
- Existing values are preserved by default.
- Use --replace-existing to replace existing values when a verified value is found.
- Missing or ambiguous online matches never erase an existing value.
- --in-place creates a timestamped backup before replacing HostSpecies.csv.
- A TSV audit report records every decision.

Only Python 3 standard-library modules are required.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import os
import re
import shutil
import socket
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable, Optional, Sequence


SCRIPT_VERSION = "1.2.0"
ENSEMBL_REST = "https://rest.ensembl.org"
KEGG_REST = "https://rest.kegg.jp"

ENSEMBL_DIVISIONS = (
    "EnsemblVertebrates",
    "EnsemblMetazoa",
    "EnsemblFungi",
    "EnsemblPlants",
    "EnsemblProtists",
)

REQUIRED_COLUMNS = {
    "Taxon_ID",
    "Scientific_name",
    "Common_name",
    "kegg",
}

REPORT_COLUMNS = (
    "Taxon_ID",
    "Scientific_name",
    "Common_name_old",
    "Common_name_verified",
    "Common_name_final",
    "Common_name_status",
    "Common_name_source",
    "kegg_old",
    "kegg_verified",
    "kegg_final",
    "kegg_status",
    "kegg_candidates",
    "changed",
    "notes",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Verify and update Common_name and kegg in MTD Explorer "
            "HostSpecies.csv using Ensembl and KEGG TaxID mappings."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "-i",
        "--hostspecies",
        default="HostSpecies.csv",
        help="Input HostSpecies.csv.",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="HostSpecies.updated.csv",
        help="Output CSV when --in-place is not used.",
    )
    parser.add_argument(
        "--report",
        default="HostSpecies.common_kegg_report.tsv",
        help="Audit report path.",
    )
    parser.add_argument(
        "--cache-dir",
        default=str(Path.home() / ".cache" / "mtd_explorer" / "hostspecies_metadata"),
        help="Directory used to cache Ensembl and KEGG responses.",
    )
    parser.add_argument(
        "--cache-max-age-days",
        type=int,
        default=7,
        help="Reuse cached responses up to this age. Use 0 to refresh.",
    )
    parser.add_argument(
        "--replace-existing",
        action="store_true",
        help=(
            "Replace non-empty Common_name/kegg values when an exact, verified "
            "TaxID match is available. Without this option, only empty fields "
            "are filled."
        ),
    )
    parser.add_argument(
        "--in-place",
        action="store_true",
        help="Replace the input CSV after creating a timestamped backup.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Fetch, compare, and write the report, but do not write an updated CSV.",
    )
    parser.add_argument(
        "--taxid",
        action="append",
        default=[],
        help="Process only this TaxID. May be repeated.",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=90,
        help="Network timeout per request, in seconds.",
    )
    parser.add_argument(
        "--retries",
        type=int,
        default=5,
        help="Number of attempts for each request.",
    )
    parser.add_argument(
        "--request-delay",
        type=float,
        default=0.40,
        help="Delay between uncached requests, in seconds (KEGG limit: at most 3 requests/second).",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=SCRIPT_VERSION,
    )

    args = parser.parse_args()

    if args.timeout < 1:
        parser.error("--timeout must be at least 1")
    if args.retries < 1:
        parser.error("--retries must be at least 1")
    if args.cache_max_age_days < 0:
        parser.error("--cache-max-age-days cannot be negative")
    if args.request_delay < 0:
        parser.error("--request-delay cannot be negative")
    if args.in_place and args.output != "HostSpecies.updated.csv":
        parser.error("--output cannot be combined with --in-place")

    return args


def say(message: str) -> None:
    print(message, flush=True)


def normalize_spaces(value: Any) -> str:
    return re.sub(r"\s+", " ", str(value or "").strip())


def normalize_scientific_name(value: Any) -> str:
    text = normalize_spaces(value).lower()
    text = re.sub(r"[^a-z0-9]+", " ", text)
    return normalize_spaces(text)


def display_common_name(value: Any) -> str:
    """
    Convert the verified Ensembl common name into the whitespace-free format
    required by HostSpecies.csv.

    Examples:
        "house mouse"       -> "House_mouse"
        "little brown bat"  -> "Little_brown_bat"
        "Mallard"           -> "Mallard"

    All whitespace and punctuation separators are normalized to underscores.
    """
    text = normalize_spaces(str(value or "").replace("_", " "))
    if not text:
        return ""

    tokens = re.findall(r"[A-Za-z0-9]+", text)
    if not tokens:
        return ""

    result = "_".join(tokens)
    return result[:1].upper() + result[1:]


def cache_is_fresh(path: Path, max_age_days: int) -> bool:
    if not path.is_file() or path.stat().st_size == 0:
        return False
    if max_age_days == 0:
        return False
    age_seconds = time.time() - path.stat().st_mtime
    return age_seconds <= max_age_days * 86400


class Downloader:
    def __init__(
        self,
        cache_dir: Path,
        timeout: int,
        retries: int,
        request_delay: float,
        cache_max_age_days: int,
    ) -> None:
        self.cache_dir = cache_dir
        self.timeout = timeout
        self.retries = retries
        self.request_delay = request_delay
        self.cache_max_age_days = cache_max_age_days
        self.user_agent = f"MTD-Explorer-HostSpecies-metadata/{SCRIPT_VERSION}"
        self.cache_dir.mkdir(parents=True, exist_ok=True)

    def get_bytes(
        self,
        url: str,
        cache_name: str,
        *,
        accept: str = "*/*",
    ) -> bytes:
        cache_path = self.cache_dir / cache_name

        if cache_is_fresh(cache_path, self.cache_max_age_days):
            return cache_path.read_bytes()

        headers = {
            "User-Agent": self.user_agent,
            "Accept": accept,
        }

        last_error: Optional[BaseException] = None
        for attempt in range(1, self.retries + 1):
            try:
                request = urllib.request.Request(url, headers=headers)
                with urllib.request.urlopen(
                    request,
                    timeout=self.timeout,
                ) as response:
                    payload = response.read()

                if not payload:
                    raise RuntimeError("server returned an empty response")

                temporary = cache_path.with_suffix(cache_path.suffix + ".tmp")
                temporary.write_bytes(payload)
                os.replace(temporary, cache_path)

                if self.request_delay:
                    time.sleep(self.request_delay)

                return payload

            except urllib.error.HTTPError as exc:
                last_error = exc

                # Do not retry permanent request errors.
                if exc.code in {400, 401, 403, 404, 405, 410}:
                    break

                if attempt < self.retries:
                    retry_after = exc.headers.get("Retry-After")
                    if retry_after and retry_after.isdigit():
                        delay = min(120, int(retry_after))
                    else:
                        delay = min(30, 2 ** (attempt - 1))
                    time.sleep(delay)

            except (
                urllib.error.URLError,
                TimeoutError,
                socket.timeout,
                OSError,
            ) as exc:
                last_error = exc
                if attempt < self.retries:
                    time.sleep(min(30, 2 ** (attempt - 1)))

        if cache_path.is_file() and cache_path.stat().st_size > 0:
            say(
                f"[WARN] Network refresh failed for {url}; "
                f"using stale cache {cache_path}"
            )
            return cache_path.read_bytes()

        raise RuntimeError(
            f"Unable to retrieve {url}: "
            f"{last_error}"
        )

    def get_json(self, url: str, cache_name: str) -> dict[str, Any]:
        payload = self.get_bytes(
            url,
            cache_name,
            accept="application/json",
        )
        data = json.loads(payload.decode("utf-8"))
        if not isinstance(data, dict):
            raise RuntimeError(f"Expected a JSON object from {url}")
        return data

    def get_text(self, url: str, cache_name: str) -> str:
        payload = self.get_bytes(
            url,
            cache_name,
            accept="text/plain",
        )
        return payload.decode("utf-8", errors="replace")


def read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file():
        raise RuntimeError(f"HostSpecies.csv not found: {path}")

    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = list(reader.fieldnames or [])
        missing = REQUIRED_COLUMNS.difference(fieldnames)
        if missing:
            raise RuntimeError(
                "HostSpecies.csv is missing required column(s): "
                + ", ".join(sorted(missing))
            )
        rows = [dict(row) for row in reader]

    seen: set[str] = set()
    duplicates: list[str] = []

    for row_number, row in enumerate(rows, start=2):
        taxid = normalize_spaces(row.get("Taxon_ID"))
        scientific_name = normalize_spaces(row.get("Scientific_name"))

        if not taxid:
            raise RuntimeError(f"Missing Taxon_ID on CSV line {row_number}")
        if not taxid.isdigit():
            raise RuntimeError(
                f"Non-numeric Taxon_ID {taxid!r} on CSV line {row_number}"
            )
        if not scientific_name:
            raise RuntimeError(
                f"Missing Scientific_name on CSV line {row_number}"
            )
        if taxid in seen:
            duplicates.append(taxid)
        seen.add(taxid)

    if duplicates:
        raise RuntimeError(
            "Duplicate Taxon_ID values in HostSpecies.csv: "
            + ", ".join(sorted(set(duplicates), key=int))
        )

    return fieldnames, rows


def write_csv(
    path: Path,
    fieldnames: Sequence[str],
    rows: Iterable[dict[str, str]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")

    with temporary.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)

    os.replace(temporary, path)


def write_report(path: Path, rows: Sequence[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")

    with temporary.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=REPORT_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)

    os.replace(temporary, path)


def fetch_ensembl_species(
    downloader: Downloader,
) -> dict[str, list[dict[str, Any]]]:
    by_taxid: dict[str, list[dict[str, Any]]] = defaultdict(list)
    errors: list[str] = []

    for division in ENSEMBL_DIVISIONS:
        query = urllib.parse.urlencode({"division": division})
        url = f"{ENSEMBL_REST}/info/species?{query}"
        cache_name = f"ensembl_info_species_{division}.json"

        try:
            payload = downloader.get_json(url, cache_name)
        except Exception as exc:
            errors.append(f"{division}: {exc}")
            continue

        for entry in payload.get("species", []):
            if not isinstance(entry, dict):
                continue
            taxid = normalize_spaces(entry.get("taxon_id"))
            if not taxid:
                continue
            item = dict(entry)
            item.setdefault("division", division)
            by_taxid[taxid].append(item)

    if not by_taxid:
        raise RuntimeError(
            "No Ensembl species metadata could be retrieved. "
            + " | ".join(errors)
        )

    if errors:
        say("[WARN] Some Ensembl divisions could not be retrieved:")
        for error in errors:
            say(f"       {error}")

    return dict(by_taxid)


def choose_ensembl_entry(
    scientific_name: str,
    mart_dataset: str,
    entries: Sequence[dict[str, Any]],
) -> tuple[Optional[dict[str, Any]], str]:
    if not entries:
        return None, "NO_TAXID_MATCH"

    scientific_norm = normalize_scientific_name(scientific_name)
    dataset_prefix = (
        normalize_spaces(mart_dataset)
        .lower()
        .split("_gene_ensembl", 1)[0]
    )

    scored: list[tuple[int, str, dict[str, Any]]] = []

    for entry in entries:
        name = normalize_scientific_name(entry.get("name"))
        display = normalize_scientific_name(entry.get("display_name"))
        aliases = {
            normalize_scientific_name(alias)
            for alias in (entry.get("aliases") or [])
            if normalize_scientific_name(alias)
        }

        score = 0
        if display == scientific_norm:
            score += 100
        if name == scientific_norm:
            score += 95
        if scientific_norm in aliases:
            score += 90

        compact_name = name.replace(" ", "")
        if dataset_prefix and dataset_prefix in compact_name:
            score += 20

        if normalize_spaces(entry.get("common_name")):
            score += 5
        if normalize_spaces(entry.get("assembly")):
            score += 2

        scored.append((score, name, entry))

    scored.sort(key=lambda item: (-item[0], item[1]))
    best_score, _, best = scored[0]

    if len(scored) > 1 and scored[1][0] == best_score:
        best_common = normalize_spaces(best.get("common_name"))
        second_common = normalize_spaces(scored[1][2].get("common_name"))
        if best_common != second_common:
            return None, "AMBIGUOUS_TAXID_MATCH"

    if best_score == 0 and len(entries) > 1:
        return None, "AMBIGUOUS_TAXID_MATCH"

    return best, "MATCHED_BY_TAXID"


def strip_kegg_name_qualifiers(name: str) -> str:
    """
    Remove KEGG common-name/strain qualifiers while keeping the leading
    scientific name.

    Examples:
        "Homo sapiens (human)" -> "Homo sapiens"
        "Escherichia coli K-12 MG1655" remains unchanged here; exact-name
        matching is intentionally conservative.
    """
    value = normalize_spaces(name)
    value = re.split(r"\s+\(", value, maxsplit=1)[0]
    return normalize_spaces(value)


def parse_kegg_genome_list(
    genome_text: str,
    *,
    source: str = "KEGG_GENOME_LIST",
) -> list[dict[str, str]]:
    """
    Parse the official KEGG /list/genome response.

    Typical records:
        T01001<TAB>hsa; Homo sapiens (human)
        gn:T01001<TAB>hsa; Homo sapiens (human)

    Returned organism codes are copied exactly from KEGG and normalized only
    to lowercase.
    """
    records: list[dict[str, str]] = []

    for raw_line in genome_text.splitlines():
        line = raw_line.strip()
        if not line:
            continue

        fields = line.split("\t", 1)
        if len(fields) != 2:
            continue

        genome_field, description = fields
        genome_match = re.search(r"\bT\d+\b", genome_field)
        description_match = re.match(
            r"\s*([a-z][a-z0-9]{1,4})\s*;\s*(.+?)\s*$",
            description,
            flags=re.IGNORECASE,
        )

        if not genome_match or not description_match:
            continue

        genome_id = genome_match.group(0)
        organism_code = description_match.group(1).lower()
        kegg_name = normalize_spaces(description_match.group(2))
        base_name = strip_kegg_name_qualifiers(kegg_name)

        records.append(
            {
                "genome_id": genome_id,
                "code": organism_code,
                "scientific_name": kegg_name,
                "base_scientific_name": base_name,
                "source": source,
            }
        )

    return records


def index_kegg_by_exact_name(
    records: Sequence[dict[str, str]],
) -> dict[str, list[dict[str, str]]]:
    by_name: dict[str, list[dict[str, str]]] = defaultdict(list)

    for record in records:
        normalized = normalize_scientific_name(
            record.get("base_scientific_name", "")
        )
        if normalized:
            by_name[normalized].append(dict(record))

    return dict(by_name)


def fetch_kegg_mapping(
    downloader: Downloader,
    species_by_taxid: dict[str, str],
) -> dict[str, list[dict[str, str]]]:
    """
    Resolve exact KEGG organism codes.

    Primary method:
        GET /list/genome/<NCBI TaxID>

    KEGG added Taxonomy-ID filtering to the genome list operation in 2026.
    Some CSV entries may represent subspecies or another rank not accepted by
    that endpoint. For those cases only, use an exact scientific-name match
    against the complete official /list/genome catalog.

    A failure for one species does not abort all remaining species. A total
    failure to retrieve the global KEGG catalog remains fatal.
    """
    global_text = downloader.get_text(
        f"{KEGG_REST}/list/genome",
        "kegg_list_genome.tsv",
    )
    global_records = parse_kegg_genome_list(
        global_text,
        source="KEGG_EXACT_NAME",
    )

    if not global_records:
        raise RuntimeError("KEGG /list/genome returned no usable entries")

    by_exact_name = index_kegg_by_exact_name(global_records)
    result: dict[str, list[dict[str, str]]] = {}
    direct_successes = 0
    direct_no_matches = 0
    direct_errors: list[str] = []

    for taxid in sorted(species_by_taxid, key=int):
        scientific_name = species_by_taxid[taxid]
        url = f"{KEGG_REST}/list/genome/{taxid}"
        cache_name = f"kegg_list_genome_taxid_{taxid}.tsv"

        direct_records: list[dict[str, str]] = []

        try:
            payload = downloader.get_text(url, cache_name)
            direct_records = parse_kegg_genome_list(
                payload,
                source="KEGG_TAXID",
            )
        except RuntimeError as exc:
            message = str(exc)

            # KEGG may return 400/404 or an empty body when a TaxID is not a
            # supported taxonomy rank or is absent from KEGG.
            expected_no_match = (
                "HTTP Error 400" in message
                or "HTTP Error 404" in message
                or "empty response" in message.lower()
            )

            if expected_no_match:
                direct_no_matches += 1
            else:
                direct_errors.append(f"{taxid}: {message}")

            # Keep public-API traffic under the documented request limit even
            # when a permanent error returns immediately.
            if downloader.request_delay:
                time.sleep(downloader.request_delay)

        if direct_records:
            result[taxid] = direct_records
            direct_successes += 1
            continue

        normalized_name = normalize_scientific_name(scientific_name)
        fallback_records = [
            {
                **record,
                "source": "KEGG_EXACT_NAME",
            }
            for record in by_exact_name.get(normalized_name, [])
        ]

        if fallback_records:
            result[taxid] = fallback_records
        else:
            direct_no_matches += 1

    say(
        "[INFO] KEGG resolution: "
        f"{direct_successes} TaxID match(es), "
        f"{sum(1 for values in result.values() if values and values[0].get('source') == 'KEGG_EXACT_NAME')} "
        "exact-name fallback match(es), "
        f"{len(species_by_taxid) - len(result)} unresolved."
    )

    if direct_errors:
        say(
            "[WARN] Some KEGG TaxID requests had network/server errors; "
            "exact-name fallback was attempted:"
        )
        for error in direct_errors[:10]:
            say(f"       {error}")
        if len(direct_errors) > 10:
            say(f"       ... and {len(direct_errors) - 10} more")

    return result

def strip_kegg_qualifiers(name: str) -> str:
    """
    KEGG names may include a strain/assembly description after a comma,
    semicolon, or opening parenthesis. Keep the leading taxonomic name for
    conservative matching.
    """
    value = normalize_spaces(name)
    value = re.split(r"\s*[;,]\s*|\s+\(", value, maxsplit=1)[0]
    return normalize_spaces(value)


def choose_kegg_candidate(
    scientific_name: str,
    candidates: Sequence[dict[str, str]],
) -> tuple[str, str, str]:
    if not candidates:
        return "", "NO_KEGG_MATCH", ""

    descriptions = "; ".join(
        (
            f"{item['code']}|{item['genome_id']}|"
            f"{item['scientific_name']}|{item.get('source', '')}"
        )
        for item in candidates
    )

    target = normalize_scientific_name(scientific_name)

    exact = [
        item
        for item in candidates
        if normalize_scientific_name(
            item.get("base_scientific_name")
            or strip_kegg_name_qualifiers(item["scientific_name"])
        )
        == target
    ]

    if len(exact) == 1:
        selected = exact[0]
    elif len(candidates) == 1:
        selected = candidates[0]
    else:
        # Multiple genomes/codes may represent strains or assemblies under one
        # TaxID. Do not guess which code the MTD analysis should use.
        return "", "AMBIGUOUS_KEGG_MATCH", descriptions

    code = normalize_spaces(selected.get("code")).lower()

    # KEGG cellular-organism codes are short lower-case alphanumeric prefixes.
    # Reject anything outside that syntax instead of writing malformed values.
    if not re.fullmatch(r"[a-z][a-z0-9]{1,4}", code):
        return "", "INVALID_KEGG_CODE_REJECTED", descriptions

    source = selected.get("source", "")
    if source == "KEGG_TAXID":
        status = "MATCHED_BY_KEGG_TAXID"
    elif source == "KEGG_EXACT_NAME":
        status = "MATCHED_BY_EXACT_KEGG_NAME"
    else:
        status = "MATCHED_BY_KEGG_CATALOG"

    return code, status, descriptions

def choose_final_value(
    old: str,
    verified: str,
    status: str,
    replace_existing: bool,
) -> tuple[str, str]:
    old = normalize_spaces(old)
    verified = normalize_spaces(verified)

    if not verified:
        if old:
            return old, f"{status}_PRESERVED_EXISTING"
        return "", status

    if not old:
        return verified, "FILLED_VERIFIED"

    if old == verified:
        return old, "ALREADY_VERIFIED"

    if replace_existing:
        return verified, "REPLACED_WITH_VERIFIED"

    return old, "DIFFERS_PRESERVED_EXISTING"


def make_backup(path: Path) -> Path:
    timestamp = dt.datetime.now().strftime("%Y%m%d-%H%M%S")
    backup = path.with_name(f"{path.name}.bak.{timestamp}")
    shutil.copy2(path, backup)
    return backup


def main() -> int:
    args = parse_args()

    input_path = Path(args.hostspecies).expanduser().resolve()
    output_path = Path(args.output).expanduser()
    report_path = Path(args.report).expanduser()
    cache_dir = Path(args.cache_dir).expanduser()

    requested_taxids = {
        normalize_spaces(value)
        for value in args.taxid
        if normalize_spaces(value)
    }
    invalid_requested = sorted(
        value for value in requested_taxids if not value.isdigit()
    )
    if invalid_requested:
        raise RuntimeError(
            "Invalid --taxid value(s): " + ", ".join(invalid_requested)
        )

    fieldnames, rows = read_csv(input_path)
    csv_taxids = {normalize_spaces(row["Taxon_ID"]) for row in rows}

    missing_requested = requested_taxids.difference(csv_taxids)
    if missing_requested:
        raise RuntimeError(
            "Requested TaxID(s) not found in HostSpecies.csv: "
            + ", ".join(sorted(missing_requested, key=int))
        )

    downloader = Downloader(
        cache_dir=cache_dir,
        timeout=args.timeout,
        retries=args.retries,
        request_delay=args.request_delay,
        cache_max_age_days=args.cache_max_age_days,
    )

    say("============================================================")
    say("MTD Explorer — HostSpecies common-name/KEGG updater")
    say(f"Input:             {input_path}")
    say(
        "Update mode:       "
        + ("replace verified values" if args.replace_existing else "fill empty fields")
    )
    say(f"In-place:          {'yes' if args.in_place else 'no'}")
    say(f"Dry run:           {'yes' if args.dry_run else 'no'}")
    say(f"Cache directory:   {cache_dir}")
    say("============================================================")

    say("[INFO] Retrieving Ensembl species metadata...")
    ensembl_by_taxid = fetch_ensembl_species(downloader)

    say("[INFO] Retrieving exact KEGG organism codes...")
    species_for_kegg = {
        normalize_spaces(row["Taxon_ID"]): normalize_spaces(
            row["Scientific_name"]
        )
        for row in rows
        if (
            not requested_taxids
            or normalize_spaces(row["Taxon_ID"]) in requested_taxids
        )
    }
    kegg_by_taxid = fetch_kegg_mapping(
        downloader,
        species_for_kegg,
    )

    reports: list[dict[str, str]] = []
    changed_rows = 0
    common_verified_count = 0
    kegg_verified_count = 0
    common_unresolved = 0
    kegg_unresolved = 0

    for row in rows:
        taxid = normalize_spaces(row.get("Taxon_ID"))
        scientific_name = normalize_spaces(row.get("Scientific_name"))
        mart_dataset = normalize_spaces(row.get("MartDatasets"))

        if requested_taxids and taxid not in requested_taxids:
            continue

        notes: list[str] = []

        ensembl_entry, common_match_status = choose_ensembl_entry(
            scientific_name,
            mart_dataset,
            ensembl_by_taxid.get(taxid, []),
        )

        common_verified = ""
        common_source = ""

        if ensembl_entry is not None:
            common_verified = display_common_name(
                ensembl_entry.get("common_name")
            )
            if common_verified:
                common_source = (
                    "Ensembl:"
                    + normalize_spaces(ensembl_entry.get("division"))
                )
                common_verified_count += 1
            else:
                common_match_status = "MATCHED_TAXID_BUT_NO_COMMON_NAME"

        if not common_verified:
            common_unresolved += 1

        common_old = normalize_spaces(row.get("Common_name"))
        common_final, common_status = choose_final_value(
            common_old,
            common_verified,
            common_match_status,
            args.replace_existing,
        )

        kegg_verified, kegg_match_status, kegg_candidates = (
            choose_kegg_candidate(
                scientific_name,
                kegg_by_taxid.get(taxid, []),
            )
        )

        if kegg_verified:
            kegg_verified_count += 1
        else:
            kegg_unresolved += 1

        kegg_old = normalize_spaces(row.get("kegg"))
        kegg_final, kegg_status = choose_final_value(
            kegg_old,
            kegg_verified,
            kegg_match_status,
            args.replace_existing,
        )

        changed = (
            common_final != common_old
            or kegg_final != kegg_old
        )

        if changed:
            changed_rows += 1

        row["Common_name"] = common_final
        row["kegg"] = kegg_final

        if common_verified and common_old and common_old != common_verified:
            notes.append("Common_name differs from verified Ensembl value")
        if kegg_verified and kegg_old and kegg_old != kegg_verified:
            notes.append("kegg differs from verified KEGG TaxID mapping")

        reports.append(
            {
                "Taxon_ID": taxid,
                "Scientific_name": scientific_name,
                "Common_name_old": common_old,
                "Common_name_verified": common_verified,
                "Common_name_final": common_final,
                "Common_name_status": common_status,
                "Common_name_source": common_source,
                "kegg_old": kegg_old,
                "kegg_verified": kegg_verified,
                "kegg_final": kegg_final,
                "kegg_status": kegg_status,
                "kegg_candidates": kegg_candidates,
                "changed": "YES" if changed else "NO",
                "notes": "; ".join(notes),
            }
        )

        say(
            f"[{taxid}] {scientific_name}: "
            f"common={common_final or '-'} ({common_status}); "
            f"kegg={kegg_final or '-'} ({kegg_status})"
        )

    write_report(report_path.resolve(), reports)

    output_written: Optional[Path] = None
    backup_path: Optional[Path] = None

    if not args.dry_run:
        if args.in_place:
            backup_path = make_backup(input_path)
            write_csv(input_path, fieldnames, rows)
            output_written = input_path
        else:
            output_resolved = output_path.resolve()
            if output_resolved == input_path:
                raise RuntimeError(
                    "Output resolves to the input file. Use --in-place to "
                    "replace HostSpecies.csv safely."
                )
            write_csv(output_resolved, fieldnames, rows)
            output_written = output_resolved

    say("============================================================")
    say(f"Rows evaluated:          {len(reports)}")
    say(f"Rows changed:            {changed_rows}")
    say(f"Common names verified:   {common_verified_count}")
    say(f"Common names unresolved: {common_unresolved}")
    say(f"KEGG codes verified:     {kegg_verified_count}")
    say(f"KEGG codes unresolved:   {kegg_unresolved}")
    say(f"Audit report:            {report_path.resolve()}")

    if output_written is not None:
        say(f"Updated CSV:             {output_written}")
    else:
        say("Updated CSV:             not written (--dry-run)")

    if backup_path is not None:
        say(f"Backup:                  {backup_path}")

    say("============================================================")

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        say("\n[ERROR] Interrupted by user.")
        raise SystemExit(130)
    except Exception as exc:
        say(f"[ERROR] {exc}")
        raise SystemExit(1)
