#!/usr/bin/env python3
"""Build a sequence accession-to-TaxID map for the cached RefSeq viral FASTAs."""

from __future__ import annotations

import argparse
import csv
import gzip
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterator, Set, Tuple

ACCESSION_PATTERN = (
    r"(?:"
    r"[A-Z]{1,4}_\d{5,12}"
    r"|"
    r"[A-Z]{1,6}\d{5,12}"
    r")"
    r"(?:\.\d+)?"
)
ACCESSION_RE = re.compile(
    rf"(?<![A-Za-z0-9_])({ACCESSION_PATTERN})(?![A-Za-z0-9])"
)
BRACKETED_ACCESSION_RE = re.compile(rf"\[({ACCESSION_PATTERN})\]")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--assembly-summary", required=True, type=Path)
    parser.add_argument("--fasta-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--conflicts", required=True, type=Path)
    return parser.parse_args()


def filename_from_ftp_path(ftp_path: str) -> str:
    ftp_path = ftp_path.rstrip("/")
    assembly = ftp_path.rsplit("/", 1)[-1]
    for suffix in ("_genomic.fna.gz", ".fna.gz", "_genomic"):
        if assembly.endswith(suffix):
            assembly = assembly[: -len(suffix)]
    return f"{assembly}_genomic.fna.gz" if assembly else ""


def extract_accession(header: str) -> str:
    match = BRACKETED_ACCESSION_RE.search(header)
    if match is None:
        token = header.split(None, 1)[0].strip("|")
        match = ACCESSION_RE.fullmatch(token)
    if match is None:
        match = ACCESSION_RE.search(header)
    return "" if match is None else match.group(1)


def fasta_headers(path: Path) -> Iterator[str]:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                yield line[1:].strip()


def main() -> int:
    args = parse_args()
    if not args.assembly_summary.is_file() or args.assembly_summary.stat().st_size == 0:
        raise FileNotFoundError(f"Missing assembly summary: {args.assembly_summary}")
    if not args.fasta_dir.is_dir():
        raise FileNotFoundError(f"Missing viral FASTA directory: {args.fasta_dir}")

    file_taxids: Dict[str, int] = {}
    with args.assembly_summary.open("rt", encoding="utf-8", errors="replace") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if not row or row[0].startswith("#") or len(row) < 20:
                continue
            taxid_text = row[5].strip()
            ftp_path = row[19].strip()
            if not taxid_text.isdigit() or ftp_path == "na":
                continue
            filename = filename_from_ftp_path(ftp_path)
            if filename:
                file_taxids[filename] = int(taxid_text)

    mapping: Dict[str, int] = {}
    conflicts: Dict[str, Set[int]] = defaultdict(set)
    files_seen = 0
    records_seen = 0
    missing_summary_files = []
    records_without_accession = 0

    for fasta_path in sorted(args.fasta_dir.glob("*.gz")):
        taxid = file_taxids.get(fasta_path.name)
        if taxid is None:
            missing_summary_files.append(fasta_path.name)
            continue
        files_seen += 1
        for header in fasta_headers(fasta_path):
            records_seen += 1
            accession = extract_accession(header)
            if not accession:
                records_without_accession += 1
                continue
            base = accession.split(".", 1)[0]
            for key in {accession, base}:
                previous = mapping.get(key)
                if previous is None:
                    mapping[key] = taxid
                elif previous != taxid:
                    conflicts[key].update({previous, taxid})

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.conflicts.parent.mkdir(parents=True, exist_ok=True)

    # Ambiguous assembly-derived mappings are not safe fallbacks. Remove them
    # from the usable map; the main helper may still resolve them through the
    # official NCBI accession2taxid files.
    for accession in conflicts:
        mapping.pop(accession, None)

    with args.output.open("wt", encoding="utf-8") as handle:
        handle.write("accession\ttaxid\n")
        for accession in sorted(mapping):
            handle.write(f"{accession}\t{mapping[accession]}\n")

    with args.conflicts.open("wt", encoding="utf-8") as handle:
        handle.write("accession\ttaxids\n")
        for accession in sorted(conflicts):
            values = ",".join(str(value) for value in sorted(conflicts[accession]))
            handle.write(f"{accession}\t{values}\n")

    print(f"Assembly entries indexed: {len(file_taxids)}")
    print(f"Cached FASTA files scanned: {files_seen}")
    print(f"Sequence records scanned: {records_seen}")
    print(f"Mappings written: {len(mapping)}")
    print(f"Conflicts found: {len(conflicts)}")
    print(f"Records without accession: {records_without_accession}")
    print(f"Cached files absent from summary: {len(missing_summary_files)}")
    print(f"Map:       {args.output}")
    print(f"Conflicts: {args.conflicts}")

    if missing_summary_files:
        for filename in missing_summary_files[:20]:
            print(f"[ERROR] Cached FASTA absent from current assembly summary: {filename}", file=sys.stderr)
        return 1
    if records_without_accession:
        print("[ERROR] Some RefSeq FASTA headers had no recognizable accession.", file=sys.stderr)
        return 1
    if conflicts:
        print(
            "[WARNING] Ambiguous assembly-derived mappings were excluded; "
            "the main helper will fall back to NCBI accession2taxid.",
            file=sys.stderr,
        )
    if not mapping:
        print("[ERROR] No RefSeq viral accession mappings were generated.", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        raise SystemExit(1)
