#!/usr/bin/env python3
"""Build an accession-to-TaxID map from Virus-Host DB release metadata."""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, Set

ACCESSION_RE = re.compile(
    r"^(?:[A-Z]{1,4}_\d{5,12}|[A-Z]{1,6}\d{5,12})(?:\.\d+)?$"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--non-segmented", required=True, type=Path)
    parser.add_argument("--segmented", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--conflicts", required=True, type=Path)
    return parser.parse_args()


def read_rows(paths: Iterable[Path]):
    for path in paths:
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(f"Missing or empty Virus-Host metadata file: {path}")
        with path.open("rt", encoding="utf-8", errors="replace", newline="") as handle:
            yield from csv.reader(handle, delimiter="\t")


def main() -> int:
    args = parse_args()
    mapping: Dict[str, int] = {}
    conflicts: Dict[str, Set[int]] = defaultdict(set)
    skipped_invalid = 0

    for row in read_rows([args.non_segmented, args.segmented]):
        if len(row) < 3:
            continue
        taxid_text = row[1].strip()
        if not taxid_text.isdigit():
            continue
        taxid = int(taxid_text)
        if taxid <= 1:
            continue

        for raw_accession in row[2].split(","):
            accession = raw_accession.strip().upper()
            if not accession:
                continue
            if ACCESSION_RE.fullmatch(accession) is None:
                skipped_invalid += 1
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

    with args.output.open("wt", encoding="utf-8") as handle:
        handle.write("accession\ttaxid\n")
        for accession in sorted(mapping):
            handle.write(f"{accession}\t{mapping[accession]}\n")

    with args.conflicts.open("wt", encoding="utf-8") as handle:
        handle.write("accession\ttaxids\n")
        for accession in sorted(conflicts):
            values = ",".join(str(value) for value in sorted(conflicts[accession]))
            handle.write(f"{accession}\t{values}\n")

    print(f"Mappings written: {len(mapping)}")
    print(f"Conflicts found:  {len(conflicts)}")
    print(f"Invalid accessions skipped: {skipped_invalid}")
    print(f"Map:       {args.output}")
    print(f"Conflicts: {args.conflicts}")

    if not mapping:
        print("[ERROR] No Virus-Host accession mappings were generated.", file=sys.stderr)
        return 1
    if conflicts:
        print("[ERROR] Conflicting Virus-Host accession mappings were detected.", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        raise SystemExit(1)
