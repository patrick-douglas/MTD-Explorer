#!/usr/bin/env python3
"""
Ensure that one TaxID exists in MTD Explorer HostSpecies.csv.

Existing rows are never rewritten by this helper. If the TaxID is absent, a
minimal MANUAL_CUSTOM row is appended using the supplied scientific name.
A timestamped backup is created before modifying the CSV.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import os
import re
import shutil
import sys
from pathlib import Path


VERSION = "1.0.0"


def clean(value: object) -> str:
    return re.sub(r"\s+", " ", str(value or "").strip())


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Check or append a manual custom-host entry in HostSpecies.csv."
        )
    )
    parser.add_argument("--csv", required=True)
    parser.add_argument("--taxid", required=True)
    parser.add_argument("--scientific-name")
    parser.add_argument(
        "--add-if-missing",
        action="store_true",
        help="Append a minimal MANUAL_CUSTOM row when TaxID is absent.",
    )
    parser.add_argument("--version", action="version", version=VERSION)
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    taxid = clean(args.taxid)
    scientific_name = clean(args.scientific_name)

    if not re.fullmatch(r"[1-9][0-9]*", taxid):
        print(f"[ERROR] Invalid TaxID: {taxid!r}", file=sys.stderr)
        return 2

    path = Path(args.csv).expanduser().resolve()
    if not path.is_file():
        print(f"[ERROR] HostSpecies.csv not found: {path}", file=sys.stderr)
        return 2

    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = list(reader.fieldnames or [])
        rows = list(reader)

    required = {"Taxon_ID", "Scientific_name"}
    missing = sorted(required - set(fieldnames))
    if missing:
        print(
            "[ERROR] HostSpecies.csv is missing required column(s): "
            + ", ".join(missing),
            file=sys.stderr,
        )
        return 2

    matches = [
        row
        for row in rows
        if clean(row.get("Taxon_ID")) == taxid
    ]

    if len(matches) > 1:
        print(
            f"[ERROR] TaxID {taxid} occurs more than once in {path}",
            file=sys.stderr,
        )
        return 2

    if matches:
        existing_name = clean(matches[0].get("Scientific_name"))
        if (
            scientific_name
            and existing_name
            and scientific_name.casefold() != existing_name.casefold()
        ):
            print(
                f"[ERROR] TaxID {taxid} is already registered as "
                f"{existing_name!r}, not {scientific_name!r}.",
                file=sys.stderr,
            )
            return 2

        print(f"HOSTSPECIES_ENTRY_STATUS=existing")
        print(f"HOSTSPECIES_ENTRY_NAME={existing_name!r}")
        return 0

    if not args.add_if_missing:
        return 3

    if not scientific_name or len(scientific_name.split()) < 2:
        print(
            "[ERROR] A binomial or trinomial --scientific-name is required "
            "to add a TaxID absent from HostSpecies.csv.",
            file=sys.stderr,
        )
        return 2

    new_row = {column: "" for column in fieldnames}
    new_row["Taxon_ID"] = taxid
    new_row["Scientific_name"] = scientific_name

    optional_values = {
        "Reference_Taxon_ID": taxid,
        "Reference_Scientific_name": scientific_name,
        "Reference_status": "MANUAL_CUSTOM",
        "MartDatasets_status": "NOT_CHECKED",
    }
    for column, value in optional_values.items():
        if column in new_row:
            new_row[column] = value

    rows.append(new_row)

    timestamp = dt.datetime.now().strftime("%Y%m%d-%H%M%S")
    backup = path.with_name(f"{path.name}.bak.manual-{timestamp}")
    shutil.copy2(path, backup)

    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)

    os.chmod(temporary, path.stat().st_mode)
    os.replace(temporary, path)

    print("HOSTSPECIES_ENTRY_STATUS=added")
    print(f"HOSTSPECIES_ENTRY_NAME={scientific_name!r}")
    print(f"HOSTSPECIES_ENTRY_BACKUP={str(backup)!r}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
