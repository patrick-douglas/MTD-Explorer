#!/usr/bin/env python3
"""
Resolve an MTD Explorer host reference from the curated HostSpecies.csv.

The program emits shell-safe variable assignments intended to be sourced by
Create_custom_host.sh.

Exit codes
----------
0: TaxID found and metadata emitted
3: TaxID not found
4: TaxID found, but required curated reference metadata is incomplete/invalid
"""

from __future__ import annotations

import argparse
import csv
import re
import shlex
import sys
from pathlib import Path
from urllib.parse import urlparse


VERSION = "1.0.0"

REQUIRED_COLUMNS = {
    "Taxon_ID",
    "Scientific_name",
    "Reference_Taxon_ID",
    "Reference_Scientific_name",
    "Ensembl_name",
    "Ensembl_division",
    "Ensembl_release",
    "Assembly",
    "Genome_URL",
    "GTF_URL",
    "Pep_URL",
    "Reference_status",
}

SHELL_FIELDS = {
    "CSV_REFERENCE_FOUND": "1",
    "CSV_TAXID": "",
    "CSV_SCIENTIFIC_NAME": "",
    "CSV_REFERENCE_TAXID": "",
    "CSV_REFERENCE_SCIENTIFIC_NAME": "",
    "CSV_ENSEMBL_NAME": "",
    "CSV_ENSEMBL_DIVISION": "",
    "CSV_ENSEMBL_RELEASE": "",
    "CSV_ASSEMBLY": "",
    "CSV_GENOME_URL": "",
    "CSV_GTF_URL": "",
    "CSV_PEP_URL": "",
    "CSV_REFERENCE_STATUS": "",
}


def clean(value: object) -> str:
    return re.sub(r"\s+", " ", str(value or "").strip())


def is_http_or_ftp_url(value: str) -> bool:
    parsed = urlparse(value)
    return parsed.scheme in {"http", "https", "ftp"} and bool(parsed.netloc)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Resolve fixed genome/GTF/protein URLs for one TaxID from a "
            "curated MTD Explorer HostSpecies.csv."
        )
    )
    parser.add_argument("--csv", required=True, help="Curated HostSpecies.csv")
    parser.add_argument("--taxid", required=True, help="Requested NCBI TaxID")
    parser.add_argument(
        "--allow-noncomplete",
        action="store_true",
        help="Allow a Reference_status that does not begin with COMPLETE",
    )
    parser.add_argument("--version", action="version", version=VERSION)
    return parser.parse_args()


def shell_assignment(name: str, value: str) -> str:
    return f"{name}={shlex.quote(value)}"


def main() -> int:
    args = parse_args()

    taxid = clean(args.taxid)
    if not re.fullmatch(r"[1-9][0-9]*", taxid):
        print(f"[ERROR] Invalid TaxID: {taxid!r}", file=sys.stderr)
        return 4

    csv_path = Path(args.csv).expanduser().resolve()
    if not csv_path.is_file():
        print(f"[ERROR] HostSpecies.csv not found: {csv_path}", file=sys.stderr)
        return 4

    with csv_path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = set(reader.fieldnames or [])

        missing_columns = sorted(REQUIRED_COLUMNS - fieldnames)
        if missing_columns:
            print(
                "[ERROR] HostSpecies.csv is not the curated version; "
                "missing column(s): " + ", ".join(missing_columns),
                file=sys.stderr,
            )
            return 4

        matches = [
            {key: clean(value) for key, value in row.items() if key is not None}
            for row in reader
            if clean(row.get("Taxon_ID")) == taxid
        ]

    if not matches:
        return 3

    if len(matches) != 1:
        print(
            f"[ERROR] TaxID {taxid} occurs {len(matches)} times in "
            f"{csv_path}; expected exactly one row.",
            file=sys.stderr,
        )
        return 4

    row = matches[0]
    status = row["Reference_status"]

    if not args.allow_noncomplete and not status.startswith("COMPLETE"):
        print(
            f"[ERROR] TaxID {taxid} has Reference_status={status!r}; "
            "automatic reference resolution requires COMPLETE or "
            "COMPLETE_OVERRIDE.",
            file=sys.stderr,
        )
        return 4

    required_values = {
        "Scientific_name": row["Scientific_name"],
        "Reference_Taxon_ID": row["Reference_Taxon_ID"],
        "Reference_Scientific_name": row["Reference_Scientific_name"],
        "Ensembl_division": row["Ensembl_division"],
        "Ensembl_release": row["Ensembl_release"],
        "Assembly": row["Assembly"],
        "Genome_URL": row["Genome_URL"],
        "GTF_URL": row["GTF_URL"],
        "Pep_URL": row["Pep_URL"],
    }

    empty = sorted(name for name, value in required_values.items() if not value)
    if empty:
        print(
            f"[ERROR] TaxID {taxid} has incomplete automatic reference "
            "metadata; empty field(s): " + ", ".join(empty),
            file=sys.stderr,
        )
        return 4

    for field in ("Genome_URL", "GTF_URL", "Pep_URL"):
        if not is_http_or_ftp_url(row[field]):
            print(
                f"[ERROR] TaxID {taxid} has an invalid {field}: {row[field]}",
                file=sys.stderr,
            )
            return 4

    reference_taxid = row["Reference_Taxon_ID"]
    if not re.fullmatch(r"[1-9][0-9]*", reference_taxid):
        print(
            f"[ERROR] TaxID {taxid} has an invalid Reference_Taxon_ID: "
            f"{reference_taxid!r}",
            file=sys.stderr,
        )
        return 4

    values = dict(SHELL_FIELDS)
    values.update(
        {
            "CSV_TAXID": taxid,
            "CSV_SCIENTIFIC_NAME": row["Scientific_name"],
            "CSV_REFERENCE_TAXID": reference_taxid,
            "CSV_REFERENCE_SCIENTIFIC_NAME": row[
                "Reference_Scientific_name"
            ],
            "CSV_ENSEMBL_NAME": row["Ensembl_name"],
            "CSV_ENSEMBL_DIVISION": row["Ensembl_division"],
            "CSV_ENSEMBL_RELEASE": row["Ensembl_release"],
            "CSV_ASSEMBLY": row["Assembly"],
            "CSV_GENOME_URL": row["Genome_URL"],
            "CSV_GTF_URL": row["GTF_URL"],
            "CSV_PEP_URL": row["Pep_URL"],
            "CSV_REFERENCE_STATUS": status,
        }
    )

    for name, value in values.items():
        print(shell_assignment(name, value))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
