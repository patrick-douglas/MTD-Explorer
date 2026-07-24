#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
from collections import defaultdict
from pathlib import Path


def normalize(token: str) -> str:
    token = token.strip().lstrip(">")
    return re.sub(r"\.\d+$", "", token)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Report the Kraken2 libraries containing unmapped accessions."
    )
    parser.add_argument("--database", required=True, type=Path)
    parser.add_argument("--unmapped", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    database = args.database.resolve()
    unmapped_path = args.unmapped.resolve()
    output = args.output.resolve()

    if not unmapped_path.is_file():
        raise SystemExit(f"[ERROR] Unmapped file not found: {unmapped_path}")

    unmapped = {
        normalize(line)
        for line in unmapped_path.read_text(
            encoding="utf-8", errors="replace"
        ).splitlines()
        if line.strip()
    }

    memberships: dict[str, set[str]] = defaultdict(set)
    counts: dict[str, set[str]] = defaultdict(set)

    for prelim in sorted((database / "library").glob("*/prelim_map.txt")):
        library = prelim.parent.name

        with prelim.open(encoding="utf-8", errors="replace") as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")

                for field in fields:
                    for token in re.split(r"[|\s,;]+", field):
                        candidate = normalize(token)

                        if candidate in unmapped:
                            memberships[candidate].add(library)
                            counts[library].add(candidate)

    output.parent.mkdir(parents=True, exist_ok=True)

    with output.open("w", encoding="utf-8", newline="") as handle:
        handle.write("accession\tlibraries\n")

        for accession in sorted(unmapped):
            libraries = ",".join(sorted(memberships.get(accession, set())))
            handle.write(f"{accession}\t{libraries}\n")

    print(f"[MTD-KRAKEN2] Unmapped accessions: {len(unmapped):,}")

    for library in sorted(counts):
        print(
            f"[MTD-KRAKEN2]   {library}: "
            f"{len(counts[library]):,}"
        )

    not_located = sum(1 for accession in unmapped if accession not in memberships)
    print(f"[MTD-KRAKEN2] Not located in prelim_map files: {not_located:,}")
    print(f"[MTD-KRAKEN2] Origin report: {output}")


if __name__ == "__main__":
    main()
