#!/usr/bin/env python3
"""Build a taxid-aware, nonredundant viral FASTA for Kraken2."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import os
import re
import sys
import tempfile
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Set, TextIO, Tuple

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
BRACKETED_ACCESSION_RE = re.compile(
    rf"\[({ACCESSION_PATTERN})\]"
)
KRAKEN_TAXID_RE = re.compile(r"(?:^|\|)kraken:taxid\|(\d+)(?:\||$)")
COMPLEMENT = str.maketrans(
    "ACGTRYKMSWBDHVN",
    "TGCAYRMKSWVHDBN",
)


@dataclass(frozen=True)
class RecordIdentity:
    source: str
    accession: str
    accession_base: str
    taxid: Optional[int]
    header: str


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("rt", encoding="utf-8", errors="replace")


def read_fasta(path: Path) -> Iterator[Tuple[str, str]]:
    header: Optional[str] = None
    sequence_parts: List[str] = []

    with open_text(path) as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, normalize_sequence("".join(sequence_parts))
                header = line[1:].strip()
                sequence_parts = []
            else:
                if header is None:
                    raise ValueError(f"Sequence found before first FASTA header in {path}")
                sequence_parts.append(line)

    if header is not None:
        yield header, normalize_sequence("".join(sequence_parts))


def normalize_sequence(sequence: str) -> str:
    return sequence.upper().replace("U", "T")


def canonical_sha256(sequence: str) -> str:
    forward = hashlib.sha256(sequence.encode("ascii")).hexdigest()
    reverse_complement = sequence.translate(COMPLEMENT)[::-1]
    reverse = hashlib.sha256(reverse_complement.encode("ascii")).hexdigest()
    return min(forward, reverse)


def extract_accession(header: str) -> Tuple[str, str]:
    """Extract the sequence accession from a FASTA header.

    Virus-Host DB may use an old GenBank locus name as the first token and put
    the real accession in brackets, for example ``PSU05771 [U05771]``.
    Prefer the bracketed accession when present.
    """
    match = BRACKETED_ACCESSION_RE.search(header)
    if match is None:
        first_token = header.split(None, 1)[0].strip("|")
        match = ACCESSION_RE.fullmatch(first_token)
    if match is None:
        match = ACCESSION_RE.search(header)
    if match is None:
        return "", ""
    versioned = match.group(1)
    base = versioned.split(".", 1)[0]
    return versioned, base


def embedded_taxid(header: str) -> Optional[int]:
    match = KRAKEN_TAXID_RE.search(header)
    if match is None:
        return None
    return int(match.group(1))


def collect_accessions(paths: Sequence[Path]) -> Tuple[Set[str], Set[str]]:
    versioned: Set[str] = set()
    base: Set[str] = set()
    for path in paths:
        for header, _sequence in read_fasta(path):
            accession, accession_base = extract_accession(header)
            if accession:
                versioned.add(accession)
                base.add(accession_base)
    return versioned, base


def load_accession_taxids(
    taxonomy_dir: Path,
    wanted_versioned: Set[str],
    wanted_base: Set[str],
) -> Dict[str, int]:
    mapping: Dict[str, int] = {}
    mapping_files = [
        taxonomy_dir / "nucl_gb.accession2taxid",
        taxonomy_dir / "nucl_wgs.accession2taxid",
    ]

    for mapping_file in mapping_files:
        if not mapping_file.is_file() or mapping_file.stat().st_size == 0:
            raise FileNotFoundError(f"Missing or empty taxonomy mapping: {mapping_file}")

        with mapping_file.open("rt", encoding="utf-8", errors="replace") as handle:
            for line_number, line in enumerate(handle, start=1):
                if line_number == 1 and line.startswith("accession\t"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 3:
                    continue
                accession_base, accession_version, taxid_text = fields[:3]
                if (
                    accession_version not in wanted_versioned
                    and accession_base not in wanted_base
                ):
                    continue
                try:
                    taxid = int(taxid_text)
                except ValueError:
                    continue
                mapping[accession_version] = taxid
                mapping.setdefault(accession_base, taxid)

    return mapping


def load_simple_accession_taxid_map(
    path: Optional[Path],
    wanted_versioned: Set[str],
    wanted_base: Set[str],
) -> Dict[str, int]:
    mapping: Dict[str, int] = {}
    if path is None:
        return mapping
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Missing or empty custom accession map: {path}")

    with path.open("rt", encoding="utf-8", errors="replace") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if line_number == 1 and fields[:2] == ["accession", "taxid"]:
                continue
            if len(fields) < 2:
                continue
            accession = fields[0].strip()
            taxid_text = fields[1].strip()
            if not accession:
                continue
            accession_base = accession.split(".", 1)[0]
            if accession not in wanted_versioned and accession_base not in wanted_base:
                continue
            try:
                taxid = int(taxid_text)
            except ValueError:
                continue
            mapping[accession] = taxid
            mapping.setdefault(accession_base, taxid)

    return mapping


def resolve_taxid(
    header: str,
    accession: str,
    accession_base: str,
    mapping: Dict[str, int],
    source_mapping: Optional[Dict[str, int]] = None,
    prefer_source_mapping: bool = False,
) -> Optional[int]:
    taxid = embedded_taxid(header)
    if taxid is not None:
        return taxid
    if not accession:
        return None

    if prefer_source_mapping and source_mapping:
        taxid = source_mapping.get(accession) or source_mapping.get(accession_base)
        if taxid is not None:
            return taxid

    taxid = mapping.get(accession) or mapping.get(accession_base)
    if taxid is not None:
        return taxid

    if source_mapping:
        return source_mapping.get(accession) or source_mapping.get(accession_base)
    return None


def kraken_header(original_header: str, accession: str, taxid: Optional[int]) -> str:
    if taxid is None or embedded_taxid(original_header) is not None:
        return original_header

    first_token, separator, remainder = original_header.partition(" ")
    sequence_id = accession or first_token
    rewritten = f"kraken:taxid|{taxid}|{sequence_id}"
    if separator and remainder:
        rewritten += f" {remainder}"
    return rewritten


def write_record(handle: TextIO, header: str, sequence: str) -> None:
    handle.write(f">{header}\n")
    for start in range(0, len(sequence), 80):
        handle.write(sequence[start : start + 80] + "\n")


def duplicate_decision(
    existing_records: Sequence[RecordIdentity],
    current: RecordIdentity,
) -> Tuple[bool, str, Optional[RecordIdentity]]:
    """Return keep, reason, matching record.

    Exact sequences are removed only when they resolve to the same TaxID. When
    TaxIDs conflict, or cannot be resolved for different accessions, both records
    are preserved so Kraken2 can retain conservative LCA behavior.
    """
    for existing in existing_records:
        if current.taxid is not None and existing.taxid is not None:
            if current.taxid == existing.taxid:
                return False, "duplicate_same_taxid", existing
            continue

        if (
            current.accession_base
            and existing.accession_base
            and current.accession_base == existing.accession_base
        ):
            return False, "duplicate_same_accession_unresolved_taxid", existing

    if existing_records:
        return True, "retained_taxid_conflict_or_unresolved", existing_records[0]
    return True, "unique_sequence", None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Combine RefSeq viral, Virus-Host DB, and optional extra FASTAs while "
            "removing only exact sequence duplicates that share a TaxID."
        )
    )
    parser.add_argument("--primary", required=True, type=Path, help="Primary viral FASTA")
    parser.add_argument("--virushost", required=True, type=Path, help="Virus-Host DB FASTA")
    parser.add_argument(
        "--extra",
        action="append",
        default=[],
        type=Path,
        help="Additional FASTA; may be supplied multiple times",
    )
    parser.add_argument("--taxonomy-dir", required=True, type=Path)
    parser.add_argument(
        "--primary-taxid-map",
        type=Path,
        help="Optional two-column accession-to-TaxID map for the primary RefSeq FASTA",
    )
    parser.add_argument(
        "--virushost-taxid-map",
        type=Path,
        help="Optional two-column accession-to-TaxID map for Virus-Host DB",
    )
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--details", required=True, type=Path)
    return parser.parse_args()


def validate_input(path: Path, label: str) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"{label} FASTA is missing or empty: {path}")


def main() -> int:
    args = parse_args()

    validate_input(args.primary, "Primary")
    validate_input(args.virushost, "Virus-Host DB")
    for extra in args.extra:
        validate_input(extra, "Extra")

    sources: List[Tuple[str, Path]] = [
        ("refseq_viral", args.primary),
        ("virushostdb", args.virushost),
    ]
    sources.extend((f"extra_{index}", path) for index, path in enumerate(args.extra, start=1))

    all_paths = [path for _label, path in sources]
    wanted_versioned, wanted_base = collect_accessions(all_paths)
    accession_taxids = load_accession_taxids(
        args.taxonomy_dir,
        wanted_versioned,
        wanted_base,
    )
    primary_taxids = load_simple_accession_taxid_map(
        args.primary_taxid_map,
        wanted_versioned,
        wanted_base,
    )
    virushost_taxids = load_simple_accession_taxid_map(
        args.virushost_taxid_map,
        wanted_versioned,
        wanted_base,
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    args.details.parent.mkdir(parents=True, exist_ok=True)

    counts: Counter[str] = Counter()
    seen: Dict[str, List[RecordIdentity]] = defaultdict(list)

    output_fd, output_tmp_name = tempfile.mkstemp(
        prefix=f".{args.output.name}.", suffix=".tmp", dir=args.output.parent
    )
    os.close(output_fd)
    output_tmp = Path(output_tmp_name)

    details_fd, details_tmp_name = tempfile.mkstemp(
        prefix=f".{args.details.name}.", suffix=".tmp", dir=args.details.parent
    )
    os.close(details_fd)
    details_tmp = Path(details_tmp_name)

    try:
        with output_tmp.open("wt", encoding="utf-8") as output_handle, details_tmp.open(
            "wt", encoding="utf-8"
        ) as details_handle:
            details_handle.write(
                "source\tstatus\tlength\taccession\ttaxid\tcanonical_sha256\t"
                "header\tmatched_source\tmatched_accession\tmatched_taxid\tmatched_header\n"
            )

            for source_label, fasta_path in sources:
                for header, sequence in read_fasta(fasta_path):
                    counts["records_seen"] += 1
                    counts[f"{source_label}_records_seen"] += 1

                    if not sequence:
                        counts["empty_sequences_skipped"] += 1
                        continue

                    accession, accession_base = extract_accession(header)
                    if source_label == "refseq_viral":
                        source_mapping = primary_taxids
                        prefer_source_mapping = True
                    elif source_label == "virushostdb":
                        source_mapping = virushost_taxids
                        prefer_source_mapping = True
                    else:
                        source_mapping = None
                        prefer_source_mapping = False
                    taxid = resolve_taxid(
                        header,
                        accession,
                        accession_base,
                        accession_taxids,
                        source_mapping,
                        prefer_source_mapping,
                    )
                    if not accession:
                        counts["records_without_accession"] += 1
                    if taxid is None:
                        counts["records_without_taxid"] += 1

                    digest = canonical_sha256(sequence)
                    identity = RecordIdentity(
                        source=source_label,
                        accession=accession,
                        accession_base=accession_base,
                        taxid=taxid,
                        header=header,
                    )
                    keep, status, matched = duplicate_decision(seen[digest], identity)
                    counts[status] += 1
                    counts[f"{source_label}_{status}"] += 1

                    if keep:
                        output_header = kraken_header(header, accession, taxid)
                        write_record(output_handle, output_header, sequence)
                        seen[digest].append(identity)
                        counts["records_written"] += 1
                        counts[f"{source_label}_records_written"] += 1
                    else:
                        counts["records_removed"] += 1
                        counts[f"{source_label}_records_removed"] += 1

                    details_handle.write(
                        "\t".join(
                            [
                                source_label,
                                status,
                                str(len(sequence)),
                                accession,
                                "" if taxid is None else str(taxid),
                                digest,
                                header.replace("\t", " "),
                                "" if matched is None else matched.source,
                                "" if matched is None else matched.accession,
                                "" if matched is None or matched.taxid is None else str(matched.taxid),
                                "" if matched is None else matched.header.replace("\t", " "),
                            ]
                        )
                        + "\n"
                    )

        if counts["records_written"] == 0 or output_tmp.stat().st_size == 0:
            raise RuntimeError("No viral sequences were written to the combined FASTA")

        os.replace(output_tmp, args.output)
        os.replace(details_tmp, args.details)

        with args.summary.open("wt", encoding="utf-8") as summary_handle:
            summary_handle.write("metric\tvalue\n")
            ordered_metrics = [
                "records_seen",
                "records_written",
                "records_removed",
                "unique_sequence",
                "duplicate_same_taxid",
                "duplicate_same_accession_unresolved_taxid",
                "retained_taxid_conflict_or_unresolved",
                "records_without_accession",
                "records_without_taxid",
                "empty_sequences_skipped",
            ]
            source_metrics: List[str] = []
            for source_label, _path in sources:
                source_metrics.extend(
                    [
                        f"{source_label}_records_seen",
                        f"{source_label}_records_written",
                        f"{source_label}_records_removed",
                        f"{source_label}_duplicate_same_taxid",
                        f"{source_label}_duplicate_same_accession_unresolved_taxid",
                        f"{source_label}_retained_taxid_conflict_or_unresolved",
                    ]
                )
            for metric in ordered_metrics + source_metrics:
                summary_handle.write(f"{metric}\t{counts[metric]}\n")

        print("============================================================")
        print("NONREDUNDANT VIRAL FASTA")
        print("============================================================")
        print(f"Records examined:  {counts['records_seen']}")
        print(f"Records written:   {counts['records_written']}")
        print(f"Duplicates removed:{counts['records_removed']}")
        print(
            "Conflicting/unresolved exact duplicates retained: "
            f"{counts['retained_taxid_conflict_or_unresolved']}"
        )
        print(f"Records without TaxID: {counts['records_without_taxid']}")
        print(f"Output:  {args.output}")
        print(f"Summary: {args.summary}")
        print(f"Details: {args.details}")

        if counts["records_without_taxid"]:
            print(
                "[ERROR] Some records could not be assigned a TaxID. "
                "The combined FASTA must not be used for Kraken2 until this is resolved.",
                file=sys.stderr,
            )
            return 2

        return 0
    except Exception:
        output_tmp.unlink(missing_ok=True)
        details_tmp.unlink(missing_ok=True)
        raise


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        raise SystemExit(1)
