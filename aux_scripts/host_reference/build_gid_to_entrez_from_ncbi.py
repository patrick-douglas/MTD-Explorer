#!/usr/bin/env python3

"""Build a persistent Ensembl gene ID to NCBI GeneID mapping."""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import re
import sys
from pathlib import Path


def normalize_gene_id(value: str) -> str:
    value = str(value).strip()
    return re.sub(r"\.\d+$", "", value)


def read_gtf_gene_ids(path: Path) -> set[str]:
    genes: set[str] = set()

    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            gene = normalize_gene_id(line)

            if gene:
                genes.add(gene)

    return genes


def main() -> None:
    parser = argparse.ArgumentParser()

    parser.add_argument("--gene2ensembl", required=True)
    parser.add_argument("--reference-taxid", required=True)
    parser.add_argument("--requested-taxid", required=True)
    parser.add_argument("--gtf-genes", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--summary", required=True)

    args = parser.parse_args()

    source = Path(args.gene2ensembl)
    gtf_file = Path(args.gtf_genes)
    output = Path(args.output)
    summary = Path(args.summary)

    if not source.is_file() or source.stat().st_size == 0:
        raise SystemExit(
            f"[ERROR] gene2ensembl cache is missing or empty: {source}"
        )

    if not gtf_file.is_file() or gtf_file.stat().st_size == 0:
        raise SystemExit(
            f"[ERROR] GTF gene-ID list is missing or empty: {gtf_file}"
        )

    gtf_genes = read_gtf_gene_ids(gtf_file)

    if not gtf_genes:
        raise SystemExit("[ERROR] No gene IDs were read from the GTF list.")

    all_reference_pairs: set[tuple[str, str]] = set()
    reference_matched_pairs: set[tuple[str, str]] = set()

    with gzip.open(
        source,
        mode="rt",
        encoding="utf-8",
        errors="replace",
        newline="",
    ) as handle:
        reader = csv.reader(handle, delimiter="\t")

        try:
            raw_header = next(reader)
        except StopIteration:
            raise SystemExit("[ERROR] gene2ensembl.gz is empty.")

        header = [
            column.lstrip("#").strip()
            for column in raw_header
        ]

        required = {
            "tax_id",
            "GeneID",
            "Ensembl_gene_identifier",
        }

        missing = required.difference(header)

        if missing:
            raise SystemExit(
                "[ERROR] Missing gene2ensembl columns: "
                + ", ".join(sorted(missing))
            )

        tax_index = header.index("tax_id")
        entrez_index = header.index("GeneID")
        ensembl_index = header.index("Ensembl_gene_identifier")

        maximum_index = max(
            tax_index,
            entrez_index,
            ensembl_index,
        )

        for row in reader:
            if len(row) <= maximum_index:
                continue

            if row[tax_index].strip() != args.reference_taxid:
                continue

            ensembl = normalize_gene_id(row[ensembl_index])
            entrez = row[entrez_index].strip()

            if not ensembl or ensembl == "-":
                continue

            if not entrez or entrez == "-":
                continue

            pair = (ensembl, entrez)
            all_reference_pairs.add(pair)

            if ensembl in gtf_genes:
                reference_matched_pairs.add(pair)

    if not all_reference_pairs:
        raise SystemExit(
            "[WARNING] No Ensembl/GeneID mappings were found for "
            f"reference TaxID {args.reference_taxid}."
        )

    if not reference_matched_pairs:
        raise SystemExit(
            "[WARNING] NCBI mappings were found, but none matched "
            "the gene IDs from the selected GTF."
        )

    mapped_gtf_genes = {
        ensembl
        for ensembl, _ in reference_matched_pairs
    }

    coverage = 100.0 * len(mapped_gtf_genes) / len(gtf_genes)

    output.parent.mkdir(parents=True, exist_ok=True)
    summary.parent.mkdir(parents=True, exist_ok=True)

    output_tmp = output.with_suffix(output.suffix + ".tmp")
    summary_tmp = summary.with_suffix(summary.suffix + ".tmp")

    with output_tmp.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.writer(handle, delimiter="\t")

        writer.writerow(
            [
                "ensembl_gene_id",
                "external_gene_name",
                "entrezgene_id",
            ]
        )

        for ensembl, entrez in sorted(reference_matched_pairs):
            writer.writerow([ensembl, "", entrez])

    with summary_tmp.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.writer(handle, delimiter="\t")

        writer.writerow(["parameter", "value"])
        writer.writerow(["requested_taxid", args.requested_taxid])
        writer.writerow(["reference_taxid", args.reference_taxid])
        writer.writerow(["gtf_unique_genes", len(gtf_genes)])
        writer.writerow(
            ["ncbi_reference_pairs", len(all_reference_pairs)]
        )
        writer.writerow(
            ["reference_matched_pairs", len(reference_matched_pairs)]
        )
        writer.writerow(
            ["mapped_gtf_genes", len(mapped_gtf_genes)]
        )
        writer.writerow(
            ["gtf_mapping_coverage_pct", f"{coverage:.2f}"]
        )
        writer.writerow(["source", str(source)])

    os.replace(output_tmp, output)
    os.replace(summary_tmp, summary)

    print("[OK] Reference-matched Ensembl -> Entrez mapping created")
    print(f"  Requested TaxID: {args.requested_taxid}")
    print(f"  Reference TaxID: {args.reference_taxid}")
    print(f"  GTF genes:       {len(gtf_genes)}")
    print(f"  Mapped genes:    {len(mapped_gtf_genes)}")
    print(f"  Coverage:        {coverage:.2f}%")
    print(f"  Mapping:         {output}")
    print(f"  Summary:         {summary}")


if __name__ == "__main__":
    main()
