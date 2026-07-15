#!/usr/bin/env python3

from __future__ import annotations

import argparse
import re
import sys
from collections import defaultdict
from pathlib import Path


GO_PATTERN = re.compile(r"^GO:\d{7}$")
GO_COLUMN_CANDIDATES = ("GOs", "GO_terms", "GO")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build a persistent master GO GMT from reference-matched "
            "eggNOG-mapper annotations."
        )
    )

    parser.add_argument(
        "--taxid",
        required=True,
        help="NCBI Taxon ID used by the MTD custom host.",
    )

    parser.add_argument(
        "--eggnog",
        required=True,
        type=Path,
        help="Reference-matched eggNOG .emapper.annotations file.",
    )

    parser.add_argument(
        "--gtf-genes",
        required=True,
        type=Path,
        help="Text file containing one reference GTF gene ID per line.",
    )

    parser.add_argument(
        "--outdir",
        required=True,
        type=Path,
        help="Persistent functional annotation output directory.",
    )

    return parser.parse_args()


def require_nonempty_file(path: Path, label: str) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"{label} was not found or is empty: {path}")


def read_gtf_genes(path: Path) -> set[str]:
    genes: set[str] = set()

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            gene = line.strip()

            if gene:
                genes.add(gene)

    if not genes:
        raise RuntimeError(f"No gene IDs were read from: {path}")

    return genes


def parse_go_terms(raw_value: str) -> set[str]:
    raw_value = raw_value.strip()

    if not raw_value or raw_value == "-":
        return set()

    terms: set[str] = set()

    for item in re.split(r"[,;]", raw_value):
        go_id = item.strip()

        if GO_PATTERN.fullmatch(go_id):
            terms.add(go_id)

    return terms


def read_eggnog(
    path: Path,
    gtf_genes: set[str],
) -> tuple[dict[str, set[str]], dict[str, set[str]], int, int]:
    header: list[str] | None = None
    go_column: str | None = None

    gene_to_go: dict[str, set[str]] = defaultdict(set)
    go_to_gene: dict[str, set[str]] = defaultdict(set)

    data_rows = 0
    rows_with_go = 0

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\r\n")

            if line.startswith("#query"):
                header = line[1:].split("\t")

                for candidate in GO_COLUMN_CANDIDATES:
                    if candidate in header:
                        go_column = candidate
                        break

                if go_column is None:
                    raise RuntimeError(
                        "No GO annotation column was found in the eggNOG header. "
                        f"Checked: {', '.join(GO_COLUMN_CANDIDATES)}"
                    )

                continue

            if not line or line.startswith("#"):
                continue

            if header is None or go_column is None:
                continue

            data_rows += 1

            fields = line.split("\t")

            if len(fields) < len(header):
                fields.extend([""] * (len(header) - len(fields)))
            elif len(fields) > len(header):
                fields = fields[: len(header)]

            record = dict(zip(header, fields))

            gene_id = record.get("query", "").strip()

            if not gene_id or gene_id not in gtf_genes:
                continue

            go_terms = parse_go_terms(record.get(go_column, ""))

            if not go_terms:
                continue

            rows_with_go += 1
            gene_to_go[gene_id].update(go_terms)

            for go_id in go_terms:
                go_to_gene[go_id].add(gene_id)

    if not gene_to_go:
        raise RuntimeError(
            "No eggNOG genes with GO terms matched the reference GTF gene IDs."
        )

    return gene_to_go, go_to_gene, data_rows, rows_with_go


def main() -> int:
    args = parse_args()

    require_nonempty_file(args.eggnog, "eggNOG annotations")
    require_nonempty_file(args.gtf_genes, "GTF gene list")

    args.outdir.mkdir(parents=True, exist_ok=True)

    gtf_genes = read_gtf_genes(args.gtf_genes)

    gene_to_go, go_to_gene, data_rows, rows_with_go = read_eggnog(
        args.eggnog,
        gtf_genes,
    )

    prefix = f"custom_taxid_{args.taxid}_eggNOG_GO"

    out_gmt = args.outdir / f"{prefix}_master.gmt"
    out_gene_table = args.outdir / f"{prefix}_master_gene_table.tsv"
    out_summary = args.outdir / f"{prefix}_master_summary.tsv"

    # Do not apply the ssGSEA 5–500 set-size filter here.
    # A large master set may become valid after intersection with host.gct.
    with out_gmt.open("w", encoding="utf-8", newline="\n") as handle:
        for go_id in sorted(go_to_gene):
            genes = sorted(go_to_gene[go_id])

            handle.write(
                "\t".join(
                    [
                        go_id,
                        "eggNOG_GO_reference_master",
                        *genes,
                    ]
                )
                + "\n"
            )

    with out_gene_table.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("gene_id\tgo_id\n")

        for gene_id in sorted(gene_to_go):
            for go_id in sorted(gene_to_go[gene_id]):
                handle.write(f"{gene_id}\t{go_id}\n")

    annotated_gtf_genes = set(gene_to_go)
    annotated_pct = (
        100.0 * len(annotated_gtf_genes) / len(gtf_genes)
        if gtf_genes
        else 0.0
    )

    with out_summary.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("metric\tvalue\n")
        handle.write(f"taxid\t{args.taxid}\n")
        handle.write(f"eggnog_data_rows\t{data_rows}\n")
        handle.write(f"eggnog_rows_with_go\t{rows_with_go}\n")
        handle.write(f"gtf_gene_count\t{len(gtf_genes)}\n")
        handle.write(
            f"gtf_genes_with_eggnog_go\t{len(annotated_gtf_genes)}\n"
        )
        handle.write(
            f"gtf_genes_with_eggnog_go_pct\t{annotated_pct:.4f}\n"
        )
        handle.write(f"unique_go_terms\t{len(go_to_gene)}\n")
        handle.write(f"master_gmt_sets\t{len(go_to_gene)}\n")
        handle.write(
            "master_set_size_filter\t"
            "none; filtering is performed against host.gct during MTD Explorer\n"
        )

    for output in (out_gmt, out_gene_table, out_summary):
        require_nonempty_file(output, "Generated functional resource")

    print("[OK] Persistent master eggNOG/GO resources created")
    print(f"[OK] Master GMT:       {out_gmt}")
    print(f"[OK] Gene table:       {out_gene_table}")
    print(f"[OK] Summary:          {out_summary}")
    print(f"[INFO] GTF genes:      {len(gtf_genes)}")
    print(f"[INFO] Annotated genes:{len(annotated_gtf_genes)}")
    print(f"[INFO] GO terms:       {len(go_to_gene)}")

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        raise SystemExit(1)
