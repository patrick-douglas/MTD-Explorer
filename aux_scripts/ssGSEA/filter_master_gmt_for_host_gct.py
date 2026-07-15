#!/usr/bin/env python3

from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path
from typing import TextIO


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Filter a persistent reference master GMT against an MTD host.gct."
        )
    )

    parser.add_argument("--taxid", required=True)
    parser.add_argument("--master-gmt", required=True, type=Path)
    parser.add_argument("--host-gct", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)

    parser.add_argument("--min-genes", type=int, default=5)
    parser.add_argument("--max-genes", type=int, default=500)
    parser.add_argument("--min-overlap-genes", type=int, default=100)
    parser.add_argument("--min-overlap-pct", type=float, default=1.0)

    return parser.parse_args()


def require_nonempty_file(path: Path, label: str) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"{label} was not found or is empty: {path}")


def open_text(path: Path) -> TextIO:
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")

    return path.open("r", encoding="utf-8", errors="replace")


def read_host_gct_genes(path: Path) -> set[str]:
    genes: set[str] = set()

    with open_text(path) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            # Standard GCT:
            # line 1 = #1.2
            # line 2 = dimensions
            # line 3 = header
            if line_number <= 3:
                continue

            line = raw_line.rstrip("\r\n")

            if not line:
                continue

            gene_id = line.split("\t", 1)[0].strip()

            if gene_id:
                genes.add(gene_id)

    if not genes:
        raise RuntimeError(f"No genes were extracted from host.gct: {path}")

    return genes


def main() -> int:
    args = parse_args()

    require_nonempty_file(args.master_gmt, "Master GMT")
    require_nonempty_file(args.host_gct, "host.gct")

    if args.min_genes < 1:
        raise RuntimeError("--min-genes must be at least 1")

    if args.max_genes < args.min_genes:
        raise RuntimeError("--max-genes must be >= --min-genes")

    if args.min_overlap_genes < 1:
        raise RuntimeError("--min-overlap-genes must be at least 1")

    if args.min_overlap_pct < 0:
        raise RuntimeError("--min-overlap-pct must be >= 0")

    args.outdir.mkdir(parents=True, exist_ok=True)

    host_genes = read_host_gct_genes(args.host_gct)

    master_rows: list[tuple[str, str, list[str]]] = []
    master_gene_union: set[str] = set()

    with args.master_gmt.open(
        "r",
        encoding="utf-8",
        errors="replace",
    ) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\r\n")

            if not line:
                continue

            fields = line.split("\t")

            if len(fields) < 3:
                raise RuntimeError(
                    f"Invalid GMT row at line {line_number}: fewer than 3 columns"
                )

            term_id = fields[0].strip()
            description = fields[1].strip() or "eggNOG_GO"
            genes = sorted({gene.strip() for gene in fields[2:] if gene.strip()})

            if not term_id or not genes:
                continue

            master_rows.append((term_id, description, genes))
            master_gene_union.update(genes)

    if not master_rows:
        raise RuntimeError("No valid gene sets were read from the master GMT")

    overlapping_genes = host_genes.intersection(master_gene_union)
    overlap_count = len(overlapping_genes)
    overlap_pct = 100.0 * overlap_count / len(host_genes)

    if (
        overlap_count < args.min_overlap_genes
        or overlap_pct < args.min_overlap_pct
    ):
        raise RuntimeError(
            "Master GMT and host.gct compatibility is too low. "
            f"host.gct genes={len(host_genes)}, "
            f"overlap={overlap_count}, "
            f"overlap_pct={overlap_pct:.4f}, "
            f"minimum_genes={args.min_overlap_genes}, "
            f"minimum_pct={args.min_overlap_pct}"
        )

    prefix = f"custom_taxid_{args.taxid}_eggNOG_GO"

    out_gmt = args.outdir / f"{prefix}.gmt"
    out_summary = args.outdir / f"{prefix}_analysis_summary.tsv"
    out_set_table = args.outdir / f"{prefix}_analysis_set_table.tsv"

    kept_rows: list[tuple[str, str, list[str]]] = []
    set_table_rows: list[tuple[str, int, int, str]] = []

    dropped_too_small = 0
    dropped_too_large = 0

    for term_id, description, master_genes in master_rows:
        analysis_genes = sorted(set(master_genes).intersection(host_genes))
        analysis_size = len(analysis_genes)

        if analysis_size < args.min_genes:
            status = "dropped_too_small"
            dropped_too_small += 1
        elif analysis_size > args.max_genes:
            status = "dropped_too_large"
            dropped_too_large += 1
        else:
            status = "kept"
            kept_rows.append((term_id, description, analysis_genes))

        set_table_rows.append(
            (
                term_id,
                len(master_genes),
                analysis_size,
                status,
            )
        )

    if not kept_rows:
        raise RuntimeError(
            "No gene sets remained after host.gct intersection and "
            f"{args.min_genes}–{args.max_genes} gene filtering"
        )

    with out_gmt.open("w", encoding="utf-8", newline="\n") as handle:
        for term_id, description, genes in kept_rows:
            handle.write(
                "\t".join([term_id, description, *genes]) + "\n"
            )

    with out_set_table.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(
            "gene_set\tmaster_gene_count\t"
            "analysis_gene_count\tstatus\n"
        )

        for term_id, master_size, analysis_size, status in set_table_rows:
            handle.write(
                f"{term_id}\t{master_size}\t{analysis_size}\t{status}\n"
            )

    with out_summary.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("metric\tvalue\n")
        handle.write(f"taxid\t{args.taxid}\n")
        handle.write(f"host_gct_genes\t{len(host_genes)}\n")
        handle.write(f"master_gmt_genes\t{len(master_gene_union)}\n")
        handle.write(f"overlap_genes\t{overlap_count}\n")
        handle.write(f"overlap_pct\t{overlap_pct:.4f}\n")
        handle.write(f"master_gene_sets\t{len(master_rows)}\n")
        handle.write(f"kept_gene_sets\t{len(kept_rows)}\n")
        handle.write(f"dropped_too_small\t{dropped_too_small}\n")
        handle.write(f"dropped_too_large\t{dropped_too_large}\n")
        handle.write(f"minimum_set_size\t{args.min_genes}\n")
        handle.write(f"maximum_set_size\t{args.max_genes}\n")

    for output in (out_gmt, out_summary, out_set_table):
        require_nonempty_file(output, "Filtered ssGSEA resource")

    print("[OK] Analysis-specific ssGSEA GMT created")
    print(f"[OK] GMT:             {out_gmt}")
    print(f"[OK] Summary:         {out_summary}")
    print(f"[OK] Set table:       {out_set_table}")
    print(f"[INFO] host.gct genes:{len(host_genes)}")
    print(f"[INFO] GMT overlap:   {overlap_count} ({overlap_pct:.4f}%)")
    print(f"[INFO] Kept sets:     {len(kept_rows)}")

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        raise SystemExit(1)
