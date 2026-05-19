#!/usr/bin/env python3
"""
bracken_to_heatmap_table_nopandas_EN.py

Convert Bracken/MTD/DEG tables into a heatmap-ready matrix with this format:

UNIQID    NAME    sample1    sample2    sample3 ...

No pandas/numpy version. Uses only the Python standard library.

It works with:
  1) Wide Bracken/MTD matrices:
       first column with feature/taxon name + numeric sample columns
       example: bracken_species_all_normalized.csv

  2) DEG tables:
       Name + raw counts + .norm + .normtrans + DESeq2 statistics
       With --prefer auto, the script prefers .normtrans columns for heatmaps.

  3) Multiple single-sample Bracken reports:
       name/taxonomy_id + new_est_reads/fraction_total_reads etc.
       The script merges one value column from each file.

Examples:
  python3 bracken_to_heatmap_table_nopandas_EN.py \
    --input bracken_species_all_normalized.csv \
    --output heatmap_table.txt

  python3 bracken_to_heatmap_table_nopandas_EN.py \
    --input bracken_species_all_DEG.csv \
    --output heatmap_table_normtrans.txt \
    --prefer normtrans

  python3 bracken_to_heatmap_table_nopandas_EN.py \
    --input sample1.bracken sample2.bracken sample3.bracken \
    --output heatmap_from_reports.txt \
    --value-col new_est_reads

  python3 bracken_to_heatmap_table_nopandas_EN.py \
    --input bracken_species_all_DEG.csv \
    --output top50_heatmap.txt \
    --prefer normtrans \
    --top 50 \
    --transform row_zscore
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import sys
from collections import Counter, OrderedDict
from typing import Dict, List, Optional, Sequence, Tuple

ID_CANDIDATES = [
    "UNIQID", "uniqid", "id", "ID", "feature_id", "FeatureID",
    "taxid", "tax_id", "taxonomy_id", "taxonomyID", "NCBI_tax_id",
]

NAME_CANDIDATES = [
    "", "NAME", "Name", "name", "taxon", "Taxon", "taxonomy", "species",
    "Species", "clade_name", "Unnamed: 0", "#NAME", "# Name",
]

BRACKEN_VALUE_CANDIDATES = [
    "new_est_reads",
    "fraction_total_reads",
    "kraken_assigned_reads",
    "added_reads",
    "reads",
    "read_count",
    "count",
    "abundance",
]

DEFAULT_EXCLUDE_REGEX = (
    r"(^$|baseMean|log2FoldChange|lfcSE|(^|[._-])stat([._-]|$)|"
    r"pvalue|padj|qvalue|FDR|PValue|P.Value|"
    r"significant|regulation|comparison|contrast|"
    r"taxonomy_lvl|taxonomy_level|taxonomic_level|"
    r"kraken_assigned_reads|added_reads|new_est_reads|fraction_total_reads)"
)


def die(msg: str, exit_code: int = 1) -> None:
    print(f"[ERROR] {msg}", file=sys.stderr)
    sys.exit(exit_code)


def warn(msg: str) -> None:
    print(f"[WARNING] {msg}", file=sys.stderr)


def info(msg: str) -> None:
    print(f"[INFO] {msg}", file=sys.stderr)


def detect_sep(path: str, user_sep: str = "auto") -> str:
    if user_sep != "auto":
        return "\t" if user_sep == r"\t" else user_sep

    with open(path, "r", encoding="utf-8", errors="replace", newline="") as handle:
        sample = handle.read(8192)

    lines = sample.splitlines()
    first_line = lines[0] if lines else ""
    if not first_line:
        die(f"Empty file: {path}")

    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=",\t;")
        return dialect.delimiter
    except Exception:
        counts = {sep: first_line.count(sep) for sep in [",", "\t", ";"]}
        return max(counts, key=counts.get)


def dedupe_headers(headers: Sequence[str]) -> List[str]:
    counts: Counter[str] = Counter()
    out: List[str] = []
    for h in headers:
        base = str(h).strip()
        if counts[base] == 0:
            out.append(base)
        else:
            out.append(f"{base}.{counts[base]}")
        counts[base] += 1
    return out


def read_any_table(path: str, sep: str = "auto") -> Tuple[List[str], List[Dict[str, str]]]:
    real_sep = detect_sep(path, sep)
    try:
        with open(path, "r", encoding="utf-8", errors="replace", newline="") as handle:
            reader = csv.reader(handle, delimiter=real_sep)
            try:
                raw_header = next(reader)
            except StopIteration:
                die(f"Empty file: {path}")

            header = dedupe_headers(raw_header)
            rows: List[Dict[str, str]] = []
            for line_number, fields in enumerate(reader, start=2):
                if not fields or all(str(x).strip() == "" for x in fields):
                    continue
                # Pad short lines or trim long lines so the script does not crash.
                if len(fields) < len(header):
                    fields = fields + [""] * (len(header) - len(fields))
                elif len(fields) > len(header):
                    fields = fields[: len(header)]
                rows.append({header[i]: fields[i] for i in range(len(header))})
    except Exception as e:
        die(f"Could not read {path}: {e}")

    return header, rows


def parse_float(value) -> Optional[float]:
    if value is None:
        return None
    s = str(value).strip().strip('"').strip("'")
    if s == "" or s.lower() in {"na", "nan", "none", "null", "inf", "-inf"}:
        return None

    # Basic support for comma decimals: 1,23 -> 1.23.
    # Values that look like thousands separators, such as 1,234.56, are not changed.
    if re.match(r"^-?\d+,\d+(e[+-]?\d+)?$", s, flags=re.IGNORECASE):
        s = s.replace(",", ".")

    try:
        x = float(s)
        if math.isfinite(x):
            return x
        return None
    except Exception:
        return None


def format_number(x: float) -> str:
    # Compact number format without unnecessary scientific notation weirdness.
    if x is None:
        return "0"
    if abs(x) < 1e-12:
        x = 0.0
    return f"{x:.10g}"


def find_column(columns: Sequence[str], candidates: Sequence[str]) -> Optional[str]:
    lower_map = {c.lower(): c for c in columns}
    for cand in candidates:
        if cand in columns:
            return cand
        if cand.lower() in lower_map:
            return lower_map[cand.lower()]
    return None


def column_numeric_ratio(rows: Sequence[Dict[str, str]], col: str) -> Tuple[int, int, float]:
    non_empty = 0
    numeric = 0
    for row in rows:
        val = row.get(col, "")
        if str(val).strip() == "":
            continue
        non_empty += 1
        if parse_float(val) is not None:
            numeric += 1
    ratio = numeric / non_empty if non_empty else 0.0
    return non_empty, numeric, ratio


def find_first_non_numeric_col(columns: Sequence[str], rows: Sequence[Dict[str, str]]) -> Optional[str]:
    for col in columns:
        non_empty, numeric, ratio = column_numeric_ratio(rows, col)
        if non_empty > 0 and ratio < 0.5:
            return col
    return None


def choose_name_id_cols(
    columns: Sequence[str],
    rows: Sequence[Dict[str, str]],
    id_col: Optional[str] = None,
    name_col: Optional[str] = None,
) -> Tuple[Optional[str], Optional[str]]:
    if name_col:
        if name_col not in columns:
            die(f"--name-col '{name_col}' does not exist. Available columns: {', '.join(columns)}")
        detected_name = name_col
    else:
        detected_name = find_column(columns, NAME_CANDIDATES) or find_first_non_numeric_col(columns, rows)

    if id_col:
        if id_col not in columns:
            die(f"--id-col '{id_col}' does not exist. Available columns: {', '.join(columns)}")
        detected_id = id_col
    else:
        detected_id = find_column(columns, ID_CANDIDATES)

    return detected_id, detected_name


def make_unique(names: Sequence[str]) -> List[str]:
    counts: Counter[str] = Counter()
    out: List[str] = []
    for n in names:
        base = str(n)
        if counts[base] == 0:
            out.append(base)
        else:
            out.append(f"{base}.{counts[base]}")
        counts[base] += 1
    return out


def clean_sample_name(name: str, strip_suffix: bool = True) -> str:
    n = str(name)
    if strip_suffix:
        n = re.sub(r"\.normtrans$", "", n)
        n = re.sub(r"\.norm$", "", n)
    return n


def condition_name_from_sample(name: str) -> str:
    n = str(name)
    upper = n.upper()
    if "LIVER" in upper:
        return "Liver"
    if re.search(r"(^|[_\-.])TEL([_\-.]|$)", upper) or "TELENCEPHALON" in upper:
        return "Telencephalon"
    return clean_sample_name(n, strip_suffix=True)


def detect_value_col(columns: Sequence[str], value_col: str = "auto") -> str:
    lower_map = {c.lower(): c for c in columns}
    if value_col != "auto":
        if value_col in columns:
            return value_col
        if value_col.lower() in lower_map:
            return lower_map[value_col.lower()]
        die(f"--value-col '{value_col}' does not exist. Columns: {', '.join(columns)}")

    for cand in BRACKEN_VALUE_CANDIDATES:
        if cand.lower() in lower_map:
            return lower_map[cand.lower()]

    die(
        "Could not find a value column in the Bracken report. "
        f"Try --value-col. Candidate columns: {', '.join(BRACKEN_VALUE_CANDIDATES)}"
    )


def is_bracken_single_report(
    columns: Sequence[str],
    rows: Sequence[Dict[str, str]],
    name_col: Optional[str],
    value_col: str,
) -> bool:
    lower_cols = {c.lower() for c in columns}
    has_name = name_col is not None or any(c.lower() in {x.lower() for x in NAME_CANDIDATES} for c in columns)
    if value_col != "auto":
        has_value = value_col in columns or value_col.lower() in lower_cols
    else:
        has_value = any(v.lower() in lower_cols for v in BRACKEN_VALUE_CANDIDATES)

    numeric_cols = 0
    for col in columns:
        non_empty, numeric, ratio = column_numeric_ratio(rows, col)
        if numeric > 0:
            numeric_cols += 1
    return bool(has_name and has_value and numeric_cols <= 8)


def sample_name_from_path(path: str) -> str:
    base = os.path.basename(path)
    base = re.sub(r"\.(bracken|tsv|txt|csv|report)$", "", base, flags=re.IGNORECASE)
    base = re.sub(r"([._-]?bracken.*)$", "", base, flags=re.IGNORECASE)
    base = re.sub(r"([._-]?kraken.*)$", "", base, flags=re.IGNORECASE)
    return base or os.path.basename(path)


def numeric_candidate_columns(
    columns: Sequence[str],
    rows: Sequence[Dict[str, str]],
    ignore_cols: Sequence[Optional[str]],
    exclude_regex: str,
    sample_regex: Optional[str],
) -> List[str]:
    ignore = {c for c in ignore_cols if c is not None}
    exclude_re = re.compile(exclude_regex, re.IGNORECASE) if exclude_regex else None
    sample_re = re.compile(sample_regex) if sample_regex else None
    candidates: List[str] = []

    for col in columns:
        if col in ignore:
            continue
        if exclude_re and exclude_re.search(col):
            continue
        if sample_re and not sample_re.search(col):
            continue
        non_empty, numeric, ratio = column_numeric_ratio(rows, col)
        if non_empty > 0 and ratio >= 0.8:
            candidates.append(col)

    return candidates


def prefer_columns(cols: Sequence[str], prefer: str) -> List[str]:
    cols = list(cols)
    normtrans = [c for c in cols if re.search(r"\.normtrans$", c)]
    norm = [c for c in cols if re.search(r"\.norm$", c)]
    raw = [c for c in cols if not re.search(r"\.norm(trans)?$", c)]

    if prefer == "normtrans":
        return normtrans
    if prefer == "norm":
        return norm
    if prefer == "raw":
        return raw

    if normtrans:
        return normtrans
    if norm:
        return norm
    return raw


def build_from_wide_matrix(
    path: str,
    sep: str,
    id_col: Optional[str],
    name_col: Optional[str],
    prefer: str,
    sample_regex: Optional[str],
    exclude_regex: str,
    fill_na: float,
    strip_suffix: bool,
    sample_name_mode: str,
) -> Tuple[List[str], List[List[object]]]:
    columns, rows = read_any_table(path, sep)
    detected_id, detected_name = choose_name_id_cols(columns, rows, id_col=id_col, name_col=name_col)

    if detected_name is None and detected_id is None:
        warn("No ID/name column was found. UNIQID/NAME will be created from row numbers.")

    ignore_cols = [detected_id, detected_name]
    candidates = numeric_candidate_columns(
        columns,
        rows,
        ignore_cols=ignore_cols,
        exclude_regex=exclude_regex,
        sample_regex=sample_regex,
    )
    chosen = prefer_columns(candidates, prefer=prefer)

    if not chosen:
        die(
            "Could not find numeric sample columns. "
            "Try --sample-regex, --prefer raw/norm/normtrans, or adjust --exclude-regex."
        )

    if sample_name_mode == "condition":
        sample_cols = make_unique([condition_name_from_sample(c) for c in chosen])
    else:
        sample_cols = make_unique([clean_sample_name(c, strip_suffix=strip_suffix) for c in chosen])

    header = ["UNIQID", "NAME"] + sample_cols
    out_rows: List[List[object]] = []

    for i, row in enumerate(rows, start=1):
        if detected_id is not None:
            uniqid = str(row.get(detected_id, "")).strip()
        elif detected_name is not None:
            uniqid = str(row.get(detected_name, "")).strip()
        else:
            uniqid = f"feature_{i}"

        if detected_name is not None:
            name = str(row.get(detected_name, "")).strip()
        elif detected_id is not None:
            name = str(row.get(detected_id, "")).strip()
        else:
            name = uniqid

        values = []
        for col in chosen:
            x = parse_float(row.get(col, ""))
            values.append(fill_na if x is None else x)

        if str(uniqid).strip() and str(name).strip():
            out_rows.append([uniqid, name] + values)

    return header, out_rows


def build_from_single_bracken_reports(
    paths: Sequence[str],
    sep: str,
    id_col: Optional[str],
    name_col: Optional[str],
    value_col: str,
    fill_na: float,
) -> Tuple[List[str], List[List[object]]]:
    sample_names: List[str] = []
    sample_values: Dict[Tuple[str, str], Dict[str, float]] = OrderedDict()

    for path in paths:
        columns, rows = read_any_table(path, sep)
        detected_id, detected_name = choose_name_id_cols(columns, rows, id_col=id_col, name_col=name_col)
        if detected_name is None:
            die(f"Could not find a name/taxon column in {path}. Use --name-col.")

        val_col = detect_value_col(columns, value_col)
        sample = sample_name_from_path(path)
        sample_names.append(sample)

        for row in rows:
            name = str(row.get(detected_name, "")).strip()
            if not name:
                continue
            if detected_id is not None and detected_id != detected_name:
                uniqid = str(row.get(detected_id, "")).strip() or name
            else:
                uniqid = name
            key = (uniqid, name)
            if key not in sample_values:
                sample_values[key] = {}
            x = parse_float(row.get(val_col, ""))
            sample_values[key][sample] = fill_na if x is None else x

    sample_names = make_unique(sample_names)
    header = ["UNIQID", "NAME"] + sample_names
    out_rows: List[List[object]] = []

    for key, values_by_sample in sample_values.items():
        uniqid, name = key
        vals = [values_by_sample.get(s, fill_na) for s in sample_names]
        out_rows.append([uniqid, name] + vals)

    return header, out_rows


def transform_matrix(header: List[str], rows: List[List[object]], transform: str) -> List[List[object]]:
    if transform == "none":
        return rows

    out: List[List[object]] = []
    for row in rows:
        prefix = row[:2]
        vals = [parse_float(v) for v in row[2:]]
        vals = [0.0 if v is None else v for v in vals]

        if transform == "log2p1":
            new_vals = [math.log2(max(0.0, v) + 1.0) for v in vals]
        elif transform == "row_zscore":
            n = len(vals)
            mean = sum(vals) / n if n else 0.0
            if n > 1:
                var = sum((v - mean) ** 2 for v in vals) / (n - 1)
                sd = math.sqrt(var)
            else:
                sd = 0.0
            if sd == 0.0:
                new_vals = [0.0 for _ in vals]
            else:
                new_vals = [(v - mean) / sd for v in vals]
        else:
            die(f"Unknown transformation: {transform}")

        out.append(prefix + new_vals)

    return out


def filter_and_sort(
    rows: List[List[object]],
    min_sum: Optional[float],
    top: Optional[int],
) -> List[List[object]]:
    out = rows

    if min_sum is not None:
        tmp = []
        for row in out:
            vals = [parse_float(v) or 0.0 for v in row[2:]]
            if sum(vals) >= min_sum:
                tmp.append(row)
        out = tmp

    if top is not None and top > 0 and len(out) > top:
        def score(row: List[object]) -> float:
            vals = [parse_float(v) or 0.0 for v in row[2:]]
            if not vals:
                return 0.0
            return sum(abs(v) for v in vals) / len(vals)

        out = sorted(out, key=score, reverse=True)[:top]

    return out


def write_table(path: str, header: List[str], rows: List[List[object]], out_sep: str) -> None:
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter=out_sep, lineterminator="\n")
        writer.writerow(header)
        for row in rows:
            clean = [str(row[0]), str(row[1])] + [format_number(parse_float(v) or 0.0) for v in row[2:]]
            writer.writerow(clean)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Convert Bracken/MTD/DEG tables into a heatmap matrix with this format: "
            "UNIQID, NAME, numeric sample columns. This version does not require pandas."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("-i", "--input", nargs="+", required=True,
                        help="Input file(s). Can be one wide matrix or multiple single-sample Bracken reports.")
    parser.add_argument("-o", "--output", default="heatmap_table.txt", help="Output file.")
    parser.add_argument("--sep", default="auto", help="Input separator: auto, ',', ';', or '\\t'.")
    parser.add_argument("--out-sep", default=r"\t", help="Output separator. Use '\\t' for tab.")
    parser.add_argument("--id-col", default=None, help="ID column name. If missing, NAME is copied into UNIQID.")
    parser.add_argument("--name-col", default=None, help="Name/taxon column name. If missing, UNIQID is copied into NAME.")
    parser.add_argument("--value-col", default="auto", help="For Bracken reports: column used as the abundance/value column.")
    parser.add_argument("--prefer", choices=["auto", "raw", "norm", "normtrans"], default="auto",
                        help="For DEG/wide tables, choose raw, .norm, or .normtrans columns. Auto prefers .normtrans.")
    parser.add_argument("--sample-regex", default=None, help="Regex to keep only selected sample columns. Example: 'LIVER|TEL'.")
    parser.add_argument("--exclude-regex", default=DEFAULT_EXCLUDE_REGEX, help="Regex for columns that should be excluded automatically.")
    parser.add_argument("--fill-na", type=float, default=0.0, help="Value used to replace NA or non-numeric values.")
    parser.add_argument("--transform", choices=["none", "log2p1", "row_zscore"], default="none",
                        help="Optional transformation applied to sample columns.")
    parser.add_argument("--min-sum", type=float, default=None, help="Remove features whose sample-value sum is lower than this value.")
    parser.add_argument("--top", type=int, default=None, help="Keep only the top N features by mean absolute value.")
    parser.add_argument("--sample-name-mode", choices=["keep", "condition"], default="keep",
                        help="keep preserves sample names; condition renames LIVER to Liver and TEL to Telencephalon.")
    parser.add_argument("--keep-suffix", action="store_true", help="Keep .norm and .normtrans suffixes in sample column names.")

    args = parser.parse_args()
    out_sep = "\t" if args.out_sep == r"\t" else args.out_sep

    missing = [p for p in args.input if not os.path.exists(p)]
    if missing:
        die("File(s) not found: " + ", ".join(missing))

    if len(args.input) > 1:
        first_cols, first_rows = read_any_table(args.input[0], args.sep)
        _, first_name = choose_name_id_cols(first_cols, first_rows, id_col=args.id_col, name_col=args.name_col)
        use_bracken_merge = is_bracken_single_report(first_cols, first_rows, first_name, args.value_col)
    else:
        use_bracken_merge = False

    if use_bracken_merge:
        info("Detected mode: multiple single-sample Bracken reports; merging by UNIQID/NAME.")
        header, rows = build_from_single_bracken_reports(
            paths=args.input,
            sep=args.sep,
            id_col=args.id_col,
            name_col=args.name_col,
            value_col=args.value_col,
            fill_na=args.fill_na,
        )
    else:
        if len(args.input) > 1:
            warn("More than one input was provided, but they do not look like single-sample Bracken reports. Only the first file will be used.")
        info("Detected mode: wide/DEG matrix; selecting numeric sample columns.")
        header, rows = build_from_wide_matrix(
            path=args.input[0],
            sep=args.sep,
            id_col=args.id_col,
            name_col=args.name_col,
            prefer=args.prefer,
            sample_regex=args.sample_regex,
            exclude_regex=args.exclude_regex,
            fill_na=args.fill_na,
            strip_suffix=not args.keep_suffix,
            sample_name_mode=args.sample_name_mode,
        )

    rows = transform_matrix(header, rows, args.transform)
    rows = filter_and_sort(rows, min_sum=args.min_sum, top=args.top)

    write_table(args.output, header, rows, out_sep)

    sample_cols = header[2:]
    info(f"Saved file: {args.output}")
    info(f"Features/rows: {len(rows)}")
    info(f"Samples/numeric columns: {len(sample_cols)}")
    info("First sample columns: " + ", ".join(sample_cols[:8]) + (" ..." if len(sample_cols) > 8 else ""))


if __name__ == "__main__":
    main()

