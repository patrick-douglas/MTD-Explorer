#!/usr/bin/env python3

import argparse
import csv
import os
import re
from pathlib import Path
from collections import defaultdict

def parse_args():
    p = argparse.ArgumentParser(
        description="Search MTD outputs for target parasite/taxa signals such as Schistosoma."
    )

    p.add_argument(
        "--mtd_res",
        required=True,
        help="MTD result directory, e.g. /home/me/projeto_caramujo/MTD_res"
    )

    p.add_argument(
        "--terms",
        default="Schistosoma,Schistosomatidae,Trematoda,Digenea,Platyhelminthes",
        help="Comma-separated search terms. Default: Schistosoma,Schistosomatidae,Trematoda,Digenea,Platyhelminthes"
    )

    p.add_argument(
        "--outdir",
        default=None,
        help="Output directory. Default: <mtd_res>/figures/Figure4_target_detection"
    )

    p.add_argument(
        "--extra_search_dir",
        default=None,
        help="Optional extra directory to search, e.g. Trinity/BLAST COI results"
    )

    p.add_argument(
        "--max_line_length",
        type=int,
        default=800,
        help="Maximum characters saved per raw matching line"
    )

    return p.parse_args()

def safe_open(path):
    try:
        return open(path, "r", encoding="utf-8", errors="replace")
    except Exception:
        return None

def detect_sep(path):
    with safe_open(path) as f:
        if f is None:
            return "\t"
        first = f.readline()
    if first.count("\t") >= first.count(","):
        return "\t"
    return ","

def norm_name(x):
    x = str(x)
    x = x.strip()
    x = x.replace("_", " ")
    x = re.sub(r"\s+", " ", x)
    return x

def lower_terms(terms):
    return [t.strip().lower() for t in terms if t.strip()]

def find_term(text, terms_lc):
    tx = str(text).lower()
    for term in terms_lc:
        if term in tx:
            return term
    return ""

def detect_taxon_col(header):
    lower = [h.strip().lower() for h in header]
    candidates = [
        "name", "taxon", "taxonomy", "classification", "#classification",
        "taxonomy_name", "scientific_name"
    ]
    for c in candidates:
        if c in lower:
            return lower.index(c)
    return 0

def clean_sample_from_col(col):
    x = os.path.basename(str(col))

    patterns = [
        r"\.phylum\.bracken(_num|_frac)?$",
        r"\.genus\.bracken(_num|_frac)?$",
        r"\.species\.bracken(_num|_frac)?$",
        r"_phylum\.bracken(_num|_frac)?$",
        r"_genus\.bracken(_num|_frac)?$",
        r"_species\.bracken(_num|_frac)?$",
        r"\.bracken(_num|_frac)?$",
        r"_bracken(_num|_frac)?$",
        r"-bracken(_num|_frac)?$",
        r"(_num|_frac)$",
        r"(\.num|\.frac)$",
    ]

    x = re.sub(r"^Report_", "", x)

    for pat in patterns:
        x = re.sub(pat, "", x)

    return x

def to_float(x):
    try:
        x = str(x).strip().replace(",", "").replace("%", "")
        if x == "":
            return 0.0
        return float(x)
    except Exception:
        return 0.0

def file_rank_from_name(path):
    name = path.name.lower()
    if "species" in name:
        return "species"
    if "genus" in name:
        return "genus"
    if "phylum" in name:
        return "phylum"
    if "family" in name:
        return "family"
    if "order" in name:
        return "order"
    if "class" in name:
        return "class"
    return "unknown"

def is_probably_text(path):
    suffix = path.suffix.lower()
    if suffix in [".png", ".pdf", ".jpg", ".jpeg", ".svg", ".html", ".gz", ".bz2", ".zip", ".bam", ".sam"]:
        return False
    return True

def is_blast_like(path):
    name = path.name.lower()
    suffix = path.suffix.lower()
    return (
        "blast" in name or
        "coi" in name or
        "cox" in name or
        suffix in [".out", ".tab", ".tsv", ".csv"]
    )

def is_taxonomy_like(path):
    name = path.name.lower()
    suffix = path.suffix.lower()

    keywords = [
        "bracken",
        "kraken",
        "krona",
        "mpa",
        "combined",
        "report",
        "taxonomy",
        "graphlan"
    ]

    return any(k in name for k in keywords) or suffix in [".mpa", ".tsv", ".csv", ".txt"]

def parse_combined_bracken(path, terms_lc):
    hits = []

    sep = detect_sep(path)

    with safe_open(path) as f:
        if f is None:
            return hits

        reader = csv.reader(f, delimiter=sep)
        rows = [r for r in reader if r and any(str(c).strip() for c in r)]

    if not rows or len(rows[0]) < 2:
        return hits

    header = rows[0]
    taxon_idx = detect_taxon_col(header)

    value_cols = []
    for idx, col in enumerate(header):
        if idx == taxon_idx:
            continue

        lc = col.lower()
        if lc.endswith("_num") or lc.endswith(".num"):
            value_cols.append((idx, clean_sample_from_col(col), "count"))
        elif lc.endswith("_frac") or lc.endswith(".frac"):
            value_cols.append((idx, clean_sample_from_col(col), "fraction"))

    # fallback: numeric columns
    if not value_cols:
        for idx, col in enumerate(header):
            if idx == taxon_idx:
                continue
            value_cols.append((idx, clean_sample_from_col(col), "value"))

    rank = file_rank_from_name(path)

    for r in rows[1:]:
        if taxon_idx >= len(r):
            continue

        taxon = norm_name(r[taxon_idx])
        matched = find_term(taxon, terms_lc)

        if not matched:
            continue

        for idx, sample, vtype in value_cols:
            val = to_float(r[idx]) if idx < len(r) else 0.0

            hits.append({
                "file": str(path),
                "rank": rank,
                "matched_term": matched,
                "taxon": taxon,
                "sample": sample,
                "value_type": vtype,
                "value": val
            })

    return hits

def parse_mpa(path, terms_lc):
    hits = []

    sep = detect_sep(path)

    with safe_open(path) as f:
        if f is None:
            return hits

        reader = csv.reader(f, delimiter=sep)
        rows = [r for r in reader if r and any(str(c).strip() for c in r)]

    if not rows or len(rows[0]) < 2:
        return hits

    header = rows[0]
    first = header[0].strip().lower()

    if first.startswith("#") or "classification" in first or "clade" in first:
        sample_names = [h.strip() for h in header[1:]]
        data_rows = rows[1:]
    else:
        sample_names = ["sample_%d" % i for i in range(1, len(header))]
        data_rows = rows

    for r in data_rows:
        if len(r) < 2:
            continue

        taxonomy = r[0]
        taxonomy_clean = taxonomy.replace("|", "; ")
        taxonomy_clean = taxonomy_clean.replace("__", "__")
        taxonomy_clean = norm_name(taxonomy_clean)

        matched = find_term(taxonomy_clean, terms_lc)

        if not matched:
            continue

        for i, sample in enumerate(sample_names, start=1):
            val = to_float(r[i]) if i < len(r) else 0.0
            hits.append({
                "file": str(path),
                "matched_term": matched,
                "taxonomy": taxonomy,
                "sample": sample,
                "value": val
            })

    return hits

def scan_raw_lines(paths, terms_lc, max_len):
    raw_hits = []

    for path in paths:
        if not path.exists() or not path.is_file():
            continue

        if not is_probably_text(path):
            continue

        with safe_open(path) as f:
            if f is None:
                continue

            for n, line in enumerate(f, start=1):
                matched = find_term(line, terms_lc)

                if matched:
                    raw_hits.append({
                        "file": str(path),
                        "line_no": n,
                        "matched_term": matched,
                        "line": line.strip()[:max_len]
                    })

    return raw_hits

def write_tsv(path, rows, fieldnames):
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        for r in rows:
            writer.writerow({k: r.get(k, "") for k in fieldnames})

def main():
    args = parse_args()

    mtd_res = Path(args.mtd_res).expanduser().resolve()
    terms = [t.strip() for t in args.terms.split(",") if t.strip()]
    terms_lc = lower_terms(terms)

    if args.outdir:
        outdir = Path(args.outdir).expanduser().resolve()
    else:
        outdir = mtd_res / "figures" / "Figure4_target_detection"

    outdir.mkdir(parents=True, exist_ok=True)

    search_roots = [mtd_res]

    if args.extra_search_dir:
        extra = Path(args.extra_search_dir).expanduser().resolve()
        if extra.exists():
            search_roots.append(extra)

    print("============================================================")
    print("MTD target taxon search")
    print("MTD result dir:", mtd_res)
    print("Search terms:", ", ".join(terms))
    print("Output dir:", outdir)
    print("============================================================")

    all_files = []
    for root in search_roots:
        for p in root.rglob("*"):
            if p.is_file():
                all_files.append(p)

    # Combined Bracken files
    bracken_combined = [
        p for p in all_files
        if p.name in ["bracken_species_all", "bracken_genus_all", "bracken_phylum_all"] or
           re.match(r"bracken_.*_all$", p.name)
    ]

    # MPA files
    mpa_files = [
        p for p in all_files
        if p.name in ["Combined.mpa", "Combined_for_graphlan.mpa"] or p.suffix.lower() == ".mpa"
    ]

    # Raw taxonomy-like files
    taxonomy_like_files = [
        p for p in all_files
        if is_taxonomy_like(p)
    ]

    # BLAST/COI-like files
    blast_like_files = [
        p for p in all_files
        if is_blast_like(p)
    ]

    print("[INFO] Combined Bracken files:", len(bracken_combined))
    for p in bracken_combined:
        print("  -", p)

    print("[INFO] MPA files:", len(mpa_files))
    for p in mpa_files:
        print("  -", p)

    bracken_hits = []
    for p in bracken_combined:
        bracken_hits.extend(parse_combined_bracken(p, terms_lc))

    mpa_hits = []
    for p in mpa_files:
        mpa_hits.extend(parse_mpa(p, terms_lc))

    raw_hits = scan_raw_lines(taxonomy_like_files, terms_lc, args.max_line_length)
    blast_hits = scan_raw_lines(blast_like_files, terms_lc, args.max_line_length)

    write_tsv(
        outdir / "target_bracken_hits.tsv",
        bracken_hits,
        ["file", "rank", "matched_term", "taxon", "sample", "value_type", "value"]
    )

    write_tsv(
        outdir / "target_mpa_hits.tsv",
        mpa_hits,
        ["file", "matched_term", "taxonomy", "sample", "value"]
    )

    write_tsv(
        outdir / "target_raw_line_hits.tsv",
        raw_hits,
        ["file", "line_no", "matched_term", "line"]
    )

    write_tsv(
        outdir / "target_blast_like_line_hits.tsv",
        blast_hits,
        ["file", "line_no", "matched_term", "line"]
    )

    # Summary by sample
    summary = defaultdict(lambda: {
        "sample": "",
        "bracken_count_sum": 0.0,
        "bracken_fraction_sum": 0.0,
        "bracken_value_sum": 0.0,
        "bracken_hit_taxa": set(),
        "mpa_value_sum": 0.0,
        "mpa_hit_taxa": set(),
        "raw_line_hits": 0,
        "blast_like_line_hits": 0
    })

    for h in bracken_hits:
        sample = h["sample"]
        summary[sample]["sample"] = sample
        v = float(h["value"])

        if h["value_type"] == "count":
            summary[sample]["bracken_count_sum"] += v
        elif h["value_type"] == "fraction":
            summary[sample]["bracken_fraction_sum"] += v
        else:
            summary[sample]["bracken_value_sum"] += v

        summary[sample]["bracken_hit_taxa"].add(h["taxon"])

    for h in mpa_hits:
        sample = h["sample"]
        summary[sample]["sample"] = sample
        summary[sample]["mpa_value_sum"] += float(h["value"])
        summary[sample]["mpa_hit_taxa"].add(h["taxonomy"])

    # raw/blast line hits are not always sample-specific, but try to infer sample from path/line
    all_samples = set(summary.keys())

    for h in raw_hits:
        text = h["file"] + " " + h["line"]
        matched_samples = [s for s in all_samples if s and s in text]
        if matched_samples:
            for s in matched_samples:
                summary[s]["raw_line_hits"] += 1
        else:
            summary["UNKNOWN"]["sample"] = "UNKNOWN"
            summary["UNKNOWN"]["raw_line_hits"] += 1

    for h in blast_hits:
        text = h["file"] + " " + h["line"]
        matched_samples = [s for s in all_samples if s and s in text]
        if matched_samples:
            for s in matched_samples:
                summary[s]["blast_like_line_hits"] += 1
        else:
            summary["UNKNOWN"]["sample"] = "UNKNOWN"
            summary["UNKNOWN"]["blast_like_line_hits"] += 1

    summary_rows = []

    for sample in sorted(summary.keys()):
        s = summary[sample]
        detected = (
            s["bracken_count_sum"] > 0 or
            s["bracken_fraction_sum"] > 0 or
            s["bracken_value_sum"] > 0 or
            s["mpa_value_sum"] > 0 or
            s["raw_line_hits"] > 0 or
            s["blast_like_line_hits"] > 0
        )

        summary_rows.append({
            "sample": s["sample"],
            "detected_any": "yes" if detected else "no",
            "bracken_count_sum": "%.10g" % s["bracken_count_sum"],
            "bracken_fraction_sum": "%.10g" % s["bracken_fraction_sum"],
            "bracken_value_sum": "%.10g" % s["bracken_value_sum"],
            "bracken_hit_taxa_n": len(s["bracken_hit_taxa"]),
            "bracken_hit_taxa": "; ".join(sorted(s["bracken_hit_taxa"])),
            "mpa_value_sum": "%.10g" % s["mpa_value_sum"],
            "mpa_hit_taxa_n": len(s["mpa_hit_taxa"]),
            "mpa_hit_taxa": "; ".join(sorted(s["mpa_hit_taxa"])),
            "raw_line_hits": s["raw_line_hits"],
            "blast_like_line_hits": s["blast_like_line_hits"]
        })

    write_tsv(
        outdir / "target_detection_summary_by_sample.tsv",
        summary_rows,
        [
            "sample",
            "detected_any",
            "bracken_count_sum",
            "bracken_fraction_sum",
            "bracken_value_sum",
            "bracken_hit_taxa_n",
            "bracken_hit_taxa",
            "mpa_value_sum",
            "mpa_hit_taxa_n",
            "mpa_hit_taxa",
            "raw_line_hits",
            "blast_like_line_hits"
        ]
    )

    print()
    print("[OK] Search finished.")
    print("Outputs:")
    print(" ", outdir / "target_bracken_hits.tsv")
    print(" ", outdir / "target_mpa_hits.tsv")
    print(" ", outdir / "target_raw_line_hits.tsv")
    print(" ", outdir / "target_blast_like_line_hits.tsv")
    print(" ", outdir / "target_detection_summary_by_sample.tsv")
    print()
    print("Quick summary:")
    print("  Bracken hits:", len(bracken_hits))
    print("  MPA hits:", len(mpa_hits))
    print("  Raw taxonomy line hits:", len(raw_hits))
    print("  BLAST-like line hits:", len(blast_hits))

    if summary_rows:
        print()
        print("Samples with any signal:")
        for r in summary_rows:
            if r["detected_any"] == "yes":
                print("  -", r["sample"])
    else:
        print()
        print("[INFO] No sample-level signal found.")

if __name__ == "__main__":
    main()
