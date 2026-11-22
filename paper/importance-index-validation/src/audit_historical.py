"""Independent reproduction of the two frozen taxon-importance rankings.

Uses standard-library arithmetic and TSV parsing, without importing the main
robustness analysis or executing the historical R script. Source files are read
only. No bootstrap, subsampling, alternate score or pipeline run is performed.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import statistics
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "provenance/independent"
WEIGHTS = (0.45, 0.30, 0.15, 0.10)
R_CODE = ROOT / "provenance/historical_importance_index.R"


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t", quoting=csv.QUOTE_NONE)
        return list(reader.fieldnames or []), list(reader)


def identity(path: Path) -> dict:
    return {
        "path": str(path.relative_to(ROOT)),
        "bytes": path.stat().st_size,
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
    }


def minmax_log(values: list[float]) -> list[float]:
    logged = [math.log10(value + 1) for value in values]
    low, high = min(logged), max(logged)
    if not math.isfinite(low) or not math.isfinite(high) or high == low:
        return [0.0] * len(values)
    return [(value - low) / (high - low) for value in logged]


def reproduce(mode: str) -> dict:
    representation = "relative_percent" if mode == "relative" else "estimated_reads"
    wide_path = ROOT / "data" / f"{representation}_matrix.tsv"
    ranked_path = ROOT / "data" / f"{representation}_historical_ranking.tsv"
    summary_path = ROOT / "data" / f"{representation}_historical_summary.tsv"
    headers, wide = read_tsv(wide_path)
    _, historical = read_tsv(ranked_path)
    _, summary = read_tsv(summary_path)
    names = [row["taxon"] for row in wide]
    samples = headers[1:]
    if len(names) != len(set(names)):
        raise ValueError("Duplicate taxon labels in exported sample matrix")
    hist_by_name = {row["taxon"]: row for row in historical}
    summary_by_name = {row["taxon"]: row for row in summary}
    if set(names) != set(hist_by_name) or set(names) != set(summary_by_name):
        raise ValueError("Different taxon universes in matrix and frozen ranking/summary")
    values = [[float(row[sample]) for sample in samples] for row in wide]
    if not all(math.isfinite(value) and value >= 0 for row in values for value in row):
        raise ValueError("Exported source matrix contains negative/nonfinite values")
    if any(not any(value > 0 for value in row) for row in values):
        raise ValueError("Exported matrix unexpectedly contains all-zero taxa")

    n = len(samples)
    calculated = []
    for taxon, row in zip(names, values):
        total = math.fsum(row)
        presence = sum(value > 0 for value in row)
        calculated.append({
            "taxon": taxon,
            "n_samples": n,
            "prevalence_n": presence,
            "prevalence_percent": 100 * presence / n,
            "mean_abundance": total / n,
            "median_abundance": statistics.median(row),
            "max_abundance": max(row),
            "total_abundance": total,
            "prevalence_score": presence / n,
        })
    for metric in ("mean", "max", "total"):
        component = minmax_log([row[f"{metric}_abundance"] for row in calculated])
        for row, value in zip(calculated, component):
            row[f"{metric}_abundance_score"] = value
    for row in calculated:
        full_components = (
            row["prevalence_score"], row["mean_abundance_score"],
            row["max_abundance_score"], row["total_abundance_score"],
        )
        score = sum(weight * value for weight, value in zip(WEIGHTS, full_components))
        row["importance_score_unrounded"] = score
        row["importance_score"] = round(score, 6)
        for key in ("prevalence_score", "mean_abundance_score", "max_abundance_score", "total_abundance_score"):
            row[key] = round(row[key], 6)
    ordered = sorted(calculated, key=lambda row: (
        -row["importance_score"], -row["prevalence_percent"],
        -row["mean_abundance"], -row["max_abundance"], row["taxon"],
    ))
    for position, row in enumerate(ordered, 1):
        row["rank_position"] = position

    numeric_fields = (
        "n_samples", "prevalence_n", "prevalence_percent", "mean_abundance",
        "median_abundance", "max_abundance", "total_abundance", "importance_score",
        "prevalence_score", "mean_abundance_score", "max_abundance_score", "total_abundance_score",
    )
    errors = {}
    mismatches = {}
    for key in numeric_fields:
        differences = [abs(row[key] - float(hist_by_name[row["taxon"]][key])) for row in calculated]
        errors[key] = max(differences)
        tolerance = 1e-12 if mode == "relative" else 1e-9
        if key.endswith("score") or key in ("n_samples", "prevalence_n"):
            tolerance = 0
        mismatches[key] = sum(difference > tolerance for difference in differences)
    order_mismatches = [
        {"recomputed_rank": i + 1, "recomputed_taxon": new["taxon"], "historical_taxon": old["taxon"]}
        for i, (new, old) in enumerate(zip(ordered, historical)) if new["taxon"] != old["taxon"]
    ]
    summary_errors = {
        key: max(abs(row[key] - float(summary_by_name[row["taxon"]][key])) for row in calculated)
        for key in numeric_fields if key in summary[0]
    }
    score_groups = {}
    complete_ties = {}
    for row in ordered:
        score_groups.setdefault(row["importance_score"], []).append(row["taxon"])
        key = (row["importance_score"], row["prevalence_percent"], row["mean_abundance"], row["max_abundance"])
        complete_ties.setdefault(key, []).append(row["taxon"])
    ties = [group for group in score_groups.values() if len(group) > 1]
    tie_break_by_label = [group for group in complete_ties.values() if len(group) > 1]
    column_totals = [math.fsum(row[j] for row in values) for j in range(n)]
    output_path = OUT / f"independently_recomputed_{mode}_ranking.tsv"
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        fields = ["rank_position", "taxon", *numeric_fields, "importance_score_unrounded"]
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(ordered)
    return {
        "representation": mode,
        "units": "within-library Bracken percentages" if mode == "relative" else "Bracken estimated read counts",
        "taxon_identity_key": "taxon label; numeric TaxIDs are not supplied in these exported matrices",
        "n_taxa": len(names),
        "n_samples": n,
        "samples": samples,
        "source_files": [identity(path) for path in (wide_path, ranked_path, summary_path)],
        "matrix_column_totals": dict(zip(samples, column_totals)),
        "all_values_finite_nonnegative": True,
        "duplicate_taxa": 0,
        "zero_total_taxa": 0,
        "max_abs_error_against_frozen_ranking": errors,
        "mismatch_counts_at_tolerance": mismatches,
        "max_abs_error_against_original_summary": summary_errors,
        "complete_rank_order_exact": not order_mismatches,
        "rank_order_mismatches": order_mismatches,
        "rounded_score_tie_groups": ties,
        "all_numeric_ranking_keys_tied_groups": tie_break_by_label,
        "dependency_total_equals_n_times_mean_max_abs_error": max(abs(row["total_abundance"] - n * row["mean_abundance"]) for row in calculated),
        "output_file": str(output_path.relative_to(ROOT)),
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    reports = [reproduce(mode) for mode in ("relative", "absolute")]
    modes = {}
    for mode in ("relative", "absolute"):
        representation = "relative_percent" if mode == "relative" else "estimated_reads"
        _, rows = read_tsv(ROOT / "data" / f"{representation}_matrix.tsv")
        modes[mode] = {row["taxon"]: row for row in rows}
    samples = reports[0]["samples"]
    relative_from_counts_errors = []
    for sample in samples:
        denom = math.fsum(float(row[sample]) for row in modes["absolute"].values())
        for taxon, row in modes["absolute"].items():
            relative_from_counts_errors.append(abs(100 * float(row[sample]) / denom - float(modes["relative"][taxon][sample])))
    report = {
        "task": "Independent exact reproduction of historical taxon importance rankings",
        "date": "2026-11-20",
        "historical_R_code": identity(R_CODE),
        "method": "Independent standard-library Python arithmetic using exported full per-library matrices; no root analysis imported; no statistical resampling performed",
        "weights_prevalence_mean_max_total": WEIGHTS,
        "presence_rule": "abundance > 0",
        "scaling": "For each abundance summary across nonzero-total taxa, log10(v+1), then min-max. Constant/nonfinite extrema produce zero component.",
        "rounding_and_sort": "Score calculated with unrounded components, then score/components rounded to6decimals; sorting by score, prevalence_percent, raw mean, raw max descending, then taxon ascending. Python lexicographic order reproduced the entire frozen ordering here; R character collation can be locale-specific in a different dataset.",
        "source_preprocessing_from_R": "Reader converts missing and negative abundance to0; aggregates duplicate taxon names before removing all-zero rows; relative mode divides each library by sum over remaining taxa and multiplies100. Full ranked table scores precomputed summary input; wide per-library export is a separate subsequent step and is independently checked here for consistency.",
        "representations": reports,
        "relative_matrix_from_absolute_max_abs_error_percentage_points": max(relative_from_counts_errors),
        "source_files_modified": False,
        "bootstrap_or_sensitivity_analysis_run": False,
        "interpretation": [
            "Total equals n times arithmetic mean within a fixed sample universe; they are not independent underlying abundance summaries.",
            "The log1p and separate min-max transforms mean their scaled score vectors are not numerically identical; removing total is still a meaningful sensitivity scenario, not a change to the frozen score definition.",
            "Maximum is sensitive to high abundance in a single sampled library; bootstrap copies of a library do not create independent new biological samples.",
            "For resampling, recompute prevalence, mean, maximum, total, min-max and rounded tie ordering; taxa absent from all selected libraries need an explicit ineligibility/rank policy.",
            "Percentages have denominator all425taxa per library before any display subset; do not renormalize only the top20.",
        ],
    }
    (OUT / "historical_reproduction_report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps({
        "n_taxa_per_mode": [r["n_taxa"] for r in reports],
        "rank_order_exact_per_mode": [r["complete_rank_order_exact"] for r in reports],
        "score_max_abs_error_per_mode": [r["max_abs_error_against_frozen_ranking"]["importance_score"] for r in reports],
        "numeric_mismatches_per_mode": [r["mismatch_counts_at_tolerance"] for r in reports],
        "relative_from_counts_max_abs_error": max(relative_from_counts_errors),
    }, indent=2))
    if any(not r["complete_rank_order_exact"] or any(r["mismatch_counts_at_tolerance"].values()) for r in reports):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
