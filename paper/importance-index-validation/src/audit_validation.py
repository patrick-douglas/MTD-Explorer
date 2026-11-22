"""Independent audit of the saved importance-index robustness outputs.

Never imports or executes the main analysis script. Recomputes all fixed
scenarios and all leave-one-out rankings with independent Python arithmetic,
plus selected saved bootstrap draws. Audits every summary from saved arrays.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import itertools
import json
import math
import random
import statistics
from fractions import Fraction
from pathlib import Path

import numpy as np

from audit_historical import minmax_log, read_tsv

PACKAGE = Path(__file__).resolve().parents[1]
OUT = PACKAGE / "provenance/independent"
DEFAULT = (.45, .30, .15, .10)
NO_TOTAL = tuple(w / math.fsum(DEFAULT[:3]) for w in DEFAULT[:3]) + (0.,)
REPS = ("relative_percent", "estimated_reads")
MODELS = {"default": DEFAULT, "no_total": NO_TOTAL}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def table(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def rank_independent(names, matrix, weights, presence=None, eligible=None):
    n = len(matrix[0])
    ids = [i for i, row in enumerate(matrix)
           if math.fsum(row) > 0 and (eligible is None or eligible[i])]
    state = []
    for i in ids:
        row = matrix[i]
        total = math.fsum(row)
        prev = sum(v > 0 for v in row) / n if presence is None else sum(presence[i]) / n
        state.append({"id": i, "taxon": names[i], "total": total,
                      "mean": total / n, "maximum": max(row), "prevalence": prev})
    comps = [[s["prevalence"] for s in state]]
    for metric in ("mean", "maximum", "total"):
        comps.append(minmax_log([s[metric] for s in state]))
    for j, s in enumerate(state):
        s["score"] = round(sum(weights[k] * comps[k][j] for k in range(4)), 6)
    ordered = sorted(state, key=lambda s: (-s["score"], -s["prevalence"],
                                         -s["mean"], -s["maximum"], s["taxon"]))
    ranks = np.full(len(names), len(names) + 1, dtype=np.int16)
    scores = np.full(len(names), np.nan)
    for rank, s in enumerate(ordered, 1):
        ranks[s["id"]] = rank
        scores[s["id"]] = s["score"]
    return ranks, scores


def quantile_independent(values, p):
    ordered = sorted(float(v) for v in values)
    if not ordered:
        return None
    h = (len(ordered) - 1) * p
    lo = math.floor(h)
    hi = math.ceil(h)
    return ordered[lo] * (hi - h) + ordered[hi] * (h - lo) if lo != hi else ordered[lo]


def metrics_independent(reference, ranks):
    n = len(ranks)
    metrics = {"eligible_taxa": int(sum(r <= n for r in ranks))}
    for k in (10, 20, 50):
        a = {i for i, r in enumerate(reference) if r <= k}
        b = {i for i, r in enumerate(ranks) if r <= k}
        metrics[f"top{k}_retained"] = len(a & b)
        metrics[f"top{k}_jaccard"] = len(a & b) / len(a | b)
    metrics["median_abs_shift_original_top20"] = statistics.median(
        abs(float(a) - float(b)) for a, b in zip(reference, ranks) if a <= 20)
    def average_bottom_ties(r):
        count = sum(x <= n for x in r)
        return [float(x) if x <= n else (count + 1 + n) / 2 for x in r]
    a, b = average_bottom_ties(reference), average_bottom_ties(ranks)
    am, bm = math.fsum(a) / n, math.fsum(b) / n
    cov = math.fsum((x - am) * (y - bm) for x, y in zip(a, b))
    var_a = math.fsum((x - am) ** 2 for x in a)
    var_b = math.fsum((y - bm) ** 2 for y in b)
    metrics["spearman_all"] = cov / math.sqrt(var_a * var_b)
    return metrics


def check_numeric(actual, expected, location, failures, tolerance=1e-11):
    if expected is None:
        equal = actual in ("", None)
    else:
        equal = abs(float(actual) - float(expected)) <= tolerance
    if not equal:
        failures.append({"location": location, "saved": actual, "independent": expected})


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    failures = []
    cfg = json.loads((PACKAGE / "protocol.json").read_text())
    summary = json.loads((PACKAGE / "results/summary.json").read_text())
    before = json.loads((PACKAGE / "provenance/protocol_identity_before_execution.json").read_text())
    hashes_match = (sha(PACKAGE / "protocol.json") == summary["protocol_sha256"] == before["protocol_sha256"]
                    and sha(PACKAGE / "src/run_validation.py") == summary["analysis_script_sha256"] == before["analysis_script_sha256"])
    if not hashes_match:
        failures.append({"location": "protocol/script identities", "error": "frozen hashes differ"})
    inputs = {}
    for rep in REPS:
        fields, rows = read_tsv(PACKAGE / "data" / f"{rep}_matrix.tsv")
        inputs[rep] = (fields[1:], {r["taxon"]: [float(r[s]) for s in fields[1:]] for r in rows})
    draws = np.load(PACKAGE / "results/bootstrap_draw_indices.npy")
    samples = inputs["estimated_reads"][0]
    groups = [s.rsplit("_", 1)[0] for s in samples]
    strata = cfg["sampling_stability"]["strata"]
    group_counts_ok = draws.shape == (5000, 15) and all(
        [sum(groups[int(i)] == g for i in row) for g in strata] == [5, 5, 5] for row in draws)
    if not group_counts_ok:
        failures.append({"location": "bootstrap draws", "error": "shape/strata mismatch"})
    rng = np.random.default_rng(cfg["sampling_stability"]["bootstrap_seed"])
    regenerated = np.hstack([rng.choice([i for i, g in enumerate(groups) if g == group], size=(5000, 5), replace=True)
                             for group in strata])
    seed_match = np.array_equal(draws, regenerated)
    if not seed_match:
        failures.append({"location": "bootstrap generator", "error": "seed regeneration mismatch"})
    selected = sorted(set([0, 1, 2, 2499, 4997, 4998, 4999] + random.Random(20261121).sample(range(5000), 17)))
    saved_design = table(PACKAGE / "results/weight_design.csv")
    grid_rows = [r for r in saved_design if r["family"] == "joint_weight_grid"]
    exact_weights = set()
    for multipliers in itertools.product((Fraction(1, 2), Fraction(1), Fraction(3, 2)), repeat=4):
        unnormalized = [Fraction(int(w * 100), 100) * m for w, m in zip(DEFAULT, multipliers)]
        exact_weights.add(tuple(v / sum(unnormalized) for v in unnormalized))
    design_weights = {
        tuple(round(float(r[f"weight_{c}"]), 14) for c in cfg["component_order"]) for r in grid_rows}
    expected_weights = {tuple(round(float(x), 14) for x in row) for row in exact_weights}
    grid_ok = len(exact_weights) == len(grid_rows) == 79 and design_weights == expected_weights
    if not grid_ok:
        failures.append({"location": "weight grid", "error": "exact79-vector design mismatch"})
    scenario_metrics = {(r["representation"], r["scenario_id"]): r for r in table(PACKAGE / "results/scenario_metrics.csv")}
    loo_metrics = {(r["representation"], r["model"], r["omitted_sample"]): r for r in table(PACKAGE / "results/loo_metrics.csv")}
    boot_metrics = {(r["representation"], r["model"], int(r["replicate"])): r for r in table(PACKAGE / "results/bootstrap_metrics.csv.gz")}
    boot_taxon = {(r["representation"], r["model"], r["taxon"]): r for r in table(PACKAGE / "results/bootstrap_taxon_summary.csv")}
    loo_taxon = {(r["representation"], r["model"], r["taxon"]): r for r in table(PACKAGE / "results/loo_taxon_summary.csv")}
    grid_taxon = {(r["representation"], r["taxon"]): r for r in table(PACKAGE / "results/weight_grid_taxon_summary.csv")}
    audit_runs = []
    for rep in REPS:
        scenario = np.load(PACKAGE / f"results/scenario_arrays_{rep}.npz")
        names = scenario["taxa"].tolist()
        matrix = [inputs[rep][1][name] for name in names]
        count_matrix = [inputs["estimated_reads"][1][name] for name in names]
        ref, _ = rank_independent(names, matrix, DEFAULT)
        design_by_id = {r["scenario_id"]: r for r in saved_design}
        scenario_checked = 0
        for index, sid in enumerate(scenario["scenario_ids"]):
            key = str(sid)
            presence, eligible = None, None
            if key in design_by_id:
                weights = tuple(float(design_by_id[key][f"weight_{c}"]) for c in cfg["component_order"])
            else:
                weights = DEFAULT
                if key.startswith("presence_reads_ge_"):
                    threshold = int(key.split("_")[-1])
                    presence = [[v >= threshold for v in row] for row in count_matrix]
                elif key.startswith("universe_prevalence_ge_"):
                    minimum = int(key.split("_")[-1])
                    eligible = [sum(v > 0 for v in row) >= minimum for row in count_matrix]
                else:
                    raise ValueError(key)
            ranks, scores = rank_independent(names, matrix, weights, presence, eligible)
            if not np.array_equal(ranks, scenario["ranks"][index]):
                failures.append({"location": f"scenario ranks {rep}/{key}", "different_taxa": int((ranks != scenario['ranks'][index]).sum())})
            if not np.array_equal(scores, scenario["scores"][index], equal_nan=True):
                failures.append({"location": f"scenario scores {rep}/{key}", "max_abs_difference": float(np.nanmax(np.abs(scores - scenario['scores'][index])))})
            actual = scenario_metrics[rep, key]
            for metric, expected in metrics_independent(ref, ranks).items():
                check_numeric(actual[metric], expected, f"scenario metric {rep}/{key}/{metric}", failures)
            scenario_checked += 1
        grid_ids = [i for i, sid in enumerate(scenario["scenario_ids"]) if str(sid) in {r["scenario_id"] for r in grid_rows}]
        grid_ranks = scenario["ranks"][grid_ids]
        for j, name in enumerate(names):
            actual = grid_taxon[rep, name]
            for metric, expected in {
                "weight_grid_top20_frequency": sum(v <= 20 for v in grid_ranks[:, j]) / len(grid_ranks),
                "weight_grid_rank_min": min(grid_ranks[:, j]),
                "weight_grid_rank_max": max(grid_ranks[:, j]),
            }.items():
                check_numeric(actual[metric], expected, f"grid taxon {rep}/{name}/{metric}", failures)
        for model, weights in MODELS.items():
            model_ref, _ = rank_independent(names, matrix, weights)
            boot = np.load(PACKAGE / f"results/bootstrap_arrays_{rep}_{model}.npz")
            loo = np.load(PACKAGE / f"results/loo_arrays_{rep}_{model}.npz")
            if not np.array_equal(boot["reference_ranks"], model_ref) or not np.array_equal(loo["reference_ranks"], model_ref):
                failures.append({"location": f"model reference {rep}/{model}", "error": "not its own model reference"})
            for b in selected:
                chosen = [int(x) for x in draws[b]]
                subsample = [[row[j] for j in chosen] for row in matrix]
                ranks, scores = rank_independent(names, subsample, weights)
                if not np.array_equal(ranks, boot["ranks"][b]):
                    failures.append({"location": f"bootstrap ranks {rep}/{model}/{b+1}", "different_taxa": int((ranks != boot['ranks'][b]).sum())})
                if not np.array_equal(scores, boot["scores"][b], equal_nan=True):
                    failures.append({"location": f"bootstrap scores {rep}/{model}/{b+1}", "max_abs_difference": float(np.nanmax(np.abs(scores - boot['scores'][b])))})
                actual = boot_metrics[rep, model, b + 1]
                for metric, expected in metrics_independent(model_ref, ranks).items():
                    check_numeric(actual[metric], expected, f"bootstrap metric {rep}/{model}/{b+1}/{metric}", failures)
                check_numeric(actual["historical_default_top20_retained"], int(((ref <= 20) & (ranks <= 20)).sum()),
                              f"bootstrap historical-default retention {rep}/{model}/{b+1}", failures)
            for omission, name in enumerate(samples):
                subsample = [[v for i, v in enumerate(row) if i != omission] for row in matrix]
                ranks, _ = rank_independent(names, subsample, weights)
                if not np.array_equal(ranks, loo["ranks"][omission]):
                    failures.append({"location": f"LOO ranks {rep}/{model}/{name}", "different_taxa": int((ranks != loo['ranks'][omission]).sum())})
                actual = loo_metrics[rep, model, name]
                for metric, expected in metrics_independent(model_ref, ranks).items():
                    check_numeric(actual[metric], expected, f"LOO metric {rep}/{model}/{name}/{metric}", failures)
            for j, name in enumerate(names):
                ranks = boot["ranks"][:, j]
                eligible = ranks <= len(names)
                values = ranks[eligible]
                actual = boot_taxon[rep, model, name]
                expected_metrics = {
                    "reference_rank": model_ref[j], "historical_default_rank": ref[j],
                    "eligibility_frequency": sum(eligible) / len(ranks),
                    "top10_frequency": sum(ranks <= 10) / len(ranks),
                    "top20_frequency": sum(ranks <= 20) / len(ranks),
                    "top50_frequency": sum(ranks <= 50) / len(ranks),
                    "rank_p025_given_eligible": quantile_independent(values, .025),
                    "rank_median_given_eligible": quantile_independent(values, .5),
                    "rank_p975_given_eligible": quantile_independent(values, .975),
                    "score_median_given_eligible": statistics.median(boot["scores"][eligible, j]) if any(eligible) else None,
                }
                for metric, expected in expected_metrics.items():
                    check_numeric(actual[metric], expected, f"bootstrap taxon {rep}/{model}/{name}/{metric}", failures)
                actual_loo = loo_taxon[rep, model, name]
                v = loo["ranks"][:, j]
                good = v <= len(names)
                for metric, expected in {
                    "eligible_omissions": sum(good), "top20_omissions": sum(v <= 20),
                    "rank_min_given_eligible": min(v[good]), "rank_max_given_eligible": max(v[good]),
                }.items():
                    check_numeric(actual_loo[metric], expected, f"LOO taxon {rep}/{model}/{name}/{metric}", failures)
            # All replicate headline metrics are recalculated from saved ranks,
            # irrespective of the subset whose numerical ranks were recomputed.
            for b, ranks in enumerate(boot["ranks"]):
                actual = boot_metrics[rep, model, b + 1]
                retained = int(((model_ref <= 20) & (ranks <= 20)).sum())
                union = int(((model_ref <= 20) | (ranks <= 20)).sum())
                check_numeric(actual["top20_retained"], retained, f"all-bootstrap retention {rep}/{model}/{b+1}", failures)
                check_numeric(actual["top20_jaccard"], retained / union, f"all-bootstrap Jaccard {rep}/{model}/{b+1}", failures)
            audit_runs.append({"representation": rep, "model": model,
                               "fixed_scenarios_recomputed": scenario_checked,
                               "LOO_rankings_recomputed": 15,
                               "bootstrap_rankings_recomputed": len(selected),
                               "bootstrap_draws_headline_metrics_checked": 5000,
                               "taxon_bootstrap_and_LOO_summaries_checked": len(names)})
    # Independently verify every headline quantile and packaged input identity.
    summary_failures = []
    for group, saved in summary["sampling_summaries"].items():
        rep, model = group.split("/")
        subset = [r for (r_rep, r_model, _), r in boot_metrics.items()
                  if r_rep == rep and r_model == model]
        for metric, values in saved.items():
            x = [float(r[metric]) for r in subset]
            expected = {"min": min(x), "p025": quantile_independent(x, .025),
                        "median": statistics.median(x),
                        "p975": quantile_independent(x, .975), "max": max(x)}
            for label, value in expected.items():
                check_numeric(values[label], value,
                              f"headline summary {group}/{metric}/{label}", summary_failures)
    for item in summary["inputs"]:
        if sha(PACKAGE / item["file"]) != item["sha256"]:
            summary_failures.append({"location": "input identity " + item["file"],
                                     "error": "SHA-256 mismatch"})
    summary_report = {
        "checks": "Independent headline summary quantiles from all 20,000 CSV rows and all four packaged input SHA-256 identities",
        "all_checks_passed": not summary_failures, "failures": summary_failures}
    (OUT / "summary_json_audit.json").write_text(json.dumps(summary_report, indent=2) + "\n")
    failures.extend(summary_failures)
    report = {
        "independent_auditor": "importance_independent_reproduction",
        "method": "Independent standard-library ranking and quantile arithmetic; no root analysis code imported or executed",
        "protocol_and_script_frozen_hashes_match": hashes_match,
        "all_5000_draws_preserve_five_libraries_per_group": group_counts_ok,
        "seed_20261120_draws_reproduced": seed_match,
        "exact_81_multiplier_products_reduce_to_79_distinct_weight_vectors": grid_ok,
        "selected_bootstrap_replicates_one_based": [i + 1 for i in selected],
        "audit_runs": audit_runs,
        "audit_totals": {"fixed_scenario_rank_and_score_vectors_independently_recomputed": 186,
                         "leave_one_out_rank_vectors_independently_recomputed": 60,
                         "bootstrap_rank_and_score_vectors_independently_recomputed": 96,
                         "bootstrap_replicate_retention_and_Jaccard_pairs_checked": 20000,
                         "bootstrap_taxon_summaries_checked": 1700,
                         "leave_one_out_taxon_summaries_checked": 1700,
                         "weight_grid_taxon_summaries_checked": 850},
        "headline_summary_quantiles_and_four_input_hashes_match": not summary_failures,
        "failures": failures,
        "all_checks_passed": len(failures) == 0,
        "interpretation_limits": [
            "Protocol is a post hoc internal robustness plan, not preregistration or independent biological validation.",
            "No-total sampling metrics correctly compare to each no-total model's own full-data reference; historical-default retention is a separately named metric.",
            "Conditional percentile rank ranges omit bootstrap draws in which a taxon is ineligible; eligibility and unconditional top20 frequencies are supplied separately.",
            "The 79 weight-grid frequencies describe those scenarios, not a sampling probability or confidence interval.",
            "Ranks are deterministic positions including secondary tie-breaks; Spearman of positions is a secondary summary and should not replace top-set metrics.",
        ],
    }
    (OUT / "new_analysis_audit.json").write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n")
    text = "Independent audit of the importance-index robustness analysis\n\n"
    text += f"All checks passed: {report['all_checks_passed']}\nFailures: {len(failures)}\n"
    text += "Checked protocol/script identities, the 81-to-79 weight grid, all 5,000 stratified draws and their seed, 186 fixed scenarios, 60 leave-one-out rankings, 96 selected bootstrap rankings, all 20,000 top-20 retention/Jaccard pairs, all taxon summaries, headline quantiles and four input hashes.\n"
    text += "Full findings are in new_analysis_audit.json. Scientific interpretation remains internal robustness only.\n"
    (OUT / "new_analysis_audit.txt").write_text(text)
    print(json.dumps({"all_checks_passed": report["all_checks_passed"], "failures": len(failures),
                      "selected_draws": len(selected), "audit_runs": audit_runs,
                      "first_failures": failures[:5]}, indent=2))
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
