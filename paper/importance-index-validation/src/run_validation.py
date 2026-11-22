"""Reproduce the frozen index and evaluate its internal ranking robustness.

Run from a published checkout with: python src/run_validation.py
Only NumPy and the Python standard library are required by this analysis.
The original matrices, rankings and historical R source are never modified.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import itertools
import json
import math
import platform
import sys
import time
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DEFAULT = np.array([.45, .30, .15, .10])
REPS = ("relative_percent", "estimated_reads")


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_csv(path, rows, fields=None):
    rows = list(rows)
    fields = fields or list(rows[0])
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "wt", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: ("" if v is None else v) for k, v in row.items()})


def read_matrix(path):
    with path.open(encoding="utf-8", newline="") as f:
        rows = list(csv.reader(f, delimiter="\t"))
    taxa = [r[0] for r in rows[1:]]
    x = np.array([list(map(float, r[1:])) for r in rows[1:]], dtype=float)
    if len(taxa) != len(set(taxa)) or not np.isfinite(x).all() or (x < 0).any():
        raise ValueError(f"Invalid source matrix: {path}")
    return taxa, rows[0][1:], x


def components(x, presence=None, eligible=None):
    # fsum reproduces the saved R summaries without row-summation artifacts.
    total = np.fromiter((math.fsum(row) for row in x), dtype=float, count=len(x))
    mean = total / x.shape[1]
    maximum = x.max(axis=1)
    prevalence = (x > 0 if presence is None else presence).mean(axis=1)
    active = total > 0
    if eligible is not None:
        active &= eligible
    if not active.any():
        raise ValueError("No eligible features")
    c = np.zeros((len(x), 4))
    c[:, 0] = prevalence
    for j, v in enumerate((mean, maximum, total), 1):
        logv = np.log10(v[active] + 1.)
        lo, hi = logv.min(), logv.max()
        c[active, j] = (logv - lo) / (hi - lo) if hi > lo else 0.
    return {"components": c, "active": active, "prevalence": prevalence,
            "mean": mean, "maximum": maximum, "total": total}


def ranked(state, weights, label_order):
    c = state["components"]
    # Match the historical left-to-right weighted sum and 6-decimal rounding.
    score = np.round(c[:, 0] * weights[0] + c[:, 1] * weights[1]
                     + c[:, 2] * weights[2] + c[:, 3] * weights[3], 6)
    ids = np.flatnonzero(state["active"])
    order = ids[np.lexsort((label_order[ids], -state["maximum"][ids],
                           -state["mean"][ids], -state["prevalence"][ids],
                           -score[ids]))]
    ranks = np.full(len(c), len(c) + 1, dtype=np.int16)
    ranks[order] = np.arange(1, len(order) + 1)
    score[~state["active"]] = np.nan
    return order, ranks, score


def spearman_positions(reference, test):
    def bottom_ties(r):
        r = r.astype(float).copy()
        absent = r > len(r)
        if absent.any():
            n_present = (~absent).sum()
            r[absent] = (n_present + 1 + len(r)) / 2
        return r - r.mean()
    a, b = bottom_ties(reference), bottom_ties(test)
    denominator = np.sqrt(np.dot(a, a) * np.dot(b, b))
    return float(np.dot(a, b) / denominator) if denominator else None


def compare(reference, test):
    m = {"eligible_taxa": int((test <= len(test)).sum()),
         "spearman_all": spearman_positions(reference, test)}
    for k in (10, 20, 50):
        a, b = reference <= k, test <= k
        overlap, union = int((a & b).sum()), int((a | b).sum())
        m[f"top{k}_retained"] = overlap
        m[f"top{k}_jaccard"] = overlap / union if union else None
    m["median_abs_shift_original_top20"] = float(np.median(np.abs(
        test[reference <= 20].astype(float) - reference[reference <= 20])))
    return m


def scenarios():
    baseline = [Fraction(x, 100) for x in (45, 30, 15, 10)]
    vectors = {}
    for mult in itertools.product((Fraction(1, 2), Fraction(1), Fraction(3, 2)), repeat=4):
        raw = [w * a for w, a in zip(baseline, mult)]
        normalized = tuple(w / sum(raw) for w in raw)
        vectors.setdefault(normalized, []).append([float(v) for v in mult])
    grid = []
    for i, (w, mults) in enumerate(vectors.items(), 1):
        grid.append({"scenario_id": f"grid_{i:03d}", "family": "joint_weight_grid",
                     "weights": np.array([float(v) for v in w]),
                     "multipliers": json.dumps(mults), "note": "Joint 0.5/1/1.5 multipliers; exact deduplication"})
    named = {"default": DEFAULT, "equal_weights": np.ones(4) / 4,
             "total_weight_added_to_mean": np.array([.45, .40, .15, 0]),
             "prevalence_mean_equal": np.array([.5, .5, 0, 0]),
             "prevalence_only": np.array([1., 0, 0, 0]),
             "mean_only": np.array([0., 1., 0, 0])}
    for i, name in enumerate(("prevalence", "mean", "maximum", "total")):
        w = DEFAULT.copy(); w[i] = 0; w /= w.sum()
        named[f"no_{name}"] = w
    comparisons = [{"scenario_id": key, "family": "named_weight_comparison", "weights": w,
                    "multipliers": "", "note": "Prespecified baseline or component ablation"}
                   for key, w in named.items()]
    return grid, comparisons, named


def quantiles(v):
    v = np.asarray(v, dtype=float)
    return {"min": float(v.min()), "p025": float(np.quantile(v, .025)),
            "median": float(np.median(v)), "p975": float(np.quantile(v, .975)),
            "max": float(v.max())}


def main(root=ROOT):
    start = time.monotonic()
    cfg = json.loads((root / "protocol.json").read_text())
    out = root / "results"; out.mkdir(exist_ok=True)
    sources, inputs = [], {}
    for rep in REPS:
        p = root / "data" / f"{rep}_matrix.tsv"
        names, samples, x = read_matrix(p)
        sources.append({"file": str(p.relative_to(root)), "sha256": sha(p), "bytes": p.stat().st_size})
        inputs[rep] = (names, samples, x)
    taxa, samples, counts = inputs["estimated_reads"]
    for rep, (names, columns, _) in inputs.items():
        if set(names) != set(taxa) or columns != samples:
            raise ValueError("Mismatched source labels")
        rowmap = {name: i for i, name in enumerate(names)}
        inputs[rep] = (taxa, samples, inputs[rep][2][[rowmap[t] for t in taxa]])
    n_taxa = len(taxa)
    label_order = np.empty(n_taxa, dtype=int)
    label_order[np.argsort(np.array(taxa), kind="stable")] = np.arange(n_taxa)
    groups = [s.rsplit("_", 1)[0] for s in samples]
    strata = [np.array([i for i, g in enumerate(groups) if g == group]) for group in cfg["sampling_stability"]["strata"]]
    if n_taxa != 425 or len(samples) != 15 or any(len(s) != 5 for s in strata):
        raise ValueError("Unexpected study dimensions or strata")
    write_csv(root / "data" / "sample_metadata.csv", [{"sample": s, "group": g} for s, g in zip(samples, groups)])
    grid, named_scenarios, named = scenarios()
    write_csv(out / "weight_design.csv", [{"scenario_id": s["scenario_id"], "family": s["family"],
        **{f"weight_{k}": float(w) for k, w in zip(cfg["component_order"], s["weights"])},
        "multipliers": s["multipliers"], "note": s["note"]} for s in grid + named_scenarios])
    rng = np.random.default_rng(cfg["sampling_stability"]["bootstrap_seed"])
    B = cfg["sampling_stability"]["bootstrap_replicates"]
    draws = np.hstack([rng.choice(ids, size=(B, 5), replace=True) for ids in strata])
    np.save(out / "bootstrap_draw_indices.npy", draws)
    write_csv(out / "bootstrap_draws.csv.gz", [{"replicate": i + 1, **{f"draw_{j+1:02d}": samples[k] for j, k in enumerate(row)}} for i, row in enumerate(draws)])
    all_scenarios, all_boot_metrics, all_boot_taxa, all_loo, all_loo_taxa = [], [], [], [], []
    baseline_rows, scenario_taxa, run_summaries = [], [], {}
    references, reference_states = {}, {}
    for rep in REPS:
        _, _, x = inputs[rep]
        state = components(x)
        reference_states[rep] = state
        ref_order, ref_ranks, ref_scores = ranked(state, DEFAULT, label_order)
        references[rep] = ref_ranks
        historical_path = root / "data" / f"{rep}_historical_ranking.tsv"
        with historical_path.open(encoding="utf-8", newline="") as f:
            historical = list(csv.DictReader(f, delimiter="\t"))
        hist_names = [r["taxon"] for r in historical]
        if hist_names != [taxa[i] for i in ref_order]:
            raise AssertionError(f"Historical ranking mismatch: {rep}")
        if any(ref_scores[i] != float(row["importance_score"]) for i, row in zip(ref_order, historical)):
            raise AssertionError(f"Historical score mismatch: {rep}")
        sources.append({"file": str(historical_path.relative_to(root)), "sha256": sha(historical_path), "bytes": historical_path.stat().st_size})
        for i in ref_order:
            baseline_rows.append({"representation": rep, "taxon": taxa[i], "rank": int(ref_ranks[i]),
                "score": float(ref_scores[i]), "prevalence_n": int((x[i] > 0).sum()),
                "mean": float(state["mean"][i]), "maximum": float(state["maximum"][i]),
                "total": float(state["total"][i]), **{f"component_{name}": float(state["components"][i, j]) for j, name in enumerate(cfg["component_order"])}})
        grid_ranks = []
        case_ranks, case_scores, case_names = [], [], []
        for scene in grid + named_scenarios:
            _, ranks, scores = ranked(state, scene["weights"], label_order)
            all_scenarios.append({"representation": rep, "scenario_id": scene["scenario_id"], "family": scene["family"],
                **{f"weight_{k}": float(w) for k, w in zip(cfg["component_order"], scene["weights"])},
                **compare(ref_ranks, ranks), "detail": scene["note"]})
            case_ranks.append(ranks); case_scores.append(scores); case_names.append(scene["scenario_id"])
            if scene["family"] == "joint_weight_grid": grid_ranks.append(ranks)
        for threshold in cfg["additional_sensitivity"]["presence_thresholds_estimated_reads"]:
            st = components(x, presence=counts >= threshold)
            _, ranks, scores = ranked(st, DEFAULT, label_order)
            key = f"presence_reads_ge_{threshold}"
            all_scenarios.append({"representation": rep, "scenario_id": key, "family": "presence_definition",
                **{f"weight_{k}": float(w) for k, w in zip(cfg["component_order"], DEFAULT)},
                **compare(ref_ranks, ranks), "detail": f"Presence only: estimated reads >= {threshold}; abundance values unchanged"})
            case_ranks.append(ranks); case_scores.append(scores); case_names.append(key)
        for minimum in cfg["additional_sensitivity"]["universe_minimum_observed_libraries"]:
            st = components(x, eligible=(counts > 0).sum(axis=1) >= minimum)
            _, ranks, scores = ranked(st, DEFAULT, label_order)
            key = f"universe_prevalence_ge_{minimum}"
            all_scenarios.append({"representation": rep, "scenario_id": key, "family": "taxon_universe",
                **{f"weight_{k}": float(w) for k, w in zip(cfg["component_order"], DEFAULT)},
                **compare(ref_ranks, ranks), "detail": f"Min-max over taxa detected in >= {minimum} libraries; no percentage reclosure"})
            case_ranks.append(ranks); case_scores.append(scores); case_names.append(key)
        case_ranks = np.array(case_ranks)
        np.savez_compressed(out / f"scenario_arrays_{rep}.npz", ranks=case_ranks,
                            scores=np.array(case_scores), scenario_ids=np.array(case_names), taxa=np.array(taxa))
        grid_ranks = np.array(grid_ranks)
        for i in ref_order:
            scenario_taxa.append({"representation": rep, "taxon": taxa[i], "reference_rank": int(ref_ranks[i]),
                "weight_grid_top20_frequency": float((grid_ranks[:, i] <= 20).mean()),
                "weight_grid_rank_min": int(grid_ranks[:, i].min()), "weight_grid_rank_max": int(grid_ranks[:, i].max())})
        for model in ("default", "no_total"):
            weights = named[model]
            _, model_ref, _ = ranked(state, weights, label_order)
            loo_ranks = []
            for omitted in range(len(samples)):
                keep = np.delete(np.arange(len(samples)), omitted)
                st = components(x[:, keep])
                _, ranks, scores = ranked(st, weights, label_order)
                loo_ranks.append(ranks)
                all_loo.append({"representation": rep, "model": model, "omitted_sample": samples[omitted],
                    "omitted_group": groups[omitted], **compare(model_ref, ranks),
                    "historical_default_top20_retained": int(((ref_ranks <= 20) & (ranks <= 20)).sum())})
            loo_ranks = np.array(loo_ranks)
            boot_ranks = np.zeros((B, n_taxa), dtype=np.int16)
            boot_scores = np.zeros((B, n_taxa), dtype=np.float64)
            for b, picked in enumerate(draws):
                st = components(x[:, picked])
                _, ranks, scores = ranked(st, weights, label_order)
                boot_ranks[b] = ranks; boot_scores[b] = scores
                all_boot_metrics.append({"representation": rep, "model": model, "replicate": b + 1,
                    **compare(model_ref, ranks), "historical_default_top20_retained": int(((ref_ranks <= 20) & (ranks <= 20)).sum())})
            np.savez_compressed(out / f"bootstrap_arrays_{rep}_{model}.npz", ranks=boot_ranks,
                scores=boot_scores, taxa=np.array(taxa), reference_ranks=model_ref)
            np.savez_compressed(out / f"loo_arrays_{rep}_{model}.npz", ranks=loo_ranks,
                omitted_samples=np.array(samples), taxa=np.array(taxa), reference_ranks=model_ref)
            for i in np.argsort(model_ref):
                eligible = boot_ranks[:, i] <= n_taxa
                v = boot_ranks[eligible, i]
                q = np.quantile(v, [.025, .5, .975]) if len(v) else [None] * 3
                all_boot_taxa.append({"representation": rep, "model": model, "taxon": taxa[i],
                    "reference_rank": int(model_ref[i]), "historical_default_rank": int(ref_ranks[i]),
                    "eligibility_frequency": float(eligible.mean()),
                    "top10_frequency": float((boot_ranks[:, i] <= 10).mean()),
                    "top20_frequency": float((boot_ranks[:, i] <= 20).mean()),
                    "top50_frequency": float((boot_ranks[:, i] <= 50).mean()),
                    "rank_p025_given_eligible": float(q[0]) if q[0] is not None else None,
                    "rank_median_given_eligible": float(q[1]) if q[1] is not None else None,
                    "rank_p975_given_eligible": float(q[2]) if q[2] is not None else None,
                    "score_median_given_eligible": float(np.median(boot_scores[eligible, i])) if eligible.any() else None})
                all_loo_taxa.append({"representation": rep, "model": model, "taxon": taxa[i],
                    "reference_rank": int(model_ref[i]), "historical_default_rank": int(ref_ranks[i]),
                    "eligible_omissions": int((loo_ranks[:, i] <= n_taxa).sum()),
                    "top20_omissions": int((loo_ranks[:, i] <= 20).sum()),
                    "rank_min_given_eligible": int(loo_ranks[loo_ranks[:, i] <= n_taxa, i].min()),
                    "rank_max_given_eligible": int(loo_ranks[loo_ranks[:, i] <= n_taxa, i].max())})
            model_metrics = [r for r in all_boot_metrics if r["representation"] == rep and r["model"] == model]
            run_summaries[f"{rep}/{model}"] = {key: quantiles([r[key] for r in model_metrics]) for key in
                ("top20_retained", "top20_jaccard", "spearman_all", "median_abs_shift_original_top20", "eligible_taxa")}
            print(f"Completed {rep}/{model}: {B} bootstrap draws and 15 omissions", flush=True)
    write_csv(out / "baseline_reproduction.csv", baseline_rows)
    write_csv(out / "scenario_metrics.csv", all_scenarios)
    write_csv(out / "weight_grid_taxon_summary.csv", scenario_taxa)
    write_csv(out / "loo_metrics.csv", all_loo)
    write_csv(out / "loo_taxon_summary.csv", all_loo_taxa)
    write_csv(out / "bootstrap_metrics.csv.gz", all_boot_metrics)
    write_csv(out / "bootstrap_taxon_summary.csv", all_boot_taxa)
    rep_comparison = compare(references["relative_percent"], references["estimated_reads"])
    state_relative = reference_states["relative_percent"]["components"]
    state_counts = reference_states["estimated_reads"]["components"]
    diag = {}
    for rep, state in reference_states.items():
        diag[rep] = {"max_abs_total_minus_15mean": float(np.max(np.abs(state["total"] - 15 * state["mean"]))),
            "mean_total_component_pearson": float(np.corrcoef(state["components"][:, 1], state["components"][:, 3])[0, 1])}
    summary = {"protocol_sha256": sha(root / "protocol.json"), "analysis_script_sha256": sha(Path(__file__)),
        "n_taxa": n_taxa, "n_samples": len(samples), "groups": dict(zip(cfg["sampling_stability"]["strata"], [len(s) for s in strata])),
        "historical_reproduction": "All 425 scores and exact ranks match in both representations",
        "joint_weight_multiplier_configurations": 81, "unique_joint_weight_vectors": len(grid),
        "named_weight_comparisons": len(named_scenarios), "bootstrap_replicates_per_model": B,
        "bootstrap_seed": cfg["sampling_stability"]["bootstrap_seed"], "sampling_summaries": run_summaries,
        "representation_comparison": rep_comparison, "component_dependence": diag,
        "software": {"python": platform.python_version(), "numpy": np.__version__},
        "elapsed_seconds": time.monotonic() - start, "inputs": sources}
    (out / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps({"unique_grid_vectors": len(grid), "summary": str(out / "summary.json"),
                      "seconds": round(summary["elapsed_seconds"], 2)}, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=ROOT)
    args = parser.parse_args()
    main(args.root.resolve())
