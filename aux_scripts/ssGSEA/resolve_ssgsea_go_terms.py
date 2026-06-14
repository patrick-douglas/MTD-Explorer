#!/usr/bin/env python3
"""
resolve_ssgsea_go_terms.py

Create a GO-name version of an ssGSEA GCT score matrix.

Resolution order:
  1) manual map, if provided
  2) local GO.db through Rscript, if available
  3) QuickGO direct term lookup, if enabled
  4) QuickGO search/secondary-ID lookup, if enabled
  5) built-in small replacement map
  6) unresolved: keep original GO ID as plot label

Outputs:
  - corrected GCT with the feature/Name column replaced by readable labels
  - TSV map documenting original GO ID, corrected GO ID, name, source, etc.
"""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import subprocess
import sys
import tempfile
import time
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple


BUILTIN_REPLACEMENTS = {
    # Observed in old eggNOG/GO-derived GMTs.
    # Keep this tiny and auditable. Prefer --manual-map for project-specific fixes.
    "GO:0017144": {
        "corrected_go_id": "GO:0006805",
        "go_name": "xenobiotic metabolic process",
        "ontology": "biological_process",
        "label": "obsolete: xenobiotic metabolic process [GO:0006805]",
        "source": "builtin_replaced_by",
        "is_obsolete": "TRUE",
        "replacement_note": "obsolete/secondary ID mapped to GO:0006805",
    }
}


@dataclass
class Resolution:
    original_go_id: str
    corrected_go_id: str
    go_name: str
    heatmap_label_base: str
    source: str
    ontology: str = ""
    is_obsolete: str = ""
    replacement_note: str = ""
    quickgo_name_raw: str = ""


@dataclass
class QuickGORecord:
    go_id: str
    name: str = ""
    aspect: str = ""
    is_obsolete: str = ""
    found: bool = False
    corrected_go_id: str = ""


def log(msg: str) -> None:
    print(msg, file=sys.stderr)


def yes_no(value: str) -> bool:
    return str(value).strip().lower() in {"1", "yes", "y", "true", "on"}


def valid_go_id(x: str) -> bool:
    if not x:
        return False
    x = str(x).strip()
    return len(x) == 10 and x.startswith("GO:") and x[3:].isdigit()


def clean_text(x: Optional[str]) -> str:
    if x is None:
        return ""
    return str(x).replace("\t", " ").replace("\r", " ").replace("\n", " ").strip()


def normalize_quickgo_label(name: str, is_obsolete: object) -> str:
    name = clean_text(name)
    obs = str(is_obsolete).strip().lower() == "true"

    if not name:
        return ""

    if not obs:
        return name

    # QuickGO often returns names like "obsolete cofactor binding".
    if name.startswith("obsolete: "):
        return name
    if name.startswith("obsolete "):
        return "obsolete: " + name.replace("obsolete ", "", 1)
    return "obsolete: " + name


def read_manual_map(path: Optional[str]) -> Dict[str, Resolution]:
    if not path:
        return {}

    manual_path = Path(path)

    if not manual_path.exists() or manual_path.stat().st_size == 0:
        log(f"[WARNING] Manual map was provided but is missing/empty: {manual_path}")
        return {}

    out: Dict[str, Resolution] = {}

    with manual_path.open(newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")

        if reader.fieldnames is None:
            log(f"[WARNING] Manual map has no header: {manual_path}")
            return {}

        for row in reader:
            old = clean_text(
                row.get("old_GO_ID")
                or row.get("original_GO_ID")
                or row.get("old")
                or row.get("from")
            )
            new_id = clean_text(
                row.get("new_GO_ID")
                or row.get("corrected_GO_ID")
                or row.get("new")
                or row.get("to")
                or old
            )
            new_name = clean_text(
                row.get("new_GO_name")
                or row.get("GO_name")
                or row.get("name")
                or row.get("term")
            )
            ontology = clean_text(row.get("ontology") or row.get("aspect"))
            note = clean_text(row.get("note") or row.get("replacement_note"))

            if not valid_go_id(old):
                continue

            if not new_id:
                new_id = old

            if not valid_go_id(new_id):
                new_id = old

            if not new_name:
                new_name = new_id

            if old != new_id:
                label = f"obsolete: {new_name} [{new_id}]"
                source = "manual_replaced_by"
                is_obsolete = "TRUE"
            else:
                label = new_name
                source = "manual_map"
                is_obsolete = ""

            out[old] = Resolution(
                original_go_id=old,
                corrected_go_id=new_id,
                go_name=new_name,
                heatmap_label_base=label,
                source=source,
                ontology=ontology,
                is_obsolete=is_obsolete,
                replacement_note=note,
            )

    log(f"[INFO] Manual GO replacements loaded: {len(out)}")
    return out


def run_go_db_resolver(go_ids: List[str], enabled: bool = True) -> Dict[str, Resolution]:
    if not enabled:
        log("[INFO] GO.db resolution disabled.")
        return {}

    rscript = shutil.which("Rscript")

    if not rscript:
        log("[WARNING] Rscript not found. Skipping GO.db resolution.")
        return {}

    unique_ids = sorted({x for x in go_ids if valid_go_id(x)})

    if not unique_ids:
        return {}

    r_code = r'''
args <- commandArgs(trailingOnly = TRUE)
ids_file <- args[1]
out_file <- args[2]

write_empty <- function(path) {
    out <- data.frame(
        GOID = character(),
        TERM = character(),
        ONTOLOGY = character(),
        stringsAsFactors = FALSE
    )
    write.table(out, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

ids <- unique(trimws(readLines(ids_file, warn = FALSE)))
ids <- ids[grepl("^GO:[0-9]{7}$", ids)]

if (!requireNamespace("AnnotationDbi", quietly = TRUE) ||
    !requireNamespace("GO.db", quietly = TRUE)) {
    message("[WARNING] AnnotationDbi or GO.db is not installed. Skipping GO.db.")
    write_empty(out_file)
    quit(status = 0)
}

all_keys <- AnnotationDbi::keys(GO.db::GO.db, keytype = "GOID")
valid <- intersect(ids, all_keys)

if (length(valid) == 0) {
    write_empty(out_file)
    quit(status = 0)
}

res <- suppressMessages(
    AnnotationDbi::select(
        GO.db::GO.db,
        keys = valid,
        keytype = "GOID",
        columns = c("TERM", "ONTOLOGY")
    )
)

res <- res[!is.na(res$GOID), , drop = FALSE]
res <- res[!duplicated(res$GOID), , drop = FALSE]

write.table(res, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
'''

    with tempfile.TemporaryDirectory(prefix="mtd_go_db_") as tmpdir:
        tmpdir_path = Path(tmpdir)
        ids_file = tmpdir_path / "go_ids.txt"
        r_file = tmpdir_path / "resolve_go_db.R"
        out_file = tmpdir_path / "go_db_map.tsv"

        ids_file.write_text("\n".join(unique_ids) + "\n")
        r_file.write_text(r_code)

        cmd = [rscript, str(r_file), str(ids_file), str(out_file)]
        log("[RUN] " + " ".join(cmd))

        proc = subprocess.run(cmd, text=True, capture_output=True)

        if proc.stdout.strip():
            log(proc.stdout.strip())

        if proc.stderr.strip():
            log(proc.stderr.strip())

        if proc.returncode != 0:
            log(f"[WARNING] GO.db resolver failed with status {proc.returncode}. Skipping GO.db.")
            return {}

        if not out_file.exists() or out_file.stat().st_size == 0:
            return {}

        resolved: Dict[str, Resolution] = {}

        with out_file.open(newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")

            for row in reader:
                goid = clean_text(row.get("GOID"))
                term = clean_text(row.get("TERM"))
                ontology = clean_text(row.get("ONTOLOGY"))

                if not valid_go_id(goid) or not term:
                    continue

                resolved[goid] = Resolution(
                    original_go_id=goid,
                    corrected_go_id=goid,
                    go_name=term,
                    heatmap_label_base=term,
                    source="GO.db",
                    ontology=ontology,
                    is_obsolete="FALSE",
                )

        log(f"[INFO] GO.db resolved: {len(resolved)} / {len(unique_ids)}")
        return resolved


def chunks(items: List[str], n: int) -> Iterable[List[str]]:
    for i in range(0, len(items), n):
        yield items[i:i + n]


def read_quickgo_cache(path: Path) -> Dict[str, QuickGORecord]:
    if not path.exists() or path.stat().st_size == 0:
        return {}

    out: Dict[str, QuickGORecord] = {}

    with path.open(newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")

        for row in reader:
            goid = clean_text(row.get("GOID"))

            if not valid_go_id(goid):
                continue

            found = str(row.get("found", "")).strip().lower() == "true"

            out[goid] = QuickGORecord(
                go_id=goid,
                name=clean_text(row.get("name")),
                aspect=clean_text(row.get("aspect")),
                is_obsolete=clean_text(row.get("is_obsolete")),
                found=found,
                corrected_go_id=clean_text(row.get("corrected_GO_ID")) or goid,
            )

    return out


def write_quickgo_cache(path: Path, cache: Dict[str, QuickGORecord]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", newline="") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(["GOID", "found", "name", "aspect", "is_obsolete", "corrected_GO_ID"])

        for goid in sorted(cache):
            rec = cache[goid]
            writer.writerow([
                goid,
                "TRUE" if rec.found else "FALSE",
                rec.name,
                rec.aspect,
                rec.is_obsolete,
                rec.corrected_go_id or goid,
            ])


def quickgo_fetch_batch(batch: List[str], timeout: int) -> Dict[str, QuickGORecord]:
    encoded_ids = ",".join(urllib.parse.quote(x, safe=":") for x in batch)
    url = f"https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/{encoded_ids}"

    req = urllib.request.Request(
        url,
        headers={
            "Accept": "application/json",
            "User-Agent": "MTD-GO-resolver/1.0",
        },
    )

    with urllib.request.urlopen(req, timeout=timeout) as response:
        data = json.loads(response.read().decode("utf-8"))

    out: Dict[str, QuickGORecord] = {}

    for item in data.get("results", []):
        goid = clean_text(item.get("id"))

        if not valid_go_id(goid):
            continue

        out[goid] = QuickGORecord(
            go_id=goid,
            name=clean_text(item.get("name")),
            aspect=clean_text(item.get("aspect")),
            is_obsolete=str(item.get("isObsolete", item.get("obsolete", ""))).upper(),
            found=True,
            corrected_go_id=goid,
        )

    # Cache negative results too.
    for goid in batch:
        if goid not in out:
            out[goid] = QuickGORecord(
                go_id=goid,
                found=False,
                corrected_go_id=goid,
            )

    return out


def run_quickgo_direct(
    go_ids: List[str],
    enabled: bool,
    cache_path: Path,
    batch_size: int,
    sleep_seconds: float,
    timeout: int,
) -> Dict[str, Resolution]:
    if not enabled:
        log("[INFO] QuickGO direct resolution disabled.")
        return {}

    unique_ids = sorted({x for x in go_ids if valid_go_id(x)})
    cache = read_quickgo_cache(cache_path)
    missing = [x for x in unique_ids if x not in cache]

    if missing:
        log(f"[INFO] QuickGO cache: {len(cache)} records")
        log(f"[INFO] QuickGO IDs to query: {len(missing)}")

        for idx, batch in enumerate(chunks(missing, batch_size), start=1):
            log(f"[INFO] Querying QuickGO batch {idx}: {len(batch)} GO IDs")

            try:
                batch_records = quickgo_fetch_batch(batch, timeout=timeout)
            except Exception as e:
                log(f"[WARNING] QuickGO batch failed: {e}")
                time.sleep(2)

                try:
                    batch_records = quickgo_fetch_batch(batch, timeout=timeout)
                except Exception as e2:
                    log(f"[ERROR] QuickGO batch failed again: {e2}")
                    batch_records = {
                        goid: QuickGORecord(
                            go_id=goid,
                            found=False,
                            corrected_go_id=goid,
                        )
                        for goid in batch
                    }

            cache.update(batch_records)
            write_quickgo_cache(cache_path, cache)
            time.sleep(sleep_seconds)

    else:
        log(f"[INFO] QuickGO cache already covers requested IDs: {len(cache)} records")

    resolved: Dict[str, Resolution] = {}

    for goid in unique_ids:
        rec = cache.get(goid)

        if rec is None or not rec.found or not rec.name:
            continue

        label = normalize_quickgo_label(rec.name, rec.is_obsolete)
        source = "QuickGO_obsolete" if str(rec.is_obsolete).lower() == "true" else "QuickGO"

        resolved[goid] = Resolution(
            original_go_id=goid,
            corrected_go_id=rec.corrected_go_id or goid,
            go_name=label,
            heatmap_label_base=label,
            source=source,
            ontology=rec.aspect,
            is_obsolete=rec.is_obsolete,
            quickgo_name_raw=rec.name,
        )

    log(f"[INFO] QuickGO direct resolved: {len(resolved)} / {len(unique_ids)}")
    return resolved


def quickgo_search_secondary(goid: str, timeout: int) -> Optional[Resolution]:
    # Best-effort only. QuickGO search schema can evolve, so this is deliberately defensive.
    query = urllib.parse.quote(goid, safe="")
    url = f"https://www.ebi.ac.uk/QuickGO/services/ontology/go/search?query={query}&limit=25"

    req = urllib.request.Request(
        url,
        headers={
            "Accept": "application/json",
            "User-Agent": "MTD-GO-resolver/1.0",
        },
    )

    with urllib.request.urlopen(req, timeout=timeout) as response:
        data = json.loads(response.read().decode("utf-8"))

    for item in data.get("results", []):
        corrected_id = clean_text(item.get("id"))

        if not valid_go_id(corrected_id):
            continue

        secondary = (
            item.get("secondaryIds")
            or item.get("secondary_ids")
            or item.get("secondaryId")
            or []
        )

        if isinstance(secondary, str):
            secondary_ids = [secondary]
        elif isinstance(secondary, list):
            secondary_ids = [clean_text(x) for x in secondary]
        else:
            secondary_ids = []

        if goid not in secondary_ids:
            continue

        name = clean_text(item.get("name"))
        aspect = clean_text(item.get("aspect"))

        if not name:
            continue

        label = f"obsolete: {name} [{corrected_id}]"

        return Resolution(
            original_go_id=goid,
            corrected_go_id=corrected_id,
            go_name=name,
            heatmap_label_base=label,
            source="QuickGO_secondary_id",
            ontology=aspect,
            is_obsolete="TRUE",
            replacement_note=f"original ID found as secondary ID of {corrected_id}",
            quickgo_name_raw=name,
        )

    return None


def run_quickgo_search(
    go_ids: List[str],
    enabled: bool,
    sleep_seconds: float,
    timeout: int,
) -> Dict[str, Resolution]:
    if not enabled:
        log("[INFO] QuickGO secondary-ID search disabled.")
        return {}

    unique_ids = sorted({x for x in go_ids if valid_go_id(x)})
    out: Dict[str, Resolution] = {}

    for idx, goid in enumerate(unique_ids, start=1):
        if idx % 25 == 1:
            log(f"[INFO] QuickGO secondary-ID search progress: {idx}/{len(unique_ids)}")

        try:
            res = quickgo_search_secondary(goid, timeout=timeout)
        except Exception:
            res = None

        if res is not None:
            out[goid] = res

        time.sleep(sleep_seconds)

    log(f"[INFO] QuickGO secondary-ID search resolved: {len(out)} / {len(unique_ids)}")
    return out


def resolve_builtin(go_ids: List[str]) -> Dict[str, Resolution]:
    out: Dict[str, Resolution] = {}

    for goid in go_ids:
        b = BUILTIN_REPLACEMENTS.get(goid)

        if not b:
            continue

        out[goid] = Resolution(
            original_go_id=goid,
            corrected_go_id=b["corrected_go_id"],
            go_name=b["go_name"],
            heatmap_label_base=b["label"],
            source=b["source"],
            ontology=b.get("ontology", ""),
            is_obsolete=b.get("is_obsolete", ""),
            replacement_note=b.get("replacement_note", ""),
        )

    if out:
        log(f"[INFO] Built-in replacements applied: {len(out)}")

    return out


def read_gct(path: Path) -> Tuple[str, str, List[str], List[List[str]]]:
    with path.open(newline="") as f:
        header1 = f.readline().rstrip("\n")
        header2 = f.readline().rstrip("\n")

        if not header1 or not header2:
            raise SystemExit(f"[ERROR] GCT file appears too short: {path}")

        reader = csv.reader(f, delimiter="\t")

        try:
            columns = next(reader)
        except StopIteration:
            raise SystemExit(f"[ERROR] GCT file has no table header: {path}")

        rows = [row for row in reader if row and any(cell.strip() for cell in row)]

    return header1, header2, columns, rows


def find_name_column(columns: List[str]) -> int:
    candidates = ["Name", "NAME", "name", "id", "ID", "feature", "Feature"]

    for c in candidates:
        if c in columns:
            return columns.index(c)

    log(f"[WARNING] No standard GCT feature-name column found. Using first column: {columns[0]}")
    return 0


def find_description_column(columns: List[str]) -> Optional[int]:
    candidates = ["Description", "DESCRIPTION", "description", "Desc", "desc"]

    for c in candidates:
        if c in columns:
            return columns.index(c)

    return None


def make_unique_labels(records: List[Resolution]) -> List[str]:
    base_counts = Counter(r.heatmap_label_base for r in records)
    seen: Dict[str, int] = defaultdict(int)
    labels: List[str] = []

    for r in records:
        label = r.heatmap_label_base

        if base_counts[label] > 1 and valid_go_id(r.original_go_id):
            label = f"{label} [{r.original_go_id}]"

        seen[label] += 1

        if seen[label] > 1:
            label = f"{label} | duplicate {seen[label]}"

        labels.append(label)

    return labels


def write_outputs(
    output_gct: Path,
    output_map: Path,
    header1: str,
    header2: str,
    columns: List[str],
    rows: List[List[str]],
    name_col_idx: int,
    description_col_idx: Optional[int],
    resolutions_for_rows: List[Resolution],
) -> None:
    labels = make_unique_labels(resolutions_for_rows)

    output_gct.parent.mkdir(parents=True, exist_ok=True)
    output_map.parent.mkdir(parents=True, exist_ok=True)

    new_rows: List[List[str]] = []

    for row, res, label in zip(rows, resolutions_for_rows, labels):
        new_row = list(row)

        if len(new_row) <= name_col_idx:
            continue

        new_row[name_col_idx] = label

        if description_col_idx is not None and len(new_row) > description_col_idx:
            new_row[description_col_idx] = res.original_go_id

        new_rows.append(new_row)

    with output_gct.open("w", newline="") as f:
        f.write(header1 + "\n")
        f.write(header2 + "\n")
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        writer.writerows(new_rows)

    with output_map.open("w", newline="") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "original_GO_ID",
            "corrected_GO_ID",
            "GO_name",
            "heatmap_label",
            "source",
            "ontology",
            "is_obsolete",
            "replacement_note",
            "quickgo_name_raw",
        ])

        for res, label in zip(resolutions_for_rows, labels):
            writer.writerow([
                res.original_go_id,
                res.corrected_go_id,
                res.go_name,
                label,
                res.source,
                res.ontology,
                res.is_obsolete,
                res.replacement_note,
                res.quickgo_name_raw,
            ])


def summarize(resolutions: List[Resolution]) -> None:
    total = len(resolutions)
    counts = Counter(r.source for r in resolutions)
    unresolved = counts.get("unresolved", 0)
    resolved = total - unresolved

    log(f"[SUMMARY] Total rows: {total}")

    for source, count in sorted(counts.items(), key=lambda x: (-x[1], x[0])):
        log(f"[SUMMARY] {source}: {count}")

    if total:
        log("[SUMMARY] Resolved percent: %.2f%%" % (resolved / total * 100))


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Resolve GO IDs in an ssGSEA GCT file to readable GO names."
    )

    parser.add_argument("--scores", required=True, help="Input ssGSEA scores GCT.")
    parser.add_argument("--out-gct", required=True, help="Output GCT with GO names.")
    parser.add_argument("--out-map", required=True, help="Output GO resolution TSV map.")
    parser.add_argument("--manual-map", default="", help="Optional TSV manual map.")
    parser.add_argument("--quickgo", default="yes", choices=["yes", "no"], help="Use QuickGO fallback.")
    parser.add_argument(
        "--quickgo-search",
        default="yes",
        choices=["yes", "no"],
        help="Try QuickGO search for secondary IDs after direct lookup fails.",
    )
    parser.add_argument(
        "--go-db",
        default="yes",
        choices=["yes", "no"],
        help="Use local Bioconductor GO.db through Rscript if available.",
    )
    parser.add_argument(
        "--quickgo-cache",
        default="",
        help="QuickGO cache TSV. Default: beside --out-map as quickgo_cache.tsv.",
    )
    parser.add_argument("--batch-size", type=int, default=100, help="QuickGO batch size.")
    parser.add_argument("--sleep", type=float, default=0.2, help="Sleep between QuickGO requests.")
    parser.add_argument("--timeout", type=int, default=60, help="HTTP timeout in seconds.")

    args = parser.parse_args()

    scores = Path(args.scores)
    output_gct = Path(args.out_gct)
    output_map = Path(args.out_map)

    if not scores.exists() or scores.stat().st_size == 0:
        raise SystemExit(f"[ERROR] Missing input GCT: {scores}")

    quickgo_cache = (
        Path(args.quickgo_cache)
        if args.quickgo_cache
        else output_map.parent / "quickgo_cache.tsv"
    )

    header1, header2, columns, rows = read_gct(scores)
    name_col_idx = find_name_column(columns)
    description_col_idx = find_description_column(columns)

    log(f"[INFO] Input GCT: {scores}")
    log(f"[INFO] Output GCT: {output_gct}")
    log(f"[INFO] Output map: {output_map}")
    log(f"[INFO] Feature column: {columns[name_col_idx]}")
    log(f"[INFO] Rows: {len(rows)}")

    original_ids: List[str] = []

    for row in rows:
        if len(row) <= name_col_idx:
            original_ids.append("")
        else:
            original_ids.append(clean_text(row[name_col_idx]).replace("\r", ""))

    unique_go_ids = sorted({x for x in original_ids if valid_go_id(x)})
    log(f"[INFO] Unique valid GO IDs detected: {len(unique_go_ids)}")

    manual_map = read_manual_map(args.manual_map)
    go_db_map = run_go_db_resolver(unique_go_ids, enabled=yes_no(args.go_db))

    remaining_after_manual_godb = [
        x for x in unique_go_ids
        if x not in manual_map and x not in go_db_map
    ]

    quickgo_map = run_quickgo_direct(
        remaining_after_manual_godb,
        enabled=yes_no(args.quickgo),
        cache_path=quickgo_cache,
        batch_size=args.batch_size,
        sleep_seconds=args.sleep,
        timeout=args.timeout,
    )

    remaining_after_quickgo = [
        x for x in remaining_after_manual_godb
        if x not in quickgo_map
    ]

    quickgo_secondary_map = run_quickgo_search(
        remaining_after_quickgo,
        enabled=yes_no(args.quickgo) and yes_no(args.quickgo_search),
        sleep_seconds=args.sleep,
        timeout=args.timeout,
    )

    remaining_after_secondary = [
        x for x in remaining_after_quickgo
        if x not in quickgo_secondary_map
    ]

    builtin_map = resolve_builtin(remaining_after_secondary)

    resolutions_for_rows: List[Resolution] = []

    for original in original_ids:
        if not valid_go_id(original):
            # Non-GO feature names are preserved.
            resolutions_for_rows.append(
                Resolution(
                    original_go_id=original,
                    corrected_go_id=original,
                    go_name=original,
                    heatmap_label_base=original,
                    source="non_GO_feature",
                )
            )
            continue

        if original in manual_map:
            res = manual_map[original]
        elif original in go_db_map:
            res = go_db_map[original]
        elif original in quickgo_map:
            res = quickgo_map[original]
        elif original in quickgo_secondary_map:
            res = quickgo_secondary_map[original]
        elif original in builtin_map:
            res = builtin_map[original]
        else:
            res = Resolution(
                original_go_id=original,
                corrected_go_id=original,
                go_name=original,
                heatmap_label_base=original,
                source="unresolved",
            )

        resolutions_for_rows.append(res)

    write_outputs(
        output_gct=output_gct,
        output_map=output_map,
        header1=header1,
        header2=header2,
        columns=columns,
        rows=rows,
        name_col_idx=name_col_idx,
        description_col_idx=description_col_idx,
        resolutions_for_rows=resolutions_for_rows,
    )

    summarize(resolutions_for_rows)

    log(f"[OK] GO-name GCT written to: {output_gct}")
    log(f"[OK] GO resolution map written to: {output_map}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
