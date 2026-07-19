#!/usr/bin/env python3

"""
Apply the MTD Explorer compatibility patch to HAllA 0.8.20.

HAllA uses a seaborn whitegrid style before generating a heatmap.
With Matplotlib 3.5.x, the remaining grid may trigger a
MatplotlibDeprecationWarning during pcolormesh creation.

The patch adds:

    ax.grid(False)

before the HAllA heatmap is generated.

The operation is idempotent and validates the patched Python module.
"""

from __future__ import annotations

import argparse
import py_compile
import shutil
import sys
from importlib import metadata
from pathlib import Path


MARKER = (
    "# MTD Explorer compatibility patch: "
    "disable grid before HAllA heatmap"
)

PATCH_LINE = "ax.grid(False)"


def fail(message: str) -> None:
    print(f"[ERROR] {message}", file=sys.stderr)
    raise SystemExit(1)


def locate_halla_report() -> tuple[Path, str]:
    try:
        installed_version = metadata.version("halla")
    except metadata.PackageNotFoundError:
        fail("The halla Python package is not installed.")

    if installed_version != "0.8.20":
        fail(
            "This compatibility patch is intended for HAllA 0.8.20, "
            f"but version {installed_version} was found."
        )

    try:
        import halla.utils.report
    except Exception as error:
        fail(f"Could not import halla.utils.report: {error}")

    report_file = Path(halla.utils.report.__file__).resolve()

    if not report_file.is_file():
        fail(f"HAllA report.py was not found: {report_file}")

    return report_file, installed_version


def compile_module(report_file: Path) -> None:
    try:
        py_compile.compile(
            str(report_file),
            doraise=True,
        )
    except py_compile.PyCompileError as error:
        fail(f"Patched HAllA module failed syntax validation: {error}")


def validate_patch(report_file: Path) -> None:
    text = report_file.read_text(encoding="utf-8")
    lines = text.splitlines()

    marker_indexes = [
        index
        for index, line in enumerate(lines)
        if MARKER in line
    ]

    if len(marker_indexes) != 1:
        fail(
            "Expected exactly one MTD Explorer patch marker in "
            f"{report_file}, but found {len(marker_indexes)}."
        )

    marker_index = marker_indexes[0]

    nearby_lines = lines[
        marker_index + 1:
        marker_index + 8
    ]

    if not any(
        line.strip() == PATCH_LINE
        for line in nearby_lines
    ):
        fail(
            "The patch marker was found, but ax.grid(False) "
            "was not found immediately after it."
        )

    compile_module(report_file)

    print("[OK] HAllA Matplotlib patch validation passed.")
    print(f"[INFO] HAllA module: {report_file}")


def apply_patch(report_file: Path) -> None:
    original_text = report_file.read_text(encoding="utf-8")

    if MARKER in original_text:
        print("[OK] HAllA Matplotlib patch is already installed.")
        validate_patch(report_file)
        return

    lines = original_text.splitlines(keepends=True)

    candidates: list[int] = []

    for index, line in enumerate(lines):
        if "_, ax = plt.subplots(figsize=figsize)" not in line:
            continue

        following_lines = lines[
            index + 1:
            index + 30
        ]

        if any(
            "sns.heatmap(" in following_line
            for following_line in following_lines
        ):
            candidates.append(index)

    if len(candidates) != 1:
        fail(
            "Expected exactly one HAllA heatmap plotting block, "
            f"but found {len(candidates)}. No changes were made."
        )

    subplot_index = candidates[0]

    indent = lines[subplot_index][
        : len(lines[subplot_index])
        - len(lines[subplot_index].lstrip())
    ]

    nearby_start = subplot_index + 1
    nearby_end = min(subplot_index + 10, len(lines))

    existing_grid_index = None

    for index in range(nearby_start, nearby_end):
        if lines[index].strip() == PATCH_LINE:
            existing_grid_index = index
            break

    backup_file = report_file.with_name(
        report_file.name + ".mtd_original"
    )

    if not backup_file.exists():
        shutil.copy2(report_file, backup_file)
        print(f"[INFO] Original HAllA module backed up to: {backup_file}")

    marker_line = f"{indent}{MARKER}\n"

    if existing_grid_index is not None:
        lines.insert(existing_grid_index, marker_line)
        print(
            "[INFO] ax.grid(False) was already present; "
            "the MTD Explorer marker was added."
        )
    else:
        insertion = [
            "\n",
            marker_line,
            f"{indent}{PATCH_LINE}\n",
        ]

        lines[
            subplot_index + 1:
            subplot_index + 1
        ] = insertion

    report_file.write_text(
        "".join(lines),
        encoding="utf-8",
    )

    validate_patch(report_file)

    print("[OK] HAllA Matplotlib compatibility patch applied.")
    print(f"[INFO] HAllA module: {report_file}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Apply or validate the MTD Explorer HAllA "
            "Matplotlib compatibility patch."
        )
    )

    parser.add_argument(
        "--check",
        action="store_true",
        help="Validate the patch without modifying HAllA.",
    )

    args = parser.parse_args()

    report_file, installed_version = locate_halla_report()

    print(f"[INFO] HAllA version: {installed_version}")
    print(f"[INFO] HAllA report module: {report_file}")

    if args.check:
        validate_patch(report_file)
    else:
        apply_patch(report_file)


if __name__ == "__main__":
    main()
