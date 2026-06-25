#!/usr/bin/env python3
"""Compare two GENIE cross-section spline XML files.

This is intentionally dependency-free so it can run on a laptop or lxplus
without a GENIE environment. It treats the old spline file as the reference,
samples both splines on the union of their knot grids, and reports pointwise
and trapezoid-area differences.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import math
import sys
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path


@dataclass
class Comparison:
    key: str
    n_ref: int
    n_cand: int
    max_abs: float
    max_rel: float
    max_rel_energy: float
    area_ref: float
    area_cand: float
    area_rel: float


def parse_splines(path: Path) -> dict[str, list[tuple[float, float]]]:
    root = ET.parse(path).getroot()
    splines: dict[str, list[tuple[float, float]]] = {}

    tunes = list(root.findall(".//genie_tune"))
    if not tunes:
        tunes = [root]

    for tune in tunes:
        tune_name = tune.attrib.get("name", "<no_tune>")
        for spline in tune.findall(".//spline"):
            name = spline.attrib["name"]
            points: list[tuple[float, float]] = []
            for knot in spline.findall("knot"):
                e_node = knot.find("E")
                xsec_node = knot.find("xsec")
                if e_node is None or xsec_node is None:
                    continue
                points.append((float(e_node.text), float(xsec_node.text)))
            if len(points) < 2:
                raise ValueError(f"{path}: spline {name} has fewer than 2 knots")
            points.sort()
            splines[f"{tune_name}/{name}"] = points

    return splines


def interp(points: list[tuple[float, float]], x: float) -> float:
    xs = [p[0] for p in points]
    i = bisect.bisect_left(xs, x)
    if i <= 0:
        return points[0][1]
    if i >= len(points):
        return points[-1][1]
    x0, y0 = points[i - 1]
    x1, y1 = points[i]
    if x1 == x0:
        return y0
    t = (x - x0) / (x1 - x0)
    return y0 + t * (y1 - y0)


def trapezoid_area(points: list[tuple[float, float]], grid: list[float]) -> float:
    area = 0.0
    last_x = grid[0]
    last_y = interp(points, last_x)
    for x in grid[1:]:
        y = interp(points, x)
        area += 0.5 * (last_y + y) * (x - last_x)
        last_x, last_y = x, y
    return area


def rel_delta(a: float, b: float, denom_floor: float) -> float:
    denom = max(abs(a), abs(b), denom_floor)
    return abs(a - b) / denom


def compare_one(
    key: str,
    ref: list[tuple[float, float]],
    cand: list[tuple[float, float]],
    rel_denom_floor: float,
    ignore_point_below: float,
) -> Comparison:
    grid = sorted({x for x, _ in ref} | {x for x, _ in cand})
    max_abs = 0.0
    max_rel = 0.0
    max_rel_energy = grid[0]

    for x in grid:
        y_ref = interp(ref, x)
        y_cand = interp(cand, x)
        abs_diff = abs(y_ref - y_cand)
        if max(abs(y_ref), abs(y_cand)) < ignore_point_below:
            rel_diff = 0.0
        else:
            rel_diff = rel_delta(y_ref, y_cand, rel_denom_floor)
        if rel_diff > max_rel:
            max_rel = rel_diff
            max_rel_energy = x
        max_abs = max(max_abs, abs_diff)

    area_ref = trapezoid_area(ref, grid)
    area_cand = trapezoid_area(cand, grid)
    area_rel = rel_delta(area_ref, area_cand, rel_denom_floor)
    return Comparison(
        key=key,
        n_ref=len(ref),
        n_cand=len(cand),
        max_abs=max_abs,
        max_rel=max_rel,
        max_rel_energy=max_rel_energy,
        area_ref=area_ref,
        area_cand=area_cand,
        area_rel=area_rel,
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("reference", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--top", type=int, default=25)
    parser.add_argument("--point-rel-warn", type=float, default=0.25)
    parser.add_argument("--area-rel-warn", type=float, default=0.05)
    parser.add_argument("--rel-denom-floor", type=float, default=1.0e-300)
    parser.add_argument("--ignore-point-below", type=float, default=0.0)
    parser.add_argument(
        "--key-contains",
        action="append",
        default=[],
        help="only compare spline keys containing this substring; repeatable",
    )
    parser.add_argument(
        "--sort-by",
        choices=("area_rel", "max_rel", "max_abs"),
        default="area_rel",
    )
    parser.add_argument("--csv", type=Path)
    args = parser.parse_args()
    if args.rel_denom_floor <= 0.0:
        parser.error("--rel-denom-floor must be positive")
    if args.ignore_point_below < 0.0:
        parser.error("--ignore-point-below must be non-negative")

    ref = parse_splines(args.reference)
    cand = parse_splines(args.candidate)
    if args.key_contains:
        def keep(key: str) -> bool:
            return all(fragment in key for fragment in args.key_contains)

        ref = {key: value for key, value in ref.items() if keep(key)}
        cand = {key: value for key, value in cand.items() if keep(key)}

    missing = sorted(set(ref) - set(cand))
    extra = sorted(set(cand) - set(ref))
    common = sorted(set(ref) & set(cand))
    if missing:
        print(f"missing candidate splines: {len(missing)}", file=sys.stderr)
        for key in missing[:10]:
            print(f"  {key}", file=sys.stderr)
    if extra:
        print(f"extra candidate splines: {len(extra)}", file=sys.stderr)
        for key in extra[:10]:
            print(f"  {key}", file=sys.stderr)
    if not common:
        print("no common splines to compare", file=sys.stderr)
        return 2

    rows = [
        compare_one(
            key,
            ref[key],
            cand[key],
            args.rel_denom_floor,
            args.ignore_point_below,
        )
        for key in common
    ]
    if args.sort_by == "max_rel":
        rows.sort(key=lambda row: (row.max_rel, row.area_rel), reverse=True)
    elif args.sort_by == "max_abs":
        rows.sort(key=lambda row: (row.max_abs, row.area_rel), reverse=True)
    else:
        rows.sort(key=lambda row: (row.area_rel, row.max_rel), reverse=True)

    if args.csv:
        with args.csv.open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(
                [
                    "key",
                    "n_ref",
                    "n_cand",
                    "max_abs",
                    "max_rel",
                    "max_rel_energy",
                    "area_ref",
                    "area_cand",
                    "area_rel",
                ]
            )
            for row in rows:
                writer.writerow(
                    [
                        row.key,
                        row.n_ref,
                        row.n_cand,
                        row.max_abs,
                        row.max_rel,
                        row.max_rel_energy,
                        row.area_ref,
                        row.area_cand,
                        row.area_rel,
                    ]
                )

    stdout = csv.writer(sys.stdout)
    stdout.writerow(
        [
            "key",
            "n_ref",
            "n_cand",
            "max_abs",
            "max_rel",
            "max_rel_energy",
            "area_rel",
            "area_ref",
            "area_cand",
        ]
    )
    for row in rows[: args.top]:
        stdout.writerow(
            [
                row.key,
                row.n_ref,
                row.n_cand,
                f"{row.max_abs:.6g}",
                f"{row.max_rel:.6g}",
                f"{row.max_rel_energy:.6g}",
                f"{row.area_rel:.6g}",
                f"{row.area_ref:.12g}",
                f"{row.area_cand:.12g}",
            ]
        )

    failed = bool(missing)
    failed = failed or any(
        row.max_rel > args.point_rel_warn or row.area_rel > args.area_rel_warn
        for row in rows
    )
    if failed:
        print(
            "comparison failed configured warning thresholds",
            file=sys.stderr,
        )
        return 2
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(1)
