#!/usr/bin/env python
"""Aggregate per-iteration MAG statistics into a summary table.

Usage:
    python summarize_bins.py <output_dir> [--completed N]
"""

import argparse
from pathlib import Path

import pandas as pd


def count_mags(output, iteration):
    binned_filtered = Path(output) / f"iter_{iteration}" / "binning_filtered"
    if not binned_filtered.is_dir():
        return 0
    return sum(1 for f in binned_filtered.glob("*.fa"))


def count_reads(path):
    if not path or not Path(path).exists():
        return None
    return Path(path).stat().st_size > 0


def build_summary(output, completed_iterations, report_tsv="quality_report.filtered.tsv"):
    rows = []
    for iteration in range(1, completed_iterations + 1):
        report = (
            Path(output)
            / f"iter_{iteration}"
            / "binning_filtered"
            / report_tsv
        )
        n_mags = count_mags(output, iteration)
        reads_1 = Path(output) / f"iter_{iteration}" / "bin_mapping" / "reads.1.fq.gz"
        has_reads = reads_1.exists() and reads_1.stat().st_size > 0
        if n_mags == 0:
            has_reads = True

        avg_comp = None
        avg_cont = None
        if report.exists() and report.stat().st_size > 0:
            q = pd.read_csv(report, sep="\t")
            if not q.empty and {"Completeness", "Contamination"}.issubset(q.columns):
                avg_comp = round(float(q["Completeness"].mean()), 2)
                avg_cont = round(float(q["Contamination"].mean()), 2)

        rows.append(
            {
                "iteration": iteration,
                "high_quality_mags": n_mags,
                "unmapped_reads_present": has_reads,
                "avg_completeness": avg_comp if avg_comp is not None else "",
                "avg_contamination": avg_cont if avg_cont is not None else "",
            }
        )
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(description="Summarize iterMAG results.")
    parser.add_argument("output_dir", help="iterMAG output directory.")
    parser.add_argument(
        "--completed", type=int, default=None,
        help="Number of completed iterations (default: infer from iter_* dirs).",
    )
    parser.add_argument(
        "--out", default=None, help="Summary output path (default: <output_dir>/summary.tsv)."
    )
    args = parser.parse_args()

    output = Path(args.output_dir)
    if not output.is_dir():
        parser.error(f"output_dir does not exist: {output}")

    if args.completed is None:
        iters = sorted(
            int(p.name.split("_", 1)[1])
            for p in output.glob("iter_*")
            if p.is_dir() and p.name.split("_", 1)[1].isdigit()
        )
        completed = max(iters) if iters else 0
    else:
        completed = args.completed

    df = build_summary(output, completed)
    out_path = Path(args.out) if args.out else output / "summary.tsv"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, sep="\t", index=False)
    print(f"Wrote summary to {out_path}")


if __name__ == "__main__":
    main()