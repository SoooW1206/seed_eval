#!/usr/bin/env python3
"""Speed ratios over many repeated timing jobs (seed_eval.py timing --no-exclusive --jobs N).

Usage: python3 -I timing_repeats.py TIMING_DIR JOB_ID ... [--ref mini_k31w20] [--out PREFIX]

TIMING_DIR: data_dir/results/<reads>/timing; each JOB_ID has job<ID>/runs.tsv. Only the given jobs are used, so
the exclusive-node jobs in the same directory stay apart. Within each job, every condition's load_seconds,
map_seconds and reads_per_cpu_second are the median over its rounds and are divided by the reference condition's
value in the same job. Across jobs: median, quartiles, range of these ratios, and the fraction of jobs with a
ratio above 1 (for seconds: slower than the reference; for reads_per_cpu_second: faster).
Real reads (runs.tsv with aligned_fraction): also aligned fraction, MAPQ 0/60 fractions and identity over all
runs, the MAPQ histogram pooled over each condition's runs, and placement.tsv summed over jobs.
Writes PREFIX_ratios.tsv (and PREFIX_real.tsv, PREFIX_mapq.tsv, PREFIX_placement.tsv for real reads).
Prints numbers only.
"""
import argparse
import collections
import csv
import statistics
import sys
from pathlib import Path

METRICS = ["load_seconds", "map_seconds", "reads_per_cpu_second"]
REAL_METRICS = ["aligned_fraction", "mapq0_fraction", "mapq60_fraction", "identity_mean"]
MAPQ_BINS = [(0, 0), (1, 9), (10, 29), (30, 59), (60, 60)]


def read_tsv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def spread(values):
    """median, Q1, Q3, min, max (quartiles as statistics.quantiles, inclusive method)."""
    values = sorted(values)
    q1, _, q3 = statistics.quantiles(values, n=4, method="inclusive") if len(values) > 1 else values * 3
    return statistics.median(values), q1, q3, values[0], values[-1]


def write(path, columns, rows):
    lines = ["\t".join(columns)] + ["\t".join(str(r[c]) for c in columns) for r in rows]
    print("\n".join(lines) + "\n")
    if path:
        Path(path).write_text("\n".join(lines) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("timing_dir")
    parser.add_argument("job_ids", nargs="+")
    parser.add_argument("--ref", default="mini_k31w20")
    parser.add_argument("--out", help="output prefix")
    args = parser.parse_args()
    out = lambda suffix: f"{args.out}_{suffix}.tsv" if args.out else None

    runs = {}
    for job in args.job_ids:
        path = Path(args.timing_dir) / f"job{job}" / "runs.tsv"
        if not path.is_file():
            sys.exit(f"no {path}")
        runs[job] = read_tsv(path)
    conditions = list(dict.fromkeys(r["condition"] for rows in runs.values() for r in rows))
    if args.ref not in conditions:
        sys.exit(f"reference {args.ref} not in the runs")

    # Per job: condition -> metric -> median over rounds
    medians = {}
    for job, rows in runs.items():
        by = collections.defaultdict(lambda: collections.defaultdict(list))
        for r in rows:
            for m in METRICS:
                if r[m] not in ("", "NA"):
                    by[r["condition"]][m].append(float(r[m]))
        medians[job] = {c: {m: statistics.median(v) for m, v in ms.items()} for c, ms in by.items()}
    print("jobs:", len(runs), "nodes:", " ".join(sorted({r["node"] for rows in runs.values() for r in rows})),
          "runs per condition:", {c: sum(r["condition"] == c for rows in runs.values() for r in rows) for c in conditions})
    print()

    columns = ["condition", "metric", "jobs", "ratio_median", "ratio_q1", "ratio_q3", "ratio_min", "ratio_max",
               "jobs_ratio_above_1", "value_median"]
    rows = []
    for m in METRICS:
        for c in conditions:
            ratios = [medians[j][c][m] / medians[j][args.ref][m] for j in runs
                      if m in medians[j].get(c, {}) and m in medians[j].get(args.ref, {})]
            values = [medians[j][c][m] for j in runs if m in medians[j].get(c, {})]
            if not ratios:
                continue
            med, q1, q3, lo, hi = spread(ratios)
            rows.append({"condition": c, "metric": m, "jobs": len(ratios), "ratio_median": f"{med:.4f}",
                         "ratio_q1": f"{q1:.4f}", "ratio_q3": f"{q3:.4f}", "ratio_min": f"{lo:.4f}",
                         "ratio_max": f"{hi:.4f}", "jobs_ratio_above_1": f"{sum(x > 1 for x in ratios)}/{len(ratios)}",
                         "value_median": f"{statistics.median(values):.4g}"})
    write(out("ratios"), columns, rows)

    if "aligned_fraction" not in next(iter(runs.values()))[0]:
        return
    # Real reads: statistics over every run, MAPQ histograms pooled per condition, placement summed over jobs
    columns = ["condition", "runs"] + [f"{m}_{s}" for m in REAL_METRICS for s in ("median", "min", "max")]
    rows = []
    for c in conditions:
        row = {"condition": c}
        for m in REAL_METRICS:
            values = [float(r[m]) for rows_ in runs.values() for r in rows_ if r["condition"] == c and r[m] not in ("", "NA")]
            row["runs"] = len(values)
            row[f"{m}_median"], row[f"{m}_min"], row[f"{m}_max"] = (f"{statistics.median(values):.6g}",
                                                                    f"{min(values):.6g}", f"{max(values):.6g}")
        rows.append(row)
    write(out("real"), columns, rows)

    columns = ["condition", "runs", "aligned"] + [f"mapq_{a}" if a == b else f"mapq_{a}_{b}" for a, b in MAPQ_BINS]
    rows = []
    for c in conditions:
        hist, n = collections.Counter(), 0
        for job, rows_ in runs.items():
            for r in rows_:
                if r["condition"] == c:
                    n += 1
                    for h in read_tsv(Path(args.timing_dir) / f"job{job}" / f"round{r['round']}_{c}" / "mapq_hist.tsv"):
                        hist[int(h["mapq"])] += int(h["count"])
        total = sum(hist.values())
        row = {"condition": c, "runs": n, "aligned": total}
        for a, b in MAPQ_BINS:
            key = f"mapq_{a}" if a == b else f"mapq_{a}_{b}"
            row[key] = f"{sum(v for q, v in hist.items() if a <= q <= b) / total:.5f}"
        rows.append(row)
    write(out("mapq"), columns, rows)

    numeric = ["reads", "aligned_both", "aligned_only_condition", "aligned_only_against", "different_start",
               "higher_score_condition", "higher_score_against", "equal_score"]
    sums = collections.defaultdict(lambda: collections.Counter())
    jobs = collections.Counter()
    for job in runs:
        path = Path(args.timing_dir) / f"job{job}" / "placement.tsv"
        if not path.is_file():
            print(f"no {path}", file=sys.stderr)
            continue
        for r in read_tsv(path):
            key = (r["condition"], r["against"])
            jobs[key] += 1
            sums[key].update({k: int(r[k]) for k in numeric})
    columns = ["condition", "against", "jobs"] + numeric
    write(out("placement"), columns, [{"condition": a, "against": b, "jobs": jobs[(a, b)], **sums[(a, b)]}
                                      for a, b in sums])


if __name__ == "__main__":
    main()
