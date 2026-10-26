#!/usr/bin/env python3
"""Per-read work in each giraffe stage, compared between conditions on the same reads.

Usage: python3 -I stage_work.py RESULTS_DIR CONDITION ... [--against NAME] [--out TABLE.tsv]

Reads RESULTS_DIR/<condition>/density/stage_counts.jsonl (one JSON object per read: read name, every
annotation.stage.<path>.results count and the aligner's per-read totals, written by
jobs/stage_counts_10k.sbatch). For every metric: median, mean, nearest-rank p99 and total per condition,
and, against --against (default: the first condition), how many reads have less / the same / more work,
with a two-sided exact sign test. aligner.bases is middle + tail bases aligned. Prints numbers only.
"""
import argparse
import json
import math
import statistics


def load(path):
    reads = {}
    for line in open(path):
        r = json.loads(line)
        name = r.pop("read")
        r["aligner.bases"] = r.get("aligner.middle.bases", 0) + r.get("aligner.tail.bases", 0)
        reads[name] = r
    return reads


def sign_p(b, c):
    """Two-sided exact sign test, summed in log space (exact big-integer sums are too slow for n in the
    thousands)."""
    n = b + c
    if n == 0:
        return 1.0
    logs = [math.lgamma(n + 1) - math.lgamma(i + 1) - math.lgamma(n - i + 1) - n * math.log(2)
            for i in range(min(b, c) + 1)]
    top = max(logs)
    return min(1.0, 2 * math.exp(top) * sum(math.exp(x - top) for x in logs))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("results_dir")
    parser.add_argument("conditions", nargs="+")
    parser.add_argument("--against")
    parser.add_argument("--out")
    args = parser.parse_args()
    data = {c: load(f"{args.results_dir}/{c}/density/stage_counts.jsonl") for c in args.conditions}
    base_name = args.against or args.conditions[0]
    base = data[base_name]
    metrics = sorted({m for reads in data.values() for r in reads.values() for m in r})
    columns = ["metric", "condition", "reads", "median", "mean", "p99", "total", "total_vs_" + base_name,
               "fewer_than_" + base_name, "same", "more_than_" + base_name, "sign_test_p"]
    rows = []
    for m in metrics:
        base_total = sum(r.get(m, 0) for r in base.values())
        for c in args.conditions:
            reads = data[c]
            v = sorted(r.get(m, 0) for r in reads.values())
            row = {"metric": m, "condition": c, "reads": len(v), "median": statistics.median(v),
                   "mean": f"{statistics.fmean(v):.4g}", "p99": v[math.ceil(0.99 * len(v)) - 1], "total": sum(v),
                   "total_vs_" + base_name: f"{sum(v) / base_total:.4g}" if base_total else "NA"}
            if c != base_name:
                shared = [n for n in reads if n in base]
                fewer = sum(reads[n].get(m, 0) < base[n].get(m, 0) for n in shared)
                more = sum(reads[n].get(m, 0) > base[n].get(m, 0) for n in shared)
                row.update({"fewer_than_" + base_name: fewer, "same": len(shared) - fewer - more,
                            "more_than_" + base_name: more, "sign_test_p": f"{sign_p(fewer, more):.2g}"})
            rows.append(row)
    lines = ["\t".join(columns)] + ["\t".join(str(r.get(col, "")) for col in columns) for r in rows]
    print("\n".join(lines))
    if args.out:
        with open(args.out, "w") as f:
            f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
